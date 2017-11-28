#' Patch size distribution extractor
#'
#' This function allows you to extract the size and type of all the patches present in an image using the GDAL library. WARNING: You need the Python library GDAL to use this function. Go to: \code{https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/} for more informations on how to proceed.
#' @param X the path of the file to analyses
#' @param Y the path of the folder containing the named pic which unique color fitting a specific type of organism potentially present on the plot
#' @param Z a chain of character to name the *.csv and *.jpg files (defaults to \code{X})
#' @return A dataframe (both in R and as a *.csv in the working directory) containing the type and size of each patch present on the analyzed image, also prints a picture with the labeled patches (same resolution as original)
#' @keywords Mosaic, image analysis, patch, size
#' @export
#' @examples #working on it

#X='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/MAI_2016.png'
#Y='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/Legend/'


IndPatchSizePoly=function(X,Y,Z=X){
  zepath=Z
  require(rgeos)
  require(raster)
  require(foreach)
  gdal_polygonizeR <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile',
                               pypath=NULL, readpoly=TRUE, quiet=TRUE) {
    cat('I am using gdal_polygonizer, a function written by John Baumgartner see: https://github.com/johnbaums',fill=T)
    if (isTRUE(readpoly)) require(rgdal)
    if (is.null(pypath)) {
      pypath <- Sys.which('gdal_polygonize.py')
    }
    if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system. You need to install the GDAL library, see:
                                   https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/
                                   for more infos.")
    owd <- getwd()
    on.exit(setwd(owd))
    setwd(dirname(pypath))
    if (!is.null(outshape)) {
      outshape <- sub('\\.shp$', '', outshape)
      f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
      if (any(f.exists))
        stop(sprintf('File already exists: %s',
                     toString(paste(outshape, c('shp', 'shx', 'dbf'),
                                    sep='.')[f.exists])), call.=FALSE)
    } else outshape <- tempfile()
    if (is(x, 'Raster')) {
      require(raster)
      writeRaster(x, {f <- tempfile(fileext='.tif')})
      rastpath <- normalizePath(f)
    } else if (is.character(x)) {
      rastpath <- normalizePath(x)
    } else stop('x must be a file path (character string), or a Raster object.')
    system2('python', args=(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"',
                                    pypath, rastpath, gdalformat, outshape)))
    if (isTRUE(readpoly)) {
      shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quiet)
      return(shp)
    }
    return(NULL)
  }

  X=FromPictoRdata(X,Y)
  legend=X[[2]]
  X=X[[1]]
  to_bind=unique(array(X))[which(unique(array(X))!=0)]
  foreach(i=1:length(to_bind))%do%{
    cat(paste('binding ',legend$`Organism type`[i],sep=''),fill=T)
    data1=matrix(0,nrow(X),ncol(X))
    data1[which(X==to_bind[i])]=1
    data1=clump(raster(data1,xmn=1,xmx=ncol(X),ymn=1,ymx=nrow(X)))
    gc()
    datap1=gdal_polygonizeR(data1)
    datap1@bbox=matrix(c(0,0,ncol(X),nrow(X)),2,2)
    #datap1 <- subset(datap1, gArea(datap1,byid = T)>=0.000001)
    datap1$DN=paste(legend$`Organism type`[to_bind[i]],datap1$DN,sep='_')
    if(i==1){
      datap=datap1
    }else{
      datap=rbind(datap,datap1,makeUniqueIDs = T)}
    rm(data1)
    rm(datap1)
    gc()
  }
  save(datap,file = paste(zepath,'Polygons.Rdata',sep='_'))

  pdf(paste(zepath,'label_and_shape.pdf',sep="_"),height = 15,width = 15)
  couleur=sample(c('blue','darkgreen','dodgerblue','darkorange','red','black','purple'),length(datap),replace=TRUE)
  couleur=sample(rainbow(length(datap)),length(datap))
  plot(datap,col='grey',border=couleur,bg='white')
  for(iii in 1:length(datap)){legend(coordinates(datap)[iii,1],coordinates(datap)[iii,2], legend=as.character(datap$DN)[iii],cex=.2,fill=couleur[iii],border = couleur[iii],bg = 'white')}
  dev.off()

  resultat=cbind(suppressWarnings(sapply(strsplit(datap$DN,'_'),function(x) paste(x[-which(!is.na(as.numeric(x)==TRUE))],collapse = '_'))),datap@data,gArea(datap,byid = T),coordinates(datap))
  colnames(resultat)=c('Group','patch_ID','x','y','Area_in_pix_num')
  write.csv(resultat,file = paste(zepath,'_results.csv',sep='_'),row.names = FALSE)
  return(resultat)
}


