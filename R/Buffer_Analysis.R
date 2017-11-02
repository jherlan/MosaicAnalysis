#' Patch Buffer Analysis
#'
#' This function allows you to measure what is within the buffer of a given size of a distribution of contiguous patch. WARNING: You need the Python library GDAL to use this function. Go to: \code{https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/} for more informations on how to proceed.
#' @param X wether (i) the path of the file to analyse, (ii) a matrix where each cell contains a number representing the cell type (as outputed by the function \code{FromPictoRdata}).
#' @param Y if \code{X} is the path of the file to analyse, \code{Y} is the path of the folder containing the named pic which unique color fitting a specific type of organism potentially present on the plot. If \code{Y=NA}, then \code{X} is considered being a matrix where each cell contains a number representing the cell type (as outputed by the function \code{FromPictoRdata}).
#' @param d the distance around each patch that define the buffer.
#' @param scale the length of one pixel edge on the image (default to 1).
#' @param minimum_size the minimum size of the patches to consider in the analysis (default to 1).
#' @return A dataframe containing the size and relative buffer composition of each patch considered in the analysis.
#' @keywords Mosaic, image analysis, buffer, cover, species, patch size, size distributions.
#' @export
#' @examples #working on it

#setwd('/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/')
#X='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/MAI_2016.png'
#Y='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/Legend/'

Buffer_Analysis=function(X,Y=NA,d,scale=1,minimum_size=1){

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




  if(is.na(Y)){
    X=X
    legend=data.frame("Organism type"=paste('Organism',seq(1:length(X)),sep='_'),"ID"=seq(1:length(X)))
    colnames(legend)[1]="Organism type"
  }else{X=FromPictoRdata(X,Y)
  legend=X[[2]]
  X=X[[1]]
  }


  require(rgeos)
  require(raster)
  require(foreach)

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
  save(datap,file = 'Polygons.Rdata')
  #load('/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/New_buffer_analysis/try_one.Rdata')
  datap <- subset(datap, gArea(datap,byid = T)>scale^2*minimum_size*1/length(X))

  tailles=gArea(datap,byid = T)
  names(tailles)=datap$DN



  buffer_size=d/scale
  tokeep=foreach(i=1:length(datap),.combine=c)%do%{
    to_test=datap@polygons[[i]]@Polygons[[1]]@coords
    a=length(which(to_test[,1]<buffer_size|to_test[,1]>ncol(X)-buffer_size))
    b=length(which(to_test[,2]<buffer_size|to_test[,2]>nrow(X)-buffer_size))
    if(a!=0|b!=0){
      todo=FALSE
    }else{todo=TRUE}
    todo
  }
  shp.sub <- subset(datap, tokeep)
  #plot(shp.sub)
  #plot(shp.sub)


  test=gBuffer(shp.sub,byid = T,width = buffer_size)
  #par(mfrow=c(1,1))
  #pdf(paste('Polygons_with_buffer_',buffer_size*scale,'.pdf',sep=''),height=50*nrow(X)/ncol(X),width=50*ncol(X)/nrow(X))
  #plot(datap,col='dodgerblue')
  #plot(test,border = 'red',add=T)
  #dev.off()
  #plot(datap, border="orange",col='dodgerblue')
  #plot(test, border="blue", col='red',add=TRUE,lty=1,lwd=3)

  gc()
  save(test,file = paste('Polygon_buffer_',buffer_size*scale,sep=''))
  test=gBuffer(test,byid = T,width = 0) # anti bug device
  datap=gBuffer(datap,byid = T,width = 0) #anti bug device
  kk=gIntersection(test,datap,byid = T,drop_lower_td = T)
  #plot(kk,col='red')
  #gArea(kk,byid = T)
  save(kk,file =  paste('Polygon_interaction_buffer_',buffer_size,sep=''))


  splitted_names=t(matrix(unlist(strsplit(names(gArea(kk,byid = T)),' ')),nrow=2))
  broken_res=split(gArea(kk,byid = T)[-which(splitted_names[,1]==splitted_names[,2])],splitted_names[-which(splitted_names[,1]==splitted_names[,2]),1])

  size_buffer=gArea(test,byid = T)-gArea(shp.sub,byid=T)

  species_names=suppressWarnings(sapply(strsplit(datap$DN,'_'),function(x) paste(x[-which(!is.na(as.numeric(x)==TRUE))],collapse = '_')))
  colony_names=datap$DN

  buffer_analysis=foreach(p=1:length(colony_names),.combine=rbind)%do%{


    final_data=matrix(0,1,length(unique(species_names)))
    j=which(names(broken_res)==rownames(data.frame(datap))[p])
    if(length(j)){
      colnames(final_data)=unique(species_names)
      pure_data=broken_res[[j]]/size_buffer[which(names(size_buffer)==names(broken_res)[j])]
      for(jj in 1:length(pure_data)){
        names(pure_data)[jj]=species_names[which(rownames(data.frame(datap))==unlist(strsplit(names(pure_data)[jj],' '))[2])]
      }
      agg_f_data=aggregate(pure_data,list(names(pure_data)),sum)
      pure_data=agg_f_data[,2]
      names(pure_data)=agg_f_data[,1]
      for(jj in 1:ncol(final_data)){
        if(length(which(names(pure_data)==colnames(final_data)[jj]))!=0){
          final_data[jj]=pure_data[which(names(pure_data)==colnames(final_data)[jj])]}
      }
    }
    final_data
  }
  buffer_analysis[which(tokeep==FALSE),]=NA
  #buffer_analysis
  size_colony=gArea(datap,byid = T)
  buffer_analysis=data.frame(species_names,colony_names,size_colony,buffer_analysis)
  colnames(buffer_analysis)=paste('Ratio of buffer containing ',colnames(buffer_analysis),sep='')
  colnames(buffer_analysis)[1:3]=c('Group name','Patch ID','Patch size')
  write.csv(buffer_analysis, paste('Table_interaction_',buffer_size,'.csv',sep=''))
  la_couleur=species_names
  #la_couleur[which(tokeep==FALSE)]='out of buffer'
  plot(datap,col=factor(la_couleur))
  plot(test,border = 'dodgerblue',add=T)
  return(buffer_analysis)

}


