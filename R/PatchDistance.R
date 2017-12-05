#' Distance bewteen patch
#'
#' This function allows you to extract the distance bewteen patch as well as their sizes, excluding colonies which are too close from the edges (i.e. wether directly touching an edge or for which an edge is closer than the identified closest neighbor).
#' @param X wether (i) the path of the file to analyse, (ii) a matrix where each cell contains a number representing the cell type (as outputed by the function \code{FromPictoRdata}) or (iii) a SpatialPolygonDataFrame.
#' @param Y if \code{X} is the path of the file to analyse, \code{Y} is the path of the folder containing the named pic which unique color fitting a specific type of organism potentially present on the plot. If \code{Y=NA}, then \code{X} is considered being a matrix where each cell contains a number representing the cell type (as outputed by the function \code{FromPictoRdata}).
#' @param Z a chain of character to name the *.csv and *.jpg files (defaults to \code{X})
#' @return A dataframe (both in R and as a *.csv in the working directory) containing the type and size of each patch present on the analyzed image, also prints a picture with the labeled patches (same resolution as original)
#' @keywords Mosaic, image analysis, patch, size
#' @export
#' @examples #working on it

#X='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/MAI_2016.png'
#Y='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/Legend/'


PatchDistance=function(X,Y=NA,Z=X){
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




  if(is.na(Y)){
    if(mode(X)!="S4"){
      X=X
      legend=data.frame("Organism type"=paste('Organism',seq(1:length(X)),sep='_'),"ID"=seq(1:length(X)))
      colnames(legend)[1]="Organism type"
      to_bind=unique(array(X))[which(unique(array(X))!=0)]

      foreach(i=1:length(to_bind))%do%{
        cat(paste('binding ',legend$`Organism type`[i],sep=''),fill=T)
        data1=matrix(0,nrow(X),ncol(X))
        data1[which(X==to_bind[i])]=1
        data1=clump(raster(data1,xmn=1,xmx=ncol(X),ymn=1,ymx=nrow(X)))
        gc()
        datap1=gdal_polygonizeR(data1)
        datap1@bbox=matrix(c(0,0,ncol(X),nrow(X)),2,2)
        datap1$DN=as.factor(seq(length(datap1$DN)))
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
      X=datap
      ladistance=gDistance(X,X,byid = T)
      diag(ladistance)=NA



      tokeep=foreach(i=1:length(X),.combine=c)%do%{
        poly_neigh=ladistance[which(rownames(ladistance)==as.numeric(X@polygons[[i]]@ID)),]
        to_test=X@polygons[[i]]@Polygons[[1]]@coords
        a=length(which(to_test[,1]-min(poly_neigh,na.rm=T)<0|((X@bbox[1,2]-to_test[,1])-min(poly_neigh,na.rm=T))<0))
        b=length(which(to_test[,2]-min(poly_neigh,na.rm=T)<0|((X@bbox[2,2]-to_test[,2])-min(poly_neigh,na.rm=T))<0))
        #cat(c(i,min(poly_neigh,na.rm=T)),fill=T)
        #cat(c(i,min(to_test,na.rm=T)),fill=T)
        #cat(min(to_test,na.rm=T)==1,fill=T)
        if(a!=0|b!=0|min(to_test,na.rm=T)==1){
          todo=FALSE
        }else{todo=TRUE}
        todo
      }
      shp.sub <- subset(X, tokeep)
      plot(X,col=tokeep)
      lavraiedistance=gDistance(X,shp.sub,byid = T)
      lavraiedistance[which(lavraiedistance==0)]=NA
      Area=gArea(X,byid = T)[tokeep]
      lavraiedistance=cbind(Area,lavraiedistance)
      colnames(lavraiedistance)=c('Area',X$DN)
      rownames(lavraiedistance)=X$DN[tokeep]
      write.csv(lavraiedistance,file='Patches_distance_matrix.csv')
      return(lavraiedistance)



    }else{
      X=X
      ladistance=gDistance(X,X,byid = T)
      diag(ladistance)=NA

      tokeep=foreach(i=1:length(X),.combine=c)%do%{
        poly_neigh=ladistance[which(rownames(ladistance)==as.numeric(X@polygons[[i]]@ID)),]
        to_test=X@polygons[[i]]@Polygons[[1]]@coords
        a=length(which(to_test[,1]-min(poly_neigh,na.rm=T)<0|((X@bbox[1,2]-to_test[,1])-min(poly_neigh,na.rm=T))<0))
        b=length(which(to_test[,2]-min(poly_neigh,na.rm=T)<0|((X@bbox[2,2]-to_test[,2])-min(poly_neigh,na.rm=T))<0))
        #cat(c(i,min(poly_neigh,na.rm=T)),fill=T)
        #cat(c(i,min(to_test,na.rm=T)),fill=T)
        #cat(min(to_test,na.rm=T)==1,fill=T)
        if(a!=0|b!=0|min(to_test,na.rm=T)==1){
          todo=FALSE
        }else{todo=TRUE}
        todo
      }
      shp.sub <- subset(X, tokeep)
      plot(X,col=tokeep)
      lavraiedistance=gDistance(X,shp.sub,byid = T)
      lavraiedistance[which(lavraiedistance==0)]=NA
      Area=gArea(X,byid = T)[tokeep]
      lavraiedistance=cbind(Area,lavraiedistance)
      colnames(lavraiedistance)=c('Area',X$DN)
      rownames(lavraiedistance)=X$DN[tokeep]
      write.csv(lavraiedistance,file='Patches_distance_matrice.csv')
      return(lavraiedistance)
    }

  }else{
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
      datap1$DN=as.factor(seq(length(datap1$DN)))
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
    X=datap
    ladistance=gDistance(X,X,byid = T)
    diag(ladistance)=NA

    tokeep=foreach(i=1:length(X),.combine=c)%do%{
      poly_neigh=ladistance[which(rownames(ladistance)==as.numeric(X@polygons[[i]]@ID)),]
      to_test=X@polygons[[i]]@Polygons[[1]]@coords
      a=length(which(to_test[,1]-min(poly_neigh,na.rm=T)<0|((X@bbox[1,2]-to_test[,1])-min(poly_neigh,na.rm=T))<0))
      b=length(which(to_test[,2]-min(poly_neigh,na.rm=T)<0|((X@bbox[2,2]-to_test[,2])-min(poly_neigh,na.rm=T))<0))
      #cat(c(i,min(poly_neigh,na.rm=T)),fill=T)
      #cat(c(i,min(to_test,na.rm=T)),fill=T)
      #cat(min(to_test,na.rm=T)==1,fill=T)
      if(a!=0|b!=0|min(to_test,na.rm=T)==1){
        todo=FALSE
      }else{todo=TRUE}
      todo
    }
    shp.sub <- subset(X, tokeep)
    plot(X,col=tokeep)
    lavraiedistance=gDistance(X,shp.sub,byid = T)
    lavraiedistance[which(lavraiedistance==0)]=NA
    Area=gArea(X,byid = T)[tokeep]
    lavraiedistance=cbind(Area,lavraiedistance)
    colnames(lavraiedistance)=c('Area',X$DN)
    rownames(lavraiedistance)=X$DN[tokeep]
    write.csv(lavraiedistance,file=paste(zepath,'_distance_matrice.csv',sep=''))
    return(lavraiedistance)



  }


}







