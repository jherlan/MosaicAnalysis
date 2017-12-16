#' Distance bewteen polygons
#'
#' None of what follows make sense
#' @param X wether (i) the path of the file to analyse, (ii) a matrix where each cell contains a number representing the cell type (as outputed by the function \code{FromPictoRdata}) or (iii) a SpatialPolygonDataFrame.
#' @param Y if \code{X} is the path of the file to analyse, \code{Y} is the path of the folder containing the named pic which unique color fitting a specific type of organism potentially present on the plot. If \code{Y=NA}, then \code{X} is considered being a matrix where each cell contains a number representing the cell type (as outputed by the function \code{FromPictoRdata}).
#' @param Z a chain of character to name the *.csv and *.jpg files (defaults to \code{X})
#' @return A dataframe (both in R and as a *.csv in the working directory) containing the type and size of each patch present on the analyzed image, also prints a picture with the labeled patches (same resolution as original)
#' @keywords Mosaic, image analysis, patch, size
#' @export
#' @examples #working on it

#X='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/MAI_2016.png'
#Y='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/Legend/'





PolyDistanceParr=function(X,cores=(detectCores()/2)){


      cat(paste('I am using here ',cores,' cores',sep=''),fill=T)
      require(doParallel)
      require(foreach)
      cores=makeCluster(cores)
      registerDoParallel(cores)
      eloign=matrix(0,length(X@polygons),length(X@polygons))
      pb <- txtProgressBar(min = 1, max = (length(X@polygons)-1), style = 3)
      for(j in 1:(length(X@polygons)-1)){
        setTxtProgressBar(pb,j)
        eloign[j,(j+1):length(X@polygons)]=foreach(i=(j+1):length(X@polygons),.combine=c,.packages='pdist')%dopar%{min(pdist(X@polygons[[j]]@Polygons[[1]]@coords,X@polygons[[i]]@Polygons[[1]]@coords)@dist)}
      }
      close(pb)
      eloign=eloign+t(eloign)
      eloign[which(eloign==0)]=NA
      stopCluster(cores)
      return(eloign)

}





