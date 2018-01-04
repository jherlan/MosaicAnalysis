#' Distance matrix bewteen polygons
#'
#' Generates the distance matrix between a set of polygons
#' @param X a list of polygon shape object.
#' @param cores the number of CPU cores to use for the analysis, defaults to 1 (CAUTION: each core will have a specific RAM requirement, take that in account very carefully).
#' @return The distance matrix.
#' @keywords Mosaic, image analysis, patch, size, distance
#' @export
#' @examples #working on it

#X='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/MAI_2016.png'
#Y='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/Legend/'





PolyDistanceParr=function(X,cores=1){


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





