#' Distance matrix bewteen polygons
#'
#' Generates the distance matrix between a set of polygons
#' @param X a list of polygon shape object.
#' @param cores the number of CPU cores to use for the analysis, defaults to 1 (CAUTION: each core will have a specific RAM requirement, take that in account very carefully).
#' @param smoother a number between 1 and 100 defining the degree of polygonal smoothing to apply, with 1 meaning that 1% of the polygons node are kept, 100 meaning that 100% are kept. (a smaller number will compute faster and/or prevent bugs).
#' @return The distance matrix.
#' @keywords Mosaic, image analysis, patch, size, distance
#' @export
#' @examples #working on it

#X='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/MAI_2016.png'
#Y='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/Legend/'





PolyDistanceParr=function(X,cores=1,smoother=1){


      cat(paste('I am using here ',cores,' cores',sep=''),fill=T)
      require(doParallel)
      require(foreach)
      cores=makeCluster(cores)
      registerDoParallel(cores)
      eloign=matrix(0,length(X@polygons),length(X@polygons))
      pb <- txtProgressBar(min = 1, max = (length(X@polygons)-1), style = 3)
      for(j in 1:(length(X@polygons)-1)){
        setTxtProgressBar(pb,j)
        eloign[j,(j+1):length(X@polygons)]=foreach(i=(j+1):length(X@polygons),.combine=c,.packages='pdist')%dopar%{
          p1=X@polygons[[j]]@Polygons[[1]]@coords
            p2=X@polygons[[i]]@Polygons[[1]]@coords

          min(pdist(p1[seq(1,nrow(p1),length.out=ceiling((smoother/100)*nrow(p1))),],p2[seq(1,nrow(p2),length.out=ceiling((smoother/100)*nrow(p2))),])@dist)

          }
      }
      close(pb)
      eloign=eloign+t(eloign)
      eloign[which(eloign==0)]=NA
      stopCluster(cores)
      return(eloign)

}





