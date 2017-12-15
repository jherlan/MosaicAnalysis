PolyDistanceParr=function(X,Y=NULL,cores=(detectCores()/2)){
  if(is.null(Y)){
    cat(paste('I am using here ',cores,' cores',sep=''),fill=T)
    require(doMC)
    registerDoMC()
    eloign=matrix(0,length(X@polygons),length(X@polygons))
    for(j in 1:(length(X@polygons)-1)){
      eloign[j,(j+1):length(X@polygons)]=foreach(i=(j+1):length(X@polygons),.combine=c)%dopar%{min(pdist(X@polygons[[j]]@Polygons[[1]]@coords,X@polygons[[i]]@Polygons[[1]]@coords)@dist)}
    }
    eloign=eloign+t(eloign)
    eloign[which(eloign==0)]=NA
    return(eloign)
  }

}





