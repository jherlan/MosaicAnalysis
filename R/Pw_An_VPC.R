#' Power Analysis for Point Counts
#'
#' This function allows you to proceed with a power analysis of Point Counts on a annotated image.
#' @param X the path of the file to analyses
#' @param Y the path of the folder containing the named pic which unique color fitting a specific type of organism potentially present on the plot. If \code{Y=NA}, then \code{X} is considered being a matrix where each cell contains a number representing the cell type (as outputed by the function \code{FromPictoRdata}).
#' @param Cover if \code{TRUE} then power analysis is runned for the cover estimation (default to \code{TRUE})
#' @param Species if \code{TRUE} then power analysis is runned for the number of species estimation (default to \code{TRUE})
#' @param Copower the desired power for the power analysis for the cover estimation (default to \code{0.9})
#' @param Cotol the tolerance for the cover estimation (default to \code{0.1})
#' @param Sppower the desired power for the power analysis for the number of species estimation (default to \code{0.8})
#' @param Spgoal the minimum ratio of the number of species to detect (default to \code{0.5})
#' @param speed the speed of the estimation, faster is less precise, minimum is 1, default to 100.
#' @param max_n the maximum number of point allowed, default to 5000.
#' @return A dataframe with the results of the power analysis.
#' @keywords Mosaic, image analysis, power analysis, cover, species, patch
#' @export
#' @examples #working on it

#X='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/MAI_2016.png'
#Y='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/Legend/'



Pw_An_VPC=function(X,Y=NA,Cover=TRUE,Species=TRUE,Copower=.9,Cotol=.1,Sppower=.8,SpGoal=.5,speed=100,max_n=5000){

  gc()
  require(png)
  require(pracma)
  require(raster)
  require(foreach)
  require(Matrix)
  #X='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/MAI_2016.png'
  #Y='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/Legend/'
  #Z='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/RESULT'

  cat('- Power Analysis started - v1.0',fill=T)

  if(is.na(Y)){X=X}else{X=FromPictoRdata(X,Y)[[1]]}

  if(Cover==TRUE){
    cat('Calculating the # of points required to measure the cover ',fill=T)
    out=0
    power=Copower  #to replace
    tol=Cotol      # to replace
    real=length(which(X!=0))/length(X)
    real_up=real+tol
    real_down=real-tol
    np=2
    res=replicate(5000,length(which(X[sample(length(X),np)]!=0))/np)
    pwd_q=quantile(res,probs=c((1-power)/2,1-(1-power)/2))
    gradiente=10
    #start_dist=dist(rbind(c(real_down,real_up),pwd_q))
    #pb <- txtProgressBar(min = -start_dist, max = 0, style = 3)
    while(out!=1){
      if(pwd_q[1]>real_down&pwd_q[2]<real_up){out=1}else{

        np=np+gradiente
        if(np>max_n){out=1}
        pwd_qp=pwd_q
        res=replicate(5000,length(which(X[sample(length(X),np)]!=0))/np)
        pwd_q=quantile(res,probs=c((1-power)/2,1-(1-power)/2))
        gradiente=max(c(floor(speed*dist(rbind(c(real_down,real_up),pwd_q))),1))
        #setTxtProgressBar(pb, -dist(rbind(c(real_down,real_up),pwd_q)))
        cat(np,'points, ')
      }
    }
    #close(pb)
    #targeted=c(real_down,real_up)
    #names(targeted)=c('2.5%','97.5%')
    cat(np,'Done!',fill=T)
    cover=c(round(tol,2),round(power,2),real,round(as.numeric(np)))
    names(cover)=c('Estimation precision','Estimation power','Actual value','Minimum # of points required')
    RES1=as.data.frame(cover)
  }


  #########################################################################################################

  if(Species==TRUE){
    cat('Calculating the # of points required to measure the # of species ',fill=T)
    out=0
    real=length(unique(array(X)[which(array(X)!=0)]))
    goal=floor(real*SpGoal)
    np=2
    res=foreach(boot=1:5000,.combine=c)%do%{ sampled=X[sample(length(X),np)]
    length(unique(sampled[which(sampled!=0)]))
    }
    pwd_q=quantile(res,probs=c((1-Sppower)/2,1-(1-Sppower)/2))
    gradiente=10
    #pb <- txtProgressBar(min = min(pwd_q)-goal, max = 0, style = 3)
    while(out!=1){
      if(min(pwd_q)>=goal){out=1}else{

        np=np+gradiente
        if(np>max_n){out=1}
        pwd_qp=pwd_q
        res=apply(replicate(5000,X[sample(length(X),np)]),2,function(X) length(unique(X[which(X!=0)])))
        pwd_q=quantile(res,probs=c((1-Sppower)/2,1-(1-Sppower)/2))
        gradiente=max(floor(abs(min(pwd_q)-goal)*speed/100),1)
        cat(np,'points, ')
      }
    }
    #close(pb)
    cat(np,'Done!',fill=T)
    targeted=goal
    resultat=list(as.numeric(np),pwd_q)
    names(resultat)=c(paste('Number of points needed to detect a minimum of',goal,' species with a power of ',power),'obtained percentiles')
    RES2=resultat

    species=c(round(goal,2),round(power,2),real,round(as.numeric(np)))
    names(species)=c('Estimation precision','Estimation power','Actual value','Minimum # of points required')
    RES2=as.data.frame(species)
  }

  if(Species==TRUE&Cover==TRUE){
    RES=cbind(RES1,RES2)
    return(RES)
  }
  if(Species==FALSE&Cover==TRUE){

    return(RES1)
  }
  if(Species==TRUE&Cover==FALSE){

    return(RES2)
  }

}


