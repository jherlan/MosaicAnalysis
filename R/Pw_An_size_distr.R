#' Power Analysis for Size Distribution
#'
#' This function allows you to proceed with a power analysis for the 5 most commons metrics used to describe a size distribution: mean, geometrical mean, kurtosis, skewness, coefficient of variation.
#' @param X wether (i) the path of the file to analyse, (ii) a matrix where each cell contains a number representing the cell type (as outputed by the function \code{FromPictoRdata}), (iii) a list containing the sizes of each organisms belonging to each group.
#' @param Y if \code{X} the path of the file to analyse, \code{Y} the path to the folder containing named images of unique colors corresponding to a specific type of organism potentially present on the analysed image. If \code{Y=NA}, then \code{X} is considered being a matrix where each cell contains a number representing the cell type (as outputed by the function \code{FromPictoRdata}).
#' @param Scale the size of the area covered by one pixel on the image (default to 1).
#' @param log10 if \code{TRUE} the distributions are logged in base 10 (default to \code{TRUE})
#' @return A dataframe with the results of the power analysis for the 5 most commons metric used to describe a size distribution: mean, geometrical mean, kurtosis, skewness, coefficient of variation. The 2.5 and 97.5 percentile of the estimation are returned for each metric and for 10 increasing number of points.
#' @keywords Mosaic, image analysis, power analysis, cover, species, patch size, size distributions.
#' @examples #working on it

Pw_An_size_distr=function(X,Y=NA,scale=1,log10=TRUE){
  require(raster)
  require(moments)

  gmean <- function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }

  if(mode(X)!='list'){
    if(is.na(Y)){X=X}else{
      raw=FromPictoRdata(X,Y)
      X=raw[[1]]
      legend=raw[[2]]
    }

    size_distribution=foreach(n=1:nrow(legend))%do%{
      if(n>1) rm(Xsub)
      gc()
      Xsub=X
      Xsub[which(Xsub!=legend$ID[n])]=0
      Xsub=as.matrix(clump(raster(Xsub)))
      if(log10==FALSE){as.numeric(table(Xsub))*scale}else{log10(as.numeric(table(Xsub))*scale)}
    }
  }else{
    size_distribution=X
    if(is.null(names(X))){
      legend=data.frame("Organism type"=paste('Organism',seq(1:length(X)),sep='_'),"ID"=seq(1:length(X)))
      colnames(legend)[1]="Organism type"
    }else{
      legend=data.frame("Organism type"=names(X),"ID"=seq(1:length(X)))
      colnames(legend)[1]="Organism type"
    }}



  Result=foreach(sp=1:nrow(legend))%do%{
    X=size_distribution[[sp]]
    gc()
    if(length(X)<10){'Not enough patches to proceed'}else{


      result=list()
      #mean
      cat(paste('working on the mean size of ',legend$`Organism type`[sp],sep=''),fill=T)
      res=cbind(round(seq(10,length(X),length.out = 10)),foreach(size_sample=round(seq(10,length(X),length.out = 10)),.combine=rbind)%do%{
        quantile(replicate(5000,mean(sample(X,size_sample))),probs=c(.025,.975))})
      colnames(res)[1]='number of colonies'
      res=t(as.data.frame(res))
      colnames(res)=NULL
      res
      result[[1]]=res


      #gmean
      cat(paste('working on the geometrical mean size of ',legend$`Organism type`[sp],sep=''),fill=T)
      res=cbind(round(seq(10,length(X),length.out = 10)),foreach(size_sample=round(seq(10,length(X),length.out = 10)),.combine=rbind)%do%{
        quantile(replicate(5000,gmean(sample(X,size_sample))),probs=c(.025,.975))})
      colnames(res)[1]='number of colonies'
      res=t(as.data.frame(res))
      colnames(res)=NULL
      result[[2]]=res


      #kurtosis
      cat(paste('working on the kurtosis of ',legend$`Organism type`[sp],sep=''),fill=T)
      res=cbind(round(seq(10,length(X),length.out = 10)),foreach(size_sample=round(seq(10,length(X),length.out = 10)),.combine=rbind)%do%{
        quantile(replicate(5000,kurtosis(sample(X,size_sample),na.rm=T)),probs=c(.025,.975))
        })
      colnames(res)[1]='number of colonies'
      res=t(as.data.frame(res))
      colnames(res)=NULL
      result[[3]]=res

      #skewness
      cat(paste('working on the skewness of ',legend$`Organism type`[sp],sep=''),fill=T)
      res=cbind(round(seq(10,length(X),length.out = 10)),foreach(size_sample=round(seq(10,length(X),length.out = 10)),.combine=rbind)%do%{
        quantile(replicate(5000,skewness(sample(X,size_sample),na.rm=T)),probs=c(.025,.975))
        })
      colnames(res)[1]='number of colonies'
      res=t(as.data.frame(res))
      colnames(res)=NULL
      result[[4]]=res

      #cv
      cat(paste('working on the coefficient of variation of ',legend$`Organism type`[sp],sep=''),fill=T)
      res=cbind(round(seq(10,length(X),length.out = 10)),foreach(size_sample=round(seq(10,length(X),length.out = 10)),.combine=rbind)%do%{
        quantile(replicate(5000,cv(sample(X,size_sample),na.rm=T)),probs=c(.025,.975))})
      colnames(res)[1]='number of colonies'
      res=t(as.data.frame(res))
      colnames(res)=NULL
      result[[5]]=res


      names(result)=c('mean','geometrical mean', 'kurtosis', 'skewness','coefficient of variation')
      result
    }


    }


  names(Result)=legend$`Organism type`
  names(size_distribution)=legend$`Organism type`
  R=list(size_distribution,Result)
  names(R)=c('size distribution','power analysis')
  return(R)

}





