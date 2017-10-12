#' Patch size distribution extractor
#'
#' This function allows you to extract the size and type of all the patches present in an image
#' @param X the path of the file to analyses
#' @param Y the path of the folder containing the named pic which unique color fitting a specific type of organism potentially present on the plot
#' @param Z a chain of character to name the *.csv and *.jpg files
#' @return A dataframe (both in R and as a *.csv in the working directory) containing the type and size of each patch present on the analyzed image, also prints a picture with the labeled patches (same resolution as original)
#' @keywords Mosaic, image analysis, patch, size
#' @export
#' @examples #working on it


IndPatchSize=function(X,Y,Z){
  require(jpeg)
  require(png)
  require(raster)
  require(foreach)

  #X='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/MAI_2016_cropped.png'
  #Y='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/Legend/'
  #Z='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/RESULT'



  list_species=dir(Y,pattern='png')
  if(length(list_species)==1){
    species=t(as.data.frame(readPNG(paste(Y,list_species,sep=''))[1,1,1:3]))
    colnames(species)=c('red','green','blue')
    rownames(species)=unlist(strsplit(list_species,'.png'))
  }else{
    species=foreach(i=1:length(list_species),.combine=rbind)%do%{
      profile=readPNG(paste(Y,list_species[i],sep=''))
      profile[1,1,1:3]
    }
    colnames(species)=c('red','green','blue')
    rownames(species)=unlist(strsplit(list_species,'.png'))}
  cat('I will be able to identify the following groups: ',paste(rownames(species),collapse=', '),fill=T)
  # analysis

  masterfile=readPNG(X)
  cat(paste('I am currenlty working on this image: ',X,sep=''),fill=T)
  step=1
  base_num=0
  quad_data=data.frame('Group'=NA,'patch_ID'=NA,'x'=NA,'y'=NA,'area_in_pix_num'=NA)
  #quad_data=foreach(ii=1:nrow(species),.combine=rbind)%do%{
  for(ii in 1:nrow(species)){
    cat(paste('I am currenlty working on this group: ',rownames(species)[ii],sep=''),fill=T)
    base_num=base_num+1
    extract=matrix(0,dim(masterfile)[1],dim(masterfile)[2])   # IS THIS USELESS ?
    adonde_coral=which(masterfile[,,1]==species[ii,1]&masterfile[,,2]==species[ii,2]&masterfile[,,3]==species[ii,3])
    if(length(adonde_coral)>1){
      extract[adonde_coral]=1    # IS THIS USELESS ?
      AA=raster(extract)
      BB=clump(AA)
      named_layer=t(matrix(BB,ncol(BB),nrow(BB)))
      remove(AA,BB)
      gc()
      named_layer=named_layer+base_num
      colony_names=unique(array(named_layer))
      colony_names=colony_names[!is.na(colony_names)]
      Height=nrow(named_layer)
      Width=ncol(named_layer)
      if(step==1){test=named_layer
      step=step+1 }else{
        test=test
        step=step+1
      }
      test[!is.na(named_layer)]=0.8
      named_layer[is.na(named_layer)]=0

      if(dir.exists('name/')==FALSE)  dir.create('name/')
      # stamping picture process
      pb <- txtProgressBar(min = 0, max = length(colony_names), style = 3)
      for(stamp in 1:length(colony_names)){
        setTxtProgressBar(pb, stamp)
        #stamping the picture
        png(file=paste('name/',colony_names[stamp],'.png',sep=''),width = 500, height = 250)
        plot(1,1,axes=F,'n',xaxt = "n" ,yaxt = "n",xlab='',ylab='')
        text(c(1,1),labels=as.character(colony_names[stamp]),cex=2)
        dev.off()
        ID=readPNG(paste('name/',colony_names[stamp],'.png',sep=''),F)
        position=which(ID[,,1]==0,arr.ind=T)
        linepos_min=min(position[,1])
        colpos_min=min(position[,2])
        linepos_max=max(position[,1])
        colpos_max=max(position[,2])
        ID=ID[(linepos_min-1):(linepos_max+1),(colpos_min-1):(colpos_max+1),]


        centroid=round(apply(which(named_layer==colony_names[stamp],arr.ind=T),2,mean))
        if(named_layer[centroid[1],centroid[2]]!=colony_names[stamp]){
          whereisit=which(named_layer==colony_names[stamp],arr.ind=T)
          centroid= whereisit[sample(dim(which(named_layer==colony_names[stamp],arr.ind=T))[1],1),]
        }

        line_max=min(c(Height,centroid[1]+dim(ID)[1]))
        col_max=min(c(Width,centroid[2]+dim(ID)[2]))

        test[(line_max-dim(ID)[1]+1):line_max,(col_max-dim(ID)[2]+1):col_max]=ID[,,1]
      }
      close(pb)
      unlink('name/', recursive = TRUE)

      species_result=as.data.frame(foreach(j=1:length(colony_names),.combine=rbind)%do%{
        pos_2D=as.numeric(round(apply(which(named_layer==colony_names[j],arr.ind=T),2,mean)))
        cbind(rownames(species)[ii],colony_names[j],pos_2D[2],pos_2D[1],length(which(array(extract[which(named_layer==colony_names[j])])==1)))
      })
      colnames(species_result)=c('Group','patch_ID','x','y','area_in_pix_num')
      rownames(species_result)=NULL
      base_num=max(colony_names)
      quad_data=rbind(quad_data,species_result)
    }

  }

  quad_data=quad_data[-1,]
  writeJPEG(test,paste(Z,'_with_ID.jpeg',sep=''))

  #mat_dist=as.matrix(dist(cbind(as.numeric(as.character(quad_data$x)),as.numeric(as.character(quad_data$y)))))
  #colnames(mat_dist)=quad_data$patch_ID
  #diag(mat_dist)=NA
  #quad_data=cbind(quad_data,mat_dist)
  write.csv(quad_data,paste(Z,'_analysis_output.csv',sep=''))
  #file.copy(from, to, overwrite = recursive, recursive = FALSE,copy.mode = TRUE, copy.date = FALSE)
  return(quad_data)
}
















