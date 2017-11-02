#' Patch size distribution extractor
#'
#' This function allows you to extract the size and type of all the patches present in an image without using the GDAL library.
#' @param X the path of the file to analyses
#' @param Y the path of the folder containing the named pic which unique color fitting a specific type of organism potentially present on the plot
#' @param Z a chain of character to name the *.csv and *.jpg files (defaults to \code{X})
#' @return A dataframe (both in R and as a *.csv in the working directory) containing the type and size of each patch present on the analyzed image, also prints a picture with the labeled patches (same resolution as original)
#' @keywords Mosaic, image analysis, patch, size
#' @export
#' @examples #working on it

#X='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/MAI_2016.png'
#Y='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/Legend/'


IndPatchSize=function(X,Y,Z=X){
  gc()
  require(jpeg)
  require(png)
  require(raster)
  require(foreach)
  require(Matrix)
  #X='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/MAI_2016.png'
  #Y='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/Legend/'
  #Z='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/RESULT'

  cat('- Analysis started - v3.0',fill=T)
  im=FromPictoRdata(X,Y,save=FALSE)
  extract=Matrix(im[[1]])
  species=as.character(im[[2]][,1])
  present_species=unique(array(extract))
  present_species=present_species[-which(present_species==0)]
  base_num=0
  step=1
  copy_master=readPNG(X)
  quad_data=data.frame('Group'=NA,'patch_ID'=NA,'x'=NA,'y'=NA,'area_in_pix_num'=NA)
  for(ii in present_species){
    cat(paste('I am currenlty working on this group: ',species[ii],sep=''),fill=T)
    sp_mat=extract
    sp_mat[which(as.matrix(sp_mat)!=ii)]=0
    AA=raster(as.matrix(sp_mat))
    remove(sp_mat)
    gc()
    BB=clump(AA)
    remove(AA)
    gc()
    named_layer=t(matrix(BB,ncol(BB),nrow(BB)))
    remove(BB)
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
      copy_master[(line_max-dim(ID)[1]+1):line_max,(col_max-dim(ID)[2]+1):col_max,1]=ID[,,1]
      copy_master[(line_max-dim(ID)[1]+1):line_max,(col_max-dim(ID)[2]+1):col_max,2]=ID[,,2]
      copy_master[(line_max-dim(ID)[1]+1):line_max,(col_max-dim(ID)[2]+1):col_max,3]=ID[,,3]
      test[(line_max-dim(ID)[1]+1):line_max,(col_max-dim(ID)[2]+1):col_max]=ID[,,1]
    }
    close(pb)
    unlink('name/', recursive = TRUE)


    species_result=as.data.frame(foreach(j=1:length(colony_names),.combine=rbind)%do%{
      pos_2D=as.numeric(round(apply(which(named_layer==colony_names[j],arr.ind=T),2,mean)))
      cbind(species[ii],colony_names[j],pos_2D[2],pos_2D[1],length(which(named_layer==colony_names[j])))
    })


    colnames(species_result)=c('Group','patch_ID','x','y','area_in_pix_num')
    rownames(species_result)=NULL
    base_num=max(colony_names)
    quad_data=rbind(quad_data,species_result)
    remove(named_layer)
    gc()
  }



  quad_data=quad_data[-1,]
  writeJPEG(test,paste(Z,'_BW_with_ID.jpeg',sep=''))
  writeJPEG(copy_master,paste(Z,'_color_with_ID.jpeg',sep=''))
  write.csv(quad_data,paste(Z,'_analysis_output.csv',sep=''))
  gc()
  return(quad_data)
}
















