#' A function to create fake Length/Age relationship for fish
#'
#' This function allows you create fake data, with added noise.
#' @param n Number of fish captured
#' @param maxA Age of the oldest fish captured (years)
#' @param K  Mean growth rate (mm/years)
#' @param Linf  Mean maximum size (mm)
#' @param noiseK  Standard deviation of the growth rate (mm/years)
#' @param noiseL  Standard deviation of the maximum size (mm)
#' @return A dataframe containing the simulated Age(year) and Length(mm) along with a plot presenting Length(mm) as a function of age(year).
#' @keywords Fish, Von Bertalanffy, simulation
#' @export
#' @examples
#' fake_fish_creator(12,32,0,01,10,0.5,1200)


require(jpeg)
require(png)
require(raster)
require(foreach)

# X is the path of the file to analyses
# Y is the path of the legend


  list_species=dir(Y,pattern='png')
  species=foreach(i=1:length(list_species),.combine=rbind)%do%{
    profile=readPNG(paste(Y,list_species[i],sep=''))
    profile[1,1,1:3]
  }
  colnames(species)=c('red','green','blue')
  rownames(species)=unlist(strsplit(list_species,'.png'))

  # analysis

    masterfile=readPNG(paste('to_be_analyzed/',to_analyze[i],sep=''))
    cat(paste('I am currenlty working on this mosaic: ',to_analyze[i],sep=''),fill=T)
    step=1
    quad_data=foreach(ii=1:nrow(species),.combine=rbind)%do%{
      cat(paste('I am currenlty working on this group: ',rownames(species)[ii],sep=''),fill=T)
      if(ii==1){base_num=0}else{base_num=base_num+1}
      extract=matrix(0,dim(masterfile)[1],dim(masterfile)[2])   # IS THIS USELESS ?
      adonde_coral=which(masterfile[,,1]==species[ii,1]&masterfile[,,2]==species[ii,2]&masterfile[,,3]==species[ii,3])
      if(length(adonde_coral)>1){
      extract[adonde_coral]=1    # IS THIS USELESS ?
      named_layer=as.matrix(clump(raster(extract)))
      named_layer=named_layer+base_num
      colony_names=unique(array(named_layer))
      colony_names=colony_names[!is.na(colony_names)]
      Height=nrow(named_layer)
      Width=ncol(named_layer)
      if(step==1){test=named_layer}else{
        test=test
      step=step+1
      }
      test[!is.na(named_layer)]=0.8
      named_layer[is.na(named_layer)]=0

      # stamping picture process
      for(stamp in 1:length(colony_names)){

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
      #if(length(colony_names)>1){
      species_result=as.data.frame(foreach(j=1:length(colony_names),.combine=rbind)%do%{
        pos_2D=as.numeric(round(apply(which(named_layer==colony_names[j],arr.ind=T),2,mean)))
        cbind(rownames(species)[ii],colony_names[j],pos_2D[2],pos_2D[1],length(which(array(extract[which(named_layer==colony_names[j])])==1)))
      })
      colnames(species_result)=c('Group','patch_ID','x','y',' area in # pix')
      rownames(species_result)=NULL
      base_num=max(colony_names)
      species_result
    }
    }
    writeJPEG(test,paste('analyzed/',unlist(strsplit(to_analyze[i],'.png')),'_identified.jpeg',sep=''))

    #mat_dist=as.matrix(dist(cbind(as.numeric(as.character(quad_data$x)),as.numeric(as.character(quad_data$y)))))
    #colnames(mat_dist)=quad_data$patch_ID
    #diag(mat_dist)=NA
    #quad_data=cbind(quad_data,mat_dist)
    write.csv(quad_data,paste('analyzed/',unlist(strsplit(to_analyze[i],'.png')),'_analysis.csv',sep=''))
    #file.copy(from, to, overwrite = recursive, recursive = FALSE,copy.mode = TRUE, copy.date = FALSE)

















