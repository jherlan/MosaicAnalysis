#' Transform a PNG into a matrix in R
#'
#' This function allows you to transform an annotated PNG file into a matrixwith every cell value corresponding to a a specific type of organism
#' @param X the path of the file to analyses
#' @param Y the path of the folder containing the named pic which unique color fitting a specific type of organism potentially present on the plot
#' @param save \code{TRUE} to save an *.Rdata file containing the obtained matrix (default to \code{FALSE})
#' @return A list containig a matrix with every cell value corresponding to a a specific type of organism present on the plot and a table providing a legend
#' @keywords Mosaic, image analysis, patch, size
#' @export
#' @examples #working on it

FromPictoRdata=function(X,Y,save=FALSE){
  gc()
  require(png)
  require(pracma)
  require(raster)
  require(foreach)
  require(Matrix)
  #X='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/MAI_2016.png'
  #Y='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/Legend/'
  #Z='/Users/yoaneynaud/Desktop/Travail/Post_doc_scripps/Mosaic/test_for_package/RESULT'

  cat('- Image Analysis - v1.0',fill=T)

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
  copy_master=masterfile
  cat(paste('I am currenlty working on this image: ',X,sep=''),fill=T)
  step=1
  base_num=0
  quad_data=data.frame('Group'=NA,'patch_ID'=NA,'x'=NA,'y'=NA,'area_in_pix_num'=NA)
  #quad_data=foreach(ii=1:nrow(species),.combine=rbind)%do%{
  #extract=simple_triplet_zero_matrix(dim(masterfile)[1],dim(masterfile)[2])
  data=Matrix(0,dim(masterfile)[1],dim(masterfile)[2],sparse=T)

  cat('I am identifying the different groups',fill=T)
  pb <- txtProgressBar(min = 0, max = nrow(species), style = 3)
  for(ii in 1:nrow(species)){
    setTxtProgressBar(pb,ii)
    base_num=base_num+1
    # IS THIS USELESS ?
    adonde_coral=which(masterfile[,,1]==species[ii,1]&masterfile[,,2]==species[ii,2]&masterfile[,,3]==species[ii,3])
    if(length(adonde_coral)>1){
      data[adonde_coral]=ii
    }}
  close(pb)
  remove(masterfile,adonde_coral)
  X=as.matrix(data)
  gc()
  Legend=as.data.frame(t(rbind(unlist(strsplit(list_species,'.png')),seq(1:length(list_species)))))
  colnames(Legend)=c('Organism type','ID')
  if(save==TRUE){save(X,file=paste(unlist(strsplit(X,'.png')),'.Rdata',sep=''))}
  resultat=list(X,Legend)
  names(resultat)=c('resulting matrix','Legend')
  return(resultat)

}


