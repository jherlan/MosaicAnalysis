#' Sorting and renaming pictures
#'
#' This function renames and sorts pictures in different folders based on the focal length provided by the EXIF file.
#' @param X The path to the folder containing the pictures. Note: The name of this folder will be used to name the newly created fodlers and in the new pictures' name.
#' @param email the email to witch updates will be sent.
#' @param erase if \code{TRUE} the original folders cointaining the pictures are erased.
#' @return
#' @keywords Mosaic, image analysis, patch, size
#' @export
#' @examples #working on it


Sorting_raw=function(X,email=NA,erase=TRUE){
  require(exif)
  require(sendmailR)





      old=dir(X)
      site=unlist(strsplit(X,'/'))[length(unlist(strsplit(X,'/')))]

      to_move=list.files(X,pattern = ".JPG$", recursive = TRUE)



      for(i in 1:length(to_move)){
        cat('\014')
        cat(paste('Working on ', site,', ', round(i / length(to_move) * 100,2), '% completed'))
        f=read_exif(paste(X,sep='/',to_move[i]))$focal_length
        new_name=paste(site,'_',f,'mm_',i,'.JPG',sep='')
        if(dir.exists(paste(X,paste(site,'_',f,'mm',sep=''),sep='/'))!=TRUE){dir.create(paste(X,paste(site,'_',f,'mm',sep=''),sep='/'))}
        file.copy(paste(X,sep='/',to_move[i]),
                  paste(X,paste(site,'_',f,'mm',sep=''),new_name,sep='/'))
        unlink(paste(X,sep='/',to_move[i]))
        #cat('Just transfered and renamed',paste(X,sep='/',to_move[i]),' to ',paste(X,paste(site,'_',f,'mm',sep=''),new_name,sep='/'),fill=T)

        if (i == length(to_move)) cat(': Done')
      }

      if(!is.na(email)){
      from <- "<ssautomaticprocessor@SSlab.com>"
      to <- paste('<',email,'>',sep='')
      subject <-lesfichiers[fi]
      body <- list(paste("I am done processing ",X,sep=''))
      sendmail(from, to, subject, body, control=list(smtpServer="ASPMX.L.GOOGLE.COM"))
      }
      if(erase){unlink(paste(X,old,sep='/'), recursive = TRUE, force = TRUE)}
}


