################################################################################
# Make bingene.df
### FUNCTION ###################################################################
readMARKFILE <- function(uploadstatus=upload$status, markfiledatapath=input$markfile$datapath,
                         header.bed=input$header.bed){
  
  if( is.null(uploadstatus) ){
    return(NULL)
  } else if( uploadstatus == "uploaded" & !is.null(markfiledatapath) ){
    return( read.delim(file=markfiledatapath, header=header.bed, stringsAsFactors=FALSE) )
  } else {
    return(NULL)
  }
  
}

################################################################################

# rm(list=ls()); gc()