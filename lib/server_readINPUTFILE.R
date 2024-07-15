################################################################################
# Read input file containing bins, genes or regions of interest
### FUNCTION ###################################################################
readINPUTFILE <- function(uploadbgstatus, inputfiledatapath, res){
  
  if( is.null(uploadbgstatus) ){
    return(NULL)
  } else if(uploadbgstatus == "uploaded" & !is.null(inputfiledatapath) ){
    
    x <- read.table(file=inputfiledatapath, header=FALSE, stringsAsFactors=FALSE)
  
    col.len <- length(x[1,])
    if(col.len >= 3){
      x <- x[,1:3]
      x <- cbind(x, ceiling( x[,2] / res ), ceiling( x[,3] / res ))
      colnames(x) <- c("chr", "start", "end", "start.bin", "end.bin")
    } else if(col.len == 1){
    } else {
      x <- NULL
    }
    
    return(x)
    
  } else {
    return(NULL)
  } 
  
}
################################################################################

# rm(list=ls()); gc()

