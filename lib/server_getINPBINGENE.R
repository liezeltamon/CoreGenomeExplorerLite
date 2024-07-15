################################################################################
# Get input bin/s or gene/s
### FUNCTION ###################################################################
getINPBINGENE <- function(inpbingene, inputfile, chr){
  
  if( is.null(inputfile) ){ x <- inpbingene } 
  else {
    
    col.len <- ncol(inputfile)
    if(col.len == 5){
      inputfile <- inputfile[ inputfile[,1] == chr, ]
      inputfile <- unlist(mapply(function(start.bin,end.bin) start.bin:end.bin, start.bin=inputfile[,"start.bin"], end.bin=inputfile[,"end.bin"],
                                 SIMPLIFY=FALSE))
    } else if(col.len == 1){
      inputfile <- inputfile[,1]
    }
    x <- paste(x=inputfile, collapse=",") 
    
  } 
  
  return(x)
  
}
################################################################################

# rm(list=ls()); gc()