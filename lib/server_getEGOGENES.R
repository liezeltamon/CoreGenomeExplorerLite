################################################################################
# Get genes overlapping with contacts for functional term enrichment
### FUNCTION ###################################################################
getEGOGENES <- function(dnetij){

  genes <- unique( strsplit(x=paste(dnetij$GENE, collapse=";"), split=";")[[1]] )
  genes <- genes[!is.na(genes) & genes != "NA"]
  return(genes)
  
}

################################################################################

# rm(list=ls()); gc()