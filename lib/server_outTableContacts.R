################################################################################
# Table of network contacts
### FUNCTION ###################################################################
outTableContacts <- function(dnetij){
  
  dnetij$bin <- as.character(dnetij$bin)
  col.nmes <- colnames(dnetij)
  col.disp <- c("bin", "gene", "accession", "partner", "Nij", "Cp")
  add.col <- c("feature", "score")
  col.disp <- c(col.disp, add.col[add.col %in% col.nmes])
  dnetij <- dnetij[,col.disp]
  return(dnetij)
  #return( datatable(dnetij, options=list(pageLength=5), rownames=F) )
  
}

################################################################################

# rm(list=ls()); gc()