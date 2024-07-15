################################################################################
# Table of enriched function terms
### FUNCTION ###################################################################
outTableFunxTerm <- function(ego){
  
  ego <- ego[,c("ONTOLOGY", "Description", "geneID", "GeneRatio", "p.adjust")]
  ego$p.adjust <- format(ego$p.adjust, digits=4, scipen=T)
  return(ego)
  #return( datatable(ego, options=list(pageLength=5), rownames=F) )
  
}

################################################################################

# rm(list=ls()); gc()