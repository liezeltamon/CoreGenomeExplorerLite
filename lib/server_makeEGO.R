################################################################################
# Make ego
### FUNCTION ###################################################################
makeEGO <- function(genes, bgrgenes, hugoEntrez.file){

  withProgress(message="GO enrichment analysis...",{
    ego <- funxAnnoWrapper(input=list(genes, bgrgenes), filePath=NULL, 
                           hugoEntrez.file=hugoEntrez.file)
  })
  return(ego)
  
}

################################################################################

# rm(list=ls()); gc()