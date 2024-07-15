################################################################################
# Sync numericInput and sliderInput for gap threshold 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# source(paste0(wk.dir, "/lib/makeNetworkData.R"))
### FUNCTION ###################################################################
makeBaseNetwork <- function(ij.df, ij.incl.ind, chr, bin.last, res, bingene.df,
                            centrobed.file, edge.res){
  
  bnet <- makeNetworkData(contact.mx=ij.df[ij.incl.ind,c("i", "j")], chr=chr,
                          cp.v=ij.df[ij.incl.ind,"Cp"], bin.last=bin.last, 
                          edge.res=edge.res, bin.len=res, mark.df=NULL, 
                          centrobed.file=centrobed.file)
  
  rownames(bnet$nodes) <- bnet$nodes$id
  
  # Append genes and transcript uniqueID on each node
  bnet$nodes <- cbind.data.frame(bnet$nodes, bingene.df[rownames(bnet$nodes),])
  # <br> html new line tag
  
  hover.text <- cbind(nodes.title=bnet$nodes$title, nodes.gene=bnet$nodes$gene)
  hover.text[is.na(hover.text)] <- ""
  bnet$nodes$title <- paste0("<p style='text-align:left;'>", 
                             "Bin: ", hover.text[,"nodes.title"], "<br>", 
                             "Gene: ", hover.text[,"nodes.gene"], "</p>")
  
  print("Base full network generated...") # REMOVE
  
  return(bnet)
  
}

################################################################################

# rm(list=ls()); gc()