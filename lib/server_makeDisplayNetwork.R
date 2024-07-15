################################################################################
# Sync numericInput and sliderInput for gap threshold 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# source(paste0(wk.dir, "/lib/makeNetworkData.R"))
# source(paste0(wk.dir, "/lib/server_makeBaseNetwork.R"))
### FUNCTION ###################################################################
makeDisplayNetwork <- function(BNET=BNET(), selnode=selnode(), bin.last=bin.last(),
                               ij.df=ij.df(), ij.incl.ind=ij.incl.ind(), res=res(),
                               bingene.df=bingene.df(), chr=chr(), edge.res=input$edge.res,
                               centrobed.file){
  
  # Display whole network - lower resolution
  if( paste(selnode, collapse=",") == "NA" ){
    
    #dnet <- makeNetworkData(contact.mx=ij.df[ij.incl.ind,c("i", "j")], chr=chr,
    #                        cp.v=ij.df[ij.incl.ind,"Cp"], bin.last=bin.last, 
    #                        edge.res=edge.res, bin.len=res,
    #                        mark.df=NULL, centrobed.file=centrobed.file)
    
    #rownames(dnet$nodes) <- dnet$nodes$id
    
    ## Append genes and transcript uniqueID on each node
    #dnet$nodes <- cbind.data.frame(dnet$nodes, bingene.df[rownames(dnet$nodes),])
    ## <br> html new line tag
    #dnet$nodes$title <- paste0(dnet$nodes$title, "<br>", dnet$nodes$gene)
    
    dnet <- makeBaseNetwork(ij.df=ij.df, ij.incl.ind=ij.incl.ind, chr=chr, bin.last=bin.last,
                            res=res, bingene.df=bingene.df, centrobed.file=centrobed.file, 
                            edge.res=edge.res)
    
    # Display part of network (default edge.res=0)
  } else {
    
    edges.TF <- BNET$edges$from%in%selnode | BNET$edges$to%in%selnode
    #edges.TF <- BNET$edges$from%in%selnode & BNET$edges$to%in%selnode
    #validate( need(sum(edges.TF)>0, "edges: Input bins/genes not in defined base full network.") )
    edges <- BNET$edges[edges.TF,]
    
    nodes.TF <- BNET$nodes$id%in%c(selnode, edges$from, edges$to)
    #validate( need(sum(nodes.TF)>0, "nodes: Input bins/genes not in defined base full network.") )
    nodes <- BNET$nodes[nodes.TF,]
    
    dnet <- list(nodes=BNET$nodes[nodes.TF,], edges=edges)
    
  }
  
  print("Network for display generated...", quote=F) # REMOVE
  
  #})
  
  return(dnet)
  
}

################################################################################

# rm(list=ls()); gc()