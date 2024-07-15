################################################################################
# Sync numericInput and sliderInput for gap threshold 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# source(paste0(wk.dir, "/lib/makeNetworkData.R"))
### FUNCTION ###################################################################
makeDisplayNetworkData <- function(noi.df, anno.df, nodes, edges, res, bin.last, chr.len){
  
  if( length(noi.df[,1]) > 0 ){
    
    noi <- sort(as.numeric( unique(c(noi.df$from, noi.df$to)) ), decreasing=F)
    
    # Get partner and Cp data per node
    temp <- sapply(X=noi, simplify=F, FUN=function(node){
      
      # Rows with noi
      p.TF <- noi.df$from==node | noi.df$to==node
      
      temp <- noi.df[p.TF,c("from", "to")]
      # Transpose so partner will be extracted per row 
      # of original dataframe and partner will still match 
      # order of Cp (label column)
      temp <- t(temp)
      
      p.v <- as.vector(temp)
      rm(temp)
      p.v <- p.v[p.v!=node]
      
      Cp.v <- noi.df$label[p.TF]
      
      if( length(p.v)!=length(Cp.v) ){
        stop("Lengths of partner nodes and Cp not equal.")
      } 
      
      pCp.df <- data.frame(partner=paste(x=p.v, collapse=";"), 
                           Cp=paste(x=Cp.v, collapse=";"),
                           stringsAsFactors=F)
      
      return(pCp.df)
      
    })
    
    df <- cbind.data.frame(bin=as.numeric(nodes[as.character(noi),"id"]), 
                           start=NA, end=NA,
                           # label column has the Cp
                           nodes[as.character(noi), colnames(anno.df)])
    
    df$end <- df$bin * res
    df$end[ df$bin == as.numeric(bin.last) ] <- chr.len
    df$start <- df$end - res + 1L
    
    df <- data.frame(df, 
                     do.call("rbind", temp),
                     Nij=sapply(X=noi, simplify=T, FUN=function(node){
                       return( sum(c(noi.df$from, noi.df$to)==node) )
                     }), 
                     stringsAsFactors=F)
    
    print("Contact data prepared.", quote=FALSE)
    
    ## Reorder columns
    #df <- df[, c("bin", "start", "end",
    #             "chr", "txStart", "txEnd", "gene", "accession",
    #             "partner", "Cp", "Nij") ]
    
  } else {
    df <- NULL
  }
  
  return(df)
  
}

################################################################################

# rm(list=ls()); gc()