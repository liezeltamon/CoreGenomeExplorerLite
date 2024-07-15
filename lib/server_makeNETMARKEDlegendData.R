################################################################################
# Make marked network legend dataframe for visLegend()
### FUNCTION ###################################################################
makeNETMARKEDlegendData <- function(MNET, colstylenet, binary.col){
  
  NETWRK <- MNET; rm(MNET)
  
  nonij.TF <- is.na(NETWRK$edges$label)
  edgewithbin.TF <- (NETWRK$edges$to - NETWRK$edges$from) > 1
  edgeToColor.TF <- nonij.TF & edgewithbin.TF
  
  fnx <- ifelse(colstylenet == "discrete", '!is.na', 'is.finite')
  col.nme <- ifelse(colstylenet == "discrete", "name.feat.olap", "score.olap.mean")
  
  eval(parse(text=paste0(
    'nonijfinite.all.TF <- c(', fnx, '(NETWRK$nodes$', col.nme, '),', 
    fnx, '(NETWRK$edges$', col.nme, ') & edgeToColor.TF)'
  )))
  
  fixed.col <- c("#0000FFB3", "#FF1493B3") # adjustcolor(col=c("blue", "deeppink"), alpha.f=0.7)
  if( sum(nonijfinite.all.TF) > 0 ){
    
    vals <- c(NETWRK$nodes[[col.nme]], NETWRK$edges[[col.nme]])[nonijfinite.all.TF]
    cols <- c(NETWRK$nodes$color, NETWRK$edges$color)[nonijfinite.all.TF]
    
    if(colstylenet == "binary"){
      
      labels.fin  <- "Overlap"
      cols.fin <- adjustcolor(col=binary.col, alpha.f=0.7)
      fixed.leg <- NULL
      
    } else if(colstylenet == "discrete"){
      
      cols.fin <- setdiff(unique(cols), fixed.col) 
      names(vals) <- cols
      labels.fin <- unname(vals[cols.fin])
      
      labels.order <- order(labels.fin)
      labels.fin <- labels.fin[labels.order] 
      cols.fin <- cols.fin[labels.order]
      
      # Fix labels
      labels.fin <- strsplit(x=labels.fin, split=";")
      labels.fin <- unlist(lapply(X=labels.fin, FUN=unique))
      
      if( length(labels.fin) != length(cols.fin) ){
        rm(labels.fin, cols.fin)
        stop("makeNETMARKEDlegendData(): Discrete colouring error.")
      }
      
      fixed.leg <- data.frame(shape="box", font.size=30, stringsAsFactors=FALSE,
                              label=c("Mult group", "No group"), color=fixed.col)
    
    } else if( grepl(x=colstylenet, pattern="gradient", fixed=TRUE) ){
      
      vals.len <- length(vals)
      #rank.toGet <- unique( c(1, seq(1, vals.len, length.out=10), vals.len) )
      rank.toGet <- unique(floor(seq(1, vals.len, length.out=10)))
      ind.ordered <- order(vals, decreasing=TRUE)[rank.toGet]
      #vals <- round(x=vals, digits=6)
      
      #labels.fin <- format(vals[ind.ordered], nsmall=3)
      labels.fin <- format(vals[ind.ordered], scientific=T, digits=3)
      cols.fin <- cols[ind.ordered]
      
      fixed.leg <- data.frame(shape="box", font.size=30, stringsAsFactors=FALSE,
                              label=c("No score"), color="#FF1493B3")
      
    }
    
    df <- data.frame(shape="box", label=labels.fin, font.size=30, color=cols.fin, stringsAsFactors=FALSE)
    df <- df[!duplicated(df), ]
    df <- rbind.data.frame(df, fixed.leg, stringsAsFactors=FALSE)
    
  } else {
    df <- data.frame(shape="box", label="No overlap", font.size=30, color="white", stringsAsFactors=FALSE)
  }

  return(df)
  
}

################################################################################

# rm(list=ls()); gc()