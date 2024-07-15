################################################################################
# Select nodes to be displayed
### FUNCTION ###################################################################
selectNodes <- function(inpbingene, inptypnet, buildtab, anno.df, anno.TF, res, bin.last){

  if(inptypnet == "genes"){
    
    g <- unique(strsplit(x=as.character(inpbingene), split=",")[[1]])
    g <- as.character(g)
    g.abs <- g[ !toupper(g) %in% anno.df$GENE ]
    validate( need(length(g.abs) == 0, 
                   paste0("Warning: ", paste(g.abs, collapse=","), 
                          " not in UCSC hg19 transcript annotation list.") ))
    
    # Identify bins the genes overlap with
    ainchr.TF <- anno.df$GENE %in% toupper(g) & anno.TF
    eval.str <- paste(ceiling(anno.df$txStart[ainchr.TF]/res), 
                      ceiling(anno.df$txEnd[ainchr.TF]/res), sep=":")
    eval.str <- paste(eval.str, collapse=",")
    eval(parse(text=paste0( 'b <- c(', eval.str, ')' )))
    
  } else if(inptypnet == "bins|regions"){
    
    # Check for unexpected characters
    nonum.inp <- gsub(x=inpbingene, pattern="[0-9]+",  replacement="")
    nonum.inp <- strsplit(x=nonum.inp, split="")[[1]]
    is.invalidchar <- any(!nonum.inp %in% c(",", ":"))
      
    if(is.invalidchar){
      b <- NULL
    } else {
      
      b <- tryCatch(
        eval(parse(text=paste0(
          'c(', inpbingene, ')'
        ))),
        error=function(e) NULL
      )
      
    }
    
    if( !is.numeric(b) ){ b <- NULL }
    
  }
  
  b <- b[ b >= 1 & b <= bin.last ]
  return( sort(unique(b), decreasing=F) )
  
}
################################################################################

# rm(list=ls()); gc()