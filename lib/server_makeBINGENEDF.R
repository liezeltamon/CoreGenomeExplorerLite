################################################################################
# Make bingene.df
### FUNCTION ###################################################################
makeBINGENEDF <- function(res, anno.TF, bin.last, chr, anno.df){
  
  bin <- 1:bin.last
  bin.end <- bin*res
  
  hits.df <- WhichOverlap(start.query=bin.end - res + 1L, 
                          end.query=bin.end, 
                          space.query=rep(x=chr, times=length(bin)),
                          start.subject=anno.df$txStart[anno.TF],
                          end.subject=anno.df$txEnd[anno.TF],
                          space.subject=anno.df$chr[anno.TF], 
                          maxgap=-1L, minoverlap=1L, type="any")
  
  # No need for this because all chr have genes
  #validate( need(nrow(hits.df) != 0, "No genes overlapping with any bin of chr.") )
  
  hits.df <- by(data=anno.df[anno.TF,][hits.df[,"subject"],],
                INDICES=bin[hits.df[,"query"]],
                FUN=function(df){
                  apply(X=df[order(df[,"gene"]),], MARGIN=2, FUN=paste, collapse=";")
                })
  hits.df <- cbind.data.frame(bin=as.numeric(names(hits.df)), do.call("rbind", hits.df), 
                              stringsAsFactors=F)
  
  rownames(hits.df) <- hits.df$bin
  
  print("bingene.df generated.")
  
  return(hits.df)
  
}

################################################################################

# rm(list=ls()); gc()