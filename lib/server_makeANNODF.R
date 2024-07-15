################################################################################
# Table of gene annotations
### FUNCTION ###################################################################
makeANNODF <- function(anno.file){
  
  anno.df <- read.delim(file=anno.file, stringsAsFactors=F, header=T)
  anno.df <- data.frame(chr=anno.df$chrom, txStart=anno.df$txStart + 1L, txEnd=anno.df$txEnd,
                        gene=anno.df$name2, accession=anno.df$name, uniqueID=anno.df$uniqueID,
                        stringsAsFactors=F, row.names=anno.df$uniqueID)
  anno.df$GENE <- toupper(anno.df$gene)
  return(anno.df)
  
}
################################################################################

# rm(list=ls()); gc()