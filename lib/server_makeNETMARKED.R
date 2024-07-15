################################################################################
# Make marked network
### FUNCTION ###################################################################
makeNETMARKED <- function(markfile, DNET, res, chr, bed.CS, header.bed, centrobed.file,
                          colouring.style, colour.nme){
  
  NETWRK <- DNET
  
  chr.TF <- markfile[,1] == chr
  if( sum(chr.TF) > 0 ){
    
    NETWRK <- modifyNetworkData(NETWRK=NETWRK, bin.len=res, chr=chr,
                                mark.df=markfile[chr.TF,], 
                                bed.CS=bed.CS, 
                                header.bed=header.bed,
                                olap.col="#2f496e", addNodes=F,
                                centrobed.file=centrobed.file,
                                colouring.style=colouring.style,
                                colour.nme=colour.nme)
    
  } 
  
  return(NETWRK)
  
}

################################################################################

# rm(list=ls()); gc()