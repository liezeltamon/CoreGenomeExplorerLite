################################################################################
# Sync numericInput and sliderInput for gap threshold 
### FUNCTION ###################################################################
filterContacts <- function(Cp.v, ij.df, topCP, gap.bin.v, gap.bin.range, ct){
  
  cp.v <- rev(Cp.v) #21:1
  if(topCP < 0){ cp.v <- rev(cp.v) }
  ij.incl.TF <- (gap.bin.v >= gap.bin.range[1]) & (gap.bin.v <= gap.bin.range[2]) 
  ij.incl.TF <- ij.incl.TF & (ij.df[,"Cp"] %in% cp.v[1:abs(topCP)])
  if(ct%in%ct.v){
    ij.incl.TF <- ij.incl.TF & ij.df[,ct]>0
  }
  
  print("Contacts filtered.")
  
  return( which(ij.incl.TF) )
  
}

################################################################################

# rm(list=ls()); gc()