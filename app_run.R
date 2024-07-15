################################################################################
# Shiny app to view and explore our core, persistent, genome organisation.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer/A2_app_cge"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
data.dir = paste0(wk.dir, "/data")
persist.dir = paste0(data.dir, "/out_basePersist")
anno.file = paste0(data.dir, "/hg19anno_ALL")
hugoEntrez.file = paste0(data.dir, "/hg19anno_SYMBOLtoENTREZID_052020")
centrobed.file = paste0(data.dir, "/ct_hg19_foi_centromoreonly_desc_DNA")
chrLen.file = paste0(data.dir, "/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
options(shiny.maxRequestSize = 50 * 1024^2) # 5 Mb (default) -> 10 Mb maximum 
# file size for file upload
gcb = "min2Mb"
#ct.v <- c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
#          "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
ct.v = c("ESC", "FC", "LC") # Abridged PERSIST.MX only contain this
Cp.v = 1:21
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(clusterProfiler)
library(compiler)
library(data.table)
library(dplyr)
library(DT) # for dataTableOutput()
library(GenomicRanges)
library(ggplot2)
library(org.Hs.eg.db)
library(RColorBrewer)
library(shiny)
library(shinydashboard)
library(shinyWidgets) # pickerInput() and updatePickerInput()
library(visNetwork)
source(paste0(wk.dir, "/lib/funxAnno.R"))
source(paste0(wk.dir, "/lib/funxAnnoWrapper.R"))
source(paste0(wk.dir, "/lib/GEN_WhichOverlap.R"))
source(paste0(wk.dir, "/lib/markCentromereOnNetworkData.R"))
source(paste0(wk.dir, "/lib/makeNetworkData.R"))
source(paste0(wk.dir, "/lib/addNodesToNetwork.R"))
source(paste0(wk.dir, "/lib/modifyNetworkData.R"))
source(paste0(wk.dir, "/app_ui.R"))
source(paste0(wk.dir, "/app_server.R"))
source(paste0(wk.dir, "/lib/server_makeANNODF.R"))
source(paste0(wk.dir, "/lib/server_makeBINGENEDF.R"))
source(paste0(wk.dir, "/lib/server_filterContacts.R"))
source(paste0(wk.dir, "/lib/server_makeBaseNetwork.R"))
source(paste0(wk.dir, "/lib/server_readINPUTFILE.R"))
source(paste0(wk.dir, "/lib/server_getINPBINGENE.R"))
source(paste0(wk.dir, "/lib/server_selectNodes.R"))
source(paste0(wk.dir, "/lib/server_makeDisplayNetwork.R"))
source(paste0(wk.dir, "/lib/server_makeDisplayNetworkData.R"))
source(paste0(wk.dir, "/lib/server_makeEGO.R"))
source(paste0(wk.dir, "/lib/server_getEGOGENES.R"))
source(paste0(wk.dir, "/lib/server_makeNETMARKED.R"))
source(paste0(wk.dir, "/lib/server_makeNETMARKEDlegendData.R"))
source(paste0(wk.dir, "/lib/server_outTableContacts.R"))
source(paste0(wk.dir, "/lib/server_outTableFunxTerm.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
shinyApp(ui=ui, server=server)
################################################################################
# rm(list=ls()); gc()


