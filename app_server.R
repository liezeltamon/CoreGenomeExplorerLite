################################################################################
# server() for app_run.R
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# Refer to app_run.R
### FUNCTION ###################################################################
server <- function(input, output, session){
  
  chrLen.df <- reactive({
    chrLen.df <- try( read.delim(file=chrLen.file, stringsAsFactors=F, header=T) )
    validate(need( chrLen.df, paste0("Error: Reading chromosome length file.") )) 
    return(chrLen.df)
  })
  
  anno.df <- reactive({
    anno.df <- try( makeANNODF(anno.file=anno.file) )
    validate(need( anno.df, paste0("Error: Generating transcript annotation object.") )) 
    return(anno.df)
  })
  
  #
  output$hexpalette <- renderUI({
    if (input$colstylenet == "binary") {
      textInput(inputId="colour.nme", label="Hex code:", value="#11DBCE", placeholder="")
    } else if(input$colstylenet == "discrete"){
      selectInput(inputId="colour.nme", label="Palette:", choices=rownames(brewer.pal.info), 
                  multiple=F, selected="Dark2")
    } else {
      selectInput(inputId="colour.nme", label="Palette:", choices=rownames(brewer.pal.info), 
                  multiple=F, selected="RdYlBu")
    }
  })
  
  #-------------------Dependency on LOAD (to urge user to think twice before loading a dataset)
  
  chr <- eventReactive(input$loaddata, {
    return( strsplit(x=input$chrRes, split="-res:")[[1]][1] )
  })
  
  res <- eventReactive(input$loaddata, {
    return( as.numeric(strsplit(x=input$chrRes, split="-res:")[[1]][2]) )
  })
  
  bin.last <- eventReactive(input$loaddata, {
    chrLen.df <- chrLen.df()
    return( ceiling(chrLen.df[chrLen.df[,"chromosome"] == chr(),"length.bp"] / res()) )
  })
  
  max.gap <- eventReactive(input$loaddata, {
    return( bin.last() - 1 - 1 )
  })
  
  observe({
    bin.last <- bin.last(); max.gap <- max.gap()
    updateSliderInput(session, inputId="gap.bin", value=c(input$gap.bin.min.text, input$gap.bin.max.text), min=50, max=max.gap, step=1)
    updateNumericInput(session, inputId="edge.res", value=50, min=50, max=max.gap, step=1)
  })
  
  # Sync numericInput and sliderInput for contact gap 
  observeEvent(input$gap.bin.text, {
    
    validate(need( try(is.numeric(input$gap.bin.text)), "Error: Invalid contact gap." )) 
    max.gap <- max.gap()
    vals <- c(input$gap.bin.min.text, input$gap.bin.max.text)
    if( any(vals == "") ){ # While user is typing
      updateSliderInput(session, inputId="gap.bin", value=c(0,0), min=50, max=max.gap, step=1)
    } else if ( all(vals != "") & !identical(vals, input$gap.bin) & diff(vals) >= 0 ){ 
      updateSliderInput(session, inputId="gap.bin", value=vals, min=50, max=max.gap, step=1)
    }
    
  }) 
  observeEvent(input$gap.bin, {
    updateNumericInput(session, inputId="gap.bin.min.text", value=input$gap.bin[1], min=50, max=max.gap(), step=1)
    updateNumericInput(session, inputId="gap.bin.max.text", value=input$gap.bin[2], min=50, max=max.gap(), step=1)
  })
  
  ij.df <- eventReactive(input$loaddata, {
    persist.mx <- try( load(file=paste0(persist.dir, "/", chr(), "_Persist_", gcb, "_topCP4.RData")) )
    validate(need( persist.mx, paste0("Error: Loading ", chr(), " contact data.") ))
    return( cbind(PERSIST.MX$hits, Cp=PERSIST.MX$ntis) )
  })
  
  gap.bin.v <- eventReactive(input$loaddata, {
    ij.df <- ij.df()
    return( ij.df[,"j"] - ij.df[,"i"] - 1 )
  })
  
  anno.TF <- eventReactive(input$loaddata, {
    return( anno.df()$chr == chr() ) 
  })
  
  bingene.df <- eventReactive(input$loaddata, {
    anno.df <- anno.df(); bin.last <- bin.last(); anno.TF <- anno.TF()
    bingene.df <- try( makeBINGENEDF(res=res(), anno.TF=anno.TF, bin.last=bin.last, chr=chr(), anno.df=anno.df) )
    validate(need( bingene.df, paste0("Error: Generating bin-gene mapping for node hover panel.") )) 
    return(bingene.df)
  })
  
  #-------------------Dependency on SELECT
  
  # Select contacts to be displayed based on Cp, gap and cell/tissue
  ij.incl.ind <- eventReactive(input$setparam, {
    ij.df <- ij.df(); gap.bin.v <- gap.bin.v() # define first to avoid error from function in case one of these encounter validate
    return( filterContacts(Cp.v=Cp.v, ij.df=ij.df, topCP=input$topCP, gap.bin.v=gap.bin.v, 
                           gap.bin.range=input$gap.bin, ct="All") ) #input$ct) ) # Can't use try() cause object may be empty
  }) 
  
  # Print number of contacts (whole chr)
  output$Nij <- renderText({
    ij.incl.ind <- ij.incl.ind()
    paste0( "Whole ", chr(), " has ", length(ij.incl.ind), " contacts based on set parameters." )
  })
  
  # Base whole network (structure resolution = 1, used for quick isolation of part of chromosome structure)
  BNET <- eventReactive(input$setparam, {
    bin.last <- bin.last(); ij.df <- ij.df(); ij.incl.ind <- ij.incl.ind(); bingene.df <- bingene.df()
    BNET <- try( makeBaseNetwork(ij.df=ij.df, ij.incl.ind=ij.incl.ind, chr=chr(), bin.last=bin.last, 
                                 res=res(), bingene.df=bingene.df, centrobed.file=centrobed.file, edge.res=0) )
    validate(need(BNET, paste0("Error: Generating base whole network.") ))
    return(BNET)
  })
  
  #-------------------Dependency on BUILD
  
  uploadbg <- reactiveValues(status=NULL)
  observeEvent(input$addfeatbg, {uploadbg$status <- "uploaded"} )
  observeEvent(input$removefeatbg, {uploadbg$status <- "removed"} )
  
  inputfile <- reactive({
    
    if( is.null(uploadbg$status) ){
      return(NULL)
    } else if(uploadbg$status == "uploaded" & !is.null(input$inputfile$datapath) ){
     
      #bed.df <- try(read.table(file=input$inputfile$datapath, header=FALSE, na.strings=c(""), stringsAsFactors=FALSE))
      bed.df <- try(data.table::fread(file=input$inputfile$datapath, header=FALSE, na.strings="", stringsAsFactors=FALSE, data.table=FALSE))
      validate(need( bed.df, "Error in reading input file." ))
      col.len <- ncol(bed.df)
      validate(need( col.len != 2, "Error: Required number of columns in input file not met." ))
      validate(need( try( na.fail(bed.df[,1:min(c(col.len, 3))])), "Error: Missing value/s in input file required columns." ))
      if(col.len >= 3){
        #validate(need( all(bed.df[,1] %in% paste0("chr", c(1:22, "X", "Y", "MT"))) , "Error: Invalid input BED chromosome name/s." ))
        validate(need( try(is.numeric(c(bed.df[,2], bed.df[,3]))), "Error: Non-numeric input BED coordinate/s." ))
        bed.df <- bed.df[,1:3]
        res <- res()
        bed.df[,2] <- bed.df[,2] + 1  # Convert to 1-based
        bed.df <- data.frame(bed.df, ceiling( bed.df[,2] / res ), ceiling( bed.df[,3] / res ), stringsAsFactors=FALSE)
        bed.df[,2] <- bed.df[,2] - 1  # Revert to 0-based
        colnames(bed.df) <- c("chr", "start", "end", "start.bin", "end.bin")
      }
      
      return(bed.df)
      
    } else {
      return(NULL)
    } 
   
  })
  
  output$downloadregbin <- downloadHandler(
    
    filename=function(){ paste0(out.id2() ,"_regionToBin.csv") }, 
    content=function(file){
      inputfile <- inputfile()
      if( !is.null(inputfile) ){
        write.csv(inputfile, file, row.names=F)
      }
    }
    
  )
  
  inpbingene <- eventReactive(input$build, {
    
    if(input$buildtab == "Part"){ 
      validate(need( input$inptypnet, "Error: Specify input type." )) 
      inputfile <- inputfile()
      validate(need( !is.null(inputfile) | input$inpbingene != "", "Error: No file uploaded/applied or no input given." )) 
      validate(need( !(!is.null(inputfile) & input$inpbingene != ""), "Error: Both file and input given. Choose one." )) 
      return( getINPBINGENE(inpbingene=input$inpbingene, inputfile=inputfile, chr=chr()) )
    } else {
      return(NULL)
    }
  
  })
  
  selnode <- eventReactive(input$build, {
    
    inpbingene <- inpbingene()
    if( input$buildtab == "Part" & !is.null(inpbingene) ){
      bin.last <- bin.last(); inptypnet <- input$inptypnet; anno.df=anno.df(); anno.TF=anno.TF()
      selnode <- selectNodes(inpbingene=inpbingene, inptypnet=inptypnet, buildtab=input$buildtab, 
                             anno.df=anno.df, anno.TF=anno.TF, res=res(), bin.last=bin.last)
      validate(need( selnode, paste0("Error: Invalid input ", inptypnet, " for ", chr(), ".") ))
    } else {
      return(NA)
    }
  
    return(selnode)
    
  })
  
  DNET <- eventReactive(input$build, {
  
    bin.last <- bin.last(); edge.res <- input$edge.res
    validate(need( try(is.numeric(edge.res)), "Error: Invalid structure resolution." ))
    validate(need( edge.res >= 0 & edge.res <= bin.last, "Error: Invalid structure resolution." ))
    showModal(modalDialog("Building structure...", footer=NULL, fade=F, size="s", easyClose=T))
    ij.df <- ij.df(); ij.incl.ind <- ij.incl.ind(); bingene.df <- bingene.df(); BNET <- BNET(); selnode <- selnode()
    DNET <- try( makeDisplayNetwork(BNET=BNET, selnode=selnode, bin.last=bin.last, ij.df=ij.df, ij.incl.ind=ij.incl.ind, 
                                    res=res(), bingene.df=bingene.df, chr=chr(), edge.res=edge.res, centrobed.file=centrobed.file) )
    validate(need( DNET, "Error: Generating display structure object." )) 
    removeModal()
    return(DNET)
    
  })
  
  # Print number of contacts (subset chr)
  output$dnetNij <- renderText({
    DNET <- DNET()
    paste0( chr(), " structure has ", sum(!is.na(DNET$edges$label)), " contacts. Check contact data below." )
  })
  
  # Contact data of displayed structure
  dnetij0 <- eventReactive(input$build, {
    
    chrLen.df <- chrLen.df(); anno.df <- anno.df(); bin.last <- bin.last()
    DNET <- DNET(); edges <- DNET$edges
    chr.len <- chrLen.df[chrLen.df[,"chromosome"] == chr(),"length.bp"]
    dnetij0 <- try( makeDisplayNetworkData(noi.df=edges[!is.na(edges$label), c("from", "to", "label")], 
                                           anno.df=anno.df, nodes=DNET$nodes, edges=edges,
                                           res=res(), bin.last=bin.last, chr.len=chr.len) )
    validate(need(dnetij0, "Error: Generating contact data table or No contacts displayed.")) 
    dnetij0 <- dnetij0[,setdiff(colnames(dnetij0), "uniqueID")]
    return(dnetij0)
    
  })

  dnetij <- reactive({
    
    dnetij <- dnetij0(); MNET <- MNET()
    if( !is.null(MNET) & !is.null(dnetij) ){
      match.ind <- match(as.numeric(dnetij$bin), as.numeric(MNET$nodes$id))
      dnetij$feature <- MNET$nodes$name.feat.olap[match.ind]
      dnetij$score <- MNET$nodes$score.olap.mean[match.ind]
    }
    return(dnetij)
    
  })
  
  #- Dependency on VIEW
  network <- eventReactive(input$visual, {
    DNET <- DNET()
    return( visNetwork(nodes=DNET$nodes, edges=DNET$edges, width=1200, height=1000) %>%
              visLayout(randomSeed = 278) )
  })
  
  #-------------------STUDY structure
  
  #- Analyse term enrichment
  
  egogenes <- eventReactive(input$build, {
    dnetij <- dnetij()
    return( getEGOGENES(dnetij=dnetij) )
  })
  
  observeEvent(input$build, {
    egogenes <- egogenes()
    updatePickerInput(session=getDefaultReactiveDomain(), selected=egogenes, inputId="GOinput", choices=egogenes)
  })
  
  ego <- eventReactive(input$GOterm, {
    validate(need(input$GOinput != "Not yet available", "No selected genes for enrichment analysis.") )
    ego <- try( makeEGO(genes=input$GOinput, bgrgenes=unique(anno.df()$gene), hugoEntrez.file=hugoEntrez.file) )
    validate(need( ego, "Error: Enrichment analysis." ))
    return(ego)
  })
  
  # Print number when enrichment analysis done.
  NGOterm <- eventReactive(input$GOterm, {
    paste0("Enrichment analysis done. Check results below." )
  })
    
  output$NGOterm <- renderText({
    NGOterm()
  })
  
  #- Mark features
  
  upload <- reactiveValues(status=NULL)
  observeEvent(input$addfeat, {
    upload$status <- "uploaded"
  })
  observeEvent(input$removefeat, {
    upload$status <- "removed"
  })
  
  markfile <- reactive({

    if( is.null(upload$status) ){
      return(NULL)
    } else if( upload$status == "uploaded" & !is.null(input$markfile$datapath) ){
      
      #bed.df <- try(read.delim(file=input$markfile$datapath, header=FALSE, na.strings=c(""), stringsAsFactors=FALSE))
      bed.df <- try(data.table::fread(file=input$markfile$datapath, header=FALSE, na.strings="", stringsAsFactors=FALSE, data.table=FALSE))
      validate(need( bed.df, "Error: Reading marker BED file." ))
      validate(need( ncol(bed.df) >= 3 , "Error: Required number of columns in marker BED not met." ))
      validate(need( try(na.fail(bed.df[,1:3])) , "Error: Missing value/s in marker BED required columns." ))
      validate(need( is.numeric(c(bed.df[,2], bed.df[,3])) , "Error: Non-numeric marker BED coordinates." ))
      return(bed.df)
      
    } else {
      return(NULL)
    }
    
  })
  
  MNET <- reactive({
    
    if( input$visualnetmarked != 0 & !is.null(input$colstylenet) & !is.null(markfile()) ){
      
      markfile <- markfile() 
      DNET <- DNET()
      MNET <- try( makeNETMARKED(markfile=markfile, DNET=DNET, res=res(), chr=chr(), bed.CS="0-based", #input$bed.CS, 
                                 centrobed.file=NULL, header.bed=FALSE, colouring.style=input$colstylenet, colour.nme=input$colour.nme) )
      validate(need( MNET, "Error: Generating marked structure object." ))
      return(MNET) #input$header.bed) )
    } else {
      return(NULL)
    }
    
  })
  
  MNET.leg <- eventReactive(input$visualnetmarked, {
    
    MNET.leg <- list()
    MNET.leg$nodes <- makeNETMARKEDlegendData(MNET=MNET(), colstylenet=input$colstylenet, binary.col=input$colour.nme)
    MNET.leg$edges <- data.frame(label=paste0("No \n bin/s"), color="black", width=3, font.size=25, 
                                 font.color="darkred", length=700, arrow="to", arrowStrikethrough=FALSE)
    return(MNET.leg)
    
  })
  
  net.marked <- eventReactive(input$visualnetmarked, {
    
    MNET <- MNET()
    validate(need( !is.null(MNET), "Error: No BED file uploaded or applied." )) 
    MNET.leg <- MNET.leg()
    return( visNetwork(nodes=MNET$nodes, edges=MNET$edges, width=1200, height=1000) %>%
            visNodes(physics=F) %>%
            visLayout(randomSeed = 278) %>%
            visLegend(addNodes=MNET.leg$nodes, useGroups=FALSE, width=0.07, addEdges=MNET.leg$edges) 
            ) 

  })
  
  net.marked.dld <- eventReactive(input$visualnetmarked, {
    
    net.marked <- net.marked()
    MNET.leg <- MNET.leg()
    return( net.marked %>%
              visNodes(physics=T)  %>%
              visLayout(randomSeed = 278) )
            #%>%
            #  visLayout(randomSeed = 278) %>%
            #  visLegend(addNodes=MNET.leg$nodes, useGroups=FALSE, width=0.07, addEdges=MNET.leg$edges) )  
    
  })
  
  # Marked network mimics displayed network and this part adjusts it when the latter is adjusted manually
  vals <- reactiveValues(coords=NULL)
  
  observe({
    
    invalidateLater(1000)
    visNetworkProxy("network") %>% visGetNodes()
    network_nodes <- input$network_nodes
    if( !is.null(network_nodes) ){
      network_nodes <- lapply(X=network_nodes, FUN=function(node){
        unlist(node[c("id", "x", "y")])
      })
      get.nodes <- do.call(rbind.data.frame, network_nodes)
      rm(network_nodes)
      colnames(get.nodes) <- c("id", "x", "y")
      vals$coords <- get.nodes
    }
    
  })
  
  observeEvent(input$matchnodes,{
    if( !is.null(net.marked()) ){
      visNetworkProxy("net.marked") %>% visUpdateNodes(nodes=vals$coords) %>%
        visLayout(randomSeed = 278)
    }
  })
  
  #---------------------------------------OUTPUTS AND DOWNLOADS
  
  # Output file ids
  
  ## Output id for downloaded genes and GO term results
  out.id0 <- reactive({
    #return( paste0(gcb, "_", chr(), "_", input$ct, "_topCP", input$topCP, "_gapBin", input$gap.bin) )
    return( paste0(gcb, "_", chr(), "_All_topCP", input$topCP, "_gapBinMin", input$gap.bin[1], "_Max", input$gap.bin[2]) )
  })
  
  ## For downloaded network
  out.id <- reactive({
    return( paste0(out.id0(), "_edgeRes", input$edge.res, "_edgelenmultDEF") )
  })
  
  ## For input bin/gene specific output
  out.id1 <- reactive({
    inpbingene <- inpbingene()
    inpbingene.id <- ifelse(nchar(inpbingene) < 75, inpbingene, paste0(substr(inpbingene, start=1, stop=75), "ETC"))
    return( paste0(out.id(), "_input", inpbingene.id) )
  })
  
  out.id2 <- eventReactive(input$build, { 
    tme <- gsub(x=format(Sys.time(), "%X"), pattern=":", replacement="") 
    return( paste0(Sys.Date(), "_time", tme) )
  })
  
  #-------------------
  
  # Displayed network
  output$network <- renderVisNetwork({
    return( network() )
  })
  
  ## Download displayed network
  output$downloadstr <- downloadHandler(
    
    filename=function(){ paste0(out.id(), ".html") }, 
    content=function(file){
      network() %>%
        visExport(type="png", name=out.id(), label=paste0("Export as png"), float="right") %>%
        visSave(file, selfcontained=T, background="white")
    }
    
  )
  
  # Displayed network - marked
  output$net.marked <- renderVisNetwork({
    net.marked()
  })
  
  ## Download displayed network - marked
  output$downloadstrM <- downloadHandler(
    
    filename=function(){ paste0(out.id(), "_marked.html") }, 
    content=function(file){
      net.marked.dld() %>%
        visExport(type="png", name=paste0(out.id(), "_marked"), label=paste0("Export as png"), float="right") %>%
        visSave(file, selfcontained=T, background="white")
    }
    
  )
  
  # Table of contacts in displayed structure
  output$dnetij <- DT::renderDataTable({
    dnetij <- dnetij()
    validate(need(!is.null(dnetij), "No contact data to display or download.") )
    return( datatable(outTableContacts(dnetij=dnetij), options=list(pageLength=5), rownames=F) )
  })
  
  ## Download contact data table
  output$downloaddnetij <- downloadHandler(
    
    filename=function(){
      paste0(out.id1(), ".csv")
    }, 
    content=function(file){
      write.csv(dnetij(), file, row.names=F)
    }
    
  )
  
  # Table of enriched GO/KEGG terms in genes overlapping contacts of displayed structure
  output$GOtable <- DT::renderDataTable({
    ego <- ego()
    validate(need(nrow(ego) > 0, "No term/s returned."))
    return( datatable(outTableFunxTerm(ego=ego), options=list(pageLength=5), rownames=F) )
  })
  
  ## Download GO/KEGG table
  output$downloadego <- downloadHandler(
    filename=function(){ paste0(out.id1(), "_GOKEGG-", out.id2(), ".csv") }, 
    content=function(file){
      write.csv(ego(), file, row.names=F)
    }
  )
  
}

# rm(list=ls()); gc()



