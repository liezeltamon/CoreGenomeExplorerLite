################################################################################
# ui() for app_run.R
### OTHER SETTINGS #############################################################
#ct.v <- c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
#          "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
#ct.v = c("ESC", "FC", "LC")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(shiny)
# library(shinydashboard)
# library(shinyWidgets)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Tab border colour
js <- '.nav-tabs-custom .nav-tabs li.active {
    border-top-color: #571666;
}"'

# Header
header <- dashboardHeader(title="Core Genome Explorer")

# Body
body <- dashboardBody(
  
  #tags$head(tags$style(HTML(
  #  '.box {margin-left: 20px; margin-right: -20px;}'
  #  ))),
  
  # Apply tab border colour
  tags$style(js),
  
  # Row for settings
  fluidRow(
    
    tabBox( 
      # Use id of tabBox to find the current tab
      title="Define contacts", id="seltab", width = 4, height = 325,
      tabPanel(title="Load",
               tags$div(tags$h4(
                 "Load contacts"
               ))
               ,
               tags$hr(),
               selectInput(inputId="chrRes", label="Dataset:", 
                           choices=paste0("chr", c(1:22, "X"), "-res:40000", sep=""), 
                           multiple=F, selected="chr1-res:40000")
               ),
      tabPanel(title="Select", 
               splitLayout(
                 tags$div(tags$h5( strong("Contact gap (in bins):") )),
                 tags$div(tags$h5( strong("      ") )),
                 tags$div(tags$h5( strong("Min:") )),
                 tags$div(tags$h5( strong("Max:") )),
                 cellWidths=c(220,20,125,125)
               ),
               splitLayout(
                 #numericInput(inputId="gap.bin.text", label=NULL, value=1000, min=0, step=1),
                 sliderInput(inputId="gap.bin", label=NULL, step=1, min=50, max=6230, value=c(1000,6230)),
                 tags$div(tags$h5( strong("      ") )),
                 #numericInput(inputId="gap.bin.text", label=NULL, value=1000, min=50, step=1),
                 numericInput(inputId="gap.bin.min.text", label=NULL, value=1000, min=50, step=1),
                 numericInput(inputId="gap.bin.max.text", label=NULL, value=6230, min=50, step=1),
                 cellWidths=c(220,20,125,125)
               ),
               
               splitLayout(
                 #selectInput(inputId="ct", label="Cell or tissue:", choices=c("All", sort(ct.v)), 
                 #            multiple=F, selected="All"),
                 #tags$div(tags$h5( strong("      ") )),
                 sliderInput(inputId="topCP", 
                             HTML(paste0("Contact persistence (c", tags$sub("p"), ") stringency:")),
                             #label=paste0("Contact persistence (Cp) stringency:", 
                             step=1, min=1, max=4, value=3),
                 cellWidths=c(220,20,250)
               )
               )
    ),
    
    tabBox( 
      title="Build structure", id="buildtab", width = 4, height = 325,
      tabPanel(title="Whole",
               tags$div(tags$h4(
                 "Build structure for whole chromosome"
               ))
               ,
               tags$hr(),
               numericInput(inputId="edge.res", label="Structure resolution (in bins)*", 
                            min=0, value=50, step=1),
              
               "*Resolution of 2 means that bins 1, 4, 7... are represented (when set to 0, 
               all bins are represented). Lower resolution to speed up display of 
               structure."
               ),
      tabPanel(title="Part", 
               tags$div(tags$h4(
                 "Build structure using input bins, genes or regions"
               ))
               ,
               tags$hr(),
               tags$div(tags$h5(
                 strong("Upload file (newline-sep. or 0-based BED):")
               )),
               splitLayout(
                 fileInput(inputId="inputfile", label=NULL, #"Upload list (line break-sep):",
                           multiple=F, accept=c("text", "text/plain")),
                 actionButton(inputId="addfeatbg", label="APPLY"),
                 actionButton(inputId="removefeatbg", label="RESET"),
                 downloadButton(outputId="downloadregbin", label="Region to bin file"),
                 cellWidths=c(200, 63, 100, 180)
               ),
               splitLayout(
                 tags$div(tags$h5( strong("OR Input (comma-sep.):") )),
                 tags$div(tags$h5( strong("") )),
                 tags$div(tags$h5( strong("") )),
                 tags$div(tags$h5( strong("Input type:") )),
                 cellWidths=c(200, 63, 100, 180)
               ),
               splitLayout(
                 textInput(inputId="inpbingene", label=NULL, value=NULL, 
                           placeholder="eg. 1,2,4:10 or RUNX1,TIAM1", width=335),
                 tags$div(tags$h5( strong("") )),
                 tags$div(tags$h5( strong("") )),
                 radioButtons(inputId="inptypnet", label=NULL, inline=T,
                              choices=c(`bins|regions`="bins|regions", genes="genes"),
                              selected=character(0)),
                 cellWidths=c(200, 63, 80, 180)
               )
               )
      ),
    
    tabBox( 
      # The id lets us use input$tabset1 on the server to find the current tab
      title="Study structure", id="studytab", width = 4, height = 325,
      tabPanel(title="GO|KEGG",
               tags$div(tags$h4(
                 "Analyse term enrichment"
               ))
               ,
               tags$hr(),
               pickerInput(inputId="GOinput", label="Genes on contacts in structure:",
                           choices="Not yet available", selected="Not yet available",
                           options=list(`actions-box`=T), multiple=T),
               actionButton(inputId="GOterm", label="Analyse"),
               downloadButton(outputId="downloadego", label="Result")
               ),
      tabPanel(title="Mark", 
               tags$div(tags$h4(
                 "Mark features on structure"
               ))
               ,
               tags$hr(),
               tags$div(tags$h5(
                 strong("Upload 0-based BED:")
               )),
               splitLayout(
                 fileInput(inputId="markfile", label=NULL, multiple=F,
                           accept=c("text/csv", "text/comma-separated-values,text/plain",".csv")),
                 #checkboxInput(inputId="header.bed", label="Header", value=F)
                 actionButton(inputId="addfeat", label="APPLY"),
                 actionButton(inputId="removefeat", label="RESET"),
                 actionButton(inputId="visualnetmarked", label="Mark"),
                 cellWidths=c(300, 63, 80, 65)
               ),
               #splitLayout(tags$div(tags$h5(
               #  strong("Coordinate system:")
               #)),
               #actionButton(inputId="visualnetmarked", label="Mark"),
               #cellWidths=c(300, 63)),
               #radioButtons(inputId="bed.CS", label=NULL,
               #             choices=c(`0-based`="0-based", `1-based`="1-based"), 
               #             selected="1-based", inline=T)
               splitLayout(
                 radioButtons(inputId="colstylenet", label="Style:", inline=T,
                              choices=c(`binary`="binary", discrete="discrete", continuous="gradient", 
                                        `continuous.#overlaps`="gradient.numolap"), selected="binary"),
                 uiOutput("hexpalette"),
                 cellWidths=c(425, 80)
               )
               )
      ),
    fluidRow(
      box(
        width=4,
        actionButton(inputId="loaddata", label="LOAD"),
        " then ",
        actionButton(inputId="setparam", label="SELECT"),
        tags$hr(),
        textOutput(outputId="Nij")
        ),
      box(
        width=4,
        actionButton(inputId="build", label="BUILD"),
        " then ",
        downloadButton(outputId="downloaddnetij", label="Contact data"),
        #actionButton(inputId="visual", label="VIEW"),
        tags$hr(),
        textOutput(outputId="dnetNij")
      ),
      box(
        width=4, 
        actionButton(inputId="filler", label="FROM THE SAHAKYAN LAB"),
        tags$hr(),
        textOutput(outputId="NGOterm")
      )
    ),
    
    # Row for networks
    fluidRow(
      column(
        width=6,
        #splitLayout(
          actionButton(inputId="visual", label="View"),
          downloadButton(outputId="downloadstr", label="Download as html"),
          "Note: Thicker edges mark centromere."
          #visNetworkOutput(outputId="network", width="800px", height="550px")
        #)
      )#,
      #column(
      #  width=6,
      #  actionButton(inputId="matchnodes", label="Match"),
      #  downloadButton(outputId="downloadstrM", label="Download as html")#,
      #  #visNetworkOutput(outputId="net.marked", width="800px", height="550px")
      #)
    ),
    
    #fluidRow(
    #  align="center",
    #  box(
    #    height="600px", 
    #    visNetworkOutput(outputId="network", width="800px", height="550px")
    #  ),
    #  box(
    #    height="600px", 
    #    visNetworkOutput(outputId="net.marked", width="800px", height="550px")
    #  )
    #),
    
    # For screenshot of structures
    column(
      width=11, style = "background-color:#FFFFFF;",
      visNetworkOutput(outputId="network", height="1100px")
    ), 
    fluidRow(
      column(
        width=6,
        actionButton(inputId="matchnodes", label="Match"),
        downloadButton(outputId="downloadstrM", label="Download as html")
      )
    ),
    column(
      width=11, style = "background-color:#FFFFFF;",
      visNetworkOutput(outputId="net.marked", height="1100px")
    ),
    ## For screenshot of structures

    # Row for contact data and GO/KEGG tables
    fluidRow(
      column(
        width=6,
        DT::dataTableOutput(outputId="dnetij", width="100%"),
        br(),
        br()
      ),
      column(
        width=6,
        DT::dataTableOutput(outputId="GOtable", width="100%"),
        br(),
        br()
      )
    )
    
    
  )
  
)

ui <- dashboardPage(header=header, sidebar=dashboardSidebar(disable=T), body=body, skin="purple")

# rm(list=ls()); gc()






