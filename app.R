# load packages
library(shiny)
library(shinyjs)
library(shinythemes)
library(DT)
library(plotly)
library(reshape2)
library(RColorBrewer)

# laod functions
source("function.R")

library(data.table)
setDTthreads(18)

# load data
library(fst)

gene_symbol <- read.fst("data/gene_symbol.fst")
gnomad_exome <- read.fst("data/gnomad_exome.fst")
row.names(gnomad_exome)<-gnomad_exome$V1
summary <- read.fst("data/summary.fst")
gtrM <- read.fst("data/gtrM.fst")
gtrS <- read.fst("data/gtrS.fst")
gpt <- read.fst("data/gpt.fst")
# genes_by_ccds_id <- read.fst("data/genes_by_ccds_id.fst")
tP <- read.fst("data/test_to_pheno.fst")

# tP2 <- tP[tP$phenotype_name == "Tuberous sclerosis 1",]
# unique(gpt[gpt$GTR_accession %in% tP2$AccessionVersion,"gene_symbol"])

tP$AccessionVersion <- as.character(tP$AccessionVersion)
tP$test_name <- as.character(tP$test_name)
tP$phenotype_name <- as.character(tP$phenotype_name)
# length(unique(as.character(tP$phenotype_name[tP$AccessionVersion %in% gpt$GTR_accession[ gpt$gene_symbol == "SIK1" ]])))
## ccds id -> ensembl id (for coverage plot & link to gnomAD browser)
ccds2ens<-readRDS("data/ccds_ens_map.rds")


# Define UI ----
ui <- fluidPage(
  useShinyjs(),
  #theme = shinytheme("cerulean"),
  #titlePanel("WEScover"),
  tags$head(tags$style(".modal-dialog{min-width:1200px}")),
  navbarPage("WEScover",
    tabPanel("Report",
      sidebarLayout(position = "left",
      sidebarPanel(
              tags$h2("User input"),
              
              div(style="display: inline-block;vertical-align:baseline; width: 70%;",
                  selectizeInput("phen",
                                 label="Phenotype",
                                 choices = NULL,
                                 multiple = TRUE)
              ),
              div(style="display: inline-block;vertical-align:baseline; width: 25%;",
                  actionButton("fGPT", "Filter")
              ),
              
              div(style="display: inline-block;vertical-align:baseline; width: 70%;",
                  selectizeInput("gpt",
                             label="GPT name",
                             choices = NULL,
                             multiple = TRUE)
                ),
              div(style="display: inline-block;vertical-align:baseline; width: 25%;",
                  actionButton("fGenes", "Filter")
                ),
              
              # selectizeInput("gene_symbol",
              #                label = "Gene symbol",
              #                choices = NULL,
              #                multiple = TRUE),
              selectizeInput("gene_symbol", "Gene symbol",
                             choices = NULL,
                             multiple = TRUE,
                             options = list(
                               splitOn = I("(function() { return /[,; ]/; })()"),
                               create = I("function(input, callback){
                                return {
                                  value: input,
                                  text: input
                                  };
                                }")
                             ),
                             width = '70%'
                            ),
               selectInput("depth_of_coverage", 
                           label = "Depth of coverage",
                           choices = c("10x", "20x", "30x"),
                           selected = "20x",
                           width = '70%'),
        #       selectizeInput("x", "Choose:",
        #                      letters,
        #                      multiple = TRUE,
        #                      options = list(
        #                        splitOn = I("(function() { return /[, ]/; })()"),
        #                        create = I("function(input, callback){
        #   return {
        #     value: input,
        #     text: input
        #    };
        # }")
        #                      )
        #       ),
        actionButton("clear", "Clear inputs", class = "btn-secondary"),
        actionButton("update", "Submit query", class = "btn-primary")
        ),
        mainPanel(
               dataTableOutput('tableMain')  
        )
      )
    ),
    tabPanel("Help"
    )
  )#,
  # tags$footer(title="© Boston Children's Hospital. All Rights Reserved.", 
  #   align = "right", 
  #   style = "
  #     position:absolute;
  #     bottom:0;
  #     width:100%;
  #     height:50px; /* Height of the footer */
  #     color: white;
  #     padding: 10px;
  #     background-color: #DCDCDC;
  #     z-index: 1000;"
  # )
)

# Define server logic ----
server <- function(input, output, session) {
  
  # updateSelectizeInput(session, 'x', choices = letters, server = TRUE)
  updateSelectizeInput(session, 'phen', choices = sort(unique(tP$phenotype_name)), 
                       server = TRUE)
    updateSelectizeInput(session, 'gene_symbol', 
                       choices = gene_symbol$gene_symbol,
                       server = TRUE)
  updateSelectizeInput(session, 'gpt', choices = sort(unique(gpt$test_name)), server = TRUE)
  observeEvent (input$clear,{
    updateSelectizeInput(session, 'phen', choices = sort(unique(tP$phenotype_name)), server = TRUE)
    updateSelectizeInput(session, 'gene_symbol', 
                         choices = gene_symbol$gene_symbol, server = TRUE,
                         label = "Gene symbol")
    updateSelectizeInput(session, 'gpt', choices = sort(unique(gpt$test_name)), server = TRUE,
                         label = "GPT name")
    updateSelectInput(session, "depth_of_coverage", choices = c("10x", "20x", "30x"), selected = "20x")
  })
  
  observeEvent (input$fGPT,{
    listPhe <- unique(gpt$test_name[ gpt$GTR_accession %in% tP$AccessionVersion[ tP$phenotype_name %in% as.character(input$phen) ] ])
    updateSelectizeInput(session, 'gpt', choices = sort(listPhe), server = TRUE,
                         label = "GPT name (filtered)")
    listGenes <- unique(gpt$gene_symbol[ gpt$GTR_accession %in% tP$AccessionVersion[ tP$phenotype_name %in% as.character(input$phen) ]])
    updateSelectizeInput(session, 'gene_symbol', choices = sort(listGenes), server = TRUE,
                         label = "Gene symbol (filtered)")
  })
  
  observeEvent (input$fGenes,{
    listGenes <- gpt$gene_symbol[ gpt$test_name %in% as.character(input$gpt) ]
    updateSelectizeInput(session, 'gene_symbol', choices = sort(listGenes), server = TRUE,
                         label = "Gene symbol (filtered)")
  })
  
  # generator of buttons for the main table
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
  # list of global values
  myValue <- reactiveValues(summary_table = NA, GPT_table = NA, 
    violon_population = NA, gene = "", ccds="", gnomAD_plot="")
  
  # observer to capture detail buttons in the main table
  observeEvent(input$detail_button, {
    selectedRow <- as.numeric(strsplit(input$detail_button, "_")[[1]][2])
    #message(paste('click on detail ', selectedRow))
    geneS <- as.character(main_table()[selectedRow, 1])
    
    myValue$summary_table <<- 
      createMainTable2(geneS, input$depth_of_coverage, summary, gtrM, gtrS)[,c(1,2,6:10)]
    myValue$GPT_table <<- 
      createGPT(selectedRow, main_table(), gtrM)
    myValue$violon_population <<- createPlot(geneS, main_table()[selectedRow, 2])
    myValue$gene <<- geneS
    myValue$ccds <<- main_table()[selectedRow, 2]
    fileAD <- 
      paste0("coverage_plots/", myValue$ccds, ".png")
    if(file.exists(fileAD)) {
      myValue$gnomAD_plot <<- list(
        src = fileAD,
        contentType = 'image/png',
        width = "100%",
        #height = 300,
        alt = paste0("gnomAD coverage for ", myValue$ccds))
    } else {
      myValue$gnomAD_plot <<- list(
        src = "coverage_plots/NO_CCDS_GAD.png",
        contentType = 'image/png',
        width = "100%",
        #height = 300,
        alt = paste0("No coverage from gnomAD for ", myValue$ccds))
    }
    showModal(modal_main())
    
  })
  
  # create reactive main table
  main_table <- reactive ({
    geneS <- c()
    
    if(length(input$phen) != 0) {

    if (length(input$gene_symbol) != 0) {
      geneS <- c(geneS, input$gene_symbol)
    } else if (length(input$gpt) != 0) {
      geneG <- as.character(unique(gpt[gpt$test_name %in% input$gpt, "gene_symbol"]))
      geneS <- c(geneS, geneG)
    } else if(length(input$phen) != 0) {

      message("OK")
      #input <- list(pehn = "Tuberous sclerosis 1")
      geneP <- unique(gpt$gene_symbol[ gpt$GTR_accession %in% tP$AccessionVersion[ tP$phenotype_name %in% as.character(input$phen) ]])
      message(length(geneP))
      geneS <- c(geneS, geneP)
    }
    }

    if (length(input$gpt) != 0) {
      geneG <- as.character(unique(gpt[gpt$test_name %in% input$gpt, "gene_symbol"]))
      geneS <- c(geneS, geneG)
    }
    
    if (length(input$gene_symbol) != 0) {
      geneS <- c(geneS, input$gene_symbol)
    }
    
    geneS <- unique(geneS)
    tbl <- createMainTable2(geneS, input$depth_of_coverage, summary, gtrM, gtrS)
    if(ncol(tbl) > 0) {
      tbl <- tbl[ , c(1:3,11:12)]
      tbl$Action <- shinyInput(actionButton, nrow(tbl), 'button_', 
                               label = "Detail", 
                               onclick = 'Shiny.onInputChange(\"detail_button\",  this.id)' )
    }
    tbl
    })
  
  # display reactive main table
  output$tableMain <- renderDataTable( {
    input$update
    isolate({
      if(ncol(main_table()) > 0) {
        formatStyle(datatable(main_table(), escape = FALSE, selection = 'none'), 
                  columns = "Global coverage (mean, %)", target = 'row', backgroundColor = 
                    styleInterval(cuts = c(95, 98), values=c("#FFE4E1", "#FFFFE0", "#F0FFF0")))
      }
    })
  })
  
  # modal window to show tables and plots
  modal_main <- function(failed = FALSE){
    modalDialog(size="l",
      tabsetPanel(type="tabs",
        tabPanel("Population summary",
          dataTableOutput('summary_table')),
        tabPanel("Coverage plots",
                 fluidRow(align="center",
                   column(4, plotOutput("violin_population"), div(style="max-height:100px;")),
<<<<<<< HEAD
                   column(8, imageOutput("gnomAD_plot"))
                 )
                 # fluidRow(align="center",
                 #          textOutput("gnomAD_address"))),
        ),
=======
                   column(8, tags$a(imageOutput("gnomAD_plot"),href=paste0("http://gnomad.broadinstitute.org/gene/", ccds2ens[as.character(myValue$ccds), 2]), target="_blank"))
                 )),
>>>>>>> db201e22460e7ca882a0716f96afeb7ae0da3be5
#        tabPanel("Violin per population",
#                 plotOutput("violin_population")),
#        tabPanel("gnomAD coverage",
#                 imageOutput("gnomAD_plot")),
        tabPanel("Gene panels",
                 dataTableOutput('GPT_table'))
      ),
      easyClose = TRUE
    )
  }

  # table displayed in our modal
  output$summary_table <- renderDataTable({
    myValue$summary_table
  }, escape = FALSE)
  
  output$violin_population <- renderPlot({
    createPlot(myValue$gene, myValue$ccds)
  })
  
  output$GPT_table <- renderDataTable({
    myValue$GPT_table
  }, escape = FALSE)
  
  output$gnomAD_plot <- renderImage({
    myValue$gnomAD_plot
  }, deleteFile = FALSE)
  
  createPlot <- function(gene, selCCDS) {
    message("[PLOT] Number of CCDS: ", selCCDS, "; Number of genes:", gene )
    selCCDS<-as.character(selCCDS)
    
    idx <- 0 
    if (input$depth_of_coverage == "10x") {
      load10x(summary)
      idx <- 2
      dta <- do.call(rbind, list(melt(AFR_10x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(AMR_10x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(EAS_10x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(EUR_10x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(SAS_10x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol"))))
    }
    
    if (input$depth_of_coverage == "20x") {
      load20x(summary)
      idx <- 3
      dta <- do.call(rbind, list(melt(AFR_20x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(AMR_20x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(EAS_20x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(EUR_20x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(SAS_20x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol"))))
      
    }
    
    if (input$depth_of_coverage == "30x") {
      load30x(summary)
      idx <- 4
      dta <- do.call(rbind, list(melt(AFR_30x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(AMR_30x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(EAS_30x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(EUR_30x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(SAS_30x[selCCDS, ], id.vars = c("CCDS", "Population", "GeneSymbol"))))
      
    }
  
    dta$CCDS <- as.character(dta$CCDS)
    dta$Population <- as.character(dta$Population)
    dta$variable <- as.character(dta$variable)
    dta$value <- as.numeric(dta$value)
    dta$GeneSymbol <- myValue$gene
    
    gAD <- gnomad_exome[selCCDS, idx]
    p1 <- ggplot(dta, aes(x=Population, y=value, color=Population)) + 
      theme_bw() + geom_violin()  +
      facet_wrap(GeneSymbol~CCDS) +
      geom_hline(yintercept=gAD, colour="black", show.legend = T) +
      scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
      ylab("Breadth of coverage") + 
      theme(legend.position="bottom", strip.text.x = element_text(size=12), axis.title=element_text(size=12))
<<<<<<< HEAD
    message("[PLOT] ", unique(dta$GeneSymbol), ", ", class(dta$GeneSymbol))
=======
    message("[PLOT] ", ccds2ens[ selCCDS, 2])

>>>>>>> db201e22460e7ca882a0716f96afeb7ae0da3be5
    #ggplotly(p1)
    p1
  }
}




# Run the app ----
shinyApp(ui = ui, server = server)