# load packages
library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(shinyjs)
library(reshape2)
library(RColorBrewer)
library(fst)
library(data.table)
setDTthreads(18)

# laod functions
source("function.R")

# load data
gene_symbol <- read.fst("data/gene_symbol.fst")
gnomad_exome <- read.fst("data/gnomad_exome.fst")
row.names(gnomad_exome) <- gnomad_exome$V1
summary <- read.fst("data/summary.fst")
gtrM <- read.fst("data/gtrM.fst")
gtrS <- read.fst("data/gtrS.fst")
gpt_tP_tG <- read.fst("data/gpt_tP_tG.fst")
ccds2ens <- readRDS("data/ccds_ens_map.rds")

# format data
gpt_tP_tG$GTR_accession <- as.character(gpt_tP_tG$GTR_accession)
gpt_tP_tG$test_name <- as.character(gpt_tP_tG$test_name)
gpt_tP_tG$phenotype_name <- as.character(gpt_tP_tG$phenotype_name)

setDTthreads(18)

# Define UI ----
ui <- fluidPage(
  useShinyjs(),
  theme = shinytheme("flatly"),
  tags$head(tags$style(".modal-dialog{min-width:1200px}")),
  ## javascript required to clean the 'check' of detail buttons in main table
  tags$script("
      Shiny.addCustomMessageHandler('resetInputValue', function(variableName){
        Shiny.onInputChange(variableName, null);
      });
    "),
  tags$style(type="text/css", "body {padding-top: 80px;} .selectize-input {height: 45px;} .action-button {height:45px; width:100%;} .center {display: block; margin-left: auto; margin-right: auto}"),
  navbarPage("WEScover", id="mainNav", windowTitle = "WEScover", position = "fixed-top", fluid = TRUE,
    tabPanel("Home",
     absolutePanel( width = "70%", left = "15%", right = "15%",
       wellPanel(
         em(h1("WEScover")),
         hr(),
         p(em('WEScover'), 'helps users to check whether genes of interest could be sufficiently covered in terms of breadth and depth by whole exome sequencing (WES). For each transcript, breadth of coverage data was calculated at 10x, 20x, and 30x read depth from the ', 
           a("1000 Genomes Project (1KGP)", href = "http://www.internationalgenome.org/", target="_blank"), 
           '(N=2,692). A user will be able to minimize the chance of false negatives by selecting a targeted gene panel test for the genes that WES cannot cover well.'),
         p('Breadth and depth of coverage for ', a(em('NOTCH1'), href = "http://gnomad.broadinstitute.org/gene/ENSG00000148400", target="_blank"),
           ' are illustrated below. For some of the exons, breadth of coverage seems to be sub-optimal that could result in false negative results with WES.'),
         tags$img(src="gnomAD_notch1.png", alt = "Coverage from gnomAD project for NOTCH1", style="width:650px;height:300px", class="center"),
         p(em('WEScover'), ' provides detailed coverage information including difference in breadth of coverage between continent-level populatios.'),
         tags$img(src="violin_notch1.png", alt = "Contintental population breath of coverage violin plot for CCDS43905.1/NOTCH1", style="width:650px;height:300px", class="center"),
         p('Phenotype, genetic test names, or gene symbols can be used to retrieve coverage information in the query window. The output summary helps users to choose WES vs. targeted gene panel testing.')
       )
     )
    ),
    tabPanel("Query",
      sidebarLayout(position = "left",
      sidebarPanel(
              tags$h2("User input"),
              fluidRow(
                column(8,
                  selectizeInput("phen",
                    label="Phenotype",
                    choices = NULL,
                    multiple = TRUE)
                ),
                column(4,
                  HTML("<label>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</label>"),
                  tags$br(),
                  actionButton("fGPT", "Filter")#, icon = icon("filter", lib = "glyphicon"))
                )
              ),
              fluidRow(
                column(8,
                  selectizeInput("gpt",
                    label="GPT name",
                    choices = NULL,
                    multiple = TRUE)
                ),
                column(4,
                  HTML("<label>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</label>"),
                  tags$br(),
                  actionButton("fGenes", "Filter")#, icon = icon("filter", lib = "glyphicon"))
                )
              ),
              fluidRow(
                column(12, 
                  selectizeInput("gene_symbol", label = "Gene symbol",
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
                    )
                  )
                )
              ),
              fluidRow(
                column(12, 
                 selectInput("depth_of_coverage",
                             label = "Depth of coverage",
                             choices = c("10x", "20x", "30x"),
                             selected = "20x")
                             #width = '70%'),
                )
              ),
              fluidRow(
                column(2),
                column(5, actionButton("clear", "Clear inputs", class = "btn-secondary")),
                column(5, actionButton("update", "Submit query", class = "btn-primary"))
              )
        ),
        mainPanel(
          dataTableOutput('tableMain')
        )
      )
    ),
    tabPanel("Help",
      includeHTML("help.html")
    ),
    tabPanel("Data",
      includeHTML("data.html"),
      hidden( 
        numericInput( inputId = 'refresh_helper', label = 'refresh_helper', value = 0 ) 
      )
    )
  )
)


# Define server logic ----
server <- function(input, output, session) {
  # fill with default values
  updateSelectizeInput(session, 'phen', choices = sort(unique(gpt_tP_tG$phenotype_name)), server = TRUE)
  updateSelectizeInput(session, 'gene_symbol', choices = gene_symbol$gene_symbol, server = TRUE)
  updateSelectizeInput(session, 'gpt', choices = sort(unique(gpt_tP_tG$test_name)), server = TRUE)
  
  # if a gene symbol is provided by url
  observe({
    query <- parseQueryString(session$clientData$url_search)
    message("query")
    if (!is.null(query[['gene']]) & length(input$gene_symbol) == 0) {
      geneS <- strsplit(query[['gene']], ",")[[1]]
      updateNavbarPage(session, "mainNav", "Query")
      updateSelectizeInput(session, 'gene_symbol', choices = gene_symbol$gene_symbol, selected = geneS, server = TRUE)
    #}
    #if(length(input$gene_symbol) != 0 & input$refresh_helper == 0) {
      message("update! ", input$refresh_helper, " ", myValue$inc)
      updateNumericInput( session = session, inputId = 'refresh_helper', value = input$refresh_helper + 1 )
    }
  })
  
  # if clear button is pushed
  observeEvent (input$clear,{
    updateSelectizeInput(session, 'phen', choices = sort(unique(gpt_tP_tG$phenotype_name)), server = TRUE)
    updateSelectizeInput(session, 'gene_symbol', choices = gene_symbol$gene_symbol, server = TRUE, label = "Gene symbol")
    updateSelectizeInput(session, 'gpt', choices = sort(unique(gpt_tP_tG$test_name)), server = TRUE, label = "GPT name")
    updateSelectInput(session, "depth_of_coverage", choices = c("10x", "20x", "30x"), selected = "20x")
  })
  
  # if filter GPT button is pushed
  observeEvent (input$fGPT,{
    listPhe <- unique(
      gpt_tP_tG$test_name[ gpt_tP_tG$GTR_accession %in% gpt_tP_tG$GTR_accession[ gpt_tP_tG$phenotype_name %in% as.character(input$phen) ] ])
    updateSelectizeInput(session, 'gpt', choices = sort(listPhe), server = TRUE, label = "GPT name (filtered)")
    listGenes <- unique(
      gpt_tP_tG$gene_symbol[ gpt_tP_tG$GTR_accession %in% gpt_tP_tG$GTR_accession[ gpt_tP_tG$phenotype_name %in% as.character(input$phen) ]])
    updateSelectizeInput(session, 'gene_symbol', choices = sort(listGenes), server = TRUE, label = "Gene symbol (filtered)")
  })
  
  # if filter gene symbol butto is pushed
  observeEvent (input$fGenes,{
    listGenes <- gpt_tP_tG$gene_symbol[ gpt_tP_tG$test_name %in% as.character(input$gpt) ]
    updateSelectizeInput(session, 'gene_symbol', choices = sort(listGenes), server = TRUE, label = "Gene symbol (filtered)")
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
      geneS <- as.character(main_table()[selectedRow, 1])
      ccds <-  as.character(main_table()[selectedRow, 2])
      withProgress(message = paste0('Query for ', ccds, '/', geneS), value = 0.1, {
      
      incProgress(0.2, detail = "(Obtaining continental population)")
      myValue$summary_table <<- 
        createMainTable2(geneS, input$depth_of_coverage, summary, gtrM, gtrS)[,c(1,2,6:10)]
      
      incProgress(0.2, detail = "(Obtaining gene panel tests)")
      myValue$GPT_table <<- createGPT(selectedRow, main_table(), gtrM)
      
      incProgress(0.2, detail = "(Creating violin plot)")
      myValue$violon_population <<- createPlot(geneS, main_table()[selectedRow, 2])
      
      incProgress(0.2, detail = "(Saving details)")
      myValue$gene <<- geneS
      myValue$ccds <<- ccds
      fileAD <- paste0("coverage_plots/", ccds2ens[as.character(myValue$ccds), 2], ".png")
      
      if(file.exists(fileAD)) {
        myValue$gnomAD_plot <<- list(
          src = fileAD,
          contentType = 'image/png',
          width = "100%",
          height = 300,
          alt = paste0("gnomAD coverage for ", ccds2ens[as.character(myValue$ccds), 2]))
      } else {
        myValue$gnomAD_plot <<- list(
          src = "coverage_plots/Dummy_coverage_plot.png",
          contentType = 'image/png',
          width = "100%",
          height = 300,
          alt = paste0("No coverage from gnomAD for ", ccds2ens[as.character(myValue$ccds), 2]))
      }
      
      setProgress(1)
    })
    showModal(modal_main())
    # reset the "check" on the button in the main table
    session$sendCustomMessage(type = 'resetInputValue', message =  "detail_button")
  })
  
  # create reactive main table
  main_table <- reactive ({
    geneS <- c()
    if(length(input$phen) != 0) {
      if (length(input$gene_symbol) != 0) {
        geneS <- c(geneS, input$gene_symbol)
      } else if (length(input$gpt) != 0) {
        geneG <- as.character(unique(gpt_tP_tG[gpt_tP_tG$test_name %in% input$gpt, "gene_symbol"]))
        geneS <- c(geneS, geneG)
      } else if(length(input$phen) != 0) {
        geneP <- unique(gpt_tP_tG$gene_symbol[ 
          gpt_tP_tG$GTR_accession %in% gpt_tP_tG$GTR_accession[ gpt_tP_tG$phenotype_name %in% as.character(input$phen) ]])
        message(length(geneP))
        geneS <- c(geneS, geneP)
      }
    }

    if (length(input$gpt) != 0) {
      geneG <- as.character(unique(gpt_tP_tG[gpt_tP_tG$test_name %in% input$gpt, "gene_symbol"]))
      geneS <- c(geneS, geneG)
    }
    
    if (length(input$gene_symbol) != 0) {
      geneS <- c(geneS, input$gene_symbol)
    }
    
    geneS <- unique(geneS)
    tbl <- createMainTable2(geneS, input$depth_of_coverage, summary, gtrM, gtrS)
    if(ncol(tbl) > 0) {
      tbl <- tbl[ , c(1:5, 11:13)]
      tbl$Action <- shinyInput(actionButton, nrow(tbl), 'button_', 
                               label = "Detail", 
                               onclick = 'Shiny.onInputChange(\"detail_button\",  this.id)')
    }
    tbl
  })
  
  
  
  observeEvent(input$update, {
    updateNumericInput( session = session, inputId = 'refresh_helper', value = input$refresh_helper + 1 )
  })
  
  # display reactive main table
  output$tableMain <- renderDataTable( {
    #input$update
    message(input$refresh_helper)
    t = input$refresh_helper
    if(t > 0 ) {
      isolate({
        if(ncol(main_table()) > 0) {
          formatStyle(datatable(main_table(), escape = FALSE, selection = 'none'), 
                    columns = "Global coverage (mean, %)", target = 'row', backgroundColor = 
                      styleInterval(cuts = c(95, 98), values=c("#FFE4E1", "#FFFFE0", "#F0FFF0")))
        }
      })
    } else {
      data.frame()
    }
  })
  
  # modal window to show tables and plots
  modal_main <- function(failed = FALSE){
    modalDialog(size="l",
      tabsetPanel(type="tabs",
        tabPanel("Population summary",
          dataTableOutput('summary_table')),
        tabPanel("Coverage plots",
          fluidRow(
            column(4, align="center", plotOutput("violin_population")),
            column(8, align="center",
              tags$label("Click the plot to go to the gnomAD server."),
              tags$a(imageOutput("gnomAD_plot"),
                     href=paste0("http://gnomad.broadinstitute.org/gene/", ccds2ens[as.character(myValue$ccds), 2]), target="_blank")
          ))
        ),
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
    dta$GeneSymbol <- gene#myValue$gene
    
    gAD <- gnomad_exome[selCCDS, idx]
    p1 <- ggplot(dta, aes(x=Population, y=value, color=Population)) + 
      theme_bw() + geom_violin()  +
      facet_wrap(GeneSymbol~CCDS) +
      geom_hline(yintercept=gAD, colour="black", show.legend = T) +
      scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
      ylab("Breadth of coverage") + 
      theme(legend.position="bottom", strip.text.x = element_text(size=12), axis.title=element_text(size=12))
    p1
  }
}




# Run the app ----
shinyApp(ui = ui, server = server)