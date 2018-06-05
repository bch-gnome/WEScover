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
summary <- read.fst("data/summary.fst")
gtrM <- read.fst("data/gtrM.fst")
gtrS <- read.fst("data/gtrS.fst")
gpt <- read.fst("data/gpt.fst")
#genes_by_ccds_id <- read.fst("data/genes_by_ccds_id.fst")

# Define UI ----
ui <- fluidPage(
  useShinyjs(),
  #theme = shinytheme("cerulean"),
  #titlePanel("WEScover"),
  navbarPage("WEScover",
    tabPanel("Report",
      sidebarLayout(position = "left",
      sidebarPanel(
              tags$h2("User input"),
              selectizeInput("gene_symbol",
                              label = "Gene symbol",
                              choices = NULL,
                              multiple = TRUE),
              selectizeInput("gpt",
                              label="Test name",
                              choices = NULL,
                              multiple = TRUE),
               selectInput("depth_of_coverage", 
                           label = "Depth of coverage",
                           choices = c("10x", "20x", "30x"),
                           selected = "20x")
               ,
               actionButton("update", "Submit query"),
               actionButton("resetInputs", "Clear inputs")
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
  
  updateSelectizeInput(session, 'gene_symbol', choices = gene_symbol$gene_symbol, server = TRUE)
  updateSelectizeInput(session, 'gpt', choices =sort(unique(gpt$test_name)), server = TRUE)
  
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
    message(paste('click on detail ', selectedRow))
    geneS <- as.character(main_table()[selectedRow, 1])
    
    myValue$summary_table <<- 
      createMainTable2(geneS, input$depth_of_coverage, summary, gtrM, gtrS)[,c(1,2,6:10)]
    myValue$GPT_table <<- 
      createGPT(selectedRow, main_table(), gtrM)
    myValue$violon_population <<- createPlot(geneS, main_table()[selectedRow, 2])
    myValue$gene <<- geneS
    myValue$ccds <<- main_table()[selectedRow, 2]
    fileAD <- 
      paste0("/opt/shiny-server/samples/sample-apps/WEScover/coverage_plots/", myValue$ccds, ".png")
    if(file.exists(fileAD)) {
      myValue$gnomAD_plot <<- list(
        src = fileAD,
        contentType = 'image/png',
        width = "100%",
        #height = 300,
        alt = paste0("gnomAD coverage for ", myValue$ccds))
    } else {
      myValue$gnomAD_plot <<- list(
        src = "/opt/shiny-server/samples/sample-apps/WEScover/coverage_plots/NO_CCDS_GAD.png",
        contentType = 'image/png',
        width = "100%",
        #height = 300,
        alt = paste0("No coverage from gnomAD for ", myValue$ccds))
    }
    showModal(modal_main())
    
  })
  
  # create reactive main table
  main_table <- reactive({
    geneS <- c()
    if (length(input$gene_symbol) != 0) {
      geneS <- c(geneS, input$gene_symbol)
    }
    
    if (length(input$gpt) != 0) {
      geneP <- as.character(unique(gpt[gpt$test_name == input$gpt, "gene_symbol"]))
      geneS <- c(geneS, geneP)
    }
    
    geneS <- unique(geneS)
    tbl <- createMainTable2(geneS, input$depth_of_coverage, summary, gtrM, gtrS)
    if(ncol(tbl) > 0) {
      tbl <- tbl[ , seq(5)]
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
        tabPanel("Violin per population",
                 plotOutput("violin_population")),
        tabPanel("gnomAD coverage",
                 imageOutput("gnomAD_plot")),
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
    
    gAD <- gnomad_exome[selCCDS, idx]
    p1 <- ggplot(dta, aes(x=Population, y=value, color=Population)) + 
      theme_bw() + geom_violin(outlier.shape = NA)  +
      facet_wrap(GeneSymbol~CCDS) +
      geom_hline(yintercept=gAD, colour="black", show.legend = T) +
      scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
      ylab("Breadth of coverage")
    
    #ggplotly(p1)
    p1
  }
}

# Run the app ----
shinyApp(ui = ui, server = server)