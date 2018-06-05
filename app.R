# load packages
library(shiny)
library(shinyjs)
library(shinythemes)
library(DT)
library(plotly)
library(reshape2)
library(RColorBrewer)
#install.packages("shinyBS")
# laod functions
source("function.R")

# load data
gene_symbol <- readRDS("../data/gene_symbol.rds")
gnomad_exome <- read.delim("gnomad_exome.txt", header = F, stringsAsFactors = FALSE)
rownames(gnomad_exome) <- gnomad_exome$V1
summary <- readRDS("../data/summary.rds")
gtrM <- read.delim("GTR_table.txt")
gtrS <- read.delim("GTR_summary.txt")
rownames(gtrS) <- gtrS$ccds_id
pheno <- readRDS("../data/phenotypes.rds")
genes_by_ccds_id <- readRDS('../data/genes_by_ccds_id.rds')

#master <- readRDS("../data/master_table.rds")
#rownames(master) <- master$gene_symbol

# input <- list()
# input$gene_symbol <- c("A2ML1", "A2M") #c("ASTN1", "A1CF")

# Define UI ----
ui <- fluidPage(
  #theme = shinytheme("cerulean"),
  # useShinyjs(),
  titlePanel("WEScover"),
  fluidRow(
    column(2,
        wellPanel(
           h2("User input"),
           # radioButtons("type_input", 
           #              label = "Type of input", 
           #              choices = c("Gene symbol", "Phenotype"),
           #              selected = "Gene symbol",
           #              inline = TRUE),
           selectizeInput("gene_symbol", 
                          label = "Gene symbol",
                          choices = gene_symbol,
                          multiple = TRUE ),
           selectizeInput("phenotype",
                          label="Phenotype",
                          choices = sort(unique(pheno$phenotype_name)),
                          multiple = TRUE),
           selectInput("depth_of_coverage", 
                       label = "Depth of coverage",
                       choices = c("10x", "20x", "30x"),
                       selected = "20x"),
           #submitButton("Submit query")
           actionButton("update", "Submit query")
    )),
    column(10,
           dataTableOutput('tableMain')  
    )#,
    # column(5, 
    #        plotlyOutput("plot")#, height = 650)
    #        #plotOutput("plot")
    # )
  ),
  hr(),
  fluidRow(
    column(3),
    column(6),
    column(3)
  )
)


# Define server logic ----
server <- function(input, output) {
  
  # generator of buttons for the main table
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
  
  # list of global values
  myValue <- reactiveValues(summary_table = NA, GPT_table = NA, violon_population = NA, gene = "", ccds="")
  
  
  # observer to capture detail buttons in the main table
  observeEvent(input$detail_button, {
    selectedRow <- as.numeric(strsplit(input$detail_button, "_")[[1]][2])
    message(paste('click on detail ', selectedRow))
    geneS <- as.character(main_table()[selectedRow, 1])
    
    myValue$summary_table <<- 
      createMainTable2(geneS, input$depth_of_coverage, summary, gtrM, gtrS)[,c(1,2,6:10)]
    message("MODAL:",colnames(myValue$summary_table))
    myValue$GPT_table <<- 
      createGPT(selectedRow, main_table(), gtrM)
    message("MODAL:",colnames(myValue$GPT_table))
    myValue$violon_population <<- createPlot(geneS, main_table()[selectedRow, 2])
    myValue$gene <<- geneS
    myValue$ccds <<- main_table()[selectedRow, 2]
    message("FN")
    showModal(modal_main())
    
  })
  
  # create reactive main table
  main_table <- reactive({
    geneS <- c()
    if (length(input$gene_symbol) != 0) {
      geneS <- c(geneS, input$gene_symbol)
    }
    
    if (length(input$phenotype) != 0) {
      geneP <- as.character(unique(pheno[pheno$phenotype_name == input$phenotype, "gene_symbol"]))
      geneS <- c(geneS, geneP)
    }
    
    geneS <- unique(geneS)
    tbl <- createMainTable2(geneS, input$depth_of_coverage, summary, gtrM, gtrS)
    if(ncol(tbl) > 0) {
      tbl <- tbl[ , seq(5)]
      tbl$Action <- shinyInput(actionButton, nrow(tbl), 'button_', label = "Detail", onclick = 'Shiny.onInputChange(\"detail_button\",  this.id)' )
    }
    tbl
  })
  
  # display reactive main table
  output$tableMain <- renderDataTable( {
    input$update
    isolate({
      if(ncol(main_table()) > 0) {
        formatStyle(datatable(main_table(), escape = FALSE, selection = 'none'), 
                  columns = "Global mean", target = 'row', backgroundColor = 
                    styleInterval(cuts = c(95, 98), values=c("#FFE4E1", "#FFFFE0", "#F0FFF0")))
      }
    })
  })
  
  # modal window to show tables and plots
  modal_main <- function(failed = FALSE){
    modalDialog(
      tabsetPanel(type="tabs",
        tabPanel("Population summary",
          dataTableOutput('summary_table')),
        tabPanel("Gene panels",
                 dataTableOutput('GPT_table')),
        tabPanel("Violin per population",
                 plotOutput("violin_population"))
      ),
      easyClose = TRUE
    )
  }

  # table displayed in our modal
  output$summary_table <- renderDataTable({
    myValue$summary_table
  }, escape = FALSE)
  
  output$violin_population <- renderPlot({
    message("RENDER")
    #myValue$violin_population
    createPlot(myValue$gene, myValue$ccds)
  })
  
  output$GPT_table <- renderDataTable({
    message("TBL")
    myValue$GPT_table
  }, escape = FALSE)
  
  
  
  
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
    message("AB")
    p1 <<- ggplot(dta, aes(x=Population, y=value, color=Population)) + 
      theme_bw() + geom_violin(outlier.shape = NA)  +
      facet_wrap(GeneSymbol~CCDS) +
      geom_hline(yintercept=gAD, colour="black", show.legend = T) +
      scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
      ylab("Breadth of coverage")
    
    message("CD")
    #ggplotly(p1)
    p1
  }
}

# Run the app ----
shinyApp(ui = ui, server = server)