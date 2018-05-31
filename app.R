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

# load data
gene_symbol <- readRDS("../data/gene_symbol.rds")
gnomad_exome <- read.delim("gnomad_exome.txt", header = F, stringsAsFactors = FALSE)
rownames(gnomad_exome) <- gnomad_exome$V1
summary <- readRDS("../data/summary.rds")
gtrM <- read.delim("GTR_table.txt")
gtrS <- read.delim("GTR_summary.txt")
rownames(gtrS) <- gtrS$ccds_id
pheno <- read.delim("test_condition_gene.txt", stringsAsFactors = FALSE)
pheno <- pheno[pheno$object == "gene", ]
pheno$condition <- sapply(strsplit(pheno$object_name, ":"), function(x) {
  if(length(x) == 1) {
    ""
  } else if(length(x) == 2) {
    x[2]
  } else {
    paste0(x[2], ":", x[3])
  }
})

# input <- list()
# input$gene_symbol <- c("A2ML1", "A2M") #c("ASTN1", "A1CF")

# Define UI ----
ui <- fluidPage(
  #theme = shinytheme("cerulean"),
  useShinyjs(),
  titlePanel("WEScover"),
  fluidRow(
    column(3,
           h2("User input"),
           radioButtons("type_input", 
                        label = "Type of input", 
                        choices = c("Gene symbol", "Phenotype"),
                        selected = "Gene symbol",
                        inline = TRUE),
           selectizeInput("gene_symbol", 
                          label = "Gene symbol",
                          choices = gene_symbol,
                          selected = gene_symbol[1],
                          multiple = TRUE ),
           selectizeInput("phenotype",
                          label="Phenotype",
                          choices = unique(pheno$condition),
                          multiple = FALSE),
           #textInput("phenotype", h3("Phenotype"), value = "Enter phenotype..."),
           selectInput("depth_of_coverage", 
                       label = "Depth of coverage",
                       choices = c("10x", "20x", "30x"),
                       selected = "20x"),
           
           numericInput("breadth_of_coverage", 
                        h3("Maximum breadth of coverage"), 
                        value = 0.95)
    ),
    hr(),
    column(6,
           dataTableOutput('tableMain')  
    ),
    column(3, 
           plotlyOutput("plot", height = 650)
           #plotOutput("plot", height = 650)
    )
  ),
  fluidRow(
    column(3),
    column(6),
    column(3)
  )
)


# Define server logic ----
server <- function(input, output) {
  # disable/enable input according to radio buttons
  observe({
    if (input$type_input == "Gene symbol") {
      shinyjs::enable("gene_symbol")
      shinyjs::disable("phenotype")
    } else {
      shinyjs::disable("gene_symbol")
      shinyjs::enable("phenotype")
    }
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
  myValue <- reactiveValues(modal_table = NA, selected_gene = "")
  
  
  # observer to capture detail buttons in the main table
  observeEvent(input$detail_button, {
    selectedRow <- as.numeric(strsplit(input$detail_button, "_")[[1]][2])
    message(paste('click on detail ', selectedRow))
   
    myValue$modal_table <<- createGPT(selectedRow, main_table(), gtrM)
    showModal(modal_main(1))
  })
  
  # observer to capture GPT buttons in the main table
  observeEvent(input$population_button, {
    selectedRow <- as.numeric(strsplit(input$population_button, "_")[[1]][2])
    message(paste('click on population ', selectedRow))
    geneS <- as.character(main_table()[selectedRow, 1])
    message('\t', geneS)
    
    myValue$modal_table <<- createMainTable(geneS, input$depth_of_coverage, summary, gtrM, gtrS)[ , c(1:2, 5:9)]
    showModal(modal_main(2))
  })
  
  # observer to capture violin buttons in the main table
  observeEvent(input$violin_button, {
    selectedRow <- as.numeric(strsplit(input$violin_button, "_")[[1]][2])
    message(paste('click on vilin ', selectedRow))
    xx <- paste0(main_table()[selectedRow, 1], "-", main_table()[selectedRow, 2])
    myValue$violin <<- xx
  })
  
  
  # create reactive main table
  main_table <- reactive({
    if(input$type_input == "Gene symbol") {
      geneS <- input$gene_symbol
    } else {
      geneS <- as.character(unique(pheno[pheno$condition == input$phenotype, "gene_symbol"]))
    }
    tbl <- createMainTable(geneS, input$depth_of_coverage, summary, gtrM, gtrS)
    tbl <- tbl[ , seq(4)]
    tbl$B1 <- shinyInput(actionButton, nrow(tbl), 'button_', label = "See detail", onclick = 'Shiny.onInputChange(\"detail_button\",  this.id)' )
    tbl$B2 <- shinyInput(actionButton, nrow(tbl), 'button_', label = "See GTP", onclick = 'Shiny.onInputChange(\"population_button\",  this.id)' )
    tbl$B3 <- shinyInput(actionButton, nrow(tbl), 'button_', label = "See as violin", onclick = 'Shiny.onInputChange(\"violin_button\",  this.id)' )
    tbl
  })
  # display reactive main table
  output$tableMain <- renderDataTable(main_table(), server = FALSE, escape = FALSE, selection = 'none')
  
  
  # # rective table for list of GTP
  # gpt_reactive <- reactive({
  #   createGPT(selectedRow, main_table(), gtrM)
  # })
  
  # modal dialog box
  modal_main <- function(type, failed = FALSE){
    modalDialog(
      dataTableOutput('modal_table'),
      easyClose = TRUE
    )
  }
  
  # table displayed in our modal
  output$modal_table <- renderDataTable({
    myValue$modal_table
  }, escape = FALSE)

  
  # #event to trigger the modal box to appear
  # observeEvent(input$tableMain_rows_selected,{
  #   showModal(modal_main())
  # })
  
  # observeEvent(input$gene_symbol, {
  #   if(input$type_input == "Gene symbol") {
  #     geneS <- input$gene_symbol
  #   } else {
  #     geneS <- as.character(unique(pheno[pheno$condition == input$phenotype, "gene_symbol"]))
  #   }
  #   selCCDS <- getCCDS(geneS, summary)
  #   if(length(selCCDS) <= 9) {
  #     show("plot")
  #     hide("tableMain")
  #   } else {
  #     hide("plot") 
  #     show("tableMain")
  #   }
  # })
  
  # create the main plot if less than 9 CCDS are selected
  output$plot <- renderPlotly({ #renderPlot({
    if(input$type_input == "Gene symbol") {
      geneS <- input$gene_symbol
    } else {
      message("hi!")
      message(input$phenotype)
      geneS <- as.character(unique(pheno[pheno$condition == input$phenotype, "gene_symbol"]))
    }
    selCCDS <- getCCDS(geneS, summary)
    if(length(selCCDS) <= 9) {
      message("[PLOT] Number of CCDS: ", length(selCCDS), "; Number of genes:", length(input$gene_symbol) )
      
      idx <- 0 
      if (input$depth_of_coverage == "10x") {
        load10x(summary)
        idx <- 2
         
        dta <- do.call(rbind, list(melt(AFR_10x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                   melt(AMR_10x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                   melt(EAS_10x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                   melt(EUR_10x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                   melt(SAS_10x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol"))))
      }
      
      if (input$depth_of_coverage == "20x") {
        load20x(summary)
        idx <- 3
        dta <- do.call(rbind, list(melt(AFR_20x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                   melt(AMR_20x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                   melt(EAS_20x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                   melt(EUR_20x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                   melt(SAS_20x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol"))))
  
      }
      
      if (input$depth_of_coverage == "30x") {
        load30x(summary)
        idx <- 4
        dta <- do.call(rbind, list(melt(AFR_30x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                   melt(AMR_30x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                   melt(EAS_30x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                   melt(EUR_30x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                   melt(SAS_30x[getCCDS(geneS, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol"))))
  
       }
      
      
      dta$CCDS <- as.character(dta$CCDS)
      dta$Population <- as.character(dta$Population)
      dta$variable <- as.character(dta$variable)
      dta$value <- as.numeric(dta$value)
      
      # gAD <- gnomad_exome[getCCDS(geneS, summary), idx]
      dta2 <- rbind(dta, data.frame(
        CCDS = getCCDS(geneS, summary), 
        Population = "gAD", 
        GeneSymbol = unlist(lapply(geneS, function(x) { rep(x, length(getCCDS(x, summary))) })),
        variable = 1, 
        value = unlist(lapply(geneS, function(x) { gnomad_exome[getCCDS(x, summary), idx] }))
      ))
      
      dta2$Population <- factor(dta2$Population, levels = c("AFR", "AMR", "EAS", "EUR", "SAS", "gAD"))
      colnames(dta2) <- c("CCDS", "Population", "GeneSymbol", "variable", "Coverage")
      colorsPopulation <- c("AFR" = "#191970", #MidnightBlue
                       "AMR" = "#4169E1", #RoyalBlue
                       "EAS" = "#008B8B", #DarkCyan
                       "EUR" = "#228B22", #ForestGreen
                       "SAS" = "#9ACD32", #YellowGreen
                       "gAD" = "#FF4500"  #OrangeRed
                       )
      colorsPopulation <- c("AFR" = "#FA8072", #LightSalmon
                            "AMR" = "#90EE90", #LightGreen
                            "EAS" = "#FFA500", #Orange
                            "EUR" = "#AFEEEE", #PaleTurquoise
                            "SAS" = "#EE82EE", #Violet
                            "gAD" = "#000000"  #Black
      )
      
      
      nc <- ifelse(length(unique(dta2$CCDS)) %% 3 == 0, 3, 2)
      
      plot_ly(dta2, x=~Population, y=~Coverage, split=~CCDS, color=~Population, type = 'violin',
              box = list(visible = T), meanline = list(visible = T))
      
      # p1 <- ggplot(dta2, aes(x=Population, y=coverage, color=Population)) + 
      #   theme_bw() + geom_boxplot(outlier.shape = NA)  + 
      #   facet_wrap(GeneSymbol~CCDS, ncol = nc) + geom_jitter(alpha=0.55) +
      #   #geom_hline(yintercept=gAD, colour="black", show.legend = T) + 
      #   scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
      #   scale_color_manual(values=colorsPopulation) +
      #   ylab("Breadth of coverage")
      # 
      # #ggplotly(p1)
      # p1
    }
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)