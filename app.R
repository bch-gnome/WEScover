# load packages

library(RColorBrewer)
library(plotly)
library(shiny)
library(reshape2)

source("function.R")

# get data

gene_symbol <- readRDS("../data/gene_symbol.rds")
gnomad_exome <- readRDS("../data/gnomad_exome.rds")
summary <- readRDS("../data/summary.rds")

# Define UI ----
ui <- fluidPage(
  titlePanel("WEScover"),
  sidebarLayout(
    sidebarPanel(
      h2("User input"),
      selectizeInput("gene_symbol", 
                  label = "Gene symbol",
                  choices = gene_symbol,
                  selected = gene_symbol[1],
                  multiple = TRUE ),
      textInput("phenotype", h3("Phenotype"), value = "Enter phenotype..."),
      selectInput("depth_of_coverage", 
            label = "Depth of coverage",
            choices = c("10x", "20x", "30x"),
            selected = "20x"),

      numericInput("breadth_of_coverage", 
                   h3("Maximum breadth of coverage"), 
                   value = 0.95)
      ),
    
    # create div for plot output
    mainPanel(
      plotlyOutput("plot", height = 650)
      )
    )
  )

# Define server logic ----
server <- function(input, output) {
  
  # call plotly function using inputs
  output$plot <- renderPlotly({

    idx <- 0 
    if (input$depth_of_coverage == "10x") {
      load10x(summary)
      idx <- 2
      
      dta <- do.call(rbind, list(melt(AFR_10x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population")),
                                 melt(AMR_10x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population")),
                                 melt(EAS_10x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population")),
                                 melt(EUR_10x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population")),
                                 melt(SAS_10x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population"))))
    }
    
    if (input$depth_of_coverage == "20x") {
      load20x(summary)
      idx <- 3
      dta <- do.call(rbind, list(melt(AFR_20x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population")),
                                 melt(AMR_20x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population")),
                                 melt(EAS_20x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population")),
                                 melt(EUR_20x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population")),
                                 melt(SAS_20x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population"))))

    }
    
    if (input$depth_of_coverage == "30x") {
      load30x(summary)
      idx <- 4
      dta <- do.call(rbind, list(melt(AFR_30x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population")),
                                 melt(AMR_30x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population")),
                                 melt(EAS_30x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population")),
                                 melt(EUR_30x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population")),
                                 melt(SAS_30x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population"))))

     }
    
    gAD <- gnomad_exome[getCCDS(input$gene_symbol, summary), idx]

    p1 <- ggplot(dta, aes(x=Population, y=value, color=Population)) + 
      theme_bw() + geom_boxplot(outlier.shape = NA)  + 
      facet_wrap(~CCDS, ncol = 2) + geom_jitter(alpha=0.55) +
      geom_hline(yintercept=gAD, colour="black", show.legend = T) + 
      scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
      ylab("Breadth of coverage")

 #   ggplotly(p1)
    
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)