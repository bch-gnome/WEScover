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

# input <- list()
# input$gene_symbol <- c("ASTN1", "A1CF")

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
      #plotOutput("plot", height = 650)
      )
    )
  )

# Define server logic ----
server <- function(input, output) {
  
  # call plotly function using inputs
  output$plot <-  renderPlotly({ #renderPlot({

    idx <- 0 
    if (input$depth_of_coverage == "10x") {
      load10x(summary)
      idx <- 2
       
      dta <- do.call(rbind, list(melt(AFR_10x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(AMR_10x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(EAS_10x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(EUR_10x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(SAS_10x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol"))))
    }
    
    if (input$depth_of_coverage == "20x") {
      load20x(summary)
      idx <- 3
      dta <- do.call(rbind, list(melt(AFR_20x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(AMR_20x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(EAS_20x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(EUR_20x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(SAS_20x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol"))))

    }
    
    if (input$depth_of_coverage == "30x") {
      load30x(summary)
      idx <- 4
      dta <- do.call(rbind, list(melt(AFR_30x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(AMR_30x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(EAS_30x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(EUR_30x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol")),
                                 melt(SAS_30x[getCCDS(input$gene_symbol, summary), ], id.vars = c("CCDS", "Population", "GeneSymbol"))))

     }
    
    gAD <- gnomad_exome[getCCDS(input$gene_symbol, summary), idx]
    dta$CCDS <- as.character(dta$CCDS)
    dta$Population <- as.character(dta$Population)
    dta$variable <- as.character(dta$variable)
    dta$value <- as.numeric(dta$value)
    
    dta2 <- rbind(dta, data.frame(CCDS = getCCDS(input$gene_symbol, summary), 
      Population = "gAD", 
      GeneSymbol = unlist(lapply(input$gene_symbol, function(x) { rep(x, length(getCCDS(x, summary))) })),
      variable = 1, 
      value = gAD))
    
    dta2$Population <- factor(dta2$Population, levels = c("AFR", "AMR", "EAS", "EUR", "SAS", "gAD"))
    colnames(dta2) <- c("CCDS", "Population", "GeneSymbol", "variable", "coverage")
    colorsPopulation <- c("AFR" = "#191970", #MidnightBlue
                     "AMR" = "#4169E1", #RoyalBlue
                     "EAS" = "#008B8B", #DarkCyan
                     "EUR" = "#228B22", #ForestGreen
                     "SAS" = "#9ACD32", #YellowGreen
                     "gAD" = "#FF4500"  #OrangeRed
                     )
    
    
    nc <- ifelse(length(unique(dta2$CCDS)) %% 3, 3, 2)
    
    p1 <- ggplot(dta2, aes(x=Population, y=coverage, color=Population)) + 
      theme_bw() + geom_boxplot(outlier.shape = NA)  + 
      facet_wrap(GeneSymbol~CCDS, ncol = nc) + geom_jitter(alpha=0.55) +
      #geom_hline(yintercept=gAD, colour="black", show.legend = T) + 
      scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
      scale_fill_manual(values=colorsPopulation) +
      ylab("Breadth of coverage")

    ggplotly(p1)
    #p1
    
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)