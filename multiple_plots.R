# load packages

library(RColorBrewer)
library(plotly)
library(shiny)

# get data

AFR_20x <- readRDS("data/AFR_20x.rds")
AMR_20x <- readRDS("data/AMR_20x.rds")
EAS_20x <- readRDS("data/EAS_20x.rds")
EUR_20x <- readRDS("data/EUR_20x.rds")
SAS_20x <- readRDS("data/SAS_20x.rds")

gene_symbol <- readRDS("data/gene_symbol.rds")
gnomad_exome <- readRDS("data/gnomad_exome.rds")
summary <- readRDS("data/summary.rds")

ui <- shinyUI(pageWithSidebar(
  
  headerPanel("Dynamic number of plots"),
  
  sidebarPanel(
    selectInput("gene_symbol", 
                label = "Gene symbol",
                choices = gene_symbol,
                selected = gene_symbol[1])
  ),
  
  mainPanel(
    # This is the dynamic UI for the plots
    plotlyOutput("plot")
  )
))

server <- function(input, output) {
  
  ccds_ids <- reactiveValues()
  
  observeEvent(input$gene_symbol, {
    ccds_ids$ccds_ids <- rownames(subset(summary, gene == input$gene_symbol))     # if the add button is clicked, increment the value by 1 and update it
  })

  # call plotly function using inputs
  output$plot <- renderPlotly({
    
    data <- list(AFR = NULL, AMR = NULL, EAS = NULL, EUR = NULL, SAS = NULL, gnomAD = NULL, ccds_id = NULL)
    for (i in 1:length(ccds_ids$ccds_ids)) {
        data[[1]] <- c(data[[1]],as.numeric(unlist(AFR_20x[c(ccds_ids$ccds_ids[i]),])))
        data[[2]] <- c(data[[2]],as.numeric(unlist(AMR_20x[c(ccds_ids$ccds_ids[i]),])))
        data[[3]] <- c(data[[3]],as.numeric(unlist(EAS_20x[c(ccds_ids$ccds_ids[i]),])))
        data[[4]] <- c(data[[4]],as.numeric(unlist(EUR_20x[c(ccds_ids$ccds_ids[i]),])))
        data[[5]] <- c(data[[5]],as.numeric(unlist(SAS_20x[c(ccds_ids$ccds_ids[i]),])))
        data[[6]] <- c(data[[6]],unlist(gnomad_exome[c(ccds_ids$ccds_ids[i]),3]))
        data[[7]] <- c(data[[7]],rep(ccds_ids$ccds_ids[i], time = 2692))
    }
    
    gAD <- data$gnomAD
    dta <- data.frame(
      variable = c(
        rep("AFR", time = length(data$AFR)),
        rep("AMR", time = length(data$AMR)),
        rep("EAS", time = length(data$EAS)),
        rep("EUR", time = length(data$EUR)),
        rep("SAS", time = length(data$SAS))
      ),
      value = c(data$AFR, data$AMR, data$EAS, data$EUR, data$SAS)
    )
    p1 <- ggplot(dta, aes(x=variable, y=value, colour=variable)) + geom_boxplot()  + 
      theme_bw() + facet_wrap(~unique(data$ccds_id), ncol = 1) + geom_jitter() +
      geom_hline(yintercept=gAD, colour="black", show.legend = T)
    
    ggplotly(p1)
    
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)