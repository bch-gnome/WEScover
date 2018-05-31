choices <- list("10x","20x","30x")
choice <- sample(choices)
data <- switch(choices, 
               "10x" = list(AFR = as.numeric(AFR_10x[input$gene_symbol,]),
                            AMR = as.numeric(AMR_10x[input$gene_symbol,]),
                            EAS = as.numeric(EAS_10x[input$gene_symbol,]),
                            EUR = as.numeric(EUR_10x[input$gene_symbol,]),
                            SAS = as.numeric(SAS_10x[input$gene_symbol,])),
               "20x" = list(AFR = as.numeric(AFR_20x[input$gene_symbol,]),
                            AMR = as.numeric(AMR_20x[input$gene_symbol,]),
                            EAS = as.numeric(EAS_20x[input$gene_symbol,]),
                            EUR = as.numeric(EUR_20x[input$gene_symbol,]),
                            SAS = as.numeric(SAS_20x[input$gene_symbol,])),
               "30x" = list(AFR = as.numeric(AFR_30x[input$gene_symbol,]),
                            AMR = as.numeric(AMR_30x[input$gene_symbol,]),
                            EAS = as.numeric(EAS_30x[input$gene_symbol,]),
                            EUR = as.numeric(EUR_30x[input$gene_symbol,]),
                            SAS = as.numeric(SAS_30x[input$gene_symbol,]))
)
