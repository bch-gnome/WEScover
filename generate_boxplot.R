library(RColorBrewer)
library(plotly)

AFR_30x <- read.delim('coverage_tables/GTR/AFR_30x.txt', header = T, sep = '\t', 
                      stringsAsFactors = F, row.names = 1)
AMR_30x <- read.delim('coverage_tables/GTR/AMR_30x.txt', header = T, sep = '\t', 
                      stringsAsFactors = F, row.names = 1)
EAS_30x <- read.delim('coverage_tables/GTR/EAS_30x.txt', header = T, sep = '\t', 
                      stringsAsFactors = F, row.names = 1)
EUR_30x <- read.delim('coverage_tables/GTR/EUR_30x.txt', header = T, sep = '\t', 
                      stringsAsFactors = F, row.names = 1)
SAS_30x <- read.delim('coverage_tables/GTR/SAS_30x.txt', header = T, sep = '\t', 
                      stringsAsFactors = F, row.names = 1)

saveRDS(AFR_30x, "data/AFR_30x.rds")
saveRDS(AMR_30x, "data/AMR_30x.rds")
saveRDS(EAS_30x, "data/EAS_30x.rds")
saveRDS(EUR_30x, "data/EUR_30x.rds")
saveRDS(SAS_30x, "data/SAS_30x.rds")

gnomad_exome <- read.delim('gnomad_exome.txt', header = T, sep = '\t', 
                           stringsAsFactors = F)
rownames(gnomad_exome) <- gnomad_exome[,1]
summary <- read.delim('GTR_summary.txt', header = T, sep = '\t', 
                           stringsAsFactors = F)
rownames(summary) <- summary[,1]

ccds_id <- rownames(subset(summary, gene =="A1BG"))
head(AFR_10x[ccds_id,])

# get ccds_ids
paste(rownames(subset(summary, gene == "A1CF")), collapse = " ")
mean(subset(summary, gene == "A1CF")[,3])
subset(summary, gene == "A1CF")

cov_20x <- list(AFR = as.numeric(AFR_20x["CCDS10.1",]),
     AMR = as.numeric(AMR_20x["CCDS10.1",]),
     EAS = as.numeric(EAS_20x["CCDS10.1",]),
     EUR = as.numeric(EUR_20x["CCDS10.1",]),
     SAS = as.numeric(SAS_20x["CCDS10.1",]))


plot_ly(type = 'box') %>%
  layout(yaxis = list(tickformat = "%")) %>%
  add_boxplot(y = cov_20x$AFR, boxpoints = 'all',
              marker = list(color = brewer.pal(6, "Set3")[1]),
              line = list(color = brewer.pal(6, "Set3")[1]),
              name = "AFR", jitter = 0.3, pointpos = -0.9) %>%
  add_boxplot(y = cov_20x$AMR, boxpoints = 'all',
              marker = list(color = brewer.pal(6, "Set3")[3]),
              line = list(color = brewer.pal(6, "Set3")[3]),
              name = "AMR", jitter = 0.3, pointpos = -0.9) %>%
  add_boxplot(y = cov_20x$EAS, boxpoints = 'all',
              marker = list(color = brewer.pal(6, "Set3")[4]),
              line = list(color = brewer.pal(6, "Set3")[4]),
              name = "EAS", jitter = 0.3, pointpos = -0.9) %>%
  add_boxplot(y = cov_20x$EUR, boxpoints = 'all',
              marker = list(color = brewer.pal(6, "Set3")[5]),
              line = list(color = brewer.pal(6, "Set3")[5]),
              name = "EUR", jitter = 0.3, pointpos = -0.9) %>%
  add_boxplot(y = cov_20x$SAS, boxpoints = 'all',
              marker = list(color = brewer.pal(6, "Set3")[6]),
              line = list(color = brewer.pal(6, "Set3")[6]),
              name = "SAS", jitter = 0.3, pointpos = -0.9) %>%
  add_segments(x = "AFR", xend = "SAS", y = gnomad_exome["CCDS10.1",3], 
               yend = gnomad_exome["CCDS10.1",3], name = "gnomAD", 
               line = list(color = "black")) %>%
  layout(title = "CCDS10.1 depth of coverage at 20x")