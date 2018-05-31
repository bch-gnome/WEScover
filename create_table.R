GTR <- read.delim("../GTR_table.txt", header = T, stringsAsFactors = F)
gene_symbol <- readRDS("data/gene_symbol.rds")
gnomad_exome <- readRDS("data/gnomad_exome.rds")

summary <- readRDS("data/summary.rds")

gene_symbol <- unique(summary$gene)
ccds_ids <- c()
global_mean <- c()
GTR_tests <- c()

for (i in 1:length(gene_symbol)) {
  ccds_ids[i] <- paste(subset(summary, gene == gene_symbol[i])$ccds_id, collapse = " ")
  global_mean[i] <- mean(subset(summary, gene == gene_symbol[i])$GTR_average_10x)
  GTR_tests[i] <- paste(subset(GTR, GeneSymbol == gene_symbol[i])$AccessionVersion, collapse = " ")
}

table <- data.frame(gene_symbol, ccds_ids, global_mean, GTR_tests)
