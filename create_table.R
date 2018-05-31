GTR <- read.delim("../../GTR_table.txt", header = T, stringsAsFactors = F)
GTR_10x <- read.delim("../../coverage_tables/GTR_10x.txt", header = T, stringsAsFactors = F)
GTR_20x <- read.delim("../../coverage_tables/GTR_20x.txt", header = T, stringsAsFactors = F)
GTR_30x <- read.delim("../../coverage_tables/GTR_30x.txt", header = T, stringsAsFactors = F)

source("function.R")

gene_symbol <- readRDS("data/gene_symbol.rds")
gnomad_exome <- readRDS("data/gnomad_exome.rds")

summary <- readRDS("../data/summary.rds")
x <- as.data.frame(table(summary$gene), stringsAsFactors = F)
y <- as.data.frame(table(GTR$GeneSymbol), stringsAsFactors = F)
y <- subset(y, Var1 %in% x$Var1)
x <- x[order(x[1,])]
y <- y[order(y[1,])]

load10x(summary)
load20x(summary)
load30x(summary)

nums <- unlist(lapply(AFR_10x, is.numeric))  
AFR_10x <- AFR_10x[,nums]
nums <- unlist(lapply(AFR_20x, is.numeric))  
AFR_20x <- AFR_20x[,nums]
nums <- unlist(lapply(AFR_30x, is.numeric))  
AFR_30x <- AFR_30x[,nums]
nums <- unlist(lapply(AMR_10x, is.numeric))  
AMR_10x <- AMR_10x[,nums]
nums <- unlist(lapply(AMR_20x, is.numeric))  
AMR_20x <- AMR_20x[,nums]
nums <- unlist(lapply(AMR_30x, is.numeric))  
AMR_30x <- AMR_30x[,nums]
nums <- unlist(lapply(EAS_10x, is.numeric))  
EAS_10x <- EAS_10x[,nums]
nums <- unlist(lapply(EAS_20x, is.numeric))  
EAS_20x <- EAS_20x[,nums]
nums <- unlist(lapply(EAS_30x, is.numeric))  
EAS_30x <- EAS_30x[,nums]
nums <- unlist(lapply(EUR_10x, is.numeric))  
EUR_10x <- EUR_10x[,nums]
nums <- unlist(lapply(EUR_20x, is.numeric))  
EUR_20x <- EUR_20x[,nums]
nums <- unlist(lapply(EUR_30x, is.numeric))  
EUR_30x <- EUR_30x[,nums]
nums <- unlist(lapply(SAS_10x, is.numeric))  
SAS_10x <- SAS_10x[,nums]
nums <- unlist(lapply(SAS_20x, is.numeric))  
SAS_20x <- SAS_20x[,nums]
nums <- unlist(lapply(SAS_30x, is.numeric))  
SAS_30x <- SAS_30x[,nums]

global_mean_10x <- rowMeans(GTR_10x)
global_mean_20x <- rowMeans(GTR_20x)
global_mean_30x <- rowMeans(GTR_30x)
AFR_mean_10x <- rowMeans(AFR_10x)
AFR_mean_20x <- rowMeans(AFR_20x)
AFR_mean_30x <- rowMeans(AFR_30x)
AMR_mean_10x <- rowMeans(AMR_10x)
AMR_mean_20x <- rowMeans(AMR_20x)
AMR_mean_30x <- rowMeans(AMR_30x)
EAS_mean_10x <- rowMeans(EAS_10x)
EAS_mean_20x <- rowMeans(EAS_20x)
EAS_mean_30x <- rowMeans(EAS_30x)
EUR_mean_10x <- rowMeans(EUR_10x)
EUR_mean_20x <- rowMeans(EUR_20x)
EUR_mean_30x <- rowMeans(EUR_30x)
SAS_mean_10x <- rowMeans(SAS_10x)
SAS_mean_20x <- rowMeans(SAS_20x)
SAS_mean_30x <- rowMeans(SAS_30x)

global_min_10x <- apply(GTR_10x, 1, FUN=min)
global_min_20x <- apply(GTR_20x, 1, FUN=min)
global_min_30x <- apply(GTR_30x, 1, FUN=min)
AFR_min_10x <- apply(AFR_10x, 1, FUN=min)
AFR_min_20x <- apply(AFR_20x, 1, FUN=min)
AFR_min_30x <- apply(AFR_30x, 1, FUN=min)
AMR_min_10x <- apply(AMR_10x, 1, FUN=min)
AMR_min_20x <- apply(AMR_20x, 1, FUN=min)
AMR_min_30x <- apply(AMR_30x, 1, FUN=min)
EAS_min_10x <- apply(EAS_10x, 1, FUN=min)
EAS_min_20x <- apply(EAS_20x, 1, FUN=min)
EAS_min_30x <- apply(EAS_30x, 1, FUN=min)
EUR_min_10x <- apply(EUR_10x, 1, FUN=min)
EUR_min_20x <- apply(EUR_20x, 1, FUN=min)
EUR_min_30x <- apply(EUR_30x, 1, FUN=min)
SAS_min_10x <- apply(SAS_10x, 1, FUN=min)
SAS_min_20x <- apply(SAS_20x, 1, FUN=min)
SAS_min_30x <- apply(SAS_30x, 1, FUN=min)

global_max_10x <- apply(GTR_10x, 1, FUN=max)
global_max_20x <- apply(GTR_20x, 1, FUN=max)
global_max_30x <- apply(GTR_30x, 1, FUN=max)
AFR_max_10x <- apply(AFR_10x, 1, FUN=max)
AFR_max_20x <- apply(AFR_20x, 1, FUN=max)
AFR_max_30x <- apply(AFR_30x, 1, FUN=max)
AMR_max_10x <- apply(AMR_10x, 1, FUN=max)
AMR_max_20x <- apply(AMR_20x, 1, FUN=max)
AMR_max_30x <- apply(AMR_30x, 1, FUN=max)
EAS_max_10x <- apply(EAS_10x, 1, FUN=max)
EAS_max_20x <- apply(EAS_20x, 1, FUN=max)
EAS_max_30x <- apply(EAS_30x, 1, FUN=max)
EUR_max_10x <- apply(EUR_10x, 1, FUN=max)
EUR_max_20x <- apply(EUR_20x, 1, FUN=max)
EUR_max_30x <- apply(EUR_30x, 1, FUN=max)
SAS_max_10x <- apply(SAS_10x, 1, FUN=max)
SAS_max_20x <- apply(SAS_20x, 1, FUN=max)
SAS_max_30x <- apply(SAS_30x, 1, FUN=max)

master <- data.frame(summary$ccds_id, summary$gene,
                     global_mean_10x, global_min_10x, global_max_10x,
                     AFR_mean_10x, AFR_min_10x, AFR_max_10x,
                     AMR_mean_10x, AMR_min_10x, AMR_max_10x,
                     EAS_mean_10x, EAS_min_10x, EAS_max_10x,
                     EUR_mean_10x, EUR_min_10x, EUR_max_10x,
                     SAS_mean_10x, SAS_min_10x, SAS_max_10x,
                     global_mean_20x, global_min_20x, global_max_20x,
                     AFR_mean_20x, AFR_min_20x, AFR_max_20x,
                     AMR_mean_20x, AMR_min_20x, AMR_max_20x,
                     EAS_mean_20x, EAS_min_20x, EAS_max_20x,
                     EUR_mean_20x, EUR_min_20x, EUR_max_20x,
                     SAS_mean_20x, SAS_min_20x, SAS_max_20x,
                     global_mean_30x, global_min_30x, global_max_30x,
                     AFR_mean_30x, AFR_min_30x, AFR_max_30x,
                     AMR_mean_30x, AMR_min_30x, AMR_max_30x,
                     EAS_mean_30x, EAS_min_30x, EAS_max_30x,
                     EUR_mean_30x, EUR_min_30x, EUR_max_30x,
                     SAS_mean_30x, SAS_min_30x, SAS_max_30x)

colnames(summary)
gene_symbol <- unique(summary$gene)

ccds_ids <- c()
n_ccds_ids <- c()
test_ids <- c()
n_test_ids <- c()

new_global_mean_10x <- c()
new_global_mean_20x <- c()
new_global_mean_30x <- c()
new_AFR_mean_10x <- c()
new_AFR_mean_20x <- c()
new_AFR_mean_30x <- c()
new_AMR_mean_10x <- c()
new_AMR_mean_20x <- c()
new_AMR_mean_30x <- c()
new_EAS_mean_10x <- c()
new_EAS_mean_20x <- c()
new_EAS_mean_30x <- c()
new_EUR_mean_10x <- c()
new_EUR_mean_20x <- c()
new_EUR_mean_30x <- c()
new_SAS_mean_10x <- c()
new_SAS_mean_20x <- c()
new_SAS_mean_30x <- c()

new_global_min_10x <- c()
new_global_min_20x <- c()
new_global_min_30x <- c()
new_AFR_min_10x <- c()
new_AFR_min_20x <- c()
new_AFR_min_30x <- c()
new_AMR_min_10x <- c()
new_AMR_min_20x <- c()
new_AMR_min_30x <- c()
new_EAS_min_10x <- c()
new_EAS_min_20x <- c()
new_EAS_min_30x <- c()
new_EUR_min_10x <- c()
new_EUR_min_20x <- c()
new_EUR_min_30x <- c()
new_SAS_min_10x <- c()
new_SAS_min_20x <- c()
new_SAS_min_30x <- c()

new_global_max_10x <- c()
new_global_max_20x <- c()
new_global_max_30x <- c()
new_AFR_max_10x <- c()
new_AFR_max_20x <- c()
new_AFR_max_30x <- c()
new_AMR_max_10x <- c()
new_AMR_max_20x <- c()
new_AMR_max_30x <- c()
new_EAS_max_10x <- c()
new_EAS_max_20x <- c()
new_EAS_max_30x <- c()
new_EUR_max_10x <- c()
new_EUR_max_20x <- c()
new_EUR_max_30x <- c()
new_SAS_max_10x <- c()
new_SAS_max_20x <- c()
new_SAS_max_30x <- c()

for (i in 1:length(gene_symbol)) {
  ccds_ids[i] <- paste(subset(summary, gene == gene_symbol[i])$ccds_id, collapse = "|||")
  n_ccds_ids[i] <- length(subset(summary, gene == gene_symbol[i])$ccds_id)
  test_ids[i] <- paste(subset(GTR, GeneSymbol == gene_symbol[i])$AccessionVersion, collapse = "|||")
  n_test_ids[i] <- length(subset(GTR, GeneSymbol == gene_symbol[i])$AccessionVersion)
  new_global_mean_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$global_mean_10x)
  new_global_mean_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$global_mean_20x)
  new_global_mean_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$global_mean_30x)
  new_AFR_mean_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AFR_mean_10x)
  new_AFR_mean_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AFR_mean_20x)
  new_AFR_mean_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AFR_mean_30x)
  new_AMR_mean_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AMR_mean_10x)
  new_AMR_mean_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AMR_mean_20x)
  new_AMR_mean_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AMR_mean_30x)
  new_EAS_mean_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EAS_mean_10x)
  new_EAS_mean_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EAS_mean_20x)
  new_EAS_mean_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EAS_mean_30x)
  new_EUR_mean_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EUR_mean_10x)
  new_EUR_mean_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EUR_mean_20x)
  new_EUR_mean_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EUR_mean_30x)
  new_SAS_mean_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$SAS_mean_10x)
  new_SAS_mean_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$SAS_mean_20x)
  new_SAS_mean_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$SAS_mean_30x)
  
  new_AFR_min_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AFR_min_10x)
  new_AFR_min_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AFR_min_20x)
  new_AFR_min_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AFR_min_30x)
  new_AMR_min_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AMR_min_10x)
  new_AMR_min_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AMR_min_20x)
  new_AMR_min_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AMR_min_30x)
  new_EAS_min_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EAS_min_10x)
  new_EAS_min_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EAS_min_20x)
  new_EAS_min_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EAS_min_30x)
  new_EUR_min_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EUR_min_10x)
  new_EUR_min_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EUR_min_20x)
  new_EUR_min_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EUR_min_30x)
  new_SAS_min_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$SAS_min_10x)
  new_SAS_min_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$SAS_min_20x)
  new_SAS_min_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$SAS_min_30x)
  
  new_AFR_max_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AFR_max_10x)
  new_AFR_max_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AFR_max_20x)
  new_AFR_max_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AFR_max_30x)
  new_AMR_max_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AMR_max_10x)
  new_AMR_max_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AMR_max_20x)
  new_AMR_max_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$AMR_max_30x)
  new_EAS_max_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EAS_max_10x)
  new_EAS_max_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EAS_max_20x)
  new_EAS_max_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EAS_max_30x)
  new_EUR_max_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EUR_max_10x)
  new_EUR_max_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EUR_max_20x)
  new_EUR_max_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$EUR_max_30x)
  new_SAS_max_10x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$SAS_max_10x)
  new_SAS_max_20x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$SAS_max_20x)
  new_SAS_max_30x[i] <- mean(subset(master, summary.gene == gene_symbol[i])$SAS_max_30x)
}
head(new)
new <- data.frame(gene_symbol, ccds_ids, n_ccds_ids, test_ids, n_test_ids,
                     new_AFR_mean_10x, new_AFR_min_10x, new_AFR_max_10x,
                     new_AMR_mean_10x, new_AMR_min_10x, new_AMR_max_10x,
                     new_EAS_mean_10x, new_EAS_min_10x, new_EAS_max_10x,
                     new_EUR_mean_10x, new_EUR_min_10x, new_EUR_max_10x,
                     new_SAS_mean_10x, new_SAS_min_10x, new_SAS_max_10x,
                     new_AFR_mean_20x, new_AFR_min_20x, new_AFR_max_20x,
                     new_AMR_mean_20x, new_AMR_min_20x, new_AMR_max_20x,
                     new_EAS_mean_20x, new_EAS_min_20x, new_EAS_max_20x,
                     new_EUR_mean_20x, new_EUR_min_20x, new_EUR_max_20x,
                     new_SAS_mean_20x, new_SAS_min_20x, new_SAS_max_20x,
                     new_AFR_mean_30x, new_AFR_min_30x, new_AFR_max_30x,
                     new_AMR_mean_30x, new_AMR_min_30x, new_AMR_max_30x,
                     new_EAS_mean_30x, new_EAS_min_30x, new_EAS_max_30x,
                     new_EUR_mean_30x, new_EUR_min_30x, new_EUR_max_30x,
                     new_SAS_mean_30x, new_SAS_min_30x, new_SAS_max_30x)

saveRDS(new, "../data/master_table.rds")
