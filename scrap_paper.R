summary <- readRDS("../data/summary.rds")
identical(summary$ccds_id, rownames(global_10x))
AFR_10x <- readRDS("../data/AFR_10x.rds")
AFR_20x <- readRDS("../data/AFR_20x.rds")
AFR_30x <- readRDS("../data/AFR_30x.rds")
AMR_10x <- readRDS("../data/AMR_10x.rds")
AMR_20x <- readRDS("../data/AMR_20x.rds")
AMR_30x <- readRDS("../data/AMR_30x.rds")
EAS_10x <- readRDS("../data/EAS_10x.rds")
EAS_20x <- readRDS("../data/EAS_20x.rds")
EAS_30x <- readRDS("../data/EAS_30x.rds")
EUR_10x <- readRDS("../data/EUR_10x.rds")
EUR_20x <- readRDS("../data/EUR_20x.rds")
EUR_30x <- readRDS("../data/EUR_30x.rds")
SAS_10x <- readRDS("../data/SAS_10x.rds")
SAS_20x <- readRDS("../data/SAS_20x.rds")
SAS_30x <- readRDS("../data/SAS_30x.rds")

global_10x <- read.delim("../../coverage_tables/GTR_10x.txt", header = T, stringsAsFactors = F)
global_20x <- read.delim("../../coverage_tables/GTR_20x.txt", header = T, stringsAsFactors = F)
global_30x <- read.delim("../../coverage_tables/GTR_30x.txt", header = T, stringsAsFactors = F)

pops <- list(global_10x, global_20x, global_30x,
             AFR_10x, AFR_20x, AFR_30x,
             AMR_10x, AMR_20x, AMR_30x,
             EAS_10x, EAS_20x, EAS_30x,
             EUR_10x, EUR_20x, EUR_30x,
             SAS_10x, SAS_20x, SAS_30x)

global_mean_10x <- c()
global_mean_20x <- c()
global_mean_30x <- c()
AFR_mean_10x <- c()
AFR_mean_20x <- c()
AFR_mean_30x <- c()
AMR_mean_10x <- c()
AMR_mean_20x <- c()
AMR_mean_30x <- c()
EAS_mean_10x <- c()
EAS_mean_20x <- c()
EAS_mean_30x <- c()
EUR_mean_10x <- c()
EUR_mean_20x <- c()
EUR_mean_30x <- c()
SAS_mean_10x <- c()
SAS_mean_20x <- c()
SAS_mean_30x <- c()

global_min_10x <- c()
global_min_20x <- c()
global_min_30x <- c()
AFR_min_10x <- c()
AFR_min_20x <- c()
AFR_min_30x <- c()
AMR_min_10x <- c()
AMR_min_20x <- c()
AMR_min_30x <- c()
EAS_min_10x <- c()
EAS_min_20x <- c()
EAS_min_30x <- c()
EUR_min_10x <- c()
EUR_min_20x <- c()
EUR_min_30x <- c()
SAS_min_10x <- c()
SAS_min_20x <- c()
SAS_min_30x <- c()

global_max_10x <- c()
global_max_20x <- c()
global_max_30x <- c()
AFR_max_10x <- c()
AFR_max_20x <- c()
AFR_max_30x <- c()
AMR_max_10x <- c()
AMR_max_20x <- c()
AMR_max_30x <- c()
EAS_max_10x <- c()
EAS_max_20x <- c()
EAS_max_30x <- c()
EUR_max_10x <- c()
EUR_max_20x <- c()
EUR_max_30x <- c()
SAS_max_10x <- c()
SAS_max_20x <- c()
SAS_max_30x <- c()

pops <- list(global_10x, global_20x, global_30x,
             AFR_10x, AFR_20x, AFR_30x,
             AMR_10x, AMR_20x, AMR_30x,
             EAS_10x, EAS_20x, EAS_30x,
             EUR_10x, EUR_20x, EUR_30x,
             SAS_10x, SAS_20x, SAS_30x)

means <- list(global_mean_10x, global_mean_20x, global_mean_30x,
              AFR_mean_10x, AFR_mean_20x, AFR_mean_30x,
              AMR_mean_10x, AMR_mean_20x, AMR_mean_30x,
              EAS_mean_10x, EAS_mean_20x, EAS_mean_30x,
              EUR_mean_10x, EUR_mean_20x, EUR_mean_30x,
              SAS_mean_10x, SAS_mean_20x, SAS_mean_30x)

mins <- list(global_min_10x, global_min_20x, global_min_30x,
             AFR_min_10x, AFR_min_20x, AFR_min_30x,
             AMR_min_10x, AMR_min_20x, AMR_min_30x,
             EAS_min_10x, EAS_min_20x, EAS_min_30x,
             EUR_min_10x, EUR_min_20x, EUR_min_30x,
             SAS_min_10x, SAS_min_20x, SAS_min_30x)

maxs <- list(global_max_10x, global_max_20x, global_max_30x,
             AFR_max_10x, AFR_max_20x, AFR_max_30x,
             AMR_max_10x, AMR_max_20x, AMR_max_30x,
             EAS_max_10x, EAS_max_20x, EAS_max_30x,
             EUR_max_10x, EUR_max_20x, EUR_max_30x,
             SAS_max_10x, SAS_max_20x, SAS_max_30x)

for (i in 1:length(pops)) {
  means[[i]] <- round(rowMeans(pops[[i]]*100),0)
  mins[[i]] <- round(apply(pops[[i]]*100, 1, FUN=min),0)
  maxs[[i]] <- round(apply(pops[[i]]*100, 1, FUN=max),0)
}

genes_by_ccds_id <- data.frame(summary$ccds_id, summary$gene, means, mins, maxs)
colnames(genes_by_ccds_id) <- c("ccds_id", "gene_symbol",
                                "global_mean_10x", "global_mean_20x", "global_mean_30x",
                                "AFR_mean_10x", "AFR_mean_20x", "AFR_mean_30x",
                                "AMR_mean_10x", "AMR_mean_20x", "AMR_mean_30x",
                                "EAS_mean_10x", "EAS_mean_20x", "EAS_mean_30x",
                                "EUR_mean_10x", "EUR_mean_20x", "EUR_mean_30x",
                                "SAS_mean_10x", "SAS_mean_20x", "SAS_mean_30x",
                                "global_min_10x", "global_min_20x", "global_min_30x",
                                "AFR_min_10x", "AFR_min_20x", "AFR_min_30x",
                                "AMR_min_10x", "AMR_min_20x", "AMR_min_30x",
                                "EAS_min_10x", "EAS_min_20x", "EAS_min_30x",
                                "EUR_min_10x", "EUR_min_20x", "EUR_min_30x",
                                "SAS_min_10x", "SAS_min_20x", "SAS_min_30x",
                                "global_max_10x", "global_max_20x", "global_max_30x",
                                "AFR_max_10x", "AFR_max_20x", "AFR_max_30x",
                                "AMR_max_10x", "AMR_max_20x", "AMR_max_30x",
                                "EAS_max_10x", "EAS_max_20x", "EAS_max_30x",
                                "EUR_max_10x", "EUR_max_20x", "EUR_max_30x",
                                "SAS_max_10x", "SAS_max_20x", "SAS_max_30x")

# genes_by_ccds_id <- apply(genes_by_ccds_id[,-1:-2], 1, function(x) round(x*100, 0))
#saveRDS(genes_by_ccds_id, file = "genes_by_ccds_id.rds")
#write.csv(genes_by_ccds_id, file = "../../genes_by_ccds_id.csv", row.names = F)

###### genes by gene symbol ######

gene_symbol <- unique(summary$gene)

