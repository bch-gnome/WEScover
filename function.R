getCCDS <- function(geneN, db) {
  message("getCCDS - ", paste(geneN, collapse = ", "))
  db$ccds_id[db$gene %in% geneN]
}

load10x <- function(summary) {
  if(!"AFR_10x" %in% ls(envir = .GlobalEnv)) {
    AFR_10x <- readRDS("../data/AFR_10x.rds")
    AFR_10x$CCDS <- rownames(AFR_10x)
    AFR_10x$Population <- "AFR"
    AFR_10x$GeneSymbol <- summary[AFR_10x$CCDS, "gene"]
    AMR_10x <- readRDS("../data/AMR_10x.rds")
    AMR_10x$CCDS <- rownames(AMR_10x)
    AMR_10x$Population <- "AMR"
    AMR_10x$GeneSymbol <- summary[AMR_10x$CCDS, "gene"]
    EAS_10x <- readRDS("../data/EAS_10x.rds")
    EAS_10x$CCDS <- rownames(EAS_10x)
    EAS_10x$Population <- "EAS"
    EAS_10x$GeneSymbol <- summary[EAS_10x$CCDS, "gene"]
    EUR_10x <- readRDS("../data/EUR_10x.rds")
    EUR_10x$CCDS <- rownames(EUR_10x)
    EUR_10x$Population <- "EUR"
    EUR_10x$GeneSymbol <- summary[EUR_10x$CCDS, "gene"]
    SAS_10x <- readRDS("../data/SAS_10x.rds")
    SAS_10x$CCDS <- rownames(SAS_10x)
    SAS_10x$Population <- "SAS"
    SAS_10x$GeneSymbol <- summary[SAS_10x$CCDS, "gene"]
    
    assign("AFR_10x", AFR_10x, envir = .GlobalEnv)
    assign("AMR_10x", AMR_10x, envir = .GlobalEnv)
    assign("EAS_10x", EAS_10x, envir = .GlobalEnv)
    assign("EUR_10x", EUR_10x, envir = .GlobalEnv)
    assign("SAS_10x", SAS_10x, envir = .GlobalEnv)
  }
}

load20x <- function(summary) {
  if(!"AFR_20x" %in% ls(envir = .GlobalEnv)) {
    AFR_20x <- readRDS("../data/AFR_20x.rds")
    AFR_20x$CCDS <- rownames(AFR_20x)
    AFR_20x$Population <- "AFR"
    AFR_20x$GeneSymbol <- summary[AFR_20x$CCDS, "gene"]
    AMR_20x <- readRDS("../data/AMR_20x.rds")
    AMR_20x$CCDS <- rownames(AMR_20x)
    AMR_20x$Population <- "AMR"
    AMR_20x$GeneSymbol <- summary[AMR_20x$CCDS, "gene"]
    EAS_20x <- readRDS("../data/EAS_20x.rds")
    EAS_20x$CCDS <- rownames(EAS_20x)
    EAS_20x$Population <- "EAS"
    EAS_20x$GeneSymbol <- summary[EAS_20x$CCDS, "gene"]
    EUR_20x <- readRDS("../data/EUR_20x.rds")
    EUR_20x$CCDS <- rownames(EUR_20x)
    EUR_20x$Population <- "EUR"
    EUR_20x$GeneSymbol <- summary[EUR_20x$CCDS, "gene"]
    SAS_20x <- readRDS("../data/SAS_20x.rds")
    SAS_20x$CCDS <- rownames(SAS_20x)
    SAS_20x$Population <- "SAS"
    SAS_20x$GeneSymbol <- summary[SAS_20x$CCDS, "gene"]
    
    assign("AFR_20x", AFR_20x, envir = .GlobalEnv)
    assign("AMR_20x", AMR_20x, envir = .GlobalEnv)
    assign("EAS_20x", EAS_20x, envir = .GlobalEnv)
    assign("EUR_20x", EUR_20x, envir = .GlobalEnv)
    assign("SAS_20x", SAS_20x, envir = .GlobalEnv)
  }
}

load30x <- function(summary) {
  if(!"AFR_30x" %in% ls(envir = .GlobalEnv)) {
    AFR_30x <- readRDS("../data/AFR_30x.rds")
    AFR_30x$CCDS <- rownames(AFR_30x)
    AFR_30x$Population <- "AFR"
    AFR_30x$GeneSymbol <- summary[AFR_30x$CCDS, "gene"]
    AMR_30x <- readRDS("../data/AMR_30x.rds")
    AMR_30x$CCDS <- rownames(AMR_30x)
    AMR_30x$Population <- "AMR"
    AMR_30x$GeneSymbol <- summary[AMR_30x$CCDS, "gene"]
    EAS_30x <- readRDS("../data/EAS_30x.rds")
    EAS_30x$CCDS <- rownames(EAS_30x)
    EAS_30x$Population <- "EAS"
    EAS_30x$GeneSymbol <- summary[EAS_30x$CCDS, "gene"]
    EUR_30x <- readRDS("../data/EUR_30x.rds")
    EUR_30x$CCDS <- rownames(EUR_30x)
    EUR_30x$Population <- "EUR"
    EUR_30x$GeneSymbol <- summary[EUR_30x$CCDS, "gene"]
    SAS_30x <- readRDS("../data/SAS_30x.rds")
    SAS_30x$CCDS <- rownames(SAS_30x)
    SAS_30x$Population <- "SAS"
    SAS_30x$GeneSymbol <- summary[SAS_30x$CCDS, "gene"]
    
    assign("AFR_30x", AFR_30x, envir = .GlobalEnv)
    assign("AMR_30x", AMR_30x, envir = .GlobalEnv)
    assign("EAS_30x", EAS_30x, envir = .GlobalEnv)
    assign("EUR_30x", EUR_30x, envir = .GlobalEnv)
    assign("SAS_30x", SAS_30x, envir = .GlobalEnv)
  }
}

