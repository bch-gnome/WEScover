#genes_by_ccds_id <- read_fst('data/genes_by_ccds_id.fst')

#summary <- read_fst("data/summary.fst")

getCCDS <- function(geneN, sm) {
  #message("getCCDS - ", paste(geneN, collapse = ", "))
  unlist(lapply(geneN, function(x) {
    sm$ccds_id[sm$gene == x]
  }))
}

load10x <- function(summary) {
  if(!"AFR_10x" %in% ls(envir = .GlobalEnv)) {
    AFR_10x <- read_fst("data/AFR_10x.fst")
    rownames(AFR_10x) <- summary$ccds_id
    AFR_10x$CCDS <- rownames(AFR_10x)
    AFR_10x$Population <- "AFR"
    AFR_10x$GeneSymbol <- summary[AFR_10x$CCDS, "gene"]
    AMR_10x <- read_fst("data/AMR_10x.fst")
    rownames(AMR_10x) <- summary$ccds_id
    AMR_10x$CCDS <- rownames(AMR_10x)
    AMR_10x$Population <- "AMR"
    AMR_10x$GeneSymbol <- summary[AMR_10x$CCDS, "gene"]
    EAS_10x <- read_fst("data/EAS_10x.fst")
    rownames(EAS_10x) <- summary$ccds_id
    EAS_10x$CCDS <- rownames(EAS_10x)
    EAS_10x$Population <- "EAS"
    EAS_10x$GeneSymbol <- summary[EAS_10x$CCDS, "gene"]
    EUR_10x <- read_fst("data/EUR_10x.fst")
    rownames(EUR_10x) <- summary$ccds_id
    EUR_10x$CCDS <- rownames(EUR_10x)
    EUR_10x$Population <- "EUR"
    EUR_10x$GeneSymbol <- summary[EUR_10x$CCDS, "gene"]
    SAS_10x <- read_fst("data/SAS_10x.fst")
    rownames(SAS_10x) <- summary$ccds_id
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
    AFR_20x <- read_fst("data/AFR_20x.fst")
    rownames(AFR_20x) <- summary$ccds_id
    AFR_20x$CCDS <- rownames(AFR_20x)
    AFR_20x$Population <- "AFR"
    AFR_20x$GeneSymbol <- summary[AFR_20x$CCDS, "gene"]
    AMR_20x <- read_fst("data/AMR_20x.fst")
    rownames(AMR_20x) <- summary$ccds_id
    AMR_20x$CCDS <- rownames(AMR_20x)
    AMR_20x$Population <- "AMR"
    AMR_20x$GeneSymbol <- summary[AMR_20x$CCDS, "gene"]
    EAS_20x <- read_fst("data/EAS_20x.fst")
    rownames(EAS_20x) <- summary$ccds_id
    EAS_20x$CCDS <- rownames(EAS_20x)
    EAS_20x$Population <- "EAS"
    EAS_20x$GeneSymbol <- summary[EAS_20x$CCDS, "gene"]
    EUR_20x <- read_fst("data/EUR_20x.fst")
    rownames(EUR_20x) <- summary$ccds_id
    EUR_20x$CCDS <- rownames(EUR_20x)
    EUR_20x$Population <- "EUR"
    EUR_20x$GeneSymbol <- summary[EUR_20x$CCDS, "gene"]
    SAS_20x <- read_fst("data/SAS_20x.fst")
    rownames(SAS_20x) <- summary$ccds_id
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
    AFR_30x <- read_fst("data/AFR_30x.fst")
    rownames(AFR_30x) <- summary$ccds_id
    AFR_30x$CCDS <- rownames(AFR_30x)
    AFR_30x$Population <- "AFR"
    AFR_30x$GeneSymbol <- summary[AFR_30x$CCDS, "gene"]
    AMR_30x <- read_fst("data/AMR_30x.fst")
    rownames(AMR_30x) <- summary$ccds_id
    AMR_30x$CCDS <- rownames(AMR_30x)
    AMR_30x$Population <- "AMR"
    AMR_30x$GeneSymbol <- summary[AMR_30x$CCDS, "gene"]
    EAS_30x <- read_fst("data/EAS_30x.fst")
    rownames(EAS_30x) <- summary$ccds_id
    EAS_30x$CCDS <- rownames(EAS_30x)
    EAS_30x$Population <- "EAS"
    EAS_30x$GeneSymbol <- summary[EAS_30x$CCDS, "gene"]
    EUR_30x <- read_fst("data/EUR_30x.fst")
    rownames(EUR_30x) <- summary$ccds_id
    EUR_30x$CCDS <- rownames(EUR_30x)
    EUR_30x$Population <- "EUR"
    EUR_30x$GeneSymbol <- summary[EUR_30x$CCDS, "gene"]
    SAS_30x <- read_fst("data/SAS_30x.fst")
    rownames(SAS_30x) <- summary$ccds_id
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

createGPT <- function(row, main_table, gpt_tP_tG) {
# createGPT <- function(row, main_table, gtrM) {
  if (as.character(main_table[row, 1]) %in% as.character(gpt_tP_tG$gene_symbol)) {
    tbl <- data.frame(
      A = as.character(main_table[row, 1]),
      B = gpt_tP_tG[gpt_tP_tG$gene_symbol == as.character(main_table[row, 1]), "GTR_accession"],
      C = gpt_tP_tG[gpt_tP_tG$gene_symbol == as.character(main_table[row, 1]), "test_name"],
      # B = unique(gtrM[gtrM$GeneSymbol == as.character(main_table[row, 1]), "AccessionVersion"]),
      # C = unique(gtrM[gtrM$GeneSymbol == as.character(main_table[row, 1]), "ObjectName"]),
      stringsAsFactors = FALSE
    )
    createLink <- function(val) {
      sprintf('<a href="https://www.ncbi.nlm.nih.gov/gtr/tests/%s" target="_blank" class="btn btn-primary">%s</a>',val, val)
    }
    tbl$B <- sapply(tbl$B, createLink)
    
    colnames(tbl) <- c("Gene Symbol", "GTR Test ID", "Test Name")
    rownames(tbl) <- seq(nrow(tbl))
    tbl <- tbl[!duplicated(tbl), ]
    tbl
  } else {
    tbl <- data.frame()
  }
}

# createMainTable2 <- function(geneS, depth, gpt_tP_tG) {
# createMainTable2 <- function(geneS, depth, summary, gtrM, gtrS) {
createMainTable2 <- function(geneS, dp, summary, assembly) {
  library(dplyr)
  
  if (length(geneS) > 0) {
    tmp <- filter(summary, gene_symbol %in% geneS, depth == dp, ver==assembly)
    
    xx <- reshape2::dcast(tmp, ccds_id + gene_symbol ~ variable, value.var = 'value')
    xx <- xx[,c("gene_symbol","ccds_id","global_mean","global_min","global_max","AFR_mean","AMR_mean","EAS_mean","EUR_mean","SAS_mean","F_statistic","p_unadj","p_adj")]
    xx <- xx[order(xx[, 3]), ]
    
    colnames(xx) <- c("Gene Symbol", "CCDS ID", "Global coverage (mean, %)", "Global coverage (min, %)", "Global coverage (max, %)", "AFR (mean, %)", 
                      "AMR (%)", "EAS (%)", "EUR (%)", "SAS (%)", "ANOVA F-statistic",
                      "Raw P-Value","Adj. P-Value")
    rownames(xx) <- seq(nrow(xx))
    #colnames(genes_by_ccds_id)
    for (ii in 3:10) {
      xx[,ii] <- round(xx[,ii] * 100.0, digits=2)
    }
    xx[,11] <- round(xx[,11], digits=2)
    xx[,12] <- format.pval(xx[,12], digits = 3, eps = 0.001)
    xx[,13] <- format.pval(xx[,13], digits = 3, eps = 0.001)
    xx
  } else {
    data.frame()
  }
  
}


globalMean2 <- function(idx, ccds, summary) {
  if(idx==3) {
    load10x(summary)
    list(
      "AFR" = genes_by_ccds_id[ccds,"AFR_mean_10x"],
      "AFR.MI" = genes_by_ccds_id[ccds,"AFR_min_10x"],
      "AFR.MA" = genes_by_ccds_id[ccds,"AFR_max_10x"],
      "AMR" = genes_by_ccds_id[ccds,"AMR_mean_10x"],
      "AMR.MI" = genes_by_ccds_id[ccds,"AMR_min_10x"],
      "AMR.MA" = genes_by_ccds_id[ccds,"AMR_max_10x"],
      "EAS" = genes_by_ccds_id[ccds,"EAS_mean_10x"],
      "EAS.MI" = genes_by_ccds_id[ccds,"EAS_min_10x"],
      "EAS.MA" = genes_by_ccds_id[ccds,"EAS_max_10x"],
      "EUR" = genes_by_ccds_id[ccds,"EUR_mean_10x"],
      "EUR.MI" = genes_by_ccds_id[ccds,"EUR_min_10x"],
      "EUR.MA" = genes_by_ccds_id[ccds,"EUR_max_10x"],
      "SAS" = genes_by_ccds_id[ccds,"SAS_mean_10x"],
      "SAS.MI" = genes_by_ccds_id[ccds,"SAS_min_10x"],
      "SAS.MA" = genes_by_ccds_id[ccds,"SAS_max_10x"],
      "GLOBAL"= genes_by_ccds_id[ccds,"global_mean_10x"],
      "GLOBAL.MI"= genes_by_ccds_id[ccds,"global_min_10x"],
      "GLOBAL.MA"= genes_by_ccds_id[ccds,"global_max_10x"],
      "F.10X" = genes_by_ccds_id[ccds,"F_statistic_10x"],
      "p.10x" = genes_by_ccds_id[ccds, "p_value_10x"]
    )
  } else if(idx == 4) {
    load20x(summary)
    list(
      "AFR" = genes_by_ccds_id[ccds,"AFR_mean_20x"],
      "AFR.MI" = genes_by_ccds_id[ccds,"AFR_min_20x"],
      "AFR.MA" = genes_by_ccds_id[ccds,"AFR_max_20x"],
      "AMR" = genes_by_ccds_id[ccds,"AMR_mean_20x"],
      "AMR.MI" = genes_by_ccds_id[ccds,"AMR_min_20x"],
      "AMR.MA" = genes_by_ccds_id[ccds,"AMR_max_20x"],
      "EAS" = genes_by_ccds_id[ccds,"EAS_mean_20x"],
      "EAS.MI" = genes_by_ccds_id[ccds,"EAS_min_20x"],
      "EAS.MA" = genes_by_ccds_id[ccds,"EAS_max_20x"],
      "EUR" = genes_by_ccds_id[ccds,"EUR_mean_20x"],
      "EUR.MI" = genes_by_ccds_id[ccds,"EUR_min_20x"],
      "EUR.MA" = genes_by_ccds_id[ccds,"EUR_max_20x"],
      "SAS" = genes_by_ccds_id[ccds,"SAS_mean_20x"],
      "SAS.MI" = genes_by_ccds_id[ccds,"SAS_min_20x"],
      "SAS.MA" = genes_by_ccds_id[ccds,"SAS_max_20x"],
      "GLOBAL"= genes_by_ccds_id[ccds,"global_mean_20x"],
      "GLOBAL.MI"= genes_by_ccds_id[ccds,"global_min_20x"],
      "GLOBAL.MA"= genes_by_ccds_id[ccds,"global_max_20x"],
      "F.20X" = genes_by_ccds_id[ccds,"F_statistic_20x"],
      "p.20x" = genes_by_ccds_id[ccds, "p_value_20x"]
    )
  } else if(idx == 5) {
    load30x(summary)
    list(
      "AFR" = genes_by_ccds_id[ccds,"AFR_mean_30x"],
      "AFR.MI" = genes_by_ccds_id[ccds,"AFR_min_30x"],
      "AFR.MA" = genes_by_ccds_id[ccds,"AFR_max_30x"],
      "AMR" = genes_by_ccds_id[ccds,"AMR_mean_30x"],
      "AMR.MI" = genes_by_ccds_id[ccds,"AMR_min_30x"],
      "AMR.MA" = genes_by_ccds_id[ccds,"AMR_max_30x"],
      "EAS" = genes_by_ccds_id[ccds,"EAS_mean_30x"],
      "EAS.MI" = genes_by_ccds_id[ccds,"EAS_min_30x"],
      "EAS.MA" = genes_by_ccds_id[ccds,"EAS_max_30x"],
      "EUR" = genes_by_ccds_id[ccds,"EUR_mean_30x"],
      "EUR.MI" = genes_by_ccds_id[ccds,"EUR_min_30x"],
      "EUR.MA" = genes_by_ccds_id[ccds,"EUR_max_30x"],
      "SAS" = genes_by_ccds_id[ccds,"SAS_mean_30x"],
      "SAS.MI" = genes_by_ccds_id[ccds,"SAS_min_30x"],
      "SAS.MA" = genes_by_ccds_id[ccds,"SAS_max_30x"],
      "GLOBAL"= genes_by_ccds_id[ccds,"global_mean_30x"],
      "GLOBAL.MI"= genes_by_ccds_id[ccds,"global_min_30x"],
      "GLOBAL.MA"= genes_by_ccds_id[ccds,"global_max_30x"],
      "F.30X" = genes_by_ccds_id[ccds,"F_statistic_30x"],
      "p.30x" = genes_by_ccds_id[ccds, "p_value_30x"]
    )
  }
}
