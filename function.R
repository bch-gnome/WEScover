genes_by_ccds_id <- readRDS('../data/genes_by_ccds_id.rds')

getCCDS <- function(geneN, sm) {
  #message("getCCDS - ", paste(geneN, collapse = ", "))
  unlist(lapply(geneN, function(x) {
    sm$ccds_id[sm$gene == x]
  }))
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


globalMean <- function(idx, ccds, summary) {
  if(idx==3) {
    load10x(summary)
    list(
      "AFR" = apply(AFR_10x[ccds, seq(ncol(AFR_10x) - 3)], MARGIN = 1, mean, na.rm = T),
      "AFR.MI" = apply(AFR_10x[ccds, seq(ncol(AFR_10x) - 3)], MARGIN = 1, min, na.rm = T),
      "AFR.MA" = apply(AFR_10x[ccds, seq(ncol(AFR_10x) - 3)], MARGIN = 1, max, na.rm = T),
      "AMR" = apply(AMR_10x[ccds, seq(ncol(AMR_10x) - 3)], MARGIN = 1, mean, na.rm = T),
      "AMR.MI" = apply(AMR_10x[ccds, seq(ncol(AMR_10x) - 3)], MARGIN = 1, min, na.rm = T),
      "AMR.MA" = apply(AMR_10x[ccds, seq(ncol(AMR_10x) - 3)], MARGIN = 1, max, na.rm = T),
      "EAS" = apply(EAS_10x[ccds, seq(ncol(EAS_10x) - 3)], MARGIN = 1, mean, na.rm = T),
      "EAS.MI" = apply(EAS_10x[ccds, seq(ncol(EAS_10x) - 3)], MARGIN = 1, min, na.rm = T),
      "EAS.MA" = apply(EAS_10x[ccds, seq(ncol(EAS_10x) - 3)], MARGIN = 1, max, na.rm = T),
      "EUR" = apply(EUR_10x[ccds, seq(ncol(EUR_10x) - 3)], MARGIN = 1, mean, na.rm = T),
      "EUR.MI" = apply(EUR_10x[ccds, seq(ncol(EUR_10x) - 3)], MARGIN = 1, min, na.rm = T),
      "EUR.MA" = apply(EUR_10x[ccds, seq(ncol(EUR_10x) - 3)], MARGIN = 1, max, na.rm = T),
      "SAS" = apply(SAS_10x[ccds, seq(ncol(SAS_10x) - 3)], MARGIN = 1, mean, na.rm = T),
      "SAS.MI" = apply(SAS_10x[ccds, seq(ncol(SAS_10x) - 3)], MARGIN = 1, min, na.rm = T),
      "SAS.MA" = apply(SAS_10x[ccds, seq(ncol(SAS_10x) - 3)], MARGIN = 1, max, na.rm = T),
      "GLOBAL"= apply(
        do.call(cbind, list(AFR_10x[ccds, seq(ncol(AFR_10x) - 3)], AMR_10x[ccds, seq(ncol(AMR_10x) - 3)],
                            EAS_10x[ccds, seq(ncol(EAS_10x) - 3)], EUR_10x[ccds, seq(ncol(EUR_10x) - 3)],
                            SAS_10x[ccds, seq(ncol(SAS_10x) - 3)])), 
        MARGIN = 1, mean, na.rm = T),
      "GLOBAL.MI"= apply(
        do.call(cbind, list(AFR_10x[ccds, seq(ncol(AFR_10x) - 3)], AMR_10x[ccds, seq(ncol(AMR_10x) - 3)],
                            EAS_10x[ccds, seq(ncol(EAS_10x) - 3)], EUR_10x[ccds, seq(ncol(EUR_10x) - 3)],
                            SAS_10x[ccds, seq(ncol(SAS_10x) - 3)])), 
        MARGIN = 1, min, na.rm = T),
      "GLOBAL.MA"= apply(
        do.call(cbind, list(AFR_10x[ccds, seq(ncol(AFR_10x) - 3)], AMR_10x[ccds, seq(ncol(AMR_10x) - 3)],
                            EAS_10x[ccds, seq(ncol(EAS_10x) - 3)], EUR_10x[ccds, seq(ncol(EUR_10x) - 3)],
                            SAS_10x[ccds, seq(ncol(SAS_10x) - 3)])), 
        MARGIN = 1, max, na.rm = T)
    )
  } else if(idx == 4) {
    load20x(summary)
    list(
      "AFR" = apply(AFR_20x[ccds, seq(ncol(AFR_20x) - 3)], MARGIN = 1, mean, na.rm = T),
      "AFR.MI" = apply(AFR_20x[ccds, seq(ncol(AFR_20x) - 3)], MARGIN = 1, min, na.rm = T),
      "AFR.MA" = apply(AFR_20x[ccds, seq(ncol(AFR_20x) - 3)], MARGIN = 1, max, na.rm = T),
      "AMR" = apply(AMR_20x[ccds, seq(ncol(AMR_20x) - 3)], MARGIN = 1, mean, na.rm = T),
      "AMR.MI" = apply(AMR_20x[ccds, seq(ncol(AMR_20x) - 3)], MARGIN = 1, min, na.rm = T),
      "AMR.MA" = apply(AMR_20x[ccds, seq(ncol(AMR_20x) - 3)], MARGIN = 1, max, na.rm = T),
      "EAS" = apply(EAS_20x[ccds, seq(ncol(EAS_20x) - 3)], MARGIN = 1, mean, na.rm = T),
      "EAS.MI" = apply(EAS_20x[ccds, seq(ncol(EAS_20x) - 3)], MARGIN = 1, min, na.rm = T),
      "EAS.MA" = apply(EAS_20x[ccds, seq(ncol(EAS_20x) - 3)], MARGIN = 1, max, na.rm = T),
      "EUR" = apply(EUR_20x[ccds, seq(ncol(EUR_20x) - 3)], MARGIN = 1, mean, na.rm = T),
      "EUR.MI" = apply(EUR_20x[ccds, seq(ncol(EUR_20x) - 3)], MARGIN = 1, min, na.rm = T),
      "EUR.MA" = apply(EUR_20x[ccds, seq(ncol(EUR_20x) - 3)], MARGIN = 1, max, na.rm = T),
      "SAS" = apply(SAS_20x[ccds, seq(ncol(SAS_20x) - 3)], MARGIN = 1, mean, na.rm = T),
      "SAS.MI" = apply(SAS_20x[ccds, seq(ncol(SAS_20x) - 3)], MARGIN = 1, min, na.rm = T),
      "SAS.MA" = apply(SAS_20x[ccds, seq(ncol(SAS_20x) - 3)], MARGIN = 1, max, na.rm = T),
      "GLOBAL"= apply(
        do.call(cbind, list(AFR_20x[ccds, seq(ncol(AFR_20x) - 3)], AMR_20x[ccds, seq(ncol(AMR_20x) - 3)],
                            EAS_20x[ccds, seq(ncol(EAS_20x) - 3)], EUR_20x[ccds, seq(ncol(EUR_20x) - 3)],
                            SAS_20x[ccds, seq(ncol(SAS_20x) - 3)])), 
        MARGIN = 1, mean, na.rm = T),
      "GLOBAL.MI"= apply(
        do.call(cbind, list(AFR_20x[ccds, seq(ncol(AFR_20x) - 3)], AMR_20x[ccds, seq(ncol(AMR_20x) - 3)],
                            EAS_20x[ccds, seq(ncol(EAS_20x) - 3)], EUR_20x[ccds, seq(ncol(EUR_20x) - 3)],
                            SAS_20x[ccds, seq(ncol(SAS_20x) - 3)])), 
        MARGIN = 1, min, na.rm = T),
      "GLOBAL.MA"= apply(
        do.call(cbind, list(AFR_20x[ccds, seq(ncol(AFR_20x) - 3)], AMR_20x[ccds, seq(ncol(AMR_20x) - 3)],
                            EAS_20x[ccds, seq(ncol(EAS_20x) - 3)], EUR_20x[ccds, seq(ncol(EUR_20x) - 3)],
                            SAS_20x[ccds, seq(ncol(SAS_20x) - 3)])), 
        MARGIN = 1, max, na.rm = T)
    )
  } else if(idx == 5) {
    load30x(summary)
    list(
      "AFR" = apply(AFR_30x[ccds, seq(ncol(AFR_30x) - 3)], MARGIN = 1, mean, na.rm = T),
      "AFR.MI" = apply(AFR_30x[ccds, seq(ncol(AFR_30x) - 3)], MARGIN = 1, min, na.rm = T),
      "AFR.MA" = apply(AFR_30x[ccds, seq(ncol(AFR_30x) - 3)], MARGIN = 1, max, na.rm = T),
      "AMR" = apply(AMR_30x[ccds, seq(ncol(AMR_30x) - 3)], MARGIN = 1, mean, na.rm = T),
      "AMR.MI" = apply(AMR_30x[ccds, seq(ncol(AMR_30x) - 3)], MARGIN = 1, min, na.rm = T),
      "AMR.MA" = apply(AMR_30x[ccds, seq(ncol(AMR_30x) - 3)], MARGIN = 1, max, na.rm = T),
      "EAS" = apply(EAS_30x[ccds, seq(ncol(EAS_30x) - 3)], MARGIN = 1, mean, na.rm = T),
      "EAS.MI" = apply(EAS_30x[ccds, seq(ncol(EAS_30x) - 3)], MARGIN = 1, min, na.rm = T),
      "EAS.MA" = apply(EAS_30x[ccds, seq(ncol(EAS_30x) - 3)], MARGIN = 1, max, na.rm = T),
      "EUR" = apply(EUR_30x[ccds, seq(ncol(EUR_30x) - 3)], MARGIN = 1, mean, na.rm = T),
      "EUR.MI" = apply(EUR_30x[ccds, seq(ncol(EUR_30x) - 3)], MARGIN = 1, min, na.rm = T),
      "EUR.MA" = apply(EUR_30x[ccds, seq(ncol(EUR_30x) - 3)], MARGIN = 1, max, na.rm = T),
      "SAS" = apply(SAS_30x[ccds, seq(ncol(SAS_30x) - 3)], MARGIN = 1, mean, na.rm = T),
      "SAS.MI" = apply(SAS_30x[ccds, seq(ncol(SAS_30x) - 3)], MARGIN = 1, min, na.rm = T),
      "SAS.MA" = apply(SAS_30x[ccds, seq(ncol(SAS_30x) - 3)], MARGIN = 1, max, na.rm = T),
      "GLOBAL"= apply(
        do.call(cbind, list(AFR_30x[ccds, seq(ncol(AFR_30x) - 3)], AMR_30x[ccds, seq(ncol(AMR_30x) - 3)],
                            EAS_30x[ccds, seq(ncol(EAS_30x) - 3)], EUR_30x[ccds, seq(ncol(EUR_30x) - 3)],
                            SAS_30x[ccds, seq(ncol(SAS_30x) - 3)])), 
        MARGIN = 1, mean, na.rm = T),
      "GLOBAL.MI"= apply(
        do.call(cbind, list(AFR_30x[ccds, seq(ncol(AFR_30x) - 3)], AMR_30x[ccds, seq(ncol(AMR_30x) - 3)],
                            EAS_30x[ccds, seq(ncol(EAS_30x) - 3)], EUR_30x[ccds, seq(ncol(EUR_30x) - 3)],
                            SAS_30x[ccds, seq(ncol(SAS_30x) - 3)])), 
        MARGIN = 1, min, na.rm = T),
      "GLOBAL.MA"= apply(
        do.call(cbind, list(AFR_30x[ccds, seq(ncol(AFR_30x) - 3)], AMR_30x[ccds, seq(ncol(AMR_30x) - 3)],
                            EAS_30x[ccds, seq(ncol(EAS_30x) - 3)], EUR_30x[ccds, seq(ncol(EUR_30x) - 3)],
                            SAS_30x[ccds, seq(ncol(SAS_30x) - 3)])), 
        MARGIN = 1, max, na.rm = T)
    )
  }
}


createGPT <- function(row, main_table, gtrM) {
  tbl <- data.frame(
    A = as.character(main_table[row, 1]),
    B = unique(gtrM[gtrM$GeneSymbol == as.character(main_table[row, 1]), "AccessionVersion"]),
    C = unique(gtrM[gtrM$GeneSymbol == as.character(main_table[row, 1]), "ObjectName"]),
    stringsAsFactors = FALSE
  )
  createLink <- function(val) {
    sprintf('<a href="https://www.ncbi.nlm.nih.gov/gtr/tests/%s" target="_blank" class="btn btn-primary">%s</a>',val, val)
  }
  tbl$B <- sapply(tbl$B, createLink)
  
  colnames(tbl) <- c("Gene Symbol", "Version", "Gene Panel Testing")
  rownames(tbl) <- seq(nrow(tbl))
  tbl
}

createMainTable <- function(geneS, depth, summary, gtrM, gtrS) {
  selCCDS <- getCCDS(geneS, summary)

  message("[TABLE] Number of CCDS: ", length(selCCDS), "; Number of genes:", length(geneS) )
  idx <- ifelse(depth == "10x", 3, ifelse(depth == "20x", 4, 5))
  xx <- do.call(rbind, lapply(geneS, function(gene) {
    message(gene)
    selCCDS <- getCCDS(gene, summary)
    if(length(selCCDS) > 0) {
      GS <- gtrS[selCCDS, "gene"]
      GPT <- length(unique(gtrM[gtrM$GeneSymbol == gene, "AccessionVersion"]))
      MM <- globalMean(idx, selCCDS, summary)
      data.frame(
        A = GS,
        B = selCCDS,
        C = GPT,
        #D = paste(round(MM$GLOBAL * 100, 0), " (", round(MM$GLOBAL.MI * 100, 0), "-", round(MM$GLOBAL.MA * 100, 0), ")", sep=""),
        D = round(MM$GLOBAL * 100, 0),
        E = paste(round(MM$AFR * 100, 0), " (", round(MM$AFR.MI * 100, 0), "-", round(MM$AFR.MA * 100, 0), ")", sep=""),
        F = paste(round(MM$AMR * 100, 0), " (", round(MM$AMR.MI * 100, 0), "-", round(MM$AMR.MA * 100, 0), ")", sep=""),
        G = paste(round(MM$EAS * 100, 0), " (", round(MM$EAS.MI * 100, 0), "-", round(MM$EAS.MA * 100, 0), ")", sep=""),
        H = paste(round(MM$EUR * 100, 0), " (", round(MM$EUR.MI * 100, 0), "-", round(MM$EUR.MA * 100, 0), ")", sep=""),
        I = paste(round(MM$SAS * 100, 0), " (", round(MM$SAS.MI * 100, 0), "-", round(MM$SAS.MA * 100, 0), ")", sep=""),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        A = gene, B = NA, C = NA, D = NA,
        E = NA, F = NA, G = NA, H = NA,
        I = NA, stringsAsFactors = FALSE
      )
    }
  }))
  colnames(xx) <- c("Gene Symbol", "CCDS", "Gene Panel Testing", "GM", "AFR (%)", 
                    "AMR (%)", "EAS (%)", "EUR (%)", "SAS (%)")
  rownames(xx) <- seq(nrow(xx))
  xx

}

# Jeff's versions

createMainTable2 <- function(geneS, depth, summary, gtrM, gtrS) {
  # selCCDS <- getCCDS(geneS, summary)
  # 
  # message("[TABLE] Number of CCDS: ", length(selCCDS), "; Number of genes:", length(geneS) )
  # idx <- ifelse(depth == "10x", 3, ifelse(depth == "20x", 4, 5))
  # xx <- do.call(rbind, lapply(geneS, function(gene) {
  #   message(gene)
  #   selCCDS <- getCCDS(gene, summary)
  #   if(length(selCCDS) > 0) {
  #     GS <- gtrS[selCCDS, "gene"]
  #     GPT <- length(unique(gtrM[gtrM$GeneSymbol == gene, "AccessionVersion"]))
  #     MM <- globalMean2(idx, selCCDS, summary)
  #     data.frame(
  #       A = GS,
  #       B = selCCDS,
  #       C = GPT,
  #       #D = paste(round(MM$GLOBAL * 100, 0), " (", round(MM$GLOBAL.MI * 100, 0), "-", round(MM$GLOBAL.MA * 100, 0), ")", sep=""),
  #       D = MM$GLOBAL,
  #       E = paste(MM$AFR, " (", MM$AFR.MI, "-", MM$AFR.MA, ")", sep=""),
  #       F = paste(MM$AMR, " (", MM$AMR.MI, "-", MM$AMR.MA, ")", sep=""),
  #       G = paste(MM$EAS, " (", MM$EAS.MI, "-", MM$EAS.MA, ")", sep=""),
  #       H = paste(MM$EUR, " (", MM$EUR.MI, "-", MM$EUR.MA, ")", sep=""),
  #       I = paste(MM$SAS, " (", MM$SAS.MI, "-", MM$SAS.MA, ")", sep=""),
  #       stringsAsFactors = FALSE
  #     )
  #   } else {
  #     data.frame(
  #       A = gene, B = NA, C = NA, D = NA,
  #       E = NA, F = NA, G = NA, H = NA,
  #       I = NA, stringsAsFactors = FALSE
  #     )
  #   }
  # }))
  if (length(geneS) > 0) {
    if(depth == "10x") {
      colS <- c("gene_symbol", "ccds_id", "global_mean_10x", "global_min_10x", "global_max_10x", "AFR_mean_10x", "AMR_mean_10x", "EAS_mean_10x",
                "EUR_mean_10x", "SAS_mean_10x")
    } else if(depth == "20x") {
      colS <- c("gene_symbol", "ccds_id", "global_mean_20x", "global_min_20x", "global_max_20x","AFR_mean_20x", "AMR_mean_20x", "EAS_mean_20x",
              "EUR_mean_20x", "SAS_mean_20x")
    } else {
      colS <- c("gene_symbol", "ccds_id", "global_mean_30x","global_min_30x", "global_max_30x", "AFR_mean_30x", "AMR_mean_30x", "EAS_mean_30x",
                "EUR_mean_30x", "SAS_mean_30x")
    }
    xx <- genes_by_ccds_id[genes_by_ccds_id$gene_symbol %in% geneS, colS]
    xx <- xx[order(xx[, 3]), ]
    colnames(xx) <- c("Gene Symbol", "CCDS", "Global mean", "Global min", "Global max", "AFR (%)", 
                      "AMR (%)", "EAS (%)", "EUR (%)", "SAS (%)")
    rownames(xx) <- seq(nrow(xx))
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
      "GLOBAL.MA"= genes_by_ccds_id[ccds,"global_max_10x"]
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
      "GLOBAL.MA"= genes_by_ccds_id[ccds,"global_max_20x"]
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
      "GLOBAL.MA"= genes_by_ccds_id[ccds,"global_max_30x"]
    )
  }
}
