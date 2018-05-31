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
  # if(length(selCCDS) <= 9) {
  #   data.frame()
  # } else {
    message("[TABLE] Number of CCDS: ", length(selCCDS), "; Number of genes:", length(geneS) )
    idx <- ifelse(depth == "10x", 3, ifelse(depth == "20x", 4, 5))
    xx <- do.call(rbind, lapply(geneS, function(geneS) {
      selCCDS <- getCCDS(geneS, summary)
      GS <- gtrS[selCCDS, "gene"]
      GPT <- length(unique(gtrM[gtrM$GeneSymbol == geneS, "AccessionVersion"]))
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
    }))
    colnames(xx) <- c("Gene Symbol", "CCDS", "Gene Panel Testing", "GM", "AFR (min-max)", 
                      "AMR (min-max)", "EAS (min-max)", "EUR (min-max)", "SAS (min-max)")
    rownames(xx) <- seq(nrow(xx))
    xx
  # }
}
