## summary table per gene
## - global mean, min, max
## - population mean, min, max
## - ANOVA F statistic
## - unadjusted P, adjusted P (FDR)
library(dplyr)
library(fst)

kgp.p3<-read.table("data/integrated_call_samples_v3.20130502.ALL.panel", header = T, sep = "\t", as.is = T, col.names = c("ID","sub_pop","pop","gender"))

collect_values<-function(dp, ref) {
  fst_file=sprintf("data/%s_%sx.fst", ref, dp)
  raw_data<-read_fst(fst_file)
  
  idset<-intersect(colnames(raw_data), kgp.p3$ID)
  kgp.map.1<-filter(kgp.p3, ID %in% idset)
  raw_data.1<-raw_data[, kgp.map.1$ID]
  
  mat<-data.frame(
    ccds_id=raw_data$ccds_id,
    gene_symbol=raw_data$gene,
    
    global_mean=apply(raw_data.1, MARGIN = 1, FUN = mean),
    global_min=apply(raw_data.1, MARGIN = 1, FUN = min),
    global_max=apply(raw_data.1, MARGIN = 1, FUN = max),
    
    AFR_mean=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "AFR" ]], MARGIN = 1, FUN = mean),
    AFR_min=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "AFR" ]], MARGIN = 1, FUN = min),
    AFR_max=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "AFR" ]], MARGIN = 1, FUN = max),
    
    AMR_mean=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "AMR" ]], MARGIN = 1, FUN = mean),
    AMR_min=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "AMR" ]], MARGIN = 1, FUN = min),
    AMR_max=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "AMR" ]], MARGIN = 1, FUN = max),
    
    EUR_mean=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "EUR" ]], MARGIN = 1, FUN = mean),
    EUR_min=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "EUR" ]], MARGIN = 1, FUN = min),
    EUR_max=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "EUR" ]], MARGIN = 1, FUN = max),
    
    EAS_mean=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "EAS" ]], MARGIN = 1, FUN = mean),
    EAS_min=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "EAS" ]], MARGIN = 1, FUN = min),
    EAS_max=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "EAS" ]], MARGIN = 1, FUN = max),
    
    SAS_mean=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "SAS" ]], MARGIN = 1, FUN = mean),
    SAS_min=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "SAS" ]], MARGIN = 1, FUN = min),
    SAS_max=apply(raw_data.1[, kgp.map.1$ID[ kgp.map.1$pop == "SAS" ]], MARGIN = 1, FUN = max),
    
    F_statistic=apply(raw_data.1, MARGIN=1, FUN=function(x){ dat=data.frame(val=x, pop=kgp.map.1$pop, stringsAsFactors = F); return(summary(aov(data=dat, formula=val ~ pop))[[1]][["F value"]][1]) }),
    p_unadj=apply(raw_data.1, MARGIN=1, FUN=function(x){ dat=data.frame(val=x, pop=kgp.map.1$pop, stringsAsFactors = F); return(summary(aov(data=dat, formula=val ~ pop))[[1]][["Pr(>F)"]][1]) }),
    
    stringsAsFactors = F
  )
  mat$p_adj=p.adjust(mat$p_unadj, method = "fdr")
  mat$depth=dp
  mat$ver=ref
  
  return(mat)
}

summary.hg38<-bind_rows(lapply(X=c(5,10,15,20,25,30,50,75,100), FUN = collect_values, ref="hg38"))
summary.b37<-bind_rows(lapply(X=c(5,10,15,20,25,30,50,75,100), FUN = collect_values, ref="b37"))

write_fst(summary.b37, path = "data/summary.b37.fst")
write_fst(summary.hg38, path = "data/summary.hg38.fst")

summary.b37.long<-melt(summary.b37, id.vars = c("ccds_id","gene_symbol","depth","ver"))
summary.hg38.long<-melt(summary.hg38, id.vars = c("ccds_id","gene_symbol","depth","ver"))

## finally merge all into long data frame
summary.long<-bind_rows(summary.b37.long, summary.hg38.long)
summary.long$depth<-paste0(summary.long$depth, "x")
write_fst(summary.long, path="data/summary.long.fst")


## prepare .fst data file for gnomAD exomes
gnomad_exome<-read.table("data/gnomad_exome.r2.1.txt", header = T, sep = "\t", as.is = T, check.names = F)
row.names(gnomad_exome)<-gnomad_exome$ccds_id
write_fst(gnomad_exome, path="data/gnomad_exome.r2.1.fst")

## transcript metadata from ensembl
library(ensembldb)
library(wiggleplotr)
library(biomaRt)
txdb = makeTxDbFromBiomart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
tx_exons = exonsBy(txdb, by="tx", use.names = T)
tx_cdss = cdsBy(txdb, by="tx", use.names = T)

ensembl_mart = useMart("ensembl", host = "grch37.ensembl.org")
ensembl_dataset = useDataset("hsapiens_gene_ensembl", mart = ensembl_mart)
tx_metadata = getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id", "external_gene_name", "strand", "gene_biotype","transcript_biotype","ccds", "hgnc_symbol"),
                    mart = ensembl_dataset)
tx_metadata = rename(tx_metadata, transcript_id = ensembl_transcript_id, gene_id = ensembl_gene_id, gene_name = external_gene_name)
save(list=c("tx_exons","tx_cdss", "tx_metadata"), file="data/tx_data.RData")

