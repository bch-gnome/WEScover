library(dplyr)
library(tidyr)
library(fst)

## prepare GTR table (gpt_tP_tG) from latest file
download.file("https://ftp.ncbi.nlm.nih.gov/pub/GTR/data/test_version.gz", destfile = "data/test_version.gz", method = "wget")

tv <- read.delim("data/test_version.gz", stringsAsFactors = F, header = T, sep = "\t")
tv <- subset(tv, now_current == 1 & test_currStat == 1 & test_pubStat ==3
             & genes != "" & lab_test_name != genes,
             select = c(test_accession_ver,
                        genes,
                        lab_test_name,
                        condition_identifiers,
                        CLIA_number))
tv <- tv[!grepl("|", tv$condition_identifiers, fixed = T),]


split_genes<-function(row) {
  ngenes <- unlist(strsplit(row[2], "|", fixed = T))
  
  res<-data.frame(
    GTR_accession = rep(row[1], times = length(ngenes)),
    test_name = rep(row[3], times = length(ngenes)),
    phenotype_name = rep(row[4], times = length(ngenes)),
    gene_symbol = ngenes,
    stringsAsFactors = F
  )
  
  return(res)
}

new_tv<- bind_rows(apply(X=tv, MARGIN=1, FUN=split_genes))
new_tv<- new_tv[!duplicated(new_tv), ]

write_fst(new_tv, "data/gpt_tP_tG.fst")
