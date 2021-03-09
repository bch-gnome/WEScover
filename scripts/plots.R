## 1. compare gene-level summary from 1000 genomes vs. gnomad exomes (v2.1)

per_gene_summary<-read.fst("data/summary.long.fst")
gnomad_exome <- read.fst("data/gnomad_exome.r2.1.fst")

filter(per_gene_summary, ver == "b37", variable == "global_mean") -> aa
colnames(aa)[6]<-"kgp"

reshape2::melt(gnomad_exome, id.vars = "ccds_id") -> bb
colnames(bb)[3]<-"gnomad"
bb$variable<-as.character(bb$variable)

xx<-inner_join(x=aa[,c("ccds_id","gene_symbol","depth","kgp")], y=bb, by=c("ccds_id"="ccds_id", "depth" = "variable"))
xx$depth<-factor(xx$depth, levels = c("5x","10x","15x","20x","25x","30x","50x","100x"), labels = c("5x","10x","15x","20x","25x","30x","50x","100x"))


library(gridExtra)
library(patchwork)
png(file="data/kgp_vs_gnomad.png", width = 6.5, height = 9, units = "in", res = 300)
grid.arrange(
  ggplot(filter(xx, depth %in% c("5x","10x"))) + geom_point(aes(x=kgp, y=gnomad), size=1) + geom_abline(slope = 1, intercept = 0, col="red") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + facet_grid(.~depth) + xlab('') + ylab(''),
  ggplot(filter(xx, depth %in% c("15x","20x"))) + geom_point(aes(x=kgp, y=gnomad), size=1) + geom_abline(slope = 1, intercept = 0, col="red") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + facet_grid(.~depth) + xlab('') + ylab(''),
  ggplot(filter(xx, depth %in% c("25x","30x"))) + geom_point(aes(x=kgp, y=gnomad), size=1) + geom_abline(slope = 1, intercept = 0, col="red") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + facet_grid(.~depth) + xlab('') + ylab(''),
  ggplot(filter(xx, depth %in% c("50x","100x"))) + geom_point(aes(x=kgp, y=gnomad), size=1) + geom_abline(slope = 1, intercept = 0, col="red") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + facet_grid(.~depth) + xlab('') + ylab(''),
  nrow=4,
  left = "Mean fraction of gnomAD exomes over X",
  bottom = "Mean breadth of coverage across 1000 Genomes exomes"
)
dev.off()

cor.test(xx$kgp[ xx$depth=="5x" ], xx$gnomad[ xx$depth == "5x"], method = "spearman") # 0.79
cor.test(xx$kgp[ xx$depth=="10x" ], xx$gnomad[ xx$depth == "10x"], method = "spearman") # 0.78
cor.test(xx$kgp[ xx$depth=="15x" ], xx$gnomad[ xx$depth == "15x"], method = "spearman") # 0.76
cor.test(xx$kgp[ xx$depth=="20x" ], xx$gnomad[ xx$depth == "20x"], method = "spearman") # 0.72
cor.test(xx$kgp[ xx$depth=="25x" ], xx$gnomad[ xx$depth == "25x"], method = "spearman") # 0.68
cor.test(xx$kgp[ xx$depth=="30x" ], xx$gnomad[ xx$depth == "30x"], method = "spearman") # 0.64
cor.test(xx$kgp[ xx$depth=="50x" ], xx$gnomad[ xx$depth == "50x"], method = "spearman") # 0.52
cor.test(xx$kgp[ xx$depth=="100x" ], xx$gnomad[ xx$depth == "100x"], method = "spearman") # 0.36

cor.test(xx$kgp[ xx$depth=="5x" ], xx$gnomad[ xx$depth == "5x"], method = "pearson") # 0.86
cor.test(xx$kgp[ xx$depth=="10x" ], xx$gnomad[ xx$depth == "10x"], method = "pearson") # 0.87
cor.test(xx$kgp[ xx$depth=="15x" ], xx$gnomad[ xx$depth == "15x"], method = "pearson") # 0.87
cor.test(xx$kgp[ xx$depth=="20x" ], xx$gnomad[ xx$depth == "20x"], method = "pearson") # 0.86
cor.test(xx$kgp[ xx$depth=="25x" ], xx$gnomad[ xx$depth == "25x"], method = "pearson") # 0.83
cor.test(xx$kgp[ xx$depth=="30x" ], xx$gnomad[ xx$depth == "30x"], method = "pearson") # 0.79
cor.test(xx$kgp[ xx$depth=="50x" ], xx$gnomad[ xx$depth == "50x"], method = "pearson") # 0.63
cor.test(xx$kgp[ xx$depth=="100x" ], xx$gnomad[ xx$depth == "100x"], method = "pearson") # 0.36


## 2. Distribution of global mean breadth of coverage 

tmp<-filter(per_gene_summary, variable %in% c("global_mean"), depth != "75x")
tmp$ver<-factor(tmp$ver, levels = c("b37","hg38"), labels = c("GRCh37 (b37/hg19)", "GRCh38 (hg38)"))
tmp$depth<-factor(tmp$depth, levels = c("5x","10x","15x","20x","25x","30x","50x", "100x"), labels = c("5x","10x","15x","20x","25x","30x","50x", "100x"))
pdf(file="data/hg38_vs_b37.pdf", width=10, height=10)
ggplot(tmp, aes(x=depth, y=value)) + geom_boxplot(aes(fill=ver), notch = T) + xlab('Threshold for target read depth') + ylab('Global mean for breadth of coverage') + labs(fill="") + theme_bw()
dev.off()

