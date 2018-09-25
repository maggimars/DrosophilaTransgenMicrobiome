library("phyloseq")
library("ggplot2")
theme_set(theme_bw())
library("dplyr")
library("tidyr")
library("RColorBrewer", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
colors<- c("grey", "#66CCFF", "#00cc66" )

P <- brewer.pal(12, "Paired")
S2 <- brewer.pal(8, "Set2")
S1 <- brewer.pal(8, "Set1")
D2 <- brewer.pal(8, "Dark2")
colors<- c(P, S2, S1,D2)
colors2 <- rep(colors, 20)

#import taxonomy (saved as tsv from taxonomy.qzv)
taxonomy<- read.delim("~/desktop/carmen/Taxonomy.tsv") 
names(taxonomy) <- c("row", "tax", "Confidence")
row.names(taxonomy) <-taxonomy[[1]]
taxonomy <- taxonomy[,(-1)]

taxonomy <-  separate(taxonomy, tax, c("D0","D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12", "D13", "D14"), sep = ";", fill = "right")
taxonomy <- taxonomy[,c(1:6)]
taxmat <- as.matrix(taxonomy)

TAX = tax_table(taxmat)

#metadata
metatable <- read.delim("~/desktop/carmen/Metadata2.txt")
#View(metatable)
row.names(metatable) <- metatable[[1]]
metatable<- metatable[,(-1)]

META <- sample_data(metatable)

#Import Phylogenetic tree (rooted, exported from qiime2)
TREE<- read_tree("~/desktop/carmen/exported-tree/tree.nwk")

#make PhyloseqObject
physeq <-  phyloseq(OTU, TAX, META, TREE)

#normalize for 1000 reads
physeqT = transform_sample_counts(physeq, function(x) 1E6 * x/sum(x))

#calculate weighted unifrac distances between samples and do PCoA
ordu = ordinate(physeqT, "PCoA", "unifrac", weighted=TRUE)

#plot PCoA
p0<-plot_ordination(physeqT, ordu, color="Treatment")+scale_color_manual(values=colors)+ggtitle("F0")

#import F1-Gut feature table 
F1table <- read.delim("~/desktop/carmen/extracted-feature-table/F1-table-bact_only/data/F1-feature-table.txt")

row.names(F1table)<-F1table[[1]]
F1table <- F1table[,(-1)]
f1mat <- as.matrix(F1table)

OTU1 = otu_table(f1mat, taxa_are_rows = TRUE)

physeq1 <-  phyloseq(OTU1, TAX, META, TREE)

physeqT1 = transform_sample_counts(physeq1, function(x) 1E6 * x/sum(x))

ordu1 = ordinate(physeqT1, "PCoA", "unifrac", weighted=TRUE)
#ordu1 = ordinate(physeqT1, "NMDS", "unifrac", weighted=FALSE)
#ordu1 = ordinate(physeqT1, "NMDS", "bray")
p1<-plot_ordination(physeqT1, ordu1, color="Treatment")+scale_color_manual(values=colors)+ ggtitle("F1")

library("DESeq2", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")

physeq1 <-  phyloseq(OTU1, TAX, META)

ddse <- phyloseq_to_deseq2(physeq1, ~Treatment)
ddse2 <- DESeq(ddse, test="Wald", fitType="parametric")
res<- results(ddse2,  cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
head(sigtab)
#NONE !

#make figure for all results

restab = cbind(as(res, "data.frame"), as(tax_table(physeq1)[rownames(res), ], "matrix"))

x = tapply(restab$log2FoldChange, restab$D4, function(x) max(x))
x = sort(x, TRUE)
restab$D4 = factor(as.character(restab$D4), levels=names(x))
x = tapply(restab$log2FoldChange, restab$D3, function(x) max(x))
x = sort(x, TRUE)
restab$D3 = factor(as.character(restab$D3), levels=names(x))
x = tapply(restab$log2FoldChange, restab$D2, function(x) max(x))
x = sort(x, TRUE)
restab$D2 = factor(as.character(restab$D2), levels=names(x))
ggplot(restab, aes(x=D4, y=log2FoldChange)) + geom_point(size=1) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +ggtitle("F1")

#theme(legend.position="none")



#import F2/F3 (T)-Gut feature table 
FTtable <- read.delim("~/desktop/carmen/extracted-feature-table/T-table-bact_only/data/T-feature-table.txt")

row.names(FTtable)<-FTtable[[1]]
FTtable <- FTtable[,(-1)]
fTmat <- as.matrix(FTtable)

OTUT = otu_table(fTmat, taxa_are_rows = TRUE)

physeqTT <-  phyloseq(OTUT, TAX, META, TREE)

physeqTTT = transform_sample_counts(physeqTT, function(x) 1E6 * x/sum(x))

orduT = ordinate(physeqTTT, "PCoA", "unifrac", weighted=TRUE)
pT<-plot_ordination(physeqTTT, orduT, color="Treatment")+scale_color_manual(values=colors)+ggtitle("F2 & F3")

#F0batch1
F01table <- read.delim("~/desktop/carmen/extracted-feature-table/1F0-table-bact_only/data/1F0-feature-table.txt")

row.names(F01table)<-F01table[[1]]
F01table <- F01table[,(-1)]
f01mat <- as.matrix(F01table)

OTU01 = otu_table(f01mat, taxa_are_rows = TRUE)

physeq01 <-  phyloseq(OTU01, TAX, META, TREE)

physeqT01 = transform_sample_counts(physeq01, function(x) 1E6 * x/sum(x))

ordu01 = ordinate(physeqT01, "PCoA", "unifrac", weighted=TRUE)

p01<-plot_ordination(physeqT01, ordu01, color="Treatment")+scale_color_manual(values=colors)+ ggtitle("F0")

#Deseq2 
physeq1 <-  phyloseq(OTU01, TAX, META)

ddse <- phyloseq_to_deseq2(physeq1, ~Treatment)
ddse2 <- DESeq(ddse, test="Wald", fitType="parametric")
res<- results(ddse2,  cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq1)[rownames(sigtab), ], "matrix"))
restab = cbind(as(res, "data.frame"), as(tax_table(physeq1)[rownames(res), ], "matrix"))
head(sigtab)

dim(sigtab)

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$D4, function(x) max(x))
x = sort(x, TRUE)
sigtab$D4 = factor(as.character(sigtab$D4), levels=names(x))
ggplot(sigtab, aes(x=D4, y=log2FoldChange, color=D4)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#Figure with all 
x = tapply(restab$log2FoldChange, restab$D4, function(x) max(x))
x = sort(x, TRUE)
restab$D4 = factor(as.character(restab$D4), levels=names(x))
restab$D3 = factor(as.character(restab$D3), levels=names(x))
ggplot(restab, aes(x=D4, y=log2FoldChange)) + geom_point(size=1) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ theme(legend.position="none") +ggtitle("F0")

#convert to relative abundance F0

sumcols <- colSums(F01table)
for(i in 1:ncol(F01table)) {
  F01table[[i]]<- (F01table[[i]]/sumcols[[i]])*100
}

f01mat <- as.matrix(F01table)
OTU01 = otu_table(f01mat, taxa_are_rows = TRUE)
physeq01 <-  phyloseq(OTU01, TAX, META, TREE)

#taxa bar plots
colors2<- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00" ,"#FFFF33" ,"#A65628" ,"#F781BF", "#999999", "#8DD3C7" ,"#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3" ,"#FDB462" ,"#B3DE69" ,"#FCCDE5", "#D9D9D9", "#BC80BD",
            "#CCEBC5" ,"#FFED6F", "#000000")

D2F0taxa<-plot_bar(physeq01, fill = "D2") + scale_y_continuous(expand = c(0, 0)) + ggtitle("F0 - Level 2") + scale_fill_manual(values=colors2) + theme(legend.title=element_blank()) 
D4F0taxa<-plot_bar(physeq01, fill = "D4") + scale_y_continuous(expand = c(0, 0)) + ggtitle("F0 - Level 4") +theme(legend.position="none") 
  
#convert to relative abundance F1

sumcols <- colSums(F1table)
for(i in 1:ncol(F1table)) {
  F1table[[i]]<- (F1table[[i]]/sumcols[[i]])*100
}

f1mat <- as.matrix(F1table)
OTU1 = otu_table(f1mat, taxa_are_rows = TRUE)
physeq1 <-  phyloseq(OTU1, TAX, META, TREE)

#taxa bar plots
D2F1taxa<-plot_bar(physeq1, fill = "D2") + scale_y_continuous(expand = c(0, 0)) + ggtitle("F1 - Level 2") + scale_fill_manual(values=colors2) + theme(legend.title=element_blank()) + geom_bar(aes(color=D2, fill=D2), stat="identity", position="stack")
D4F1taxa<-plot_bar(physeq1, fill = "D4") + scale_y_continuous(expand = c(0, 0)) + ggtitle("F1 - Level 4") +theme(legend.position="none") + geom_bar(width = 0.75, aes(color=D4, fill=D4), stat="identity", position="stack") + scale_fill_manual(values=colors2) + scale_color_manual(values=colors2)

D4F1taxa<-plot_bar(physeq1, fill = "D4") + scale_y_continuous(expand = c(0, 0)) + ggtitle("F1 - Level 4") + geom_bar(width = 0.75, aes(color=D4, fill=D4), stat="identity", position="stack") + scale_fill_manual(values=colors2) + scale_color_manual(values=colors2)


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
legend <- g_legend(D4F1taxa)
Prettylegend<- grid.arrange(legend)


#convert to relative abundance F2/3

sumcols <- colSums(FTtable)
for(i in 1:ncol(FTtable)) {
  FTtable[[i]]<- (FTtable[[i]]/sumcols[[i]])*100
}

Tmat <- as.matrix(FTtable)
OTUT = otu_table(Tmat, taxa_are_rows = TRUE)
physeqT <-  phyloseq(OTUT, TAX, META, TREE)

#taxa bar plots
D2FTtaxa<-plot_bar(physeqT, fill = "D2") + scale_y_continuous(expand = c(0, 0)) + ggtitle("F2&3 - Level 2") + theme(legend.title=element_blank()) +theme(legend.position="none") 
D4FTtaxa<-plot_bar(physeqT, fill = "D4") + scale_y_continuous(expand = c(0, 0)) + ggtitle("F2&3 - Level 4") +theme(legend.position="none") 


