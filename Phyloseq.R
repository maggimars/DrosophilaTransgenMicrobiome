library("phyloseq")
library("ggplot2")
theme_set(theme_bw())
library("dplyr")
library("tidyr")
library("RColorBrewer")
library(gridExtra)
library(vegan)

P <- brewer.pal(12, "Paired")
S2 <- brewer.pal(8, "Set2")
S1 <- brewer.pal(8, "Set1")
D2 <- brewer.pal(8, "Dark2")
colors<- c(P, S2, S1,D2)
colors2<- c(P, S2, S1,D2)

#import taxonomy (saved as tsv from taxonomy.qzv)
taxonomy<- read.delim("Taxonomy.tsv")
names(taxonomy) <- c("row", "tax", "Confidence")
row.names(taxonomy) <-taxonomy[[1]]
taxonomy <- taxonomy[,(-1)]

taxonomy <-  separate(taxonomy, tax, c("D0","D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12", "D13", "D14"), sep = ";", fill = "right")
taxonomy <- taxonomy[,c(1:6)]
taxmat <- as.matrix(taxonomy)

TAX = tax_table(taxmat)

#import metadata
metatable <- read.delim("Metadata2.txt")
#View(metatable)
row.names(metatable) <- metatable[[1]]
metatable<- metatable[,(-1)]

META <- sample_data(metatable)

#Import Phylogenetic tree (rooted, exported from qiime2)
TREE<- read_tree("tree.nwk")

#Import Feature Table
Svtab <- read.delim("feature-table.txt")
row.names(Svtab)<-Svtab[[1]] #make OTU ID the row names
Svtab<-Svtab[,-(1)] # remove OTU ID column
fmat <- as.matrix(Svtab)
OTU = otu_table(fmat, taxa_are_rows = TRUE)

#create phyloseq object
ps<- phyloseq(OTU, TAX, META, TREE)

#subset to only include GUT samples (exclude media samples)

ps<- subset_samples(ps, SampleType=="Gut")

rank_names(ps)
table(tax_table(ps)[, "D0"], exclude = NULL)
#remove sequences that have D0 unassigned
# Filter Physeq object - Bacteria ONLY
ps <- subset_taxa(ps, D0=="D_0__Bacteria")
table(tax_table(ps)[, "D2"], exclude = NULL)
#remove phyla that cannot be part of microbiome
filterPhyla = c("D_2__Chloroplast" , " D_2__Cyanobacteria")
ps = subset_taxa(ps, !D2 %in% filterPhyla)

prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

plyr::ddply(prevdf, "D2", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(ps),color=D2)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~D2) + theme(legend.position="none")

prevalenceThreshold = 0.1 * nsamples(ps) # remove OTS found in less than 10% of the samples
prevalenceThreshold

keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)

table(tax_table(ps2)[, "D2"], exclude = NULL)

prevdf2 = apply(X = otu_table(ps2),
               MARGIN = ifelse(taxa_are_rows(ps2), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf2 = data.frame(Prevalence = prevdf2,
                    TotalAbundance = taxa_sums(ps2),
                    tax_table(ps2))

plyr::ddply(prevdf2, "D1", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

ggplot(prevdf2, aes(TotalAbundance, Prevalence / nsamples(ps),color=D2)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~D1) + theme(legend.position="none")

#transform to relative abundance
ps2ra = transform_sample_counts(ps2, function(x){x / sum(x)})

#subset by generation
ps2raF0 <- subset_samples(ps2ra, Generation=="F0")
ps2raF1 <- subset_samples(ps2ra, Generation=="F1")

ordu1F0 = ordinate(ps2raF0, "PCoA", "unifrac", weighted=TRUE)
ordu1F1 = ordinate(ps2raF1, "PCoA", "unifrac", weighted=TRUE)

p1F0<-plot_ordination(ps2raF0, ordu1F0, color="Treatment")+scale_color_manual(values=colors)+ ggtitle("F0 - Weighted Unifrac")
p1F1<-plot_ordination(ps2raF1, ordu1F1, color="Treatment")+scale_color_manual(values=colors)+ ggtitle("F1 - Weighted Unifrac")

ordu2F0 = ordinate(ps2raF0, "PCoA", "bray")
ordu2F1 = ordinate(ps2raF1, "PCoA", "bray")

p2F0<- plot_ordination(ps2raF0, ordu2F0, color="Treatment")+scale_color_manual(values=colors)+ ggtitle("F0 - Bray-Curtis")
p2F1<- plot_ordination(ps2raF1, ordu2F1, color="Treatment")+scale_color_manual(values=colors)+ ggtitle("F1 - Bray-Curtis")

#Relabundance Plots
taxabarplot<-plot_bar(ps2ra, fill = "D2") + scale_y_continuous(expand = c(0, 0)) + ggtitle("Level 2") + scale_fill_manual(values=colors2) + theme(legend.title=element_blank()) + geom_bar(aes( fill=D2), stat="identity", position="stack")

taxabarplot + theme(legend.position="none")

#legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
legend <- g_legend(taxabarplot)
SVlegend <- grid.arrange(legend)

##################################
#Remove all taxa that are D1 == NA

filterPhyla = c("D_1__Acidobacteria", "D_1__Actinobacteria", "D_1__Bacteroidetes", "D_1__Firmicutes", " D_1__Parcubacteria", "D_1__Proteobacteria", "D_1__Verrucomicrobia")
ps3 = subset_taxa(ps2, D1 %in% filterPhyla)
table(tax_table(ps3)[, "D1"], exclude = NULL)

#Final Prevalence Plot
prevdf3 = apply(X = otu_table(ps3),
                MARGIN = ifelse(taxa_are_rows(ps3), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})

prevdf3 = data.frame(Prevalence = prevdf3,
                     TotalAbundance = taxa_sums(ps3),
                     tax_table(ps3))

plyr::ddply(prevdf3, "D1", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

ggplot(prevdf3, aes(TotalAbundance, Prevalence / nsamples(ps),color=D2)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~D2) + theme(legend.position="none")

#transform to relative abundance
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})

#subset by generation
ps3raF0 <- subset_samples(ps3ra, Generation=="F0")
ps3raF1 <- subset_samples(ps3ra, Generation=="F1")

ordu1F0 = ordinate(ps3raF0, "PCoA", "unifrac", weighted=TRUE)
ordu1F1 = ordinate(ps3raF1, "PCoA", "unifrac", weighted=TRUE)

p1F0<-plot_ordination(ps3raF0, ordu1F0, color="Treatment")+scale_color_manual(values=colors)+ ggtitle("F0 - Weighted Unifrac")
p1F1<-plot_ordination(ps3raF1, ordu1F1, color="Treatment")+scale_color_manual(values=colors)+ ggtitle("F1 - Weighted Unifrac")

ordu2F0 = ordinate(ps3raF0, "PCoA", "bray")
ordu2F1 = ordinate(ps3raF1, "PCoA", "bray")

p2F0<- plot_ordination(ps3raF0, ordu2F0, color="Treatment")+scale_color_manual(values=colors)+ ggtitle("F0 - Bray-Curtis")
p2F1<- plot_ordination(ps3raF1, ordu2F1, color="Treatment")+scale_color_manual(values=colors)+ ggtitle("F1 - Bray-Curtis")

#Relabundance Plots
filterGens = c("F0", "F1")
ps4 <- subset_samples(ps3ra, Generation %in% filterGens)
taxabarplot<-plot_bar(ps4, fill = "D2") + scale_y_continuous(expand = c(0, 0)) + ggtitle("F0 & F1 - Level 2 Taxonomy") + scale_fill_manual(values=colors2) + theme(legend.title=element_blank()) + geom_bar(aes( fill=D2), stat="identity", position="stack")

taxabarplot + theme(legend.position="none")

#legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
legend <- g_legend(taxabarplot)
SVlegend <- grid.arrange(legend)


#PERMANOVA with Vegan

#F0
OTUs <- t(data.frame(otu_table(ps3raF0))) #get data frame of OTUs from phyloseq object
meta <- metatable[metatable$Generation=="F0",]
meta <- meta[meta$SampleType=="Gut",]

adonis(vegdist(OTUs, method = "bray") ~ Treatment, data = meta)

#RESULTS = SIGNIFICANT
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Treatment  2   0.57473 0.28737  6.9801 0.39932  0.001 ***
#Residuals 21   0.86456 0.04117         0.60068
#Total     23   1.43929                 1.00000

#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#pairwise CD ~ NS F0
CDNS0 <- subset_samples(ps3raF0, Treatment %in% c("CD" , "NSD"))
OTUs <- t(data.frame(otu_table(CDNS0))) #get data frame of OTUs from phyloseq object
metaPW <- meta[meta$Treatment %in% c("NSD", "CD"),]
adonis(vegdist(OTUs, method = "bray") ~ Treatment, data = metaPW)

#Results = Significant
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#Treatment  1   0.36891 0.36891  9.0982 0.39389  0.001 ***
#Residuals 14   0.56767 0.04055         0.60611
#Total     15   0.93658                 1.00000

#F1
OTUs <- t(data.frame(otu_table(ps3raF1))) #get data frame of OTUs from phyloseq object
meta <- metatable[metatable$Generation=="F1",]
meta <- meta[meta$SampleType=="Gut",]
adonis(vegdist(OTUs, method = "bray") ~ Treatment, data = meta)

#RESULTS = Not Significant
#          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
#Treatment  2    0.1999 0.099948  1.0886 0.04928  0.336
#Residuals 42    3.8561 0.091812         0.95072
#Total     44    4.0560                  1.00000

#pairwise CD ~ NS F1
CDNS1 <- subset_samples(ps3raF1, Treatment %in% c("CD" , "NSD"))
OTUs <- t(data.frame(otu_table(CDNS1))) #get data frame of OTUs from phyloseq object
metaPW <- meta[meta$Treatment %in% c("NSD", "CD"),]
adonis(vegdist(OTUs, method = "bray") ~ Treatment, data = metaPW)

#RESULTS = NOT Significant
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
#Treatment  1    0.1456 0.145605  1.6786 0.05656   0.17
##Total     29    2.5744                  1.00000

## DESeq2 ##
### Differential abundance of features ####
library("DESeq2", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
#fully filtered (prevalence, taxonomy) physeq object  = ps3

#subset by Diet
testdiets <- c("CD", "NSD")
ps3cn <- subset_samples(ps3, Treatment %in% testdiets)

#subset by generation
ps3F0 <- subset_samples(ps3cn, Generation=="F0")
ps3F1 <- subset_samples(ps3cn, Generation=="F1")

#F0
ddse <- phyloseq_to_deseq2(ps3F0, ~Treatment)
ddse2 <- DESeq(ddse, test="Wald", fitType="parametric")
res<- results(ddse2,  cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
head(sigtab)
#15 Svs are Significantly differentially abundant between NSD and CD in the F0
# log2 fold change (MLE): Treatment NSD vs CD
#NSD was compared to CD
# + meand more in NSD and - means less in NSD

#add taxonomy to significant features
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps3F0)[rownames(sigtab), ], "matrix"))

#figure for F0 results
#D2 level
x = tapply(sigtab$log2FoldChange, sigtab$D2, function(x) max(x))
x = sort(x, TRUE)
sigtab$D2 = factor(as.character(sigtab$D2), levels=names(x))
# D5 level
x = tapply(sigtab$log2FoldChange, sigtab$D5, function(x) max(x))
x = sort(x, TRUE)
sigtab$D5 = factor(as.character(sigtab$D5), levels=names(x))
#plot F0
F0sigplot<- ggplot(sigtab, aes(x=D2, y=log2FoldChange, color = D5)) + geom_point(size=3, alpha = 0.7) + theme(legend.title=element_blank()) +
   ggtitle("F0 - log2 fold change (MLE): Treatment NSD vs CD ") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))# +scale_color_manual(values=P)

F0sigplot

#make fold change table for all results (not significant)
restab = cbind(as(res, "data.frame"), as(tax_table(ps3F0)[rownames(res), ], "matrix"))
restabCC <- restab[complete.cases(restab),]
#figure for all F0 results
F0resplot<- ggplot(restabCC, aes(x=D2, y=log2FoldChange)) + geom_point(size=3, alpha = 0.7) + theme(legend.title=element_blank()) +
  ggtitle("F1 - log2 fold change (MLE): Treatment NSD vs CD ") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))# +scale_color_manual(values=P)

F0resplot

#F1
ddse <- phyloseq_to_deseq2(ps3F1, ~Treatment)
ddse2 <- DESeq(ddse, test="Wald", fitType="parametric")
res<- results(ddse2,  cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
head(sigtab)
# No significantly different Svs
#log2 fold change (MLE): Treatment NSD vs CD

#make fold change table for all results (not significant)
restab = cbind(as(res, "data.frame"), as(tax_table(ps3F1)[rownames(res), ], "matrix"))

restabCC <- restab[complete.cases(restab),]
#figure for F1 results

F1resplot<- ggplot(restabCC, aes(x=D2, y=log2FoldChange)) + geom_point(size=3, alpha = 0.7) + theme(legend.title=element_blank()) +
  ggtitle("F1 - log2 fold change (MLE): Treatment NSD vs CD ") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))# +scale_color_manual(values=P)

F1resplot



