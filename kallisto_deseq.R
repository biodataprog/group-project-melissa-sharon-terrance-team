#!/usr/bin/env Rscript
library(DESeq2)
library(tximport)
library(dplyr)
library(ggplot2)
library(magrittr)
library(Biobase)
library(pheatmap)
library(RColorBrewer)

samples <- read.table("samplerun.csv",header=TRUE,sep=",")
samples$Name = sprintf("%s.%s.%s",samples$Genotype,samples$Treatment,samples$Replicate)
samples$Name
files <- file.path("output",samples$Name,"abundance.tsv")
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)
colnames(txi.kallisto$counts) <- samples$Name
colnames(txi.kallisto$abundance) <- samples$Name
write.csv(txi.kallisto$abundance,"reports/kallisto.TPM.csv")
write.csv(txi.kallisto$counts,"reports/kallisto.counts.csv")

# DEseq2 analyses
geno = samples$Genotype
treatment = samples$Treatment

sampleTable <- data.frame(condition=treatment,
                         genotype = geno)

#sampleTable <- data.frame(condition=treatment)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$genotype <- factor(sampleTable$genotype)
rownames(sampleTable) = samples$Name


dds <- DESeqDataSetFromTximport(txi.kallisto,sampleTable,~ condition + genotype)
#dds <- DESeqDataSetFromTximport(txi.kallisto,sampleTable, ~ condition)

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

dds <- estimateSizeFactors(dds)
rld <- rlog(dds, blind = TRUE)
#vsd <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "parametric")
#vsd <- vst(dds)

head(assay(rld), 3)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
  #as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "varianceStabilizingTransformation"))

colnames(df)[1:2] <- c("x", "y")

#pdf("plots/RNASeq_kallisto.pdf")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:30]
df <- as.data.frame(colData(dds)[,c("condition","genotype")])
#df <- as.data.frame(colData(dds)[,c("condition")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE,annotation_col=sampleTable)

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$genotype, sep="-")
#rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

pcaData <- plotPCA(rld, intgroup=c("condition", "genotype"), returnData=TRUE)
#pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype, shape=condition)) +
#ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
