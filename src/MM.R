library(edgeR)
library(cowplot)
library(ggplot2)

# Read in Counts
counts <- read.csv("../data/counts/MM_CountMatrix.csv",row.names=1)
pDat <- read.csv("../data/meta/SampleInfo.csv")
pDat <- pDat[pDat$Species=="Mouse",]
rownames(pDat) <- pDat$Sample
stopifnot(identical(rownames(pDat), colnames(counts)))


# Move stats to PD
pDat$Unmapped <- t(counts)[,1]
pDat$MultiMap <- t(counts)[,2]
pDat$NoFeat <- t(counts)[,3]
pDat$Ambiguous <- t(counts)[,4]
counts <- counts[-c(1:4),]
pDat$TotalReads <- colSums(counts)

# Add Gene Annotation
library(biomaRt)
ensembl  <- useMart("ensembl",dataset="mmusculus_gene_ensembl",host="useast.ensembl.org") # at the time of writing the ensembl.org had issues
gene.data <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), mart=ensembl)
gene.data <- gene.data[!is.na(gene.data[,1]),]
gene.data$uniq <- scater::uniquifyFeatureNames(gene.data[,1],gene.data[,2])

# this is a fuck up because they are non-unqiue fucking ENSMUSG blablablbal
gene.data <- gene.data[rownames(counts)[-((nrow(counts)-5):nrow(counts))],]
rownames(counts) <- c(gene.data$uniq,rownames(counts)[(nrow(counts)-5):nrow(counts)])

# Set up DGEList
keep <- rowSums(counts)>0
counts <- counts[keep,]
y <- DGEList(counts=counts,samples=pDat)
# Normalize
y <- calcNormFactors(y)

# Get Log-Cpm
logcpm <- cpm(y,log=TRUE)

out <- list()
for (gen in c("SV40LT","largeT","Puro","RASG12V")) {
    tmp <- pDat 
    tmp$Cond <- paste0(tmp$Transformation,tmp$Source)
    tmp$Gene <- t(logcpm)[,gen]
    out[[gen]] <- ggplot(tmp, aes(x=Cond, y=Gene)) +
		    geom_point() +
		    geom_boxplot() +
		    ggtitle(gen) +
		    ylab("LogCPM") +
		    ylim(c(-4,15)) +
		    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
plot_grid(plotlist=out)

# Compute PCA and rmv transgenes
logcpm <- logcpm[-((nrow(logcpm)-5):nrow(logcpm)),]
pcs <- prcomp(t(logcpm))
logcpm.pc <- pcs$x
vars <- apply(logcpm.pc,2,var)
rel.vars <- round(vars/sum(vars)*100,1)

pDat <- y$samples
pDat$PC1 <- logcpm.pc[,1]
pDat$PC2 <- logcpm.pc[,2]


# Plot
ggplot(pDat, aes(x=PC1, y=PC2, color=Transformation, shape=Source)) +
    geom_point(size=4) +
    theme_bw() +
    xlab(paste0("PC1(",rel.vars["PC1"],"%)")) +
    ylab(paste0("PC2(",rel.vars["PC2"],"%)")) +
    scale_color_brewer(palette="Dark2") +
    coord_fixed(1)

#Loadings
ld.pc1 <- pcs$rotation[,1]
plot(ld.pc1[order(ld.pc1,decreasing=TRUE)])

pos <- ld.pc1[ld.pc1 > 0.02]
pos <- pos[order(pos,decreasing=TRUE)]
neg <- ld.pc1[ld.pc1 < -0.02][order(ld.pc1)]
neg <- neg[order(neg)]
# tbc
