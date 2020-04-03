# Script to compute differential expression of different NMR cell lines
library(edgeR)
library(cowplot)
library(ggplot2)
library(viridis)
library(biomaRt)
library(ggrepel)
library(dplyr)
theme_set(theme_cowplot())

# Read in Counts
counts <- read.csv("../data/counts/NMR_Counts_Clean.csv",row.names=1)
pDat <- read.csv("../data/meta/NMR_ColData.csv")
orthologs <- read.csv("../data/meta/media-28.csv")

# Add Gene Annotation from the orthologs
orthologs <- read.csv("../data/meta/media-28.csv",stringsAsFactor=FALSE)
allGenes <- data.frame("nmr.gene_id"=rownames(counts),stringsAsFactors=FALSE)
allGenes <- left_join(allGenes,orthologs)
rownames(allGenes) <- allGenes$nmr.gene_id
allGenes$uniq <- scater::uniquifyFeatureNames(allGenes$nmr.gene_id,allGenes$nmr.symbol)
rownames(allGenes) <- rownames(counts) <- allGenes$uniq

# extraGenes to highlight in the volcano
intgenes <- c("RASG12V","Puro","RasV12","largeT","SV40LT")

# ---- UntransformedVsTransformed(HADI) ----

# DGEL
y <- DGEList(counts=counts,samples=pDat)
y <- y[,y$samples$Source=="Hadi" & y$samples$Transformation!="SV40LT"]
y$samples$Source <- factor(y$samples$Source)
y$samples$Transformation <- factor(y$samples$Transformation)

# Normalize
y <- calcNormFactors(y)

# Design
mdl <- model.matrix(~0+Transformation,data=y$samples)

# estimate Disp
y <- estimateDisp(y,design=mdl)

# fit
fit <- glmQLFit(y,mdl)

cons <- makeContrasts(UnvTr=TransformationSV40LTRasHadi-TransformationUntransformed,levels=mdl)

qlf <- glmTreat(fit, contrast=cons, lfc=1)
tp <- topTags(qlf,n=nrow(y))$table
tp$Gene <- rownames(tp)
int <- tp[1:20,]
int <- rbind(int,tp[setdiff(intgenes,rownames(int)),])
ngenes <- sum(tp$FDR<0.01)

p0 <- ggplot(tp, aes(x=logFC, y=-log10(FDR))) + 
    geom_point(color="grey50",size=0.5) +
    geom_point(data=tp[tp$FDR<0.01,],color="black",size=0.5) +
    geom_hline(yintercept=2,color="red",lty="dashed") +
    geom_vline(xintercept=c(1,-1)) +
    geom_text_repel(data=int,aes(label=Gene),color="blue",force=1) +
    ggtitle("Transformed versus Untransformed (Hadi)")

# ---- Transformed(Hadi)VsTransformed(Tian) ----

# DGEL
y <- DGEList(counts=counts,samples=pDat)
y <- y[,c("H2","A3","B3","A4","B4","C4")]
y$samples$Source <- factor(y$samples$Source)
y$samples$Transformation <- factor(y$samples$Transformation)

# Normalize
y <- calcNormFactors(y)

# Design
mdl <- model.matrix(~0+Source,data=y$samples)

# estimate Disp
y <- estimateDisp(y,design=mdl)

# fit
fit <- glmQLFit(y,mdl)

cons <- makeContrasts(SourceHadi-SourceTian,levels=mdl)

qlf <- glmTreat(fit, contrast=cons, lfc=1)
tp <- topTags(qlf,n=nrow(y))$table
tp$Gene <- rownames(tp)
int <- tp[1:20,]
int <- rbind(int,tp[setdiff(intgenes,rownames(int)),])
ngenes <- c(ngenes,sum(tp$FDR<0.01))

p1 <- ggplot(tp, aes(x=logFC, y=-log10(FDR))) + 
    geom_point(color="grey50",size=0.5) +
    geom_point(data=tp[tp$FDR<0.01,],color="black",size=0.5) +
    geom_hline(yintercept=2,color="red",lty="dashed") +
    geom_vline(xintercept=c(1,-1)) +
    geom_text_repel(data=int,aes(label=Gene),color="blue") +
    ggtitle("Transformed (Hadi) versus Transformed (Tian)")


library(cowplot)
p0 <- p0 + xlim(c(-17,17)) + ylim(c(0,7)) + labs(subtitle=paste0(ngenes[1],"DE Genes"))
p1 <- p1 + xlim(c(-17,17)) + ylim(c(0,7))+ labs(subtitle=paste0(ngenes[2],"DE Genes"))
plot_grid(p0,p1,nrow=1)

