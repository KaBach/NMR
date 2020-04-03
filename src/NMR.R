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
counts <- read.csv("../data/counts/NMR_Counts.csv",row.names=1)
pDat <- read.csv("../data/meta/SampleInfo.csv")
orthologs <- read.csv("../data/meta/media-28.csv")
pDat <- pDat[pDat$Species=="NMR",]
rownames(pDat) <- pDat$Sample
stopifnot(identical(rownames(pDat), colnames(counts)))

# Move count stats to PD
pDat$Unmapped <- t(counts)[,1]
pDat$MultiMap <- t(counts)[,2]
pDat$NoFeat <- t(counts)[,3]
pDat$Ambiguous <- t(counts)[,4]
counts <- counts[-c(1:4),]
pDat$TotalReads <- colSums(counts)

# Remove lowly expressed
keep <- rowSums(counts)>0
counts <- counts[keep,]

# Convert NMR to Entrez Gene
rmgenes <- read.table("../data/meta/rm_genes.tsv") # These are 24 tRNAs and a few ribosomal RNAs that don't have entrez geneID mapping in the gtf file
counts <- counts[!rownames(counts) %in% rmgenes$V1,]
conv <- read.csv("../data/meta/gene_id_mapping_clean.csv",stringsAsFactors=FALSE,row.names=1) # Conv has more genes, some of them were not mapped to in star
conv <- conv[rownames(counts),]
rownames(counts) <- conv

# Add Gene Annotation from the orthologs
orthologs <- read.csv("../data/meta/media-28.csv",stringsAsFactor=FALSE)
allGenes <- data.frame("nmr.gene_id"=rownames(counts),stringsAsFactors=FALSE)
allGenes <- left_join(allGenes,orthologs)
rownames(allGenes) <- allGenes$nmr.gene_id
allGenes$uniq <- scater::uniquifyFeatureNames(allGenes$nmr.gene_id,allGenes$nmr.symbol)
rownames(allGenes) <- rownames(counts) <- allGenes$uniq

# extraGenes to highlight in the volcano
intgenes <- c("Has2","Hyal2","RASG12V","Puro","RasV12","largeT","SV40LT")

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


# ---- UntransformedVsTransformed(Tian) ----

# DGEL
y <- DGEList(counts=counts,samples=pDat)
y <- y[,y$samples$Source=="Tian" & y$samples$Transformation!="SV40LTRasHadi"]
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

cons <- makeContrasts(UnvTr=TransformationSV40LTRasTian-TransformationUntransformed,levels=mdl)

qlf <- glmTreat(fit, contrast=cons, lfc=1)
tp <- topTags(qlf,n=nrow(y))$table
tp$Gene <- rownames(tp)
int <- tp[1:20,]
int <- rbind(int,tp[setdiff(intgenes,rownames(int)),])
ngenes <- c(ngenes,sum(tp$FDR<0.01))

p2 <- ggplot(tp, aes(x=logFC, y=-log10(FDR))) + 
    geom_point(color="grey50",size=0.5) +
    geom_point(data=tp[tp$FDR<0.01,],color="black",size=0.5) +
    geom_hline(yintercept=2,color="red",lty="dashed") +
    geom_vline(xintercept=c(1,-1)) +
    geom_text_repel(data=int,aes(label=Gene),color="blue") +
    ggtitle("Transformed versus Untransformed (Tian)")


# ---- UntransformedHadivsUntransformedTian
# DGEL
y <- DGEList(counts=counts,samples=pDat)
y <- y[,c("B2","C2","D2","C3","D3","E3")]
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

p3 <- ggplot(tp, aes(x=logFC, y=-log10(FDR))) + 
    geom_point(color="grey50",size=0.5) +
    geom_point(data=tp[tp$FDR<0.01,],color="black",size=0.5) +
    geom_hline(yintercept=2,color="red",lty="dashed") +
    geom_vline(xintercept=c(1,-1)) +
    geom_text_repel(data=int,aes(label=Gene),color="blue") +
    ggtitle("Untransformed (Hadi) versus Untransformed (Tian)")


### Don't we love GO


genesWithOrth <- allGenes[rownames(tp),c("uniq","mouse.gene_id")]
tp.sub <- tp[!is.na(genesWithOrth$mouse.gene_id),]
genesWithOrth <- genesWithOrth[!is.na(genesWithOrth$mouse.gene_id),]
stopifnot(identical(genesWithOrth$uniq,rownames(tp.sub)))
rownames(tp.sub) <- genesWithOrth$mouse.gene_id

sig <- tp.sub[tp.sub$FDR < 0.01,]
deList <- list("Up"=rownames(sig)[sig$logFC>0], "Down"=rownames(sig)[sig$logFC<0])# topGO
univrs <- rownames(tp.sub)
out <- list()

for (i in seq_along(deList)) {
    gens <- deList[[i]]
    subst <- names(deList)[i]
    # Gene universe
    alG <- factor(as.numeric(univrs %in% gens))
    names(alG) <- univrs

    # ---- GOanalysis ----

    library(topGO)
    #prepare Data for topGO
    GO.data <- new("topGOdata", description="Lib GO",ontology="BP", allGenes=alG, 
		   annot=annFUN.org, mapping="org.Mm.eg.db",
		   nodeSize=20, ID="ensembl")
    result.classic <- runTest(GO.data, statistic="fisher")
    res <- GenTable(GO.data, Fisher.classic=result.classic, orderBy="topgoFisher", topNodes=50, numChar=300)

    res$Fisher.classic <- as.numeric(res$Fisher.classic)
    res$Term <- factor(res$Term, levels=res$Term[order(res$Fisher.classic,decreasing=TRUE)])

    out[[subst]] <- res[1:10,] %>%  # this only works because they are ordered by fisher classic
			mutate(hitsPerc=Significant*100/nrow(sig)) %>% 
			ggplot(aes(y=hitsPerc, 
				   x=Term, 
				   fill=Fisher.classic
				   )) +
			    geom_col() +
			    expand_limits(x=0) +
			    labs(y="Hits (%)", x="GO term", colour="p value", size="Count") +
			    scale_fill_viridis(begin=1,end=0) +
			    coord_flip()
}
o0 <- out[[1]] %+% ggtitle(paste(names(out)[1],"in Hadi"))
o1 <- out[[2]] %+% ggtitle(paste(names(out)[2],"in Hadi"))
ggsave("GEOHadi.svg",width=12)

library(cowplot)
p0 <- p0 + xlim(c(-17,17)) + ylim(c(0,7)) + labs(subtitle=paste0(ngenes[1],"DE Genes"))
p0
ggsave("Volcanoe0.svg")
p1 <- p1 + xlim(c(-17,17)) + ylim(c(0,7))+ labs(subtitle=paste0(ngenes[2],"DE Genes"))
p1
ggsave("Volcanoe1.svg")
p2 <- p2 + xlim(c(-17,17)) + ylim(c(0,7))+ labs(subtitle=paste0(ngenes[3],"DE Genes"))
p2
ggsave("Volcanoe2.svg")
p3 <- p3 + xlim(c(-17,17)) + ylim(c(0,7))+ labs(subtitle=paste0(ngenes[4],"DE Genes"))
p3
ggsave("Volcanoe3.svg")

plot_grid(p0,p2,p1,p3,nrow=2)

