# Script to clean up the NMR data

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

# Remove not expressed genes
keep <- rowSums(counts)>0
counts <- counts[keep,]

# Convert NMR to Entrez Gene
rmgenes <- read.table("../data/meta/rm_genes.tsv") # These are 24 tRNAs and a few ribosomal RNAs that don't have entrez geneID mapping in the gtf file
counts <- counts[!rownames(counts) %in% rmgenes$V1,]
conv <- read.csv("../data/meta/gene_id_mapping_clean.csv",stringsAsFactors=FALSE,row.names=1) # Conv has more genes, some of them were not mapped to in star
conv <- conv[rownames(counts),]
rownames(counts) <- conv

write.csv(counts,"../data/counts/NMR_Counts_Clean.csv")
write.csv(pDat,"../data/meta/NMR_ColData.csv")
