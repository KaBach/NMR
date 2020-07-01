# this file only exists because I am too lazy to use bash in the first place
map <- read.table(file="../data/meta/gene_id_mapping.tsv",sep=" ",stringsAsFactors=FALSE)
map$V1 <- gsub(";","",map$V1) # remove trailing ;
map$V2 <- unlist(lapply(strsplit(map$V2,","), function(x) x[[1]])) # remove GenBank bit
map$V2 <- gsub("GeneID:|Genbank:|;","",map$V2)
#Drop the 13 mitochondrial genes that exist twice
map <- map[!grepl("YP_",map$V2),]
# remove the duplicated entries
map <- map[!duplicated(map$V1),]
write.csv(map,file="../data/meta/gene_id_mapping_clean.csv",row.names=FALSE,quote=FALSE)
