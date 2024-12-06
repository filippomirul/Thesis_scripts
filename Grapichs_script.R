library(rtracklayer)
library(dplyr)
library(cluster)
library(factoextra)
library(dbscan)
library("DESeq2")
library(ggplot2)
library(biomaRt)
library(tidyr)
library(AnnotationHub)
library(clusterProfiler)
library(enrichplot)


#### Gene DESeq2 ####

gene_matrix <- "C:\\Users\\filoa\\Desktop\\Tesi_Trento\\Scripts_tesi\\File\\Rkat_gene_matrix.tsv"

gene_matrix <- as.data.frame(read.csv(gene_matrix, sep = "\t", na.strings = "NaN"))
gene_condition <- data.frame(Haplotype = c(rep("Haplotype_1", 9),rep("Haplotype_2", 9)),
                             Sample = c(rep("Berry",3),rep("Leaf",3),rep("Tends",3),
                                        rep("Berry",3),rep("Leaf",3),rep("Tends",3)),
                             row.names = c(colnames(gene_matrix)[2:10],colnames(gene_matrix)[13:21]))
gene_count <- gene_matrix[,c(2:10, 13:21)]
rownames(gene_count) <- gene_matrix$Gene_haplotype_1

dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = gene_condition, design = ~ Haplotype)

keep <- rowSums(counts(dds)) >= 30
summary(keep)
dds <- dds[keep,]

dds$Haplotype <- relevel(dds$Haplotype, ref="Haplotype_1")

dds <- DESeq(dds)

res <- results(dds, alpha = 0.01)
r <- as.data.frame(res)
r$Haplotype_2 <- gene_matrix$Gene_haplotype_2[gene_matrix$Gene_haplotype_1 %in% rownames(r)]
r <- r[r$padj <= 0.01,]

#write.table(r[order(r$padj, decreasing = FALSE),],
#           file="C:\\Users\\filoa\\Desktop\\Tesi_Trento\\Scripts_tesi\\File\\DGE.tsv",
#          row.names=T,col.names=T,sep="\t",quote=F)

res
summary(res)
plotMA(res, main= "Differential expressed gene")

rdl <- rlogTransformation(dds, blind = FALSE)
r <- assay(rdl)
hist(assay(rdl))
PCAA <- plotPCA(rdl, intgroup = "Haplotype")
PCAA + geom_text(aes(label = name), size = 3.5, check_overlap = FALSE) + ggtitle("PCA gene expression")

PCAA + 
  geom_point() +  # Add points for PCA if not already included
  geom_text_repel(aes(label = name), size = 3.5) +  # Use geom_text_repel to avoid overlap
  ggtitle("PCA Gene Expression") +
  theme_minimal()


a <- rowSums(gene_count[,1:9])
b <- rowSums(gene_count[,10:18])
a <- data.frame(Haplotype_1 = a, Haplotype_2 =b)
t <- gather(a, key = "Haplotype", value = "Reads_per_gene")

ggplot(data=t,aes(Haplotype, Reads_per_gene)) +
  geom_boxplot(colour="indianred",fill="indianred",alpha=0.7) +
  theme_bw() +
  scale_y_log10()


r <- as.data.frame(res[res$padj < 0.1 & abs(res$log2FoldChange) > 2,])
#&( res$log2FoldChange > 2 | res$log2FoldChange < -2)
names <- rownames(r)
sel_gene_mat <- gene_matrix[gene_matrix$Gene_haplotype_1 %in% names,]

pca <- prcomp(sel_gene_mat[,c(2:10,13:21)], scale = TRUE)
table_pca <- as.data.frame(pca$x)
screeplot(pca)
plot(table_pca$PC1, table_pca$PC2, pch=19)

fviz_nbclust(sel_gene_mat[,c(2:10,13:21)], kmeans, method = "wss")
km <- kmeans(sel_gene_mat[,c(2:10,13:21)], centers = 4, nstart = 10)
fviz_cluster(km, data = sel_gene_mat[,c(2:10,13:21)])


t <- gene_count %>%
  mutate(Leaf_hap_1 = Leaf1.1 + Leaf2.1 + Leaf3.1) %>% 
  mutate(Berry_hap_1 = Berry1.1+ Berry2.1 + Berry3.1) %>% 
  mutate(Tends_hap_1 = Tend1.1 + Tend2.1 + Tend3.1) %>% 
  mutate(Leaf_hap_2 = Leaf1.2 + Leaf2.2 + Leaf3.2) %>% 
  mutate(Berry_hap_2 = Berry1.2+ Berry2.2 + Berry3.2) %>% 
  mutate(Tends_hap_2 = Tend1.2 + Tend2.2 + Tend3.2)
t <- t[,19:24]
t <- gather(t, key = "Sample", value = "Reads_distribution")

ggplot(data=t,aes(Sample,Reads_distribution)) +
  geom_boxplot(colour="indianred",fill="indianred",alpha=0.7) +
  theme_bw() +
  scale_y_log10()

#### Transcrit with singleltons ####

m <- "C:\\Users\\filoa\\Desktop\\Tesi_Trento\\Scripts_tesi\\File\\Matrix_transcript.tsv"
tr <- as.data.frame(read.csv(m, sep = "\t", na.strings = c("NA", "NaN")))

tr$X <- NULL

count_mat <- tr[, 2:19]
count_mat[is.na(count_mat)] <- 0
count_mat <- round(count_mat)
rownames(count_mat) <- ifelse(is.na(tr$Transcript_id_hap1), tr$Transcript_id_hap2, tr$Transcript_id_hap1)
tr[, 20:22]<- NULL


cond <- data.frame(Haplotype = c(rep("Haplotype_1", 9),rep("Haplotype_2", 9)),
                   Sample = c(rep("Berry",3),rep("Leaf",3),rep("Tends",3),
                              rep("Berry",3),rep("Leaf",3),rep("Tends",3)),
                   row.names = c(colnames(tr)[2:19]))

dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = cond, design = ~ Haplotype)

keep <- rowSums(counts(dds)) >= 10
summary(keep)
dds <- dds[keep,]

dds$Haplotype <- relevel(dds$Haplotype, ref="Haplotype_1")

dds <- DESeq(dds)

res <- results(dds, alpha = 0.01)
res_tr_full <- as.data.frame(res)
res_tr_full <- res_tr_full[res_tr_full$padj <= 0.01,]

#write.table(res_tr_full[order(res_tr_full$padj, decreasing = FALSE),],
#            file="C:\\Users\\filoa\\Desktop\\Tesi_Trento\\Scripts_tesi\\File\\DTE.tsv",
#            row.names=T,col.names=T,sep="\t",quote=F)

summary(res)
plotMA(res, main= "Differential expressed transcripts")

rdl <- rlogTransformation(dds, blind = FALSE)
head(assay(rdl))
hist(assay(rdl))
PCAA <- plotPCA(rdl, intgroup = "Haplotype")
PCAA + geom_text(aes(label = name), size = 4.5) + ggtitle("PCA with allele specific transcripts")

PCAA + 
  geom_point() +  # Add points for PCA if not already included
  geom_text_repel(aes(label = name), size = 3.5) +  # Use geom_text_repel to avoid overlap
  ggtitle("PCA transcripts") +
  theme_minimal()

summary(rownames(res_tr_full) %in% rownames(res_tr_link))
# Mode   FALSE    TRUE 
# logical   46942     555 

t <- count_mat %>%
  mutate(Leaf_hap_1 = Leaf1_hap1 + Leaf2_hap1 + Leaf3_hap1) %>% 
  mutate(Berry_hap_1 = Berry1_hap1+ Berry2_hap1 + Berry3_hap1) %>% 
  mutate(Tends_hap_1 = Tend1_hap1 + Tend2_hap1 + Tend3_hap1) %>% 
  mutate(Leaf_hap_2 = Leaf1_hap2 + Leaf2_hap2 + Leaf3_hap2) %>% 
  mutate(Berry_hap_2 = Berry1_hap2+ Berry2_hap2 + Berry3_hap2) %>% 
  mutate(Tends_hap_2 = Tend1_hap2 + Tend2_hap2 + Tend3_hap2)
t <- t[,19:24]
t <- gather(t, key = "sample", value = "read_number")

ggplot(data=t,aes(sample,read_number + 1)) +
  geom_boxplot(colour="indianred",fill="indianred",alpha=0.7) +
  theme_bw() +
  scale_y_log10()


d <-as.data.frame(read.csv("C:\\Users\\filoa\\Desktop\\gene_predictions.gff3", sep = "\t", na.strings = c("NA", "NaN")))
