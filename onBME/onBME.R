library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("EnhancedVolcano")
library("magrittr")
library("ggpubr")
library("ggrepel")
library("org.Hs.eg.db")

################ Load files ##################
t387 <- read.csv("./T387_counts.txt", header = TRUE, skip = 1, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
mgg31 <- read.csv("./MGG31_counts.txt", header = TRUE, skip = 1, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
all_counts_full <- as.data.frame(c(t387,mgg31[,6:11]))
simple_counts <- as.data.frame(c(t387[,6:11],mgg31[,6:11]))
row.names(simple_counts) <- rownames(t387)
sampleTable_all <- read.csv("sample_names_all_counts.txt", row.names = 1)
colnames(simple_counts) <- rownames(sampleTable_all)

################ Modular analysis ##################
modules <- read.csv("../GBM_modules.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
gene_length <- all_counts_full[,"Length"]
names(gene_length) <- rownames(all_counts_full)
tpm <- simple_counts/gene_length
tpm <- t( t(tpm) * 1e6 / colSums(tpm) )
tpm_log <- log2(tpm+1)
module_score <- matrix(nrow = 12, ncol = 4)
rownames(module_score) <- colnames(tpm_log)
colnames(module_score) <- colnames(modules)[1:4]
for (i in 1:4){
  print(i)
  subtype <- modules[,i][modules[,i] != ""]
  module_score[,i] <- tpm_log[subtype,] %>% colMeans()
}
# 2D representation score calculation
score <- function(mat){
  d <- pmax(mat[,"OPC"], mat[,"NPC"]) - pmax(mat[,"AC"],mat[,"MES"])
  x <- mat[,"OPC"]-mat[,"NPC"]
  x[d<0] <- (mat[,"MES"]-mat[,"AC"])[d<0]
  xy <- data.frame(Column1 = x,
                   Column2 = d)
  colnames(xy) <- c("x","y")
  return(xy)
}

xy_coord <- score(module_score)
xy_coord[,"group"] <- sampleTable_all$condition
ggscatter(
  xy_coord, x = "x", y = "y",
  title = "Scatter Plot",
  xlab = "Relative Meta-module Score", ylab = "Relative Meta-module Score",
  color = "group",
  size = 3,
  xlim = c(-3,3),
  ylim = c(-3.5,3.5)
) +
  scale_x_continuous(limits = c(-3, 3), sec.axis = dup_axis(),breaks = seq(-3, 3, by = 1)) + # Duplicate axes
  scale_y_continuous(limits = c(-4, 4), sec.axis = dup_axis(), breaks = seq(-4, 4, by = 1)) +
  coord_cartesian(clip = "off") + # Prevent clipping
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 1) +  # Horizontal line at y=0
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 1) +
  geom_text_repel(
    data = aggregate(cbind(x, y) ~ group, xy_coord, mean), # Group centroids
    aes(x = x, y = y, label = group),
    color = "black", size = 5, fontface = "bold"
  )+
  theme_minimal()+
  theme(
    axis.ticks.length = unit(-0.1, "cm"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 14, color = "black", face = "bold"),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    axis.title.y = element_text(size = 16, color = "black", face = "bold"),
    axis.title.x.top = element_blank(),
    axis.title.y.right = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )


################ DE analysis #################
## Create DESeq class object
dds <- DESeqDataSetFromMatrix(countData = simple_counts, colData = sampleTable_all, design = ~ condition)
# Only save non-empty rows
dds <- dds[ rowSums(counts(dds)) > 10, ]
# Subset cell line and define control group
dds <- subset(dds,1:6)
dds$condition <- relevel(dds$condition, ref = "T387-2D")
# Run Differential Expression
dds <- DESeq(dds)
bme <- results(dds, name = "condition_T387.onBME_vs_T387.2D")
bmeShrunk <- lfcShrink(dds, coef = "condition_T387.onBME_vs_T387.2D")
bmeSig <- subset(bmeShrunk, padj < 0.05)

#rlog transformation to make the data homoskedastic for PCA, may also use variance stabilizing transformation 
rld <- rlog(dds)
head(assay(rld),3)

# Plot Volcano with EnhancedVolcano" package
EnhancedVolcano(bmeShrunk,
                lab = rownames(bmeShrunk),
                x = 'log2FoldChange',
                y = 'pvalue',
                )

# Plot PCA
plotPCA(rld, intgroup = "condition", ntop = 10000)

# Save csv files
bmeSigByP <- bmeSig[order(bmeSig$padj),]
bmeSigByPDF <- as.data.frame(bmeSigByP)
bmeSigByLFC <- bmeSig[order(bmeSig$log2FoldChange),]
bmeSigByLFC <- as.data.frame(bmeSigByLFC)
bmeSigName <- as.data.frame(row.names(bmeSig))
write.csv(bmeSigName, file = "bmeDESigName.csv", row.names = F, quote = F, col.names = F)
write.csv(bmeSigByP, file = "bmeDESigByP.csv")
write.csv(bmeSigByLFC, file = "bmeDESigbyLFC.csv")

################## GSEA and GO enrichment ###################
library(clusterProfiler)
library(GSEABase)
library(fgsea)
library(tibble)
library(dplyr)
# Prepare gene set objects
hallmarks <- gmtPathways("../h.all.v7.4.symbols.gmt")
mes_pro <- read.csv("../MES_vs_proN.csv", header = TRUE, skip = 16)
mes_pro <- mes_pro[mes_pro$log2fc < -2 & mes_pro$corrected.pvalue < 0.001,]
mes_pro <- list(c(mes_pro$hugo))
names(mes_pro) <- "MES vs Pro-Neural"

# Prepare GSEA input
grpgsea <- data.frame(
  Symbol <- rownames(bmeSigByLFC),
  log2FoldChange <- bmeSigByLFC[,"log2FoldChange"]
)
ranks <- deframe(grpgsea)

# GSEA analysis and plotting
gsea_res <- fgsea(mes_pro[1],ranks)
p_value <- gsea_res[gsea_res$pathway == "MES vs Pro-Neural", "padj"]
plotEnrichment(mes_pro[[1]],ranks, ticksSize = 0.5)+
  ggtitle("MES vs Proneural") +
  xlab("Ranked Genes") +
  ylab("Enrichment Score (ES)") +
  theme_minimal(base_size = 28) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(face = "bold"))+
  annotate("text", x = 9000, y = 0.5, label = paste0("padj = ", signif(p_value, 3)), size = 10, color = "black")

grpGSEAres <- fgsea(pathways = hallmarks, stats = ranks)
pathway_name <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
p_value <- grpGSEAres[grpGSEAres$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "padj"]
plotEnrichment(hallmarks$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,ranks,ticksSize = 0.5)+
  ggtitle("HALLMARK_EMT") +
  xlab("Ranked Genes") +
  ylab("Enrichment Score (ES)") +
  theme_minimal(base_size = 28) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
  annotate("text", x = 9000, y = 0.4, label = paste0("padj = ", signif(p_value, 3)), size = 10, color = "black")


grpGSEAres <- fgsea(pathways = hallmarks, stats = ranks)
grpGSEAresTidy <- grpGSEAres %>% as_tibble() %>% arrange(desc(NES))
grpGSEAresTidy %>% 
  dplyr::select(-leadingEdge) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(grpGSEAresTidy[grpGSEAresTidy$padj <0.05,], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  theme(axis.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size=13,face = "bold"),
        legend.position = "none") +
  scale_fill_manual(values = c("#00CCFF", "#FFCC66")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA")

# GO enrichment using clusterProfiler package

de_names <- read.csv("./bmeDESigName.csv")
de_names <- as.vector(de_names[,1])
de_names <- trimws(de_names)

ego <- enrichGO(gene = de_names, 
                OrgDb = org.Hs.eg.db, 
                keyType = "SYMBOL", 
                ont = "BP", 
                pAdjustMethod = "BH",
                readable = TRUE)
ego
goplot(ego)
barplot(ego, showCategory = 10)
categorys <- c("response to oxygen levels", "response to decreased oxygen levels","response to hypoxia",
              "T cell differentiation", "extracellular matrix organization", "extracellular structure organization",
              "T cell differentiation in thymus")
dotplot(ego,showCategory = categorys, font.size = 18, ) +
  theme(
    plot.margin = margin(10, 10, 10, 10), 
    aspect.ratio = 1.6,
    legend.position = "right",
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black", face = "bold"),
    axis.text.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.title.y = element_text(size = 16, color = "black", face = "bold"),
  )



