library("GenomeInfoDb")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("EnhancedVolcano")

################# Load Files #################
csvfile <- file.path(".", "sample_names_MGG31.txt")
sampleTable <- read.csv(csvfile, row.names = 1)

# Generate summarizedExperiment
count_path <-file.path(".","MGG31_simple_counts.txt")
simple_counts <- read.csv(count_path, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
colnames(simple_counts) <- rownames(sampleTable)

# Create DESeq class object
dds <- DESeqDataSetFromMatrix(countData = simple_counts, colData = sampleTable, design = ~ condition)
# Only save non-empty rows
dds <- dds[ rowSums(counts(dds)) > 10, ]
dds$condition <- relevel(dds$condition, ref = "DMSO")

# Run Differential Expression
dds <- DESeq(dds)
go289 <- results(dds, name = "condition_GO289_vs_DMSO")
cx <- results(dds, name = "condition_CX4945_vs_DMSO")
sgc <-results(dds, name = "condition_SGC.CK2.1_vs_DMSO")

goShrunk <- lfcShrink(dds, coef = "condition_GO289_vs_DMSO")
cxShrunk <- lfcShrink(dds, coef = "condition_CX4945_vs_DMSO")
sgcShrunk <- lfcShrink(dds, coef = "condition_SGC.CK2.1_vs_DMSO")
goSig <- subset(goShrunk, padj < 0.05)
cxSig <- subset(cxShrunk, padj < 0.05)
sgcSig <- subset(sgcShrunk, padj < 0.05)

################## Export CSV ################## 
goSigByLFC <- goSig[order(goSig$log2FoldChange),] %>% as.data.frame()
write.csv(goSigByLFC, file = "goDESigbyLFC.csv", quote = F)
goSigName <- as.data.frame(row.names(goSig))
write.csv(goSigName, file = "goDESigName.csv", row.names = F, quote = F)

cxSigByLFC <- cxSig[order(cxSig$log2FoldChange),] %>% as.data.frame()
write.csv(cxSigByLFC, file = "cxDESigByLFC.csv", quote = F)
cxSigName <- as.data.frame(row.names(cxSig))
write.csv(cxSigName, file = "gcxDESigName.csv", row.names = F, quote = F)

sgcSigByLFC <- sgcSig[order(sgcSig$log2FoldChange),] %>% as.data.frame()
write.csv(sgcSigByLFC, file = "sgcDESigByLFC.csv", quote = F)
sgcSigName <- as.data.frame(row.names(sgcSig))
write.csv(sgcSigName, file = "sgcDESigName.csv", row.names = F, quote = F)


# Plot Volcano with EnhancedVolcano" package
EnhancedVolcano(sgcShrunk,
                lab = rownames(sgcShrunk),
                x = 'log2FoldChange',
                y = 'pvalue')

# GSEA Analysis
library(GSEABase)
library(fgsea)
library(tibble)
mes_pro <- read.csv("../MES_vs_proN.csv", header = TRUE, skip = 16)
mes_pro <- mes_pro[mes_pro$log2fc < -2 & mes_pro$corrected.pvalue < 0.001,]
mes_pro <- list(c(mes_pro$hugo))
names(mes_pro) <- "MES vs Pro-Neural"
hallmarks <- gmtPathways("../h.all.v7.4.symbols.gmt")

grpgsea <- data.frame(
  Symbol <- rownames(sgcSigByLFC),
  log2FoldChange <- sgcSigByLFC[,"log2FoldChange"]
)
ranks <- deframe(grpgsea)

gsea_res <- fgsea(mes_pro[1],ranks)
p_value <- gsea_res[gsea_res$pathway == "MES vs Pro-Neural", "padj"]
plotEnrichment(mes_pro[[1]],ranks, ticksSize = 0.5)+
  ggtitle("MES vs Proneural") +
  xlab("Ranked Genes") +
  ylab("Enrichment Score (ES)") +
  theme_minimal(base_size = 28) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(face = "bold"))+
  annotate("text", x = 6000, y = 0.2, label = paste0("padj = ", signif(p_value, 3)), size = 10, color = "black")

gsea_res <- fgsea(hallmarks,ranks)
p_value <- gsea_res[gsea_res$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "padj"]
plotEnrichment(hallmarks$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,ranks, ticksSize = 0.5)+
  ggtitle("HALLMARK_EMT") +
  xlab("Ranked Genes") +
  ylab("Enrichment Score (ES)") +
  theme_minimal(base_size = 28) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(face = "bold"))+
  annotate("text", x = 6000, y = 0.2, label = paste0("padj = ", signif(p_value, 3)), size = 10, color = "black")

# GO enrichment
library(clusterProfiler)

go_names <- rownames(goSig[abs(goSig$log2FoldChange) >0.5,])
cx_names <- rownames(cxSig[abs(cxSig$log2FoldChange) >0.5,])
sgc_names <- rownames(sgcSig[abs(sgcSig$log2FoldChange) >0.5,])

de_names <- sgc_names
ego <- enrichGO(gene = de_names, 
                OrgDb = org.Hs.eg.db, 
                keyType = "SYMBOL", 
                ont = "BP", 
                pAdjustMethod = "BH",
                readable = TRUE)
ego
goplot(ego)

catgorys <- c("response to oxygen levels", "response to decreased oxygen levels","response to hypoxia",
              "T cell differentiation", "extracellular matrix organization", "extracellular structure organization",
              "T cell differentiation in thymus")
dotplot(ego,showCategory = catgorys, font.size = 18, ) +
  theme(
    plot.margin = margin(10, 10, 10, 10), 
    #aspect.ratio = 1.6,
    legend.position = "right",
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black", face = "bold"),
    axis.text.y = element_text(size = 20, color = "black", face = "bold"),
    axis.title.x = element_text(size = 20, color = "black", face = "bold"),
    axis.title.y = element_text(size = 16, color = "black", face = "bold"),
  )


