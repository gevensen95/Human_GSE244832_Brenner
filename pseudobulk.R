library(SingleCellTools)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

rna <- readRDS('data/hLIVER_seurat_231212.RDS')

ts_onco <- read_csv('../../Gene_Lists/updated_tumor_oncogene_human.csv') %>%
  filter(Gene %in% rownames(rna))
ts_onco <- ts_onco[!duplicated(ts_onco$Gene),]
ts_onco[which(ts_onco$Gene=='NFE2L2'),]$Category <- 'Oncogene'
ts_onco[which(ts_onco$Gene=='NOTCH1'),]$Category <- 'Oncogene'
ts_onco <- ts_onco %>%
  arrange(Category, Gene)
isgs <- read_csv('../../Gene_Lists/human_isgs.csv') %>%
  filter(Gene %in% rownames(rna))
isgs <- isgs[!duplicated(isgs$Gene),]

Idents(rna) <- paste(rna$condition, rna$celltype)

mash_hep <- FindMarkers(rna, ident.1 = 'MASH Hepatocyte',
                        ident.2 = 'NORMAL Hepatocyte', logfc.threshold = 0,
                        features = unique(c(ts_onco$Gene, isgs$Gene)))
masl_hep <- FindMarkers(rna, ident.1 = 'MASL Hepatocyte',
                        ident.2 = 'NORMAL Hepatocyte', logfc.threshold = 0,
                        features = unique(c(ts_onco$Gene, isgs$Gene)))
mash_hep$gene <- rownames(mash_hep)
mash_hep <- mash_hep %>%
  filter(p_val_adj < 0.05)
masl_hep$gene <- rownames(masl_hep)
masl_hep <- masl_hep %>%
  filter(p_val_adj < 0.05)

ts_onco_sig <- ts_onco %>%
  filter(Gene %in% unique(c(mash_hep$gene, masl_hep$gene)))

isgs_sig <- isgs %>%
  filter(Gene %in% unique(c(mash_hep$gene, masl_hep$gene)))


cts <- GetAssayData(rna, assay = "SCT", slot = "counts") # genes x cells (dgCMatrix)
meta <- rna@meta.data

# group id per cell = Sample x CellType
meta$Group <- paste(meta$Patient, 
                    paste(meta$condition, meta$celltype,
                          sep = '_'), sep = '_')

# rowsum over cells per group (on transpose for speed-friendly rowsum)
# result: genes x groups (sums of UMI counts)
pb_mat <- rowsum(t(cts), group = meta$Group) # groups x genes
pb_mat <- t(pb_mat)                           # genes x groups

# Optional: drop empty groups (can happen if some labels are unused)
pb_mat <- pb_mat[, colSums(pb_mat) > 0, drop = FALSE]

lib_size <- colSums(pb_mat)
pb_cpm   <- t(t(pb_mat) / lib_size * 1e6)
pb_log   <- log1p(pb_cpm)  # log(1 + CPM)

# If you have condition/sample metadata, you can carry it into a colData (optional)
suppressPackageStartupMessages(library(DESeq2))
coldata <- data.frame(Group = colnames(pb_mat), row.names = colnames(pb_mat))
dds <- DESeqDataSetFromMatrix(countData = round(pb_mat), colData = coldata, design = ~ 1)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
pb_log <- assay(vsd)   # variance-stabilized counts

genes_to_plot <- isgs$Gene

present <- intersect(rownames(pb_log), genes_to_plot)
missing <- setdiff(genes_to_plot, present)
if (length(missing)) message("Missing genes (not in counts): ", paste(missing, collapse = ", "))

mat_plot <- pb_log[present, , drop = FALSE]
mat_plot <- mat_plot[,grep('Hepatocyte', colnames(mat_plot))]
# center/scale each gene (row)
mat_scaled <- t(scale(t(mat_plot)))  # z-score per row
# keep NAs (flat genes) as 0 for plotting
mat_scaled[is.na(mat_scaled)] <- 0

# Column annotations (optional): split Group back into Sample & CellType for better labels
annot <- NULL
if (grepl("_", colnames(mat_scaled)[1])) {
  annot <- data.frame(
    Condition = str_split_fixed(colnames(mat_scaled), '_', 4)[,3],
    row.names = colnames(mat_scaled)
  )
}

annot$Condition <- factor(annot$Condition, levels = c('NORMAL', 'MASL', 'MASH'))
annot <- annot %>% arrange(Condition)
mat_scaled2 <- mat_scaled[,rownames(annot)]
annot_colors <- setNames(c('#91D6E4', '#F37C79', '#860F0C'),
                         c('NORMAL', 'MASL', 'MASH'))
column_ha <- columnAnnotation(Condition = annot$Condition,
                              col = list(Condition = annot_colors), 
                              show_annotation_name = FALSE,
                              gp = gpar(col = "black"))
#cols <- rev(brewer.pal(7, 'RdBu'))
#col_fun <- colorRamp2(c(-3, -2, -1, 0, 1, 2, 3), cols)

cols <- rev(brewer.pal(5, 'RdBu'))
col_fun <- colorRamp2(c(-2, -1, 0, 1, 2), cols)

pathway.heat <- Heatmap(mat_scaled2,
                        col = col_fun,
                        clustering_distance_rows = "manhattan",
                        clustering_distance_columns = 'manhattan',
                        clustering_method_rows = "ward.D2",
                        clustering_method_columns = 'ward.D2',
                        cluster_rows = T,
                        cluster_columns = F,
                        show_column_names = F,
                        #height = nrow(as.matrix(scaled_mat))*unit(2, "mm"),
                        show_row_names = TRUE,
                        show_row_dend = FALSE,
                        top_annotation = column_ha,
                        #rect_gp = gpar(col = "black", lwd = 1),
                        show_column_dend = FALSE,
                        #left_annotation = row_ha,
                        #column_dend_reorder = c(rep(5,8),rep(1,8)),
                        #heatmap_legend_param = list(direction = 'horizontal'),
                        row_names_gp = gpar(fontsize = 5),
                        name = 'Z-Score')
pdf(file = 'plots/isgs_significant_hepatocytes.pdf')
draw(pathway.heat, heatmap_legend_side = 'right', annotation_legend_side = 'right',
     merge_legend = TRUE)
dev.off()
