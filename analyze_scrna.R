library(SingleCellTools)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

rna <- readRDS('data/hLIVER_seurat_231212.RDS')

DimPlot(rna, label = T, raster = F)

markers <- FindMarkers(rna, ident.1 = 'Schwann', logfc.threshold = 0.5, only.pos = T)
head(markers %>% arrange(-avg_log2FC))
table(str_split_fixed(colnames(rna), '-', 3)[,3])

Idents(rna) <- paste(rna$condition, rna$celltype)

genes <- read_csv('../../Gene_Lists/updated_tumor_oncogene_human.csv') %>%
  filter(Gene %in% rownames(rna))
genes <- genes[!duplicated(genes$Gene),]
genes[which(genes$Gene=='NFE2L2'),]$Category <- 'Oncogene'
genes[which(genes$Gene=='NOTCH1'),]$Category <- 'Oncogene'
genes <- genes %>%
  arrange(Category, Gene)

DotPlot(rna, features = genes$Gene, #col.min = 0,
        idents = c('NORMAL Hepatocyte', 'MASL Hepatocyte',
                   'MASH Hepatocyte')) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 5, 
                                   color = ifelse(genes$Category == 'Oncogene',
                                                  '#D95F02', '#7570B3'))) +
  labs(x = '', y = '') +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  )


mash_hep <- FindMarkers(rna, ident.1 = 'MASH Hepatocyte',
                        ident.2 = 'NORMAL Hepatocyte', logfc.threshold = 0,
                        features = genes$Gene)
masl_hep <- FindMarkers(rna, ident.1 = 'MASL Hepatocyte',
                        ident.2 = 'NORMAL Hepatocyte', logfc.threshold = 0,
                        features = genes$Gene)
mash_hep$gene <- rownames(mash_hep)
mash_hep <- mash_hep %>%
  filter(p_val_adj < 0.05)
masl_hep$gene <- rownames(masl_hep)
masl_hep <- masl_hep %>%
  filter(p_val_adj < 0.05)

genes_sig <- genes %>%
  filter(Gene %in% unique(c(mash_hep$gene, masl_hep$gene)))
DotPlot(rna, features = genes_sig$Gene, #col.min = 0,
        idents = c('NORMAL Hepatocyte', 'MASL Hepatocyte',
                   'MASH Hepatocyte')) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 5, 
                                   color = ifelse(genes_sig$Category == 'Oncogene',
                                                  '#D95F02', '#7570B3'))) +
  labs(x = '', y = '') +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  )
table(str_split_fixed(colnames(rna), '-', 3)[,3])

rna@meta.data %>% group_by(Patient, condition) %>%
  tally() %>%
  group_by(condition) %>% tally()


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

Idents(rna) <- paste(rna$condition, 
                     paste(rna$Patient, rna$celltype, sep = '_'),
                     sep = '_')
Idents(rna) <- paste(rna$condition, rna$celltype)
counts <- AggregateExpression(rna, assays = 'SCT', return.seurat = T)
sct <- as.matrix(counts@assays$SCT@data)
sct <- sct[c(isgs$Gene[1:10], 'XIST'),]
sct <- sct[,grep('Hepatocyte', colnames(sct))]
sct <- sct[,sort(colnames(sct))]

annot_colors <- setNames(c('#91D6E4', '#F37C79', '#860F0C'),
                         c('NORMAL', 'MASL', 'MASH'))
column_ha <- columnAnnotation(Condition = str_split_fixed(colnames(sct), '-', 2)[,1],
                              col = list(Condition = annot_colors), 
                              show_annotation_name = FALSE,
                              gp = gpar(col = "black"))
cols <- rev(brewer.pal(7, 'RdBu'))
col_fun <- colorRamp2(c(-3, -2, -1, 0, 1, 2, 3), cols)

scaled_mat <- t(scale(t(sct)))
pathway.heat <- Heatmap(scaled_mat,
                        col = col_fun,
                        clustering_distance_rows = "manhattan",
                        clustering_distance_columns = 'manhattan',
                        clustering_method_rows = "ward.D2",
                        clustering_method_columns = 'ward.D2',
                        cluster_rows = T,
                        cluster_columns = T,
                        show_column_names = TRUE,
                        #height = nrow(as.matrix(scaled_mat))*unit(2, "mm"),
                        show_row_names = TRUE,
                        show_row_dend = FALSE,
                        top_annotation = column_ha,
                        #rect_gp = gpar(col = "black", lwd = 1),
                        show_column_dend = FALSE,
                        #left_annotation = row_ha,
                        #column_dend_reorder = c(rep(5,8),rep(1,8)),
                        #heatmap_legend_param = list(direction = 'horizontal'),
                        row_names_gp = gpar(fontsize = 7),
                        name = 'Z-Score')
#pdf(file = '../../plots_for_aaron/snRNAseq_hep_sub_tumor_suppressors_oncogenes_extended2.pdf',
 #   height = 12)
draw(pathway.heat, heatmap_legend_side = 'right', annotation_legend_side = 'right',
     merge_legend = TRUE)
dev.off()



row_colors <- brewer.pal(3, 'Dark2')
row_colors <- setNames(row_colors,
                       c('Interferon', 'Oncogene', 'Tumor Suppressor'))
row_ha <- rowAnnotation(Pathway = some_genes.df$Category,
                        col = list(Pathway = row_colors),
                        show_annotation_name = FALSE,
                        gp = gpar(col = "black"))