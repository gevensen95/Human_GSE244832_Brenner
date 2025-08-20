library(SingleCellTools)

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


head(rna)
#AAACTACCAGAAACCCGAGATA