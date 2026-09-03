library(Seurat)
library(dplyr)
library(data.table)
library(DoubletFinder)
library(BuenColors)

process_1 <- function(id){
  
  # import first pass
  meta_df <- fread(paste0("../output/",id,"_orf_df.tsv")) %>%
    filter(assignable & is_cell) %>%
    mutate(barcodeNew = paste0(barcode, id)) %>% data.frame()
  rownames(meta_df) <- meta_df$barcodeNew
  
  # do an initial seurat stuff for clustering
  gex <- Read10X_h5(paste0("../data/gex/",id,"_counts.h5"))
  colnames(gex) <- paste0(substr(colnames(gex), 1, 16),"-1", id )
  cells_v1 <- intersect(colnames(gex), meta_df$barcodeNew)
  meta_df <- meta_df[cells_v1,]
  
  # Do initial QC and doublet filtering
  so <- CreateSeuratObject(gex[,cells_v1])
  so <- SCTransform(so) %>% RunPCA() %>% RunUMAP(dims = 1:20) 
  so[["SCT"]] <- as(object = so[["SCT"]], Class = "Assay")
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  so_sweep_list <- paramSweep(so, PCs = 1:10, sct = TRUE)
  sweep.stats_so <- summarizeSweep(so_sweep_list, GT = FALSE)
  so_doubs <- find.pK(sweep.stats_so)
  pK <- as.numeric(as.character(so_doubs$pK[which.max(so_doubs$BCmetric)]))
  
  so <- FindNeighbors(so, dims = 1:10)
  so <- FindClusters(so, resolution = 0.1)
  
  # look at the embedding briefly
  annotations   <- so@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  DimPlot(so)
  
  # Do the doublet finder stuff
  doublet_rate <- 0.075
  nExp_poi     <- round(doublet_rate * nrow(so@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  so <- doubletFinder(so, PCs = 1:10, pN = 0.25, pK = pK,
                      nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)
  
  pANN_col <- paste("pANN", 0.25, pK, nExp_poi, sep = "_")
  so <- doubletFinder(so, PCs = 1:10, pN = 0.25, pK = pK,
                      nExp = nExp_poi.adj, reuse.pANN = pANN_col, sct = TRUE)
  colnames(so@meta.data)[dim(so@meta.data)[2]] <- "doublet_finder"
  DimPlot(so, group.by = "doublet_finder")
  cells_v2 <- rownames(so@meta.data %>% filter(doublet_finder == "Singlet"))
  
  # annotate cell types based on marker genes
  anno_df <- data.frame(t(GetAssayData(so)[c("CD3E", "CD3D",   "MPO","CHI3L1", "CD93", "CLEC2B"),])) %>%
    mutate(cluster = so$seurat_clusters) %>%
    group_by(cluster) %>% summarize(mCD3E = mean(CD3E), mCD3D = mean(CD3D),
                                    mMPO = mean(MPO), mCHI = mean(CHI3L1),
                                    mCD93 = mean(CD93), mCLEC2B = mean(CLEC2B)) %>%
    mutate(celltype = case_when(
      mCD3E + mCD3D > 1.5 ~ "Jurkat", 
      mCD93 + mCLEC2B > 2.9 ~ "U937", 
      mMPO + mCHI > 3.5 ~ "HL60",
      TRUE ~ "unknown"
    ))
  print(id)
  print(dim(anno_df))
  print(anno_df)
  
  # Add cell type and finalize
  vec <- anno_df$celltype; names(vec) <- as.character(anno_df$cluster)
  meta_df$celltype <- unname(vec[as.character(so@meta.data$seurat_clusters)])
  so$celltype <- meta_df$celltype
  meta_df2 <- meta_df[cells_v2,]
  
  # final cells
  meta_final <- meta_df2 %>% filter(celltype != "unknown")
  return(list(gex[,rownames(meta_final)], meta_final))
}
list_of_stuff <- lapply(c("anello_unstim", "anello_TNFa", "anello_IFNg", "anello_both"), process_1)

# make the data frame of full metadata
full_meta <- lapply(list_of_stuff, function(i) i[[2]]) %>% rbindlist()
rownames(full_meta) <- full_meta$barcodeNew
full_meta$condition <- stringr::str_split_fixed(full_meta$barcodeNew, "_", 2)[,2]

# make the gex
full_gex <- lapply(list_of_stuff, function(i) i[[1]]) %>% do.call(what = "cbind",.)

# full seurat object
so <- CreateSeuratObject(counts = full_gex, meta.data = full_meta)

# Do standard dimensionality reduction
so <- NormalizeData(so) %>% FindVariableFeatures() %>% ScaleData() %>%
  RunPCA() %>% 
  FindNeighbors() 
so <- so %>% RunUMAP(dims = 1:20) 

p1 <- DimPlot(so, group.by = c("celltype"), shuffle = TRUE) + theme_void() +
  scale_color_manual(values = c("#B87333", "#74c476", "#6F4E37")) +
  theme(legend.position = "none") + ggtitle("")

p2 <- DimPlot(so, group.by = c("condition"), shuffle = TRUE) + theme_void() +
  scale_color_manual(values = c("#3B5B8C", "#2A788E", "#7AD151", "lightgrey"))+ 
  theme(legend.position = "none") + ggtitle("")
p2
p3 <- DimPlot(so, group.by = c("max_id"), shuffle = TRUE) + theme_void() +
  scale_color_manual(values = c("#74c476", jdb_palette("corona"))) +
  theme(legend.position = "none") + ggtitle("")


cowplot::ggsave2(p1, file = "../plots/celltypes_umap.png", width = 4, height = 4, dpi = 400)
cowplot::ggsave2(p2, file = "../plots/stims_umap.png", width = 4, height = 4, dpi = 400)
cowplot::ggsave2(p3, file = "../plots/orfs_umap.png", width = 4, height = 4, dpi = 400)

# Loop over all the things
orfs <- sort(unique(so$max_id))[-1]
stims <- c("both", "IFNg", "TNFa", "unstim")
celltypes <- c("HL60", "Jurkat", "U937")

lapply(orfs, function(orf_1){
  print("----")
  print(orf_1)
  print("----")
  lapply(stims, function(stim_1){
    lapply(celltypes, function(celltype_1){
      so_subset <- subset(so, max_id %in% c(orf_1,"meGFP") & condition == stim_1 & celltype == celltype_1)
      n_orf <- sum(so_subset$max_id != "meGFP")
      diff_df <- so_subset %>% NormalizeData() %>% FindMarkers(ident.1 = orf_1, ident.2 = "meGFP", group.by = "max_id")
      diff_df$orf_1 <- orf_1
      diff_df$stim_1 <- stim_1
      diff_df$celltype_1 <- celltype_1
      diff_df$n_orf <- n_orf
      diff_df$gene <- rownames(diff_df)
      diff_df
    }) %>% rbindlist() 
  })%>% rbindlist() 
})%>% rbindlist() -> full_diff_df

# count number of DEGs
full_diff_df %>% 
  filter(stim_1 == "TNFa") %>% 
  group_by(orf_1, celltype_1) %>% summarize(nDEGs = sum(p_val_adj < 0.1)) %>% arrange(desc(nDEGs)) %>% head(20)


p_ndeg <- full_diff_df %>% 
  group_by( celltype_1, stim_1) %>% summarize(nDEGs = sum(p_val_adj < 0.05)) %>% arrange(desc(nDEGs)) %>%
  ggplot(aes(x = celltype_1, y = stim_1, fill = nDEGs, label = nDEGs)) +
  geom_tile() + 
  geom_text() + coord_flip() + pretty_plot(fontsize = 7) + 
  scale_x_discrete(limits = rev, expand = c(0,0)) + scale_y_discrete(limits = rev, expand = c(0,0)) +
  scale_fill_gradientn(colors = c("white", jdb_palette("solar_rojos")[c(2:9)]))
cowplot::ggsave2(p_ndeg, file = "../plots/grid_count.pdf", width = 2.3, height = 1.5)

order <- c("MZ286102_1_ORF1", "MZ286188_1_ORF1", "OZ258728_1_ORF1", "MZ286258_1_ORF1", "MZ286197_1_ORF1", "OP549859_1_ORF1",
           "MZ286102_1_ORF2", "MZ286188_1_ORF2", "OZ258728_1_ORF2", "MZ286258_1_ORF2", "MZ286197_1_ORF2", "OP549859_1_ORF2")

p1 <- full_diff_df %>% group_by(orf_1) %>%
  filter(stim_1 == "both" & celltype_1 == "U937") %>% summarize(nDEGs = sum(p_val_adj < 0.05)) %>%
  mutate(orf_1 = factor(as.character(orf_1), levels = rev(order))) %>%
  arrange(orf_1) %>%
  mutate(lineage = rep(c("g", "g", "b", "b", "a", "a"), 2)) %>%
  ggplot(aes(x = nDEGs, y = orf_1, fill = lineage))+ 
  geom_bar(stat="identity",color = "black") + pretty_plot(fontsize = 5) + L_border() +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values =c("#984EA3","#7A0403FF","#4777EFFF")) +
  theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/bars_go.pdf", width = 1.6, height = 1.4)


full_diff_df %>%
  filter(orf_1 == "MZ286197_1_ORF2") %>% arrange(p_val_adj) %>% data.frame() %>% head(20)

full_diff_df %>%
  arrange(p_val) %>%
  filter( p_val_adj < 0.05) %>% data.frame()

# make volcano plots
pvolc <- full_diff_df %>%
  filter(pct.1 > 0.05 & pct.2 > 0.01) %>%
  mutate(color = case_when(
    p_val_adj < 0.1 & avg_log2FC > 0 ~ "a_up",
    p_val_adj < 0.1 & avg_log2FC < 0 ~ "b_down",
    TRUE ~ "ns"
  )) %>% 
  filter(orf_1 == "MZ286188_1_ORF2" & celltype_1 == "U937" & stim_1 == "both") %>% 
  ggplot(aes(x = avg_log2FC, y = -1*log10(p_val_adj), label = gene, color = color)) +
  geom_point(size = 0.5) + scale_y_continuous(expand = c(0,0)) + pretty_plot(fontsize = 7) + 
  L_border() + labs(x = "log2FC", y = "-log10p") + 
  scale_color_manual(values = c("dodgerblue2", "firebrick", "black")) +
  theme(legend.position = "none")
cowplot::ggsave2(pvolc, file = "../plots/volcano_MZ286188.pdf", width = 1.7, height = 1.4)

pvolcO <- full_diff_df %>%
  filter(pct.1 > 0.05 & pct.2 > 0.01) %>%
  mutate(color = case_when(
    p_val_adj < 0.1 & avg_log2FC > 0 ~ "a_up",
    p_val_adj < 0.1 & avg_log2FC < 0 ~ "b_down",
    TRUE ~ "ns"
  )) %>% 
  filter(orf_1 == "MZ286102_1_ORF2" & celltype_1 == "U937" & stim_1 == "both") %>% 
  ggplot(aes(x = avg_log2FC, y = -1*log10(p_val_adj), label = gene, color = color)) +
  geom_point(size = 1) + scale_y_continuous(expand = c(0,0)) + pretty_plot(fontsize = 7) + 
  L_border() + labs(x = "log2FC", y = "-log10p") + 
  scale_color_manual(values = c("dodgerblue2", "firebrick", "black")) +
  theme(legend.position = "none")
cowplot::ggsave2(pvolcO, file = "../plots/volcano_MZ286102.pdf", width = 2, height = 1.7)


pvolc2 <- full_diff_df %>%
  filter(pct.1 > 0.05 & pct.2 > 0.01) %>%
  mutate(color = case_when(
    p_val_adj < 0.1 & avg_log2FC > 0 ~ "a_up",
    p_val_adj < 0.1 & avg_log2FC < 0 ~ "b_down",
    TRUE ~ "ns"
  )) %>% 
  filter(orf_1 == "MZ286197_1_ORF2" & celltype_1 == "U937" & stim_1 == "both") %>% 
  ggplot(aes(x = avg_log2FC, y = -1*log10(p_val_adj), label = gene, color = color)) +
  geom_point(size = 1) + scale_y_continuous(expand = c(0,0)) + pretty_plot(fontsize = 7) + 
  L_border() + labs(x = "log2FC", y = "-log10p") + 
  scale_color_manual(values = c("dodgerblue2", "firebrick", "black")) +
  theme(legend.position = "none")
cowplot::ggsave2(pvolc2, file = "../plots/volcano_MZ286197-gamma.pdf", width = 2, height = 1.7)


# now look at pathways
library(msigdbr)
library(fgsea)
library(clusterProfiler)

hm <- msigdbr(species = "Homo sapiens", collection = "H")
t2g <- hm %>% select(gs_name, gene_symbol)

hallmark <- split(hm$gene_symbol, hm$gs_name)
stopifnot(length(hallmark) == 50)

down <- full_diff_df %>% filter(orf_1 == "MZ286188_1_ORF2" & celltype_1 == "U937" & stim_1 == "both") %>%
  filter(avg_log2FC < 0, p_val_adj < 0.1, pct.1 > 0.05, pct.2 > 0.01) %>%
  pull(gene)

enricher(gene = down, TERM2GENE = t2g, universe = full_diff_df$gene,
         pvalueCutoff = 1, qvalueCutoff = 1) %>% data.frame() %>%
  head(5) %>% 
  mutate(ID = gsub("HALLMARK_", "", ID)) %>%
  mutate(ID = factor(as.character(ID), levels = rev(as.character(ID)))) %>%
  ggplot(aes(y = ID, x = -1*log10(qvalue), color = FoldEnrichment)) +
  geom_point(size = 4) + pretty_plot(fontsize = 6) + L_border() +
  scale_x_continuous(limits = c(0, 14)) +
  scale_color_gradientn(colors = c("grey", jdb_palette("solar_rojos")[c(2:9)]), limits = c(0,17)) -> p1
cowplot::ggsave2(p1, file = "../plots/pathways_ORA_188_ORF2-U937-both.pdf", width = 3.4, height = 1.4)
  

#----
### Now do module scores
#----
library(UCell)
tnfa <- intersect(unique(hm$gene_symbol[hm$gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]), rownames(full_gex))
so <- AddModuleScore_UCell(so, features = list(TNFA = tnfa), name = "", ncores = 4)

# handy annotation vectors
genus <- c(MZ286102 = "alpha", MZ286188 = "alpha",
           MZ286258 = "beta",  OZ258728 = "beta",
           MZ286197 = "gamma", OP549859 = "gamma")

cols <- c(meGFP = "black", ORF1 = "grey70",
          alpha = "#984EA3", beta = "#7A0403FF", gamma = "#4777EFFF")

pdf_df <- so@meta.data %>%
  mutate(virus = ifelse(max_id == "meGFP", "meGFP", gsub("_1_ORF[12]", "", max_id)),
         orf   = ifelse(max_id == "meGFP", "meGFP", gsub(".*_1_", "", max_id)),
         fill  = case_when(orf == "meGFP" ~ "meGFP", orf == "ORF1" ~ "ORF1",
                           TRUE ~ unname(genus[virus])),
         lab   = ifelse(max_id == "meGFP", "meGFP", paste(virus, orf)),
         condition = factor(condition, levels = c("unstim", "TNFa", "IFNg", "both")))

lev <- c("meGFP",
         paste(names(genus), "ORF1"),
         paste(names(genus), "ORF2"))
pdf_df$lab <- factor(pdf_df$lab, levels = lev)

violin_1 <- function(df){
  ref <- df %>% filter(max_id == "meGFP") %>%
    group_by(celltype, condition) %>% summarize(m = median(TNFA), .groups = "drop")
  
  ggplot(df, aes(x = lab, y = TNFA, fill = fill)) +
    geom_violin(scale = "width", width = 0.85, linewidth = 0.2, color = "black") +
    geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", linewidth = 0.2, color = "black") +
    scale_fill_manual(values = cols) +
    labs(x = "", y = "UCell TNFA") +
    pretty_plot(fontsize = 7) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(data = ref, aes(yintercept = m), linetype = 2, linewidth = 0.2,
               color = "black", inherit.aes = FALSE) 
}

# all cell types / stims
pall <- violin_1(pdf_df) + facet_grid(celltype ~ condition, scales = "free_y")
cowplot::ggsave2(pall, file = "../plots/violin_ucell_tnfa_all.pdf", width = 7.8, height = 3.9)

# U937 + TNFa alone
pu <- violin_1(pdf_df %>% filter(celltype == "U937" & condition == "TNFa"))+ L_border()
cowplot::ggsave2(pu, file = "../plots/violin_ucell_tnfa_U937_TNFa.pdf", width = 2.6, height = 1.5)

lapply(c("HL60", "U937", "Jurkat"), function(celltype_1){
  lapply(levels(pdf_df$condition), function(stim_1){
    
    ref <- pdf_df %>% filter(celltype == celltype_1, condition == stim_1, max_id == "meGFP") %>% pull(TNFA)
    
    lapply(setdiff(unique(pdf_df$max_id), "meGFP"), function(strain_1){
      x <- pdf_df %>% filter(celltype == celltype_1, condition == stim_1, max_id == strain_1) %>% pull(TNFA)
      if (length(x) < 15) return(NULL)
      
      data.frame(celltype = celltype_1, stim = stim_1, strain = strain_1,
                 virus = gsub("_1_ORF[12]", "", strain_1),
                 orf   = gsub(".*_1_", "", strain_1),
                 n = length(x),
                 delta = mean(x) - mean(ref),
                 p = wilcox.test(x, ref)$p.value)
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist() -> wilcox_df

wilcox_df$fdr <- p.adjust(wilcox_df$p, "BH")

wilcox_df <- wilcox_df %>% group_by(celltype, stim) %>%
  mutate(z_orf1 = (delta - mean(delta[orf == "ORF1"])) / sd(delta[orf == "ORF1"])) %>% ungroup()

wilcox_df %>%
  filter(celltype == "U937" & stim == "TNFa") %>% arrange(z_orf1)
