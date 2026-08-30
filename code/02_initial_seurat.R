library(Seurat)
library(dplyr)
library(data.table)
library(DoubletFinder)

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
  so <- SCTransform(so) %>% RunPCA() %>% RunUMAP(dims = 1:10) 
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
  doublet_rate <- 0.05 
  nExp_poi     <- round(doublet_rate * nrow(so@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  so <- doubletFinder(so, PCs = 1:10, pN = 0.25, pK = pK,
                      nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)
  
  pANN_col <- paste("pANN", 0.25, pK, nExp_poi, sep = "_")
  so <- doubletFinder(so, PCs = 1:10, pN = 0.25, pK = pK,
                      nExp = nExp_poi.adj, reuse.pANN = pANN_col, sct = TRUE)
  colnames(so@meta.data)[dim(so@meta.data)[2]] <- "doublet_finder"
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
so <- so %>% RunUMAP(dims = 1:30) 

DimPlot(so, group.by = c("condition", "celltype", "max_id"), shuffle = TRUE)

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
full_diff_df %>% group_by(orf_1, stim_1, celltype_1) %>% summarize(nDEGs = sum(p_val_adj < 0.01)) %>% arrange(desc(nDEGs)) %>% head(20)

full_diff_df %>% group_by(orf_1) %>% summarize(nDEGs = sum(p_val_adj < 0.1)) %>% arrange(desc(nDEGs)) %>% head(20)


full_diff_df %>%
  filter(celltype_1 == "HL60" & p_val_adj < 0.01) %>% arrange(p_val_adj) %>% data.frame()

full_diff_df %>%
  filter(celltype_1 == "Jurkat" & p_val_adj < 0.01) %>% arrange(avg_log2FC)

