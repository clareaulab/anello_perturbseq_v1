library(data.table)
library(dplyr)
library(Matrix)
library(matrixStats)
library(Seurat)
library(BuenColors)

import_kb_probes <- function(i){
  genes <- fread(paste0("../data/orf_probes/",i,".genes.txt"), header = FALSE)[[1]]
  barcodes <- fread(paste0("../data/orf_probes/",i,".barcodes.txt"), header = FALSE)[[1]]
  mtx <- fread(paste0("../data/orf_probes/",i,".mtx"), skip = 3)
  
  # create matrix
  sm <- sparseMatrix(
    i = c(mtx[[1]],length(barcodes)),
    j = c(mtx[[2]],length(genes)),
    x = c(mtx[[3]],0)
  )
  rownames(sm) <- paste0(barcodes, "-1")
  colnames(sm) <- genes
  return(data.matrix(sm))
}




# Function to assign each barcode
assign_per_barcode <- function(mat){
 
  # strip trailing "_probeN" to get the parent feature/gene name
  feature <- sub("_genscript", "", sub("_probe[0-9]+$", "", colnames(mat)))
  
  # sum columns sharing the same feature name
  feature_mat <- sapply(unique(feature), function(f) {
    cols <- which(feature == f)
    if (length(cols) == 1) mat[, cols] else rowSums(mat[, cols, drop = FALSE])
  })
  rownames(feature_mat) <- rownames(mat)
  
  summary_df <- data.frame(
    barcode = rownames(feature_mat),
    max_id = colnames(feature_mat)[max.col(feature_mat)],
    n_4 = rowSums(feature_mat >= 4),
    n_10 = rowSums(feature_mat >= 10),
    max_umi = rowMaxs(feature_mat),
    total_umi = rowSums(feature_mat)
  ) %>%
    mutate(assignable = n_4 == 1 | ((max_umi/total_umi) >= 0.75 & max_umi >= 2))
}

process_sample <- function(id){
  sm_orfs <- import_kb_probes(id)
  barcodes <- colnames(Read10X_h5(paste0("../data/gex/",id,"_counts.h5"))) %>% substr(., 1, 16) %>% paste0(., "-1")
  assign_df <- assign_per_barcode(sm_orfs) %>% mutate(is_cell = barcode %in% barcodes)
  
  write.table(assign_df, file = paste0("../output/",id,"_orf_df.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}
process_sample(id = "anello_both")
process_sample(id = "anello_IFNg")
process_sample(id = "anello_TNFa")
process_sample(id = "anello_unstim")
