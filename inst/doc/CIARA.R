## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7, 
  fig.height=5
)



## ----setup--------------------------------------------------------------------
library(CIARA)

required <- c("Seurat")
if (!all(unlist(lapply(required, function(pkg) requireNamespace(pkg, quietly = TRUE)))))
  knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
umap_elmir <- readRDS(system.file("extdata", "annot_umap.rds", package = "CIARA"))

coordinate_umap_human <- umap_elmir[, 2:3]


## ---- eval = FALSE------------------------------------------------------------
#    human_data_seurat <- cluster_analysis_integrate_rare(raw_counts_human_data, "Human_data", 0.1, 5, 30)
#  
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  norm_human_data <- as.matrix(Seurat::GetAssayData(human_data_seurat, slot = "data", assay = "RNA"))
#  
#  knn_human_data <- as.matrix(human_data_seurat@graphs$RNA_nn)
#  
#  
#  
#  

## -----------------------------------------------------------------------------

original_cluster_human <- as.vector(umap_elmir$cluster_id)
names(original_cluster_human) <- names(umap_elmir$cell_name)
plot_umap(coordinate_umap_human, original_cluster_human)

## ---- eval = FALSE------------------------------------------------------------
#  background <- get_background_full(norm_human_data, threshold = 1, n_cells_low = 3, n_cells_high = 20)
#  result <- CIARA(norm_human_data, knn_human_data, background, cores_number = 1, p_value = 0.001, odds_ratio = 2, local_region = 1, approximation = FALSE)

## -----------------------------------------------------------------------------
load(system.file("extdata", "result.Rda", package = "CIARA"))

## -----------------------------------------------------------------------------
ciara_genes <- row.names(result)[result[, 1]<1]
ciara_genes_top <- row.names(result)[order(as.numeric(result[, 1]))]

## ---- eval = FALSE------------------------------------------------------------
#  norm_human_data_ciara <- norm_human_data[ciara_genes, ]

## -----------------------------------------------------------------------------
load(system.file("extdata", "norm_human_data_ciara.Rda", package = "CIARA"))

## -----------------------------------------------------------------------------

p <- list()
for(i in (ciara_genes_top)[1:5]) {
  q <- plot_gene(norm_human_data_ciara, coordinate_umap_human, i, i)
  p <- list(p, q)
}
p

## -----------------------------------------------------------------------------
background <- row.names(result)

## -----------------------------------------------------------------------------
plot_genes_sum(coordinate_umap_human, norm_human_data_ciara, (ciara_genes), "Sum from top CIARA genes")


## -----------------------------------------------------------------------------
norm_counts_small <- apply(norm_human_data_ciara, 1, function(x) {
  y <- x/sum(x)
  return(y)
  })
gene_sum <- apply(norm_counts_small, 1, sum)
    
genes_name_text <- selection_localized_genes(norm_human_data_ciara, ciara_genes, min_number_cells = 4, max_number_genes = 4)
colnames(coordinate_umap_human) <- c("UMAP_1", "UMAP_2")

if ((requireNamespace("plotly", quietly = TRUE))) {
plot_interactive(coordinate_umap_human, gene_sum, genes_name_text, min_x = NULL, max_x = NULL, min_y = NULL, max_y = NULL)
  }

## ---- eval = FALSE------------------------------------------------------------
#  human_data_ciara <- cluster_analysis_integrate_rare(raw_counts_human_data, "Elmir data", 0.01, 5, 30, (ciara_genes))
#  

## ---- eval = FALSE------------------------------------------------------------
#  if ((requireNamespace("clustree", quietly = TRUE))) {
#    find_resolution(human_data_ciara, seq(0.01, 1, 0.1))
#    }

## ---- eval = FALSE------------------------------------------------------------
#  ciara_cluster_human <- human_data_ciara$RNA_snn_res.0.01

## -----------------------------------------------------------------------------
load(system.file("extdata", "ciara_cluster_human.Rda", package = "CIARA"))

## -----------------------------------------------------------------------------
plot_umap(coordinate_umap_human, ciara_cluster_human)
plot_gene(norm_human_data_ciara, coordinate_umap_human, "NANOS3", "NANOS3")
plot_gene(norm_human_data_ciara, coordinate_umap_human, "NANOG", "NANOG")
plot_gene(norm_human_data_ciara, coordinate_umap_human, "DPPA5", "DPPA5")






## -----------------------------------------------------------------------------
final_cluster_human <- merge_cluster(original_cluster_human, ciara_cluster_human, 10)

## -----------------------------------------------------------------------------
final_cluster_human[final_cluster_human == "2-step_2"] <- "PGC"

## -----------------------------------------------------------------------------
plot_umap(coordinate_umap_human, final_cluster_human)

## ---- eval = FALSE------------------------------------------------------------
#  
#  result_test <- test_hvg(raw_counts_human_data, final_cluster_human, (ciara_genes), background, 100, 0.05)
#  

## ---- eval = FALSE------------------------------------------------------------
#  result_test[[2]]

## ---- eval = FALSE------------------------------------------------------------
#  
#  raw_endoderm <- raw_counts_human_data[, as.vector(umap_elmir$cluster_id) == "Endoderm"]
#  raw_hemo <- raw_counts_human_data[, as.vector(umap_elmir$cluster_id) == "Hemogenic Endothelial Progenitors"]
#  raw_exe_meso <- raw_counts_human_data[, as.vector(umap_elmir$cluster_id) == "ExE Mesoderm"]
#  
#  
#  
#  combined_endoderm <- cluster_analysis_sub(raw_endoderm, 0.2, 5, 30, "Endoderm")
#  
#  combined_hemo <- cluster_analysis_sub(raw_hemo, 0.6, 5, 30, "Hemogenic Endothelial Progenitors")
#  
#  combined_exe_meso <- cluster_analysis_sub(raw_exe_meso, 0.5, 5, 30, "ExE Mesoderm")
#  
#  
#  
#  
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  all_sub_cluster <- c(combined_endoderm$seurat_clusters, combined_hemo$seurat_clusters, combined_exe_meso$seurat_clusters)
#  final_cluster_human_version_sub <- merge_cluster(final_cluster_human, all_sub_cluster)

## -----------------------------------------------------------------------------
load(system.file("extdata", "final_cluster_human_version_sub.Rda", package = "CIARA"))

## -----------------------------------------------------------------------------
plot_umap(coordinate_umap_human, final_cluster_human_version_sub)

## -----------------------------------------------------------------------------
table(as.vector(final_cluster_human_version_sub))

## ---- eval = FALSE------------------------------------------------------------
#  Seurat::DefaultAssay(human_data_seurat) <- "RNA"
#  markers_human_final <- markers_cluster_seurat(human_data_seurat, final_cluster_human_version_sub, names(human_data_seurat$RNA_snn_res.0.1), 5)
#  
#  markers_human_top_final <- markers_human_final[[1]]
#  markers_human_all_final <- markers_human_final[[3]]

## ---- eval = FALSE------------------------------------------------------------
#  white_black_markers <- white_black_markers(final_cluster_human_version_sub, "Hemogenic Endothelial Progenitors_4", norm_human_data, markers_human_all_final, 0)
#  sum(white_black_markers)
#  

## ---- eval = FALSE------------------------------------------------------------
#  white_black_markers <- white_black_markers(final_cluster_human_version_sub, "Endoderm_2", norm_human_data, markers_human_all_final, 0)
#  sum(white_black_markers)

## ---- eval = FALSE------------------------------------------------------------
#  white_black_markers <- white_black_markers(final_cluster_human_version_sub, "ExE Mesoderm_0", norm_human_data, markers_human_all_final, 0)
#  sum(white_black_markers)

## ---- eval = FALSE------------------------------------------------------------
#  
#  top_endo <- white_black_markers(final_cluster_human_version_sub, "Endoderm_2", norm_human_data, markers_human_all_final, 0)
#  top_endo <- names(top_endo)[top_endo]
#  
#  
#  mean_top_endo <- apply(norm_human_data[top_endo, final_cluster_human_version_sub == "Endoderm_2"], 1, mean)
#  mean_top_endo <- sort(mean_top_endo, decreasing = T)
#  
#  top_endo <- names(mean_top_endo)
#  names(top_endo) <- rep("Endoderm_2", length(top_endo))
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  top_hemo <- white_black_markers(final_cluster_human_version_sub, "Hemogenic Endothelial Progenitors_4", norm_human_data, markers_human_all_final, 0)
#  top_hemo <- names(top_hemo)[top_hemo]
#  
#  
#  mean_top_hemo <- apply(norm_human_data[top_hemo, final_cluster_human_version_sub == "Hemogenic Endothelial Progenitors_4"], 1, mean)
#  mean_top_hemo <- sort(mean_top_hemo, decreasing = T)
#  
#  top_hemo <- names(mean_top_hemo)
#  names(top_hemo) <- rep("Hemogenic Endothelial Progenitors_4", length(top_hemo))

## ---- eval = FALSE------------------------------------------------------------
#  top_meso <- white_black_markers(final_cluster_human_version_sub, "ExE Mesoderm_1", norm_human_data, markers_human_all_final, 0)
#  top_meso <- names(top_meso)[top_meso]
#  
#  
#  mean_top_meso <- apply(norm_human_data[top_meso, final_cluster_human_version_sub == "ExE Mesoderm_1"], 1, mean)
#  mean_top_meso <- sort(mean_top_meso, decreasing = T)
#  
#  top_meso <- names(mean_top_meso)
#  names(top_meso) <- rep("ExE Mesoderm_1", length(top_meso))
#  
#  
#  

## ---- eval = FALSE------------------------------------------------------------
#  norm_human_data_plot <- norm_human_data[c(top_endo, top_hemo, top_meso), ]

## -----------------------------------------------------------------------------
load(system.file("extdata", "norm_human_data_plot.Rda", package = "CIARA"))
load(system.file("extdata", "top_meso.Rda", package = "CIARA"))
load(system.file("extdata", "top_endo.Rda", package = "CIARA"))
load(system.file("extdata", "top_hemo.Rda", package = "CIARA"))


## -----------------------------------------------------------------------------



toMatch <- c("Endoderm")

plot_balloon_marker(norm_human_data_plot[, grep(paste(toMatch, collapse="|"), final_cluster_human)], final_cluster_human_version_sub[grep(paste(toMatch, collapse="|"), final_cluster_human)], top_endo, 20, max_size=5, text_size=10)


toMatch <- c("Hemogenic Endothelial Progenitors")
plot_balloon_marker(norm_human_data_plot[, grep(paste(toMatch, collapse = "|"), final_cluster_human)], final_cluster_human_version_sub[grep(paste(toMatch, collapse = "|"), final_cluster_human)], top_hemo, 20, max_size = 5, text_size = 8)


toMatch <- c("ExE Mesoderm")
plot_balloon_marker(norm_human_data_plot[, grep(paste(toMatch, collapse = "|"), final_cluster_human)], final_cluster_human_version_sub[grep(paste(toMatch, collapse = "|"), final_cluster_human)], top_meso, length(top_meso), max_size = 5, text_size = 8)




## -----------------------------------------------------------------------------

plot_gene(norm_human_data_ciara, coordinate_umap_human, "C1R", "C1R")
plot_gene(norm_human_data_ciara, coordinate_umap_human, "NANOS3", "NANOS3")
plot_gene(norm_human_data_ciara, coordinate_umap_human, "NANOG", "NANOG")
plot_gene(norm_human_data_ciara, coordinate_umap_human, "SOX17", "SOX17")
plot_gene(norm_human_data_ciara, coordinate_umap_human, "DPPA5", "DPPA5")
plot_gene(norm_human_data_ciara, coordinate_umap_human, "CSF1", "CSF1")

## -----------------------------------------------------------------------------
utils::sessionInfo()

