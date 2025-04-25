.runBBKNN <- function(pca, block, confounder = NULL,
                      neighbors_within_batch = 3,
                      n_pcs = 50, trim = NULL,
                      computation = "annoy",
                      annoy_n_trees = 10,
                      pynndescent_n_neighbors = 30,
                      pynndescent_random_state = 0,
                      metric = "euclidean",
                      set_op_mix_ratio = 1,
                      local_connectivity = 1,
                      approx = NULL,
                      use_annoy = NULL, use_faiss = NULL,
                      scanpy_logging = False) {
    anndata <- reticulate::import("anndata", convert = FALSE)
    bbknn <- reticulate::import("bbknn", convert = FALSE)
    sc <- reticulate::import("scanpy", convert = FALSE)
    new_dimred(
        t(mnn_res$corrected),
        merge_order = mnn_res$merge.order,
        n_pairs = mnn_res$num.pairs
    )
}
