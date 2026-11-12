compat_data <- function() {
    old <- get_seed()
    on.exit(restore_seed(old))
    set.seed(101)
    counts <- matrix(rpois(400L * 60L, 5), nrow = 400L,
        dimnames = list(paste0("gene", seq_len(400)), paste0("cell", seq_len(60))))
    list(
        counts = counts,
        logcounts = log2(counts + 1),
        groups = rep(c("A", "B", "C"), each = 20L),
        block = rep(c("batch1", "batch2"), each = 30L),
        pcs = matrix(rnorm(10L * 60L), nrow = 10L,
            dimnames = list(NULL, colnames(counts)))
    )
}
