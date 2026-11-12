# scend 0.0.0.9000

* Expose options added since scrapper 1.6.0 for variance modelling, PCA,
  feature scoring, and UMAP. Existing handling of older options is preserved.
  `scoreMarkers()` retains its existing arguments and output structure.
* Select the PCA score-centering and UMAP optimization-thread argument names
  using the installed scrapper version (1.7.1 or later), and use
  `rlang::inject()` for dynamic calls.
* Restore ordinary data frames for QC, variance modelling (including per-block
  results), and marker tables with scrapper's S4 DataFrame outputs.
* Fix marker scoring without AUC, numeric group labels, unused factor levels,
  unnamed genes in pairwise tables, and thread forwarding to effect summaries.
* Adapt UMAP seeds and automatic epoch selection, and use the current PCA and
  UMAP argument names when the installed scrapper supports them. Explicit
  t-SNE neighbor counts are passed to scrapper for handling.
* Draw missing seeds from R's RNG so that `set.seed()` controls default runs.
  Explicit seeds now preserve the caller's RNG state.
* Fix feature scoring with default seeds and differently sized feature sets.
  The `weights` attribute is now a list of data frames, one per set,
  containing the input feature row indices and their weights.
* Fix `quickMNN(threads = ...)` and expose the current MNN `n_steps` and
  `merge_policy` options. The obsolete `n_mads`, `robust_iterations`,
  `robust_trim`, `mass_cap`, and `order` arguments have been removed.
  The old `reference_policy` argument has also been removed in favor of
  `merge_policy`.
