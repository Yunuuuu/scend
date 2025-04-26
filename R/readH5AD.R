#' Read an H5AD File
#'
#' @param file A string specifying the path to the `.h5ad` file.
#' @param delayed A logical value indicating whether the assay data should be
#' loaded as [`HDF5Array`][HDF5Array::HDF5ArraySeed()] or
#' [`H5SparseMatrix`][HDF5Array::H5SparseMatrix] objects from the **HDF5Array**
#' package.
#' @return A [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment]
#' object.
#' @export
readH5AD <- function(file, delayed = TRUE) {
    # https://github.com/scverse/anndata/blob/main/src/anndata/_io/h5ad.py#L177-L274
    contents <- rhdf5::h5dump(file, load = FALSE)

    # Let's read in the X matrix first... if it's there.
    assays <- list()
    if (!is.null(contents$X)) {
        assays[["X"]] <- h5_read_matrix(
            file,
            "X",
            contents[["X"]],
            delayed = delayed
        )
    }
    for (layer in names(contents[["layers"]])) {
        rlang::try_fetch(
            {
                assays[[layer]] <- h5_read_matrix(
                    file,
                    file.path("layers", layer),
                    contents[["layers"]][[layer]],
                    delayed = delayed
                )
            },
            error = function(cnd) {
                cli::cli_warn(
                    "Setting additional assay from 'layers' failed",
                    parent = cnd
                )
            }
        )
    }
    sce <- SingleCellExperiment::SingleCellExperiment(assays)
    # Adding the various pieces of data.
    rlang::try_fetch(
        {
            coldata <- h5_read_dim_data(file, "obs", contents[["obs"]])
            if (!is.null(coldata)) SummarizedExperiment::colData(sce) <- coldata
        },
        error = function(cnd) {
            cli::cli_warn("Setting 'colData' failed", parent = cnd)
        }
    )

    rlang::try_fetch(
        {
            rowdata <- h5_read_dim_data(file, "var", contents[["var"]])
            if (!is.null(rowdata)) {
                SummarizedExperiment::rowData(sce) <- rowdata
                # Manually set SCE rownames, because setting rowData
                # doesn't seem to set them. (Even tho setting colData
                # does set the colnames)
                rownames(sce) <- rownames(rowdata)
            }
        },
        error = function(cnd) {
            cli::cli_warn("Setting 'rowData' failed", parent = cnd)
        }
    )

    # Adding the reduced dimensions and other bits and pieces.
    rlang::try_fetch(
        {
            SingleCellExperiment::reducedDims(sce) <- h5_read_dim_mats(
                file, "obsm", contents[["obsm"]]
            )
        },
        error = function(cnd) {
            cli::cli_warn("Setting 'reducedDims' failed", parent = cnd)
        }
    )

    rlang::try_fetch(
        {
            rowmat <- h5_read_dim_mats(file, "varm", contents[["varm"]])
            if (length(rowmat)) {
                SummarizedExperiment::rowData(sce) <- cbind(
                    SummarizedExperiment::rowData(sce),
                    rlang::inject(S4Vectors::DataFrame(!!!lapply(rowmat, I)))
                )
            }
        },
        error = function(cnd) {
            cli::cli_warn("Setting 'varm' failed", parent = cnd)
        }
    )

    # Adding pairings, if any exist.
    rlang::try_fetch(
        {
            SingleCellExperiment::rowPairs(sce) <- h5_read_dim_pairs(
                file, "varp", contents[["varp"]]
            )
        },
        error = function(cnd) {
            cli::cli_warn("Setting 'rowPairs' failed", parent = cnd)
        }
    )

    rlang::try_fetch(
        {
            SingleCellExperiment::colPairs(sce) <- h5_read_dim_pairs(
                file, "obsp", contents[["obsp"]]
            )
        },
        error = function(cnd) {
            cli::cli_warn("Setting 'colPairs' failed", parent = cnd)
        }
    )

    if (!is.null(contents$uns)) {
        rlang::try_fetch(
            {
                uns <- h5_read_vec(file, "uns")
                if (!is.null(uns$X_name) &&
                    !is.null(contents$X) &&
                    rlang::is_string(uns$X_name)) {
                    if (SummarizedExperiment::assayNames(sce)[1L] == "X") {
                        SummarizedExperiment::assayNames(sce)[1L] <- uns$X_name
                    }
                    uns$X_name <- NULL
                }
                S4Vectors::metadata(sce) <- uns
            },
            error = function(cnd) {
                cli::cli_warn("Setting 'metadata' failed", parent = cnd)
            }
        )
    }
    sce
}

h5_read_matrix <- function(file, path, fields, delayed) {
    # How to determine if a dataset is a sparse matrix
    if (is.data.frame(fields)) {
        mat <- HDF5Array::HDF5Array(file, path)
    } else {
        mat <- HDF5Array::H5SparseMatrix(file, path)
    }

    if (!delayed) {
        if (DelayedArray::is_sparse(mat)) {
            mat <- methods::as(mat, "sparseMatrix")
        } else {
            mat <- as.matrix(mat)
        }
    }
    mat
}

h5_read_dim_data <- function(file, path, fields) {
    cat_fields <- names(fields[["__categories"]])
    fields <- setdiff(names(fields), "__categories")
    data <- lapply(fields, function(field) {
        o <- h5_read_vec(file, file.path(path, field))
        if (!is.factor(o)) o <- as.vector(o)
        o
    })
    names(data) <- fields

    # for AnnData versions <= 0.7
    for (field in cat_fields) {
        levels <- rhdf5::h5read(file, file.path(path, "__categories", field))
        levels <- as.vector(levels)
        codes <- data[[field]] + 1L
        data[[field]] <- factor(levels[codes], levels = levels)
    }

    ## rhdf5::h5readAttributes(file, "var") |> str()
    ## List of 4
    ##  $ _index          : chr "feature_id"
    ##  $ column-order    : chr [1:4(1d)] "feature_is_filtered" "feature_name" "feature_reference" "feature_biotype"
    ##  $ encoding-type   : chr "dataframe"
    ##  $ encoding-version: chr "0.2.0"
    attributes <- rhdf5::h5readAttributes(file, path)
    index <- attributes[["_index"]]
    if (!is.null(index)) {
        indices <- data[[index]]
    } else {
        indices <- NULL
    }

    colorder <- attributes[["column-order"]]
    if (!is.null(colorder)) data <- data[colorder]

    if (length(data)) {
        out <- rlang::inject(S4Vectors::DataFrame(!!!data))
        rownames(out) <- indices
    } else if (!is.null(indices)) {
        out <- S4Vectors::DataFrame(row.names = indices)
    } else {
        out <- NULL
    }
    out
}

h5_read_vec <- function(file, path) {
    out <- rhdf5::h5read(file, path)
    attrs <- rhdf5::h5readAttributes(file, path)

    # Convert categorical element for AnnData v0.8+
    if (identical(attrs[["encoding-type"]], "categorical") &&
        all(c("codes", "categories") %in% names(out))) {
        codes <- out[["codes"]] + 1
        codes[codes == 0] <- NA
        levels <- out[["categories"]]

        ord <- as.logical(attrs[["ordered"]])

        out <- factor(levels[codes], levels = levels, ordered = ord)
        return(out)
    }

    # Handle booleans. Non-nullable booleans have encoding-type
    # "array", so we have to infer the type from the enum levels
    if (is.factor(out) && identical(levels(out), c("FALSE", "TRUE"))) {
        out <- as.logical(out)
        return(out)
    }

    # Recursively convert element members
    if (is.list(out) && !is.null(names(out))) {
        for (field in names(out)) {
            out[[field]] <- h5_read_vec(file, file.path(path, field))
        }
        names(out) <- make.names(names(out))
    }
    out
}

h5_read_dim_mats <- function(file, path, fields) {
    out <- lapply(fields, function(field) {
        # Because everything's transposed.
        t(rhdf5::h5read(file, file.path(path, field)))
    })
    names(out) <- fields
    out
}

h5_read_dim_pairs <- function(file, path, fields) {
    out <- lapply(fields, function(field) {
        out <- HDF5Array::H5SparseMatrix(file, file.path(path, field))
        as(out, "sparseMatrix")
    })
    names(out) <- fields
    out
}
