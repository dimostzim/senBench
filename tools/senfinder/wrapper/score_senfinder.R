# Classifies the samples of one count table with the released SenFinder script and model.
# Usage: Rscript --vanilla score_senfinder.R <counts.tsv.gz> <annotation.tsv.gz> <SenFinder repository> <its commit>
#            <model_genes.tsv> <set absent genes to zero: TRUE or FALSE> <out dir>

.libPaths(.Library)  # no library outside the environment

TOOL <- "SenFinder"
STATUS_SCORE <- "Senescence_status"
PROBABILITY_SCORE <- "Senescence_probability"
# senFinder.R writes Yes for the class X1 of the model and No for the class X0.
SENESCENT_CLASS <- "X1"
NOT_SENESCENT_CLASS <- "X0"
STATUS_SCORES <- c(No = 0, Yes = 1)
STATUS_CALLS <- c(No = "not senescent", Yes = "senescent")
RELEASED_SCRIPT <- "senFinder.R"
RELEASED_MODEL <- "senFinder.rds"
RELEASED_OUTPUT <- "predicted_results.csv"
MODEL_INPUT <- "senfinder_input.csv"
SAMPLE_COLUMN <- "Sample_id"

read_counts <- function(counts_path) {
    table <- read.delim(counts_path, check.names = FALSE, quote = "")
    if (names(table)[1] != "GeneID") stop("the first column of the count table is not GeneID")
    if (anyDuplicated(table$GeneID)) stop("the count table repeats a GeneID")
    if (anyDuplicated(names(table))) stop("the count table repeats a sample name")
    if (ncol(table) < 2) stop("the count table has no sample")
    counts <- as.matrix(table[, -1, drop = FALSE])
    rownames(counts) <- table$GeneID
    if (!is.numeric(counts) || anyNA(counts) || any(counts < 0) || any(counts != round(counts))) {
        stop("the count table holds values that are not whole numbers of zero or more")
    }
    counts
}

read_annotation <- function(annotation_path) {
    columns <- c("GeneID", "Symbol", "EnsemblGeneID", "Length")
    annotation <- read.delim(annotation_path, check.names = FALSE, quote = "")[, columns]
    if (anyDuplicated(annotation$GeneID)) stop("the annotation table repeats a GeneID")
    annotation
}

# Returns TPM with GeneID as row names. Genes without an annotation row have no length and are left out.
convert_counts_to_tpm <- function(counts, annotation) {
    annotation <- annotation[match(rownames(counts), annotation$GeneID), ]
    keep <- !is.na(annotation$GeneID)
    counts <- counts[keep, , drop = FALSE]
    annotation <- annotation[keep, ]
    if (anyNA(annotation$Length) || any(annotation$Length <= 0)) stop("a gene has no length in the annotation table")
    if (any(colSums(counts) == 0)) stop("a sample has no reads: ", toString(colnames(counts)[colSums(counts) == 0]))
    rate <- counts / annotation$Length
    t(t(rate) / colSums(rate)) * 1e6
}

read_model_genes <- function(model_genes_path, model) {
    model_genes <- read.delim(model_genes_path, quote = "")
    if (!identical(model_genes$model_column, model$coefnames)) {
        stop(model_genes_path, " does not list the columns of the model in the order of the model")
    }
    model_genes
}

# Adds GeneID and matched_by to the table of model genes. A model gene is matched to the annotation row with its
# Ensembl identifier, or, when no row has that identifier, to the row with its symbol. GeneID is NA when the
# annotation table has no such row or the count table has no such GeneID.
find_model_genes <- function(model_genes, annotation, count_gene_ids) {
    find_one <- function(ensembl_gene_id, symbol) {
        by_ensembl <- annotation$GeneID[annotation$EnsemblGeneID == ensembl_gene_id]
        by_symbol <- annotation$GeneID[annotation$Symbol == symbol]
        gene_ids <- if (length(by_ensembl) > 0) by_ensembl else by_symbol
        if (length(gene_ids) > 1) stop("the annotation table has several genes for ", ensembl_gene_id)
        found <- length(gene_ids) == 1 && gene_ids %in% count_gene_ids
        data.frame(
            GeneID = if (found) gene_ids else NA,
            matched_by = if (!found) NA else if (length(by_ensembl) > 0) "EnsemblGeneID" else "Symbol"
        )
    }
    matches <- unname(Map(find_one, model_genes$ensembl_gene_id, model_genes$symbol))
    cbind(model_genes, do.call(rbind, matches))
}

build_model_input <- function(tpm, model_genes) {
    values <- lapply(model_genes$GeneID, function(gene_id) {
        if (is.na(gene_id)) rep(0, ncol(tpm)) else unname(tpm[as.character(gene_id), ])
    })
    names(values) <- model_genes$model_column
    cbind(setNames(data.frame(colnames(tpm)), SAMPLE_COLUMN), as.data.frame(values))
}

# Runs senFinder.R of the repository as its README does. The script writes its result to the working directory.
run_released_script <- function(tool_repository, work_dir) {
    previous_directory <- setwd(work_dir)
    on.exit(setwd(previous_directory))
    status <- system2(
        file.path(R.home("bin"), "Rscript"),
        c("--vanilla", shQuote(file.path(tool_repository, RELEASED_SCRIPT)),
          shQuote(file.path(tool_repository, RELEASED_MODEL)), shQuote(MODEL_INPUT))
    )
    if (status != 0) stop(RELEASED_SCRIPT, " ended with exit code ", status)
    read.csv(RELEASED_OUTPUT, colClasses = "character")
}

# The class probabilities of the model, for the features as senFinder.R reads them.
predict_class_probabilities <- function(model, model_input_path) {
    new_data <- read.csv(model_input_path)
    features <- new_data[, -which(names(new_data) == SAMPLE_COLUMN)]
    predict(model, newdata = features, type = "prob")
}

score_samples <- function(counts_path, annotation_path, tool_repository, tool_commit, model_genes_path,
                          set_absent_genes_to_zero, out_dir) {
    start_time <- proc.time()[["elapsed"]]
    counts <- read_counts(counts_path)
    annotation <- read_annotation(annotation_path)
    tpm <- convert_counts_to_tpm(counts, annotation)

    tool_repository <- normalizePath(tool_repository)
    model_path <- file.path(tool_repository, RELEASED_MODEL)
    model <- readRDS(model_path)
    model_genes <- find_model_genes(read_model_genes(model_genes_path, model), annotation, rownames(tpm))
    absent <- is.na(model_genes$GeneID)
    if (any(absent) && !set_absent_genes_to_zero) {
        stop("the count table has no gene for ", sum(absent), " of the ", nrow(model_genes),
             " genes of the model, and the model needs a value for each of them: ",
             toString(sprintf("%s (%s)", model_genes$model_column[absent], model_genes$ensembl_gene_id[absent])))
    }

    work_dir <- tempfile("senfinder_")
    stopifnot(dir.create(work_dir))
    on.exit(unlink(work_dir, recursive = TRUE))
    model_input_path <- file.path(work_dir, MODEL_INPUT)
    write.csv(build_model_input(tpm, model_genes), model_input_path, row.names = FALSE)

    status <- run_released_script(tool_repository, work_dir)
    probabilities <- predict_class_probabilities(model, model_input_path)
    samples <- colnames(counts)
    stopifnot(
        identical(names(status), c(SAMPLE_COLUMN, STATUS_SCORE)),
        identical(status[[SAMPLE_COLUMN]], samples),
        status[[STATUS_SCORE]] %in% names(STATUS_SCORES),
        nrow(probabilities) == length(samples),
        all(is.finite(probabilities[[SENESCENT_CLASS]])),
        # caret calls the class with the higher probability, so the two outputs must agree.
        identical(status[[STATUS_SCORE]] == "Yes",
                  probabilities[[SENESCENT_CLASS]] > probabilities[[NOT_SENESCENT_CLASS]])
    )

    scores <- rbind(
        data.frame(sample = samples, tool = TOOL, score_name = STATUS_SCORE,
                   score = unname(STATUS_SCORES[status[[STATUS_SCORE]]]),
                   call = unname(STATUS_CALLS[status[[STATUS_SCORE]]])),
        data.frame(sample = samples, tool = TOOL, score_name = PROBABILITY_SCORE,
                   score = probabilities[[SENESCENT_CLASS]], call = "")
    )
    packages <- c("caret", "kernlab", "e1071", "dplyr", "tidyr", "ggplot2", "lattice", "jsonlite")
    run_info <- list(
        tool = TOOL,
        tool_repository = tool_repository,
        tool_commit = tool_commit,
        model = list(path = model_path, md5 = unname(tools::md5sum(model_path))),
        main_score = STATUS_SCORE,
        counts = list(path = normalizePath(counts_path), md5 = unname(tools::md5sum(counts_path))),
        annotation = list(path = normalizePath(annotation_path), md5 = unname(tools::md5sum(annotation_path))),
        samples = ncol(counts),
        genes_in_counts = nrow(counts),
        genes_in_counts_without_annotation = nrow(counts) - nrow(tpm),
        genes_in_tpm = nrow(tpm),
        tool_genes = nrow(model_genes),
        tool_genes_found = sum(!absent),
        tool_genes_missing = sum(absent),
        tool_genes_set_to_zero = model_genes$model_column[absent],
        tool_gene_table = model_genes,
        r_version = R.version.string,
        package_versions = lapply(setNames(packages, packages), function(package) as.character(packageVersion(package))),
        library_paths = .libPaths(),
        seconds = round(proc.time()[["elapsed"]] - start_time, 1)
    )

    stopifnot(dir.create(out_dir, recursive = TRUE))
    write.table(scores, file.path(out_dir, "scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    stopifnot(file.copy(file.path(work_dir, c(MODEL_INPUT, RELEASED_OUTPUT)), out_dir))
    jsonlite::write_json(run_info, file.path(out_dir, "run_info.json"), auto_unbox = TRUE, pretty = TRUE,
                         na = "null")
}

suppressPackageStartupMessages(library(caret))  # the class of the model and its predict method

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 7, arguments[6] %in% c("TRUE", "FALSE"))
score_samples(
    counts_path = arguments[1], annotation_path = arguments[2], tool_repository = arguments[3],
    tool_commit = arguments[4], model_genes_path = arguments[5],
    set_absent_genes_to_zero = as.logical(arguments[6]), out_dir = arguments[7]
)
