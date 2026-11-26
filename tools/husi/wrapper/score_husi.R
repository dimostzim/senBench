# Scores the samples of one count table with hUSI, using the functions and gene weights of the hUSI repository.
# Usage: Rscript --vanilla score_husi.R <counts.tsv.gz> <annotation.tsv.gz> <hUSI repository> <its commit> <out dir>

.libPaths(.Library)  # no library outside the environment

TOOL <- "hUSI"
SCALED_SCORE <- "hUSI"
UNSCALED_SCORE <- "hUSI_unscaled"
# SSE_hUSI returns 1 for the samples above its threshold and 0 for the others.
SSE_CLASS_NAMES <- c("0" = "not senescent", "1" = "senescent")

read_counts <- function(counts_path) {
    table <- read.delim(counts_path, check.names = FALSE, quote = "")
    if (names(table)[1] != "GeneID") stop("the first column of the count table is not GeneID")
    if (anyDuplicated(table$GeneID)) stop("the count table repeats a GeneID")
    if (anyDuplicated(names(table))) stop("the count table repeats a sample name")
    counts <- as.matrix(table[, -1, drop = FALSE])
    rownames(counts) <- table$GeneID
    if (!is.numeric(counts) || anyNA(counts) || any(counts < 0) || any(counts != round(counts))) {
        stop("the count table holds values that are not whole numbers of zero or more")
    }
    if (ncol(counts) < 2) stop("hUSI is scaled between the lowest and highest sample: at least two samples are needed")
    counts
}

read_annotation <- function(annotation_path) {
    annotation <- read.delim(annotation_path, check.names = FALSE, quote = "")[, c("GeneID", "Symbol", "Length")]
    if (anyDuplicated(annotation$GeneID)) stop("the annotation table repeats a GeneID")
    annotation
}

# Returns TPM with gene symbols as row names. Genes without an annotation row are left out, and of the genes
# that share a symbol the first is kept, as ID2symbol of the hUSI repository does (R/functions.R).
convert_counts_to_symbol_tpm <- function(counts, annotation) {
    annotation <- annotation[match(rownames(counts), annotation$GeneID), ]
    keep <- !is.na(annotation$GeneID) & !duplicated(annotation$Symbol)
    counts <- counts[keep, , drop = FALSE]
    annotation <- annotation[keep, ]
    if (anyNA(annotation$Length) || any(annotation$Length <= 0)) stop("a gene has no length in the annotation table")
    if (any(colSums(counts) == 0)) stop("a sample has no reads: ", toString(colnames(counts)[colSums(counts) == 0]))
    rate <- counts / annotation$Length
    tpm <- t(t(rate) / colSums(rate)) * 1e6
    rownames(tpm) <- annotation$Symbol
    tpm
}

# sc/hUSI.R loads the gene weights from a path relative to the repository and defines cal_hUSI, minmax and SSE_hUSI.
source_husi_functions <- function(tool_repository) {
    previous_directory <- setwd(tool_repository)
    on.exit(setwd(previous_directory))
    source("sc/hUSI.R")
}

# scoreOCLR is the function of the authors' analysis scripts that returns the correlation before scaling. Only
# its definition is evaluated, because the rest of R/functions_new.R loads packages that scoring does not need.
define_scoreOCLR <- function(tool_repository) {
    expressions <- parse(file.path(tool_repository, "R", "functions_new.R"), keep.source = FALSE)
    is_definition <- vapply(expressions, function(expression) {
        is.call(expression) && identical(expression[[1]], as.name("<-")) &&
            identical(expression[[2]], as.name("scoreOCLR"))
    }, logical(1))
    stopifnot(sum(is_definition) == 1)
    eval(expressions[[which(is_definition)]], globalenv())
}

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 5)
counts_path <- arguments[1]
annotation_path <- arguments[2]
tool_repository <- arguments[3]
tool_commit <- arguments[4]
out_dir <- arguments[5]

start_time <- proc.time()[["elapsed"]]
counts <- read_counts(counts_path)
annotation <- read_annotation(annotation_path)
tpm <- convert_counts_to_symbol_tpm(counts, annotation)

source_husi_functions(tool_repository)
define_scoreOCLR(tool_repository)

husi <- cal_hUSI(tpm)
husi_unscaled <- scoreOCLR(profile = tpm, model = SenOCLR, m = "spearman")[[1]]
if (!all(is.finite(husi_unscaled))) stop("hUSI returned a correlation that is not a number")
if (!all(is.finite(husi))) stop("hUSI cannot be scaled: all samples have the same correlation")
stopifnot(
    identical(names(husi), colnames(counts)),
    identical(names(husi_unscaled), colnames(counts)),
    isTRUE(all.equal(minmax(husi_unscaled), husi))
)
sse_class <- SSE_hUSI(husi)
stopifnot(identical(names(sse_class), colnames(counts)), as.character(sse_class) %in% names(SSE_CLASS_NAMES))

scores <- rbind(
    data.frame(sample = names(husi), tool = TOOL, score_name = SCALED_SCORE, score = unname(husi),
               call = unname(SSE_CLASS_NAMES[as.character(sse_class)])),
    data.frame(sample = names(husi_unscaled), tool = TOOL, score_name = UNSCALED_SCORE,
               score = unname(husi_unscaled), call = "")
)

tool_genes <- names(SenOCLR$w)
genes_without_annotation <- sum(!rownames(counts) %in% annotation$GeneID)
run_info <- list(
    tool = TOOL,
    tool_repository = tool_repository,
    tool_commit = tool_commit,
    main_score = SCALED_SCORE,
    counts = list(path = normalizePath(counts_path), md5 = unname(tools::md5sum(counts_path))),
    annotation = list(path = normalizePath(annotation_path), md5 = unname(tools::md5sum(annotation_path))),
    samples = ncol(counts),
    genes_in_counts = nrow(counts),
    genes_in_counts_without_annotation = genes_without_annotation,
    genes_in_counts_with_repeated_symbol = nrow(counts) - genes_without_annotation - nrow(tpm),
    genes_in_tpm = nrow(tpm),
    tool_genes = length(tool_genes),
    tool_genes_found = sum(tool_genes %in% rownames(tpm)),
    tool_genes_missing = sum(!tool_genes %in% rownames(tpm)),
    r_version = R.version.string,
    package_versions = lapply(
        c(dplyr = "dplyr", mclust = "mclust", jsonlite = "jsonlite"),
        function(package) as.character(packageVersion(package))
    ),
    library_paths = .libPaths(),
    seconds = round(proc.time()[["elapsed"]] - start_time, 1)
)

stopifnot(dir.create(out_dir, recursive = TRUE))
write.table(scores, file.path(out_dir, "scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
jsonlite::write_json(run_info, file.path(out_dir, "run_info.json"), auto_unbox = TRUE, pretty = TRUE)
