.libPaths(.Library)

TOOL <- "hUSI"
SCALED_SCORE <- "hUSI"
UNSCALED_SCORE <- "hUSI_unscaled"
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

source_husi_functions <- function(tool_repository) {
    previous_directory <- setwd(tool_repository)
    on.exit(setwd(previous_directory))
    source("sc/hUSI.R")
}

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
stopifnot(length(arguments) == 4)
counts_path <- arguments[1]
annotation_path <- arguments[2]
tool_repository <- arguments[3]
out_dir <- arguments[4]

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

stopifnot(dir.create(out_dir, recursive = TRUE))
write.table(scores, file.path(out_dir, "scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
