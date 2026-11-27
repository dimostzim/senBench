.libPaths(.Library)

TOOL <- "SENCAN"
SCORE_NAME <- "SENCAN_score"
CLASSIFIER_FILE <- file.path("Data", "SENCAN_classifier.rda")
CLASS_NAMES <- c("FALSE" = "not senescent", "TRUE" = "senescent")

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
    counts
}

read_annotation <- function(annotation_path) {
    annotation <- read.delim(annotation_path, check.names = FALSE, quote = "")[, c("GeneID", "EnsemblGeneID")]
    if (anyDuplicated(annotation$GeneID)) stop("the annotation table repeats a GeneID")
    annotation
}

load_classifier <- function(tool_repository) {
    classifier <- new.env()
    load(file.path(tool_repository, CLASSIFIER_FILE), envir = classifier)
    environment(classifier$calcNormFactors_WSPref) <- classifier
    classifier
}

sum_counts_by_ensembl_gene <- function(counts, ensembl_gene) {
    has_ensembl_gene <- !is.na(ensembl_gene) & ensembl_gene != ""
    rowsum(counts[has_ensembl_gene, , drop = FALSE], ensembl_gene[has_ensembl_gene])
}

impute_missing_genes <- function(classifier_counts, reference_counts) {
    missing <- is.na(classifier_counts[, 1])
    if (!any(missing)) return(classifier_counts)
    for (sample in colnames(classifier_counts)) {
        correlations <- apply(reference_counts[!missing, ], 2, pcaPP::cor.fk, x = classifier_counts[!missing, sample])
        nearest <- tail(order(correlations), n = 2)
        classifier_counts[missing, sample] <- rowMeans(reference_counts[missing, nearest, drop = FALSE])
    }
    classifier_counts
}

standardize_expression <- function(classifier_counts, classifier) {
    dge <- edgeR::DGEList(counts = classifier_counts, genes = data.frame(ensembl_id = rownames(classifier_counts)))
    dge <- classifier$calcNormFactors_WSPref(dge, classifier$edger_ref_column)
    log_cpm <- edgeR::cpm(dge, log = TRUE)
    lowest_reference_value <- min(classifier$dmCPM)
    log_cpm[log_cpm < lowest_reference_value] <- lowest_reference_value
    standardized <- (log_cpm - classifier$gene_mean) / classifier$gene_sd
    t(standardized[classifier$use_gene_ix, , drop = FALSE])
}

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 4)
counts_path <- arguments[1]
annotation_path <- arguments[2]
tool_repository <- arguments[3]
out_dir <- arguments[4]

counts <- read_counts(counts_path)
annotation <- read_annotation(annotation_path)
classifier <- load_classifier(tool_repository)
suppressPackageStartupMessages(library(glmnet))

ensembl_gene <- annotation$EnsemblGeneID[match(rownames(counts), annotation$GeneID)]
ensembl_counts <- sum_counts_by_ensembl_gene(counts, ensembl_gene)
classifier_genes <- classifier$dge$genes$ensembl_id
stopifnot(
    identical(rownames(classifier$dmCPM), classifier_genes),
    identical(names(classifier$gene_mean), classifier_genes),
    identical(names(classifier$gene_sd), classifier_genes)
)
classifier_counts <- ensembl_counts[match(classifier_genes, rownames(ensembl_counts)), , drop = FALSE]
rownames(classifier_counts) <- classifier_genes
gene_found <- !is.na(classifier_counts[, 1])
if (any(colSums(classifier_counts[gene_found, , drop = FALSE]) == 0)) {
    stop("a sample has no reads in the genes of the classifier")
}
classifier_counts <- impute_missing_genes(classifier_counts, classifier$dge$counts)

standardized <- standardize_expression(classifier_counts, classifier)
score <- as.vector(predict(classifier$full_classifier, standardized, s = "lambda.1se", type = "response"))
class <- as.vector(predict(classifier$full_classifier, standardized, s = "lambda.1se", type = "class"))
if (!all(is.finite(score))) stop("SENCAN returned a score that is not a number")
stopifnot(identical(rownames(standardized), colnames(counts)), class %in% names(CLASS_NAMES))

scores <- data.frame(sample = rownames(standardized), tool = TOOL, score_name = SCORE_NAME, score = score,
                     call = unname(CLASS_NAMES[class]))

stopifnot(dir.create(out_dir, recursive = TRUE))
write.table(scores, file.path(out_dir, "scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
