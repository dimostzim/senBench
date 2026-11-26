# Scores the samples of one count table with the SENCAN classifier (Jochems 2021, Cell Rep), using the copy of
# the classifier object in the hUSI repository.
# Usage: Rscript --vanilla score_sencan.R <counts.tsv.gz> <annotation.tsv.gz> <hUSI repository> <its commit> <out dir>

.libPaths(.Library)  # no library outside the environment

TOOL <- "SENCAN"
SCORE_NAME <- "SENCAN_score"
CLASSIFIER_FILE <- file.path("Data", "SENCAN_classifier.rda")
# The classifier was fitted with the classes FALSE and TRUE; TRUE is "senescent".
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

# Returns the objects of the classifier file in an environment: the reference samples (dge, dmCPM,
# edger_ref_column), the mean and standard deviation of every gene (gene_mean, gene_sd), the genes given to the
# elastic net (use_gene_ix), the fitted elastic net (full_classifier) and the normalisation function.
load_classifier <- function(tool_repository) {
    classifier <- new.env()
    load(file.path(tool_repository, CLASSIFIER_FILE), envir = classifier)
    # The function calls .calcFactorTMMwsp, which is stored in the same file, so it must look it up there.
    environment(classifier$calcNormFactors_WSPref) <- classifier
    classifier
}

# Returns counts with Ensembl gene identifiers as row names. ensembl_gene holds the identifier of every row of
# counts: NA for a gene without a row in the annotation table, "" for a gene without an Ensembl gene identifier.
# Both kinds are left out; the counts of genes that share an identifier are summed.
sum_counts_by_ensembl_gene <- function(counts, ensembl_gene) {
    has_ensembl_gene <- !is.na(ensembl_gene) & ensembl_gene != ""
    rowsum(counts[has_ensembl_gene, , drop = FALSE], ensembl_gene[has_ensembl_gene])
}

# Fills the genes of the classifier that the input lacks: per sample, with the mean count of the two reference
# samples whose counts correlate best (Kendall) with the sample over the genes present. This is the rule of
# cal_SENCAN in R/functions_new.R of the hUSI repository.
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

# Returns the standardised expression of the genes given to the elastic net, samples in rows: normalisation
# against the reference sample, log2 counts per million, the floor of the reference samples, then centring and
# scaling with the mean and standard deviation of the reference samples.
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
stopifnot(length(arguments) == 5)
counts_path <- arguments[1]
annotation_path <- arguments[2]
tool_repository <- arguments[3]
tool_commit <- arguments[4]
out_dir <- arguments[5]

start_time <- proc.time()[["elapsed"]]
counts <- read_counts(counts_path)
annotation <- read_annotation(annotation_path)
classifier <- load_classifier(tool_repository)
suppressPackageStartupMessages(library(glmnet))  # predict and coef of the fitted elastic net

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

coefficients <- coef(classifier$full_classifier, s = "lambda.1se")[-1, 1]
tool_genes <- names(coefficients)[coefficients != 0]
elastic_net_genes <- names(classifier$use_gene_ix)
found_genes <- classifier_genes[gene_found]
run_info <- list(
    tool = TOOL,
    tool_repository = tool_repository,
    tool_commit = tool_commit,
    classifier_file = CLASSIFIER_FILE,
    main_score = SCORE_NAME,
    counts = list(path = normalizePath(counts_path), md5 = unname(tools::md5sum(counts_path))),
    annotation = list(path = normalizePath(annotation_path), md5 = unname(tools::md5sum(annotation_path))),
    samples = ncol(counts),
    genes_in_counts = nrow(counts),
    genes_in_counts_without_annotation = sum(is.na(ensembl_gene)),
    genes_in_counts_without_ensembl_gene = sum(ensembl_gene == "", na.rm = TRUE),
    ensembl_genes_in_counts = nrow(ensembl_counts),
    tool_genes = length(tool_genes),
    tool_genes_found = sum(tool_genes %in% found_genes),
    tool_genes_missing = sum(!tool_genes %in% found_genes),
    tool_genes_missing_ensembl_ids = tool_genes[!tool_genes %in% found_genes],
    tool_genes_missing_share_of_absolute_coefficients = round(
        sum(abs(coefficients[tool_genes[!tool_genes %in% found_genes]])) / sum(abs(coefficients)), 3),
    elastic_net_input_genes = length(elastic_net_genes),
    elastic_net_input_genes_found = sum(elastic_net_genes %in% found_genes),
    elastic_net_input_genes_missing = sum(!elastic_net_genes %in% found_genes),
    normalisation_genes = length(classifier_genes),
    normalisation_genes_found = sum(gene_found),
    normalisation_genes_missing = sum(!gene_found),
    r_version = R.version.string,
    package_versions = lapply(
        c(edgeR = "edgeR", limma = "limma", glmnet = "glmnet", Matrix = "Matrix", pcaPP = "pcaPP",
          jsonlite = "jsonlite"),
        function(package) as.character(packageVersion(package))
    ),
    library_paths = .libPaths(),
    seconds = round(proc.time()[["elapsed"]] - start_time, 1)
)

stopifnot(dir.create(out_dir, recursive = TRUE))
write.table(scores, file.path(out_dir, "scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
jsonlite::write_json(run_info, file.path(out_dir, "run_info.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
