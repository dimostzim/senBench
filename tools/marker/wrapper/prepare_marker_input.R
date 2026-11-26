# Preparation of the expression matrix that markeR scores: raw counts, genes named by symbol, filtered and normalised.

# A gene is kept when its mean raw count over the samples is above this value (markeR, ?counts_example).
MIN_MEAN_COUNT <- 70

read_counts <- function(counts_path) {
  table <- read.delim(counts_path, check.names = FALSE, quote = "")
  if (colnames(table)[1] != "GeneID") stop("the first column of the count table is not GeneID")
  counts <- as.matrix(table[, -1, drop = FALSE])
  rownames(counts) <- table$GeneID
  if (!is.numeric(counts) || anyNA(counts) || any(counts < 0) || any(counts != round(counts))) {
    stop("the count table holds values that are not whole numbers of reads")
  }
  if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop("the count table holds a repeated GeneID or a repeated sample")
  }
  counts
}

read_symbols_of_gene_ids <- function(gene_ids, annotation_path) {
  annotation <- read.delim(annotation_path, quote = "", colClasses = "character")
  symbols <- annotation$Symbol[match(gene_ids, annotation$GeneID)]
  if (anyNA(symbols) || any(symbols == "")) stop("a GeneID of the count table has no Symbol in the annotation table")
  setNames(symbols, gene_ids)
}

# Filtering and normalisation of the markeR paper (branch paper, Preprocessing_SenescenceDatasets.Rmd),
# without its split of the samples by condition and without its batch correction across studies.
normalise_counts <- function(counts) {
  filtered_counts <- counts[rowMeans(counts) > MIN_MEAN_COUNT, , drop = FALSE]
  dge_list <- edgeR::calcNormFactors(edgeR::DGEList(filtered_counts), method = "TMM")
  as.data.frame(2^limma::voom(dge_list, plot = FALSE)$E)
}
