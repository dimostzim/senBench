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

normalise_counts <- function(counts) {
  filtered_counts <- counts[rowMeans(counts) > MIN_MEAN_COUNT, , drop = FALSE]
  dge_list <- edgeR::calcNormFactors(edgeR::DGEList(filtered_counts), method = "TMM")
  as.data.frame(2^limma::voom(dge_list, plot = FALSE)$E)
}

TOOL <- "markeR"

tabulate_scores <- function(scores_by_method, samples) {
  score_tables <- list()
  for (method in names(scores_by_method)) {
    for (gene_set in names(scores_by_method[[method]])) {
      score_name <- paste(gene_set, method, sep = ":")
      scores <- scores_by_method[[method]][[gene_set]]
      scores <- scores[match(samples, as.character(scores$sample)), ]
      if (anyNA(scores$sample) || nrow(scores) != length(samples) || !all(is.finite(scores$score))) {
        stop("markeR returned no finite score for at least one sample: ", score_name)
      }
      score_tables[[score_name]] <- data.frame(
        sample = samples, tool = TOOL, score_name = score_name, score = scores$score, call = ""
      )
    }
  }
  do.call(rbind, score_tables)
}

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 3)
counts_path <- arguments[1]
annotation_path <- arguments[2]
out_dir <- arguments[3]

counts <- read_counts(counts_path)
symbol_counts <- rowsum(counts, read_symbols_of_gene_ids(rownames(counts), annotation_path))
normalised_counts <- normalise_counts(symbol_counts)

utils::data("genesets_example", package = "markeR")
scores_by_method <- markeR::CalculateScores(
  data = normalised_counts, metadata = NULL, gene_sets = genesets_example, method = "all"
)
scores <- tabulate_scores(scores_by_method, colnames(counts))

stopifnot(dir.create(out_dir))
write.table(scores, file.path(out_dir, "scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
