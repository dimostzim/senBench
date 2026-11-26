# Score the samples of one study with the senescence gene sets that the markeR package ships.
# Arguments: <counts.tsv.gz> <annotation.tsv.gz> <new output folder> <markeR_source.tsv> <prepare_marker_input.R>

TOOL <- "markeR"

list_gene_set_genes <- function(gene_set) {
  unique(if (is.data.frame(gene_set)) gene_set[[1]] else gene_set)
}

describe_gene_set_coverage <- function(gene_set, input_genes, scored_genes) {
  genes <- list_gene_set_genes(gene_set)
  missing_genes <- setdiff(genes, input_genes)
  removed_genes <- setdiff(intersect(genes, input_genes), scored_genes)
  list(
    genes_in_set = length(genes),
    genes_found_in_input = length(genes) - length(missing_genes),
    genes_missing_from_input = length(missing_genes),
    genes_removed_by_expression_filter = length(removed_genes),
    genes_scored = length(intersect(genes, scored_genes)),
    names_of_genes_missing_from_input = I(missing_genes),
    names_of_genes_removed_by_expression_filter = I(removed_genes)
  )
}

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
stopifnot(length(arguments) == 5)
counts_path <- arguments[1]
annotation_path <- arguments[2]
out_dir <- arguments[3]
tool_source <- read.delim(arguments[4])
source(arguments[5])

counts <- read_counts(counts_path)
# markeR gene sets name genes by symbol; the counts of GeneIDs that share a symbol are added.
symbol_counts <- rowsum(counts, read_symbols_of_gene_ids(rownames(counts), annotation_path))
normalised_counts <- normalise_counts(symbol_counts)

utils::data("genesets_example", package = "markeR")
scores_by_method <- markeR::CalculateScores(
  data = normalised_counts, metadata = NULL, gene_sets = genesets_example, method = "all"
)
scores <- tabulate_scores(scores_by_method, colnames(counts))

loaded_packages <- sort(loadedNamespaces())
run_info <- list(
  tool = TOOL,
  tool_url = tool_source$url,
  tool_commit = tool_source$commit,
  samples = ncol(counts),
  gene_ids_in_input = nrow(counts),
  gene_symbols_in_input = nrow(symbol_counts),
  min_mean_count = MIN_MEAN_COUNT,
  gene_symbols_after_expression_filter = nrow(normalised_counts),
  gene_sets = lapply(
    genesets_example, describe_gene_set_coverage,
    input_genes = rownames(symbol_counts), scored_genes = rownames(normalised_counts)
  ),
  r_version = as.character(getRversion()),
  package_versions = as.list(setNames(
    vapply(loaded_packages, function(package) as.character(packageVersion(package)), character(1)),
    loaded_packages
  ))
)

stopifnot(dir.create(out_dir))
write.table(scores, file.path(out_dir, "scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
jsonlite::write_json(run_info, file.path(out_dir, "run_info.json"), auto_unbox = TRUE, pretty = TRUE)
