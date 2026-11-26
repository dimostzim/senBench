# markeR: the senescence gene sets that the package ships, scored on bulk RNA-seq counts

markeR (Martins-Silva, Kaizeler and Barbosa-Morais 2026, NAR Genom Bioinform, DOI 10.1093/nargab/lqag057) is an R package that scores gene sets in bulk RNA-seq samples. It is a toolkit for evaluating gene sets, with no trained model and no rule that calls a sample senescent. This wrapper computes, for every sample of one study, the three scores of markeR for each of the three senescence gene sets that the package ships. That gives nine scores per sample.

Tool version: commit `c6528f0fa1570602b077e4688865dcd68af3d058` (2026-07-08, package version 1.3.1) of https://github.com/DiseaseTranscriptomicsLab/markeR, installed unchanged from the clone `/SCRATCH/dtzim01/senescence/tools/markeR/` on biogemt-1. Paths of the repository below refer to this commit unless the branch `paper` (commit `81793db1c2e0c487a1f1bbe3b72ae6491c95e96a`) is named.

## Files

| File | Content |
| --- | --- |
| `run_tool.sh` | The entry point with the fixed interface. |
| `score_gene_sets.R` | The call of markeR and the two output files. |
| `prepare_marker_input.R` | Conversion of the input: raw counts, genes named by symbol, filtered and normalised. The gene sets wrapper (`scripts/tools/genesets/`) uses it too. |
| `environment.yml` | The environment with pinned versions (R 4.5.3 and the packages that markeR imports). |
| `setup_env.sh` | Creates the environment and installs markeR from the clone, through `git archive`, so that the clone is not written to. |
| `conda_explicit_lock.txt`, `r_package_versions.tsv` | The lock: the exact conda packages, and the version of every R package in the environment. |

## How to run

`bash run_tool.sh --counts <counts.tsv.gz> --annotation <annotation.tsv.gz> --out-dir <new folder>`

The wrapper refuses an `--out-dir` that exists. It uses the environment `/SCRATCH/dtzim01/senescence/tools/envs/marker/` on biogemt-1; `--env-prefix <folder>` names another environment made by `setup_env.sh`. The environment was made with:

`bash setup_env.sh /SCRATCH/dtzim01/senescence/tools/envs/env.sh /SCRATCH/dtzim01/senescence/tools/markeR /SCRATCH/dtzim01/senescence/tools/envs/marker`

Output: `scores.tsv` (columns `sample`, `tool`, `score_name`, `score`, `call`; `call` is empty because markeR assigns no class) and `run_info.json` (commit, number of genes at each step, per gene set the genes found, missing and removed by the expression filter, with their names, and the version of every loaded package). The wrapper exits with an error and writes nothing when the count table holds a value that is not a whole number, when a `GeneID` has no symbol in the annotation table, or when markeR returns a score that is not a finite number.

## What markeR needs, with the source of each point

| Need | What the authors state | Source |
| --- | --- | --- |
| Expression matrix | Genes by samples, "filtered and normalised, non log-transformed". markeR applies the log2 itself. | `README.md`, section Input Requirements; `vignettes/markeR.Rmd` |
| Filter | Genes with a mean count of 70 or less in every condition (quiescent, proliferative, senescent) are removed. | `R/data.R` (`?counts_example`); paper, Methods, Preprocessing; branch `paper`, `preprocessing/Filtering_Normalisation_Batch/Preprocessing_SenescenceDatasets.Rmd` |
| Normalisation | Normalisation factors with `edgeR::calcNormFactors` (TMM), then `limma::voom`; the log2 values of `voom` are turned back with `2^`. | Same three sources; the `2^` step is in the file of the branch `paper` |
| Batch correction | In the paper, the dataset was regressed out across the 25 datasets while cell type and cellular state were kept. The package documentation does not ask for it. | Paper, Methods, Batch Effect Correction |
| Gene identifiers | Gene symbols as row names. | `R/data.R`; the gene sets of the package hold gene symbols |
| Gene sets | A named list. A character vector is a gene set without direction; a data frame with a gene column and a column of 1 and -1 is a gene set with direction. | `README.md`; `R/CalculateScores.R` |
| Scoring methods | `logmedian`, `ssGSEA` and `ranking`; the vignette and the paper compute all three (`method = "all"`). No method has a parameter that the user sets. | `R/CalculateScores.R`; `vignettes/markeR.Rmd`; paper, Figure 3 and Supplementary Figure S3 |

The wrapper does the following, per study:

1. It reads the raw counts and names every gene by the `Symbol` of the annotation table.
2. It keeps the genes with a mean raw count above 70 over all samples of the study.
3. It computes TMM normalisation factors (`edgeR::calcNormFactors`, method `TMM`), takes the log2 counts per million of `limma::voom`, and turns them back with `2^`. These are the calls of the file of the branch `paper`.
4. It calls `markeR::CalculateScores(data, metadata = NULL, gene_sets = genesets_example, method = "all")`, where `genesets_example` is the data object of the package.

## Gene sets

The package ships one data object with gene sets, `genesets_example` (`data/genesets_example.rda`), with three gene sets. `score_name` uses the names of the data object.

| Gene set | Genes | Direction | Origin, as the repository states it |
| --- | --- | --- | --- |
| `HernandezSegura` | 55 | 32 up, 23 down | Hernandez-Segura et al. 2017 (Curr Biol, DOI 10.1016/j.cub.2017.07.033) |
| `REACTOME_CELLULAR_SENESCENCE` | 195 | none | The gene set of this name in MSigDB |
| `LiteratureMarkers` | 7 | 5 up, 2 down | "A small, curated gene set of commonly reported senescence markers" (`R/data.R`) |

The paper's case study used nine gene sets (among them SenMayo, CellAge, SeneQuest and CSGene). Seven of them are not in the package, so they are not scored here. The branch `paper` holds a gene set file (`data/SenescenceSigntures_Bidirectional.rds`), which I did not open.

## Scores

`score_name` is `<gene set>:<scoring method>`. The method names are those of the package.

| Method | What the score is (from the code) | Depends on the other samples of the study |
| --- | --- | --- |
| `logmedian` | log2(value + 1); every gene is centred on its median over the samples of the input; the score is the mean over the genes of the set. | Yes, through the median |
| `ssGSEA` | The package's own implementation (`R/ssGSEA_alternative.R`, alpha 0.25), on log2(value): genes are ranked within the sample and a running sum over the ranked list is summed. | Only through the filter |
| `ranking` | Genes are ordered by expression within the sample; the score is the sum of the positions of the genes of the set, divided by the number of genes in the matrix. | Only through the filter |

For a gene set with direction, every method computes the score of the up genes and the score of the down genes separately and subtracts the second from the first (`ssGSEA` weights each by its share of the genes of the set).

Direction: a higher score means more senescent. For `HernandezSegura` and `LiteratureMarkers` this follows from the subtraction. For `REACTOME_CELLULAR_SENESCENCE` the authors state that genes without direction are "assumed to be upregulated upon phenotype" (paper, legend of Figure 1); the vignette reports that this gene set separates senescent from proliferative samples weakly.

Main score: `HernandezSegura:logmedian`. The authors name no single main score. This choice is mine: the paper's main figure of scores (Figure 3A) shows log2-median scores, and of the three shipped gene sets the paper reports the Hernandez-Segura signature as the one that separates senescent samples from both proliferative and quiescent samples.

## Deviations from the authors' use

1. **Filter without conditions.** The authors keep a gene when its mean count is above 70 in at least one condition. The wrapper receives no labels, so it keeps a gene when its mean count is above 70 over all samples of the study. Every gene that the wrapper keeps would also be kept by the authors' rule on the same samples. The reverse does not hold: a gene that is expressed in one condition only can be removed.
2. **One study at a time, no batch correction.** The paper filtered, normalised and batch-corrected 545 samples of 25 datasets together. The wrapper filters and normalises the samples of one study and applies no batch correction (the paper's correction needs the cellular state of every sample).
3. **Read counts.** The paper used kallisto against the RefSeq transcriptome (NCBI release 109). The input here is GEO's uniformly processed count table. The threshold of 70 is applied to these counts unchanged.
4. **Gene symbols.** The symbol is the `Symbol` column of GEO's annotation table. Where several `GeneID`s share one symbol, their counts are added (in the table of 2026-11-15 this concerns one symbol, `TRNAV-CAC`, with three `GeneID`s). Gene set symbols are matched as they are; old symbols are not translated.
5. **Software versions.** The paper used markeR 1.0.0 with limma 3.62.2. The wrapper uses the pinned commit (markeR 1.3.1) with limma 3.68.5 and edgeR 4.10.7 on R 4.5.3. The five files with the scoring code (`R/CalculateScores*.R`, `R/ssGSEA_alternative.R`) are identical between the tag `v1.0.0` and the pinned commit (`git diff`, checked on 2026-11-15).
6. **One thread.** The wrapper sets `OMP_NUM_THREADS` and `OPENBLAS_NUM_THREADS` to 1. The scores of the test study are identical with and without this setting.

## Released behaviour that differs from the documentation (kept as released)

- **`ranking`.** The README and the paper say that the sum of ranks is divided by the number of genes in the set. The code divides by the number of genes in the expression matrix (`R/CalculateScores_Ranking.R`, `length(universe_genes)`). Within one study this changes every sample by the same factor.
- **`ssGSEA`.** The paper describes a modified version of the `gsva` function of GSVA. The code calls its own function (`R/ssGSEA_alternative.R`); the calls of GSVA are commented out. For a gene set with direction, the paper describes the difference of the two scores; the code multiplies each by the share of up or down genes in the whole gene set first.
- **Default method.** The documentation of `CalculateScores` gives `logmedian` as the default; in the code the default is `ssGSEA`. The wrapper passes `method = "all"`, so the default is not used.
- **Names.** `R/data.R` describes the gene sets as `Literature_Senescence` and `REACTOME_Senescence`; in the data object they are `LiteratureMarkers` and `REACTOME_CELLULAR_SENESCENCE`.
- **Missing genes.** Genes of a set that are absent from the matrix are left out of the score without a message.

## Data behind the tool

markeR has no fitted model, so no training data in the strict sense exist. Two links between the tool and public datasets matter for a benchmark:

- **Hernandez-Segura signature.** It was derived from expression data of senescent cells. The markeR repository and paper give only the reference. `kb/identification/tools.tsv` lists the datasets from the Hernandez-Segura paper: E-MTAB-5403, GSE56293, GSE61130, GSE64553, GSE63577, GSE53356 and GSE58910.
- **Datasets of the markeR case study.** The authors report the performance of the gene sets on E-MTAB-5403, E-MTAB-9714, E-MTAB-10969, GSE175533, GSE60340, GSE63577, GSE64553, GSE214410, GSE222400, GSE224070, GSE230181, GSE235768, GSE247831, GSE74324, GSE75643, GSE94928, GSE130727, GSE250224, GSE196724, GSE206677, GSE217718, GSE110268 and GSE213323 (README of the branch `paper`, where GSE213323 is written `200213323`). Nothing was fitted on them. The test study of this wrapper, GSE130727, is one of them.

## Limits

- **No class.** markeR calls no sample senescent; `call` is empty.
- **Scores are relative to the study.** `logmedian` is centred on the median of the samples in the input file. Its value for a sample changes with the other samples of the study, and values of different studies are not on one scale. An input with one sample gives a `logmedian` of exactly 0 for every gene set (checked); a study in which all samples have the same state gives scores that vary around 0.
- **`ranking` across gene sets.** The package documentation states that `ranking` scores of different gene sets cannot be compared.
- **Genes not found.** In GEO's annotation table, two of the 55 genes of `HernandezSegura` have no row under the symbol of the gene set: `FAM214B` and `MT-CYB`. They are left out of the score. `run_info.json` lists the missing genes of every run.
- **The filter depends on sequencing depth.** The threshold is a raw count, so a study with fewer reads per sample keeps fewer genes, and genes of a set can be removed. In the test study, 15,757 of 39,374 gene symbols passed; 7 of the 195 genes of `REACTOME_CELLULAR_SENESCENCE` were removed, none of the other two gene sets.

## Test run

GSE130727 (34 samples), on biogemt-1, 2026-11-15: `/SCRATCH/dtzim01/senescence/runs/20261115-0011-tool-test-marker/` (`output/scores.tsv`, `output/run_info.json`). The run took 7 seconds on one thread. `scores.tsv` has 306 rows: one row for each of the 34 samples and 9 scores; all scores are finite, and every score takes 34 different values. An earlier run with the wrapper before the one-thread setting, `/SCRATCH/dtzim01/senescence/runs/20261115-0009-tool-test-marker/`, gave an identical `scores.tsv`.

Checked in addition: `LiteratureMarkers:logmedian` recomputed by hand from the normalised values agrees with `scores.tsv` (largest difference 5e-15); scoring the samples in reversed column order gives the same score for every sample and score (largest difference 5e-13). The scores were not compared with sample titles or labels.
