# hUSI wrapper

hUSI (human universal senescence index) is the tool of Wang et al., Nat Aging 2025 (DOI 10.1038/s43587-025-00886-2, PMC12176644). This wrapper runs the released hUSI code and gene weights, unchanged, on the raw counts of one study of bulk RNA-seq samples.

- Repository: https://github.com/WJPina/HUSI, commit `300da32f111f8c8603e822eca8e84c9b5a76db48` (2025-09-11). The wrapper reads the unchanged clone `/SCRATCH/dtzim01/senescence/tools/HUSI` on biogemt-1 and stops if the clone is at another commit or holds a changed file.
- The same repository holds a copy of the SENCAN classifier. This wrapper does not use it.

## How to run

```
bash run_tool.sh --counts <counts.tsv.gz> --annotation <annotation.tsv.gz> --out-dir <new folder>
```

`--out-dir` must not exist. Two optional arguments name the environment and the clone when they are not at their places on biogemt-1: `--env-prefix` (default `/SCRATCH/dtzim01/senescence/tools/envs/husi`) and `--tool-repository` (default `/SCRATCH/dtzim01/senescence/tools/HUSI`). The folder that holds the environment must also hold `env.sh`.

Files of this folder:

| File | Content |
| --- | --- |
| `run_tool.sh` | Entry point. Checks the arguments and the commit of the clone, then calls `score_husi.R`. |
| `score_husi.R` | Reads the two input files, converts the counts, calls the hUSI functions, writes the output. |
| `setup_env.sh`, `environment.yml`, `install_r_packages.R` | Create the environment. |
| `conda_explicit_lock.txt`, `r_packages_lock.tsv` | The installed conda packages and R packages with their versions. |

## What hUSI needs

| Need | Value | Source |
| --- | --- | --- |
| Input | A matrix of genes by samples with normalised expression values. For bulk RNA-seq the authors use TPM (transcripts per million). | README of the repository ("input normalized gene expression matrix"); paper, Methods, "Validation of hUSI on in vivo datasets" ("using TPM normalized gene expression matrix"); `R/Validation.R` (GTEx TPM, TCGA TPM, GSE113957 FPKM converted to TPM, the authors' own RNA-seq GSE282274 as StringTie TPM) |
| Gene identifiers | Human gene symbols as row names | `sc/hUSI.R`, function `cal_hUSI` (`intersect(names(w), rownames(profile))`) |
| Gene weights | 13,346 genes with one weight each, in `sc/SenOCLR_features.rdata` (object `SenOCLR`, element `w`). The file is identical to `Data/SenOCLR_l2=1_drop.rdata`, the final model of `R/Model.R`. | Repository; md5 of both files compared |
| Genes used for a sample | The genes that are both in the weights and in the input. Genes of the weights that the input lacks are left out; nothing is filled in. | `sc/hUSI.R`, `cal_hUSI` |
| Score | Spearman correlation, over these genes, between the weights and the expression values of the sample. The correlations of all samples of the input are then scaled so that the lowest becomes 0 and the highest 1. | `sc/hUSI.R`, `cal_hUSI` and `minmax`; paper, Methods, "Model training and hyperparameter optimization" |
| Parameters | None at scoring. The L2 penalty of 1 and the choice of Spearman correlation are fixed in the released weights and code. | Paper, same section |
| Class | Optional. `SSE_hUSI` splits the samples of the input into two groups at the hUSI value that gives the smallest sum of squared errors within the two groups. `GMM_hUSI` fits a Gaussian mixture and returns numbered states. The authors recommend `SSE_hUSI` when the data hold two states, senescent and non-senescent. | `sc/hUSI.R`; README, section "Note"; paper, Methods, "Classification of senescent cells based on hUSI" |
| Software | R 4.0.5, Seurat 4.2.1, mclust 6.0.1, dplyr 1.1.4, ggplot2 3.5.0, pROC 1.18.0 | README, section "Usage" |

## What the wrapper does

1. It reads the count table and stops if a value is not a whole number of zero or more, if a `GeneID` or a sample name is repeated, or if the table has fewer than two samples.
2. It gives each gene its `Symbol` and `Length` from the annotation table. Genes of the count table without a row in the annotation table are left out. Where several genes have the same symbol, the first is kept. Both rules are those of `ID2symbol` in `R/functions.R` of the repository. The numbers of genes left out are written to `run_info.json`.
3. It converts counts to TPM: count divided by `Length`, then scaled so that the values of a sample sum to one million. This is the formula of `Counts2TPM` in `R/functions.R`. It stops if a sample has no reads.
4. It sources `sc/hUSI.R` from the clone and calls `cal_hUSI` on the TPM matrix. The result is the score `hUSI`.
5. It calls `scoreOCLR(profile, SenOCLR, 'spearman')`, the function of `R/functions_new.R` that the authors' analysis scripts use, which returns the correlation before scaling. The result is the score `hUSI_unscaled`. The wrapper evaluates only the definition of `scoreOCLR`, because the rest of that file loads packages that scoring does not need. It stops unless scaling `hUSI_unscaled` gives `hUSI`.
6. It calls `SSE_hUSI` on `hUSI` and writes the class as `call`.
7. It writes `scores.tsv` and `run_info.json` after all steps succeeded.

The wrapper receives no labels. A failure in any step ends the run with a non-zero exit code and without an output folder.

## Scores

| `score_name` | Meaning | Direction | `call` |
| --- | --- | --- | --- |
| `hUSI` | The output of the released function `cal_hUSI`: the Spearman correlation between the gene weights and the TPM values of the sample, scaled over the samples of the count table to the range 0 to 1. This is the score the authors present as the main result. | Higher means more senescent. The authors compute the AUROC with senescent samples expected above non-senescent samples (`R/Validation.R`, `pROC::roc(..., levels=c('non-senescent','senescent'), direction='<')`). | The class of `SSE_hUSI`: `senescent` for the samples above the threshold (value 1 of the function), `not senescent` for the others (value 0) |
| `hUSI_unscaled` | The same Spearman correlation before scaling. Its possible range is -1 to 1. | Higher means more senescent. | Empty |

`hUSI` and `hUSI_unscaled` rank the samples of one count table in the same order, so an AUROC within one study is the same for both. `hUSI` of a sample changes when other samples are added to or removed from the count table, because the lowest and the highest sample define 0 and 1. `hUSI_unscaled` of a sample depends on that sample only. Only `hUSI_unscaled` can be compared between studies that were scored in separate runs.

## Deviations from the authors' use

1. **Source of the TPM values.** The authors' bulk RNA-seq inputs are TPM values from the pipelines of the data sources (GTEx, TCGA) or from StringTie (their own data and the training data). The wrapper computes TPM from the GEO counts and the column `Length` of the GEO annotation table, as the interface requires. Checked on GSE130727: this conversion reproduces GEO's own TPM table (`GSE130727_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz`) to the four digits that GEO prints.
2. **No further transformation of the TPM values.** In one analysis, the test of robustness to batch effects in `R/Model.R`, the authors apply `normalize_exp` before scoring: log2, scaling of each sample, and centring of each gene over the samples of the pooled data. The same transformation is applied to the training data before the weights are learned. In all applications to bulk RNA-seq in `R/Validation.R` and in the README the TPM values are scored as they are. The wrapper follows these applications. Centring each gene over the samples of a study would change the ranks of the genes within a sample, and therefore the score.
3. **Gene symbols.** The authors prescribe no translation of older or newer gene symbols, and the wrapper makes none. With the GEO annotation table, 13,337 of the 13,346 genes of the weights are found and 9 are missing (SMIM39, OR2I1P, FKBP1C, RAMACL, TCAF2C, NPIPA9, PNMA8C, TAF11L11, GAGE2E).
4. **Software.** R 4.0.5, dplyr 1.1.4 and mclust 6.0.1 are the versions of the README. Seurat, ggplot2 and pROC are not installed: the README uses them to read its example data, to plot and to compute an AUROC, and the scoring functions do not call them. jsonlite 1.8.8 is added to write `run_info.json`. The R packages are built from source from the CRAN snapshot of 2023-12-15 (`https://packagemanager.posit.co/cran/2023-12-15`), because conda-forge has no build of dplyr 1.1.4 or mclust 6.0.1 for R 4.0.
5. **`hUSI_unscaled` is an addition.** The released user function `cal_hUSI` returns only the scaled score. The unscaled correlation comes from the authors' function `scoreOCLR`, which their scripts call before `minmax`.
6. **`GMM_hUSI` is not run.** It returns numbered states without a rule for which states are senescent, and `scores.tsv` has one `call` per row.

## Data the released weights were fitted on

- 843 bulk RNA-seq samples of cultured human cells from 73 studies were collected (`Data/trainSet_meta843.csv`). After the filters of the paper, 770 samples remained: 385 senescent and 385 non-senescent, from 54 cell lines and 34 cell types (paper, Methods, "Data collection and pre-processing").
- The weights were learned with one-class logistic regression (`gelnet`, L2 penalty 1) from the 385 senescent samples only (paper, Methods, "Model training and hyperparameter optimization"; `R/Model.R`).
- Expression values: StringTie TPM of protein-coding genes (GENCODE version 31), without mitochondrial and ribosomal protein genes, without genes with low expression, and without the 85 most specific genes of each cell type. Before the fit the values were log-transformed, scaled per sample and centred per gene (`normalize_exp` in `R/functions_new.R`).
- The test study of this wrapper, GSE130727, is in `Data/trainSet_meta843.csv` (37 rows). Scores of hUSI on GSE130727 are therefore scores on data that the tool was collected from. The studies of this table are listed in `kb/identification/tool_training_datasets.tsv`.

## Limits

- `hUSI` and the class of `SSE_hUSI` are relative to the samples of one count table. `SSE_hUSI` calls at least one sample `senescent` and at least one `not senescent` in every count table with two or more different scores, whatever the samples are.
- A count table needs at least two samples. With one sample the scaling divides zero by zero.
- The weights are gene symbols of human genes. The wrapper does not accept other species.
- The released code was not changed, and no bug in it was found that affects scoring. One property is kept as released: `cal_hUSI` scales the correlations inside the function, so the released user function cannot return a score that is comparable between runs.
- Not checked: the example of the README (`Data/Aarts2017.rdata`) was not reproduced, because it needs Seurat, which is not installed.

## Test run

GSE130727 (34 samples), on biogemt-1, folder `/SCRATCH/dtzim01/senescence/runs/20261115-0012-tool-test-husi/` (`log.txt`, `GSE130727/scores.tsv`, `GSE130727/run_info.json`). The run took 2 seconds. `scores.tsv` has 68 rows: one row per sample for each of the two scores. All scores are finite, and the 34 values of each score are all different. A recomputation that uses neither the functions of the tool nor the wrapper (weights from `sc/SenOCLR_features.csv`, ranks and Pearson correlation in base R) agrees with `hUSI_unscaled` to 3e-8 and with `hUSI` to 2e-7. The scores were not compared with sample titles or labels.

## Environment

`/SCRATCH/dtzim01/senescence/tools/envs/husi` on biogemt-1 (1.4 GiB), created with:

```
source /SCRATCH/dtzim01/senescence/tools/envs/env.sh
bash /SCRATCH/dtzim01/senescence/code/tools/husi/setup_env.sh /SCRATCH/dtzim01/senescence/tools/envs/husi
```
