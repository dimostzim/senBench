# SenFinder wrapper (Mahmud 2024 classifier)

SenFinder is the released model of Mahmud et al., GeroScience 2024 (DOI 10.1007/s11357-024-01485-6, PMID 39738795, PMC12181610): a support vector machine that classifies a bulk RNA-seq sample of cultured human cells as senescent or not senescent from the expression of 20 genes. The repository calls the model `senFinder`. In `kb/identification/tools.tsv` the tool is the row `Mahmud 2024 classifier (called SenFinder by Bridge 2026)`. This wrapper runs the released script and the released model, unchanged, on the raw counts of one study.

- Repository: https://github.com/XiaoDongLab/SenFinder, commit `5f56506300baab2ed4c20af5a52c6252c4578b50` (2024-12-06). The wrapper reads the unchanged clone `/SCRATCH/dtzim01/senescence/tools/SenFinder` on biogemt-1 and stops if the clone is at another commit or holds a changed file.
- The repository releases a fitted model (`senFinder.rds`) and the script that applies it (`senFinder.R`). No model is trained by this wrapper.

## Result: one of the 20 genes of the model cannot be taken from the count table

The model needs a value for each of its 20 genes. One of them, `Lnc_PLOD2_2` (Ensembl gene ENSG00000261051, a long non-coding RNA that overlaps PLOD2 and has no gene symbol), is not a gene of the NCBI annotation. A count table with NCBI gene ids as `GeneID` therefore never holds it, and the wrapper cannot compute its TPM from the two input files.

- **As the interface requires, the wrapper stops in this case** (exit code 1, no output folder) and names the genes without a value. With the fixed input this happens for every study, so with the fixed interface SenFinder is not runnable.
- **The option `--set-absent-genes-to-zero` runs the model with 0 for these genes.** The README of the repository tells the user that the input file "should not have NA, Inf, use 0 instead through basic cleaning process". The option follows this sentence. It is a deviation from the interface, which allows no made-up value for a missing gene unless the tool's own code makes it. Whether the benchmark uses the option waits for the user (see Limits for what the option changes).

How I checked that the gene is absent: Table S3 of the paper's supplement gives the Ensembl identifier; the annotation table has no row with it; the Ensembl record (GRCh38, chromosome 3, 146,059,585 to 146,061,679, minus strand) has no cross-reference to an NCBI gene; and no gene of the annotation table overlaps this interval on NC_000003.12.

## How to run

```
bash run_tool.sh --counts <counts.tsv.gz> --annotation <annotation.tsv.gz> --out-dir <new folder> [--set-absent-genes-to-zero]
```

`--out-dir` must not exist. Two optional arguments name the environment and the clone when they are not at their places on biogemt-1: `--env-prefix` (default `/SCRATCH/dtzim01/senescence/tools/envs/senfinder`) and `--tool-repository` (default `/SCRATCH/dtzim01/senescence/tools/SenFinder`). The folder that holds the environment must also hold `env.sh`.

Files of this folder:

| File | Content |
| --- | --- |
| `run_tool.sh` | Entry point. Checks the arguments and the commit of the clone, then calls `score_senfinder.R`. |
| `score_senfinder.R` | Reads the two input files, converts the counts, writes the input file of the model, runs `senFinder.R` of the clone, writes the output. |
| `model_genes.tsv` | The 20 columns of the model with the Ensembl identifier and the gene symbol of Table S3. |
| `setup_env.sh`, `environment.yml` | Create the environment. |
| `conda_explicit_lock.txt` | The installed conda packages (R and all R packages) with their builds. |

## What SenFinder needs

| Need | Value | Source |
| --- | --- | --- |
| Input | A CSV file with one row per sample: the column `Sample_id` and one column per gene with the TPM (transcripts per million) of the gene. No NA and no infinite value. | README of the repository, "Step 1. Prepare input file"; `test_file.csv` |
| Genes | 20 genes, chosen by the authors with the mRMR algorithm (minimum redundancy, maximum relevance). The columns of the model are `RRM2B`, `CPXM1`, `PHGDH`, `PIK3IP1`, `YIU2`, `WASF1`, `MRPS9`, `TCF12`, `DNAI3`, `CAMK2N2`, `NSD1`, `MRPL11`, `F8`, `RPL10A`, `ISY1`, `LUCAT1`, `PURPL`, `CAPN10_DT`, `Lnc_PLOD2_2`, `LINC00663`. | `senFinder.rds` (element `coefnames`); README |
| Gene identifiers | Table S3 of the supplement gives an Ensembl gene identifier and a symbol for each of the 20 genes. The authors selected the genes by Ensembl identifier. | Supplement file `11357_2024_1485_MOESM2_ESM.xlsx`, Table S3; `pap_data_creation_unix_10292024.R` at commit `9e6ab08` |
| Normalisation | TPM as they are: no log transformation. The authors' TPM values come from RSEM after Trim Galore, alignment with STAR to GRCh37 with the Ensembl annotation of release 87 plus the ERCC spike-in sequences, and removal of duplicate reads with Picard. The training data stored in the model are TPM values with two decimals. | `rm_pipeline_gse56293.sh` and the three other pipeline scripts at commit `9e6ab08` (removed from the repository in commit `db10414`); element `trainingData` of `senFinder.rds` |
| Model | An object of the R package caret (`train`, method `svmLinear`): a support vector machine of the package kernlab with a linear kernel and cost C = 1. kernlab centres and scales each gene with the mean and standard deviation of the training data, which are stored in the model. The model also holds a function that turns the distance of a sample from the boundary into a class probability. | `senFinder.rds`, read with caret 6.0.94 and kernlab 0.9.33 |
| Parameters | None when the model is applied. | `senFinder.R` |
| How to apply | `Rscript senFinder.R "./senFinder.rds" "./test_file.csv"`. The script writes `predicted_results.csv` with the columns `Sample_id` and `Senescence_status` to the working directory. | README, "Step 2" |
| Software | R 4.4.1, dplyr 1.1.4, tidyr 1.3.1, ggplot2 3.5.1, lattice 0.22.6, caret 6.0.94, e1071 1.7.14 | README, "Dependencies" |

I could not read the Methods of the paper. PubMed Central and Europe PMC answer automated requests with a browser check, the publisher's page shows the abstract only, and I did not work around either. What this file says about the paper comes from the abstract, the supplement (open at the publisher) and the repository with its history.

## What the wrapper does

1. It reads the count table and stops if a value is not a whole number of zero or more, if a `GeneID` or a sample name is repeated, or if the table has no sample.
2. It converts counts to TPM: count divided by `Length` of the annotation table, then scaled so that the values of a sample sum to one million. All genes of the count table with a row in the annotation table enter the sum. It stops if a sample has no reads.
3. It finds each of the 20 genes of the model in the annotation table: the row whose `EnsemblGeneID` is the Ensembl identifier of Table S3, or, when no row has that identifier, the row whose `Symbol` is the symbol of Table S3. With the GEO annotation table, 18 genes are found by Ensembl identifier, `LINC00663` by symbol (its row has no Ensembl identifier), and `Lnc_PLOD2_2` is not found.
4. If a gene of the model has no row in the annotation table or no row in the count table, it stops and names these genes. With `--set-absent-genes-to-zero` it gives them the value 0 for every sample and goes on.
5. It writes the input file of the model (`senfinder_input.csv`: `Sample_id` and the 20 columns of the model) and runs `senFinder.R` of the clone with `senFinder.rds` of the clone on it, in a temporary working directory. The result, `predicted_results.csv`, is `Senescence_status`.
6. It asks the same model object for its class probabilities (`predict(svm_model, newdata = features, type = "prob")`, with `features` read from the same file as `senFinder.R` reads it). The probability of the class `X1` is `Senescence_probability`. It stops unless `Senescence_status` is `Yes` exactly for the samples whose probability of `X1` is the higher of the two.
7. It writes `scores.tsv` and `run_info.json` after all steps succeeded, and copies `senfinder_input.csv` and `predicted_results.csv` next to them.

The wrapper receives no labels. A failure in any step ends the run with a non-zero exit code and without an output folder.

## Scores

| `score_name` | Meaning | Direction | `call` |
| --- | --- | --- | --- |
| `Senescence_status` | The output of the released script: 1 for `Yes`, 0 for `No`. This is the result the authors present; they report accuracy of the class. | 1 means senescent ("'Yes' means predicted as senescence", README). | `senescent` for `Yes`, `not senescent` for `No` |
| `Senescence_probability` | The probability of the senescent class (`X1`) that the released model returns. Range 0 to 1. | Higher means more senescent. | Empty |

`Senescence_status` is `Yes` exactly when `Senescence_probability` is above 0.5, because caret calls the class with the higher probability (function `predict` of the method `svmLinear` in caret 6.0.94). Both scores of a sample depend on that sample only, so they can be compared between studies. For a measure that needs a ranking of the samples, such as an AUROC, `Senescence_probability` is the score to use (my suggestion; the authors report only the class).

## Deviations from the authors' use

1. **Source of the TPM values.** The authors' TPM values come from RSEM with the Ensembl annotation of release 87 on GRCh37, after removal of duplicate reads. The wrapper computes TPM from the GEO counts (reads counted on GRCh38 with the NCBI annotation) and the column `Length`, as the interface requires. I did not check whether GEO's counts include duplicate reads. The model applies fixed means and standard deviations of the authors' TPM values to these values and has no step that adjusts for a different pipeline.
2. **No rounding.** The TPM values of the training data and of `test_file.csv` have two decimals, as RSEM prints them. The authors do not prescribe rounding, and the wrapper does not round.
3. **Gene identifiers.** The model's columns are names. The wrapper finds the genes by the Ensembl identifiers of Table S3 and, for `LINC00663`, by symbol. That the column `YIU2` of the model is the gene YJU2 (ENSG00000105248) is my inference: Table S3 and Figure S2 of the supplement list YJU2, the other 19 columns match the other 19 genes of Table S3, and neither the HGNC gene name database nor the GEO annotation table has a gene with the symbol or synonym YIU2.
4. **`--set-absent-genes-to-zero`** (only when the option is given): 0 for genes of the model that the count table lacks, after the sentence of the README on NA values. Without the option nothing is filled in.
5. **`Senescence_probability` is an addition.** The released script writes only the class. The probability comes from the same model object through caret's own `predict`.
6. **Software.** R and the six packages are the versions of the README. kernlab 0.9.33 is added: the model is an object of a kernlab class and cannot be applied without kernlab, the README does not list it, and 0.9.33 was the current version when the model was released (my choice). jsonlite 1.8.9 is added to write `run_info.json`. The packages that the solver added are in `conda_explicit_lock.txt`; conda-forge built caret, dplyr, tidyr and kernlab under R 4.4.3, and R 4.4.1 prints a warning about this when it loads them.

## Data the released model was fitted on

- 114 labelled bulk RNA-seq samples of cultured human cells: 71 not senescent and 43 senescent (element `trainingData` of `senFinder.rds`; README: "we used all the 114 samples"; Table S2 of the supplement has the same numbers).
- The samples come from four studies (Table S1): GSE63577 and GSE64553 (replicative senescence), GSE56293 (bleomycin) and E-MTAB-5403 (irradiation). Table S2 lists fibroblasts (BJ, MRC-5, HFF, IMR-90, WI-38, foreskin fibroblasts), keratinocytes and melanocytes. The six quiescent fibroblast samples of E-MTAB-5403 are labelled not senescent. A further 48 samples of these studies (intermediate passages, and 4 and 20 days after irradiation) are labelled unknown and were not used for the fit.
- The 20 genes were selected on the same 114 samples.
- The test study of this wrapper, GSE130727, is not among the four studies. The studies are listed in `kb/identification/tool_training_datasets.tsv`.

## Limits

- **One gene has no value with the fixed input** (see Result). With `--set-absent-genes-to-zero` the model receives 0 for `Lnc_PLOD2_2` in every sample. The model is linear, so this moves every sample by the same amount towards `not senescent`: 0.43 units of the distance from the boundary, compared with a sample at the training mean of this gene (0.38 TPM); the distances of the training samples range from -3.6 to 6.2. The order of the samples by `Senescence_probability` is then the order given by the other 19 genes, and what this gene would have added is lost. In the training data the gene is 0 in 4 of 114 samples, and the calls of the 114 training samples do not change when it is set to 0.
- **Other genes can be absent** in count tables that come from authors' files. Of the 291 count tables in `tables/benchmark_bulk/read975/counts` on biogemt-1 (read on 2026-11-15), 212 hold the 19 other genes, 57 hold 18 of them, and 22 hold fewer. With the option, each absent gene is 0, and `run_info.json` names them (`tool_genes_set_to_zero`).
- **The scale of the input matters.** The model compares TPM values with fixed means of the authors' TPM values. In the test input, 14 of the 20 genes have at least one sample outside the range of the training data.
- The released code was not changed. Three properties are kept as released: the README does not list kernlab, which the model needs; a column of the model is named `YIU2`; and `senFinder.R` loads ggplot2 and lattice, which it does not use.
- Not checked: the repository gives no expected output for `test_file.csv`, so the test shows that the released script runs on it, and nothing more. The accuracy that the paper reports was not reproduced.

## Test run

On biogemt-1, folder `/SCRATCH/dtzim01/senescence/runs/20261115-0213-tool-test-senfinder/` (`command.sh`, `log.txt`).

1. `senFinder.R` of the clone on `test_file.csv` of the clone: exit code 0, five samples classified (`repository_test_file/predicted_results.csv`), 2 seconds.
2. The wrapper on GSE130727 (34 samples): it stops with exit code 1 and the message that the count table has no gene for `Lnc_PLOD2_2` (ENSG00000261051); no output folder. 4 seconds.
3. The wrapper on GSE130727 with `--set-absent-genes-to-zero`: exit code 0, 5 seconds, folder `GSE130727_absent_genes_set_to_zero/`. `scores.tsv` has 68 rows: one row per sample for each of the two scores. All scores are finite; the 34 values of `Senescence_probability` are all different, and `Senescence_status` takes both values. 19 genes of the model were found and 1 was set to 0.
4. The wrapper refuses an output folder that exists (exit code 1).

Two recomputations that use neither the wrapper's functions nor `predict` agree with the output of step 3: TPM of the 19 genes from the two input files (largest difference 5e-12), and the probability from the support vectors, the scaling values and the probability function stored in the model (largest difference 1e-15). The scores were not compared with sample titles or labels.

## Environment

`/SCRATCH/dtzim01/senescence/tools/envs/senfinder` on biogemt-1 (1.5 GiB), created with:

```
source /SCRATCH/dtzim01/senescence/tools/envs/env.sh
bash /SCRATCH/dtzim01/senescence/code/tools/senfinder/setup_env.sh /SCRATCH/dtzim01/senescence/tools/envs/senfinder
```
