# SENCAN wrapper

SENCAN is the classifier of Jochems et al., Cell Rep 2021 (DOI 10.1016/j.celrep.2021.109441, PMC8333195): an elastic net that calls a bulk RNA-seq sample of cancer cells senescent or not senescent. This wrapper applies the released classifier object, unchanged, to the raw counts of one study of bulk RNA-seq samples.

- The paper states that it "does not report original code". SENCAN has no repository of its own.
- Classifier object: the file `Data/SENCAN_classifier.rda` of https://github.com/WJPina/HUSI (the hUSI repository), commit `300da32f111f8c8603e822eca8e84c9b5a76db48` (2025-09-11). The wrapper reads the unchanged clone `/SCRATCH/dtzim01/senescence/tools/HUSI` on biogemt-1 and stops if the clone is at another commit or holds a changed file.
- The hUSI authors apply the classifier object in the function `cal_SENCAN` of `R/functions_new.R`. The wrapper makes the same steps as this function, with one difference in how genes are matched (deviation 2).
- The paper says that a script for calculating scores with the classifier is at the SENESCopedia website (Method details, "Senescence classifier"). The website was not used and the script was not read, by instruction. Whether `cal_SENCAN` is a copy of that script is therefore not known.

Three sources are named below: **paper** (the SENCAN paper), **object** (the content of the classifier object), and **hUSI use only** (the step is in `cal_SENCAN` and neither in the paper nor fixed by the object).

## How to run

```
bash run_tool.sh --counts <counts.tsv.gz> --annotation <annotation.tsv.gz> --out-dir <new folder>
```

`--out-dir` must not exist. Two optional arguments name the environment and the clone when they are not at their places on biogemt-1: `--env-prefix` (default `/SCRATCH/dtzim01/senescence/tools/envs/sencan`) and `--tool-repository` (default `/SCRATCH/dtzim01/senescence/tools/HUSI`). The folder that holds the environment must also hold `env.sh`.

Files of this folder:

| File | Content |
| --- | --- |
| `run_tool.sh` | Entry point. Checks the arguments and the commit of the clone, then calls `score_sencan.R`. |
| `score_sencan.R` | Reads the two input files, prepares the counts, applies the classifier object, writes the output. |
| `setup_env.sh`, `environment.yml` | Create the environment. |
| `conda_explicit_lock.txt` | The installed conda packages with their versions. All R packages are conda packages. |

## What the classifier object holds

| Object | Content |
| --- | --- |
| `dge` | Estimated counts of 20,029 genes in 60 reference samples: the authors' cancer cell lines, untreated or treated with alisertib or etoposide. Per gene: `ensembl_id`, `gene_name`, `hgnc_symbol`. |
| `dmCPM` | log2 counts per million of the 60 reference samples. Its lowest value is -2.99466. |
| `edger_ref_column` | The counts of one reference sample, `HCT116_P_batch3`, against which samples are normalised. |
| `calcNormFactors_WSPref`, `.calcFactorTMMwsp` | The authors' normalisation function and the edgeR function it calls. |
| `gene_mean`, `gene_sd` | Mean and standard deviation of every gene, used to centre and scale. |
| `use_gene_ix` | The 17,405 genes given to the elastic net. |
| `full_classifier` | The fitted elastic net (`cv.glmnet`, `alpha = 0.7`, classes `FALSE` and `TRUE`). At `lambda.1se`, 137 genes have a coefficient other than zero. |

Checks of the object against the paper, all made in the environment of this wrapper:

- The 137 genes, their coefficients and the intercept are identical to Table S4 of the paper (largest difference 0).
- The mixing parameter 0.7, the 17,405 input genes and the lowest value -2.99 are those of the paper (Method details, "Senescence classifier").
- `edgeR::cpm(dge, log = TRUE)` reproduces `dmCPM` (largest difference 4e-15), and `calcNormFactors_WSPref` on the 60 reference samples reproduces the normalisation factors stored in `dge` (difference 0). The log2 counts per million of the classifier are therefore those of `edgeR::cpm` with its default prior count.

## What SENCAN needs

| Need | Value | Source |
| --- | --- | --- |
| Input | Gene-level estimated counts of the sample, re-quantified from the raw reads with kallisto 0.46 against GENCODE version 34 (without unexpressed pseudogenes and without transcripts of support level 5 that have no consensus coding sequence), summed per gene | Paper, Method details, "RNA-sequencing" and "Senescence classifier". The wrapper cannot follow this step (deviation 1). |
| Gene identifiers | The genes of GENCODE version 34. The object names them by Ensembl gene identifier (row names of `dmCPM`, names of `gene_mean`, `gene_sd` and `use_gene_ix`) and also holds `gene_name` and `hgnc_symbol`. The paper gives no rule for an input with other gene identifiers. `cal_SENCAN` matches the row names of its input to `gene_name`. | Paper; object; matching by gene name: hUSI use only |
| Genes that the input lacks | No rule in the paper: an input quantified as the paper prescribes has every gene. `cal_SENCAN` fills them per sample with the mean count of the two reference samples whose counts correlate best with the sample (Kendall correlation over the genes present, `pcaPP::cor.fk`). | hUSI use only |
| Normalisation between samples | edgeR, trimmed mean of M-values with singleton pairing, against the reference sample of the authors' data | Paper, "Senescence classifier"; object (`calcNormFactors_WSPref` with `edger_ref_column`) |
| Samples normalised together | `calcNormFactors_WSPref` divides the factors of the samples it is given by their geometric mean. The paper does not say which samples are normalised together. `cal_SENCAN` passes all samples of a data set in one call. | Object (the division); hUSI use only (all samples of a data set together) |
| Expression values | log2 counts per million; values below the lowest value of the authors' data (-2.99) are set to that value | Paper, "Senescence classifier"; object (`dmCPM`) |
| Centring and scaling | Per gene, with the mean and standard deviation of the authors' data | Paper, "Senescence classifier"; object (`gene_mean`, `gene_sd`) |
| Parameters | The penalty chosen by the one-standard-error rule (`s = "lambda.1se"`). Nothing is set at scoring. | Paper, "Senescence classifier"; object |
| Score | The probability of the class senescent from the logistic elastic net | Paper, Figure 6 ("SENCAN scores", "senescence score"); object |
| Class | Senescent above a score of 0.5. The text of the paper states no threshold. Figure 6B to 6D draw a dashed line at 0.5, and with 0.5 the scores of Table S5 give the counts of the text (44 of 49 correct, precision 0.94, sensitivity 0.8, specificity 0.97) when the five irradiated U2OS samples are counted as senescent, as Figure 6C draws them. Table S5 marks these five samples as not senescent. | Paper, Figure 6 and Table S5 |
| Software | R 3.6.0, edgeR 3.26.8, glmnet (version not given) | Paper, "RNA-sequencing", "Senescence classifier", "Quantification and statistical analysis" |

## What the wrapper does

1. It reads the count table and stops if a value is not a whole number of zero or more, or if a `GeneID` or a sample name is repeated.
2. It gives each gene its `EnsemblGeneID` from the annotation table. Genes without a row in the annotation table and genes without an Ensembl gene identifier are left out. The counts of genes that share an Ensembl gene identifier are summed. The numbers are written to `run_info.json`.
3. It places the counts on the 20,029 genes of the classifier object by Ensembl gene identifier.
4. It fills the genes that the input lacks by the rule of `cal_SENCAN` (table above).
5. It calls `calcNormFactors_WSPref` of the object on all samples of the count table together, computes log2 counts per million with `edgeR::cpm(log = TRUE)`, sets values below the lowest value of `dmCPM` to that value, and centres and scales with `gene_mean` and `gene_sd`.
6. It calls `predict` on `full_classifier` with `s = "lambda.1se"` for the 17,405 genes of `use_gene_ix`: once for the probability (the score) and once for the class (the call).
7. It writes `scores.tsv` and `run_info.json` after all steps succeeded.

The wrapper receives no labels. A failure in any step ends the run with a non-zero exit code and without an output folder.

## Score

| `score_name` | Meaning | Direction | `call` |
| --- | --- | --- | --- |
| `SENCAN_score` | The probability of the class senescent from the elastic net, between 0 and 1. It is the only score, and the one the authors present. | Higher means more senescent (paper, Figure 6: senescent samples near 1, untreated samples near 0). | The class that `predict(type = "class")` returns: `senescent` for a score above 0.5, `not senescent` otherwise |

## Deviations from the authors' use

1. **Quantification.** The paper re-quantifies every sample from the raw reads with kallisto against GENCODE version 34. The wrapper receives GEO's uniformly processed raw counts per NCBI gene, as the interface requires. The counts differ in the read assignment (GEO counts aligned reads per gene, kallisto estimates transcript abundances) and in the set of genes.
2. **Gene identifiers.** The wrapper matches by Ensembl gene identifier, through the column `EnsemblGeneID` of the GEO annotation table. `cal_SENCAN` matches gene symbols to `gene_name`. Reason: the object names its genes by Ensembl gene identifier, and gene names of GENCODE version 34 have since changed (for example DDX58 is now RIGI). Genes found in the GEO table, by Ensembl gene identifier and by gene name: 119 and 114 of the 137 genes with a coefficient; 15,774 and 15,729 of the 17,405 genes given to the elastic net; 17,372 and 17,346 of all 20,029 genes. For 41 of the 16,735 genes found both ways, the two ways point to different NCBI genes.
3. **Genes that share an Ensembl gene identifier.** In the GEO annotation table, 103 genes carry an Ensembl gene identifier that another gene also carries. The wrapper sums their counts. The paper has no such case. 74 of the 20,029 genes of the object are affected, none of the 137 genes with a coefficient.
4. **Genes that the input lacks are filled in.** This follows `cal_SENCAN` and is not in the paper. With the GEO table, 2,657 of the 20,029 genes are missing, 1,631 of the 17,405 genes given to the elastic net, and 18 of the 137 genes with a coefficient. These 18 genes carry 21.6% of the sum of the absolute coefficients. Among them are the gene with the largest coefficient (ENSG00000240219, AL512306.2, 0.879) and the gene with the most negative one (ENSG00000276517, AL133243.3, -0.801). 13 of the 18 have a clone-based name such as AL512306.2 and no gene symbol. Their Ensembl gene identifiers are listed in `run_info.json`.
5. **Samples normalised together.** The wrapper normalises all samples of one count table in one call, as `cal_SENCAN` does with a data set.
6. **Software.** R 4.0.5, edgeR 3.32.1, glmnet 4.1-2, pcaPP 2.0-2. R 4.0.5 is the version that the hUSI README states; the classifier object was written by R 4.0.3. The versions of the paper (R 3.6.0, edgeR 3.26.8) were not installed. The wrapper uses the `.calcFactorTMMwsp` stored in the object. It differs from the function of edgeR 3.32.1 only in using `rep.int` where edgeR 3.32.1 uses `rep_len`. jsonlite 1.8.0 is added to write `run_info.json`.
7. **Class.** `cal_SENCAN` returns the score only. The wrapper adds the class from `predict(type = "class")`.

The steps of `cal_SENCAN` are written out in `score_sencan.R` and not called from the clone, because `cal_SENCAN` takes a Seurat object and matches by gene name. Checked on GSE130727: the unchanged `cal_SENCAN`, evaluated from the clone and given the same counts with Ensembl gene identifiers as gene names, returns the scores of the wrapper (largest difference 4e-16).

## Data the released classifier was fitted on

- The authors' RNA-seq of 13 cancer cell lines (HepG2, Huh7, Hep3B, MCF-7, T47D, MDA-MB-231, SUM159, HCT116, LoVo, RKO, A549, H358, PC9), untreated and made senescent with alisertib or etoposide; ArrayExpress E-MTAB-9970 (paper, "Data and code availability" and key resources table). The 60 reference samples of the object are of these cell lines.
- 45 untreated samples of the Cancer Cell Line Encyclopedia, three cell lines of each of 15 cancer types (paper, "Senescence classifier", Table S3). They are not in the object.
- The elastic net was fitted on all these samples after a cross-validation that leaves out one cell line (paper, "Senescence classifier").
- The test study of this wrapper, GSE130727, is not among the fitted data. The authors used 25 of its samples as an external validation set (paper, Figure 6D). The paper gives no table of the scores of these samples.

## Limits

- **The score of a sample depends on the other samples of the count table.** `calcNormFactors_WSPref` divides the normalisation factors by their geometric mean over the samples it is given. Checked on GSE130727: GSM3752524 scores 0.463 when all 34 samples are scored together and 0.271 when it is scored with two other samples. A sample scored alone gets the factor 1, so it is not adjusted to the reference sample at all. Scores are comparable within one count table. Between count tables they are comparable only as far as this dependence allows.
- **Filled-in genes carry values of the fitted data.** The filled-in counts are those of two reference samples, which are samples the classifier was fitted on. They are copied as counts, without scaling to the number of reads of the scored sample. The reference samples have 12 to 23 million estimated counts; the samples of GSE130727 have 29 to 153 million reads in the genes found. This is the behaviour of `cal_SENCAN` and was kept.
- The classifier was fitted on cancer cell lines. The authors state that a low score does not rule out senescence, that a sample is called senescent only if most of its cells are senescent, and that the accuracy in vivo is unclear (paper, Results, "SENCAN classifier for cancer cell senescence").
- Human genes only.
- Not checked: the authors' script at the SENESCopedia website; the scores of the paper for GSE130727 (the paper prints none); a run with the software versions of the paper. Read from the code and not run: `cal_SENCAN` fails on a count table with one sample or with exactly one missing gene; `score_sencan.R` handles both.

## Test run

GSE130727 (34 samples), on biogemt-1, folder `/SCRATCH/dtzim01/senescence/runs/20261115-0019-tool-test-sencan/` (`command.sh`, `log.txt`, `code/`, `GSE130727/scores.tsv`, `GSE130727/run_info.json`). The run took 15 seconds. `scores.tsv` has 34 rows, one per sample. All scores are finite, and the 34 values are all different (0.0009 to 1.0000). Both classes occur among the calls. A second run gave an identical `scores.tsv`. The scores were not compared with sample titles or labels.

Also checked: the wrapper refuses an existing `--out-dir`; it stops on a count that is not a whole number and on an annotation table without the column `EnsemblGeneID`; it scores a count table with one sample.

## Environment

`/SCRATCH/dtzim01/senescence/tools/envs/sencan` on biogemt-1 (1.4 GiB), created with:

```
source /SCRATCH/dtzim01/senescence/tools/envs/env.sh
bash /SCRATCH/dtzim01/senescence/code/tools/sencan/setup_env.sh /SCRATCH/dtzim01/senescence/tools/envs/sencan
```
