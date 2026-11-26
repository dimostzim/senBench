# SenCat RNA score: wrapper for bulk RNA-seq raw counts

The SenCat RNA score is the score from RNA features of the row `SenCat ML score` in `kb/identification/tools.tsv` (Anerillas 2026, Mol Cell, DOI 10.1016/j.molcel.2026.05.017). This folder holds a wrapper that applies the released weights to the raw counts of one study with the tool's own scripts. The score from protein features is not covered.

## Is a fitted model released?

Yes. The repository `wf-ml-markers-senescence` releases fixed weights, and its README describes how to apply them to new samples. No model has to be trained, and the wrapper trains none.

Released:

- The RNA marker file `analysis/transcriptomics_loose5000f_tuned_common_features.csv`: 5,000 Ensembl gene ids, each with a gene symbol and one coefficient (2,543 positive, 2,457 negative, none zero; md5 `672d8d8850f489da161af3b1d1627434`).
- The normalisation script `workflow/scripts/data/normalize_counts.py` and the scoring script `workflow/scripts/cls/marker_classifier.py`.
- The training script `workflow/scripts/cls/ml_classifier.py` and the configuration `config/config.yml`.

Not released:

- A threshold or a class call. The scoring script writes one number per sample.
- An intercept. The 14 fitted logistic regression models are not in the repository; only the mean of their coefficients is.
- The training table (`Supplementary table 1 transcriptomic results.xlsx`, Table S1 of the paper) and the external IMR-90 dataset of Figure 3E (`23_all_sample_gene_counts.txt`, `23_coldata.xlsx`). The paper gives no accession for the IMR-90 dataset, so the authors' external result cannot be repeated here.
- Package versions, apart from Python 3.11.10 (`workflow/envs/data.yaml`, `workflow/envs/pydeseq2.yaml`).

The three released copies differ:

| Copy | RNA marker file | Scoring script |
| --- | --- | --- |
| Zenodo DOI 10.5281/zenodo.19333691 (version 1.0.0, tag `v1.0.0`) | absent; the archive holds only tables of scores per sample | needs the label column `is_sen` |
| Zenodo DOI 10.5281/zenodo.20008425 (version 1.0.1, tag `v1.0.1`, commit `6c16d47`), the archive cited in the paper | present, same md5 as in the clone | needs the label column `is_sen` in the input (`y = adata.obs['is_sen']`), so it stops on samples without labels |
| Clone `/SCRATCH/dtzim01/senescence/tools/wf-ml-markers-senescence`, commit `0e23fcc` (2026-08-19) | present | scores samples without labels (commit `5d78c0b`, 2026-07-28) |

The wrapper uses the clone at commit `0e23fcc`. The score formula is the same in versions 1.0.1 and `0e23fcc`. `normalize_counts.py` is identical in both.

The second repository, `SenCat-Companion` (commit `a30a162`), holds R scripts for elastic-net feature selection and cross-validation on tables of clinical traits. It holds no RNA weights and no fitted model, and the wrapper does not use it.

## What the tool needs

| Need | What the authors prescribe | Source |
| --- | --- | --- |
| Input | An `h5ad` file: `adata.X` is samples by genes with normalised counts, `adata.var_names` are the gene identifiers of the marker file | Repository `README.md`, "Prepare data" |
| Gene identifiers | Ensembl gene ids without version | Index column `converted_alias` of the RNA marker file; paper, Method details, "Machine learning (ML) determination of senescence markers" |
| Normalisation of raw counts | `normalize_counts.py`. It removes genes with a count in at most 5% of the samples and genes with a total count of at most 380. If fewer than 100 genes have a count in every sample, it adds 1 to every count. It then divides each sample by its DESeq2 size factor (pydeseq2, `fit_size_factors`) | Repository `README.md`, "If your data are not normalized"; rule `normalize_counts` in `workflow/rules/data.smk`; paper: "We applied the same normalization and scoring strategy to this dataset without retraining the model" |
| Design | `normalize_counts.py` needs a design formula. For the external IMR-90 dataset the configuration gives `sample`, a column that holds the sample name | `config/config.yml`, `NORMALIZATION_DESIGN`; `workflow/scripts/data/anndata_from_txt.py` line 24 |
| Transformation | `log1p`, applied inside `marker_classifier.py` | `marker_classifier.py` line 71; repository `README.md`, Notes |
| Score | The sum over the markers found in the input of coefficient times `log1p` of the normalised count. Markers that are not in the input are skipped | `marker_classifier.py` lines 25 to 53 |
| Parameters | None | |

## What the wrapper does

`bash run_tool.sh --counts <counts.tsv.gz> --annotation <annotation.tsv.gz> --out-dir <new folder>`

1. It reads the raw counts and stops if a value is not a whole number of zero or more, or if a `GeneID` or a sample occurs twice.
2. It converts NCBI gene ids to Ensembl gene ids with the column `EnsemblGeneID` of the annotation table (section Deviations, point 1) and writes the counts as `h5ad`, with the sample name in the column `sample`.
3. It runs the tool's `normalize_counts.py` with `--design sample`.
4. It runs the tool's `marker_classifier.py` with the RNA marker file.
5. It stops if the tool's result does not hold exactly one finite score for every input sample, in the input order. It then writes the output.

Output in `--out-dir`:

- `scores.tsv`: one row per sample, `tool` = `sencat`, `score_name` = `rna_score`, `call` empty.
- `run_info.json`: tool commit, md5 of the RNA marker file, the counts of the gene identifier conversion, the number of markers found and missing with the reason, and package versions.
- `sencat_results.csv`: the file that `marker_classifier.py` wrote, unchanged.
- `marker_classifier.log` and `normalize_counts.log`: the tool's own logs. `normalize_counts.log` is empty, because that script sends its messages to a logger without an output (released behaviour).

The wrapper receives no labels. It refuses an `--out-dir` that exists. The two `h5ad` files are written to a temporary folder inside `--out-dir` and removed at the end.

`run_tool.sh` expects the folder with the clones and `envs/` at `/SCRATCH/dtzim01/senescence/tools`; the variable `SENESCENCE_TOOLS_DIR` names another place. It sources `envs/env.sh`, uses the Python of `envs/sencat`, limits numeric libraries to 8 threads, and sets `PYTHONDONTWRITEBYTECODE=1` so that Python writes nothing into the clone.

## The score

`rna_score` is the only score and the one the authors present. The repository README states that "higher values indicate stronger similarity to the senescence-associated signature". The paper describes it as "a linear combination of log-transformed RNA or protein abundance weighted by their ML-derived coefficients". It is not a probability: it has no intercept and no fixed range, and it can be negative. The tool assigns no class.

## Deviations from the authors' use

1. **Gene identifier conversion.** The authors' input is a count table indexed by Ensembl gene id. The benchmark input is indexed by NCBI gene id, so the wrapper converts with the column `EnsemblGeneID` of GEO's annotation table. Genes without an Ensembl gene id are dropped before the normalisation (12,299 of 39,376 rows of the annotation). Where several NCBI genes have the same Ensembl gene id, their counts are summed (100 Ensembl gene ids, 203 NCBI genes). The authors prescribe nothing for this case; the sum is the wrapper's choice. 120 of the 5,000 markers have no row in the annotation table and can therefore never be found with this input; 17 markers are among the summed Ensembl gene ids.
2. **Design `sample`.** The repository README leaves the design to the user. The wrapper uses the design that the authors' configuration gives for their own external dataset. It names no condition, so no label enters. In the test run the scores equal a calculation by hand with median-of-ratios size factors, which do not use the design.
3. **One environment in place of two, with versions chosen at install.** The authors give two conda environments (`data.yaml` for scoring, `pydeseq2.yaml` for normalisation) and fix only Python 3.11.10. The wrapper uses one environment with Python 3.11.10, scanpy, anndata and pandas from conda-forge and pydeseq2 from pip, at the versions available on 2026-11-15 (`environment.yml`, `conda_explicit_lock.txt`, `pip_freeze.txt`). The channel `defaults` of the authors' files is not used. `openpyxl`, which the scoring does not need, is not installed.
4. **Clone in place of the archived release.** The paper cites Zenodo version 1.0.1, whose scoring script stops on samples without labels (table above). The wrapper uses commit `0e23fcc`, with the same marker file and the same formula.

The tool's code and the RNA marker file are used unchanged. Nothing is tuned or fitted.

## Data the released weights were fitted on

From the paper and the repository; not rerun here.

- SenCat bulk RNA-seq of 14 primary human cell types under more than 30 senescence models: GEO GSE284666 (190 samples), part of SuperSeries GSE285826. A study that holds these samples is training data of this tool.
- Label: taken from the sample name. Treatment `P` (proliferating) or `EV` (empty vector) is not senescent; every other treatment is senescent (`workflow/scripts/data/anndata_from_excel.py`).
- Values: the workflow trains on Table S1 of the paper as it is and does not normalise it. The paper states that counts were normalised with DESeq2 (v1.32.0). Table S1 was not opened here.
- The released weights come from the profile `loose5000f` of `config/config.yml`: 5,000 genes selected with `SelectKBest` on all samples, 14 logistic regression models that each leave one cell type out, and the mean coefficient of a second set of 14 models fitted on the 5,000 genes. The Methods of the paper describe 1,000 genes, an L1 penalty, a weight threshold of 0.05 and markers that pass in all 14 models. `kb/identification/benchmarks.md`, section "Differences between the SenCat Methods and the released code", lists these points.

## Limits

- **The score of a sample depends on the other samples in the counts file.** The gene filter and the size factors of `normalize_counts.py` are computed over all samples of the file. The same sample scored with other samples can get another score, and the number of markers that enter the score differs between studies. Scores of different studies are therefore not on one scale. `run_info.json` records the number of markers found for each run.
- **Missing markers are skipped, not replaced.** The sum then runs over fewer markers. In the test run 4,873 of 5,000 markers were found (96.3% of the summed absolute coefficients): 120 have no Ensembl gene id in the annotation table and 7 were removed by the tool's gene filter.
- The wrapper was run on one study of 34 samples. Files with very few samples were not tried.
- pydeseq2 0.5.4 prints two `ImplicitModificationWarning` lines when `normalize_counts.py` runs. The run ends with exit code 0.
- The repository has no licence file. The Zenodo record is CC BY 4.0.
- Not checked: whether the released weights are the weights behind the figures of the paper, and whether they equal Table S4 of the paper.

## Environment

On biogemt-1: `/SCRATCH/dtzim01/senescence/tools/envs/sencat` (1.1 GiB). To create it again:

```
source /SCRATCH/dtzim01/senescence/tools/envs/env.sh
export PYTHONNOUSERSITE=1
micromamba create -y -p /SCRATCH/dtzim01/senescence/tools/envs/sencat --file /SCRATCH/dtzim01/senescence/code/tools/sencat/conda_explicit_lock.txt
/SCRATCH/dtzim01/senescence/tools/envs/sencat/bin/python -m pip install --no-user --no-deps pydeseq2==0.5.4 formulaic-contrasts==1.0.0 session-info==1.0.1 stdlib-list==0.12.0
```

`pip_freeze.txt` is the output of `pip list --format=freeze` and lists every Python package of the environment with its version.

## Test run

GEO's uniformly processed raw counts of GSE130727 (34 samples), on biogemt-1, started through the wrapper:

- Run folder: `/SCRATCH/dtzim01/senescence/runs/20261115-0009-tool-test-sencat/` (output in `out/`), marker `DONE`, 5 seconds.
- `scores.tsv` has 34 rows, one per input sample in the input order; all scores are finite and all 34 are distinct.
- The scores equal a calculation by hand from the two input files and the RNA marker file (largest absolute difference 1.5e-14).
- After the run the clone has no changed or new file.
- The scores were not compared with sample titles or labels.

The folder `/SCRATCH/dtzim01/senescence/runs/20261115-0008-tool-test-sencat/` holds an earlier attempt with the marker `FAILED`. The wrapper ended with exit code 0 and the same scores; the marker comes from an expected-output path that was given relative to the wrong folder when the run was started.
