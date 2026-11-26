# SenCID wrapper

SenCID (Tao, Yu and Han, Cell Metab 2024, DOI 10.1016/j.cmet.2024.03.009; repository https://github.com/JackieHanLab/SenCID, commit `f93cba19b9d44e84f89504ecc6aaeb108ce8e508` of 2024-02-22, version 1.0.0) gives each sample a senescence score from each of six models, SID1 to SID6, and has a recommendation model that names one of the six for the sample.

This wrapper runs the released code and the released model files, unchanged, on the raw counts of one study of bulk RNA-seq. State on 2026-11-15: the environment builds from public package channels with the versions the repository pins, and the wrapper runs.

```
bash run_tool.sh --counts <counts.tsv.gz> --annotation <annotation.tsv.gz> --out-dir <new folder>
```

- `--counts`: tab-separated, gzip. First column `GeneID` (NCBI gene id), then one column of raw read counts per sample.
- `--annotation`: GEO's `Human.GRCh38.p13.annot.tsv.gz`. The columns `GeneID` and `Symbol` are used.
- `--out-dir`: a new folder. The wrapper stops if it exists.

The wrapper exits with a non-zero code when the counts are not whole numbers of zero or more, when SenCID writes no output, when the rows of the output are not the samples of the input, or when a score is not finite.

## Files

| File | Content |
| --- | --- |
| `run_tool.sh` | Entry point. Sets the environment and starts `run_sencid.py`. |
| `run_sencid.py` | Converts gene ids to symbols, runs SenCID's command line, writes `scores.tsv` and `run_info.json`. |
| `config.sh` | The paths on biogemt-1, the commit of the clone and the limit of 8 threads. |
| `install_env.sh` | Builds the environment and installs SenCID from the clone. |
| `environment.yml` | The versions. A copy of the repository's `SenCID.yaml` with every pin unchanged. |
| `environment.explicit.lock`, `environment.pip.lock` | What was installed (`micromamba env export --explicit`, `pip freeze --all`). |

## What SenCID needs

The Methods of the paper were not read: the publisher's article page refuses automated requests (error 403), and Europe PMC and Unpaywall list no other open copy. Every entry below comes from the repository or from the replies of the first author (GitHub user SherryMiyano) in the issues of the repository, read on 2026-11-15.

| Need | Value | Source |
| --- | --- | --- |
| Input | Raw counts, genes by samples. Bulk samples are accepted. | `SenCID/api.py`, docstring of `SenCID` ("SenCID inputs the raw counts of data"). `SenCID/demo/SenCID_tutorial.ipynb` ("for each cell (or sample for bulk RNA-seq)"). Issue 12: "The input for SenCID is only a count matrix." |
| Gene identifiers | HGNC gene symbols. | `SenCID/DataPro.py` line 100 ("Features of count data should be hgnc symbol: BIRC5, BRCA1, etc"). `SenCID/resource/seneset.txt` holds the 1,290 symbols that SenCID reads; this file is called the gene list below. |
| Normalisation | None before SenCID. SenCID normalises the counts itself (next section). | `SenCID/DataPro.py`, functions `ScaleData`, `Cpm_added` and `Zcol`. |
| Denoising | Off for bulk RNA-seq. Denoising fits an autoencoder (the package `dca`) to the input and is meant for single cells. | Issue 6: "always turn off DCA imputation (by setting parameter denoising = False) when you perform SenCID in bulk RNA-seq data". The same in issues 1, 9 and 13. `SenCID/api.py`: "Recommended for single cell data. True by default." |
| Models | All six (`sidnums` 1 to 6). | Default of `SenCID/api.py` and of `SenCID/__main__.py`. |
| Calls | `binarize` on. | Default of `SenCID/__main__.py`; the tutorial sets it on. |
| Species | Human. | Issue 8: "currently the models are all trained by human datasets". |
| Software | Python 3.6.15, TensorFlow 2.2.0, scikit-learn 0.24.2 and the other pins. | `SenCID.yaml`. The model files were saved with scikit-learn 0.24.2 (read from the files). |

### What SenCID does with the counts

Read from `SenCID/DataPro.py` and `SenCID/Pred.py`:

1. It keeps the genes of the gene list that are in the input and have a total count of at least 1 over the samples of the input file.
2. It adds every other gene of the gene list with the value 0.01 in every sample.
3. In each sample it divides by the sum over the 1,290 genes and multiplies by 1,000,000. The sum is over the gene list, not over all genes of the sample.
4. It takes log2(value + 1) and, within each sample, subtracts the mean and divides by the standard deviation over the 1,290 genes.
5. Each of the six models reads its own genes from this table (19, 13, 186, 87, 83 and 224 genes for SID1 to SID6). The recommendation model reads all 1,290 genes.

With denoising off, nothing is fitted to the input. Step 1 is the only step in which the samples of a file affect each other. On the test input, scoring the first 5 samples alone gave the same scores as scoring all 34 (largest difference 0).

## Output

`scores.tsv` has the columns `sample`, `tool`, `score_name`, `score`, `call`, with 20 rows per sample. `<k>` is 1 to 6.

| `score_name` | Value | Direction | `call` |
| --- | --- | --- | --- |
| `SID<k>_Score` | SenCID's `SID_Score` of model SID`<k>`: the probability of class 1 from a logistic regression on the genes of the model. 0 to 1. | Higher means more senescent. | SenCID's `Binarization` of model SID`<k>`: `senescent` for 1, `not senescent` for 0. |
| `SID<k>_Decision` | SenCID's `Decision` of model SID`<k>`: the decision value of a linear support vector machine on the same genes. Centred at 0, no fixed range. | Higher means more senescent. | Empty. |
| `rec_SID<k>` | SenCID's `rec_SID<k>`: the decision value of the recommendation model for SID`<k>`. | It says how well the sample fits SID`<k>`. It does not measure senescence. | Empty. |
| `RecSID` | The number `<k>` of the model that SenCID's recommendation model names for the sample (SenCID's column `RecSID`). | A class, not a measure. | Empty. |
| `RecSID_Score` | The `SID<k>_Score` of the model named in `RecSID`. **This is the main result.** | Higher means more senescent. | The `call` of that `SID<k>_Score`. |

Sources for the meaning and the direction: the docstring of `SenCID/api.py` ("SID scores ranging from 0 to 1", "Decision Function values similar to SID scores, but without the ranges and centered by 0", "crude Senescent(1)-Or-Not(0) judgements by the threshold 0.5"), and issue 1, where the author writes that the logistic regression was trained "with the aim of normalizing the SID scores to all fall within a range of 0 to 1, enabling interpretation as probabilities of senescence". Checked in the model files: every SID model is a scikit-learn `SVC` with a linear kernel plus a `LogisticRegression`, both with the classes 0 and 1; the recommendation model is an `SVC` with a linear kernel and the classes `SID1` to `SID6`.

Why `RecSID_Score` is the main result: the tutorial says "the SID with the highest recommendation score is considered to be most suitable for senescence evaluation". In issue 4 the author advises to run all six models and to "separate the different cell types into different SIDs according to the 'recSID' results" and to "refer their respective SID scores". In issue 1 the author writes that `rec_SID6` only says whether a sample belongs to SID6 and that the `SID_Score` of that model is the measure of senescence. On the test input `RecSID` was the model with the highest `rec_SID<k>` in 34 of 34 samples.

Other files in `--out-dir`:

- `run_info.json`: commit of SenCID, the command that was run, the number of genes of the gene list found in the input and the number missing with their names, the missing genes of each model, the genes with a count of zero in all samples, the seconds SenCID took, and the package versions.
- `sencid_output/`: the seven files that SenCID wrote, unchanged.
- `sencid_log.txt`: what SenCID printed.

## Deviations from the authors' use

1. **Denoising is off.** The default of the code is on. The author prescribes off for bulk RNA-seq (see the table above), so the wrapper passes `--denoising f`. With denoising on, the scores would depend on a random fit to the input (issue 13).
2. **The recommendation is applied to each sample.** The tutorial chooses one model for "one single cell type" from the recommendation scores of all its cells. The wrapper is given no cell type, so each sample gets the score of the model named for that sample. Samples of one study can therefore be scored by different models: on the test input 24 samples were assigned SID5 and 10 samples SID4. The repository does not say whether the scores of two different models can be compared with each other. All six scores and `RecSID` are in `scores.tsv`, so a rule with one model per study or per cell type can be computed from this output without a new run.
3. **Gene symbols come from the column `Symbol` of the GEO annotation, and nothing else.** Ten of the 1,290 genes of the gene list are absent from that column, because the gene has a newer symbol there: ARNTL (now BMAL1), C12orf45 (NOPCHAP1), CCDC58 (MIX23), DDX58 (RIGI), FAM214A (ATOSA), FAM214B (ATOSB), GBA (GBA1), PHB (PHB1), TMEM251 (LYSET) and WDR63 (DNAI3). SenCID's own code gives a missing gene the value 0.01 (step 2). Affected model genes: CCDC58 and FAM214B in SID3, WDR63 in SID5, and DDX58, TMEM251 and WDR63 in SID6; all ten are read by the recommendation model. The wrapper does not search the column `Synonyms`. Size of the effect, measured once on the test input without labels: with the nine symbols that have exactly one match in `Synonyms` restored, `RecSID_Score` changed by at most 0.044, `RecSID` and the `call` of `RecSID_Score` changed in no sample, the `call` of `SID2_Score` changed in 2 of 34 samples and the `call` of `SID3_Score` in 1 of 34.
4. **Three gene ids share the symbol TRNAV-CAC.** The wrapper passes all three rows on. SenCID's own reader renames repeated names. The symbol is not in the gene list.
5. **The channels of `SenCID.yaml` are reduced to conda-forge.** Every conda package pinned there is a conda-forge build. The pin `ca-certificates=2022.6.15` resolved to version 2022.6.15.2 with the same build string.
6. **The home folder is replaced for the SenCID process.** On import, the packages `kopt` and Keras read and write settings in the home folder. The wrapper points `HOME` to a temporary folder inside `--out-dir` and removes it afterwards, so that nothing is written to the home folder and no existing Keras setting is read.
7. **At most 8 threads** for the numeric libraries (`config.sh`).

## Released behaviour that the wrapper keeps

- SenCID's command line catches every error of a file, prints it and exits with code 0 (`SenCID/__main__.py`, lines 92 to 102). The wrapper does not change this. It checks that all seven output files exist and stops otherwise.
- Counts per million are computed over the 1,290 genes of the gene list (step 3), not over all genes.

## Data the released models were fitted on

- The repository holds no list of training data and no training code (issues 3 and 23 ask for it).
- The abstract of the paper says the tool was "trained on 602 samples from 52 senescence transcriptome datasets spanning 30 cell types". The list is Table S1 of the paper; `kb/identification/tool_training_datasets.md` and `.tsv` hold it with the role of each dataset (250 samples of 20 datasets are marked as used for training).
- Read from the released model files: the six support vector machines were fitted on 39, 6, 45, 21, 83 and 18 samples (SID1 to SID6; 212 in total), and the recommendation model on 572 samples and 1,290 genes. I did not reconcile 212 with the 250 samples that Table S1 marks. The recommendation model was fitted on more samples than Table S1 marks as used for training, so it has probably seen samples of datasets that the table marks as not used for training. This is my inference from the two numbers.
- The test input of this wrapper, GSE130727, is one of the training datasets (32 of its 37 samples are marked in Table S1). The test run shows that the wrapper works and nothing about how well SenCID separates senescent samples.

## Limits

- The Methods of the paper were not read. If they prescribe something that the repository and the issues do not state, the wrapper does not follow it.
- Human bulk RNA-seq raw counts only. Normalised values are refused by the wrapper.
- The author advises against comparing scores of bulk samples and of single cells directly, and refers to the paper (Figure S6D to H) for the choice of the threshold 0.5 (issue 1).
- Whether SenCID separates quiescent from senescent cells is asked in issue 25 and not answered.
- The training data are cultured cells (see `kb/identification/tools.tsv`).
- The wrapper was tested on one study (34 samples).

## Environment

`bash install_env.sh` creates `/SCRATCH/dtzim01/senescence/tools/envs/sencid` on biogemt-1 with micromamba 2.3.0 from `environment.yml`, copies the clone `/SCRATCH/dtzim01/senescence/tools/SenCID` at its commit to a temporary folder, installs SenCID from that copy with `pip install --no-deps --no-build-isolation`, runs `pip check`, compares the installed code and model files with the clone, writes the commit to `sencid_source.json` in the environment, and writes the two locks. The clone is not changed.

- Built on 2026-11-15 on biogemt-1 from conda-forge and PyPI. `pip check` reported no broken requirement, and the installed files equal the clone. Size 2.3 GiB.
- `environment.pip.lock` leaves out SenCID itself, because pip would list the temporary folder; SenCID is pinned by the commit in `config.sh`.
- During the export micromamba printed `ERROR: unknown command "inspect"` (pip 21.2.4 has no `inspect` command). The lock was written.
- SenCID's own tutorial command with denoising on (`SenCID --filepath . --denoising t --fileclass txt --binarize t`) ran in this environment on the demo file of the repository (200 cells, 16 s). Its scores were not compared with anything.

## Test

`/SCRATCH/dtzim01/senescence/runs/20261115-0010-tool-test-sencid/` on biogemt-1, started through the wrapper on GSE130727 (34 samples): 7 s in total, of which SenCID took 3 s. `scores.tsv` has 680 rows, one for each of the 34 samples and 20 score names; every score is finite; every score name has 34 different values except `RecSID`, which has two. 1,280 of the 1,290 genes of the gene list were found. A second run gave an identical `scores.tsv`. The scores were not compared with sample titles or labels.
