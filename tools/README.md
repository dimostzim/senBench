# The published tools

This folder runs each published tool that senBench compares with Agentomics, in the way the benchmark ran it. Each tool has one folder, which holds:

- `setup.sh`: clones the tool's own repository at the commit the benchmark used and builds the wrapper's environment.
- The benchmark's wrapper of the tool: `run_tool.sh`, the scripts that convert the counts and call the tool's unchanged code, and the lock files of the environment. Every wrapper file equals the copy that the benchmark ran; `SHA256SUMS` gives their sha256 (`sha256sum --check SHA256SUMS` in this folder).
- `predictions.tsv`: the tool's score of every labelled `test` sample that it scored (columns `dataset`, `sample`, `score`), the scores of the results table.

The tools' own code is not in this repository; `setup.sh` downloads it, and each tool keeps its own licence.

| Tool | Folder | Paper | Repository, commit | Licence of the tool | Score in senBench |
| --- | --- | --- | --- | --- | --- |
| hUSI | `husi/` | Wang 2025, Nat Aging | `WJPina/HUSI`, `300da32` | Apache-2.0 | `hUSI_unscaled` (not a released score) |
| SenCID | `sencid/` | Tao 2024, Cell Metab | `JackieHanLab/SenCID`, `f93cba1` | MIT | `RecSID_Score` |
| SenCat | `sencat/` | Anerillas 2026, Mol Cell | `maragkakislab/wf-ml-markers-senescence`, `0e23fcc` | none in the repository | `rna_score` |
| SENCAN | `sencan/` | Jochems 2021, Cell Rep | the classifier file in `WJPina/HUSI`, `300da32` | none given for SENCAN; the hUSI repository is Apache-2.0 | `SENCAN_score` |
| markeR | `marker/` | Martins-Silva 2026, NAR Genom Bioinform | `DiseaseTranscriptomicsLab/markeR`, `c6528f0` | Artistic-2.0 | `HernandezSegura:ssGSEA` |
| SenFinder | `senfinder/` | Mahmud 2024, GeroScience | `XiaoDongLab/SenFinder`, `5f56506` | AGPL-3.0 | `Senescence_probability` |

For every score in the table, a higher value means more senescent.

## Scoring `test` with a tool

After `bash build.sh`, from the repository's folder, for SenCID:

```
bash tools/sencid/setup.sh
python3 pipeline/predict.py sencid
```

The same two commands, with the tool's folder name, run each of the other tools. `python3 pipeline/score.py` then scores the new predictions with the others.

`setup.sh` installs into `build/tools/`, or into the folder given as its argument; the six tools share one install folder. It needs Linux on x86-64, bash, git, curl and tar, and access to GitHub, conda-forge, bioconda and PyPI (hUSI also to a dated CRAN snapshot of Posit). Proxy settings are taken from the shell. It writes:

- `envs/env.sh` and `envs/bin/micromamba` (micromamba 2.3.0, the version that built the benchmark's environments, checked by sha256), the package cache and temporary files under `envs/`;
- the clone of the tool (`HUSI`, `SenCID`, `SenFinder`, `markeR` or `wf-ml-markers-senescence`), holding only the pinned commit;
- the environment `envs/<tool>`, built from the wrapper's explicit lock, which names the exact build of every conda package, followed by the steps the benchmark's own setup took (R packages from source for hUSI, markeR from its clone, the pip packages of SenCID and SenCat, SenCID from its clone);
- for SenCID only, `wrappers/sencid/`: a copy of the wrapper's three scripts whose `config.sh` names this install folder (the wrapper reads its paths from that file; the other two scripts are unchanged).

`setup.sh` stops when a clone is at another commit or holds a changed file, when the installed R or pip packages differ from the wrapper's lock, or when the environment exists (delete it to build it again). All six together take about 9 GB, 5 GB of it the package cache.

`pipeline/predict.py <tool>` calls the wrapper once per labelled sample of `build/all980_v10/split/test_samples.tsv`, on a count table that holds only that sample, cut from the round's count table of its dataset (`build/all980_v10/counts/`), so that every score is computed from one sample alone. hUSI is the exception (below). It runs 8 wrapper calls at a time (`--workers`), each with one thread for the numeric libraries, as the benchmark did, and writes the tool's `predictions.tsv`, replacing the shipped file; the inputs, the outputs and the log of every call are in `build/predict/<tool>/`, which must not exist yet. A sample on which the wrapper stops gets no row. The wrappers' default paths of their environments and clones name the folders of our build; `predict.py` overrides them with `--env-prefix` and `--tool-repository` (hUSI, SENCAN, SenFinder), with `--env-prefix` (markeR), with the variable `SENESCENCE_TOOLS_DIR` (SenCat), and with the copy of the wrapper that `setup.sh` writes (SenCID).

## The wrappers' input and output

Every wrapper takes the same arguments:

```
bash run_tool.sh --counts <counts.tsv.gz> --annotation <Human.GRCh38.p13.annot.tsv.gz> --out-dir <new folder>
```

- `--counts`: a tab-separated table, gzip or plain: the first column `GeneID` (NCBI gene id), then one column of raw read counts (whole numbers, zero or more) per sample.
- `--annotation`: GEO's gene annotation table, `https://ftp.ncbi.nlm.nih.gov/geo/rnaseq_counts/annot/Human.GRCh38.p13.annot.tsv.gz`; `download.sh` fetches it with the benchmark's data, into `downloads/geo_uniform_counts/annotation/`.
- `--out-dir`: a folder that does not exist yet.

The wrapper writes `scores.tsv` (columns `sample`, `tool`, `score_name`, `score`, `call`) and `run_info.json` (the tool's commit, the genes found and missing, the package versions). It receives no labels. It stops with a non-zero exit code and writes no output when the counts are not whole numbers of zero or more, or when the tool returns a score that is not a finite number. In the benchmark such a stop is a result: the sample has no score from that tool.

## hUSI

hUSI, the human universal senescence index, is the tool of Wang 2025 (Nat Aging, DOI 10.1038/s43587-025-00886-2). Its weights come from a one-class logistic regression fitted on 385 senescent bulk RNA-seq samples of cultured human cells: one weight for each of 13,346 genes. The score of a sample is the Spearman correlation, over the genes it shares with the weights, between the weights and the sample's TPM values. The released function then rescales the correlations of the samples of its input from 0 (the lowest) to 1 (the highest).

- Repository: `https://github.com/WJPina/HUSI`, commit `300da32f111f8c8603e822eca8e84c9b5a76db48` (2025-09-11). Licence: Apache-2.0.
- Needs: TPM values of genes named by symbol. The wrapper computes TPM from the counts and the gene lengths of GEO's annotation table, and reads the weights and functions of the clone (`sc/hUSI.R`).
- Score in senBench: `hUSI_unscaled`, the correlation before hUSI's rescaling. It is not a released score. The released function `cal_hUSI` returns only the rescaled score `hUSI`; the wrapper also writes the correlation, computed with the authors' function `scoreOCLR` (`R/functions_new.R`), which their analysis scripts call before rescaling. `hUSI_unscaled` of a sample depends on that sample alone; `hUSI` depends on the other samples of the input, because the lowest and the highest sample define 0 and 1.
- One sample: the wrapper refuses a count table with one sample, because the released score cannot be rescaled over one sample. The benchmark therefore took `hUSI_unscaled` from the wrapper's runs on the count table of each dataset with the labelled samples of that dataset, and `predict.py` calls it in the same way, on the round's count table of each `test` dataset. A sample's `hUSI_unscaled` is the same whichever other samples are given (checked on GSE130727: 0.04712 for its first sample with either of two other samples).

What the wrapper changes against the tool's published use:

1. TPM is computed in the wrapper from the counts and the gene lengths of the annotation table. The authors used TPM from StringTie or from the data sources. On GSE130727 the wrapper's TPM equals GEO's own TPM table to the four digits GEO prints.
2. TPM values are scored as they are. In one analysis the authors first log-transform, scale and centre each gene over pooled samples; their applications to bulk data do not.
3. 9 of the 13,346 genes with a weight are not in the annotation table under their symbol.
4. The wrapper writes the second score, `hUSI_unscaled`; the released function returns only `hUSI`.
5. `call` is the class of the tool's function `SSE_hUSI`. The tool's other classifier, `GMM_hUSI`, is not run.
6. R packages are built at the versions the repository states (R 4.0.5, dplyr 1.1.4, mclust 6.0.1). Seurat, ggplot2 and pROC are not installed; the scoring functions do not call them.

## SenCID

SenCID is the tool of Tao 2024 (Cell Metab, DOI 10.1016/j.cmet.2024.03.009). It holds six models, SID1 to SID6, one for each of six senescence identities found in the authors' data; each scores a sample from 0 to 1 with a logistic regression on 13 to 224 genes. A seventh model recommends which of the six to use for the sample, from 1,290 genes.

- Repository: `https://github.com/JackieHanLab/SenCID`, commit `f93cba19b9d44e84f89504ecc6aaeb108ce8e508` (2024-02-22, version 1.0.0). Licence: MIT.
- Needs: raw counts of genes named by HGNC symbol, with denoising switched off for bulk RNA-seq (the first author's advice in the repository's issues). SenCID normalises the counts itself: it keeps the genes of its list of 1,290 that have a total count of at least 1 over the samples of the input, gives every other gene of the list the value 0.01, divides by the sum over the 1,290 genes, takes log2(value + 1) and turns the values of each sample into z-scores. The environment keeps every version pin of the repository (Python 3.6.15, TensorFlow 2.2.0, scikit-learn 0.24.2).
- Score in senBench: `RecSID_Score`, the score of the model that SenCID recommends for the sample, as the tool's tutorial uses it; `call` is SenCID's own call at 0.5. The wrapper also writes the scores of all six models.
- One sample: with one sample, every gene of the list with no read in that sample gets the value 0.01; with several samples, only genes with no read in any of them do.

What the wrapper changes against the tool's published use:

1. Denoising is off. The default in the code is on; the first author prescribes off for bulk RNA-seq in the repository's issues.
2. The tool recommends one of its six models. The wrapper takes the recommendation per sample, because it is given no cell type; the tutorial picks one model per cell type.
3. 10 of the tool's 1,290 genes are absent under their current symbol. SenCID's own code fills them with 0.01.
4. The environment keeps every version pin of the repository; packages come from conda-forge and PyPI only.
5. `HOME` points to a temporary folder for the SenCID process, because two of its libraries write settings to the home folder.

The paper's STAR Methods agree with these choices: for bulk samples, each sample gets the model with the highest recommendation index, and 0.5 is the threshold for a call.

## SenCat

The SenCat RNA score is the score of Anerillas 2026 (Mol Cell, DOI 10.1016/j.molcel.2026.05.017). Its repository releases fixed weights: one coefficient for each of 5,000 Ensembl genes. The score of a sample is the sum, over the genes found in the input, of the coefficient times the logarithm of the gene's normalised count plus one. It has no threshold and no class call.

- Repository: `https://github.com/maragkakislab/wf-ml-markers-senescence`, commit `0e23fcc926f31fc2206631ff9b64a852fb81fe89` (2026-08-19). The archived release cited by the paper is on Zenodo (DOI 10.5281/zenodo.20008425). Licence: the repository has no licence file; the Zenodo record is CC BY 4.0.
- Needs: counts by Ensembl gene id, normalised with the tool's `normalize_counts.py`, which removes genes with a count in at most 5% of the samples or a total count of at most 380 and divides each sample by its DESeq2 size factor. The wrapper converts NCBI gene ids to Ensembl gene ids, then runs the tool's own `normalize_counts.py` and `marker_classifier.py`.
- Score in senBench: `rna_score`, the only score. A higher value means more similar to the senescence signature (the repository's README).
- One sample: the DESeq2 size factor is 1, and the gene filter keeps the genes with a count above 380 in that sample.

What the wrapper changes against the tool's published use:

1. NCBI gene ids are converted to Ensembl gene ids with the annotation table. 120 of the 5,000 genes can never be found with this table. Counts of NCBI genes that share an Ensembl id are added.
2. The repository clone replaces the archived release that the paper cites, because the archived script stops on samples without labels. The weights file is the same in both.
3. One environment at current package versions replaces the authors' two.
4. The tool's gene filter and size factors are computed over the samples of the input.

## SENCAN

SENCAN is the classifier of Jochems 2021 (Cell Rep, DOI 10.1016/j.celrep.2021.109441): an elastic net on 137 genes that gives each sample a probability of being senescent. The paper reports no code. The classifier file in the hUSI repository, `Data/SENCAN_classifier.rda`, equals the paper's Table S4 and the file the authors' website offers (same sha256), and the wrapper reads it from there.

- Classifier file: `https://github.com/WJPina/HUSI`, commit `300da32f111f8c8603e822eca8e84c9b5a76db48`. Licence: none given for SENCAN; the hUSI repository is Apache-2.0.
- Needs: read counts per Ensembl gene id. The classifier normalises the input samples against its reference samples, takes log counts per million, centres and scales them, and fills a gene that the input lacks with the mean count of the two reference samples that correlate best with the sample. The authors' scoring script (`SENCAN_example.r` on their website) takes these steps, and the wrapper follows it.
- Score in senBench: `SENCAN_score`, the probability of senescent from the elastic net, from 0 to 1; `call` is senescent above 0.5.
- One sample: all samples of a table are normalised together, so with one sample the normalisation factor is 1 and the sample is not adjusted towards the authors' reference sample.

What the wrapper changes against the tool's published use:

1. The paper quantifies raw reads with kallisto on GENCODE 34. The wrapper starts from read counts per NCBI gene, so this step cannot be followed.
2. Genes that the input lacks are filled per sample with the mean count of the two reference samples that correlate best with it, as the authors' script does. GEO's count tables lack 18 of the 137 genes, which carry 21.6% of the summed absolute coefficients, among them the largest positive and the most negative coefficient.
3. The filled-in counts are not rescaled to the depth of the sample, as in the authors' script.
4. All samples of a count table are normalised together, as in the authors' script. A sample's score therefore changes with the other samples of the table.
5. Genes are matched by Ensembl gene id (119 of 137 found), as in the authors' script.
6. `call` is senescent above 0.5, the dashed line of the paper's Figure 6. The text of the paper states no threshold.
7. R 4.0.5 and edgeR 3.32.1 replace the paper's R 3.6.0 and edgeR 3.26.8.

## markeR

markeR is an R package of Martins-Silva 2026 (NAR Genom Bioinform, DOI 10.1093/nargab/lqag057) that scores gene sets in bulk RNA-seq samples. It has no trained model and no rule that calls a sample senescent. It ships three senescence gene sets (`HernandezSegura`, `REACTOME_CELLULAR_SENESCENCE` and `LiteratureMarkers`) and three scoring methods (`logmedian`, `ssGSEA` and `ranking`); the wrapper writes all nine scores of a sample.

- Repository: `https://github.com/DiseaseTranscriptomicsLab/markeR`, commit `c6528f0fa1570602b077e4688865dcd68af3d058` (2026-07-08, package version 1.3.1). Licence: Artistic-2.0.
- Needs: filtered and normalised counts of genes named by symbol. The wrapper names the genes by the symbol of GEO's annotation table, keeps the genes with a mean raw count above 70 over the samples of the input, normalises with edgeR's TMM factors and limma's `voom`, as the authors' preprocessing does, and calls `markeR::CalculateScores` with the package's gene sets.
- Score in senBench: `HernandezSegura:ssGSEA`, the gene set of Hernandez-Segura 2017 (Curr Biol; 55 genes, 32 up and 23 down in senescence; 53 of them are in the annotation table) scored with markeR's own ssGSEA, which ranks the genes within the sample. The score of the down genes is subtracted from that of the up genes, so a higher score means more senescent. The wrapper's main score is `HernandezSegura:logmedian`, but `logmedian` centres every gene on its median over the samples of the input and is 0 for every gene set when the input holds one sample; the benchmark therefore scores `ssGSEA`.
- One sample: the expression filter keeps the genes with a count above 70 in that sample. When one direction of a gene set keeps no gene after the filter, its `ssGSEA` score is not a finite number, and the wrapper stops without output, because it stops when any of its nine scores is not finite. On the 687 `test` samples it stopped 42 times, so markeR is scored on 645 samples.

What the wrapper changes against the tool's published use:

1. The authors keep a gene when its mean count is above 70 in at least one condition. The wrapper sees no conditions, so it keeps a gene when its mean count over all samples of the input is above 70. This is stricter.
2. Each count table is scored alone. The paper processed 545 samples of 25 datasets together and corrected for the dataset.
3. Counts of gene ids that share a symbol are added (one case).
4. Two genes of `HernandezSegura` (`FAM214B`, `MT-CYB`) are not in the annotation table under these symbols.
5. markeR 1.3.1 replaces the paper's 1.0.0; the scoring files are the same in both.
6. The authors name no main score; the wrapper takes `HernandezSegura:logmedian` (above for the benchmark's choice).

Kept as released, although the documentation says otherwise: `ranking` divides by the number of genes in the matrix; `ssGSEA` is the package's own function; the default method in the code is `ssGSEA`. The wrapper stops for the whole input when any of its nine scores is not a finite number for any sample.

## SenFinder

SenFinder is the released model of Mahmud 2024 (GeroScience, DOI 10.1007/s11357-024-01485-6): a linear support vector machine on the TPM values of 20 genes, which scales each gene with the means and standard deviations of the authors' training samples stored in the model.

- Repository: `https://github.com/XiaoDongLab/SenFinder`, commit `5f56506300baab2ed4c20af5a52c6252c4578b50` (2024-12-07). Licence: AGPL-3.0.
- Needs: TPM values, not log-transformed, of the 20 genes, with no missing value (the README: use 0 instead). The wrapper computes TPM from the counts and the gene lengths of GEO's annotation table.
- Score in senBench: `Senescence_probability`, the probability of the class senescent that the released model returns, from 0 to 1. The released script writes only the class; the wrapper writes both.
- One sample: SenFinder's score of a sample does not depend on the other samples: on the `test` samples of an earlier round, the scores of one sample per call equalled those of the samples scored together (largest difference 1e-15).

One of the 20 genes, `Lnc_PLOD2_2` (ENSG00000261051), has no NCBI gene id, so a count table of the benchmark can never supply it. The benchmark runs the wrapper with `--set-absent-genes-to-zero`, which sets it to 0, as the repository's README prescribes for a missing value ("19 of 20 genes"). Every sample then moves by the same amount towards not senescent, so the ranking of the samples comes from the other 19 genes, and the calls are biased towards not senescent.

What the wrapper changes against the tool's published use:

1. The absent gene is set to 0 (above).
2. TPM is computed in the wrapper from the counts and gene lengths. The authors used RSEM after STAR on GRCh37.
3. The model's column `YIU2` is read as the gene YJU2.
4. The wrapper writes the probability of the class senescent in addition to the class.

## Sources

The benchmark documented each wrapper in a README of its own, with the source of every statement above and the data each tool was fitted on; this repository no longer ships these READMEs, nor the wrappers' `environment.yml` and setup scripts, which `setup.sh` replaces. Commit `6aece08` holds them, for example `git show 6aece08:tools/sencid/wrapper/README.md`. Those files, and the comments in some wrapper files, name the folders of our build host and the setup scripts by which the benchmark built its environments.
