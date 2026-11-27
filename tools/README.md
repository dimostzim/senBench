# The published tools

`bash install.sh <tool>` installs one tool and `bash run.sh <tool>` scores `test` with it; the tool is `hUSI`, `SenCID`, `SenCat`, `SENCAN`, `markeR` or `SenFinder`. The tools' own code is not in this repository.

## Install

`tools/<tool>/setup.sh` clones the tool's repository at a fixed commit into `build/tools/` and builds its environment in `build/envs/<tool>` with micromamba 2.3.0, from the lock files of the tool's folder. It needs git, curl and tar, and access to GitHub, conda-forge, bioconda and PyPI (hUSI also to Posit's CRAN snapshot of 2023-12-15). It stops when the clone is at another commit or holds a changed file, or when the installed packages differ from the lock files. The six tools take about 9 GB.

| Tool | Paper | Repository, commit | Licence of the tool | Installed after the conda packages |
| --- | --- | --- | --- | --- |
| hUSI | Wang 2025, Nat Aging | `WJPina/HUSI`, `300da32` | Apache-2.0 | dplyr, mclust and jsonlite from source |
| SenCID | Tao 2024, Cell Metab | `JackieHanLab/SenCID`, `f93cba1` | MIT | the pip packages of `environment.pip.lock`, then SenCID from the clone |
| SenCat | Anerillas 2026, Mol Cell | `maragkakislab/wf-ml-markers-senescence`, `0e23fcc` | none in the repository | four pip packages (pydeseq2 0.5.4 and three it needs) |
| SENCAN | Jochems 2021, Cell Rep | the classifier file `Data/SENCAN_classifier.rda` of `WJPina/HUSI`, `300da32` | none given for SENCAN | nothing |
| markeR | Martins-Silva 2026, NAR Genom Bioinform | `DiseaseTranscriptomicsLab/markeR`, `c6528f0` | Artistic-2.0 | markeR from the clone |
| SenFinder | Mahmud 2024, GeroScience | `XiaoDongLab/SenFinder`, `5f56506` | AGPL-3.0 | nothing |

## Run

```
bash tools/<tool>/run_tool.sh <counts.tsv.gz> <annotation.tsv.gz> <output folder>
```

- **Counts:** a tab-separated table with the column `GeneID` (NCBI gene id) and one column of raw read counts per sample.
- **Annotation:** GEO's gene annotation table, `downloads/geo_uniform_counts/annotation/Human.GRCh38.p13.annot.tsv.gz`.
- **Output folder:** it must not exist. The script writes `scores.tsv` into it, with the columns `sample`, `tool`, `score_name`, `score` and `call`. It stops without `scores.tsv` when a count is not a whole number of zero or more, or when a score is not a finite number.

`pipeline/predict.py <tool>` calls `run_tool.sh` once per `test` sample, on a count table that holds that sample only, and writes the score named below to `predictions/<tool>.tsv`. A sample on which `run_tool.sh` stops gets no row. hUSI is called once per dataset, because its released score needs two samples; the score taken from it depends on the sample alone. The files of every call are in `build/predict/<tool>/`.

## The score of each tool, and what differs from the tool's published use

| Tool | Score in `predictions/` (higher is more senescent) | What `run_tool.sh` changes against the tool's published use |
| --- | --- | --- |
| hUSI | `hUSI_unscaled`: the correlation of the sample's TPM values with the tool's gene weights, before the tool rescales it over the samples of the input | TPM is computed from the counts and GEO's gene lengths; the tool releases only the rescaled score; 9 of its 13,346 genes are not in the annotation table |
| SenCID | `RecSID_Score`: the score of the one of the tool's six models that the tool recommends for the sample | Denoising is off, as the first author prescribes for bulk RNA-seq; 10 of the tool's 1,290 genes are absent and get the tool's fill value |
| SenCat | `rna_score`: the weighted sum of the tool's 5,000 genes | NCBI gene ids are converted to Ensembl gene ids (120 genes are never found); on one sample the tool's filter keeps only genes with a count above 380 |
| SENCAN | `SENCAN_score`: the probability of senescent | The counts come from NCBI genes, not from kallisto; 18 of the 137 genes are absent and are filled from the tool's reference samples; one sample is not adjusted to the reference sample |
| markeR | `HernandezSegura:ssGSEA`: the gene set of Hernandez-Segura 2017 scored with the tool's ssGSEA | Each sample is scored alone, without the paper's correction across datasets; a gene is kept when its count is above 70; the script stops when any of the tool's nine scores is not finite (42 `test` samples) |
| SenFinder | `Senescence_probability`: the probability of senescent of the released model | TPM is computed from the counts and GEO's gene lengths; the one of the 20 genes that has no NCBI gene id is set to 0, as the tool's README prescribes for a missing value |
