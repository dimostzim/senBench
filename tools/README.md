# The published tools

This folder runs each published tool that senBench compares with Agentomics, in the way the benchmark ran it. Each tool has its own folder:

- `wrapper/`: the benchmark's wrapper of the tool: `run_tool.sh`, the scripts that convert the counts and call the tool's unchanged code, the environment files and the locks. Every file equals the copy that the benchmark ran; `SHA256SUMS` gives their sha256 (`sha256sum --check SHA256SUMS` in this folder). The wrapper's own `README.md` documents it as written for the benchmark, with the source of every statement.
- `setup.sh`: clones the tool's own repository at the commit the benchmark used and builds the wrapper's environment.
- `README.md`: the tool, what it needs, how to score one sample, what the wrapper changes against the tool's published use, and the score that senBench reports.

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

## Setting up

```
bash tools/<tool>/setup.sh /path/to/tools_install
```

The six tools share one install folder. `setup.sh` needs Linux on x86-64, bash, git, curl and tar, and access to GitHub, conda-forge, bioconda and PyPI (hUSI also to a dated CRAN snapshot of Posit). Proxy settings are taken from the shell. It writes:

- `envs/env.sh` and `envs/bin/micromamba` (micromamba 2.3.0, the version that built the benchmark's environments, checked by sha256), the package cache and temporary files under `envs/`;
- the clone of the tool (`HUSI`, `SenCID`, `SenFinder`, `markeR` or `wf-ml-markers-senescence`), holding only the pinned commit;
- the environment `envs/<tool>`, built from the wrapper's explicit lock, which names the exact build of every conda package, followed by the steps the wrapper's own setup takes (R packages from source for hUSI, markeR from its clone, the pip packages of SenCID and SenCat, SenCID from its clone);
- for SenCID only, `wrappers/sencid/`: a copy of the wrapper whose `config.sh` names this install folder (the wrapper reads its paths from that file; every other file of the copy is unchanged).

The wrappers' default paths of their environments and clones, and the commands in the wrappers' own READMEs, name the folders of our build. The commands in each tool's README override them: with `--env-prefix` and `--tool-repository` (hUSI, SENCAN, SenFinder), with `--env-prefix` (markeR), with the variable `SENESCENCE_TOOLS_DIR` (SenCat), and with the copy of the wrapper that `setup.sh` writes (SenCID). `env.sh` does not change these paths; it places micromamba, the package cache and the temporary files in the install folder.

`setup.sh` stops when a clone is at another commit or holds a changed file, when the installed R or pip packages differ from the wrapper's lock, or when the environment exists (delete it to build it again). All six together take about 9 GB, 5 GB of it the package cache. With 16 CPUs the six were set up in 6 minutes, SenCID taking the longest (3.4 minutes).

## Input and output

Every wrapper takes the same arguments:

```
bash run_tool.sh --counts <counts.tsv.gz> --annotation <Human.GRCh38.p13.annot.tsv.gz> --out-dir <new folder>
```

- `--counts`: a tab-separated table, gzip or plain: the first column `GeneID` (NCBI gene id), then one column of raw read counts (whole numbers, zero or more) per sample. The benchmark gave each tool one sample per call, a table with the columns `GeneID` and that sample, so that each score depends on its sample alone (hUSI is the exception, see `husi/README.md`). One sample's table can be cut out of a dataset's count table, for example `zcat GSE130727.tsv.gz | cut -f1,2 | gzip > one_sample.tsv.gz` for its first sample.
- `--annotation`: GEO's gene annotation table, `https://ftp.ncbi.nlm.nih.gov/geo/rnaseq_counts/annot/Human.GRCh38.p13.annot.tsv.gz`; `download.sh` fetches it with the benchmark's data.
- `--out-dir`: a folder that does not exist yet.

The wrapper writes `scores.tsv` (columns `sample`, `tool`, `score_name`, `score`, `call`) and `run_info.json` (the tool's commit, the genes found and missing, the package versions). It receives no labels. It stops with a non-zero exit code and writes no output when the counts are not whole numbers of zero or more, or when the tool returns a score that is not a finite number. In the benchmark such a stop is a result: the sample has no score from that tool.

## Smoke test

On a Linux server with 16 CPUs, each `setup.sh` was run into a new install folder and each wrapper on one sample of GSE130727, a dataset of the training data (2026-11-26). Each `scores.tsv` was identical byte for byte to the one that the benchmark's own install gives on the same table, and the clones were unchanged afterwards.
