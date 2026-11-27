# senBench

senBench is a benchmark of methods that call a bulk RNA-seq sample of cultured human cells senescent or not senescent, from that sample alone. It tests them on studies from laboratories that no method was trained on. It holds public RNA-seq datasets with a label for every sample, a split of the studies by laboratory, and the scores of seven models on the held-out part `test`: Agentomics, a model that the agent system Agentomics-ML built on senBench's training data, and six published tools.

## Quick start

On Linux (x86-64), from the repository's folder:

```
bash download.sh
bash build.sh
python3 pipeline/predict.py agentomics
bash tools/sencid/setup.sh
python3 pipeline/predict.py sencid
python3 pipeline/score.py
```

1. `download.sh` downloads every public file the build needs into `downloads/` (18.76 GiB; python3 and curl).
2. `build.sh` builds the count tables, the sample table, the split, the training data, `test` and the dataset in Agentomics-ML's layout into `build/`, and compares every output with the sha256 of our own build (it creates its conda environment first; `--conda mamba` or `--conda micromamba` if conda is not installed, `--threads` for more than 8 CPUs).
3. `predict.py agentomics` scores `test` with the Agentomics model and writes `agentomics_model/predictions.tsv` (it creates the model's environment, 7.2 GB, the first time).
4. `setup.sh` installs a published tool into `build/tools/`, and `predict.py` with the tool's folder name (`husi`, `sencid`, `sencat`, `sencan`, `marker`, `senfinder`) scores `test` with it and writes `tools/<tool>/predictions.tsv`.
5. `score.py` scores the seven prediction files and writes `results.tsv` (it needs numpy and pandas; `build/environments/count-tables/bin/python`, which `build.sh` creates, has them).

Each step can run alone: `score.py` on the shipped predictions gives the shipped `results.tsv` byte for byte, without any download.

## The task

A method receives the raw read counts of one sample and returns a score that is higher for senescent samples. A sample is `senescent` when the study made the cells senescent, and `not senescent` when it is a control of the study: proliferating, quiescent, early passage, untreated, vehicle or empty vector cells. The score must be computed from the sample alone, so that it does not depend on which other samples a study holds. Each tool was therefore given a count table with one sample, except hUSI, whose score is a correlation computed on each sample alone but whose wrapper needs at least two samples (`tools/README.md`).

## The data

The benchmark holds 3,573 labelled samples (1,860 senescent, 1,713 not senescent) of 318 datasets from GEO and ArrayExpress. Their counts are GEO's uniformly processed counts, the authors' own count tables, or, for 24 datasets, counts that we made from the raw reads. Language models read the labels from the public records of each dataset. The labels and our counts of the raw reads are shipped in `pipeline/data/` and not made again; everything else is built from the downloaded records and count files. `pipeline/README.md` describes how the data are made, how we counted the raw reads, the label rules, the split and the sizes.

## The split

Studies are grouped into laboratory groups: two studies are in one group when their records share the contact person and institute, the last contributor, a PubMed id, a BioProject, a SuperSeries, or two contributors. Four more links were made by hand from the records, three of them found by a check for people shared between the parts. Set 1 holds the laboratory groups with no dataset in the list of data that any published method was trained, validated or tested on, or made from. Its groups were divided at random into `test` and `development`, within strata of inducer, cell type and quiescent controls. The training data are `development` and every study outside set 1. For Agentomics-ML they were divided into train (2,107 samples, 92 laboratory groups) and validation (688 samples, half of the laboratory groups of `development`, 37). `test` holds 687 labelled samples (342 senescent, 345 not senescent) in 78 datasets from 71 laboratory groups. No laboratory group is in two parts. With the links made, checks for shared people, shared sample accessions and copied expression found none across the parts.

## Results

Agentomics and the six published tools on `test`. No sample of these laboratories was in the data that Agentomics was trained or selected on, and none is in the lists of the data that the tools' authors used. Each sample was scored alone: Agentomics computes a sample's probability from that sample's counts only, and each tool was given a count table that held that one sample (hUSI: below). Each measure is given with its 95% interval; values are rounded to three decimals, and `results.tsv` holds them to four, with the counts behind each value.

| Method | Score | AUROC across studies | Pooled AUROC | Pooled AUPRC | Samples |
| --- | --- | --- | --- | --- | --- |
| Agentomics | probability of senescent | 0.834 (0.775 to 0.883) | 0.864 (0.797 to 0.913) | 0.878 (0.819 to 0.923) | 687 |
| hUSI | `hUSI_unscaled` | 0.772 (0.716 to 0.824) | 0.800 (0.728 to 0.862) | 0.822 (0.750 to 0.883) | 687 |
| SenCID | `RecSID_Score` | 0.752 (0.691 to 0.808) | 0.777 (0.708 to 0.837) | 0.795 (0.727 to 0.857) | 687 |
| SenCat | `rna_score` | 0.750 (0.699 to 0.800) | 0.761 (0.697 to 0.821) | 0.772 (0.715 to 0.832) | 687 |
| SENCAN | `SENCAN_score` | 0.721 (0.669 to 0.768) | 0.722 (0.678 to 0.771) | 0.710 (0.628 to 0.784) | 687 |
| markeR | `HernandezSegura:ssGSEA` | 0.663 (0.615 to 0.711) | 0.702 (0.635 to 0.766) | 0.725 (0.636 to 0.799) | 645 |
| SenFinder | `Senescence_probability` | 0.621 (0.589 to 0.656) | 0.636 (0.593 to 0.689) | 0.609 (0.572 to 0.662) | 687 |

Agentomics minus each tool, on the same resamples and on the samples both have a score for:

| Method | AUROC across studies | Pooled AUROC | Pooled AUPRC | Samples |
| --- | --- | --- | --- | --- |
| hUSI | +0.062 (+0.026 to +0.099) | +0.065 (+0.023 to +0.103) | +0.056 (+0.024 to +0.088) | 687 |
| SenCID | +0.082 (+0.040 to +0.124) | +0.087 (+0.045 to +0.125) | +0.084 (+0.044 to +0.119) | 687 |
| SenCat | +0.084 (+0.045 to +0.119) | +0.104 (+0.057 to +0.142) | +0.107 (+0.063 to +0.137) | 687 |
| SENCAN | +0.113 (+0.077 to +0.150) | +0.143 (+0.084 to +0.179) | +0.169 (+0.114 to +0.221) | 687 |
| markeR | +0.169 (+0.123 to +0.215) | +0.160 (+0.110 to +0.201) | +0.141 (+0.085 to +0.207) | 645 |
| SenFinder | +0.213 (+0.168 to +0.255) | +0.228 (+0.168 to +0.275) | +0.269 (+0.215 to +0.306) | 687 |

Every interval of the difference lies above 0. On the 645 samples that markeR scores, Agentomics gives 0.832, 0.862 and 0.866.

- **AUROC across studies.** For two studies from different laboratory groups, the probability that a senescent sample of one study scores higher than a not senescent sample of the other, averaged over all such ordered pairs of studies. It asks whether a score means the same in different laboratories. A random score gives 0.5.
- **Pooled AUROC.** All samples are ranked together by their score; the share of all pairs of a senescent and a not senescent sample in which the senescent sample scores higher (ties count half). A random score gives 0.5.
- **Pooled AUPRC.** In the same single ranking, the average precision with which the senescent samples are found (scikit-learn's average precision, tied scores taken as one threshold). A random score gives about the share of senescent samples: 0.498 on the 687 samples, 0.470 on markeR's 645.
- **95% interval.** The 2.5th and 97.5th percentiles over 2,000 resamples of the laboratory groups, drawn with replacement (seed 20261116). A resample keeps all samples of the groups it draws, so the interval shows how much the result depends on which laboratories are in `test`. Every method uses the same resamples, and the difference between Agentomics and a method is computed in each of them.
- **Samples.** The samples a method has a score for. markeR's wrapper stopped on 42 samples (a score that was not a finite number), so markeR is scored, and compared with Agentomics, on the other 645.

Notes on the methods. Agentomics is scored by its probability of senescent from `agentomics_model/inference.py` with the shipped model file. hUSI is scored by `hUSI_unscaled`, the correlation before hUSI's rescaling, which the wrapper writes and the tool does not release; hUSI's wrapper refuses a single sample, so its scores come from runs on each dataset's labelled samples together, and the unscaled correlation of a sample depends on that sample alone. markeR is scored by the Hernandez-Segura gene set with `ssGSEA`, because markeR's `logmedian` is 0 for every sample scored alone. SenFinder runs with 19 of its 20 genes; the gene without an NCBI gene id is set to 0, as its README prescribes for a missing value. Our scoring also covered 19 published gene sets scored with markeR's `ssGSEA`; senBench reports the models only.

## Scoring a model

`pipeline/score.py` computes the results table from the prediction files of the seven models (`agentomics_model/predictions.tsv`, `tools/<tool>/predictions.tsv`: columns `dataset`, `sample` and `score`, one row per scored `test` sample) and the labels of `test` (`pipeline/test_samples.tsv`), with the measure code and resampling of our scoring. A new model is scored against the others by writing its predictions in the same form and naming it as the reference:

```
python3 pipeline/score.py --add MyModel=build/MyModel.tsv --reference MyModel --out build/results_MyModel.tsv
```

`pipeline/README.md` (Scoring) gives the rules and the format.

## Agentomics

`agentomics_model/` holds the model's code (`train.py`, `inference.py` and the modules they import), the trained model file, the environment, and the model's predictions of `test`. The model compares each gene's count with the counts of reference genes of similar abundance in the same sample, adds products of pairs of selected features, and fits one logistic regression. `agentomics_model/README.md` describes the model, its input and output, how to train it again, the one line changed in the agent's code, and the checks against the run.

## The published tools

`tools/` has one folder per tool (`husi`, `sencid`, `sencat`, `sencan`, `marker`, `senfinder`), each with the benchmark's wrapper, a `setup.sh` that clones the tool at the commit the benchmark used and builds its environment, and the tool's predictions of `test`. `tools/README.md` describes the shared input and output, the install folder, and for each tool its score and what its wrapper changes against the tool's published use.

## Folders

- `pipeline/`: the code that downloads, builds, checks, predicts and scores, its environment, the input tables (`pipeline/data/`), the sha256 of our build (`checksums.tsv`) and the labels of `test`.
- `agentomics_model/`: the Agentomics model and its predictions.
- `tools/`: the six published tools and their predictions.
- `results.tsv`: the results table.

## Limitations

- **Labels.** The sample labels were read by language models from the public records, and checked by further readings by language models. No person checked them.
- **The split was not made blind.** This split was drawn after results on earlier splits of the same datasets had been seen, and about half of the `test` studies were in the part `development` of earlier splits, where the tools' results on them were read. The round is a test bed, not the final benchmark.
- **Few quiescent controls.** `test` holds quiescent controls in 2 studies (7 samples), so it says little about telling senescent cells from quiescent ones.
- **The tools are scored on one sample each, which some were not designed for.** SENCAN and SenCat normalise over the samples they are given, markeR's `logmedian` centres over them, and hUSI rescales over them. On one sample, SENCAN's sample is not adjusted to its reference, SenCat keeps only genes with a count above 380, markeR keeps only genes with a count above 70 and gave no score for 42 samples, and hUSI is scored by its correlation before rescaling, which the tool does not release. Each tool also receives counts made by another pipeline than its authors' (`tools/README.md`).
- **Agentomics was trained on this benchmark's own label definition**, on data from the same kind of public records, while each tool was built on its authors' data and definitions. A reference model shows how much training on these data and labels can give. A default logistic regression on the standardised log2 counts per million of 20,117 genes, trained on the same train split, reaches an AUROC across studies of 0.801 (0.738 to 0.854) on `test`, higher than every published tool. No difference between it and Agentomics was computed.
- **The lists of the tools' training data are incomplete for SenCat**, so a `test` laboratory may have contributed to SenCat's training without being in its list.
- **Our counts of raw reads are inputs of the build.** `pipeline/README.md` gives the method of the 25 tables counted from raw reads (the 24 datasets above and E-MTAB-14464, which the build reads only to leave it out by the rule for genomic DNA); making them again needs the raw reads (1,135 GiB) and the aligner and counter at the stated versions.

## Licence

The code of this repository is under the MIT licence (`LICENSE`). The tools' own code is not part of this repository; `setup.sh` downloads it, and each tool keeps its own licence (`tools/README.md`).
