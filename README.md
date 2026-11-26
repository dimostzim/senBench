# senBench

senBench is a benchmark of methods that call a bulk RNA-seq sample of cultured human cells senescent or not senescent, from that sample alone. It tests them on studies from laboratories that no method was trained on. It holds public RNA-seq datasets with a label for every sample, a split of the studies by laboratory, and the scores of seven models on the held-out part `test`: Agentomics, a model that the agent system Agentomics-ML built on senBench's training data, and six published tools.

## The task

A method receives the raw read counts of one sample and returns a score that is higher for senescent samples. A sample is `senescent` when the study made the cells senescent, and `not senescent` when it is a control of the study: proliferating, quiescent, early passage, untreated, vehicle or empty vector cells. The score must be computed from the sample alone, so that it does not depend on which other samples a study holds. Each tool was therefore given a count table with one sample, except hUSI, whose score is a correlation computed on each sample alone but whose wrapper needs at least two samples (`tools/husi/README.md`).

## The data

The benchmark holds 3,573 labelled samples (1,860 senescent, 1,713 not senescent) of 318 datasets from GEO and ArrayExpress. Their counts are GEO's uniformly processed counts, the authors' own count tables, or, for 24 datasets, counts made by us from the raw reads. Language models read the labels from the public records of each dataset; the labels are shipped in `data/` and not made again. `docs/data.md` describes how the data are made, the label rules, the split and the sizes.

`download.sh` downloads every public file the benchmark needs into a folder (default `downloads/`), 1,158 GiB in all, almost all of it raw reads; with `--no-raw-reads` it downloads 18.76 GiB, and `build.sh` then takes the counts of the datasets counted from raw reads through `--counted-reads`. `build.sh` builds from these files the count tables, the sample table, the split, the final sets and the dataset in Agentomics-ML's layout (default `build/`), and compares each output with the sha256 of our own build in `checksums/`:

```
bash download.sh downloads
bash build.sh --threads 16 downloads build
```

## The split

Studies are grouped into laboratory groups: two studies are in one group when their records share the contact person and institute, the last contributor, a PubMed id, a BioProject, a SuperSeries, or two contributors. Four more links were made by hand from the records, three of them found by a check for people shared between the parts. Set 1 holds the laboratory groups with no dataset in the list of data that any published method was trained, validated or tested on, or made from. Its groups were divided at random into `test` and `development`, within strata of inducer, cell type and quiescent controls. The training data are `development` and every study outside set 1. For Agentomics-ML they were divided into train (2,107 samples, 92 laboratory groups) and validation (688 samples, half of the laboratory groups of `development`, 37). `test` holds 687 labelled samples (342 senescent, 345 not senescent) in 78 datasets from 71 laboratory groups. No laboratory group is in two parts. With the links made, checks for shared people, shared sample accessions and copied expression found none across the parts.

## Results

Each measure is given with its 95% interval from 2,000 resamples of the laboratory groups. The AUROC across studies is the probability that a senescent sample of one study scores higher than a not senescent sample of a study from another laboratory group, averaged over all such pairs of studies. The pooled AUROC and the pooled AUPRC rank all samples together. `results.md` defines each measure and gives the differences between Agentomics and each tool; `results.tsv` holds the values with their counts, and `predictions/` the scores behind them.

| Method | Score | AUROC across studies | Pooled AUROC | Pooled AUPRC | Samples |
| --- | --- | --- | --- | --- | --- |
| Agentomics | probability of senescent | 0.834 (0.775 to 0.883) | 0.864 (0.797 to 0.913) | 0.878 (0.819 to 0.923) | 687 |
| hUSI | `hUSI_unscaled` | 0.772 (0.716 to 0.824) | 0.800 (0.728 to 0.862) | 0.822 (0.750 to 0.883) | 687 |
| SenCID | `RecSID_Score` | 0.752 (0.691 to 0.808) | 0.777 (0.708 to 0.837) | 0.795 (0.727 to 0.857) | 687 |
| SenCat | `rna_score` | 0.750 (0.699 to 0.800) | 0.761 (0.697 to 0.821) | 0.772 (0.715 to 0.832) | 687 |
| SENCAN | `SENCAN_score` | 0.721 (0.669 to 0.768) | 0.722 (0.678 to 0.771) | 0.710 (0.628 to 0.784) | 687 |
| markeR | `HernandezSegura:ssGSEA` | 0.663 (0.615 to 0.711) | 0.702 (0.635 to 0.766) | 0.725 (0.636 to 0.799) | 645 |
| SenFinder | `Senescence_probability` | 0.621 (0.589 to 0.656) | 0.636 (0.593 to 0.689) | 0.609 (0.572 to 0.662) | 687 |

A random score gives 0.5 for both AUROCs and 0.498 for the pooled AUPRC. For each tool and each measure, the interval of the difference between Agentomics and the tool lies above 0.

## Scoring a model

`predictions/` holds the score of every labelled `test` sample by each of the seven models, and `scoring/score.py` computes the results table from them, with the same measure code and resampling as our scoring; on the shipped predictions it gives `results.tsv` byte for byte:

```
python scoring/score.py --predictions predictions/*.tsv --reference Agentomics --out build/results.tsv
```

A new model is scored against the others by adding its prediction file (columns `dataset`, `sample`, `score`) and naming it as the reference. `scoring/README.md` describes the format, the environment and the tests.

## Agentomics

`agentomics_model/` holds the model's code (`train.py`, `inference.py` and the modules they import), the trained model file and the environment. The model compares each gene's count with the counts of reference genes of similar abundance in the same sample, adds products of pairs of selected features, and fits one logistic regression. To score the samples of the dataset's `test` part:

```
conda env create --prefix build/agentomics-env --file agentomics_model/environment.yml
build/agentomics-env/bin/python agentomics_model/inference.py --input build/agentomics/senescence_bulk_all980_v10/test/input --output build/agentomics_test_predictions.csv
```

`agentomics_model/README.md` describes the model, its input and output, how to train it again, the one line changed in the agent's code, and the checks against the run.

## The published tools

`tools/` has one folder per tool (`husi`, `sencid`, `sencat`, `sencan`, `marker`, `senfinder`), each with the benchmark's wrapper, a `setup.sh` that clones the tool at the commit the benchmark used and builds its environment, and a README. To score the first sample of a dataset with SenCID:

```
bash tools/sencid/setup.sh build/tools_install
zcat build/all980_v10/counts/GSE130727.tsv.gz | cut -f1,2 | gzip > build/one_sample.tsv.gz
bash build/tools_install/wrappers/sencid/run_tool.sh --counts build/one_sample.tsv.gz --annotation downloads/geo_uniform_counts/annotation/Human.GRCh38.p13.annot.tsv.gz --out-dir build/sencid_scores
```

`tools/README.md` describes the shared input and output and the install folder; each tool's README gives its command, its score and what its wrapper changes against the tool's published use.

## Limitations

- **Labels.** The sample labels were read by language models from the public records, and checked by further readings by language models. No person checked them.
- **The split was not made blind.** This split was drawn after results on earlier splits of the same datasets had been seen, and about half of the `test` studies were in the part `development` of earlier splits, where the tools' results on them were read. The round is a test bed, not the final benchmark.
- **Few quiescent controls.** `test` holds quiescent controls in 2 studies (7 samples), so it says little about telling senescent cells from quiescent ones.
- **The tools are scored on one sample each, which some were not designed for.** SENCAN and SenCat normalise over the samples they are given, markeR's `logmedian` centres over them, and hUSI rescales over them. On one sample, SENCAN's sample is not adjusted to its reference, SenCat keeps only genes with a count above 380, markeR keeps only genes with a count above 70 and gave no score for 42 samples, and hUSI is scored by its correlation before rescaling, which the tool does not release. Each tool also receives counts made by another pipeline than its authors' (`tools/*/README.md`).
- **Agentomics was trained on this benchmark's own label definition**, on data from the same kind of public records, while each tool was built on its authors' data and definitions. A reference model shows how much training on these data and labels can give. A default logistic regression on the standardised log2 counts per million of 20,117 genes, trained on the same train split, reaches an AUROC across studies of 0.801 (0.738 to 0.854) on `test`, higher than every published tool. No difference between it and Agentomics was computed.
- **The lists of the tools' training data are incomplete for SenCat**, so a `test` laboratory may have contributed to SenCat's training without being in its list.

## Licence

The code of this repository is under the MIT licence (`LICENSE`). The tools' own code is not part of this repository; `setup.sh` downloads it, and each tool keeps its own licence (`tools/README.md`).
