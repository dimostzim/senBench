# The data, the build and the scorer

This folder holds the code of the three commands of the repository (`download.sh`, `build.sh`, `run.sh`) and the tables they read (`data/`).

## Terms

- **Dataset:** one GEO series (GSE) or ArrayExpress experiment (E-MTAB).
- **Study:** the datasets that share samples, such as a GEO SuperSeries and its SubSeries. Each sample is kept once.
- **Laboratory group:** the datasets whose records share the contact person and institute, the last contributor, two contributors, a PubMed id, a BioProject or a SuperSeries, directly or through other datasets, or that one of four links made by hand joins.
- **Set 1:** the laboratory groups of which no dataset is among the data that a published method was trained, validated or tested on, or made from.
- **`test` and `development`:** the two parts into which the laboratory groups of set 1 are divided at random, within strata of inducer, cell type and quiescent controls.
- **Training data:** the samples of `development` and of the studies outside set 1.

## Inputs

`download.sh` downloads the 1,374 files of `data/downloads.tsv` into `downloads/` (18.76 GiB): the records and the submitters' count files of the GEO series, GEO's uniform count tables with its gene annotation table, and the files of the ArrayExpress experiments. A file downloaded before is not downloaded again. `downloads/download_summary.tsv` says for every file whether its md5 equals that of the file we built from.

The tables of `data/` are shipped, and `build.sh` does not make them again:

| File | Content |
| --- | --- |
| `sample_labels.tsv.gz` | The label of every sample (`senescent`, `not senescent`, `excluded` or `unclear`), its inducer, kind of control, cell type, cell line, time after induction and other treatment. Two readings of the public records by a language model gave each label, with a third reading where the two differed |
| `dataset_labels.tsv` | Whether each dataset holds cultured cells |
| `label_corrections.tsv` | The labels that later readings of the records and papers changed. `test` takes the rows with `selection` `test reading`; the training data take every row |
| `raw_read_counts/` | Our counts of the raw reads of 25 datasets that have no raw counts on GEO or ArrayExpress, made as NCBI makes GEO's uniform counts: HISAT2 2.1.0 on GRCh38, then featureCounts (Subread 2.0.6) on NCBI's Annotation Release 109.20190905 |
| `rebuilt_datasets.tsv` | The 325 datasets whose count table the build makes |
| `geo_uniform_counts_state.tsv` | The GEO datasets that have a uniform count table, with its address |
| `count_layouts.json` | For 104 datasets with count files of the submitters: the file, the gene column and the column of each sample, as a language model read them from the files |
| `count_column_reassignments.tsv` | The GEO column that each sample of three datasets takes, and the datasets built from the submitters' table instead |
| `count_records_not_rebuilt.json` | The rows of the table of count tables for 153 datasets that have no raw counts or no kept sample |
| `geo_series_contacts.tsv` | The fields of the records that join datasets into laboratory groups, with every name as a sha256 |
| `laboratory_group_links.tsv` | The four links made by hand |
| `method_source_datasets.tsv`, `list_methods.tsv`, `accession_pubmed_ids.tsv` | The datasets that each published method was built or tested on, the methods, and the PubMed ids of listed accessions outside GEO |

## Outputs

`build.sh` uses the environment that `install.sh` makes from `environment.yml` in `build/envs/pipeline`, and makes four folders in `build/`, in this order. A folder that exists is not made again; to make one again, delete it and the folders after it.

| Folder | Content |
| --- | --- |
| `sample_table/` | The fields of the GEO records and the ArrayExpress sample sheets, the organism and assay of every sample (`samples.tsv.gz`), and the fields that the label readings read (`label_inputs.jsonl.gz`) |
| `counts/` | One count table per dataset (`<dataset>.tsv.gz`), from GEO's uniform counts, the submitters' raw or expected counts, or our counts of the raw reads; `count_tables.tsv` and `count_table_columns.tsv.gz` describe them |
| `all980_v10/` | The benchmark: the labelled human samples with raw counts (`benchmark_bulk_samples.tsv`), their count tables by NCBI gene id (`counts/`, `benchmark_bulk_counts.tsv.gz`), and in `split/` the laboratory groups, the split (`study_split.tsv`), `training_data_samples.tsv` and `test_samples.tsv` |
| `agentomics/` | The training data and `test` in the layout of Agentomics-ML (`senescence_bulk_all980_v10/`, with `train/`, `validation/` and `test/`), and the key from each anonymous `id` to its sample (`senescence_bulk_all980_v10_key/id_key.tsv`) |

The benchmark holds 3,573 labelled samples (1,860 senescent, 1,713 not senescent) of 318 datasets. The training data hold 2,795 labelled samples of 129 laboratory groups, and `test` 687 (342 senescent, 345 not senescent) of 71 laboratory groups. Validation is half of the laboratory groups of `development` (688 samples); train is the rest of the training data (2,107 samples).

At the end `build.sh` writes `build/checksum_comparison.tsv`: for each of the 667 outputs listed in `checksums.tsv`, whether it equals the output of our build (`equal`, `differs` or `missing`). It exits with 1 if an output differs or is missing. A gzip file is compared by its uncompressed content.

## Predictions and scores

`install.sh` installs each method into `build/envs/` (the Agentomics model from `agentomics_model/environment.yml`, a tool with its `setup.sh`). `run.sh` runs `predict.py <method>` for each method, and at the end `score.py`.

- `predict.py <method>` scores the samples of `build/all980_v10/split/test_samples.tsv` and writes `predictions/<method>.tsv` (columns `dataset`, `sample`, `score`), replacing the shipped file. For Agentomics it runs `agentomics_model/inference.py` on `build/agentomics/senescence_bulk_all980_v10/test/input`; for a tool, `tools/README.md` describes the call. On another processor, some probabilities of Agentomics differ from the shipped ones after the 15th decimal.
- `score.py` scores every file of `predictions/` on the labels of `test_samples.tsv` and writes `results.tsv`: per method the AUROC within studies, the AUROC across studies, the pooled AUROC and the pooled AUPRC with their 95% intervals from 2,000 resamples of the laboratory groups (seed 20261116), and the Agentomics model minus the method on the same resamples and on the samples both have a score for. It needs numpy and pandas (`build/envs/pipeline/bin/python pipeline/score.py`). On the shipped predictions it gives the shipped `results.tsv`.

To score another model, write its score of every `test` sample to `predictions/<name>.tsv` in the same form and run `score.py`. The score of a sample must come from that sample alone, and the model must not be trained or selected on `test`.
