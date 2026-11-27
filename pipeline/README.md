# How the benchmark data are made

The benchmark is a table of raw RNA-seq read counts of cultured human cells, in which every sample is labelled `senescent` or `not senescent`. It holds 3,573 labelled samples (1,860 `senescent`, 1,713 `not senescent`) of 318 public datasets from GEO and ArrayExpress, which form 310 studies. This version of the benchmark is our round `all980_v10`. `download.sh` fetches every public file the round is made from, and `build.sh` turns these files and the tables in `pipeline/data/` into the round's count tables, the benchmark table, the split, the training data and `test`, and the dataset in the layout of Agentomics-ML; it then compares each output with the sha256 of our own build in `pipeline/checksums.tsv`. The labels and our counts of the raw reads of 25 datasets are shipped in `pipeline/data/` and are not made again: language models made the labels by reading the public records, and `build.sh` makes no model call and reads no raw reads.

## Terms

- **Dataset:** one GEO series (GSE) or ArrayExpress experiment (E-MTAB).
- **Study:** the datasets that share samples, such as a GEO SuperSeries and its SubSeries; a study is named after its dataset with the most samples. Each sample is kept once, in one dataset.
- **Count table:** one raw read count per gene and sample of a dataset, with NCBI gene ids as rows.
- **Laboratory group:** the datasets that come from one laboratory, found from the records of the datasets (below). The split and the resampling of the scores use laboratory groups as units.
- **Set 1:** the laboratory groups of which no dataset is among the data that a scored method was trained, validated or tested on, or made from. Studies outside set 1 are called `not in set 1`.
- **`test` and `development`:** the two parts into which the laboratory groups of set 1 are divided.
- **Training data:** the samples of `development` and of the studies `not in set 1`.

## What the benchmark holds

Every count of this section is of labelled samples (`senescent` or `not senescent`).

| Table (in the build folder) | Studies | Laboratory groups | Datasets | Samples | `senescent` | `not senescent` |
| --- | --- | --- | --- | --- | --- | --- |
| The round, `all980_v10/benchmark_bulk_samples.tsv` | 310 | | 318 | 3,573 | 1,860 | 1,713 |
| Training data, `development` part (`split/training_data_samples.tsv`) | 98 | 74 | 99 | 1,138 | 548 | 590 |
| Training data, `not in set 1` | 125 | 55 | 132 | 1,657 | 922 | 735 |
| Training data, together | 223 | 129 | 231 | 2,795 | 1,470 | 1,325 |
| `test` (`split/test_samples.tsv`) | 78 | 71 | 78 | 687 | 342 | 345 |

The training data use the column `state_for_training`, the round's label with every recorded correction applied; 33 of the 2,828 samples of the training data are `excluded` by a correction. `test` has one label column, `state`. The two tables share no sample, study or laboratory group.

The dataset in the layout of Agentomics-ML (`agentomics/senescence_bulk_all980_v10/`) holds the same samples with an anonymous `id`, the raw counts of all 39,376 genes of GEO's annotation table (`NA` where a dataset's table has no row for a gene) and `label`. Validation is half of the laboratory groups of `development` that hold a labelled sample, drawn with seed 20261124; train is the rest of the training data. The key from `id` to sample is written apart, in `agentomics/senescence_bulk_all980_v10_key/id_key.tsv`.

| Split | Samples | `senescent` | `not_senescent` | Studies | Laboratory groups | `NA` cells |
| --- | --- | --- | --- | --- | --- | --- |
| train | 2,107 | 1,150 | 957 | 171 | 92 | 14.2% |
| validation | 688 | 320 | 368 | 52 | 37 | 14.9% |
| test | 687 | 342 | 345 | 78 | 71 | 19.6% |

## Where each kind of count comes from

Each dataset of the round has one count table, of one of four kinds (`count_table_state` in `benchmark_bulk_datasets.tsv`):

| Kind | Datasets in the round | Source | How the table is made |
| --- | --- | --- | --- |
| `GEO uniform` | 196 | GEO's raw counts per NCBI gene id, which NCBI computes for human RNA-seq series (`<GSE>_raw_counts_GRCh38.p13_NCBI.tsv.gz`) | The column of each sample is found by its GEO sample accession. In three datasets (GSE76125, GSE156761, GSE125632) GEO attaches reads to other samples than the submitters' own files do; `data/count_column_reassignments.tsv` gives the column that each of their samples takes, from a rule fixed before it was applied (a correlation with the submitters' counts of at least 0.5, as a mutual best match) and from checks by hand of evidence that does not rest on senescence genes, and a sample without a column is left out |
| `authors raw counts` | 83 | The supplementary files of the submitters on GEO or ArrayExpress | A stored layout (the dataset's entry in `data/count_layouts.json`) names the file, the gene column and the column of each sample; it was written once by Claude Sonnet 5.5 from a description of the files and is applied as stored. The values must be whole numbers that are not negative, rows of reads not assigned to a gene are left out, and at most 5% of the rows may lack a gene identifier or repeat one. GSE159551 is built from the submitters' table because GEO's reads are exchanged between its states |
| `authors expected counts, rounded` | 15 | Supplementary files of expected counts per gene (as RSEM writes them) | As for raw counts, with each value rounded half to even; a file whose values are scaled or per transcript is refused |
| `counted from raw reads` | 24 | Our counts of the raw reads (below), shipped as `data/raw_read_counts/<dataset>.tsv.gz` | Taken as they are |

Every table is mapped to NCBI gene ids with GEO's annotation table `Human.GRCh38.p13.annot.tsv.gz` (39,376 genes), and the counts of rows that map to one gene are added. The round leaves out a dataset whose table maps to fewer than 10,000 NCBI gene ids (6 datasets) and a dataset whose median sample has reads on more than 35,000 genes, as genomic DNA in the libraries gives (E-MTAB-14464, also counted from raw reads). `build.sh` makes the tables of these 325 datasets (`data/rebuilt_datasets.tsv`). The table of count tables also lists 153 labelled datasets whose tables the round never reads (their tables hold no raw counts, or none of their samples is kept); their records are shipped in `data/count_records_not_rebuilt.json`, so that the reasons in `benchmark_bulk_left_out.tsv` are those of our build.

The column of a sample is found in the table by its accession, by its title, or by the stored layout (`column_match` in `benchmark_bulk_samples.tsv`); the sample table that `build.sh` makes from the GEO records and the ArrayExpress sample sheets gives the titles, the organisms and the assays.

## How we counted the raw reads of 25 datasets

For 25 datasets that have no raw counts on GEO or ArrayExpress (24 of the round, and E-MTAB-14464, which the round leaves out by the rule for genomic DNA), we counted the raw reads ourselves in the way NCBI counts GEO's uniform tables, from 2026-11-18 to 2026-11-24. The 25 tables are shipped in `data/raw_read_counts/`, 20.2 MB in all; the reads themselves (1,135 GiB, 927 files of 536 runs of 476 samples) are not part of the download.

1. **Reads.** The FASTQ files that ENA serves for each run, with ENA's md5 checked after the download; for a paired run the two files of the pair. For the 6 runs of GSE352985, of which ENA had no FASTQ files yet, the SRA files of NCBI, written as FASTQ files by `fasterq-dump` (sra-tools 3.2.1) with `--split-3`. A GEO sample's runs are the runs of the SRA experiment that its record names.
2. **Alignment.** HISAT2 2.1.0 with its default options, against NCBI's HISAT2 index of the GRCh38 analysis set without alternate loci (`GCA_000001405.15_GRCh38_no_alt_analysis_set`, md5 `de3609995cd7498322bf2f7b9c57d302`). A run with 50% or fewer of its reads aligned is left out, as NCBI leaves it out, and a sample without a run left has no column: 2 samples of GSE146585 and 9 of GSE290151, whose libraries hold short inserts. The tables therefore hold 465 of the 476 samples.
3. **Counting.** featureCounts (Subread 2.0.6) on the exons of NCBI's Annotation Release 109.20190905, one feature per NCBI gene id: the genes GEO counts, with lengths equal to those of GEO's annotation table for all 39,376 genes. A read that aligns to several places, or that overlaps several genes, is counted at each of them with an equal fraction (`-M -O --fraction`; for paired reads also `-p --countReadPairs`, so that a pair counts once). The counts of the runs of a sample are summed and then rounded half up.
4. **Threads.** HISAT2 chooses among equally good alignments of a read differently with another number of threads. One sample of GSE250218 aligned twice with 5 threads gave our counts exactly; with 6 threads 7 of its 39,376 genes differed (by at most 2 reads), with 1 thread 18 genes (by at most 5). featureCounts gives the same counts on the same alignments with 1, 5 or 6 threads. Our tables were counted with 5 threads per sample (GSE146585, GSE216964, GSE233531, GSE233770, GSE250218, GSE252219, GSE288138, GSE290151, GSE291046, GSE301995, GSE302820, GSE336332, GSE345082), 10 threads (GSE352985, E-MTAB-15005) or 6 threads (the other 10 datasets). Counted again with these thread counts in an earlier version of `build.sh`, GSE250218 (8 samples, 3 GiB of reads) and GSE352985 (6 samples, 13 GiB) gave our tables exactly.
5. **Check against GEO.** The method was chosen on two datasets that have a GEO uniform table (GSE160442, 6 samples of paired reads; GSE155903, 8 samples of single reads with 2 runs each): HISAT2 2.1.0, the featureCounts options and rounding half up after summing agreed best with GEO's counts among the versions and options tried, and were fixed before any of the 25 datasets was counted. With them, the correlation of log counts with GEO's is 0.9995 to 0.9997, the total counts are 100.06% to 100.15% of GEO's, and 97.0% to 98.6% of the genes are within 5% of GEO's count (78.4% to 94.0% identical). The counts are therefore close to GEO's, not identical.

| Dataset | Samples counted | Runs | GiB of reads | Threads | Dataset | Samples counted | Runs | GiB of reads | Threads |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| E-MTAB-10415 | 20 | 20 | 60.09 | 6 | GSE233770 | 4 | 4 | 11.02 | 5 |
| E-MTAB-13786 | 32 | 32 | 184.76 | 6 | GSE250218 | 8 | 8 | 3.02 | 5 |
| E-MTAB-14201 | 65 | 65 | 142.11 | 6 | GSE252219 | 4 | 4 | 6.50 | 5 |
| E-MTAB-14464 | 18 | 18 | 66.62 | 6 | GSE288138 | 20 | 80 | 9.03 | 5 |
| E-MTAB-15005 | 24 | 24 | 117.42 | 10 | GSE290151 | 5 of 14 | 14 | 10.98 | 5 |
| E-MTAB-16030 | 12 | 12 | 18.07 | 6 | GSE291046 | 12 | 12 | 16.59 | 5 |
| E-MTAB-17783 | 6 | 6 | 10.06 | 6 | GSE301995 | 9 | 9 | 20.36 | 5 |
| E-MTAB-3101 | 9 | 9 | 17.20 | 6 | GSE302820 | 8 | 8 | 6.99 | 5 |
| E-MTAB-5403 | 78 | 78 | 76.29 | 6 | GSE328437 | 71 | 71 | 231.06 | 6 |
| E-MTAB-8430 | 20 | 20 | 59.95 | 6 | GSE336332 | 6 | 6 | 15.96 | 5 |
| GSE146585 | 10 of 12 | 12 | 10.21 | 5 | GSE345082 | 6 | 6 | 8.77 | 5 |
| GSE216964 | 8 | 8 | 7.88 | 5 | GSE352985 | 6 | 6 | 13.11 | 10 |
| GSE233531 | 4 | 4 | 11.00 | 5 | | | | | |

The table of count tables (`counts/count_tables.tsv`, compared by `build.sh`) holds the sha256 of each of these 25 files, so the build checks that the shipped files are ours.

## How the labels were made

The labels are in `data/`: `sample_labels.tsv.gz` (32,462 samples of 1,007 datasets), `dataset_labels.tsv`, and `label_corrections.tsv` (339 rows). They are shipped as our build used them.

1. **Readings of the records.** Two readings by Claude Sonnet 5.5 (effort `high`) of the fields of each dataset that GEO or ArrayExpress holds: the title, summary and overall design of the series, and the fields of every sample. A reading saw no paper and no data. Where the two readings gave another `state` for any sample, Claude Fable 5.2 (effort `high`) read the whole dataset a third time. The label of a sample is the state that two readings agree on, otherwise `unclear`. The columns `state_reading_1` to `state_reading_3` hold the three readings.
2. **Corrections** (`label_corrections.tsv`; the column `source` names the reading). `L.4`: the corrections that a reading of the papers recommended for 13 datasets, in which at least four of six methods ranked the samples against the labels. `M.2`: Claude Fable 5.2 (effort `xhigh`) read the records and papers of 33 datasets (these 13 and the 20 that were compared by hand in the trials of the readings). Their rows have `selection` `label check` or `trial`. `M.4`, `M.7` and `M.8`: Claude Fable 5.2 (effort `xhigh`) read the GEO or ArrayExpress record and the paper of every dataset of `test`, without any score of a method, any earlier check of the labels or any expression value; their changes are rows with `selection` `test reading`.
3. **Rules of use.** `test` takes the round's label with only the `test reading` rows applied, and a `test` dataset left with samples of one state or none leaves `test`. The training data take every correction that names a sample (`label check`, `trial` and `test reading`); a correction to any state other than `senescent` or `not senescent` makes the sample `excluded`.

## The split

1. **Laboratory groups.** Two datasets are in one laboratory group when their GEO records share the contact person and institute, the last contributor, a PubMed id, a BioProject, a SuperSeries, or two contributors, directly or through other datasets. Names are compared as sha256 hashes of the names written without accents, with German umlauts as two letters and with the words of a name in sorted order, so that a name with its first and last name exchanged gives the same hash. `data/geo_series_contacts.tsv` holds these hashes and no name. Four links were made by hand, each with its evidence in `data/laboratory_group_links.tsv`.
2. **Set 1.** A laboratory group with a dataset in the list of any method (`data/method_source_datasets.tsv`: the datasets each method was trained, validated or tested on, or made from, with the role; a dataset that a method was only applied to does not count) is `not in set 1`; a list entry outside GEO is matched by its PubMed id (`data/accession_pubmed_ids.tsv`).
3. **`test` and `development`.** The laboratory groups of set 1 are divided once, stratified by the class of the inducer, the class of the cell type and whether the study has quiescent controls, by `numpy.random.default_rng(20261116).permutation` within each stratum; the first part drawn is `test`. The split is made from scratch, with no amendment of parts (`data/split_amendments.tsv` is empty). `split/study_split.tsv` gives the laboratory group, stratum, part and set of every study.

## What `download.sh` downloads

`download.sh` downloads the files listed in `data/downloads.tsv` into `downloads/`, or into the folder given as its argument. Measured on our own downloads:

| Folder of the download | Files | Datasets | GiB | What the source states and the downloader checks |
| --- | --- | --- | --- | --- |
| `geo` (series records, supplementary count files) | 1,107 | 980 | 18.24 | nothing |
| `geo_uniform_counts` (and GEO's annotation table) | 197 | 196 | 0.13 | nothing |
| `biostudies_arrayexpress` (IDF, SDRF, processed files) | 70 | 27 | 0.38 | size |
| Together | 1,374 | | 18.76 | |

A test of an earlier version of `download.sh` on 2026-11-26 reached 2.8 MiB/s, at which the download takes about 2 hours. `download.sh` downloads 3 files at a time, pauses one second after each file, waits a minute after a refused or cut transfer and five minutes after an answer of HTTP 429, gives a file up after 8 attempts in a row without new bytes, and stops when less than 10 GB of the disk stay free. A file that holds an error document of the server is set aside under the name `<file>.error_document` and counts as failed; the next run of `download.sh` downloads it again. GEO states neither the size nor a checksum of its files, so `download_summary.tsv` tells whether the md5 of each file equals that of the file our build used, and a record that GEO has changed since our download (2026-11-10 to 2026-11-18) is seen before the build.

## What `build.sh` does and checks

`build.sh` creates the conda environment of `pipeline/environment.yml` (Python 3.12.7, numpy 2.5.4, pandas 2.2.3, openpyxl 3.1.5 and xlrd 2.0.1, every package pinned to the build of our environment, for Linux x86-64) in `build/environments/count-tables`, checks that the data folder holds every listed file, and makes, in order:

1. `sample_table/`: the fields of the GEO records and ArrayExpress sample sheets, the sample table and the fields that the label readings read (`label_inputs.jsonl.gz`);
2. `counts/`: the count tables, `count_tables.tsv` and `count_table_columns.tsv.gz`;
3. `all980_v10/`: the round (`benchmark_bulk_samples.tsv`, `benchmark_bulk_datasets.tsv`, `benchmark_bulk_left_out.tsv`, `counts/`, `benchmark_bulk_counts.tsv.gz`), the class tables, and `split/` with the laboratory groups, the split, `training_data_samples.tsv` and `test_samples.tsv`;
4. `agentomics/`: the dataset in the layout of Agentomics-ML and its key.

A step whose folder exists is not run again, so a build that stopped continues with the step that did not finish; to run a step again, delete its folder and the folders of the steps after it. On 16 CPUs the build took 16 minutes.

At the end, `checksum_comparison.tsv` gives, for each of the 670 outputs listed in `pipeline/checksums.tsv`, whether it equals the output of our build (`equal`, `differs` or `missing`), and `build.sh` exits with 1 if an output differs or is missing. A gzip file is compared by its uncompressed content, because gzip stores the time of writing. Two outputs are not compared, because they differ between two builds from the same files:

- `all980_v10/split/study_split_record.json` holds the time and the folders of the build;
- `sample_table/samples.tsv.gz` holds in `senescence_evidence` a text that the rules quote from fields joined in the order of a Python set, which changes with Python's hash seed in some rows (51 rows between two seeds).

The round reads neither of these values; the tables themselves, their columns, the labels and every output after them are compared.

## Scoring

`score.py` scores the predictions of methods on the 687 labelled samples of `test`. For each method it computes the AUROC across studies, the pooled AUROC and the pooled AUPRC with their 95% intervals from 2,000 resamples of the laboratory groups (seed 20261116), and the difference between a reference method and each other method on the same resamples, on the samples both have scores for. The README of the repository defines the measures. It needs Python 3 with numpy and pandas; `build/environments/count-tables/bin/python`, which `build.sh` creates, has them.

Without options, `score.py` reads the predictions of the seven models (`agentomics_model/predictions.tsv` and `tools/<tool>/predictions.tsv`), takes Agentomics as the reference, and writes `results.tsv` in the repository's folder, replacing it; on the shipped predictions it gives the shipped `results.tsv` byte for byte (in 2 seconds, with pandas 2.2.3, 2.3.1 and 3.0.3). `test_samples.tsv` holds the study, dataset, sample, laboratory group and state of every labelled sample of `test`: the file that `build.sh` writes as `build/all980_v10/split/test_samples.tsv`.

To score a new model, write its predictions as a tab-separated file with the columns `dataset`, `sample` and `score`, one row per sample of `test_samples.tsv`, a higher score meaning more senescent, and score it with the seven models, as the reference:

```
python3 pipeline/score.py --add MyModel=build/MyModel.tsv --reference MyModel --out build/results_MyModel.tsv
```

A model that scores `test/input` of the senBench dataset gives each sample an anonymous `id`; the dataset's key (`build/agentomics/senescence_bulk_all980_v10_key/id_key.tsv`, columns `id`, `split`, `study`, `dataset`, `sample`, `laboratory_group`, `part`, `label`) maps it to `dataset` and `sample`. The benchmark's rules apply: the score of a sample comes from that sample alone, and the model is trained and selected on the training data only, never on `test`. A method with no score for some samples is scored on the samples it has, and compared with the reference on the samples both have. Every method scored on the same samples uses the same resamples, because the resamples depend only on the laboratory groups. `score.py` stops on a prediction file whose columns differ, that names a sample twice or one that is not in `test_samples.tsv`, or that holds a score that is not a finite number.

`score.py` holds our scorer's measure code unchanged (`compute_average_precision`, `resample_measures`, `describe_random_value`, `score_methods`, `compare_methods`, and the functions of the measures across studies that they call); only the reading of the scores differs. Our scorer read Agentomics' predictions through the dataset's key and the tools' scores from the folders of their runs, and found the study and laboratory group of each sample through the sample table and the split; `score.py` reads prediction files and `test_samples.tsv` instead, and writes one table in place of two, with the same values to four decimals.

## Predictions

`predict.py <method>` scores the `test` samples of a build (`build/`, or the folder given with `--build-dir`) with one of the seven models and writes its `predictions.tsv`, replacing the shipped file. For Agentomics (`python3 pipeline/predict.py agentomics`), it creates the model's environment from `agentomics_model/environment.yml` in `build/environments/agentomics` when it does not exist yet (with conda, or the program given with `--conda`; 7.2 GB, most of it torch), runs `agentomics_model/inference.py` on `agentomics/senescence_bulk_all980_v10/test/input`, and maps each anonymous id to its dataset and sample through the dataset's key; the score is the probability of senescent. On another processor than the run's, some probabilities differ from the shipped ones in their last digits, by at most 4e-16 (`agentomics_model/README.md`, Checks). For a tool, `tools/README.md` describes the call. The work files are in `build/predict/<method>/`, which must not exist yet. `predict.py` needs only Python 3.

## The files of this folder

| File | What it does or holds |
| --- | --- |
| `downloads.py` | Downloads the files of `data/downloads.tsv` (`download`, which `download.sh` runs) and checks a data folder (`check`) |
| `build_sample_table.py` | Copies the fields of the GEO records and the ArrayExpress sample sheets into tables, builds the sample table with its rules, and writes the label inputs |
| `build_count_tables.py` | Builds and checks the count table of every dataset of `data/rebuilt_datasets.tsv`, and writes `count_tables.tsv` and `count_table_columns.tsv.gz` |
| `build_benchmark_bulk.py` | Builds the round: the labelled samples with raw counts, one dataset per sample, the studies, the count tables by NCBI gene id |
| `build_study_split.py` | Builds the laboratory groups, set 1 and the split, and the tables of the classes of cell type and time |
| `write_final_sets.py` | Writes the training data and `test` with the label corrections |
| `build_agentomics_dataset.py` | Writes the dataset in the layout of Agentomics-ML and its key |
| `compare_checksums.py` | Compares the outputs of a build with `checksums.tsv` |
| `score.py`, `test_samples.tsv` | The scorer and the labels of `test` |
| `predict.py` | Scores `test` with Agentomics or a tool |
| `environment.yml` | The environment of the build and the scorer |
| `checksums.tsv` | 670 rows: path in the build folder, `file` or `uncompressed content`, sha256 of our output |
| `data/downloads.tsv` | The 1,374 files to download: folder, accession, file, address, the size the source states, and the size and md5 of our file |
| `data/sample_labels.tsv.gz`, `data/dataset_labels.tsv`, `data/label_corrections.tsv` | The labels (above) |
| `data/inventory.tsv` | Repository, relevance, title, organism, assay and number of downloaded files of the 1,007 labelled datasets, for the sample table |
| `data/rebuilt_datasets.tsv` | The 325 datasets whose count table `build.sh` makes, with the folder of the download that holds each |
| `data/geo_uniform_counts_state.tsv` | Whether GEO has a uniform table for each of the 313 GEO datasets, and its address |
| `data/count_column_reassignments.tsv` | The GEO columns of the samples of three datasets, and the datasets built from the submitters' table instead |
| `data/count_layouts.json` | The stored layout of the submitters' files of 104 datasets: the answer of the model, the sha256 of the request, the model and the effort |
| `data/count_records_not_rebuilt.json` | The records of the 153 datasets whose tables the round never reads |
| `data/raw_read_counts/` | Our count tables of the raw reads of 25 datasets (above) |
| `data/geo_series_contacts.tsv` | The fields of the GEO and ArrayExpress records that link datasets, with every name as a sha256 |
| `data/method_source_datasets.tsv`, `data/main_scores.tsv`, `data/accession_pubmed_ids.tsv` | The datasets of the methods' lists with their roles, the methods, and the PubMed ids of listed accessions outside GEO |
| `data/laboratory_group_links.tsv`, `data/split_amendments.tsv` | The four links made by hand, and the amendments of parts (none) |
