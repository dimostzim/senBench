# How the benchmark data are made

The benchmark is a table of raw RNA-seq read counts of cultured human cells, in which every sample is labelled `senescent` or `not senescent`. It holds 3,573 labelled samples (1,860 `senescent`, 1,713 `not senescent`) of 318 public datasets from GEO and ArrayExpress, which form 310 studies. This version of the benchmark is our round `all980_v10`. `download.sh` fetches every public file the round is made from, and `build.sh` turns these files and the tables in `data/` into the round's count tables, the benchmark table, the split, the training data and `test`, and the dataset in the layout of Agentomics-ML; it then compares each output with the sha256 of our own build in `checksums/`. The labels are shipped in `data/` and are not made again: language models made them by reading the public records, and `build.sh` makes no model call.

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
| `GEO uniform` | 196 | GEO's raw counts per NCBI gene id, which NCBI computes for human RNA-seq series (`<GSE>_raw_counts_GRCh38.p13_NCBI.tsv.gz`) | The column of each sample is found by its GEO sample accession. In three datasets (GSE76125, GSE156761, GSE125632) GEO attaches reads to other samples than the submitters' own files do; `data/counts/count_column_reassignments.tsv` gives the column that each of their samples takes, from a rule fixed before it was applied (a correlation with the submitters' counts of at least 0.5, as a mutual best match) and from checks by hand of evidence that does not rest on senescence genes, and a sample without a column is left out |
| `authors raw counts` | 83 | The supplementary files of the submitters on GEO or ArrayExpress | A stored layout (`data/counts/layouts/<dataset>.json`) names the file, the gene column and the column of each sample; it was written once by Claude Sonnet 5.5 from a description of the files and is applied as stored. The values must be whole numbers that are not negative, rows of reads not assigned to a gene are left out, and at most 5% of the rows may lack a gene identifier or repeat one. GSE159551 is built from the submitters' table because GEO's reads are exchanged between its states |
| `authors expected counts, rounded` | 15 | Supplementary files of expected counts per gene (as RSEM writes them) | As for raw counts, with each value rounded half to even; a file whose values are scaled or per transcript is refused |
| `counted from raw reads` | 24 | The raw reads at ENA, and for GSE352985 NCBI's SRA files | HISAT2 2.1.0 aligns each run to NCBI's index of the GRCh38 analysis set without alternate loci; a run with 50% or fewer reads aligned is left out; featureCounts (Subread 2.0.6) counts on the exons of NCBI's Annotation Release 109.20190905, the genes GEO counts, with every alignment of a read counted and a share of a read given to each gene it overlaps (`-M -O --fraction`; `-p --countReadPairs` for paired reads); the runs of a sample are summed and rounded half up. Each sample is aligned with the number of threads that our count used (5, 6 or 10, the column `threads` of `data/counts/raw_read_runs.tsv`), because HISAT2 chooses among equally good alignments of a read differently with another number of threads: one sample of GSE250218 aligned with 6 threads in place of 5 gave other counts for 7 of 39,376 genes, by at most 2 |

Every table is mapped to NCBI gene ids with GEO's annotation table `Human.GRCh38.p13.annot.tsv.gz` (39,376 genes), and the counts of rows that map to one gene are added. The round leaves out a dataset whose table maps to fewer than 10,000 NCBI gene ids (6 datasets) and a dataset whose median sample has reads on more than 35,000 genes, as genomic DNA in the libraries gives (E-MTAB-14464, also counted from raw reads). `build.sh` makes the tables of these 325 datasets. The table of count tables also lists 153 labelled datasets whose tables the round never reads (their tables hold no raw counts, or none of their samples is kept); their rows are shipped in `data/counts/records_not_rebuilt/`, so that the reasons in `benchmark_bulk_left_out.tsv` are those of our build.

The column of a sample is found in the table by its accession, by its title, or by the stored layout (`column_match` in `benchmark_bulk_samples.tsv`); the sample table that `build.sh` makes from the GEO records and the ArrayExpress sample sheets gives the titles, the organisms and the assays.

## How the labels were made

The labels are in `data/labels/`: `sample_labels.tsv.gz` (32,462 samples of 1,007 datasets), `dataset_labels.tsv`, and `label_corrections.tsv` (339 rows). They are shipped as our build used them.

1. **Readings of the records.** Two readings by Claude Sonnet 5.5 (effort `high`) of the fields of each dataset that GEO or ArrayExpress holds: the title, summary and overall design of the series, and the fields of every sample. A reading saw no paper and no data. Where the two readings gave another `state` for any sample, Claude Fable 5.2 (effort `high`) read the whole dataset a third time. The label of a sample is the state that two readings agree on, otherwise `unclear`. The columns `state_reading_1` to `state_reading_3` hold the three readings.
2. **Corrections** (`label_corrections.tsv`; the column `source` names the reading). `L.4`: the corrections that a reading of the papers recommended for 13 datasets, in which at least four of six methods ranked the samples against the labels. `M.2`: Claude Fable 5.2 (effort `xhigh`) read the records and papers of 33 datasets (these 13 and the 20 that were compared by hand in the trials of the readings). Their rows have `selection` `label check` or `trial`. `M.4`, `M.7` and `M.8`: Claude Fable 5.2 (effort `xhigh`) read the GEO or ArrayExpress record and the paper of every dataset of `test`, without any score of a method, any earlier check of the labels or any expression value; their changes are rows with `selection` `test reading`.
3. **Rules of use.** `test` takes the round's label with only the `test reading` rows applied, and a `test` dataset left with samples of one state or none leaves `test`. The training data take every correction that names a sample (`label check`, `trial` and `test reading`); a correction to any state other than `senescent` or `not senescent` makes the sample `excluded`.

## The split

1. **Laboratory groups.** Two datasets are in one laboratory group when their GEO records share the contact person and institute, the last contributor, a PubMed id, a BioProject, a SuperSeries, or two contributors, directly or through other datasets. Names are compared as sha256 hashes of the names written without accents, with German umlauts as two letters and with the words of a name in sorted order, so that a name with its first and last name exchanged gives the same hash. `data/split/geo_series_contacts.tsv` holds these hashes and no name. Four links were made by hand, each with its evidence in `data/split/laboratory_group_links.tsv`.
2. **Set 1.** A laboratory group with a dataset in the list of any method (`data/split/method_source_datasets.tsv`: the datasets each method was trained, validated or tested on, or made from, with the role; a dataset that a method was only applied to does not count) is `not in set 1`; a list entry outside GEO is matched by its PubMed id (`data/split/accession_pubmed_ids.tsv`).
3. **`test` and `development`.** The laboratory groups of set 1 are divided once, stratified by the class of the inducer, the class of the cell type and whether the study has quiescent controls, by `numpy.random.default_rng(20261116).permutation` within each stratum; the first part drawn is `test`. The split is made from scratch, with no amendment of parts (`data/split/split_amendments.tsv` is empty). `split/study_split.tsv` gives the laboratory group, stratum, part and set of every study.

## What `download.sh` needs

Measured on our own downloads:

| Folder of the download | Files | Datasets | GiB | What the source states and the downloader checks |
| --- | --- | --- | --- | --- |
| `geo` (series records, supplementary count files) | 1,107 | 980 | 18.24 | nothing |
| `geo_uniform_counts` (and GEO's annotation table) | 197 | 196 | 0.13 | nothing |
| `biostudies_arrayexpress` (IDF, SDRF, processed files) | 70 | 27 | 0.38 | size |
| `reference` (HISAT2 index, annotation, assembly report) | 4 | | 3.95 | size and md5 |
| `raw_reads` (536 runs of 476 samples) | 927 | 25 | 1,135.04 | size and md5 |
| Together | 2,305 | | 1,157.75 | |

At the 14.9 MiB/s that our download of raw reads from ENA reached on 2026-11-18, the whole download takes about 22 hours; a test of `download.sh` on 2026-11-26 reached 2.8 MiB/s (3.08 GiB in 19 minutes), at which it would take about 5 days. `download.sh` downloads 3 files at a time, pauses one second after each file, waits a minute after a refused or cut transfer and five minutes after an answer of HTTP 429, and stops when less than 10 GB of the disk stay free. `download.sh --no-raw-reads` leaves out the raw reads and the reference (1,138.99 GiB; 18.76 GiB remain); `build.sh` then needs the count tables of the 25 datasets counted from raw reads through `--counted-reads`. A file whose size or checksum differs from what the source states is set aside under the name `<file>.checksum_mismatch` and counts as failed; the next run of `download.sh` downloads it again. Where the source states neither (GEO), `download_summary.tsv` tells whether the md5 of the file equals that of the file our build used, so that a record that GEO has changed since our download (2026-11-10 to 2026-11-18) is seen before the build.

## What `build.sh` does and checks

`build.sh` creates the three conda environments of `pipeline/environments/` (pinned to the packages of our build, for Linux x86-64), checks that the data folder holds every file it needs, and makes, in order:

1. `reference/`: the HISAT2 index and the gene features that GEO counts (`make_geo_gene_features.py`), only when reads are to be counted;
2. `raw_read_counts/<dataset>/`: the counts of one dataset's raw reads (`count_raw_reads.py`), and `raw_read_counts/counted/`, one table per dataset;
3. `sample_table/`: the fields of the GEO records and ArrayExpress sample sheets, the sample table and the fields that the label readings read (`label_inputs.jsonl.gz`);
4. `counts/`: the count tables, their records, `count_tables.tsv` and `count_table_columns.tsv.gz`;
5. `all980_v10/`: the round (`benchmark_bulk_samples.tsv`, `benchmark_bulk_datasets.tsv`, `benchmark_bulk_left_out.tsv`, `counts/`, `benchmark_bulk_counts.tsv.gz`), the class tables, and `split/` with the laboratory groups, the split, `training_data_samples.tsv` and `test_samples.tsv`;
6. `agentomics/`: the dataset in the layout of Agentomics-ML and its key.

A step whose folder exists is not run again, so a build that stopped continues with the step that did not finish; to run a step again, delete its folder and the folders of the steps after it. Counting the raw reads takes most of the time: our counts of these datasets took about 40 hours on 16 CPUs (in the check of `build.sh`, GSE250218, 8 samples and 3 GiB of reads, took 6 minutes, and GSE352985, 6 samples and 13 GiB, 29 minutes), and every other step together took 16 minutes on 16 CPUs.

At the end, `checksum_comparison.tsv` gives, for each of the 695 outputs listed in `checksums/all980_v10.tsv`, whether it equals the output of our build (`equal`, `differs`, `missing`, or `not built` for the reference when no reads were counted), and `build.sh` exits with 1 if an output differs or is missing. A gzip file is compared by its uncompressed content, because gzip stores the time of writing. Three outputs are not compared, because they differ between two builds from the same files:

- `all980_v10/split/study_split_record.json` holds the time and the folders of the build;
- `counts/count_tables.tsv` holds, for a dataset counted from raw reads, the sha256 of its counted table as a gzip file;
- `sample_table/samples.tsv.gz` holds in `senescence_evidence` a text that the rules quote from fields joined in the order of a Python set, which changes with Python's hash seed in some rows (51 rows between two seeds).

The round reads none of these values; the tables themselves, their columns, the labels and every output after them are compared.
