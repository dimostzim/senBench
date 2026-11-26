# Scoring

`score.py` scores the predictions of methods on the 687 labelled samples of `test`. For each method it computes the AUROC across studies, the pooled AUROC and the pooled AUPRC with their 95% intervals from 2,000 resamples of the laboratory groups (seed 20261116), and the difference between a reference method and each other method on the same resamples, on the samples both have scores for. `results.md` defines the measures.

## Files

| File | Content |
| --- | --- |
| `score.py` | The scorer. |
| `across_studies_measures.py`, `benchmark_bulk_classes.py`, `score_benchmark_bulk.py` | The modules of the measures, the names of the two states, and the AUROC that `across_studies_measures.py` imports. |
| `test_samples.tsv` | Study, dataset, sample, laboratory group and state of every labelled sample of `test`: the file that `build.sh` writes as `build/all980_v10/split/test_samples.tsv` (sha256 in `checksums/all980_v10.tsv`), shipped here because `data/` does not hold it. |
| `tests/` | Tests of the measures and of `score.py`. |

## Environment

`score.py` needs Python 3 with numpy and pandas. Our scoring ran in the environment of `pipeline/environments/count-tables.yml` (Python 3.12.7, numpy 2.5.4, pandas 2.2.3, pytest 9.2.0):

```
conda env create --prefix build/scoring-env --file pipeline/environments/count-tables.yml
```

`build.sh` creates the same environment as `build/environments/count-tables`, which can be used in its place.

## Scoring the shipped predictions

From the repository's folder:

```
build/scoring-env/bin/python scoring/score.py --predictions predictions/*.tsv --reference Agentomics --out build/results.tsv
```

`build/results.tsv` then equals `results.tsv` byte for byte. It took 2 seconds. A method is named after its prediction file (`predictions/SenCID.tsv` is `SenCID`), and the methods are ordered by their AUROC across studies.

## Scoring a new model

1. Write the model's predictions as a tab-separated file with the columns `dataset`, `sample` and `score`, one row per sample of `test_samples.tsv`, a higher score meaning more senescent. A model that scores `test/input` of the senBench dataset gives each sample an anonymous `id`; the dataset's key (`build/agentomics/senescence_bulk_all980_v10_key/id_key.tsv`, columns `id`, `split`, `study`, `dataset`, `sample`, `laboratory_group`, `part`, `label`) maps it to `dataset` and `sample`. The benchmark's rules apply: the score of a sample comes from that sample alone, and the model is trained and selected on the training data only, never on `test`.
2. Score it against the shipped predictions, with the new model as reference:

```
build/scoring-env/bin/python scoring/score.py --predictions predictions/*.tsv build/MyModel.tsv --reference MyModel --out build/results_MyModel.tsv
```

The output has the columns of `results.tsv`, with `reference_minus_method` the new model minus each other method. A method with no score for some samples is scored on the samples it has, and compared with the reference on the samples both have. Every method scored on the same samples uses the same resamples, because the resamples depend only on the laboratory groups. `score.py` stops on a prediction file whose columns differ, that names a sample twice or one that is not in `test_samples.tsv`, or that holds a score that is not a finite number.

## What differs from our scorer

`score.py` is our scorer, `score_single_sample_methods.py` (sha256 `08ff8c862155fc97569702b9cea038a3b007edd1bc6a82e02cb13b7194c39605`), with its reading of the scores replaced:

- **Unchanged:** `read_table`, `compute_average_precision`, `resample_measures`, `describe_random_value`, `score_methods` and `compare_methods`, byte for byte; `count_entered` and the name `AUROC_ACROSS`, copied unchanged from `score_benchmark_bulk_across_studies.py`.
- **Replaced:** our scorer read Agentomics' predictions through the dataset's key and the tools' scores from the folders of their runs, and found the study and laboratory group of each sample through the sample table and the split (`read_model_scores`, `find_main_method`, `read_tool_scores`, `read_one_sample_scores`, `select_test_rows`, `keep_test_set`, `select_test_set_rows`). `score.py` reads prediction files (`read_predictions`) and `test_samples.tsv` (`read_test_samples`, `join_test_labels`) instead.
- **Output:** one table (`tabulate_results`) in place of `measures.tsv` and `differences.tsv`, with the same values to four decimals; its command line takes prediction files and a reference.

`across_studies_measures.py`, `benchmark_bulk_classes.py` and `score_benchmark_bulk.py` are our modules, unchanged (`benchmark_bulk_classes.py` is also in `pipeline/`, the same file). `score_benchmark_bulk.py` is needed only for the `compute_auroc` that `across_studies_measures.py` imports. `tests/test_across_studies_measures.py` is our test file, unchanged; `tests/test_score.py` holds the three tests of the measures of our scorer's test file, unchanged but for the import, tests of the reading of prediction files, and the test that the shipped predictions give `results.tsv`.

## Tests

```
cd scoring/tests
../../build/scoring-env/bin/python -m pytest
```

18 tests pass.
