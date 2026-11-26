# Predictions

One file per model of the results table, with the score of every labelled `test` sample that the model scored. Each file has three columns: `dataset`, `sample` (the names of the benchmark's sample table and of `scoring/test_samples.tsv`) and `score`; a higher score means more senescent. Each score is written with the shortest decimal that gives back the same number. A sample on which a tool's call failed has no row: `markeR.tsv` has 645 rows, every other file 687.

| File | Model | Score | Source |
| --- | --- | --- | --- |
| `Agentomics.tsv` | Agentomics | `probability_1`, the probability of senescent | The Agentomics-ML run's `eval_predictions_test.csv` (sha256 `1de842f0f4c369dfeccb9321d2b255c4a53ebf51c17315fa373ccec86ac01944`), its prediction of `test/input` of the senBench dataset, with each anonymous `id` mapped to its dataset and sample through the dataset's key (`build.sh` writes it as `build/agentomics/senescence_bulk_all980_v10_key/id_key.tsv`). `agentomics_model/inference.py` gives these probabilities again (`agentomics_model/README.md`, Checks). |
| `hUSI.tsv` | hUSI | `hUSI_unscaled` | The hUSI wrapper's runs on the count table of each dataset with its labelled samples (our run `20261124-2242-tool-runs-all980_v9`), because the wrapper refuses a single sample; the unscaled correlation of a sample depends on that sample alone. |
| `SenCID.tsv` | SenCID | `RecSID_Score` | The wrapper called on one sample at a time (our run `20261126-1252-one-sample-tool-runs-all980_v10-test`), each sample's table cut from its dataset's count table of the round (`build/all980_v10/counts/`). |
| `SenCat.tsv` | SenCat | `rna_score` | As SenCID. |
| `SENCAN.tsv` | SENCAN | `SENCAN_score` | As SenCID. |
| `markeR.tsv` | markeR | `HernandezSegura:ssGSEA` | As SenCID; the wrapper stopped on 42 samples, because one of its nine scores was not a finite number. |
| `SenFinder.tsv` | SenFinder | `Senescence_probability` | As SenCID, with `--set-absent-genes-to-zero`. |

The files were written with the reading functions of our scorer, with the arguments of our scoring run, so they hold exactly the scores it used: `scoring/score.py` on these files gives `results.tsv` byte for byte. The wrappers and the commands that make such scores are in `tools/`.
