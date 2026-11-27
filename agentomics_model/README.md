# Agentomics

Agentomics is the model that Agentomics-ML built for senBench: iteration 19 of the Agentomics-ML run on the senBench dataset (`senescence_bulk_all980_v10`). An agent (gpt-6.2-sol through Codex, Agentomics-ML commit `bfce4b6`) wrote, trained and checked one model per iteration for 20 iterations, from the train and validation splits only; `test` was hidden from it. Iteration 19 had the highest validation AUROC (0.862), and the rule fixed before the run made the iteration with the highest validation AUROC the model to report. This folder holds the agent's code, the trained model file and the environment.

## What the model does

The model gives each sample the probability that it is senescent. It uses the sample's own counts only: no statistic is computed over the other samples of the input, so a sample gets the same probability whatever else is scored with it.

1. **Ordinal features of 8,804 genes.** Each scored gene is compared with 128 reference genes of similar typical abundance in the same sample. Its feature is the share of these reference genes with a lower count than the gene (ties count half), turned into a normal score and standardised with the means and standard deviations of the training samples. Because a gene is compared only with other genes of its own sample, the feature needs no normalisation for sequencing depth. The 4,860 reference genes were chosen on the training samples: genes detected in most samples whose within-sample rank is weakly related to the label. A feature is set to the training mean when the gene is `NA` in the sample or fewer than 32 of its reference genes have a count.
2. **Quadratic terms of 128 selected features.** 128 ordinal features were chosen on the training samples for a difference between senescent and not senescent samples that is the same across groups of samples with different shares of missing and zero counts, leaving out features that correlate above 0.85 with one already chosen. Each of their 8,128 pairs gives one more feature, the product of the two, standardised on the training samples.
3. **One logistic regression** on the 16,932 features (8,804 ordinal and 8,128 quadratic), with one weight per feature and an intercept. It was fitted on the train split only. Besides the log loss, the fit has a ranking loss on pairs of a senescent and a not senescent training sample that resemble each other, a penalty on the size of the weights, and penalties that keep the score from following the share of missing and zero counts, or the main directions of variation of 1,024 genes weakly related to the label among the not senescent training samples.

## Files

| File | Content |
| --- | --- |
| `inference.py` | Scores a table of samples with a trained model file. |
| `train.py` | Fits the representation and the logistic regression on a train split, saves the model file and scores the validation split once. |
| `representation.py`, `pipeline.py`, `matching.py`, `helpers/` | The modules that `train.py` and `inference.py` import. |
| `artifacts/model.joblib` | The trained model of iteration 19: the gene lists, the reference genes, the training means and scales, the pairs and the weights. |
| `prepare_splits.py` | Writes the train and validation splits of the senBench dataset in the form `train.py` reads (added for this repository; see below). |
| `environment.yml` | The packages of the agent's environment, at its versions. |

## Input and output

The input is a folder holding `data.csv`: one row per sample, a column `id`, and one column per gene with the sample's raw read count (whole numbers, not normalised). Genes are named by their symbol in GEO's annotation table `Human.GRCh38.p13.annot.tsv.gz`; a symbol shared by several genes becomes `<symbol>_<GeneID>`, and a gene without a symbol `GeneID_<id>`. A cell is `NA` where the count table of the sample's study has no row for the gene. The model reads 9,804 genes; other columns are ignored, and a gene missing from the header counts as `NA`. The folders `train/input`, `validation/input` and `test/input` of the senBench dataset (`pipeline/README.md`) have this form.

The output is a CSV file with one row per input row, in the input order: `id`, `prediction` (1 for senescent, when the probability is at least 0.5, otherwise 0), `probability_0` (not senescent) and `probability_1` (senescent). The benchmark scores `probability_1`. `agentomics_model/predictions.tsv` holds `probability_1` of every labelled `test` sample, by dataset and sample.

## How to use it

Build the environment (conda or micromamba; 7.2 GB, most of it torch with its CUDA libraries):

```
conda env create --prefix /path/to/agentomics-env --file agentomics_model/environment.yml
```

Score samples with the trained model (`--artifacts-dir` defaults to `agentomics_model/artifacts`):

```
/path/to/agentomics-env/bin/python agentomics_model/inference.py --input /path/to/senescence_bulk_all980_v10/test/input --output predictions.csv
```

Train the model again. `train.py` reads labels in the form that Agentomics-ML gives its agent (`labels.csv` with `id` and `numeric_label`, 1 for senescent), so `prepare_splits.py` writes the train and validation splits in that form first, linking their `input` folders:

```
/path/to/agentomics-env/bin/python agentomics_model/prepare_splits.py --dataset-dir /path/to/senescence_bulk_all980_v10 --out-dir /path/to/prepared
/path/to/agentomics-env/bin/python agentomics_model/train.py --train-data /path/to/prepared/train --validation-data /path/to/prepared/validation --artifacts-dir /path/to/new_artifacts
```

`train.py` writes `model.joblib` into `--artifacts-dir`, replacing a file of that name, so give it a new folder rather than `agentomics_model/artifacts`. It also writes the validation predictions to `agentomics_model/validation_predictions.csv`, next to itself, as it did in the agent's folder. Scoring the 687 `test` samples took 10 to 12 seconds and training 75 seconds, on a CPU with 8 threads or on an NVIDIA A10. When torch finds a CUDA GPU, the code fits and scores on it, as in the agent's run; otherwise it computes on the CPU.

## Changes to the agent's code

The agent's files come from `workspace/best_iteration_snapshot/` of the run's folder (run `20261126-1416-senescence_bulk_all980_v10-gpt-6.2-sol-20`). One line was changed, because it named a folder of the agent's container:

- `inference.py` (the file of `model_inference/`): the default model folder `Path('/workspace/run/current_iteration/model_training/training_artifacts')` became `Path(__file__).resolve().parent / 'artifacts'`.

Every other file is the agent's, byte for byte: `train.py`, `representation.py`, `pipeline.py`, `matching.py` and `helpers/` from `model_training/`, `model.joblib` from `model_training/training_artifacts/`. `model_training/` held a second `inference.py`, which differs from the one of `model_inference/` only in its first line and in its default model folder; `train.py` imports two names from it that are the same in both, so only one is kept. `environment.yml` holds every conda package of the agent's environment (`runtime_info/environment.yml`) and, of its pip packages, torch 2.9.0 and the packages torch needs, all at the agent's versions; 26 pip packages that neither the code nor torch imports (transformers, huggingface-hub and others) are left out. `prepare_splits.py` is new.

## Checks

- **Predictions.** `inference.py` with `artifacts/model.joblib` on `test/input` of the senBench dataset, compared with the run's `eval_predictions_test.csv` (sha256 `1de842f0…`, 687 rows). The run scored on an NVIDIA RTX A4000. On an NVIDIA A10, 679 of the 687 rows are identical; 8 rows of the last batch of 15 rows (the code scores 96 rows at a time) differ in the last digit of a probability (at most 1.1e-16). On the CPU, 473 rows differ in the last digits (at most 3.9e-16). The ids and the class calls are identical in both.
- **Training.** `train.py` on the train and validation splits of the senBench dataset gives the same representation as the model file: the same genes, reference genes, means, scales, pairs and product moments. The weights differ by at most 3.9e-7 on the A10 and 3.0e-7 on the CPU (the largest weight is 0.012), the intercept by 6.6e-5 and 3.2e-5. Two fits on the CPU are identical, so the difference comes from the arithmetic of the hardware, which changes the path of the optimiser. The validation probabilities differ from the agent's by at most 1.7e-5, no validation call changes, and the validation AUROC is 0.862 as in the run.
