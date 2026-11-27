# Agentomics

Agentomics is the model that an agent of Agentomics-ML wrote and trained on the train and validation splits of senBench; `test` was hidden from the agent. The model gives a sample the probability that it is senescent, from the sample's own counts: each gene's count is compared with the counts of reference genes of similar abundance in the same sample, products of pairs of 128 selected features are added, and one logistic regression is fitted on all features.

## Run

`bash install.sh Agentomics` builds the environment of `environment.yml` in `build/envs/Agentomics` (7.2 GB, most of it torch), and `bash run.sh Agentomics` scores `test` and writes `predictions/Agentomics.tsv`. To score other samples:

```
build/envs/Agentomics/bin/python agentomics_model/inference.py --input <input folder> --output <predictions.csv>
```

- **Input:** a folder with `data.csv`: one row per sample, the column `id`, and one column of raw read counts per gene. Genes are named as in `build/agentomics/senescence_bulk_all980_v10/test/input/data.csv`; a cell is `NA` where the count is not known.
- **Output:** one row per sample with `id`, `prediction` (1 for senescent), `probability_0` and `probability_1`. senBench scores `probability_1`.

The code computes on a CUDA GPU when torch finds one, otherwise on the CPU; some probabilities then differ after the 15th decimal.

## Train again

`train.py` reads labels as Agentomics-ML gives them to its agent, so `prepare_splits.py` writes the train and validation splits in that form first:

```
build/envs/Agentomics/bin/python agentomics_model/prepare_splits.py --dataset-dir build/agentomics/senescence_bulk_all980_v10 --out-dir build/agentomics_prepared
build/envs/Agentomics/bin/python agentomics_model/train.py --train-data build/agentomics_prepared/train --validation-data build/agentomics_prepared/validation --artifacts-dir build/agentomics_artifacts
```

`train.py` writes `model.joblib` into `--artifacts-dir` and the validation predictions to `agentomics_model/validation_predictions.csv`. The weights of a new fit differ from those of `artifacts/model.joblib` by at most 4e-7, because the arithmetic of the hardware changes the path of the optimiser.

## What changed in the agent's code

`train.py`, `inference.py`, `representation.py`, `pipeline.py`, `matching.py`, `training_reporter.py` and `artifacts/model.joblib` are the agent's files, except for two lines:

- `inference.py`: the default model folder, which named a folder of the agent's container, became `Path(__file__).resolve().parent / 'artifacts'`.
- `train.py`: `from helpers.training_reporter import TrainingReporter` became `from training_reporter import TrainingReporter`, because `training_reporter.py` moved out of the folder `helpers/`.

`prepare_splits.py` and `environment.yml` are ours. `environment.yml` holds the conda packages of the agent's environment and, from pip, torch 2.9.0 with the packages it needs, at the agent's versions.
