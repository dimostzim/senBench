# senBench

senBench is a benchmark of methods that call a bulk RNA-seq sample of cultured human cells senescent or not senescent. A method receives the raw read counts of one sample and returns a score that is higher for senescent samples. senBench holds public datasets from GEO and ArrayExpress with a label for every sample and a split of the datasets by laboratory. It compares Agentomics, a model that the agent system Agentomics-ML built on the training data, with six published tools on the held-out part `test` (687 samples of 78 datasets from 71 laboratory groups).

## Quick start

```
bash download.sh
bash build.sh
bash run.sh
```

The three commands need Linux (x86-64) with bash, python3, git and curl. `download.sh` downloads the public files into `downloads/` (18.76 GiB). `build.sh` builds the count tables, the split, the training data and `test` into `build/`, and compares each output with the sha256 of our build. `run.sh` installs our model and the six tools, scores `test` with each into `predictions/`, and writes `results.tsv`; `bash run.sh SenCID` runs one method.

## Results

| Method | AUROC across studies | Pooled AUROC | Pooled AUPRC |
| --- | --- | --- | --- |
| Agentomics | 0.834 (0.775 to 0.883) | 0.864 (0.797 to 0.913) | 0.878 (0.819 to 0.923) |
| hUSI | 0.772 (0.716 to 0.824) | 0.800 (0.728 to 0.862) | 0.822 (0.750 to 0.883) |
| SenCID | 0.752 (0.691 to 0.808) | 0.777 (0.708 to 0.837) | 0.795 (0.727 to 0.857) |
| SenCat | 0.750 (0.699 to 0.800) | 0.761 (0.697 to 0.821) | 0.772 (0.715 to 0.832) |
| SENCAN | 0.721 (0.669 to 0.768) | 0.722 (0.678 to 0.771) | 0.710 (0.628 to 0.784) |
| markeR | 0.663 (0.615 to 0.711) | 0.702 (0.635 to 0.766) | 0.725 (0.636 to 0.799) |
| SenFinder | 0.621 (0.589 to 0.656) | 0.636 (0.593 to 0.689) | 0.609 (0.572 to 0.662) |

- **AUROC across studies:** the probability that a senescent sample of one study scores higher than a not senescent sample of a study from another laboratory group, averaged over all such pairs of studies.
- **Pooled AUROC:** the probability that a senescent sample scores higher than a not senescent sample, with all samples in one ranking.
- **Pooled AUPRC:** the average precision with which that ranking finds the senescent samples; a random score gives about 0.5.
- **Brackets:** the 95% interval over 2,000 resamples of the laboratory groups. markeR gave no score for 42 samples and is scored on the other 645. `results.tsv` also holds Agentomics minus each tool.

## Folders

- `pipeline/`: the code that downloads, builds, predicts and scores, with the labels and the other input tables.
- `agentomics_model/`: the code and the trained model of Agentomics.
- `tools/`: for each published tool, its install script and the script that runs it on a count table.
- `predictions/`: the score that each method gave each `test` sample.

## Limitations

- Language models read the labels from the public records; no person checked them.
- The split was made after results on earlier splits of the same datasets had been seen, so senBench is a test bed.
- `test` holds quiescent controls in 2 studies only (7 samples).
- Each tool is run on one sample at a time, which some tools were not designed for.
- A plain logistic regression trained on the same data reaches an AUROC across studies of 0.801, above every tool.

The code of this repository is under the MIT licence (`LICENSE`); each tool keeps its own licence.
