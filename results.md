# Results

Agentomics and the six published tools on `test`: 687 labelled samples of cultured human cells (342 senescent, 345 not senescent) in 78 datasets from 71 laboratory groups. No sample of these laboratories was in the data that Agentomics was trained or selected on, and none is in the lists of the data that the tools' authors used. Each sample was scored alone: Agentomics computes a sample's probability from that sample's counts only, and each tool was given a count table that held that one sample. The exception is hUSI (below). Values are rounded to three decimals; `results.tsv` holds them to four, with the counts behind each value, and `predictions/` holds the score of every sample.

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

## The measures

- **AUROC across studies.** For two studies from different laboratory groups, the probability that a senescent sample of one study scores higher than a not senescent sample of the other; the mean over all such ordered pairs of studies. It asks whether a score means the same in different laboratories. A random score gives 0.5.
- **Pooled AUROC.** All samples are ranked together by their score. It is the share of all pairs of a senescent and a not senescent sample in which the senescent sample scores higher (ties count half). A random score gives 0.5.
- **Pooled AUPRC.** In the same single ranking, the average precision with which the senescent samples are found (scikit-learn's average precision, tied scores taken as one threshold). A random score gives about the share of senescent samples: 0.498 on the 687 samples, 0.470 on markeR's 645.
- **95% interval.** The 2.5th and 97.5th percentiles over 2,000 resamples of the laboratory groups, drawn with replacement (seed 20261116). A resample keeps all samples of the groups it draws, so the interval reflects how much the result depends on which laboratories are in `test`. Every method uses the same resamples.
- **Agentomics minus a method.** The difference of a measure between Agentomics and the method, computed in each of the same 2,000 resamples; its interval is the 2.5th and 97.5th percentiles of these differences.
- **Samples.** The samples a method has a score for. markeR's wrapper stopped on 42 samples (a score that was not a finite number), so markeR is scored, and compared with Agentomics, on the other 645.

## Notes on the methods

- **Agentomics** is scored by its `probability_1` (the probability of senescent) from `agentomics_model/inference.py` with the shipped model file. The run's prediction file has sha256 `1de842f0f4c369dfeccb9321d2b255c4a53ebf51c17315fa373ccec86ac01944`.
- **hUSI** is scored by `hUSI_unscaled`, the correlation before hUSI's rescaling, which the wrapper writes and the tool does not release. hUSI's wrapper refuses a single sample, so its scores come from runs on each dataset's labelled samples together. The unscaled correlation of a sample depends on that sample alone.
- **markeR** is scored by the Hernandez-Segura gene set with `ssGSEA`, because markeR's `logmedian` is 0 for every sample scored alone.
- **SenFinder** runs with 19 of its 20 genes; the gene without an NCBI gene id is set to 0, as its README prescribes for a missing value.
- The scoring also covered 19 published gene sets scored with markeR's `ssGSEA`. senBench reports the models only.

## Where the numbers come from

`predictions/` holds the score of every labelled `test` sample by each of the seven models, one file per model (`predictions/README.md` names the source of each). `scoring/score.py` computes the measures, the intervals and the differences from these files and the labels of `test` (`scoring/test_samples.tsv`); from the repository's folder,

```
python scoring/score.py --predictions predictions/*.tsv --reference Agentomics --out build/results.tsv
```

gives `results.tsv` byte for byte. `results.tsv` was written from the tables of our scoring run, whose scorer `score.py` copies: the measure code is the same, and only the reading of the scores differs (`scoring/README.md`). Its columns are the method, the measure, the value with its interval, the value of a random score, the samples (with the senescent and not senescent ones), the datasets and laboratory groups, and the reference minus the method with its interval and the samples both have scores for.
