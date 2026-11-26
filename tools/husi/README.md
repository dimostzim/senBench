# hUSI

hUSI, the human universal senescence index, is the tool of Wang 2025 (Nat Aging, DOI 10.1038/s43587-025-00886-2). Its weights come from a one-class logistic regression fitted on 385 senescent bulk RNA-seq samples of cultured human cells: one weight for each of 13,346 genes. The score of a sample is the Spearman correlation, over the genes it shares with the weights, between the weights and the sample's TPM values. The released function then rescales the correlations of the samples of its input from 0 (the lowest) to 1 (the highest).

- Repository: `https://github.com/WJPina/HUSI`, commit `300da32f111f8c8603e822eca8e84c9b5a76db48` (2025-09-11). Licence: Apache-2.0.
- Needs: TPM values of genes named by symbol. The wrapper computes TPM from the counts and the gene lengths of GEO's annotation table, and reads the weights and functions of the clone (`sc/hUSI.R`).

## Score in senBench

`hUSI_unscaled`: the correlation before hUSI's rescaling. It is not a released score. The released function `cal_hUSI` returns only the rescaled score `hUSI`; the wrapper also writes the correlation, computed with the authors' function `scoreOCLR` (`R/functions_new.R`), which their analysis scripts call before rescaling. A higher value means more senescent. `hUSI_unscaled` of a sample depends on that sample alone; `hUSI` depends on the other samples of the input, because the lowest and the highest sample define 0 and 1.

## Scoring one sample

```
bash tools/husi/setup.sh /path/to/tools_install
bash tools/husi/wrapper/run_tool.sh --counts <counts.tsv.gz> --annotation <Human.GRCh38.p13.annot.tsv.gz> --out-dir <new folder> --env-prefix /path/to/tools_install/envs/husi --tool-repository /path/to/tools_install/HUSI
```

The wrapper refuses a count table with one sample, because the released score cannot be rescaled over one sample. The benchmark therefore did not call hUSI on one sample: it took `hUSI_unscaled` from the wrapper's runs on the count table of each dataset with the labelled samples of that dataset. To score one sample, give it together with at least one other sample whose correlation differs; its `hUSI_unscaled` is the same whichever other sample is given (checked on GSE130727: 0.04712 for its first sample with either of two other samples).

## What the wrapper changes against the tool's published use

1. TPM is computed in the wrapper from the counts and the gene lengths of the annotation table. The authors used TPM from StringTie or from the data sources. On GSE130727 the wrapper's TPM equals GEO's own TPM table to the four digits GEO prints.
2. TPM values are scored as they are. In one analysis the authors first log-transform, scale and centre each gene over pooled samples; their applications to bulk data do not.
3. 9 of the 13,346 genes with a weight are not in the annotation table under their symbol.
4. The wrapper writes the second score, `hUSI_unscaled`; the released function returns only `hUSI`.
5. `call` is the class of the tool's function `SSE_hUSI`. The tool's other classifier, `GMM_hUSI`, is not run.
6. R packages are built at the versions the repository states (R 4.0.5, dplyr 1.1.4, mclust 6.0.1). Seurat, ggplot2 and pROC are not installed; the scoring functions do not call them.

The wrapper's `README.md` gives the source of each point and the data the weights were fitted on.
