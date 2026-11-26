# SenCID

SenCID is the tool of Tao 2024 (Cell Metab, DOI 10.1016/j.cmet.2024.03.009). It holds six models, SID1 to SID6, one for each of six senescence identities found in the authors' data; each scores a sample from 0 to 1 with a logistic regression on 13 to 224 genes. A seventh model recommends which of the six to use for the sample, from 1,290 genes.

- Repository: `https://github.com/JackieHanLab/SenCID`, commit `f93cba19b9d44e84f89504ecc6aaeb108ce8e508` (2024-02-22, version 1.0.0). Licence: MIT.
- Needs: raw counts of genes named by HGNC symbol, with denoising switched off for bulk RNA-seq (the first author's advice in the repository's issues). SenCID normalises the counts itself: it keeps the genes of its list of 1,290 that have a total count of at least 1 over the samples of the input, gives every other gene of the list the value 0.01, divides by the sum over the 1,290 genes, takes log2(value + 1) and turns the values of each sample into z-scores. The environment keeps every version pin of the repository (Python 3.6.15, TensorFlow 2.2.0, scikit-learn 0.24.2).

## Score in senBench

`RecSID_Score`: the score of the model that SenCID recommends for the sample, as the tool's tutorial uses it. A higher value means more senescent; `call` is SenCID's own call at 0.5. The wrapper also writes the scores of all six models.

## Scoring one sample

```
bash tools/sencid/setup.sh /path/to/tools_install
bash /path/to/tools_install/wrappers/sencid/run_tool.sh --counts <one_sample.tsv.gz> --annotation <Human.GRCh38.p13.annot.tsv.gz> --out-dir <new folder>
```

SenCID's wrapper reads the paths of its environment and clone from `config.sh`, so it runs from the copy that `setup.sh` writes into the install folder, whose `config.sh` names that folder. The benchmark called it with one sample per count table. With one sample, every gene of the list with no read in that sample gets the value 0.01; with several samples, only genes with no read in any of them do.

## What the wrapper changes against the tool's published use

1. Denoising is off. The default in the code is on; the first author prescribes off for bulk RNA-seq in the repository's issues.
2. The tool recommends one of its six models. The wrapper takes the recommendation per sample, because it is given no cell type; the tutorial picks one model per cell type.
3. 10 of the tool's 1,290 genes are absent under their current symbol. SenCID's own code fills them with 0.01.
4. The environment keeps every version pin of the repository; packages come from conda-forge and PyPI only.
5. `HOME` points to a temporary folder for the SenCID process, because two of its libraries write settings to the home folder.

The paper's STAR Methods agree with these choices: for bulk samples, each sample gets the model with the highest recommendation index, and 0.5 is the threshold for a call.

The wrapper's `README.md` gives the source of each point and the data the models were fitted on; its statement that the paper's Methods were not read predates their reading.
