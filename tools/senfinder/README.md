# SenFinder

SenFinder is the released model of Mahmud 2024 (GeroScience, DOI 10.1007/s11357-024-01485-6): a linear support vector machine on the TPM values of 20 genes, which scales each gene with the means and standard deviations of the authors' training samples stored in the model.

- Repository: `https://github.com/XiaoDongLab/SenFinder`, commit `5f56506300baab2ed4c20af5a52c6252c4578b50` (2024-12-07). Licence: AGPL-3.0.
- Needs: TPM values, not log-transformed, of the 20 genes, with no missing value (the README: use 0 instead). The wrapper computes TPM from the counts and the gene lengths of GEO's annotation table.

## Score in senBench

`Senescence_probability`: the probability of the class senescent that the released model returns, from 0 to 1. The released script writes only the class; the wrapper writes both. A higher value means more senescent.

One of the 20 genes, `Lnc_PLOD2_2` (ENSG00000261051), has no NCBI gene id, so a count table of the benchmark can never supply it. The benchmark runs the wrapper with `--set-absent-genes-to-zero`, which sets it to 0, as the repository's README prescribes for a missing value ("19 of 20 genes"). Every sample then moves by the same amount towards not senescent, so the ranking of the samples comes from the other 19 genes, and the calls are biased towards not senescent.

## Scoring one sample

```
bash tools/senfinder/setup.sh /path/to/tools_install
bash tools/senfinder/wrapper/run_tool.sh --counts <one_sample.tsv.gz> --annotation <Human.GRCh38.p13.annot.tsv.gz> --out-dir <new folder> --set-absent-genes-to-zero --env-prefix /path/to/tools_install/envs/senfinder --tool-repository /path/to/tools_install/SenFinder
```

The benchmark called the wrapper with one sample per count table. SenFinder's score of a sample does not depend on the other samples: on the `test` samples of an earlier round, the scores of one sample per call equalled those of the samples scored together (largest difference 1e-15).

## What the wrapper changes against the tool's published use

1. The absent gene is set to 0 (above).
2. TPM is computed in the wrapper from the counts and gene lengths. The authors used RSEM after STAR on GRCh37.
3. The model's column `YIU2` is read as the gene YJU2.
4. The wrapper writes the probability of the class senescent in addition to the class.

The wrapper's `README.md` gives the source of each point and the data the model was fitted on.
