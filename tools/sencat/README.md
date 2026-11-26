# SenCat

The SenCat RNA score is the score of Anerillas 2026 (Mol Cell, DOI 10.1016/j.molcel.2026.05.017). Its repository releases fixed weights: one coefficient for each of 5,000 Ensembl genes. The score of a sample is the sum, over the genes found in the input, of the coefficient times the logarithm of the gene's normalised count plus one. It has no threshold and no class call.

- Repository: `https://github.com/maragkakislab/wf-ml-markers-senescence`, commit `0e23fcc926f31fc2206631ff9b64a852fb81fe89` (2026-08-19). The archived release cited by the paper is on Zenodo (DOI 10.5281/zenodo.20008425). Licence: the repository has no licence file; the Zenodo record is CC BY 4.0.
- Needs: counts by Ensembl gene id, normalised with the tool's `normalize_counts.py`, which removes genes with a count in at most 5% of the samples or a total count of at most 380 and divides each sample by its DESeq2 size factor. The wrapper converts NCBI gene ids to Ensembl gene ids, then runs the tool's own `normalize_counts.py` and `marker_classifier.py`.

## Score in senBench

`rna_score`, the only score. A higher value means more similar to the senescence signature (the repository's README).

## Scoring one sample

```
bash tools/sencat/setup.sh /path/to/tools_install
SENESCENCE_TOOLS_DIR=/path/to/tools_install bash tools/sencat/wrapper/run_tool.sh --counts <one_sample.tsv.gz> --annotation <Human.GRCh38.p13.annot.tsv.gz> --out-dir <new folder>
```

`SENESCENCE_TOOLS_DIR` names the install folder, which holds the clone and `envs/`. The benchmark called the wrapper with one sample per count table. With one sample, the DESeq2 size factor is 1, and the gene filter keeps the genes with a count above 380 in that sample.

## What the wrapper changes against the tool's published use

1. NCBI gene ids are converted to Ensembl gene ids with the annotation table. 120 of the 5,000 genes can never be found with this table. Counts of NCBI genes that share an Ensembl id are added.
2. The repository clone replaces the archived release that the paper cites, because the archived script stops on samples without labels. The weights file is the same in both.
3. One environment at current package versions replaces the authors' two.
4. The tool's gene filter and size factors are computed over the samples of the input.

The wrapper's `README.md` gives the source of each point and the data the weights were fitted on.
