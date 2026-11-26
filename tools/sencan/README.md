# SENCAN

SENCAN is the classifier of Jochems 2021 (Cell Rep, DOI 10.1016/j.celrep.2021.109441): an elastic net on 137 genes that gives each sample a probability of being senescent. The paper reports no code. The classifier file in the hUSI repository, `Data/SENCAN_classifier.rda`, equals the paper's Table S4 and the file the authors' website offers (same sha256), and the wrapper reads it from there.

- Classifier file: `https://github.com/WJPina/HUSI`, commit `300da32f111f8c8603e822eca8e84c9b5a76db48`. Licence: none given for SENCAN; the hUSI repository is Apache-2.0.
- Needs: read counts per Ensembl gene id. The classifier normalises the input samples against its reference samples, takes log counts per million, centres and scales them, and fills a gene that the input lacks with the mean count of the two reference samples that correlate best with the sample. The authors' scoring script (`SENCAN_example.r` on their website) takes these steps, and the wrapper follows it.

## Score in senBench

`SENCAN_score`, the probability of senescent from the elastic net, from 0 to 1; `call` is senescent above 0.5.

## Scoring one sample

```
bash tools/sencan/setup.sh /path/to/tools_install
bash tools/sencan/wrapper/run_tool.sh --counts <one_sample.tsv.gz> --annotation <Human.GRCh38.p13.annot.tsv.gz> --out-dir <new folder> --env-prefix /path/to/tools_install/envs/sencan --tool-repository /path/to/tools_install/HUSI
```

The benchmark called the wrapper with one sample per count table. All samples of a table are normalised together, so with one sample the normalisation factor is 1 and the sample is not adjusted towards the authors' reference sample.

## What the wrapper changes against the tool's published use

1. The paper quantifies raw reads with kallisto on GENCODE 34. The wrapper starts from read counts per NCBI gene, so this step cannot be followed.
2. Genes that the input lacks are filled per sample with the mean count of the two reference samples that correlate best with it, as the authors' script does. GEO's count tables lack 18 of the 137 genes, which carry 21.6% of the summed absolute coefficients, among them the largest positive and the most negative coefficient.
3. The filled-in counts are not rescaled to the depth of the sample, as in the authors' script.
4. All samples of a count table are normalised together, as in the authors' script. A sample's score therefore changes with the other samples of the table.
5. Genes are matched by Ensembl gene id (119 of 137 found), as in the authors' script.
6. `call` is senescent above 0.5, the dashed line of the paper's Figure 6. The text of the paper states no threshold.
7. R 4.0.5 and edgeR 3.32.1 replace the paper's R 3.6.0 and edgeR 3.26.8.

The wrapper's `README.md` gives the source of each point; its sentence that the authors' script was not read predates the reading.
