# markeR

markeR is an R package of Martins-Silva 2026 (NAR Genom Bioinform, DOI 10.1093/nargab/lqag057) that scores gene sets in bulk RNA-seq samples. It has no trained model and no rule that calls a sample senescent. It ships three senescence gene sets (`HernandezSegura`, `REACTOME_CELLULAR_SENESCENCE` and `LiteratureMarkers`) and three scoring methods (`logmedian`, `ssGSEA` and `ranking`); the wrapper writes all nine scores of a sample.

- Repository: `https://github.com/DiseaseTranscriptomicsLab/markeR`, commit `c6528f0fa1570602b077e4688865dcd68af3d058` (2026-07-08, package version 1.3.1). Licence: Artistic-2.0.
- Needs: filtered and normalised counts of genes named by symbol. The wrapper names the genes by the symbol of GEO's annotation table, keeps the genes with a mean raw count above 70 over the samples of the input, normalises with edgeR's TMM factors and limma's `voom`, as the authors' preprocessing does, and calls `markeR::CalculateScores` with the package's gene sets.

## Score in senBench

`HernandezSegura:ssGSEA`: the gene set of Hernandez-Segura 2017 (Curr Biol; 55 genes, 32 up and 23 down in senescence; 53 of them are in the annotation table) scored with markeR's own ssGSEA, which ranks the genes within the sample. The score of the down genes is subtracted from that of the up genes, so a higher score means more senescent. The wrapper's main score is `HernandezSegura:logmedian`, but `logmedian` centres every gene on its median over the samples of the input and is 0 for every gene set when the input holds one sample; the benchmark therefore scores `ssGSEA`.

## Scoring one sample

```
bash tools/marker/setup.sh /path/to/tools_install
bash tools/marker/wrapper/run_tool.sh --counts <one_sample.tsv.gz> --annotation <Human.GRCh38.p13.annot.tsv.gz> --out-dir <new folder> --env-prefix /path/to/tools_install/envs/marker
```

The benchmark called the wrapper with one sample per count table. With one sample, the expression filter keeps the genes with a count above 70 in that sample. When one direction of a gene set keeps no gene after the filter, its `ssGSEA` score is not a finite number, and the wrapper stops without output, because it stops when any of its nine scores is not finite. On the 687 `test` samples it stopped 42 times, so markeR is scored on 645 samples.

## What the wrapper changes against the tool's published use

1. The authors keep a gene when its mean count is above 70 in at least one condition. The wrapper sees no conditions, so it keeps a gene when its mean count over all samples of the input is above 70. This is stricter.
2. Each count table is scored alone. The paper processed 545 samples of 25 datasets together and corrected for the dataset.
3. Counts of gene ids that share a symbol are added (one case).
4. Two genes of `HernandezSegura` (`FAM214B`, `MT-CYB`) are not in the annotation table under these symbols.
5. markeR 1.3.1 replaces the paper's 1.0.0; the scoring files are the same in both.
6. The authors name no main score; the wrapper takes `HernandezSegura:logmedian` (above for the benchmark's choice).

Kept as released, although the documentation says otherwise: `ranking` divides by the number of genes in the matrix; `ssGSEA` is the package's own function; the default method in the code is `ssGSEA`. The wrapper stops for the whole input when any of its nine scores is not a finite number for any sample.

The wrapper's `README.md` gives the source of each point and the origin of the three gene sets.
