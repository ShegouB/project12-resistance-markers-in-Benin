# Week 13, Day 3: differential expression with DESeq2 (GSE62132, K1 strain)

Block 6 (Transcriptomics and Functional Genomics), Week 13. Topic: negative
binomial model, size factors, dispersion shrinkage, Wald test, LFC shrinkage,
and robustness checks at very small sample size.

## Objective

Test DHA vs vehicle (0.1% DMSO) differential expression in P. falciparum,
separately for ring, trophozoite and schizont, using the Day 2 count matrix.
With only two replicates per condition, a large part of the work was checking
which results survive scrutiny.

## Inputs

- `data/counts_all.txt`: featureCounts matrix from Day 2 (5,611 genes x 12
  samples: DHA and vehicle, three stages, two replicates each).
- `data/samples.tsv`, `data/annotation.gff` (NCBI RefSeq GCF_000002765.6; used
  for gene names and products).
- For scripts 07 and 08 only: the six untreated runs, aligned and counted with
  the Day 2 pipeline (`counts_untreated.txt`, `samples_untreated.tsv`,
  `scripts/02b_align_count_untreated.sh`, which live in the Day 2 folder).

## Method (primary analysis, `scripts/03_deseq2.R`)

- One DESeq2 model per stage, design `~ treatment`, reference level `vehicle`.
  Stages were analysed separately because they differ far more than any drug
  effect: in the all-sample PCA, PC1 (70.4%) separates ring from the other two
  stages and PC2 (19.8%) separates trophozoite from schizont.
- Genes kept with at least 10 counts in at least 2 samples of the stage.
- Wald test, BH adjustment, `padj < 0.05`. apeglm shrinkage for ranking and
  plots; significance comes from the unshrunk test.
- VST (blind) only for PCA and sample distances, never for testing.

## Primary results (2-group model, DHA vs vehicle)

| Stage | Genes tested | padj < 0.05 | Up / down | \|LFC\| >= 1 |
|---|---|---|---|---|
| Ring | 3,090 of 3,497 | 413 | 113 / 300 | 101 |
| Trophozoite | 4,497 of 4,497 | 1,703 | 923 / 780 | 676 |
| Schizont | 4,301 of 4,301 | 0 | 0 / 0 | 0 |

Size factors (DHA rep1, DHA rep2, vehicle rep1, vehicle rep2): ring 1.066,
1.086, 0.996, 0.882; trophozoite 1.056, 1.211, 0.890, 0.893; schizont 1.098,
1.079, 0.926, 0.903. DHA libraries get larger size factors than vehicle
libraries at similar depth, which reflects composition (see below).

## Checks and what they showed

1. **Within-stage PCA.** In ring and trophozoite, PC1 (91% and 98% of
   within-stage variance) separates DHA from vehicle, and the two replicates of
   each condition sit close together. In schizont, PC1 separates replicate 1
   from replicate 2 in both conditions, not DHA from vehicle.
2. **Composition explains the size-factor split.** Ribosomal protein genes are
   a large share of each library and fall in DHA (trophozoite: about 10.4% in
   vehicle vs about 6.9% in DHA; ring: about 17.5% vs about 16%). When a large
   block shrinks, reads redistribute to other genes, raising the DHA size
   factors. Re-estimating size factors without ribosomal protein and rRNA genes
   changed DEG counts by 1% or less (ring 413 to 423, trophozoite 1,703 to
   1,712, schizont 0).
3. **Replicate pairing.** `~ replicate + treatment` gave fold changes with
   Spearman >= 0.99 against the primary model (ring 381 DEGs, Jaccard 0.81;
   trophozoite 1,278 DEGs, Jaccard 0.74; schizont 0). Whether replicates were
   actually paired by experiment is not stated in GEO.
4. **Ribosomal protein genes, tested properly.** A Wilcoxon test on shrunken
   LFCs was biased (apeglm shrinks low-expression genes harder, so highly
   expressed gene sets differ by expression alone, and schizont gave p = 7.6e-13
   with zero DEGs). The replacement compares the Wald statistic against random
   gene sets matched on expression decile (5,000 draws; p hits the 1/5001
   floor, and the null ignores co-regulation, so only effect sizes are
   reported): ring median Wald -2.74 vs null -1.19 (81% of genes negative),
   trophozoite -5.19 vs -0.58 (69% negative), schizont +0.69 vs +0.20 (26%
   negative).
5. **rRNA outlier.** Five rRNA loci appeared about 6 to 8 times lower in ring
   DHA, but this is a library composition difference: rRNA is 0.46% and 0.71% of
   counted fragments in ring vehicle, 0.11% and 0.43% in ring untreated, 0.12%
   and 0.10% in ring DHA, and 0.02% to 0.16% in the other eight samples. GEO
   lists a single extraction protocol and library selection (cDNA) for all 18
   samples, so the cause is undocumented. rRNA loci were excluded from
   interpretation.
6. **Retracted analysis.** Scripts 05 (part 1) and 06 (part B) regressed DHA
   fold changes on stage axes built from the vehicle samples. That is circular,
   because the vehicle samples are also the denominator of the DHA fold change,
   so those R2 values and sign-agreement percentages are not used. The stage
   test was redone with axes from the independent untreated samples (script 07):
   R2 0.115 in ring and 0.201 in trophozoite (0.101 and 0.209 without ribosomal
   and rRNA genes). Slopes were small, and the two trophozoite slopes
   (-0.137 ring-to-trophozoite, +0.127 trophozoite-to-schizont) point in
   opposite directions, which does not describe a single developmental shift.
   There is a modest stage-related component and no uniform delay.

## Three-group model (DHA, vehicle, untreated; scripts 07 and 08)

| Stage | Vehicle vs untreated | DHA vs untreated | DHA vs vehicle | Jaccard of the two DHA lists | LFC Spearman |
|---|---|---|---|---|---|
| Ring | 1,357 | 1,256 | 81 | 0.02 | 0.067 |
| Trophozoite | 317 | 1,566 | 1,312 | 0.51 | 0.783 |
| Schizont | 225 | 340 | 0 | 0.00 | 0.321 |

Distances between VST group means (top 500 variable genes) and between the two
replicates of each group:

| Stage | Replicate distance (untreated / vehicle / DHA) | DHA to vehicle | DHA to untreated | Vehicle to untreated |
|---|---|---|---|---|
| Ring | 23.2 / 10.5 / 8.0 | 18.8 | 36.8 | 37.1 |
| Trophozoite | 10.8 / 4.5 / 4.4 | 30.4 | 34.0 | 10.5 |
| Schizont | 10.7 / 12.5 / 15.1 | 8.9 | 20.4 | 20.6 |

Consensus genes (significant and in the same direction against both vehicle and
untreated): ring 7 (4 up, 3 down), trophozoite 973 (593 up, 380 down),
schizont 0. Files: `results/consensus_<stage>.csv`.

## Conclusions

**Trophozoite: supported.** DHA is 30.4 units from vehicle and 34.0 from
untreated, while vehicle is 10.5 from untreated (the same as the distance between two untreated
replicates, 10.8). The DHA vs vehicle and DHA vs untreated lists overlap well
(Jaccard 0.51, LFC Spearman 0.78). Ribosomal protein genes are lower in DHA
against both controls (median Wald -3.43 vs vehicle, -3.50 vs untreated;
69% and 75% of genes negative). Exported-protein genes are enriched among the
consensus up genes (59 of 593, 9.9%, vs 2.2% of tested genes; OR about 10).
Top consensus genes include HAD3 (PF3D7_1226100), SBP1 (PF3D7_0501300), SEMP1
(PF3D7_0702400), FIKK4.1 (PF3D7_0424500), HSP101 (PF3D7_1116800), a
heat shock protein 40 type II (PF3D7_0113700) and exported protein family 1
(PF3D7_0114000), with nucleoside phosphorylase (PF3D7_0513300) among the down
genes.

**Schizont: negative control.** DHA vs vehicle has no DEGs, and DHA and vehicle
are close (8.9) while both sit about 20 from untreated. The 225 and 340 DEGs
against untreated therefore reflect a DMSO or handling difference shared by
both treated groups.

**Ring: unresolved.** The two untreated ring replicates are 23.2 apart, more
than DHA is from vehicle (18.8), so the effect is smaller than between-replicate
disagreement in one control group. The DEG count is unstable (413 in the
2-group model, 81 in the 3-group model) and only 7 genes are consensus (EVP1,
HSP70, exported lipase 1 and YHM2 up; M1AAP, TKL4 and falcilysin down).
Ribosomal protein genes look lower than vehicle (median Wald -1.43) but higher
than untreated (+2.85), so the ribosomal drop in ring is mostly vehicle being
high. A real but weak effect is likely, because the ring fold changes correlate
with trophozoite fold changes (Spearman 0.47 vs trophozoite-vehicle, 0.42 vs
trophozoite-untreated), but individual ring genes are not trustworthy.

**Reproducibility across stages.** Ring and trophozoite DHA fold changes
correlate (Spearman 0.49; 0.45 without ribosomal protein and rRNA genes). Of
268 genes significant in both, 90% change in the same direction. Of 78 ring-up
genes, 61 are also up in trophozoite (about 20 expected by chance). A handling
artifact could in principle reproduce this, but the schizont null argues
against it.

## Limitations

- Two replicates per condition: 2 residual degrees of freedom per stage.
  DESeq2 runs at 2 vs 2 because dispersions are shrunk toward a trend fitted
  across all genes, which gives stability at the cost of power for modest
  changes and sensitivity to the fitted trend.
- The vehicle itself differs from untreated in every stage (225 to 1,357 DEGs),
  so DMSO or handling matters and the choice of control changes the ring result.
- The exported-protein keyword list (PHIST, EMP1, FIKK, PRESAN, RESA, exported,
  GEXP, ETRAMP) was chosen after looking at top-gene lists, so its p-values are
  optimistic. The consistent direction in two stages is the informative part.
- The consensus definition (significant against both controls) is conservative
  but not ideal: a gene whose DMSO response opposes its DHA response can fail
  against untreated. Vehicle is the matched control and is the primary
  comparison.
- The cause of the ring library composition differences (rRNA 0.10% to 0.71%,
  ribosomal share 9.7% to 18.3%) is unknown.
- Shallow depth (1 to 2 million fragments) limits testing of weakly expressed
  genes.
- No functional enrichment was run on Day 3, and no biological interpretation of
  the trophozoite consensus list beyond the gene families named above.

## Folder contents

```
week13/day3_deseq2/
  data/        counts_all.txt, samples.tsv, annotation.gff
  results/
    counts_matrix.csv, coldata.csv
    pca_all_samples.png, sample_distance_heatmap.png
    deseq2_<stage>_all_genes.csv, deseq2_<stage>_significant.csv
    deseq2_summary_by_stage.csv, sessionInfo_day3.txt
    consensus_<stage>.csv, stage_axis_residuals_<stage>.csv
    deseq2_console.txt, qc_checks_console.txt, stage_shift_console.txt,
    confounds_console.txt, untreated_reference_console.txt, consensus_console.txt
  scripts/
    03_deseq2.R              per-stage DESeq2, PCA, annotation
    04_qc_checks.R           PCA by stage, composition, matched-null ribosomal test
    05_stage_shift.R         stage regression (part 1 retracted), paired design
    06_confounds.R           rRNA share, exported-protein enrichment, cross-stage
    07_untreated_reference.R untreated-based stage axes, 3-group model, GEO fields
    08_consensus.R           consensus genes, replicate distances, cross-stage matrix
```

## Reproduce

```bash
# inside the Kaggle notebook, R environment rnaseq_r (R 4.5.3, DESeq2 1.50.2, apeglm 1.32.0)
micromamba run -n rnaseq_r Rscript scripts/03_deseq2.R
for s in 04_qc_checks 05_stage_shift 06_confounds 07_untreated_reference 08_consensus; do
  micromamba run -n rnaseq_r Rscript scripts/$s.R
done
```

Scripts 07 and 08 need the Day 2 untreated counts (see Inputs).

## Next

Day 4: edgeR (TMM normalisation, exact test, quasi-likelihood fit) on the same
data, with the Jaccard overlap against these DESeq2 lists.
