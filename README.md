# RNA Structure Engine

Testing different RNA secondary structure prediction tools based on probing wet-lab results (e.g. SHAPE, DMS, etc.). For UCA1.

Generic RNA secondary structure prediction engine. Stage 0 currently runs one RNAstructure/Deigan line on the `test-1/xist` case slot, populated with Busan et al. 2019 *E. coli* 16S rRNA cell-free SHAPE-MaP positive-control data.

## Run

Stage 0:

```bash
make build
make run CASE=xist TEST=test-1 SCHEME=rnastr-deigan REAGENT=1M7
```

On this workstation, if Docker Hub is temporarily unavailable but the previously built local RNAstructure image exists, the same pipeline can be run with:

```bash
make run-local CASE=xist TEST=test-1 SCHEME=rnastr-deigan REAGENT=1M7
```

Stage 1, four-scheme horizontal expansion:

```bash
make stage1 CASE=xist TEST=test-1 REAGENT=1M7
```

Local fallback:

```bash
make stage1-local CASE=xist TEST=test-1 REAGENT=1M7 WASHIETL_SAMPLE_SIZE=5
```

Results are written to:

```text
outputs/test-1/xist/
```

## What Stage 0 Does

```text
inputs/test-1/xist/raw/Ecoli_cellfree_16S_1M7_profile.txt
  -> adapter: ShapeMapper profile to FASTA + RNAstructure .shape
  -> RNAstructure Fold with Deigan SHAPE pseudoenergies, m=1.8 and b=-0.6
  -> RNAstructure partition
  -> ProbabilityPlot text output
  -> per-position Shannon entropy
  -> reactivity + entropy profile plot
  -> RNAstructure draw structure figure from the Fold CT
```

Stage 0 does not run RNAfold, Zarringhalam, Washietl, MEA, efn2, ProbKnot, UCA1, DMS, or reference base-pair scoring.

## What Stage 1 Does

Stage 1 runs the same 16S / 1M7 input across four schemes:

```text
rnastr-deigan
rnafold-deigan
rnafold-zarringhalam
rnafold-washietl
```

The outputs are organized by scheme:

```text
outputs/test-1/xist/stage1/schemes/<scheme>/
outputs/test-1/xist/stage1/comparisons/
```

Each scheme emits MFE structure, pairing probabilities, Shannon entropy, MEA structure when numerically valid, and an energy summary. RNAstructure uses `Fold`, `partition`, `MaxExpect`, and `efn2`; RNAfold schemes use `RNAfold -p --MEA`.

Washietl requires a perturbation vector from `RNApvmin`. The default `WASHIETL_SAMPLE_SIZE=5` is intentionally bounded so the 16S benchmark finishes during engineering runs. With this bounded setting, ViennaRNA can emit `NaN` ensemble/MEA quantities for Washietl. The pipeline keeps the Washietl MFE structure and marks partition/MEA/ensemble/entropy as `done_with_numeric_warning` in:

```text
outputs/test-1/xist/stage1/algorithm_applicability.tsv
```

The pairwise consistency report excludes zero-pair or numerically invalid structures from comparison.

## Inputs

The current default input is a ShapeMapper tab-delimited profile:

```text
inputs/test-1/xist/raw/Ecoli_cellfree_16S_1M7_profile.txt
```

Required columns:

```text
Nucleotide    Sequence    Norm_profile
```

The adapter writes:

```text
outputs/test-1/xist/prepared/ecoli_16s_cellfree_1m7_rnastr-deigan.fa
outputs/test-1/xist/prepared/ecoli_16s_cellfree_1m7_rnastr-deigan.shape
outputs/test-1/xist/prepared/ecoli_16s_cellfree_1m7_rnastr-deigan.converted.tsv
```

`Norm_profile` is used as the SHAPE reactivity column. `nan` values are written as RNAstructure missing values (`-999`). Negative normalized values are clamped to `0` for folding and preserved in the converted TSV as raw values.

## Outputs

The main output files are:

```text
outputs/test-1/xist/ecoli_16s_cellfree_1m7_rnastr-deigan.ct
outputs/test-1/xist/ecoli_16s_cellfree_1m7_rnastr-deigan.dbn
outputs/test-1/xist/ecoli_16s_cellfree_1m7_rnastr-deigan.pfs
outputs/test-1/xist/ecoli_16s_cellfree_1m7_rnastr-deigan.pairing_probability.txt
outputs/test-1/xist/ecoli_16s_cellfree_1m7_rnastr-deigan.shannon_entropy.tsv
outputs/test-1/xist/ecoli_16s_cellfree_1m7_rnastr-deigan.reactivity_entropy.png
outputs/test-1/xist/ecoli_16s_cellfree_1m7_rnastr-deigan.reactivity_entropy.svg
outputs/test-1/xist/ecoli_16s_cellfree_1m7_rnastr-deigan.mfe_structure.png
outputs/test-1/xist/ecoli_16s_cellfree_1m7_rnastr-deigan.mfe_structure.svg
outputs/test-1/xist/manifest.json
outputs/test-1/xist/run.log
```

Stage 1 outputs include:

```text
outputs/test-1/xist/stage1/schemes/rnastr-deigan/
outputs/test-1/xist/stage1/schemes/rnafold-deigan/
outputs/test-1/xist/stage1/schemes/rnafold-zarringhalam/
outputs/test-1/xist/stage1/schemes/rnafold-washietl/
outputs/test-1/xist/stage1/comparisons/pairwise_structure_consistency.tsv
outputs/test-1/xist/stage1/comparisons/pairwise_mfe_jaccard.png
outputs/test-1/xist/stage1/comparisons/pairwise_mfe_jaccard.svg
outputs/test-1/xist/stage1/algorithm_applicability.tsv
outputs/test-1/xist/stage1/manifest.json
```

## How To Read The Outputs

Lower per-position Shannon entropy indicates a more concentrated pairing ensemble and therefore a more determined local structural model. Pairing probabilities report ensemble support for base pairs, not a crystallographic truth label.

Stage 0 proves the engine can run the RNAstructure/Deigan line end to end on a positive-control rRNA profile. It does not yet prove benchmark accuracy against reference base pairs. Reference sensitivity, PPV, and FDR scoring should be added in the benchmark phase.

## Algorithm Matrix

| Algorithm / Scheme | RNAstr-Deigan | RNAfold-Deigan | RNAfold-Zarringhalam | RNAfold-Washietl |
|---|:---:|:---:|:---:|:---:|
| MFE / Fold | Stage 0 | Stage 1 | Stage 1 | Stage 1 |
| partition / pairing probability | Stage 0 | Stage 1 | Stage 1 | Stage 1 |
| MEA | Stage 1, `MaxExpect` | Stage 1, `RNAfold --MEA` | Stage 1 | Stage 1 |
| ProbKnot | Stage 3 | N/A | N/A | N/A |
| efn2 / ensemble energy | Stage 1 | Stage 1 | Stage 1 | Stage 1 |
| Shannon entropy | Stage 0 | Stage 1 | Stage 1 | Stage 1 |

N/A means: `N/A - RNAfold cannot represent pseudoknots`. RNAfold also has no command named `MaxExpect`; MEA structures are produced with `RNAfold --MEA`.

RNAstructure and ViennaRNA ensemble energy and Shannon entropy values must not be directly compared numerically because Turner parameter versions and partition implementations can differ. Cross-tool checks are implementation consistency checks at the structure/base-pair/domain level, not cross-validation.

ProbKnot is an RNAstructure subcommand that consumes RNAstructure partition pairing probabilities:

```text
partition -> ProbKnot -> pseudoknot-capable structure
```

It is not a preprocessing step before both tools.

## Versions

The Stage 0 Docker environment pins:

| Component | Version |
|---|---|
| Base image | `mambaorg/micromamba:1.5.10` |
| Python | `3.12.4` |
| RNAstructure | `6.5` |
| ViennaRNA | `2.7.2` |
| ImageMagick | `7.1.1_39` |

Stage 0 uses RNAstructure `draw` for MFE structure images. VARNA is not part of the Stage 0 dependency set.

## Citations

Busan, S., Weidmann, C. A., Sengupta, A., & Weeks, K. M. (2019). Guidelines for SHAPE reagent choice and detection strategy for RNA structure probing studies. *Biochemistry*, 58(21), 2655-2664.

Reuter, J. S., & Mathews, D. H. (2010). RNAstructure: Software for RNA secondary structure prediction and analysis. *BMC Bioinformatics*, 11, 129.

Deigan, K. E., Li, T. W., Mathews, D. H., & Weeks, K. M. (2009). Accurate SHAPE-directed RNA structure determination. *PNAS*, 106(1), 97-102.

Lorenz, R., et al. (2011). ViennaRNA Package 2.0. *Algorithms for Molecular Biology*, 6, 26.

Zarringhalam, K., et al. (2012). Integrating chemical footprinting data into RNA secondary structure prediction. *PLoS ONE*, 7(10), e45160.

Washietl, S., et al. (2012). RNA folding with soft constraints. *Nucleic Acids Research*, 40(10), 4261-4272.

Lorenz, R., et al. (2016). RNA folding with hard and soft constraints. *Algorithms for Molecular Biology*, 11, 8.

Bellaousov, S., & Mathews, D. H. (2010). ProbKnot: Fast prediction of RNA secondary structure including pseudoknots. *RNA*, 16(10), 1870-1880.
