# Single-Nucleus Muscle NFR and Nucleosome QTL Analysis

This repository contains the code and workflow definitions used for single-nucleus ATAC-seq analyses of chromatin accessibility, nucleosome positioning, and downstream QTL interpretation in human skeletal muscle.

The repository centers on a set of Snakemake workflows in `smk/` plus one Nextflow workflow for LDSC analyses in `nextflow/`.

## Paper

This repository accompanies:

> **Wang X, Robertson C.C. et al.**
> **bioRxiv** (2025)

## What This Repository Does

The main analyses in this repository identify and interpret:

- `nfrQTLs`: QTLs for nucleosome-free regions
- `nucQTLs`: QTLs for nucleosomal regions
- `caQTLs`: QTLs for accessibility peaks derived from MACS2 summits and narrowpeaks

Downstream workflows then connect those QTLs to:

- NFR versus nucleosomal colocalization
- causal inference testing
- GWAS colocalization
- eQTL overlap
- motif enrichment
- multi-track locus plots
- LDSC annotation enrichment

## Repository Layout

```text
muscle_nfr_nucQTL/
├── bin/                    # R and Python analysis scripts
├── env/                    # Conda environments
├── nextflow/               # LDSC workflow
├── run_snakemake/          # Helper scripts to launch Snakemake
├── smk/                    # Snakemake workflows and YAML configs
└── README.md
```

## Requirements

This codebase was developed with:

- `R >= 4.2.3`
- `Python >= 3.8`
- `Snakemake >= 7.30.1`
- `Nextflow >= 24.10.6`

External tools used across workflows include:

- [QTLtools](https://qtltools.github.io/qtltools/pages/basic_install.html)
- [samtools](https://www.htslib.org/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [MACS2 / MACS3](https://github.com/macs3-project/MACS)
- UCSC utilities such as `bedGraphToBigWig`
- [deepTools](https://deeptools.readthedocs.io/)
- [HOMER](http://homer.ucsd.edu/homer/)
- featureCounts

Conda environment files are provided in `env/` and `smk/envs/`, but some workflows also depend on cluster modules, Singularity, and external reference files specified in the YAML configs.

## Before Running Anything

Every workflow is configured through files in `smk/`. Before launching a pipeline, inspect:

- `smk/<workflow>.yaml`
- `smk/<workflow>_cluster.yaml`
- the corresponding `smk/<workflow>.smk`

In particular, confirm:

- input BAM and VCF paths
- sample metadata paths
- chromosome size and genome reference paths
- output directories
- cluster account, time, memory, and CPU settings
- trait-specific metadata for GWAS workflows

Many configs contain absolute or environment-specific paths. They are not plug-and-play across systems.

## How To Run A Workflow

All Snakemake workflows are launched from `run_snakemake/`.

From the repository root:

```bash
cd run_snakemake
bash dryrun.sh <workflow>
bash run.sh <workflow>
```

Example:

```bash
cd run_snakemake
bash dryrun.sh peakcall
bash run.sh peakcall
```

This launches:

- Snakefile: `smk/<workflow>.smk`
- config: `smk/<workflow>.yaml`
- cluster config: `smk/<workflow>_cluster.yaml`

Useful helpers:

- `bash dryrun.sh <workflow>`: inspect the DAG before submitting jobs
- `bash run.sh <workflow>`: submit the workflow to Slurm
- `bash unlock.sh <workflow>`: unlock a working directory after interrupted runs

## Recommended Execution Order

The repository has many workflows, but not all are meant to be run in a single linear chain. The safest way to think about it is:

1. Run the core NFR and nucleosome pipeline first.
2. Run downstream interpretation workflows from those outputs.
3. Run optional caQTL, GWAS, motif, plotting, or LDSC workflows only if you need those analyses.

### Core Workflow Order

Run these first if you want to reproduce the main NFR and nucQTL analyses.

| Step | Workflow | Run when | Purpose |
| --- | --- | --- | --- |
| 1 | `peakcall` | first | Calls peaks and produces peak-level signal tracks used downstream |
| 2 | `sample_prep` | after `peakcall` inputs are available | Splits BAMs into NFR and nucleosomal fractions and prepares sample-level tracks |
| 3 | `QTLprep` | after `peakcall` and `sample_prep` | Builds peak SAF/BED files, counts reads, tidies counts, filters features, normalizes counts, and writes QTL BED inputs |
| 4 | `QTLscan` | after `QTLprep` | Computes covariates and runs the main QTL scan for NFR and nucleosomal phenotypes |
| 5 | `coloc_qtls` | after `QTLscan` | Fine-maps QTLs with SuSiE and tests colocalization between NFR and nucleosomal QTL signals |
| 6 | `cit` | after `coloc_qtls` | Runs causal inference testing on colocalized NFR/nucQTL pairs |

### Core Dependency Summary

If you only need the main manuscript-style analysis, the order is:

```text
peakcall
  -> sample_prep
  -> QTLprep
  -> QTLscan
  -> coloc_qtls
  -> cit
```

### What Each Core Step Produces

- `peakcall`
  - calls NFR or nucleosomal peak-related files and signal tracks
  - provides peak definitions used later by `QTLprep`
- `sample_prep`
  - derives NFR and nucleosomal read partitions from sample BAMs
  - generates sample-level BAM and bigWig-style inputs used for counts and visualization
- `QTLprep`
  - converts peaks to SAF/BED
  - counts reads with `featureCounts`
  - tidies and filters raw counts
  - normalizes counts with `bin/normalization.R`
  - writes QTL phenotype BED files
- `QTLscan`
  - computes genotype and phenotype covariates
  - compresses/indexes phenotype BED files
  - runs `QTLtools cis`
  - merges and summarizes QTL scan results
- `coloc_qtls`
  - corrects phenotype BEDs
  - runs SuSiE fine-mapping for NFR and nuc traits
  - performs coloc between NFR and nucleosomal signals
  - writes coloc summary tables and regression summaries
- `cit`
  - reads the coloc summaries
  - filters by `PPH4`
  - runs CIT and FDR summarization

## Optional Workflow Branches

Use these only after the required upstream outputs exist.

### Peak Comparison and Peak-Level Summaries

| Workflow | Run after | Purpose |
| --- | --- | --- |
| `peak_compare` | `peakcall` | Compares MACS2 peak calls to nucleoATAC-derived regions and computes summit distances |

### Variant-in-Peak and Layered Plotting

| Workflow | Run after | Purpose |
| --- | --- | --- |
| `snps_inpeak` | `coloc_qtls` | Extracts and annotates variants within colocalized peak sets |
| `layer_plot` | `snps_inpeak` | Builds motif summaries and layered regional plots for selected GWAS traits |
| `plot_multitrack` | GWAS coloc outputs plus QTL SuSiE outputs | Generates locus-level multi-track plots for selected GWAS peaks |

### eQTL Overlap

| Workflow | Run after | Purpose |
| --- | --- | --- |
| `eqtl_overlap` | `coloc_qtls` recommended | Summarizes eQTL overlap with QTL signals and fits enrichment models |

This workflow relies on QTL overlap outputs and references the NFR SuSiE directory in its logistic modeling step, so it is best treated as downstream of the main QTL discovery and fine-mapping steps.

### Motif Analysis

| Workflow | Run after | Purpose |
| --- | --- | --- |
| `motif_homer` | `coloc_qtls` | HOMER enrichment comparing colocated versus background NFR QTL peaks |
| `motif_meme` | `QTLscan`, `coloc_qtls`, and optionally `cit` | MEME/SEA analyses on QTL peaks, coloc-derived sets, SNP-aware FASTA sets, and CIT-derived sets |

`motif_meme` pulls from several result directories, so the exact prerequisites depend on which target block is active in `rule all`.

### GWAS Colocalization

These workflows are separate branches for different GWAS inputs. Use the one that matches your summary statistics source.

| Workflow | Run after | Purpose |
| --- | --- | --- |
| `coloc_gwas` | `QTLscan` and QTL SuSiE outputs for the selected region | General GWAS colocalization workflow using lead variants and SuSiE |
| `coloc_gwas_ukbb` | before UKBB muscle-specific coloc workflows | Lifts over UKBB GWAS to hg38, identifies lead variants, and fine-maps GWAS loci |
| `coloc_gwas_nonukbb` | `QTLscan` and QTL SuSiE outputs for NFR | Colocalization against a predefined set of non-UKBB GWAS traits |
| `coloc_gwas_muscle` | trait-specific | Muscle disease or muscle-focused GWAS coloc workflow |
| `coloc_gwas_height` | trait-specific | Height GWAS coloc workflow |
| `coloc_gwas_grip_strength` | trait-specific | Grip strength GWAS coloc workflow |
| `coloc_gwas_metrate` | trait-specific | Metabolic-rate GWAS coloc workflow |
| `coloc_gwas_mv` | trait-specific | Multivariate GWAS coloc workflow |

Practical rule:

- run the appropriate `coloc_gwas*` workflow only after the relevant QTL phenotype, QTL scan, and QTL SuSiE outputs already exist
- check the active `rule all` section in the target Snakefile, because several of these workflows have been repointed to different output generations over time

### caQTL Branch

These workflows are separate from the core NFR/nucQTL chain.

| Workflow | Run after | Purpose |
| --- | --- | --- |
| `caQTL` | `peakcall` and sample-level BAM availability | Builds accessibility peaks from MACS2 summits, quantifies them, runs caQTL scans, and fine-maps caQTLs |
| `caQTL_np` | `peakcall` and sample-level BAM availability | Alternate caQTL workflow using MACS2 narrowPeak definitions |
| `ca_partition` | `caQTL` plus NFR/nucQTL peak BEDs | Computes overlap partitions between caQTL, nfrQTL, and nucQTL peak sets |
| `ca_partition_np` | `caQTL_np` plus NFR/nucQTL peak BEDs | Partition analysis for the narrowPeak-based caQTL branch |

Recommended order for the accessibility branch:

```text
peakcall
  -> caQTL or caQTL_np
  -> ca_partition or ca_partition_np
```

### CIT Matching and Downsampling Variants

| Workflow | Run after | Purpose |
| --- | --- | --- |
| `cit_match` | `coloc_qtls` or CIT-related inputs | Matching-based CIT follow-up workflow |
| `cit_downsample` | `cit` inputs plus downsampling setup | Downsampled CIT sensitivity analysis |

These are follow-up or sensitivity-analysis workflows rather than first-pass pipelines.

### MR / Steiger Split Workflow

| Workflow | Run after | Purpose |
| --- | --- | --- |
| `MR_split` | normalized counts, sample metadata, coloc peak list, and genotype resources are prepared | Splits samples, reruns QTL scans on split datasets, harmonizes exposure/outcome effects, and runs Steiger directionality analyses |

`MR_split` is self-contained but expects many explicit config entries such as:

- raw NFR counts
- raw nucleosomal counts
- sample info
- peak list to keep
- VCF and eigenvector files
- output root

Treat this as a specialized downstream analysis, not part of the default run order.

## Minimal Run Recipes

### Reproduce The Main NFR/NucQTL Analysis

```bash
cd run_snakemake
bash dryrun.sh peakcall
bash run.sh peakcall

bash dryrun.sh sample_prep
bash run.sh sample_prep

bash dryrun.sh QTLprep
bash run.sh QTLprep

bash dryrun.sh QTLscan
bash run.sh QTLscan

bash dryrun.sh coloc_qtls
bash run.sh coloc_qtls

bash dryrun.sh cit
bash run.sh cit
```

### Add GWAS Colocalization

After the core QTL workflows are complete, choose the relevant GWAS branch:

```bash
cd run_snakemake
bash dryrun.sh coloc_gwas
bash run.sh coloc_gwas
```

or:

```bash
cd run_snakemake
bash dryrun.sh coloc_gwas_nonukbb
bash run.sh coloc_gwas_nonukbb
```

or one of the other specialized `coloc_gwas*` workflows.

### Add Motif Analysis

After `coloc_qtls` and any required upstream QTL outputs:

```bash
cd run_snakemake
bash dryrun.sh motif_homer
bash run.sh motif_homer

bash dryrun.sh motif_meme
bash run.sh motif_meme
```

## Workflow Selection Guide

Use this quick guide if you are unsure what to run.

- If you want the main discovery pipeline: run `peakcall -> sample_prep -> QTLprep -> QTLscan -> coloc_qtls -> cit`
- If you want NFR versus nucleosomal peak comparison only: run `peakcall`, then `peak_compare`
- If you want GWAS colocalization: finish QTL discovery first, then run the appropriate `coloc_gwas*` workflow
- If you want eQTL overlap: finish QTL discovery and coloc/fine-mapping first, then run `eqtl_overlap`
- If you want motif enrichment on colocated signals: finish `coloc_qtls`, then run `motif_homer` and/or `motif_meme`
- If you want accessibility QTLs from MACS2 summits: run `caQTL` or `caQTL_np`, then the matching partition workflow
- If you want Steiger or split-sample directionality analyses: run `MR_split`

## Nextflow Workflow: LDSC

The `nextflow/ldsc_allgwas_consensus-cluster/` directory contains a separate Nextflow workflow for LDSC annotation enrichment analyses.

Entry point:

- `nextflow/ldsc_allgwas_consensus-cluster/scripts/main.nf`

This workflow:

- builds LDSC annotation files from BED annotations
- calculates LD scores
- runs single-annotation enrichment against GWAS summary statistics
- runs joint annotation enrichment if provided with joint annotation files

This is independent of the main Snakemake execution chain and should be treated as a separate downstream analysis.

## Important Notes

- Many paths in this repository are hard-coded to a specific computing environment.
- Some `rule all` targets have been changed over time and often point to the latest analysis generation only.
- Several Snakefiles contain historical or commented-out targets. Always inspect the active `rule all` section before launching a workflow.
- Some raw data remain controlled access, so not every analysis can be rerun from public inputs alone.

## Citation

If you use this repository or adapt these workflows, please cite:

```text
Wang, X., Robertson, C. C., Varshney, A., Manickam, N., Orchard, P., Laakso, M., Tuomilehto, J., Lakka, T. A., Mohlke, K. L., Boehnke, M., Scott, L. J., Koistinen, H. A., Collins, F. S., & Parker, S. C. J. (2025). Genetic integration with cell-specific nucleosome positioning resolves causal relationships underlying chromatin accessibility profiles (p. 2025.09.09.674883). bioRxiv. https://doi.org/10.1101/2025.09.09.674883

```

## Contact

- Alice Wang: `xiaoouw [at] umich.edu`
- Cassie Robertson: `ccrober [at] umich.edu`
- Stephen Parker: `scjp [at] umich.edu`
