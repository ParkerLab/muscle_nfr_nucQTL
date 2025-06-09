# Single Nucleus Muscle NFR and Nucleosome QTL Analysis

This repository contains the code and data analysis pipelines for the paper:

> **Wang X** **Robertson C.C. et al.**,.*
> **bioRxiv** (2025).

## Overview

This project investigates genetic variation affecting chromatin accessibility and nucleosome positioning in human skeletal muscle using snATAC-seq data. We identify:

* **nfrQTLs** (nucleosome-free region QTLs)
* **nucQTLs** (nucleosomal region QTLs)

We also analyze their relationship to regulatory annotations, fine-mapped GWAS variants, eQTLs, and transcription factor binding.

## Repository Structure

```
muscle_nfr_nucQTL/
├── bin/                    # All R and Python based codes we used
├── env/                    # yml files for the conda environments
├── nextflow/               # LDSC analysis from Varshney et al. 2024 paper
├── run_snakemake/          # Bash commands that help you run the snakemake files
├── smk/                    # Snakemake files that manage each step of our analysis
└── README.md               # Project overview and setup instructions
```

## Requirements

This repository was developed using:

* R (≥ 4.2.3)
* Python (≥ 3.8)
* Snakemake (≥ 7.30.1)
* Nextflow (≥ 24.10.6)
* Common bioinformatics tools:

  * [QTLtools] (https://qtltools.github.io/qtltools/pages/basic_install.html)
  * [Samtools](https://www.htslib.org/)
  * [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)
  * [MACS2](https://github.com/macs3-project/MACS)
  * [UCSC.utils] (https://genome.axolotl-omics.org/util.html)
  * [deepTools] (https://deeptools.readthedocs.io/en/develop/index.html)

Additional dependencies can be found in the respective `bin/` and `smk/` folders.

## Usage

### Running a Snakemake Pipeline

To execute a Snakemake workflow, open and use the `run.sh` script from the `run_snakemake` directory:

```bash
bash run.sh <pipeline_name>
```

Where `<pipeline_name>` corresponds to one of the workflows in the `smk/` directory. For example:

```bash
bash run.sh peakcall
```

This will run the pipeline defined in `smk/peak_calling.smk`.

Make sure you have configured paths and resources in the corresponding `.yaml` config file. You can also test your scripts by running:
```bash
bash dryrun.sh <pipeline_name>
```

### Main Pipelines

* `peakcall` – Call NFR and nucleosomal peaks from snATAC-seq BAM files
* `sample_prep` – Partition the files for anlaysis
* `QTLprep` – Prepare files for the QTL scan
* `QTLscan` – nfrQTL and nucQTL scanning
* `coloc_qtls` – Colocalization analysis between the nfrQTLs and the nucQTLs
* `cit` – Causal inference analysis between nfrQTLs and nucQTLs
* `coloc_gwas*` – Colocalization analysis between the QTLs and the GWAS summary statistics
* `motif_meme` – Motif enrichment analysis for different stages

## Notes
Some of the raw data used in this study are currently under controlled access. We are making our best effort to make all processed data and summary statistics publicly available in the near future.

Certain directories and file paths in the code may include date information. Please adjust paths and config settings based on your own computing environment.

## Citation

If you use this code or refer to the analyses in this repository, please cite:

```
Wang X and Robertson C.C et al., ...
```

## Contact

For questions or feedback, please contact:
**Alice Wang** – xiaoouw \[at] umich.edu
**Cassie Robertson** - ccrober \[at] umich.edu
**Stephen Parker** – scjp \[at] umich.edu

---

