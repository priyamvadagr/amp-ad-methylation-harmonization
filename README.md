# AMP-AD Methylation Harmonization

Harmonization project for Illumina DNA methylation array data across AMP-AD cohorts — **MSBB**, **ROSMAP**, **ROSMAP_APOE4**, and **MOA-PAD**. Takes each cohort from raw IDATs through a common QC, normalization, and probe-filtering process, and applies the Shireby CorticalClock DNAm-age predictor for downstream age-acceleration analysis.

## Pipeline

```
raw IDATs -> compiled RGChannelSet -> sample QC -> noob normalization -> probe filtering -> PCA QC
                                                  \-> CorticalClock age prediction (pre-filter beta)
```

1. **Metadata reconciliation** (`Scripts/Misc/get_sample_info.Rmd`) — merges assay/biospecimen/individual metadata per cohort, harmonizes to the ADKP data dictionary where needed, cross-checks metadata against available IDAT files, and builds the specimen↔IDAT manifest.
2. **Compile IDATs** (`Scripts/Data_processing/compile_idat_files.R`) — reads raw `.idat` pairs into one minfi `RGChannelSet` per cohort.
3. **Sample QC** (`Scripts/Data_processing/sample_level_QC.R`) — `ewastools` control-probe QC + minfi detection-p/intensity/sex/bead-count QC. Flags only; nothing is dropped at this stage.
4. **Normalization** (`Scripts/Data_processing/normalize_data.R`) — noob (within-sample) normalization, dropping flagged samples and masking high detection-p probe values.
5. **Probe filtering** (`Scripts/Data_processing/filter_probes_per_cohort.R`) — Zhou-mask filtering (low call rate, cross-reactive, SNP artefact, non-autosomal).
6. **PCA QC** (`Scripts/Data_processing/pca_qc_per_cohort.R`) — PCA on filtered beta values, testing PCs against technical and biological covariates to catch batch effects.
7. **CorticalClock** (`Scripts/Cortical_clock/`) — DNAm-age prediction on the pre-probe-filter noob beta matrix, with per-tissue accuracy analysis and age-acceleration comparisons.

Each stage's outputs and file-naming conventions are documented in the `README.md` within its `Scripts/` subdirectory.

## Repository layout

- `Scripts/functions/` — shared R functions sourced by the driver scripts below; not run directly.
- `Scripts/Data_processing/` — main pipeline drivers (steps 2–6 above).
- `Scripts/Misc/` — metadata reconciliation and manifest building.
- `Scripts/Cortical_clock/` — CorticalClock age prediction and analysis.
- `Results/` — QC figures, PCA diagnostics, and clock accuracy plots/statistics per cohort.
- `Logs/`, `*.csv`, `*.rds` — pipeline logs and intermediate/processed data artifacts; gitignored, not tracked here.

Cohort-specific file paths (raw data, manifests, QC/normalized output directories) are centralized in `Scripts/functions/cohort_config.R` — this repo contains no raw or processed methylation data itself.

## Setup

R dependencies are managed with [`renv`](https://rstudio.github.io/renv/); `renv::restore()` will install the pinned package versions from `renv.lock`. Running the CorticalClock scripts additionally requires the [CorticalClock](https://github.com/gemmashireby/CorticalClock) repository to be downloaded separately.
Zhou mask tables can be downloaded from:
[HM450.hg19.manifest.pop.tsv.gz](https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotationData/main/Anno/HM450/HM450.hg19.manifest.pop.tsv.gz) and [EPIC.hg19.manifest.pop.tsv.gz](https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotationData/main/Anno/EPIC/EPIC.hg19.manifest.pop.tsv.gz)