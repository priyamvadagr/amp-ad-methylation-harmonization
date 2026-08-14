# functions

Shared R functions sourced by the driver scripts in `../Data_processing`, `../Cortical_clock`, and `../Misc`. Not run directly.

- **cohort_config.R** — Canonical per-cohort file locations (`raw_dir`, `manifest`, `combined_dir`, `individual_csv`, `qc_dir`, `normalized_dir`) shared across the pipeline so paths can't drift between stages. Also defines `select_cohorts()`, which resolves command-line arguments (built-in cohort names, or ad-hoc `NAME:...` specs) into cohort configs.

- **build_idat_manifest.R** — `build_idat_manifest()`: maps each specimen to its green/red IDAT file, producing one row per specimenID/individualID with the matching Grn/Red file paths.

- **harmonize_metadata.R** — `harmonize_adkp()`: generic harmonizer that maps a study's raw individual metadata onto the ADKP harmonized data dictionary via a per-study mapping table (source column, type, value maps, constants, derived fields). Includes `validate_harmonized()` to check conformance and `combine_studies()` to stack harmonized cohorts. Based on [syn73713784](https://www.synapse.org/Synapse:syn73713784)

- **read_cohort_idats.R** — `read_cohort_idats()`: reads a cohort's raw IDAT pairs into one minfi object (`RGChannelSet`/`RGChannelSetExtended`), optionally restricted to a sample sheet, skipping corrupt files (bad magic header) and reading in memory-bounded chunks to avoid OOM kills on large cohorts.

- **sample_qc.R** — Sample-level QC: `ewastools_qc()` (control-probe QC — bisulfite conversion + 17-metric sample failure — from raw idats, chunked), `sample_qc()` (minfi detection-p/intensity/sex/bead-count QC, merged with ewastools flags; flags only, drops nothing), and `render_qc_plots()` (diagnostic figures + a control-probe report restricted to flagged samples). Also defines `beadcount_from_rgset()` for per-CpG bead-count QC.

- **noob_normalize.R** — `normalize_noob()`: memory-conscious noob (within-sample) normalization of a QC'd `RGChannelSet` — drops flagged samples, masks high-detection-p probe values, and writes the `GenomicMethylSet`, beta matrix, and (optionally) M-values, freeing large objects as soon as each is no longer needed.

- **filter_probes.R** — Zhou-mask probe filtering per Lundin et al. ([2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10996836/#s0002:~:text=PMCID%3A%20PMC10996836%C2%A0%C2%A0PMID%3A-,38571307,-ABSTRACT) & [2025](https://pubmed.ncbi.nlm.nih.gov/40480020/#:~:text=DOI%3A-,10.1016/j.atherosclerosis.2025.120219,-Abstract)): `drop_low_callrate()`, `drop_crossreactive()`, `drop_snp()`, `drop_nonautosomal()`, each composable individually, plus `filter_probes()` which chains all four and returns a per-step summary. Supports both 450K and EPIC platforms via separate Zhou InfiniumAnnotation mapping/SNP mask tables.

- **pca_qc.R** — `pca_qc()`: PCA on the most-variable, complete-case M-values of a filtered beta matrix, testing each PC for association (R²-primary, p-value secondary) with technical and biological covariates to detect batch effects before modelling. Writes scores, variance-explained, association tables, and elbow/scatter/heatmap figures.

- **cortical_clock.R** — Helpers for the Shireby CorticalClock predictor: `clock_cpgs_from_coef()` / `subset_clock_cpgs()` to read the clock's 347 CpGs and subset a beta matrix to them, `plot_clock_missingness()` for pre-clock coverage diagnostics (per-sample/per-CpG missingness, detection-p heatmap, UpSet co-occurrence), and `plot_clock_by_tissue()` for per-tissue predicted-vs-chronological-age accuracy plots.
