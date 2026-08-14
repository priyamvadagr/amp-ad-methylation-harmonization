# Cortical_clock

Runs the Shireby et al. CorticalClock DNAm-age predictor ([CorticalClock github](https://github.com/gemmashireby/CorticalClock)) on each cohort's pre-probe-filter noob beta matrix, and analyzes the resulting age predictions. Shared helper functions live in [`../functions/cortical_clock.R`](../functions/cortical_clock.R); cohort paths come from `../functions/cohort_config.R`. Need to have CorticalClock downloaded.

- **check_cpG_coverage.R** — For each cohort, checks how many of the 347 clock CpGs are present/covered in the pre-filter noob beta matrix before running the clock, and plots per-sample/per-CpG missingness plus detection-p and masking co-occurrence. Writes coverage CSVs and a `_clock_missingness.pdf` per cohort.

- **run_cortical_clock.R** — Runs `CorticalClock()` function from Shireby et al. on each cohort's pre-filter noob beta matrix (probe filtering may drop clock CpGs), isolating each cohort's fixed-name outputs in its own directory and renaming them to `<cohort>_*`. Also produces per-tissue accuracy plots, since a clock trained on cortex can be systematically off in non-cortical tissue. Writes `<cohort>_CorticalPred.csv`, `_accuracy_statistics.csv`, `_CorticalClockplot.pdf`, and per-tissue accuracy PDFs.

- **compare_ROSMAP_age_acc.R** — Ad-hoc analysis script: merges ROSMAP's clock age-acceleration output with individual metadata, derives AD vs. control diagnosis labels from Braak/CERAD/cogdx (and a pathology-only variant harmonizable across cohorts), and plots/tests age-acceleration differences between groups (Wilcoxon + Welch t-test, violin/box/jitter plot).
