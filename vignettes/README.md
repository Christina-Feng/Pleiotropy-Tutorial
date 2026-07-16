# Vignettes Folder

This folder contains the step-by-step markdown walkthroughs and the primary R Markdown template for the tutorial workflow.

## Tutorial Walkthroughs

* **`00_setup.md`**: Environment setup and package verification.
* **`01_prepare_data.md`**: Summary statistics loading and alignment.
* **`02_asset.md`**: Subset-based meta-analysis using`ASSET` package.
* **`03_placo.md`**: Composite null hypothesis testing using `PLACO`.
* **`04_gpa.md`**: Genetic Pleiotropy Analysis with annotations using `GPA`.
* **`05_cpbayes.md`**: SNP-level Bayesian cross-phenotype analysis using `CPBayes`.
* **`06_gcpbayes.md`**: Gene-level Bayesian cross-phenotype analysis using `GCPBayes` on the *BCAT1* gene.
* **`07_run_all.md`**: Executing the sequential pipeline and compiling results.

---

## R Markdown Template

* **`pleiotropy_workflow.Rmd`**: The master R Markdown document that imports scripts from the `R/` folder to compile steps 01–06 into a single HTML report (`results/pleiotropy_workflow.html`).
