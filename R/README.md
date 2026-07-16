# R Scripts Folder

This folder contains the raw R scripts for each step of the pleiotropy analysis workflow tutorial. These scripts can be run directly in R or RStudio.

## Analysis Scripts

* **`00_setup.R`**: Environment setup and package installation.
* **`01_prepare_data.R`**: Combine summary statistics into `prepared_summary_statistics.rds`.
* **`02_asset.R`**: Subset-based meta-analysis using the `ASSET` package.
* **`03_placo.R`**: Composite null hypothesis testing using the `PLACO` package.
* **`04_gpa.R`**: Genetic Pleiotropy Analysis with pathway annotation using the `GPA` package.
* **`05_cpbayes.R`**: SNP-level Bayesian cross-phenotype analysis using the `CPBayes` package.
* **`06_gcpbayes.R`**: Gene or pathway-level Bayesian cross-phenotype analysis using the `GCPBayes` package on the *BCAT1* gene as an example.
* **`07_run_all.R`**: Pipeline runner to execute steps 01–06 and compile the HTML report.

## Running the Scripts

To run method by method individually:
```R
source("R/02_asset.R")
```

Alternatively, to run the complete workflow and compile all the results in one go:
```R
source("R/07_run_all.R")
```
