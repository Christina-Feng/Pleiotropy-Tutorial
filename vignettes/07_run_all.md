# Step 07: Running the Entire Workflow

This script runs the entire pleiotropy tutorial analysis pipeline in sequence, reproducing all analysis steps and generating a comprehensive HTML workflow report.

## Step-by-Step Analysis & Results

### 1. Render Workflow Report
Render the companion R Markdown report (`vignettes/pleiotropy_workflow.Rmd`), compiling all analysis results, figures, and outputs:

```R
rm(list = ls())

# Load setup first to check directories and packages
source("R/00_setup.R")

if (requireNamespace("rmarkdown", quietly = TRUE)) {
  cat("Rendering the workflow R Markdown report...\n")
  rmarkdown::render(
    input = "vignettes/pleiotropy_workflow.Rmd",
    output_dir = "results",
    output_file = "pleiotropy_workflow.html",
    quiet = FALSE
  )
  cat("Workflow completed successfully. HTML report saved to results/pleiotropy_workflow.html\n")
} else {
  cat("rmarkdown package not found. Running the workflow sequentially...\n")
  source("R/01_prepare_data.R")
  source("R/02_asset.R")
  source("R/03_placo.R")
  source("R/04_gpa.R")
  source("R/05_cpbayes.R")
  source("R/06_gcpbayes.R")
  cat("Workflow completed successfully.\n")
}
```

**Output:**
```text
Rendering the workflow R Markdown report...

processing file: pleiotropy_workflow.Rmd
  |.......................                               |  40% [prepare-data]   
  |................................................      |  86% [gcpbayes]       
  |......................................................| 100% [session-info]   

output file: pleiotropy_workflow.knit.md

/Applications/Pandoc/pandoc +RTS -K512m -RTS pleiotropy_workflow.knit.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output pleiotropy_workflow.html --embed-resources --standalone --toc --toc-depth 3 --number-sections --mathjax 

Output created: results/pleiotropy_workflow.html

Workflow completed successfully. HTML report saved to results/pleiotropy_workflow.html
```

You can view the rendered HTML report here: [pleiotropy_workflow.html](../results/pleiotropy_workflow.html)

## Description & Key Tasks

1. **Environment Checking**: Calls `00_setup.R` to verify required CRAN and Bioconductor packages.
2. **Reproducible Report Compilation**: Knits `vignettes/pleiotropy_workflow.Rmd` using the `rmarkdown` library, executing each step of the pipeline (data prep, ASSET, PLACO, GPA, CPBayes, and GCPBayes) and documenting the outputs in a single HTML file.
3. **Sequential Fallback**: If the `rmarkdown` package is not available locally, it gracefully falls back to executing the raw `.R` scripts sequentially.
4. **Compile RDS Datasets**: Regardless of the rendering mode, all analysis outputs are compiled and written as `.rds` dataset files in the `results/` folder.

## Running this script

You can execute this script from the root directory of the project in R/RStudio using:

```R
source("R/07_run_all.R")
```
