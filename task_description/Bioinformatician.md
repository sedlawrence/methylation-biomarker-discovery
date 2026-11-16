# Bioinformatician Take-Home Task

## Dataset

[GSE59685](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59685) (publicly available via GEO)

This dataset contains 450K DNA methylation array data from 531 samples across multiple tissues 
(entorhinal cortex, frontal cortex, superior temporal gyrus, cerebellum, and whole blood) from 122 donors. 
Samples include both Alzheimer's disease cases and controls, with Braak stage recorded as a 
continuous measure of neuropathology severity.

## Task

Using the whole blood samples only from GSE59685, perform a differential methylation analysis against cases vs controls 
by converting Braak stage to a binary variable. Please use Python or R for this task. 

### Step 1 — Data Loading and QC

- Download and load the data (note the compressed beta file is about 2GB)
- Perform QC: flag and remove poor quality samples and probes, with clear justification for each exclusion decision
- Filter non-variable and cross-reactive probes and other necessary pre-processing steps

### Step 2 — Normalisation

- Apply an appropriate normalisation method and justify your choice
- Assess and adjust for cell-type composition if appropriate

### Step 3 — Differential Methylation Analysis

- Identify differentially methylated positions (DMPs) using the appropriate libraries
- Apply appropriate covariate adjustment (age, sex, cell type estimates, etc)
- Report results with effect sizes, standard errors, and corrected p-values

### Step 4 — Output

- A ranked table of top DMPs with genomic annotation (gene name, CpG island context)
- Visualisation (e.g. volcano plot, top CpG beta distribution)
- A brief methods summary covering all analytical decisions and their justifications

**Feel free to use AI tools or not for this task. 
If you do, please explain what you used them for, how you used them, and what you were able to achieve with them.**