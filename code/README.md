

# Code Directory

This directory contains all R Quarto (`.qmd`) documents for data processing, statistical analysis, and visualization of *Montipora capitata* coral embryo survival/timing/abnormality based on annotations from scope images.

## Data Structure

### Sample Naming Convention
Samples follow the format: `{cross}{treatment}{hpf}`
- **Cross**: 1-10 (parent colony crosses)
- **Treatment**: C (control), L (low), M (mid), H (high) PVC leachate exposure
- **HPF**: 4, 9, 14 (hours post-fertilization)
- Example: `1C4`, `2L9`, `10H14`

### Embryo Classification
Each embryo is classified by:
- **Stage**: egg, cleavage, morula, prawnchip, earlygastrula
- **Status**: typical, uncertain, malformed


## File Descriptions

- **`01_pull_data.qmd`**: Imports annotation data from Google Sheets 
- **`02_tidy_data.qmd`**: Creates tidy dataframes from raw annotation data
  - makes: 
  `tidy_bros.csv` (each row = individual embryo)
  `tidy_vials.csv` (each row = sample from a single experimental unit - 20mL scintillation vial)
- **`03_survival.qmd`**: Survival rate analysis using negative binomial GLM, and ANOVA
- **`04_timing.qmd`**: Analyzes developmental stage counts using negative binomial GLM
- **`05_abnormality.qmd`**: Analyzes abnormality counts using using negative binomial GLM
- **`06_figure.qmd`**: Creates a figure with `patchwork` and `cowplot` showing survival timing and abnormality simultaneously

   

### Statistical Analysis Methods

The analysis employs two main statistical approaches depending on the type of data:

#### Count Data Analysis (ANOVA)

**Files**: `anova.qmd`, `survival.qmd`, `count_viz.qmd`, `kaplan_meier.qmd`

**Method**: These analyses use **count data** (number of viable embryos) and apply **one-way ANOVA** followed by post-hoc tests to compare embryo survival across treatment groups.

- **`anova.qmd`**: One-way ANOVA on embryo counts with normality tests (Shapiro-Wilk) and homogeneity of variance checks
- **`survival.qmd`**: Survival rate analysis using boxplots, violin plots, GLM, and ANOVA on count data
- **`count_viz.qmd`**: Visualizations of embryo count data
- **`kaplan_meier.qmd`**: Kaplan-Meier survival analysis using time-to-event data

**Data type**: Discrete counts of embryos (e.g., number of viable embryos per sample)

#### Proportional Data Analysis (Dirichlet Regression)

**Method**: These analyses use **proportional data** and apply **Dirichlet regression** and **beta regression** to model continuous proportions that sum to 1.

- **`timing.qmd`**: Analyzes developmental stage proportions (egg, cleavage, morula, prawnchip, gastrula) using Dirichlet regression and beta regression
- **`abnormality.qmd`**: Analyzes abnormality status proportions (typical, uncertain, malformed) using Dirichlet regression

**Data type**: Continuous proportions representing the fraction of embryos in each category (e.g., proportion reaching each developmental stage)

**Why Dirichlet regression?** Unlike ANOVA, Dirichlet regression is specifically designed for compositional data where multiple proportions must sum to 1. This is appropriate for analyzing how embryos are distributed across multiple developmental stages or status categories.

**Reference**: Douma, J. C., & Weedon, J. T. (2019). Analysing continuous proportions in ecology and evolution: A practical introduction to beta and Dirichlet regression. *Methods in Ecology and Evolution*, 10(9), 1412-1430. ([see `douma_weedon_2019_fig1.jpg`](douma_weedon_2019_fig1.jpg))



### trial_n_error
- **`count_viz.qmd`**: Visualizes embryo count data
- **`kaplan_meier.qmd`**: Kaplan-Meier like survival analysis
- **`douma_weedon_2019_fig1.jpg`**: Reference figure for analysis methods




