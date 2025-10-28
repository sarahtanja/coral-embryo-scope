# Code Directory

This directory contains all R Quarto (`.qmd`) documents for data processing, statistical analysis, and visualization of *Montipora capitata* coral embryo survival data.

## Analysis Workflow Overview

The analysis pipeline follows a structured approach where data is imported, tidied, and then analyzed through multiple complementary statistical methods. The workflow is designed to assess the impact of PVC leachate exposure on embryo survival, development, and abnormality rates.

### Data Processing Pipeline

1. **`pull_data.qmd`**: Import annotation data from Google Sheets and copy microscopy images
2. **`tidy.qmd`**: Create tidy dataframes from raw annotations
   - Outputs: `tidy_bros.csv` (each row = individual embryo), `tidy_vials.csv` (each row = sample)

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

**Files**: `timing.qmd`, `abnormality.qmd`

**Method**: These analyses use **proportional data** and apply **Dirichlet regression** and **beta regression** to model continuous proportions that sum to 1.

- **`timing.qmd`**: Analyzes developmental stage proportions (egg, cleavage, morula, prawnchip, gastrula) using Dirichlet regression and beta regression
- **`abnormality.qmd`**: Analyzes abnormality status proportions (typical, uncertain, malformed) using Dirichlet regression

**Data type**: Continuous proportions representing the fraction of embryos in each category (e.g., proportion reaching each developmental stage)

**Why Dirichlet regression?** Unlike ANOVA, Dirichlet regression is specifically designed for compositional data where multiple proportions must sum to 1. This is appropriate for analyzing how embryos are distributed across multiple developmental stages or status categories.

**Reference**: Douma, J. C., & Weedon, J. T. (2019). Analysing continuous proportions in ecology and evolution: A practical introduction to beta and Dirichlet regression. *Methods in Ecology and Evolution*, 10(9), 1412-1430. ([see `douma_weedon_2019_fig1.jpg`](douma_weedon_2019_fig1.jpg))

## File Descriptions

- **`pull_data.qmd`**: Imports annotation data from Google Sheets and copies microscopy images
- **`tidy.qmd`**: Creates tidy dataframes from raw annotation data
- **`anova.qmd`**: One-way ANOVA on embryo counts with normality and variance tests
- **`survival.qmd`**: Survival rate analysis using count data, boxplots, GLM, and ANOVA
- **`abnormality.qmd`**: Analyzes abnormality proportions using Dirichlet regression
- **`timing.qmd`**: Analyzes developmental stage proportions using Dirichlet and beta regression
- **`count_viz.qmd`**: Visualizes embryo count data
- **`kaplan_meier.qmd`**: Kaplan-Meier survival analysis
- **`douma_weedon_2019_fig1.jpg`**: Reference figure for regression analysis methods

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

## Running the Analysis

All `.qmd` files can be rendered using Quarto:

```bash
# From the /code directory
quarto render filename.qmd
```

Or in RStudio by clicking "Render" or using Ctrl+Shift+K (Cmd+Shift+K on Mac).

**Note**: Run `pull_data.qmd` and `tidy.qmd` first to generate the required tidy datasets before running the analysis scripts.

---

> The egg–sperm bundles released by M. capitata measured approximately 1 mm and contained around 15 ± 5.1 oocytes (mean ± SD, n = 214, from 26 colonies), surrounding a central mass of spermatozoa (Fig. 4). — Padilla-Gamino et al., 2011
