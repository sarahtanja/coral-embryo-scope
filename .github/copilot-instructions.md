# Copilot Instructions for coral-embryo-scope

## Quick Reference

**Language**: R with Quarto (`.qmd`)  
**Primary Libraries**: tidyverse, dplyr, ggplot2  
**Data Format**: CSV files with tidy data structure  
**Working Directory**: `/code` for running scripts  
**Key Datasets**: 
- Raw: `/data/scope_annotation_data/*.csv` (120 files)
- Processed: `/data/output/tidy_bros.csv`, `/data/output/tidy_vials.csv`

## Project Overview

This repository contains code and data for a survival analysis of *Montipora capitata* coral embryos exposed to polyvinyl chloride (PVC) leachate pollution. The project uses R and Quarto for data processing, statistical analysis, and visualization.

## Repository Structure

- **`/code`**: R Quarto (`.qmd`) documents for analysis pipeline
- **`/data`**: Raw and processed data files
  - `/data/metadata`: Experiment metadata
  - `/data/scope_annotation_data`: Individual sample annotation CSV files (120 samples)
  - `/data/output`: Processed tidy datasets (`tidy_bros.csv`, `tidy_vials.csv`)
- **`/images`**: Microscopy images (120 subdirectories, one per sample)
- **`/plots`**: Generated visualizations

## Data Analysis Workflow

The analysis follows a structured pipeline:

1. **`pull_data.qmd`**: Import annotation data from Google Sheets and copy microscopy images
2. **`tidy.qmd`**: Create tidy dataframes from raw annotations
   - Outputs: `tidy_bros.csv` (each row = embryo), `tidy_vials.csv` (each row = sample)
3. **`anova.qmd`**: One-way ANOVA on embryo counts with normality tests
4. **`survival.qmd`**: Survival rate analysis using boxplots, violin plots, GLM, and ANOVA
5. **`abnormality.qmd`**: Abnormality proportion analysis using Dirichlet regression
6. **`timing.qmd`**: Developmental stage proportion analysis using Dirichlet and beta regression
7. **`count_viz.qmd`**: Embryo count data visualizations
8. **`kaplan_meier.qmd`**: Kaplan-Meier survival analysis

## Key Data Structures

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

### Factor Levels (Ordered)
When working with factors, maintain these orderings:
```r
# embryo_phase
bro_levels <- c("egg", "cleavage", "morula", "prawnchip", "earlygastrula")

# treatment
treat_levels <- c("control", "low", "mid", "high")

# status
status_levels <- c("typical", "uncertain", "malformed")
```

## Coding Standards

### R/Quarto Documents
- Use Quarto YAML header with both `gfm` and `html` output formats
- Include standard `knitr::opts_chunk$set()` in setup chunk:
  - `echo = TRUE`, `eval = FALSE`, `warning = FALSE`, `message = FALSE`, `comment = ""`
- Structure documents with clear section headers using `#` markdown syntax

### Tidyverse Style
- Prefer `tidyverse` functions (`dplyr`, `ggplot2`, etc.)
- Use pipe operator `%>%` for chaining operations
- Use `read_csv()` for reading data files
- Use `str()` to inspect data structure after transformations

### Data Loading
Standard paths for data files:
```r
# Metadata
metadata <- read_csv("../data/metadata/scope-metadata.csv")

# Tidy datasets
tidy_bros <- read_csv("../data/output/tidy_bros.csv")
tidy_vials <- read_csv("../data/output/tidy_vials.csv")
```

### Visualization
- Save plots to `../plots/` directory
- Use consistent naming: `{analysis_type}_{plot_type}.png`
- Set resolution: `dpi = 300`
- Use `theme_minimal()` or `theme_journal` for consistency

### Statistical Analysis
- Check normality assumptions before parametric tests
- Use Q-Q plots, density plots, and Shapiro-Wilk tests
- Document non-normal distributions and appropriate alternative tests
- Filter data by `hpf_factor` for time-specific analyses

## Important Notes

- Total of 120 samples (10 crosses × 4 treatments × 3 time points)
- Some samples may have missing or NA values
- 4 hpf data may show binomial distribution (fertilization success vs. failure)
- Always verify factor ordering when creating categorical variables
- Use `list.files()` with `pattern = "\\.csv$"` to read all CSV files from annotation directory
- Images follow same naming convention as data files

## Common Tasks

### Reading Multiple CSV Files
```r
csv <- list.files("../data/scope_annotation_data", pattern = "\\.csv$", full.names = TRUE)
data <- map_dfr(csv, read_csv)
```

### Creating Tidy Dataset
- Join annotation data with metadata using sample identifiers
- Calculate proportions and counts per sample
- Set factor levels with `mutate(factor(..., levels = ..., ordered = TRUE))`

### Filtering by Time Point
```r
data_14hpf <- tidy_vials %>% filter(hpf_factor == "14")
```

## Related Repositories

- **Gene expression analysis**: [coral-embryo-RNAseq](https://github.com/sarahtanja/coral-embryo-RNAseq)
- **Microbiome analysis**: [coral-embryo-microbiome](https://github.com/sarahtanja/coral-embryo-microbiome)

## Development Environment Setup

### Prerequisites
- **R**: Version 4.0 or higher
- **RStudio**: Recommended IDE for working with `.qmd` files
- **Quarto**: CLI tool for rendering Quarto documents
  - Install from: https://quarto.org/docs/get-started/

### Required R Packages
Install the following packages before running analyses:
```r
install.packages(c(
  "tidyverse",    # Data manipulation and visualization
  "dplyr",        # Data manipulation
  "ggplot2",      # Data visualization
  "googlesheets4" # Google Sheets integration (for pull_data.qmd)
))
```

### Workspace Setup
- Open `coral-embryo-scope.Rproj` in RStudio to set the correct working directory
- Quarto documents use relative paths (e.g., `../data/`, `../plots/`)
- Ensure you're in the `/code` directory when running `.qmd` files

## Running/Rendering Quarto Documents

### Rendering Individual Files
```bash
# From the /code directory
quarto render filename.qmd
```

### Running in RStudio
- Open the `.qmd` file in RStudio
- Click "Render" button or use Ctrl+Shift+K (Cmd+Shift+K on Mac)
- Output will generate both `.md` and `.html` files

### Analysis Pipeline Order
Execute documents in this sequence:
1. `pull_data.qmd` - Pull data from Google Sheets (one-time setup)
2. `tidy.qmd` - Create tidy datasets (run when raw data changes)
3. All other `.qmd` files can be run in any order after tidy datasets exist

### Important Notes
- Set `eval = TRUE` in chunk options to actually run code (default is `eval = FALSE`)
- Code chunks are set to not evaluate by default for safety
- Always check that required data files exist before running analyses

## Common Pitfalls and Troubleshooting

### Path Issues
- ❌ **Problem**: "Cannot find file" errors
- ✅ **Solution**: Ensure you're in the `/code` directory when running scripts. All paths are relative to this location.

### Factor Ordering
- ❌ **Problem**: Plots or analyses show incorrect ordering of treatments or stages
- ✅ **Solution**: Always set factor levels explicitly using the ordered levels defined in the coding standards section

### Missing Data
- ❌ **Problem**: Cross 2 data causes errors
- ✅ **Solution**: Cross 2 has been removed from analyses. Use `tidy_bros %>% filter()` to exclude it if needed

### Google Sheets Authentication
- ❌ **Problem**: `pull_data.qmd` fails with authentication errors
- ✅ **Solution**: Run `googlesheets4::gs4_auth()` interactively to authenticate with Google

### Rendering Issues
- ❌ **Problem**: Quarto rendering fails with package errors
- ✅ **Solution**: Check that all required packages are installed. Use `install.packages()` for missing packages.

### Data Loading Errors
- ❌ **Problem**: "File not found" when loading tidy datasets
- ✅ **Solution**: Run `tidy.qmd` first to generate `tidy_bros.csv` and `tidy_vials.csv` in `/data/output`

### Empty or NA Values
- ❌ **Problem**: Unexpected NA values in data
- ✅ **Solution**: This is expected for some samples. The code includes QAQC sections to handle zeros and NAs appropriately.

## File Modification Guidelines

### Files You Should Modify
- ✅ `.qmd` files in `/code` - Analysis scripts and documentation
- ✅ `/plots` - Generated visualizations (created by scripts)
- ✅ `/data/output` - Processed datasets (created by `tidy.qmd`)

### Files You Should NOT Modify
- ❌ `/data/scope_annotation_data/*.csv` - Raw annotation data from Google Sheets
- ❌ `/data/metadata/scope-metadata.csv` - Original metadata
- ❌ `/images` - Original microscopy images
- ❌ `.Rproj.user/` - RStudio project files (in .gitignore)

### Generated Files
The following files are auto-generated and should not be manually edited:
- `.html` files in `/code` (generated from `.qmd` files)
- `.md` files in `/code` (generated from `.qmd` files)  
- Files in `/code/tidy_files/` (Quarto rendering artifacts)

## Getting Help

When working with this repository:
- Check existing `.qmd` files for patterns and examples
- Refer to `douma_weedon_2019_fig1.jpg` in `/code` for regression analysis methods
- Consult README.md for detailed folder descriptions
- Review this file for coding standards and common patterns
