# coral-embryo-scope

A repo for the microscopy imaging annotations and survival analysis of *Montipora capitata* embryos exposed to increasing levels of PVC leachate pollution.

This repository contains code and output for a **survival** analysis of coral embryos exposed to 3 increasing levels of polyvinyl chloride (PVC) leachate. *Montipora capitata* rice coral embryos were subjected to a single-stressor ecotoxicology assay designed to test the impact of PVC leachate on their early development. You can check out the complimentary repositories for a **gene expression** [coral-embryo-RNAseq](https://github.com/sarahtanja/coral-embryo-RNAseq) and **microbiome** (coral-embryo-microbiome) [coral-embryo-microbiome](https://github.com/sarahtanja/coral-embryo-microbiome) analysis from the same ecotoxicological assay.

## Repository Structure

```mermaid
graph TD
    A[coral-embryo-scope] --> B[code/]
    A --> C[data/]
    A --> D[images/]
    A --> E[plots/]
    
    B --> B1[pull_data.qmd]
    B --> B2[tidy.qmd]
    B --> B3[anova.qmd]
    B --> B4[survival.qmd]
    B --> B5[abnormality.qmd]
    B --> B6[timing.qmd]
    B --> B7[count_viz.qmd]
    B --> B8[kaplan_meier.qmd]
    B --> B9[douma_weedon_2019_fig1.jpg]
    
    C --> C1[metadata/]
    C --> C2[output/]
    C --> C3[scope_annotation_data/]
    
    C1 --> C1A[scope-metadata.csv]
    C2 --> C2A[tidy_bros.csv]
    C2 --> C2B[tidy_vials.csv]
    C2 --> C2C[supertidy_bros.csv]
    C3 --> C3A[120 sample CSV files]
    
    D --> D1[Sample folders<br/>e.g., 1C4, 1C9, 1C14, etc.]
    D1 --> D1A[Annotated microscopy images]
    
    E --> E1[Generated plot images]
```

## Analysis Workflow

The analysis follows a structured pipeline where data is pulled, tidied, and then analyzed through multiple complementary approaches:

```mermaid
flowchart TB
    Start([Raw Data]) --> PullData[pull_data.qmd<br/>Import from Google Sheets]
    PullData --> Tidy[tidy.qmd<br/>Create tidy datasets]
    
    Tidy --> TidyBros[(tidy_bros.csv<br/>Each row = embryo)]
    Tidy --> TidyVials[(tidy_vials.csv<br/>Each row = sample)]
    
    TidyVials --> Anova[anova.qmd<br/>One-way ANOVA<br/>Embryo counts]
    TidyVials --> Survival[survival.qmd<br/>Survival rate analysis]
    TidyVials --> Abnormality[abnormality.qmd<br/>Abnormality proportions]
    TidyVials --> Timing[timing.qmd<br/>Developmental stages]
    TidyVials --> CountViz[count_viz.qmd<br/>Count visualizations]
    
    TidyBros --> KaplanMeier[kaplan_meier.qmd<br/>Survival analysis]
    
    Anova --> Plots1[Plots: ANOVA results]
    Survival --> Plots2[Plots: Survival boxplots]
    Abnormality --> Plots3[Plots: Status proportions]
    Timing --> Plots4[Plots: Stage proportions]
    CountViz --> Plots5[Plots: Count visualizations]
    KaplanMeier --> Plots6[Plots: Survival curves]
    
    style PullData fill:#e1f5ff
    style Tidy fill:#e1f5ff
    style TidyBros fill:#ffe1e1
    style TidyVials fill:#ffe1e1
    style Anova fill:#fff4e1
    style Survival fill:#fff4e1
    style Abnormality fill:#fff4e1
    style Timing fill:#fff4e1
    style CountViz fill:#fff4e1
    style KaplanMeier fill:#fff4e1
```

## Folder Descriptions

### `/code`
Contains all R Quarto (`.qmd`) documents for data processing, analysis, and visualization:

- **`pull_data.qmd`**: Imports annotation data from Google Sheets and copies microscopy images from local desktop to the repository
- **`tidy.qmd`**: Creates tidy dataframes from raw annotation data. Outputs `tidy_bros.csv` (each row is an embryo) and `tidy_vials.csv` (each row is a sample)
- **`anova.qmd`**: Performs one-way ANOVA on embryo counts, tests for normality, and compares viable embryo counts across treatments
- **`survival.qmd`**: Analyzes survival rates of embryos across different treatments and time points using boxplots, violin plots, GLM, and ANOVA
- **`abnormality.qmd`**: Analyzes abnormality proportions (typical, uncertain, malformed) using Dirichlet regression
- **`timing.qmd`**: Analyzes developmental stage proportions (egg, cleavage, morula, prawnchip, gastrula) using Dirichlet and beta regression
- **`count_viz.qmd`**: Visualizes embryo count data
- **`kaplan_meier.qmd`**: Performs Kaplan-Meier survival analysis
- **`douma_weedon_2019_fig1.jpg`**: Reference figure for regression analysis methods

### `/data`
Contains all data files organized into three subdirectories:

#### `/data/metadata`
- **`scope-metadata.csv`**: Metadata for microscopy samples including cross IDs, parent information, treatments, and time points

#### `/data/output`
Processed tidy datasets ready for analysis:
- **`tidy_bros.csv`**: Tidy dataframe where each row represents an individual embryo with its stage and status classifications
- **`tidy_vials.csv`**: Tidy dataframe where each row represents a sample (microscopy slide) with counts and proportions of embryos by stage and status
- **`supertidy_bros.csv`**: Additional processed embryo-level data

#### `/data/scope_annotation_data`
Contains 120 CSV files (one per sample) with individual embryo annotations. Each file is named using the convention:
- Cross ID (1-10)
- Treatment (C=control, L=low, M=mid, H=high)
- Hours post-fertilization (4, 9, 14)
- Example: `1C4.csv`, `2L9.csv`, `10H14.csv`

### `/images`
Contains 120 subdirectories (one per sample), each holding annotated microscopy images captured on a Nikon DS-Fi 3 camera and annotated using Nikon NIS Elements BR 4.6.00 64-bit software. Subdirectories follow the same naming convention as the data files (e.g., `1C4/`, `2L9/`, `10H14/`).

### `/plots`
Contains generated visualization outputs from the analysis scripts:
- **`counts_survival_boxplot.png`**: Boxplot of embryo survival counts
- **`counts_viable_boxplot.png`**: Boxplot of viable embryo counts
- **`density_viable.png`**: Density plot of viable embryos
- **`proportion_stage_stackedbar.png`**: Stacked bar chart of developmental stage proportions
- **`proportion_stagexstatus_stackedbar.png`**: Stacked bar chart showing stage by status proportions
- **`proportion_status_stackedbar.png`**: Stacked bar chart of embryo status proportions
- **`viablecounts_survival_boxplot.png`**: Boxplot of viable embryo survival counts
