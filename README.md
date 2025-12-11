# coral-embryo-scope

> The egg–sperm bundles released by *M. capitata* measured approximately 1 mm and contained around 15 ± 5.1 oocytes (mean ± SD, n = 214, from 26 colonies), surrounding a central mass of spermatozoa (Fig. 4). — Padilla-Gamino et al., 2011


A repo for the microscopy imaging annotations and survival analysis of *Montipora capitata* embryos exposed to increasing levels of PVC leachate pollution.

This repository contains code and output for a **survival** analysis of coral embryos exposed to 3 increasing levels of polyvinyl chloride (PVC) leachate. *Montipora capitata* rice coral embryos were subjected to a single-stressor ecotoxicology assay designed to test the impact of PVC leachate on their early development. You can check out the complimentary repositories for a **gene expression** [coral-embryo-RNAseq](https://github.com/sarahtanja/coral-embryo-RNAseq) and **microbiome** (coral-embryo-microbiome) [coral-embryo-microbiome](https://github.com/sarahtanja/coral-embryo-microbiome) analysis from the same ecotoxicological assay.

## Repository Structure

```mermaid
graph TD
    A[coral-embryo-scope] --> B[code/]
    A --> C[data/]
    A --> D[scope-images/]
    A --> E[output/]
    A --> F[annotation-guide/]
    A --> G[svg/]
    
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
    C --> C2[scope_annotation_data/]
    
    C1 --> C1A[scope-metadata.csv]
    C2 --> C2A[120 sample CSV files]
    
    D --> D1[Sample folders<br/>e.g., 1C4, 1C9, 1C14, etc.]
    D1 --> D1A[Annotated microscopy images]
    
    E --> E1[dataframes/]
    E --> E2[figs/]
    E --> E3[tables/]
    
    E1 --> E1A[tidy_bros.csv<br/>tidy_vials.csv<br/>supertidy_bros.csv]
    E2 --> E2A[Generated plot images]
    E3 --> E3A[Generated tables]
    
    F --> F1[Image annotation examples<br/>by developmental stage]
    
    G --> G1[SVG/PNG graphics<br/>of developmental stages]
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
    
    Anova --> Plots1[output/figs/<br/>ANOVA results]
    Survival --> Plots2[output/figs/<br/>Survival boxplots]
    Abnormality --> Plots3[output/figs/<br/>Status proportions]
    Timing --> Plots4[output/figs/<br/>Stage proportions]
    CountViz --> Plots5[output/figs/<br/>Count visualizations]
    KaplanMeier --> Plots6[output/figs/<br/>Survival curves]
    
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
Contains raw data files organized into two subdirectories:

#### `/data/metadata`
- **`scope-metadata.csv`**: Metadata for microscopy samples including cross IDs, parent information, treatments, and time points

#### `/data/scope_annotation_data`
Contains 120 CSV files (one per sample) with individual embryo annotations. Each file is named using the convention:
- Cross ID (1-10)
- Treatment (C=control, L=low, M=mid, H=high)
- Hours post-fertilization (4, 9, 14)
- Example: `1C4.csv`, `2L9.csv`, `10H14.csv`

### `/scope-images`
Contains 120 subdirectories (one per sample), each holding annotated microscopy images captured on a Nikon DS-Fi 3 camera and annotated using Nikon NIS Elements BR 4.6.00 64-bit software. Subdirectories follow the same naming convention as the data files (e.g., `1C4/`, `2L9/`, `10H14/`).

### `/output`
Contains all processed data and analysis outputs organized into three subdirectories:

#### `/output/dataframes`
Processed tidy datasets ready for analysis:
- **`tidy_bros.csv`**: Tidy dataframe where each row represents an individual embryo with its stage and status classifications
- **`tidy_vials.csv`**: Tidy dataframe where each row represents a sample (microscopy slide) with counts and proportions of embryos by stage and status
- **`supertidy_bros.csv`**: Additional processed embryo-level data
- **`tidy_status.csv`**: Status-level summary data
- **`tidy_timing.csv`**: Timing/stage-level summary data
- **`prop_summary.csv`**: Proportion summary statistics
- **`status_summary.csv`**: Status summary statistics

#### `/output/figs`
Contains generated visualization outputs from the analysis scripts:
- **`embryo_survival_box.png`**: Boxplot of embryo survival counts
- **`embryo_stage_stackedbar.png`**: Stacked bar chart of developmental stage proportions
- **`embryo_status_stackedbar.png`**: Stacked bar chart of embryo status proportions
- **`figure_survival_timing_morphology.png`**: Combined figure showing survival, timing, and morphology

#### `/output/tables`
Contains generated tables and statistical outputs:
- **`table_stage_composition.html`**: HTML table of stage composition
- **`table_survival_stargazer.html`**: HTML survival analysis table
- **`.doc`** files: Word-compatible table outputs

### `/annotation-guide`
Contains reference images and guides for annotating embryo microscopy images. Includes example images organized by developmental stage (egg, cleavage, morula, prawnchip, early gastrula) and status (typical, malformed). Also contains Quarto documents for annotation protocols and scope statistics.

### `/svg`
Contains graphical representations (SVG/PNG format) of coral embryo developmental stages:
- **`egg.png`**: Egg stage illustration
- **`cleavage.png`**: Cleavage stage illustration
- **`morula.png`**: Morula stage illustration
- **`prawnchip.png`**: Prawnchip stage illustration
- **`earlygastrula.png`**: Early gastrula stage illustration
