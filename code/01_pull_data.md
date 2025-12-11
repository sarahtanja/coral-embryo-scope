# Download .csv Data from Google Sheet
Sarah Tanja
2025-01-20

- [<span class="toc-section-number">1</span> Background](#background)
  - [<span class="toc-section-number">1.1</span> Inputs](#inputs)
  - [<span class="toc-section-number">1.2</span> Outputs](#outputs)
- [<span class="toc-section-number">2</span> Install & load
  packages](#install--load-packages)
- [<span class="toc-section-number">3</span> Import dataframes to a
  list](#import-dataframes-to-a-list)
- [<span class="toc-section-number">4</span> Save each sheet as
  .csv](#save-each-sheet-as-csv)
- [<span class="toc-section-number">5</span> Copy in all images from
  local desktop](#copy-in-all-images-from-local-desktop)
- [<span class="toc-section-number">6</span> Summary & Next
  Steps](#summary--next-steps)

------------------------------------------------------------------------

# Background

In August 2024, we made 120 permanent slides of *Montipora capitata*
embryos at 3 developmental stages: (4 hours post-fertilization;
cleavage) (9 hours post-fertilization; prawn chip) (14 hours
post-fertilization; early gastrula)

These embryos were exposed to increasing levels of polyvinyl chloride
(PVC) leachate: 0 mg/L (control group: 0.22$\mu m$ filtered seawater
from Kaneohe Bay, Oahu) 0.01 mg/L (low level exposure group, nominally
equivalent to 0.01 mg of PVC microplastics (\<500$\mu m$ in particle
size) per liter of seawater) 0.1 mg/L (mid level exposure group) 1 mg/L
(high level exposure group)

The embryos were fixed in 4% Z-fix, mounted with glycerol, and
photographed at 40X using a Nikon DS-Fi 3 camera. Embryos in images were
then annotated manually via Nikon’s NIS Elements BR 4.6.00 64 bit
software.

Here we pull in the images and annotation data files for analysis in R.

## Inputs

- annotation data files populated in a google sheet
- annotated images originally taken on Nikon DS-Fi 3 camera
- annotations added via Nikon’s NIS Elements BR 4.6.00 64 bit software

## Outputs

- data sheets pulled from Google Sheet to git repo in
  `data/scope_annotation_data`
- original annotated images archived in repo’s `images` folder

Here we are pulling our data from a Google Sheet, with a single ‘sheet’
containing data for a single sample.

# Install & load packages

``` r
library(tidyverse)
library(dplyr)
library(stringr)
library(ggimage)
library(googlesheets4)
```

# Import dataframes to a list

``` r
# Authenticate your Google account
gs4_auth()
```

``` r
# This is the Google Sheet URL
sheet_url <- "https://docs.google.com/spreadsheets/d/1Gj7aiv3p68cmMo6j-T_U_r57rQo8CoIWyrbhQE7Jt0g/edit?gid=1617870343#gid=1617870343"
```

Download all sheets in the workbook. Each sheet is data from a single
microscopy slide.

There are 120 slides (i.e. samples).

All of these data frames will be imported into a list named
`sheets_list` .

``` r
# Get all sheet names (these are the sample_id)
all_sheets <- sheet_names(sheet_url)

# Read each sheet into a named list of data frames
sheets_list <- lapply(all_sheets, function(sheet) {
  read_sheet(sheet_url, sheet = sheet)
})

names(sheets_list) <- all_sheets
```

# Save each sheet as .csv

``` r
# Make sure the folder exists (creates it if it doesn't)
if (!dir.exists("../data/scope_annotation_data")) dir.create("../data/scope_annotation_data", recursive = TRUE)

# Save each sheet as a CSV in the ../data folder
invisible(lapply(all_sheets, function(sheet) {
  write.csv(
    sheets_list[[sheet]],
    file = file.path("../data/scope_annotation_data", paste0(sheet, ".csv")),
    row.names = FALSE
  )
})
)
```

# Copy in all images from local desktop

-a: archive mode (preserves structure, permissions, etc.)

-v: verbose (shows you what’s happening)

–include ‘\*/’: always include folders (so structure is kept)

–include ‘\*.jpg’: include jpg files

–exclude ‘\*’: exclude everything else

/path/to/source/: your starting directory (note the trailing /)

/path/to/destination/: where you want the .jpg files and folder
structure copied

``` bash
rsync -av --include '*/' --include '*.jpg' --exclude '*' "/c/Users/Minerva/Desktop/coral-embryo-leachate-scope-pics/" "/c/Users/Minerva/Documents/GitHub/coral-embryo-scope/images/"
```

# Summary & Next Steps

We now have annotated scope images located in `images` and one .csv file
per sample (n=120 samples, or 120 microscopy slides) in
`data/scope_annotation_data` .

Next steps are to compile the scope annotation data in `code/tidy.qmd`
and visualize it for initial trends in `code/counts.qmd` and
`code/proportions.qmd`
