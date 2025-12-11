# Make tidy dataframes

2024-04-01

- [<span class="toc-section-number">1</span> Background](#background)
  - [<span class="toc-section-number">1.1</span> Input](#input)
- [<span class="toc-section-number">2</span> Install & load
  packages](#install--load-packages)
- [<span class="toc-section-number">3</span> tidy_bros : each row is an
  embryo](#tidy_bros--each-row-is-an-embryo)
  - [<span class="toc-section-number">3.1</span> Add
    metadata](#add-metadata)
  - [<span class="toc-section-number">3.2</span> Collapse to
    morula](#collapse-to-morula)
  - [<span class="toc-section-number">3.3</span> Order
    factors](#order-factors)
  - [<span class="toc-section-number">3.4</span> QAQC zeros &
    NAs](#qaqc-zeros--nas)
  - [<span class="toc-section-number">3.5</span> Drop cross
    2](#drop-cross-2)
  - [<span class="toc-section-number">3.6</span> Problems with
    eggs](#problems-with-eggs)
  - [<span class="toc-section-number">3.7</span> Save
    dataframe](#save-dataframe)
- [<span class="toc-section-number">4</span> tidy_vials : each row is a
  sample](#tidy_vials--each-row-is-a-sample)
  - [<span class="toc-section-number">4.1</span> make a column for
    viable embryo counts](#make-a-column-for-viable-embryo-counts)
  - [<span class="toc-section-number">4.2</span> QAQC zeros &
    NAs](#qaqc-zeros--nas-1)
  - [<span class="toc-section-number">4.3</span> Order
    factors](#order-factors-1)
  - [<span class="toc-section-number">4.4</span> Save
    dataframe](#save-dataframe-1)
- [<span class="toc-section-number">5</span> Summary & Next
  Steps](#summary--next-steps)
  - [<span class="toc-section-number">5.1</span> embryos_summary
    :](#embryos_summary-)
  - [<span class="toc-section-number">5.2</span> stage_totals
    :](#stage_totals-)

# Background

Each sample (aka each microscopy slide) resulted in a total embryo
count.

For each sample, each embryo counted was classified by: - stage - status

## Input

- each sample’s own annotation .csv found in
  `data/scope_annotation_data` \## Output
- `data/output/tidy_bros` each row is an embryo
- `data/output/tidy_vials` each row is a sample

# Install & load packages

``` r
library(tidyverse)
```

# tidy_bros : each row is an embryo

There should be 120 csv files, representing each sample, named by the
convention of: - 10 ‘crosses’ of somewhat random parentage **(1-10)** -
4 ‘treatments’ encompassing a control, low, mid, and high pvc leachate
exposure **(C, L, M, H)** - 3 ‘embryonic stages’ or hours
post-fertilization **(4, 9, 14)**

The columns for each file are: - sample_name = col_character() -
embryo_no = col_integer() - embryo_phase = col_factor() - status =
col_factor - notes = col_factor ::: callout-note `list.files()` finds
all .csv files in the directory. `map_dfr(read_csv)` reads each file and
binds rows together (\_dfr = data frame row-bind). :::

``` r
# list all csv files
csv <- list.files("../data/scope_annotation_data", pattern = "\\.csv$", full.names = TRUE)

# Combine all csv's into one tidy df and keep sheet name as row name
tidy_embryos <- csv %>% 
    map_dfr(~ read_csv(.x, col_types = 
                         cols(sample_name = col_character()))) 

# Change sample_name column to sample_id to match metadata column name  
tidy_embryos <- tidy_embryos %>% 
  mutate(sample_id = sample_name,
         stage = embryo_phase) %>%
  select(sample_id, stage, status) 
  

summary(tidy_embryos)
```

      sample_id            stage              status         
     Length:2065        Length:2065        Length:2065       
     Class :character   Class :character   Class :character  
     Mode  :character   Mode  :character   Mode  :character  

``` r
unique(tidy_embryos$stage)
```

    [1] "earlygastrula" "cleavage"      "4to16cell"     "egg"          
    [5] "prawnchip"     "32to64cell"    NA             

``` r
unique(tidy_embryos$status)
```

    [1] "typical"   "malformed" "uncertain" NA         

## Add metadata

``` r
metadata <- read_csv("../data/metadata/scope-metadata.csv")
```

``` r
metadata_clean <- metadata %>%
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>%
  select(
    sample_id,
    parent_1 = parent_colony_a,
    parent_2 = parent_colony_b,
    treatment = leachate,
    hpf,
    date
  )

tidy_embryos <- tidy_embryos %>%
  left_join(metadata_clean, by = "sample_id")

str(tidy_embryos)
```

    tibble [2,065 × 8] (S3: tbl_df/tbl/data.frame)
     $ sample_id: chr [1:2065] "10C14" "10C14" "10C14" "10C14" ...
     $ stage    : chr [1:2065] "earlygastrula" "earlygastrula" "earlygastrula" "earlygastrula" ...
     $ status   : chr [1:2065] "typical" "typical" "typical" "typical" ...
     $ parent_1 : chr [1:2065] "B" "B" "B" "B" ...
     $ parent_2 : chr [1:2065] "J" "J" "J" "J" ...
     $ treatment: chr [1:2065] "control" "control" "control" "control" ...
     $ hpf      : num [1:2065] 14 14 14 14 14 14 14 14 14 14 ...
     $ date     : Date[1:2065], format: "2024-06-07" "2024-06-07" ...

## Collapse to morula

``` r
tidy_bros <- tidy_embryos %>% 
  mutate(stage = if_else(stage %in% c("4to16cell", "32to64cell"), "morula", stage))

unique(tidy_bros$stage)
```

    [1] "earlygastrula" "cleavage"      "morula"        "egg"          
    [5] "prawnchip"     NA             

## Order factors

``` r
# Set your levels in order
# embryo_phase
bro_levels <- c("egg", "cleavage", "morula", "prawnchip", "earlygastrula")

# treatment
treat_levels <- c("control", "low", "mid", "high")

# status
status_levels <- c("typical", "uncertain", "malformed")


# Mutate dataframe to set factor levels
tidy_bros <- tidy_bros %>%
  mutate(stage = factor(stage, 
                        levels = bro_levels,
                        ordered = TRUE),
         treatment = factor(treatment, 
                            levels = treat_levels,
                            ordered = TRUE),
         status = factor(status, 
                         levels = status_levels,
                         ordered = TRUE))

str(tidy_bros)
```

    tibble [2,065 × 8] (S3: tbl_df/tbl/data.frame)
     $ sample_id: chr [1:2065] "10C14" "10C14" "10C14" "10C14" ...
     $ stage    : Ord.factor w/ 5 levels "egg"<"cleavage"<..: 5 5 5 5 5 5 5 5 5 5 ...
     $ status   : Ord.factor w/ 3 levels "typical"<"uncertain"<..: 1 1 1 1 1 1 1 1 1 1 ...
     $ parent_1 : chr [1:2065] "B" "B" "B" "B" ...
     $ parent_2 : chr [1:2065] "J" "J" "J" "J" ...
     $ treatment: Ord.factor w/ 4 levels "control"<"low"<..: 1 1 1 1 1 1 1 1 1 1 ...
     $ hpf      : num [1:2065] 14 14 14 14 14 14 14 14 14 14 ...
     $ date     : Date[1:2065], format: "2024-06-07" "2024-06-07" ...

What colony crosses were used in spawning for microscopy?

``` r
tidy_bros %>% distinct(parent_1, parent_2)
```

    # A tibble: 10 × 2
       parent_1  parent_2
       <chr>     <chr>   
     1 B         J       
     2 H         M       
     3 30-Orange 1-White 
     4 D         N       
     5 M         B       
     6 1-White   J/N     
     7 1-White   M       
     8 1-White   J       
     9 1-White   E       
    10 J         G       

## QAQC zeros & NAs

Which samples had ZERO embryos to count?

``` r
tidy_bros %>% 
  filter(if_any(everything(), is.na))
```

    # A tibble: 9 × 8
      sample_id stage status parent_1  parent_2 treatment   hpf date      
      <chr>     <ord> <ord>  <chr>     <chr>    <ord>     <dbl> <date>    
    1 2C14      <NA>  <NA>   30-Orange 1-White  control      14 2024-05-07
    2 2C9       <NA>  <NA>   30-Orange 1-White  control       9 2024-05-07
    3 2H14      <NA>  <NA>   30-Orange 1-White  high         14 2024-05-07
    4 2H9       <NA>  <NA>   30-Orange 1-White  high          9 2024-05-07
    5 2L14      <NA>  <NA>   30-Orange 1-White  low          14 2024-05-07
    6 2L9       <NA>  <NA>   30-Orange 1-White  low           9 2024-05-07
    7 2M14      <NA>  <NA>   30-Orange 1-White  mid          14 2024-05-07
    8 3L9       <NA>  <NA>   D         N        low           9 2024-05-07
    9 7H9       <NA>  <NA>   H         M        high          9 2024-06-07

> [!WARNING]
>
> Here we have 9 samples that will need some attention… they did not
> contain any (0) viable embryos… only disintegrating embryos or embryo
> debris. How do we handle these zeros? How do we count a ‘zero’? Of
> note cross 2 (from the rack corals; 30-Orange (O) and 1-White (W) )
> had no fertilization and it should be removed entirely from the
> analysis. Note that at 4 hpf cross 2 eggs had not yet dissolved
> despite having zero fertilization! This is good evidence that the
> unfertilized eggs are still visible at 4 hpf.

## Drop cross 2

Filter out (that is exclude) cross 2

``` r
tidy_bros <- tidy_bros %>% 
  filter(!grepl("2", sample_id))
```

There should now be no `NAs` in `status` except for samples 3L9 & 7H9..
which had zero viable embryos! **…But what do we do with that?**

``` r
tidy_bros %>% 
  filter(if_any(everything(), is.na))
```

    # A tibble: 2 × 8
      sample_id stage status parent_1 parent_2 treatment   hpf date      
      <chr>     <ord> <ord>  <chr>    <chr>    <ord>     <dbl> <date>    
    1 3L9       <NA>  <NA>   D        N        low           9 2024-05-07
    2 7H9       <NA>  <NA>   H        M        high          9 2024-06-07

![Sample 3L9 only contained two dissolving debris
fragments](../images/3L9/3L9_4xstitch_anno.jpg) ![Sample 7H9 contained
many dissolving debris fragments, but none appearing
viable](../images/7H9/7H9_manualstitch_anno.jpg)

## Problems with eggs

- At 4 hpf eggs are still intact and visible even if unfertilized
  (e.g. cross 2), however they MAY BE FERTILIZED and developing normally
  but just not yet cleaved. So we cannot assume that all eggs at 4 hpf
  are unfertilized.

- At 9 and 14 hpf unfertilized eggs may have dissolved and no longer be
  visible.

- So are eggs visible at 9 or 14 hpf fertilized but developmentally
  delayed? OR are they unfertilized and just lingering for an unknown
  reason?

- Can you visually tell if an egg is fertilized prior to cellular
  division? Australian Intitute of Marine Science (AIMS) spawning
  research links:

- https://www.facebook.com/australianmarinescience/posts/take-a-look-at-the-development-of-a-coral-embryo-through-the-lens-of-our-microsc/262783479224765/

- Unfertilized eggs typically show no signs of cellular division. If
  they turn ‘fuzzy’, they are not viable and may be decomposing.

## Save dataframe

``` r
write_csv(tidy_bros, "../output/dataframes/tidy_bros.csv")
```

# tidy_vials : each row is a sample

Counts and proportions of combined and distinct stage and status values
in each vial

``` r
# Count + proportion per group
vials <- tidy_bros %>%
  count(sample_id, stage, status, name = "n") %>%
  group_by(sample_id) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = c(stage, status),
    values_from = c(n, prop),
    values_fill = 0
  ) %>% 
  select(!matches("NA"))%>%
  left_join(metadata_clean, by = "sample_id")

str(vials)
```

    tibble [108 × 36] (S3: tbl_df/tbl/data.frame)
     $ sample_id                   : chr [1:108] "10C14" "10C4" "10C9" "10H14" ...
     $ n_earlygastrula_typical     : int [1:108] 20 0 0 8 0 0 12 0 0 19 ...
     $ n_earlygastrula_uncertain   : int [1:108] 1 0 0 3 0 0 0 0 0 0 ...
     $ n_earlygastrula_malformed   : int [1:108] 2 0 0 4 0 0 1 0 0 2 ...
     $ n_egg_typical               : int [1:108] 0 9 0 0 10 0 0 9 0 0 ...
     $ n_cleavage_typical          : int [1:108] 0 12 0 0 6 0 0 7 0 0 ...
     $ n_morula_typical            : int [1:108] 0 8 1 0 4 0 0 10 0 0 ...
     $ n_prawnchip_typical         : int [1:108] 0 0 13 0 0 17 1 0 18 0 ...
     $ n_prawnchip_uncertain       : int [1:108] 0 0 2 0 0 2 0 0 7 0 ...
     $ n_prawnchip_malformed       : int [1:108] 0 0 3 0 0 1 0 0 3 0 ...
     $ n_cleavage_malformed        : int [1:108] 0 0 0 0 1 0 0 1 0 0 ...
     $ n_egg_malformed             : int [1:108] 0 0 0 0 0 0 1 0 0 0 ...
     $ n_morula_malformed          : int [1:108] 0 0 0 0 0 0 0 0 0 0 ...
     $ n_egg_uncertain             : int [1:108] 0 0 0 0 0 0 0 0 0 0 ...
     $ n_cleavage_uncertain        : int [1:108] 0 0 0 0 0 0 0 0 0 0 ...
     $ n_morula_uncertain          : int [1:108] 0 0 0 0 0 0 0 0 0 0 ...
     $ prop_earlygastrula_typical  : num [1:108] 0.87 0 0 0.533 0 ...
     $ prop_earlygastrula_uncertain: num [1:108] 0.0435 0 0 0.2 0 ...
     $ prop_earlygastrula_malformed: num [1:108] 0.087 0 0 0.267 0 ...
     $ prop_egg_typical            : num [1:108] 0 0.31 0 0 0.476 ...
     $ prop_cleavage_typical       : num [1:108] 0 0.414 0 0 0.286 ...
     $ prop_morula_typical         : num [1:108] 0 0.2759 0.0526 0 0.1905 ...
     $ prop_prawnchip_typical      : num [1:108] 0 0 0.684 0 0 ...
     $ prop_prawnchip_uncertain    : num [1:108] 0 0 0.105 0 0 ...
     $ prop_prawnchip_malformed    : num [1:108] 0 0 0.158 0 0 ...
     $ prop_cleavage_malformed     : num [1:108] 0 0 0 0 0.0476 ...
     $ prop_egg_malformed          : num [1:108] 0 0 0 0 0 ...
     $ prop_morula_malformed       : num [1:108] 0 0 0 0 0 0 0 0 0 0 ...
     $ prop_egg_uncertain          : num [1:108] 0 0 0 0 0 0 0 0 0 0 ...
     $ prop_cleavage_uncertain     : num [1:108] 0 0 0 0 0 0 0 0 0 0 ...
     $ prop_morula_uncertain       : num [1:108] 0 0 0 0 0 0 0 0 0 0 ...
     $ parent_1                    : chr [1:108] "B" "B" "B" "B" ...
     $ parent_2                    : chr [1:108] "J" "J" "J" "J" ...
     $ treatment                   : chr [1:108] "control" "control" "control" "high" ...
     $ hpf                         : num [1:108] 14 4 9 14 4 9 14 4 9 14 ...
     $ date                        : Date[1:108], format: "2024-06-07" "2024-06-07" ...

``` r
tidy_vials <- tidy_bros %>%
  group_by(sample_id) %>%
  summarise(
    n_embryos = n(),
    n_typical = sum(status == "typical"),
    n_malformed = sum(status == "malformed"),
    n_uncertain = sum(status == "uncertain"),
    prop_typical = mean(status == "typical"),
    prop_malformed = mean(status == "malformed"),
    prop_uncertain = mean(status == "uncertain"),
    n_egg = sum(stage == "egg"),
    prop_egg = mean(stage == "egg"),
    n_cleavage = sum(stage == "cleavage"),
    prop_cleavage = mean(stage == "cleavage"),
    n_morula = sum(stage == "morula"),
    prop_morula = mean(stage == "morula"),
    n_prawnchip = sum(stage == "prawnchip"),
    prop_prawnchip = mean(stage == "prawnchip"),
    n_earlygastrula = sum(stage == "earlygastrula"),
    prop_earlygastrula = mean(stage == "earlygastrula"),
    .groups = "drop"
  )

str(tidy_vials)
```

    tibble [108 × 18] (S3: tbl_df/tbl/data.frame)
     $ sample_id         : chr [1:108] "10C14" "10C4" "10C9" "10H14" ...
     $ n_embryos         : int [1:108] 23 29 19 15 21 20 15 27 28 21 ...
     $ n_typical         : int [1:108] 20 29 14 8 20 17 13 26 18 19 ...
     $ n_malformed       : int [1:108] 2 0 3 4 1 1 2 1 3 2 ...
     $ n_uncertain       : int [1:108] 1 0 2 3 0 2 0 0 7 0 ...
     $ prop_typical      : num [1:108] 0.87 1 0.737 0.533 0.952 ...
     $ prop_malformed    : num [1:108] 0.087 0 0.1579 0.2667 0.0476 ...
     $ prop_uncertain    : num [1:108] 0.0435 0 0.1053 0.2 0 ...
     $ n_egg             : int [1:108] 0 9 0 0 10 0 1 9 0 0 ...
     $ prop_egg          : num [1:108] 0 0.31 0 0 0.476 ...
     $ n_cleavage        : int [1:108] 0 12 0 0 7 0 0 8 0 0 ...
     $ prop_cleavage     : num [1:108] 0 0.414 0 0 0.333 ...
     $ n_morula          : int [1:108] 0 8 1 0 4 0 0 10 0 0 ...
     $ prop_morula       : num [1:108] 0 0.2759 0.0526 0 0.1905 ...
     $ n_prawnchip       : int [1:108] 0 0 18 0 0 20 1 0 28 0 ...
     $ prop_prawnchip    : num [1:108] 0 0 0.947 0 0 ...
     $ n_earlygastrula   : int [1:108] 23 0 0 15 0 0 13 0 0 21 ...
     $ prop_earlygastrula: num [1:108] 1 0 0 1 0 ...

``` r
tidy_vials <- tidy_vials %>% 
  left_join(vials, by = "sample_id") %>% 
  relocate(c(treatment, hpf, parent_1, parent_2, date), .after = sample_id)
```

## make a column for viable embryo counts

This means malformed = dead, uncertain(torn) = viable. Ecologically, a
malformed embryo isn’t progressing.. but a torn embryo could continue to
develop!

So we’re only counting the total number of embryos that are typical or
uncertain. To do that let’s make new columns called n_viable and
prop_viable for the number of viable embryos and proportion of viable
embryos.

> [!IMPORTANT]
>
> We say that eggs at 4 hpf are viable unless they are malformed.
>
> At 9 and 14 hpf we say that all eggs are non-viable
> (i.e. unfertilized) so we subtract those from the viable count as
> well.

``` r
tidy_vials <- tidy_vials %>%
  mutate(
    n_viable = case_when(
      hpf == 4 ~ n_embryos - (n_egg_malformed + 
                              n_cleavage_malformed + 
                              n_morula_malformed +
                              n_prawnchip_malformed +
                              n_earlygastrula_malformed),
      hpf %in% c(9, 14) ~ n_embryos - (n_egg + 
                                        n_egg_malformed + 
                                        n_cleavage_malformed + 
                                        n_morula_malformed +
                                        n_prawnchip_malformed +
                                        n_earlygastrula_malformed)
    )
  ) %>%
  relocate(n_viable, .before = n_embryos) %>%
  mutate(prop_viable = n_viable / n_embryos) %>%
  relocate(prop_viable, .before = n_embryos)
```

## QAQC zeros & NAs

``` r
tidy_vials %>% 
  filter(if_any(everything(), is.na))
```

    # A tibble: 2 × 55
      sample_id treatment   hpf parent_1 parent_2 date       n_viable prop_viable
      <chr>     <chr>     <dbl> <chr>    <chr>    <date>        <int>       <dbl>
    1 3L9       low           9 D        N        2024-05-07       NA          NA
    2 7H9       high          9 H        M        2024-06-07       NA          NA
    # ℹ 47 more variables: n_embryos <int>, n_typical <int>, n_malformed <int>,
    #   n_uncertain <int>, prop_typical <dbl>, prop_malformed <dbl>,
    #   prop_uncertain <dbl>, n_egg <int>, prop_egg <dbl>, n_cleavage <int>,
    #   prop_cleavage <dbl>, n_morula <int>, prop_morula <dbl>, n_prawnchip <int>,
    #   prop_prawnchip <dbl>, n_earlygastrula <int>, prop_earlygastrula <dbl>,
    #   n_earlygastrula_typical <int>, n_earlygastrula_uncertain <int>,
    #   n_earlygastrula_malformed <int>, n_egg_typical <int>, …

3L9 had zero embryos & 7H9 had zero embryos These samples show ‘1’
embryo, but the row was just a placeholder and really there were no
viable embryos found in these sample (see above). Here we simply replace
the ‘1’ with a ‘0’ and all the NA values in the columns with ‘0’.

``` r
tidy_vials <- tidy_vials %>%
  mutate(across(
    c(n_viable, prop_viable, n_embryos),
    ~if_else(sample_id %in% c("3L9", "7H9"), 0, .)
  ))
```

Should no longer be NA… n_viable and prop_viable and n_embryos should be
0

``` r
tidy_vials %>% 
  filter(sample_id == c("3L9", "7H9"))
```

    # A tibble: 2 × 55
      sample_id treatment   hpf parent_1 parent_2 date       n_viable prop_viable
      <chr>     <chr>     <dbl> <chr>    <chr>    <date>        <dbl>       <dbl>
    1 3L9       low           9 D        N        2024-05-07        0           0
    2 7H9       high          9 H        M        2024-06-07        0           0
    # ℹ 47 more variables: n_embryos <dbl>, n_typical <int>, n_malformed <int>,
    #   n_uncertain <int>, prop_typical <dbl>, prop_malformed <dbl>,
    #   prop_uncertain <dbl>, n_egg <int>, prop_egg <dbl>, n_cleavage <int>,
    #   prop_cleavage <dbl>, n_morula <int>, prop_morula <dbl>, n_prawnchip <int>,
    #   prop_prawnchip <dbl>, n_earlygastrula <int>, prop_earlygastrula <dbl>,
    #   n_earlygastrula_typical <int>, n_earlygastrula_uncertain <int>,
    #   n_earlygastrula_malformed <int>, n_egg_typical <int>, …

## Order factors

Sometimes for plotting purposes we need hpf to be ordered as a factor
and not a numerical data

``` r
tidy_vials <- tidy_vials %>% 
  mutate(hpf = factor(hpf, levels = c("4","9","14"), ordered = TRUE ),
         treatment = factor(treatment, levels = c("control", "low", "mid", "high"), ordered = TRUE),
         )
```

## Save dataframe

``` r
write_csv(tidy_vials, "../output/dataframes/tidy_vials.csv")
```

# Summary & Next Steps

Ok! We have two tidy dataframes in long format to work with that should
allow us to answer all our questions.

The following are thought-experiment dataframes… not fully flushed out
yet

## embryos_summary :

``` r
# Count embryos per (sample_id, stage, status)
embryos_summary <- tidy_bros %>%
  count(sample_id, stage, status, name = "n_embryos") %>%
  group_by(sample_id) %>%
  mutate(proportion = n_embryos / sum(n_embryos)) %>%
  ungroup()
```

## stage_totals :

``` r
# Total counts per stage per sample_id
stage_totals <- tidy_bros %>%
  count(sample_id, stage, name = "total_n")
```

``` r
sessionInfo()
```

    R version 4.2.3 (2023-03-15)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Ubuntu 24.04.3 LTS

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so

    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
     [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] lubridate_1.9.4 forcats_1.0.0   stringr_1.6.0   dplyr_1.1.4    
     [5] purrr_1.2.0     readr_2.1.5     tidyr_1.3.1     tibble_3.3.0   
     [9] ggplot2_4.0.1   tidyverse_2.0.0

    loaded via a namespace (and not attached):
     [1] compiler_4.2.3     pillar_1.11.1      RColorBrewer_1.1-3 tools_4.2.3       
     [5] bit_4.6.0          digest_0.6.37      timechange_0.3.0   jsonlite_2.0.0    
     [9] evaluate_1.0.5     lifecycle_1.0.4    gtable_0.3.6       pkgconfig_2.0.3   
    [13] rlang_1.1.6        cli_3.6.5          rstudioapi_0.17.1  parallel_4.2.3    
    [17] yaml_2.3.10        xfun_0.54          fastmap_1.2.0      withr_3.0.2       
    [21] knitr_1.50         generics_0.1.4     vctrs_0.6.5        hms_1.1.3         
    [25] bit64_4.6.0-1      grid_4.2.3         tidyselect_1.2.1   glue_1.8.0        
    [29] R6_2.6.1           vroom_1.6.5        rmarkdown_2.29     farver_2.1.2      
    [33] tzdb_0.5.0         magrittr_2.0.4     scales_1.4.0       htmltools_0.5.8.1 
    [37] S7_0.2.1           utf8_1.2.6         stringi_1.8.7      crayon_1.5.3      
