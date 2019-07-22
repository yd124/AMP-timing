---
title: "RV217 - Data Cleaning"
author: "Bryan Mayer"
date: "7/22/2019"
output: 
  html_document:
    keep_md: true
    toc: true
---




The following script cleans and tidies the raw RV217 data.

**Some notes:**

  - `days` was recalculated to be from the time of the first positive viral load. 
    - Original days variable is now called `days_dx` for days from diagnosis (APTIMA  detection).
  - Used Dan's variable name conventions so our datasets have similar headers
  - No long removing missing IDs, updated processing to handle this better
  - Added a variable for arv. arv day is relative to first positive viral load.

**Some issues:**

- Not every participant has site information
- Not sure what many of the variables mean
- APTIMA needs an LLoQ? Currently: 
    - NR, R or R/NR = NA
- Date variable was used as needed but not generally cleaned for IDs without issues and without ARV (`draw_date`)

# Load data and clean names


```r
library(tidyverse)
```

```
## Registered S3 methods overwritten by 'ggplot2':
##   method         from 
##   [.quosures     rlang
##   c.quosures     rlang
##   print.quosures rlang
```

```
## ── Attaching packages ─────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──
```

```
## ✔ ggplot2 3.1.1     ✔ purrr   0.3.2
## ✔ tibble  2.1.3     ✔ dplyr   0.8.1
## ✔ tidyr   0.8.3     ✔ stringr 1.4.0
## ✔ readr   1.3.1     ✔ forcats 0.4.0
```

```
## ── Conflicts ────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
library(kableExtra)
```

```
## 
## Attaching package: 'kableExtra'
```

```
## The following object is masked from 'package:dplyr':
## 
##     group_rows
```

```r
fix_names = function(x){
  #this loads the header exactly
  suppressMessages({
    x = read_csv("../../data/RV217Master.csv", n_max = 0) %>%
      names()
    })    
  
  x %>%
    str_replace_all("#", "") %>%
    str_replace_all("%", "pct") %>%
    tidy_names(syntactic = T, quiet = T) 
}


# read.csv loads the data than read_csv
rv217_raw = read.csv("../../data/RV217Master.csv", stringsAsFactors = F,
                 na.strings = "") %>%
  as_tibble(.name_repair = fix_names) %>%
  #using Dan's naming conventions
  rename(
    ID = PIN,
    draw_date = Draw.Date,
    VL = Abbott.RealTime.HIV.1.RNA..Copies.ml.,
    log10VL = Abbott.RealTime.HIV.1.RNA.Log..Copies.ml.,
    VL_site = VL.run.by.MHRP..M..or.Site..S.,
    visit_code = Visit
  ) %>%
  mutate(
    VLunit = "copies/mL"
  )

glimpse(rv217_raw)
```

```
## Observations: 1,489
## Variables: 23
## $ Notes             <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ Priority          <chr> "Priority 1", NA, NA, NA, NA, NA, NA, NA, NA, …
## $ Site              <chr> "Ug", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
## $ ID                <int> 10066, 10066, 10066, 10066, 10066, 10066, 1006…
## $ visit_code        <chr> "A", "B", "D", "F", "H", "SBV", "SBV", "1", "2…
## $ draw_date         <chr> "4-Feb-10", "18-Feb-10", "22-Jul-10", "6-Jan-1…
## $ days.from.1st.pos <int> -588, -574, -420, -252, -84, -21, 0, 1, 4, 8, …
## $ APTIMA            <chr> NA, NA, NA, NA, NA, "NR", "12.38", "15.12", "1…
## $ log10VL           <chr> NA, NA, NA, NA, NA, NA, NA, "4.76", "5.83", "6…
## $ VL                <int> NA, NA, NA, NA, NA, NA, NA, 57544, 676083, 181…
## $ VL_site           <chr> NA, NA, NA, NA, NA, NA, NA, "M", "M", "M", "M"…
## $ EIA               <chr> "NR", "NR", "NR", "NR", "NR", NA, NA, "NR", NA…
## $ WB                <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "I…
## $ Fiebig.Stage      <chr> NA, NA, NA, NA, NA, NA, NA, "I-II", "I-II", "I…
## $ CD4.pct           <dbl> NA, NA, NA, NA, NA, NA, NA, 47, NA, NA, NA, 25…
## $ CD8.pct           <dbl> NA, NA, NA, NA, NA, NA, NA, 25, NA, NA, NA, 59…
## $ NK.pct            <dbl> NA, NA, NA, NA, NA, NA, NA, 13, NA, NA, NA, 7,…
## $ B.pct             <dbl> NA, NA, NA, NA, NA, NA, NA, 13, NA, NA, NA, 8,…
## $ CD4.              <dbl> NA, NA, NA, NA, NA, NA, NA, 824, NA, NA, NA, 3…
## $ CD8.              <dbl> NA, NA, NA, NA, NA, NA, NA, 442, NA, NA, NA, 9…
## $ NK.               <dbl> NA, NA, NA, NA, NA, NA, NA, 212, NA, NA, NA, 1…
## $ B.                <dbl> NA, NA, NA, NA, NA, NA, NA, 212, NA, NA, NA, 1…
## $ VLunit            <chr> "copies/mL", "copies/mL", "copies/mL", "copies…
```

# Error checking

## Dealing with missing IDS

Missing PIN seems to come from populated notes columns

- no cases of viral load or cell data


```r
missing_id = rv217_raw %>%
  subset(is.na(ID))

nrow(missing_id)
```

```
## [1] 18
```

```r
select(missing_id, contains("VL"), contains("CD"), contains("B.")) %>%
  na.omit %>%
  nrow()
```

```
## [1] 0
```

- There are **three** rows without any data so removing them.


```r
full_missing_rows = map(1:nrow(rv217_raw), function(i){
  if(all(is.na(select(rv217_raw, -VLunit)[i, ]))) i
}) %>%
  flatten_dbl()

full_missing_rows
```

```
## [1] 843 844 863
```

```r
x = nrow(rv217_raw)
x
```

```
## [1] 1489
```

```r
rv217_raw = rv217_raw[-c(full_missing_rows), ]
nrow(rv217_raw)
```

```
## [1] 1486
```

```r
if(x - nrow(rv217_raw) != 3) stop("total missing rows has changed") else rm(x)

# confirm
map(1:nrow(rv217_raw), function(i){
  if(all(is.na(select(rv217_raw, -VLunit)[i, ]))) i
}) %>%
  flatten_dbl()
```

```
## numeric(0)
```

- Identifying which missing IDs are within a series of the same ID

If they are, we can just use locf (assuming master is sorted), otherwise it has to be dealt with manually. There were three cases and in each of these cases, the relevant note/visit code aligns with using locf.


```r
tmp_missing_id = map_df(which(is.na(rv217_raw$ID)), function(i){
  rows = c(i - 1, i, i + 1)
  tibble(
    missing_row = i,
    check_rows = str_c(rows, collapse = ","),
    check_id = paste(map_dbl(rows, ~(rv217_raw$ID[.x])), collapse = ","),
    same = str_split_fixed(check_id, ",", n = 3)[,1] == str_split_fixed(check_id, ",", n = 3)[,3] 
  )
})

tmp_missing_id
```

```
## # A tibble: 15 x 4
##    missing_row check_rows     check_id       same 
##          <int> <chr>          <chr>          <lgl>
##  1         350 349,350,351    40100,NA,40100 TRUE 
##  2         427 426,427,428    40168,NA,40168 TRUE 
##  3         519 518,519,520    40257,NA,40265 FALSE
##  4         672 671,672,673    40436,NA,40436 TRUE 
##  5         798 797,798,799    40577,NA,40577 TRUE 
##  6         813 812,813,814    10753,NA,10753 TRUE 
##  7         817 816,817,818    10753,NA,10753 TRUE 
##  8         885 884,885,886    40503,NA,40503 TRUE 
##  9         930 929,930,931    40737,NA,40737 TRUE 
## 10        1000 999,1000,1001  10374,NA,10502 FALSE
## 11        1089 1088,1089,1090 40134,NA,40134 TRUE 
## 12        1168 1167,1168,1169 40283,NA,40283 TRUE 
## 13        1222 1221,1222,1223 40320,NA,40435 FALSE
## 14        1362 1361,1362,1363 40652,NA,40652 TRUE 
## 15        1387 1386,1387,1388 40700,NA,40700 TRUE
```

```r
tmp_missing_id %>%
  subset(!same) %>%
  select(check_rows) %>%
  unlist() %>%
  paste(collapse = ",") %>%
  str_split(",") %>%
  unlist() %>%
  as.numeric %>%
  slice(rownames_to_column(rv217_raw), ., .preserve = T) %>%
  select(rowname, ID, draw_date, Notes, visit_code, days.from.1st.pos)
```

```
## # A tibble: 9 x 6
##   rowname    ID draw_date Notes              visit_code   days.from.1st.pos
##   <chr>   <int> <chr>     <chr>              <chr>                    <int>
## 1 518     40257 15-Jun-11 <NA>               17                         439
## 2 519        NA 1-Aug-12  ARV started 8-1-12 <NA>                       852
## 3 520     40265 25-Feb-10 <NA>               A                         -265
## 4 999     10374 21-Feb-12 <NA>               17                         407
## 5 1000       NA <NA>      died 3-8-12        <NA>                        NA
## 6 1001    10502 14-Jun-11 <NA>               A                         -448
## 7 1221    40320 17-Nov-10 <NA>               B                           NA
## 8 1222       NA <NA>      <NA>               discontinued                NA
## 9 1223    40435 8-Feb-11  <NA>               A                         -624
```


```r
rv217_raw = rv217_raw %>%
  mutate(ID = zoo::na.locf(ID))
  

rv217_raw %>%
  subset(is.na(ID)) %>%
  nrow()
```

```
## [1] 0
```

## Variable types

Checking types of our important endpoints (viral load and cell counts). Random character string in log10VL


```r
rv217_raw %>%
  select(contains("VL"), contains("CD"), contains("B.")) %>%
  summarize_all(class) %>%
  gather()
```

```
## # A tibble: 10 x 2
##    key     value    
##    <chr>   <chr>    
##  1 log10VL character
##  2 VL      integer  
##  3 VL_site character
##  4 VLunit  character
##  5 CD4.pct numeric  
##  6 CD8.pct numeric  
##  7 CD4.    numeric  
##  8 CD8.    numeric  
##  9 B.pct   numeric  
## 10 B.      numeric
```

```r
map(rv217_raw$log10VL, function(i){
  if(!is.na(i) & is.na(as.numeric(i))) i
}) %>%
  flatten_chr()
```

```
## Warning in .f(.x[[i]], ...): NAs introduced by coercion
```

```
## [1] "Aptima Reactive"
```

```r
rv217_raw$log10VL[which(rv217_raw$log10VL == "Aptima Reactive")] = NA
rv217_raw$log10VL = as.numeric(rv217_raw$log10VL)

rv217_raw %>%
  select(contains("VL"), contains("CD"), contains("B.")) %>%
  summarize_all(class) %>%
  gather()
```

```
## # A tibble: 10 x 2
##    key     value    
##    <chr>   <chr>    
##  1 log10VL numeric  
##  2 VL      integer  
##  3 VL_site character
##  4 VLunit  character
##  5 CD4.pct numeric  
##  6 CD8.pct numeric  
##  7 CD4.    numeric  
##  8 CD8.    numeric  
##  9 B.pct   numeric  
## 10 B.      numeric
```

# Variable munging

##  Priority, Site

There is one unique value per participant. Not sure what prioirty is but site is important so will clean that into Uganda or Thailand unless it is missing.


```r
priority_site_dat = rv217_raw %>% group_by(ID) %>%
  summarize(
    total_priority = n_distinct(na.omit(Priority)),
    priority = paste0(na.omit(Priority), collapse = ""),
    total_site = n_distinct(na.omit(Site)),
    site_raw = paste0(na.omit(Site), collapse = "")
  )

priority_site_dat %>%
  subset(total_priority > 1 | total_site > 1) %>%
  nrow()
```

```
## [1] 0
```

```r
unique(priority_site_dat$priority)
```

```
## [1] "Priority 1" "Priority 4" "Priority 3" ""           "RV398"
```

```r
unique(priority_site_dat$site_raw)
```

```
## [1] "Ug"       "Uganda"   "UG"       "TH"       ""         "Thailand"
```

```r
priority_site_dat = priority_site_dat %>%
  mutate(site = case_when(
    str_detect(tolower(site_raw), "th") ~ "Thailand",
    str_detect(tolower(site_raw), "ug") ~ "Uganda",
    TRUE ~ NA_character_
  ))

unique(priority_site_dat$site)
```

```
## [1] "Uganda"   "Thailand" NA
```

```r
with(priority_site_dat, table(site_raw, site, useNA = "ifany"))
```

```
##           site
## site_raw   Thailand Uganda <NA>
##                   0      0    8
##   TH             25      0    0
##   Thailand        9      0    0
##   Ug              0      8    0
##   UG              0      4    0
##   Uganda          0      1    0
```

```r
priority_site_dat %>%
  subset(is.na(site)) %>%
  select(ID, site) 
```

```
## # A tibble: 8 x 2
##      ID site 
##   <int> <chr>
## 1 40067 <NA> 
## 2 40320 <NA> 
## 3 40814 <NA> 
## 4 41002 <NA> 
## 5 41074 <NA> 
## 6 41086 <NA> 
## 7 41136 <NA> 
## 8 41146 <NA>
```



```r
rv217_raw = rv217_raw %>%
  select(-Priority, -Site) %>%
  left_join(select(priority_site_dat, ID, priority, site), by = "ID") %>%
  select(-visit_code, -Notes, -EIA, -WB, everything()) # moves misc. variables to the end
```

## APTIMA cleaning

Make a variable that is measurement only truncated at minimum (for R and R/NR), truncated at 0 for NR and left missing for weird inputs


```r
aptima_vals = map_chr(rv217_raw$APTIMA, function(i){
  if(!is.na(i) & is.na(suppressWarnings(as.numeric(i)))) return(i) else return(NA_character_)
}) 

aptima_vals %>% unique()
```

```
## [1] NA                      "NR"                    "29,33"                
## [4] "25,27"                 "R"                     "R/NR"                 
## [7] "IND"                   "no follow up 4 months"
```

```r
rv217_raw = rv217_raw %>%
  mutate(
    APTIMA_num = suppressWarnings(as.numeric(APTIMA))
    )
```

# days variable

## Days variable error checking

There are several cases with missing day variables. Can these be calculated?


```r
missing_time_all = rv217_raw %>%
  group_by(ID) %>%
  mutate(
    any_missing_day = any(is.na(days.from.1st.pos))
  ) %>%
  subset(any_missing_day)

# missing times but has date
date_missing_time = missing_time_all %>%
  group_by(ID) %>%
  mutate(
   vl_missing_time = any(is.na(days.from.1st.pos) & !is.na(draw_date))
  ) %>%
  subset(vl_missing_time)

ids_missing_days = unique(date_missing_time$ID)

date_missing_time %>%
  group_by(ID) %>%
  summarize(
    any(days.from.1st.pos == 0)
    )
```

```
## # A tibble: 5 x 2
##      ID `any(days.from.1st.pos == 0)`
##   <int> <lgl>                        
## 1 40320 NA                           
## 2 40436 TRUE                         
## 3 40491 TRUE                         
## 4 40737 TRUE                         
## 5 41146 TRUE
```


```r
# recode a weird entry
date_missing_time$draw_date[date_missing_time$draw_date == "NO SAMPLES"] <- NA

clean_vl_missing_time = date_missing_time %>%
  group_by(ID) %>%
  mutate(
    infection_time = if(any(!is.na(days.from.1st.pos))) na.omit(draw_date[days.from.1st.pos == 0]) else{
                             na.omit(draw_date[APTIMA == "R"])},
    days.from.1st.pos = as.numeric(difftime(as.Date(draw_date, format = "%d-%b-%y"),
                                 as.Date(infection_time, format = "%d-%b-%y"), 
                                 units = "days"))
  ) %>%
  select(-infection_time, -any_missing_day, -vl_missing_time)
```

## Days variable creation


```r
nrow(rv217_raw)
```

```
## [1] 1486
```

```r
rv217_raw = rv217_raw %>%
  subset(!ID %in% ids_missing_days) %>% # these are processed separately
  bind_rows(clean_vl_missing_time) %>% # and added back here
  rename(days_dx = days.from.1st.pos) %>%
  group_by(ID) %>%
  mutate(
    posVL = !is.na(VL),
    days = days_dx - min(days_dx[posVL])
  ) 

nrow(rv217_raw)
```

```
## [1] 1486
```

# ARV and Primary Kinetics


## ARV

Defining primary kinetics as first pos VL -> min(first ARV, 365 days).

All ARV information is in the Notes. There were three cases where the ARV note text date did not match the ARV note date. I found this through manual inspection and denoted them in the outputted sheet.



```r
ARV_notes = rv217_raw %>% 
  subset(str_detect(Notes, "ARV")) %>%
  select(ID, draw_date, days, Notes) %>%
  add_count() %>% 
  group_by(ID) %>%
  mutate(first_note_day = min(days)) %>%
  subset(days == first_note_day) %>%
  rename(arv_note = Notes) %>%
  select(ID, draw_date, days, arv_note) %>%
  mutate(
    arv_matches_date = !ID %in% c(10742, 10204, 40700),
    arv_date = case_when(
      arv_matches_date ~ draw_date,
      ID == 10742 ~ "12-Mar-14",
      ID == 10204 ~ "03-Aug-12",
      ID == 40700 ~ "29-Apr-14",
      TRUE ~ NA_character_
      ),
    adj_time = as.numeric(difftime(as.Date(draw_date, format = "%d-%b-%y"),
                                 as.Date(arv_date, format = "%d-%b-%y"), 
                                 units = "days")),
    arv_day = days - adj_time
    )

ARV_notes
```

```
## # A tibble: 25 x 8
## # Groups:   ID [25]
##       ID draw_date  days arv_note arv_matches_date arv_date adj_time
##    <int> <chr>     <dbl> <chr>    <lgl>            <chr>       <dbl>
##  1 10066 10-Mar-14   906 ARV?     TRUE             10-Mar-…        0
##  2 10203 27-Mar-13   530 ARV sta… TRUE             27-Mar-…        0
##  3 10428 15-May-13   846 ARV sta… TRUE             15-May-…        0
##  4 10435 15-May-12   505 ARV sta… TRUE             15-May-…        0
##  5 10742 25-Feb-14   251 "ARV st… FALSE            12-Mar-…      -15
##  6 40061 21-May-10   188 ARV sta… TRUE             21-May-…        0
##  7 40094 19-Mar-13  1019 ARV sta… TRUE             19-Mar-…        0
##  8 40100 3-Nov-10    301 ARV STA… TRUE             3-Nov-10        0
##  9 40168 5-Sep-12    754 ARV sta… TRUE             5-Sep-12        0
## 10 40250 13-Jul-12   709 ARV sta… TRUE             13-Jul-…        0
## # … with 15 more rows, and 1 more variable: arv_day <dbl>
```

```r
write_csv(ARV_notes, "ARV_notes.csv")
```

## Primary kinetics


```r
rv217_raw = ARV_notes %>%
  select(ID, arv_day) %>%
  right_join(rv217_raw, by = "ID") %>%
  mutate(
    arv_note = str_detect(Notes, "ARV"),
    any_arv = any(is.na(arv_day)),
    primary_kinetics = days >= 0 & days <= 365 & (is.na(arv_day) | days < arv_day)
  )
```

# Finalize

Order and save


```r
names(rv217_raw) = str_remove_all(names(rv217_raw), "\\.")
glimpse(rv217_raw)
```

```
## Observations: 1,486
## Variables: 30
## Groups: ID [55]
## $ ID               <int> 10066, 10066, 10066, 10066, 10066, 10066, 10066…
## $ arv_day          <dbl> 906, 906, 906, 906, 906, 906, 906, 906, 906, 90…
## $ draw_date        <chr> "4-Feb-10", "18-Feb-10", "22-Jul-10", "6-Jan-11…
## $ days_dx          <dbl> -588, -574, -420, -252, -84, -21, 0, 1, 4, 8, 1…
## $ APTIMA           <chr> NA, NA, NA, NA, NA, "NR", "12.38", "15.12", "15…
## $ log10VL          <dbl> NA, NA, NA, NA, NA, NA, NA, 4.76, 5.83, 6.26, 7…
## $ VL               <int> NA, NA, NA, NA, NA, NA, NA, 57544, 676083, 1819…
## $ VL_site          <chr> NA, NA, NA, NA, NA, NA, NA, "M", "M", "M", "M",…
## $ FiebigStage      <chr> NA, NA, NA, NA, NA, NA, NA, "I-II", "I-II", "I-…
## $ CD4pct           <dbl> NA, NA, NA, NA, NA, NA, NA, 47, NA, NA, NA, 25,…
## $ CD8pct           <dbl> NA, NA, NA, NA, NA, NA, NA, 25, NA, NA, NA, 59,…
## $ NKpct            <dbl> NA, NA, NA, NA, NA, NA, NA, 13, NA, NA, NA, 7, …
## $ Bpct             <dbl> NA, NA, NA, NA, NA, NA, NA, 13, NA, NA, NA, 8, …
## $ CD4              <dbl> NA, NA, NA, NA, NA, NA, NA, 824, NA, NA, NA, 39…
## $ CD8              <dbl> NA, NA, NA, NA, NA, NA, NA, 442, NA, NA, NA, 92…
## $ NK               <dbl> NA, NA, NA, NA, NA, NA, NA, 212, NA, NA, NA, 12…
## $ B                <dbl> NA, NA, NA, NA, NA, NA, NA, 212, NA, NA, NA, 12…
## $ VLunit           <chr> "copies/mL", "copies/mL", "copies/mL", "copies/…
## $ priority         <chr> "Priority 1", "Priority 1", "Priority 1", "Prio…
## $ site             <chr> "Uganda", "Uganda", "Uganda", "Uganda", "Uganda…
## $ Notes            <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
## $ visit_code       <chr> "A", "B", "D", "F", "H", "SBV", "SBV", "1", "2"…
## $ EIA              <chr> "NR", "NR", "NR", "NR", "NR", NA, NA, "NR", NA,…
## $ WB               <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "IN…
## $ APTIMA_num       <dbl> NA, NA, NA, NA, NA, NA, 12.38, 15.12, 15.54, 14…
## $ posVL            <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE…
## $ days             <dbl> -589, -575, -421, -253, -85, -22, -1, 0, 3, 7, …
## $ arv_note         <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
## $ any_arv          <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE…
## $ primary_kinetics <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE…
```

```r
rv217_raw %>%
  select(ID, draw_date, days, VL, log10VL, VLunit, posVL, primary_kinetics, CD4, CD8, NK, B,
         APTIMA, APTIMA_num, days_dx, contains("pct"), everything()) %>%
  write_csv("../../data/RV217Clean.csv")
```
