---
title: "RV217 - Data Cleaning"
author: "Bryan Mayer"
date: "7/10/2019"
output: 
  html_document:
    keep_md: true
    toc: true
---




The following script cleans and tidies the raw RV217 data. It also cleans and saves supplementary table 4 from the paper at the end. 

**Some notes:**

  - `days` was recalculated to be from the time of the first positive viral load. 
    - Original days variable is now called `days_dx` for days from diagnosis (APTIMA  detection).
  - Used Dan's variable name conventions so our datasets have similar headers
  - Removed any row with missing ID and combined missing day/date: there was never any output on those rows.
    - This drops a lot of the notes columns

**Some issues:**

- Not every participant has site information
- Not sure what many of the variables mean
- APTIMA needs an LLoQ? Currently: 
    - NR = 0
    - R or R/NR = min(APTIMA)
- Do we need to worry about ARV, some of the dates in the notes are not exact
- Did not clean the date variable (`draw_date`)

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
## ── Attaching packages ────────────────────────────────────────────────────────── tidyverse 1.2.1 ──
```

```
## ✔ ggplot2 3.1.1     ✔ purrr   0.3.2
## ✔ tibble  2.1.3     ✔ dplyr   0.8.1
## ✔ tidyr   0.8.3     ✔ stringr 1.4.0
## ✔ readr   1.3.1     ✔ forcats 0.4.0
```

```
## ── Conflicts ───────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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
rv217 = read.csv("../../data/RV217Master.csv", stringsAsFactors = F,
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

glimpse(rv217)
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

Missing PIN seems to come from populated notes columns
- no cases of viral load or cell data in them so removing


```r
missing_id = rv217 %>%
  subset(is.na(ID))

select(missing_id, contains("VL"), contains("CD"), contains("B.")) %>%
  na.omit %>%
  nrow()
```

```
## [1] 0
```

```r
rv217 = subset(rv217, !is.na(ID))
```

Checking types of our important endpoints (viral load and cell counts). Random character string in log10VL


```r
rv217 %>%
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
walk(rv217$log10VL, function(i){
  if(!is.na(i) & is.na(as.numeric(i))) print(i)
})
```

```
## Warning in .f(.x[[i]], ...): NAs introduced by coercion
```

```
## [1] "Aptima Reactive"
```

```r
rv217$log10VL[which(rv217$log10VL == "Aptima Reactive")] = NA
rv217$log10VL = as.numeric(rv217$log10VL)

rv217 %>%
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
priority_site_dat = rv217 %>% group_by(ID) %>%
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
rv217 = rv217 %>%
  select(-Priority, -Site) %>%
  left_join(select(priority_site_dat, ID, priority, site), by = "ID") %>%
  select(-visit_code, -Notes, -EIA, -WB, everything()) # moves misc. variables to the end
```

## APTIMA cleaning

Make a variable that is measurement only truncated at minimum (for R and R/NR), truncated at 0 for NR and left missing for weird inputs


```r
aptima_vals = map_chr(rv217$APTIMA, function(i){
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
rv217 = rv217 %>%
  mutate(
    APTIMA_num_temp = suppressWarnings(as.numeric(APTIMA)),
    APTIMA_num = case_when(
      APTIMA == "NR" ~ 0,
      APTIMA %in% c("R", "R/NR") ~ min(APTIMA_num_temp, na.rm = T),
      TRUE ~ APTIMA_num_temp
    )
  ) %>%
  select(-APTIMA_num_temp)
```


# days variable

## Days variable error checking

There are several cases with missing day variables. Many are due to lack of a date variable in combination with no viral load or cell counts, so we can remove those rows.


```r
missing_time_all = rv217 %>%
  group_by(ID) %>%
  mutate(
    any_missing_day = any(is.na(days.from.1st.pos))
  ) %>%
  subset(any_missing_day)

# missing times of + VL collection
vl_missing_time = missing_time_all %>%
  group_by(ID) %>%
  mutate(
   vl_missing_time = any(is.na(days.from.1st.pos) & !is.na(VL))
  ) %>%
  subset(vl_missing_time)

ids_missing_days = unique(vl_missing_time$ID)

#missing both (checking if we can just exclude these)
missing_time_vl = missing_time_all %>%
  subset(!ID %in% ids_missing_days) %>%
  select(ID, days.from.1st.pos, draw_date, contains("VL"), contains("CD"), contains("B."))
  
# one weird entry
missing_time_vl %>%
  subset(is.na(days.from.1st.pos) & !is.na(draw_date))
```

```
## # A tibble: 1 x 13
## # Groups:   ID [1]
##      ID days.from.1st.p… draw_date log10VL    VL VL_site VLunit CD4.pct
##   <int>            <int> <chr>       <dbl> <int> <chr>   <chr>    <dbl>
## 1 40491               NA NO SAMPL…      NA    NA <NA>    copie…      NA
## # … with 5 more variables: CD8.pct <dbl>, CD4. <dbl>, CD8. <dbl>,
## #   B.pct <dbl>, B. <dbl>
```


```r
clean_vl_missing_time = vl_missing_time %>%
  group_by(ID) %>%
  mutate(
    infection_time = na.omit(draw_date[APTIMA == "R"]),
    days.from.1st.pos = as.numeric(difftime(as.Date(draw_date, format = "%d-%b-%y"),
                                 as.Date(infection_time, format = "%d-%b-%y"), 
                                 units = "days"))
  ) %>%
  select(-infection_time, -any_missing_day, -vl_missing_time)
```

## Days variable creation


```r
rv217 = rv217 %>%
  subset(!is.na(days.from.1st.pos) & !is.na(draw_date)) %>%
  subset(!ID %in% ids_missing_days) %>% # these are processed separately
  bind_rows(clean_vl_missing_time) %>% # and added back here
  rename(days_dx = days.from.1st.pos) %>%
  group_by(ID) %>%
  mutate(
    posVL = !is.na(VL),
    days = days_dx - min(days_dx[posVL])
  ) 

any(is.na(rv217$days))
```

```
## [1] FALSE
```


# Finalize

Order and save


```r
names(rv217) = str_remove_all(names(rv217), "\\.")
glimpse(rv217)
```

```
## Observations: 1,464
## Variables: 26
## Groups: ID [55]
## $ ID          <int> 10066, 10066, 10066, 10066, 10066, 10066, 10066, 100…
## $ draw_date   <chr> "4-Feb-10", "18-Feb-10", "22-Jul-10", "6-Jan-11", "2…
## $ days_dx     <dbl> -588, -574, -420, -252, -84, -21, 0, 1, 4, 8, 11, 15…
## $ APTIMA      <chr> NA, NA, NA, NA, NA, "NR", "12.38", "15.12", "15.54",…
## $ log10VL     <dbl> NA, NA, NA, NA, NA, NA, NA, 4.76, 5.83, 6.26, 7.51, …
## $ VL          <int> NA, NA, NA, NA, NA, NA, NA, 57544, 676083, 1819701, …
## $ VL_site     <chr> NA, NA, NA, NA, NA, NA, NA, "M", "M", "M", "M", "M",…
## $ FiebigStage <chr> NA, NA, NA, NA, NA, NA, NA, "I-II", "I-II", "I-II", …
## $ CD4pct      <dbl> NA, NA, NA, NA, NA, NA, NA, 47, NA, NA, NA, 25, NA, …
## $ CD8pct      <dbl> NA, NA, NA, NA, NA, NA, NA, 25, NA, NA, NA, 59, NA, …
## $ NKpct       <dbl> NA, NA, NA, NA, NA, NA, NA, 13, NA, NA, NA, 7, NA, N…
## $ Bpct        <dbl> NA, NA, NA, NA, NA, NA, NA, 13, NA, NA, NA, 8, NA, N…
## $ CD4         <dbl> NA, NA, NA, NA, NA, NA, NA, 824, NA, NA, NA, 395, NA…
## $ CD8         <dbl> NA, NA, NA, NA, NA, NA, NA, 442, NA, NA, NA, 925, NA…
## $ NK          <dbl> NA, NA, NA, NA, NA, NA, NA, 212, NA, NA, NA, 122, NA…
## $ B           <dbl> NA, NA, NA, NA, NA, NA, NA, 212, NA, NA, NA, 127, NA…
## $ VLunit      <chr> "copies/mL", "copies/mL", "copies/mL", "copies/mL", …
## $ priority    <chr> "Priority 1", "Priority 1", "Priority 1", "Priority …
## $ site        <chr> "Uganda", "Uganda", "Uganda", "Uganda", "Uganda", "U…
## $ Notes       <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
## $ visit_code  <chr> "A", "B", "D", "F", "H", "SBV", "SBV", "1", "2", "3"…
## $ EIA         <chr> "NR", "NR", "NR", "NR", "NR", NA, NA, "NR", NA, "NR"…
## $ WB          <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "IND", "…
## $ APTIMA_num  <dbl> NA, NA, NA, NA, NA, 0.00, 12.38, 15.12, 15.54, 14.30…
## $ posVL       <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRU…
## $ days        <dbl> -589, -575, -421, -253, -85, -22, -1, 0, 3, 7, 10, 1…
```

```r
rv217 %>%
  select(ID, draw_date, days, VL, log10VL, VLunit, posVL, CD4, CD8, NK, B,
         APTIMA, APTIMA_num, days_dx, contains("pct"), everything()) %>%
  write_csv("../../data/RV217Clean.csv")
```
