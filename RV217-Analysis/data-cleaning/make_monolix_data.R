library(tidyverse)

tmp = read_csv("../../data/RV217Clean.csv") %>%
  subset(!is.na(log10VL) & days >= 0) %>%
  group_by(ID) %>%
  mutate(total = n()) %>%
  subset(total > 2) 

tmp %>%
  select(ID, days, log10VL) %>%
  write_csv("../../data/RV217Mono.csv")
  
tmp %>%
  select(ID, days, log10VL, CD4, CD8) %>%
  gather(ytype_def, y, log10VL, CD4, CD8) %>%
  mutate(
    ytype = case_when(
      ytype_def == "log10VL" ~ 1,
      ytype_def == "CD4" ~ 2,
      ytype_def == "CD8" ~ 3
    )
  ) %>%
  write_csv("../../data/RV217MonoCells.csv")



