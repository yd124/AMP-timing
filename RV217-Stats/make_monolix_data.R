library(tidyverse)

read_csv("../data/RV217Clean.csv") %>%
  subset(!is.na(log10VL) & days >= 0) %>%
  select(ID, days, log10VL) %>%
  write_csv("../data/RV217Mono.csv")
  