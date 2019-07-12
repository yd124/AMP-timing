library(tidyverse)

read_csv("../data/RV217Clean.csv") %>%
  subset(!is.na(log10VL) & days >= 0) %>%
  group_by(ID) %>%
  mutate(total = n()) %>%
  subset(total > 2) %>%
  select(ID, days, log10VL) %>%
  write_csv("../data/RV217Mono.csv")
  