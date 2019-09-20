library(tidyverse)

exclude_ids = c("10502", "40640", "41002")

tmp = read_csv("../../data/RV217Clean.csv") %>%
  subset(!is.na(log10VL) & primary_kinetics) %>%
  group_by(ID) %>%
  mutate(total = n()) %>%
  subset(total > 2)

tmp %>%
  select(ID, days, log10VL) %>%
  subset(!ID %in% exclude_ids) %>%
  write_csv("../../data/RV217Mono.csv")

tmp %>%
  select(ID, days, log10VL, CD4, CD8) %>%
  subset(!ID %in% exclude_ids) %>%
  gather(ytype_def, y, log10VL, CD4, CD8) %>%
  mutate(
    ytype = case_when(
      ytype_def == "log10VL" ~ 1,
      ytype_def == "CD4" ~ 2,
      ytype_def == "CD8" ~ 3
    )
  ) %>%
  write_csv("../../data/RV217MonoCells.csv")


## compare to Dan's data
dan_data = read_csv("../../Monolix/Data/DBRout_RV217.csv")
test = tmp %>%
  subset(days_dx <= 100) %>%
  select(ID, days_dx, log10VL) %>%
  left_join(dan_data, by = c("ID", "days_dx" = "days"))

#all(test$log10VL.x == test$log10VL.y, na.rm = T)
