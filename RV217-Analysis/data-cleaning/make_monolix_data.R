library(tidyverse)

exclude_ids = c("10502", "40640", "41002")

raw_dat = read_csv("../../data/RV217Clean.csv") 


# ----- Censoring ------------

tmp_cens = raw_dat %>%
  subset(primary_kinetics | days < 0) %>%
  group_by(ID) %>%
  mutate(total = n()) %>%
  subset(total > 2)

tmp_cens %>%
  mutate(cens = days < 0) %>%
  select(ID, days, log10VL, cens) %>%
  subset(!ID %in% exclude_ids) %>%
  write_csv("../../data/RV217MonoCens.csv")


# ----- No Censoring ------------

tmp = raw_dat %>%
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

# ----- APTIMA prediction + censoring ------------
# model comes from the APTIMA-VL analysis

mod_data = raw_dat %>%
  subset(ID != 40320) %>% # poor data
  subset(days == 0) %>%
  select(ID, VL, log10VL, APTIMA_num) %>%
  na.omit() %>% # complete case subsetting but not actually necessary
  mutate(
    ID = factor(ID), 
    log10_APTIMA = log10(APTIMA_num)
  )

pred_mod = lm(log10VL ~ APTIMA_num, 
              data = subset(mod_data, APTIMA_num >= 9 & APTIMA_num <= 34))

raw_dat %>%
  subset(primary_kinetics | days < 0) %>%
  group_by(ID) %>%
  mutate(total = n()) %>%
  subset(total > 2) %>%
  mutate(
    predicted_VL = !is.na(APTIMA_num) & (days < 0 & (APTIMA_num >= 9 & APTIMA_num <= 34)),
    log10VL_pred = if_else(predicted_VL, 
                           predict(pred_mod, newdata = data.frame(APTIMA_num = APTIMA_num)),
                           log10VL),
    cens = days < 0,
    cens_pred = is.na(log10VL_pred)
    ) %>%
  group_by(ID) %>%
  mutate(
    first_pred_day = if(any(predicted_VL)) -(min(-days[days < 0])) else 0,
    cens_firstpred = days < first_pred_day
  ) %>%
  select(ID, days, log10VL, cens, APTIMA_num, predicted_VL, log10VL_pred,
         cens_pred, first_pred_day, cens_firstpred) %>%
  write_csv("../../data/RV217MonoVLPred.csv")





# ----- Compare to Dan's data ------------

dan_data = read_csv("../../Monolix/Data/DBRout_RV217.csv")
test = tmp %>%
  subset(days_dx <= 100) %>%
  select(ID, days_dx, log10VL) %>%
  left_join(dan_data, by = c("ID", "days_dx" = "days"))

#all(test$log10VL.x == test$log10VL.y, na.rm = T)
