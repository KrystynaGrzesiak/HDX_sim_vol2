
library(dplyr)

all_params = readRDS("./all_params.RDS")


nrow(all_params)


all_params = all_params %>% 
  select(-protection_factor, -fraction) %>% 
  unique() 



all_params1 = all_params %>% 
  mutate(protection_factor = 10)

all_params2 = all_params %>% 
  mutate(protection_factor = 100)

all_params3 = all_params %>% 
  mutate(protection_factor = 100000)


all_params_new_pf = rbind(all_params1, all_params2, all_params3)


nrow(all_params_new_pf)


saveRDS(all_params_new_pf, file = "all_params_new_pf.RDS")
