
library(dplyr)

spectra_by_ph_seq <- do.call(rbind, spectra_by_ph_seq)

spectra_by_ph_seq %>% 
  group_by(PH , Sequence, PF, Charge) %>% 
  summarise(is_like_ppvq = max(Exposure) < 43200) %>% 
  filter(is_like_ppvq == TRUE)

