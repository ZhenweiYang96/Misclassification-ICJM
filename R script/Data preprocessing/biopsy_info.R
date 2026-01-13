
# Packages ----------------------------------------------------------------
library(tidyverse)

# Load data ---------------------------------------------------------------
load("Data/Biopsy.RData")
load("Data/pass_id.RData")

ids <- biopsy$CISNET_ID[biopsy$CISNET_ID %in% pass.id$CISNET_ID]
biopsy <- biopsy %>% 
  filter(CISNET_ID %in% pass.id$CISNET_ID) %>% 
  mutate(time.cmp2 = rep(pass.id$time.cmp2, table(ids))) %>% 
  filter(TimeSince_Dx <= time.cmp2)

save(biopsy, file="Data/biopsy_sub.RData")

