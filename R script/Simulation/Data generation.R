library(PCaASSim)

set.seed(2024)
seed <- sample(1:1e5, 200)
save(seed, file="Output/Simulation/Seeds/DGM_seed.RData")

for (i in 1:200) {
  sim <- mcicjmsim(n = 1000, seed = seed[i],
                   sens_fit = 0.75, sens_sim = 0.75, keep_complete = T)
  train.data <- sim$dat[sim$dat$CISNET_ID %in% 1:500,]
  train.biopsy <- sim$biopsy[sim$biopsy$CISNET_ID%in% 1:500,]
  save(train.data,
       file = paste0("Output/Simulation/Training sets/traindata_",
                     i,
                     ".RData"))
  save(train.biopsy,
       file = paste0("Output/Simulation/Training sets/trainbiopsy_",
                     i,
                     ".RData"))
}
