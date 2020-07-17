rm(list=ls())
libs <- c("rstan", "readr")
invisible(lapply(libs, library, character.only = TRUE))

dir <- file.path("~/bayes_cpm_paper")

# read in and compile model
ord_mod_file1 <- read_file(file.path(dir,"ordinal_model_1.stan"))
ord_mod1 <- stan_model(model_code = ord_mod_file1)
saveRDS(ord_mod1, file = file.path(dir,"ordinal_model_1.rds"))

# make simarray
nsamps <- c(25,50,100,200,400)

simarray0 <- expand.grid(nsamps,1:5)
simarray <- cbind(simarray0,round(runif(25,0,1e5)))
names(simarray) <- c("nsamps","rep","seed")
simarray <- simarray[order(simarray$seed),]
saveRDS(simarray, file = file.path(dir,"bayes_cpm_simarray.rds"))
