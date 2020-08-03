rm(list=ls())

dir<-file.path("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper/sims")

set.seed(77384)
# make simarray
nsamps <- c(25,50,100,200,400)

simarray0 <- expand.grid(nsamps,1:10)
simarray <- cbind(simarray0,round(runif(50,0,1e5)))
names(simarray) <- c("nsamps","rep","seed")
simarray <- simarray[order(simarray$seed),]
saveRDS(simarray, file = file.path(dir,"bayes_cpm_simarray2.rds"))

# in accre slurm scripts make array over 1-50, time=16:00:00
# in sim files use bayes_cpm_simarray2.rds and sim=100 for sim function
