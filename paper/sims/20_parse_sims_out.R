rm(list=ls())
libs <- c("dplyr", "stringr", "readr", "tidyr", "purrr", "ggplot2")
invisible(lapply(libs, library, character.only = TRUE))

#seminar directory (sims a for conc=1/J & probit link)
sdir <- file.path("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/biostat_sem")
ssdir <- file.path(sdir,"sims","out")

#paper directory (sims for other functions)
pdir <- file.path("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper")
psdir <- file.path(pdir,"sims","out")

figdir <- file.path(pdir,"fig")

#source(file.path(dir,"cpm_functions.r"))

# load sim arrays
simarray1 <- readRDS(file.path(pdir,"sims","bayes_cpm_simarray.rds")) # for a,b,c
simarray2 <- readRDS(file.path(pdir,"sims","bayes_cpm_simarray2.rds")) # for d

# function to parse simulation .out files to get processor info
parse_lines_pr <- function(lines_in){
  mn<-grep("model.name", lines_in)
  lisub <- lines_in[mn] #subset of lines_in

return(lisub)
}


# function to parse simulation .out files to get sampling times
parse_lines_tm <- function(lines_in, simarrnum, simarr=simarray1, nsims=200){
  iternum0 <- grep("\\[1\\] [1-9]", lines_in) # beginning of iter
  ctot0 <- grep("Chain [1-2]:  Elapsed Time:", lines_in) + 2 # total time
  lisub <- lines_in[sort(c(iternum0,ctot0))] #subset of lines_in
  
  iternum <- grep("\\[1\\] [1-9]", lisub) 
  
  listc <- listb <- lista <- vector("list", nsims)
  
  for (i in 1:length(lisub)){
    if(i %in% iternum){
      iter <- as.numeric(substr(lisub[i],5,10))
    } else {
      lista[[iter]] <- rbind(lista[[iter]], lisub[i])
    }
  }
  
  listb <- lapply(lista, data.frame)
  
  for(n in 1:nsims){
    df <- listb[[n]]
    if (nrow(df)) {
      listc[[n]] <- separate(df, col=X..i.., into=c("chn","secs","type"), sep=":|seconds") %>% 
        mutate(iter=n)
    }
  }
  
  dd <- bind_rows(listc) %>% 
    group_by(iter) %>% 
    mutate(pos=1:n(), nc=n(), 
           secs=as.numeric(secs),
           mod=if_else(pos<3,"full","cens"),
           mod=if_else(nc==4,mod,"unk"),
           nsamps=simarr[simarrnum,"nsamps"],
           rep=simarr[simarrnum,"rep"]) %>% 
    select(-c(pos,nc)) %>% ungroup()
  
  return(dd)
}

#old parse funs and test/comparison
if(0){
parse_lines_tm00 <- function(lines_in, simarrnum, simarr=simarray1, nsims=200){
  iternum <- grep("\\[1\\] [1-9]", lines_in) 
  c1tot <- grep("Chain 1:  Elapsed Time:", lines_in)+2 # total time
  c2tot <- grep("Chain 2:  Elapsed Time:", lines_in)+2 # total time
  
  listc <- lista <- vector("list", nsims)
  
  for (i in 1:length(lines_in)){
    if(i %in% iternum){
      iter <- as.numeric(substr(lines_in[i],5,10))
    }
    if(i %in% c(c1tot,c2tot)){
      lista[[iter]]<-rbind(lista[[iter]], lines_in[i])
    }
  }
  
  blist <- lapply(alist,data.frame)
  
  for(n in 1:nsims){
    df <- blist[[n]]
    if (nrow(df)) {
      listc[[n]] <- separate(df, col=X..i.. , into=c("chn","time"), sep=":") %>% 
        separate(col=time, into=c("secs","type"), sep="seconds") %>% 
        mutate(iter=n, secs=as.numeric(secs))
    }
  }
  
  dd <- bind_rows(listc) %>% group_by(iter) %>% 
    mutate(pos=1:n(), nc=n(), 
           mod=if_else(pos<3,'full','cens'),
           mod=if_else(nc==4,mod,'unk')) %>% 
    select(-c(pos,nc)) %>% ungroup() %>% 
    mutate(nsamps=simarr[simarrnum,"nsamps"], rep=simarr[simarrnum,"rep"])
  
  return(dd)
}

parse_lines_tm0 <- function(lines_in, simarrnum, simarr=simarray1, nsims=200){
  iternum0 <- grep("\\[1\\] [1-9]", lines_in) 
  c1tot0 <- grep("Chain 1:  Elapsed Time:", lines_in)+2 # total time
  c2tot0 <- grep("Chain 2:  Elapsed Time:", lines_in)+2 # total time
  
  lisub <- lines_in[sort(c(iternum0,c1tot0,c2tot0))] #subset of lines_in
  
  iternum <- grep("\\[1\\] [1-9]", lisub) 
  c1tot <- grep("Chain 1:", lisub) # total time
  c2tot <- grep("Chain 2:", lisub) # total time
  
  listc <- blist <- lista <- vector("list", nsims)
  
  for (i in 1:length(lisub)){
    if(i %in% iternum){
      iter <- as.numeric(substr(lisub[i],5,10))
    }
    if(i %in% c(c1tot,c2tot)){
      lista[[iter]]<-rbind(lista[[iter]], lisub[i])
    }
  }
  
  blist <- lapply(alist,data.frame)
  
  for(n in 1:nsims){
    df <- blist[[n]]
    if (nrow(df)) {
      listc[[n]] <- separate(df, col=X..i.. , into=c("chn","time"), sep=":") %>% 
        separate(col=time, into=c("secs","type"), sep="seconds") %>% 
        mutate(iter=n, secs=as.numeric(secs))
    }
  }
  
  dd <- bind_rows(listc) %>% group_by(iter) %>% 
    mutate(pos=1:n(), nc=n(), 
           mod=if_else(pos<3,'full','cens'),
           mod=if_else(nc==4,mod,'unk')) %>% 
    select(-c(pos,nc)) %>% ungroup() %>% 
    mutate(nsamps=simarr[simarrnum,"nsamps"], rep=simarr[simarrnum,"rep"])
  
  return(dd)
}

li1 <- read_lines("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper/sims/out/sim_c3_24.out")
li2 <- read_lines("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper/sims/out/sim_a2_12.out")

mb1 <- microbenchmark(parse_lines_tm00(li1,simarrnum=24),
                      parse_lines_tm0(li1,simarrnum=24),
                      parse_lines_tm(li1,simarrnum=24),
                      times=25,
                      check='equal')

ggplot2::autoplot(mb1)


mb2 <- microbenchmark(parse_lines_tm00(li2,simarrnum=12),
                      parse_lines_tm0(li2,simarrnum=12),
                      parse_lines_tm(li2,simarrnum=12),
                      times=25,
                      check='equal')

ggplot2::autoplot(mb2)

}

#!paste0("sim_",i,j,"_",k)
#!foo<-read_lines(file.path(ssdir,"sim_1.out"))

# load sim data
# sims a, b
for (i in letters[1:2]){
  for (j in 0:2){
    for (k in 1:25){
        nm <- paste0("sim_",i,j,"_",k)
        nm_prc <- paste0(nm,"_proc")
        if (i=="a" & j=="0"){
          dir <- ssdir
          fn <- paste0("sim_",k,".out")
        } else {
          dir <- psdir
          fn <- paste0(nm,".out")
        }
        #dir <- ifelse(i=="a" & j=="0",ssdir,psdir)
        fp<-file.path(dir,fn)
        print(fp); print(nm)
        try(assign(nm, parse_lines_tm(read_lines(fp),k)))
        try(assign(nm_prc, parse_lines_pr(read_lines(fp))))
    }
  }
}

kk<-25
for (i in letters[1:2]){
  for (j in 0:2){
     # paste0("sim_",i,j,"_",k) 
     sim_setting <- paste0("sim_",i,j)
    #! print(sim_setting)
     sim_nm <- paste0(sim_setting,"_",1:kk)
    #! print(sim_nm)
    simnmlist <- map(sim_nm,get) 
    try(assign(sim_setting, bind_rows(simnmlist) %>% mutate(scenario=sim_setting) ))
    rm(list=sim_nm)
  }
}


#ls()[grep("_proc",ls())]


pltw<-10; plth<-5

# summarize sim scenario a 
sim_a<-bind_rows(sim_a0,sim_a1,sim_a2)

sim_a %>% group_by(scenario,nsamps,mod) %>% 
  summarize(mn_time=mean(secs,na.rm=TRUE),md_time=median(secs,na.rm=TRUE))

sim_a %>% mutate(nsamps=factor(nsamps,levels=c(25,50,100,200,400)),
                 scenario=factor(scenario,labels=c("alpha=1/J",
                                                   "alpha=1/(0.8+0.35*J)",
                                                   "alpha=1/(2+(J/3))"))) %>% 
  filter(mod!='unk') %>% 
  ggplot(aes(x=nsamps,y=secs)) +
  geom_boxplot(aes(color=mod)) +
  scale_y_log10() +
  facet_wrap(~scenario)

ggsave(file.path(figdir,"sim_a_MCMC_sample_times.png"),width=pltw,height=plth)

# summarize sim scenario b - model is same as sim scenario a
sim_b<-bind_rows(sim_b0,sim_b1,sim_b2)

sim_b %>% group_by(scenario,nsamps,mod) %>% 
  summarize(mn_time=mean(secs,na.rm=TRUE),md_time=median(secs,na.rm=TRUE))

sim_b %>% mutate(nsamps=factor(nsamps,levels=c(25,50,100,200,400)),
                 scenario=factor(scenario,labels=c("alpha=1/J",
                                                   "alpha=1/(0.8+0.35*J)",
                                                   "alpha=1/(2+(J/3))"))) %>% 
  filter(mod!='unk') %>% 
  ggplot(aes(x=nsamps,y=secs)) +
    geom_boxplot(aes(color=mod)) +
    scale_y_log10() +
    facet_wrap(~scenario)

ggsave(file.path(figdir,"sim_b_MCMC_sample_times.png"),width=pltw,height=plth)

# sims c
# sim_c0_n - logistic link, log transform, conc=1/ncats
# sim_c1_n - logistic link, log transform, conc=1/(0.8 + 0.35*max(ncats, 3))
# sim_c2_n - logistic link, log transform, conc=1/(2+(ncats/3))
# sim_c3_n - logistic link, log transform, conc=1/2

for (i in letters[3]){
  for (j in 0:3){
    for (k in 1:25){
      nm <- paste0("sim_",i,j,"_",k)
      nm_prc <- paste0(nm,"_proc")
      if (i=="a" & j=="0"){
        dir <- ssdir
        fn <- paste0("sim_",k,".out")
      } else {
        dir <- psdir
        fn <- paste0(nm,".out")
      }
      #dir <- ifelse(i=="a" & j=="0",ssdir,psdir)
      fp<-file.path(dir,fn)
      print(fp); print(nm)
      try(assign(nm, parse_lines_tm(read_lines(fp),k)))
      try(assign(nm_prc, parse_lines_pr(read_lines(fp))))
    }
  }
}

kk<-25
for (i in letters[3]){
  for (j in 0:3){
    # paste0("sim_",i,j,"_",k) 
    sim_setting <- paste0("sim_",i,j) 
    #! print(sim_setting)
    sim_nm <- paste0(sim_setting,"_",1:kk)
    #! print(sim_nm)
    simnmlist <- map(sim_nm,get) 
    try(assign(sim_setting, bind_rows(simnmlist) %>% mutate(scenario=sim_setting) ))
    rm(list=sim_nm)
  }
}

# summarize sim scenario c
sim_c<-bind_rows(sim_c0,sim_c1,sim_c2,sim_c3)

sim_c %>% group_by(scenario,nsamps,mod) %>% 
  summarize(mn_time=mean(secs,na.rm=TRUE),md_time=median(secs,na.rm=TRUE))

sim_c %>% filter(scenario!='sim_c3') %>% 
  mutate(nsamps=factor(nsamps,levels=c(25,50,100,200,400)),
         scenario=factor(scenario,labels=c("alpha=1/J",
                                           "alpha=1/(0.8+0.35*J)",
                                           "alpha=1/(2+(J/3))"))) %>% 
  filter(mod!='unk') %>% 
  ggplot(aes(x=nsamps,y=secs)) +
  geom_boxplot(aes(color=mod)) +
  scale_y_log10() +
  facet_grid(~scenario)

ggsave(file.path(figdir,"sim_c_MCMC_sample_times.png"),width=pltw,height=plth)

# sims d
# sim_d0_n - loglog link, no transform, conc=1/ncats
# sim_d1_n - loglog link, no transform, conc=1/(0.8 + 0.35*max(ncats, 3))
# sim_d2_n - loglog link, no transform, conc=1/(2+(ncats/3))
# sim_d3_n - loglog link, no transform, conc=1/2

# check
#li3 <- read_lines("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper/sims/out/sim_d2_12.out")
#parse_lines_tm(li3,12,simarr=simarray2,nsims=100)

kk<-50
for (i in letters[4]){
  for (j in 0:3){
    for (k in 1:kk){
      nm <- paste0("sim_",i,j,"_",k)
      nm_prc <- paste0(nm,"_proc")
      if (i=="a" & j=="0"){
        dir <- ssdir
        fn <- paste0("sim_",k,".out")
      } else {
        dir <- psdir
        fn <- paste0(nm,".out")
      }
      #dir <- ifelse(i=="a" & j=="0",ssdir,psdir)
      fp<-file.path(dir,fn)
      print(fp); print(nm)
      try(assign(nm, parse_lines_tm(read_lines(fp),k,simarr=simarray2,nsims=100)))
      try(assign(nm_prc, parse_lines_pr(read_lines(fp))))
    }
  }
}

for (i in letters[4]){
  for (j in 0:3){
    # paste0("sim_",i,j,"_",k) 
    sim_setting <- paste0("sim_",i,j) 
    #! print(sim_setting)
    sim_nm <- paste0(sim_setting,"_",1:kk)
    #! print(sim_nm)
    simnmlist <- map(sim_nm,get) 
    try(assign(sim_setting, bind_rows(simnmlist) %>% mutate(scenario=sim_setting) ))
    rm(list=sim_nm)
  }
}

# summarize sim scenario d
sim_d<-bind_rows(sim_d0,sim_d1,sim_d2,sim_d3)

sim_d %>% group_by(scenario,nsamps,mod) %>% 
  summarize(mn_time=mean(secs,na.rm=TRUE),md_time=median(secs,na.rm=TRUE))

sim_d %>% filter(scenario!='sim_d3') %>%
 mutate(nsamps=factor(nsamps,levels=c(25,50,100,200,400)),
        scenario=factor(scenario,labels=c("alpha=1/J",
                                          "alpha=1/(0.8+0.35*J)",
                                          "alpha=1/(2+(J/3))"))) %>% 
  filter(mod!='unk') %>% 
  ggplot(aes(x=nsamps,y=secs)) +
  geom_boxplot(aes(color=mod)) +
  scale_y_log10() +
  facet_grid(~scenario)

ggsave(file.path(figdir,"sim_d_MCMC_sample_times.png"),width=pltw,height=plth)



## combined plot
pltw<-10; plth<-8

sim_a_mod <- sim_a %>% mutate(scn='scenario 1', prior=substr(scenario,6,7))
# drop alpha=1/2 priors for logit and loglog
sim_c_mod <- sim_c %>% filter(scenario!='sim_c3') %>% mutate(scn='scenario 2', prior=substr(scenario,6,7))
sim_d_mod <- sim_d %>% filter(scenario!='sim_d3') %>% mutate(scn='scenario 3', prior=substr(scenario,6,7))

rbind(sim_a_mod,sim_c_mod,sim_d_mod) %>%
  mutate(nsamps=factor(nsamps,levels=c(25,50,100,200,400)),
         prior=factor(prior,labels=c("alpha==1/J",
                                           "alpha==1/(0.8+0.35*J)",
                                           "alpha==1/(2+(J/3))"))) %>% 
  filter(mod!='unk') %>% 
  mutate(mod=factor(mod,labels=c("censored Y","uncensored Y"))) %>% 
  ggplot(aes(x=nsamps,y=secs)) +
  geom_boxplot(aes(color=mod)) +
  scale_y_log10() + xlab("sample size") + ylab("seconds") +
  scale_color_discrete(name="outcome") +
  facet_grid(scn~prior, labeller = labeller(prior=label_parsed), switch="y")


ggsave(file.path(figdir,"sim_MCMC_sample_times.png"),width=pltw,height=plth)
ggsave(file.path(figdir,"sim_MCMC_sample_times.tiff"),width=pltw,height=plth,dpi=600, compression = "lzw")


## summarize processors used 
processor_list<-map(ls()[grep("_proc",ls())],get)
processor_info<-do.call(rbind,processor_list)

table(processor_info)

### scratch 
if (0){
la<-read_lines("/home/nathan/Dropbox/njames/school/PhD/orise_ra/bayes_cpm/paper/sims/out/sim_c3_23.out")

a<-grep("R version", la)
b<-grep("Platform", la)
c<-grep("Running under:", la)
d<-grep("model.name", la)
e<-grep("\\[1\\] [1-9]", la)
f<-grep("Chain 1:  Elapsed Time:", la)+2 # total time
g<-grep("Chain 2:  Elapsed Time:", la)+2 # total time

la2<-la[sort(c(a,b,c,d,e,f,g))]

# la.small<-head(la2,16)
# data.frame(iter=la[e],chns=la[sort(c(f,g))])
# foo<-strsplit(la2,la[e])

alist<-vector("list", 200)
for (i in 1:length(la)){
  if(i %in% e){
    iter <- as.numeric(substr(la[i],5,10))
  }
  if(i %in% c(f,g)){
    alist[[iter]]<-rbind(alist[[iter]], la[i])
  }
}

blist<-lapply(1:200, function(x) data.frame(alist[[x]]))

# separate(blist[[1]], col=alist..x.. ,into=c("chn","time"),sep=":")
# separate(blist[[2]], col=alist..x.. ,into=c("chn","time"),sep=":")

clist<-vector("list", 200)
for(n in 1:200){
  df<-blist[[n]]
  if (nrow(df)) {
    df2 <- separate(df, col=alist..x.. , into=c("chn","time"),sep=":")
    clist[[n]] <- separate(df2, col=time, into=c("secs","type"),sep="seconds") %>%
      mutate(iter=n,secs=as.numeric(secs))
    }
}

dd<-bind_rows(clist) %>% group_by(iter) %>%
  mutate(pos=1:n(), nc=n(),
         mod=if_else(pos<3,'full','cens'),
         mod=if_else(nc==4,mod,'unk')) %>%
  select(-c(pos,nc)) %>% ungroup()

#group_by(mod) %>% summarize(mn_time=mean(secs),md_time=median(secs))

# doesn't work because of failed sims for some iters
#clist <- lapply(blist, function(x) separate(x,col="alist..x..", into=c("chn","time"),sep=":"))

}
