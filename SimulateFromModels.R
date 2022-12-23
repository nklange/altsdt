# General information -----------------------

# Simulation of data under UVGaussian, Gumbel, ExGaussian model
# LargeN - for mimicry and correlation analysis
# Realistic item types, mimicking MWW2007 set up for confidence ratings

  library("cowplot")
  library("tidyverse")


  source("FitSimple.R")


  # Simulate data for Mimicry and correlations analysis ---------------------

# check out a reasonable range of parameter estimates

fits <- readRDS("EmpiricalFits/SB2021_e12_fits.rds")


genmodel <- "ExGaussNormEVSDT"
bestfit <- fits %>%
 # mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id,condition) %>%
  arrange(objective) %>%
  dplyr::slice(1)


datas <- bestfit %>%
  filter(model %in% genmodel) %>% group_by(condition) %>%
  mutate(strength = ifelse((condition == "A" | condition == "B"), "high","low"),
         variability = ifelse((condition == "A" | condition == "C"), "high","low")) %>%
  group_by(exp,condition,model) %>%
  select(c(model,id,condition,muo,betao,c1,dc1,dc2,dc3,dc4)) %>% ungroup() %>%
  pivot_longer(!c(exp,model,id,condition),names_to = "parameter",values_to="value") %>%
  group_by(parameter) %>%
  filter(value < quantile(value,.975) & value > quantile(value,.025)) %>%
  spread(parameter,value) %>%
  drop_na() %>%
  select(-c(exp,model,id,condition)) %>%
  select(names(get_start_par(genmodel,6)))

# for each model: identify mean/sd of paramaters and correlations between them

correlations <- cor(datas)
means <- as.numeric(t(datas %>% mutate_all(mean) %>% .[1,])[,1])
sds <- as.numeric(t(datas %>% mutate_all(sd) %>% .[1,])[,1])
sigmas <- sds %*% t(sds) * correlations

# Use multivariate to generate parameters on basis of multivariate distribution
# Sample 1000 possible sets of generating parameters per model

# Genpar: mvn from par estimates

parameters <- tmvtnorm::rtmvnorm(n=10000, mean = means, sigma=sigmas,
                                 lower=c(-Inf,0,-Inf,0,0,0,0),
                                 upper=rep(Inf,length(means)),
                                 algorithm="gibbs",
                                 burn.in.samples=1000,
                                 thinning = 100)
saveRDS(parameters,file="LargeNSimulation/simulate_genExGaussNormEVSDT_parametervalues.rds")


genmodel <- "GaussianUVSDT"
simulation <- NULL
parametertibble <-NULL
parametersext <- readRDS(paste0("LargeNSimulation/simulate_gen",genmodel,"_parametervalues.rds"))


for(i in c(1:dim(parametersext)[[1]])){

# for(j in c(2:51)){ # for each set of values, make 50 small-N datasets

id <- paste0("gen",genmodel,"_par",i)
#id <- paste0(parametersext[i,]$parid,"_set1")
parid <-paste0("gen",genmodel,"_par",i)# parametersext[i,]$parid#
pars <- as.numeric(parametersext[i,])#parameters[i,]


sim <- PredictSDT(data = NULL, model = genmodel,par = pars,
           itemspertype = c(50000,50000))

# itemspertype choice as a LargeN (50000) and standard-experiment (100)

simtibble <- sim %>% mutate(id = id,
                            parid = parid,
               genmodel = genmodel)

partibble <- tibble(value = pars,
                    parameter = names(get_start_par(genmodel))) %>%
  mutate(parid = parid,
         genmodel = genmodel)

  # partibble <- pars %>%
  #   set_colnames(names(get_start_par(genmodel))) %>%
  # mutate(parid = parid,
  #        genmodel = genmodel)


simulation <- simulation %>% bind_rows(simtibble)
parametertibble <- parametertibble %>% bind_rows(partibble)



}

#saveRDS(simulation,file=paste0("Simulations/LargeNSimulation/simulate_gen",genmodel,"_data_LN.rds"))
#saveRDS(parametertibble,file=paste0("Simulations/LargeNSimulation/simulate_gen",genmodel,"_parametervalues_LN.rds"))

# Simulated data for confidence ratings analysis-----------------------------


fits <- readRDS("EmpiricalFits/SB2021_e12_fits.rds")


bestfit <- fits %>%
  # mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id,condition) %>%
  arrange(objective) %>%
  dplyr::slice(1)


 genmodel <- "ExGaussNormEVSDT"

datas <- bestfit %>%
  filter(model %in% genmodel) %>% group_by(condition) %>%
  mutate(strength = ifelse((condition == "A" | condition == "B"), "high","low"),
         variability = ifelse((condition == "A" | condition == "C"), "high","low")) %>%
  group_by(exp,condition,model) %>%
  select(c(model,id,condition,names(get_start_par(genmodel,numcategories=6)))) %>%
  group_by(exp,model,id,condition) %>%
 # mutate(c5 = sum(c(c1,dc1,dc2,dc3,dc4))) %>% ungroup() %>%
  #select(-dc1,-dc2,-dc3,-dc4) %>%
  pivot_longer(!c(exp,model,id,condition),names_to = "parameter",values_to="value") %>%
  group_by(parameter) %>%
  filter(value < quantile(value,.975) & value > quantile(value,.025)) %>%
  spread(parameter,value) %>%
  drop_na() %>%
  select(-c(exp,model,id,condition)) %>%
  select(names(get_start_par(genmodel,6)))

#
# # for each model: identify mean/sd of paramaters and correlations between them
#
correlations <- cor(datas)
means <- as.numeric(t(datas %>% mutate_all(mean) %>% .[1,])[,1])
sds <- as.numeric(t(datas %>% mutate_all(sd) %>% .[1,])[,1])
sigmas <- sds %*% t(sds) * correlations
#
# # Use multivariate to generate parameters on basis of multivariate distribution
# # Sample 1000 possible sets of generating parameters per model
#
## Genpar: mvn from par estimates

parameters <- tmvtnorm::rtmvnorm(n=10000, mean = means, sigma=sigmas,
                                 lower=c(-Inf,0,-Inf,rep(0,4)),
                                 upper=rep(Inf,7),
                                 algorithm="gibbs",
                                 burn.in.samples=1000,
                                 thinning = 100)
# saveRDS(parameters,file="ConfidenceSimulation/simulate_genExGaussNormEVSDT_parametervalues4_unif.rds")
#
#
#
genmodel <- "ExGaussNormEVSDT"
simulation <- NULL
parametertibble <-NULL
parametersext <- readRDS(paste0("Simulations/ConfidenceSimulation/simulate_gen",genmodel,"_parametervalues4_unif.rds"))


for(i in c(1:dim(parametersext)[[1]])){


# Extend 5 thresholds to 23 thresholds ...

  # Create increments that are based on distances between points on
  # original 5 thresholds

      incr1 <- runif(4)
      incr1 <- cumsum(incr1/sum(incr1) * parametersext[i,4])
      incr2 <- runif(4)
      incr2 <- cumsum(incr2/sum(incr2) * parametersext[i,5])

      incr3 <- runif(4)
      incr3 <- cumsum(incr3/sum(incr3) * parametersext[i,6])
      incr4 <- runif(4)
      incr4 <- cumsum(incr4/sum(incr4) * parametersext[i,7])

      # add intermediate ratings with lowest and highest are extended beyond
      # the extremes of the original 5 thresholds
    prec1 <- c(parametersext[i,3] - incr1[[3]],
               parametersext[i,3] - incr1[[2]],
               parametersext[i,3] - incr1[[1]])
    c1 <- c(parametersext[i,3])
    prec2 <- c(parametersext[i,3] + incr1[[1]],
               parametersext[i,3] + incr1[[2]],
               parametersext[i,3] + incr1[[3]])
    c2 <- c(parametersext[i,3] + incr1[[4]])
    prec3 <- c(c2 + incr2[[1]],
               c2 + incr2[[2]],
               c2 + incr2[[3]])
    c3 <- c(c2 +incr2[[4]])
    prec4 <- c(c3 + incr3[[1]],
               c3 + incr3[[2]],
               c3 + incr3[[3]])
    c4 <- c(c3 + incr3[[4]])
    prec5 <- c(c4 + incr4[[1]],
               c4 + incr4[[2]],
               c4 + incr4[[3]])
    c5 <- c(c4 + incr4[[4]])
    postc5 <- c(c5 + incr1[[1]],
                c5 + incr1[[2]],
                c5 + incr1[[3]])

  critsort <- c(prec1,c1,prec2,c2,prec3,c3,prec4,c4,prec5,c5,postc5)


  numcategories <- 24
  crits <- cumsum(parametersext[i,c(3:length(parametersext[i,]))])

  # UNIFORM CRIT WITH SOME RESTRICTION ON EXTREMES!
  # critsort<- sort(runif(numcategories-1,
  #                       min=crits[[1]]-0.5,
  #                       max=crits[[5]]+0.5))
  # crit<-c(critsort[1],diff(critsort))

  # EQUIDISTANT CRIT WITH SOME RESTRICTION ON EXTREMES!
 # even<-abs(diff(c((crits[[1]]-0.5), (crits[[5]]+0.5))))/(numcategories-2)
 #
 # crit <- c(crits[[1]]-0.5,rep(even,22))
 #

  id <- paste0("gen",genmodel,"_par",i)
  #id <- paste0(parametersext[i,]$parid,"_set1")
  parid <-paste0("gen",genmodel,"_par",i)# parametersext[i,]$parid#

  if(genmodel=="GumbelEVSDT"){
    pars <- c(parametersext[i,c(1)],crit)
  } else {
    pars <- c(parametersext[i,c(1:2)],crit)
  }


  sim <- PredictSDT(data = NULL, model = genmodel,par = pars,
                    itemspertype = c(300,300),numcategories=numcategories)

  # itemspertype & numcategories based on Exp 1 in Mickes, Wais & Wixted

  simtibble <- sim %>% mutate(id = id,
                              parid = parid,
                              genmodel = genmodel)

  partibble <- tibble(value = pars,
                      parameter = names(get_start_par(genmodel,numcategories))) %>%
    mutate(parid = parid,
           genmodel = genmodel)

  # partibble <- pars %>%
  #   set_colnames(names(get_start_par(genmodel))) %>%
  # mutate(parid = parid,
  #        genmodel = genmodel)


  simulation <- simulation %>% bind_rows(simtibble)
  parametertibble <- parametertibble %>% bind_rows(partibble)



}

#saveRDS(simulation,file=paste0("Simulations/ConfidenceSimulation/simulate_gen",genmodel,"_data4_unif.rds"))
#saveRDS(parametertibble,file=paste0("Simulations/ConfidenceSimulation/simulate_gen",genmodel,"_parametervalues4_unif.rds"))


# recode 24point to 6 point scale and fit

genmodel <- "ExGaussNormEVSDT"
sim <- readRDS(paste0("Simulations/ConfidenceSimulation/simulate_gen",genmodel,"_data5_equi.rds"))

sim1 <- sim %>%
  mutate(response = case_when(response %in% c(1:4) ~ 1,
                              response %in% c(5:8) ~ 2,
                              response %in% c(9:12) ~ 3,
                              response %in% c(13:16) ~ 4,
                              response %in% c(17:20)~5,
                              response %in% c(21:24) ~6)) %>%
  mutate(response = factor(response)) %>%
  group_by(id,oldnew,response) %>%
  mutate(Freq = sum(Freq)) %>%
  distinct()

#saveRDS(sim1,file=paste0("ConfidenceSimulation/simulate_gen",genmodel,"_data5_equi_sixpoint.rds"))

