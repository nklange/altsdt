library(tidyverse)
library(waffle)
library(cowplot)

source("FitSimple.R")

models <- c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT","baselineUVSDT")
Gsquaredsingle <- function(modelname,dp,expdata,totitems){

  obs <- c(dp$datalist$New,
           dp$datalist$Old)
  exp <- expdata*totitems
  exp[exp==0]<-1/6

  f <- tibble( GSq = 2*sum(obs*log(obs/exp))) # Calculate G^2
  f$df     <- length(obs)-length(names(get_start_par(modelname,numcategories=6)))-1
  f$pGSq   <- 1-pchisq(f$GSq,f$df)       # p value of GSq

  return(f)
}
Gsquaredmult <- function(modelname,dp,expdata,totitems){

  obs <- c(dp$datalist$high_ev_1,
           dp$datalist$low_ev_1,
           dp$datalist$high_ev_0,
           dp$datalist$low_ev_0)
  exp <- expdata*totitems
  exp[exp==0]<-1/6

  f <- tibble( GSq = 2*sum(obs*log(obs/exp))) # Calculate G^2
  f$df     <- length(obs)-length(names(get_start_par(modelname,numcategories=6)))-1
  f$pGSq   <- 1-pchisq(f$GSq,f$df)       # p value of GSq

  return(f)
}



allfitsSB2020 <- readRDS("EmpiricalFits/SB2020_e123_fits.rds")
dataSB2020 <- bind_rows(readRDS("EmpiricalData/SB2020_data/SB2020_e1.rds") %>% mutate(item = as.character(item)),
                    readRDS("EmpiricalData/SB2020_data/SB2020_e2.rds") %>% mutate(item = as.character(item)),
                    readRDS("EmpiricalData/SB2020_data/SB2020_e3.rds") %>% mutate(item = as.character(item))) %>%
  mutate(oldnew=isold) %>%
  mutate(oldnew = ifelse(oldnew==1,"Old","New")) %>%
  mutate(response = rating)


SB2020Gsq <- NULL
SB2020respfreq <- NULL
for(index in c(1:dim(allfitsSB2020)[1])){

  idn<-allfitsSB2020[index,] %>% .$id
  condn <- allfitsSB2020[index,] %>% .$condition
  modelname <- allfitsSB2020[index,] %>% .$model

  dp<-prep_data(data = dataSB2020 %>%
                  filter(id == idn) %>%
                  filter(condition == condn),freqdat =F,numcategories=6)

  pars <- as.numeric(allfitsSB2020[index,] %>% ungroup() %>%
                       select(names(get_start_par(modelname,numcategories=6))))
  expdata <- predict_frequencies(data_list=dp$datalist,modelname,pars)

  Gsquared <- Gsquaredsingle(modelname,dp,expdata,totitems=130) %>%
    mutate(id = idn,
           condition= condn,
           model = modelname,
           exp = allfitsSB2020[index,] %>% .$exp)

  SB2020Gsq <- SB2020Gsq %>% bind_rows(Gsquared)
  freq <- bind_cols(t(expdata),
                    id = idn,
                    condition = condn,
                    model = modelname,
                    exp = allfitsSB2020[index,] %>% .$exp)

  SB2020respfreq <- SB2020respfreq %>% bind_rows(freq)

}



allfitsSB <- readRDS("EmpiricalFits/SB2021_e12_fits.rds")


dataSB12 <- bind_rows(readRDS("../EmpiricalData/data_SB2022_e1.rds"),
                      readRDS("../EmpiricalData/data_SB2022_e2.rds")) %>%
  mutate(oldnew = ifelse(oldnew==0,"New","Old"))

SB12Gsq <- NULL
SB12respfreq <- NULL
for(index in c(1:dim(allfitsSB)[1])){

  idn<-allfitsSB[index,] %>% .$id
  condn <- allfitsSB[index,] %>% .$condition
  modelname <- allfitsSB[index,] %>% .$model

  dp<-prep_data(data = dataSB12 %>%
                  filter(id == idn) %>%
                  filter(condition == condn),freqdat =F,numcategories=6)

  pars <- as.numeric(allfitsSB[index,] %>% ungroup() %>%
    select(names(get_start_par(modelname,numcategories=6))))
  expdata <- predict_frequencies(data_list=dp$datalist,modelname,pars)

  Gsquared <- Gsquaredsingle(modelname,dp,expdata,totitems=60) %>%
    mutate(id = idn,
           condition= condn,
           model = modelname,
           exp = allfitsSB[index,] %>% .$exp)

  SB12Gsq <- SB12Gsq %>% bind_rows(Gsquared)
  freq <- bind_cols(t(expdata),
                            id = idn,
                              condition = condn,
                              model = modelname,
                              exp = allfitsSB[index,] %>% .$exp)

  SB12respfreq <- SB12respfreq %>% bind_rows(freq)

}


fitsMWW <- list.files("../EmpiricalFits/MWW2007/","MWW",full.names=T)
fitsMWW <- fitsMWW[grepl(paste(models,collapse="|"),fitsMWW)]

allfitsMWW <-  purrr::map(fitsMWW, readRDS) %>%
  bind_rows() %>%
  group_by(exp,id,condition,model) %>%
  filter(!model %in% "GaussianUVSDT_equi") %>%
  filter(objective == min(objective))
dataMWW <- readRDS("../EmpiricalData/MWWData/MWW2007_sixpoint_preprocesseddata.rds")

MWWGsq <- NULL
MWWrespfreq <- NULL

for(index in c(1:dim(allfitsMWW)[1])){

  idn<-allfitsMWW[index,] %>% .$id
  #condn <- allfitsSB[index,] %>% .$condition
  modelname <- allfitsMWW[index,] %>% .$model

  dp<-prep_data(data = dataMWW %>%
                  filter(id == idn) ,freqdat =F,numcategories=6)

  pars <- as.numeric(allfitsMWW[index,] %>% ungroup() %>%
                       select(names(get_start_par(modelname,numcategories=6))))
  expdata <- predict_frequencies(data_list=dp$datalist,modelname,pars)

  Gsquared <- Gsquaredsingle(modelname,dp,expdata,totitems=round(dp$trialsNew)) %>%
    mutate(id = idn,
          # condition= condn,
           model = modelname,
           exp = allfitsMWW[index,] %>% .$exp)

  MWWGsq <- MWWGsq %>% bind_rows(Gsquared)
  freq <- bind_cols(t(expdata),
                           id = idn,
                                              condition = condn,
                                              model = modelname,
                                              exp = allfitsMWW[index,] %>% .$exp)

  MWWrespfreq <- MWWrespfreq %>% bind_rows(freq)


}


bestfitsSB3 <- bind_rows(readRDS("../EmpiricalFits/BestFits/SB2021_e3_modelfits_GumbelEVSDT.rds"),
                         readRDS("../EmpiricalFits/BestFits/SB2021_e3_modelfits_baselineUVSDT.rds"),
                         readRDS("../EmpiricalFits/BestFits/SB2021_e3_modelfits_ExGaussNormEVSDT.rds"))

dataSB3 <- SB2021_e3
SB3Gsq <- NULL
SB3respfreq <- NULL

for(index in c(1:dim(bestfitsSB3)[1])){

  idn<-bestfitsSB3[index,] %>% .$id
  #condn <- allfitsSB[index,] %>% .$condition
  modelname <- bestfitsSB3[index,] %>% .$model

  dp<-prep_data(data = dataSB3 %>%
                  filter(id == idn) ,freqdat =F,numcategories=6)

  pars <- as.numeric(bestfitsSB3[index,] %>% ungroup() %>%
                       select(names(get_start_par(modelname,numcategories=6))))
  expdata <- predict_frequencies(data_list=dp$datalist,modelname,pars)

  Gsquared <- Gsquaredmult(modelname,dp,expdata,totitems=round(dp$trials_LN)) %>%
    mutate(id = idn,
           # condition= condn,
           model = modelname,
           exp = bestfitsSB3[index,] %>% .$exp)

  SB3Gsq <-SB3Gsq %>% bind_rows(Gsquared)
  freq <- bind_cols(t(expdata),
                           id = as.character(idn),
                                            condition = condn,
                                            model = modelname,
                                            exp = bestfitsSB3[index,] %>% .$exp)
  SB3respfreq <- SB3respfreq %>% bind_rows(freq)
}

Gsqs <- bind_rows(SB2020Gsq,SB12Gsq,SB3Gsq %>% mutate(id=as.character(id)),MWWGsq)
Preds <- list(SB2020_e123 = SB2020respfreq,
              SB2021_e12 = SB12respfreq,
              SB2021_e3  = SB3respfreq ,
              MWW2007 = MWWrespfreq)
#saveRDS(Gsqs,file="Gsquared_allexps.rds")
#saveRDS(Preds,file="Preds_allexps.rds")



