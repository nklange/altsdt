
library("tidyverse")

SB2020 <- bind_rows(readRDS("EmpiricalData/SB2020_data/SB2020_e1.rds") %>% mutate(item = as.character(item)),
                    readRDS("EmpiricalData/SB2020_data/SB2020_e2.rds") %>% mutate(item = as.character(item)),
                    readRDS("EmpiricalData/SB2020_data/SB2020_e3.rds") %>% mutate(item = as.character(item))) %>%
  mutate(oldnew=isold) %>%
  mutate(oldnew = ifelse(oldnew==1,"Old","New")) %>%
  mutate(response = rating)



# Model fitting ----------------------------------------------------------------
#
source("FitSimple.R")
models <- c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT")

for (model in models){

  for(subjid in unique(SB2020$id)){


    fullsubj <- NULL
    for(cond in unique(SB2020 %>% filter(id == subjid) %>% .$condition)){

      data <- SB2020 %>% filter(id == subjid) %>% filter(condition==cond)

      fit <- FitSDT(data = data, model = model, rep=50, freqdat = F,numcategories=6) %>%
        mutate(condition = cond)

      fullsubj <- fullsubj %>% bind_rows(fit)
    }

    saveRDS(fullsubj,file = paste0(""))

  }

}

source("FitMultiple.R")



SB2021 <- readRDS("EmpiricalData/data_SB2022_e3.rds") %>%
  mutate(item = as.character(item)) %>%
  mutate(oldnew=isold) %>%
  mutate(oldnew = ifelse(oldnew==1,"Old","New")) %>%
  mutate(response = rating)

# data already has hardcorded:
# oldnew -> only 1 set of new items (low_ev_0), the other set of new items (high_ev_0) coded as 'old'
# condition -> turn 2 actual conditions (low_ev, high_ev) into four by attaching the original old/new
# oldnew_true -> actual old/new coding
# condition_true -> actual condition coding
# old/new and condition coding gets used in FitMultiple.R/prep_data function

for(modelname in c("baselineUVSDT","GumbelEVSDT","ExGaussNormEVSDT")){
  allfits2 <- NULL
  for(idn in unique(SB2021$id)){



    fit <- FitSDT(data = SB2021_e1 %>% filter(id == idn),
                  rep = 25,
                  model = modelname,
                  idn = idn)

    best <- fit %>% filter(objective == min(objective)) %>% .[1,] %>%
      mutate(model = modelname,
             id = idn)

    fitted <- best %>%
      mutate(AIC = 2 * objective + 2 * npar)

    allfits2 <- bind_rows(allfits2,fitted)

  }
  saveRDS(allfits2,file="")
}




