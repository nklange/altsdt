# General information -----------------------------------------

# Correlation UVSDT sigma ratio with empirical SDratio of confidence rating
# Simulations described in SimulateFromModels.R

library(tidyverse)

# Load simulated data -------------------------

# 3 ways of setting decision thresholds for simulation
# for paper: 1. empirical criteria
# 2. uniform distribution of criteria where extremes are given-ish by data (not reported)
# 3. equi-distant criteria where extremes are given-ish by data (not reported)

# Simulation of this is in SimulateFromModels_MWW.R (same principle as
# SimulateFromModels.R but using MWW fits as basis, larger confidence scale,
# standard number of items per experiment)

# _data3: simulated 24point confidence rating scale from empirical data
# _data4_unif: uniform distribution of confidence ratings
# _data5_equi: equidistant criteria
# for paper: _sixpoint: after simulation of _data3, binning to six confidence rating bins (MWW approach)

GetGsquared <- function(LL,observedData){


  temp <- c(observedData[1:6] * log(observedData[1:6]/sum(observedData[1:6])),
            observedData[7:12] * log(observedData[7:12]/sum(observedData[7:12])))
  temp[observedData == 0] <- 0

  GSq = 2*(LL - -sum(temp) ) # Calculate G^2 (better approximator to Chi-sq dist than Pearson Chi-Sq)
  # f$df     <- length(observedData)-length(names(get_start_par(model)))-1            # df model (L&F,p.63) n_observations - n_freeparas - 1
  #f$pGSq   <- 1-pchisq(f$GSq,f$df)       # p value of GSq



  return(GSq)
}

load_files <- function(path,pattern) {

  files <- dir(path, pattern = pattern,full.names = TRUE,recursive = FALSE)
  # search for test data (dataTEST) file in all subdirectories (recursive) of path
  tables <- lapply(files, readRDS)
  do.call(bind_rows, tables)
  # collate all of them into 1 dataframe
}

# Empirical Threshold: Confidence3_bestfits.rds (24pt) "EmpThresh24"
# Uniform Threshold: ConfUnif_bestfits.rds (24pt) "Unif24"
# Equidistant Threshold: ConfEqui_bestfits.rds (24pt) "Equi24"


# FOR PAPER: Empirical Threshold: Confidence3_sixpoint_bestfits.rds (24 -> 6pt) "EmpThresh6"
# Uniform Threshold: ConfUnif_bestfits_sixpoint.rds (24 -> 6pt) "Unif6"
# Equidistant Threshold: ConfEqui_bestfits_sixpoint.rds (24 -> 6pt) "Equi6"


extractSDRatio<-function(ThresholdString){

  if (ThresholdString == "Equi6"){

    datasim <- "data5_equi_sixpoint"
    parsim <- "parametervalues5_equi"
    fitsim <- "ConfEqui_bestfits_sixpoint"
  } else if (ThresholdString == "Equi24"){

    datasim <- "data5_equi"
    parsim <- "parametervalues5_equi"
    fitsim <- "ConfEqui_bestfits"

  } else if (ThresholdString == "Unif6"){

    datasim <- "data4_unif_sixpoint"
    parsim <- "parametervalues4_unif"
    fitsim <- "ConfUnif_bestfits_sixpoint"

  } else if (ThresholdString == "Unif24"){

    datasim <- "data4_unif"
    parsim <- "parametervalues4_unif"
    fitsim <- "ConfUnif_bestfits"

  } else if (ThresholdString == "EmpThresh6"){

    datasim <- "data3_sixpoint"
    parsim <- "parametervalues3"
    fitsim <- "Confidence3_sixpoint_bestfits"

  } else if (ThresholdString == "EmpThresh24"){

    datasim <- "data3"
    parsim <- "parametervalues3"
    fitsim <- "Confidence3_bestfits"

  }

MimicryFits <- NULL
for (gen in c("genGumbelEVSDT","genGaussianUVSDT","genExGaussNormEVSDT")){

    fits <- readRDS(paste0("Simulations/ConfidenceSimulation/",fitsim,".rds")) %>%
    filter(genmodel == substr(gen,start=4,stop=nchar(gen))) %>% arrange(id,model)

  if(ThresholdString == c("EmpThresh6")){
  data <- readRDS(paste0("Simulations/ConfidenceSimulation/simulate_",datasim,".rds")) %>%
    arrange(id) %>%  filter(genmodel == substr(gen,start=4,stop=nchar(gen)))
  } else {
    data <- readRDS(paste0("Simulations/ConfidenceSimulation/simulate_",gen,"_",datasim,".rds")) %>%
      arrange(id) %>% mutate(response = as.numeric(as.character(response)))
  }

  LL <- fits %>% mutate(genmodel = substr(gen,start=4,stop=nchar(gen))) %>%
    group_by(genmodel,id,model) %>%
    group_nest(keep=T,.key="fit")

  observedData <- data %>% group_by(genmodel,id)  %>% group_nest(keep=T)%>%
    dplyr::slice(rep(1:n(),each=1))

  parameters <- readRDS(paste0("Simulations/ConfidenceSimulation/simulate_",gen,"_",parsim,".rds")) %>%
    pivot_wider(values_from="value",names_from="parameter") %>%  group_by(genmodel,parid)  %>%
    group_nest(keep=T)%>%dplyr::slice(rep(1:n(),each=1)) %>% rename("genparameters" = data)

  getGsquaredvals <- LL %>% bind_cols(observedData %>% ungroup() %>%  dplyr::select(data)) %>%
    bind_cols(parameters %>% ungroup() %>% dplyr::select(genparameters)) %>%
    mutate(Gsq = map2(.x = fit,.y = data, .f = ~GetGsquared(.x$objective,.y$Freq))) %>%
    mutate(objective = map(.x = fit,.f = ~mean(.x$objective))) %>%
    mutate(Dev = map(.x = fit,.f = ~mean(2*.x$objective))) %>%
    mutate(AIC = map(.x = fit,.f = ~mean(.x$AIC)))


  MimicryFits <- MimicryFits %>% bind_rows(getGsquaredvals)

}

MimicryFits

}

choosedata <- function(ThresholdString){


  MimicryFits <- extractSDRatio(ThresholdString)

  detSDratio <- MimicryFits %>% .$data %>% bind_rows()

  Data_SDratio <- NULL
  for (i in unique(detSDratio$id)){

    une <- detSDratio %>% filter(id == i)

    new <- sd(rep(une %>% filter(oldnew == "New") %>% .$response,
                  une %>% filter(oldnew == "New") %>% .$Freq))
    old <- sd(rep(une %>% filter(oldnew == "Old") %>% .$response,
                  une %>% filter(oldnew == "Old") %>% .$Freq))

    out <- tibble(id = i,
                  genmodel = unique(une$genmodel),
                  SDlure = new,
                  SDtarget = old) %>%
      mutate(SDratio = new/old)

    Data_SDratio <- Data_SDratio %>% bind_rows(out)

  }

  detsigmaratio <- MimicryFits %>% .$fit %>% bind_rows() %>% filter(model == "GaussianUVSDT")


  SD_Ratio <- Data_SDratio %>% arrange(genmodel,id)
  PlotRatios_sim <- detsigmaratio %>% dplyr::select(genmodel,id,sigo)  %>% arrange(genmodel,id) %>%
    mutate(sigmaratio = 1/sigo,
           stdratio = SD_Ratio$SDratio) %>%
    filter(stdratio != "Inf")



  PlotRatios <- bind_rows(PlotRatios_sim) %>%
    mutate(genmodel = factor(genmodel, levels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                             labels=c("Gen: U-V Gaussian","Gen: Gumbel","Gen: Ex-Gaussian"))) %>%
    mutate(designatealpha = "low")

  lim <- 1.5
  cortext <- tibble(genmodel=c("Gen: U-V Gaussian","Gen: Gumbel","Gen: Ex-Gaussian"),
                    corrlim = c(
                                cor(PlotRatios %>% filter(stdratio < lim & sigmaratio < lim) %>%
                                      filter(genmodel=="Gen: U-V Gaussian") %>%
                                      .$stdratio,
                                    PlotRatios %>% filter(stdratio < lim & sigmaratio < lim) %>%
                                      filter(genmodel=="Gen: U-V Gaussian") %>%
                                      .$sigmaratio),
                                cor(PlotRatios %>% filter(stdratio < lim & sigmaratio < lim) %>%
                                      filter(genmodel=="Gen: Gumbel") %>%
                                      .$stdratio,
                                    PlotRatios %>% filter(stdratio < lim & sigmaratio < lim) %>%
                                      filter(genmodel=="Gen: Gumbel") %>%
                                      .$sigmaratio),
                                cor(PlotRatios %>% filter(stdratio < lim & sigmaratio < lim) %>%
                                      filter(genmodel=="Gen: Ex-Gaussian") %>%
                                      .$stdratio,
                                    PlotRatios %>% filter(stdratio < lim & sigmaratio < lim) %>%
                                      filter(genmodel=="Gen: Ex-Gaussian") %>%
                                      .$sigmaratio)),
                    corr = c(
                             cor(PlotRatios %>%
                                   filter(genmodel=="Gen: U-V Gaussian") %>%
                                   .$stdratio,
                                 PlotRatios %>%
                                   filter(genmodel=="Gen: U-V Gaussian") %>%
                                   .$sigmaratio),
                             cor(PlotRatios %>%
                                   filter(genmodel=="Gen: Gumbel") %>%
                                   .$stdratio,
                                 PlotRatios %>%
                                   filter(genmodel=="Gen: Gumbel") %>%
                                   .$sigmaratio),
                             cor(PlotRatios %>%
                                   filter(genmodel=="Gen: Ex-Gaussian") %>%
                                   .$stdratio,
                                 PlotRatios %>%
                                   filter(genmodel=="Gen: Ex-Gaussian") %>%
                                   .$sigmaratio))) %>%
    mutate(designatealpha = NA) %>%
    mutate(genmodel = factor(genmodel,
                             levels=c("Gen: U-V Gaussian","Gen: Gumbel","Gen: Ex-Gaussian")))




  meantext <- tibble(genmodel=c("Gen: U-V Gaussian","Gen: Gumbel","Gen: Ex-Gaussian"),
                     meanstd = c(
                                 mean(PlotRatios %>%
                                        filter(genmodel=="Gen: U-V Gaussian") %>%
                                        .$stdratio),
                                 mean(PlotRatios %>%
                                        filter(genmodel=="Gen: Gumbel") %>%
                                        .$stdratio),
                                 mean(PlotRatios %>%
                                        filter(genmodel=="Gen: Ex-Gaussian") %>%
                                        .$stdratio)),
                     meansigma = c(
                                   mean(PlotRatios %>%
                                          filter(genmodel=="Gen: U-V Gaussian") %>%
                                          .$sigmaratio),
                                   mean(PlotRatios %>%
                                          filter(genmodel=="Gen: Gumbel") %>%
                                          .$sigmaratio),
                                   mean(PlotRatios %>%
                                          filter(genmodel=="Gen: Ex-Gaussian") %>%
                                          .$sigmaratio))) %>%
    mutate(designatealpha = NA) %>%
    mutate(genmodel = factor(genmodel,
                             levels=c("Gen: U-V Gaussian","Gen: Gumbel","Gen: Ex-Gaussian")))

  return(list(PlotRatios,cortext,meantext))
}

PlotRatios_EmpThresh <- choosedata("EmpThresh6")
#PlotRatios_Unif<-choosedata("Unif6")
#PlotRatios_Equi<-choosedata("Equi6")


# MWW2007 ------------

dat_MWW <- readRDS("EmpiricalData/MWWData/MWW2007_sixpoint_preprocesseddata.rds")

datSD_MWW <- dat_MWW %>% group_by(exp,id,oldnew) %>% summarize(sdr = sd(response)) %>%
  group_by(exp,id) %>% mutate(stdratio = sdr[1]/sdr[2]) %>% filter(oldnew=="Old") %>%
  select(-oldnew,-sdr)


## U-V Gaussian Fits to empirical data

fits_MWW <- load_files("EmpiricalFits/MWW2007/","MWW2007")

bestfit_MWW <- fits_MWW %>% filter(model=="GaussianUVSDT") %>%
  # mutate(objective = ifelse(objective < 0, objective * -1, objective)) %>%
  mutate(AIC = 2 * npar + 2 * objective) %>%
  group_by(model,id,condition) %>%
  arrange(objective) %>%
  dplyr::slice(1) %>% ungroup() %>%
  select(exp,id,sigo) %>%
  mutate(sigmaratio= 1/sigo)

empSD_MWW <- inner_join(bestfit_MWW,datSD_MWW)%>%
  filter(stdratio != "Inf") %>%
  mutate(genmodel = "MWW2007",
         designatealpha="VHigh")



corrMWW <- tibble(genmodel = "MWW2007",
       corrlim = cor(empSD_MWW %>% filter(sigmaratio < 1.5 & stdratio < 1.5) %>% .$sigmaratio,
                     empSD_MWW %>% filter(sigmaratio < 1.5 & stdratio < 1.5) %>% .$stdratio),
       corr = cor(empSD_MWW %>% .$sigmaratio,
                  empSD_MWW %>% .$stdratio),
       designatealpha = NA)

meanMWW <- tibble(genmodel = "MWW2007",
                     meanstd = mean(empSD_MWW %>% .$stdratio),
                     meansigma = mean(empSD_MWW %>% .$sigmaratio),
                     designatealpha=NA)

## U-V Gaussian fits to simulated data  (only show for MWW not SB2021 for paper)

datacorr <-bind_rows(PlotRatios_EmpThresh[[1]] %>% filter(genmodel != "SB2021"),
                     empSD_MWW) %>%
  mutate(genmodel = factor(genmodel,levels= c("MWW2007","Gen: U-V Gaussian", "Gen: Gumbel", "Gen: Ex-Gaussian"),
                           labels = c("MWW (2007)","Generating model\nU-V Gaussian", "Generating model\nGumbel", "Generating model\nEx-Gaussian")))

datacorr_text <- bind_rows(corrMWW,
                          PlotRatios_EmpThresh[[2]] %>% filter(genmodel != "SB2021"),
) %>%
  mutate(genmodel = factor(genmodel,levels= c("MWW2007","Gen: U-V Gaussian", "Gen: Gumbel", "Gen: Ex-Gaussian"),
                           labels = c("MWW (2007)","Generating model\nU-V Gaussian", "Generating model\nGumbel", "Generating model\nEx-Gaussian")))


datamean_text <- bind_rows(PlotRatios_EmpThresh[[3]] %>% filter(genmodel != "SB2021"),
                           meanMWW)%>%
  mutate(genmodel = factor(genmodel,levels= c("MWW2007","Gen: U-V Gaussian", "Gen: Gumbel", "Gen: Ex-Gaussian"),
                           labels = c("MWW (2007)","Generating model\nU-V Gaussian", "Generating model\nGumbel", "Generating model\nEx-Gaussian")))

# Make graph: Figure for  paper -----------------

ggplot(datacorr,
       aes(x=stdratio,y=sigmaratio,color=genmodel,alpha=designatealpha,size=designatealpha))+
  # annotate("text", x = 0, y = 4, aes(label = lm_eqn_poly1(testcorr %>% filter(model == "GaussianUVSDT"))),
  #          parse = TRUE, color="black",hjust = 0)+
  geom_point() +
 # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")+
  #geom_hex(bins=70) +
 # scale_fill_continuous(type = "viridis") +
  #geom_smooth()+

  scale_alpha_discrete(range=c(.1,.5))+
   scale_size_discrete(range=c(0.5,2))+
  geom_text(data = datacorr_text,aes(x = 0,y=1.5,label=paste0("graph: r = ", stringr::str_remove(round(corrlim,2), "^0+")),
                                                                                  color=genmodel),hjust = 0)+
  geom_text(data = datacorr_text,aes(x = 0,y=1.35,label=paste0("total: r = ", stringr::str_remove(round(corr,2), "^0+")),
                                                                                  color=genmodel),hjust = 0)+
  geom_text(data = datamean_text,aes(x = 0.75,y=0,label=paste0("M = ", round(meanstd,2)),
                                                                                  color=genmodel),hjust = 0.5)+
  geom_text(data = datamean_text,aes(x = 0,y=0.75,label=paste0("M = ", round(meansigma,2)),
                                                                                  color=genmodel),angle = 90,hjust=0.5)+
  scale_color_manual(values = c("#616161","#E69F00","#009E73","#0072B2"))+
  scale_x_continuous(name=expression(paste(SD[n]/SD[o]," (Confidence Ratings)")),limits=c(0,1.5))+
  scale_y_continuous(name=expression(paste(sigma[n]/sigma[o], " (U-V Gaussian estimates)")),limits=c(0,1.5))+

  facet_grid(.~genmodel)+
  theme(
    axis.text = element_text(size=12),

    #legend.position = legendpos,
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = "transparent",colour="black"),
    strip.placement = "outside",
    legend.position = "none",
    strip.background = element_rect(fill = "transparent",
                                    colour="transparent"),
    strip.text = element_text(size = 12))

ggsave(paste0("Figures/StandardDeviations_sixpoint.png"), units="cm", width=22, height=7, dpi=600)

