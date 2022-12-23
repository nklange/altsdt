# General information -----------------------------------------------------------

# Identify best fits of all fits
# integrate Gsquared: Gsquared is calculated in AnalyseGsquared, collected in:
# EmpiricalFits/Gsquared_allexps.rds
# make visualizations


library(tidyverse)
library(waffle)
library(cowplot)


models <- c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT","baselineUVSDT")

Gsquared <- readRDS("EmpiricalFits/Gsquared_allexps.rds")

# SB2020, E1, E2, E3 ----------------------------------------------------------

allfitsSB2020 <- readRDS("EmpiricalFits/SB2020_e123_fits.rds")
bestfitSB2020 <- allfitsSB2020 %>%
  mutate(AIC = 2 * objective + 2 *npar) %>%
  group_by(exp,id,condition) %>%
  mutate(deltaAIC = AIC - min(AIC)) %>%
  mutate(exp = case_when(exp=="SB2020_e1"~"Exp 1",
                         exp=="SB2020_e2" ~"Exp 2",
                         exp=="SB2020_e3"~"Exp 3")) %>%
  mutate(condition = case_when(condition == "FixedStudy" ~ "Fixed Study",
                               condition == "VariableStudy" ~ "Variable Study",
                               condition == "FixedIntervalNback" ~ "Fixed N-back",
                               condition == "VariableIntervalNback" ~ "Variable N-back",
                               condition == "FixedFrequency" ~ "Fixed Word Frequency",
                               condition == "VariableFrequency" ~ "Variable Word Frequency")) %>%
  #mutate(exp = ifelse(exp == "SB2021_e1","Exp 1", "Exp 2")) %>%
  mutate(exp = paste0(exp,", ",condition)) %>%
  mutate(model = factor(model,levels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                        labels = c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"))) %>%
  select(exp,id,model,AIC,deltaAIC)
# mutate(condition = factor(condition,




GsqSBB2020 <- Gsquared %>%
  filter(exp %in% c("SB2020_e1","SB2020_e2","SB2020_e3")) %>%

  mutate(exp = case_when(exp=="SB2020_e1"~"Exp 1",
                         exp=="SB2020_e2" ~"Exp 2",
                         exp=="SB2020_e3"~"Exp 3")) %>%
  mutate(condition = case_when(condition == "FixedStudy" ~ "Fixed Study",
                               condition == "VariableStudy" ~ "Variable Study",
                               condition == "FixedIntervalNback" ~ "Fixed N-back",
                               condition == "VariableIntervalNback" ~ "Variable N-back",
                               condition == "FixedFrequency" ~ "Fixed Word Frequency",
                               condition == "VariableFrequency" ~ "Variable Word Frequency")) %>%
  #mutate(exp = ifelse(exp == "SB2021_e1","Exp 1", "Exp 2")) %>%
  mutate(exp = paste0(exp,", ",condition)) %>%
  group_by(exp,model) %>%
  summarize(values= length(pGSq[pGSq > .05])) %>%
  mutate(percentage = paste0(round(values/40 * 100,0),"%")) %>% #hardcoded number ppts
  mutate(printval = paste0(values,"/",40,"\n",percentage))

bestfitSB2020_best <- bestfitSB2020 %>%
  group_by(exp,id) %>%
  filter(AIC == min(AIC)) %>%
  # mutate(exp = paste0(exp,"_",condition)) %>%
  group_by(exp,model) %>%
  summarize(values=length(model)) %>%
  # mutate(exp = factor(exp,labels=c("D (2007), Exp 1a","D (2007), Exp 1b"),
  #                     levels= c("Exp 1a","Exp 1b"))) %>%
  group_by(exp) %>%
  mutate(total = sum(values)) %>%
  mutate(percentage = paste0(round(values/total * 100,0),"%")) %>%
  mutate(printval = paste0(values,"/",total,"\n",percentage))


meanAICSB2020 <- bestfitSB2020 %>%
  #mutate(exp = paste0(exp,"_",condition)) %>%
  group_by(exp,model) %>% summarize(meanAIC = mean(AIC)) %>%
  arrange(meanAIC) %>%
  group_by(exp) %>%
  mutate(deltaAIC = meanAIC - min(meanAIC))


# SB2021 E1, E2 -----------------------------------------------------------------

allfitsSB <- readRDS("EmpiricalFits/SB2021_e12_fits.rds")


bestfitSB12 <- allfitsSB %>%
  mutate(AIC = 2 * objective + 2 *npar) %>%
  group_by(exp,id,condition) %>%
  mutate(deltaAIC = AIC - min(AIC)) %>%
  mutate(exp = ifelse(exp == "SB2021_e1","Exp 1", "Exp 2")) %>%
  mutate(exp = paste0(exp,", ",condition)) %>%
  mutate(model = factor(model,levels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                        labels = c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"))) %>%
  select(exp,id,model,AIC,deltaAIC)
  # mutate(condition = factor(condition,
  #                           levels = c("A","B","C","D"),
  #                           labels = c("High strength\nHigh variability\n",
  #                                      "High strength\nLow variability\n",
  #                                      "Low strength\nHigh variability\n",
  #                                      "Low strength\nLow variability\n")))

GsqSBB <- Gsquared %>%
  filter(exp %in% c("SB2021_e1","SB2021_e2")) %>%

  mutate(exp = ifelse(exp == "SB2021_e1","Exp 1", "Exp 2")) %>%
  mutate(exp = paste0(exp,", ",condition)) %>%
  group_by(exp,model) %>%
  summarize(values= length(pGSq[pGSq > .05])) %>%
  mutate(percentage = paste0(round(values/64 * 100,0),"%")) %>%
  mutate(printval = paste0(values,"/",64,"\n",percentage))

bestfitSB12_best <- bestfitSB12 %>%
  group_by(exp,id) %>%
  filter(AIC == min(AIC)) %>%
 # mutate(exp = paste0(exp,"_",condition)) %>%
  group_by(exp,model) %>%
  summarize(values=length(model)) %>%
  # mutate(exp = factor(exp,labels=c("D (2007), Exp 1a","D (2007), Exp 1b"),
  #                     levels= c("Exp 1a","Exp 1b"))) %>%
  group_by(exp) %>%
  mutate(total = sum(values)) %>%
  mutate(percentage = paste0(round(values/total * 100,0),"%")) %>%
  mutate(printval = paste0(values,"/",total,"\n",percentage))


meanAICSB12 <- bestfitSB12 %>%
  #mutate(exp = paste0(exp,"_",condition)) %>%
  group_by(exp,model) %>% summarize(meanAIC = mean(AIC)) %>%
  arrange(meanAIC) %>%
  group_by(exp) %>%
  mutate(deltaAIC = meanAIC - min(meanAIC))


# MWW (2007), E1, E2, E3 -----------------------------------------

fitsMWW <- list.files("EmpiricalFits/MWW2007/","MWW",full.names=T)

bestfitsMWW <- purrr::map(fitsMWW, readRDS) %>%

  bind_rows() %>%
  filter(model %in% c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT")) %>%
  group_by(exp,id,model) %>%
  filter(objective == min(objective)) %>%
  mutate(AIC = 2 * objective + 2 *npar) %>%
  group_by(exp,id) %>%
  filter(!model %in% "GaussianUVSDT_equi") %>%
  mutate(model = factor(model,levels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                        labels = c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"))) %>%
  mutate(exp = factor(exp,levels=c("MWW2007_e1","MWW2007_e2","MWW2007_e3"),
                      labels= c("Exp 1","Exp 2","Exp 3"))) %>%
  group_by(exp,id) %>%
  mutate(deltaAIC = AIC - min(AIC))

GsqMWW  <- Gsquared %>%
  filter(exp %in% c("MWW2007_e1","MWW2007_e2","MWW2007_e3")) %>%

  mutate(exp = factor(exp,levels=c("MWW2007_e1","MWW2007_e2","MWW2007_e3"),
                      labels= c("Exp 1","Exp 2","Exp 3"))) %>%
  group_by(exp,model) %>%
  summarize(values= length(pGSq[pGSq > .05]),
            total = length(pGSq)) %>%
  mutate(percentage = paste0(round(values/total * 100,0),"%")) %>%
  mutate(printval = paste0(values,"/",total,"\n",percentage))


bestfitsMWW_best <- bestfitsMWW %>%
  group_by(exp,id) %>%
  filter(AIC == min(AIC)) %>%
  group_by(exp,model) %>%
  summarize(values=length(model)) %>%
  # mutate(exp = factor(exp,labels=c("D (2007), Exp 1a","D (2007), Exp 1b"),
  #                     levels= c("Exp 1a","Exp 1b"))) %>%
  group_by(exp) %>%
  mutate(total = sum(values)) %>%
  mutate(percentage = paste0(round(values/total * 100,0),"%")) %>%
  mutate(printval = paste0(values,"/",total,"\n",percentage))

meanAICMWW <- bestfitsMWW %>%  group_by(exp,model) %>% summarize(meanAIC = mean(AIC)) %>%
  arrange(exp,meanAIC) %>%
  group_by(exp) %>%
  mutate(deltaAIC = meanAIC - min(meanAIC))

# SB2021 E3 -------------------------------------------------------------------

bestfitsSB3 <-  bind_rows(readRDS("EmpiricalFits/SB2021_e3_modelfits_GumbelEVSDT.rds"),
                          readRDS("EmpiricalFits/SB2021_e3_modelfits_baselineUVSDT.rds"),
                          readRDS("EmpiricalFits/SB2021_e3_modelfits_ExGaussNormEVSDT.rds")) %>%
  mutate(model = factor(model,levels=c("baselineUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                        labels = c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"))) %>%
  mutate(AIC = 2 * objective + 2 *npar) %>%
  mutate(exp = "Exp 3") %>%
  group_by(exp,id) %>%
  mutate(deltaAIC = AIC - min(AIC))



GsqSB3 <-  Gsquared %>%
  filter(exp %in% c("SB2022_e3")) %>%

  mutate(model = factor(model,levels=c("baselineUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                        labels = c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"))) %>%
  mutate(exp = "Exp 3") %>%
  group_by(exp,model) %>%
  summarize(values= length(pGSq[pGSq > .05]),
            total = length(pGSq)) %>%
  mutate(percentage = paste0(round(values/total * 100,0),"%")) %>%
  mutate(printval = paste0(values,"/",total,"\n",percentage))


bestfitsSB3_best <- bestfitsSB3 %>%
  group_by(exp,id) %>%
  filter(AIC == min(AIC)) %>%
  group_by(exp,model) %>%
  summarize(values=length(model)) %>%
  # mutate(exp = factor(exp,labels=c("D (2007), Exp 1a","D (2007), Exp 1b"),
  #                     levels= c("Exp 1a","Exp 1b"))) %>%
  group_by(exp) %>%
  mutate(total = sum(values)) %>%
  mutate(percentage = paste0(round(values/total * 100,0),"%")) %>%
  mutate(printval = paste0(values,"/",total,"\n",percentage))

meanAICSB3<- bestfitsSB3 %>%  group_by(exp,model) %>% summarize(meanAIC = mean(AIC)) %>%
  arrange(meanAIC) %>%
  group_by(exp) %>%
  mutate(deltaAIC = meanAIC - min(meanAIC))

# Make Empirical Fits plots ---------------------------------------------------
# Make waffle plot part of plot


waffleSB3 <- ggplot(bestfitsSB3_best,
       aes(fill=model,values=values))+
  geom_waffle(n_rows = 5, size = 0.33, colour = "white",keep=T)+
  # scale_fill_manual(values=c("#E69F00","#009E73","#56B4E9","#CC79A7"))+
  scale_fill_manual(values=c("#E69F00","#009E73","#56B4E9"))+
  # scale_fill_manual(values=c("#E69F00","#CC79A7"))+
  scale_y_continuous(name="Winning")+
  # facet_grid(condition~exp)+
  ggtitle("Spanton & Berry (2022)")+
  coord_equal()+
  facet_wrap(exp~.)+
  theme( axis.ticks= element_blank(),
         axis.text =element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_text(size=12,face="plain",hjust=0.5),
         panel.background = element_blank(),
         legend.position = "none",
         title = element_text(size=14,hjust=0.5,face="bold"),
         #legend.title=element_blank(),
         strip.background.x = element_rect(fill = "white"),
         strip.background.y = element_rect(fill = "white"),
         strip.placement = "outside",
         strip.text = element_text(size=12))

# Make deltaAIC part of plot



deltaAICSB3 <- ggplot(bestfitsSB3)+
  coord_cartesian(ylim=c(-3,9))+
  geom_hline(yintercept=0,size=0.5)+

  geom_jitter(data=bestfitsSB3,
              aes(x=model,y=deltaAIC,color=model),size=3,alpha=0.3,height = 0,width=0.2)+
  scale_color_manual(values=c("#E69F00","#009E73","#56B4E9"))+
  scale_y_continuous(name=expression(paste(Delta, " AIC")),
                     breaks = c(-2,0,2,4,6,8),
                     labels = c(expression(paste(Delta,"AIC = 0")),
                                0,2,4,6,
                                expression(paste(G^2,">.05"))))+
  geom_point(data = meanAICSB3,
             aes(x=model,y=deltaAIC),
             size=5,shape=15,color="white",alpha=.7)+
  geom_point(data = meanAICSB3,aes(x=model,y=deltaAIC),color="black",
             size=4,shape=7)+
  #  annotate("text",x=1,y=-1,label=expression(paste(Delta,"AIC = 0")))+
  geom_text(data = bestfitsSB3_best,size=3,
            aes(x=model,y=-2,label=printval))+
  #  annotate("text",x=1,y=20,label=expression(paste(G^2,">.05")))+
  geom_text(data = GsqSB3,size=3,
            aes(x=model,y=8,label=printval))+
  scale_x_discrete(labels=c("U-V Gaussian","Gumbel","Ex-Gaussian"))+

  facet_wrap(exp~.,nrow=1)+
  ggtitle("Spanton & Berry (2022)")+
  theme(axis.text.x = element_text(size=9),
        axis.text.y = element_text(size = c(10,9,9,9,9,10)),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_line(color=c("transparent","black","black","black","black","transparent")),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill="white"),
        strip.placement = "outside")

# Make deltaAIC part of plot


  deltaAICSB12 <- ggplot(bestfitSB12)+
  coord_cartesian(ylim=c(-3,9))+
  geom_hline(yintercept=0,size=0.5)+
  geom_jitter(data=bestfitSB12,
              aes(x=model,y=deltaAIC,color=model),size=3,alpha=0.3,height = 0,width=0.2)+
  scale_color_manual(values=c("#E69F00","#009E73","#56B4E9"))+
  scale_y_continuous(name=expression(paste(Delta, " AIC")),
                     breaks = c(-2,0,2,4,6,8),
                     labels = c(expression(paste(Delta,"AIC = 0")),
                                0,2,4,6,
                                expression(paste(G^2,">.05"))))+
  geom_point(data = meanAICSB12,
             aes(x=model,y=deltaAIC),
             size=5,shape=15,color="white",alpha=.5)+
  geom_point(data = meanAICSB12,
             aes(x=model,y=deltaAIC),color="black",
             size=4,shape=7)+
  #  annotate("text",x=1,y=-1,label=expression(paste(Delta,"AIC = 0")))+
  geom_text(data = bestfitSB12_best,size=3,
            aes(x=model,y=-2,label=printval))+
  #  annotate("text",x=1,y=20,label=expression(paste(G^2,">.05")))+
  geom_text(data = GsqSBB,size=3,
            aes(x=model,y=8,label=printval))+
  scale_x_discrete(labels=c("U-V Gaussian","Gumbel","Ex-Gaussian"))+

  facet_wrap(exp~.,nrow=2,dir="h")+
  ggtitle("Spanton & Berry (2022)")+
  theme(axis.text.x = element_text(size=9),
        axis.text.y = element_text(size = c(10,9,9,9,9,10)),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_line(color=c("transparent","black","black","black","black","transparent")),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill="white"),
        strip.placement = "outside")


# Make deltaAIC part of plot

deltaAICMWW <- ggplot(bestfitsMWW)+

  coord_cartesian(ylim=c(-3,9))+
  geom_hline(yintercept=0,size=0.5)+
  geom_jitter(data=bestfitsMWW,
              aes(x=model,y=deltaAIC,color=model),size=3,alpha=0.3,height = 0,width=0.2)+
  scale_color_manual(values=c("#E69F00","#009E73","#56B4E9"))+
  scale_y_continuous(name=expression(paste(Delta, " AIC")),
                     breaks = c(-2,0,2,4,6,8),
                     labels = c(expression(paste(Delta,"AIC = 0")),
                                0,2,4,6,
                                expression(paste(G^2,">.05"))))+
  geom_point(data = meanAICMWW ,
             aes(x=model,y=deltaAIC),
             size=5,shape=15,color="white",alpha=.7)+
  geom_point(data = meanAICMWW,aes(x=model,y=deltaAIC),color="black",
             size=4,shape=7)+

  #  annotate("text",x=1,y=-1,label=expression(paste(Delta,"AIC = 0")))+
  geom_text(data = bestfitsMWW_best,size=3,
            aes(x=model,y=-2,label=printval))+
  #  annotate("text",x=1,y=20,label=expression(paste(G^2,">.05")))+
  geom_text(data = GsqMWW,size=3,
            aes(x=model,y=8,label=printval))+
  scale_x_discrete(labels=c("U-V Gaussian","Gumbel","Ex-Gaussian"))+

  facet_wrap(exp~.,nrow=1)+
  ggtitle("Mickes, Wais & Wixted (2007)")+
  theme(axis.text.x = element_text(size=9),
        axis.text.y = element_text(size = c(10,9,9,9,9,10)),
        axis.ticks.x = element_blank(),
        axis.ticks.y=element_line(color=c("transparent","black","black","black","black","transparent")),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill="white"),
        strip.placement = "outside")






SB3MWW <- plot_grid(deltaAICSB3, deltaAICMWW,rel_widths=c(.3,.7),nrow=1,align="v")

plot_grid(deltaAICSB12,SB3MWW,nrow=2,rel_heights=c(1.8,1))



