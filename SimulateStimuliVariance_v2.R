# General informations -------------------

# Simulate effect of changing variance of item effects on ROCs
# for data in SB2021_e1,_e2, ie. between-subject conditions

library(tidyverse)

sim <- function(m_high,m_low,sd_r2_high,sd_r2_low){

  mu_low = m_low
  mu_high = m_high
 # r1 = rnorm(1,0,sd_r1)

# r2_ohh = rnorm(1000,0,sd_r2_high)
# r2_ohl = rnorm(1000,0,sd_r2_low)
# r2_ohh = rnorm(1000,0,sd_r2_high)
# r2_ohl = rnorm(1000,0,sd_r2_low)
# r2_new_high = rnorm(1000,0,sd_r2_high)
# r2_new_low = rnorm(1000,0,sd_r2_low)


r2_ohh = rnorm(1000,m_high,sd_r2_high)
r2_ohl = rnorm(1000,m_low,sd_r2_low)
r2_ohh = rnorm(1000,m_high,sd_r2_high)
r2_ohl = rnorm(1000,m_low,sd_r2_low)
r2_new_high = rnorm(1000,0,sd_r2_high)
r2_new_low = rnorm(1000,0,sd_r2_low)

baseline <- rnorm(1000,0,1)

  crit_high <- c(m_high - 1,m_high-0.5,m_high, m_high+0.5, m_high+1)
  crit_low <- c(m_low - 1.5 ,m_low-0.75,m_low, m_low+0.75, m_low+1.5)

  return(list(item_ohh = baseline + rnorm(1000,m_high,sd_r2_high),
              item_ohl = baseline + rnorm(1000,m_high,sd_r2_low),
              item_olh = baseline + rnorm(1000,m_low,sd_r2_high),
              item_oll = baseline + rnorm(1000,m_low,sd_r2_low),
              item_nhh = baseline + rnorm(1000,0,sd_r2_high),
              item_nhl = baseline + rnorm(1000,0,sd_r2_low),
              item_nlh = baseline + rnorm(1000,0,sd_r2_high),
              item_nll = baseline + rnorm(1000,0,sd_r2_low),
              crit_high  = crit_low,
              crit_low = crit_low))
}
makerating <- function(crit,testitem){

  rating<-NULL
  for(i in seq_along(testitem)){
    rating[i] <- which(sort(c(crit,testitem[i])) ==testitem[i])
  }

  as_tibble(rating) %>%
    group_by(value) %>%
    summarize(Freq = length(value)) %>%
    mutate(value = factor(value,levels = c(1:6))) %>%
    complete(value,fill = list(Freq = 0)) %>%
    .$Freq

}

makeFreq <- function(test){

    # pO_hh <- pO_hl <-
    # pO_lh <- pO_ll <-
    # pN_hh <- pN_hl <-
    # pN_lh <- pN_ll <- matrix(ncol=6,nrow=length(test$item_ohh))

  #   crits <- c(-Inf,test$crit,Inf)
  # for (i in 1:6){
  #
  #   # pO_hh[,i] <- pnorm(crits[i+1],mean=test$item_ohh,sd=1)-
  #   #   pnorm(crits[i],mean=test$item_ohh,sd=1)
  #   pO_hl[,i] <- pnorm(crits[i+1],mean=test$item_ohl,sd=1)-
  #     pnorm(crits[i],mean=test$item_ohl,sd=1)
    # pO_lh[,i] <- pnorm(test$crit[i+1],mean=test$item_olh,sd=1)-
    #   pnorm(test$crit[i],mean=test$item_olh,sd=1)
    # pO_ll[,i] <- pnorm(test$crit[i+1],mean=test$item_oll,sd=1)-
    #   pnorm(test$crit[i],mean=test$item_oll,sd=1)
    # pN_hh[,i] <- pnorm(test$crit[i+1],mean=test$item_nhh,sd=1)-
    #   pnorm(test$crit[i],mean=test$item_nhh,sd=1)
    # pN_hl[,i] <- pnorm(test$crit[i+1],mean=test$item_nhl,sd=1)-
    #   pnorm(test$crit[i],mean=test$item_nhl,sd=1)
    # pN_lh[,i] <- pnorm(test$crit[i+1],mean=test$item_nlh,sd=1)-
    #   pnorm(test$crit[i],mean=test$item_nlh,sd=1)
    # pN_ll[,i] <- pnorm(test$crit[i+1],mean=test$item_nll,sd=1)-
    #   pnorm(test$crit[i],mean=test$item_nll,sd=1)
#   }



  data.frame(Freq=c(makerating(test$crit_high,test$item_ohh),
                    makerating(test$crit_high,test$item_ohl),
                    makerating(test$crit_low,test$item_olh),
                    makerating(test$crit_high,test$item_nhh),
                    makerating(test$crit_low,test$item_oll),
                    makerating(test$crit_high,test$item_nhl),
                    makerating(test$crit_low,test$item_nlh),
                    makerating(test$crit_low,test$item_nll)),
             condition_strength = rep(rep(c("high","high","low","low"),each=6),2),
             condition_var = rep(rep(c("high","low"),each=6),4),
             condition = rep(rep(c("A","B","C","D"),each=6),2),
             isold = rep(c(1,0),each=24),
             rating = rep(c(1:6),8)
  )

}



exp <- NULL

for(i in c(1:200)){
  simdata <- sim(2,1,2.5,1)
  dfFreq <- makeFreq(simdata) %>% mutate(id = paste0("simSB2022_e1_",i))
  exp <- bind_rows(exp,dfFreq)

}





makeROC <- function(exp){




  makerocdata <- exp %>%
    group_by(id,condition,isold,rating) %>%
    mutate(cumresp = Freq) %>%
    mutate(rating = factor(rating, levels = c(max(exp$rating):1))) %>%
    arrange(id,condition,isold,rating) %>%
    mutate(isold = factor(isold,levels=c(0,1),labels=c("INew","IOld"))) %>%
    group_by(id,condition,isold) %>%
    complete(rating,fill = list(cumresp = 0)) %>%
    summarize(cumprop = cumsum(cumresp/sum(cumresp)),
              rating = rating) %>%
    ungroup() %>%
    pivot_wider(names_from="isold",values_from="cumprop",
                id_cols = c("id","condition","rating")) %>%
    filter(rating != 1)


  roc <- ggplot(makerocdata,aes(x = INew,y=IOld,group=condition,shape=condition,
                                text=condition,fill=condition,color=condition))+
    #facet_wrap(id~.,ncol=6)+
    geom_line(size=1)+
    geom_point(size=3)+

    geom_abline(intercept=0,slope=1,color="black",linetype="dashed")+
    theme_bw()+
    coord_fixed(xlim=c(0,1),ylim=c(0,1))+
    #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+
    scale_shape_manual(values=c(21,22,24,25),
                       labels = c("High Strength, High Variance",
                                  "High Strength, Low Variance",
                                  "Low Strength, High Variance",
                                  "Low Strength, Low Variance"))+
    scale_y_continuous(name="p('old'|old)")+
    scale_x_continuous(name = expression(paste("p('old'|new)"^.)))+
    theme_bw()+
    scale_fill_manual(values = c("#0072B2","#E69F00","#CC79A7","#009E73"),
                       labels = c("High Strength, High Variance",
                                  "High Strength, Low Variance",
                                  "Low Strength, High Variance",
                                  "Low Strength, Low Variance"))+
    scale_color_manual(values = c("#0072B2","#E69F00","#CC79A7","#009E73"),

                       labels = c("High Strength, High Variance",
                                  "High Strength, Low Variance",
                                  "Low Strength, High Variance",
                                  "Low Strength, Low Variance"))+
    theme(axis.text = element_text(size=8),
          axis.title = element_text(size=10),
          #legend.position = c(0.7,0.2),
          legend.title = element_blank(),
          legend.background = element_rect(fill="transparent",color="transparent"),
          panel.grid = element_blank(),
          plot.background = element_rect(fill="white"),
          strip.background = element_rect(fill="black"),
          strip.text = element_text(color="white",face = "bold",size=13))

  list(roc)
}

library(cowplot)
roczroc <- plot_grid(makeROC(exp %>% filter(id == "simSB2022_e1_1"))[[1]],
          makeROC(exp %>% filter(id == "simSB2022_e1_1"))[[2]],rel_widths=c(0.5,0.5),
          align='hv')


exp <- exp %>% mutate(oldnew = isold) %>%
  mutate(oldnew = ifelse(oldnew==1,"Old","New")) %>%
  mutate(response=rating)

source("FitSimple.R")
for(modelname in c("GaussianUVSDT")){
  allfits2 <- NULL
  for(idn in unique(exp$id)){
    for(cond in c("A","B","C","D")){

    fit <- FitSDT(data = exp %>% filter(id == idn) %>% filter(condition == cond),
                  rep = 5,
                  model = modelname,
                  #idn = idn,
                  freqdat= T)

    best <- fit %>% filter(objective == min(objective)) %>% .[1,] %>%
      mutate(model = modelname,
             id = idn,
             condition = cond)
    #dp <- prep_data(SB2021_e3 %>% filter(id == idn),freqdat =F,numcategories = 6)

    # pars <- tibble::rownames_to_column(as.data.frame(t(best %>% select(any_of(c("d_ho","d_lo","d_hn","sig_ho","sig_lo","sig_hn",
    #                                                                             "beta_ho","beta_lo","beta_hn",
    #                                                                     "c1","dc1","dc2","dc3","dc4"))))), "parameter") %>%
    #   mutate(parameter = factor(parameter,levels=c("d_ho","d_lo","d_hn",
    #                                                "sig_ho","sig_lo","sig_hn",
    #                                                "beta_ho","beta_lo","beta_hn",
    #                                                "c1","dc1","dc2","dc3","dc4"))) %>%
    #   .$V1

    # expdata <- predict_frequencies(data_list=dp$datalist,modelname,pars)
    #
    # Gsquared <- GetGsquared(modelname,dp,expdata)

    fitted <- best %>%
      mutate(AIC = 2 * objective + 2 * npar)

    allfits2 <- bind_rows(allfits2,fitted)

    }
  }
  saveRDS(allfits2,file=paste0("SimulateVarianceItemEffects2_modelfits_",modelname,".rds"))

}

modelname<-"GaussianUVSDT"

allfits2 <- readRDS(paste0("Simulations/ItemVariance/SimulateVarianceItemEffects2_modelfits_",modelname,".rds"))

all <-allfits2


testsALL <- all %>% select(id,condition,muo,sigo) %>%
  pivot_longer(values_to="value",names_to="par",
               cols = c(muo,sigo)) %>%
  mutate(condition = factor(condition,levels=c("A","B","C","D"),
                            labels=c("HS\nHV",
                                     "HS\nLV",
                                     "LS\nHV",
                                     "LS\nLV"))) %>%
  mutate(par= factor(par,levels=c("muo","sigo"),
                     labels=c(expression(paste("U-V Gaussian: ", mu[o] - mu[n])),expression(paste("U-V Gaussian: ",sigma[o])))))



Exp12pars<-ggplot(testsALL,
       aes(x=condition,y=value,fill = condition)) +
  #geom_flat_violin(position = position_nudge(x = .2, y = 0), colour = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  #geom_jitter(width = 0.1, alpha = 0.05, shape = 16, size = 1) +
  #geom_point(data = forgraph_means, aes(x = Parameter, y = mean_estimate), size = 3, colour = "black", shape = 16) +
  scale_fill_manual(values = c("#0072B2","#E69F00","#CC79A7","#009E73")) +
  # scale_x_discrete(labels = c(expression(sigma[nh]),
  #                             expression(sigma[oh]),
  #                             expression(sigma[ol]),
  #                             expression(mu[oh]),
  #                             expression(mu[ol])))+

  scale_x_discrete(labels=c("H/H","H/L",
                            "L/H","L/L"))+

  facet_grid(.~par,labeller=label_parsed)+
  theme_bw() +
  theme(axis.text.x = element_text(size=10),
        legend.position = "none",
        strip.text = element_text(size=14),
        panel.grid = element_blank(),
        axis.title.x = element_blank())+
  scale_y_continuous(limits=c(0,1.5),breaks=c(0,1,2),name="parameter value")+

  #scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
  scale_size(guide = "none")



fitsSB <- list.files("EmpiricalFits/Fits2/","fit",full.names = T)
fitsSB <- fitsSB[grepl(paste(models,collapse="|"),fitsSB)]
allfitsSB <-  purrr::map(fitsSB, readRDS) %>%
  bind_rows() %>%
  group_by(exp,id,condition,model) %>%
  filter(objective == min(objective))



all <-allfitsSB %>% filter(model=="GaussianUVSDT")


testsEmp <- all %>% select(id,condition,muo,sigo) %>%
  pivot_longer(values_to="value",names_to="par",
               cols = c(muo,sigo)) %>%
  mutate(condition = factor(condition,levels=c("A","B","C","D"),
                            labels=c("HS\nHV",
                              "HS\nLV",
                              "LS\nHV",
                              "LS\nLV"))) %>%
  mutate(par= factor(par,levels=c("muo","sigo"),
                     labels=c(expression(paste("U-V Gaussian: ", mu[o] - mu[n])),expression(paste("U-V Gaussian: ",sigma[o])))))

tests <- bind_rows(testsALL %>% mutate(type="aSimulation"),testsEmp %>% mutate(type="bemp")) %>%
  mutate(type = factor(type,
                       labels=c("simulation","'S&B' (2022)")))

parests <- ggplot(tests,
       aes(x=condition,y=value,fill = condition)) +
  #geom_flat_violin(position = position_nudge(x = .2, y = 0), colour = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
 # geom_jitter(width = 0.1, alpha = 0.1, shape = 16, size = 1) +
  #geom_point(data = forgraph_means, aes(x = Parameter, y = mean_estimate), size = 3, colour = "black", shape = 16) +
  scale_fill_manual(values = c("#0072B2","#E69F00","#CC79A7","#000000")) +
  # scale_x_discrete(labels = c(expression(sigma[nh]),
  #                             expression(sigma[oh]),
  #                             expression(sigma[ol]),
  #                             expression(mu[oh]),
  #                             expression(mu[ol])))+
  facet_grid(par~type,scales="free",labeller=label_parsed)+
  scale_y_continuous(limits=c(0,4),name="estimate")+
  theme_bw() +
  theme(legend.position = "none",
        axis.text= element_text(size=8),
        strip.background = element_rect(color="transparent",fill="transparent"))+
  scale_x_discrete(name="condition")+
  #scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
  scale_size(guide = "none")


plot_grid(makeROC(exp %>% filter(id == "simSB2022_e1_1"))[[1]],Exp12pars,nrow=2,
          rel_heights=c(0.5,0.5),labels=c("A","B"))

ggsave("Figures/SimStimVar_v2.png", units="cm", width=15, height=15, dpi=600)

#
# muo_anova <- afex::aov_ez(id = "id",
#                       dv = "muo",
#                       data = allfits2,
#                       within = "condition",
#                       anova_table = list(es = "pes"))
#
# muo_margins <- emmeans::emmeans(muo_anova, ~ condition)
# # Pairwise comparisons
# muo_pairs <- pairs(sigma_margins, adjust = "bonf")
#
#
# sigo_anova <- afex::aov_ez(id = "id",
#                           dv = "sigo",
#                           data = allfits2,
#                           within = "condition",
#                           anova_table = list(es = "pes"))
#
# sigo_margins <- emmeans::emmeans(sigo_anova, ~ condition)
# # Pairwise comparisons
# sigo_pairs <- pairs(sigo_margins, adjust = "bonf")
#
#




