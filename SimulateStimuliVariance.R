# General informations -------------------

# Simulate effect of changing variance of item effects on ROCs
# for data in SB2021_e3, ie. within-subject conditions

library(tidyverse)

sim <- function(basestrength,sd_r2_high,sd_r2_low){

  mu = basestrength
 # r1 = rnorm(1,0,sd_r1)

  r2_old_high = rnorm(1000,mu,sd_r2_high)
  r2_old_low = rnorm(1000,mu,sd_r2_low)
  r2_new_high = rnorm(1000,0,sd_r2_high)
  r2_new_low = rnorm(1000,0,sd_r2_low)

  # r2_old_high = rnorm(1000,0,sd_r2_high)
  # r2_old_low = rnorm(1000,0,sd_r2_low)
  # r2_new_high = rnorm(1000,0,sd_r2_high)
  # r2_new_low = rnorm(1000,0,sd_r2_low)


  munew = rnorm(1000,0,1)

  m_low<-1

  crit<- c(m_low - 1.5 ,m_low-0.75,m_low, m_low+0.75, m_low+1.5)

  return(list(item_oh = munew + r2_old_high,
              item_ol = munew + r2_old_low,
              item_nh = munew + r2_new_high,
              item_nl = munew + r2_new_low,
              crit  = crit))
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

  data.frame(Freq=c(makerating(test$crit,test$item_oh),
                    makerating(test$crit,test$item_ol),
                    makerating(test$crit,test$item_nh),
                    makerating(test$crit,test$item_nl)),
             condition = rep(c("high_ev_1","low_ev_1","high_ev_0","low_ev_0"),each=6),
             isold = rep(c(1,0),each=12),
             rating = rep(c(1:6),4)
  )

}





exp <- NULL

for(i in c(1:200)){
simdata <- sim(1,2.5,1)
dfFreq <- makeFreq(simdata) %>% mutate(id = paste0("simSB2022_e3_",i))
exp <- bind_rows(exp,dfFreq)

}

exp <- exp %>% mutate(isold = ifelse(condition=="low_ev_0",0,1))


example_exp <- exp %>% filter(id == "simSB2022_e3_1")


makeROC <- function(exp){




      expOld <- exp %>% filter(isold == 1)
      expNew <- exp %>% filter(isold == 0)

      expInd <- list()
      for(i in seq_along(unique(expOld$condition))){

        expInd[[i]] <- expOld %>% filter(condition == unique(expOld$condition)[[i]]) %>%
          bind_rows(expNew) %>% mutate(condition = unique(expOld$condition)[[i]])
      }

      ex <- bind_rows(expInd) %>%
        ungroup() %>%
        select(id,condition,isold,rating,Freq) %>%
        mutate(condition = factor(condition,levels=c("low_ev_0","high_ev_0","low_ev_1","high_ev_1"),
                                  labels=c("New, low V","New, high V","Old, low V", "Old, high V")))


  makerocdata <- ex %>%
    group_by(id,condition,isold,rating) %>%
    mutate(cumresp = Freq) %>%
    mutate(rating = factor(rating, levels = c(max(ex$rating):1))) %>%
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
#
  # makezrocdata <- ex %>%
  #   group_by(id,condition,isold,rating) %>%
  #   mutate(cumresp = Freq) %>%
  #   mutate(rating = factor(rating, levels = c(max(ex$rating):1))) %>%
  #   arrange(id,condition,isold,rating) %>%
  #   mutate(isold = factor(isold,levels=c(0,1),labels=c("INew","IOld"))) %>%
  #   group_by(id,condition,isold) %>%
  #   complete(rating,fill = list(cumresp = 0)) %>%
  #   summarize(cumprop = qnorm(cumsum(cumresp/sum(cumresp))),
  #             rating = rating) %>%
  #   ungroup() %>%
  #   pivot_wider(names_from="isold",values_from="cumprop",
  #               id_cols = c("id","condition","rating")) %>%
  #   filter(rating != 1)

  makerocdata


}

makerocdata <- makeROC(example_exp)

Exp3simrocall <- ggplot(makerocdata %>%  mutate(id = ". vs New-Low-V items"),aes(x = INew,y=IOld,group=condition,shape=condition,
                              text=condition,color=condition))+
  facet_wrap(id~.,ncol=6)+

  geom_point(size=3)+
  geom_line(size=1)+
  geom_abline(intercept=0,slope=1,color="black",linetype="dashed")+
  theme_bw()+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+
  scale_shape_manual(values=c(16,17,18),
                     guide = guide_legend(reverse = TRUE))+
  scale_x_continuous(name="p('old'|new, low V)")+
  scale_y_continuous(name = "p('old'|condition)")+
  theme_bw()+
  scale_color_manual(values = c("#0072B2","#E69F00","#CC79A7"),
                     guide = guide_legend(reverse = TRUE))+

  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        legend.position = c(0.7,0.2),
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank())

makerocdata_new <- makeROC(example_exp %>% filter(condition %in% c("high_ev_0","low_ev_0")) %>%
                             mutate(isold = ifelse(condition == "low_ev_0",0,1)))

rocnewonly <- ggplot(makerocdata_new[[1]] %>%  mutate(id = "New-High-V vs New-Low-V items"),aes(x = INew,y=IOld,group=condition,shape=condition,
                                                                                   text=condition,color=condition))+
  facet_wrap(id~.,ncol=6)+

  geom_point(size=3)+
  geom_line(size=1)+
  geom_abline(intercept=0,slope=1,color="black",linetype="dashed")+
  theme_bw()+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+
  scale_shape_manual(values=c(16,17,18))+
  scale_x_continuous(name="p('old'|new, low V)")+
  scale_y_continuous(name = "p('old'|new, high V)")+
  theme_bw()+
  scale_color_manual(values = c("#0072B2","#E69F00","#CC79A7"))+

  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        legend.position = "none",
        strip.background = element_rect(fill="black"),
        strip.text = element_text(color="white",face = "bold",size=13))


makerocdata_old <- makeROC(example_exp %>% filter(condition %in% c("high_ev_1","low_ev_1")) %>%
                             mutate(isold = ifelse(condition == "low_ev_1",0,1)))

modelname <- "baselineUVSDT"

bestfitsSB3 <- readRDS(paste0("Simulations/ItemVariance/SimulateVariance_modelfits_",modelname,"_SB2020e3.rds"))

bestfits <- bestfitsSB3 %>% filter(model == "baselineUVSDT")

simROC <- NULL
simEST <- NULL
for(i in unique(bestfits$id)[1]){

  best <- bestfits %>% filter(id == i)


  rocdatUV <-  tibble(crit = seq(-8,8,.01)) %>%
    mutate(hit = pnorm(crit,mean=best$d_ho,sd=best$sig_ho,lower.tail=F),
           fa = pnorm(crit,mean=best$d_lo,sd=best$sig_lo,lower.tail = F)) %>%
    mutate(id = best$id)

  two <- tibble(sig_ho = best$sig_ho,
                sig_lo = best$sig_lo) %>%
    mutate(ratio = sig_ho/sig_lo) %>%
    mutate(ratiodesc = ifelse(ratio>1,1,0)) %>%
    mutate(id = best$id)

  simROC <- simROC %>% bind_rows(rocdatUV)
  simEST <- simEST %>% bind_rows(two)

}



Exp3rocoldonly <- ggplot(makerocdata_old %>%  mutate(id = "Old-High-V vs Old-Low-V items"),aes(x = INew,y=IOld,group=condition,shape=condition,
                                                                                                text=condition,fill=condition,color=condition))+
  facet_wrap(id~.,ncol=6)+
  geom_text(data = simEST %>% mutate(condition = "Old, high V",
                                     id = "Old-High-V vs Old-Low-V items"),
            aes(x=0.2,y=0.9,
                label=paste("frac(sigma[oh],sigma[ol])","==",round(ratio,2))),parse=T,
            color="black",size=5)+
  geom_path(data = simROC %>% mutate(condition="Old, high V",
                                     id = "Old-High-V vs Old-Low-V items"),
            aes(x=fa,y=hit,group=condition),color="black",size=1)+
  geom_line(size=1)+
  geom_point(size=5,color="black")+

  geom_abline(intercept=0,slope=1,color="black",linetype="dashed")+
  theme_bw()+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+
  scale_shape_manual(values=c(21))+
  scale_x_continuous(name="p('old'|old, low V)")+
  scale_y_continuous(name = "p('old'|old, high V)")+

  scale_fill_manual(values = c("#CC79A7"))+
  scale_color_manual(values = c("#CC79A7"))+
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        legend.position = "none")

# zroc <- ggplot(makezrocdata,aes(x = INew,y=IOld,group=condition,
#                                 text=condition,color=condition,shape=condition))+
#   facet_wrap(id~.,ncol=6)+
#
#   geom_point(size=3)+
#   geom_line(size=1)+
#   geom_abline(intercept=0,slope=1,color="black",linetype="dashed")+
#   theme_bw()+
#   coord_fixed(xlim=c(-2.5,2.5),ylim=c(-2.5,2.5))+
#   #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+
#   scale_shape_manual(values=c(16,17,18))+
#   scale_x_continuous(name=expression(paste(phi^-1,"(p('old'|new))")))+
#   scale_y_continuous(name=expression(paste(phi^-1,"(p('old'|old))")))+
#   theme_bw()+
#   scale_color_manual(values = c("#0072B2","#E69F00","#CC79A7"))+
#   theme(axis.text = element_text(size=13),
#         axis.title = element_text(size=13),
#         strip.background = element_rect(fill="black"),
#         strip.text = element_text(color="white",face = "bold",size=12))



modelname <-"baselineUVSDT"

allfits2 <- readRDS(paste0("SimulateVariance_modelfits_",modelname,"_SB2020e3.rds"))
#
# all <- allfits2 %>% mutate(d_hn = 0) %>%
#   mutate(sig_hn = log(sig_hn),
#          sig_ho = log(sig_ho),
#          sig_lo = log(sig_lo))
# meand <- all %>%
#   ungroup() %>%
#   summarise(across(c(d_lo,d_ho,d_hn,sig_ho,sig_lo,sig_hn), mean))
# sdd <-  all %>%
#   ungroup() %>%
#   summarise(across(c(d_lo,d_ho,d_hn,sig_ho,sig_lo,sig_hn), sd))
#
# to_excl <- all %>% filter(d_lo > (meand$d_lo + 3 * sdd$d_lo) |
#
#                             d_ho > (meand$d_ho + 3 * sdd$d_ho) |
#
#                             d_hn >  (meand$d_hn + 3 * sdd$d_hn)|
#
#                             sig_ho > (meand$sig_ho + 3 * sdd$sig_ho) |
#                             sig_lo >(meand$sig_lo + 3 * sdd$sig_lo) |
#                             sig_hn >(meand$sig_hn + 3 * sdd$sig_hn))
#
# to_excl2 <- all %>% filter(d_lo < (meand$d_lo - 3 * sdd$d_lo) |
#
#                              d_ho < (meand$d_ho - 3 * sdd$d_ho) |
#
#                              d_hn <  (meand$d_hn - 3 * sdd$d_hn))
#



#all <- all %>% filter(!(id %in% to_excl$id)) %>% filter(!(id %in% to_excl2$id))

testsALL <- allfits2 %>% select(id,d_lo,d_ho,d_hn,sig_ho,sig_lo,sig_hn) %>%
  pivot_longer(values_to="value",names_to="par",
               cols = c(d_hn,d_lo,d_ho,sig_ho,sig_lo,sig_hn)) %>%
  mutate(condition = ifelse(grepl("h",par),"high","low")) %>%
  mutate(par= factor(par,levels=c("sig_hn","sig_ho","sig_lo","d_hn","d_ho","d_lo"))
  )




ttestBF(testsALL %>% filter(par == "d_ho") %>% .$value,
        mu=0)

parest <- ggplot(testsALL,
       aes(x=par,y=value,fill = condition)) +
  #geom_flat_violin(position = position_nudge(x = .2, y = 0), colour = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.1, shape = 16, size = 1) +
  #geom_point(data = forgraph_means, aes(x = Parameter, y = mean_estimate), size = 3, colour = "black", shape = 16) +
  scale_fill_manual(values = c("#0072B2","#E69F00","#CC79A7","#000000")) +
  # scale_x_discrete(labels = c(expression(sigma[nh]),
  #                             expression(sigma[oh]),
  #                             expression(sigma[ol]),
  #                             expression(mu[oh]),
  #                             expression(mu[ol])))+
  scale_x_discrete(labels=c(expression(paste(sigma[nh])),
                            expression(paste(sigma[oh])),
                            expression(paste(sigma[ol])),
                            expression(paste(mu[nh])),
                            expression(paste(mu[oh])),
                            expression(paste(mu[ol]))))+

  #facet_wrap(.~partype,labeller=label_parsed,scales="free_x")+
  theme_bw() +
  scale_y_continuous(limits=c(0,2.5),breaks=c(0,1,2,3),name="estimate")+
  theme(axis.text.x = element_text(size = 15),
        legend.position = c(0.8,0.8),
        axis.title.x = element_blank()) +
  #scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
  scale_size(guide = "none")




# Exp 3 data:

exp3 <- readRDS("EmpiricalData/data_SB2022_e3.rds")
exp3f <- exp3 %>% mutate(isold = oldnew,condition) %>%
  group_by(exp,id,condition,isold,rating) %>% summarize(Freq=length(rating))


# Exp 3 fit


allfits <- readRDS("SB2021_e3_modelfits_baselineUVSDT.rds")

all <- allfits %>%
  group_by(id,model) %>%  filter(objective == min(objective)) %>%
  distinct() %>%
  mutate(sig_ho = sig_ho,
         sig_hn = sig_hn,
         sig_lo = sig_lo) %>%
  arrange(id)


simROC <- NULL
simEST <- NULL
for(i in unique(all$id)){

  best <- all %>% filter(id == i)


rocdatUV <-  tibble(crit = seq(-8,8,.01)) %>%
  mutate(hit = pnorm(crit,mean=best$d_ho,sd=best$sig_ho,lower.tail=F),
         fa = pnorm(crit,mean=best$d_lo,sd=best$sig_lo,lower.tail = F)) %>%
  mutate(id = best$id)

two <- tibble(sig_ho = best$sig_ho,
                 sig_lo = best$sig_lo) %>%
  mutate(ratio = sig_ho/sig_lo) %>%
  mutate(ratiodesc = ifelse(ratio>1,1,0)) %>%
  mutate(id = best$id)

simROC <- simROC %>% bind_rows(rocdatUV)
simEST <- simEST %>% bind_rows(two)

}

rocalldat <- NULL
for(idn in unique(exp3f$id)){

  ones <- makeROC( exp3f %>% filter(id == idn))
  rocalldat <- rocalldat %>% bind_rows(ones)
}
ggplot(rocalldat,aes(x = INew,y=IOld,group=interaction(id,condition),shape=condition,
                                                                text=condition,color=condition))+
  #facet_wrap(id~.,ncol=6)+

  #geom_point(size=3,alpha=0.2)+
  geom_line(size=1,alpha=0.2)+
  geom_abline(intercept=0,slope=1,color="black",linetype="dashed")+
  theme_bw()+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+
  scale_shape_manual(values=c(16,17,18))+
  scale_x_continuous(name="p('old'|new, low V)")+
  scale_y_continuous(name = "p('old'|condition)")+
  theme_bw()+
  scale_color_manual(values = c("#0072B2","#E69F00","#CC79A7"))+

  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        #legend.position = "none",
        strip.background = element_rect(fill="black"),
        strip.text = element_text(color="white",face = "bold",size=13))



makerocdata_old <- makeROC(example_exp %>% filter(condition %in% c("high_ev_1","low_ev_1")) %>%
                             mutate(isold = ifelse(condition == "low_ev_1",0,1)))



rocolddat <- NULL
for(idn in unique(exp3f$id)){

  ones <- makeROC( exp3f %>% filter(id == idn) %>% filter(condition %in% c("high_ev_1","low_ev_1")) %>%
                     mutate(isold = ifelse(condition == "low_ev_1",0,1)))
  rocolddat <- rocolddat %>% bind_rows(ones)
}



rocoldonly <- ggplot(rocolddat,aes(x = INew,y=IOld,group=interaction(id,condition),shape=condition,
                                   text=condition,color=condition))+
  facet_wrap(id~.,ncol=15)+



  geom_path(data = simROC %>% mutate(condition = "Old, high V"),aes(x=fa,y=hit,group=id),size=1,color="black")+
  geom_text(data = simEST %>% mutate(condition = "Old, high V"),
            aes(x=0.3,y=0.8,
                label=paste("frac(sigma[oh],sigma[ol])","==",round(ratio,2))),parse=T,
            color="black",size=4)+

 # paste("'Mean petal' ==", round(Petal.Width, digits=2), "* Omega"))


 # geom_line(size=1)+
  geom_point(size=4,color="black",fill="#CC79A7")+
  geom_abline(intercept=0,slope=1,color="black",linetype="dashed")+
  theme_bw()+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+
  scale_shape_manual(values=c(21))+
  scale_x_continuous(name="p('old'|old, low V)")+
  scale_y_continuous(name = "p('old'|old, high V)")+
  theme_bw()+
  scale_color_manual(values = c("#CC79A7"))+

  theme(axis.text = element_text(size=8),
        axis.title = element_text(size=8),
        legend.position = "none",
        strip.background = element_rect(fill="transparent"),
        strip.text = element_text(color="black",size=8))


uniqueid <- sample(unique(rocolddat %>% .$id),size=4,replace=F)

rocoldonly <- ggplot(rocolddat %>% filter(id %in% uniqueid),aes(x = INew,y=IOld,group=interaction(id,condition),shape=condition,
                                   text=condition,color=condition))+
  facet_wrap(id~.,ncol=4)+



  geom_path(data = simROC  %>% filter(id %in% uniqueid)%>% mutate(condition = "Old, high V"),aes(x=fa,y=hit,group=id),size=1,color="black")+
  geom_text(data = simEST %>% filter(id %in% uniqueid) %>% mutate(condition = "Old, high V"),
            aes(x=0.3,y=0.8,
                label=paste("frac(sigma[oh],sigma[ol])","==",round(ratio,2))),parse=T,
            color="black",size=4)+

  # paste("'Mean petal' ==", round(Petal.Width, digits=2), "* Omega"))


  # geom_line(size=1)+
  geom_point(size=4,color="black",fill="#CC79A7")+
  geom_abline(intercept=0,slope=1,color="black",linetype="dashed")+
  theme_bw()+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+
  scale_shape_manual(values=c(21))+
  scale_x_continuous(name="p('old'|old, low V)")+
  scale_y_continuous(name = "p('old'|old, high V)")+
  theme_bw()+
  scale_color_manual(values = c("#CC79A7"))+

  theme(axis.text = element_text(size=8),
        axis.title = element_text(size=8),
        legend.position = "none",
        strip.background = element_rect(fill="transparent"),
        strip.text = element_text(color="black",size=8))


source("FitMultiple.R")


  allfits2 <- NULL
for(idn in unique(exp$id)){
  for(modelname in c("baselineUVSDT")){


    fit <- FitSDT(data = exp %>% filter(id == idn),
       rep = 5,
       model = modelname,
       idn = idn,
       freqdat= T)

best <- fit %>% filter(objective == min(objective)) %>% .[1,] %>%
  mutate(model = modelname,
         id = idn)

fitted <- best %>%
  mutate(AIC = 2 * objective + 2 * npar)

allfits2 <- bind_rows(allfits2,fitted)

  }
}

#saveRDS(allfits2,file=paste0("SimulateVariance_modelfits_",modelname,"_SB2020e3.rds"))







allfits2<-readRDS("SimulateVariance_modelfits_baselineUVSDT_SB2020e3.rds")

library(afex)
library(emmeans)
library(BayesFactor)
# Run ANOVA
sigma_anova <- aov_ez(id = "id",
                      dv = "value",
                      data = tests,
                      anova_table = list(es = "pes"))

# Get BF
lmBF(value ~ par, data = tests, whichRandom = "id")

# Run post hoc tests
# Marginal means
sigma_margins <- emmeans(sigma_anova, ~ par)
# Pairwise comparisons
sigma_pairs <- pairs(sigma_margins, adjust = "bonf")
sigma_pairs

all <- allfits2 %>% filter(model == "baselineUVSDT")
tests2 <- all %>%

  select(id,d_lo,d_ho,d_hn,sig_ho,sig_lo,sig_hn) %>%
  pivot_longer(values_to="value",names_to="par",
               cols = c(d_lo,d_ho,d_hn,sig_ho,sig_lo,sig_hn)) %>%
  mutate(type = ifelse(par %in% c("d_lo","d_ho","d_hn"),"mu","sigma")) %>%
  mutate(type = factor(type)) %>%
  mutate(par = case_when(par == "d_lo" ~ "Old\nLow Var",
                         par == "d_ho" ~ "Old\nHigh Var",
                         par == "d_hn" ~ "New\nHigh Var",
                         par == "sig_lo" ~ "Old\nLow Var",
                         par == "sig_ho" ~ "Old\nHigh Var",
                         par == "sig_hn" ~ "New\nHigh Var")) %>%
  mutate(par = factor(par, levels = c( "Old\nHigh Var","Old\nLow Var","New\nHigh Var"))) %>%
  mutate(type= factor(type,levels=c("mu","sigma"),
                     labels=c(expression(paste("U-V Gaussian: ", mu)),expression(paste("U-V Gaussian: ",sigma)))))

t.test(all$d_ho - all$d_hn,all$d_lo,paired=T)

Exp3par <- ggplot(tests2,
       aes(x=par,y=value,fill = par)) +
  #geom_flat_violin(position = position_nudge(x = .2, y = 0), colour = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.05, shape = 16, size = 1) +
  #geom_point(data = forgraph_means, aes(x = Parameter, y = mean_estimate), size = 3, colour = "black", shape = 16) +
  #scale_fill_manual(values = c("#0072B2","#E69F00","#CC79A7","#009E73")) +
  # scale_x_discrete(labels = c(expression(sigma[nh]),
  #                             expression(sigma[oh]),
  #                             expression(sigma[ol]),
  #                             expression(mu[oh]),
  #                             expression(mu[ol])))+



  facet_grid(.~type,labeller=label_parsed, scales = "free_x")+
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.title.x = element_blank())+
  scale_y_continuous(limits=c(-0.5,2.5),breaks=c(0,1,2),name="parameter value")+

  #scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
  scale_size(guide = "none")

# Make plot for paper

library(cowplot)
ROCS <- plot_grid(Exp3simrocall,Exp3rocoldonly,rel_widths=c(0.55,0.45),nrow=1,labels=c("A","B"))
Exp3parp <- plot_grid(Exp3par,labels=c("C"))
comb1 <- plot_grid(ROCS,Exp3parp,ncol=1,rel_heights = c(0.6,0.4))
Exp3oldonly <- plot_grid(rocoldonly,labels="D")
comb2 <- plot_grid(comb1,Exp3oldonly,ncol=1,rel_heights=c(2,1))


mu_anova <- aov_ez(id = "id",
                   dv = "value",
                   data = tests2,
                   within = "par",
                   anova_table = list(es = "pes"))

# Run post hoc tests
# Marginal means
mu_margins <- emmeans(mu_anova, ~ par)
# Pairwise comparisons
mu_pairs <- pairs(mu_margins, adjust = "bonf")
mu_pairs



testsALL <- all %>% select(id,d_lo,d_ho,sig_ho,sig_lo,sig_hn) %>%
  pivot_longer(values_to="value",names_to="par",
               cols = c(d_lo,d_ho,sig_ho,sig_lo,sig_hn)) %>%
  mutate(par = factor(par,levels=c("sig_hn","sig_ho","sig_lo","d_ho","d_lo"))) %>%
  mutate(value = ifelse(grepl("sig",par),exp(value),value)) %>%
  mutate(variability = ifelse(par %in% c("sig_hn","sig_ho","d_ho"),"High","Low"))





ggplot(testsALL,
       aes(x=par,y=value,fill = variability)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), colour = NA) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.1, shape = 16, size = 1) +
  #geom_point(data = forgraph_means, aes(x = Parameter, y = mean_estimate), size = 3, colour = "black", shape = 16) +
  scale_fill_manual(values = c("#00798c", "#edae49")) +
  scale_x_discrete(labels = c(expression(sigma[nh]),
                              expression(sigma[oh]),
                              expression(sigma[ol]),
                              expression(mu[oh]),
                              expression(mu[ol])))+
  theme_classic() +
  theme(axis.text.x = element_text(size = 15)) +
  #scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)) +
  scale_size(guide = "none")


# UV at different levels of ------------------------------



simROC <- NULL

  rocdatUV <-  tibble(crit = seq(-8,8,.01)) %>%
      mutate(hit = pnorm(crit,mean=1,sd=2,lower.tail=F),
           fa = pnorm(crit,mean=2,sd=1,lower.tail = F)) %>%
    mutate(type="m3")

  simROC <- simROC %>% bind_rows(rocdatUV)

  parse.labels <- function(x) parse(text = x)

  simROC <- simROC %>% mutate(type = factor(type,levels=c("m1","m2","m3")))#,
                                            #labels=expression(paste("mu[oh] > mu[ol]","mu[oh] = mu[ol]","mu[oh] < mu[ol]"))))

ggplot(simROC,aes(x = fa,y=hit,color=type,group=type))+
  #facet_wrap(id~.,ncol=6)+
  geom_path(aes(color=type))+

  geom_abline(intercept=0,slope=1,color="black",linetype="dashed")+
  theme_bw()+
  coord_fixed(xlim=c(0,1),ylim=c(0,1))+
  #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+
  scale_shape_manual(values=c(21))+
  scale_x_continuous(name="p('old'|old, low V)")+
  scale_y_continuous(name = "p('old'|old, high V)")+
  theme_bw()+
  scale_color_manual(values = c("#CC79A7","#000000","#0072B2"))+
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        #legend.position = "none",
        strip.background = element_rect(fill="black"),
        strip.text = element_text(color="white",face = "bold",size=13))



