# General information ---------------------------------------
# Illustrate U-V Gaussian, Gumbel and ExGaussian models

source("FitSimple.R")
library(cowplot)
models <- c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT")

# density graphs, based on median par estimates from fits to Spanton & Berry

generatedensities <- function(model,par){

  nitems <- 1e6
if(model == "GaussianUVSDT"){
      strength_old <- rnorm(nitems,mean = 3,sd = 1.5)
      strength_new <- rnorm(nitems,mean = 0, sd = 1)
  #     crits <- c(-Inf,cumsum(as.numeric(par[c(3:7)])),Inf)
      strength_old2 <- NA
   } else if(model == "GumbelEVSDT"){

       # switch the sign on muo to get equivalent results to the pgumbel fitting
       strength_old <- ordinal::rgumbel(nitems,location = 3,scale = 1,max=F)
       strength_new <- ordinal::rgumbel(nitems,location = 0,scale = 1,max=F)
  #    crits <- c(-Inf,cumsum(as.numeric(par[c(2:6)])),Inf)
       strength_old2 <- NA
   } else if(model == "ExGaussNormEVSDT"){

     strength_old <- brms::rexgaussian(nitems,mu = 3,sigma = 1,beta=1.5)
     strength_new <- rnorm(nitems,mean = 0, sd = 1)
     strength_old2 <- NA
   } else if(model=="GaussianEVSDT"){

       strength_old <- rnorm(nitems,mean = 3,sd = 1)
       strength_new <- rnorm(nitems,mean = 0, sd = 1)
       strength_old2 <- NA
   } else if(model=="Mix0"){

     strength_old <- rnorm(nitems,mean = 3,sd = 1)
     strength_new <- rnorm(nitems,mean = 0, sd = 1)
     strength_old2 <- rnorm(nitems,mean = 0.2,sd = 1)
   } else if(model=="MixA"){

     strength_old <- rnorm(nitems,mean = 3,sd = 1)
     strength_new <- rnorm(nitems,mean = 0, sd = 1)
     strength_old2 <- rnorm(nitems,mean = 1,sd = 1)
   }

  return(list(old = strength_old,new=strength_new,old2 = strength_old2))
}

densityplot <- function(modelname,par){

genden <- generatedensities(modelname)


dats <- tibble(type = rep(c("old","new"),each=1e6),
               value = c(genden$old,genden$new)) %>%
  mutate(model = modelname)

}

# Make ROC and zROC
# values the same as the density graphs

rocdatExG1 <- tibble(crit = seq(-8,15,.01)) %>%
  mutate(hit = brms::pexgaussian(crit,mu=2+1.5,sigma=1,beta=1.5,lower.tail=F),
         fa = pnorm(crit,lower.tail = F)) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa)) %>%
  mutate(model="ExG") %>%
  mutate(type="a")%>%
  mutate(type2 = "muo")

rocdatExG2 <- tibble(crit = seq(-8,15,.01)) %>%
  mutate(hit = brms::pexgaussian(crit,mu=1 + 1.5,sigma=1,beta=1.5,lower.tail=F),
         fa = pnorm(crit,lower.tail = F)) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa)) %>%
  mutate(model="ExG") %>%
  mutate(type="b")%>%
  mutate(type2 = "muo")

rocdatExG3 <- tibble(crit = seq(-8,15,.01)) %>%
  mutate(hit = brms::pexgaussian(crit,mu=0+1.5,sigma=1,beta=1.5,lower.tail=F),
         fa = pnorm(crit,lower.tail = F)) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa)) %>%
  mutate(model="ExG") %>%
  mutate(type="c")%>%
  mutate(type2 = "muo")



rocdatExG4 <- tibble(crit = seq(-8,15,.01)) %>%
  mutate(hit = brms::pexgaussian(crit,mu=0+1.5,sigma=1,beta=1.5,lower.tail=F),
         fa = pnorm(crit,lower.tail = F)) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa)) %>%
  mutate(model="ExG") %>%
  mutate(type="a")%>%
  mutate(type2 = "sigma")

rocdatExG5 <- tibble(crit = seq(-8,15,.01)) %>%
  mutate(hit = brms::pexgaussian(crit,mu=0 + 2.5,sigma=1,beta=2.5,lower.tail=F),
         fa = pnorm(crit,lower.tail = F)) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa)) %>%
  mutate(model="ExG") %>%
  mutate(type="b")%>%
  mutate(type2 = "sigma")

rocdatExG6 <- tibble(crit = seq(-8,15,.01)) %>%
  mutate(hit = brms::pexgaussian(crit,mu=0+3.5,sigma=1,beta=3.5,lower.tail=F),
         fa = pnorm(crit,lower.tail = F)) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa)) %>%
  mutate(model="ExG") %>%
  mutate(type="c")%>%
  mutate(type2 = "sigma")

rocdatExG <- bind_rows(rocdatExG1,rocdatExG2,rocdatExG3,
                       rocdatExG4,rocdatExG5,rocdatExG6)



rocdatUV1 <-  tibble(crit = seq(-8,8,.01)) %>%
  mutate(hit = pnorm(crit,mean=1.5,sd=1.5,lower.tail=F),
         fa = pnorm(crit,lower.tail = F)) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa)) %>%
  mutate(model="UV") %>%
  mutate(type="a")%>%
  mutate(type2 = "sigma")


rocdatUV2 <-  tibble(crit = seq(-8,8,.01)) %>%
  mutate(hit = pnorm(crit,mean=3,sd=2.5,lower.tail=F),
         fa = pnorm(crit,lower.tail = F)) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa)) %>%
  mutate(model="UV") %>%
  mutate(type="b")%>%
  mutate(type2 = "sigma")


rocdatUV3 <-  tibble(crit = seq(-8,8,.01)) %>%
  mutate(hit = pnorm(crit,mean=3,sd=3.5,lower.tail=F),
         fa = pnorm(crit,lower.tail = F)) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa)) %>%
  mutate(model="UV") %>%
  mutate(type="c")%>%
  mutate(type2 = "sigma")


rocdatUV4 <-  tibble(crit = seq(-8,8,.01)) %>%
  mutate(hit = pnorm(crit,mean=3,sd=1.5,lower.tail=F),
         fa = pnorm(crit,lower.tail = F)) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa)) %>%
  mutate(model="UV") %>%
  mutate(type="a")%>%
  mutate(type2 = "muo")

rocdatUV5<-  tibble(crit = seq(-8,8,.01)) %>%
  mutate(hit = pnorm(crit,mean=2,sd=1.5,lower.tail=F),
         fa = pnorm(crit,lower.tail = F)) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa)) %>%
  mutate(model="UV") %>%
  mutate(type="b")%>%
  mutate(type2 = "muo")


rocdatUV6 <-  tibble(crit = seq(-8,8,.01)) %>%
  mutate(hit = pnorm(crit,mean=1,sd=1.5,lower.tail=F),
         fa = pnorm(crit,lower.tail = F)) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa)) %>%
  mutate(model="UV") %>%
  mutate(type="c") %>%
  mutate(type2 = "muo")

rocdatUV <- bind_rows(rocdatUV1,rocdatUV2,rocdatUV3,
                      rocdatUV4,rocdatUV5,rocdatUV6)

rocdatGumbel1 <- tibble(crit = seq(-8,15,.01)) %>%
  mutate(hit = exp(-exp(crit - 3)),
         fa =  exp(-exp(crit))) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa))%>%
  mutate(model="Gumbel") %>%
  mutate(type="a") %>%
  mutate(type2 = "muo")

rocdatGumbel2 <- tibble(crit = seq(-8,15,.01)) %>%
  mutate(hit = exp(-exp(crit - 2)),
         fa =  exp(-exp(crit))) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa))%>%
  mutate(model="Gumbel") %>%
  mutate(type="b") %>%
  mutate(type2 = "muo")


rocdatGumbel3 <- tibble(crit = seq(-8,15,.01)) %>%
  mutate(hit = exp(-exp(crit - 1)),
         fa =  exp(-exp(crit))) %>%
  mutate(threshcolor=NA,
         shapecolor=NA) %>%
  mutate(crit = as.character(round(crit,2))) %>%
  mutate(zhit = qnorm(hit),
         zfa = qnorm(fa))%>%
  mutate(model="Gumbel") %>%
  mutate(type="c") %>%
  mutate(type2 = "muo")


rocdatGumbel <- bind_rows(rocdatGumbel1,rocdatGumbel2,rocdatGumbel3)



rocdat <- bind_rows(rocdatGumbel,rocdatExG,rocdatUV)

makezROC <- function(rocdat,type,model){
if(type=="muo") {

  namelegend<-expression(paste(mu[o]))
  namelabels <- c(3,2,1)
} else {

  namelegend <- ifelse(model=="ExG",expression(paste(beta[o])),expression(paste(sigma[o])))
  namelabels <- c(1.5,2.5,3.5)
}



  rocdatEV <-  tibble(crit = seq(-8,8,.01)) %>%
    mutate(hit = pnorm(crit,mean=1.5,sd=1,lower.tail=F),
           fa = pnorm(crit,lower.tail = F)) %>%
    mutate(threshcolor=NA,
           shapecolor=NA) %>%
    #mutate(crit = as.character(round(crit,2))) %>%
    mutate(zhit = qnorm(hit),
           zfa = qnorm(fa)) %>%
    mutate(model="UV") %>%
    mutate(type="a")%>%
    mutate(type2 = "sigma")


ggplot(rocdatEV ,
       aes(x = fa,y = hit)) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  #coord_cartesian(xlim=c(-20,20),ylim=c(-20,20))+
  scale_y_continuous(name=expression(paste(phi^-1,"[p('old'|old)]")))+
  scale_x_continuous(name=expression(paste(phi^-1,"[p('old'|new)]")))+
 # geom_abline(intercept=parEVSDT[1],slope=1,color="black",size=1)+
  geom_abline(intercept=0,slope=1,color="grey",linetype="dotted")+
  geom_path()+

  theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title=element_blank(),
        #legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))


ggplot(rocdatEV ,
       aes(x = zfa,y = zhit)) +
  coord_cartesian(xlim=c(-2.5,2.5),ylim=c(-2.5,2.5))+
  #coord_cartesian(xlim=c(-20,20),ylim=c(-20,20))+
  scale_y_continuous(name=expression(paste(phi^-1,"[p('old'|old)]")))+
  scale_x_continuous(name=expression(paste(phi^-1,"[p('old'|new)]")))+
  # geom_abline(intercept=parEVSDT[1],slope=1,color="black",size=1)+
  geom_abline(intercept=0,slope=1,color="black")+
  geom_path(size=1)+

  theme( aspect.ratio=1,
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    #legend.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = "transparent"),
    strip.background = element_rect(fill = "white"),
    strip.placement = "outside",
    strip.text = element_text(size=12))


ggplot(rocdatEV ,
       aes(x = fa,y = hit)) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  #coord_cartesian(xlim=c(-20,20),ylim=c(-20,20))+
  scale_y_continuous(name=expression(paste(phi^-1,"[p('old'|old)]")))+
  scale_x_continuous(name=expression(paste(phi^-1,"[p('old'|new)]")))+
   geom_abline(intercept=0,slope=1,color="black")+
  #geom_abline(intercept=0,slope=1,color="grey",linetype="dotted")+
  geom_path(color="red",linetype="dashed",size=1)+
  geom_path(data = rocdatUV1,aes(x = fa,y = hit),size=1)+
  theme(
    aspect.ratio=1,
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank(),
    #legend.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = "transparent"),
    strip.background = element_rect(fill = "white"),
    strip.placement = "outside",
    strip.text = element_text(size=12))

}
makeROC <- function(rocdat,type,model) {

  if(type=="muo") {

    namelegend<-ifelse(model=="ExG",expression(paste(mu[o[G]])),expression(paste(mu[o])))
    if(model == "Gumbel"){
      namelabels <- c(3,2,1)
      } else if (model == "UV"){
        namelabels <- c(3,2,1)
        } else {
          namelabels <- c(2,1,0)
        }

    ann <- ifelse(model=="ExG",expression(paste(beta[o]," = 1.5")),
                  ifelse(model=="UV",expression(paste(sigma[o]," = 1.5")),
                         ""))


  } else {

    namelegend <- ifelse(model=="ExG",expression(paste(beta[o])),expression(paste(sigma[o])))
    namelabels <- c(1.5,2.5,3.5)
    ann <- ifelse(model=="ExG",expression(paste(mu[o[G]]," = 0")),expression(paste(mu[o]," = 3")))
  }

  ggplot(rocdat ,
       aes(x = fa,y = hit,linetype=type)) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  scale_y_continuous(name="p('old'|old)")+
  scale_x_continuous(name="p('old'|new)")+
  scale_linetype_manual(values=c("solid","dashed","dotted"),name=namelegend,
                          labels=namelabels)+
  # geom_abline(intercept=parEVSDT[1],slope=1,color="black",size=1)+
  geom_abline(intercept=0,slope=1,color="grey",linetype="dotted")+
  geom_path()+
  annotate("text",label=ann,x=0.3,y=0.1,size=3.5)+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_text(size=10),
    legend.title=element_text(size=10,hjust=1),
    legend.text=element_text(size=10),
    panel.background = element_rect(fill = "white"),
    legend.position = c(0.75,0.35),
    legend.background = element_rect(fill="transparent"),
    legend.key = element_rect(fill="transparent"),
    legend.key.size = unit(0.4, 'cm'),
    panel.border = element_rect(colour = "black", fill = "transparent"),
    strip.background = element_rect(fill = "white"),
    strip.placement = "outside",
    strip.text = element_text(size=12))

}

UVROC <- plot_grid(makeROC(rocdat %>% filter(type2=="muo") %>% filter(model=="UV"),"muo","UV"),
                   makeROC(rocdat %>% filter(type2=="sigma") %>% filter(model=="UV"),"sigma","UV"),nrow=1)

GumbelROC <- plot_grid(makeROC(rocdat %>% filter(model=="Gumbel"),"muo","Gumbel"),NULL,
                   nrow=1)

ExGROC <- plot_grid(makeROC(rocdat %>% filter(type2=="muo") %>% filter(model=="ExG"),"muo","ExG"),
                   makeROC(rocdat %>% filter(type2=="sigma") %>% filter(model=="ExG"),"sigma","ExG"),
                    nrow=1)


#c("#E69F00", "#009E73", "#56B4E9")



Ill_UVSDT <- ggplot(dats %>% filter(model=="U-V Gaussian"),aes(x = value,fill=type))+
  geom_density(size=1,alpha=0.7) +
  coord_cartesian(xlim = c(-5,8),ylim=c(0,0.4))+
  scale_fill_manual(values=c("#8c6100","#E69F00"),labels=c("New","Old"))+
  scale_x_continuous(breaks=c(0,3),labels=c(expression(paste(mu[n])),expression(paste(mu[o]))))+
  annotate("segment",x=0,xend=0,y=0,yend=.39,linetype="dotted")+
  annotate("segment",x=3,xend=3,y=0,yend=.27,linetype="dotted")+
  annotate("text",x=2,y=0.35,hjust=0,size=3.5,label=expression(paste("Free: ", mu[o], ", ",sigma[o])))+
  annotate("text",x=2,y=0.39,hjust=0,size=3.5,label=expression(paste("Fixed: ", mu[n]," = 0, ", sigma[n]," = 1")))+
  # geom_vline(xintercept = critsprop[[1]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[2]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[3]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[4]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[5]], linetype="dashed",color="grey") +
  facet_grid(.~model)+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),

        panel.background = element_rect(fill = "white"),
        legend.position = c(0.2,0.6),
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))

UVSDTmodel <-plot_grid(Ill_UVSDT,UVROC,nrow=2,rel_heights=c(1.5,1))

Ill_Gum <- ggplot(dats %>% filter(model=="Gumbel"),aes(x = value,fill=type))+
  geom_density(size=1,alpha=0.7) +
  coord_cartesian(xlim = c(-7,5),ylim=c(0,0.4))+
  scale_fill_manual(values=c("#005941","#009E73"),labels=c("New","Old"))+
  scale_x_continuous(breaks=c(0,3),labels=c(expression(paste(mu[n])),expression(paste(mu[o]))))+
  annotate("segment",x=0,xend=0,y=0,yend=.37,linetype="dotted")+
  annotate("segment",x=3,xend=3,y=0,yend=.37,linetype="dotted")+
  annotate("text",x=-5,y=0.35,hjust=0,size=3.5,label=expression(paste("Free: ", mu[o])))+
  annotate("text",x=-5,y=0.39,hjust=0,size=3.5,label=expression(paste("Fixed: ", mu[n]," = 0, ", beta[n]," = ",beta[o]," = 1")))+
  # geom_vline(xintercept = critsprop[[1]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[2]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[3]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[4]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[5]], linetype="dashed",color="grey") +
  facet_grid(.~model)+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),

        panel.background = element_rect(fill = "white"),
        legend.position = c(0.2,0.6),
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))


Gumbelmodel <-plot_grid(Ill_Gum,GumbelROC,nrow=2,rel_heights=c(1.5,1))

Ill_EG <- ggplot(dats %>% filter(model=="Ex-Gaussian"),aes(x = value,fill=type))+
  geom_density(size=1,alpha=0.7) +
  coord_cartesian(ylim=c(0,0.4))+

  scale_fill_manual(values=c("#2b5e7a","#56B4E9"),labels=c("New","Old"))+

  annotate("text",x=1.5,y=0.35,hjust=0,size=3.5,label=expression(paste("Free: ", mu[o[G]], ", ",beta[o]," = 1/",lambda[o])))+
  annotate("text",x=1.5,y=0.39,hjust=0,size=3.5,label=expression(paste("Fixed: ", mu[n]," = 0, ", sigma[n]," = ",sigma[o[G]]," = 1")))+
  scale_x_continuous(breaks=c(0,3-1.5),
                     limits = c(-5,11),
                     labels=c(expression(paste(mu[n])),expression(paste(mu[o[G]]))))+
  annotate("segment",x=0,xend=0,y=0,yend=.39,linetype="dotted")+
  annotate("segment",x=3-1.5,xend=3-1.5,y=0,yend=.2,linetype="dotted")+
  # geom_vline(xintercept = critsprop[[1]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[2]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[3]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[4]], linetype="dashed",color="grey") +
  # geom_vline(xintercept = critsprop[[5]], linetype="dashed",color="grey") +
  facet_grid(.~model)+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "transparent"),
        legend.position = c(0.7,0.6),
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))

ExGmodel <- plot_grid(Ill_EG, ExGROC,nrow=2,rel_heights=c(1.5,1))

#based on median par estimates from fits to Spanton & Berry

#plot_grid(Ill_UVSDT,Ill_Gum,Ill_EG,nrow=1)
plot_grid(UVSDTmodel,Gumbelmodel,ExGmodel,nrow=1)

ggsave("Figures/IllustrateModels2.png", units="cm", width=22, height=9, dpi=600)
