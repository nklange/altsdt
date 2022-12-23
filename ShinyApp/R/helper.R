#library(shinydashboard)
#library(plotly)

library(DT)
library(tidyverse)
library(cowplot)
library(tidybayes)
library(shinythemes)
library(shinyWidgets)
library("scales")


modelpred <- function(fit,model){

  if (model == "GaussianUVSDT"){
    fitU <- fit %>% filter(model == "GaussianUVSDT")
  ROC <- tibble(crit = seq(-8,8,.01)) %>%
    mutate(hit= pnorm(crit,mean=fitU$muo,sd=fitU$sigo,lower.tail=F),

           fa = pnorm(crit,lower.tail = F)) %>%
    mutate(threshcolor=NA,
           shapecolor=NA) %>%
    mutate(crit = as.character(round(crit,2))) %>%
    mutate(zhit = qnorm(hit),

           zfa = qnorm(fa))

  } else if (model == "GumbelEVSDT"){
    fitG <- fit %>% filter(model == "GumbelEVSDT")

    ROC <- tibble(crit = seq(-8,8,.01)) %>%
      mutate(hit= ordinal::pgumbel(crit,location=fitG$muo,scale=1,max=T),

             fa = ordinal::pgumbel(crit,max=T)) %>%
      mutate(threshcolor=NA,
             shapecolor=NA) %>%
      mutate(crit = as.character(round(crit,2))) %>%
      mutate(zhit = qnorm(hit),

             zfa = qnorm(fa))

  } else if(model == "ExGaussNormEVSDT"){
    fitE <- fit %>% filter(model == "ExGaussNormEVSDT")

    ROC <- tibble(crit = seq(-8,8,.01)) %>%
      mutate(hit= brms::pexgaussian(crit,mu=fitE$muo,sigma=1,beta = fitE$betao,lower.tail = F),

             fa =pnorm(crit,lower.tail = F)) %>%
      mutate(threshcolor=NA,
             shapecolor=NA) %>%
      mutate(crit = as.character(round(crit,2))) %>%
      mutate(zhit = qnorm(hit),

             zfa = qnorm(fa))

}

  return(ROC)
}

modelpred_mult <- function(fit,model){

  if (model == "baselineUVSDT"){
    fitU <- fit %>% filter(model == "baselineUVSDT")
    ROC_high_ev1 <- tibble(crit = seq(-8,8,.01)) %>%
      mutate(hit= pnorm(crit,mean=fitU$d_ho,sd=fitU$sig_ho,lower.tail=F),

             fa = pnorm(crit,lower.tail = F)) %>%
      mutate(threshcolor=NA,
             shapecolor=NA) %>%
      mutate(crit = as.character(round(crit,2))) %>%
      mutate(zhit = qnorm(hit),

             zfa = qnorm(fa)) %>%
      mutate(condition = "high_ev_1")

    ROC_low_ev1 <- tibble(crit = seq(-8,8,.01)) %>%
      mutate(hit= pnorm(crit,mean=fitU$d_lo,sd=fitU$sig_lo,lower.tail=F),

             fa = pnorm(crit,lower.tail = F)) %>%
      mutate(threshcolor=NA,
             shapecolor=NA) %>%
      mutate(crit = as.character(round(crit,2))) %>%
      mutate(zhit = qnorm(hit),

             zfa = qnorm(fa)) %>%
      mutate(condition = "low_ev_1")

    ROC_high_ev0 <- tibble(crit = seq(-8,8,.01)) %>%
      mutate(hit= pnorm(crit,mean=fitU$d_hn,sd=fitU$sig_hn,lower.tail=F),

             fa = pnorm(crit,lower.tail = F)) %>%
      mutate(threshcolor=NA,
             shapecolor=NA) %>%
      mutate(crit = as.character(round(crit,2))) %>%
      mutate(zhit = qnorm(hit),

             zfa = qnorm(fa)) %>%
      mutate(condition = "high_ev_0")

    ROC <- bind_rows(ROC_low_ev1,ROC_high_ev1,ROC_high_ev0)

  } else if (model == "GumbelEVSDT"){
    fitG <- fit %>% filter(model == "GumbelEVSDT")

    ROC_high_ev1 <- tibble(crit = seq(-8,8,.01)) %>%
      mutate(hit= ordinal::pgumbel(crit,location=fitG$d_ho,scale=1,max=T),

             fa = ordinal::pgumbel(crit,max=T)) %>%
      mutate(threshcolor=NA,
             shapecolor=NA) %>%
      mutate(crit = as.character(round(crit,2))) %>%
      mutate(zhit = qnorm(hit),

             zfa = qnorm(fa)) %>%
      mutate(condition = "high_ev_1")

    ROC_low_ev1 <- tibble(crit = seq(-8,8,.01)) %>%
      mutate(hit= ordinal::pgumbel(crit,location=fitG$d_lo,scale=1,max=T),

             fa = ordinal::pgumbel(crit,max=T)) %>%
      mutate(threshcolor=NA,
             shapecolor=NA) %>%
      mutate(crit = as.character(round(crit,2))) %>%
      mutate(zhit = qnorm(hit),

             zfa = qnorm(fa)) %>%
      mutate(condition = "low_ev_1")

    ROC_high_ev0 <- tibble(crit = seq(-8,8,.01)) %>%
      mutate(hit= ordinal::pgumbel(crit,location=fitG$d_hn,scale=1,max=T),

             fa = ordinal::pgumbel(crit,max=T)) %>%
      mutate(threshcolor=NA,
             shapecolor=NA) %>%
      mutate(crit = as.character(round(crit,2))) %>%
      mutate(zhit = qnorm(hit),

             zfa = qnorm(fa)) %>%
      mutate(condition = "high_ev_0")

    ROC <- bind_rows(ROC_low_ev1,ROC_high_ev1,ROC_high_ev0)


  } else if(model == "ExGaussNormEVSDT"){
    fitE <- fit %>% filter(model == "ExGaussNormEVSDT")

    ROC_high_ev1 <- tibble(crit = seq(-8,8,.01)) %>%
      mutate(hit= brms::pexgaussian(crit,mu=fitE$d_ho,sigma=1,beta = fitE$beta_ho,lower.tail = F),

             fa =pnorm(crit,lower.tail = F)) %>%
      mutate(threshcolor=NA,
             shapecolor=NA) %>%
      mutate(crit = as.character(round(crit,2))) %>%
      mutate(zhit = qnorm(hit),

             zfa = qnorm(fa)) %>%
      mutate(condition = "high_ev_1")

    ROC_low_ev1 <- tibble(crit = seq(-8,8,.01)) %>%
      mutate(hit= brms::pexgaussian(crit,mu=fitE$d_lo,sigma=1,beta = fitE$beta_lo,lower.tail = F),

             fa =pnorm(crit,lower.tail = F)) %>%
      mutate(threshcolor=NA,
             shapecolor=NA) %>%
      mutate(crit = as.character(round(crit,2))) %>%
      mutate(zhit = qnorm(hit),

             zfa = qnorm(fa)) %>%
      mutate(condition = "low_ev_1")

    ROC_high_ev0 <- tibble(crit = seq(-8,8,.01)) %>%
      mutate(hit= brms::pexgaussian(crit,mu=fitE$d_hn,sigma=1,beta = fitE$beta_hn,lower.tail = F),

             fa =pnorm(crit,lower.tail = F)) %>%
      mutate(threshcolor=NA,
             shapecolor=NA) %>%
      mutate(crit = as.character(round(crit,2))) %>%
      mutate(zhit = qnorm(hit),

             zfa = qnorm(fa)) %>%
      mutate(condition = "high_ev_0")
    ROC <- bind_rows(ROC_low_ev1,ROC_high_ev1,ROC_high_ev0)


  }

  return(ROC)
}

facet_strip_bigger <- function(gp, size){
  if(missing(gp)){
    print("this function needs a facet_wrap ggplotly object")
  }
  if(missing(size)){
    print("this function needs 'size' argument to be specified as integer. 80 will be introduced as default")
    size <- 80
  }

  n_facets <- c(1:length(gp[["x"]][["layout"]][["shapes"]]))

  for(i in n_facets){
    if(n_facets[i] %% 2 == 0){
      gp[["x"]][["layout"]][["shapes"]][[i]][["y0"]] <- + as.numeric(size)
      gp[["x"]][["layout"]][["shapes"]][[i]][["y1"]] <- 0
    }
  }

  return(gp)
}

makeROC <- function(exp,fit){

  if(unique(exp$exp) %in% c("SB2020_e1","SB2020_e2","SB2020_e3")){
    ex <- exp %>%
      mutate(rating = as.numeric(as.character(rating))) %>%
      ungroup() %>%
      select(exp,id,condition,oldnew,rating)

    allpred <- NULL
    for(i in unique(ex$condition)){
      for(model in c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT")){


        preds <- modelpred(fit %>% filter(condition == i),model) %>%
          mutate(model = model,
                 condition = i)

        allpred <- bind_rows(allpred,preds)
      }

    }

  } else if(unique(exp$exp) %in% c("SB2021_e1","SB2021_e2")){
  ex <- exp %>%
        mutate(rating = as.numeric(as.character(rating))) %>%
        ungroup() %>%
        select(exp,id,condition,oldnew,rating)

  allpred <- NULL
  for(i in c("A","B","C","D")){
    for(model in c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT")){


      preds <- modelpred(fit %>% filter(condition == i),model) %>%
        mutate(model = model,
               condition = i)

      allpred <- bind_rows(allpred,preds)
    }

  }

  } else if(unique(exp$exp) %in% c( "MWW2007_e1","MWW2007_e2","MWW2007_e3")){

    ex <- exp %>%
      mutate(rating = as.numeric(as.character(rating))) %>%
      ungroup() %>%
      select(exp,id,oldnew,rating) %>%
      mutate(condition = "-")

    allpred <- NULL

      for(model in c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT")){


        preds <- modelpred(fit,model) %>%
          mutate(model = model,
                 condition = "-")


        allpred <- bind_rows(allpred,preds)


    }

  } else if(unique(exp$exp) %in% c("SB2021_e3")){

  #   }
  #
  # } else {
  #
  #
  #
    expOld <- exp %>% filter(oldnew == 1)
    expNew <- exp %>% filter(oldnew == 0)
  #
    expInd <- list()
    for(i in seq_along(unique(expOld$condition))){

      expInd[[i]] <- expOld %>% filter(condition == unique(expOld$condition)[[i]]) %>%
        bind_rows(expNew) %>% mutate(condition = unique(expOld$condition)[[i]])
    }

    ex <- bind_rows(expInd) %>%
      ungroup() %>%
      select(exp,id,condition,oldnew,rating) %>%
      mutate(condition = factor(condition, levels=c("high_ev_1","low_ev_1","high_ev_0"),
                                labels = c("High Var, Old","Low Var, Old","High Var, New")))
  #
  # }

    allpred <- NULL

      for(model in c("baselineUVSDT","GumbelEVSDT","ExGaussNormEVSDT")){


        preds <- modelpred_mult(fit,model) %>%
          mutate(model = ifelse(model == "baselineUVSDT","GaussianUVSDT",model))

        allpred <- bind_rows(allpred,preds)


      }

    allpred <- allpred %>%
      mutate(condition = factor(condition, levels=c("high_ev_1","low_ev_1","high_ev_0"),
                                labels = c("High Var, Old","Low Var, Old","High Var, New")))

  } else if(unique(exp$exp) %in% c("D2007_e1a","D2007_e1b")) {

    expOld <- exp %>% filter(isold == 1)
    expNew <- exp %>% filter(isold == 0)
    #
    expInd <- list()
    for(i in seq_along(unique(expOld$condition))){

      expInd[[i]] <- expOld %>% filter(condition == unique(expOld$condition)[[i]]) %>%
        bind_rows(expNew) %>% mutate(condition = unique(expOld$condition)[[i]])
    }

    ex <- bind_rows(expInd) %>%
      ungroup() %>%
      select(exp,id,condition,isold,rating) %>%
      mutate(oldnew = isold) %>%
      mutate(condition = factor(condition, levels=c("high_ev_1","low_ev_1","high_ev_0"),
                                labels = c("High Freq, Old","Low Freq, Old","Low Freq, New")))
    #
    # }

    allpred <- NULL

    for(model in c("baselineUVSDT","GumbelEVSDT","ExGaussNormEVSDT")){


      preds <- modelpred_mult(fit,model) %>%
        mutate(model = ifelse(model == "baselineUVSDT","GaussianUVSDT",model))

      allpred <- bind_rows(allpred,preds)


    }

    allpred <- allpred %>%
      mutate(condition = factor(condition, levels=c("high_ev_1","low_ev_1","high_ev_0"),
                                labels = c("High Freq, Old","Low Freq, Old","Low Freq, New")))

}
  makerocdata <- ex %>%
    group_by(exp,id,condition,oldnew,rating) %>%
    summarize(cumresp = length(rating)) %>%
    mutate(rating = factor(rating, levels = c(max(ex$rating):1))) %>%
    arrange(exp,id,condition,oldnew,rating) %>%
    mutate(oldnew = factor(oldnew,levels=c(0,1),labels=c("INew","IOld"))) %>%
    group_by(exp,id,condition,oldnew) %>%
    complete(rating,fill = list(cumresp = 0)) %>%
    summarize(cumprop = cumsum(cumresp/sum(cumresp)),
              rating = rating) %>%
    ungroup() %>%
    pivot_wider(names_from="oldnew",values_from="cumprop",
                id_cols = c("exp","id","condition","rating")) %>%
    filter(rating != 1)

  makezrocdata <- ex %>%
    group_by(exp,id,condition,oldnew,rating) %>%
    summarize(cumresp = length(rating)) %>%
    mutate(rating = factor(rating, levels = c(max(ex$rating):1))) %>%
    arrange(exp,id,condition,oldnew,rating) %>%
    mutate(oldnew = factor(oldnew,levels=c(0,1),labels=c("INew","IOld"))) %>%
    group_by(exp,id,condition,oldnew) %>%
    complete(rating,fill = list(cumresp = 0)) %>%
    summarize(cumprop = qnorm(cumsum(cumresp/sum(cumresp))),
              rating = rating) %>%
    ungroup() %>%
    pivot_wider(names_from="oldnew",values_from="cumprop",
                id_cols = c("exp","id","condition","rating")) %>%
    filter(rating != 1)





  ROC <- ggplot(makerocdata,aes(x = INew,y=IOld))+
    facet_wrap(.~condition,ncol=4)+


    geom_line(data = allpred %>% filter(model == "GaussianUVSDT"),
              aes(x=fa,y=hit),color="#E69F00",size=1)+

    geom_line(data = allpred %>% filter(model == "GumbelEVSDT"),
              aes(x=fa,y=hit),color="#009E73",size=1)+
    geom_line(data = allpred %>% filter(model == "ExGaussNormEVSDT"),
              aes(x=fa,y=hit),color="#CC79A7",size=1)+
    geom_point(data=makerocdata,aes(x = INew,y=IOld),size=5)+
    geom_abline(intercept=0,slope=1,color="blue",linetype="dashed")+
    theme_bw()+
    coord_fixed(xlim=c(0,1),ylim=c(0,1))+
    #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+
    scale_shape_manual(values=c(16,17),name="Data")+
    scale_x_continuous(name="p('old'|new)")+
    scale_y_continuous(name = "p('old'|old)")+
    theme_bw()+
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14),
          strip.background = element_rect(fill="black"),
          strip.text = element_text(color="white",face = "bold",size=16))


  zROC <- ggplot(makezrocdata,aes(x = INew,y=IOld))+
    facet_wrap(.~condition,ncol=4)+


    geom_line(data = allpred %>% filter(model == "GaussianUVSDT"),
              aes(x=zfa,y=zhit),color="#E69F00",size=1)+

    geom_line(data = allpred %>% filter(model == "GumbelEVSDT"),
              aes(x=zfa,y=zhit),color="#009E73",size=1)+
    geom_line(data = allpred %>% filter(model == "ExGaussNormEVSDT"),
              aes(x=zfa,y=zhit),color="#CC79A7",size=1)+
    geom_point(data=makezrocdata,aes(x = INew,y=IOld),size=5)+
    geom_abline(intercept=0,slope=1,color="blue",linetype="dashed")+
    theme_bw()+
    coord_fixed(xlim=c(-2.5,2.5),ylim=c(-2.5,2.5))+
    #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+
    scale_shape_manual(values=c(16,17),name="Data")+
    scale_x_continuous(name=expression(paste(phi^-1,"[p('old'|new)]")))+
    scale_y_continuous(name=expression(paste(phi^-1,"[p('old'|old)]")))+
    theme_bw()+
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14),
          strip.background = element_rect(fill="black"),
          strip.text = element_text(color="white",face = "bold",size=16))

  plot_grid(ROC,zROC,nrow=2,rel_heights=c(0.5,0.5),align="hv")
}
makeDeviance <- function(fit,Gsq){

  if(unique(fit$exp) %in% c("SB2021_e1","SB2021_e2","SB2020_e1","SB2020_e2","SB2020_e3")){

  deltas <- fit %>%
    mutate(AIC = 2 * objective + 2 * npar) %>%
    group_by(condition) %>%
    mutate(deltaAICt = AIC - min(AIC)) %>%
    mutate(model = factor(model,levels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT")))%>%
    mutate(deltaAIC = ifelse(deltaAICt >= 20, 20,deltaAICt))


ggplot(deltas,aes(x = model,y=deltaAIC,fill=model))+
    facet_wrap(.~condition,ncol=4)+

  geom_bar(stat = "identity")+
  annotate("text",x=2,y=19,size=5,label=expression(paste(italic(p)[G^2])))+
  scale_x_discrete(labels=c("U-V Gaussian","Gumbel","Ex-Gaussian"))+
  scale_y_continuous(limits=c(-2,20),name=expression(paste(Delta,"AIC")),breaks=seq(0,20,5))+
  theme_bw()+
  scale_fill_manual(values = c("#E69F00","#009E73","#CC79A7"))+
  #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+

  geom_text(data = Gsq, aes(x = model, y = 18,
                            label=ifelse(pGSq < .01,"0",
                                         paste0(sub("^0+", "", round(pGSq,2))))),
            size=5)+
  geom_text(data=deltas,aes(x=model,y=-2,label=round(deltaAICt,0)))+

    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=14),
          legend.position = "none",
          axis.title.y = element_text(size=14),
          axis.title.x = element_blank(),
          strip.background = element_rect(fill="black"),
          strip.text = element_text(color="white",face = "bold",size=16))

} else if(unique(fit$exp) %in% c("MWW2007_e1","MWW2007_e2","MWW2007_e3")){

  deltas <- fit %>%
    mutate(AIC = 2 * objective + 2 * npar) %>%
    mutate(deltaAICt = AIC - min(AIC)) %>%
    mutate(model = factor(model,levels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"))) %>%
    mutate(condition = "-") %>%
    mutate(deltaAIC = ifelse(deltaAICt >= 20, 20,deltaAICt))


  ggplot(deltas,aes(x = model,y=deltaAIC,fill=model))+
    facet_wrap(.~condition,ncol=4)+


    geom_bar(stat = "identity")+
    annotate("text",x=2,y=19,size=5,label=expression(paste(italic(p)[G^2])))+
    scale_x_discrete(labels=c("U-V Gaussian","Gumbel","Ex-Gaussian"))+
    scale_y_continuous(limits=c(-2,20),name=expression(paste(Delta,"AIC")),breaks=seq(0,20,5))+
    theme_bw()+
    scale_fill_manual(values = c("#E69F00","#009E73","#CC79A7"))+
    #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+

    geom_text(data = Gsq %>% mutate(condition = "-"), aes(x = model, y = 18,
                              label=ifelse(pGSq < .01,"0",
                                           paste0(sub("^0+", "", round(pGSq,2))))),
              size=5)+
    geom_text(data=deltas,aes(x=model,y=-2,label=round(deltaAICt,0)))+

    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=14),
          legend.position = "none",
          axis.title.y = element_text(size=14),
          axis.title.x = element_blank(),
          strip.background = element_rect(fill="black"),
          strip.text = element_text(color="white",face = "bold",size=16))

  } else {

    deltas <- fit %>%
      mutate(AIC = 2 * objective + 2 * npar) %>%
      #group_by(condition) %>%
      mutate(deltaAICt = AIC - min(AIC)) %>%
      mutate(model = factor(model,levels=c("baselineUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                            labels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"))) %>%
      mutate(deltaAIC = ifelse(deltaAICt >= 20, 20,deltaAICt))


    ggplot(deltas,aes(x = model,y=deltaAIC,fill=model))+
      #facet_wrap(.~condition,ncol=4)+


      geom_bar(stat = "identity")+
      annotate("text",x=2,y=19,size=5,label=expression(paste(italic(p)[G^2])))+
      scale_x_discrete(labels=c("U-V Gaussian","Gumbel","Ex-Gaussian"))+
      scale_y_continuous(limits=c(-2,20),name=expression(paste(Delta,"AIC")),breaks=seq(0,20,5))+
      theme_bw()+
      scale_fill_manual(values = c("#E69F00","#009E73","#CC79A7"))+
      #ggtitle(paste0(unique(exp$exp)," (Manipulation: ",unique(exp$manipulation),"; N = ", length(unique(exp$id)),")"))+

      geom_text(data = Gsq, aes(x = model, y = 18,
                                label=ifelse(pGSq < .01,"0",
                                             paste0(sub("^0+", "", round(pGSq,2))))),
                size=5)+
      geom_text(data=deltas,aes(x=model,y=-2,label=round(deltaAICt,0)))+
      theme_bw()+
      theme(axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=14),
            legend.position = "none",
            axis.title.y = element_text(size=14),
            axis.title.x = element_blank(),
            strip.background = element_rect(fill="black"),
            strip.text = element_text(color="white",face = "bold",size=16))

}

}
makeFreqprediction <- function(exp,Pred){

  if(unique(exp$exp) %in% c("SB2020_e1","SB2020_e2","SB2020_e3")){


    ex <- exp %>%
      mutate(rating = as.numeric(as.character(rating))) %>%
      ungroup() %>%
      select(exp,id,condition,oldnew,rating) %>%
      group_by(condition,oldnew,rating) %>%
      summarize(Freq = length(rating)) %>%
      ungroup() %>%
      mutate(condition = factor(condition),
             oldnew = factor(oldnew,levels = c(0,1)),
             rating = factor(rating, levels = c(1:6)))

    fullresp <- ex %>% expand(condition,oldnew,rating)

    prepdat <- ex %>% right_join(fullresp) %>%
      #mutate(Freq = ifelse(is.na(Freq),1/6,Freq)) %>%
      mutate(Freq = ifelse(is.na(Freq),0,Freq)) %>%
      arrange(condition,oldnew,rating) %>%
      group_by(condition,oldnew) %>%
      mutate(Freq  = Freq/sum(Freq))

    colnames(Pred) <- c(paste0("V",1:12),"id","condition","model","exp")

    residuals <- Pred %>%
      # mutate(across(contains("V"), ~. * 60)) %>%
      pivot_longer(cols = contains("V"),names_to = "response",values_to="Freq") %>%
      mutate(oldnew = rep(rep(rep(c(0,1),each=6),2),3),
             rating = rep(rep(rep(c(1:6),2),2),3)) %>%
      group_by(model) %>%
      nest() %>%
      mutate(relfreq = map(data, ~.$Freq - prepdat$Freq)) %>%
      unnest() %>%
      mutate(model = factor(model,levels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"))) %>%
      mutate(oldnew = factor(oldnew,levels=c(0,1),
                             labels=c("New","Old")))

    propdat <- prepdat %>%
      mutate(oldnew = factor(oldnew,levels=c(0,1),
                             labels=c("New","Old")))

  } else if(unique(exp$exp) %in% c("SB2021_e1","SB2021_e2")){


  ex <- exp %>%
    mutate(rating = as.numeric(as.character(rating))) %>%
    ungroup() %>%
    select(exp,id,condition,oldnew,rating) %>%
    group_by(condition,oldnew,rating) %>%
    summarize(Freq = length(rating)) %>%
    ungroup() %>%
    mutate(condition = factor(condition, levels = c("A","B","C","D")),
           oldnew = factor(oldnew,levels = c(0,1)),
           rating = factor(rating, levels = c(1:6)))

  fullresp <- ex %>% expand(condition,oldnew,rating)

  prepdat <- ex %>% right_join(fullresp) %>%
    #mutate(Freq = ifelse(is.na(Freq),1/6,Freq)) %>%
    mutate(Freq = ifelse(is.na(Freq),0,Freq)) %>%
    arrange(condition,oldnew,rating) %>%
    group_by(condition,oldnew) %>%
    mutate(Freq  = Freq/sum(Freq))

  colnames(Pred) <- c(paste0("V",1:12),"id","condition","model","exp")

  residuals <- Pred %>%
   # mutate(across(contains("V"), ~. * 60)) %>%
    pivot_longer(cols = contains("V"),names_to = "response",values_to="Freq") %>%
    mutate(oldnew = rep(rep(rep(c(0,1),each=6),4),3),
           rating = rep(rep(rep(c(1:6),2),4),3)) %>%
    group_by(model) %>%
    nest() %>%
    mutate(relfreq = map(data, ~.$Freq - prepdat$Freq)) %>%
    unnest() %>%
    mutate(model = factor(model,levels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"))) %>%
    mutate(oldnew = factor(oldnew,levels=c(0,1),
                           labels=c("New","Old")))

  propdat <- prepdat %>%
    mutate(oldnew = factor(oldnew,levels=c(0,1),
                           labels=c("New","Old")))

  } else if(unique(exp$exp) == "SB2021_e3") {

    ex <- exp %>%
      mutate(rating = as.numeric(as.character(rating))) %>%
      ungroup() %>%
      select(exp,id,condition,oldnew,rating) %>%
      group_by(condition,oldnew,rating) %>%
      summarize(Freq = length(rating)) %>%
      ungroup() %>%
      mutate(condition = factor(condition, levels = c("high_ev_1","low_ev_1","high_ev_0","low_ev_0")),
             oldnew = factor(oldnew,levels = c(0,1)),
             rating = factor(rating, levels = c(1:6)))

    fullresp <- ex %>% expand(condition,rating)

    prepdat <- ex %>% right_join(fullresp) %>%
      mutate(oldnew = ifelse(condition %in% c("high_ev_1","low_ev_1"),1,0))  %>%
      #mutate(Freq = ifelse(is.na(Freq),1/6,Freq)) %>%
      mutate(Freq = ifelse(is.na(Freq),0,Freq)) %>%
      mutate(Freq = ifelse(is.nan(Freq),0,Freq)) %>%
      arrange(condition,oldnew,rating) %>%
      group_by(condition,oldnew) %>%
      mutate(Freq  = Freq/sum(Freq))

    colnames(Pred) <- c(paste0("V",1:24),"id","condition","model","exp")

    residuals <- Pred %>%
      # mutate(across(contains("V"), ~. * 60)) %>%
      pivot_longer(cols = contains("V"),names_to = "response",values_to="Freq") %>%
      mutate(oldnew = rep(rep(c(1,1,0,0),each=6),3),
             condition = rep(rep(c("high_ev_1","low_ev_1","high_ev_0","low_ev_0"),each=6),3),
             rating = rep(rep(c(1:6),4),3)) %>%
      group_by(model) %>%
      nest() %>%
      mutate(relfreq = map(data, ~.$Freq - prepdat$Freq)) %>%
      unnest() %>%
      mutate(model = factor(model,levels=c("baselineUVSDT","GumbelEVSDT","ExGaussNormEVSDT"))) %>%
      mutate(oldnew = ifelse(condition=="high_ev_0",0,oldnew)) %>%
      mutate(oldnew = factor(oldnew,levels=c(0,1),
                             labels=c("New","Old"))) %>%
      mutate(conditiont = condition) %>%
      mutate(condition = ifelse(grepl("high",condition),"High Var","Low Var"))

    propdat <- prepdat %>%
      mutate(oldnew = as.character(oldnew)) %>%
      mutate(oldnew = ifelse(condition=="high_ev_0","0",oldnew)) %>%
      mutate(oldnew = factor(oldnew,levels=c(0,1),
                             labels=c("New","Old"))) %>%
      mutate(conditiont = condition) %>%
      mutate(condition = ifelse(grepl("high",condition),"High Var","Low Var"))

  }else if(unique(exp$exp) %in% c("MWW2007_e1","MWW2007_e2","MWW2007_e3")){


    ex <- exp %>%
      mutate(rating = as.numeric(as.character(rating))) %>%
      ungroup() %>%
      select(exp,id,oldnew,rating) %>%
      group_by(oldnew,rating) %>%
      summarize(Freq = length(rating)) %>%
      ungroup() %>%
      mutate(oldnew = factor(oldnew,levels = c(0,1)),
             rating = factor(rating, levels = c(1:6)))

    fullresp <- ex %>% expand(oldnew,rating)



    prepdat <- ex %>% right_join(fullresp) %>%
      #mutate(Freq = ifelse(is.na(Freq),1/6,Freq)) %>%
      mutate(Freq = ifelse(is.na(Freq),0,Freq)) %>%
      arrange(oldnew,rating) %>%
      group_by(oldnew) %>%
      mutate(Freq  = Freq/sum(Freq)) %>%
      mutate(condition = "-")

    colnames(Pred) <- c(paste0("V",1:12),"id","condition","model","exp")

    residuals <- Pred %>%

      # mutate(across(contains("V"), ~. * 60)) %>%
      pivot_longer(cols = contains("V"),names_to = "response",values_to="Freq") %>%
      mutate(oldnew = rep(rep(c(0,1),each=6),3),
             rating = rep(rep(c(1:6),2),3)) %>%
      group_by(model) %>%
      nest() %>%
      mutate(relfreq = map(data, ~.$Freq - prepdat$Freq)) %>%
      unnest() %>%
      mutate(model = factor(model,levels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"))) %>%
      mutate(oldnew = factor(oldnew,levels=c(0,1),
                             labels=c("New","Old"))) %>%
        mutate(condition = "-")

    propdat <- prepdat %>%
      mutate(oldnew = factor(oldnew,levels=c(0,1),
                             labels=c("New","Old"))) %>%
      mutate(condition = "-")

    } else if(unique(exp$exp) %in% c("D2007_e1a","D2007_e1b")) {

    ex <- exp %>%
      mutate(rating = as.numeric(as.character(rating))) %>%
      ungroup() %>%
      mutate(oldnew = isold) %>%
      select(exp,id,condition,oldnew,rating) %>%
      group_by(condition,oldnew,rating) %>%
      summarize(Freq = length(rating)) %>%
      ungroup() %>%
      mutate(condition = factor(condition, levels = c("high_ev_1","low_ev_1","high_ev_0","low_ev_0")),
             oldnew = factor(oldnew,levels = c(0,1)),
             rating = factor(rating, levels = c(1:6)))

    fullresp <- ex %>% expand(condition,rating)

    prepdat <- ex %>% right_join(fullresp) %>%
      mutate(oldnew = ifelse(condition %in% c("high_ev_1","low_ev_1"),1,0))  %>%
      #mutate(Freq = ifelse(is.na(Freq),1/6,Freq)) %>%
      mutate(Freq = ifelse(is.na(Freq),0,Freq)) %>%
      mutate(Freq = ifelse(is.nan(Freq),0,Freq)) %>%
      arrange(condition,oldnew,rating) %>%
      group_by(condition,oldnew) %>%
      mutate(Freq  = Freq/sum(Freq))

    colnames(Pred) <- c(paste0("V",1:24),"id","condition","model","exp")

    residuals <- Pred %>%
      # mutate(across(contains("V"), ~. * 60)) %>%
      pivot_longer(cols = contains("V"),names_to = "response",values_to="Freq") %>%
      mutate(oldnew = rep(rep(c(1,1,0,0),each=6),3),
             condition = rep(rep(c("high_ev_1","low_ev_1","high_ev_0","low_ev_0"),each=6),3),
             rating = rep(rep(c(1:6),4),3)) %>%
      group_by(model) %>%
      nest() %>%
      mutate(relfreq = map(data, ~.$Freq - prepdat$Freq)) %>%
      unnest() %>%
      mutate(model = factor(model,levels=c("baselineUVSDT","GumbelEVSDT","ExGaussNormEVSDT"))) %>%
      mutate(oldnew = ifelse(condition=="high_ev_0",0,oldnew)) %>%
      mutate(oldnew = factor(oldnew,levels=c(0,1),
                             labels=c("New","Old"))) %>%
      mutate(conditiont = case_when(condition == "high_ev_0" ~ "low_ev_0",
                                    condition == "low_ev_0" ~ "high_ev_0",
                                    condition == "high_ev_1" ~ "high_ev_1",
                                    condition == "low_ev_1" ~ "low_ev_1",
                                    TRUE~ "NA")) %>%
      mutate(condition = ifelse(grepl("high",conditiont),"High Freq","Low Freq"))

    propdat <- prepdat %>%
      mutate(oldnew = as.character(oldnew)) %>%
      mutate(oldnew = ifelse(condition=="high_ev_0","0",oldnew)) %>%
      mutate(oldnew = factor(oldnew,levels=c(0,1),
                             labels=c("New","Old"))) %>%
      mutate(conditiont = case_when(condition == "high_ev_0" ~ "low_ev_0",
                                    condition == "low_ev_0" ~ "high_ev_0",
                                    condition == "high_ev_1" ~ "high_ev_1",
                                    condition == "low_ev_1" ~ "low_ev_1",
                                    TRUE~ "NA")) %>%
      mutate(condition = ifelse(grepl("high",conditiont),"High Freq","Low Freq"))

  }


  propresp <- ggplot(propdat
             ,aes(x=rating,y=Freq))+
    facet_grid(oldnew ~ condition)+
    #geom_hline(yintercept=0)+
    geom_bar(stat = "identity")+
    geom_point(data = residuals,aes(x = rating,y=Freq,color=model),size=5,
               position = position_dodge(.5))+
    scale_y_continuous(name="Proportion of responses")+
    scale_x_discrete(labels=c(1:6))+
    scale_color_manual(values = c("#E69F00","#009E73","#CC79A7"),
                       labels=c("U-V Gaussian","Gumbel","Ex-Gaussian"))+
    #
    theme_bw()+
    theme(axis.text = element_text(size=14),
          #legend.position = "none",
          axis.title.y = element_text(size=14),
          axis.title.x = element_blank(),
          strip.background = element_rect(fill="black"),
          strip.text = element_text(color="white",face = "bold",size=16))

  res <- ggplot(residuals,aes(x=rating,y=relfreq,color=model))+
    facet_grid(oldnew ~ condition)+
    geom_hline(yintercept=0)+
    geom_line(size=1)+
    geom_point(size=5)+
    scale_y_continuous(name="Residuals")+
    scale_x_continuous(breaks=c(1:6))+
    scale_color_manual(values = c("#E69F00","#009E73","#CC79A7"),
                       labels=c("U-V Gaussian","Gumbel","Ex-Gaussian"))+

    theme_bw()+
    theme(axis.text = element_text(size=14),
          #legend.position = "none",
          axis.title.y = element_text(size=14),
          axis.title.x = element_blank(),
          strip.background = element_rect(fill="black"),
          strip.text = element_text(color="white",face = "bold",size=16))

  plot_grid(propresp,res,rel_heights=c(0.5,0.5),align="hv",ncol=1)

}

makePar <- function(fit){

if(unique(fit$exp) %in% c("MWW2007_e1","MWW2007_e2","MWW2007_e3")){
tabs <- fit %>%
  select(model,muo,sigo,betao,c1,dc1,dc2,dc3,dc4) %>%
  mutate(across(where(is.numeric),~round(.,2))) %>%
  mutate(model= factor(model,levels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                       labels=c("U-V Gaussian","Gumbel","Ex-Gaussian"))) %>%
  arrange(model)
colnames(tabs) <- c("Model", "d",
                    "sig_o","beta_o", "c1","dc1","dc2","dc3","dc4")

} else if(unique(fit$exp) %in% c("SB2021_e1","SB2021_e2","SB2020_e1","SB2020_e2","SB2020_e3")){

  tabs <- fit %>% ungroup() %>%
    select(condition,model,muo,sigo,betao,c1,dc1,dc2,dc3,dc4) %>%
    mutate(across(where(is.numeric),~round(.,2))) %>%
    mutate(model= factor(model,levels=c("GaussianUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                         labels=c("U-V Gaussian","Gumbel","Ex-Gaussian"))) %>%
    arrange(condition,model)
  colnames(tabs) <- c("Condition","Model", "d",
                      "sig_o","beta_o", "c1","dc1","dc2","dc3","dc4")


} else if(unique(fit$exp) %in% c("SB2021_e3")){

  tabs <- fit %>% ungroup() %>%
    select(model,d_ho,d_lo,d_hn,sig_ho,sig_lo,sig_hn,
           beta_ho, beta_lo, beta_hn,c1,dc1,dc2,dc3,dc4) %>%
    mutate(across(where(is.numeric),~round(.,2))) %>%
    mutate(model= factor(model,levels=c("baselineUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                         labels=c("U-V Gaussian","Gumbel","Ex-Gaussian"))) %>%
    arrange(model)

}  else if(unique(fit$exp) %in% c("D2007_e1a","D2007_e1b")){

  tabs <- fit %>% ungroup() %>%
    select(model,d_ho,d_lo,d_hn,sig_ho,sig_lo,sig_hn,
           beta_ho, beta_lo, beta_hn,c1,dc1,dc2,dc3,dc4) %>%
    mutate(across(where(is.numeric),~round(.,2))) %>%
    mutate(model= factor(model,levels=c("baselineUVSDT","GumbelEVSDT","ExGaussNormEVSDT"),
                         labels=c("U-V Gaussian","Gumbel","Ex-Gaussian"))) %>%
    arrange(model) %>%
    rename("d_ln" = d_hn,
           "sig_ln" = sig_hn,
           "beta_ln" = beta_hn)



}


return(tabs)
}
