# General information --------------------------

# Functions to fit SDT models to data sets with 1 set of new items
# and multiple sets of old items

library(tidyverse)


prep_data <- function(data,freqdat,numcategories){

  if(freqdat){

    prepdat <- data %>%
      mutate(Freq = ifelse(is.na(Freq) | Freq==0,1/numcategories,Freq)) %>%
      arrange(condition,rating)

  } else {

    preprepdat <- data %>% group_by(condition,rating) %>%
      summarize(Freq = length(rating)) %>%
      ungroup() %>%
      mutate(condition = factor(condition,levels = c("high_ev_1","low_ev_1",
                                                     "high_ev_0" ,"low_ev_0")),
             rating = factor(rating, levels = c(1:numcategories)))

    fullresp <- preprepdat %>% tidyr::expand(condition,rating)

    prepdat <- preprepdat %>% right_join(fullresp) %>%
      mutate(Freq = ifelse(is.na(Freq),1/numcategories,Freq)) %>%
     # mutate(Freq = ifelse(is.na(Freq),0,Freq)) %>%
      arrange(condition,rating)

  }

  dl <- split(prepdat$Freq, f = prepdat$condition)
  trials_LN <- sum(prepdat %>% filter(condition == "low_ev_0") %>% .$Freq)
  trials_HN <- sum(prepdat %>% filter(condition == "high_ev_0") %>% .$Freq)
  trials_LO <- sum(prepdat %>% filter(condition == "low_ev_1") %>% .$Freq)
  trials_HO <- sum(prepdat %>% filter(condition == "high_ev_1") %>% .$Freq)

  id <- unique(data$id)
  confidence <- sort(unique(prepdat$rating))
  if (is.null(unique(data$cvid))) {
    cvid <- id
  } else {
    cvid <- unique(data$cvid)
  }
  if (is.null(unique(data$exp))) {
    exp <- NA
  } else {
    exp <- unique(data$exp)
  }
  if (is.null(unique(data$leftout))) {
    leftout <- 0
  }else {
    leftout <- unique(data$leftout)
  }


  return(list(datalist = dl, confidence = confidence,
              trials_LN = trials_LN,
              trials_HN = trials_HN,
              trials_LO = trials_LO,
              trials_HO = trials_HO,
              id = id,
              cvid = cvid, exp = exp, leftout = leftout))

}

get_par_limits <- function(model,numcategories){

  if (model %in% c("baselineUVSDT")){
    # "d"    "sigo" "c1"   "dc1"  "dc2"  "dc3"  "dc4"

    lower <- c(-Inf,-Inf,-Inf,
               .Machine$double.eps,.Machine$double.eps,.Machine$double.eps,
               -Inf, rep(.Machine$double.eps,(numcategories-2)))
    #lower <- c(-Inf, .Machine$double.eps,rep(-Inf,(numcategories-1)))
    upper <- Inf

  } else if (model %in% "GumbelEVSDT"){

    lower <- c(-Inf,-Inf,-Inf,

               -Inf, rep(.Machine$double.eps,(numcategories-2)))
    #lower <- c(-Inf, .Machine$double.eps,rep(-Inf,(numcategories-1)))
    upper <- Inf

  } else if (model %in% c("ExGaussNormEVSDT")){

    lower <- c(-Inf,-Inf,-Inf,
               .Machine$double.eps,.Machine$double.eps,.Machine$double.eps,
               -Inf, rep(.Machine$double.eps,(numcategories-2)))
    #lower <- c(-Inf, .Machine$double.eps,rep(-Inf,(numcategories-1)))
    upper <- Inf

  }
  list(lower,upper)
}

get_start_par <- function(model,numcategories){


  #crit <- as_tibble(t(sort(runif(numcategories-1, min= -3, max=3)))) %>% set_names(paste0("c",c(1:(numcategories-1))))

  if (model %in% c("baselineUVSDT")){

    dc <- as_tibble(t(runif(numcategories-2, min= .Machine$double.eps, max=0.3))) %>% set_names(paste0("dc",c(1:(numcategories-2))))
    pstartunif <- tibble(
      d_ho   = runif(1, min=0, max=3),
      d_lo = runif(1, min=0, max=3),
      d_hn = runif(1, min=0, max=1),
     # d_ln = runif(1, min=0, max=1),
      sig_ho = runif(1, min=1, max=4),
     sig_lo = runif(1, min=1, max=4),
     sig_hn = runif(1, min=1, max=4),
      #crit
      c1   = runif(1, min=-3, max=0),
      dc
    )

    pstart <- bind_rows(pstartunif)

  } else if (model %in% "GumbelEVSDT"){

    dc <- as_tibble(t(runif(numcategories-2, min= .Machine$double.eps, max=0.3))) %>% set_names(paste0("dc",c(1:(numcategories-2))))
    pstartunif <- tibble(
      d_ho   = runif(1, min=-3, max=0),
      d_lo = runif(1, min=-3, max=0),
       d_hn = runif(1, min=-3, max=0),
      # d_ln = runif(1, min=0, max=1), =d_hn = 0
      #sig_ho = runif(1, min=1, max=4),
      #sig_lo = runif(1, min=1, max=4), = sig_ho
     # sig_hn = runif(1, min=1, max=4), = sig_ln = 1
      #crit
      c1   = runif(1, min=-3, max=0),
      dc
    )

    pstart <- bind_rows(pstartunif)

  }else if (model %in% "ExGaussNormEVSDT"){

    dc <- as_tibble(t(runif(numcategories-2, min= .Machine$double.eps, max=0.3))) %>% set_names(paste0("dc",c(1:(numcategories-2))))
    pstartunif <- tibble(
      d_ho   = runif(1, min=0, max=3),
      d_lo = runif(1, min=0, max=3),
      d_hn = runif(1, min=0, max=3),
      beta_ho =  runif(1, min=0.1, max=3),
      beta_lo =  runif(1, min=0.1, max=3),
      beta_hn =  runif(1, min=0.1, max=3),
      # d_ln = runif(1, min=0, max=1), =d_hn = 0
      #sig_ho = runif(1, min=1, max=4),
      #sig_lo = runif(1, min=1, max=4), = sig_ho
      # sig_hn = runif(1, min=1, max=4), = sig_ln = 1
      #crit
      c1   = runif(1, min=-3, max=0),
      dc
    )

    pstart <- bind_rows(pstartunif)

  }
  return(pstartunif)
}
get_start_par2 <- function(model,datalist){


  H_hi <-  sum(datalist$high_ev_1[4:6])/ sum(datalist$high_ev_1)
  H_lo <-  sum(datalist$low_ev_1[4:6]) / sum(datalist$low_ev_1)
  Hn_hi <-  sum(datalist$high_ev_0[4:6])/ sum(datalist$high_ev_0)
  FA_lo <- sum(datalist$low_ev_0[4:6]) / sum(datalist$low_ev_0)

  # Estimate mu old for both types of old item distribution
  est_muo_hi <- qnorm(H_hi) - qnorm(FA_lo)
  est_muo_lo <- qnorm(H_lo) - qnorm(FA_lo)
  est_mun_hi <- qnorm(Hn_hi) - qnorm(FA_lo)

  # Take zROC slope as estimate of sigo low and high
  cProbO_hi <- rev(cumsum(rev( datalist$high_ev_1 / sum(datalist$high_ev_1) )))[2:6]
  cProbO_lo <- rev(cumsum(rev(datalist$low_ev_1 / sum(datalist$low_ev_1)  )))[2:6]
  cProbN_lo <- rev(cumsum(rev( datalist$low_ev_0 / sum(datalist$low_ev_0)  )))[2:6]
  cProbN_hi <- rev(cumsum(rev( datalist$high_ev_0  / sum(datalist$high_ev_0)  )))[2:6]

  # Estimate sigmas from zROC, or use a fixed estimate if zROC calculation isn't possible
  if (any(c(cProbO_hi, cProbN_lo) == 0) | any(c(cProbO_hi, cProbN_lo) == 1)) {
    est_sigo_hi <- 1/0.8
  } else {
    est_sigo_hi <- 1 / round(lm(qnorm(cProbO_hi) ~ qnorm(cProbN_lo))$coef[2], digits = 2)
  }

  if (any(c(cProbO_lo, cProbN_lo) == 0) | any(c(cProbO_lo, cProbN_lo) == 1)) {
    est_sigo_lo <- 1/0.8
  } else {
    est_sigo_lo <- 1 / round(lm(qnorm(cProbO_lo) ~ qnorm(cProbN_lo))$coef[2], digits = 2)
  }

  if (any(c(cProbN_hi, cProbN_lo) == 0) | any(c(cProbN_hi, cProbN_lo) == 1)) {
    est_sign_hi <- 1/0.8
  } else {
    est_sign_hi <- 1 / round(lm(qnorm(cProbN_hi) ~ qnorm(cProbN_lo))$coef[2], digits = 2)
  }

  # Ensure no 1s and 0s are in cProb vectors, else these will throw errors when fed to qnorm
  cProbO_hi[cProbO_hi == 0] <- 1e-10
  cProbO_lo[cProbO_lo == 0] <- 1e-10
  cProbN_lo[cProbN_lo == 0] <- 1e-10
  cProbN_hi[cProbN_hi == 0] <- 1e-10

  cProbO_hi[cProbO_hi == 1] <- 1 - 1e-10
  cProbO_lo[cProbO_lo == 1] <- 1 - 1e-10
  cProbN_lo[cProbN_lo == 1] <- 1 - 1e-10
  cProbN_hi[cProbN_hi == 1] <- 1 - 1e-10


  est_c1  <- qnorm(1-cProbN_lo[1]) # c values are relative to new
  est_dc2 <- qnorm(1-cProbN_lo[2]) - qnorm(1-cProbN_lo[1]) # diff - c2 will = c1 + dc2. dc2 = c2 - c1
  est_dc3 <- qnorm(1-cProbN_lo[3]) - qnorm(1-cProbN_lo[2])
  est_dc4 <- qnorm(1-cProbN_lo[4]) - qnorm(1-cProbN_lo[3])
  est_dc5 <- qnorm(1-cProbN_lo[5]) - qnorm(1-cProbN_lo[4])


      tmpvals <- tibble(
        d_ho  = abs(rnorm(1,mean=est_muo_hi,sd=0.1)),
       sig_ho = abs(rnorm(1,mean=est_sigo_hi,sd=0.1)),
        d_lo   = abs(rnorm(1,mean=est_muo_lo,sd=0.1)),
        sig_lo  = abs(rnorm(1,mean=est_sigo_lo,sd=0.1)),
        d_hn = abs(rnorm(1,mean=est_mun_hi,sd=0.1)),
        sig_hn= abs(rnorm(1, mean=est_sign_hi,sd=0.1)),
       beta_ho = runif(1, min=0.1, max=3),
       beta_lo = runif(1,min=0.1,max=3),
       beta_hn = runif(1,min=0.1,max=3),
        c1     = rnorm(1,mean=est_c1,sd=0.1),
        dc1    = abs(rnorm(1,mean=est_dc2,sd=0.1)),
        dc2    = abs(rnorm(1,mean=est_dc3,sd=0.1)),
        dc3    = abs(rnorm(1,mean=est_dc4,sd=0.1)),
        dc4    = abs(rnorm(1,mean=est_dc5,sd=0.1))
      )


    if(model=="baselineUVSDT"){
      tmpvals <- tmpvals %>% select(d_ho,d_lo,d_hn,
                                    sig_ho,sig_lo,sig_hn,
                                    c1,dc1,dc2,dc3,dc4)
    } else if(model =="GumbelEVSDT"){

      tmpvals <- tmpvals %>% mutate(
        d_ho = -d_ho,
        d_lo = -d_lo,
        d_hn = -d_hn,
      c1 = -(c1+dc1+dc2+dc3+dc4)
      ) %>%
        select(d_ho,d_lo,d_hn,

              c1,dc1,dc2,dc3,dc4)

    } else if(model=="ExGaussNormEVSDT"){

      tmpvals <- tmpvals %>% select(d_ho,d_lo,d_hn,
                                    beta_ho,beta_lo,beta_hn,
                                    c1,dc1,dc2,dc3,dc4)
    }
      return(tmpvals)
}


optfunction <- function(model, data_list, par){

  if(model %in% c("baselineUVSDT")){
    out <- baseline_gaussian_uvsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  } else if (model %in% "GumbelEVSDT"){
    out <- gumbel_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")

  } else if (model %in% "ExGaussNormEVSDT"){
    out <- exGaussNorm_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
}

  return(-out)
}


predict_frequencies <- function(data_list,model,par){

  out <- vector()
  if(model %in% c("baselineUVSDT")){
    out <- baseline_gaussian_uvsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  }  else if (model %in% "GumbelEVSDT"){
    out <- gumbel_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")

  } else if (model %in% "ExGaussNormEVSDT"){
    out <- exGaussNorm_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  }


  return(out)

}

baseline_gaussian_uvsdt_opt <- function(data_list,par,predictorLL){


  d_ho     <- as.numeric(par[1])
  d_lo     <- as.numeric(par[2])
  d_hn     <- as.numeric(par[3])
  d_ln     <- 0
  sig_ho  <- as.numeric(par[4])
  sig_lo  <- as.numeric(par[5])
  sig_hn  <- as.numeric(par[6])
  sig_ln  <- 1
  c     <- cumsum(as.numeric(c(par[7:length(par)])))
 # sign  <- 1
  I <- c(-Inf,c,Inf) # put criteria into larger array

  # Likelihood of every trial
  # New items
  pHO_likJ <-  pLO_likJ <-  pHN_likJ <-  pLN_likJ <- vector()
  for (i in 1:(length(c)+1)){

    pHO_likJ[i] <- pnorm(I[i+1],mean=d_ho,sd=sig_ho)-pnorm(I[i],mean=d_ho,sd=sig_ho)
    pLO_likJ[i] <- pnorm(I[i+1],mean=d_lo,sd=sig_lo)-pnorm(I[i],mean=d_lo,sd=sig_lo)
    pHN_likJ[i] <- pnorm(I[i+1],mean=d_hn,sd=sig_hn)-pnorm(I[i],mean=d_hn,sd=sig_hn)
    pLN_likJ[i] <- pnorm(I[i+1],mean=d_ln,sd=sig_ln)-pnorm(I[i],mean=d_ln,sd=sig_ln)



  }



  if(predictorLL == "LL"){


    LN_likJ <-  data_list$low_ev_0 * log(pLN_likJ)
    HN_likJ <-  data_list$high_ev_0 * log(pHN_likJ)
    LO_likJ <-  data_list$low_ev_1 * log(pLO_likJ)
    HO_likJ <-  data_list$high_ev_1 * log(pHO_likJ)

    llk <- sum(c(LN_likJ,HN_likJ,LO_likJ,HO_likJ))

    if (is.na(llk)) llk <- -1e10
    if (llk == -Inf) llk <- -1e10

    return(llk)

  } else {
    return(c(pHO_likJ,pLO_likJ,pHN_likJ,pLN_likJ))
  }


}

gumbel_evsdt_opt <- function(data_list,par,predictorLL){

  d_ho     <- as.numeric(par[1])
  d_lo     <- as.numeric(par[2])
  d_hn     <- as.numeric(par[3])
  d_ln <- 0

  c     <- cumsum(as.numeric(c(par[4:length(par)])))

  # Constraints

  I <- c(-Inf,c,Inf) # put criteria into larger array


  # Likelihood of every trial
  # New items
  pHO_likJ <-  pLO_likJ <-  pHN_likJ <-  pLN_likJ <- vector()
  for (i in 1:(length(c)+1)){

    pHO_likJ[i] <- ordinal::pgumbel(I[i+1],location=d_ho,scale=1,max=F)-
      ordinal::pgumbel(I[i],location=d_ho,scale=1,max=F)
    pLO_likJ[i] <-  ordinal::pgumbel(I[i+1],location=d_lo,scale=1,max=F)-
      ordinal::pgumbel(I[i],location=d_lo,scale=1,max=F)
    pHN_likJ[i] <-  ordinal::pgumbel(I[i+1],location=d_hn,scale=1,max=F)-
      ordinal::pgumbel(I[i],location=d_hn,scale=1,max=F)
    pLN_likJ[i] <-  ordinal::pgumbel(I[i+1],location=d_ln,scale=1,max=F)-
      ordinal::pgumbel(I[i],location=d_ln,scale=1,max=F)



  }

  if(predictorLL == "LL"){


    LN_likJ <-  data_list$low_ev_0 * log(pLN_likJ)
    HN_likJ <-  data_list$high_ev_0 * log(pHN_likJ)
    LO_likJ <-  data_list$low_ev_1 * log(pLO_likJ)
    HO_likJ <-  data_list$high_ev_1 * log(pHO_likJ)

    llk <- sum(c(LN_likJ,HN_likJ,LO_likJ,HO_likJ))

    if (is.na(llk)) llk <- -1e10
    if (llk == -Inf) llk <- -1e10

    return(llk)

  } else {
    return(c(pHO_likJ,pLO_likJ,pHN_likJ,pLN_likJ))
  }

}
exGaussNorm_evsdt_opt <- function(data_list,par,predictorLL){

  d_ho     <- as.numeric(par[1])
  d_lo     <- as.numeric(par[2])
  d_hn     <- as.numeric(par[3])
  d_ln <- 0
 # sigo <- 1
  # lambda <- 1/par[2]
  # sigo  <- sqrt(sigma +  1/lambda^2)
  sig_ln<-1

  beta_ho <- as.numeric(par[4])
  beta_lo <- as.numeric(par[5])
  beta_hn <- as.numeric(par[6])
  c     <- cumsum(as.numeric(c(par[7:length(par)])))


  # Constraints for equal-variance!

  I <- c(-Inf,c,Inf) # put criteria into larger array

  # Likelihood of every trial
  # New items


  pHO_likJ <-  pLO_likJ <-  pHN_likJ <-  pLN_likJ <- vector()
  for (i in 1:(length(c)+1)){

    evi_ho <- brms::pexgaussian(I[i],mu=d_ho,sigma=1,beta=beta_ho)
    eviplus_ho <- brms::pexgaussian(I[i+1],mu=d_ho,sigma=1,beta=beta_ho)
    evi_lo <- brms::pexgaussian(I[i],mu=d_lo,sigma=1,beta=beta_lo)
    eviplus_lo <- brms::pexgaussian(I[i+1],mu=d_lo,sigma=1,beta=beta_lo)
    evi_hn <- brms::pexgaussian(I[i],mu=d_hn,sigma=1,beta=beta_hn)
    eviplus_hn <- brms::pexgaussian(I[i+1],mu=d_hn,sigma=1,beta=beta_hn)

    pHO_likJ[i] <- ifelse(is.na(eviplus_ho),0,eviplus_ho) -
      ifelse(is.na(evi_ho),0,evi_ho)
    pLO_likJ[i] <- ifelse(is.na(eviplus_lo),0,eviplus_lo) -
      ifelse(is.na(evi_lo),0,evi_lo)
    pHN_likJ[i] <- ifelse(is.na(eviplus_hn),0,eviplus_hn) -
      ifelse(is.na(evi_hn),0,evi_hn)

    pLN_likJ[i] <-  pnorm(I[i+1],mean=d_ln,sd=sig_ln)-
      pnorm(I[i],mean=d_ln,sd=sig_ln)

  }

  if(predictorLL == "LL"){


    LN_likJ <-  data_list$low_ev_0 * log(pLN_likJ)
    HN_likJ <-  data_list$high_ev_0 * log(pHN_likJ)
    LO_likJ <-  data_list$low_ev_1 * log(pLO_likJ)
    HO_likJ <-  data_list$high_ev_1 * log(pHO_likJ)

    llk <- sum(c(LN_likJ,HN_likJ,LO_likJ,HO_likJ))

    if (is.na(llk)) llk <- -1e10
    if (llk == -Inf) llk <- -1e10

    return(llk)

  } else {
    return(c(pHO_likJ,pLO_likJ,pHN_likJ,pLN_likJ))
  }



} #new: pnorm, old: exgaussian


fit_nlminb <- function(data, model, rep, startpar, freqdat, numcategories,idn){

  dp <- prep_data(data,freqdat,numcategories)
  out_list <- vector("list",rep)


  if (is.null(startpar)) {
    startpar1 <- NULL
    for (reps in c(1:rep)) {
      startpar1[[reps]] <- get_start_par(model,length(dp$confidence))
    }
    startpar2 <- NULL
    for (reps in c(1:rep)){
      startpar2[[reps]] <- get_start_par2(model,dp$datalist)
    }

    startpar <- bind_rows(startpar1,startpar2)
  }

  for (i in c(1:(rep*2))){

    start <- startpar %>% dplyr::slice(i) %>% unlist()
    tic <- Sys.time()


    tmp <- tryCatch(nlminb(start,objective = optfunction,
                           data_list = dp$datalist,
                           lower = get_par_limits(model,length(dp$confidence))[[1]],
                           upper = get_par_limits(model,length(dp$confidence))[[2]],
                           model = model,
                           control = list(eval.max = 300, iter.max = 300, trace = 1,
                                          rel.tol = 1e-5, x.tol = 1.5e-8)), error = function(e) NA)

    if (is.list(tmp)) {
      out_list[[i]] <- bind_cols(tidyr::spread(tibble::enframe(tmp$par),
                                               name, value) %>%
                                   setNames(sort(names(get_start_par(model,length(dp$confidence))))),
                                 tibble::tibble(npar = length(names(get_start_par(model,length(dp$confidence))))),
                                 tibble::as_tibble(tmp[c(2, 3, 4, 6)]), tibble::tibble(time = Sys.time() -
                                                                                         tic), tibble::tibble(model = model), tibble::tibble(id = dp$id),
                                 tibble::tibble(leftout = dp$leftout),
                                 tibble::tibble(rep = i),
                                 tibble::tibble(exp = dp$exp), tibble::tibble(cvid = dp$cvid),
                                 tidyr::spread(tibble::enframe(start, "start"),
                                               start, value) %>% setNames(paste0("s_", sort(names(get_start_par(model,length(dp$confidence)))))))
    }
  }
  dplyr::bind_rows(out_list)
}

FitSDT <- function (data, model, rep = rep,
                    startpar = NULL, freqdat = F,
                    numcategories=6,idn) {

  res <- fit_nlminb(data = data,
                    rep = rep,
                    model = model,
                    startpar = startpar,
                    freqdat = freqdat,
                    numcategories = numcategories,
                    idn=idn)

  return(res)
}

#
#
# GetGsquared <- function(model,obsdata,expdata){
#
#   obs <- c(obsdata$datalist$high_ev_1,
#            obsdata$datalist$low_ev_1,
#            obsdata$datalist$high_ev_0,
#            obsdata$datalist$low_ev_0)
#
#   f <- tibble( GSq = 2*sum(obs*log(obs/(expdata*60)))) # Calculate G^2
#   f$df     <- length(obs)-length(names(get_start_par(model,numcategories=6)))-1
#   f$pGSq   <- 1-pchisq(f$GSq,f$df)       # p value of GSq
#
#   return(f)
# }
