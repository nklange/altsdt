# General information --------------------------

# Functions to fit SDT models to data sets with 1 set of new
# items and 1 set of old items

library(tidyverse)

prep_data <- function(data,freqdat,numcategories){

  if(freqdat){

    prepdat <- data %>%
      mutate(Freq = ifelse(is.na(Freq) | Freq==0,1/numcategories,Freq)) %>%
      arrange(oldnew,response)

  } else {

    preprepdat <- data %>% group_by(oldnew,response) %>%
      summarize(Freq = length(response)) %>%
      ungroup() %>%
      mutate(oldnew = factor(oldnew,levels = c("New","Old")),
             response = factor(response, levels = c(1:numcategories)))

    fullresp <- preprepdat %>% tidyr::expand(oldnew,response)

    prepdat <- preprepdat %>% right_join(fullresp) %>%
      mutate(Freq = ifelse(is.na(Freq),1/numcategories,Freq)) %>%
      #mutate(Freq = ifelse(is.na(Freq),0,Freq)) %>%
      arrange(oldnew,response)

  }

  dl <- split(prepdat$Freq, f = prepdat$oldnew)
  trialsNew <- sum(prepdat %>% filter(oldnew == "New") %>% .$Freq)
  trialsOld <- sum(prepdat %>% filter(oldnew == "Old") %>% .$Freq)

  id <- unique(data$id)
  confidence <- sort(unique(prepdat$response))
  if (is.null(unique(data$cvid))) {
    cvid <- id
  }
  else {
    cvid <- unique(data$cvid)
  }
  if (is.null(unique(data$exp))) {
    exp <- NA
  }
  else {
    exp <- unique(data$exp)
  }
  if (is.null(unique(data$leftout))) {
    leftout <- 0
  }
  else {
    leftout <- unique(data$leftout)
  }


  return(list(datalist = dl, confidence = confidence, trialsNew = trialsNew, trialsOld = trialsOld, id = id,
              cvid = cvid, exp = exp, leftout = leftout))

}

get_start_par <- function(model,numcategories){


  #crit <- as_tibble(t(sort(runif(numcategories-1, min= -3, max=3)))) %>% set_names(paste0("c",c(1:(numcategories-1))))

  if (model %in% c("GaussianUVSDT")){
    dc <- as_tibble(t(runif(numcategories-2, min= .Machine$double.eps, max=0.3))) %>% set_names(paste0("dc",c(1:(numcategories-2))))
    pstartunif <- tibble(
      muo   = runif(1, min=0, max=3),
      sigo = runif(1, min=1, max=4),
      #crit
      c1   = runif(1, min=-3, max=0),
      dc
      )

    pstart <- bind_rows(pstartunif)

  }  else if (model %in% c("GumbelEVSDT")){
    dc <- as_tibble(t(runif(numcategories-2, min= .Machine$double.eps, max=0.3))) %>% set_names(paste0("dc",c(1:(numcategories-2))))
    pstartunif <- tibble(
      muo    = runif(1, min=-3, max=0),
      c1   = runif(1, min=-3, max=0),
      dc
    )

    pstart <- bind_rows(pstartunif)

  }   else if (model %in% c("ExGaussNormEVSDT")){
    dc <- as_tibble(t(runif(numcategories-2, min= .Machine$double.eps, max=0.3))) %>% set_names(paste0("dc",c(1:(numcategories-2))))
    pstartunif <- tibble(
      muo   = runif(1, min=0, max=3),
      betao = runif(1, min=0.1, max=3),
      c1   = runif(1, min=-3, max=0),
      dc
    )

    pstart <- bind_rows(pstartunif)

  }
  return(pstart)
}
get_par_limits <- function(model,numcategories){

  if (model %in% c("GaussianUVSDT")){
    # "d"    "sigo" "c1"   "dc1"  "dc2"  "dc3"  "dc4"

    lower <- c(-Inf, .Machine$double.eps,-Inf, rep(.Machine$double.eps,(numcategories-2)))
    #lower <- c(-Inf, .Machine$double.eps,rep(-Inf,(numcategories-1)))
    upper <- Inf

  } else if (model %in% c("GaussianUVSDT_equi")){
    # "d"    "sigo" "c1"   "dc1"  "dc2"  "dc3"  "dc4"

    lower <- c(-Inf, .Machine$double.eps,-Inf,.Machine$double.eps)
    #lower <- c(-Inf, .Machine$double.eps,rep(-Inf,(numcategories-1)))
    upper <- Inf

  } else if (model %in% c("GaussianEVSDT","GumbelEVSDT","GumbelFlipEVSDT",
                          "GumbelNormEVSDT","GumbelLargeEVSDT","GumbelLargeNormEVSDT")){
    # "d"  "c1"   "dc1"  "dc2"  "dc3"  "dc4"

    lower <- c(-Inf, -Inf, rep(.Machine$double.eps,(numcategories-1)))
    upper <- Inf

  } else if (model %in% c("ExGaussNormEVSDT","GumbelUVSDT")){
    # "d" "betao" "c1"   "dc1"  "dc2"  "dc3"  "dc4"

    lower <- c(-Inf,.Machine$double.eps,-Inf, rep(.Machine$double.eps,(numcategories-2)))
    upper <- Inf

  }else if (model %in% c("GaussianEVSDT2")){
    # "d" "betao" "c1"   "dc1"  "dc2"  "dc3"  "dc4"

    lower <- c(-Inf)
    upper <- Inf

  }


  list(lower,upper)
}


optfunction <- function(model, data_list, par){

  if(model %in% c("GaussianUVSDT")){
    out <- gaussian_uvsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  }  else if(model %in% c("GumbelEVSDT")){
    out <- gumbel_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  } else if(model %in% c("ExGaussNormEVSDT")){
    out <- exGaussNorm_evsdt_opt(data_list = data_list,par = par,predictorLL = "LL")
  }

  return(-out)
}


gaussian_uvsdt_opt <- function(data_list,par,predictorLL){


  d     <- as.numeric(par[1])
  sigo  <- as.numeric(par[2])
  c     <- cumsum(as.numeric(c(par[3:length(par)])))
  sign  <- 1
  I <- c(-Inf,c,Inf) # put criteria into larger array

  # Likelihood of every trial
  # New items
  pNlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pNlikJ[i] <- pnorm(I[i+1],mean=0,sd=sign)-pnorm(I[i],mean=0,sd=sign)
  }

  # Old items
  pOlikJ <- vector()
  for (i in 1:(length(c)+1)){

    pOlikJ[i] <- pnorm(I[i+1],mean=d,sd=sigo)-pnorm(I[i],mean=d,sd=sigo)
  }



  if(predictorLL == "LL"){

    DataN <- data_list$New
    DataO <- data_list$Old

    NlikJ <- DataN * log(pNlikJ)
    OlikJ <- DataO * log(pOlikJ)


    llk <- sum(c(NlikJ,OlikJ))

    if (is.na(llk)) llk <- -1e10
    if (llk == -Inf) llk <- -1e10

    return(llk)

  } else {
    return(c(pNlikJ,pOlikJ))
  }


}
gumbel_evsdt_opt <- function(data_list,par,predictorLL){

  d     <- as.numeric(par[1])
  sigo  <- 1
  c     <- cumsum(as.numeric(c(par[2:length(par)])))

  # Constraints
  sign  <- 1
  I <- c(-Inf,c,Inf) # put criteria into larger array

  # likelihood by response category

  pNlikJ <- vector()
  pOlikJ <- vector()

  for (i in 1:(length(c)+1)){

    pNlikJ[i] <- ordinal::pgumbel(I[i+1],location=0,scale=sign,max=F)-ordinal::pgumbel(I[i],location=0,scale=sign,max=F)
    pOlikJ[i] <- ordinal::pgumbel(I[i+1],location=d,scale=sigo,max=F)-ordinal::pgumbel(I[i],location=d,scale=sigo,max=F)

  }

  if(predictorLL == "LL"){

    Data <-  c(data_list$New,data_list$Old)
    llk <- Data * log(c(pNlikJ,pOlikJ))
    llk[Data == 0] <- 0

    llk <- sum(llk)
    if (is.na(llk)) llk <- -1e10
    if (llk == -Inf) llk <- -1e10

    return(llk)

  } else {
    return(c(pNlikJ,pOlikJ))
  }
}
exGaussNorm_evsdt_opt <- function(data_list,par,predictorLL){

  d     <- as.numeric(par[1])
  sigo <- 1
  # lambda <- 1/par[2]
  # sigo  <- sqrt(sigma +  1/lambda^2)

  betao <- as.numeric(par[2])
  c     <- cumsum(as.numeric(c(par[3:length(par)])))


  # Constraints for equal-variance!
  sign  <- sigo
  I <- c(-Inf,c,Inf) # put criteria into larger array

  # Likelihood of every trial
  # New items
  pNlikJ <- vector()
  for (i in 1:(length(c)+1)){
    pNlikJ[i] <- pnorm(I[i+1],0,sign)-pnorm(I[i],0,sign)

    # beta has to be > 0 (maybe try with super small beta -> super big lambda -> as gaussian
    # as possible)
    # evaluateboundaryi <- brms::pexgaussian(I[i],mu=0,sigma=sign,beta=0)
    # evaluateboundaryiplus <- brms::pexgaussian(I[i+1],mu=0,sigma=sign,beta=0)
    #
    # pNlikJ[i] <- ifelse(is.na(evaluateboundaryiplus),0,evaluateboundaryiplus) -
    #   ifelse(is.na(evaluateboundaryi),0,evaluateboundaryi)
  }




  # Old items
  pOlikJ <- vector()
  for (i in 1:(length(c)+1)){


    evaluateboundaryi <- brms::pexgaussian(I[i],mu=d,sigma=sigo,beta=betao)
    evaluateboundaryiplus <- brms::pexgaussian(I[i+1],mu=d,sigma=sigo,beta=betao)

    pOlikJ[i] <- ifelse(is.na(evaluateboundaryiplus),0,evaluateboundaryiplus) -
      ifelse(is.na(evaluateboundaryi),0,evaluateboundaryi)

    #pOlikJ[i] <- evaluateboundaryiplus - evaluateboundaryi

    # this is a hack for quick test. pexGaussian(-Inf,mu,sigma,beta) evaluates to NaN while
    # pnorm(-Inf,mead,sd) evaluates to 0, so does pgumbel(-Inf,location,scale,max=F)
    # pexGaussian(Inf,mu,sigma,beta) evaluates to 1


  }


  if(predictorLL == "LL"){

    DataN <- data_list$New
    DataO <- data_list$Old

    NlikJ <- DataN * log(pNlikJ)
    OlikJ <- DataO * log(pOlikJ)


    llk <- sum(c(NlikJ,OlikJ))

    if (is.na(llk)) llk <- -1e10
    if (llk == -Inf) llk <- -1e10

    return(llk)

  } else {
    return(c(pNlikJ,pOlikJ))
  }
} #new: pnorm, old: exgaussian



fit_nlminb <- function(data, model, rep, startpar, freqdat, numcategories){


  dp <- prep_data(data,freqdat,numcategories)
  out_list <- vector("list",rep)

  if (is.null(startpar)) {
    startpar <- NULL
    for (reps in c(1:rep)) {
      startpar[[reps]] <- get_start_par(model,length(dp$confidence))
    }
    startpar <- startpar %>% bind_rows()
  }

  for (i in seq_len(rep)){


    start <- startpar %>% dplyr::slice(i) %>% unlist()
    tic <- Sys.time()
    tmp <- tryCatch(nlminb(start,objective = optfunction,
                           data_list = dp$datalist,
                           lower = get_par_limits(model,length(dp$confidence))[[1]],
                           upper = get_par_limits(model,length(dp$confidence))[[2]],
                           model = model,
                           control = list(eval.max = 300, iter.max = 300, trace = 1,
                                          rel.tol = 1e-05, x.tol = 1.5e-08)), error = function(e) NA)
    if (is.list(tmp)) {
      out_list[[i]] <- bind_cols(tidyr::spread(tibble::enframe(tmp$par),
                                               name, value) %>% setNames(sort(names(get_start_par(model,
                                                                                                  length(dp$confidence))))),
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

FitSDT <- function (data, model, rep = rep, startpar = NULL, freqdat = F, numcategories=6) {

  res <- fit_nlminb(data = data,
                    rep = rep,
                    model = model,
                    startpar = startpar,
                    freqdat = freqdat,
                    numcategories = numcategories)

  return(res)
}

PredictSDT <- function(data = NULL, model, par, itemspertype = NULL,numcategories=6){



  if (is.null(data)){


    # using rmultinom
    # alternative: using sample(x,n,prob,replace=T) but unsampled categories are not listed in results

    probs <- predict_frequencies(data_list = NULL, model = model, par = as.numeric(par), numcategories=numcategories)

    predfreq <- c(rmultinom(c(1:numcategories),itemspertype[[1]],prob = probs[1:numcategories]),
                  rmultinom(c(1:numcategories),itemspertype[[2]],prob = probs[(numcategories+1):(numcategories * 2)]))

    simulated <- tibble(Freq = predfreq,
                        oldnew = rep(c("New","Old"),each = numcategories),
                        response = rep(c(1:numcategories),2))

    return(simulated)

  } else {

    dp <- prep_data(data,freqdat=F,numcategories)
    obsfreq<- c(dp$datalist$New,dp$datalist$Old)

    predfreq <- predict_frequencies(data_list = dp$datalist, model = model, par = par, numcategories = numcategories) *
      c(rep(dp$trialsNew,length(dp$confidence)),rep(dp$trialsOld,length(dp$confidence)))

    # Multiply each probability with the total number of new and old items
    # complicated here only to future-proof it for data sets with different numbers of targets/lures


    return(list(
      predicted = predfreq,
      observed = obsfreq
    ))

  }





}


predict_frequencies <- function(data_list,model,par,numcategories){

  out <- vector()
  if(model %in% c("GaussianUVSDT")){
    out <- gaussian_uvsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("GumbelEVSDT")){
    out <- gumbel_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("ExGaussNormEVSDT")){
    out <- exGaussNorm_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  } else if(model %in% c("GaussianEVSDT")){
    out <- gaussian_evsdt_opt(data_list = data_list,par = par,predictorLL = "predict")
  }

  return(out)

}

