# -----------------------------------------------------------------------------
# ------- RESOLVE Trial
# ------- Analysis of trial results for ITT and PP effects
# -----------------------------------------------------------------------------
#
# Last updated: July 30, 2025
# By: Joy Shi
#
# Description: Creating a function to analyze RESOLVE trial data under different 
# approaches.
# - Multiple imputation to handle missingness in covariates
# - ITT analyses (unadjusted, adjusted, and IP weighted for censoring)
# - PP analyses
#    - via IP weighting (unadjusted, adjustment for baseline covariates only, 
#      and adjustment for time-varying covariates)
#    - via g-formula (unadjusted, adjustment for baseline covariates only, 
#      and adjustment for time-varying covariates)
#
# Function is a more simple bootstrapping function to perform the PP analysis
# using IPW. If other analyses are of interest, then can update the
# bootstrapping function and generate separate functions for each of the 
# analyses of interest, or (ideally) wrap all of the analyses in one function
# (so that we don't have to duplicate the imputation each time)
#
# Please check "Note" for items to pay attention to
#
# -----------------------------------------------------------------------------

# ---- (0) Set-up ----

### Load libraries ###
library(tidyverse)
library(data.table)
library(mice)
library(pbapply)

### Import data ###
data <- read.csv("./data_PAIN.full.csv")
data <- data.table(data)

### Data checks ### 
### Note: reorder data here as needed for multiple imputation
data[id==4, edn:=4]


# ---- (1) Bootstrapping function ----
# Note: simple bootstrapping function to be customized for different analyses
# Here, using the time-varying PP analysis as an example

boot.ipw_pp <- function(b, outcome, standardize=T, truncate=0.999, seed=5438593){
  
  # Specifying seeds
  set.seed(seed)
  seed.list <- floor(runif(b+1)*10e8)
  
  # Looping over bootstrapped samples
  results <- pblapply(0:b, function(x){

    # First, run the analysis in the original data to get point estimate
    if (x==0){
      boot.data <- copy(data)
      boot.data[, newid:=id]
      boot.data[, id:=NULL]
    }
    # In subsequent loops, sample IDs with replacement to get new bootstrapped
    # sample. Then merge with original data to get the rest of the data
    if (x>0){
      set.seed(seed.list[x+1])
      boot.data <- data.table(id=sample(unique(data$id), replace=T))
      boot.data[, newid:=1:.N]
      boot.data <- merge(boot.data, data)
      boot.data[, id:=NULL]
    }

    # Conduct multiple imputation
    # Note: for analyses in the original data, get imputations for many copies
    # of the dataset (e.g., m=10 copies) and  average over all imputations.
    # This will provide a more reliable estimate (that isn't dependent on the
    # variability in the imputation). In the bootstrapped sampled, only get an
    # imputation for one copy (i.e., m=1) in order to appropriately reflect the
    # sample size of the original data

    if (x==0) m <- 20
    if (x>0)  m <- 1

    if (outcome=="rmdq") list.temp <- c("pain", "rmdq")
    if (outcome=="pain") list.temp <- c("rmdq", "pain")

    ## Reshape to wide ##
    boot.wide <- boot.data %>%
      filter(time<=4) %>%
      pivot_wider(
        id_cols     = c("newid", "age", "sex", "edn", "wrk.absnc", "allocn"),
        values_from = c("qolvas", "isi", "pain.bthr", "pain.prsnt",
                        "ceq", "pseq", "frbq", "tsk", "bbq", "pcs", "dass",
                        list.temp,
                        "adhred", "receiv"),  # Updated 20250815: reordered
        names_from  = "time",
        names_vary  = 'slowest'
      )

    ## Specifying predictor matrix for MICE ##
    pred.mat <- make.predictorMatrix(boot.wide)
    pred.mat["newid",] <- 0           # Do not impute ID or use as predictor
    pred.mat[,"newid"] <- 0
    pred.mat[2:19, 20:81] <- 0     # For covariates at time 0, do not use future info
    pred.mat[2:34, 35:81] <- 0     # For covariates at time 1, do not use future info
    pred.mat[2:49, 50:81] <- 0     # For covariates at time 2, do not use future info
    pred.mat[2:64, 65:81] <- 0     # For covariates at time 3, do not use future info
    pred.mat[65:81,] <- 0          # Do not impute variables at fourth time point
    pred.mat[,65:81] <- 0
    pred.mat[c(paste0("adhred_", seq(0,3))),] <- 0   # Do not impute adhred or use as predictor
    pred.mat[,c(paste0("adhred_", seq(0,3)))] <- 0

    ## Specifying method for imputation ##
    meth <- make.method(boot.wide)
    meth[grep("_4", names(meth))] <- "" # Do not impute at fourth time point
    # Note: for simplicity, using pmm to impute for each variable but
    # could use alternate methods if that is of interest

    ## MICE ##
    options(warn=-1)
    mice.results <- mice(
      data            = boot.wide,
      predictorMatrix = pred.mat,
      m               = m,
      method          = meth,
      seed            = seed.list[x+1],
      printFlag       = F,
    )
    options(warn=0)

    ## Obtaining an imputed dataset and changing back to long format ##
    boot.mice <- complete(mice.results, action="long") %>%
      mutate(bootid=seq(1:n())) %>%
      pivot_longer(
        cols          = qolvas_0:receiv_4,
        names_to      = c(".value", "time"),
        names_pattern = "([A-Za-z.]+)_([0-9])") %>%
      mutate(time=as.numeric(time)) %>%
      arrange(bootid, time) %>%
      data.table()

    # RUNNING ANALYSIS
    # Note: update below if interested in adapting this for a different
    # analysis

    # Calculating cumulative averages
    boot.mice[,`:=`(dass_cumavg       = ifelse(time<=3, cummean(dass), lag(cummean(dass))),
                    rmdq_cumavg       = ifelse(time<=3, cummean(rmdq), lag(cummean(rmdq))),
                    pcs_cumavg        = ifelse(time<=3, cummean(pcs), lag(cummean(pcs))),
                    bbq_cumavg        = ifelse(time<=3, cummean(bbq), lag(cummean(bbq))),
                    tsk_cumavg        = ifelse(time<=3, cummean(tsk), lag(cummean(tsk))),
                    frbq_cumavg       = ifelse(time<=3, cummean(frbq), lag(cummean(frbq))),
                    pseq_cumavg       = ifelse(time<=3, cummean(pseq), lag(cummean(pseq))),
                    ceq_cumavg        = ifelse(time<=3, cummean(ceq), lag(cummean(ceq))),
                    pain.prsnt_cumavg = ifelse(time<=3, cummean(pain.prsnt), lag(cummean(pain.prsnt))),
                    pain.bthr_cumavg  = ifelse(time<=3, cummean(pain.bthr), lag(cummean(pain.bthr))),
                    isi_cumavg        = ifelse(time<=3, cummean(isi), lag(cummean(isi))),
                    qolvas_cumavg     = ifelse(time<=3, cummean(qolvas), lag(cummean(qolvas))),
                    pain_cumavg       = ifelse(time<=3, cummean(pain), lag(cummean(pain)))), by=bootid]

    # Identify when each participant is censored due to treatment deviation
    boot.mice[, devtime:=min(fcase(
      allocn==0 & time>0 & receiv!=1, time,
      allocn==1 & time>0 & receiv!=2, time,
      default=99)), by=bootid]

    # Fitting treatment weights

    # Denominator of weights
    # Note: was having issues with convergence in some bootstrapped samples
    # (likely because there were few people that deviated)
    # As a result, removed squared terms and
    # modeled as a function of cumulative averages
    trt.denom0 <- withCallingHandlers(
      glm( # Fit stratified models by allocn
        receiv ~ age + sex + wrk.absnc + # as.factor(edn) +
          dass_cumavg + rmdq_cumavg + pcs_cumavg + bbq_cumavg + tsk_cumavg +
          frbq_cumavg + pseq_cumavg + ceq_cumavg + pain.prsnt_cumavg +
          pain.bthr_cumavg + isi_cumavg + qolvas_cumavg + pain_cumavg +
          as.factor(time),
        data=boot.mice[allocn==0 & time>0 & time<=3 & time<=devtime,],
        family=binomial()),
      warning = function(w){
        if (grepl("glm.fit: fitted probabilities numerically 0 or 1 occurred", conditionMessage(w))) {
          invokeRestart("muffleWarning") # Note: suppress this warning, not an issue
        }
        if (grepl("glm.fit: algorithm did not converge", conditionMessage(w))) {
          warning(paste0("glm.fit: algorithm did not converge for trt.denom0 in Iteration ", x), call.=F)
          invokeRestart("muffleWarning")
          # Note: if model did not converge, also print which model did not converge and in
          # which iteration
        }
      })


    trt.denom1 <- withCallingHandlers(
      glm( # Fit stratified models by allocn
        I(receiv/2) ~ age + sex + wrk.absnc + # as.factor(edn) +
          dass_cumavg + rmdq_cumavg + pcs_cumavg + bbq_cumavg + tsk_cumavg +
          frbq_cumavg + pseq_cumavg + ceq_cumavg + pain.prsnt_cumavg +
          pain.bthr_cumavg + isi_cumavg + qolvas_cumavg + pain_cumavg +
          as.factor(time),
        data=boot.mice[allocn==1 & time>0 & time<=3 & time<=devtime,],
        family=binomial()),
      warning = function(w){
        if (grepl("glm.fit: fitted probabilities numerically 0 or 1 occurred", conditionMessage(w))) {
          invokeRestart("muffleWarning") # Note: suppress this warning, not an issue
        }
        if (grepl("glm.fit: algorithm did not converge", conditionMessage(w))) {
          warning(paste0("glm.fit: algorithm did not converge for trt.denom1 in Iteration ", x), call.=F)
          invokeRestart("muffleWarning")
          # Note: if model did not converge, also print which model did not converge and in
          # which iteration
        }
      })

    if (standardize==F){ # If not standardizing, then the numerator of the
                         # weights will just be dependent on time

      # Numerator of weights
      trt.num0 <- glm(
        receiv ~ as.factor(time),
        data=boot.mice[allocn==0 & time>0 & time<=3 & time<=devtime,],
        family=binomial())

      trt.num1 <- glm(
        I(receiv/2) ~ as.factor(time),
        data=boot.mice[allocn==1 & time>0 & time<=3 & time<=devtime,],
        family=binomial())
    }
    if (standardize==T) { # If standardizing, then numerator of the weights
                          # will also include time-fixed covariates

      trt.num0 <- glm(
        receiv ~ as.factor(time) + age + sex + wrk.absnc,
        data=boot.mice[allocn==0 & time>0 & time<=3 & time<=devtime,],
        family=binomial())

      trt.num1 <- glm(
        I(receiv/2) ~ as.factor(time) + age + sex + wrk.absnc,
        data=boot.mice[allocn==1 & time>0 & time<=3 & time<=devtime,],
        family=binomial())

    }

    # Pulling predicted probabilities and estimating weights
    boot.mice[time>0 & time<=3, `:=`(
      trt.denom.pred0=predict(trt.denom0, newdata=boot.mice[time>0 & time<=3,], type="response"),
      trt.denom.pred1=predict(trt.denom1, newdata=boot.mice[time>0 & time<=3,], type="response"),
      trt.num.pred0=predict(trt.num0, newdata=boot.mice[time>0 & time<=3,], type="response"),
      trt.num.pred1=predict(trt.num1, newdata=boot.mice[time>0 & time<=3,], type="response"))]

    boot.mice[, sw_t:=fcase(
      time==0|time>3,                            1,   # Weights = 1 at baseline or after time 3
      time>0 & time<=3 & allocn==0 & receiv==1,  trt.num.pred0/trt.denom.pred0,
      time>0 & time<=3 & allocn==0 & receiv==0,  0,   # Weight = 0 if deviated
      time>0 & time<=3 & allocn==1 & receiv==2,  trt.num.pred1/trt.denom.pred1,
      time>0 & time<=3 & allocn==1 & receiv==0,  0    # Weight = 0 if deviated
    )]

    # Fitting censoring weights

    cens.model1 <- ~ allocn +  age + sex + wrk.absnc +
      dass_cumavg + pcs_cumavg + bbq_cumavg +
      qolvas_cumavg + pseq_cumavg + ceq_cumavg + isi_cumavg +
      rmdq_cumavg

    # Note: Because there are so few censoring events, this model
    # sometimes does not converge in bootstrapped samples. Removed:
    # - quadratic terms
    # - as.factor(edn)
    # - pain_cumavg
    # - tsk_cumavg
    # - pain.bthr_cumavg
    # - pain.prsnt_cumavg
    # - frbq_cumavg
    # Note: would typically want to model as a function of actual treatment
    # received (receiv), but here, missingness in pain is nearly perfectly
    # predicted by receiv

    # Specifying denominator model
    if (outcome=="pain") cens.model1 <- update(cens.model1, !is.na(pain) ~ .)
    if (outcome=="rmdq") cens.model1 <- update(cens.model1, !is.na(rmdq) ~ .)

    # Specifying numerator model
    if (standardize==F){ # If not standardizing, then numerator weights will
                         # only include treatment and time
      if (outcome=="pain"){ cens.model2 <- !is.na(pain) ~ allocn }
      if (outcome=="rmdq"){ cens.model2 <- !is.na(rmdq) ~ allocn }
    }
    if (standardize==T){ # If standardizing, then numerator weights will
                         # also include time-fixed covariates
      if (outcome=="pain"){ cens.model2 <- !is.na(pain) ~ allocn + age + sex + wrk.absnc }
      if (outcome=="rmdq"){ cens.model2 <- !is.na(rmdq) ~ allocn + age + sex + wrk.absnc }
    }


    # Fitting denominator of weights for pain
    cens.denom <- withCallingHandlers(
      glm(cens.model1,
          data=boot.mice[time==4,],
          family=binomial()),
      warning = function(w){
        if (grepl("glm.fit: fitted probabilities numerically 0 or 1 occurred", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
        if (grepl("glm.fit: algorithm did not converge", conditionMessage(w))) {
          warning(paste0("glm.fit: algorithm did not converge for cens.denom in Iteration ", x), call.=F)
          invokeRestart("muffleWarning")
          # Note: if model did not converge, also print which model did not converge and in
          # which iteration
        }
      })

    # Fitting numerator of weights for pain
    cens.num <- glm(cens.model2, data=boot.mice[time==4,], family=binomial())

    # Obtaining weights
    boot.mice[time==4, sw_c:=fcase(
      !is.na(boot.mice[time==4,][[outcome]]), predict(cens.num, newdata=boot.mice[time==4,], type="response")/
        predict(cens.denom, newdata=boot.mice[time==4,], type="response"),
      is.na(boot.mice[time==4,][[outcome]]),  0)]

    # Taking cumulative weights
    boot.mice[, sw:=cumprod(ifelse(time==4, sw_t*sw_c, sw_t)), by=bootid]

    # Truncate weights at specified truncation level
    boot.mice[time==4, sw:=pmin(sw, quantile(boot.mice[time==4,]$sw, truncate))]

    # Fit outcome model
    if (standardize==F){ # If not standardizing, then outcome model will only include
      # treatment assignment
      if (outcome=="pain"){ out.model <- pain ~ allocn }
      if (outcome=="rmdq"){ out.model <- rmdq ~ allocn }
    }
    if (standardize==T){ # If standardizing, outcome model will also include
      # time-fixed covariates
      if (outcome=="pain"){ out.model <- pain ~ allocn + age + sex + wrk.absnc }
      if (outcome=="rmdq"){ out.model <- rmdq ~ allocn + age + sex + wrk.absnc }
    }

    outcome.model <- glm(out.model, data=boot.mice[time==4 & time<devtime & sw!=0,], weights=sw)

    if (standardize==F){ # If not standardizing, then can use estimate from model
      estimate.temp <- unname(outcome.model$coefficients[2])
    }

    if (standardize==T){ # Otherwise, need to standardize
      boot.standardize <- rbind(
        copy(boot.mice[time==4,])[, allocn:=0], # Create copy where everyone's allocn set to 0
        copy(boot.mice[time==4,])[, allocn:=1]  # Create copy where everyone's allocn set to 1
      )
      boot.standardize[, pain.predict:=predict(outcome.model, newdata=boot.standardize)] # Obtain predicted values of outcome
      estimate.temp <- mean(boot.standardize[allocn==1,]$pain.predict) - mean(boot.standardize[allocn==0,]$pain.predict) # Obtain standardized estimate

    }

    return(c(Iteration = x,
             Estimate = estimate.temp))

  })  
  
  boot.results <- data.table(do.call(rbind, results))
  
  return(list(
    c(Estimate = boot.results[Iteration==0,]$Estimate,
      CI_LL = unname(quantile(boot.results[Iteration!=0,]$Estimate, p=0.025)),
      CI_UL = unname(quantile(boot.results[Iteration!=0,]$Estimate, p=0.975))),
    boot.results))
  
}

boot.ipw_pp(b=500, standardize=T, truncate=0.999, outcome="pain")
boot.ipw_pp(b=500, standardize=T, truncate=0.995, outcome="rmdq") 
