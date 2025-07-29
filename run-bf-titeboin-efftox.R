rm(list = ls())
library(dplyr)
library(BOIN)
library(readr)
library(boinet)
library(bfboinet)
library(filelock)
library(tibble)
# library(TITEgBOIN)

### Functions =======================
patient.enrol <- function(n, rate, unit = "month", type = "poisson") {
  type <- match.arg(type)
  if (type == "poisson") {
    inter_arrival <- rexp(n, rate = rate)
    cum_arrival <- cumsum(inter_arrival)
  }
  return(cum_arrival)
}

time.to.event <- function(n, pi, kk = 0.5, window = 1, unit = "month", type = "weibull") {
  # kk < 1: first kk month with kk probability of DLT
  if (type == "weibull") {
    ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
    alpha <- log(log(1 - pi) / log(1 - (kk * pi))) / log(2)
    lambda <- -log(1 - pi) / (window^alpha)
    time <- (-log(runif(n)) / lambda)^(1 / alpha)
  }

  return(time)
}

next.dose <- function(target, d, y, y.obs, n, n.pend, elimi, lambda.e, lambda.d) {
  # decide next dose
  conflict <- FALSE
  suspend <- FALSE
  escalate <- FALSE
  ndose <- length(n)

  phat <- y / n

  if (elimi[d] == 1) {
    next_d <- which(elimi == 1)[1] - 1
    return(list(next_d = next_d, suspend = suspend, escalate = escalate))
  }

  conflict_cond1 <- phat[seq_len(d - 1)] >= lambda.d & (y.obs[seq_len(d - 1)] / n[seq_len(d - 1)]) >= target
  conflict_cond2 <- phat[seq_len(d - 1)] <= lambda.e & (y.obs[seq_len(d - 1)] / n[seq_len(d - 1)]) <= target
  if (any(conflict_cond1)) {
    conflict_set <- seq_len(d - 1)[conflict_cond1]
    conflict <- TRUE
  } else if (phat[d] <= lambda.e & (y.obs[d] / n[d]) <= target) {
    # current dose: escalate
    if (any((!conflict_cond1) & (!conflict_cond2))) {
      # bf doses: stay
      conflict_set <- seq_len(d - 1)[(!conflict_cond1) & (!conflict_cond2)]
      conflict <- TRUE
    } else {
      # bf doses: escalate
      # no conflict
      escalate <- TRUE
      if (d != ndose & elimi[d + 1] == 0) {
        next_d <- d + 1
      } else {
        next_d <- d
      }
    }
  } else {
    # no conflict
    if (phat[d] <= lambda.e & (y.obs[d] / n[d]) <= target) {
      escalate <- TRUE
      if (d != ndose & elimi[d + 1] == 0) {
        next_d <- d + 1
      } else {
        next_d <- d
      }
    } else if (phat[d] >= lambda.d & (y.obs[d] / n[d]) >= target) {
      if (d != 1) {
        next_d <- d - 1
      } else {
        next_d <- d
        if (n.pend > 0) suspend <- TRUE
      }
    } else {
      next_d <- d
    }
  }

  if (conflict) {
    conflict_d <- max(conflict_set)
    y_pool <- sum(y[conflict_d:d])
    yobs_pool <- sum(y.obs[conflict_d:d])
    n_pool <- sum(n[conflict_d:d])

    if ((y_pool / n_pool) <= lambda.e & (yobs_pool / n_pool) <= target) {
      escalate <- TRUE
      if (d != ndose & elimi[d + 1] == 0) {
        next_d <- d + 1
      } else {
        next_d <- d
      }
    } else if ((y_pool / n_pool) >= lambda.d & (yobs_pool / n_pool) >= target) {
      if (d != 1) {
        next_d <- d - 1
      } else {
        next_d <- d
        if (n.pend > 0) suspend <- TRUE
      }
    } else {
      next_d <- d
    }
  }

  return(list(next_d = next_d, suspend = suspend, escalate = escalate))
}


### BF-TITE-BOIN =======================
run.bfboin <- function(target, 
                       pT.true, 
                       pE.true, 
                       w0 = 2 / 3, 
                       ncohort, 
                       cohortsize,
                       fixed.cohort = TRUE, 
                       accuralrate = 3,
                       DLTwindow = 1, 
                       effwindow = 3,
                       n.earlystop = 100, 
                       n.backfilling = 100, 
                       bf.type = c("tox", "utility"), 
                       bf.extended = FALSE, 
                       utility.method = c("mean", "pp"),
                       startdose = 1, 
                       init.size = 3, 
                       titration = FALSE,
                       pT.saf = 0.6 * target, 
                       pT.tox = 1.4 * target, 
                       pE.low = 0.25, 
                       cutoff.eli = 0.95, 
                       cutoff.eff = 0.95,
                       extrasafe = FALSE, 
                       offset = 0.05,
                       boundMTD = FALSE, 
                       prior.p = rep(1 / 3, 3),
                       min.complete = c(tox = 0.51, eff = 0.51),
                       min.MF = c(tox = 0.25, eff = 0),
                       ntrial = 1000, 
                       seed = 6) {
  
  bf.type <- match.arg(bf.type)
  utility.method <- match.arg(utility.method)
  
  set.seed(seed)
  # utility: (no tox, eff), (no tox, no eff), (tox, eff), (tox, no eff)
  utility_score <- c(100, 100 * w0 / (w0 + 1), 100 / (w0 + 1), 0)
  lower_u <- utility_score[1] * (1 - pT.tox) * pE.low + utility_score[2] * (1 - pT.tox) * (1 - pE.low) + utility_score[3] * pT.tox * pE.low
  ub <- lower_u + (100 - lower_u)/2
  
  
  if (cohortsize == 1) titration <- FALSE
  lambda_e <- log((1 - pT.saf) / (1 - target)) / log(target * (1 - pT.saf) / (pT.saf * (1 - target)))
  lambda_d <- log((1 - target) / (1 - pT.tox)) / log(pT.tox * (1 - target) / (target * (1 - pT.tox)))

  ndose <- length(pT.true)
  npts <- ncohort * cohortsize
  nprior <- 1
  priortox <- target / 2 # the prior mean for the Beta prior
  prioreff <- 0.50 / 2

  Ytox_s1 <- Yeff_s1 <- Ytox_s2 <- Yeff_s2 <- N_s1 <- N_s2 <- matrix(rep(0, ndose * ntrial), ncol = ndose)
  dselect <- matrix(NA, nrow = 2,ncol = ntrial) %>% {`rownames<-`(., c("MTD", "OBD"))}
  duration <- matrix(NA, nrow = 2,ncol = ntrial) %>% {`rownames<-`(., c("dose_finding", "extended_bf"))}

  for (trial in 1:ntrial) {
    # generate patient arrival
    arrival_t <- arrival_t0 <- patient.enrol(n = (10 * npts), rate = accuralrate)
    # each patient got a total of `ndose` possible DLT times (ndose x totaln matrix)
    DLT_t <- sapply(1:ndose, function(rr) time.to.event(n = length(arrival_t), pi = pT.true[rr], window = DLTwindow)) %>% t()
    eff_t <- sapply(1:ndose, function(rr) time.to.event(n = length(arrival_t), pi = pE.true[rr], window = effwindow)) %>% t()


    # initialize
    n <- elimi <- rep(0, ndose)
    earlystop <- 0 # dose 1 stopping indicator
    d <- startdose
    this_cosize <- init_cosize <- 1
    bf_dose <- dv <- c()
    DLTresult_t <- effresult_t <- assess_t <- decision_t <- c() # DLTresult_t: time of obtaining the DLT results (y/n)
    trial_log <- list(
      dose_finding = c(), backfilling = c(),
      y_tox = c(), y_eff = c(), n = c()
    )
    patient_type <- c(1) # 1 for dose-finding; 2 for backfiling
    ii <- 1
    trialcont <- TRUE
    early_susp <- suspend <- FALSE
    # ft <- TRUE
    # if (titration) {
    #   z <- (runif(ndose) < pT.true)
    #   if (sum(z) == 0) {
    #     d <- ndose
    #     n[1:ndose] <- 1
    #   } else {
    #     d <- which(z == 1)[1]
    #     n[1:d] <- 1
    #     y_tox[d] <- 1
    #   }
    # }

    while (trialcont) {
      if (!suspend) {
        # update trial timeline
        if (ii == 1) {
          assess_t[ii] <- arrival_t[ii]
        }

        # eff & tox data collection
        if (patient_type[ii] == 2) {
          current_d <- bf_d
          trial_log$backfilling[ii] <- dv[ii] <- bf_d
        } else {
          current_d <- d
          trial_log$dose_finding[ii] <- dv[ii] <- d
        }
        DLTresult_t[ii] <- assess_t[ii] + min(DLT_t[current_d, ii], DLTwindow)
        effresult_t[ii] <- assess_t[ii] + min(eff_t[current_d, ii], effwindow)

        y_tox0 <- replace(numeric(ndose), current_d, (1 * (DLT_t[current_d, ii] <= DLTwindow)))
        y_eff0 <- replace(numeric(ndose), current_d, (1 * (eff_t[current_d, ii] <= effwindow)))
        n0 <- replace(numeric(ndose), current_d, 1)

        if (ii == 1) {
          trial_log$y_tox <- matrix(y_tox0, ncol = 1)
          trial_log$y_eff <- matrix(y_eff0, ncol = 1)
          trial_log$n <- matrix(n0, ncol = 1)
        } else {
          trial_log$y_tox <- cbind(trial_log$y_tox, y_tox0)
          trial_log$y_eff <- cbind(trial_log$y_eff, y_eff0)
          trial_log$n <- cbind(trial_log$n, n0)
        }
      }

      n <- rowSums(trial_log$n)

      ## add timeline for pending pts
      if (sum(patient_type == 1) >= npts) {
        # trial stops: max patients reached
        trialcont <- FALSE
        decision_t[ii] <- assess_t[ii] + max(DLTwindow, effwindow)
      } else if (patient_type[ii] == 1 & this_cosize == cohortsize & n[d] >= n.earlystop) {
        # Stop trial if number of patients assigned to single dose reaches "n.earlystop" and the decision is to stay
        # If early stop number reached, let the last patient complete the follow-up
        decision_DFt <- assess_t[ii] + DLTwindow
        decision_t[ii] <- arrival_t[ii + 1] # for tite calculation
        early_susp <- TRUE
      } else {
        decision_t[ii] <- arrival_t[ii + 1]
      }

      finish_idx <- list(
        tox = (DLTresult_t <= decision_t[ii]), eff = (effresult_t <= decision_t[ii])
      )
      followup_t <- list(tox = numeric(length = ii), eff = numeric(length = ii))
      followup_t$tox[finish_idx$tox] <- DLTresult_t[finish_idx$tox] - assess_t[finish_idx$tox] # until the most recent decision time point
      followup_t$tox[!finish_idx$tox] <- decision_t[ii] - assess_t[!finish_idx$tox]

      followup_t$eff[finish_idx$eff] <- effresult_t[finish_idx$eff] - assess_t[finish_idx$eff]
      followup_t$eff[!finish_idx$eff] <- decision_t[ii] - assess_t[!finish_idx$eff]


      if (sum(finish_idx$tox) == 1) {
        ytox_obs <- trial_log$y_tox[, finish_idx$tox]
      } else {
        ytox_obs <- rowSums(trial_log$y_tox[, finish_idx$tox])
      }
      if (sum(finish_idx$eff) == 1) {
        yeff_obs <- trial_log$y_eff[, finish_idx$eff]
      } else {
        yeff_obs <- rowSums(trial_log$y_eff[, finish_idx$eff])
      }


      # approx & exact estimation
      totalt_v <- list(tox = c(), eff = c())
      phat_tox <- phat_eff <- matrix(NA, nrow = 2, ncol = ndose)
      phat_tox_exact <- phat_eff_exact <- matrix(NA, nrow = 2, ncol = ndose)
      for (dd in 1:ndose) {
        cset <- (dv == dd)
        # tox
        totalt <- followup_t$tox[(cset & (!finish_idx$tox))]
        totalt <- 3 * prior.p[1] * totalt * (totalt <= DLTwindow / 3) +
          ((prior.p[1] - prior.p[2]) * DLTwindow + 3 * prior.p[2] * totalt) * (DLTwindow / 3 < totalt & totalt <= 2 * DLTwindow / 3) +
          ((prior.p[1] + prior.p[2] - 2 * prior.p[3]) * DLTwindow + 3 * prior.p[3] * totalt) * (2 * DLTwindow / 3 < totalt & totalt <= DLTwindow)
        totalt_v$tox[dd] <- sum(totalt) / DLTwindow

        n_pending <- sum(cset & (!finish_idx$tox))
        if (n_pending == 0) {
          phat_tox[, dd] <- phat_tox_exact[, dd] <- c(ytox_obs[dd], n[dd])
        } else {
          # compute the approx estimated toxicity rate based on the imputed data
          phat0 <- (ytox_obs[dd] + priortox * nprior) / (n[dd] - n_pending + nprior)
          phat_tox[, dd] <- c(
            (ytox_obs[dd] + phat0 / (1 - phat0) * (n_pending - totalt_v$tox[dd])), n[dd]
            )
          # exact estimate
          E_yi <- (phat0 * (1 - totalt / DLTwindow)) / (phat0 * (1 - totalt / DLTwindow) + 1 - phat0)
          phat_tox_exact[, dd] <- c((ytox_obs[dd] + sum(E_yi)), n[dd])
        }

        # eff
        totalt <- followup_t$eff[(cset & (!finish_idx$eff))]
        totalt <- 3 * prior.p[1] * totalt * (totalt <= effwindow / 3) +
          ((prior.p[1] - prior.p[2]) * effwindow + 3 * prior.p[2] * totalt) * (effwindow / 3 < totalt & totalt <= 2 * effwindow / 3) +
          ((prior.p[1] + prior.p[2] - 2 * prior.p[3]) * effwindow + 3 * prior.p[3] * totalt) * (2 * effwindow / 3 < totalt & totalt <= effwindow)
        totalt_v$eff[dd] <- sum(totalt) / effwindow

        n_pending <- sum(cset & (!finish_idx$eff))
        if (n_pending == 0) {
          phat_eff[, dd] <- phat_eff_exact[, dd] <- c(yeff_obs[dd], n[dd])
        } else {
          # compute the estimated efficacy rate based on the imputed data
          phat0 <- (yeff_obs[dd] + prioreff * nprior) / (n[dd] - n_pending + nprior)
          phat_eff[, dd] <- c(
            (yeff_obs[dd] + phat0 / (1 - phat0) * (n_pending - totalt_v$eff[dd])), n[dd]
            )
          # exact estimate
          E_yi <- (phat0 * (1 - totalt / effwindow)) / (phat0 * (1 - totalt / effwindow) + 1 - phat0)
          phat_eff_exact[, dd] <- c((yeff_obs[dd] + sum(E_yi)), n[dd])
        }
      }

      # for n > 0, calculate utility
      uhat <- vector(mode = "numeric", length = ndose)
      if (utility.method == "mean") {
        for (dd in which(n > 0)) {
          uhat[dd] <- phat_eff_exact[1, dd] / phat_eff_exact[2, dd] - w0 * (phat_tox_exact[1, dd] / phat_tox_exact[2, dd])
        }
      } else if (utility.method == "pp") {
        quasi_x <- vector(mode = "numeric", length = ndose)
        for (dd in which(n > 0)) {
          quasi_x[dd] <-  (utility_score[3] * phat_eff_exact[1, dd] + utility_score[2] * (n[dd] - phat_tox_exact[1, dd]))/100
          uhat[dd] <- 1 - pbeta(ub/100, (quasi_x[dd] + 1), (n[dd] - quasi_x[dd] + 1))
        }
      }
      

      if (this_cosize == cohortsize) {
        for (dd in 1:ndose) {
          if (1 - pbeta(target, ytox_obs[dd] + 1, n[dd] - ytox_obs[dd] + 1) > cutoff.eli & n[dd] >= 3) {
            elimi[dd:ndose] <- 1
            break
          }
        }
        if (elimi[1] == 1) {
          earlystop <- 1
          trialcont <- FALSE
        }
      }

      # backfilling: candidate doses `bf_dose` (open & close a dose)
      # admissible set:
      if (d > 1) {
        # safety
        Aset <- seq_len(d - 1)
        
        # efficacy
        if (bf.type == "utility") {
          prob_eff <- pbeta(pE.low, (phat_eff_exact[1, Aset] + prioreff * nprior), (n[Aset] - phat_eff_exact[1, Aset] + (1 - prioreff) * nprior))
          Aset <- Aset[n[Aset] <= 3 | prob_eff < cutoff.eff]
        } else if (bf.type == "tox") {
          Aset <- Aset[cumsum(yeff_obs[seq_len(d - 1)]) >= 1]
        }
       
      } else {
        Aset <- c()
      }
      bf_dose <- Aset
     

      # close
      if (length(bf_dose) > 0) {
        dlt1 <- ((phat_tox[1, bf_dose] / phat_tox[2, bf_dose]) > lambda_d)
        dlt2 <- sapply(bf_dose, function(rr) {
          ((phat_tox[1, rr] + phat_tox[1, (rr + 1)]) / (phat_tox[2, rr] + phat_tox[2, (rr + 1)])) > lambda_d
        })

        bf_dose <- bf_dose[!(dlt1 & dlt2)]
        bf_dose <- bf_dose[n[bf_dose] < n.backfilling] # close if reaches n.backfilling
        bf_dose <- bf_dose[elimi[bf_dose] == 0] # only include doses not eliminated
      }

      suspend <- FALSE
      # during n.earlystop triggered suspension period
      if (trialcont & early_susp) {
        if (arrival_t[ii + 1] < decision_DFt) {
          if (length(bf_dose) > 0) {
            # backfillling
            bf_d <- switch(
              bf.type,
              utility = bf_dose[which.max(uhat[bf_dose])],
              tox     = max(bf_dose)
            )
            assess_t[ii + 1] <- arrival_t[ii + 1]
            patient_type[ii + 1] <- 2
            ii <- ii + 1
            next
          } else {
            # suspend
            arrival_t <- arrival_t[-(ii + 1)]
            suspend <- TRUE
            next
          }
        } else {
          # out of the suspend period
          early_susp <- FALSE
          # arrival_t[arrival_t >= decision_DFt] <- arrival_t[arrival_t >= decision_DFt] - (min(arrival_t[arrival_t >= decision_DFt]) - decision_DFt) ## ??? (to match TITE-BOIN results)
          
          # dose-finding & check whether to stop the trial using the observed tox
          ph <- phat_tox[1, d] / n[d]
          for (dd in 1:ndose) {
            ntox.curr1 <- rowSums(trial_log$y_tox)[dd]
            if (1 - pbeta(target, ntox.curr1 + 1, n[dd] - ntox.curr1 + 1) > cutoff.eli & n[dd] >= 3) {
              elimi[dd:ndose] <- 1
              break
            }
          }

          if (!(ph <= lambda_e & (ytox_obs[d] / n[d]) <= target) & !(ph >= lambda_d & (ytox_obs[d] / n[d]) >= target)) {
            # if stay then stop
            trialcont <- FALSE
          } else if (ph <= lambda_e & (ytox_obs[d] / n[d]) <= target) {
            # if escalate
            if ((d == ndose | elimi[min(d + 1, ndose)] == 1)) {
              # but next dose eliminated or current dose is highest dose
              trialcont <- FALSE
            }
          } else if (ph >= lambda_d & (ytox_obs[d] / n[d]) >= target & d == 1) {
            # if de-escalate
            trialcont <- FALSE
          }
        }
      }
      # end of earlystop


      # allocate next patient (if trial continues)
      if (trialcont) {
        # this_cosize: size of the current dose
        if (ii < init.size) {
          # initital cohort (size = 3): to dose-finding cohort
          init_cosize <- init_cosize + 1
          assess_t[ii + 1] <- arrival_t[ii + 1]
          d <- d
          patient_type[ii + 1] <- 1

        } else {
          if (ii == init.size) this_cosize <- min(init_cosize, cohortsize)
          if (fixed.cohort) {
            # fixed cohort size (e.g., = 3)
            if (this_cosize == cohortsize) {
              # Suspending rules
              #----------------------------------------------------#
              # Rule 1: Suspend the accrual to await more data if < min.complete % patients
              # have complete DLT window at the current dose d (default min.complete=51).

              # (Rule 2) Delay escalation decision if the last pending patient at the current
              # dose have completed < min.MF % of the DLT assessment window (default min.MF=25).
              #----------------------------------------------------#

              cset_susp <- (dv == d)
              npend_susp <- list(
                tox = sum(cset_susp & (!finish_idx$tox)), eff = sum(cset_susp & (!finish_idx$eff))
              )

              # obs y's trigger de-escalate (tox only)
              if ((ytox_obs[d] / n[d]) >= lambda_d) {
                # no need to suspend; directly go to the dose-finding cohort
                if (d == 1) {
                  d <- d
                  if (npend_susp$tox > 0) {
                    arrival_t <- arrival_t[-(ii + 1)]
                    suspend <- TRUE
                  } else {
                    this_cosize <- 1
                    assess_t[ii + 1] <- arrival_t[ii + 1]
                    patient_type[ii + 1] <- 1
                  }
                } else {
                  d <- d - 1
                  this_cosize <- 1
                  assess_t[ii + 1] <- arrival_t[ii + 1]
                  patient_type[ii + 1] <- 1
                }
                if (!suspend) ii <- ii + 1
                next
              }

              susp_rule1 <- (((n[d] - npend_susp$tox) < (min.complete["tox"] * n[d])) | 
                               ((n[d] - npend_susp$eff) < (min.complete["eff"] * n[d]))) %>% {`names<-`(., NULL)}

              decide_tmp <- next.dose(
                target = target, d = d, y = phat_tox[1, ], y.obs = ytox_obs, n = phat_tox[2, ],
                n.pend = npend_susp$tox, elimi = elimi, lambda.e = lambda_e, lambda.d = lambda_d
              )
              d_tmp <- decide_tmp$next_d
              susp_rule3 <- decide_tmp$suspend

              SFTL_tox <- ifelse(
                npend_susp$tox == 0,
                1,
                min(followup_t$tox[cset_susp & (!finish_idx$tox)]) / DLTwindow
              )
              SFTL_eff <- ifelse(
                npend_susp$eff == 0,
                1,
                min(followup_t$eff[cset_susp & (!finish_idx$eff)]) / effwindow
              )

              susp_rule2 <- (decide_tmp$escalate) & (
                (SFTL_tox < min.MF["tox"] | npend_susp$tox == n[d]) | (SFTL_eff < min.MF["eff"])
              ) %>% {`names<-`(., NULL)}


              if (length(bf_dose) > 0) {
                if (susp_rule1 | susp_rule2 | susp_rule3) {
                  # backfilling
                  bf_d <- switch(
                    bf.type,
                    utility = bf_dose[which.max(uhat[bf_dose])],
                    tox     = max(bf_dose)
                  )
                  assess_t[ii + 1] <- arrival_t[ii + 1]
                  patient_type[ii + 1] <- 2
                } else {
                  # dose-finding
                  d <- d_tmp
                  this_cosize <- 1
                  assess_t[ii + 1] <- arrival_t[ii + 1]
                  patient_type[ii + 1] <- 1
                }
              } else {
                if (susp_rule1 | susp_rule2 | susp_rule3) {
                  arrival_t <- arrival_t[-(ii + 1)]
                  suspend <- TRUE
                } else {
                  # start a new cohort
                  d <- d_tmp
                  this_cosize <- 1
                  assess_t[ii + 1] <- arrival_t[ii + 1]
                  patient_type[ii + 1] <- 1
                }
              }
            } else {
              # to dose-finding part & current dose cohort
              d <- d
              this_cosize <- this_cosize + 1
              assess_t[ii + 1] <- arrival_t[ii + 1]
              patient_type[ii + 1] <- 1
            }
          }
        }

        # if continue
        if (!suspend) ii <- ii + 1
      }
    }


    Ytox_s1[trial, ] <- rowSums(trial_log$y_tox)
    Yeff_s1[trial, ] <- rowSums(trial_log$y_eff)
    N_s1[trial, ] <- n
    duration["dose_finding", trial] <- (max(assess_t) - arrival_t[1]) + max(DLTwindow, effwindow)
    
    if (earlystop == 1) {
      dselect[ , trial] <- 99  # early stopping: no need to extend
      # OBD: 99 means not entering extended backfilling due to early stop or not planned
      duration["extended_bf", trial] <- duration["dose_finding", trial]
      N_s2[trial, ] <- 0
    } else {
      d_mtd <- select.mtd(
        target = target, npts = n, ntox = rowSums(trial_log$y_tox), cutoff.eli = cutoff.eli,
        extrasafe = extrasafe, offset = offset, boundMTD = boundMTD, p.tox = pT.tox
      )$MTD
      dselect["MTD", trial] <- d_mtd
      
      ### Extended backfilling
      Aset_ext <- seq_len(d_mtd)
      prob_eff_ext <- pbeta(pE.low, (Yeff_s1[trial, Aset_ext] + prioreff * nprior), (n[Aset_ext] - Yeff_s1[trial, Aset_ext] + (1 - prioreff) * nprior))
      Aset_ext <- Aset_ext[n[Aset_ext] <= 3 | prob_eff_ext < cutoff.eff]
      # top 2 utility
      if (length(Aset_ext) > 2) {
        Aset_ext <- Aset_ext[order(uhat[Aset_ext], decreasing = TRUE)][1:2] 
      }
      if (length(Aset_ext) == 0) {
        # end of trial & no OBD
        n_ext <- rep(0, ndose)
        duration["extended_bf", trial] <- duration["dose_finding", trial]
        dselect["OBD", trial] <- 99
        N_s2[trial, ] <- 0
        next
      } else {
        n_ext <- replace(numeric(ndose), Aset_ext, pmax(n.backfilling - n[Aset_ext], 0))
      }
      
      uhat_final <- vector(mode = "numeric", length = ndose)
      if (!bf.extended | all(n_ext == 0)) {
        # end of trial & select OBD
        duration["extended_bf", trial] <- duration["dose_finding", trial]
        if (utility.method == "mean") {
          for (dd in Aset_ext) {
            uhat_final[dd] <- (rowSums(trial_log$y_eff)[dd] / n[dd] - w0 * (rowSums(trial_log$y_tox)[dd] / n[dd])) %>% round(., 6)
          }
        } else if (utility.method == "pp") {
          quasi_x <- vector(mode = "numeric", length = ndose)
          for (dd in Aset_ext) {
            quasi_x[dd] <-  (utility_score[3] * rowSums(trial_log$y_eff)[dd] + utility_score[2] * (n[dd] - rowSums(trial_log$y_tox)[dd]))/100
            uhat_final[dd] <- (1 - pbeta(ub/100, (quasi_x[dd] + 1), (n[dd] - quasi_x[dd] + 1))) %>% round(., 6)
          }
        }
        dselect["OBD", trial] <- Aset_ext[which.max(uhat_final[Aset_ext])]
        N_s2[trial, ] <- 0
      } else {
        # timeline
        # completion of dose-finding stage
        stage1_endt <- max(assess_t) + max(DLTwindow, effwindow)
        # completion of last patient in extended backfilling
        stage2_endt <- arrival_t[arrival_t > stage1_endt][sum(n_ext)] + max(DLTwindow, effwindow)
        duration["extended_bf", trial] <- stage2_endt - arrival_t[1]
        
        # outcome collection
        if (sum(n_ext) == 1) {
          ext_label <- which(n_ext > 0)
        } else {
          ext_label <- rep(which(n_ext>0), times = n_ext[n_ext>0]) %>% sample()
        }
         
        ytox_ext <- sapply(seq_along(ext_label), function(rr) {
          replace(numeric(ndose), ext_label[rr], (1 * (DLT_t[ext_label[rr], (seq_len(sum(n_ext)) + sum(n))[rr]] <= DLTwindow)))
        }) %>% rowSums()
        yeff_ext <- sapply(seq_along(ext_label), function(rr) {
          replace(numeric(ndose), ext_label[rr], (1 * (eff_t[ext_label[rr], (seq_len(sum(n_ext)) + sum(n))[rr]] <= effwindow)))
        }) %>% rowSums()
        if (utility.method == "mean") {
          for (dd in Aset_ext) {
            uhat_final[dd] <- ((rowSums(trial_log$y_eff) + yeff_ext)[dd] / (n + n_ext)[dd] - w0 * ((rowSums(trial_log$y_tox) + ytox_ext)[dd] / (n + n_ext)[dd])) %>% round(., 6)
          }
        } else if (utility.method == "pp") {
          quasi_x <- vector(mode = "numeric", length = ndose)
          for (dd in Aset_ext) {
            quasi_x[dd] <-  (utility_score[3] * (rowSums(trial_log$y_eff) + yeff_ext)[dd] + utility_score[2] * ((n + n_ext)[dd] -(rowSums(trial_log$y_tox) + ytox_ext)[dd]))/100
            uhat_final[dd] <- (1 - pbeta(ub/100, (quasi_x[dd] + 1), ((n + n_ext)[dd] - quasi_x[dd] + 1))) %>% round(., 6)
          }
        }
        dselect["OBD", trial] <- Aset_ext[which.max(uhat_final[Aset_ext])]
        N_s2[trial, ] <- n_ext
      }
  
    }
    
  } # end of for loop

  sel_MTD <- sel_OBD <- rep(0, ndose)
  nptsdose <- rbind(EN.s1 = apply(N_s1, 2, mean), EN.s2 = apply(N_s2, 2, mean))
  nptsdose <- rbind(nptsdose, total = colSums(nptsdose))
  ntoxdose <- apply(Ytox_s1, 2, mean)
  for (i in 1:ndose) {
    sel_MTD[i] <- sum(dselect["MTD", ] == i) / ntrial * 100
    sel_OBD[i] <- sum(dselect["OBD", ] == i) / ntrial * 100
  }

  out <- list(
    sel_MTD = sel_MTD, sel_OBD = sel_OBD, npatients = nptsdose,
    ntox = ntoxdose, totaltox = sum(Ytox_s1) / ntrial, 
    totaln = c(s1 = sum(N_s1), s2 = sum(N_s2), total = sum(N_s1+N_s2)) / ntrial,
    percentstop = sum(dselect["MTD", ] == 99) / ntrial * 100, 
    duration = rowMeans(duration),
    simu.setup = data.frame(
      target = target, pT.true = pT.true, pE.true = pE.true, 
      ncohort = ncohort, cohortsize = cohortsize,
      startdose = startdose, pT.saf = pT.saf, pT.tox = pT.tox,
      cutoff.eli = cutoff.eli, extrasafe = extrasafe,
      offset = offset, ntrial = ntrial, dose = 1:ndose
    ),
    flowchart = TRUE, lambda_e = lambda_e, lambda_d = lambda_d
  )
  
  class(out) <- "boin"
  return(out)
}


MAIN.func <- function(target = 0.25, pT.true, pE.true, 
                      accuralrate, DLTwindow, effwindow, w0 = 2/3, 
                      eff.complete = 0.21, nsimu = 5000, seed = 123) {
  ndose <- length(pT.true)
  
  ### our method (no tite)
  result_ours0 <- run.bfboin(
    target = target, pT.true = pT.true, pE.true = pE.true, w0 = w0,
    ncohort = 10, cohortsize = 3, fixed.cohort = TRUE,
    accuralrate = accuralrate, DLTwindow = DLTwindow, effwindow = effwindow,
    n.earlystop = 12, n.backfilling = 12, bf.type = "utility", bf.extended = TRUE,
    utility.method = "mean", startdose = 1, init.size = 3, 
    titration = FALSE, pE.low = 0.20, cutoff.eff = 0.90,
    min.complete = c(tox = 1, eff = 1), min.MF = c(tox = 1, eff = 1),
    ntrial = nsimu, seed = seed
  )
  
  ### BF-BOIN
  result_bf <- run.bfboin(
    target = target, pT.true = pT.true, pE.true = pE.true, w0 = w0,
    ncohort = 10, cohortsize = 3, fixed.cohort = TRUE,
    accuralrate = accuralrate, DLTwindow = DLTwindow, effwindow = effwindow,
    n.earlystop = 12, n.backfilling = 12, bf.type = "tox", bf.extended = FALSE,
    utility.method = "mean", startdose = 1, init.size = 3, 
    titration = FALSE, pE.low = 0.20, cutoff.eff = 0.90,
    min.complete = c(tox = 1, eff = 1), min.MF = c(tox = 1, eff = 1),
    ntrial = nsimu, seed = seed
  )
  
  
  ### our method (mean utility)
  result_ours1 <- run.bfboin(
    target = target, pT.true = pT.true, pE.true = pE.true, w0 = w0,
    ncohort = 10, cohortsize = 3, fixed.cohort = TRUE,
    accuralrate = accuralrate, DLTwindow = DLTwindow, effwindow = effwindow,
    n.earlystop = 12, n.backfilling = 12, bf.type = "utility", bf.extended = TRUE,
    utility.method = "mean", startdose = 1, init.size = 3, 
    titration = FALSE, pE.low = 0.20, cutoff.eff = 0.90,
    min.complete = c(tox = 0.51, eff = eff.complete), min.MF = c(tox = 0.25, eff = 0),
    ntrial = nsimu, seed = seed
  )
  
  ### our method (utility posterior prob)
  result_ours2 <- run.bfboin(
    target = target, pT.true = pT.true, pE.true = pE.true, w0 = w0,
    ncohort = 10, cohortsize = 3, fixed.cohort = TRUE,
    accuralrate = accuralrate, DLTwindow = DLTwindow, effwindow = effwindow,
    n.earlystop = 12, n.backfilling = 12, bf.type = "utility", bf.extended = TRUE,
    utility.method = "pp", startdose = 1, init.size = 3, 
    titration = FALSE, pE.low = 0.20, cutoff.eff = 0.90,
    min.complete = c(tox = 0.51, eff = eff.complete), min.MF = c(tox = 0.25, eff = 0),
    ntrial = nsimu, seed = seed
  )
  
  ### TITE-BOIN-ET
  result_titeet <- tite.boinet(
    n.dose = 5, start.dose = 1, size.cohort = 3, n.cohort = 10, 
    toxprob = pT.true, effprob = pE.true, 
    phi = 0.25, delta = 0.50, tau.T = DLTwindow * 30, tau.E = effwindow * 30, 
    te.corr = 0, accrual = 30 / accuralrate, gen.enroll.time="exponential", 
    stopping.npts = 30, obd.method = "utility.scoring", 
    seed.sim = seed, n.sim = nsimu
    )
  
  
  ### BF-BOIN-ET
  result_bfet <- get.oc.backboinet(
    target_T = target, target_E = 0.25, toxprob = pT.true, effprob = pE.true, 
    n.dose = 5, startdose = 1, ncohort = 10, cohortsize = 3, 
    tau.T = DLTwindow, tau.E = effwindow, te.corr = 0, accrual = accuralrate, 
    gen.enroll.time = "exponential", 
    stopping.npts = 12, suspend = (1 - 0.51), seed.sim = seed, n.sim = nsimu
  )
  
  
  # overdose N
  EN_overdose <- list(ours0 = rep(-1, 3), bf = rep(-1, 3), ours1 = rep(-1, 3), 
                      ours2 = rep(-1, 3), bfet = -1, titeet = -1)
  if (any(pT.true >= target)) {
    maxtol <- which(pT.true>=target) %>% min()
    if (maxtol < (ndose - 1)) {
      EN_overdose$ours0 <- rowSums(result_ours0$npatients[, (maxtol + 1):ndose])
      EN_overdose$bf <- rowSums(result_bf$npatients[, (maxtol + 1):ndose])
      EN_overdose$ours1 <- rowSums(result_ours1$npatients[, (maxtol + 1):ndose])
      EN_overdose$ours2 <- rowSums(result_ours2$npatients[, (maxtol + 1):ndose])
      EN_overdose$bfet <- sum(result_bfet$n.patient[(maxtol + 1):ndose])
      EN_overdose$titeet <- sum(result_titeet$n.patient[(maxtol + 1):ndose])
    } else if (maxtol == (ndose - 1)) {
      EN_overdose$ours0 <- result_ours0$npatients[, (maxtol + 1):ndose]
      EN_overdose$bf <- result_bf$npatients[, (maxtol + 1):ndose]
      EN_overdose$ours1 <- result_ours1$npatients[, (maxtol + 1):ndose]
      EN_overdose$ours2 <- result_ours2$npatients[, (maxtol + 1):ndose]
      EN_overdose$bfet <- sum(result_bfet$n.patient[(maxtol + 1):ndose])
      EN_overdose$titeet <- sum(result_titeet$n.patient[(maxtol + 1):ndose])
    }
  }
  
  
  ### summary
  settings <- rbind(
    True.Tox = pT.true, True.Eff = pE.true,
    True.Utility = (pE.true - w0 * pT.true)
  ) %>% round(2)
  colnames(settings) <- paste("Dose", 1:5)
  
  selection <- rbind(
    MTD.ours0 = c(result_ours0$sel_MTD, result_ours0$percentstop),
    OBD.ours0 = c(result_ours0$sel_OBD, NA),
    MTD.bf = c(result_bf$sel_MTD, result_bf$percentstop),
    OBD.bf = c(result_bf$sel_OBD, NA),
    MTD.ours1 = c(result_ours1$sel_MTD, result_ours1$percentstop),
    OBD.ours1 = c(result_ours1$sel_OBD, NA),
    MTD.ours2 = c(result_ours2$sel_MTD, result_ours2$percentstop),
    OBD.ours2 = c(result_ours2$sel_OBD, NA),
    OBD.titeet = c(result_titeet$prop.select, result_titeet$prop.stop),
    OBD.bfet = c(result_bfet$prop.select, result_bfet$prop.stop)
  ) %>% round(2)
  colnames(selection) <- c(paste("Dose", 1:5), "early.stop")
  
  EN <- rbind(
    ours0 = cbind(result_ours0$npatients, result_ours0$totaln, EN_overdose$ours0) %>%
      { rownames(.) <- c("ours0.s1", "ours0.s2", "ours0.total"); . },
    bf = cbind(result_bf$npatients, result_bf$totaln, EN_overdose$bf) %>%
      { rownames(.) <- c("bf.s1", "bf.s2", "bf.total"); . },
    ours1 = cbind(result_ours1$npatients, result_ours1$totaln, EN_overdose$ours1) %>%
      { rownames(.) <- c("ours1.s1", "ours1.s2", "ours1.total"); . },
    ours2 = cbind(result_ours2$npatients, result_ours2$totaln, EN_overdose$ours2) %>%
      { rownames(.) <- c("ours2.s1", "ours2.s2", "ours2.total"); . },
    titeet = c(result_titeet$n.patient, sum(result_titeet$n.patient), EN_overdose$titeet),
    bfet = c(result_bfet$n.patient, result_bfet$totaln, EN_overdose$bfet)
  )  %>% round(2)
  colnames(EN) <- c(paste("Dose", 1:5), "EN", "EN.Overdose")
  
  duration <- c(
    ours0 = c(result_ours0$duration[1], diff(result_ours0$duration), result_ours0$duration[2]),
    bf = c(result_bf$duration[1], diff(result_bf$duration), result_bf$duration[2]),
    ours1 = c(result_ours1$duration[1], diff(result_ours1$duration), result_ours1$duration[2]),
    ours2 = c(result_ours2$duration[1], diff(result_ours2$duration), result_ours2$duration[2]),
    titeet = result_titeet$duration/30,
    bfet = result_bfet$duration
  )  %>% round(2)
  
  res <- list(
    settings = settings, selection = selection, EN = EN, duration = duration
  )
  
  return(res)
}


### settings ========================
pT.true <- rbind(
  c(0.1, 0.18, 0.35, 0.40, 0.50),
  c(0.05, 0.15, 0.25, 0.35, 0.50),
  c(0.02, 0.06, 0.10, 0.20, 0.35),
  c(0.01, 0.03, 0.05, 0.12, 0.22),
  c(0.10, 0.20, 0.35, 0.43, 0.50),
  c(0.02, 0.06, 0.10, 0.20, 0.35),
  c(0.05, 0.10, 0.20, 0.35, 0.40),
  c(0.01, 0.05, 0.15, 0.18, 0.35)
)
pE.true <- rbind(
  c(0.35, 0.35, 0.37, 0.39, 0.39),
  c(0.10, 0.35, 0.35, 0.38, 0.39),
  c(0.05, 0.10, 0.35, 0.35, 0.40),
  c(0.05, 0.10, 0.15, 0.35, 0.36),
  c(0.10, 0.36, 0.37, 0.40, 0.41),
  c(0.05, 0.10, 0.15, 0.35, 0.37),
  c(0.35, 0.36, 0.37, 0.40, 0.41),
  c(0.05, 0.35, 0.36, 0.37, 0.38)
)

all_config <- expand.grid(
  Scenarrio = 1:nrow(pT.true),
  DLT_window = c(1, 2, 3),
  eff_window = c(1, 2, 3), 
  accural_rate = c(1, 2, 3),
  eff_complete = c(0.21, 0.31, 0.41)
)

### RUN code ========================
output <- MAIN.func(
  pT.true = pT.true[all_config$Scenarrio[sc00], ], 
  pE.true = pE.true[all_config$Scenarrio[sc00], ], 
  accuralrate = all_config$accural_rate[sc00], 
  DLTwindow = all_config$DLT_window[sc00], 
  effwindow = all_config$eff_window[sc00], 
  eff.complete = all_config$eff_complete[sc00], 
  nsimu = 5000
)
output$settings <- cbind(setting.idx = sc00, output$settings) %>% 
  as.data.frame() %>% tibble() %>% 
  mutate(metric = c("True.Tox", "True.Eff", "True.Utility"), .after = 1)
output$selection <- cbind(setting.idx = sc00, output$selection) %>% 
  as.data.frame() %>% tibble() %>% 
  mutate(method = c("MTD.ours0", "OBD.ours0", "MTD.bf", "OBD.bf", "MTD.ours1", "OBD.ours1", "MTD.ours2", "OBD.ours2", "OBD.titeet", "OBD.bfet"), .after = 1)
output$EN <- cbind(setting.idx = sc00, output$EN, duration = output$duration) %>% 
  as.data.frame() %>% tibble() %>% 
  mutate(method = c("ours0.s1", "ours0.s2", "ours0.total", "bf.s1", "bf.s2", "bf.total", "ours1.s1", "ours1.s2", "ours1.total", "ours2.s1", "ours2.s2", "ours2.total", "titeet", "bfet"), .after = 1)



resultpath.setting <- "./results/output_settings.csv"
resultpath.sel <- "./results/output_selection.csv"
resultpath.EN <- "./results/output_EN.csv"

### write file in a parallel way
lockfile <- "./results/lockfile.lock"
lock <- lock(lockfile, timeout = Inf)
if (!file.exists(resultpath.setting)) {
  write_excel_csv(output$settings, file = resultpath.setting, append = T, col_names = TRUE)
} else {
  write_excel_csv(output$settings, file = resultpath.setting, append = T, col_names = FALSE)
}

if (!file.exists(resultpath.sel)) {
  write_excel_csv(output$selection, file = resultpath.sel, append = T, col_names = TRUE, na = "")
} else {
  write_excel_csv(output$selection, file = resultpath.sel, append = T, col_names = FALSE, na = "")
}

if (!file.exists(resultpath.EN)) {
  write_excel_csv(output$EN, file = resultpath.EN, append = T, col_names = TRUE)
} else {
  write_excel_csv(output$EN, file = resultpath.EN, append = T, col_names = FALSE)
}

unlock(lock)
