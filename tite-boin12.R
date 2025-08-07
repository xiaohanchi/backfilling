get.boundary <- function(target, targetE, ncohort, cohortsize = 3, p.saf = NA, p.tox = NA, cutoff.eli = 0.95,
                         cutoff.eli.E = 0.90) {
  # if the user does not provide p.saf and p.tox, use the default values
  if (is.na(p.saf)) p.saf <- 0.6 * target
  if (is.na(p.tox)) p.tox <- 1.4 * target

  ### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
  npts <- ncohort * cohortsize
  ntrt <- NULL
  b.e <- NULL
  b.d <- NULL
  elim <- NULL
  elimE <- NULL

  for (n in (1:ncohort) * cohortsize)
  {
    error.min <- 3
    for (m1 in 0:(n - 1))
    {
      for (m2 in (m1 + 1):n)
      {
        error1 <- pbinom(m1, n, target) + 1 - pbinom(m2 - 1, n, target)
        error2 <- 1 - pbinom(m1, n, p.saf)
        error3 <- pbinom(m2 - 1, n, p.tox)

        error <- error1 + error2 + error3
        if (error < error.min) {
          error.min <- error
          cutoff1 <- m1
          cutoff2 <- m2
        }
      }
    }
    ntrt <- c(ntrt, n)
    b.e <- c(b.e, cutoff1)
    b.d <- c(b.d, cutoff2)

    elimineed <- 0 # indicating whether elimination is needed
    elimineedE <- 0
    if (n < 3) {
      elim <- c(elim, NA)
      elimE <- c(elimE, NA)
    } # require treating at least 3 patients before eliminating a dose
    else {
      for (ntox in 3:n) # determine elimination boundary, prior beta(1,1) is used in beta-binomial model
      {
        if (1 - pbeta(target, ntox + 1, n - ntox + 1) > cutoff.eli) {
          elimineed <- 1
          break
        }
      }
      if (elimineed == 1) {
        elim <- c(elim, ntox)
      } else {
        elim <- c(elim, NA)
      } # set the elimination boundary large such that no elimination will actually occurs

      for (neff in n:0) {
        if (pbeta(targetE, neff + 1, n - neff + 1) > cutoff.eli.E) {
          elimineedE <- 1
          break
        }
      }
      if (elimineedE == 1) {
        elimE <- c(elimE, neff)
      } else {
        elimE <- c(elimE, NA)
      }
    }
  }
  for (i in 1:length(b.d)) {
    if (!is.na(elim[i]) && (b.d[i] > elim[i])) b.d[i] <- elim[i]
  }
  boundaries <- rbind(ntrt, elim, b.d, b.e, elimE)
  rownames(boundaries) <- c(
    "Number of patients treated", "Eliminate if # of DLT >=",
    "Deescalate if # of DLT >=", "Escalate if # of DLT <=", "Eliminate if # of Eff <="
  )
  colnames(boundaries) <- rep("", ncohort)

  return(boundaries)
}


get.pts.outcome <- function(d, dist, coh.size, pE.true, pT.true, r, A, alpha = c(0.5, 0.5)) {
  d.pE <- pE.true[d]
  d.pT <- pT.true[d]
  # from marginal prob to joint prob
  get.joint.p <- function(d.pE, d.pT, r) {
    p.ab <- function(d.pE, d.pT, a, b, r) {
      d.pE^a * (1 - d.pE)^(1 - a) * d.pT^b * (1 - d.pT)^(1 - b) + (-1)^(a + b) * d.pE * (1 - d.pE) * d.pT * (1 - d.pT) * ((exp(r) - 1) / (exp(r) + 1))
    }
    p <- c(p.ab(d.pE, d.pT, a = 1, b = 0, r), p.ab(d.pE, d.pT, a = 0, b = 0, r), p.ab(d.pE, d.pT, a = 1, b = 1, r), p.ab(d.pE, d.pT, a = 0, b = 1, r))
    return(p)
  }

  library(truncdist)
  joint.p <- get.joint.p(d.pE, d.pT, r = r)
  rmultinom(coh.size, size = 1, prob = joint.p) -> ys
  ys
  E.out <- unlist(lapply(1:coh.size, function(x) ifelse((ys[1, x] == 1 | ys[3, x] == 1), 1, 0)))
  E.out
  T.out <- unlist(lapply(1:coh.size, function(x) ifelse((ys[3, x] == 1 | ys[4, x] == 1), 1, 0)))
  T.out

  if (dist == "Uniform") {
    T.time <- ifelse(T.out == 1, runif(1, 0, A[1]), A[1] + 0.001)
    E.time <- ifelse(E.out == 1, runif(1, 0, A[2]), A[2] + 0.001)
  }
  ## Generate time to efficacy or toxicity
  ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
  if (dist == "Weibull") {
    # Generate time to DLT
    # T.shape=log(log(1-d.pT)/log(1-d.pT/2))/log(2)
    T.pihalft <- d.pT * alpha[1]
    T.shape <- log(log(1 - d.pT) / log(1 - T.pihalft)) / log(2)
    T.scale <- A[1] / exp(log(-log(1 - d.pT)) / T.shape)
    #  print(paste("T.scale=",T.scale, "and T.shape=",T.shape))
    T.time <- ifelse(T.out == 1, rtrunc(1, "weibull", a = 0, b = A[1], shape = T.shape, scale = T.scale), A[1] + 0.001)
    T.time

    #   #Generate time to tumor response
    E.pihalft <- d.pE * alpha[2]
    E.shape <- log(log(1 - d.pE) / log(1 - E.pihalft)) / log(2)
    E.scale <- A[2] / exp(log(-log(1 - d.pE)) / E.shape)
    # print(paste("E.scale=",E.scale, "and E.shape=",E.shape))
    E.time <- ifelse(E.out == 1, rtrunc(1, "weibull", a = 0, b = A[2], shape = E.shape, scale = E.scale), A[2] + 0.001)
    E.time
  }




  list(eff = E.out, tox = T.out, t.tox = T.time, t.eff = E.time)
}

get.tite.boin12.approx <- function(pT.true, pE.true, rho, targetT = 0.35, targetE = 0.25, ncohort = 10, cohortsize = 3, startdose = 1,
                                   accrual.rate = 2, arrive.dist = "Uniform", event.dist = "Weibull", A = c(2, 3), maxpen = 0.5, n.earlystop = ncohort * cohortsize, p.saf = 0.6 * targetT, p.tox = 1.4 * targetT, cutoff.eli.T = 0.95,
                                   N1 = 6, N2 = 9, cutoff.eli.E = 0.90, u11 = 60, u00 = 40, ntrial = 100, rseed = 200, alpha = c(0.5, 0.5), RunIn = FALSE) {
  set.seed(rseed)

  ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1) {
      return(x)
    }
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol))) {
        break
      }
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl]) / sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }

  ## Set up before simulation
  ndose <- length(pT.true)
  npts <- ncohort * cohortsize
  YT <- matrix(rep(0, ndose * ntrial), ncol = ndose) # store toxicity outcome
  YE <- matrix(rep(0, ndose * ntrial), ncol = ndose) # store efficacy outcome
  N <- matrix(rep(0, ndose * ntrial), ncol = ndose) # store the number of patients
  dselect_mtd <- dselect <- rep(0, ntrial) # store the selected dose level
  durationV <- rep(0, ntrial) # record trial duratoin
  npendV <- rep(0, ntrial) # number of times suspending the trial
  sel_mtd <- sel <- rep(0, ndose)
  pts <- rep(0, ndose)
  dlt <- rep(0, ndose)
  eff <- rep(0, ndose)
  poorall <- 0
  temp <- get.boundary(targetT, targetE, ncohort, cohortsize, cutoff.eli = cutoff.eli.T, cutoff.eli.E = cutoff.eli.E)
  b.e <- temp[4, ] # escalation boundary
  b.d <- temp[3, ] # deescalation boundary
  b.elim <- temp[2, ] # elimination boundary
  b.elimE <- temp[5, ]
  p.saf <- targetT * 0.6
  p.tox <- targetT * 1.4
  lambda1 <- log((1 - p.saf) / (1 - targetT)) / log(targetT * (1 - p.saf) / (p.saf * (1 - targetT)))
  lambda2 <- log((1 - targetT) / (1 - p.tox)) / log(p.tox * (1 - targetT) / (targetT * (1 - p.tox)))
  u01 <- 100
  u10 <- 0
  utility <- c(u11, u10, u01, u00)
  # Assume independence between toxicity and efficacy
  targetP <- c(targetE * targetT, targetT * (1 - targetE), (1 - targetT) * targetE, (1 - targetT) * (1 - targetE))

  # Calculate the benchmark utility
  uu <- sum(targetP * utility) # highest unacceptable utility
  uu <- uu + (100 - uu) / 2 # benchmark utility (i.e., desirable utility)


  # Calculate true utility
  # from marginal prob to joint prob
  get.joint.p <- function(pE, pT, r) {
    p.ab <- function(pE, pT, a, b, r) {
      pE^a * (1 - pE)^(1 - a) * pT^b * (1 - pT)^(1 - b) + (-1)^(a + b) * pE * (1 - pE) * pT * (1 - pT) * ((exp(r) - 1) / (exp(r) + 1))
    }
    p <- c(p.ab(pE, pT, a = 1, b = 0, r), p.ab(pE, pT, a = 0, b = 0, r), p.ab(pE, pT, a = 1, b = 1, r), p.ab(pE, pT, a = 0, b = 1, r))
    return(p)
  }


  p10 <- p01 <- p00 <- p11 <- rep(0, ndose)

  for (d in 1:ndose) {
    p.joint <- get.joint.p(pE = pE.true[d], pT = pT.true[d], r = rho)
    p01[d] <- p.joint[1]
    p00[d] <- p.joint[2]
    p11[d] <- p.joint[3]
    p10[d] <- p.joint[4]
  }
  u.true <- u11 * p11 + u00 * p00 + 100 * p01
  for (trial in 1:ntrial) {
    # if (trial %% 500 == 0) print(paste("ntrial=", trial))
    ## record individual patient data
    yT.true <- yE.true <- NULL # toxicity/efficacy indicator for each subject
    dv <- NULL # dose for each subject

    t.enter <- NULL # time enter the study
    t.eff.event <- t.tox.event <- NULL # time to event
    t.decision <- 0 # decision making time


    ## summarize data
    n.yT <- n.yE <- rep(0, ndose) ## number of observed DLT/efficacy at each dose level
    y01 <- y10 <- y11 <- y00 <- rep(0, ndose) ## number of different outcomes at each dose level
    n <- ess.T <- ess.E <- rep(0, ndose) ## number of patients, effective toxicity, efficacy sample size
    earlystop <- 0 ## indicate whether the trial terminates early
    # =1;ndose=5
    d <- startdose ## current dose level
    elimi <- rep(0, ndose) ## indicate whether doses are eliminated due to toxicity
    elimiE <- rep(0, ndose) ## indicate whether doses are eliminated due to efficacy
    safe <- 0
    npend <- 0
    posH <- rep(1 - uu / 100, ndose)

    for (i in 1:ncohort)
    {
      # print(paste0("ncohort=",i))
      # time of entry
      for (j in 1:cohortsize)
      {
        if (j == 1) {
          t.enter <- c(t.enter, t.decision)
        } else {
          if (arrive.dist == "Uniform") {
            t.enter <- c(t.enter, t.enter[length(t.enter)] + runif(1, 0, 2 / accrual.rate))
          }
          if (arrive.dist == "Exponential") {
            t.enter <- c(t.enter, t.enter[length(t.enter)] + rexp(1, rate = accrual.rate))
          }
        }
      }

      # generate data for the new patients from the truth

      obscohort <- get.pts.outcome(
        d = d, dist = event.dist, coh.size = cohortsize, pE.true = pE.true, pT.true = pT.true,
        r = rho, A = A, alpha = alpha
      )
      t.eff.event <- c(t.eff.event, obscohort$t.eff)
      t.tox.event <- c(t.tox.event, obscohort$t.tox)
      yE.true <- c(yE.true, obscohort$eff)
      yT.true <- c(yT.true, obscohort$tox)

      n[d] <- n[d] + cohortsize


      dv <- c(dv, rep(d, cohortsize))

      t.decision <- t.enter[length(t.enter)]

      nobs <- -1
      pending <- 1

      npend <- npend - 1

      if (n[d] >= n.earlystop) {
        break
      }
      if (earlystop == 1) {
        t.decision <- max(t.enter) + max(A)
        break
      }
      while (pending == 1) {
        npend <- npend + 1
        pending <- 0
        # update decision time
        if (i == ncohort) {
          t.decision <- t.decision + max(A)
        } else {
          if (arrive.dist == "Uniform") {
            t.decision <- t.decision + runif(1, 0, 2 / accrual.rate)
          }
          if (arrive.dist == "Exponential") {
            t.decision <- t.decision + rexp(1, rate = accrual.rate)
          }
        }




        cset <- (dv == d) # current dose indicator

        # indicator for whether or not toxicity and efficacy is observed. O means pending.
        delta.tox <- as.numeric(((t.enter + t.tox.event) < t.decision))[cset]
        delta.eff <- as.numeric(((t.enter + t.eff.event) < t.decision))[cset]

        # number of patients who have pending outcomes
        n.pend <- length(unique(c(which(delta.tox == 0), which(delta.eff == 0))))
        n.pend



        # check whether the trial should be suspended, where the number of patients on current dose is sum(cset)
        if (n.pend > sum(cset) * maxpen) {
          pending <- 1
        } else {
          # n[d]<-n[d]+cohortsize

          ## update utility values for doses that have been used to treat patients
          for (dd in which(n != 0)) {
            cset <- (dv == dd)
            # individual follow-up time
            follow.t.tox <- pmin((t.decision - t.enter)[cset], A[1])
            follow.t.eff <- pmin((t.decision - t.enter)[cset], A[2])
            # indicator for whether or not toxicity and efficacy is observed. O means pending.
            delta.tox <- as.numeric(((t.enter + t.tox.event) < t.decision))[cset]
            delta.eff <- as.numeric(((t.enter + t.eff.event) < t.decision))[cset]
            ## calculate \tilde{n_T},\tilde{n_E}: number of observed toxicity/ efficacy
            t.nT <- sum(yT.true[cset][delta.tox == 1] == 1)
            t.nE <- sum(yE.true[cset][delta.eff == 1] == 1)


            ## calculate m_T, m_E: complete assessment window, but no DLT/ efficacy
            mT <- sum(yT.true[cset][delta.tox == 1] == 0)
            mE <- sum(yE.true[cset][delta.eff == 1] == 0)

            ## calculate t_T/A_T, t_T/A_E
            ft.T <- ifelse(sum(delta.tox == 0) == 0, 0, sum(follow.t.tox[delta.tox == 0]) / A[1])
            ft.E <- ifelse(sum(delta.eff == 0) == 0, 0, sum(follow.t.eff[delta.eff == 0]) / A[2])

            ## effective sample size for toxicity and efficacy
            ess.T[dd] <- t.nT + mT + ft.T
            ess.E[dd] <- t.nE + mE + ft.E


            n.yT[dd] <- t.nT
            n.yE[dd] <- t.nE


            # Number of observations: use effective sample size
            n.currT <- ess.T[dd]
            n.currE <- ess.E[dd]

            # determine whether current dose level is toxic (don't use ESS because we don't want to accidently eliminate doses)
            if (1 - pbeta(targetT, n.yT[dd] + 1, n[dd] - n.yT[dd] + 1) > cutoff.eli.T & n[dd] > 3) {
              elimi[dd:ndose] <- 1
              if (d == 1) {
                earlystop <- 1
                break
              }
            }

            # determine whether current dose level is futile
            # print(n)
            # print(paste("dd=",dd,", n.currE=",n.currE,", n.yE[dd]=",n.yE[dd]))
            if (pbeta(targetE, n.yE[dd] + 1, n.currE - n.yE[dd] + 1) > cutoff.eli.E) {
              elimi[dd] <- 1
            }





            ## calculate utility for quasi binomial data
            # For observed set: patients with both efficacy and toxicity observed

            pE.approx <- n.yE[dd] / ess.E[dd]
            pT.approx <- n.yT[dd] / ess.T[dd]
            O.ind <- which(delta.tox == 1 & delta.eff == 1)
            u_O <- (sum(yT.true[cset][O.ind] == 0 & yE.true[cset][O.ind] == 1) * u01 +
              sum(yT.true[cset][O.ind] == 0 & yE.true[cset][O.ind] == 0) * u00 +
              sum(yT.true[cset][O.ind] == 1 & yE.true[cset][O.ind] == 1) * u11) / 100

            ## For pending set
            # (1) toxicity is observed but efficacy is pending
            tox.seen <- which(delta.tox == 1 & delta.eff == 0)
            u_tox.seen <- 0
            if (!is.null(tox.seen)) {
              temp.pE <- pE.approx * (1 - follow.t.eff[tox.seen] / A[2]) /
                (1 - pE.approx * follow.t.eff[tox.seen] / A[2])

              u_tox.seen <- sum(c(
                (yT.true[cset][tox.seen] == 0) * temp.pE * u01,
                (yT.true[cset][tox.seen] == 0) * (1 - temp.pE) * u00,
                (yT.true[cset][tox.seen] == 1) * temp.pE * u11
              )) / 100
            }
            # (2) efficacy is observed but toxicity is pending
            eff.seen <- which(delta.tox == 0 & delta.eff == 1)
            u_eff.see <- 0

            if (!is.null(eff.seen)) {
              temp.pT <- pT.approx * (1 - follow.t.tox[eff.seen] / A[1]) /
                (1 - pT.approx * follow.t.tox[eff.seen] / A[1])
              u_eff.seen <- sum(c(
                (1 - temp.pT) * (yE.true[cset][eff.seen] == 1) * u01,
                (1 - temp.pT) * (yE.true[cset][eff.seen] == 0) * u00,
                temp.pT * (yE.true[cset][eff.seen] == 1) * u11
              )) / 100
            }


            # (3) both efficacy and toxicity are pending
            both.miss <- which(delta.tox == 0 & delta.eff == 0)
            u_both.miss <- 0

            temp.pE <- pE.approx * (1 - follow.t.eff[both.miss] / A[2]) /
              (1 - pE.approx * follow.t.eff[both.miss] / A[2])

            temp.pT <- pT.approx * (1 - follow.t.tox[both.miss] / A[1]) /
              (1 - pT.approx * follow.t.tox[both.miss] / A[1])

            if (!is.null(both.miss)) {
              u_both.miss <- sum((1 - temp.pT) * temp.pE * u01 +
                (1 - temp.pT) * (1 - temp.pE) * u00 +
                temp.pT * temp.pE * u11) / 100
            }

            ## "The number of events for the quasi-binomial data
            u_curr <- u_O + u_tox.seen + u_eff.seen + u_both.miss

            posH[dd] <- 1 - pbeta(uu / 100, 1 + u_curr, n[dd] - u_curr + 1)
          }
          # stop the trial if all doses are eliminated
          if (sum(elimi == 1) == ndose) {
            earlystop <- 1
            break
          }


          posH <- posH * (1 - elimi)
          if (n[d] >= N1) {
            safe <- 1
          } else {
            safe <- 0
          }




          # print("n.yT")
          #
          # print(n.yT)
          # print("ess.T")
          # print(ess.T)
          ## BOIN12 dose assignment: DLT probability was calculated using observed #ofDLT divided by effective sample size
          # Number of observationsL: use effective sample size
          n.currT <- ess.T[d]
          ## added the 3+3 RunIn on 4/2/2021
          if (n.currT == 3 & n.yT[d] == 1 & targetT == 0.25 & RunIn) {
            ## use 3+3 Run-in
            d_opt <- d
          } else if (n.yT[d] / n.currT >= lambda2 && d != 1) {
            if (sum(elimi[1:(d - 1)] == 0) > 0) { # there is at least one dose below d still available. de-escalate to the highest dose amont those available ones.
              d_opt <- max(which(elimi[1:(d - 1)] == 0))
            } else {
              if (elimi[d] == 1) { # if current dose is eliminated, stop the trial.
                earlystop <- 1
                break
              } else {
                d_opt <- d
              }
            }
          } else if (n.yT[d] / n.currT >= lambda2 && d == 1) {
            if (elimi[d] == 0) {
              d_opt <- d
            } else {
              earlystop <- 1
              break
            }
          } else {
            admi_set <- d
            if (d > 1) {
              if (sum(elimi[1:(d - 1)] == 0) > 0) {
                admi_set <- c(admi_set, max(which(elimi[1:(d - 1)] == 0)))
              }
            }
            if (d < ndose) {
              if (safe == 0) {
                if (sum(elimi[(d + 1):ndose] == 0) > 0) {
                  admi_set <- c(admi_set, d + min(which(elimi[(d + 1):ndose] == 0)))
                }
              } else {
                if (n.yT[d] / n.currT <= lambda1 & sum(elimi[(d + 1):ndose] == 0) > 0) {
                  admi_set <- c(admi_set, d + min(which(elimi[(d + 1):ndose] == 0)))
                }
              }
            }

            temp.posH <- posH[admi_set] + runif(length(admi_set)) * (10^-15)
            d_opt <- admi_set[which.max(temp.posH)]
          }

          if (elimi[d_opt] == 1) {
            earlystop <- 1
            break
          }
          if (sum(elimi) == ndose) {
            earlystop <- 1
            break
          }

          if (d < ndose) {
            if (sum(elimi[(d + 1):ndose] == 0) > 0) {
              d_temp <- d + min(which(elimi[(d + 1):ndose] == 0))
              if (n[d] >= N2 & n[min(d_temp, ndose)] == 0 & n.yT[d] / n.currT < lambda2) {
                d_opt <- d_temp
              }
            }
          }

          d <- d_opt
        }
      }

      #  if(earlystop==1){t.decision = max(t.enter)+max(A); break;}
    }


    yT <- sapply(1:ndose, function(x) sum(yT.true[dv == x]))
    yE <- sapply(1:ndose, function(x) sum(yE.true[dv == x]))
    y01 <- y00 <- y11 <- y10 <- rep(0, ndose)
    for (d in 1:ndose) {
      y01[d] <- sum(yT.true[which(dv == d)] == 0 & yE.true[which(dv == d)] == 1)
      y00[d] <- sum(yT.true[which(dv == d)] == 0 & yE.true[which(dv == d)] == 0)
      y11[d] <- sum(yT.true[which(dv == d)] == 1 & yE.true[which(dv == d)] == 1)
      y10[d] <- sum(yT.true[which(dv == d)] == 1 & yE.true[which(dv == d)] == 0)
    }

    YT[trial, ] <- yT
    YE[trial, ] <- yE

    N[trial, ] <- n
    durationV[trial] <- t.decision
    npendV[trial] <- npend

    ## Code is the same as in BOIN12 hereafter
    if (earlystop == 0) {
      pT <- rep(targetT, ndose)
      pE <- rep(targetE, ndose)
      pT[n != 0] <- (yT[n != 0]) / (n[n != 0])
      pE[n != 0] <- (yE[n != 0]) / (n[n != 0])
      pT <- pava(pT, w = 1 / ((yT + 0.05) * (n - yT + 0.05) / ((n + 0.1)^2 * (n + 0.1 + 1)))) + 0.001 * seq(1, ndose)
      d_mtd <- which.min(abs(pT - targetT))
      pT[n != 0] <- (yT[n != 0]) / (n[n != 0])
      u <- (u01 * y01 + u10 * y10 + u11 * y11 + u00 * y00) / 100
      u <- (u + 1) / (n + 2)
      u[elimi == 1] <- -100
      u[elimiE == 1] <- -100
      u[n == 0] <- -100
      d_opt <- which.max(u[1:d_mtd])
      
      dselect_mtd[trial] <- d_mtd
      dselect[trial] <- d_opt
      sel_mtd[dselect_mtd[trial]] <- sel_mtd[dselect_mtd[trial]] + 1 / ntrial * 100
      sel[dselect[trial]] <- sel[dselect[trial]] + 1 / ntrial * 100
    } else {
      dselect_mtd[trial] <- dselect[trial] <- 99
    }
    pts <- pts + n / ntrial
    dlt <- dlt + yT / ntrial
    eff <- eff + yE / ntrial
  }
  sel_mtd <- round(sel_mtd, 1)
  sel <- round(sel, 1)
  pts <- round(pts, 1)
  dlt <- round(dlt, 1)
  u.true <- round(u.true, 1)
  earlystop <- sum(dselect == 99) / ntrial * 100
  results <- list(
    pT = pT.true, pE = pE.true, u = u.true, sel_mtd = sel_mtd, sel = sel,
    pts = pts, dlt = dlt, eff = eff, earlystop = earlystop, duration = mean(durationV)
  )
  return(results)
}
