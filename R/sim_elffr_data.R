#' Generate data to be used for longitudinal FoFR
#'
#' @param func_predictors list of strings of mathematical expressions with s and u as arguments for the response and predictor domain
#' defining the true functional predictor coefficient functions
#' @param family family of response function: "gaussian", "binomial", or "poisson"
#' @param I number of subjects
#' @param J average number of visits for each subject
#' @param L dimension of response domain
#' @param U dimension of predictor domains - list if more than 1 functional predictor, else a number
#' @param SNR_B relative importance of random effects
#' @param SNR_sigma signal-to-noise ratio
#' @param seed integer seed for generating random data
#' @param psi2 whether to include subject-visit specific random effects (otherwise only subject-specific)
#'
#' @export
#'
#' @returns Returns a dataframe with (1) ID - integer for subject identification, (2) visit - visit number for the
#' i'th subject, (3) X. scalar predictors (4) W. - functional predictors, (5) Y - response, and
#' (6) eta - true response / linear predictor value
#'
sim_elffr_data <- function(func_predictors, family = "gaussian", I = 100, J = 5, L = 100, U = 100,
                           SNR_B = 1, SNR_sigma = 1, seed = 1, psi2 = FALSE){
  set.seed(seed)

  # initialize functional predictor objects
  W <- list()
  functional_eff <- list()

  grid <- seq(0, 1, length = L)
  beta_true <- matrix(NA, 3, L)

  # intercept + 2 scalar predictors
  beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
  beta_true[2,] = dnorm(grid, .6, .15)/20
  beta_true[3,] = dnorm(grid, 0.2, 0.1)/60 + dnorm(0.35, 0.1)/200 -
                  dnorm(grid, 0.65, 0.06)/250 + dnorm(1, 0.07)/60

  rownames(beta_true) <- c("Intercept", "x1", "x2")

  psi_true <- matrix(NA, 2, L)
  psi_true[1,] <- (1.5 - sin(2*grid*pi) - cos(2*grid*pi))
  psi_true[1,] <- psi_true[1,] / sqrt(sum(psi_true[1,]^2))
  psi_true[2,] <- sin(4*grid*pi)
  psi_true[2,] <- psi_true[2,] / sqrt(sum(psi_true[2,]^2))

  ## generate number of visits for each subject from poisson distribution
  J_subj <- pmax(rpois(I, J), 1)

  ## generate fixed effects
  n <- sum(J_subj)
  X_des = cbind(1, rnorm(n, 0, 5), rnorm(n, 0, 3))
  fixef <- X_des %*% beta_true

  ## generate random effects
  subj <- as.factor(rep(1:I, J_subj))
  Z_des <- model.matrix( ~ 0 + subj)

  ## simulate score function
  c_true <- mvtnorm::rmvnorm(I, mean = rep(0, 2), sigma = diag(c(3, 1.5)))
  b_true <- c_true %*% psi_true
  ranef = Z_des %*% b_true

  ## by default do not add subject-visit random deviation
  # if(psi2){
  #   psi2_true <- matrix(NA, 2, L)
  #   psi2_true[1,] <- (cos(2*grid*pi) + sin(2*grid*pi))
  #   psi2_true[1,] <- psi2_true[1,] / sqrt(sum(psi2_true[1,]^2))
  #   psi2_true[2,] <- (cos(4*grid*pi) - sin(4*grid*pi))
  #   psi2_true[2,] <- psi2_true[2,] / sqrt(sum(psi2_true[2,]^2))
  #
  #   c2_true <- rmvnorm(n, mean = rep(0, 2), sigma = diag(c(3, 1.5)))
  #   ranef <- ranef + c2_true %*% psi2_true
  # }

  # generate functional effects
  for (i in 1:length(func_predictors)){
    ff.formula <- function(s, u, formula_str = func_predictors[[i]])
    {
      eval(parse(text = formula_str))
    }

    famm_dat <- pffrSim_modified(ff.formula, scenario = "ff",
                                 n = n, nxgrid = U[[i]], nygrid = L, SNR = 10)
    W[[i]] <- famm_dat$X1
    functional_eff[[i]] <- attributes(famm_dat)$truth$etaTerms$X1
  }

  # combination of fixed effects (scalar + functional)
  full_fixef <- fixef + Reduce("+", functional_eff)

  ## generate linear predictors
  ranef <- sd(fixef)/sd(ranef)/SNR_B*ranef ## adjust for relative importance of random effects
  eta_true <- full_fixef + ranef

  ## generate longitudinal functional data
  Y_obs <- matrix(NA, n, L)
  p_true <- plogis(eta_true)
  lam_true <- exp(eta_true)
  sd_signal <- sd(eta_true)
  for(i in 1:n){
    for(j in 1:L){
      if(family == "gaussian"){
        Y_obs[i, j] <- rnorm(1, mean = eta_true[i, j], sd = sd_signal/SNR_sigma)
      }else if(family == "binomial"){
        Y_obs[i, j] <- rbinom(1, 1, p_true[i, j])
      }else if(family == "poisson"){
        Y_obs[i, j] <- rpois(1, lam_true[i, j])
      }
    }
  }

  ## combine simulated data
  visit <- rep(1:I, J_subj)
  for(i in 1:I){
    visit[which(subj == i)] <- 1:J_subj[i]
  }

  dat.sim <- data.frame(ID = subj, visit = visit, X1 = X_des[,2], X2 = X_des[,3],
                        Y = I(Y_obs), eta = I(eta_true))

  for (i in 1:length(func_predictors)){
    ff_name <- paste0("W", i)
    dat.sim[[ff_name]] <- W[[i]]
  }

  return(dat.sim)
}

