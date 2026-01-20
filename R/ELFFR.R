#' Fit longitudinal function-on-function regression models using an efficient 3 step approach
#'
#' `elffr` (efficient longitudinal function-on-function regression) fits longitudinal function-on-function regression models by (1) fitting
#' longitudinal scalar-on-function models with a random intercept at each index of the response domain, (2) smoothing pointwise
#' fixed effect estimates by univariate/bivariate smoothers and (3) constructing confidence regions for coefficient estimates
#' by either an analytic approach in the case of Gaussian responses or by a bootstrap procedure for general response distributions.
#'
#' @usage elffr(formula, subj_var, visit_var, data, num_cores = 8, kz = 30, kb = 30, nknots_scalar = 10, nknots_fbps = 5,
#'                 bspline_p = 3, diff_penalty = 2, smooth.cov = TRUE, analytic = TRUE,
#'                 parallel = TRUE, var = TRUE, num_boot = 100, CMA = TRUE,
#'                 family = "gaussian", argvals.functional = NULL, silent = TRUE,
#'                 seed = 1)
#'
#' @param formula Formula object with ff() to specify functional predictors
#' @param subj_var Name of subject ID variable in data
#' @param visit_var Name of visit/longitudinal variable in data
#' @param data Data object generated from the `lfofr_sim_data` function
#' @param num_cores If using parallel computing, the number of cores to be used for parallel processing
#' @param kz dimension of functional principal component basis for functional predictors
#' @param kb number of knots in the spline basis for the functional coefficient
#' @param nknots_scalar number of knots for smoothing pointwise estimates for the scalar predictor coefficient
#' @param nknots_fbps number of knots for smoothing pointwise estimates for the functional predictor
#' by the sandwich smoother (refund::fbps). Can be a single number or a list of knots on the domain [0,1] for each dimension
#' @param bspline_p degree of B-splines in sandwich smoother
#' @param diff_penalty order of differencing penalty in sandwich smoother
#' @param smooth.cov logical - whether to smooth covariance matrix in pointwise longitudinal scalar-on-function model
#' @param analytic logical - whether to use analytic inference
#' @param parallel logical - whether to use parallel computing
#' @param var logical - whether to compute variance estimates (else only returns estimated coefficients)
#' @param num_boot if using bootstrap inference, number of replicates
#' @param CMA logical - whether to compute correlation and multiplicity adjusted (CMA) confidence bands
#' @param family character - set to "gaussian" for normally-distributed outcomes
#' @param argvals.functional list specifying numeric arguments of functional predictors
#' in the same order as the ff() terms supplied to `formula`
#' Used for adjusting scale of estimates. If NA, defaults to equally spaced interval 1,...,U.
#' Note, ELFFR only works for data on regularly spaced grids.
#' @param seed seed used for drawing samples in bootstrap inference
#'
#' @returns For the case of variance estimation with analytic inference:
#' \item{`scalar.est`}{estimated coefficient curves for scalar predictors}
#' \item{`functional.est`}{estimated coefficient surfaces for functional predictors}
#' \item{`var.scalar`}{variance estimates for scalar predictor coefficient curves}
#' \item{`var.functional`}{variance estimates for functional predictor coefficient surfaces}
#' \item{`time`}{time to fit in seconds}
#' \item{`CMA_quantiles`}{quantiles for CMA confidence regions}
#' @export
#' @importFrom mgcv s
elffr <- function(formula, subj_var, visit_var, data, num_cores = 8, kz = 30, kb = 30, nknots_scalar = 10, nknots_fbps = 5,
                  bspline_p = 3, diff_penalty = 2, smooth.cov = TRUE, analytic = TRUE,
                  parallel = TRUE, var = TRUE, num_boot = 100, CMA = TRUE,
                  family = "gaussian", argvals.functional = NULL, silent = TRUE,
                  seed = 1){
  #############################################
  ##### Step 0: Setup
  #############################################

  ## If parallel, get number of cores
  if (parallel & !is.integer(num_cores)){
    num_cores <- as.integer(round(parallel::detectCores() * 0.75))
    if (!silent){message(paste("Number of cores for parallelization:", num_cores))}
  }

  ## For non-Gaussian family, only do bootstrap inference
  if (family != "gaussian") analytic <- FALSE

  ## Organize objects specified in formula for use in `lpfr`
  model.char <- as.character(formula)
  Y_label <- model.char[[2]] # outcome
  varnames <- unlist(strsplit(model.char[[3]], "\\s*\\+\\s*"))

  # get functional predictor variable names
  ff_indices <- grep("^ff\\(.*?\\)$", varnames)
  ff_varnames <- sub("^ff\\((.*?)\\)$", "\\1", varnames[ff_indices])

  # get scalar predictor variable names
  scalar_varnames <- varnames[-ff_indices]

  #############################################
  ##### Step 1: Fit pointwise lpfr models
  #############################################
  ptm <- proc.time()

  # Get parameters from data
  I <- length(unique(data[[subj_var]]))
  L <- dim(data[[Y_label]])[2]
  argvals <- 1:L
  U.list <- sapply(ff_varnames, function(ff_name) dim(data[[ff_name]])[2])
  p <- length(varnames)
  p.ff <- length(ff_varnames)
  p.scalar <- length(scalar_varnames)
  n <- nrow(data)

  # get argvals and step size for functional predictors
  ff.stepsize <- list()
  if (is.null(argvals.functional)){
    argvals.functional <- list()
    for (i in 1:p.ff){
      argvals.functional[[i]] <- seq(1, U.list[[i]], length = U.list[[i]])
    }
    ff.stepsize <- rep(1, p.ff)
  }else{
    for (i in 1:p.ff){
      ff.stepsize[[i]] <- (max(argvals.functional[[i]]) - min(argvals.functional[[i]])) / (length(argvals.functional[[i]]))
    }
  }

  if (family == "gaussian" & analytic & var & any(c(L, U.list) > 200)){
    message(paste("One or more functional variables are dense! (>200)\n",
                  "Consider subsampling along the functional domains",
                  "or using bootstrap inference"))
  }


  if (analytic == TRUE){ # Get design matrices for later if using analytic inference
    # Fit temporary model to get X_star
    s_l <- 1

    # Create matrix of covariates
    covariate_matrix <- as.matrix(data[,scalar_varnames])

    # Fit LPFR model
    fit.sl <- refund::lpfr(Y = data[[Y_label]][,s_l], subj = data[[subj_var]],
                           covariates = covariate_matrix,
                           funcs = lapply(ff_varnames, function(ff_name) data[[ff_name]]),
                           method = "REML", kz = kz,
                           kb = kb, smooth.cov = smooth.cov)

    Z <- fit.sl$X[,((p.scalar+2):(I+p.scalar+1))]
    X.star <- fit.sl$X[,-((p.scalar+2):(I+p.scalar+1))]

    rm(covariate_matrix, fit.sl)
  }

  ## Build marginal model
  fit_lpfr <- function(s_l){
    # Create matrix of covariates
    covariate_matrix <- as.matrix(data[,scalar_varnames])

    # Get list of functional predictors
    func_list <- lapply(ff_varnames, function(ff_name) as.matrix(data[[ff_name]]))

    # Fit pointwise LPFR model
    fit.sl <- refund::lpfr(Y = data[[Y_label]][,s_l], subj = data[[subj_var]],
                           covariates = covariate_matrix,
                           funcs = func_list,
                           method = "REML", kz = kz,
                           kb = kb, smooth.cov = smooth.cov)


    if(analytic == TRUE){
      ## Get estimated random effect variance components + penalty tuning parameter
      vcov.obj <- mgcv::gam.vcomp(fit.sl$fit)
      if (is.list(vcov.obj)){
        ranef_sd <- vcov.obj$vc
      } else{
        ranef_sd <- vcov.obj
      }
      ranef_var <- ranef_sd^2

      sigma2_err <- ranef_var[,1][nrow(ranef_var)][[1]]
      HHat <- ranef_var[,1][1:(nrow(ranef_var)-1)]

      sigma2_subj <- HHat[[1]] # variance estimate of random intercept
      lambdas <- (1 / HHat[-c(1,nrow(ranef_var))]) # penalty tuning parameters = 1/(variance of random coefs)

      ## Get fixed effect coefficient estimates (scalar + functional predictors)
      beta.star <- fit.sl$fit$coefficients[-((p.scalar+2):(p.scalar+I+1))]

      return(list(betaHat = fit.sl$beta.covariates, gammaHat = fit.sl$BetaHat, sigma2_err = sigma2_err,
                  sigma2_subj = sigma2_subj, lambdas = lambdas, beta.star = beta.star))
    }else{
      return(list(betaHat = fit.sl$beta.covariates, gammaHat = fit.sl$BetaHat))
    }
  }

  if (!silent){cat("Step 1: Fit Massively Univariate LPFR Models\n")}

  ## Fit many lpfr models
  if (parallel == TRUE){
    if(.Platform$OS.type == "windows"){
      cl <- parallel::makePSOCKcluster(num_cores)
      # add global variables to clusters
      parallel::clusterExport(cl, varlist = c("data", "scalar_varnames", "ff_varnames",
                                    "Y_label", "subj_var", "kz", "kb",
                                    "smooth.cov", "analytic", "p.scalar", "I"),
                    envir = environment())
      # add necessary packages to clusters
      parallel::clusterEvalQ(cl, {
        library(refund)
        library(mgcv)
      })
      fits <- tryCatch(
        parallel::parLapply(cl = cl, X = 1:L,
                  fun = fit_lpfr),
        warning = function(w){
          print(paste0("Warning during model fitting:", "\n", w))
        },
        error = function(e){
          stop(paste0("Error during model fitting:", "\n", e))
        },
        finally = {if(!silent) print("Finished fitting LPFR models")}
      )
    } else {
      fits <- parallel::mclapply(1:L, fit_lpfr, mc.cores = num_cores)
    }
  } else {
    fits <- lapply(1:L, fit_lpfr)
  }

  ## Obtain pointwise estimates
  # scalar predictors
  betaHat_raw <- lapply(fits, function(x) x[[1]])
  betaHat_raw <- t(do.call(rbind, betaHat_raw))

  # functional predictors
  gamma_mats_raw <- list()
  for (i in 1:p.ff){
    gamma_list <- lapply(fits, function(x) x[[2]][[i]]) # functional predictor
    gamma_mats_raw[[i]] <- do.call(cbind, gamma_list) # get raw estimate by column binding vectors
  }

  if (analytic == TRUE){
    beta.star_raw <- lapply(fits, function(x) x[[6]])
    beta.star_raw <- (do.call(rbind, beta.star_raw))

    ## Get scalar and functional design matrices
    X.scalar <- X.star[,1:(p.scalar+1)]
    beta.scalar_raw <- beta.star_raw[,1:(p.scalar+1)]
    X.functional <- X.star[,-(1:(p.scalar+1))]
    beta.functional_raw <- beta.star_raw[,-(1:(p.scalar+1))]
  }
  #############################################
  ##### Step 2: Smoothing
  #############################################
  if (!silent){cat("Step 2: Smoothing\n")}
  ## Smoothing scalar predictors and getting smoother parameters
  betaHat <- matrix(NA, nrow = p.scalar+1, ncol = L)
  lambda.scalar <- rep(NA, p.scalar+1) # smoothing parameter
  for(r in 1:(p.scalar+1)){
    #fit_smooth <- mgcv::gam(betaHat_raw[r,] ~ s(argvals, bs = "cr", k = (nknots_scalar + 1)), method = "REML")
    fit_smooth <- mgcv::gam(formula = betaHat_raw[r,] ~ s(argvals, bs = "cr", k = (nknots_scalar + 1)), method = "REML")
    betaHat[r,] <- fit_smooth$fitted.values
    lambda.scalar[r] <- fit_smooth$sp ## get smoothing parameter
  }
  if (analytic == TRUE){ ## get smoothing matrices
    ## scalar predictors
    sm <- mgcv::smoothCon(s(argvals, bs = "cr", k = (nknots_scalar + 1)),
                    data=data.frame(argvals=argvals), absorb.cons=TRUE)
    S <- sm[[1]]$S[[1]] ## penalty matrix
    B <- sm[[1]]$X ## basis functions

    ## functional predictors
    gamma_S.list <- list()
    GammaHat.list <- list()
    for (i in 1:p.ff){
      # fits assume argvals are 1...U; if the interval is smaller (say [0,1])
      # we need to adjust the estimates accordingly by multiplying by 1/stepsize
      smooth_GammaHat <- refund::fbps(gamma_mats_raw[[i]] * 1/ff.stepsize[[i]],
                              knots = nknots_fbps, p = bspline_p, m = diff_penalty)
      GammaHat.list[[i]] <- smooth_GammaHat$Yhat

      ## Get functional predictor smoother matrix
      gamma_lambda <- smooth_GammaHat$lambda
      gamma_smoothers <- get_gamma_smoothers(gamma_mats_raw[[i]], nknots_fbps, bspline_p, diff_penalty)

      B1 <- gamma_smoothers$B1 # basis functions
      P1 <- gamma_smoothers$P1 # penalty matrix
      B2 <- gamma_smoothers$B2
      P2 <- gamma_smoothers$P2

      S1 <- B1 %*% solve((t(B1) %*% B1 + gamma_lambda[1] * P1)) %*% t(B1)
      S2 <- B2 %*% solve((t(B2) %*% B2 + gamma_lambda[2] * P2)) %*% t(B2)

      gamma_S.list[[i]] <- kronecker(S2, S1) # tensor product of univariate smoothers
    }
  } else{
    GammaHat.list <- list()
    for (i in 1:p.ff){
      smooth_GammaHat <- refund::fbps(gamma_mats_raw[[i]] * U.list[[i]],
                              knots = nknots_fbps, p = bspline_p, m = diff_penalty)
      GammaHat.list[[i]] <- smooth_GammaHat$Yhat
    }
  }

  #############################################
  ##### Step 3: Inference
  #############################################
  if (var ==  TRUE){
    if (analytic == TRUE){ ### Analytic inference (Gaussian response)
      if (!silent){cat("Step 3: Analytic Inference\n")}
      ## Collect variance estimates into rows
      sigma2_errs <- t(as.matrix(unlist(lapply(fits, function(x) x[[3]]))))
      sigma2_subjs <- t(as.matrix(unlist(lapply(fits, function(x) x[[4]]))))

      ### 1. Smooth variance estimates across time points
      HHat <- t(apply(sigma2_subjs, 1, function(b) smooth.spline(x = argvals, y = b)$y))
      HHat[which(HHat < 0)] <- 0
      RHat <- t(apply(sigma2_errs, 1, function(b) smooth.spline(x = argvals, y = b)$y))
      RHat[which(RHat < 0)] <- 0

      ### 2. Estimate covariance between random effects G(s1, s2)
      GHat.scalar <- matrix(NA, L, L)
      GHat.functional <- matrix(NA, L, L)

      ## Ghat for scalar predictors
      for(i in 1:L){
        for(j in i:L){
          GHat.scalar[i,j] = GHat.scalar[j,i] = cov(data$Y[,i], data$Y[,j], use = "pairwise.complete.obs") -
            t(matrix(betaHat[,i], ncol = 1)) %*% var(X.scalar) %*% matrix(betaHat[,j], ncol = 1) -
            t(matrix(beta.functional_raw[i,], ncol = 1)) %*% var(X.functional) %*% matrix(beta.functional_raw[j,], ncol = 1)
        }
      }
      diag(GHat.scalar) <- HHat # Set diagonal as smoothed variance estimates
      GHat.scalar_smooth <- refund::fbps(GHat.scalar)$Yhat # Smooth GHat.scalar to reduce variability
      diag(GHat.scalar_smooth)[which(diag(GHat.scalar_smooth) < 0)] <- diag(GHat.scalar)[which(diag(GHat.scalar_smooth) < 0)]

      ## GHat for functional predictors
      for(i in 1:L){
        for(j in i:L){
          GHat.functional[i,j] = GHat.functional[j,i] = cov(data$Y[,i], data$Y[,j], use = "pairwise.complete.obs")
        }
      }
      diag(GHat.functional) <- HHat # Set diagonal as smoothed variance estimates
      GHat.functional_smooth <- refund::fbps(GHat.functional)$Yhat # Smooth GHat.functional to reduce variability
      diag(GHat.functional_smooth)[which(diag(GHat.functional_smooth) < 0)] <- diag(GHat.functional)[which(diag(GHat.functional_smooth) < 0)]

      ### 3.1 Get intra-location variance estimates Var(\hat \beta(s_l))
      if (!silent){cat("Step 3.1: Estimating Intra-location Variance\n")}
      intra_vars.scalar <- array(NA, dim = c(p.scalar+1, p.scalar+1, L))
      intra_vars.functional <- array(NA, dim = c(p.ff*kb, p.ff*kb, L))

      ### Compute variance-covariance matrix at each location of the response domain
      ## Scalar
      for (i in 1:L){
        V.mat <- Z %*% diag(HHat[[i]], I) %*% t(Z) + diag(RHat[[i]], nrow(data))
        V.inv <- solve(V.mat)
        intra_vars.scalar[,,i] <- solve(t(X.scalar) %*% V.inv %*% X.scalar)
      }

      ## Functional
      lambdas <- lapply(fits, function(x) x[[5]])

      for (i in 1:L){
        # setup penalty matrix
        blocks <- list()
        for (j in 1:p.ff){
          block.ff <- list(matrix(0, 2, 2), as.matrix(lambdas[[i]][[j]]*diag(kb-2)))
          blocks <- c(blocks, block.ff)
        }

        D.lambda <- as.matrix(Matrix::bdiag(blocks))

        V.mat <- Z %*% diag(HHat[[i]], I) %*% t(Z) + diag(RHat[[i]], nrow(data))
        V.inv <- solve(V.mat)
        mat <- solve(t(X.functional) %*% V.inv %*% X.functional + D.lambda) %*% t(X.functional)
        intra_vars.functional[,,i] <- mat %*% V.inv %*% t(mat)
      }

      #### 4. Get inter-location covariance estimates Cov(\hat \beta(s_1), \hat \beta(s_2))
      if (!silent){cat("Step 3.2: Estimating Inter-location Variance\n")}
      cov.all.scalar <- array(NA, dim = c(p.scalar+1, p.scalar+1, L, L))
      cov.all.functional <- array(NA, dim = c(p.ff*kb, p.ff*kb, L, L))

      ## Scalar
      if(parallel == TRUE){
        get_cov_all_parallel.scalar <- function(i,j){
          Gmat <- diag(GHat.scalar_smooth[i,j], I)

          V.mat.i <- Z %*% diag(HHat[[i]], I) %*% t(Z) + diag(RHat[[i]], nrow(data))
          V.inv.i <- solve(V.mat.i)

          V.mat.j <- Z %*% diag(HHat[[j]], I) %*% t(Z) + diag(RHat[[j]], nrow(data))
          V.inv.j <- solve(V.mat.j)

          mat.i <- solve(t(X.scalar) %*% V.inv.i %*% X.scalar) %*% t(X.scalar) %*% V.inv.i
          mat.j <- solve(t(X.scalar) %*% V.inv.j %*% X.scalar) %*% t(X.scalar) %*% V.inv.j

          cov.val <- mat.i %*% Z %*% Gmat %*% t(Z) %*% t(mat.j)

          return(cov.val)
        }
        tmp <- list()
        for(i in 1:L){
          if(.Platform$OS.type == "windows"){
            # clear old variables
            parallel::clusterEvalQ(cl, rm(list = ls()))

            # add global variables to clusters
            parallel::clusterExport(cl, varlist = c("GHat.scalar_smooth", "I", "Z", "HHat",
                                          "RHat", "data", "X.scalar"),
                          envir = environment())
            tmp[[i]] <- tryCatch(
              parallel::parLapply(cl = cl, X = i:L,
                        fun = get_cov_all_parallel.scalar,
                        i = i),
              warning = function(w){
                print(paste0("Warning during model fitting:", "\n", w))
              },
              error = function(e){
                stop(paste0("Error during model fitting:", "\n", e))
              }
            )
          }else{
            tmp[[i]] <- parallel::mclapply(i:L, get_cov_all_parallel.scalar, i = i, mc.cores = num_cores)
          }
          for(j in 1:length(tmp[[i]])){
            j_idx <- i - 1 + j
            cov.all.scalar[,,i,j_idx] = cov.all.scalar[,,j_idx,i] = tmp[[i]][[j]]
          }
        }
      }else{
        for(i in 1:L){
          for(j in i:L){
            Gmat <- diag(GHat.scalar_smooth[i,j], I)

            V.mat.i <- Z %*% diag(HHat[[i]], I) %*% t(Z) + diag(RHat[[i]], nrow(data))
            V.inv.i <- solve(V.mat.i)

            V.mat.j <- Z %*% diag(HHat[[j]], I) %*% t(Z) + diag(RHat[[j]], nrow(data))
            V.inv.j <- solve(V.mat.j)

            mat.i <- solve(t(X.scalar) %*% V.inv.i %*% X.scalar) %*% t(X.scalar) %*% V.inv.i
            mat.j <- solve(t(X.scalar) %*% V.inv.j %*% X.scalar) %*% t(X.scalar) %*% V.inv.j

            cov.all.scalar[,,i,j] = cov.all.scalar[,,j,i] = mat.i %*% Z %*% Gmat %*% t(Z) %*% t(mat.j)
          }
        }
      }

      ## Functional
      if(parallel == TRUE){
        get_cov_all_parallel.functional <- function(i,j){
          # setup penalty matrices
          blocks <- list()
          for (k in 1:p.ff){
            block.ff <- list(matrix(0, 2, 2), as.matrix(lambdas[[i]][[k]]*diag(kb-2)))
            blocks <- c(blocks, block.ff)
          }
          D.lambda.i <- as.matrix(Matrix::bdiag(blocks))

          blocks <- list()
          for (k in 1:p.ff){
            block.ff <- list(matrix(0, 2, 2), as.matrix(lambdas[[j]][[k]]*diag(kb-2)))
            blocks <- c(blocks, block.ff)
          }
          D.lambda.j <- as.matrix(Matrix::bdiag(blocks))

          # compute covariance
          Gmat <- diag(GHat.functional_smooth[i,j], I)

          V.mat.i <- Z %*% diag(HHat[[i]], I) %*% t(Z) + diag(RHat[[i]], nrow(data))
          V.inv.i <- solve(V.mat.i)

          V.mat.j <- Z %*% diag(HHat[[j]], I) %*% t(Z) + diag(RHat[[j]], nrow(data))
          V.inv.j <- solve(V.mat.j)

          mat.i <- solve(t(X.functional) %*% V.inv.i %*% X.functional + D.lambda.i) %*% t(X.functional) %*% V.inv.i
          mat.j <- solve(t(X.functional) %*% V.inv.j %*% X.functional + D.lambda.j) %*% t(X.functional) %*% V.inv.j

          cov.val <- mat.i %*% Z %*% Gmat %*% t(Z) %*% t(mat.j)

          return(cov.val)
        }
        tmp <- list()
        for(i in 1:L){
          if(.Platform$OS.type == "windows"){
            # clear old variables
            parallel::clusterEvalQ(cl, rm(list = ls()))

            # add necessary packages to clusters
            parallel::clusterEvalQ(cl, library(Matrix))

            # add global variables to clusters
            parallel::clusterExport(cl, varlist = c("p.ff", "lambdas", "kb", "GHat.functional_smooth",
                                          "I", "Z", "HHat", "RHat", "data", "X.functional"),
                          envir = environment())
            tmp[[i]] <- tryCatch(
              parallel::parLapply(cl = cl, X = i:L,
                        fun = get_cov_all_parallel.functional,
                        i = i),
              warning = function(w){
                print(paste0("Warning during model fitting:", "\n", w))
              },
              error = function(e){
                stop(paste0("Error during model fitting:", "\n", e))
              }
            )
          } else{
            tmp[[i]] <- parallel::mclapply(i:L, get_cov_all_parallel.functional, i = i, mc.cores = num_cores)
          }
          for(j in 1:length(tmp[[i]])){
            j_idx <- i - 1 + j
            cov.all.functional[,,i,j_idx] = cov.all.functional[,,j_idx,i] = tmp[[i]][[j]]
          }
        }
      }else{
        for(i in 1:L){
          for(j in i:L){
            # setup penalty matrices
            blocks <- list()
            for (k in 1:p.ff){
              block.ff <- list(matrix(0, 2, 2), as.matrix(lambdas[[i]][[k]]*diag(kb-2)))
              blocks <- c(blocks, block.ff)
            }
            D.lambda.i <- as.matrix(Matrix::bdiag(blocks))

            blocks <- list()
            for (k in 1:p.ff){
              block.ff <- list(matrix(0, 2, 2), as.matrix(lambdas[[j]][[k]]*diag(kb-2)))
              blocks <- c(blocks, block.ff)
            }
            D.lambda.j <- as.matrix(Matrix::bdiag(blocks))

            # compute covariance
            Gmat <- diag(GHat.functional_smooth[i,j], I)

            V.mat.i <- Z %*% diag(HHat[[i]], I) %*% t(Z) + diag(RHat[[i]], nrow(data))
            V.inv.i <- solve(V.mat.i)

            V.mat.j <- Z %*% diag(HHat[[j]], I) %*% t(Z) + diag(RHat[[j]], nrow(data))
            V.inv.j <- solve(V.mat.j)

            mat.i <- solve(t(X.functional) %*% V.inv.i %*% X.functional + D.lambda.i) %*% t(X.functional) %*% V.inv.i
            mat.j <- solve(t(X.functional) %*% V.inv.j %*% X.functional + D.lambda.j) %*% t(X.functional) %*% V.inv.j

            cov.all.functional[,,i,j] = cov.all.functional[,,j,i] = mat.i %*% Z %*% Gmat %*% t(Z) %*% t(mat.j)
          }
        }
      }

      ### 5. Construct covariance matrices for scalar + functional predictors
      if (!silent){cat("Step 3.3: Finalizing Variance Estimation\n")}
      ## Scalar Predictors
      if (!silent){cat("Step 3.3.1: Scalar Predictor Variance Estimation\n")}
      cov.scalar <- replicate(p.scalar+1, matrix(NA, nrow = L, ncol = L), simplify = FALSE)
      cov.scalar.trimmed = list()
      for (i in 1:L){
        for (j in i:L){
          if(i == j){
            for(k in 1:(p.scalar+1)){
              cov.scalar[[k]][i,j] <- intra_vars.scalar[k,k,i]
            }
          } else{
            for(k in 1:(p.scalar+1)){
              cov.scalar[[k]][i,j] = cov.scalar[[k]][j,i] = cov.all.scalar[k,k,i,j]
            }
          }
        }
      }

      # smooth covariance matrices
      for (i in 1:(p.scalar+1)){
        M <- B %*% solve(t(B) %*% B + lambda.scalar[i]*S) %*% t(B) + matrix(1/L, nrow = L, ncol = L)
        cov.scalar.smoothed <- M %*% cov.scalar[[i]] %*% t(M)

        # trim non-positive eigenvalues
        edcomp <- eigen(cov.scalar.smoothed)
        eigen.positive <- which(edcomp$values > 0)
        if(length(eigen.positive) == L){
          cov.scalar.trimmed[[i]] <- cov.int.smoothed
        } else{
          cov.scalar.trimmed[[i]] <- edcomp$vectors[,eigen.positive] %*% diag(edcomp$values[eigen.positive]) %*% t(edcomp$vectors[,eigen.positive])
        }
      }

      ## Functional Predictors
      if (!silent){cat(paste0("Step 3.3.2: Functional Predictor Variance Estimation", "\n",
                              "Note, this can take some time for dense functional data\n"))}
      phi.list <- list()
      cov.functional <- list()
      for (i in 1:p.ff){
        t <- seq(0, 1, length = U.list[[i]])
        num=kb-2
        qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
        knots <- quantile(t, qtiles)
        # Basis functions
        phi.list[[i]] <- cbind(1, t, sapply(knots, function(k) ((t - k > 0) * (t - k))))

        ## Get covariance matrix for functional predictor
        cov.functional[[i]] <-  array(0, dim = c(L*U.list[[i]], L*U.list[[i]]))
      }

      get_cov_functional <- function(s,u,t,v,phi,ff.idx){
        if (s == u){
          cov.tmp <- intra_vars.functional[((ff.idx - 1)*kb+1):(ff.idx*kb),
                                           ((ff.idx - 1)*kb+1):(ff.idx*kb),
                                           s]
          cov.tmp.functional <- t(as.matrix(phi[t,])) %*% cov.tmp %*% as.matrix(phi[v,])
        }else{
          cov.tmp <- cov.all.functional[((ff.idx - 1)*kb+1):(ff.idx*kb),
                                        ((ff.idx - 1)*kb+1):(ff.idx*kb),
                                        s,u]
          cov.tmp.functional <- t(as.matrix(phi[t,])) %*% cov.tmp %*% as.matrix(phi[v,])
        }
        return(cov.tmp.functional)
      }

      if(parallel == TRUE){
        covBeta.parallel <- function(s.par, u.par, U, phi, ff.idx){
          cov.beta <- array(0, dim = c(U, U))
          for(t in 1:U){
            for(v in t:U){
              tmp <- get_cov_functional(s.par,u.par,t,v,phi,ff.idx)
              cov.beta[t, v] <- cov.beta[v, t] <- tmp
            }
          }
          return(cov.beta)
        }
        for (i in 1:p.ff){
          tmp <- list()
          for(s in argvals){
            if(.Platform$OS.type == "windows"){
              # clear old variables
              parallel::clusterEvalQ(cl, rm(list = ls()))

              # add global variables to clusters
              parallel::clusterExport(cl, varlist = c("get_cov_functional", "intra_vars.functional",
                                            "cov.all.functional", "kb"),
                            envir = environment())

              tmp[[s]] <- tryCatch(
                parallel::parLapply(cl = cl, X = s:L,
                          fun = covBeta.parallel,
                          s.par=s, U=U.list[[i]],
                          phi=phi.list[[i]], ff.idx=i),
                warning = function(w){
                  print(paste0("Warning during model fitting:", "\n", w))
                },
                error = function(e){
                  stop(paste0("Error during model fitting:", "\n", e))
                }
              )
            } else{
              tmp[[s]] <- parallel::mclapply(s:L, covBeta.parallel, s.par=s, U=U.list[[i]],
                                             phi=phi.list[[i]], ff.idx=i, mc.cores = num_cores)
            }
            for(idx in 1:length(tmp[[s]])){
              u <- s - 1 + idx
              cov.functional[[i]][((s-1)*U.list[[i]]+1):(s*U.list[[i]]), ((u-1)*U.list[[i]]+1):(u*U.list[[i]])] <-
                cov.functional[[i]][((u-1)*U.list[[i]]+1):(u*U.list[[i]]), ((s-1)*U.list[[i]]+1):(s*U.list[[i]])] <- tmp[[s]][[idx]]
            }
          }
        }
      } else{
        for(s in 1:L){
          for(u in s:L){
            for (i in 1:p.ff){
              for(t in 1:U.list[[i]]){
                for(v in t:U.list[[i]]){
                  tmp <- get_cov_functional(s,u,t,v,phi.list[[i]],i)
                  cov.functional[[i]][(u-1)*U.list[[i]]+v, (s-1)*U.list[[i]]+t] <-
                    cov.functional[[i]][(u-1)*U.list[[i]]+t, (s-1)*U.list[[i]]+v] <-
                    cov.functional[[i]][(s-1)*U.list[[i]]+v, (u-1)*U.list[[i]]+t] <-
                    cov.functional[[i]][(s-1)*U.list[[i]]+t, (u-1)*U.list[[i]]+v] <- tmp
                }
              }
            }
          }
        }
      }

      ## Smooth functional covariance matrix
      cov.functional.trimmed <- list()
      for (i in 1:p.ff){
        cov.functional[[i]] <- cov.functional[[i]] * ((1/ff.stepsize[[i]])^2) # estimate needs to be multiplied by length of domain U in our simulation setting
        cov.functional.smoothed <- gamma_S.list[[i]] %*% cov.functional[[i]] %*% t(gamma_S.list[[i]])

        # trim negative eigenvalues
        edcomp <- eigen(cov.functional.smoothed)
        eigen.positive <- which(edcomp$values > 0)
        if(length(eigen.positive) == L*U.list[[i]]){
          cov.functional.trimmed[[i]] <- cov.functional.smoothed
        }else{
          cov.functional.trimmed[[i]] <- edcomp$vectors[,eigen.positive] %*% diag(edcomp$values[eigen.positive]) %*% t(edcomp$vectors[,eigen.positive])
        }
      }

      ### 6. Organize final variance estimates
      if (!silent){cat("Step 3.4: Finalizing Variance Estimation\n")}
      ## Intercept + scalar predictors
      var.scalar <- list()
      for(i in 1:(p.scalar+1)){
        var.scalar[[i]] <- diag(cov.scalar[[i]])
      }

      ## Functional predictors
      var.ff <- list()
      for (i in 1:p.ff){
        var.ff[[i]] <- matrix(NA, U.list[[i]], L)
        vars.tmp <- diag(cov.functional.trimmed[[i]])
        for(s in 1:L){
          var.ff[[i]][,s] <- vars.tmp[((s-1)*U.list[[i]] + 1):((s-1)*U.list[[i]] + U.list[[i]])]
        }
      }

      ### 7. CMA quantiles
      if (CMA == TRUE){
        qn.scalar <- rep(0, length = p.scalar + 1)
        qn.functional <- rep(0, length = p.ff)
        N <- 10000 ## sample size in simulation-based approach

        # Scalar predictors
        for(i in 1:(p.scalar+1)){
          # ensure PD matrix
          cov.pd <- as.matrix(Matrix::nearPD(cov.scalar.trimmed[[i]], ensureSymmetry = TRUE)$mat)
          samples <- mvtnorm::rmvnorm(N, mean = betaHat[i,], sigma = cov.pd)
          un <- rep(NA, N)
          for(j in 1:N){
            un[j] <- max(abs((samples[j,] - betaHat[i,])/sqrt(diag(cov.pd))))
          }
          qn.scalar[[i]] <- quantile(un, 0.95)
        }

        ## Functional predictors
        for(i in 1:p.ff){
          cov.pd <- as.matrix(Matrix::nearPD(cov.functional.trimmed[[i]],
                                     ensureSymmetry = TRUE)$mat)

          gammaVec <- ks::vec(GammaHat.list[[i]])
          samples <- mvtnorm::rmvnorm(N, mean = gammaVec, sigma = cov.pd)
          un <- rep(NA, N)
          for(j in 1:N){
            un[j] <- max(abs((samples[j,] - gammaVec)/sqrt(diag(cov.pd))))
          }
          qn.functional[[i]] <- quantile(un, 0.95)
        }
        qn <- list(qn.scalar = qn.scalar, qn.functional = qn.functional)
      }else {
        qn <- NULL
      }
    }
    else{
      ############################
      #### Bootstrap inference
      ############################
      ## Run bootstrapped models
      if (!silent){cat("Step 3: Bootstrap Inference\n")}
      if(!(parallel & .Platform$OS.type == "windows")){cl <- NA}
      boot_fits <- run_boot(formula = formula, data = data,
                            subj_var = subj_var, visit_var = visit_var,
                            seed = seed, p = p, argvals = argvals,
                            U = U, num_boot = num_boot,
                            splines = "ps", nknots_scalar = nknots_scalar,
                            nknots_fbps = nknots_scalar,
                            bspline_p = bspline_p, diff_penalty = diff_penalty,
                            scenario = "default", num_cores = num_cores,
                            kz = kz, kb = kb, smooth.cov = smooth.cov,
                            family = family, cl = cl)

      #### Get variance estimates from bootstrap
      ### Scalar covariates
      # For each time point, get bootstrap estimates of scalar covariates
      res <- lapply(boot_fits[2], function(x) lapply(x, function(y) y[[1]]))
      res <- do.call(c, res)

      # Create a 3D array to hold all results
      res_array <- array(unlist(res), dim = c(p.scalar+1, L, num_boot))

      # Function to variance at each time point across bootstrapped samples
      get_var <- function(x) {
        apply(x, c(1, 2), function(y) var(y))
      }

      # Apply the variance function across the third dimension
      var.est <- get_var(res_array)
      var.scalar <- split(var.est, seq_len(p.scalar+1))

      ### Functional predictors
      # Getting CI by bootstrapped samples
      res.ff <- lapply(boot_fits[2], function(x) lapply(x, function(y) (y[[2]])))
      res.ff <- do.call(c, res.ff)
      var.ff <- list()

      for (i in 1:p.ff){
        # Combine the matrices into a 3D array
        res.tmp <- lapply(res.ff, function(x) x[[i]])
        res.ff_3d <- array(unlist(res.tmp), dim = c(U.list[[i]], L, num_boot))

        # Calculate the variance along the third dimension
        var.ff[[i]] <- apply(res.ff_3d, c(1, 2), var)
      }

      ### Get quantiles for CMA inference
      ## Step 1: get max order statistics across bootstrapped samples
      ## Scalar predictors
      d.b <- matrix(NA, nrow = length(res), ncol = p.scalar+1)

      for (i in 1:length(res)){
        d.vec <- abs(res[[i]] - betaHat)/sqrt(var.est)
        d.b[i,] <- apply(d.vec, 1, max)
      }

      ## Functional predictors
      db_ff <- matrix(NA, nrow = length(res), ncol = p.ff)

      for (i in 1:p.ff){
        for (j in 1:length(res.ff)){
          d.mat <- abs(res.ff[[j]][[i]] - GammaHat.list[[i]])/sqrt(var.ff[[i]])
          db_ff[j,i] <- max(d.mat)
        }
      }

      ## Step 2: Obtain empirical quantiles
      #TODO: ff cma quantiles seem large in initial tests
      cma.quantiles.scalar <- apply(d.b, 2, quantile, probs = 0.95)
      names(cma.quantiles.scalar) <- c("Intercept", scalar_varnames)
      cma.quantiles.functional <- apply(db_ff, 2, quantile, probs = 0.95)
      names(cma.quantiles.functional) <- ff_varnames
    }
  }
  time <- (proc.time() - ptm)[3]
  ## add names to return objects
  scalar.est <- split(betaHat, seq_len(p.scalar+1))
  names(scalar.est) <- c("Intercept", scalar_varnames)
  names(GammaHat.list) <- ff_varnames
  if (var == TRUE){
    names(var.scalar) <- c("Intercept", scalar_varnames)
    names(var.ff) <- ff_varnames
  }

  if(var ==  TRUE){
    if(analytic == TRUE){
      return(list(scalar.est = scalar.est, functional.est = GammaHat.list,
                  var.scalar = var.scalar, var.functional = var.ff,
                  time = time, CMA_quantiles = qn))
    } else{
      return(list(scalar.est = scalar.est, functional.est = GammaHat.list,
                  var.scalar = var.scalar, var.functional = var.ff,
                  time  = time,
                  CMA_quantiles = c(cma.quantiles.scalar, cma.quantiles.functional)))
    }
  } else{
    return(list(scalar.est = scalar.est,
                functional.est = GammaHat.list, time  = time))
  }
}
