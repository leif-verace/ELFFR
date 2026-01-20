#' Helper function to fit bootstrapped ELFFR models
#' @importFrom mgcv s
run_boot <- function(formula, data, subj_var, visit_var, seed, p, argvals, U.list, num_boot = 50,
                     splines = "ps", nknots_scalar = 5,
                     nknots_fbps = 3, bspline_p = 3, diff_penalty = 2,
                     scenario = "default", num_cores = 8,
                     kz = 30, kb = 30, smooth.cov = TRUE,
                     family = "gaussian", cl = NA){
  ptm <- proc.time()

  ## Organize objects specified in formula for use in `lpfr`
  model.char <- as.character(formula)
  Y_label <- model.char[[2]] # outcome
  varnames <- unlist(strsplit(model.char[[3]], "\\s*\\+\\s*"))

  # get functional predictor variable names
  ff_indices <- grep("^ff\\(.*?\\)$", varnames)
  ff_varnames <- sub("^ff\\((.*?)\\)$", "\\1", varnames[ff_indices])

  # get scalar predictor variable names
  scalar_varnames <- varnames[-ff_indices]

  p.ff <- length(ff_varnames)
  p.scalar <- length(scalar_varnames)

  ####################################################
  ## Create temporary directories to store temp data
  ####################################################
  # Get the path of the main temporary directory
  main_temp_dir <- tempdir()

  # Define paths for the two new temporary directories
  temp_dir.main <- file.path(main_temp_dir, "temp_dir.main")
  temp_dir.sub <- file.path(main_temp_dir, "temp_dir.sub")

  # Create the new directories
  dir.create(temp_dir.main, showWarnings = FALSE, recursive = TRUE)
  dir.create(temp_dir.sub, showWarnings = FALSE, recursive = TRUE)

  #################################################
  #### Set parameters and run bootstrap B times
  #################################################
  # Define bootstrap interval
  id_vec <- unique(data[[subj_var]], na.rm = TRUE)
  B_sample_size <- length(id_vec) # create bootstrap samples with same size as original data
  pb <- txtProgressBar(min = 0,
                       max = num_boot, # Maximum value
                       style = 3,
                       width = 50,
                       char = "=")   # Character to use for the bar

  for (b in 1:num_boot){
    #print(paste("Bootstrap: ", b, sep = ""))
    boot_sample <- sample(id_vec, B_sample_size, replace = TRUE) # sampled ID list

    boot_dat <- tibble::tibble(!!subj_var := boot_sample) |>
      dplyr::left_join(
        data,
        by = subj_var,
        relationship = "many-to-many"
      ) |>
      dplyr::group_by(.data[[subj_var]], .data[[visit_var]]) |>
      dplyr::mutate(group_num = cumsum(.data[[visit_var]] == dplyr::first(.data[[visit_var]]))) |>
      dplyr::ungroup() |>
      dplyr::mutate(new_id = paste0(.data[[subj_var]], group_num)) |>
      dplyr::select(-group_num)

    # Convert id to numeric so the random intercept works in mgcv
    boot_dat[[subj_var]] <- as.numeric(boot_dat[[subj_var]])
    boot_dat$new_id <- as.numeric(boot_dat$new_id)

    # Convert variables to factors
    boot_dat <- boot_dat |>
      dplyr::mutate(across(all_of(visit_var), factor))

    #############################################################
    #### Build model using marginal approach
    if (parallel == TRUE){
      if(.Platform$OS.type == "windows"){
        # clear old variables
        parallel::clusterEvalQ(cl, rm(list = ls()))

        # add necessary packages to clusters
        parallel::clusterEvalQ(cl, library(refund))

        # add global variables to clusters
        parallel::clusterExport(cl, varlist = c("subj_var"),
                      envir = environment())
        tryCatch(
          parallel::parLapply(cl = cl, X = 1:L,
                    fun = lpfr_boot,
                    b = b,
                    boot_dat = boot_dat, temp_dir = temp_dir.sub,
                    kz = kz, kb = kb, smooth.cov = smooth.cov,
                    family = family, scalar_varnames = scalar_varnames,
                    ff_varnames = ff_varnames, Y_label = Y_label),
          warning = function(w){
            print(paste0("Warning during model fitting:", "\n", w))
          },
          error = function(e){
            stop(paste0("Error during model fitting:", "\n", e))
          }
        )
      } else{
        parallel::mclapply(1:L, lpfr_boot, mc.cores = num_cores, b = b,
                           boot_dat = boot_dat, temp_dir = temp_dir.sub,
                           kz = kz, kb = kb, smooth.cov = smooth.cov,
                           family = family, scalar_varnames = scalar_varnames,
                           ff_varnames = ff_varnames, Y_label = Y_label)
      }
    } else{
      lapply(1:L, lpfr_boot, b = b,
             boot_dat = boot_dat, temp_dir = temp_dir.sub,
             kz = kz, kb = kb, smooth.cov = smooth.cov,
             family = family)
    }

    ############################################################
    #### Regression smoothing
    ## Load data from RDS files
    fits <- list.files(temp_dir.sub, full.names=TRUE) |> purrr::map(readRDS)

    ## Sort so time points are sequential 1, .. , L (not read sequentially by list.files)
    indices <- sapply(fits, function(x) x[[1]])
    order_indices <- order(indices)
    sorted_fits <- fits[order_indices]

    ##################################################
    #### Smooth Scalar Regression Coefficients
    # Set parameters, initialize objects to store values
    betaHat_raw <- lapply(sorted_fits, function(x) x[[2]])
    betaHat_raw <- t(do.call(rbind, betaHat_raw))

    betaHat <- matrix(NA, nrow = p.scalar+1, ncol = L)
    lambda.scalar <- rep(NA, p.scalar+1) ## smoothing parameter
    for(r in 1:(p.scalar+1)){
      fit_smooth <- mgcv::gam(formula = betaHat_raw[r,] ~ s(argvals, bs = "cr", k = (nknots_scalar + 1)), method = "REML")
      betaHat[r,] <- fit_smooth$fitted.values
      lambda.scalar[r] <- fit_smooth$sp ## get smoothing parameter
    }

    #### Smooth Functional Coefficient Estimates
    ## Get matrix of γ(s, u) estimated coefficients
    GammaHat.list <- list()
    for (i in 1:p.ff){
      gamma_list <- lapply(fits, function(x) x[[3]][[i]]) # functional predictor
      gamma_mat_raw <- do.call(cbind, gamma_list) # get raw estimate by column binding vectors
      # use sandwich smoother on gamma matrix
      smooth_GammaHat <- refund::fbps(gamma_mat_raw, knots = nknots_fbps, p = bspline_p, m = diff_penalty)
      GammaHat.list[[i]] <- smooth_GammaHat$Yhat * U.list[[i]]
    }

    saveRDS(list(scalar.est = betaHat, functional.est = GammaHat.list, b = b),
            file = paste(temp_dir.main, "/smoothed_fits_boot_", L, "_", b, ".RDS", sep=""))

    # Remove files from temp sub directory
    remove.sub <- list.files(temp_dir.sub, full.names = TRUE, recursive = TRUE)

    # Remove all files and directories
    unlink(remove.sub, recursive = TRUE)
    setTxtProgressBar(pb, b)
  }
  close(pb)
  ## Load data from RDS files
  fits_boot <- list.files(temp_dir.main, full.names=TRUE) |> purrr:::map(readRDS)

  ## Sort boot outputs
  indices_boot <- sapply(fits_boot, function(x) x[[3]])
  order_indices_boot <- order(indices_boot)
  sorted_fits_boot <- fits_boot[order_indices_boot]

  unlink(temp_dir.main, recursive = TRUE)

  return(list(time = (proc.time() - ptm)[3], res = sorted_fits_boot))
}
