#' Simulate example data for lfofr - adapted from pffrSim in refund
#'
#' @param scenario see Description
#' @param n number of observations
#' @param nxgrid number of evaluation points of functional covariates
#' @param nygrid number of evaluation points of the functional response
#' @param SNR the signal-to-noise ratio for the generated data: empirical
#'   variance of the additive predictor divided by variance of the errors.
#' @param propmissing proportion of missing data in the response, default = 0.
#'   See Details.
#' @param limits a function that defines an integration range, see
#'   \code{\link{ff}}
#' @importFrom splines spline.des
#' @importFrom stats dnorm rnorm
#' @return a named list with the simulated data, and the true components of the
#'   predictor etc as attributes.
pffrSim_modified <- function(
    func_predictor,
    scenario="all",
    n = 100,
    nxgrid = 40,
    nygrid = 60,
    SNR = 10,
    propmissing=0,
    limits = NULL){

  mc <- match.call()
  for(i in 2:length(mc)) if(is.symbol(mc[[i]]))
    mc[[i]] <- get(deparse(mc[[i]]), envir=parent.frame())

  ## generates random functions...
  rf <- function(x=seq(0,1,length=100), bs.dim=7, center=FALSE) {
    nk <- bs.dim - 2
    xu <- max(x)
    xl <- min(x)
    xr <- xu - xl
    xl <- xl - xr * 0.001
    xu <- xu + xr * 0.001
    dx <- (xu - xl)/(nk - 1)
    kn <- seq(xl - dx * 3, xu + dx * 3,
              length = nk + 4 + 2)
    X <- splines::spline.des(kn, x, 4, x * 0)$design

    drop(X %*% rnorm(bs.dim))
  }

  test1 <- function(s, t){
    s*cos(pi*abs(s-t)) - .19
  }
  test2 <- function(s, t, ss = 0.3, st = 0.4)
  {
    cos(pi * s) * sin(pi * t) + (s * t) ^ 2 - 0.11
  }
  test3 <- function(s, t, ss = 0.3, st = 0.4)
  {
    sin(2*pi*(1+s)*t)
  }
  s <- seq(0, 1, length=nxgrid)
  t <- seq(0, 1, length=nygrid)

  mu.t <- matrix(1 + dbeta(t, 2,7), nrow=n, ncol=nygrid, byrow=TRUE)

  data <- list()
  etaTerms <- list()
  etaTerms$int <- mu.t

  #functional covariates
  data$X1 <- I(t(replicate(n, rf(s))))

  L <- matrix(1/nxgrid, ncol=nxgrid, nrow=n)
  LX1 <- L*data$X1
  beta1.st <- outer(s, t, func_predictor)
  if(!is.null(limits)){
    range <- outer(s, t, limits)
    beta1.st  <- beta1.st * range
  }
  etaTerms$X1 <- LX1 %*% beta1.st

  #     data$X2 <- I(t(replicate(n, rf(s))))
  #     LX2 <- L*data$X2
  #     beta2.st <- outer(s, t, test2)
  #     etaTerms$X2 <- LX2%*%beta2.st

  #scalar covariates
  data$xlin <- I(rnorm(n))
  beta.t <- matrix(scale(-dnorm(4 * (t - .2))), nrow = n, ncol = nygrid,
                   byrow = T)
  etaTerms$xlin <- data$xlin * beta.t

  data$xsmoo <- I(rnorm(n))
  etaTerms$xsmoo <- outer(drop(scale(cos(data$xsmoo))), (t-.5), "*")

  data$xfactor <- sample(gl(3, n/3), replace = TRUE)
  etaTerms$xfactor <-  2* as.numeric(data$xfactor) +
    sin(2 * outer(as.numeric(data$xfactor), t))
  if ("2factor" %in% scenario) {
    data$x2factor <- sample(gl(3, n/3), replace = TRUE)
    etaTerms$x2factor <-  2 * as.numeric(data$x2factor) +
      cos(2 * outer(as.numeric(data$x2factor), t))
  }

  data$xte1 <- I(rnorm(n))
  data$xte2 <- I(rnorm(n))
  etaTerms$xte <- matrix(drop(scale(-data$xte1*data$xte2^2)),
                         ncol=nygrid,
                         nrow=n)

  data$xconst <- I(rnorm(n))
  etaTerms$xconst <- matrix(2*data$xconst,
                            ncol=nygrid,
                            nrow=n)

  if(length(scenario)==1){
    eta <- mu.t + switch(scenario,
                         "int" = 0,
                         "all" = Reduce("+", etaTerms),
                         "ff" =  etaTerms$X1, #+ etaTerms$X2,
                         "lin" = etaTerms$xlin,
                         "smoo" = etaTerms$xsmoo,
                         "te" = etaTerms$xte,
                         "const" = etaTerms$xconst,
                         "factor" = etaTerms$xfactor)
  } else {
    stopifnot(all(scenario %in% c("int" ,"ff", "lin", "smoo", "te", "const",
                                  "factor", "2factor")))
    eta <- 0*mu.t
    if("int" %in% scenario) eta <- eta + mu.t
    if("ff" %in% scenario) eta <- eta + etaTerms$X1 #+ etaTerms$X2
    if("lin" %in% scenario) eta <- eta + etaTerms$xlin
    if("smoo" %in% scenario) eta <- eta + etaTerms$xsmoo
    if("te" %in% scenario) eta <- eta + etaTerms$xte
    if("const" %in% scenario) eta <- eta + etaTerms$xconst
    if("factor" %in% scenario) eta <- eta + etaTerms$xfactor
    if("2factor" %in% scenario) eta <- eta + etaTerms$x2factor
  }


  eps <-  sd(as.vector(eta))/sqrt(SNR) * matrix(scale(rnorm(n*nygrid)),
                                                nrow=n)
  data$Y <- I(eta + eps)



  if(propmissing == 0){
    return(structure(as.data.frame(data, rownames=1:n), xindex=s, yindex=t,
                     truth=list(eta=eta, etaTerms=etaTerms), call=mc))
  } else {
    missing <- sample(c(rep(T, propmissing*n*nygrid),
                        rep(F, n*nygrid-propmissing*n*nygrid)))
    data <- as.data.frame(data, rownames=1:n)

    ydata <- data.frame(.obs = rep(1:n, each=nygrid)[!missing],
                        .index = rep(t, times=n)[!missing],
                        .value = as.vector(t(data$Y))[!missing])

    return(structure(list(data=data, ydata=ydata), xindex=s, yindex=t,
                     truth=list(eta=eta, etaTerms=etaTerms),
                     call=mc))
  }
}
