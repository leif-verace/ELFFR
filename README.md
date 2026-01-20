
# ELFFR R Package

<!-- badges: start -->
<!-- badges: end -->

ELFFR implements a marginal approach to fitting longitudinal function-on-function
regression models. This approach allows models to be fit quickly, even when the number
of subjects is large. The package currently only allows for a random intercept term for
random effects.

## Installation

The R package may be downloaded by running the following command within R or RStudio:

``` r
devtools::install_github("leif-verace/ELFFR")
```

## Package Usage

Below is an example demonstrating the syntax and usage of the main
`elffr` function on simulated data:

``` r
library(ELFFR)

## Set parameters for simulating data
func_predictors <- c("5*sin(0.5*pi*(s+0.5)^2)*cos(pi*u+0.5)",
                     "sin(u*pi/3)*cos(s*pi/2)")
U <- c(25, 30)

## Simulate data
data <- sim_elffr_data(func_predictors = func_predictors,
                       I = 100, L = 25,  U = U, J = 5,
                       SNR_B = 0.5, SNR_sigma = 1.5,
                       family = "gaussian")

## Run ELFFR on data
formula <- as.formula(Y ~ X1 + X2 + ff(W1) + ff(W2))
subj_var <- "ID"
visit_var <- "visit"
argvals.functional <- list(seq(0,1,length.out=U[[1]]), seq(0,1,length.out=U[[2]]))
nknots_fbps <- list(seq(0, 1, length = 10), seq(0, 1, length = 5))

res <- elffr(formula, subj_var, visit_var, data, num_cores = 8L, kz = 15, kb = 15,
             nknots_scalar = 8, nknots_fbps = nknots_fbps,
             analytic = TRUE, parallel = TRUE, var = TRUE, CMA = TRUE,
             family = "gaussian", argvals.functional = argvals.functional,
             silent = FALSE)
```

