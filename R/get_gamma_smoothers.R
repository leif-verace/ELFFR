#' Helper function to obtain sandwich smoother matrices
#' for analytic inference
get_gamma_smoothers <-
  function(data,
           knots,
           p,
           m,
           subj = NULL,
           covariates = NULL,
           knots.option = "equally-spaced",
           periodicity = c(FALSE, FALSE),
           lambda = NULL,
           selection = "GCV",
           search.grid = T,
           search.length = 100,
           method = "L-BFGS-B",
           lower = -20,
           upper = 20,
           control = NULL) {
    ## data dimension
    data_dim = dim(data)
    n1 = data_dim[1]
    n2 = data_dim[2]

    ## subject ID
    if(is.null(subj)) subj = 1:n2
    subj_unique = unique(subj)
    I = length(subj_unique)
    ## covariates for the two axis
    if(!is.list(covariates)) {

      x=(1:n1)/n1-1/2/n1; ## if NULL, assume equally distributed
      z = (1:n2)/n2-1/2/n2
    }
    if(is.list(covariates)){

      x = covariates[[1]]
      z = covariates[[2]]
    }

    ## B-spline basis setting
    p1 = rep(p,2)[1]
    p2 = rep(p,2)[2]
    m1 = rep(m,2)[1]
    m2 = rep(m,2)[2]


    select_knots <- function(t,knots=10,p=3,option="equally-spaced"){

      qs <- seq(0,1,length=knots+1)

      if(option=="equally-spaced"){
        knots <- (max(t)-min(t))*qs + min(t)
      }
      if(option=="quantile"){
        knots <- as.vector(quantile(t,qs))
      }

      K <- length(knots)
      knots_left <- 2*knots[1]-knots[p:1+1]
      knots_right <- 2*knots[K] - knots[K-(1:p)]

      return(c(knots_left,knots,knots_right))
    }
    ## knots
    if(!is.list(knots)){
      K1 = rep(knots,2)[1]
      xknots = select_knots(x,knots=K1,option=knots.option)
      K2 = rep(knots,2)[2]
      zknots = select_knots(z,knots=K2,option=knots.option)
    }

    if(is.list(knots)){

      xknots = knots[[1]]
      K1 = length(xknots)-1
      knots_left <- 2*xknots[1]-xknots[p1:1+1]
      knots_right <- 2*xknots[K1] - xknots[K1-(1:p1)]
      xknots <- c(knots_left,xknots,knots_right)

      zknots= knots[[2]]
      K2 = length(zknots)-1
      knots_left <- 2*zknots[1]- zknots[p2:1+1]
      knots_right <- 2*zknots[K2] - zknots[K2-(1:p2)]
      zknots <- c(knots_left,zknots,knots_right)
    }
    ##
    pspline.setting <- function(x,knots=select_knots(x,35),p=3,m=2,periodicity=FALSE,weight=NULL){

      # x: the marginal data points
      # knots: the list of interior knots or the numbers of interior knots
      # p: degrees for B-splines, with defaults values 3
      # m: orders of difference penalty, with default values 2
      #require(splines)
      #require(Matrix)

      ### design matrix
      K = length(knots)-2*p-1
      B = splines::spline.des(knots=knots, x=x, ord = p+1,outer.ok = TRUE)$design
      if(periodicity){
        Bint = B[,-c(1:p,K+1:p)]
        Bleft = B[,1:p]
        Bright = B[,K+1:p]
        B = cbind(Bint,Bleft+Bright)
      }


      difference.penalty <-function(m,p,K,periodicity=periodicity){

        # parameter  m: difference order
        # parameter  p: degree of B-splines
        # parameter  K: number of interior knots
        c = rep(0,m+1)

        for(i in 0:m)
          c[i+1] = (-1)^(i+1)*factorial(m)/(factorial(i)*factorial(m-i))

        if(!periodicity){

          M = matrix(0,nrow=K+p-m,ncol=K+p)
          for(i in 1:(K+p-m)) M[i,i:(i+m)] = c
        }
        if(periodicity){

          M = matrix(0,nrow=K,ncol=K)
          for(i in 1:(K-m)) M[i,i:(i+m)] = c
          for(i in (K-m+1):K) M[i,c(i:K,1:(m-K+i))] = c
        }

        return(M)
      }


      P = difference.penalty(m,p,K,periodicity)
      P1 = Matrix::Matrix(P)
      P2 = Matrix::Matrix(t(P))
      P = P2%*%P1

      List = list(
        "B" = B,
        "P" = as.matrix(P))

      return(List)
    }
    List1 <- pspline.setting(x,xknots,p1,m1,FALSE)
    B1 <- List1$B
    P1 <- List1$P

    List2 <- pspline.setting(z,zknots,p2,m2,FALSE)
    B2 <- List2$B
    P2 <- List2$P

    return(list(B1 = B1, P1 = P1, B2 = B2, P2 = P2))
  }
