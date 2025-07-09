#' @noRd
#' @keywords internal
# Ordinal probit regression
ordinal <- function(formula,data)
{
  n <- nrow(data)

  # design matrix
  D <- stats::model.matrix(formula,data)

  # response variable # assumed to be 1:K
  Y <- stats::model.frame(formula,data,na.action = "na.pass")  # na.pass ensures that Y contains data
  Y <- stats::model.extract(Y,"response")
  # if(min(Y)!=1) { stop("Response variable needs to be an integer sequence 1:K") }
  if(min(Y)<1) { stop("Response variable needs to be an integer sequence 1:K") }

  # quantile scale
  K <- max(Y)
  P <- seq(0,1,length.out=K+1)
  Q <- stats::qnorm(P)
  Q[K+1] <- Inf # fix NaN

  nloglike <- function(par,zero=0)
  {
    # linear effect
    LE <- c(D %*% par)

    NLL <- sum(-log(stats::pnorm(Q[Y+1]-LE)-stats::pnorm(Q[Y]-LE)) - zero/n)
    return(NLL)
  }

  par <- numeric(ncol(D))
  names(par) <- colnames(D)
  OPT <- ctmm::optimizer(par,nloglike)
  par <- OPT$par
  loglike <- -OPT$value
  AIC <- 2*ncol(D) + 2*OPT$value

  # safe COV calculation
  DERIV <- ctmm:::genD(par,nloglike)
  COV <- ctmm:::cov.loglike(DERIV$hess,DERIV$grad)

  R <- list(est=par,COV=COV,loglike=loglike,AIC=AIC)
  return(R)
}
