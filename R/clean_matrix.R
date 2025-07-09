#' @noRd
#' @keywords internal
# VAR-COV matrix cleaning
clean_matrix <- function(COV,precision=1/3)
{
 MIN <- .Machine$double.eps^precision
 MAX <- 1/MIN
 VAR <- diag(COV)
 cov_names <- dimnames(COV)

 # minimal bound
 VAR <- pmax(VAR,MIN)
 VAR <- pmin(VAR,MAX)
 SD <- sqrt(VAR)
 SD <- SD %o% SD

 MAX.COR <- Inf
 MAX.DVAR <- Inf
 while(TRUE)
 {
  COR <- COV / SD
  # remove NA and NaN
  NAS <- is.na(COR)
  COR[NAS] <- 0

  # clamp spurious correlations
  MAX.COR.OLD <- MAX.COR
  MAX.COR <- max(abs(COR))
  COR[] <- pmax(COR,-1)
  COR[] <- pmin(COR,+1)
  COV <- COR * SD

  # clamp variances
  EIGEN <- eigen(COV)
  EIGEN$values <- pmax(EIGEN$values,MIN)
  EIGEN$values <- pmin(EIGEN$values,MAX)
  COV <- EIGEN$vectors %*% diag(EIGEN$values) %*% t(EIGEN$vectors)
  MAX.DVAR.OLD <- MAX.DVAR
  MAX.DVAR <- max(abs((diag(COV)/VAR-1)))
  diag(COV) <- VAR

  # no longer improving
  if(MAX.COR.OLD-MAX.COR<=MIN && MAX.DVAR.OLD-MAX.DVAR<=MIN) { break }
 }

 dimnames(COV) <- cov_names

 return(COV)
}
