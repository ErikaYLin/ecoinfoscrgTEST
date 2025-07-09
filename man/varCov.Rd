\name{varCov}
\alias{varCov}
\title{Variance-Covariance matrix preparation for meta-analysis model}
\description{
Retrieves and cleans the variance-covariance matrix for the final selected models from each study site. The adjusted variances extracted from the covariance matrices and the adjusted covariance matrices are used to build the meta-analytic regression model (see \code{\link{meta_regression}}). This function includes the \code{\link{getBetas}} function internally, thus it is not necessary to use run \code{\link{getBetas}} prior to using this function.
}
\usage{
varCov(data, predictors, mods, Betas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Cleaned and standardized output following \code{\link{wranglingCleaning}} and \code{\link{standardize}}. This should be a list of the following items: \code{state}, \code{predictors}, \code{numerical_vars}, \code{categorical_vars}, and \code{binary_vars}.}
  \item{predictors}{Vector of numeric and factored predictor variables to be used in the meta-analysis model. See \code{vignette(package = ecoinfoscrg)} for examples.}
  \item{mods}{List of the same length as the number of studies, where each item is the selected model for the study site. Model selection is to be conducted using the "build-up pair-down" method, via \code{BUPD} (parallelized) or \code{BUPD_nonparallel} (for Windows operating systems).}
  \item{Betas}{List of the same length as the number of studies, where each item is the model summary of beta and standard error estimates for selected predictors. These results should be saved upon running model selection for each study. See \code{BUPD} or \code{BUPD_nonparallel} and associated vignettes at \code{vignette(package = "ecoinfoscrg")}.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Returns a list of:

  \item{VAR}{Data frame of variance values extracted from the adjusted variance-covariance matrix of each study.}
  \item{adjusted_cov_matrix}{List of cleaned variance-covariance matrices from the selected models of studies to be included in the meta-analysis model.}
  \item{combined_betas}{Data frame of the combined and standardized beta estimates for all predictors across studies, including a column that indicates the study.}
  \item{combined_betas_only}{Data frame of the combined and standardized beta estimates for all predictors across studies.}
  \item{combined_se}{Data frame of the combined and standardized standard error estimates for all predictors across studies, including a column that indicates the study.}
  \item{combined_se_only}{Data frame of the combined and standardized standdard error estimates for all predictors across studies.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{getBetas}}, \code{\link{meta_regression}}, and vignettes.
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.

%%## The function is currently defined as

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
