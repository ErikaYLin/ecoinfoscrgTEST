\name{BUPD}
\alias{BUPD}
\title{Model selection for Living Shoreline Suitability Model studies}
\description{
Conducts AIC-based model selection using the "build-up pair-down" method, in which the model is evaluated by adding predictors to the best model in a forward stepwise fashion (build-up phase), then removed in a backward stepwise fashion (pair-down phase). Ordinal probit regression is used for this data to model an ordered, categorical response for shoreline suitability (categories 1, 2, 3, with 3 being the best and 1 being the worst) as a function of the combined predictors from the Living Shoreline suitability studies. The probit link function allows for the response to be transformed from integers of 1, 2, and 3, to a probability ranging 0-1. This function has an option for parallel computing (for faster run time), however, we recommend using the non-parallel option on non-Unix-based operating systems (i.e. Windows systems).
}
\usage{
BUPD(data, predictors, parallel = FALSE)
}
\arguments{
  \item{data}{Cleaned and standardized output following \code{\link{wranglingCleaning}} and \code{\link{standardize}}. This should be a list of the following items: \code{state}, \code{predictors}, \code{numerical_vars}, \code{categorical_vars}, and \code{binary_vars}.}
  \item{predictors}{Data frame of numeric and factored predictor variables with values combined across studies. This should follow the same format as the predictors object output from \code{\link{standardize}}. See \code{vignette("Meta-Analysis_Model_Full_Workflow", package = ecoinfoscrg)} and other vignettes for examples.}
  \item{parallel}{Logical, if \code{TRUE}, the process is parallelized, which allows faster computation on systems where enabled.}
}
\details{
This function is to be used to add studies to the meta-analysis model and should not be considered a standalone function.
}
\value{
Returns a list of model selection results for each study, which include:
  \item{final_model}{Final selected model object for the study.}
  \item{final_form}{Character string, formula of the selected model.}
  \item{odds_ratios}{Odds ratios of the selected model, calculated using the beta coefficents}
  \item{coefs}{Model summary in table format. Includes the beta coefficient estimates and their standard error values.}
  \item{run_time}{Total computational time for the BUPD process of each study.}
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
\code{\link{standardize}}, \code{\link{getBetas}}, \code{\link{varCov}}, \code{\link{meta_regression}}, and vignettes.
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.

%%## The function is currently defined as
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
