\name{meta_regression}
\alias{meta_regression}
\title{SCRG Living shoreline suitability meta-analysis model}
\description{
Meta-analytic regression model built as part of the 2024 Florida Statewide Coastal Restoration Guide (SCRG) project. The model is uses data from multiple studies based on the Living Shoreline Suitability Model framework, with the aim of predicting the suitability of sites for living shoreline restoration. This function can be used to update the state-wide model with additional local studies (see \code{vignette("Meta-Analysis_Model_Full_Workflow", package = "ecoinfoscrg")}) or to make suitability predictions for specific coastlines (see \code{vignette("Meta-Analysis_Predictions", package = "ecoinfoscrg")}).
}
\usage{
meta_regression(data, varCov)
}
\arguments{
  \item{data}{Cleaned and standardized output following \code{\link{wranglingCleaning}} and \code{\link{standardize}}. This should be a list of the following items: \code{state}, \code{predictors}, \code{numerical_vars}, \code{categorical_vars}, and \code{binary_vars}.}
  \item{varCov}{Output from the \code{\link{varCov}} function. This should be a list of the following items: \code{VAR}, \code{combined_betas}, \code{combined_betas_only}, \code{combined_se}, and \code{combined_se_only}.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{Returns a meta-analytic regression model built from the input models selected for each individual study. See \code{\link[metafor]{rma.mv}} for more details.}
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
\code{\link{varCov}}, \code{\link{wranglingCleaning}}, \code{\link{standardize}}, \code{\link[metafor]{rma.mv}}, and vignettes.}
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
