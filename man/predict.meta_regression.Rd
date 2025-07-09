\name{predict.meta_regression}
\alias{predict.meta_regression}
\title{Predict from Florida Living Shoreline meta-analysis model}
\description{
Makes predictions from the Living Shoreline suitability meta-analysis model developed for the 2024 Florida Statewide Coastal Restoration Guide (SCRG). Requires the input of the meta-analysis model, updated with data from the sites for which suitability predictions will be made. Suitability for Living Shoreline restoration is represented as a probability for each location point, ranging 0-1 (0-100\%). If no data is available, the \code{\link[metafor]{predict.rma}} function can help retrieve the state-wide average suitability from the model average. See \code{vignette(package = "ecoinfoscrg")} for examples.
}
\usage{
\S3method{predict}{meta_regression}(data, meta_regression)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Cleaned and standardized output following \code{\link{wranglingCleaning}} and \code{\link{standardize}}. This should be a list of the following items: \code{state}, \code{predictors}, \code{numerical_vars}, \code{categorical_vars}, and \code{binary_vars}.}
  \item{meta_regression}{Meta-analytical regression model output from \code{\link{meta_regression}}.}
}
\details{
The predictions can be mapped to determine suitable locations for Living Shoreline restoration. See vignettes for examples.
}
\value{
Returns a vector of predicted values for the data provided.
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
\code{\link{meta_regression}}, \code{\link{wranglingCleaning}}, \code{\link{standardize}}, \code{\link[metafor]{predict.rma}}, and vignettes.
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
