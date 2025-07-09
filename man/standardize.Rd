\name{standardize}
\alias{standardize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Standardize and impute missing data for Florida Living Shoreline meta-analysis model}
\description{
Standardizes predictor data after cleaning (see \code{\link{wranglingCleaning}}) to be centered on 0 with a standard deviation of 1, ensuring that values are comparable across sites and variables. Any categorical variables are converted into dummy variable format (categories as individual indicators). This function also includes the option to impute any missing data at the study-specific scale or statewide, using the median or mean values for each predictor. Kriging, a geostatistical technique for spatial interpolation can also be used to impute locally.
}
\usage{
standardize(data, site.method = c("krige", "medianImpute", "meanImpute"),
            state.method = c("meanImpute", "medianImpute"), duplicates = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Cleaned output from \code{\link{wranglingCleaning}}. This should be a list of the following items: \code{state}, \code{predictors}, \code{numerical_vars}, \code{categorical_vars}, and \code{binary_vars}.}
  \item{site.method}{Within-study imputation method: "krige" (recommended) for auto-Kriging using \code{\link[automap]{autoKrige}}, "medianImpute" to use the site-wide median using \code{\link[caret]{preProcess}}, or "meanImpute" to use the site-wide mean. \code{NULL} to skip this step.}
  \item{state.method}{State-wide imputation method for predictor information not available at the study site level. Options include only "medianImpute" and "meanImpute" (recommended). \code{NULL} to skip this step.}
  \item{duplicates}{\emph{Only for Kriging}. Logical, \code{TRUE} removes duplicate points from the input data, however, this can be time-consuming for large datasets.}
}
\details{
Standardization uses a stored mean and standard deviation value from the studies used to build the original meta-analysis model. This ensures that any additional data is appropriately handled for the model. Note that the data sets used to build the original model were imputed using site-wide Kriging and the state-wide mean. See \code{vignettes(package = "ecoinfoscrg")} for examples of use.
}
\value{
Returns a list of the same format as \code{\link{wranglingCleaning}}, updated with a standardized (and imputed) data frame of predictors.
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
\code{\link{wranglingCleaning}} and vignettes.
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
