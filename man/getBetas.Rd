\name{getBetas}
\alias{getBetas}
\title{Organize and combine beta estimates across studies}
\description{
Retrieves beta coefficient estimates from the selected model of each study (after \code{BUPD} or \code{BUPD_nonparallel}) and produces a table of combined beta estimates for predictors across studies. Coefficient estimates are standardized for comparability between studies. Standard error estimates are combined in the same fashion.
}
\usage{
getBetas(data, predictors, Betas)
}
\arguments{
  \item{data}{Cleaned and standardized output following \code{\link{wranglingCleaning}} and \code{\link{standardize}}. This should be a list of the following items: \code{state}, \code{predictors}, \code{numerical_vars}, \code{categorical_vars}, and \code{binary_vars}.}
  \item{predictors}{Vector of numeric and factored predictor variables to be used in the meta-analysis model. See \code{vignette(package = ecoinfoscrg)} for examples.}
  \item{Betas}{List of the same length as the number of studies, where each item is the model summary of beta and standard error estimates for selected predictors. These results should be saved upon running model selection for each study. See \code{BUPD} or \code{BUPD_nonparallel} and associated vignettes at \code{vignette(package = "ecoinfoscrg")}.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Returns a list of:
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
\code{\link{standardize}}, model selection scripts (\code{BUPD} or \code{BUPD_nonparallel}, and vignettes.
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
