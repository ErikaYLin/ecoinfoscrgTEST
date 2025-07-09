\name{wranglingCleaning}
\alias{wranglingCleaning}
\title{Clean shapefile (.shp) data for Florida Living Shoreline meta-analysis model}
\description{
Automates the process of cleaning and combining data sets to be used in the Florida Living Shoreline suitability model. The output is a list that includes the combined response variables and study names, the corresponding table of predictor variables, and vectors of the names of all numeric variables, categorical variables, and binary (indicator) variables. The cleaned data is to be subsequently processed using the \code{standardize} function. See \code{\link{standardize}} for more details and \code{vignette(package = "ecoinfoscrg")} for examples of usage.
}
\usage{
wranglingCleaning(data, response)
}
\arguments{
  \item{data}{A named list of shape (.shp) files, 1 file per study to be cleaned and merged.}
  \item{response}{A vector of shoreline suitability response variable names, 1 per study (vector length should be equal to data list length. These should correspond to the column name in each shape file that describes the results of a living shoreline suitability evaluation. If the data is meant for predicting suitability at a certain site, stating \code{response = NULL} will bypass the response column cleaning step.}
}
\details{
Note that this is not intended to be a stand-alone function; \code{wranglingCleaning} should be followed up with data standardization using the \code{standardize} function, which also includes the option to impute any missing data where possible (recommended). Please check the cleaned data for any format discrepancies prior to standardization, as some additional manual cleaning may be needed for specific cases. See \code{vignette("Meta-Analysis_Model_Full_Workflow", package = "ecoinfoscrg")} for an example of how to deal with minor discrepencies.
}
\value{
Returns a list with:
%%  If it is a LIST, use
\item{state}{A two-column data frame of all response values (all location points combined) and the corresponding study names. Reponses are categorical values 1, 2, and 3, with 1 indicating poor suitability and 3 indicating good suitability for living shoreline restoration. Note that existing living shorelines may be categorized as unsuitable for further restoration. If \code{response = NULL}, a single-column data frame of study names only will be returned instead.}
\item{predictors}{A data frame of combined predictors, merging all studies into a single predictor data set.}
\item{numerical_vars}{A vector of all numerical predictor variable names.}
\item{categorical_vars}{A vector of all categorical (discrete) predictor variable names.}
\item{binary_vars}{A vector of all binary predictor variable names that indicate the presence or absence of feature.}
%% ...
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
\code{\link{ecoinfoscrg}}, \code{\link{standardize}}, and vignettes.
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
