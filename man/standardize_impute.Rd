\name{standardize_impute}
\alias{krigePredictors}
\alias{medianImpute}
\alias{meanImpute}
\title{Imputation techniques for missing predictor data in SCRG model
}
\description{
Options for imputing missing predictor data to be used in the 2024 SCRG Living Shoreline suitability meta-analysis model. Imputation can be performed within studies or at the state-wide level. \emph{Internal functions.}
}
\usage{krigePredictors(site, var, formula, duplicates = TRUE)
medianImpute(site, var, method = "medianImpute", k = 5)
meanImpute(site, var)}
\arguments{
  \item{site}{Study site.}
  \item{var}{Variable to be imputed}
  \item{formula}{Formula for Kriging model}
  \item{duplicates}{Logical for whether or not to remove duplicate entries prior to Kriging. Can be time-consuming for large data sets.}
  \item{method}{Always set to "medianImpute" to ensure that \code{\link[caret]{preProcess}} uses the correct imputation method.}
  \item{k}{Technically not needed, but would be passed on to \code{\link[caret]{preProcess}} for K-nearest neighbors imputation.}
}

%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
 \keyword{ internal }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
