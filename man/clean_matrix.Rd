\name{clean_matrix}
\alias{clean_matrix}
\title{Covariance matrix cleaning in varCov for SCRG model}
\description{
Cleans variance-covariance matrix for varCov function for the 2024 SCRG Living Shoreline suitability meta-analysis model.
}
\usage{
clean_matrix(COV, precision = 1/3)
}
\arguments{
  \item{COV}{Covariance matrix of selected model.}
  \item{precision}{Power to determine machine precision value.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Returns a cleaned covariance matrix.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
 \keyword{ internal }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
