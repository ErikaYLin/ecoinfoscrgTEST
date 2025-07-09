\name{utils}
\alias{scrg__drawlinestr}
\alias{scrg__get_features_from_sf}
\alias{scrg__get_line_from_linestrings}
\alias{scrg__get_geometry_in_bbox}
\alias{scrg__get_geometry_in_bbox2}
\alias{scrg__get_geometry_in_polygon}
\alias{scrg__make_string_distance_diff}
\alias{scrg__multiline_length}
\alias{scrg__venn_set}
\title{Various functions for handling geospatial data}
\description{
\emph{Internal functions} for handling geospatial data using \code{\link{sf}} for the 2024 SCRG project.
}
\usage{scrg__drawlinestr(coords, crs=4326)
scrg__get_features_from_sf(sf)
scrg__get_line_from_linestrings(lstrings, test_line,
                                element_number, crs="EPSG:4326")
scrg__get_geometry_in_bbox(points, sf_object)
scrg__get_geometry_in_bbox2(point_list, sf_object, byrow=TRUE)
scrg__get_geometry_in_polygon(polygon, sf_object)
scrg__make_string_distance_diff(x, y, dropMatches=NULL)
scrg__multiline_length(multiline_str)
scrg__venn_set(x, y)}
%%\arguments{
%%  \item{x}{
%%     ~~Describe \code{x} here~~
%%}
%%}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
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
