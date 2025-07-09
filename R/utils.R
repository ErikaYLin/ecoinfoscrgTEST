# library(sf)
# library(stringdist)
#' @noRd
#' @keywords internal
scrg__drawlinestr <- function(coords, crs=4326) {
  return(sf::st_sfc(sf::st_linestring(matrix(coords, ncol=2)), crs=crs))
}

#' @noRd
#' @keywords internal
scrg__get_features_from_sf <- function(sf) {
  # This function takes a geospatial object and returns the features as a data
  # frame.
  return(as.data.frame(sf::st_drop_geometry(sf)))
}

#' @noRd
#' @keywords internal
scrg__get_line_from_linestrings <- function(
      lstrings, test_line, element_number, crs="EPSG:4326"
    )
  {
  total_length <- length(test_line[[1]][[1]])
  return(
    paste(
      "LINESTRING(",
      sf::st_geometry(lstrings)[[1]][[1]][element_number],
      " ",
      sf::st_geometry(test_line)[[1]][[1]][element_number+(total_length/2)],
      ")") |> sf::st_as_sfc(crs = crs)
  )
}

#' @noRd
#' @keywords internal
scrg__get_geometry_in_bbox <- function(points, sf_object) {
  # This function takes a series of points and an sf object, and returns a
  # subset of the sf object that falls within the bounding box formed by the
  # series of points.
  p1 <- points[1:2]
  p2 <- points[3:4]
  bbox_matrix <- matrix(
    c(p1[1], p1[2], p2[1], p1[2], p2[1], p1[2], p2[1], p2[2], p2[1], p2[2],
      p1[1], p2[2], p1[1], p2[2], p1[1], p1[2]),
    ncol=2,
    byrow=TRUE
  )
  bbox <- sf::st_polygon(list(bbox_matrix))
  return(subset(sf_object, sf::st_within(sf_object, bbox, sparse = FALSE)))
}

#' @noRd
#' @keywords internal
scrg__get_geometry_in_bbox2 <- function(point_list, sf_object, byrow=TRUE) {
  # This function takes a series of points and an sf object, and returns a
  # subset of the sf object that falls within the bounding box formed by the
  # series of points.
  bbox_matrix <- matrix(point_list, ncol=2, byrow=byrow)
  bbox <- sf::st_polygon(list(bbox_matrix))
  return(subset(sf_object, sf::st_within(sf_object, bbox, sparse = FALSE)))
}

#' @noRd
#' @keywords internal
scrg__get_geometry_in_polygon <- function(polygon, sf_object) {
  # Returns an sf object that is a subset of all entries with a geometry that
  # is within the provided polygon.
  return(subset(sf_object, sf::st_within(sf_object, polygon, sparse = FALSE)))
}

#' @noRd
#' @keywords internal
scrg__make_string_distance_diff <- function(x, y, dropMatches=NULL) {
  # Takes two character vectors and compares the string distances between each
  # element of both vectors. Function returns a matrix of the results.
  x <- sort(x)
  y <- sort(y)
  str_dists <- list()
  to_drop <- list()
  for(a in x){
    for(b in y){
      dist <- stringdist::stringdist(a,b)
      str_dists <- append(str_dists, dist)
      if((dist == 0) && (dropMatches)){
        to_drop <- append(to_drop, a)
      }
    }
  }
  dists <- matrix(
    str_dists,
    nrow=length(x),
    ncol=length(y),
    dimnames=list(y,x)
  )
  if(dropMatches){
    return(dists[,x[!(x %in% to_drop)]])
  } else {
    return(dists)
  }
}

#' @noRd
#' @keywords internal
scrg__multiline_length <- function(multiline_str) {
  # Returns the number of LINESTRING in a MULTILINESTRING
  return(multiline_str[[1]][[1]] / 2)
}

#' @noRd
#' @keywords internal
scrg__venn_set <- function(x, y) {
  # Making a venn diagram of two vectors to determine what is not in the other
  # vector
  # From: https://stackoverflow.com/a/17599164
  both <- union(x,y)
  inX <- both %in% x
  inY <- both %in% y
  return(table(inX,inY))
}
