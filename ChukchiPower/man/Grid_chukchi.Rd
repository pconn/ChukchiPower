\docType{data}
\name{Grid_chukchi}
\alias{Grid_chukchi}
\title{Chukchi sea gridded data}
\description{
  A list of 4 objects:
  1) Grid - a SpatialPolygonsDataFrame (see sp package) housing location of grid cells and associated covariates
  2) Adj - an adjacency matrix describing the neighbors of each grid cell.
  3) Chukchi_centroids - a SpatialPoints object giving the centroid of each grid cell
  4) K.data - a list housing information about knots used for smoothing.  This includes
     a) K - a (520 x 20) matrix holding distances from each grid cell to each knot (20         knots)
     b) Q.knot A structure matrix suitable for intrinsic conditionally autoregressive 
     (ICAR) modelling of knots (see e.g. Rue and Held 2005)
}
\author{
  Paul Conn \email{paul.conn@noaa.gov}
}
\keyword{data}
\keyword{datasets}

