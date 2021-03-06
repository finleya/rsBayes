\name{aoi_hls_rast_info}
\docType{data}

\alias{aoi_hls_rast_info}

\title{Enhanced Vegetation Index data}

\description{
  The rows in this data frame hold the satellite, year, and day of the
  year (doy) for the bands in the \code{raster::RasterBrick}
  \code{aoi_hls_rast}.
}

\usage{
data(aoi_hls_rast_info)
}

\format{
  The data are held in a 651-by-3 data frame with the following columns:
  \describe{
  \item{sat}{L30 for Landsat 8 and S30 for Sentinel-2}
  \item{year}{of data collection}
  \item{doy}{day of the year of data collection}
}

Rows in \code{aoi_hls_rast_info} correspond to bands in
\code{aoi_hls_rast}.
}


\keyword{datasets}
