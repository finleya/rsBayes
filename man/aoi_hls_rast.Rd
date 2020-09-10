\name{aoi_hls_rast}
\docType{data}

\alias{aoi_hls_rast}

\title{Enhanced Vegetation Index data}

\description{
  This \code{raster::RasterBrick} comprises 144 (21-by-21) 30 m Enhanced Vegetation
  Index (EVI) values computed from Harmonized Landsat
  Sentinel-2 (HLS) reflectance data (the L30 and S30 products) from January 1, 2014 to December 31, 2019, for tile ID 18TYN (available via
  \url{https://hls.gsfc.nasa.gov/data/}). This small area of interest is
  forested land on the Quabbin Reservoir Watershed System in central
  Massachusetts. Note, except for the format, this is the same data as \code{aoi_hls}.
}

\usage{
data(aoi_hls_rast)
}

\format{
 The bands in this 21-by-21 pixel \code{raster::RasterBrick} correspond to
 dates that are held in the \code{aoi_hls_rast_info} data frame.
}


\keyword{datasets}
