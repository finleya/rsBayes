\name{aoi_hls}
\docType{data}

\alias{aoi_hls}

\title{Enhanced Vegetation Index data}

\description{
  The dataset comprises 144 (21-by-21) 30 m Enhanced Vegetation
  Index (EVI) values computed from Harmonized Landsat
  Sentinel-2 (HLS) reflectance data (the L30 and S30 products) from January 1, 2014 to December 31, 2019, for tile ID 18TYN (available via
  \url{https://hls.gsfc.nasa.gov/data/}). This small area of interest is
  forested land on the Quabbin Reservoir Watershed System in central
  Massachusetts. 
}

\usage{
data(aoi_hls)
}

\format{
  The data are held in a 273420-by-7 data frame with the following columns:
  \describe{
  \item{pixel}{is the pixel id}
  \item{x}{easting coordinate}
  \item{y}{northing coordinate}
  \item{sat}{L30 for Landsat 8 and S30 for Sentinel-2}
  \item{year}{of data collection}
  \item{doy}{day of the year of data collection}
  \item{evi}{computed enhanced vegetation index}
}

The proj4string for these data is  "+proj=utm +zone=18 +ellps=WGS84 +units=m +no_defs".
}


\keyword{datasets}
