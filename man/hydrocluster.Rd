\name{hydrocluster}
\alias{hydrocluster}
\title{hydrocluster: K-means clustering of stations hydroclimatic time series}
\usage{
hydrocluster(file,shp,clusters)
}
\description{
K-means clustering of a hydroclimatic time series, e.g., precipitation, streamflow, and evaluation by silhouette (Rau et al, 2017). Spatial visualization and control quality by silhouettes.
}
\arguments{
file: A dataframe object including head station names, latitude, longitude and monthly or annual data

shp: A spatial object, a polygon shapefile of study region

clusters: A numeric value, number of clusters
}
\examples{
hydrocluster(file=database4,shp=region,clusters=4)
}
