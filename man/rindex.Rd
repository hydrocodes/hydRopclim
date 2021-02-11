\name{rindex}
\alias{rindex}
\title{rindex: Estimation of runoff index in an ungauged basin}
\usage{
rindex(data,a,l,p)
}
\description{
Estimation of a monthly runoff index in an ungauged watershed through GR2M model and geomorphometric parameters (Rau et al, 2019).
}
\arguments{
data: A monthly dataframe object including Dates (Date), Precipitation (P) and  Potential Evapotranspiration (PET) in mm

a: A numeric value: Watershed area in km2

l: A numeric value: Watershed main channel lenght in km

p: A numeric value: Watershed perimeter in km2
}
\examples{
rindex(data=database7, a=2352, l=88.3, p=261)
}
