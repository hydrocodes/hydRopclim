\name{spatial_grad}
\alias{spatial_grad}
\title{spatial_grad: Mean spatial climatology corrected by gradients for a basin}
\usage{
spatial_grad(DEM, temp_stations, prec_stations, ccas, grad_temp, grad_pr)
}
\description{
Mean spatial climatology (e.g. monthly precipitation, temperature and potential evapotranspiration) corrected by gradients for a basin (Rau et al, 2018). 
}
\arguments{
DEM: A digital elevation model in masl

temp_stations: A dataframe object: head station names, long in degrees, lat in degrees, elevation in masl and monthly temperature data with dates %b-%Y

prec_stations: A dataframe object: head station names, long in degrees, lat in degrees, elevation in masl and monthly precipitation data with dates %b-%Y

ccas: A shapefile object containing polygons representing basins 

grad_temp: A numerical value of temperature gradient (C/m) 

grad_pr: A numerical value of precipitation gradiente (mm/m) 
}
\examples{
spatial_grad(DEM, temp_stations, prec_stations, ccas, grad_temp=-6.5/1000, grad_pr=4*10^(-4))
}
