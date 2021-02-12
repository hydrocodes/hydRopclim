\name{hydrochange}
\alias{hydrochange}
\title{hydrochange: Hydroclimatic change analysis at annual time step - Tm}
\usage{
hydrochange(data,lat)
}
\description{
Hydroclimatic change analysis at annual time step for a database that includes Mean Temperature. Includes: estimation of potential evapotranspiration by Oudin method, simulation of actual evapotranspiration by Budyko-Zhang model, quantifying impacts of climate and human activities on runoff change and quantifying watershed sensitivity and adaptation (Rau et al, 2018).
}
\arguments{
data: An annual dataframe object including Dates(Date), Precipitation (P in mm), Mean temperature (Tm in C) and Runoff (R in mm)

lat: A numerical value, watershed latitude
}
\examples{
hydrochange(data=database5, lat=-12)
}
