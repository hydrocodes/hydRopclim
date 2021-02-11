\name{hydrochange2}
\alias{hydrochange2}
\title{hydrochange2: Hydroclimatic change analysis at annual time step - PET}
\usage{
hydrochange2(data,lat)
}
\description{
Hydroclimatic change analysis at annual time step for a database that includes Potential Evapotranspiration. Includes: Simulation of actual evapotranspiration by Budyko-Zhang model, quantifying impacts of climate and human activities on runoff change and quantifying watershed sensitivity and adaptation (Rau et al, 2018).
}
\arguments{
data: An annual dataframe object including Dates(Date), Precipitation (P in mm), Potential Evapotranspiration (PET in mm) and Runoff (R in mm)

lat: A numerical value, watershed latitude
}
\examples{
hydrochange2(data=database6, lat=-12)
}
