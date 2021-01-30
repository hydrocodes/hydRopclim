\name{pgridcorr}
\alias{pgridcorr}
\title{pgridcorr: Correcting precipitation grid}
\usage{
pgridcorr(data)
}
\description{
Correction of monthly grid-point precipitation, e.g., TRMM (Tropical Rainfall Measuring Mission) or others gridded products in function of an in-situ station (Condom et al, 2011).
}
\arguments{
data: A monthly dataframe object including Dates (Date), in-situ data (Station) and grid data (Grid)
}
\examples{
pgridcorr(data=database1)
}
