# hydRopclim: An R package for easy hydroclimatic calculations
<img src="https://github.com/hydrocodes/hydRopclim/blob/main/tutorial/hydropclim.jpeg" width="150">
## 1. What is hydRopclim?
`hydRopclim` is an R package for automatizing hydroclimatic calculations. `hydRopclim` seeks to helping non-expert R users in the hydroclimatological field. 
Some results and interpretations require a supervised opinion before drawing conclusions.

## 2. What is hydRopclim for?
Six main functions are implemented in `hydRopclim`. Their applications cover the topics of deterministic hydroclimatology and watershed hydrology
with a focus over the Peruvian Andes and Pacific slope and coast, however they can be applied to any similar geographical context under a critical hydrological criterion. 

The main functions are:
- `pgridcorr()`: Correction of monthly grid-point precipitation, e.g., TRMM (Tropical Rainfall Measuring Mission) (Condom, Rau & Espinoza, 2011) or others gridded products in function of an in-situ station .
- `tgridcorr()`: Correction of monthly grid-point mean temperature at target elevations, e.g., NCEP NCAR or others reanalysis products in function of an in-situ station (Rau, Condom & Lavado, 2013).
- `indexcorrl()`: Estimation of hydroclimatic seasonal indexes, and their running correlations, e.g., ENSO index versus precipitation (Bourrel, Rau, Labat et al, 2015). Also, it contains 4 complementary functions.
  - `seasavg()`, `seasavg2()` and `seassum()`: Calculation of seasonal average vector index for a season of n-months; seasonal average and sum matrix indexes for a season of n-months from m-hydroclimatic variables, respectively.
  - `zscorem()`: Transformation of m-hydroclimatic monthly variables into m-zscores indexes over 12-months each one.
- `hydrocluster()`: K-means clustering of a hydroclimatic time series, e.g., precipitation, streamflow; spatial visualization and evaluation by "silhouettes" (Rau, Bourrel, Labat et al, 2017).
- `hydrochange()`: Hydroclimatic change analysis at annual time step for a database that includes mean temperature. Estimation of potential evapotranspiration by Oudin method,
simulation of actual evapotranspiration by Budyko-Zhang model, quantifying impacts of climate and human activities on runoff change and quantifying watershed sensitivity and
adaptation (Rau, Bourrel, Labat et al, 2018).
  - `hydrochange2()`: Hydroclimatic change analysis at annual time step for a database that includes potential evapotranspiration. Simulation of actual evapotranspiration by
  Budyko-Zhang model, quantifying impacts of climate and human activities on runoff change and quantifying watershed sensitivity and adaptation (Rau, Bourrel, Labat et al, 2018).
- `rindex()`: Estimation of a monthly runoff index in an ungauged watershed through GR2M model and geomorphometric parameters (Rau, Bourrel, Labat et al, 2019).


## 3. How to install hydRopclim?
The `hydRopclim` package must be installed from Github hydrocodes repository, following the next 2 steps.

**Step 1**: In Rstudio, install `devtools` package from CRAN

**Step 2**: In Rstudio console or on your script, please write 

```r
devtools::install_github(c("hydrocodes/hydRopclim"))
```
or also :

```r
library(devtools)
install_github("hydrocodes/hydRopclim")
```
During the installation, please check in R console and skip other updates with an empty line or selecting option "None".

That’s all! Finally, do not forget call the package in your script and if is necessary install and call other packages required in some functions. 
Here a list of `hydRopclim` functions that work fine with the next packages:
- For `indexcorrl()`: `reshape2`, `ggplot2`, `wesanderson`, `cowplot`
- For `hydrocluster()`: `stats`, `cluster`, `sp`, `rgdal`
- For `rindex()`: `ggplot2`

Example: `rindex()` requires ggplot2 package for plotting runoff index superavit and deficit time series:
```r
library(hydRopclim)
library(ggplot2)
rindex(data,a,l,p)
```
Please, check tutorial folder for codelines examples and more details:
https://github.com/hydrocodes/hydRopclim/tree/main/tutorial

## 4. Credits
`hydRopclim` was developed by Pedro Rau and allows to reproduce some base figures shown in Rau (2017). For any issue or suggestion please write to: pedro.rau.ing@gmail.com

`hydRopclim` could be not possible without runnning the next softwares and packages: R (R Core Team, 2020), Rstudio (RStudio Team, 2020), stats (R Core Team, 2020), cluster (Maechler et al, 2019), sp (Bivand et al, 2013), rgdal (Bivand et al, 2015), ggplot2 (Wickham, 2016), reshape2 (Wickham, 2007), wesanderson (Ram et al, 2018), cowplot (Wilke, 2020).

## 5. Versions
v 1.0 - January 30, 2021

## 6. How to cite?

Rau, P. 2021. hydRopclim: An R package for easy hydroclimatic calculations. figshare. Software. https://doi.org/10.6084/m9.figshare.13670191, 
GitHub repository: https://github.com/hydrocodes/hydRopclim

## 7. References:

Bivand RS, Pebesma E, Gomez-Rubio V, 2015. Applied spatial data analysis with R, Second edition. Springer, NY. https://asdar-book.org/

Bivand RS, Keitt T, Rowlingson B, 2015. rgdal: Bindings for the Geospatial Data Abstraction Library. R package version 1.0–4. http://CRAN.R-project.org/package=rgdal

Bourrel L, Rau P, Dewitte B, Labat D, Lavado W, Coutaud A, Vera A, Alvarado A, Ordoñez J, 2015. Low-frequency modulation and trend of the relationship between precipitation and ENSO along the Northern to Center Peruvian Pacific coast. Hydrological Processes. 29(6):1252-1266. http://dx.doi.org/10.1002/hyp.10247

Condom T, Rau P, Espinoza JC, 2011. Correction of TRMM 3B43 monthly precipitation data over the mountainous areas of Peru during the period 1998-2007. Hydrological Processes. 25(12):1924-1933. http://dx.doi.org/10.1002/hyp.7949

Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2019). cluster: Cluster Analysis Basics and Extensions. R package version 2.1.0. https://cran.r-project.org/web/packages/cluster/index.html

R Core Team, 2020. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

Ram K, Wickham H , Richards C, Baggett A, 2018. A Wes Anderson Palette Generator. https://cran.r-project.org/web/packages/wesanderson/index.html

Rau P, Condom T, Lavado, W. 2013. Spatio-temporal analysis of monthly temperature in the mountainous regions of Peru. An approach for NCEP NCAR Reanalysis data correction. Proceedings of the 35th IAHR World Congress. 12:10602-10612. https://doi.org/10.13140/2.1.4591.9522

Rau P, Bourrel L, Labat D, Melo P, Dewitte B, Frappart F, Lavado W, Felipe O, 2017. Regionalization of rainfall over the Peruvian Pacific slope and coast. International Journal of Climatology 37(1):143-158. http://dx.doi.org/10.1002/joc.4693

Rau P. 2017. Precipitation, runoff and water balance regimes variability along the Peruvian Pacific slope and coast: ENSO influence and sensitivity to hydroclimatic change (PhD thesis). Université Toulouse III Paul Sabatier. France 267pp, tel-01627597. https://tel.archives-ouvertes.fr/tel-01627597

Rau P, Bourrel L, Labat D, Frappart F, Ruelland D, Lavado W, Dewitte B, Felipe O, 2018. Hydroclimatic change disparity of Peruvian Pacific drainage catchments. Theoretical and Applied Climatology. 134(1-2):139-153. http://dx.doi.org/10.1007/s00704-017-2263-x

Rau P, Bourrel L, Labat D, Ruelland D, Frappart F, Lavado W, Dewitte B, Felipe O, 2019. Assessing multi-decadal runoff (1970‒2010) using regional hydrological modelling under data and water scarcity conditions in Peruvian Pacific catchments. Hydrological Processes. 33(1):20-35. https://doi.org/10.1002/hyp.13318

RStudio Team (2020). RStudio: Integrated Development for R. RStudio, PBC, Boston, MA URL http://www.rstudio.com/

Wickham H (2007). “Reshaping Data with the reshape Package.” Journal of Statistical Software, 21(12), 1–20. http://www.jstatsoft.org/v21/i12/

Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org

Wilke, 2020. cowplot: Streamlined Plot Theme and Plot Annotations for 'ggplot2'. https://cran.r-project.org/web/packages/cowplot/index.html
