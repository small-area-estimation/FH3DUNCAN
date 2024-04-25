# Model-based estimation of small area dissimilarity indexes: An application to sex occupational segregation in Spain

### Bugallo M., Esteban M.D., Morales D., Pagliarella M.C.

### Universidad Miguel Hernández de Elche, Spain

### Università degli Studi di Cassino e del Lazio Meridionale, Italy

## Objective

This repository contains the R scripts of the *Application to real data* from the paper entitled "Model-based estimation of small area dissimilarity indexes: An application to sex occupational segregation in Spain".

## Installation and requirements

The R code available in this repository was tested on a 3.3 GHz Intel Xeon W Mac computer with 32 GB RAM. For data manipulation, task processing and graphical results, the following R packages are required:

-   library(car)
-   library(maptools)
-   library(RColorBrewer)

In order to facilitate the manipulation of the code, we recommend the GUI, "Graphical User Interface", [R Studio](https://posit.co/downloads/). Once all the libraries are installed, simply download the project sources and run the scripts.

## Guide of use

To obtain a complete set of results of the *Application to real data* section of the paper, you need to follow the steps:

1.  *Obtain direct estimators*. Run the R script `FH3_EPA.R` which returns the output result files `YHajek_covariates.csv` and `var_covariates.csv`. The datasets required to run the script are `EPA_2019T4.csv`, `EPA_2020T1.csv`, `EPA_2020T2.csv`, ...,`EPA_2021T4.csv`, which must be allocated in the same folder as the script.

2.  *Fit a three-fold Fay-Herriot model*. Run the R script `FH3_model.R`. The data sets needed to run this script are `YHajek_covariates.csv` and `var_covariates.csv` (results of the previous step), which must be in the same folder as the script. It also needs the R scripts `REML.3foldFH.indep.R` and `BETA.U.3foldFH.indep.R` and `stderr.3foldFH.indep.R`.

3.  *Obtain Bootstrap MSE estimators*. Run the R script `FH3_BootRMSE.R`. The data sets required to run this script are `YHajek_covariates.csv` and `var_covariates.csv`, which must be in the same folder as the script. Please note that the bootstrap procedure uses a lot of computing time. Our script has been created to perform 500 bootstrap iterations. We recommend decreasing this number if you want to test the code.

4.  *Map the results*. Using the maptools package, the `Maps_DSI.R` script provides several map plots of the above results. Please note that you need to have the files corresponding to the location you want to plot. In our case we have these files in the \`spainmap' folder.

## Warning

The maptools package is now archived. If you experience problems with this package, maybe you can fix it with the following instance:

``` r
install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")
```
