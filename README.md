![alt text](https://github.com/wheelemw/ToxicRDocs/blob/main/Toxic-R_Web_Graphic_V13.jpg)

# Welcome to the Official Git Repository of ToxicR!

## What is ToxicR?

- test workflow 

ToxicR is an R package that utilizes the core functionality of the US EPA's Benchmark Dose Software and the NTP's BMDexpress, which was developed in unison with researchers from the National Institute of Environmental Health. Unlike those packages, ToxicR is fully interactive and can be used to create custom analysis pipelines within the R programming environment. 

## Current Release

Version 1.1.0 is the most recent release. 


### Installation

Depending on your system, cut and paste the following code into your R terminal. 

## Recommended Method
### Compile Yourself
If you have the package devtools, you can download and install directly from GitHub!


`library(devtools)`

`devtools::install_github("NIEHS/ToxicR")`

***Note:*** For Windows, you will need the rtools executable available at: https://cran.r-project.org/bin/windows/Rtools/

***Note:*** If you have a MacOS, you will need to download the GNU Scientific Library. 
To do this, go to a command line and type

`brew install gsl`

This assumes you have HomeBrew installed. If you do not go to https://brew.sh, which will give you the instructions on how to install. 

***Note:*** For MacOS, you also need to install xcode.  To do this go to the command line and type

`xcode-select -install`

***Note:*** For Linux, you will also need the GNU Scientific Library.  The install depends on your flavor of Linux. 
For Ubuntu, type

`sudo apt-get install libgsl-dev`

## Alternative Methodology

**First, install the required packages**
 
`install.packages(c("Rcpp","RcppEigen","RcppGSL","ggplot2","shiny","coda","scales","tidyverse","forcats","ggridges","doBy","multcomp","dplyr","rmarkdown", "actuar","ggpubr", "testthat","gridExtra","VIM","knitr", "modules", "plotly" ))`

### Windows R 4.3.0

`download.file("https://github.com/NIEHS/ToxicR/releases/download/v1.10.0/ToxicR_23.4.1.1.0.R4.3.zip", 
              "ToxicR_23.4.1.1.0.zip")`
`install.packages("ToxicR_23.4.1.1.0.zip", repos = NULL, type = "win.binary")`

### Windows R 4.2.3

`download.file("https://github.com/NIEHS/ToxicR/releases/download/v1.10.0/ToxicR_23.4.1.1.0R4.2.3.zip", 
              "ToxicR_23.4.1.1.0.zip")`
`install.packages("ToxicR_23.4.1.1.0.zip", repos = NULL, type = "win.binary")`

## Citing ToxicR

All dose-response methodologies used in ToxicR were developed in the following manuscripts: 

Wheeler, M.W., Blessinger, T., Shao, K., Allen, B.C., Olszyk, L., Davis, J.A. and Gift, J.S. (2020), Quantitative Risk Assessment: Developing a Bayesian Approach to Dichotomous Dose–Response Uncertainty. Risk Analysis, 40: 1706-1722. https://doi.org/10.1111/risa.13537

Wheeler, M. W., Cortiñas Abrahantes, J., Aerts, M., Gift, J. S., & Allen Davis, J. (2022). Continuous model averaging for benchmark dose analysis: Averaging over distributional forms. Environmetrics, e2728. https://doi.org/10.1002/env.2728

Interuniversity Institute for Biostatistics and statistical Bioinformatics, 2022. EFSA Platform for Bayesian Benchmark Dose Analysis. EFSA Supporting Publications, 19(12), p.7740E.

Wheeler, M.W., Lim, S., House, J.S., Shockley, K.R., Bailer, A.J., Fostel, J., Yang, L., Talley, D., Raghuraman, A., Gift, J.S. and Davis, J.A., 2023. ToxicR: A computational platform in R for computational toxicology and dose–response analyses. Computational Toxicology, 25, p.100259.

<img src="https://github.com/wheelemw/ToxicRDocs/blob/main/NIEHS.png" width="100" height="100"> <img src="https://github.com/wheelemw/ToxicRDocs/blob/main/NTP.gif" width="200" height="100">
  
