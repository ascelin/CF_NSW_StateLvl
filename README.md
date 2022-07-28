# CF_NSW_StateLvl

#Packages required for running the analysis
#Important: the package "mlr3spatiotempcv" has an issue with the coordinate system so need to install an older version of the package manually first

#Load required libraries
packages <- c("sf","terra","data.table","tidyverse",
              "mlr3","mlr3learners","mlr3viz","mlr3tuning",
              "iml","future","furrr","purrr","xgboost",
              "lattice","tictoc","scico","ggtext","mlr3spatiotempcv")

#Install the old version of mlr3spatitempcv manually# we need devtools to do that so first install that package
if(!require(devtools)) install.packages("devtools")
library(devtools)
if (!require("mlr3spatiotempcv")) install_version("mlr3spatiotempcv", 
                                                  version = "1.0.1", 
                                                  repos = "http://cran.us.r-project.org")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

#Load the packages
lapply(packages, require, character.only=TRUE)

