# phageELF: Estimator of Lytic Function for DLS data

This package provides functions to analyze and plot Dynamic Light Scattering data of phage particles paired with titer data, train a model of titer loss and AUC differences, and predict titer loss for new DLS data. It accompanies the paper "Monitoring and Predicting Phage Bioactivity Using Dynamic Light Scattering" by Dharmaraj et al. and the Shiny app [Phage-ELF](https://jp22.shinyapps.io/shinyapp/).

This package contains 7 functions:

- `run_all` is a wrapper function that allows you to replicate from start to end the analysis from the Phage-ELF shiny app. 
- `AUC_diff` calculates the difference in AUC between different treatments and their control from DLS data.
- `lm_train` trains a linear model with the difference in AUC as the indepedent variable and titer loss as the dependent variable.
- `lm_test` predicts titer loss using the trained model from `lm_train` and differences in AUC calculated from new DLS data.
- `plot_DLS` creates a plot of DLS data for one control and its associated treatments.
- `plot_model` creates a plot of the linear model with the training data.
- `plot_test_data` creates a plot of the linear model with both the training and test data. 

You can get more information about each function with the help function: `?run_all` for example. This requires installing the package first. 

## System Requirements

The `phageELF` package requires only a standard computer and should be compatible with Windows, Mac and Linux operating systems.

Users should have `R` version 4.2.2 installed. The package functions with all packages in their latest version as they appear on `CRAN` on March 2023: 

```
tibble_3.0.3
tidyr_1.1.0
ggplot2_3.3.2
```

## Installation

If you have not used `R` before, you can install R [here](https://www.r-project.org/) and RStudio [here](https://www.rstudio.com/products/rstudio/). 

Run the following code in your `R` console in order to install and load this package:

``` r
install.packages("devtools")
devtools::install_github("jpourtois/phageELF", build_vignettes = TRUE)
library("phageELF")
```
Package installation should only take a few seconds. 

## Demo

A demo of all the functions in the `phageELF` package is available in the demo vignette:

```
vignette('phageELF-Demo','phageELF')
```

The user must install the package with `build_vignettes = TRUE` to access the vignette. 
