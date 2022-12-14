# phageELF: Estimator of Lytic Function for DLS data

This package provides functions to analyze and plot Dynamic Light Scattering data of phage particles paired with titer data, train a model of titer loss and difference in AUC and perdict titer loss for new DLS data. It accompanies the paper "Counting Bacteriophages and Predicting Their Bioactivity Using Dynamic Light Scattering" by Dharmaraj et al. and the Shiny app [Phage-ELF](https://jp22.shinyapps.io/shinyapp/).

This package contains 7 functions:

- `run_all` is a wrapper function that allows you to replicate from start to end the analysis from the Phage-ELF shiny app. 
- `AUC_diff` calculates the difference in AUC between different treatments and their control from DLS data.
- `lm_train` trains a linear model with the difference in AUC as the indepedent variable and titer loss as the dependent variable.
- `lm_test` predicts titer loss using the trained model from `lm_train` and differences in AUC calculated from new DLS data.
- `plot_DLS` creates a plot of DLS data for one control and its associated treatments.
- `plot_model` creates a plot of the linear model with the training data.
- `plot_test_data` creates a plot of the linear model with both the training and test data. 

You can get more information about each function with the help function: `?run_all` . This requires installing the package first. 

## Installation

Run the following code in your R console in order to install this package:

``` r
install.packages("devtools")
devtools::install_github("jp22/phageELF")
```
If you have not used R before, you can install R [here](https://www.r-project.org/) and RStudio [here](https://www.rstudio.com/products/rstudio/). 
