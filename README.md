 

This repo is a fork of Perez et al 2020, with the code encapsulated into a package.

The package can be installed with

``` r
# install pak if necessary
install.packages("pak")
pak::pak("richardjtelford/thermalTolerance") 

## OR if this fails, try
install.packages("remotes")
remotes::install_github("richardjtelford/thermalTolerance") 
```

Open the package with 

``` r
library(thermalTolerance)
```

And run it with (but with many more bootstrap samples.)
``` r
# htdata is demo data
htol <- psiiht(temperature = htdata$temperature, fvfm = htdata$fvfm, id = htdata$id, 
          control_temp = 23, warming = TRUE, boots = 5)
summary(htol)
autoplot(htol)
```

The original readme follows:

This file contains from sample data from a paper (that can be found here: https://onlinelibrary.wiley.com/doi/10.1111/pce.13990), and a small piece of R code that can used to estimate the PSII heat tolerance of PSII measured with maximum quantum yield of PSII (i.e. Fv/Fm). The function is pretty straightforward - it will provide estimates of 3 different PSII heat tolerances; Tcrit, T50 & T95. For more information on these heat tolerances, you can read this paper: https://onlinelibrary.wiley.com/doi/abs/10.1111/pce.13990. The user defines the temperature treatments used, the response variable (Fv/Fm), if they want to plot their results, the number of bootstrap iteration they would like, and a unique identifier (like a species) that indicate how the data should be spilt before curves are estimated. 
 
Right now, Fv/Fm is the response variable, but other response variable like membrane leakage or cell death could be used. However,  the parameters in the nls function would have to be changed (e.g. the 0.8 value indicating the y-intercept would change).\
\
I figured different users will probably have different preference for plotting, so the plots that are made with this function are pretty basic. For example, multi-paneled plots are probably desirable for large datasets\
\
The only required package for this function is the 'car' package, which is used to in the 'logit' function, which is used to estimate starting parameters for the nls function. The following warning message is common 'In logit(FvFm) : proportions remapped to (0.025, 0.975)'

If you use this code or this data, please consider citing the papers above. Thanks.

https://onlinelibrary.wiley.com/doi/10.1111/pce.13990