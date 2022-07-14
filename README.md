# TV MSGARCH
<!-- badges: start -->
[![R-CMD-check](https://github.com/mkaywins/MSGARCH/workflows/R-CMD-check/badge.svg)](https://github.com/mkaywins/MSGARCH/actions)
<!-- badges: end -->


Please note: This is a fork of [MSGARCH package](https://github.com/keblu/MSGARCH). 
More about `MSGARCH` available at [http://keblu.github.io/MSGARCH/](http://keblu.github.io/MSGARCH/).


This fork includes time-varying transition probabilities for the methods `FitML()` and 
`CreateSpec()`. This fork was created for the purpose of my master thesis. I do not take any responsibility for the functionality of the package and advise anyone to use this fork at your own risk.


### Installation & Setup

In order to install this fork of the [MSGARCH package](https://github.com/keblu/MSGARCH) one needs to install the `devtools` package first. 
``` r 
install.packages("devtools")
require("devtools")
```
Make sure that `MSGARCH` is uninstalled on your machine. Then download the fork of the package from this repository.

``` r
devtools::install_github("mkaywins/MSGARCH", subdir="Package")

```
Load the package.
``` r
library(MSGARCH)

```

### How to use?
Before we fit a (TV-)MSGARCH model to data, we need to define the model specification through 
`CreateSpec`. This step is analogous to the basic version of the package with 
the key difference that in the `switch.spec` argument the variable `do.tvp` must be 
set to `TRUE`. 

``` r
# To fit the model with time-varying probabilities (TVP), we must set
# do.tvp = TRUE:
spec = CreateSpec(switch.spec = list(do.tvp=TRUE),
                  variance.spec = list(model = c("sGARCH", "sGARCH")),
                  distribution.spec = list(distribution = c("norm", "norm")))
``` 


``` r
# Constructing the covariate matrix Z
Z = matrix(1, nrow = length(data)-1, ncol = 2)
Z[,2] = data[1:2499]
data = data[2:2500]
```

Subsequently, we fit the model through the MLE method `FitML`. Note that the 
covariate matrix `Z` only needs to be supplied for time-varying switching. Thus, 
it only works if the object `spec` includes the correct specifications (`do.tvp = TRUE` and 
`do.mix = FALSE`). 
``` r
# We fit the model via the maximum likelihood method:
fit = FitML(spec = spec, data = data,  Z = Z)
```



### Which extensions were made?

The following extensions were made to the `MSGARCH` package:

* Incorporate time-varying transition probabilities into the ML estimation ✔️
* Enable one-step-ahead predictions for S3 method `predict()` ✔️
* Enable one-step-ahead predictions for method `PredPdf()`️ ✔️
* Enable one-step-ahead predictions for method `Risk()` ✔️
* Make ML-Fit object (estimated through time-varying switching) compatible for method `State()` ✔️
* S3 methods like `print()` and `summary()` indicate wheter time-varying MSGARCH was used ✔️
* `TransMat()` returns the initial transition matrix ✔️
* `ExtractStateFit()` is enabled for MLFit object which were estimated through time-varying switching ✔️
* S3 method `plot()` work for time-varying switching ✔️
* `Volatility()` works for time-varying switching ✔️
* `UncVol()` works for time-varying switching ✔️



The following functionalists are not yet implemented for time-varying switching:

* Time-varying Mixture modelling, i.e. `do.tvp=TRUE` and `do.mix=TRUE` is not yet supported ❌
* `PIT()` is not working for time-varying switching ❌
* `simulate.MSGARCH_SPEC()` simulate is work in progress ❌
* multi-step-ahead predictions are not supported for time-varying switching ❌

### Important Info
I want to re-emphasize that I created this fork solely for the purpose of my master's 
thesis. I don't want to make any pull-request to the original package and I am not planing 
on maintain this repo in the future. 

Contact:  

* [LinkedIn - Max Kuttner](https://www.linkedin.com/in/maximilian-kuttner/) 
* [Uni - WU Vienna](https://www.wu.ac.at/en/programs/masters-programs/quantitative-finance/overview)



-------------------------------------------------------------------------------------------------------

## Please cite the package in publications!

By using `MSGARCH` you agree to the following rules: 

1) You must cite [Ardia et al. (2019)](https://doi.org/10.18637/jss.v091.i04) in working papers and published papers that use `MSGARCH`.
2) You must place the following URL in a footnote to help others find `MSGARCH`: [https://CRAN.R-project.org/package=MSGARCH](https://CRAN.R-project.org/package=MSGARCH). 
3) You assume all risk for the use of `MSGARCH`.

**Ardia, D., Bluteau, K., Boudt, K., Catania, L., Trottier, D.-A. (2019).**    
Markov-switching GARCH models in R: The MSGARCH package.    
_Journal of Statistical Software_, 91(4), 1-38.    
[https://doi.org/10.18637/jss.v091.i04](https://doi.org/10.18637/jss.v091.i04)

## Other references

**Ardia, D., Bluteau, K., Boudt, K., Catania, L. (2018).**  
Forecasting risk with Markov-switching GARCH models: A large-scale performance study   
_International Journal of Forecasting_, 34(4), 733-747.                                               
[https://doi.org/10.1016/j.ijforecast.2018.05.004](https://doi.org/10.1016/j.ijforecast.2018.05.004)

**Ardia, D., Bluteau, K., Ruede, M. (2019).**    
Regime changes in Bitcoin GARCH volatility dynamics.    
_Finance Research Letters_, 29, 266-271.                                         
[https://doi.org/10.1016/j.frl.2018.08.009](https://doi.org/10.1016/j.frl.2018.08.009)
