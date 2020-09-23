# SLTCA

[![CRAN Status Badge](http://www.r-pkg.org/badges/version/SLTCA)](http://cran.r-project.org/web/packages/SLTCA)
[![Downloads badge](https://cranlogs.r-pkg.org/badges/SLTCA)](https://cranlogs.r-pkg.org/badges/SLTCA)
[![Total downloads](https://cranlogs.r-pkg.org/badges/grand-total/SLTCA)](https://cranlogs.r-pkg.org/badges/grand-total/SLTCA)
[![](https://img.shields.io/badge/doi-10.1111/biom.13366-blue.svg)](https://doi.org/10.1111/biom.13366)

SLTCA: Scalable and Robust Latent Trajectory Class Analysis Using Artificial Likelihood

# News

The package is on [CRAN](https://cran.r-project.org/package=SLTCA) now. This repository will mainly serve as a platform for bug reporting (at [Issues](https://github.com/tengfei-emory/SLTCA/issues)) while we will post the lastest version on CRAN. To access the latest version on GitHub, please go to branch [CRAN](https://github.com/tengfei-emory/SLTCA/tree/CRAN). Thanks for considering using our software!

# Installation Guide
```{r}
library(devtools)
install_github('tengfei-emory/SLTCA')
library(SLTCA)
```
Currently `SLTCA` supports R version >= 3.3.0.

# Example: analyze a simulated dataset

## Data simulation

By default, the function `simulation(n)` generates a dataset with n observations under the scenario 1 described by Hart, Fei and Hanfelt (2020). 
```{r}
# generate a dataset with 500 individuals
dat <- simulation(500)
```
Specifically, it returns a data frame of 2 latent classes with 6 longitudinal features `y.1` to `y.6`, including count (`y.1` and `y.2`), binary (`y.3` and `y.4`) and continuous (`y.5` and `y.6`) features. The data frame also consists of individual identifiers (`id`), corresponding time of longitudinal features (`time`) and the number of visit (`num_obs`). In addition, variable `baselinecov` is a binary baseline risk factor of latent classes. Variable `latent` is the true latent class labels.

## Model fitting

The analysis for the dataset `dat` can be conducted by running `SLTCA` function:

```{r}
res <- SLTCA(k=1,dat=dat,num_class=2,id="id",time="time",num_obs="num_obs",features=paste("y.",1:6,sep=''),
             Y_dist=c('poi','poi','bin','bin','normal','normal'),
             covx="baselinecov",ipw=1,stop="tau",tol=0.005,max=50,
             varest=T,balanced=T,MSC='EQIC',verbose=T)
```
Please refer to the function documentation for more details.

# References

Hart, Fei and Hanfelt (2020), [Scalable and Robust Latent Trajectory Class Analysis Using Artificial Likelihood](https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.13366). Biometrics, Accepted Author Manuscript. doi:10.1111/biom.13366
