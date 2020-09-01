# SLTCA
SLTCA: Scalable and Robust Latent Trajectory Class Analysis Using Artificial Likelihood

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

Hart, Fei and Hanfelt (2020), Scalable and Robust Latent Trajectory Class Analysis Using Artificial Likelihood. Biometrics, accepted.
