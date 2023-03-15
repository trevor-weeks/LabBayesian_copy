Lab 9 Bayesian RSF
================
Josh Nowak, Mark Hebblewhite, Sarah Straughan
March 15, 2023

The code in this repo was created during a lab with WILD 562. The
sim_fit script shows the user one simple way to simulate data for a
binomial regression, in this case a resource selection function (RSF).
The script also allows the user to fit a basic regression considering
covariates (model_one.txt) and a second model with an individual random
effect (model_two.txt). For the sake of creating good habits a few
functions were written in the script to help with summarizing results,
but please recognize that there are established packages that accomplish
these same tasks better than what was written here
(<https://github.com/mjskay/tidybayes>).

The purpose of these scripts is to provide a simple entry point that
allows the user to become familiar with the simulated/fit workflow and
the querks of running an analysis in R. For those interested in fitting
RSFs in R I would consider reading the ecology and spatial task views in
R to get a feel for the types of analyses that are packaged for you. In
addition, those interested in Bayesian methods should consider
alternative ways to call the models such as rjags, rstan and jagsUI.

``` r
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#load or install these packages:
packages <- c("tidyverse", "mcmcplots", "R2jags", "ResourceSelection","purrr", "glmmTMB")
#run function to install packages
ipak(packages)
```

    ## Warning in checkMatrixPackageVersion(): Package version inconsistency detected.
    ## TMB was built with Matrix version 1.5.3
    ## Current Matrix version is 1.5.1
    ## Please re-install 'TMB' from source using install.packages('TMB', type = 'source') or ask CRAN for a binary version of 'TMB' matching CRAN's 'Matrix' package

    ##         tidyverse         mcmcplots            R2jags ResourceSelection 
    ##              TRUE              TRUE              TRUE              TRUE 
    ##             purrr           glmmTMB 
    ##              TRUE              TRUE

## Model 1 - no random effects

``` r
rsf_m1dat <- tibble::tibble(
  ndvi = rnorm(length(id)),
  pres = rbinom(length(id), size = 1, prob = plogis(0.5 + 0.3 * ndvi))
)

# Exploring priors
hist(plogis(rnorm(100000, 0, sqrt(1/0.001))), breaks = 100, col = "dodgerblue")
```

![](README_files/figure-gfm/Model%201-1.png)<!-- -->

``` r
hist(plogis(rnorm(100000, 0, sqrt(1/0.5))), breaks = 100, col = "dodgerblue")
```

![](README_files/figure-gfm/Model%201-2.png)<!-- -->

``` r
hist(plogis(runif(10000, -5, 5)))
```

![](README_files/figure-gfm/Model%201-3.png)<!-- -->

``` r
# Gather data for JAGS - must be named list
jdat_m1 <- list(
  NOBS = nrow(rsf_m1dat),

  NDVI = rsf_m1dat$ndvi,
  PRES = rsf_m1dat$pres
)

# Create initial values
jinits <- function(){
  list(
    alpha = rnorm(1),
    ndvi_eff = rnorm(1)
  )
}

# Parameters to monitor
params_m1 <- c("alpha", "ndvi_eff")
```

You will need to have JAGS downloaded before you can run this code.
Download [Here](http://www.sourceforge.net/projects/mcmc-jags/files).
JAGS stands for Just Another Gibbs Sampler and quoting the program
author, Martyn Plummer, “It is a program for analysis of Bayesian
hierarchical models using Markov Chain Monte Carlo (MCMC) simulation…”

``` r
#Call JAGS
fit_m1 <- jags(
  data = jdat_m1,
  inits = jinits,
  parameters.to.save = params_m1,
  model.file = "models/model_one.txt",
  n.chains = 3,
  n.burnin = 500,
  n.iter = 600,
  n.thin = 1
)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 1
    ##    Unobserved stochastic nodes: 2
    ##    Total graph size: 10
    ## 
    ## Initializing model

Here we write a function to summarise results, we write a function
because this action willvbe repeated multiple times

``` r
summ_fun <- function(x, param){
  tibble::tibble(
    Parameter = param,
    Mean = mean(x$BUGS$sims.list[[param]]),
    SD = sd(x$BUGS$sims.list[[param]]),
    LCL = quantile(x$BUGS$sims.list[[param]], probs = .025),
    UCL = quantile(x$BUGS$sims.list[[param]], probs = .975)
  )
}
```

Next, we create a function to determine if a parameter value is greater
than 0

``` r
grtr_zero <- function(x, param){
  sum(x$BUGS$sims.list[[param]] > 0)/length(x$BUGS$sims.list[[param]])
}

summ_fun(fit_m1, "alpha")
```

    ## # A tibble: 1 × 5
    ##   Parameter  Mean    SD   LCL   UCL
    ##   <chr>     <dbl> <dbl> <dbl> <dbl>
    ## 1 alpha      11.9  28.7 -38.0  71.3

``` r
summ_fun(fit_m1, "ndvi_eff")
```

    ## # A tibble: 1 × 5
    ##   Parameter  Mean    SD   LCL   UCL
    ##   <chr>     <dbl> <dbl> <dbl> <dbl>
    ## 1 ndvi_eff   23.1  22.9 -18.2  67.1

``` r
grtr_zero(fit_m1, "alpha")
```

    ## [1] 0.65

``` r
grtr_zero(fit_m1, "ndvi_eff")
```

    ## [1] 0.8433333

The purrr package implements some clean functions aimed at the tenants
of functional programming, here we loop over a series of inputs while
calling a function

``` r
purrr::map_df(c("alpha", "ndvi_eff"), ~summ_fun(fit_m1, .x))
```

    ## # A tibble: 2 × 5
    ##   Parameter  Mean    SD   LCL   UCL
    ##   <chr>     <dbl> <dbl> <dbl> <dbl>
    ## 1 alpha      11.9  28.7 -38.0  71.3
    ## 2 ndvi_eff   23.1  22.9 -18.2  67.1

``` r
#  More generic
purrr::map_df(params_m1, ~summ_fun(fit_m1, .x))
```

    ## # A tibble: 2 × 5
    ##   Parameter  Mean    SD   LCL   UCL
    ##   <chr>     <dbl> <dbl> <dbl> <dbl>
    ## 1 alpha      11.9  28.7 -38.0  71.3
    ## 2 ndvi_eff   23.1  22.9 -18.2  67.1

The mcmcplots package has several useful utilities to help with
assessing convergence and examining model outputs

``` r
mcmcplots::mcmcplot(fit_m1)
```

    ##                                                                                 Preparing plots for alpha.  33% complete.

    ##                                                                                 Preparing plots for deviance.  67% complete.

    ##                                                                                 Preparing plots for ndvi.  100% complete.

## Model 2 - Bayesian Model with Random Intercept for Individual Elk

``` r
rsf_m2dat <- tibble::tibble(
  id = rep(1:3, each = 5),
  ndvi = rnorm(length(id)),
  pres = rbinom(length(id), size = 1, prob = plogis(0.5 + 0.3 * ndvi))
)

# Exploring priors
hist(plogis(rnorm(100000, 0, sqrt(1/0.001))), breaks = 100, col = "dodgerblue")
```

![](README_files/figure-gfm/Model%202-1.png)<!-- -->

``` r
hist(plogis(rnorm(100000, 0, sqrt(1/0.5))), breaks = 100, col = "dodgerblue")
```

![](README_files/figure-gfm/Model%202-2.png)<!-- -->

``` r
hist(plogis(runif(10000, -5, 5)))
```

![](README_files/figure-gfm/Model%202-3.png)<!-- -->

``` r
# Gather data for JAGS - must be named list
jdat_m2 <- list(
  NOBS = nrow(rsf_m2dat),
  NIND = n_distinct(rsf_m2dat$id),

  IND = rsf_m2dat$id,
  NDVI = rsf_m2dat$ndvi,
  PRES = rsf_m2dat$pres
)


# Parameters to monitor
params_m2 <- c("alpha", "ndvi_eff", "ind_eff", "sd_ind")

# Call JAGS
fit_m2 <- jags(
  data = jdat_m2,
  inits = jinits,
  parameters.to.save = params_m2,
  model.file = "models/model_two.txt",
  n.chains = 3,
  n.burnin = 500,
  n.iter = 600,
  n.thin = 1
)
```

    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 15
    ##    Unobserved stochastic nodes: 6
    ##    Total graph size: 105
    ## 
    ## Initializing model

``` r
summ_fun(fit_m2, "alpha")
```

    ## # A tibble: 1 × 5
    ##   Parameter    Mean    SD   LCL   UCL
    ##   <chr>       <dbl> <dbl> <dbl> <dbl>
    ## 1 alpha     -0.0549  1.47 -3.62  2.71

``` r
summ_fun(fit_m2, "ndvi_eff")
```

    ## # A tibble: 1 × 5
    ##   Parameter  Mean    SD    LCL   UCL
    ##   <chr>     <dbl> <dbl>  <dbl> <dbl>
    ## 1 ndvi_eff  0.860 0.797 -0.717  2.37

``` r
grtr_zero(fit_m2, "alpha")
```

    ## [1] 0.4933333

``` r
grtr_zero(fit_m2, "ndvi_eff")
```

    ## [1] 0.88

``` r
purrr::map_df(c("alpha", "ndvi_eff"), ~summ_fun(fit_m2, .x))
```

    ## # A tibble: 2 × 5
    ##   Parameter    Mean    SD    LCL   UCL
    ##   <chr>       <dbl> <dbl>  <dbl> <dbl>
    ## 1 alpha     -0.0549 1.47  -3.62   2.71
    ## 2 ndvi_eff   0.860  0.797 -0.717  2.37

``` r
#  More generic
purrr::map_df(params_m2, ~summ_fun(fit_m2, .x))
```

    ## # A tibble: 4 × 5
    ##   Parameter    Mean    SD     LCL   UCL
    ##   <chr>       <dbl> <dbl>   <dbl> <dbl>
    ## 1 alpha     -0.0549 1.47  -3.62    2.71
    ## 2 ndvi_eff   0.860  0.797 -0.717   2.37
    ## 3 ind_eff    0.0165 1.48  -3.15    3.42
    ## 4 sd_ind     1.59   1.54   0.0127  5.54

``` r
mcmcplots::mcmcplot(fit_m2)
```

    ##                                                                                 Preparing plots for alpha.  14% complete.

    ##                                                                                 Preparing plots for deviance.  29% complete.

    ##                                                                                 Preparing plots for ind.  43% complete.

    ##                                                                                 Preparing plots for ind.  57% complete.

    ##                                                                                 Preparing plots for ind.  71% complete.

    ##                                                                                 Preparing plots for ndvi.  86% complete.

    ##                                                                                 Preparing plots for sd.  100% complete.

Also check out the Bayesian task view in R and tidybayes in particular

## RSF analysis of mountain goats (Section 4.1)

Authors: S. Muff, J. Signer, J. Fieberg

\-[Link to
Paper](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13087)
-[Link to
Code](https://conservancy.umn.edu/bitstream/handle/11299/204737/Goats_RSF.R?sequence=20&isAllowed=y)

To install the INLA-package in R, you have to manually add the r-inla
repository as they are not on CRAN. You may have to restart your R
session before installation – this is the type of thing where in the
afternoon the install code did not work for me and the next morning it
worked – the only difference being that I had restarted R in between and
ran the following code before anything else:

Load libraries and data

``` r
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
```

``` r
library(INLA)
data(goats)
str(goats)
```

    ## 'data.frame':    19014 obs. of  8 variables:
    ##  $ STATUS   : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ ID       : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ ELEVATION: int  651 660 316 334 454 343 429 493 400 442 ...
    ##  $ SLOPE    : num  38.5 39.7 20.5 34.1 41.6 ...
    ##  $ ET       : num  35.4 70.7 50 35.4 25 ...
    ##  $ ASPECT   : num  243 270 279 266 258 ...
    ##  $ HLI      : num  0.918 0.884 0.713 0.864 0.935 ...
    ##  $ TASP     : num  0.947 0.699 0.575 0.745 0.829 ...

``` r
#Scale and center variables
goats$ELEVATION <- scale(goats$ELEVATION)
goats$TASP <- scale(goats$TASP)

#Use and available data by animal
with(goats, prop.table(table(ID, STATUS), 1))
```

    ##     STATUS
    ## ID           0         1
    ##   1  0.6666667 0.3333333
    ##   2  0.6666667 0.3333333
    ##   3  0.6666667 0.3333333
    ##   4  0.6666667 0.3333333
    ##   5  0.6666667 0.3333333
    ##   6  0.6666667 0.3333333
    ##   7  0.6666667 0.3333333
    ##   8  0.6666667 0.3333333
    ##   9  0.6666667 0.3333333
    ##   10 0.6666667 0.3333333

### Model M1: using `glmmTMB()`

- Fixed Effects: elevation
- Random Effects: intercept only

``` r
goats.M1 <- glmmTMB(STATUS ~ ELEVATION  + (1|ID), family=binomial(), data = goats)
summary(goats.M1)
```

    ##  Family: binomial  ( logit )
    ## Formula:          STATUS ~ ELEVATION + (1 | ID)
    ## Data: goats
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  24199.8  24223.4 -12096.9  24193.8    19011 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups Name        Variance Std.Dev.
    ##  ID     (Intercept) 0.007625 0.08732 
    ## Number of obs: 19014, groups:  ID, 10
    ## 
    ## Conditional model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.68924    0.03192 -21.592   <2e-16 ***
    ## ELEVATION    0.11844    0.04634   2.556   0.0106 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Model M2

- Fixed Effects: elevation, aspect
- Random Effects: intercept only

``` r
goats.M2 <- glmmTMB(STATUS ~ ELEVATION + TASP + (1|ID), family=binomial(), data = goats)
summary(goats.M2)
```

    ##  Family: binomial  ( logit )
    ## Formula:          STATUS ~ ELEVATION + TASP + (1 | ID)
    ## Data: goats
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  23331.1  23362.6 -11661.6  23323.1    19010 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups Name        Variance Std.Dev.
    ##  ID     (Intercept) 0.01303  0.1142  
    ## Number of obs: 19014, groups:  ID, 10
    ## 
    ## Conditional model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.74221    0.03970 -18.694  < 2e-16 ***
    ## ELEVATION    0.14362    0.02579   5.568 2.58e-08 ***
    ## TASP         0.51913    0.01883  27.566  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Model M3

- Fixed Effects: elevation, aspect
- Random Effects: intercepts and slopes (elevation, aspect)

``` r
goats.M3 <- glmmTMB(STATUS ~ TASP +  ELEVATION + (1|ID) + (0+ELEVATION |ID) + (0+TASP|ID), family=binomial(), data = goats)
summary(goats.M3)
```

    ##  Family: binomial  ( logit )
    ## Formula:          
    ## STATUS ~ TASP + ELEVATION + (1 | ID) + (0 + ELEVATION | ID) +  
    ##     (0 + TASP | ID)
    ## Data: goats
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  21978.5  22025.7 -10983.3  21966.5    19008 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups Name        Variance Std.Dev.
    ##  ID     (Intercept) 0.9572   0.9783  
    ##  ID.1   ELEVATION   1.4011   1.1837  
    ##  ID.2   TASP        0.1043   0.3229  
    ## Number of obs: 19014, groups:  ID, 10
    ## 
    ## Conditional model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.39285    0.31158  -1.261    0.207    
    ## TASP         0.66416    0.10565   6.286 3.25e-10 ***
    ## ELEVATION    0.06934    0.37622   0.184    0.854    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Model M4 (with fixed intercept variance)

- Fixed Effects: elevation
- Random Effects: intercept (with large fixed variance), elevation,
  aspect

Here, we also use a weighted likelihood. To this end, we need to create
a variable for the weights, where used points (`STATUS=1`) keep weight
1, and available points (`STATUS=0`) obtain a large weight $W$ (here
$W=1000$):

``` r
goats$weight <- 1000^(1-goats$STATUS)
```

We fit the same model as under M3, again using `glmmTMB`. Note that we
have to manually fix the variance of the intercept first. Start by
setting up the model, but do not yet fit it:

``` r
goats.M4.tmp <- glmmTMB(STATUS ~   ELEVATION + TASP + (1|ID) + (0+ELEVATION |ID) + (0+TASP|ID), family=binomial(), data = goats,doFit=F, weights = weight)
```

Then fix the standard deviation of the first random term, which is the
`(1|ID)` component in the above model equation. We use $\sigma=10^3$,
which corresponds to a variance of $10^6$:

``` r
goats.M4.tmp$parameters$theta[1] = log(1e3)
```

We need to tell `glmmTMB` not to change the first entry of the vector of
variances, and give all other variances another indicator to make sure
they can be freely estimated:

``` r
goats.M4.tmp$mapArg = list(theta=factor(c(NA,1:2)))
```

Then fit the model and look at the results:

``` r
goats.M4 <- glmmTMB:::fitTMB(goats.M4.tmp)
summary(goats.M4)
```

    ##  Family: binomial  ( logit )
    ## Formula:          
    ## STATUS ~ ELEVATION + TASP + (1 | ID) + (0 + ELEVATION | ID) +  
    ##     (0 + TASP | ID)
    ## Data: goats
    ## Weights: weight
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## 105999.6 106038.9 -52994.8 105989.6    19009 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups Name        Variance  Std.Dev. 
    ##  ID     (Intercept) 1.000e+06 1000.0000
    ##  ID.1   ELEVATION   9.257e-01    0.9621
    ##  ID.2   TASP        1.221e-01    0.3495
    ## Number of obs: 19014, groups:  ID, 10
    ## 
    ## Conditional model:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -3.633e-05  3.162e+02   0.000    1.000    
    ## ELEVATION    1.172e-01  3.055e-01   0.384    0.701    
    ## TASP         6.505e-01  1.127e-01   5.770 7.92e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Model M4 (with intercept variance estimated)

For comparison, we again fit model M4, but without fixing the intercept
variance, letting it be estimated instead. Importantly, estimating the
intercept variance is the current standard procedure. For this
particular RSF case, it does not lead to a real difference, as expected
due to the many observations per individual. This confirms that the
decision to fix or estimate the intercept variance is not critical for
RSFs, in contrast to SSFs (see Discussion in the paper).

``` r
goats.M4.2 <- glmmTMB(STATUS ~   ELEVATION + TASP + (1|ID) + (0+ELEVATION |ID) + (0+TASP|ID), family=binomial(), data = goats, weights = weight)
summary(goats.M4.2)
```

    ##  Family: binomial  ( logit )
    ## Formula:          
    ## STATUS ~ ELEVATION + TASP + (1 | ID) + (0 + ELEVATION | ID) +  
    ##     (0 + TASP | ID)
    ## Data: goats
    ## Weights: weight
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ## 105869.2 105916.3 -52928.6 105857.2    19008 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups Name        Variance Std.Dev.
    ##  ID     (Intercept) 0.6448   0.8030  
    ##  ID.1   ELEVATION   0.9151   0.9566  
    ##  ID.2   TASP        0.1208   0.3476  
    ## Number of obs: 19014, groups:  ID, 10
    ## 
    ## Conditional model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  -7.3852     0.2554 -28.918  < 2e-16 ***
    ## ELEVATION     0.1177     0.3037   0.387    0.698    
    ## TASP          0.6501     0.1121   5.797 6.73e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## INLA (only model M4)

Let us now carry the analysis of model M4 with random intercept
$\mathsf{N}(0,\sigma_{ID}^2)$ and fixed variance $\sigma_{ID}^2=10^6$
using INLA. A peculiarity of INLA is that the same variable cannot be
used more than once. So for ID we need to generate two new (but
identical) variables:

``` r
goats$ID2 <- goats$ID3 <- goats$ID
```

For the fixed effects we use the INLA (default) priors
$\beta \sim \mathsf{N}(0,\sigma_\beta^2)$ with $\sigma_\beta^2=10^4$.
The precisions of the priors are thus set to:

``` r
prec.beta.TASP  <- 1e-4
prec.beta.ELEVATION  <- 1e-4
```

The INLA formula with the fixed effects `TASP` and `ELEVATION`, plus
three random effects: one for the individual-specific intercept, and two
random slopes for `TASP` and `ELEVATION`. Note that the precision (thus
$1/\sigma^2$) for `ID` is fixed (`fixed=T`) at the value of $10^{-6}$
(thus the variance is fixed at $10^6$). The precisions for the random
slopes for `TASP` and `ELEVATION` are given PC(1,0.05) priors:

``` r
formula.inla <-STATUS ~  TASP  + ELEVATION +
  f(ID,model="iid",hyper=list(theta = list(initial=log(1e-6),fixed=T))) +
  f(ID2,TASP,values=1:10,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(1,0.05)))) +
  f(ID3,ELEVATION,values=1:10,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(1,0.05))))
```

The actual INLA call is then given as follows:

``` r
#inla.setOption(enable.inla.argument.weights=TRUE) #this line did not work but rest of code worked without it??
goats.M4.inla  <- inla(formula.inla, family ="binomial", data=goats, weights=goats$weight,
                       control.fixed = list(
                         mean = 0,
                         prec = list(TASP = prec.beta.TASP,
                                     ELEVATION = prec.beta.ELEVATION)
                       )
)
```

The summary for the posterior distribution of the fixed effects is given
as follows:

``` r
goats.M4.inla$summary.fixed
```

    ##                   mean          sd   0.025quant   0.5quant  0.975quant
    ## (Intercept) -7.3828278 316.2253265 -627.5495303 -7.3828278 612.7838747
    ## TASP         0.6485443   0.1269227    0.3979946  0.6477524   0.9032699
    ## ELEVATION    0.1169434   0.3098050   -0.4977457  0.1168881   0.7319257
    ##                   mode          kld
    ## (Intercept) -7.3828278 5.527851e-11
    ## TASP         0.6461330 1.012595e-05
    ## ELEVATION    0.1167736 1.323423e-06

Since variances are parameterized and treated as precisions, the summary
of the respective posterior distributions is given for the precisions:

``` r
goats.M4.inla$summary.hyperpar
```

    ##                       mean        sd 0.025quant 0.5quant 0.975quant      mode
    ## Precision for ID2 7.617567 3.7264495  2.5148791 6.918226  16.832045 5.5717052
    ## Precision for ID3 1.179362 0.4827767  0.4716325 1.102258   2.340417 0.9527702

Source R functions for calculating posterior means and medians of the
precisions. Not currently functioning

``` r
source("inla_emarginal.R")
source("inla_mmarginal.R")
inla_emarginal(goats.M4.inla)
inla_mmarginal(goats.M4.inla)
```
