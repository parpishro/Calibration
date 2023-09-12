
<!-- badges: start -->
[![R-CMD-check](https://github.com/parpishro/FBC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/parpishro/FBC/actions/workflows/R-CMD-check.yaml)

[![codecov](https://codecov.io/gh/parpishro/FBC/branch/main/graph/badge.svg)](https://codecov.io/gh/parpishro/FBC)

[![CRAN
status](https://www.r-pkg.org/badges/version/FBC)](https://CRAN.R-project.org/package=FBC)

[![Dependencies](https://img.shields.io/badge/dependencies-0/0-brightgreen?style=flat)](#)

[![License:
GPL-3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://choosealicense.com/licenses/gpl-3.0/)

[![LifeCycle](https://img.shields.io/badge/lifecycle-experimental-orange)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

[![Project Status:
Concept](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
<!-- badges: end -->

------------------------------------------------------------------------

# FBC: Full Bayesian Calibration

`FBC` is a package that uses the data from both field experiment and
computer simulation to calibrate the simulator model [^1]. A field
experiment relates a set of experimental variables as inputs to a
response as output. A computer simulation also relates the same
experimental variables as inputs to the response as output but also
include calibration inputs. These inputs either represent unknown but
fixed physical properties that are governed by the physical system in
field experiments (and therefore need not be specified in field
experiments) or represent various aspects of the simulation model such
as tuning hyperparameter of the model. In both cases, the calibration
parameters must be estimated for a simulator model to mimic the physical
system adequately. The goal of the calibration is to estimate
calibration parameters by fitting a statistical model to the observed
data using both field and simulator data. Estimated parameters of the
calibration model can be used for inference or calibrated prediction.

The `FBC` package is a based on the well-known Kennedy and O’Hagan (KOH)
calibration model (2001). In the KOH model, field and simulator data are
stacked in a single model, where simulator response is fitted using a
Gaussian Process (GP) that models the simulator model and field response
is fitted using a three-component model that accounts for simulator
model using simulator GP, bias-correction using another independent GP,
and a Gaussian random error term representing measurement error. Using
two independent GPs and a Gaussian random error introduces new
hyperparameters to the calibration model that must be estimated along
with calibration parameters. In its original formulation, KOH formulates
a hierarchical Bayesian framework with two sequential phases. In the
first phase, model hyperparameters are predicted using maximum
likelihood estimation (MLE) method. In the second phase, Bayesian
analysis is employed to derive the posterior distribution of calibration
parameters while fixing the hyperparameters found in the first phase.

Neither KOH’s original formulation nor later suggestions are fully
Bayesian as they use MLE methods to estimate hyperparameters and fix the
hyperparameters in the second phase which ignores the added uncertainty
of the estimated hyperparameters, largely due to computational
infeasibility (Kennedy & O’Hagan, 2001; Higdon et al., 2004; Liu et al.,
2009). On the other hand, `FBC` runs a fully Bayesian model that
includes both calibration parameters and model hyperparameters in its
Bayesian framework. Implementation of the package optimizes the
calibration process and memory management to increase computational
efficiency. Moreover, the fully implemented Bayesian framework, enables
the user to apply the expert knowledge about any of the model parameters
using prior specifications. All common prior distribution are
implemented in `FBC` to allows for high degree of flexibility in
specifying the prior information or the lack thereof. `FBC` runs a
variation of Markov Chain Monte Carlo (MCMC) algorithm called
Metropolis-Within-Gibbs algorithm to find the posterior distribution of
calibration parameters and model hyperparameters. Furthermore, `FBC`
enables calibrated prediction for new input configurations using the
calibration model.

The current vignette is structured into following sections: The first
section, Using `FBC`, explains the package functionality through a
simple pedagogic example. Also in this section, the notation for inputs
and outputs of both computer and physical experiments are introduced. In
the second section, Calibration Model, generalizes the example
introduced in the first section to build model components and shows how
a calibration model based on KOH model is built internally. Using the
general notation, while referencing the example, modelling choices are
justified. In the third section, Parameter Estimation, the theoretical
results to characterize the posterior distribution of model parameters
are presented along with procedure for MCMC-based estimation of
parameters. In the fourth section, Calibrated Prediction, the
theoretical results to derive calibrated predictions for new input
configuration using the estimated parameters are offered. In the fifth
section, Implementation, some of the implementation features and choices
of `FBC` package are explained along with its limitations. In the sixth
and last section, Applications, three more examples are presented to
demonstrate the full functionality and limitations of the `FBC` package.
And finally, the Appendix provides further details to some of the
concepts presented in the sections.

`FBC` package has two main public functions: `calibrate()` and
`predict()`. As function names suggest, `calibrate()` takes the field
and simulation training data to calibrate the simulator model and
`predict` takes a new field input configuration and the calibration
model object to predict the field response and its associated
uncertainty. In addition, `FBC` has four more public functions to help
with prior specifications, summarizing, and visualization of
calibration. Furthermore, there are three more public function that are
not directly used by user to build calibration model but provide
functionalities that are not offered by base R. This section explains
how to use the mentioned functions and what to expect in their use and
results using a simple pedagogic example. Throughout this section, the
user is expected to be well-versed with calibration model in general and
KOH model in particular. It is designed to be a self-contained section
on the usage of `FBC` package but can be augmented with section six,
Applications, to further examine the functionalities and limitation of
the package. Less experienced users are encouraged to start from next
section which aims to explain the theory behind calibration models and
calibrated predictions in more details.

### Setup

The simplest and safest method to obtain `FBC` package is through CRAN
using following command.

``` r
install.packages(FBC)
```

Alternatively, the development version of the `FBC` package can be
installed directly from Github using `devtools` package. Note that,
current vignette is based on published version of the package on CRAN
and the development version may contain further functionalities.

``` r
devtools::install_github("parpishro/FBC")
```

After installing the package, it must be loaded into the R session.

``` r
library(FBC)
```

### Data

Building a calibration model requires data from both field experiment
and computer simulation. To focus on the functionality of the package,
we use a simple pedagogic example as experimental setting. In the field
setting, a wiffle ball is dropped from different heights and the time it
takes to hit the ground is measured. The experimental input is height
($h$). and the response is time ($t$). In the simulation setting, the
physical process of a falling ball is modelled by a mathematical
equation that relates $h$ to $t$, however, one must also provide gravity
($g$) as calibration input. Note that in the field experiment, the
earth’s gravity is fixed but unknown to the experimenter [^2]. The field
experiment has been performed by Derek Bingham and Jason Loeppky and it
is provided by Robert Gramacy in his book “surrogates” [^3]. The
mathematical model used in simulation of the ball drop experiment can be
easily derived from introductory physics results.

$$
t = \sqrt{\frac{2h}{g}}
$$

The simulation code takes $h$ and $g$ are taken as experimental and
calibration inputs and returns $t$ as simulator response based on above
mathematical equation. A common method to choose the different
combination of the inputs from acceptable ranges is Latin Hypercube
Sampling (LHS) method (McKay et al., 1979) [^4]. The `FBC` package
requires the training data to be in matrix format, where for both
matrices, the first column always represents response vector and
following columns represent experimental inputs (for both data matrices)
and calibration inputs (only for simulation data matrix). The data
matrices for ball example are produced in advance and loaded into the
package environment under the name of `ballField` and `ballSim` [^5].

Below, the structures and dimensions of both data matrices can be
inspected.

``` r
head(ballField, 3)
```

    ##         t     h
    ## [1,] 0.27 0.178
    ## [2,] 0.22 0.356
    ## [3,] 0.27 0.534

``` r
dim(ballField)
```

    ## [1] 63  2

``` r
head(ballSim, 3)
```

    ##          t     h      g
    ## [1,] 0.404 0.998 12.220
    ## [2,] 0.487 0.940  7.909
    ## [3,] 0.450 0.886  8.736

``` r
dim(ballSim)
```

    ## [1] 100   3

Figure 1 shows the distribution of time versus height for both physical
and simulation experiments. Note that for higher height values, the
simulation responses (blue) underestimate their corresponding field
response (red), displaying a systemic bias for simulation model.

<!-- For HTML output -->

<div class="figure-box">

**Figure 1:** *The time versus height plot for both field experiment and
simulation.*

![](README_files/figure-gfm/ball%20plot%20html-1.png)<!-- -->

</div>

The ball example is a simple pedagogical example as it has only one
experimental input and one calibration input. Using this toy example, we
demonstrate the functionality of the package without getting to the
details of a complex mathematical model. More complex and real-world
examples are covered in the last section.

### Building a Calibration Model

The `calibrate()` function takes three sets of arguments from user:
data, MCMC, and prior specification parameters. Other than data
arguments which must be supplied by user, all other arguments have
reasonable default values for ball example [^6].

``` r
calMod <- calibrate(sim = ballSim, field = ballField,                                 # Data 
                    nMCMC = 11000, nBurn = 1000, thinning = 50,                       # MCMC
                    kappaDist = "beta", kappaP1 = 1.1, kappaP2 = 1.1, # Priors
                    hypers = set_hyperPriors(),
                    showProgress = FALSE) 
```

The first and second arguments `sim` and `field` must be supplied by
user. As mentioned earlier, both must be either in matrix format
representing the simulation and field data respectively. In both
`ballSim` and `ballField` of our ball example, the first column
represent the response $\bf{t}$ and the second column represent
experimental input $\bf{h}$. Additionally, `ballSim` has a third column
that represents calibration input $\bf{g}$. Calibration inputs are
implicit in field experiment and are absent in the data matrix
`ballField`.

$$
\begin{aligned}
&\quad \quad \text{Response}  \quad  \text{Experimental} \quad \text{Calibration} \\
\texttt{ballSim} \quad &= \
\begin{bmatrix}
\quad \quad \bf{t} & \quad \quad \quad \quad \bf{h} & \quad \quad \quad \quad \quad \bf{g} \quad \quad   
\end{bmatrix} \\ \\
\texttt{ballField} &= \
\begin{bmatrix}
\quad \quad \bf{t} & \quad \quad \quad \quad \bf{h} & \quad \quad  
\end{bmatrix}
\end{aligned}
$$

The second set of arguments consists of MCMC parameters: `nMCMC`
represents the number of total MCMC iterations , `nBurn` represent the
number of burn-in iterations to be removed from the beginning of the
chain, and `thinning` indicate the sampling rate to remove the
autocorrelation from the sampled draws. For example, when
`thinning = 50`, for every 50 draws from the result only one will be
kept in order to remove the autocorrelation between draws.

The third and last set of arguments consists of prior specification for
each parameter of the model. The main goal of calibration is to estimate
the calibration parameters ($g$ in ball example). However, employing KOH
calibration model introduces eight additional classes of hyperparameters
to the model. The parameters in all eight classes are not known in
advance and must be estimated. As `FBC` employs a full Bayesian
framework, the priors for all model parameters must be specified in
advance. Throughout the package implementation and current guide a
consistent notation is used to denote parameters: $\bf{\kappa}$,
$\bf{\theta_s}$, $\bf{\alpha_s}$, $\bf{\theta_b}$, $\bf{\alpha_b}$,
$\sigma^2_s$, $\sigma^2_b$, $\sigma^2_e$, and $\mu_b$ denote calibration
parameters, correlation scale parameters of simulator GP, correlation
smoothness parameters of simulator GP, correlation scale parameters of
bias-correction GP, correlation smoothness parameters of bias-correction
GP, marginal variance of simulator GP, marginal variance of
bias-correction GP, measurement error variance, and mean of
bias-correction GP respectively.

For each class of parameters, there are four associated arguments. 1)
distribution type, which is a character string determining the prior
distribution family (suffixed with “dist” in argument names). Currently,
`FBC` supports almost all common distributions and can be chosen from
“uniform”, “normal”, “normalTr”, “lognormal”, “gamma”, “inversegamma”,
“beta”, “betashift”, “logbeta”, “logistic”,“exponential”, “fixed” [^7].
2) Initial value, which is double representing starting point of the
parameter in MCMC algorithm . 3) First distribution parameter (suffixed
with “p1” in argument names) and 4) second distribution parameter
(suffixed with “p2” in argument names) are doubles representing the two
parameters of the chosen distribution [^8], [^9]. Table 1 summarizes all
classes of parameters along with their corresponding argument names for
prior specifications.

**Table 1:** *Argument names to specify priors for each parameter
class.*

| Parameter Class                                    | Distribution  | Initial Value | Shape Parameters          |
|:---------------------------------------------------|:--------------|:--------------|:--------------------------|
| True field calibration inputs ($\kappa$)           | `kappaDist`   | `kappaInit`   | `kappaP1` & `kappaP2`     |
| Simulation GP scale ($\theta_s$)                   | `thetaSDist`  | `thetaSInit`  | `thetaSP1` & `thetaSP2`   |
| Simulation GP smoothness ($\alpha_s$)              | `alphaSDist`  | `alphaSInit`  | `alphaSP1` & `alphaSP2`   |
| Bias-correction GP scale ($\theta_b$)              | `thetaBDist`  | `thetaBInit`  | `thetaBP1` & `thetaBP2`   |
| Bias-correction GP smoothness ($\alpha_b$)         | `alphaBDist`  | `alphaBInit`  | `alphaBP1` &`alphaBP2`    |
| Simulation marginal variance ($\sigma^2_s$)        | `sigma2SDist` | `sigma2SInit` | `sigma2SP1` & `sigma2SP2` |
| Bias-correction marginal variance ($\sigma^2_b$)   | `sigma2BDist` | `sigma2BInit` | `sigma2BP1` & `sigma2BP2` |
| Measurement error variance ($\sigma^2_{\epsilon}$) | `sigma2EDist` | `sigma2EInit` | `sigma2EP1` & `sigma2EP2` |
| Bias-correction mean ($\mu_b$)                     | `muBDist`     | `muBInit`     | `muBP1` & `muBP2`         |

From the classes of parameters mentioned, the prior for calibration
parameters should be specified by user based on field knowledge as these
parameters are problem-specific [^10]. In contrast, all other parameters
have reasonable default values based on KOH calibration model literature
(Kennedy & O’Hagan, 2001; Higdon et al. 2004; Chen et al. 2017). For
this reason and to avoid unwanted complexity, the priors for all other
classes of parameters are specified with `hypers` argument and using a
helper function `set_hyperPriors()`. In particular, the default value of
`hypers` argument is `set_hyperPriors()` without any argument.
Nevertheless, expert knowledge about one or more of these parameters can
be supplied using arguments of the `set_hyperPriors()` function, which
are defined in Table 1, to change the default prior specifications.

Note that some classes of parameters may contain more than one
parameters. In this case, the argument values can be vector instead of
default scaler. In this case all four fields of that parameter class
must be also in vector format with same length as number of parameters
in the class. For example, if there are five calibration parameters,
`kappaDist` can either be a scaler, in which case the distribution for
all calibration parameters will be set to that scaler value and
`kappaInit`, `kappaP1`, and `kappaP2` must be also scaler, or can be a
vector of length five that supplies the distribution types for all
calibration parameters and `kappaInit`, `kappaP1`, and `kappaP2` must
also be vector of length 5.

Moreover, there is a logical argument `showProgress` that indicate
whether function must show the progress in calibration on console. This
will put the `calibrate()` in interactive mode and will show the
percentage of the MCMC draws along with sample draws [^11].

The output of the `calibrate()` function is an object of class `fbc`
that contains the samples from posterior joint distribution of
parameters, along with other model information. Below, the components of
the a `fbc` object is displayed.

``` r
names(calMod)
```

    ## [1] "Phi"        "estimates"  "logPost"    "priors"     "acceptance"
    ## [6] "vars"       "data"       "scale"      "indices"

The first and main component of the `calMod` is matrix `Phi` whose
columns represent the sample of posterior densities for each unknown
parameter of the model in the same order as parameter classes in
table 1. Since $\kappa$, $\theta_s$, $\alpha_s$, $\theta_b$, $\alpha_b$
parameter classes may contain more than one parameter, they are suffixed
by a number that represents the index of the parameter in the class. For
example, if there are 3 calibration inputs, the first 3 columns of the
matrix `Phi` represent posterior density of calibration parameters and
the column headers will be `kappa1`, `kappa2`, and `kappa3`. In the ball
example, there is only one calibration parameter $g$, which is denoted
by $\kappa_1$ and represented by `kappa1` in matrix `Phi`. Each row of
the matrix `Phi` represents a MCMC draw.

``` r
head(calMod$Phi, 3)
```

    ##   kappa1 thetaS1 thetaS2 alphaS1 alphaS2 thetaB1 alphaB1 sigma2S sigma2B
    ## 1 0.5783  0.6744  0.1233  1.9794  1.9607  0.6394  1.6308  3.2139  0.1689
    ## 2 0.1822  0.6744  0.1160  1.9826  1.9743  0.6493  1.5982  4.1193  0.1662
    ## 3 0.3921  1.0425  0.1449  1.9851  1.9850  0.9704  1.5978  2.8183  0.1150
    ##   sigma2E     muB
    ## 1  0.0985 -0.2902
    ## 2  0.0984  0.9832
    ## 3  0.0984  0.2554

Table 2 provides an overview of the model parameters and the notation to
represent them in ball example. Later sections will explain in detail
why are these hyperparameters introduced, what do they represent, and
how they are estimated.

**Table 2:** Notation used in matrix `Phi` to represent parameters in
ball example.

| Column    | Notation              | Description                                                       |
|:----------|:----------------------|:------------------------------------------------------------------|
| `kappa1`  | $\kappa_1$            | Unknown value of true calibration input $g$                       |
| `thetaS1` | $\theta_{s1}$         | Scale parameter of $h$ input for simulator correlation            |
| `thetaS2` | $\theta_{s2}$         | Scale parameter of $g$ input for simulator correlation            |
| `alphaS1` | $\alpha_{s1}$         | Smoothness parameter of $h$ input for simulator correlation       |
| `alphaS2` | $\alpha_{s2}$         | Smoothness parameter of $g$ input for simulator correlation       |
| `thetaB1` | $\theta_{b1}$         | Scale parameter of $h$ input for bias-correction correlation      |
| `alphaB1` | $\alpha_{b1}$         | Smoothness parameter of $h$ input for bias-correction correlation |
| `sigma2S` | $\sigma^2_s$          | Marginal variance of simulator covariance                         |
| `sigma2B` | $\sigma^2_b$          | Marginal variance of bias-correction covariance                   |
| `sigma2E` | $\sigma^2_{\epsilon}$ | Variance of random measurement error in field                     |
| `muB`     | $\mu_b$               | bias-correction mean                                              |

There are ten more elements in `caqlMod` other than matrix `Phi`, which
are listed here along with a brief description. The data frame
`estimates`, which provides a summary table of all model parameters,
includes the mean, mode, median, standard deviation, and 50% and 80%
upper and lower quantiles of the marginal distribution of all
parameters. The vector `logPost` contains the posterior log likelihood
given a parameter draw from `Phi` matrix. The nested list `priors`
contains prior specifications for all parameters. The vector
`acceptance` represent the acceptance rate of each parameter after
running the MCMC algorithm with adaptive proposal. It is important to
note that, the current implementation of MCMC algorithm employs
Metropolis-Within-Gibbs variation, which is a one-dimensional proposal
scheme and the optimal acceptance rate must be close to 0.44. The vector
of strings `vars` contains the parameter notation used in code and as
`Phi` column headers. The rest of the components in `Phi` are not of
great importance for user, however, they are needed for calibrated
prediction. The list of matrices and vectors `data` includes the
training data in the forms that match KOH model components. The numeric
vector `scale` contains scaling factors that are used to scale the
training data during calibration. The named list `indices` contain the
indices of parameters in each row of `Phi` matrix. The function list
`priorFns` includes the prior functions that are created during
calibration based on given prior specifications. And finally,
`proposalSD` is a numeric vector that represents the final standard
deviation of proposal for each parameter.

### Calibrated Prediction

Similar to any other predict function, `predict()` requires a model
object argument called `object`, which in `FBC` package must be a `fbc`
object, along with an argument representing a new input configuration
called `newdata`. Moreover, current implementation of the `predict()`,
support two different methods of prediction: Maximum A Posteriori
(“MAP”) and MCMC-based fully Bayesian (“Bayesian”) methods. The method
can be selected using `method` argument that can take a character string
value of either “MAP” or “Bayesian”.

``` r
predsMAP   <- predict(object = calMod, newdata = matrix(c(2.2, 2.4), ncol = 1), method = "MAP")
predsBayes <- predict(object = calMod, newdata = matrix(c(2.2, 2.4), ncol = 1), method = "Bayesian")
```

The return value of `predict()` is a list consisting of two fields:
`pred`, which is a vector of the predicted response for every new input
configuration (rows of `newdata`), and `se`, which is a vector of the
predicted response’s standard errors.

``` r
predsMAP  
```

    ## $pred
    ## [1] 0.6968523 0.7319061
    ## 
    ## $se
    ## [1] 0.0801 0.0801

``` r
predsBayes
```

    ## $pred
    ## [1] 0.697 0.732
    ## 
    ## $se
    ## [1] 0.0746 0.0746

In “Bayesian” method, which is the default value of `method` argument,
MCMC draws of calibration model parameters are used to form a
distribution of each predicted value. In particular, each rows of matrix
`Phi` from the output of the calibration model, is used to predict the
response for every row of `newdata`. Therefore, a distribution of
response is created for each new input configuration. Then, the
predictive mean and variance of the resulting distribution are used to
compute the point predictions as well as standard errors for
predictions.

In “MAP” method, the row of `Phi` matrix that results in maximum log
posterior, is extracted and taken as model parameters to compute both
point estimates and standard errors. This method is much faster than the
“Bayesian” method as it only computes the prediction for each row of
observation once.

### Specifying Parameter Priors:

As mentioned hyperparameters of the calibration model can be set using
`set_hyperParameters()` function. All arguments of this function, which
collectively specify the priors for all hyperparameters, have reasonable
default values. Therefore `set_hyperParameters()` can be used without
arguments to set the hyperparameters. In fact, the default value of the
`hyper` argument in `calibrate()` is the function
`set_hyperParameters()` without any argument. Nevertheless, when there
is prior belief about structure of correlation structures (either
simulator GP or bias-correction GP), these beliefs can be applied to the
model in the form of prior specification using `set_hyperPriors()`. For
example in the following snippet, only correlation scale parameters of
the simulation GP are set to “beta” distributions and the second
parameter of beta distribution is set to 6. The first parameter of the
beta distribution, the initial value for these parameters, and all other
parameter specification remain unchanged.

``` r
priors <- set_hyperPriors(thetaSDist = "beta", thetaBP2 = 6)
```

### Summarizing the Calibration Model

Both `summary()` and `print()` generic functions are implemented to work
with the output of `calibrate()` function. In particular, given a `fbc`
object, `summary()` returns the `estimate` component of `calibrate()`
output, which is a data frame containing statistical summary of
calibration parameters. Similarly, `print()` will display the same
summary data frame in the console.

``` r
calModSum <- summary(calMod) 
print(calMod)
```

    ##       mean  median    mode   lwr50  upr50   lwr80   upr80     sd
    ## 1   0.3154  0.2803  0.2885  0.1710 0.4453  0.0930  0.6058 0.1990
    ## 2   1.3368  1.3245  1.1670  1.1384 1.5366  0.9934  1.7134 0.2925
    ## 3   0.3443  0.3019  0.2735  0.2560 0.4280  0.1930  0.4892 0.1304
    ## 4   1.9972  1.9978  1.9980  1.9971 1.9985  1.9959  1.9991 0.0030
    ## 5   1.9980  1.9995  1.9996  1.9980 2.0000  1.9957  2.0000 0.0051
    ## 6   0.8740  0.8477  0.8390  0.6629 1.0787  0.5160  1.2629 0.2801
    ## 7   1.8164  1.8625  1.9446  1.7202 1.9419  1.5982  1.9624 0.1464
    ## 8   7.2093  6.8783  4.8156  4.8918 9.1035  3.9589 11.8601 2.9346
    ## 9   0.1727  0.1704  0.1796  0.1262 0.2010  0.1048  0.2448 0.0540
    ## 10  0.1130  0.1127  0.1239  0.0997 0.1246  0.0965  0.1302 0.0133
    ## 11 -0.1727 -0.3197 -0.5427 -0.6788 0.3548 -0.8869  0.7469 0.5987

### Visualization of Calibration Model

The implementation of this generic function `plot()` in `FBC` package
enables visualization of a calibration model results. Given a
calibration model in the form of a `fbc` object using argument `x`,
`plot()` can visualize the model in three different mode, which can be
chosen by supplying the `type` argument.

In particular, `type` must be a character string form “density”,
“trace”, and “fits”. Using “density”, which is the default value of the
`type` argument, `plot()` will plot the marginal posterior density
distribution for the given parameter using `parameter` argument. It does
so by estimating a density function given the MCMC-based posterior draws
of the parameter. It also plots the prior distribution of the given
parameter in the same plot for ease of comparison and visualizes the
empirical mode of the posterior density. The `parameter` argument must
be a character string consistent with the notation used throughout the
package and current vignette , namely from “kappa”, “thetaS”, “alphaS”,
“thetaB”, “alphaB”, “sigma2S”, “sigma2B”,“sigma2E”, or “muB”. The
default value for `parameter` argument is “kappa” or calibration
parameters, which usually are the parameters of the interest. Note that
`parameter` argument characterizes the class of parameters and when
there is more than one parameter in that class,`plot()` will plot the
density distribution for all parameters in that class in separate plots.

``` r
# Note that there are two correlation scale parameters in simulator GP and there will be two plots
plot(calMod, parameter = "thetaS")
```

![](README_files/figure-gfm/plot%20density-1.png)<!-- -->![](README_files/figure-gfm/plot%20density-2.png)<!-- -->

Using “trace” as `type` argument, `plot()` will plot the progression of
the given parameter as MCMC draws are taken. It is similar to the time
series of the parameter but indexed with number of iteration in MCMC
rather than time. The trace plot can be used to determine whether there
is good mixing in MCMC draws. Using “trace” as `type` argument also
requires supplying the `parameter` from aforementioned list of possible
parameter classes. And similar to density plots, trace plots will be
plotted for all of the parameters in the given class in separate plots.

``` r
# Note that there are two correlation smoothness parameters in simulator GP and there will be 
# two plots
plot(calMod, parameter = "alphaS", type = "trace")
```

![](README_files/figure-gfm/plot%20trace-1.png)<!-- -->![](README_files/figure-gfm/plot%20trace-2.png)<!-- -->

And finally using “fits” as `type` argument, `plot()` will plot the
fitted values of the response versus all experimental variables in
separate plots. The `parameter` argument is not required for this type
and will be ignored. Fits plots can be used to visually determine the
goodness of fits versus actual response during interpolations.
Internally, `plot()` will use `predict()` function (using “MAP” method)
to compute the fitted values for the training input configurations and
will plot them in along with actual values. Furthermore, `xlab` argument
can be supplied with a character string to characterize the experimental
variables’ names. If not supplied, the x-axis will labelled by “x1”,
“x2”, until last experimental variable.

``` r
# Plots the fitted values versus all experimental inputs along with actual values in separate plots
plot(calMod, type = "fits", xlab = "height")
```

![](README_files/figure-gfm/plot%20fits-1.png)<!-- -->

There are also three more exported functions in `FBC` package that are
not required to build calibration models or to perform calibrated
prediction. These functions are used internally in the package but since
they offer functionalists that are not supported in base R, they are
exported for use.

### Computing Correlation

This function is used to compute the correlation between rows of two
given matrices assuming a power exponential correlation family
structure. This family is a generalization of correlation, where
distinct scale and smoothness hyperparameters are used for different
dimensions of the given data (columns of given matrices) when computing
the correlation. In `FBC`, `correlation()` function treats different
dimensions using separate scale and smoothness parameters. Note that
only the first matrix, characterized by argument `X` must be supplied
and the default value for the second matrix, characterized by `Y` is
`NULL`. When only a single matrix is supplied, `correlation()` will
compute the correlation of that matrix with itself. Other than the given
matrix or matrices, which must have same number of columns, user must
supply two vectors with `theta` and `alpha` arguments to characterize
scale and smoothness parameters. The length of both vectors either must
be same and equal to the number of columns in given matrices, or they
must be scaler in which case for all dimensions same values of scale or
smoothness will be used.

``` r
X     <- matrix(c(1, 3, 5,
                  2, 2, 6,
                  1, 4, 1), nrow = 3, byrow = TRUE)

Y     <- matrix(c(7, 3, 0,
                  2, 2, 4), nrow = 2, byrow = TRUE)

sc    <- c(1, 2, 3) # scale parameters of correlation structure
sm    <- c(2, 1, 2) # smoothness parameters of correlation structure

# correlation of a matrix with itself
round(correlation(X, theta = sc, alpha = sm), 5)
```

    ##         [,1]    [,2] [,3]
    ## [1,] 1.00000 0.00248    0
    ## [2,] 0.00248 1.00000    0
    ## [3,] 0.00000 0.00000    1

``` r
# correlation between two matrices
round(correlation(X, Y, theta = sc, alpha = sm), 5)
```

    ##      [,1]    [,2]
    ## [1,]    0 0.00248
    ## [2,]    0 0.00001
    ## [3,]    0 0.00000

### Building Prior Functions

This function will create a prior function based on given distribution
and distribution parameters. The function arguments are `prior`, which
characterize the distribution family of the prior and two other the
arguments, `p1` and `p2` that characterize the parameters of the chosen
distribution. The `prior` argument can take a character string value
from “uniform”, “normal”, “normalTr”, “lognormal”, “gamma”,
“inversegamma”, “beta”, “betashift”, “logbeta”,
“logistic”,“exponential”, “fixed”. For example for `prior = "gamma"`,
`p1` and `p2` determine shape and scale of the gamma distribution, or
for Gaussian distribution, `p1` and `p2` determine mean and standard
deviation of the Gaussian distribution. When “fixed” is used for
`prior`, the prior function is simply `x = 1`. The output of
`build_prior()` is a function that given a value, computes the
probability density of the chosen distribution. To reduce the load of
internal computation during creation of the calibration model,
`build_prior()` computes the log of density and if used externally must
be transformed.

``` r
# create a prior function for beta(2, 5). Note that the function compute log of priors and must be transformed
pr_fun <- build_prior(dist = "beta", p1 = 2, p2 = 5)$fun
round(exp(pr_fun(c(-1, 0, 0.1, 0.5, 0.9, 1, 2))), 3)
```

    ## [1] 0.000 0.000 1.968 0.938 0.003 0.000 0.000

``` r
# create a prior function for a Uniform distribution with lower bound of -10, and upper bound of 10
pr_fun <- build_prior(dist = "uniform", p1 = -10, p2 = 10)$fun
round(exp(pr_fun(c(-11, -5, 0, 4, 10, 12))), 3)
```

    ## [1] 0.00 0.05 0.05 0.05 0.05 0.00

``` r
# create a prior function for Gaussian distribution with mean of 1 and standard deviation of 2
pr_fun <- build_prior(dist = "normal", p1 = 1, p2 = 2)$fun
round(exp(pr_fun(c(-9, -5, -3, -1, 1, 3, 5, 7, 11))), 3)
```

    ## [1] 0.000 0.002 0.027 0.121 0.199 0.121 0.027 0.002 0.000

### Estimating the Mode of a Continious Variable

The function `pmode()` computes an estimate of the mode for a continuous
distribution. It works similar to a histogram, in which the domain of
the distribution is broken into same length bins. For each bin, the
draws of the distribution that fall within that bin is counted and at
the end the mean of the bin with highest count is returned as estimated
mode. The function also takes the number of bins through `breaks`
argument. The default value is `NULL`, which makes `pmode` to determine
the number of required bins dynamically based on number of data points
and the domain of the distribution.

``` r
# find the estimated mode of a vector
vec <- runif(100, 0, 10)
pmode(vec) 
```

    ## [1] 5.596472

``` r
pmode(vec, breaks = 10)
```

    ## [1] 7.317986

[^1]: *Field experiments, which are sometimes called physical
    experiments in other texts, will be represented by f subscript (for
    field) and computer simulations, which are sometimes called computer
    experiments in other texts, will be represented by s subscript (for
    simulation) throughout this vignette.*

[^2]: *Although there are relatively precise estimates of earth’s
    gravity, namely $9.81 m/s^2$, it is still an estimate and subject to
    some level of uncertainty.*

[^3]: *The raw data can be downloaded from
    [here](https://bobby.gramacy.com/surrogates/ball.csv)*

[^4]: *To produce the design matrix, `lhs` R package (Rob Carnell, 2022)
    is used. For more information on LHS method please click
    [here](https://en.wikipedia.org/wiki/Latin_hypercube_sampling#cite_note-C3M-1)
    and for information on `lhs` package please click
    [here](https://cran.r-project.org/web/packages/lhs/index.html).*

[^5]: *The preprocess code that produces `ballField` and `ballSim` from
    raw data can be found
    [here](https://github.com/parpishro/FBC/blob/main/data-raw/ball/preprocess.R).*

[^6]: *Running the `calibrate()` function with reasonable number of MCMC
    runs is too lengthy to run while knitting the current vignette.
    However, this command is run beforehand and the result is saved in
    package data to be loaded and used in the vignette. This `fbc`
    object, which is called `calMod`, can be found
    [here](https://github.com/parpishro/FBC/blob/main/inst/calMod.rda).*

[^7]: *“betashift” refers to a beta distribution that is shifted one
    unit to right to cover \[1 2\] interval and it is used to specify
    the priors for correlation smoothness parameters which must be
    constrained to \[1 2\] interval. Moreover, choosing “fixed” as
    distribution will exclude that class of parameters from the MCMC
    sampling. In this case, the given initial value will be used as
    fixed parameter value and p1 and p2 arguments will be used.*

[^8]: *For example, if “normal” distribution is used, p1 and p2
    represent mean and variance of the distribution and if “uniform”
    distribution is used, p1 and p2 represent lower and upper bound of
    the distribution.*

[^9]: *Not all distribution types require two arguments. In particular,
    “exponential” distribution only requires rate parameter (p1) and
    “jeffreys” requires none. In these cases the unused arguments are
    ignored.*

[^10]: *Although a vague prior is specified as default for calibration
    parameters values, the user is encouraged to specify the prior based
    on the field knowledge. The default vague prior is specified using a
    uniform distribution with lower bound of 0 and upper bound of 1,
    which characterize the lower and upper bound of parameter domain
    after standardization ($U(0, 1$).*

[^11]: *Running `calibrate()` with high number of parameters or MCMC
    runs will be a lengthy process. This argument shows the progress of
    algorithm on percentage basis, along with sample of results, which
    can be useful for debugging purposes.*
