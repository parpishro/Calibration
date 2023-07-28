---
title: "FBC: Full Bayesian Calibration"
author: "Parham Pishrobat and William J. Welch"
date: "2023-07-18"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{FBC}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




## Introduction

`FBC` is a package that uses the data from both physical experiment and computer simulation to 
calibrate the simulator model. Calibration is the process of fitting a statistical model to the 
observed data by adjusting calibration parameters. These parameters can represent various aspects of
the simulation model, such as tuning hyperparameter of the model itself or unknown but fixed 
physical properties that are governed by physical system. 

The `FBC` package is a based on the well-known Kennedy and O'Hagan (KOH) calibration model (2001). 
In its original formulation, KOH formulates a hierarchical Bayesian framework with two sequential 
phases. In the first phase, model hyperparameters are predicted using maximum likelihood estimation 
(MLE) method. In the second phase, Bayesian analysis is used to derive the posterior distribution of
calibration parameters while fixing the hyperparameters found in the first phase. Neither KOH's 
original formulation nor later suggestions are fully Bayesian as they use MLE methods to estimate 
hyperparameters, largely due to computational infeasibility (CITATION: Higdon, Berger, other). 

On the other hand, `FBC` runs a fully Bayesian model that includes both calibration parameters and 
model hyperparameters in its Bayesian framework. Implementation of the package optimizes the 
calibration process and memory management to increase computational efficiency. Moreover, the fully
implemented Bayesian framework, enables the user to input the information about all parameters and 
hyperparameters in the  form of prior specification. Common prior distribution are implemented in 
`FBC` to allows for high degree of flexibility in specifying the expert knowledge or the lack 
thereof. `FBC` runs a Markov Chain Monte Carlo (MCMC) algorithm to find the posterior distribution
of calibration parameters and model hyperparameters. Furthermore, `FBC` can be used for prediction 
of mean response and model discrepancy for a new test input configuration. 

(TODO: edit this paragraph to reflect the structure of the vignette in the end)
The current vignette is structured into following sections: The first section (`FBC` Usage) explains
the package functionality with a simple pedagogic example. Also in this section, the notation for 
inputs and outputs of both computer and physical experiments are introduced. The second section 
(Calibration Model) generalizes the example introduced in the first section to build model 
components. Using the general notation while referencing the example, we describe the KOH 
calibration model and its modelling choices. The third section (Parameter Estimation), presents the 
theoretical results that characterize the posterior distribution of model parameters and MCMC-based 
estimation of parameters. The fourth section (Prediction), provides the results to derive MCMC-based 
predictions for new input configuration using the estimated parameters and MCMC samples of model 
components. The fifth section (Implementation) explains some of the implementation features and 
choices that distinguishes `FBC` from other implementations while justifying those choices. Finally,
the last section (Application) provides two real-world and well-cited examples to demonstrate the 
full functionality and limitations of the `FBC` package. In the end, an appendix is provided to 
provide further details and a reference page to link for further readings.



## 1. `FBC` Usage

`FBC` package has two main public functions: `calibrate()` and `predict()`. As the name suggests, 
`calibrate()` takes the field and simulation data to calibrate the simulator model and `predict` 
takes a new input configuration and the calibration model (`fbc` object) to predict the field 
response at new input configuration and quantifies the uncertainty in prediction. In addition, `FBC` 
has a few helper function to aid in arguments entry and visualization.


### 1.1 Setup

Building a calibration model requires data from field experiment and computer simulation. To focus 
on the functionality of the package, we use a simple pedagogic example as experimental setting. 
In this setting, a wiffle ball is dropped from different heights and the time it takes to hit the 
ground is measured. This experiment has one experimental input height ($h$) and one calibration 
input gravity ($g$) to produce the response, time ($t$). Note that in the field experiment, the 
earth gravity is fixed but unknown to the experimenter (at least to some level of uncertainty) and
therefore it is not part of the input data. The field experiment with aforementioned specification 
has been performed by Derek Bingham and Jason Loeppky (CITATION). The data field is loaded in the 
package environment under `ballField` name. To increase robustness, `calibrate()` requires the 
training field data to include both field response vector $\bf{t}$ and input matrix $\bf{h}$ (a 
vector in ball example) in a single input matrix `ballField`, where first column always contain the 
response, followed by experimental input variables as columns in the field data matrix.



```r
library(FBC)
head(ballField, 3)
>         t     h
> [1,] 0.27 0.178
> [2,] 0.22 0.356
> [3,] 0.27 0.534
dim(ballField)
> [1] 63  2
```


Simulation of the ball drop experiment is easy to model using introductory physics results: 


$$
t = \sqrt{\frac{2h}{g}}
$$


We have implemented the above mathematical model as a code that takes $h$ and $g$ as experimental 
and calibration inputs and returns $t$ as simulator response. Latin Hypercube Design (CITATION) is
used to create the input design matrix consisting of two columns: $h$ and $g$. Similar to the field 
data, `calibrate` requires the simulation data to be packaged into a single matrix `ballSim`, where
first column always represents response vector ($\bf{t}$), followed by experimental and calibration
inputs as columns $[\bf{h} \quad \bf{g}]$.  



```r
head(ballSim, 3)
>          t     h      g
> [1,] 0.404 0.998 12.220
> [2,] 0.487 0.940  7.909
> [3,] 0.450 0.886  8.736
dim(ballSim)
> [1] 100   3
```


Figure 1 shows the distribution of time versus height for both physical and simulation experiments. Note that for higher height values, the simulation responses (blue) underestimate their corresponding field response (red), displaying a systemic bias for simulation model. 



![plot of chunk ball plot](figure/ball plot-1.png)



The ball example is a suitable pedagogical example as it has only one experimental input $\bf{h}$) 
and one calibration input $\bf{g}$ to model response $\bf{t}$. Using this toy example, we 
demonstrate the functionality of the package without getting to the details of a complex 
mathematical model. More complex and real-world examples are covered in the last section. 


### 1.2 \ \ `calibrate()` \ \ 


**Arguments:** The `calibrate()` function takes three sets of arguments from user: data, MCMC, and prior 
specification parameters. Other than data arguments which must be supplied by user, all other 
arguments have reasonable default values for ball example.


```r
output <- calibrate(sim = ballSim, field = ballField,                                 # Data 
                    nMCMC = 11000, nBurn = 1000, thinning = 50,                       # MCMC
                    kappaDist = "beta", kappaInit = NA, kappaP1 = 1.1, kappaP2 = 1.1, # Priors
                    hypers = set_hyperPriors(),
                    showProgress = FALSE) 
```


The first and second arguments `sim` and `field` must be supplied by user. Both must be either 
in matrix or dataframe format representing the simulation and field data respectively. In both 
`ballSim` and `ballField` of our ball example, the first column represent the response $\bf{t}$ and 
the second column represent experimental input $\bf{h}$. Additionally, `ballSim` has a third column 
that represents calibration input $\bf{g}$. Calibration inputs are implicit in field experiment and 
are absent in the data matrix `ballField`. 


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

The second set of arguments consists of MCMC parameters such as number of total iterations `nMCMC`, 
number of burn-in iterations to be removed from the beginning of the chain `nBurn`, and `thinning`, 
which indicate the sampling rate to remove the autocorrelation from the sampled draws. For example,
when `thinning = 50`, for every 50 draws from the result only one will be kept in order to remove 
the correlation between draws. 


The third and last set of arguments consists of prior specification for each parameter of the model. 
The main goal of calibration is to estimate the calibration parameters ($g$ in ball example). 
However, employing KOH calibration model introduces seven additional classes of hyperparameters to 
the model. All eight classes of parameters are unknown and must be estimated. As `FBC` employs a 
full Bayesian framework, the priors for all model parameters must be specified in advance. 
Throughout the package implementation and current guide a consistent notation is used to denote 
parameters: $\bf{\kappa}$ (calibration parameters), $\bf{\theta_s}$ (simulator correlation scale 
parameters), $\bf{\alpha_s}$ (simulator correlation smoothness parameters), $\bf{\theta_b}$ 
(bias-correction correlation scale parameters), $\bf{\alpha_b}$ (bias-correction correlation 
smoothness parameters), $\sigma^2_s$ (marginal simulator variance), $\sigma^2_b$ (marginal 
bias-correction variance), $\sigma^2_{\epsilon}$ (measurement error variance), and $\mu_b$ (mean of
bias-correction Gaussian Process. All calibration model components along with the introduced 
parameters are explained in more detail in next section. 

For each class of parameters, there are four associated arguments. 1) distribution (suffixed with 
"dist" in argument names) type which is a character string determining the prior distribution. 
Currently, `FBC` supports almost all common distributions and can be chosen from ("uniform", 
"gaussian", "gamma", "beta", "lognormal", "logistic", "betashift", "exponential", "inversegamma", "jeffreys", "fixed"). "betashift" refers to a beta distribution that is shifted one unit to right 
to cover [1 2] interval and it is used to specify the priors for correlation smoothness parameters 
which must be constrained to [1 2] interval. Moreover, choosing "fixed" as distribution will exclude 
that class of parameters from the MCMC sampling. In this case, the given initial value will be used 
as fixed parameter value and p1 and p2 arguments are not used. 2) Initial value (suffixed with 
"init" in argument names) which represent starting point of the parameter(s) in MCMC algorithm. 3) 
and 4) represent the two parameters (p1 and p2) of the chosen distribution. For example, if 
"gaussian" distribution is used, p1 and p2 represent mean and variance of the distribution and if 
"uniform" distribution is used, p1 and p2 represent lower and upper bound of the distribution.
Please note that not all distribution types require two arguments. In particular, 
"exponential" distribution only requires tone parameter (p1) and "jeffreys" requires none. In these 
cases the unused arguments are ignored.. 


From the classes of parameters mentioned, the prior for calibration parameters must be specified by 
user based on field knowledge as these parameters are problem-specific, even though a general vague 
prior is specified as default values. In contrast, all other parameters have reasonable default 
values based on literature and can be left unchanged. For this reason and to avoid unwanted 
complexity, all other priors are specified using a helper function set_hyperPriors()`. Using this
function without argument (which is the default value of the `hyper` argument) will specify all 
parameters based on their default values. However, expert knowledge about one or more of these 
parameters can be supplied using arguments of the `set_hyperPriors()` function to change the default
prior specifications. 


```r
priors <- set_hyperPriors(thetaSDist = "beta", thetaBP2 = 6)
```


Thge following table summarizes all prior arguments to build the calibration model. 

\rule{\linewidth}{0.75pt}
**Table 1:** *Argument names to specify priors for each parameter class.*

|Parameter Class                                   |Distribution |Initial Value|Shape Parameters         |
|:-------------                                    |:------      |:-----       |:-----------             |
|True field calibration inputs ($\kappa$)          |`kappaDist`  |`kappaInit`  |`kappaP1`  &  `kappaP2`  |
|                                                  |             |             |                         |
|Simulation GP scale ($\theta_s$)                  |`thetaSDist` |`thetaSInit` |`thetaSP1` & `thetaSP2`  |
|                                                  |             |             |                         |
|Simulation GP smoothness ($\alpha_s$)             |`alphaSDist` |`alphaSInit` |`alphaSP1` & `alphaSP2`  |
|                                                  |             |             |                         | 
|Bias-correction GP scale ($\theta_b$)             |`thetaBDist` |`thetaBInit` |`thetaBP1` & `thetaBP2`  | 
|                                                  |             |             |                         |
|Bias-correction GP smoothness ($\alpha_b$)        |`alphaBDist` |`alphaBInit` |`alphaBP1` &`alphaBP2`   | 
|                                                  |             |             |                         |
|Simulation marginal variance ($\sigma^2_s$)       |`sigma2SDist`|`sigma2SInit`|`sigma2SP1` & `sigma2SP2`|
|                                                  |             |             |                         | 
|Bias-correction marginal variance ($\sigma^2_b$)  |`sigma2BDist`|`sigma2BInit`|`sigma2BP1` & `sigma2BP2`|
|                                                  |             |             |                         | 
|Measurement error variance ($\sigma^2_{\epsilon}$)|`sigma2EDist`|`sigma2EInit`|`sigma2EP1` & `sigma2EP2`| 
|                                                  |             |             |                         | 
|Bias-correction mean ($\mu_b$)                    |`muBDist`    |`muBInit`    |`muBP1`    & `muBP2`     |


Distribution arguments must be chosen from the string list of implemented prior distributions and
initial value and the two shape parameter arguments must be numeric. Note that some classes of 
parameters may contain more than one parameters. In this case, the argument values can be vector
instead of default scaler. In this case all four fields of that parameter class must be in vector
format with same length as number of parameters in the same class. For example if there are five calibration parameters, `kappaDist` can either be a scaler, in which case the distribution for all calibration parameters will be set to that scaler value (and `kappaInit`, `kappaP1`, and `kappaP2`
must be also scaler), or a vector of length five that supplies the distribution types for all 
calibration parameters (and `kappaInit`, `kappaP1`, and `kappaP2` must also be vector of length 5).

Moreover, there is a logical argument `showProgress` that indicate whether function must show the
progress in calibration on console. This will put the `calibrate()` in interactive mode and will
show the percentage of the MCMC draws along with sample draws. 

**Output:** The output of the `calibrate()` function is a `fbc` object that contains the samples 
from posterior joint distribution of parameters, along with other model information.






```r
names(output)
>  [1] "Phi"        "estimates"  "logPost"    "priors"     "acceptance" "vars"       "data"      
>  [8] "scale"      "indices"    "priorFns"   "proposalSD"
```


The first and main component of the output is matrix `Phi` whose columns represent the sample of 
posterior densities for each unknown parameter of the model in the same order as parameter classes in
table 1. Since $\kappa$, $\theta_s$, $\alpha_s$, $\theta_b$, $\alpha_b$ parameter classes may 
contain more than one parameter, they are suffixed by a number that represents the index of the 
parameter in the class. For example, if there are 3 calibration inputs, the first $p$ columns of 
the matrix `Phi` represent posterior density of calibration parameters and the column headers will 
be `kappa1`, `kappa2`, and `kappa3`. In the ball example, there is only one calibration parameter 
$g$, which is denoted by $\kappa_1$ and represented by `kappa1` in matrix `Phi`. Each row of the 
matrix `Phi` represents a MCMC draw.



```r
head(output$Phi, 3)
>   kappa1 thetaS1 thetaS2 alphaS1 alphaS2 thetaB1 alphaB1 sigma2S sigma2B sigma2E   muB
> 1   6.16    5.14    0.69    1.94    1.84    0.68    1.90    0.09    0.30    0.09 -0.35
> 2   8.72    4.92    1.04    1.94    1.97    0.17    1.81    0.09    0.66    0.09  0.24
> 3  13.70    4.35    0.76    1.90    1.93    0.43    1.79    0.10    0.44    0.08  0.15
```


All other columns of matrix `Phi` represent the samples for other model hyperparameters. Table 2
provides an overview of the model parameters and the notation to represent them in ball example. 
Later sections will explain in detail why are these hyperparameters introduced, what do they 
represent, and how they are estimated.


**Table 2:** Notation used in matrix `Phi` to represent parameters in ball example.

|Column   |Notation             |Description                                                       |
|:---     |:---                 |:-------------------                                              |
|`kappa1` |$\kappa_1$           |Unknown value of true calibration input $g$                       |
|         |                     |                                                                  |
|`thetaS1`|$\theta_{s1}$        |Scale parameter of $h$ input for simulator correlation            |
|         |                     |                                                                  |
|`thetaS2`|$\theta_{s2}$        |Scale parameter of $g$ input for simulator correlation            |
|         |                     |                                                                  |
|`alphaS1`|$\alpha_{s1}$        |Smoothness parameter of $h$ input for simulator correlation       |
|         |                     |                                                                  |
|`alphaS2`|$\alpha_{s2}$        |Smoothness parameter of $g$ input for simulator correlation       |
|         |                     |                                                                  |
|`thetaB1`|$\theta_{b1}$        |Scale parameter of $h$ input for bias-correction correlation      |
|         |                     |                                                                  |
|`alphaB1`|$\alpha_{b1}$        |Smoothness parameter of $h$ input for bias-correction correlation |
|         |                     |                                                                  |
|`sigma2S`|$\sigma^2_s$         |Marginal variance of simulator covariance                         |
|         |                     |                                                                  |
|`sigma2B`|$\sigma^2_b$         |Marginal variance of bias-correction covariance                   |
|         |                     |                                                                  |
|`sigma2E`|$\sigma^2_{\epsilon}$|Variance of random measurement error in field                     |
|         |                     |                                                                  |
|`muB`    |$\mu_b$              |bias-correction mean                                              |



`estimates` is a data frame that provides a summary table of all model parameters. `summary()` 
function also provides the same statistics given a `fbc` object. `logPost` is vector of computed
posterior log likelihood given a parameter draw from `Phi` matrix. Therefore, the length of the
`logPost` vector is same as number of rows in `Phi`. `priors` is a nested list that contains prior specifications for all parameters (including calibration parameters). `acceptance` is a vector of same length as number of unknown parameters and represent the acceptance rate of each parameter in MCMC
algorithm. The current implementation of MCMC algorithm employs Metropolis-Within-Gibbs variation, 
which is a one-dimensional proposal scheme and the optimal acceptance rate must be close to 0.44. 
`vars` is a vector of character string with the same length as number of parameters and contains the parameter notation used in code and as column headers. `data` is a list consisting of a training 
data in the form of vectors and matrices. `scale` is a numeric vector that contains scaling factors 
used during calibration to scale the training data. `indices` is a list of vectors containing the 
indices of parameters in each row of `Phi` matrix. `priorFns` is a list of prior functions that are
created during calibration based on given prior specifications. And finally, `proposalSD` is a 
numeric vector that represents the final standard deviation of proposal for each parameter.



### 1.4 `predict()`

TODO: 1. mean simulator response, 2. mean field response, 3. mean bias

### 1.5 Helpers

#### `control()`:



#### `summary()`:



#### `plot()`:

TODO: 1. format data input arguments, 2. summary, 3. plots?


\pagebreak


## 2. Calibration Model

Calibration model is statistical model that represent both field and simulator response as function of input configuration. In this section, we explain the theory behind building a calibration model and relate the notation used in the formulation of calibration model and in the package implementation to our ball example.


### 2.1 Data

A computer experiment or simulation is simply running a computer code at different input configurations and recording the response. The code is an implementation of the mathematical model that is intended to mimic the physical experiment. In general, a simulation has $p$ experimental inputs but also has $q$ additional calibration inputs that are either tuning parameters or unknown physical properties that are not controllable by field experimenter. Table 2 describes the response and inputs of a computer experiment with $m$ observations using matrix notation. The subscript $s$ is used to denote simulation.


\rule{\linewidth}{0.75pt}
**Table 2:** Notation used to represent simulation data component of calibration model.

| ***Notation*** | ***Description***                                 | ***Ball Example***                            |
| ---            | ---------------                                   | ------                                        |
| $m$            | Number of simulation runs                         | 100                                           |
|                |                                                   |                                               |
| $p$            | Number of experimental inputs                     | 1                                             |
|                |                                                   |                                               |
| $q$            | Number of calibration inputs                      | 1                                             |
|                |                                                   |                                               |
| $\bf{x}_s$     | Simulation input vector containing $(p+q)$ inputs | ($h$, $g$) (vector of length 2)               |
|                |                                                   |                                               |
| $y_s$          | Univariate simulation response                    | $t$ (scaler)                                  |
|                |                                                   |                                               |
| $\bf{X}_s$     | $(m \times (p+q))$ simulation input matrix        | [$\bf{h}$ $\bf{g}$] ($(100 \times 2)$ matrix) |
|                |                                                   |                                               |
| $\bf{y}_s$     | Vector of $m$ univariate simulation response      | $\bf{t}$ (vector of length $100$)             |



On the other hand, a physical experiment consists of $n$ observations of a physical property, each with $p$ experimental inputs. The calibration inputs are implicit in physical experiment and their values are unknown but assumed to be fixed throughout observations. To represent both experiments in a unified structure, the unknown calibration inputs of field experiment $\kappa$ must be augmented to experimental inputs, so that both physical and computer experiments have $(p+q)$ inputs ($h$ and $g$ in ball example) and a univariate response ($t$ in ball example). Assuming vector $\bf{x}_f$ represent experimental input, the augmented input vector is denoted by $\bf{x}_{\kappa}$. Since $\bf{x}_f$ has $p$ elements and $\bf{\kappa}$ has $q$ elements, both will have $(p+q)$ elements similar to $\bf{x}_s$. Stacking the input vectors, we can represent all field and augmented input vectors using $\bf{X_f}$ and $\bf{X_{\kappa}}$. Similarly, vector $\bf{y_f}$ represents all field observations. Subscripts $f$ and $\kappa$are used to denote field and augmented data respectively. 

Elements of vector $\kappa$ are parameters of the calibration model and will be estimated by `calibrate()` [^1]. 

\

[^1]: *The first $q$ columns of matrix `Phi` in `calibrate` output represent the MCMC samples of posterior densities for calibration parameters. In our ball example, there is only one calibration parameter (gravity), which is denoted by $\kappa_1$ and represented by `kappa1` in `Phi`.*</em>.

\pagebreak

\rule{\linewidth}{0.75pt}
**Table 3:** Notation used to represent field data component of calibration model.

| ***Notation***    | ***Description***                                                | ***Ball Example***                         |
| ---               | ---------------                                                  | -------                                    |
| $n$               | Number of field observations                                     | 63                                         |
|                   |                                                                  |                                            |
| $p$               | Number of experimental inputs                                    | 1                                          |
|                   |                                                                  |                                            |
| $\bf{x}_f$        | Field input vector containing $p$ experimental inputs            | $h$ (scaler)[^2]                          |
|                   |                                                                  |                                            |
| $\kappa$          | Vector of unknown true calibration inputs in field experiment    | true value of gravity $\kappa_1$           |
|                   |                                                                  |                                            |
| $\bf{x}_{\kappa}$ | Augmented field input vector containing $(p+q)$ inputs[^3]       | ($h$, $\kappa_1$)$^2$ (vector of length 2) |
|                   |                                                                  |                                            |
| $y_f$             | Univariate field response                                        | $t$ (scaler)                               |
|                   |                                                                  |                                            |
| $\bf{X}_f$        | $(n \times p)$ field input matrix                                | $\bf{h}$ (vector$^2$ of length $63$)       |
|                   |                                                                  |                                            |
| $\bf{X}_{\kappa}$ | $(n \times (p+q))$ augmented field input matrix                  | $[\bf{h} \quad \bf{\kappa_1}]$[^4]              |
|                   |                                                                  |                                            |
| $\bf{y}_f$        | Vector of $n$ univariate field response                          | $\bf{t}$ (vector of length $63$)           |
|                   |                                                                  |                                            |

[^2]: *In ball example there is only one experimental input and therefore $\bf{x}_f$ is a vector of length one or scaler. Similarly, in its matrix notation, $\bf{X}_f$ is a $(n \times 1)$ matrix or a vector of length $n$.*</em>.
[^3]: *Note that calibration parameters $\bf{\kappa}$ are often assumed to be unchanged throughout field experiment. Therefore, same ($\kappa$) vector is augmented to all of the field input configurations.*</em>.
[^4]: *Since calibration input is fixed for all field observation, $\bf{\kappa_1} = (\kappa_1, ... , \kappa_1)$ vector of length 63} ] ($(100 \times 2)$ matrix).*</em>.



All of the data components in Table 2 and 3 are generated internally to build the calibration model. The user is only required to provide two matrices representing the field and simulation datasets in their entirety through `sim` and `field` arguments of `calibrate()`. 


### 2.2 KOH Model 

#### Random Functions:
KOH models the functional relationship between simulation input and output as a realization of a random function $\eta(\bf{x_s})$. Similarly, KOH models the functional relationship between field input and output as a realization of random function $\eta(\bf{x}_{\kappa})$ but acknowledges a systemic model discrepancy and measurement errors. As a result, KOH models the discrepancy by adding a bias-correction term as realization of another random function $\delta_{\kappa}(\bf{x}_f)$. The error term $\epsilon$ is considered to be an independent draw from a normal distribution with zero mean and unknown variance $\sigma^2_{\epsilon}$. 


$$
\begin{aligned}
\epsilon  \ \ &\sim  \mathcal{N}(0, \sigma^2_{\epsilon}) \\ \\
y_f           &=     \eta(\bf{x_\kappa}) + \delta_{\kappa}(\bf{x_f}) + \epsilon \\ \\
y_s           &=     \eta (\bf{x_s})  \\
\end{aligned}
$$

Therefore other than $\bf{\kappa}$ and $\sigma^2_{\epsilon}$ parameters, the random functions $\eta(.)$ and $\delta_{\kappa}(.)$ are also unknown and must be specified. KOH models $\eta(.)$ and $\delta_{\kappa}(.)$  by two independent Gaussian Processes (GP). 


$$
\begin{aligned}
\eta(.) \  &\sim GP \ (0, \ \sigma^2_s . R_s(., .)) \\
\delta_{\kappa}(.) &\sim GP \ (0, \ \sigma^2_b . R_b(., .)) \\ \\
\end{aligned}
$$


Where $\sigma^2_s$ and $\sigma^2_b$ are marginal variance of simulator and bias-correction GPs, and $R_s(., .))$ and $R_b(., .))$ are correlation matrix of simulator GP (using full input matrix $\bf{X_s}$ or $\bf{X_{\kappa}}$) and bias-correction GP (using field input matrix $X_f$).

Note that means of GPs are considered to be zero because `calibrate()` function first standardizes simulation response $\bf{y}_s$ (mean zero and standard deviation of one) and then scales field response according to $\bf{y}_s$'s scaling factors. Furthermore, the simulator inputs (both experimental and calibration) are scaled to span $[0, 1]$ and the scaling factors of simulation experimental inputs are used to scale field experimental inputs. As a result considering zero mean for both processes seem reasonable. (TODO: look into the effect of constant mean for discrepancy GP).

#### Correlation Structure:
`FBC` employs a power exponential correlation family to represent the correlation structure of both GPs. Assuming $\bf{x}$ and $\bf{x'}$ are two rows of full input matrix (either $\bf{X_s}$ or $\bf{X_{\kappa}}$) and $\bf{x_f}$ and $\bf{x_f'}$ are two rows of field experimental input matrix ($\bf{X_f}$), the correlation matrices $R_s(\bf{x}, \bf{x'})$ and $R_b(\bf{x_f}, \bf{x_f'})$ are defined as following:


$$
\begin{aligned}
R_s(\bf{x}, \ \bf{x}') \ &= \prod^{p+q}_{i=1} e^{-\theta_i |x_i - x_i'|^{\alpha_i}} \\ \\
R_b(\bf{x_f}, \ \bf{x_f}')       &= \prod^{p}_{j=1} e^{-\theta_j |x_j - x_j'|^{\alpha_j}} \\
\end{aligned}
$$


Using separable power exponential correlation family introduces two new hyperparameters for each input: scale ($\theta_i$) and smoothness ($\alpha_i$). Together, they flexibly determine the shape of correlation structure. Table 4 introduce the notation used for hyperparameters of $\eta(.)$ and $\delta_{\kappa}(.)$.


\rule{\linewidth}{0.75pt}
**Table 4:** *New Parameters To Build KOH Calibration Model.*

|**GP**              | **Hyperparameters**                                                          | **Ball Example (Columns of `Phi`)**                  |
| ---                | ----------                                                                   | --------------                                      |
|                    |                                                                              |                                                     |
| $\eta(.)$          |$(\theta_{s1},...,\theta_{s(p+q)},\alpha_{s1},...,\alpha_{s(p+q)},\sigma^2_s)$|(`thetaS1`, `thetaS2`,`alphaS1`,`alphaS2`, `sigma2S`)|
|                    |                                                                              |                                                     |
|$\delta_{\kappa}(.)$| $(\theta_{b1},...,\theta_{bp},\alpha_{b1},...,\alpha_{bp},\sigma^2_b)$       | (`thetaB1`, `alphaB1`, `sigma2B`)                   |


#### Full Model:
After augmentation of true calibration inputs (vector $\bf{\kappa}$ to field data, both simulation and field experiment have the same input. KOH combines both components to build a joint model. The joint vector of all parameters in the final calibration model is denoted by:


$$
\phi = (\kappa_1, \ ... , \ \kappa_q, \ 
\theta_{s1}, \ ... , \ \theta_{s(p+q)}, \
\alpha_{s1}, \ ... , \ \alpha_{s(p+q)}, \
\theta_{b1}, \ ... , \ \theta_{bp}, \
\alpha_{b1}, \ ... , \ \alpha_{bp}, \ 
\sigma^2_s, \ \sigma^2_b, \ \sigma^2_{\epsilon})
$$


Note that in the ball example, the column headers of matrix `Phi` exactly match to model parameters.



```r
#head(output$Phi, 3)
```


### 2.3 Model Parameters

Table 5 provides a general overview of all model parameters, the notations, and corresponding parameters in the ball example.


**Table 5:** General notation used to represent model parameters and an example of corresponding identifiers in `Phi` matrix.

| **Parameter**                     | **General Notation**                               | **Ball Example**      | 
| --------------                    | --------                                           | -------------         |
| True field calibration inputs     | $\kappa = (\kappa_1, ... , \kappa_q)$              | (`kappa1`)            | 
|                                   |                                                    |                       | 
| Simulation GP scale               | $\theta_s = (\theta_{s,1}, ... , \theta_{s,p+q})$  | (`thetaS1` `thetaS2`) |  
|                                   |                                                    |                       | 
| Simulation GP smoothness          | $\alpha_s = (\alpha_{s,1}, ... , \alpha_{s,p+q})$  | (`alphaS1` `alphaS2`) |  
|                                   |                                                    |                       | 
| Bias-correction GP scale          | $\theta_b = (\theta_{b,1}, ... , \theta_{b,p})$    | (`thetaB1`)           |  
|                                   |                                                    |                       | 
| Bias-correction GP smoothness     | $\alpha_b = (\alpha_{b,1}, ... , \alpha_{b,p})$    | (`alphaB1`)           | 
|                                   |                                                    |                       | 
| Simulation marginal variance      | $\sigma^2_s$                                       | `sigma2S`             |
|                                   |                                                    |                       | 
| Bias-correction marginal variance | $\sigma^2_b$                                       | `sigma2B`             |
|                                   |                                                    |                       | 
| Measurement error variance        | $\sigma^2_{\epsilon}$                              | `sigma2E`             |  



Each row of matrix `Phi` represents a draw from joint distribution of parameters (MCMC run) and each column represents a parameter in the model. User-given initial values for parameters are used to initialize the first row of the `Phi` matrix. Then, each row will be used to find another sample from joint parameter space to fill the next row of `Phi` until matrix `Phi` is complete.


## 3. Parameter Estimation

`FBC` employs a full Bayesian approach to jointly estimate all parameters. To find the marginal posterior density distribution for each parameter, we need prior specification for each parameter (prior knowledge) and the joint likelihood estimation (full data).


### 3.1 Bayesian Analysis


Because `FBC` uses a full Bayesian framework, expert knowledge or opinion can be applied to the model parameters as prior specification. Variety of common prior distributions are implemented in `FBC` and can be used to specify the priors for each parameter (see section 1.2). There are seven classes of parameters and all have been specified in `calibrate()` using default values. Of those
seven classes, calibration parameters (vector $\bf{\kappa}$) and perhaps measurement error variance (scaler $\sigma^2_{\epsilon}$) are application-dependent. It is recommended for user to specify the prior arguments for these parameters based on prior knowledge or consensus. Nevertheless, prior for calibration parameters is defaulted to Beta(1.1, 1.1) distribution [^5]. It is close to standard uniform distribution (U(0, 1)) but densities approach to zero sharply as samples approach boundaries. This default choice has been made to ensure a somewhat non-informative prior while de-emphasizing on boundary values[^6]. For all other classes of parameters reasonable priors have been specified using default values. Priors for correlation scale parameters (vectors $\theta_s$ and $\theta_b$) have been set to Gamma(1.1, 0.1) distribution. Similarly, the priors for correlation smoothness parameters (vectors $\alpha_s$ and $\alpha_b$) have been set to Beta(5, 2) distribution that is shifted one unit to right to span [1, 2] as is the acceptable range for moothness parameters. This choice emphasizes higher (closer to 2 than 1) smoothness parameters. Finally, the priors for marginal simulator and  bias-correction and measurement error variances (scalers $\sigma^2_s$, $\sigma^2_b$, and $\sigma^2_{\epsilon}$) are set to be Inverse Gamma(1.5, 1.5). This emphasizes very low variances and de-emphasizes higher values.

[^5]: *Note that the range of Beta distribution is the span of [0, 1], which is the range of calibration inputs for simulator after scaling* </em>.
[^6]: *Boundary values of [0, 1] corresponds to $-\infty$ and $\infty$ in the original scale of calibration parameters.* </em>.



By representing both data in a joint calibration model, we can compute the conditional likelihood of response given a parameter vector
(See Appendix). Therefore, given the prior specifications above, we can derive joint posterior distribution of parameters given data:


$$
\mathcal{P} [\phi | \mathcal{D}]      \ \propto \ L(\mathcal{D} | \phi) \ .  \ \mathcal{P} [\phi]     
$$

Where $\mathcal{D}$ represent full data (field and simulation). However, the above formulation is intractable and thus we need a simulation-based method to sample from joint posterior distribution. `FBC` implements a version of Markov Chain Monte Carlo (MCMC) simulation.

### 3.2 MCMC Simulation

MCMC simulation algorithm is used to draw samples from joint posterior distribution of the parameter space and build `Phi` matrix 
row by row. MCMC algorithm creates a Markov chain by updating parameters in each iteration according to a proposal scheme. Then given this parameter vector ($\phi^{(i)}$) and data ($\mathcal{D}$), the posterior likelihood can be computed. If posterior probability density of the parameters is larger than the density of a random draw from standard uniform distribution, the algorithm keeps the that parameter configuration by writing the next row of matrix `Phi`, otherwise updates the new row by last parameter vector. In either case, the new row will be used to generate the next proposal. Note that the first row of `Phi`, which is needed to start the algorithm, is supplied by user through initial values for the parameters. The detailed algorithm is presented in the Appendix.



### 3.3 Parameter Posterior Distributions

At the end of the MCMC run, sample of the joint posterior distribution for each parameter (a column in matrix `Phi`) can be used as an approximation of marginal posterior distribution. Center measures such as mean, mode, or median are provided as parameter estimate depending on the application and distribution shape. Furthermore, 50% and 80% credible sets are formed for each parameters to quantify the uncertainty in estimation. These statistics are provided in `summary` element in the output of `calibrate()` and additionally. (TODO: fix estimate->summary)



```r
#output$estimates
```

Alternatively, posterior density kernels can be visualized over their assumed prior to investigate the effect of data on priors for each parameter. 



```r
# plot(output)
```


\pagebreak

## 4. Prediction

TODO


\pagebreak

## 5. Implementation

TODO

\pagebreak

## 6. Application

### Ball Example

##### Data: 


\pagebreak

### Spot Weld Example

##### Description:

##### Experimental Input:

The physical model has three inputs: gauge ($G$), load ($L$), and current ($C$):

- Gauge ($G$):
- Load ($L$): 
- Current ($C$):


##### Calibration Input:

The simulation model has one additional input, $\tau$ that affects the amount of
heat produced in the metal sheets. $\tau$ cannot be controlled in the physical
experiment and its value is unknown. However it has to be specified for the 
simulation model as calibration input $t$:

- Heat generation factor $\tau$: Factor affecting amount of heat produced

##### Mapping Parameters:

To map the spot weld data (both field and simulation data) and parameters to 
`FBC` input configuration, we use the dagger $\dagger$ superscript to distinguish process
parameters and variables with `FBC` variables and parameters:

$$
\begin{aligned}
\text{x}_1 \quad &\longrightarrow \quad G^{\dagger} \\
\text{x}_2 \quad &\longrightarrow \quad L^{\dagger} \\
\text{x}_3 \quad &\longrightarrow \quad C^{\dagger} \\
\kappa_1 \quad &\longrightarrow \quad \tau^{\dagger}
\end{aligned}
$$

\pagebreak

### Kinetic Example

\pagebreak

## Appendix

### MCMC Algorithm

\rule{\linewidth}{1pt}
$$\textbf{Implementation of Metropolis within Gibbs algorithm}$$
\rule{\linewidth}{1pt}

Let $\Phi$ be the matrix of parameter values (columns) indexed by MCMC 
iterations. Each column represents (after MCMC completes) the posterior density
of a parameters. Since all parameters are included in $\Phi$ but have 
overlapping indices, the parameter densities (columns) are renamed to 
$(\phi_1, ... , \phi_d)$, where $d = 4p + 3q + 3$ is total number of parameters
to have a unique index:   


$$
\begin{aligned}
& \ \ \mathbf{(1)}  \quad \quad \quad  \quad  \quad \quad \quad ... \quad \quad  \quad \quad \quad \quad \quad \quad \quad  \quad \quad \quad  \quad  \quad \ ... \quad  \quad \quad \quad \quad \quad \quad  \quad \quad \quad \quad \quad \quad  \quad \quad \quad \ \ \  \mathbf{(d)} \\ \\
& \quad  \textbf{k}_1^* \quad  ...  \ \textbf{k}_q^*  \quad \quad \textbf{t}_{s1}^* \ ... \ \textbf{t}_{s(p+q)}^* \quad \textbf{a}_{s1}^* \ \ ... \ \ \textbf{a}_{s(p+q)}^* \quad \ \ \textbf{t}_{b1}^* \ \ ...  \ \textbf{t}_{bp}^* \quad \quad \ \ \ \textbf{a}_{b1}^* \ \ ... \ \textbf{a}_{bp}^*  \quad \quad \textbf{v}_s^* \quad \ \ \ \textbf{v}_b^* \quad \ \ \ \textbf{v}_e^* \\ \\
\Phi = 
&\begin{bmatrix} 
k_1^{(1)} \ ... \ k_q^{(1)} & t_{s1}^{(1)} \ ... \ t_{s(p+q)}^{(1)} & a_{s1}^{(1)} \ ... \ a_{s(p+q)}^{(1)} & t_{b1}^{(1)} \ ... \ t_{bp}^{(1)} & a_{b1}^{(1)} \ ... \ a_{bp}^{(1)}  & v_s^{(1)} & v_b^{(1)} & v_e^{(1)} \\ \\
... \  \    & \  \  & \ ... \  &  \  \  &  \ ... \ &  &  & ... \\ \\
k_1^{(N)} \ ... \ k_q^{(N)} & t_{s1}^{(N)} \ ... \ t_{s(p+q)}^{(N)} & a_{s1}^{(N)} \ ... \ a_{s(p+q)}^{(N)} & t_{b1}^{(N)} \ ... \ t_{b(p+q)}^{(N)} & a_{b1}^{(N)} \ ... \ a_{b(p+q)}^{(N)}  & v_s^{(N)} & v_b^{(N)} & v_e^{(N)} \\ 
\end{bmatrix} \\
\end{aligned} \\
$$
\linebreak
\begin{itemize}
\item Initialize the first row with user-given initial values:
$$
\begin{aligned}
\phi^{(1)} &= (\phi_1^{(1)}, \ ... , \ \phi_d^{(1)}) \\
& = (k_1, \ ... , \ k_q, \ \theta_0^{(1)}, \ ..., \ \theta_0^{(p+q)}, \ \alpha_0^{(1)}, \ ..., \ \alpha_0^{(p+q)}, \  \theta_0^{(1)}, \ ..., \ \theta_0^{(p)}, \ \alpha_0^{(p)}, \ ..., \ \alpha_0^{(p)}, \ v_{s0}, \ v_{b0}, \ v_{e0})
\end{aligned}
$$
\item At each iteration $i \in (2, ... , N)$, and to update $j$-th parameter ($j \in (1,... , d)$ and first ($j-1$) parameters are already updated):   
\begin{enumerate} 
\item Propose a new value for $\phi_j^{(i)}$ based on its last update \footnote{If it is first parameter, the last update is the last row (\(\phi^{(i-1)}\)}. $\phi_j^{(i-1)}$, called Metropolis Update (MU): 

$$\phi_j^* = \mathcal{N}(\phi_j^{(i-1)}, \sigma^2_p)$$
    
where $\sigma^2_p$ is adaptively adjusted \footnote{Every \(50\) iterations, acceptance rate (\(AR\)) is computed. If \(AR < 0.44\), proposal variance \(\sigma^2_p\) is decreased, and vice versa.It has been shown that for one-dimensional proposals used in Metropolis within Gibbs algorithm, the optimal acceptance rate is \(0.44\) (TODO: citation).} to ensure (faster) convergence.    

\item Form the parameter vector $\phi^*$ based on MU:
$$
\begin{aligned}
\phi^{(last)} &= (\phi_1^{(i)}, ... , \phi_{j-2}^{(i)}, \ \phi_{j-1}^{(i)}, \ \phi_{j}^{(i-1)}, \ ..., \ \phi_d^{(i-1)})    \\
\phi^*  \ \ \ \ \ &= (\phi_1^{(i)}, ... , \phi_{j-1}^{(i)}, \ \phi_{j}^*, \ \phi_{j+1}^{(i-1)}, \ ..., \ \phi_d^{(i-1)})    \\
\end{aligned}
$$ 
\item Draw a random sample $u$ from $U(0, 1)$ and take its log: $\ln(u)$  
\item Compute the difference between the log of joint posterior density given current   and last parameter vectors:
$$
h(\phi^*, \ \phi^{(last)}) = \ln(L(\text{y} | \phi^*)) + \ln(\mathcal{P}[\phi^*]) - \ln(L(\text{y} | \phi^{(last)})) + \ln(\mathcal{P}[\phi^{(last)}])
$$
\item If $h(\phi^*, \ \phi^{(last)}) > \ln(u)$, set:
$$
\phi_{j}^{(i)} = \phi_{j}^*
$$
otherwise, 
$$
\ \quad \phi_{j}^{(i)} = \phi_{j}^{(i-1)}
$$
\item The update vector now is:
$$
\phi^{(last)}  \ = (\phi_1^{(i)}, ... , \phi_{j-1}^{(i)}, \ \phi_{j}^{(i)}, \ \phi_{j+1}^{(i-1)}, \ ..., \ \phi_d^{(i-1)})
$$

\end{enumerate}

\item If all parameters are updated, go to next iteration of $i$
\item When iterations of $i$ is completed, return the matrix of $\Phi$ that 
contains joint posterior density distribution of all parameters. Marginal 
distribution of each parameter can be used for point prediction and uncertainty 
quantification (credible interval).
\end{itemize}

\rule{\linewidth}{2pt}


### Bayesian Analysis

In the Bayesian framework, the joint probability distribution of all parameters
and hyperparameters of the calibration model given data 
($\mathcal{P}[\phi|\text{y}]$) can be derived:

$$
\begin{aligned}
L(\text{y} | \phi)                   &=       \ |C|^{-\frac{1}{2}} \ . \ e^{- \frac{1}{2} \text{y}. C^{-1}.\text{y}^T} \\ \\
\mathcal{P} [\phi | \text{y}]      \ &\propto \ L(\text{y} | \phi) \ .  \ \mathcal{P} [\phi]                         \\ \\
\end{aligned}
$$

Taking the log from both sides will decrease computational load and increase speed:


$$
\begin{aligned}
\ln(L(\text{y} | \phi))             &= -\frac{1}{2} \ln(|C|) - \frac{1}{2} \text{y}. C^{-1}.\text{y}^T \\ \\
\ln(\mathcal{P}[\phi | \text{y}])\  &\propto   \ln(L(\text{y} | \phi)) + \ln(\mathcal{P}[\phi) \\ \\
&= -\frac{1}{2} \ln(|C|) - \frac{1}{2} \text{y}. C^{-1}.\text{y}^T   &\text{(log liklihood given full data: X, y)} \\ \\
&+ \ \sum_{i=1}^q \mathcal{P}[\kappa_i]  &\text{(priors for calibration parameters)}\\ \\ 
&+ \ \sum_{i=1}^{p+q} \mathcal{P}[\theta_{si}] + \sum_{i=1}^{p+q} \mathcal{P}[\alpha_{si}] \ + \ \mathcal{P}[\sigma_s^2] \ &\text{(priors for} \ \eta(.) \ \text{hyperparameters)} \\ \\
&+ \ \sum_{i=1}^{p} \mathcal{P}[\theta_{bi}] + \sum_{i=1}^{p} \mathcal{P}[\alpha_{bi}] \ + \ \mathcal{P}[\sigma_b^2] &\text{(priors for} \ \delta_{\kappa}(.) \ \text{hyperparameters)}\\ \\
&+ \ \mathcal{P}[\sigma_{\epsilon}^2] &\text{(prior for measurement error variance)}\\ 
\end{aligned}
$$

Above equation is intractable and thus a simulation-based method must be used to
sample from posterior distribution. `FBC` implements a version Markov Chain 
Monte Carlo (MCMC) simulation. 


### Full Model

KOH model combines simulation and field data to form a joint model
using joint dataset:



$$
\begin{aligned}
\bf{y} \  &=
\begin{bmatrix}
\bf{y}_f  \\
\bf{y}_s \\ 
\end{bmatrix} = 
(y_1, \ ... \ , \ y_{n+m})^T
\ \ &\text{(joint vector of responses)}  \\ \\
\bf{X}   &=
\begin{bmatrix}
\bf{X_{\kappa}}  \\
\bf{X_s} \\
\end{bmatrix} =
\begin{bmatrix}
\bf{x_1} \ \bf{x_2}  \ ... \ \bf{x_{p+q}} \\
\end{bmatrix}
\ \ &\text{(joint input matrix)}  \\ \\
\bf{x_i} &= (x_1, x_2, ..., x_{n+m}) \ \ \ \ \forall i \in \{1, 2, ..., p+q\} \ 
\end{aligned}
$$


Since $\bf{X_f}$ is a sub-matrix of $\bf{X}$, we can represent the functional
relationship between input and response for full model with $\zeta(.)$, which
is considered to be realization of a random function and derived from $\eta(.)$
and $\delta_{\kappa}(.)$:

$$
\bf{y} \  = \zeta(\bf{X}) 
$$



Because both $\eta(.)$ and $\delta_{\kappa}(.)$ are GPs, $\zeta(.)$ can also be
considered a zero mean GP:

$$
\zeta(.) \sim GP \ (0, \ C(., .))
$$

Where covariance matrix $C$ is dependent to full input matrix $\bf{X}$ and 
hyperparameters of $\eta(.)$ and $\delta_{\kappa}(.)$, which in turn are 
dependent to model hyperparameters.


Using matrix notation, input/output relationship of both experiments is 
presented below in matrix notation. This relations can be used to derive a 
relationship for joint data:

$$
\text{y} \  = \zeta(X) = \eta(X) + 
\begin{bmatrix}
\delta_{\kappa}(X_f) +\mathcal{E} & 0 \\
0 & 0  \\
\end{bmatrix} \\ \\
$$

Since $X_f$ is a subset of matrix $X$, the joint response $\text{y}$ can be 
modeled as a random function $\zeta(X)$. In the ball example, the full input 
matrix $X$ is a $(163 \times 2)$ matrix by stacking $X_{\kappa}$ on $X_s$ Note 
that since $\kappa$ is unknown, the initial value `c0` is used internally to 
build $X_{\kappa}$.






Where, $C_{\eta}$ and $C_{\delta}$ are covariance matrices of full input matrix 
$X$ and original field input matrix $X_f$. And $C$ is characterized as:


$$
\begin{aligned}
C(X, X) &=  C_{\eta}(X, X) + 
\begin{bmatrix} 
C_{\delta}(X_f, X_f) + \sigma^2_{\epsilon}.I_n   &   0  \\
0                               &   0 \\
\end{bmatrix} \\ \\ 
&= 
\begin{bmatrix} 
C_{\eta}(X_f, X_f) + C_{\delta} + \sigma^2_{\epsilon}.I_n & C_{\eta}(X_s, X_f)  \\ 
C_{\eta}(X_f, X_s)                                & C_{\eta}(X_s, X_s)  \\
\end{bmatrix}
\end{aligned}
$$

Note that $C_s$ covariance matrix of the full data $X$ is further divided to
components $C_{\eta}(X_f, X_f)$, $C_{\eta}(X_s, X_f)$, $C_{\eta}(X_f, X_s)$, and
$C_{\eta}(X_s, X_s)$ to optimize the computation. 






\pagebreak

## References












Vignettes are long form documentation commonly included in packages. Because they are part of the distribution of the package, they need to be as compact as possible. The `html_vignette` output type provides a custom style sheet (and tweaks some options) to ensure that the resulting html is as small as possible. The `html_vignette` format:

- Never uses retina figures
- Has a smaller default figure size
- Uses a custom CSS stylesheet instead of the default Twitter Bootstrap style

## Vignette Info

Note the various macros within the `vignette` section of the metadata block above. These are required in order to instruct R how to build the vignette. Note that you should change the `title` field and the `\VignetteIndexEntry` to match the title of your vignette.

## Styles

The `html_vignette` template includes a basic CSS theme. To override this theme you can specify your own CSS in the document metadata as follows:

    output: 
      rmarkdown::html_vignette:
        css: mystyles.css

## Figures

The figure sizes have been customised so that you can easily put two images side-by-side. 


```r
plot(1:10)
plot(10:1)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-2.png)

You can enable figure captions by `fig_caption: yes` in YAML:

    output:
      rmarkdown::html_vignette:
        fig_caption: yes

Then you can use the chunk option `fig.cap = "Your figure caption."` in **knitr**.

## More Examples

You can write math expressions, e.g. $Y = X\beta + \epsilon$, footnotes, and tables, e.g. using `knitr::kable()`.



|                  |  mpg| cyl|  disp|  hp| drat|    wt|  qsec| vs| am| gear| carb|
|:-----------------|----:|---:|-----:|---:|----:|-----:|-----:|--:|--:|----:|----:|
|Mazda RX4         | 21.0|   6| 160.0| 110| 3.90| 2.620| 16.46|  0|  1|    4|    4|
|Mazda RX4 Wag     | 21.0|   6| 160.0| 110| 3.90| 2.875| 17.02|  0|  1|    4|    4|
|Datsun 710        | 22.8|   4| 108.0|  93| 3.85| 2.320| 18.61|  1|  1|    4|    1|
|Hornet 4 Drive    | 21.4|   6| 258.0| 110| 3.08| 3.215| 19.44|  1|  0|    3|    1|
|Hornet Sportabout | 18.7|   8| 360.0| 175| 3.15| 3.440| 17.02|  0|  0|    3|    2|
|Valiant           | 18.1|   6| 225.0| 105| 2.76| 3.460| 20.22|  1|  0|    3|    1|
|Duster 360        | 14.3|   8| 360.0| 245| 3.21| 3.570| 15.84|  0|  0|    3|    4|
|Merc 240D         | 24.4|   4| 146.7|  62| 3.69| 3.190| 20.00|  1|  0|    4|    2|
|Merc 230          | 22.8|   4| 140.8|  95| 3.92| 3.150| 22.90|  1|  0|    4|    2|
|Merc 280          | 19.2|   6| 167.6| 123| 3.92| 3.440| 18.30|  1|  0|    4|    4|



Also a quote using `>`:

> "He who gives up [code] safety for [code] speed deserves neither."
([via](https://twitter.com/hadleywickham/status/504368538874703872))


#' ## Calibration Model
#'
#' This implementation is based on Kennedy-O'Hagan (KOH) calibration model. In
#' their seminal paper*, Kennedy and O'Hagan augmented simulator output with
#' field observation and fit a three-component model that accounts for simulator
#' input and bias correction using two independent Gaussian Processes (GP) and a
#' third term representing measurement error.
#' \deqn{z_i = \eta (x_i, \kappa) + \delta(x_i) + e_i}
#' In the original paper, the authors proposed a two-stage hierarchical Bayesian
#' model. In the first stage, point estimates of GP hyperparameters are computed
#' using maximum likelihood estimation (MLE) method. In the second stage, these
#' hyperparameters are fixed at their estimated value and run a MCMC algorithm
#' to sample calibration parameters.
#' In contrast, `calibrate` runs MCMC algorithm to sample the posterior
#' distribution of all parameters/hyperparameters. As a result, priors must be
#' specified carefully to reflect the prior expert belief about the distribution
#' and initial values must be chosen as close as possible to prior means.








### Introduction

### KOH Calibration Model

Kennedy and O'Hagan (KOH) proposed a Bayesian calibration model to integrate the 
information from field observations, $y_f$, (and its corresponding matrix of 
experimental input, $X_f$) and computer simulation, $y_s$ (and its corresponding 
matrix of experimental/calibration input, $X_s$). Often, in computer simulation
there are extra parameters that are either are unkwon in the field setting or
are tuning parameters specific to simulation without physical interpretation. 
These inputs are calibration parameters and tuning them for the computer model is
one of the main objective of calibration. 

Moreover, KOH assumes both field and simulation results are realization of real 
process, $Y$, while admitting to a systemic bias in simulation models, and 
measurement error in the form of independent zero-mean Gaussian error terms with
variance $\sigma^2_{\epsilon}$.:  


$$
\begin{eqnarray}
Y|_{X_f} : y(x) &=& y_f(x_f) +  \epsilon \ \ &, \text{where} \ x_f: 1 \times p \ \text{vector of experimental inputs} \\
Y|_{X_s} : y(x) &=& y_s(x_s) + y_b(x_b) + \epsilon \ \ &, \ \text{where} \ x_s: 1 \times p+q \ \text{vector of experimental + calibration inputs}
\end{eqnarray}
$$

In this formulation, $x_f$ is in fact $1 \times p+q$ vector, similar to $x_s$, 
because it is realization of the same real process. However, calibration 
parameters are set (fixed or variable) at their best value but their value is 
not known. Since all field observation have the same calibration parameters,
they are removed from notation. The current text uses $(1 \times p+q)$ vector 
$x_s$ for simulation input, $(1 \times p+q)$ vector $x_f$ for field data 
augmented with calibration parameters (at their current updated value and same 
for all field observations), and $(1 \times p)$ vector $x_b$ for field data 
without augmenting calibration parameters. The subscript $b$ is chosen to 
indicate that it is only input to bias term. Stacking the two data, we can write
the KOH calibration model in a single expression:

$$
Y|_{X} : y(x) = y_s(x) + y_b(x) + \epsilon \ \ , \ \text{where} \ x: 1 \times p+q \ \text{vector of experimental + calibration inputs}
$$

this formulation, emphasizes that model has systemic bias, even under best 
calibration parameters. However, KOH postulates that reasonable correction can
be made through $y_b(.)$ 

Above specification leads to three unknowns: calibration parameters (denoted by
$\kappa$), variance of measurement error ($\sigma^2_{\epsilon}$) and functional
form of bias correction term (denoted by $y_b(.)$ noting that the actual bias is
$-y_b(.)$). Moreover, KOH assumes a Gaussian Process (GP) for both simulation
term (often the simulation is also expensive leading to limited runs. Fitting a 
GP as a surrogate to couble bias-correction GP) and bias-correction term, which 
introduces new hyper-parameters to the model. 

Specifically, assuming a generalized correlation structure such as power 
exponential correlation family, the scale (denoted by $\theta_s$ and $\theta_b$),
smoothness (denoted by $\alpha_s$ and $\alpha_b$), and their marginal variance
(denoted by $\sigma_s^2$ and $\sigma_b^2$) are added to the model. 

$$
\Omega = (Y_s(.), \Delta_b(.), \kappa, \sigma_{\epsilon}^2 ) = (\kappa, \theta_s, \alpha_s, \sigma_s^2, \theta_b, \alpha_b, \sigma^2_b, \sigma_{\epsilon}^2)
$$

After loading of the data (either toy data or user data), the combined 
observation vector is scaled by mean and standard deviation of simulation runs 
such that $y_s$ terms have mean zero and standard deviation of 1 after scaling.
Furthermore, input matrix is scaled according to mean and standard deviation of
columns of the input matrix (consisting of experimental and calibration values)
such that $X_s$ columns span $[0, 1]$. These scaling removes the need for 
parametrization of the mean of the GPs as they will be considered zero. 
Following the results from Chen et al. (2016), the fitting a linear regression
model does not necessarily improve the performance as the GP models are flexible.
Therefore current package does not consider fitting regression model on the mean
GPs.

### Bayesian Inference

The framework of KOH model is Bayesian. Prior information and restrictions about
model parameters and hyperparameters can be specified through informative or
non-informative priors. 

For calibration parameters, often non-informative uniform prior is used over 
wide but finite domain, unless expert opinion is available in the form of prior
distributions. An alternative choice can be a regularizing prior that prevent 
from over-concentration of posterior density on boundaries of domain.

... other pror specification

The likelihood for posterior density is in following form along with its log 
form (where both sides are log transformed for computational efficiency):


$$
\begin{align}
&L(\Omega \mid Y) = 
\end{align}
$$
where augmented covariance matrix is given by:

$$
\Sigma \mid \Omega = 
$$

This expression is interactable and must be solved through optimization methods
or it can be used to sample from its joint posterior.


### Monte Carlo Markov Chain Sampling

### Metropolis-Within-Gibbs Algorithm Imlementation

### Adaptive Proposal

### Parameter Estimation

### Prediction

### Toy Example

##### Data

##### Model

##### MCMC Runs

##### Parameter Estimation

##### Prediction

##### Diagnostics

### Future Development

### Consclusion






(denoted by $\alpha$) and marginal variance












