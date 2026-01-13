# Title: A Bayesian Joint Modelling for Misclassified Interval-censoring and Competing Risks

**Author**: Zhenwei Yang 

**Affiliation**: Department of Biostatistics, Erasmus Medical Center

**DOI**: 10.48550/arXiv.2404.09362

****

## Content

* Misclassification-ICJM.Rproj
* [Rscript](#Rscript)
* [Output](#Output)

****

### Rscript

* Data analysis

|File name| Description|
|:----------:|:--------------|
|MCICJM_fixed.R| the 4 MCICJMs fitted on the Canary PASS data (assuming 60%, 75%, 80% and 100% sensitivity)|
|MCICJM_inits.R | prepare for the inits of the longitudinal submodel |
|MCICJM_prior.R | the MCICJM using uniform prior for sensitivity (0.5, 0.9) |

* Data preprocessing

|File name| Description|
|:----------:|:--------------|
|biopsy_info.R | generate a dataset recording the biopsy schedule |

* Functions

|File name| Description|
|:----------:|:--------------|
|function.R | help and auxiliary functions|
|MCICJM fit_fixed.R | a function to fit the MCICJM assuming a fixed sensitivity |
|MCICJM fit_prior.R | a function to fit the MCICJM using a prior for sensitivty |
|MCICJM_mixed.R | a function to fit a mixed model to prepare for the inits of the longitudinal submodel in the MCICJM |

* Simulation

|File name| Description|
|:----------:|:--------------|
|Dataset generation.R | training sets generation based on the MCICJM with 75% fixed sensitivity |
|MCICJM_fixed60.R | fitting the MCICJM assuming 60% fixed sensitivity based on 200 simulated sets|
|MCICJM_fixed75.R | fitting the MCICJM assuming 75% fixed sensitivity based on 200 simulated sets|
|MCICJM_fixed100.R | fitting the MCICJM assuming 100% fixed sensitivity based on 200 simulated sets|
|MCICJM_unif50_80.R | fitting the MCICJM assuming a sensitivity prior Unif(0.5, 0.8) based on 200 simulated sets|
|MCICJM_unif60_90.R | fitting the MCICJM assuming a sensitivity prior Unif(0.6, 0.9) based on 200 simulated sets|
|MM_fit.R | fitting mixed models to prepare for the inits of the longitudinal submodels in the MCICJMs |

* Tables&figures

|File name| Description|
|:----------:|:--------------|
|Fig1.ppt | for Figure 1 in the manuscript|
|Figures_main.R | code to generate figures in the manuscript|
|Figures_supplement.R | code to generate figures in the supplementary materials|


### Output

* Data analysis

|File name| Description|
|:----------:|:--------------|
|inits | inits of the longitudinal submodel in the MCICJM fitted on the Canary PASS data |
|MCICJM_0.6.R | model results (including the MCMC samples) of the MCICJM with 60% fixed sensitivity |
|MCICJM_0.8.R | model results (including the MCMC samples) of the MCICJM with 80% fixed sensitivity |
|MCICJM_0.75.R | model results (including the MCMC samples) of the MCICJM with 75% fixed sensitivity |
|MCICJM_1.R | model results (including the MCMC samples) of the MCICJM with 100% fixed sensitivity |
|MCICJM_unif5090.R | model results (including the MCMC samples) of the MCICJM with a sensitivity prior Unif(0.5, 0.9) |

* Plots

|File name| Description|
|:----------:|:--------------|
|Main text | all figures used in the manuscripts|
|summary storage| storage datasets for generating the plots|
|Supplementary| all figures used in the supplementary materials|

### Package dependencies

- rajgs
- mcmcplots
- GLMMadaptive
- ggplot2
- tidyverse
- splines
- future
- mcmcse
- mvtnorm
- JMbayes2
- JMbayes
- MASS
- doParallel
- truncnorm
- Matrix
- latex2exp
- cowplot
- PCaASSim

> [!Note]
> - Please make sure to install above-mentioned packages by `install.packages()` before running the R code
