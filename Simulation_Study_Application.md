# Simulation Study Applications
This markdown file illustrates the code and applications for the simulation studies shown in the *Zero-Inflated Negative-Binomial Stochastic Block Model* paper.
There are two simulation studies experimented in the paper.
The simulation study $1$ focuses on an artificial network randomly generated from the zero-inflated Negative-Binomial stochastic block model (ZINB-SBM) and the simulation study $2$ (SS2) focuses on an artificial network randomly generated from the zero-inflated Poissin stochastic block model (ZIP-SBM).
Both simulation studies fit both ZINB-SBM and ZIP-SBM to the artificial datasets.

The source function code for the implementations of the inference of both models is included in the file [`Functions_for_ZINB_SBM.R`].
The explanations of the code are written as comments beside the code.
Note here that the practitioners need to load the source functions everytime after cleaning the enviroment.

``` r
rm(list=ls()) # remove all in the enviroment
gc() # Free unused memory
source("Functions_for_ZINB_SBM.R") # load the functions
```

The packages required are also included in such a file.

