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

The packages required in the applications are also included in such a file.

## Simulation Study $1$

We start from the simulation process of the artificial network we focus on in the paper.
Recall here that we randomly generated a network from the ZINB-SBM with the settings: $N = 75$, $K = 3$, $p = 0.15$, $\boldsymbol{\Pi} = \left(0.45, 0.35, 0.2\right)$ and 

$$\boldsymbol{R} = \begin{pmatrix}0.10 & 0.05 & 0.05 \\\0.02 & 0.60 & 0.05 \\\0.02 & 0.10 & 2.50 \end{pmatrix}, \boldsymbol{Q} = \begin{pmatrix}0.10 & 0.20 & 0.25 \\\0.15 & 0.15 & 0.25 \\\ 0.15 & 0.20 & 0.35\end{pmatrix}.$$





$\boldsymbol{\Pi} = \left(0.35, 0.20, 0.45\right)$ and 

$$\boldsymbol{R} = \begin{pmatrix}0.60 & 0.05 & 0.02 \\\0.10 & 2.50 & 0.02 \\\0.05 & 0.05 & 0.10 \end{pmatrix}, \boldsymbol{Q} = \begin{pmatrix}0.15 & 0.25 & 0.15 \\\0.20 & 0.35 & 0.15 \\\ 0.20 & 0.25 & 0.10\end{pmatrix}.$$



