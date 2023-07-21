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
We randomly generated a network from the ZINB-SBM with the settings: $N = 75$, $K = 3$, $p = 0.15$, $\boldsymbol{\Pi} = \left(0.45, 0.35, 0.2\right)$ and 

$$\boldsymbol{R} = \begin{pmatrix}0.10 & 0.05 & 0.05 \\\0.02 & 0.60 & 0.05 \\\0.02 & 0.10 & 2.50 \end{pmatrix}, \boldsymbol{Q} = \begin{pmatrix}0.10 & 0.20 & 0.25 \\\0.15 & 0.15 & 0.25 \\\ 0.15 & 0.20 & 0.35\end{pmatrix}.$$

Note here that the label switching method we applied in the paper is the one introduced in [Rastelli and Friel (2018)](https://doi.org/10.1007/s11222-017-9786-y), so that, the uniqueness of each clustering state is assured. 
Such a label switching method can be described as: assign the first node to the cluster $1$ and then iteratively assign the next node either to a new empty cluster or an existing cluster.
We applied the label switching on the simulated data to ensure the uniqueness of the clsutering state.
This also leads to the label-switching of the clustering dependent parameters, $\boldsymbol{\Pi},\boldsymbol{R},\boldsymbol{Q}$, which are used for generating the network shown above.
This is the reason that the initial parameters are different from the ones we illustrated in the paper.
The data we show in the paper is already label-switched.

The simulation can be implemented via the code shown below:

``` r
# Simulation study 1 of ZINB-SBM N = 75, K = 3
# Simulate an artificial data from the ZINB-SBM
SS1_ZINBSBM_N75_K3 <- 
  Simulation_Directed_ZINBSBM(N = 75 ,K = 3 ,Pi = c(0.45,0.35,0.2),p = 0.15,
                              R = matrix(c(0.1,0.02,0.02,
                                           0.05,0.6,0.1,
                                           0.05,0.05,2.5),3,3),
                              Q = matrix(c(0.1,0.15,0.15,
                                           0.2,0.15,0.2,
                                           0.25,0.25,0.35),3,3))
# # Store the simulation study 1 artificial dataset
# write.csv(SS1_ZINBSBM_N75_K3$Y,"SS1_ZINBSBM_N75K3_obsY.csv", row.names = FALSE)
# write.csv(SS1_ZINBSBM_N75_K3$X,"SS1_ZINBSBM_N75K3_obsX.csv", row.names = FALSE)
# write.csv(SS1_ZINBSBM_N75_K3$nu,"SS1_ZINBSBM_N75K3_obsnu.csv", row.names = FALSE)
# write.csv(SS1_ZINBSBM_N75_K3$Z,"SS1_ZINBSBM_N75K3_obsZ.csv", row.names = FALSE)
```

The artificial network data we used in the paper are included in the files within this repository and can be loaded by:

``` r
# # Load the simulation study 1 artificial dataset
SS1_ZINBSBM_N75_K3 <- 
  list(Y = as.matrix(read.csv("SS1_ZINBSBM_N75K3_obsY.csv",header = TRUE)),
       X = as.matrix(read.csv("SS1_ZINBSBM_N75K3_obsX.csv",header = TRUE)),
       nu = as.matrix(read.csv("SS1_ZINBSBM_N75K3_obsnu.csv",header = TRUE)),
       Z = as.matrix(read.csv("SS1_ZINBSBM_N75K3_obsZ.csv",header = TRUE)))
colnames(SS1_ZINBSBM_N75_K3$Y) <- c()
colnames(SS1_ZINBSBM_N75_K3$X) <- c()
colnames(SS1_ZINBSBM_N75_K3$nu) <- c()
colnames(SS1_ZINBSBM_N75_K3$Z) <- c()
```

Then we apply label switching on the latent clustering $\boldsymbol{z}$ and those clustering dependent parameters, $\boldsymbol{\Pi},\boldsymbol{R},\boldsymbol{Q}$, of the simulated network, and obtain the initial $\boldsymbol{z},\boldsymbol{\Pi},\boldsymbol{R},\boldsymbol{Q}$.
Recall here that we treat the label-switched initial clustering and parameters as the "true" references in the experiments.

``` r
#--------------------------------------------------------------------------------------------------------------------------------------------
# Label switch Z
SS1_ZINBSBM_N75_K3_LSZ <- LabelSwitching_SG2003_ZINBSBM(Z = list(SS1_ZINBSBM_N75_K3$Z))$Z[[1]]
image(t(SS1_ZINBSBM_N75_K3$Y)[order(SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)),rev(order(SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)))])
group_counts <- (as.numeric(table(SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3))))
abline(v = -1/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1))))
abline(h = 1-(-1/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)))))
#--------------------------------------------------------------------------------------------------------------------------------------------
# Label switch initial R,Q,Pi and evaluate the initial mean and variance
res <- LabelSwitching_SG2003_ZINBSBM(Z = list(SS1_ZINBSBM_N75_K3$Z),
                                     Pi = list(c(0.45,0.35,0.2)),
                                     R = list(matrix(c(0.1,0.02,0.02,
                                                       0.05,0.6,0.1,
                                                       0.05,0.05,2.5),3,3)),
                                     Q = list(matrix(c(0.1,0.15,0.15,
                                                       0.2,0.15,0.2,
                                                       0.25,0.25,0.35),3,3)),
                                     Acceptance_count_R = list(matrix(0,3,3)))
SS1_ZINBSBM_N75_K3_obs_InitialR <- res$R[[1]]
SS1_ZINBSBM_N75_K3_obs_InitialQ <- res$Q[[1]]
SS1_ZINBSBM_N75_K3_obs_InitialPi <- res$Pi[[1]]
```

The label-switched latent parameters now agree with the ones we illustrated in the paper: $\boldsymbol{\Pi} = \left(0.35, 0.20, 0.45\right)$ and 

$$\boldsymbol{R} = \begin{pmatrix}0.60 & 0.05 & 0.02 \\\0.10 & 2.50 & 0.02 \\\0.05 & 0.05 & 0.10 \end{pmatrix}, \boldsymbol{Q} = \begin{pmatrix}0.15 & 0.25 & 0.15 \\\0.20 & 0.35 & 0.15 \\\ 0.20 & 0.25 & 0.10\end{pmatrix}.$$

Once we obtained the label-switched latent clustering and parameters, we can evaluate the initial mean and variance of the distribution assumed for the non-missing weights.

``` r
SS1_ZINBSBM_N75_K3_obs_Initialmean <- SS1_ZINBSBM_N75_K3_obs_InitialR*(1-SS1_ZINBSBM_N75_K3_obs_InitialQ)/SS1_ZINBSBM_N75_K3_obs_InitialQ
SS1_ZINBSBM_N75_K3_obs_Initialvar <- SS1_ZINBSBM_N75_K3_obs_InitialR*(1-SS1_ZINBSBM_N75_K3_obs_InitialQ)/SS1_ZINBSBM_N75_K3_obs_InitialQ^2
```



