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
Such a label switching method can be described as: assigning the first node to the cluster $1$, and then iteratively assigning the next node either to a new empty cluster or an existing cluster.
We applied the label switching on the simulated data to ensure the uniqueness of the clsutering.
This also leads to the label-switching of the above clustering dependent parameters, $\boldsymbol{\Pi},\boldsymbol{R},\boldsymbol{Q}$, for the generation of the network.

The simulation of the ZINB-SBM can be implemented via the function `Simulation_Directed_ZINBSBM()` included in the source code file [`Functions_for_ZINB_SBM.R`] and the code is shown below:

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

The artificial network data we focus on in the paper is included in the files within this repository and can be loaded by:

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

Then we apply label switching on the latent clustering $\boldsymbol{z}$ and those clustering dependent parameters, $\boldsymbol{\Pi},\boldsymbol{R},\boldsymbol{Q}$, of the simulated network. 
The label switching methods for the ZINB-SBM can be implemented by the function `LabelSwitching_SG2003_ZINBSBM()` included in the [`Functions_for_ZINB_SBM.R`].
Recall here that we treat the label-switched initial clustering and the clustering dependent parameters as the "true" references in the experiments.

``` r
#--------------------------------------------------------------------------------------------------------------------------------------------
# Label switch Z
SS1_ZINBSBM_N75_K3_LSZ <- LabelSwitching_SG2003_ZINBSBM(Z = list(SS1_ZINBSBM_N75_K3$Z))$Z[[1]]
image(t(SS1_ZINBSBM_N75_K3$Y)[order(SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)),rev(order(SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)))])
group_counts <- (as.numeric(table(SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3))))
abline(v = -1/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1))))
abline(h = 1-(-1/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)))))
#--------------------------------------------------------------------------------------------------------------------------------------------
# Label switch initial R,Q,Pi
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

It can be checked that the label-switched latent parameters are now: $\boldsymbol{\Pi} = \left(0.35, 0.20, 0.45\right)$ and 

$$\boldsymbol{R} = \begin{pmatrix}0.60 & 0.05 & 0.02 \\\0.10 & 2.50 & 0.02 \\\0.05 & 0.05 & 0.10 \end{pmatrix}, \boldsymbol{Q} = \begin{pmatrix}0.15 & 0.25 & 0.15 \\\0.20 & 0.35 & 0.15 \\\ 0.20 & 0.25 & 0.10\end{pmatrix},$$

which argee with the ones we illustrated in the paper.

Once we obtained the label-switched latent clustering and parameters, we can evaluate the initial mean and variance of the distribution assumed for the non-missing weights between each pair of clusters.

``` r
#  Evaluate the initial mean and variance
SS1_ZINBSBM_N75_K3_obs_Initialmean <- SS1_ZINBSBM_N75_K3_obs_InitialR*(1-SS1_ZINBSBM_N75_K3_obs_InitialQ)/SS1_ZINBSBM_N75_K3_obs_InitialQ
SS1_ZINBSBM_N75_K3_obs_Initialvar <- SS1_ZINBSBM_N75_K3_obs_InitialR*(1-SS1_ZINBSBM_N75_K3_obs_InitialQ)/SS1_ZINBSBM_N75_K3_obs_InitialQ^2
```

We can also evaluate the initial $\boldsymbol{P_{m0}}$ (the probability of the zero interaction being missing zero) assumed for this network.

``` r
# Evaluate the initial P_m0
SS1_ZINBSBM_N75_K3_obs_InitialProbObs0Missing0 <- matrix(0,3,3)
for (k1 in 1:3){
  for (k2 in 1:3){
    SS1_ZINBSBM_N75_K3_obs_InitialProbObs0Missing0[k1,k2] <- SS1_ZINBSBM_N75_K3_obs_Initialp/(SS1_ZINBSBM_N75_K3_obs_Initialp + (1-SS1_ZINBSBM_N75_K3_obs_Initialp)*dnbinom(0,SS1_ZINBSBM_N75_K3_obs_InitialR[k1,k2],SS1_ZINBSBM_N75_K3_obs_InitialQ[k1,k2]))
  }
}
SS1_ZINBSBM_N75_K3_obs_InitialProbObs0Missing0
```

The reference probability of missing zero, $p$, is simply the initial setting:

``` r
# Set initial p
SS1_ZINBSBM_N75_K3_obs_Initialp <- 0.15
```

The plots of the adjacency matrix $Y$, the adjacency matrix $Y$ conditional on true clustering $\boldsymbol{z}^{\*}$, i.e $\boldsymbol{Y}|\boldsymbol{z}^{*}$, and the indicator of whether $y_{ij}$ is non-zero (dark color) or not (light color) conditional on $\boldsymbol{z}^{\*}$ shown as Figure $1$ of the paper can be recovered by:

``` r
par(mfrow=c(1,3),mai = c(0.05, 0.05, 0.2, 0.05),mgp=c(0.1,0.1,0))
image(t(SS1_ZINBSBM_N75_K3$Y),axes = FALSE,xlab = "",ylab = "",main = "Adjacency Matrix Y")
image(t(SS1_ZINBSBM_N75_K3$Y)[order(SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)),rev(order(SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)))],axes = FALSE,xlab = "",ylab = "",main = TeX(r'($Y|$ True $z^*$)',bold = TRUE))
group_counts <- (as.numeric(table(SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3))))
abline(v = -1/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1))))
abline(h = 1-(-1/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)))))

image(t(1*(SS1_ZINBSBM_N75_K3$Y!=0))[order(SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)),rev(order(SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)))],axes = FALSE,xlab = "",ylab = "",main = TeX(r'($Y\neq 0|$ True $z^*$)',bold = TRUE))
group_counts <- (as.numeric(table(SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3))))
abline(v = -1/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1))))
abline(h = 1-(-1/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)))))
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42),mgp=c(3,1,0))
```

### SS1 ZINB-SBM Implementations

The implementations of applying partially collapsed Metropolis within Gibbs algorithm (PCMwG) for the ZINB-SBM on the network are based on the function `Directed_ZINBSBM_PCMwG()` included in the source code file [`Functions_for_ZINB_SBM.R`].
Note that such a function will automatically label switch the initial clustering input to the function in order to ensure the uniqueness of the input clustering.
Other steps are exactly the same as the Algorithm $2$ stated in the paper.

The function `Directed_ZINBSBM_PCMwG_FixedZ()` aims to apply further inference conditional on the fixed summarized clustering as we disucssed in the paper.
Recall here that, within such a function, the inference step of the clustering is removed and instead the clustering is fixed at the summarized clustering we obtained from the outputs of the function `Directed_ZINBSBM_PCMwG()`.

The PCMwG algorithm for the ZINB-SBM is implemented for $40,000$ iterations for each fixed $K = 2,3,4,5$.
The $p$ prior setting for the function `Directed_ZINBSBM_PCMwG()` is $p \sim Beta(1,9)$ which can be changed by inputting prior parameters.
The prior settings of other parameters are set by default as we discussed in the paper, that is, $\boldsymbol{\Pi} \sim \text{Dirichlet}(\alpha, \dots, \alpha)$, $q_{gh} \sim \text{Beta}(\beta_{q1}, \beta_{q2})$ for $g,h=1,2,\dots,K$, and the prior distribution of $\boldsymbol{R}$ is simply positive uniform $\text{U}(0,\text{UpperBound})$ where the "UpperBound" here can be a big enough value so that the $\boldsymbol{R}$ prior term can be cancelled by the fraction in the acceptance ratio of the Metropolis-Hastings step.



