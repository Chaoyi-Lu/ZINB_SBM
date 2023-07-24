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
The $p$ prior setting for the function `Directed_ZINBSBM_PCMwG()` is $p \sim \text{Beta}(1,9)$ which can be changed by inputting prior parameters.
The prior settings of other parameters are set by default as we discussed in the paper, that is, $\boldsymbol{\Pi} \sim \text{Dirichlet}(\alpha, \dots, \alpha)$, $q_{gh} \sim \text{Beta}(\beta_{q1}, \beta_{q2})$ for $g,h=1,2,\dots,K$, and the prior distribution of $\boldsymbol{R}$ is simply positive uniform $\text{U}(0,\text{UpperBound})$ where the "UpperBound" here can be a big enough value so that the $\boldsymbol{R}$ prior term can be cancelled by the fraction in the acceptance ratio of the Metropolis-Hastings (M-H) step. Recall also here that the proprosal distrbution of $\boldsymbol{R}$ in the M-H step is $r'\_{gh} \sim \text{U}(\text{max}(0,r_{gh}^{(t-1)}-\epsilon),r_{gh}^{(t-1)}+\epsilon)$ for each pair of $g,h = 1,\dots,K$ where $r_{gh}^{(t-1)}$ is the current state of the $r_{gh}$ and the proposal epsilon $\epsilon$ here is tuned to be $0.175$ where such an epsilon will also be applied in the real data application.
The implementations are applied by the code shown below.

``` r
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Simulation study 1: apply the ZINB-SBM Metropolis within Gibbs algorithm for fixed K = 2, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS1_ZINBSBM_N75_K3_Fixed_K2_T40000_1 <- Directed_ZINBSBM_PCMwG(Y = SS1_ZINBSBM_N75_K3$Y, K = 2, T = 40000, eps_R = 0.175,beta1 = 1, beta2 = 9)
end.time <- Sys.time()
SS1_ZINBSBM_N75_K3_Fixed_K2_T40000_1_time <- end.time - start.time
SS1_ZINBSBM_N75_K3_Fixed_K2_T40000_1_time # Time difference of 1.868442 hours
# save.image("SS1_ZINBSBM_N75_K3_Fixed_K2_T40000_1_prior_p_Beta_1_9.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Simulation study 1: apply the ZINB-SBM Metropolis within Gibbs algorithm for fixed K = 3, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1 <- 
  Directed_ZINBSBM_PCMwG(Y = SS1_ZINBSBM_N75_K3$Y, K = 3, T = 40000, eps_R = 0.175,beta1 = 1, beta2 = 9)
end.time <- Sys.time()
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_time <- end.time - start.time
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_time # Time difference of 2.605677 hours
# save.image("SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_prior_p_Beta_1_9.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Simulation study 1: apply the ZINB-SBM Metropolis within Gibbs algorithm for fixed K = 4, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS1_ZINBSBM_N75_K3_Fixed_K4_T40000_1 <- Directed_ZINBSBM_PCMwG(Y = SS1_ZINBSBM_N75_K3$Y, K = 4, T = 40000, eps_R = 0.175,beta1 = 1, beta2 = 9)
end.time <- Sys.time()
SS1_ZINBSBM_N75_K3_Fixed_K4_T40000_1_time <- end.time - start.time
SS1_ZINBSBM_N75_K3_Fixed_K4_T40000_1_time # Time difference of 3.291847 hours
# save.image("SS1_ZINBSBM_N75_K3_Fixed_K4_T40000_1_prior_p_Beta_1_9.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# Simulation study 1: apply the ZINB-SBM Metropolis within Gibbs algorithm for fixed K = 5, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS1_ZINBSBM_N75_K3_Fixed_K5_T40000_1 <- Directed_ZINBSBM_PCMwG(Y = SS1_ZINBSBM_N75_K3$Y, K = 5, T = 40000, eps_R = 0.175,beta1 = 1, beta2 = 9)
end.time <- Sys.time()
SS1_ZINBSBM_N75_K3_Fixed_K5_T40000_1_time <- end.time - start.time
SS1_ZINBSBM_N75_K3_Fixed_K5_T40000_1_time # Time difference of 4.01542 hours
# save.image("SS1_ZINBSBM_N75_K3_Fixed_K5_T40000_1_prior_p_Beta_1_9.RData")
```

Once we obtained the outputs from each fixed $K$ case, we can apply label switching on the clustering and check the mixing performance of the clustering.
We can first obtain the summarized missing zero probability $\tilde{p}$, and the summarized $\tilde{\boldsymbol{\nu}}$ the entry of which indicates the proportion of the times for the corresponding $y_{ij}$ being inferred as a missing zero.
We can also obtain the summarized $\tilde{\boldsymbol{z}}$ and $\tilde{\boldsymbol{P_{m0}}}$ to evaluate the model selection criterion, integrated classification log-partially-collapsed-likelihood (ICPCL), by maximizing ICPCL with respect to $\boldsymbol{R}$ as we discussed in Section $3.2$ of the paper.
Based on the model selection criterion values, we determine which $K$ case is the best one and apply further summary statistics for the analysis.

As it's shown in the Section $4.1$ of the paper, the $K=3$ case is the best pick and we take this case as an example here.
We start from the label switching process once we obatined the outputs `SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1` from the implementation.

``` r
# # Simulation study 1: apply the ZINB-SBM Metropolis within Gibbs algorithm for fixed K = 3, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
# rm(list=ls())
# gc()
# source("Functions_for_ZINB_SBM.R")
# start.time <- Sys.time()
# SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1 <- 
#   Directed_ZINBSBM_PCMwG(Y = SS1_ZINBSBM_N75_K3$Y, K = 3, T = 40000, eps_R = 0.175,beta1 = 1, beta2 = 9)
# end.time <- Sys.time()
# SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_time <- end.time - start.time
# SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_time # Time difference of 2.605677 hours
# # save.image("SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_prior_p_Beta_1_9.RData")

# Apply label switching
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS <-
  LabelSwitching_SG2003_ZINBSBM(Z = SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Z,
                                Pi = SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Pi,
                                R = SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$R,
                                Q = SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Q,
                                Acceptance_count_R = SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Acceptance_count_R)
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Z <- c() # delete the original data in order to save memory
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Pi <- c()
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$R <- c()
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Q <- c()
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Acceptance_count_R <- c()
gc() # Free unused R memory
```

After the label switching, we can check the mixing of the clustering posterior samples by evaluating the rand index between each iteration's $\boldsymbol{z}^{(t)}$ and the true clustering $\boldsymbol{z}^*$.

``` r
## check rand index for each iteration
require("fossil")
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI <- c()
for (t in 1:40001){
  SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI <-
    c(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI,
      rand.index(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]%*%c(1:ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]])),
                 SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3))) # calculate and store the rand index between the t's clustering and the reference clustering
}
plot(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI,type = "l",xlab = "",ylab = "", main = "Rand Index",cex.axis = 0.8) # obtain the rand index plot 
```

Similar rand index plots can also be obtained for $K=2,4,5$ cases and also for the ZIP-SBM cases. 
The code for the ZIP-SBM cases will follow in the next subsection of this file.
The plots of the rand index mixing performance shown as Figure $2$ of the paper can be recovered by:

``` r
par(mfrow=c(2,2),mai = c(0.3, 0.3, 0.15, 0.05), mgp=c(0.9,0.2,0))
plot(SS1_ZINBSBM_N75_K3_Fixed_K2_T40000_1_LSZ_RI,type = "l",xlab = "",ylab = "", main = "",cex.axis = 0.8, ylim = c(0.35,1),col = 2,lty = 1)
lines(SS1_ZIPSBM_N75_K3_Fixed_K2_T40000_1_LSZ_RI,type = "l", col = 1,lty = 2)
title(xlab = "Iteration",ylab = "RI", main = "K=2 Cases Rand Index", mgp=c(0.9,0.1,0),cex.main=0.8,cex.lab = 0.8)

plot(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI,type = "l",xlab = "",ylab = "", main = "",cex.axis = 0.8, ylim = c(0.35,1),col = 2,lty = 1)
lines(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI,type = "l", col = 1,lty = 2)
title(xlab = "Iteration",ylab = "RI", main = "K=3 Cases Rand Index", mgp=c(0.9,0.1,0),cex.main=0.8,cex.lab = 0.8)

plot(SS1_ZINBSBM_N75_K3_Fixed_K4_T40000_1_LSZ_RI,type = "l",xlab = "",ylab = "", main = "",cex.axis = 0.8, ylim = c(0.35,1),col = 2,lty = 1)
lines(SS1_ZIPSBM_N75_K3_Fixed_K4_T40000_1_LSZ_RI,type = "l", col = 1,lty = 2)
title(xlab = "Iteration",ylab = "RI", main = "K=4 Cases Rand Index", mgp=c(0.9,0.1,0),cex.main=0.8,cex.lab = 0.8)

plot(SS1_ZINBSBM_N75_K3_Fixed_K5_T40000_1_LSZ_RI,type = "l",xlab = "",ylab = "", main = "",cex.axis = 0.8, ylim = c(0.35,1),col = 2,lty = 1)
lines(SS1_ZIPSBM_N75_K3_Fixed_K5_T40000_1_LSZ_RI,type = "l", col = 1,lty = 2)
title(xlab = "Iteration",ylab = "RI", main = "K=5 Cases Rand Index", mgp=c(0.9,0.1,0),cex.main=0.8,cex.lab = 0.8)
legend("bottomright", legend=c("ZINB-SBM","ZIP-SBM"),
       col=1:2, lty = 2:1, cex=0.6)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42),mgp=c(3,1,0))
```

Back to the SS1 ZINB-SBM $K=3$ case, we can also check some specific clustering states of the posterior chain:

``` r
# Check some specific clustering states
table(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[1]]%*%c(1:3),SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)) # initial state
table(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[2]]%*%c(1:3),SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)) # first state
table(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[3]]%*%c(1:3),SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)) # second state
table(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[40001]]%*%c(1:3),SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)) # end state
```

Note here it's possible that the asymptotic behavior we discussed in the paper might happen, so that, the posterior clustering chain would end up with the clustering state which has less number of clusters than the fixed $K$ input to the algorithm.
Based on the reasons we provided in the paper, we propose here to rerun the code again and obtain another round of the outputs if this situation happens.
This is the multiple implementation we might need to apply as we discussed in the paper.

If the multiple implementations are not required, we can then obtain part of the summary statistics.
We can summarize missing zero probability $p$ by simply evaluating the posterior mean of the samples and we can also check the mixing of the $p$ posterior chain by the trace plot.

``` r
## Summarize p
plot(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$p[1:40001],type = "l") # trace plot of posterior samples of p from iteration 0 to 40000
plot(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001],type = "l") # trace plot of posterior samples of p from iteration 20000 to 40000
hist(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001]) # histogram plot of posterior samples of p from iteration 20000 to 40000
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp <-
  mean(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001]) # obtain the summarized p by the posterior mean after burn-in
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp
```

The acceptance rate of the $\boldsymbol{R}$ M-H step can be checked by:

``` r
## check acceptance rate
Reduce(`+`, SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Acceptance_count_R[20001:40001])/20001
```

Then we aim to summarize the posterior clustering samples by minimizing an unbaised estimator of the expected posterior VI loss as we discussed in the section $3.2$ of the paper.
We start from obtaining the marginal posterior mode which is used as the initial clustering state of the greedy algorithm for the minimization.

``` r
# Obtain the marginal posterior mode of the Z chain
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States <- list() # initialize a list which will store different clustering states; the clustering states are labeled from 1,2,3... and are put at the 1st,2nd,3rd... element of the list, respectively
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration <- list() # store all the t's (iteration number) which provides the same cluster as the clustering state 1,2,3...
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop <- 20001:40001 # the "IterationLoop" which stores all the iteration t's which we focus on, that is, all the iteration t's after burn-in
StatesLabelIndicator = 0 # initialize the label for the clustering states
while (length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop)!=0){ # if the "IterationLoop" is not empty
  StatesLabelIndicator <- StatesLabelIndicator + 1 # assign the next label to the next clustering state
  SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState <- SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop[1]]] # extract the first clustering state for the "IterationLoop"
  SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States[[StatesLabelIndicator]] <- SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState # store the first state within the "IterationLoop" with label "StatesLabelIndicator" in the list which will contain all different unique states
  SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp <- c() # create a vector to temporarily store all the iteration t's whose clustering is the same as the first clustering state within the "IterationLoop"
  for (t in SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop){ # loop over all the current existing iterations in "IterationLoop"
    if (sum(c(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]%*%1:ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]))==
            c(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState%*%1:ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState)))==nrow(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState)){ # if the t's clustering is the same as the "FirstState"
      SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp <- c(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp,t) # store the iteration t in the temporary vector
    }
  }
  SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration[[StatesLabelIndicator]] <- SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp # store all the t's as the list element
  SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop <- (20001:40001)[-(unlist(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration)-20000)] # remove all the iterations we have stored and then move to the next clustering state
}
length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States) # check the number of different clustering states
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency <- c() # check the number of times one clustering state occurs
for (t in 1:length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States)){
  SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency <- 
    c(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency,
      length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration[[t]]))
}
which.max(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency) # find the marginal posterior mode, i.e. the most frequent clustering state
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ <- SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States[[which.max(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency)]] # initialize the clustering state as the marginal posterior mode
table(c(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%1:ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ)),SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3),dnn = c("","")) # compare the marginal posterior mode with the true lcustering
```

Then we can input the marginal posterior mode to the algorithm provided by [Rastelli and Friel (2018)](https://doi.org/10.1007/s11222-017-9786-y) in the R package `GreedyEPL` and obtain the summarized clustering.

``` r
require(GreedyEPL) # obtain the summarized Z by the greedy algorithm proposed by Rastelli and Friel (2018)
Z_temp <- c()
for (t in 20001:40001){ # transform Z matrices to z vectors and put all vectors together as a T X N matrix if T is the number of iterations here
  Z_temp <- rbind(Z_temp,c(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]%*%1:ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]])))
}
output <- MinimiseEPL(Z_temp, list(Kup = 10, loss_type = "VI",
                                   decision_init = c(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%
                                                       1:ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ))))
table(output$decision,SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3),dnn = c("","")) # check whether the the output is the same as the marginal posterior mode
output$EPL # check VI posterior loss: 0
```

In most of the cases, the summarized clustering is the same as the marginal posterior mode in our experiments.
This makes sense because our ZINB-SBM is a high-dimensional complex model where the model parameters are the mixing of discrete and continuous cases, so that, the posterior clustering chain is usually stably concentrated around the end clustering state forming the stationary distribution of the posterior clustering.
However, in the case of that the summarized clustering is not the same as the marginal posterior mode, we transform the output clustering, which is in vector form of the clustering, from the function `MinimiseEPL` above to the matrix form we focus on in the experiments.

``` r
# # if the clustering obtained by MinimiseEPL() is not the same as the marginal posterior mode, we treat such a clustering as the summarized clustering
# SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ <- matrix(0,length(output$decision),max(output$decision)) # Save output as summarized Z
# for (i in 1:length(output$decision)){SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[i,output$decision[i]]<-1}
# # # if the clustering obtained by MinimiseEPL() has smaller K than the fixed K, we add the zero column for the corresponding empty clusters
# # SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ <- cbind(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ,0)
```







