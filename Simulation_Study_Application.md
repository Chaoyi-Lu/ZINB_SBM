# Simulation Study Applications
This markdown file illustrates the code and applications for the simulation studies shown in the *Zero-Inflated Negative-Binomial Stochastic Block Model* paper.
There are two simulation studies experimented in the paper.
The simulation study $1$ focuses on an artificial network randomly generated from the zero-inflated Negative-Binomial stochastic block model (ZINB-SBM) and the simulation study $2$ (SS2) focuses on an artificial network randomly generated from the zero-inflated Poissin stochastic block model (ZIP-SBM).
Both simulation studies fit both ZINB-SBM and ZIP-SBM to the artificial datasets.

The source function code for the implementations of the inference of both models is included in the file [`Functions_for_ZINB_SBM.R`].
The explanations of the code are written as comments beside the code.
Note here that the practitioners need to load the source functions everytime after cleaning the environment.

``` r
rm(list=ls()) # remove all in the environment
gc() # Free unused memory
source("Functions_for_ZINB_SBM.R") # load the functions
```

The packages required in the applications are also included in such a file.

## 1. Simulation Study $1$

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

The artificial network data we focus on in the paper is included in the files within this repository and can be loaded by the code below.
Note that the latent clustering variable is mainly in the matrix form, $\boldsymbol{Z}$, in the code. 
Recall from Section $2.1$ of the paper that the latent clustering variable in the form of $\boldsymbol{Z}$ is a $N \times K$ matrix whose $i$ th row is the vector $\boldsymbol{z_i}=(z_{i1},z_{i2},\dots,z_{iK})$ where $z_{ik}=1$ if node $i$ belongs to cluster $k$, and $z_{ik}=0$ otherwise.

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

The number or proportion of true zeros/missing zeros/non zeros can be checked by:

``` r
# Number/proportion of true 0
sum(SS1_ZINBSBM_N75_K3$nu[SS1_ZINBSBM_N75_K3$Y == 0]==0)-75
(sum(SS1_ZINBSBM_N75_K3$nu[SS1_ZINBSBM_N75_K3$Y == 0]==0)-75)/(75*74)
# Number/proportion  of missing 0
sum(SS1_ZINBSBM_N75_K3$nu[SS1_ZINBSBM_N75_K3$Y == 0]==1)
(sum(SS1_ZINBSBM_N75_K3$nu[SS1_ZINBSBM_N75_K3$Y == 0]==1))/(75*74)
# Number/proportion  of non 0
sum(SS1_ZINBSBM_N75_K3$nu[SS1_ZINBSBM_N75_K3$Y != 0]==0)
(sum(SS1_ZINBSBM_N75_K3$nu[SS1_ZINBSBM_N75_K3$Y != 0]==0))/(75*74)
```

Then we apply label switching on the latent clustering $\boldsymbol{Z}$ and those clustering dependent parameters, $\boldsymbol{\Pi},\boldsymbol{R},\boldsymbol{Q}$, of the simulated network. 
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

which argee with the ones we illustrate in the paper.

Once we obtained the label-switched latent clustering and parameters, we can evaluate the initial mean and variance of the Negative-Binomial distribution proposed for the non-missing weights between each pair of clusters.

``` r
#  Evaluate the initial mean and variance
SS1_ZINBSBM_N75_K3_obs_Initialmean <- SS1_ZINBSBM_N75_K3_obs_InitialR*(1-SS1_ZINBSBM_N75_K3_obs_InitialQ)/SS1_ZINBSBM_N75_K3_obs_InitialQ
SS1_ZINBSBM_N75_K3_obs_Initialvar <- SS1_ZINBSBM_N75_K3_obs_InitialR*(1-SS1_ZINBSBM_N75_K3_obs_InitialQ)/SS1_ZINBSBM_N75_K3_obs_InitialQ^2
```

The reference probability of missing zero, $p$, is simply the initial setting:

``` r
# Set initial p
SS1_ZINBSBM_N75_K3_obs_Initialp <- 0.15
```

We can also evaluate the initial $\boldsymbol{P_{m0}}$ (the probability of the zero interaction being missing zero) proposed for this network.

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

The plot of the adjacency matrix $Y$, the adjacency matrix $Y$ plotted based on the true clustering $\boldsymbol{Z}^{\*}$ and the indicator of whether each $y_{ij}$ is non-zero (dark color) or not (light color) plotted according to $\boldsymbol{Z}^{\*}$ shown as Figure $1$ of the paper can be recovered by:

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

### 1.1 SS1 ZINB-SBM Implementations

The implementations of applying partially collapsed Metropolis within Gibbs algorithm (PCMwG) for the ZINB-SBM on the network are based on the function `Directed_ZINBSBM_PCMwG()` included in the source code file [`Functions_for_ZINB_SBM.R`].
Note that such a function will automatically label switch the initial clustering input to the function in order to ensure the uniqueness of the input clustering.
Other steps are exactly the same as the Algorithm $2$ stated in the paper.

The function `Directed_ZINBSBM_PCMwG_FixedZ()` aims to apply further inference conditional on the fixed summarized clustering as we disucssed in the paper.
Recall here that, within such a function, the inference step of the clustering is removed and instead the clustering is fixed to be the summarized clustering we obtained from the outputs of the above function `Directed_ZINBSBM_PCMwG()`.

The PCMwG algorithm for the ZINB-SBM is implemented for $40,000$ iterations for each fixed $K = 2,3,4,5$.
The $p$ prior setting for the function `Directed_ZINBSBM_PCMwG()` is $p \sim \text{Beta}(1,9)$ which can be changed by inputting prior parameters.
The prior settings of other parameters are set by default as we discussed in the paper, that is, $\boldsymbol{\Pi} \sim \text{Dirichlet}(1, \dots, 1)$, $q_{gh} \sim \text{Beta}(1, 1)$ for $g,h=1,2,\dots,K$, and the prior distribution of $\boldsymbol{R}$ is simply positive uniform $\text{U}(0,\text{UpperBound})$ where the "UpperBound" here can be a big enough value so that the $\boldsymbol{R}$ prior term can be cancelled by the fraction in the acceptance ratio of the Metropolis-Hastings (M-H) step. Recall also here that the proprosal distrbution of $\boldsymbol{R}$ in the M-H step is $r'\_{gh} \sim \text{U}(\text{max}(0,r_{gh}^{(t-1)}-\epsilon),r_{gh}^{(t-1)}+\epsilon)$ for each pair of $g,h = 1,\dots,K$ where $r_{gh}^{(t-1)}$ is the current state of the $r_{gh}$ and the proposal epsilon $\epsilon$ here is tuned to be $0.175$ where such an epsilon will also be applied in the real data application.
The implementations are applied by the code shown below.
Note here that we also provide a reference implementation time for each case in the code below and it's shown to take significant time to finish.
So practitioners might prefer reducing the number of iterations for quicker completion because the results illustrated in the paper show that the posterior chains usually converge quickly within thousands of iterations for the simulation studies (note also that this may or may not be the case for real data applications).

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
We can also obtain the summarized $\tilde{\boldsymbol{Z}}$ and $\tilde{\boldsymbol{P_{m0}}}$ to evaluate the model selection criterion, integrated classification log-partially-collapsed-likelihood (ICPCL), by maximizing ICPCL with respect to $\boldsymbol{R}$ as we discussed in Section $3.2$ of the paper.
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

After the label switching, we can check the mixing of the clustering posterior samples by evaluating the rand index between each iteration's $\boldsymbol{Z}^{(t)}$ and the true clustering $\boldsymbol{Z}^*$.

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

Once we obtained the summarized clustering $\tilde{\boldsymbol{Z}}$, we can then obtain the summarized $\widetilde{\boldsymbol{P_{m0}}}$.
We first summarize the posterior $\boldsymbol{\nu}$ chain by applying the posterior mean.

``` r
## Summarize nu
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu <- matrix(0,nrow(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$nu[[1]]),ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$nu[[1]]))
for (t in 20001:40001){
  SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu <- SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu +
    SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$nu[[t]]
}
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu <- SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu/20001
```

Then, based on the summarized $\tilde{\boldsymbol{\nu}}$ and the clustering, we can obtain the summarized $\widetilde{\boldsymbol{P_{m0}}}$ by:

``` r
## Summarize the P_m0
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0 <- matrix(0,3,3)
for (k1 in 1:3){
  for (k2 in 1:3){
    SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0[k1,k2] <- 
      sum(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu[
        SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1,
        SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2])/
      (length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu[
        SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1,
        SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2][
          SS1_ZINBSBM_N75_K3$Y[SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1,
                               SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2]==0])-
         sum((SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1)*
               (SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2)))
  }
}
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0
```

Once we obtained the summarized $\widetilde{\boldsymbol{P_{m0}}}$ and $\tilde{\boldsymbol{Z}}$, we can approximate the ICPCL by the function `Directed_ZINBSBM_ApproximateICPCL_OptR` provided in the file [`Functions_for_ZINB_SBM.R`].
Such a function applies a greedy algorithm to search for an optimized $\boldsymbol{R}$ which maximize each partially collapsed $\boldsymbol{Y}_{gh}$ term of the ICPCL criterion as we discussed in Section $3.2$ of the paper.
The prior settings for the ICPCL are the same as those applied for the inference and the initial state of $\boldsymbol{R}$ for the greedy search is the posterior mean of $\boldsymbol{R}$ which can be obtained by:

``` r
# ## Obtain the initial R to be optimized for the ICPCL by the posterior mean of R
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR <-
  apply(array(unlist(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$R),
              dim = c(nrow(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$R[[1]]),
                      ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$R[[1]]),
                      length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$R)))[,,20001:40001],1:2,mean)
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR 
```

Then we can input all we need to the function and evaluate the ICPCL:

``` r
# Obtain the approximate ICPCL
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_ICPCLandOptR <- 
  Directed_ZINBSBM_ApproximateICPCL_OptR(Y = SS1_ZINBSBM_N75_K3$Y,
                                         ProbObs0Missing0 = SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0,
                                         Z = SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ,
                                         R_0 = SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR,
                                         alpha=1, beta1=1,beta2=9, betaq1=1,betaq2=1)
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_ICPCLandOptR$value # -6449.478
```

The concaveness of the log $\boldsymbol{Y}_{gh}$ term of the ICPCL with respect to the $r_{gh}$ can be checked by the function `Directed_AZINBSBMRghQghp_ApproximatedICL_CheckConcave()` provided in [`Functions_for_ZINB_SBM.R`].
The concaveness plots Figure $3$ shown in Section $4.1$ of the paper can be recovered by:

``` r
# Check the concaveness for each Y_gh term of the ICPCL w.r.t. r_gh
res <- Directed_ZINBSBM_ApproximateICPCL_logY_gh_term_CheckConcave(Y = SS1_ZINBSBM_N75_K3$Y,
                                                             ProbObs0Missing0 = SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0,
                                                             Z = SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ,
                                                             alpha=1, beta1=1,beta2=9, betaq1=1,betaq2=1,rUB = 6)
par(mfrow=c(3,3),mai = c(0.25, 0.15, 0.15, 0.05),mgp = c(1,0.25,0))
plot(seq(0.001,6,0.001),res[1,1,],type = "l", main = "g=1,h=1",xlab = TeX(r'($r_{11}$)'),ylab = "")
plot(seq(0.001,6,0.001),res[1,2,],type = "l", main = "g=1,h=2",xlab = TeX(r'($r_{12}$)'),ylab = "")
plot(seq(0.001,6,0.001),res[1,3,],type = "l", main = "g=1,h=3",xlab = TeX(r'($r_{13}$)'),ylab = "")
plot(seq(0.001,6,0.001),res[2,1,],type = "l", main = "g=2,h=1",xlab = TeX(r'($r_{21}$)'),ylab = "")
plot(seq(0.001,6,0.001),res[2,2,],type = "l", main = "g=2,h=2",xlab = TeX(r'($r_{22}$)'),ylab = "")
plot(seq(0.001,6,0.001),res[2,3,],type = "l", main = "g=2,h=3",xlab = TeX(r'($r_{23}$)'),ylab = "")
plot(seq(0.001,6,0.001),res[3,1,],type = "l", main = "g=3,h=1",xlab = TeX(r'($r_{31}$)'),ylab = "")
plot(seq(0.001,6,0.001),res[3,2,],type = "l", main = "g=3,h=2",xlab = TeX(r'($r_{32}$)'),ylab = "")
plot(seq(0.001,6,0.001),res[3,3,],type = "l", main = "g=3,h=3",xlab = TeX(r'($r_{33}$)'),ylab = "")
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42))
```

We apply the same process for all the fixed $K = 2,3,4,5$ cases and obtain the ICPCL table shown as Table $1$ in Section $4.1$ of the paper.
Since the code for other $K$ cases are similar as above, we propose not to provide more details here.

Once we picked the best $K$ case, we can apply further inference conditional on the summarized clustering $\tilde{\boldsymbol{Z}}$ in order to summarize those clustering dependent parameters, that is, $\boldsymbol{\Pi},\boldsymbol{R},\boldsymbol{Q}$.
The further inference is also implemented for the same number of iterations as above, that is, $40,000$ iterations and all other settings are also the same.
Note here that the `Z_s` in the code below denotes/means the summarized $\tilde{\boldsymbol{Z}}$.

``` r
# Further inference of R,Q,Pi conditional on summarized z
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s <-
  Directed_ZINBSBM_PCMwG_FixedZ(Y = SS1_ZINBSBM_N75_K3$Y,
                                K = 3, T = 40000, eps_R = 0.175,beta1 = 1, beta2 = 9,
                                Z_0 = SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ)
```

According to the further inference outputs (after $20,000$-iteration burn-in), we can summarize the clustering dependent parameters as:

``` r
## Summarize Pi
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedPi <-
  apply(array(unlist(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi),
              dim = c(nrow(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi[[1]]),
                      ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi[[1]]),
                      length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi)))[,,20001:40001],1,mean)
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedPi

## Summarize Q
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ <-
  apply(array(unlist(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Q),
              dim = c(nrow(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Q[[1]]),
                      ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Q[[1]]),
                      length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Q)))[,,20001:40001],1:2,mean)
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ
```

And also summarize $\boldsymbol{R}$ by:

```r
# ## Summarize R
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedR <-
  apply(array(unlist(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$R),
              dim = c(nrow(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$R[[1]]),
                      ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$R[[1]]),
                      length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$R)))[,,20001:40001],1:2,mean)
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedR
# SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR # Compare with posterior mean R
# SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_ICPCLandOptR$R # Compare with the optimized R used for ICPCL
# # acceptance rate for R in the further inference conditional on Z_s
# Reduce(`+`, SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Acceptance_count_R[20001:40001])/20001
```

Note here that we can also compare the summarized $\tilde{\boldsymbol{R}}$ with the posterior mean of $\boldsymbol{R}$ and with the $\boldsymbol{R}$ which optimize the ICPCL as the commentted code above.
It can be checked that all the three $\boldsymbol{R}$'s we obtained are all most the same as each other.
This result can be expected because the posterior clustering chain was shown to be very stable for this experiment.
We can also check the acceptance rate of the $\boldsymbol{R}$ M-H step for the further inference, and the acceptance rate is also shown to be similar as the one we obtained in the PCMwG implementation without fixing the clustering.

Based on the summarized $\tilde{\boldsymbol{R}}$ and $\tilde{\boldsymbol{Q}}$, we can also obtain the summarized mean and variance of the distribution assumed for the edge weights between each pair of summarized clusters:

```r
# Summarized mean 
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedR*
  (1-SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ)/
  SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ
# Summarized var
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedR*
  (1-SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ)/
  SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ^2
```

which can be checked that they agree well with the reference ones provided at the beginning of this section.

The posterior density plot of $p$ as well as the posterior density plot of $\boldsymbol{\Pi},\boldsymbol{R},\boldsymbol{Q}$ conditional on $\tilde{\boldsymbol{Z}}$ shown as Figure $4$ of the paper can be recovered via the code below.

```r
# Transform the list to the array
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredR <-
  array(unlist(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$R),
        dim = c(nrow(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$R[[1]]),
                ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$R[[1]]),
                length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$R)))[,,20001:40001]
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredQ <-
  array(unlist(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Q),
        dim = c(nrow(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Q[[1]]),
                ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Q[[1]]),
                length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Q)))[,,20001:40001]
# Make the 2 X 2 figure
par(mfrow=c(2,2),mai = c(0.3, 0.25, 0.2, 0.05), mgp=c(1.25,0.5,0))
# Histogram of the posterior p samples
hist(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001],ylab = "",xlab="", main = TeX(r'(Posterior Density of $p$)'))
abline(v=0.15,col = 2,lty=2)
par(xpd=TRUE)
text(0.14,1200, TeX(r'($p^*$)'), pos = 4,col=2)
par(xpd=FALSE)
# Posterior density plots of Pi conditional on Z_s
plot(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredPi[1,],bw=0.01),col = 2, xlim=c(0,0.7),ylim=c(0,10), ylab = "",xlab="", main = TeX(r'(Posterior Density of $\pi_k$|\widetilde{\textbf{z}})'))
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredPi[2,],bw=0.01),col = 3)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredPi[3,],bw=0.01),col = 4)
abline(v=0.35,col = 2,ylim=c(0,9),lty=2)
abline(v=0.2,col = 3,ylim=c(0,9),lty=2)
abline(v=0.45,col = 4,ylim=c(0,9),lty=2)
legend("topright", legend=c(TeX(r'($\pi_1$)'),TeX(r'($\pi_2$)'),TeX(r'($\pi_3$)')),
       col=2:4, lty = 1, cex=0.6)
# Posterior density plots of R conditional on Z_s
plot(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredR[1,1,]),col = 1, xlim=c(0,0.9),ylim=c(0,60), ylab = "",xlab="", main = TeX(r'(Posterior Density of $r_{gh}$|\widetilde{\textbf{z}})'))
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredR[1,2,],bw=0.005),col = 2)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredR[1,3,],bw=0.005),col = 3)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredR[2,1,],bw=0.005),col = 4)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredR[2,2,]),col = 5)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredR[2,3,],bw=0.005),col = 6)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredR[3,1,],bw=0.005),col = 7)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredR[3,2,],bw=0.005),col = 8)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredR[3,3,],bw=0.005),col = "rosybrown")
legend("topright", legend=c("1,1","1,2","1,3", "2,1","2,2","2,3", "3,1","3,2","3,3"),
       col=c(1:8,"rosybrown"), lty = 1, cex=0.6)
# Posterior density plots of Q conditional on Z_s
plot(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredQ[1,1,]),col = 1, xlim=c(0,0.5),ylim=c(0,35), ylab = "",xlab="", main = TeX(r'(Posterior Density of $q_{gh}$|\widetilde{\textbf{z}})'))
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredQ[1,2,]),col = 2)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredQ[1,3,]),col = 3)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredQ[2,1,]),col = 4)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredQ[2,2,]),col = 5)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredQ[2,3,]),col = 6)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredQ[3,1,]),col = 7)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredQ[3,2,]),col = 8)
lines(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredQ[3,3,]),col = "rosybrown")
legend("topright", legend=c("1,1","1,2","1,3", "2,1","2,2","2,3", "3,1","3,2","3,3"),
       col=c(1:8,"rosybrown"), lty = 1, cex=0.6)
# Add the extra posterior density plot of r_22 conditional on summarized z since it's scale is very different from others
par(new=TRUE,mfcol=c(5,5), mfg=c(4,2),mai=c(0.1,0.1,0.1,0.1),mgp=c(1,0.5,0)) # Add the r_22 posterior density additionally
plot(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredR[2,2,]),col = 5,main = "", ylab = "",xlab="")
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42),mgp=c(3,1,0))
```

Figure $5$ in the paper shows the summarized $\tilde{\boldsymbol{\nu}}$ plotted according to the summarized $\tilde{\boldsymbol{Z}}$ as well as the comparison of the bar plot of the posterior missing weight samples and the Negative-Binomial distribution with summarized and reference parameters for the interaction from node $2$ to node $1$.
The figure can be recovered by the code provided below.

```r
par(mfrow=c(1,2),mai = c(0.3, 0.15, 0.2, 0.15), mgp=c(1.25,0.5,0))
layout(matrix(c(1,2,2), nrow = 1, ncol = 3, byrow = TRUE))
## Summarized nu plot based on summarized z
image(t(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu)[order(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)),
                                                                     rev(order(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)))],axes = FALSE,main = TeX(r'(Summarized $\widetilde{\nu}$)',bold=TRUE))
group_counts <- (as.numeric(table(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3))))
abline(v = -1/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1))))
abline(h = 1-(-1/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS1_ZINBSBM_N75_K3$Y)-1)))))
## Compare posterior missing weights with NB distribution
# Take entry 2,1 as an example to check X
res_list <- c()
for (t in 20001:40001){
  if(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$nu[[t]][2,1]==1){ # if the y_21=0 is inferred as a missing zero
    res_list <- c(res_list,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1$X[[t]][2,1]) # store the sampled missing weight
  }
}
plot(table(res_list), main="Entry 2,1 Missing Weights Posterior Density Comparison with NB Distribution",xlab = "Weight",ylab = "frequency or density")
lines(0:23,table(res_list)[1]*dnbinom(0:23,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedR[SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[2,]%*%1:3,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[1,]%*%1:3],
                                      SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ[SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[2,]%*%1:3,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[1,]%*%1:3])/
        dnbinom(0,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedR[SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[2,]%*%1:3,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[1,]%*%1:3],
                SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ[SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[2,]%*%1:3,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[1,]%*%1:3]),col = 2,lty=2)
lines(0:23,table(res_list)[1]*dnbinom(0:23,SS1_ZINBSBM_N75_K3_obs_InitialR[SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[2,]%*%1:3,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[1,]%*%1:3],
                                      SS1_ZINBSBM_N75_K3_obs_InitialQ[SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[2,]%*%1:3,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[1,]%*%1:3])/
        dnbinom(0,SS1_ZINBSBM_N75_K3_obs_InitialR[SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[2,]%*%1:3,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[1,]%*%1:3],
                SS1_ZINBSBM_N75_K3_obs_InitialQ[SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[2,]%*%1:3,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ[1,]%*%1:3]),col = 3)
legend("topright", legend=c("Bar plot of missing weights", "NB distribution with summarized r,q", "NB distribution with reference r,q"),
       col=1:3, lty=c(1,2,1), cex=1)
par(mfrow=c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
```

Up to this point, we have illustrated all the outputs we show on the paper.
Here we also provide some extra summary statistics of interest which are not included in the paper but the performance could be expected.

Apart from the summarized parameters, we can also check the posterior mean of $\boldsymbol{\Pi}$ and $\boldsymbol{Q}$, and compare them with the summarized ones.

```r
## Posterior mean Pi
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanPi <-
  apply(array(unlist(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Pi),
              dim = c(nrow(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Pi[[1]]),
                      ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Pi[[1]]),
                      length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Pi)))[,,20001:40001],1,mean)
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanPi
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedPi # Compare with summarized Pi
## Posterior mean Q
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ <-
  apply(array(unlist(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Q),
              dim = c(nrow(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Q[[1]]),
                      ncol(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Q[[1]]),
                      length(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Q)))[,,20001:40001],1:2,mean)
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ # Compare with summarized Q
```

It can be expected and checked that the posterior mean of them are almost the same as the summarized ones due to the fact that the posterior clustering chain is shown to be very stable.
This also leads to the evaluation of the mean and variance of the edge weights based on the posterior mean of $\boldsymbol{\Pi}$ and $\boldsymbol{Q}$:

```r
# Check distribution mean based on PosteriorMeanR and PosteriorMeanQ
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR*
  (1-SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ)/
  SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ
# Check distribution var based on PosteriorMeanR and PosteriorMeanQ
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR*
  (1-SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ)/
  SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ^2
```

which are also agree well with the summarized ones and reference ones.

Recall here that the summarized $\tilde{p}$ is obatined directly by posterior mean, but here we also have the posterior samples of $p$ conditional on the summarized $\tilde{\boldsymbol{Z}}$ from the further inference.
Thus we can also compare the summarized $\tilde{p}$ with the posterior mean of $p|\tilde{\boldsymbol{Z}}$:

```r
# p Posterior mean conditional on summarized z
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s <-
  mean(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$p[20001:40001])
hist(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$p[20001:40001])
plot(density(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$p[20001:40001]))
plot(SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$p[20001:40001], type = "l")
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp # Compare with summarized p
```

which are also very close to each other.

Similar well agreement can also be discovered for the statistic $\boldsymbol{P_{m0}}$ evaluated by the posterior mean or the posterior mean conditional on the summarized clustering as well as the summarized $\tilde{\boldsymbol{P_{m0}}}$:

```r
# Compare P_m0 evaluated by posterior mean
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp/
  (SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp+
     (1-SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp)*
     dnbinom(0,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR,
             SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ))
# with P_m0 evaluated by posterior mean conditional on Z_s
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s/
  (SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s+
     (1-SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s)*
     dnbinom(0,SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedR,
             SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ))
# and with summarized P_m0
SS1_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0
```

All the well agreements shown above bring stronger confidence of the applications of the further inference conditional on the summarized clustering.

### 1.2 SS1 ZIP-SBM Implementations

The implementations of the Gibbs sampler for ZIP-SBM is similar to those of the ZINB-SBM.
Here we leverage the function `Directed_ZIPSBM_Gibbs()` provided in the source code file for the implementations and the function `Directed_ZIPSBM_Gibbs_FixedZ()` is the corresponding Gibbs sampler with the clustering fixed.
The prior distribution for the parameter $\lambda$ of the Poisson embeding is set as the gamma distribution $\lambda \sim \text{Ga}(1,1)$ and all other settings are the same as those of the ZINB-SBM cases.
The implementations can be applied by the code:

```r
# Simulation study 1: apply the ZIP-SBM Gibbs sampler for fixed K = 2, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS1_ZIPSBM_N75_K3_Fixed_K2_T40000_1 <- 
  Directed_ZIPSBM_Gibbs(Y = SS1_ZINBSBM_N75_K3$Y, K = 2, T = 40000, 
                        alpha=1, beta1 = 1,beta2 = 9, alpha1=1,alpha2=1, Z_0=NA)
end.time <- Sys.time()
SS1_ZIPSBM_N75_K3_Fixed_K2_T40000_1_time <- end.time - start.time
SS1_ZIPSBM_N75_K3_Fixed_K2_T40000_1_time # Time difference of 1.176335 hours
# save.image("SS1_ZIPSBM_N75_K3_Fixed_K2_T40000_1_prior_p_Beta_1_9.RData")
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Simulation study 1: apply the ZIP-SBM Gibbs sampler for fixed K = 3, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1 <- 
  Directed_ZIPSBM_Gibbs(Y = SS1_ZINBSBM_N75_K3$Y, K = 3, T = 40000, 
                        alpha=1, beta1 = 1,beta2 = 9, alpha1=1,alpha2=1, Z_0=NA)
end.time <- Sys.time()
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_time <- end.time - start.time
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_time # Time difference of 1.564769 hours
# save.image("SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_prior_p_Beta_1_9.RData")
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Simulation study 1: apply the ZIP-SBM Gibbs sampler for fixed K = 4, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS1_ZIPSBM_N75_K3_Fixed_K4_T40000_1 <- 
  Directed_ZIPSBM_Gibbs(Y = SS1_ZINBSBM_N75_K3$Y, K = 4, T = 40000, 
                        alpha=1, beta1 = 1,beta2 = 9, alpha1=1,alpha2=1, Z_0=NA)
end.time <- Sys.time()
SS1_ZIPSBM_N75_K3_Fixed_K4_T40000_1_time <- end.time - start.time
SS1_ZIPSBM_N75_K3_Fixed_K4_T40000_1_time # Time difference of 1.991501 hours
# save.image("SS1_ZIPSBM_N75_K3_Fixed_K4_T40000_1_prior_p_Beta_1_9.RData")
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Simulation study 1: apply the ZIP-SBM Gibbs sampler for fixed K = 5, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS1_ZIPSBM_N75_K3_Fixed_K5_T40000_1 <- 
  Directed_ZIPSBM_Gibbs(Y = SS1_ZINBSBM_N75_K3$Y, K = 5, T = 40000, 
                        alpha=1, beta1 = 1,beta2 = 9, alpha1=1,alpha2=1, Z_0=NA)
end.time <- Sys.time()
SS1_ZIPSBM_N75_K3_Fixed_K5_T40000_1_time <- end.time - start.time
SS1_ZIPSBM_N75_K3_Fixed_K5_T40000_1_time # Time difference of 2.428095 hours
# save.image("SS1_ZIPSBM_N75_K3_Fixed_K5_T40000_1_prior_p_Beta_1_9.RData")
```

Similar as ZINB-SBM cases, we also take $K=3$ case as an example below for the summarizing process.
Note here that the function `LabelSwitching_SG2003_ZIPSBM()` is the label-switching function coded specifically for the ZIP-SBM.

``` r
# # Simulation study 1: apply the ZIP-SBM Gibbs sampler for fixed K = 3, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
# rm(list=ls())
# gc()
# source("Functions_for_ZINB_SBM.R")
# start.time <- Sys.time()
# SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1 <- 
#   Directed_ZIPSBM_Gibbs(Y = SS1_ZINBSBM_N75_K3$Y, K = 3, T = 40000, 
#                         alpha=1, beta1 = 1,beta2 = 9, alpha1=1,alpha2=1, Z_0=NA)
# end.time <- Sys.time()
# SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_time <- end.time - start.time
# SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_time # Time difference of 1.564769 hours
# # save.image("SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_prior_p_Beta_1_9.RData")

# # Summarize the outputs
# Apply label switching
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS <-
  LabelSwitching_SG2003_ZIPSBM(Z = SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1$Z,
                               Pi = SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1$Pi,
                               Lambda = SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1$Lambda)
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1$Z <- c()
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1$Pi <- c()
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1$Lambda <- c()
gc()

## check rand index for each iteration
require("fossil")
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI <- c()
for (t in 1:40001){
  SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI <-
    c(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI,
      rand.index(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]%*%c(1:ncol(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]])),
                 SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)))
}
plot(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI,type = "l",xlab = "",ylab = "", main = "Rand Index",cex.axis = 0.8)

table(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[1]]%*%c(1:3),SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)) # initial state
table(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[2]]%*%c(1:3),SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)) # first state
table(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[3]]%*%c(1:3),SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)) # second state
table(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[40001]]%*%c(1:3),SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3)) # end state

## Summarize p
plot(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1$p[1:40001],type = "l")
hist(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001])
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp <-
  mean(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001])

# Obtain the marginal posterior mode of the Z chain
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States <- list()
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration <- list()
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop <- 20001:40001
StatesLabelIndicator = 0
while (length(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop)!=0){
  StatesLabelIndicator <- StatesLabelIndicator + 1
  SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState <- SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop[1]]]
  SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States[[StatesLabelIndicator]] <- SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState
  SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp <- c()
  for (t in SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop){
    if (sum(c(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]%*%1:ncol(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]))==
            c(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState%*%1:ncol(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState)))==nrow(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState)){
      SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp <- c(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp,t)
    }
  }
  SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration[[StatesLabelIndicator]] <- SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp
  SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop <- (20001:40001)[-(unlist(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration)-20000)]
}
length(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States)
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency <- c()
for (t in 1:length(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States)){
  SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency <- c(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency,
                                                               length(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration[[t]]))
}
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency
which.max(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency)
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ <- SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States[[which.max(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency)]]
table(c(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%1:ncol(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ)),SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3),dnn = c("",""))
require(GreedyEPL) # obtain the summarized Z by the greedy algorithm proposed by Rastelli and Friel (2018)
Z_temp <- c()
for (t in 20001:40001){
  Z_temp <- rbind(Z_temp,c(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]%*%1:ncol(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]])))
}
output <- MinimiseEPL(Z_temp, list(Kup = 10, loss_type = "VI",
                                   decision_init = c(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%
                                                       1:ncol(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ))))
table(output$decision,SS1_ZINBSBM_N75_K3_LSZ%*%c(1:3),dnn = c("",""))
output$EPL # check VI posterior loss: -1.110223e-15

## Summarize nu
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu <- matrix(0,nrow(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1$nu[[1]]),ncol(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1$nu[[1]]))
for (t in 20001:40001){
  SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu <- SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu +
    SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1$nu[[t]]
}
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu <- SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu/20001

## Summarize the P_m0
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0 <- matrix(0,3,3)
for (k1 in 1:3){
  for (k2 in 1:3){
    SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0[k1,k2] <- sum(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu[SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1,SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2])/
      (length(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu[SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1,SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2][SS1_ZINBSBM_N75_K3$Y[SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1,SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2]==0])-sum((SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1)*(SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2)))
  }
}
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0

# Check ExactICL with summarized Z and nu
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_ICPCLandOptR <- 
  Directed_ZIPSBM_ExactICL(Y = UKfaculty_directed_weighted_sparse_network_adj, 
                                       ProbObs0Missing0 = SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0,
                                       Z = SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ, 
                                       alpha1=1,alpha2=1, beta1=1,beta2=9, alpha=1)
SS1_ZIPSBM_N75_K3_Fixed_K3_T40000_1_ICPCLandOptR$value # -7678.607
```

Note here that the function `Directed_ZIPSBM_ExactICL()` above is applied for the evaluation of the exact integrated complete data log likelihood ($\text{ICL}_{\text{ex}}$).
Due to the fact that none of the ZIP-SBM cases fit well to the simulation study $1$ artificial network, that is, none of the cases provide the summarized clustering which agrees with the true clustering and the summarized $\tilde{p},\widetilde{\boldsymbol{P\_{m0}}}$ are far away from the true references, we propose not to apply further implementations for these cases.

## 2. Simulation Study $2$

The second simmulation study instead focuses on an artificial network randomly generated from the ZIP-SBM with the settings: $N = 75, K = 3, p = 0.15$ and 

$$\boldsymbol{\Pi} = \left(0.3, 0.4, 0.3\right), \boldsymbol{\lambda} = \begin{pmatrix} 2.0 & 0.7 & 0.9\\\0.1 & 2.0 & 0.3\\\0.5 & 1.1 &2.5\end{pmatrix}$$

```r
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
# Simulation study 2 of AZIPSBMp N = 75, K = 3
SS2_ZIPSBM_N75_K3 <- 
  Simulation_Directed_ZIPSBM(N = 75 ,K = 3 ,Pi = c(0.3,0.4,0.3), p = 0.15,
                             Lambda = matrix(c(2,0.1,0.5,
                                               0.7,2.0,1.1,
                                               0.9,0.3,2.5),3,3))
# write.csv(SS2_ZIPSBM_N75_K3$Y,"SS2_ZIPSBM_N75K3_obsY.csv", row.names = FALSE)
# write.csv(SS2_ZIPSBM_N75_K3$X,"SS2_ZIPSBM_N75K3_obsX.csv", row.names = FALSE)
# write.csv(SS2_ZIPSBM_N75_K3$nu,"SS2_ZIPSBM_N75K3_obsnu.csv", row.names = FALSE)
# write.csv(SS2_ZIPSBM_N75_K3$Z,"SS2_ZIPSBM_N75K3_obsZ.csv", row.names = FALSE)
```

where the function `Simulation_Directed_ZIPSBM()` simulates one sample from the ZIP-SBM.
The network data we used in the paper can be loaded by the code:

```r
SS2_ZIPSBM_N75_K3 <- 
  list(Y = as.matrix(read.csv("SS2_ZIPSBM_N75K3_obsY.csv",header = TRUE)),
       X = as.matrix(read.csv("SS2_ZIPSBM_N75K3_obsX.csv",header = TRUE)),
       nu = as.matrix(read.csv("SS2_ZIPSBM_N75K3_obsnu.csv",header = TRUE)),
       Z = as.matrix(read.csv("SS2_ZIPSBM_N75K3_obsZ.csv",header = TRUE)))
colnames(SS2_ZIPSBM_N75_K3$Y) <- c()
colnames(SS2_ZIPSBM_N75_K3$X) <- c()
colnames(SS2_ZIPSBM_N75_K3$nu) <- c()
colnames(SS2_ZIPSBM_N75_K3$Z) <- c()
```

And we can check the number or proportion of difference types of zeros via:

```r
# Number/proportion of true 0
sum(SS2_ZIPSBM_N75_K3$nu[SS2_ZIPSBM_N75_K3$Y == 0]==0)-75
(sum(SS2_ZIPSBM_N75_K3$nu[SS2_ZIPSBM_N75_K3$Y == 0]==0)-75)/(75*74)

# Number/proportion of missing 0
sum(SS2_ZIPSBM_N75_K3$nu[SS2_ZIPSBM_N75_K3$Y == 0]==1)
(sum(SS2_ZIPSBM_N75_K3$nu[SS2_ZIPSBM_N75_K3$Y == 0]==1))/(75*74)

# Number/proportion of non 0
sum(SS2_ZIPSBM_N75_K3$nu[SS2_ZIPSBM_N75_K3$Y != 0]==0)
(sum(SS2_ZIPSBM_N75_K3$nu[SS2_ZIPSBM_N75_K3$Y != 0]==0))/(75*74)
```

The label switching process follows:

```r
# Label switch Z
SS2_ZIPSBM_N75_K3_LSZ <- LabelSwitching_SG2003_ZIPSBM(Z = list(SS2_ZIPSBM_N75_K3$Z))$Z[[1]]
image(t(SS2_ZIPSBM_N75_K3$Y)[order(SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3)),rev(order(SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3)))])
group_counts <- (as.numeric(table(SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3))))
abline(v = -1/(2*(nrow(SS2_ZIPSBM_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS2_ZIPSBM_N75_K3$Y)-1))))
abline(h = 1-(-1/(2*(nrow(SS2_ZIPSBM_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS2_ZIPSBM_N75_K3$Y)-1)))))
#--------------------------------------------------------------------------------------------------------------------------------------------
# Label switch initial Lambda, Pi
res <- LabelSwitching_SG2003_ZIPSBM(Z = list(SS2_ZIPSBM_N75_K3$Z),
                                    Pi = list(c(0.3,0.4,0.3)),
                                    Lambda = list(matrix(c(2,0.1,0.5,
                                                           0.7,2.0,1.1,
                                                           0.9,0.3,2.5),3,3)))
SS2_ZIPSBM_N75_K3_obs_InitialLambda <- res$Lambda[[1]]
SS2_ZIPSBM_N75_K3_obs_InitialPi <- res$Pi[[1]]
```

which brings the label switched initial $\boldsymbol{\Pi}$ and $\boldsymbol{Z}$:

$$\boldsymbol{\Pi} = \left(0.3, 0.3, 0.4\right), \boldsymbol{\lambda} = \begin{pmatrix} 2.5 & 0.5 & 1.1\\\0.9 & 2.0 & 0.7\\\0.3 & 0.1 &2.0\end{pmatrix}$$

as shown in Section $4.2$ of the paper.

The reference missing zero probability $p^*$ and the conditional missing zero probability $\boldsymbol{P_{m0}}^\*$ for SS2 network can be set or evaluated as:

```r
# Set initial p
SS2_ZIPSBM_N75_K3_obs_Initialp <- 0.15

# Evaluate the initial P_m0
SS2_ZIPSBM_N75_K3_obs_InitialProbObs0Missing0 <- matrix(0,3,3)
for (k1 in 1:3){
  for (k2 in 1:3){
    SS2_ZIPSBM_N75_K3_obs_InitialProbObs0Missing0[k1,k2] <- SS2_ZIPSBM_N75_K3_obs_Initialp/(SS2_ZIPSBM_N75_K3_obs_Initialp + (1-SS2_ZIPSBM_N75_K3_obs_Initialp)*dpois(0,SS2_ZIPSBM_N75_K3_obs_InitialLambda[k1,k2]))
  }
}
SS2_ZIPSBM_N75_K3_obs_InitialProbObs0Missing0
```

### 2.1 SS2 ZINB-SBM Implementations

The implementations as well as all the settings of SS2 ZINB-SBM cases are similar as those of SS1 cases:

```r
# Simulation study 2: apply the ZINB-SBM Metropolis within Gibbs algorithm for fixed K = 2, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS2_ZINBSBM_N75_K3_Fixed_K2_T40000_1 <- 
  Directed_ZINBSBM_PCMwG(Y = SS2_ZIPSBM_N75_K3$Y, K = 2, T = 40000, eps_R = 0.175,beta1 = 1, beta2 = 9)
end.time <- Sys.time()
SS2_ZINBSBM_N75_K3_Fixed_K2_T40000_1_time <- end.time - start.time
SS2_ZINBSBM_N75_K3_Fixed_K2_T40000_1_time # Time difference of 1.71906 hours
# save.image("SS2_ZINBSBM_N75_K3_Fixed_K2_T40000_1_prior_p_Beta_1_9.RData")
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Simulation study 2: apply the ZINB-SBM Metropolis within Gibbs algorithm for fixed K = 4, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS2_ZINBSBM_N75_K3_Fixed_K4_T40000_1 <- 
  Directed_ZINBSBM_PCMwG(Y = SS2_ZIPSBM_N75_K3$Y, K = 4, T = 40000, eps_R = 0.175,beta1 = 1, beta2 = 9)
end.time <- Sys.time()
SS2_ZINBSBM_N75_K3_Fixed_K4_T40000_1_time <- end.time - start.time
SS2_ZINBSBM_N75_K3_Fixed_K4_T40000_1_time # Time difference of 3.013449 hours
# save.image("SS2_ZINBSBM_N75_K3_Fixed_K4_T40000_1_prior_p_Beta_1_9.RData")
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Simulation study 2: apply the ZINB-SBM Metropolis within Gibbs algorithm for fixed K = 5, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS2_ZINBSBM_N75_K3_Fixed_K5_T40000_1 <- 
  Directed_ZINBSBM_PCMwG(Y = SS2_ZIPSBM_N75_K3$Y, K = 5, T = 40000, eps_R = 0.175,beta1 = 1, beta2 = 9)
end.time <- Sys.time()
SS2_ZINBSBM_N75_K3_Fixed_K5_T40000_1_time <- end.time - start.time
SS2_ZINBSBM_N75_K3_Fixed_K5_T40000_1_time # Time difference of 3.713538 hours
# save.image("SS2_ZINBSBM_N75_K3_Fixed_K5_T40000_1_prior_p_Beta_1_9.RData")
```

And the $K=3$ case as well as the example summarizing processes are also similar, so we plan not to provide more explanations for the code below.
The code illustrated here just aims to make everything more convenient for the reader to recover the work.

``` r
# Simulation study 2: apply the ZINB-SBM Metropolis within Gibbs algorithm for fixed K = 3, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1 <- 
  Directed_ZINBSBM_PCMwG(Y = SS2_ZIPSBM_N75_K3$Y, K = 3, T = 40000, eps_R = 0.175,beta1 = 1, beta2 = 9)
end.time <- Sys.time()
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_time <- end.time - start.time
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_time # Time difference of 2.442501 hours
# save.image("SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_prior_p_Beta_1_9.RData")

# # Summarize the outputs
# Apply label switching
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS <-
  LabelSwitching_SG2003_ZINBSBM(Z = SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Z,
                                Pi = SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Pi,
                                R = SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$R,
                                Q = SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Q,
                                Acceptance_count_R = SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Acceptance_count_R)
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Z <- c()
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Pi <- c()
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$R <- c()
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Q <- c()
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$Acceptance_count_R <- c()
gc()

## check rand index for each iteration
require("fossil")
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI <- c()
for (t in 1:40001){
  SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI <-
    c(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI,
      rand.index(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]%*%c(1:ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]])),
                 SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3)))
}
plot(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI,type = "l",xlab = "",ylab = "", main = "Rand Index",cex.axis = 0.8)
# Check some specific clustering states
table(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[1]]%*%c(1:3),SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3)) # initial state
table(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[2]]%*%c(1:3),SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3)) # first state
table(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[3]]%*%c(1:3),SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3)) # second state
table(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[40001]]%*%c(1:3),SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3)) # end state

## Summarize p
plot(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$p[1:40001],type = "l")
hist(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001])
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp <-
  mean(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001])
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp

# Obtain the marginal posterior mode of the Z chain
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States <- list()
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration <- list()
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop <- 20001:40001
StatesLabelIndicator = 0
while (length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop)!=0){
  StatesLabelIndicator <- StatesLabelIndicator + 1
  SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState <- SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop[1]]]
  SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States[[StatesLabelIndicator]] <- SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState
  SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp <- c()
  for (t in SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop){
    if (sum(c(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]%*%1:ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]))==
            c(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState%*%1:ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState)))==nrow(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState)){
      SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp <- c(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp,t)
    }
  }
  SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration[[StatesLabelIndicator]] <- SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp
  SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop <- (20001:40001)[-(unlist(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration)-20000)]
}
length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States)
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency <- c()
for (t in 1:length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States)){
  SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency <- c(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency,
                                                                length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration[[t]]))
}
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency
which.max(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency)
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ <- SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States[[which.max(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency)]]
table(c(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%1:ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ)),SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3),dnn = c("",""))
require(GreedyEPL) # obtain the summarized Z by the greedy algorithm proposed by Rastelli and Friel (2018)
Z_temp <- c()
for (t in 20001:40001){
  Z_temp <- rbind(Z_temp,c(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]%*%1:ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]])))
}
output <- MinimiseEPL(Z_temp, list(Kup = 10, loss_type = "VI",
                                   decision_init = c(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%
                                                       1:ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ))))
table(output$decision,SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3),dnn = c("",""))
output$EPL # check VI posterior loss: 0

## Check R M-H step acceptance rate
Reduce(`+`, SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Acceptance_count_R[20001:40001])/20001

## Summarize nu
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu <- matrix(0,nrow(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$nu[[1]]),ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$nu[[1]]))
for (t in 20001:40001){
  SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu <- SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu +
    SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1$nu[[t]]
}
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu <- SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu/20001

## Summarize the P_m0
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0 <- matrix(0,3,3)
for (k1 in 1:3){
  for (k2 in 1:3){
    SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0[k1,k2] <- 
      sum(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu[SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1,
                                                            SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2])/
      (length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu[SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1,
                                                                SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2][
                                                                  SS2_ZIPSBM_N75_K3$Y[SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1,
                                                                                      SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2]==0])-
         sum((SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1)*(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2)))
  }
}
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0


# ## Obtain the initial R to be optimized for the ICPCL by the posterior mean of R
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR <-
  apply(array(unlist(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$R),
              dim = c(nrow(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$R[[1]]),
                      ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$R[[1]]),
                      length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$R)))[,,20001:40001],1:2,mean)

# Obtain the approximate ICPCL
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_ICPCLandOptR <- 
  Directed_ZINBSBM_ApproximateICPCL_OptR(Y = SS2_ZIPSBM_N75_K3$Y,
                                         ProbObs0Missing0 = SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0,
                                         Z = SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ,
                                         R_0 = SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR,
                                         alpha=1, beta1=1,beta2=9, betaq1=1,betaq2=1)
# -8098.871

# Further inference of R,Q,Pi conditional on summarized z
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s <-
  Directed_ZINBSBM_PCMwG_FixedZ(Y = SS2_ZIPSBM_N75_K3$Y,
                                             K = 3, T = 40000, eps_R = 0.175,beta1 = 1, beta2 = 9,
                                             Z_0 = SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ)
# Summarize Pi
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedPi <-
  apply(array(unlist(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi),
              dim = c(nrow(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi[[1]]),
                      ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi[[1]]),
                      length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi)))[,,20001:40001],1,mean)
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedPi

# Pi density plots
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredPi <- 
  array(unlist(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi),
        dim = c(nrow(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi[[1]]),
                ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi[[1]]),
                length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi)))[,,20001:40001]
plot(density(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredPi[1,]),col = 2, xlim=c(0.1,0.6),ylim=c(0,10), ylab = "",xlab="", main = TeX(r'($\pi_k$ Posterior Density)'))
lines(density(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredPi[2,]),col = 3)
lines(density(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredPi[3,]),col = 4)
legend("topright", legend=c(TeX(r'($\pi_1$)'),TeX(r'($\pi_2$)'),TeX(r'($\pi_3$)')),
       col=2:4, lty = 1, cex=0.7)

# ## Summarize R
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedR <-
  apply(array(unlist(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$R),
              dim = c(nrow(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$R[[1]]),
                      ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$R[[1]]),
                      length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$R)))[,,20001:40001],1:2,mean)
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedR
# SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR # Compare with posterior mean of R
# SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_ICPCLandOptR$R # Compare with the optimized R used for ICPCL
# # acceptance rate for R in the further inference conditional on Z_s
# Reduce(`+`, SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Acceptance_count_R[20001:40001])/20001

## Summarize Q
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ <-
  apply(array(unlist(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Q),
              dim = c(nrow(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Q[[1]]),
                      ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Q[[1]]),
                      length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Q)))[,,20001:40001],1:2,mean)
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ

# Summarized mean
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedR*
  (1-SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ)/
  SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ
# Summarized var
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedR*
  (1-SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ)/
  SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ^2
```

Similar extra summary statistics as we mentioned in SS1 can also be checked here:

``` r
## Posterior mean Pi
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanPi <-
  apply(array(unlist(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Pi),
              dim = c(nrow(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Pi[[1]]),
                      ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Pi[[1]]),
                      length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Pi)))[,,20001:40001],1,mean)
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanPi
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedPi # Compare with summarized posterior samples

## Posterior mean Q
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ <-
  apply(array(unlist(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Q),
              dim = c(nrow(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Q[[1]]),
                      ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Q[[1]]),
                      length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_LS$Q)))[,,20001:40001],1:2,mean)
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ # Compare with summarized posterior samples

# Check distribution mean based on PosteriorMeanR and PosteriorMeanQ
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR*
  (1-SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ)/
  SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ
# Check distribution var based on PosteriorMeanR and PosteriorMeanQ
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR*
  (1-SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ)/
  SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ^2

# p Posterior mean conditional on Z_s
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s <-
  mean(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$p[20001:40001])
hist(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$p[20001:40001])
plot(density(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$p[20001:40001]))
plot(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$p[20001:40001], type = "l")
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp # Compare with summarized posterior samples

# Compare P_m0 evaluated by posterior mean
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp/
  (SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp+
     (1-SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp)*
     dnbinom(0,SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanR,
             SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanQ))
# with P_m0 evaluated by posterior mean conditional on Z_s
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s/
  (SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s+
     (1-SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s)*
     dnbinom(0,SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedR,
             SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedQ))
# and with summarized P_m0
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0
```

### 2.2 SS2 ZIP-SBM Implementations

This subsection shows the code and applications of the Gibbs sampler for ZIP-SBM fit to the SS2 artificial network.

``` r
# Simulation study 2: apply the ZIP-SBM Gibbs sampler for fixed K = 2, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS2_ZIPSBM_N75_K3_Fixed_K2_T40000_1 <- 
  Directed_ZIPSBM_Gibbs(Y = SS2_ZIPSBM_N75_K3$Y, K = 2, T = 40000, 
                        alpha=1, beta1 = 1,beta2 = 9, alpha1=1,alpha2=1, Z_0=NA)
end.time <- Sys.time()
SS2_ZIPSBM_N75_K3_Fixed_K2_T40000_1_time <- end.time - start.time
SS2_ZIPSBM_N75_K3_Fixed_K2_T40000_1_time # Time difference of 1.023112 hours
# save.image("SS2_ZIPSBM_N75_K3_Fixed_K2_T40000_1_prior_p_Beta_1_9.RData")
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Simulation study 2: apply the ZIP-SBM Gibbs sampler for fixed K = 4, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS2_ZIPSBM_N75_K3_Fixed_K4_T40000_1 <- 
  Directed_ZIPSBM_Gibbs(Y = SS2_ZIPSBM_N75_K3$Y, K = 4, T = 40000, 
                        alpha=1, beta1 = 1,beta2 = 9, alpha1=1,alpha2=1, Z_0=NA)
end.time <- Sys.time()
SS2_ZIPSBM_N75_K3_Fixed_K4_T40000_1_time <- end.time - start.time
SS2_ZIPSBM_N75_K3_Fixed_K4_T40000_1_time # Time difference of 1.730288 hours
# save.image("SS2_ZIPSBM_N75_K3_Fixed_K4_T40000_1_prior_p_Beta_1_9.RData")
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Simulation study 2: apply the ZIP-SBM Gibbs sampler for fixed K = 5, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS2_ZIPSBM_N75_K3_Fixed_K5_T40000_1 <- 
  Directed_ZIPSBM_Gibbs(Y = SS2_ZIPSBM_N75_K3$Y, K = 5, T = 40000, 
                        alpha=1, beta1 = 1,beta2 = 9, alpha1=1,alpha2=1, Z_0=NA)
end.time <- Sys.time()
SS2_ZIPSBM_N75_K3_Fixed_K5_T40000_1_time <- end.time - start.time
SS2_ZIPSBM_N75_K3_Fixed_K5_T40000_1_time # Time difference of 2.157126 hours
# save.image("SS2_ZIPSBM_N75_K3_Fixed_K5_T40000_1_prior_p_Beta_1_9.RData")
```

The $K=3$ case with the summarizing process can be applied by the code:

``` r
# Simulation study 2: apply the ZIP-SBM Gibbs sampler for fixed K = 3, N = 75, T = 40000, Round 1 with p prior Beta(1,9)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1 <- 
  Directed_ZIPSBM_Gibbs(Y = SS2_ZIPSBM_N75_K3$Y, K = 3, T = 40000, 
                        alpha=1, beta1 = 1,beta2 = 9, alpha1=1,alpha2=1, Z_0=NA)
end.time <- Sys.time()
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_time <- end.time - start.time
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_time # Time difference of 1.437432 hours
# save.image("SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_prior_p_Beta_1_9.RData")

# # Summarize the outputs
# Apply label switching
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS <-
  LabelSwitching_SG2003_ZIPSBM(Z = SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$Z,
                               Pi = SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$Pi,
                               Lambda = SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$Lambda)
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$Z <- c()#in order to save memory
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$Pi <- c()
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$Lambda <- c()
gc()

## check rand index for each iteration
require("fossil")
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI <- c()
for (t in 1:40001){
  SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI <-
    c(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI,
      rand.index(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]%*%c(1:ncol(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]])),
                 SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3)))
}
plot(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_RI,type = "l",xlab = "",ylab = "", main = "Rand Index",cex.axis = 0.8)
# Check some specific clustering states
table(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[1]]%*%c(1:3),SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3)) # initial state
table(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[2]]%*%c(1:3),SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3)) # first state
table(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[3]]%*%c(1:3),SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3)) # second stat
table(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[40001]]%*%c(1:3),SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3)) # end state

## Summarize p
plot(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$p[1:40001],type = "l")
hist(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001])
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp <-
  mean(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001])
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp

# Obtain the marginal posterior mode of the Z chain
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States <- list()
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration <- list()
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop <- 20001:40001
StatesLabelIndicator = 0
while (length(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop)!=0){
  StatesLabelIndicator <- StatesLabelIndicator + 1
  SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState <- SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop[1]]]
  SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States[[StatesLabelIndicator]] <- SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState
  SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp <- c()
  for (t in SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop){
    if (sum(c(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]%*%1:ncol(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]))==
            c(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState%*%1:ncol(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState)))==nrow(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_FirstState)){
      SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp <- c(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp,t)
    }
  }
  SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration[[StatesLabelIndicator]] <- SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration_temp
  SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_IterationLoop <- (20001:40001)[-(unlist(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration)-20000)]
}
length(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States)
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency <- c()
for (t in 1:length(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States)){
  SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency <- c(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency,
                                                               length(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesIteration[[t]]))
}
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency
which.max(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency)
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ <- SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_States[[which.max(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LSZ_StatesFrequency)]]
table(c(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%1:ncol(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ)),SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3),dnn = c("",""))

require(GreedyEPL) # obtain the summarized Z by the greedy algorithm proposed by Rastelli and Friel (2018)
Z_temp <- c()
for (t in 20001:40001){
  Z_temp <- rbind(Z_temp,c(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]]%*%1:ncol(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Z[[t]])))
}
output <- MinimiseEPL(Z_temp, list(Kup = 10, loss_type = "VI",
                                   decision_init = c(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%
                                                       1:ncol(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ))))
table(output$decision,SS2_ZIPSBM_N75_K3_LSZ%*%c(1:3),dnn = c("",""))
output$EPL # VI posterior loss: 0

## Summarize p
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp <-
  mean(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001])
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp
plot(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001],type = "l")
hist(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$p[20001:40001])

## Summarize nu
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu <- matrix(0,nrow(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$nu[[1]]),ncol(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$nu[[1]]))
for (t in 20001:40001){
  SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu <- SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu +
    SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1$nu[[t]]
}
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu <- SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu/20001

## Summarize the P_m0
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0 <- matrix(0,3,3)
for (k1 in 1:3){
  for (k2 in 1:3){
    SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0[k1,k2] <-
      sum(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu[SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1,SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2])/
      (length(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizednu[SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1,SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2][
        SS2_ZIPSBM_N75_K3$Y[SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1,SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2]==0])-
         sum((SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k1)*(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ%*%c(1:3)==k2)))
  }
}
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0

# Check ExactICL with summarized Z and nu
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_ICPCLandOptR <-
  Directed_ZIPSBM_ExactICL(Y = SS2_ZIPSBM_N75_K3$Y,
                           ProbObs0Missing0 = SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0,
                           Z = SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ,
                           alpha1=1,alpha2=1, beta1=1,beta2=9, alpha=1)
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_ICPCLandOptR$value # -8120.121


# Further inference of Lambda,Pi conditional on summarized z
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s <-
  Directed_ZIPSBM_Gibbs_FixedZ(Y = SS2_ZIPSBM_N75_K3$Y,
                                      K = 3, T = 40000,beta1 = 1, beta2 = 9,
                                      Z_0 = SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedZ)
# Summarize Inferred Pi
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedPi <-
  apply(array(unlist(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi),
              dim = c(nrow(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi[[1]]),
                      ncol(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi[[1]]),
                      length(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi)))[,,20001:40001],1,mean)
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedPi

# ## Summarize Lambda
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedLambda <-
  apply(array(unlist(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Lambda),
              dim = c(nrow(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Lambda[[1]]),
                      ncol(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Lambda[[1]]),
                      length(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Lambda)))[,,20001:40001],1:2,mean)
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedLambda
```

Note here that, different from SS1 ZIP-SBM cases, the SS2 ZIP-SBM cases fit well to the network data and thus further inference conditional on the summarized clustering is implemented in this experiment.

Some extra summary statistics follow:

``` r
## Posterior mean Pi
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanPi <-
  apply(array(unlist(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Pi),
              dim = c(nrow(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Pi[[1]]),
                      ncol(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Pi[[1]]),
                      length(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Pi)))[,,20001:40001],1,mean)
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanPi
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedPi # Compare with summarized posterior samples

##  Posterior mean Lambda
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanLambda <-
  apply(array(unlist(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Lambda),
              dim = c(nrow(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Lambda[[1]]),
                      ncol(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Lambda[[1]]),
                      length(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_LS$Lambda)))[,,20001:40001],1:2,mean)
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanLambda 
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedLambda # Compare with summarized posterior samples

# p Posterior mean conditional on Z_s
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s <-
  mean(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$p[20001:40001])
hist(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$p[20001:40001])
plot(density(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$p[20001:40001]))
plot(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$p[20001:40001], type = "l")
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp # Compare with summarized posterior samples

# Compare P_m0 evaluated by posterior mean
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp/
  (SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp+
     (1-SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Summarizedp)*
     dpois(0,SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanLambda))
# with P_m0 evaluated by posterior mean conditional on Z_s
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s/
  (SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s+
     (1-SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_PosteriorMeanp_CondZ_s)*
     dpois(0,SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedLambda))
# and with summarized P_m0
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_SummarizedProbObs0Missing0
```

### 2.3 SS2 Figures Recovery

In this section, we provide the code for the figures we illustrate in the paper.
The Figure $6$ which shows the plots of the observed adjacency matrix $\boldsymbol{Y}$ in Section $4.2$ of the paper can be recovered by:

``` r
par(mfrow=c(1,3),mai = c(0.05, 0.05, 0.2, 0.05),mgp=c(0.1,0.1,0))
image(t(SS_MDAZIPSBMp_N75_K3$Y),axes = FALSE,xlab = "",ylab = "",main = "Adjacency Matrix Y")
image(t(SS_MDAZIPSBMp_N75_K3$Y)[order(SS_MDAZIPSBMp_N75_K3_LSZ%*%c(1:3)),rev(order(SS_MDAZIPSBMp_N75_K3_LSZ%*%c(1:3)))],axes = FALSE,xlab = "",ylab = "",main = TeX(r'($Y|$ True $z^*$)',bold = TRUE))
group_counts <- (as.numeric(table(SS_MDAZIPSBMp_N75_K3_LSZ%*%c(1:3))))
abline(v = -1/(2*(nrow(SS_MDAZIPSBMp_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS_MDAZIPSBMp_N75_K3$Y)-1))))
abline(h = 1-(-1/(2*(nrow(SS_MDAZIPSBMp_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS_MDAZIPSBMp_N75_K3$Y)-1)))))

image(t(1*(SS_MDAZIPSBMp_N75_K3$Y!=0))[order(SS_MDAZIPSBMp_N75_K3_LSZ%*%c(1:3)),rev(order(SS_MDAZIPSBMp_N75_K3_LSZ%*%c(1:3)))],axes = FALSE,xlab = "",ylab = "",main = TeX(r'($Y\neq 0|$ True $z^*$)',bold = TRUE))
group_counts <- (as.numeric(table(SS_MDAZIPSBMp_N75_K3_LSZ%*%c(1:3))))
abline(v = -1/(2*(nrow(SS_MDAZIPSBMp_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS_MDAZIPSBMp_N75_K3$Y)-1))))
abline(h = 1-(-1/(2*(nrow(SS_MDAZIPSBMp_N75_K3$Y)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(SS_MDAZIPSBMp_N75_K3$Y)-1)))))
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42),mgp=c(3,1,0))
```

The rand index plots in Figure $7$ can be recovered by:

``` r
par(mfrow=c(2,2),mai = c(0.3, 0.3, 0.15, 0.05), mgp=c(0.9,0.2,0))
plot(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K2_T40000_1_LSZ_RI,type = "l",xlab = "",ylab = "", main = "",cex.axis = 0.8, ylim = c(0.45,1),col = 2,lty = 1)
lines(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K2_T40000_1_LSZ_RI,type = "l", col = 1,lty = 2)
title(xlab = "Iteration",ylab = "RI", main = "K=2 Cases Rand Index", mgp=c(0.9,0.1,0),cex.main=0.8,cex.lab = 0.8)

plot(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_LSZ_RI,type = "l",xlab = "",ylab = "", main = "",cex.axis = 0.8, ylim = c(0.45,1),col = 2,lty = 1)
lines(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_LSZ_RI,type = "l", col = 1,lty = 2)
title(xlab = "Iteration",ylab = "RI", main = "K=3 Cases Rand Index", mgp=c(0.9,0.1,0),cex.main=0.8,cex.lab = 0.8)

plot(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K4_T40000_1_LSZ_RI,type = "l",xlab = "",ylab = "", main = "",cex.axis = 0.8, ylim = c(0.45,1),col = 2,lty = 1)
lines(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K4_T40000_1_LSZ_RI,type = "l", col = 1,lty = 2)
title(xlab = "Iteration",ylab = "RI", main = "K=4 Cases Rand Index", mgp=c(0.9,0.1,0),cex.main=0.8,cex.lab = 0.8)

plot(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K5_T40000_1_LSZ_RI,type = "l",xlab = "",ylab = "", main = "",cex.axis = 0.8, ylim = c(0.45,1),col = 2,lty = 1)
lines(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K5_T40000_1_LSZ_RI,type = "l", col = 1,lty = 2)
title(xlab = "Iteration",ylab = "RI", main = "K=5 Cases Rand Index", mgp=c(0.9,0.1,0),cex.main=0.8,cex.lab = 0.8)
legend("bottomright", legend=c("ZINB-SBM","ZIP-SBM"),
       col=1:2, lty = 2:1, cex=0.6)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42),mgp=c(3,1,0))
```

The posterior density plots in Figure $8$ based on the code:

``` r
SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_InferredLambda <-
  array(unlist(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Lambda),
        dim = c(nrow(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Lambda[[1]]),
                ncol(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Lambda[[1]]),
                length(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Lambda)))[,,20001:40001]
par(mfrow=c(2,2),mai = c(0.3, 0.25, 0.2, 0.05), mgp=c(1.25,0.5,0))
plot(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1$p[20001:40001],bw=0.002),ylim=c(0,35),ylab = "",xlab="", main = TeX(r'(Posterior Density of $p$)'), col = 2, xlim=c(0.09,0.2))
abline(v=0.15,col = 1,lty=2)
lines(density(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1$p[20001:40001],bw=0.002), col = 3)
par(xpd=TRUE)
text(0.1425,10, TeX(r'($p^*$)'), pos = 4,col=1)
par(xpd=FALSE)
legend("topleft", legend=c("ZINB-SBM","ZIP-SBM"),
       col=2:3, lty = 1, cex=0.6)

plot(density(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_InferredLambda[1,1,]),col = 1, xlim=c(0,2.65),ylim=c(0,30), ylab = "",xlab="", main = TeX(r'(Posterior Density of $\lambda_{gh}$|\widetilde{\textbf{z}})'))
lines(density(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_InferredLambda[1,2,]),col = 2)
lines(density(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_InferredLambda[1,3,]),col = 3)
lines(density(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_InferredLambda[2,1,]),col = 4)
lines(density(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_InferredLambda[2,2,]),col = 5)
lines(density(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_InferredLambda[2,3,]),col = 6)
lines(density(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_InferredLambda[3,1,]),col = 7)
lines(density(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_InferredLambda[3,2,]),col = 8)
lines(density(SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_InferredLambda[3,3,]),col = "rosybrown")
legend("topright", legend=c("1,1","1,2","1,3", "2,1","2,2","2,3", "3,1","3,2","3,3"),
       col=c(1:8,"rosybrown"), lty = 1, cex=0.6)

plot(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredR[1,1,],bw=1.75),col = 1, xlim=c(-3,41),ylim=c(0,0.175), ylab = "",xlab="", main = TeX(r'(Posterior Density of $r_{gh}$|\widetilde{\textbf{z}})'))
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredR[1,2,],bw=1.75),col = 2)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredR[1,3,],bw=1.75),col = 3)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredR[2,1,],bw=1.75),col = 4)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredR[2,2,],bw=1.75),col = 5)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredR[2,3,],bw=1.75),col = 6)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredR[3,1,],bw=1.75),col = 7)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredR[3,2,],bw=1.75),col = 8)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredR[3,3,],bw=1.75),col = "rosybrown")
legend("topright", legend=c("1,1","1,2","1,3", "2,1","2,2","2,3", "3,1","3,2","3,3"),
       col=c(1:8,"rosybrown"), lty = 1, cex=0.6)

plot(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredQ[1,1,],bw=0.01),col = 1, xlim=c(0.8,1.025),ylim=c(0,40), ylab = "",xlab="", main = TeX(r'(Posterior Density of $q_{gh}$|\widetilde{\textbf{z}})'))
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredQ[1,2,],bw=0.01),col = 2)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredQ[1,3,],bw=0.01),col = 3)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredQ[2,1,],bw=0.01),col = 4)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredQ[2,2,],bw=0.01),col = 5)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredQ[2,3,],bw=0.01),col = 6)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredQ[3,1,],bw=0.01),col = 7)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredQ[3,2,],bw=0.01),col = 8)
lines(density(SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_InferredQ[3,3,],bw=0.01),col = "rosybrown")
legend("topleft", legend=c("1,1","1,2","1,3", "2,1","2,2","2,3", "3,1","3,2","3,3"),
       col=c(1:8,"rosybrown"), lty = 1, cex=0.6)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42),mgp=c(3,1,0))
```

The comparison of the Poisson distribution with summarized or reference parameters and the Negative-Binomial distribution with summarized parameters shown in Figure $9$ can be recovered by:

``` r
par(mfrow=c(3,3),mai = c(0.2, 0.2, 0.2, 0.1),mgp=c(0.75,0.25,0))
plot(0:10,dpois(0:10,SS_MDAZIPSBMp_N75_K3_obs_InitialLambda[1,1]),xlab = "",ylab = "", type = "b", lty = 1, pch = 16, col = 1)
lines(0:10,dnbinom(0:10,SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredR[1,1],
                   SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredQ[1,1]),col = 2, type = "b", lty = 2, pch = 17)
lines(0:10,dpois(0:10,SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredLambda[1,1]),xlab = "",ylab = "", type = "b", lty = 3, pch = 18,col = 3)
title(xlab = "",ylab = "", main = "g=1,h=1", mgp=c(1,1.1,0),cex.main=0.8,cex.lab = 0.8)
legend("topright", legend=c("True Reference",TeX(r'($NB(\widetilde{r}_{gh},\widetilde{q}_{gh})$)'),TeX(r'($Pois(\widetilde{\lambda}_{gh})$)')),
       col=1:3, lty = 1:3, cex=0.6)

plot(0:10,dpois(0:10,SS_MDAZIPSBMp_N75_K3_obs_InitialLambda[1,2]),xlab = "",ylab = "", type = "b", lty = 1, pch = 16, col = 1)
lines(0:10,dnbinom(0:10,SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredR[1,2],
                   SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredQ[1,2]),col = 2, type = "b", lty = 2, pch = 17)
lines(0:10,dpois(0:10,SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredLambda[1,2]),xlab = "",ylab = "", type = "b", lty = 3, pch = 18,col = 3)
title(xlab = "",ylab = "", main = "g=1,h=1", mgp=c(1,1.1,0),cex.main=0.8,cex.lab = 0.8)

plot(0:10,dpois(0:10,SS_MDAZIPSBMp_N75_K3_obs_InitialLambda[1,3]),xlab = "",ylab = "", type = "b", lty = 1, pch = 16, col = 1)
lines(0:10,dnbinom(0:10,SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredR[1,3],
                   SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredQ[1,3]),col = 2, type = "b", lty = 2, pch = 17)
lines(0:10,dpois(0:10,SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredLambda[1,3]),xlab = "",ylab = "", type = "b", lty = 3, pch = 18,col = 3)
title(xlab = "",ylab = "", main = "g=1,h=1", mgp=c(1,1.1,0),cex.main=0.8,cex.lab = 0.8)

plot(0:10,dpois(0:10,SS_MDAZIPSBMp_N75_K3_obs_InitialLambda[2,1]),xlab = "",ylab = "", type = "b", lty = 1, pch = 16, col = 1)
lines(0:10,dnbinom(0:10,SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredR[2,1],
                   SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredQ[2,1]),col = 2, type = "b", lty = 2, pch = 17)
lines(0:10,dpois(0:10,SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredLambda[2,1]),xlab = "",ylab = "", type = "b", lty = 3, pch = 18,col = 3)
title(xlab = "",ylab = "", main = "g=1,h=1", mgp=c(1,1.1,0),cex.main=0.8,cex.lab = 0.8)

plot(0:10,dpois(0:10,SS_MDAZIPSBMp_N75_K3_obs_InitialLambda[2,2]),xlab = "",ylab = "", type = "b", lty = 1, pch = 16, col = 1)
lines(0:10,dnbinom(0:10,SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredR[2,2],
                   SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredQ[2,2]),col = 2, type = "b", lty = 2, pch = 17)
lines(0:10,dpois(0:10,SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredLambda[2,2]),xlab = "",ylab = "", type = "b", lty = 3, pch = 18,col = 3)
title(xlab = "",ylab = "", main = "g=1,h=1", mgp=c(1,1.1,0),cex.main=0.8,cex.lab = 0.8)

plot(0:10,dpois(0:10,SS_MDAZIPSBMp_N75_K3_obs_InitialLambda[2,3]),xlab = "",ylab = "", type = "b", lty = 1, pch = 16, col = 1)
lines(0:10,dnbinom(0:10,SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredR[2,3],
                   SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredQ[2,3]),col = 2, type = "b", lty = 2, pch = 17)
lines(0:10,dpois(0:10,SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredLambda[2,3]),xlab = "",ylab = "", type = "b", lty = 3, pch = 18,col = 3)
title(xlab = "",ylab = "", main = "g=1,h=1", mgp=c(1,1.1,0),cex.main=0.8,cex.lab = 0.8)

plot(0:10,dpois(0:10,SS_MDAZIPSBMp_N75_K3_obs_InitialLambda[3,1]),xlab = "",ylab = "", type = "b", lty = 1, pch = 16, col = 1)
lines(0:10,dnbinom(0:10,SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredR[3,1],
                   SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredQ[3,1]),col = 2, type = "b", lty = 2, pch = 17)
lines(0:10,dpois(0:10,SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredLambda[3,1]),xlab = "",ylab = "", type = "b", lty = 3, pch = 18,col = 3)
title(xlab = "",ylab = "", main = "g=1,h=1", mgp=c(1,1.1,0),cex.main=0.8,cex.lab = 0.8)

plot(0:10,dpois(0:10,SS_MDAZIPSBMp_N75_K3_obs_InitialLambda[3,2]),xlab = "",ylab = "", type = "b", lty = 1, pch = 16, col = 1)
lines(0:10,dnbinom(0:10,SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredR[3,2],
                   SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredQ[3,2]),col = 2, type = "b", lty = 2, pch = 17)
lines(0:10,dpois(0:10,SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredLambda[3,2]),xlab = "",ylab = "", type = "b", lty = 3, pch = 18,col = 3)
title(xlab = "",ylab = "", main = "g=1,h=1", mgp=c(1,1.1,0),cex.main=0.8,cex.lab = 0.8)

plot(0:10,dpois(0:10,SS_MDAZIPSBMp_N75_K3_obs_InitialLambda[3,3]),xlab = "",ylab = "", type = "b", lty = 1, pch = 16, col = 1)
lines(0:10,dnbinom(0:10,SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredR[3,3],
                   SS2_MDAZINBSBMRghQghpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredQ[3,3]),col = 2, type = "b", lty = 2, pch = 17)
lines(0:10,dpois(0:10,SS2_MDAZIPSBMpALLZ_N75_K3_Fixed_K3_T40000_1_SummarizedInferredLambda[3,3]),xlab = "",ylab = "", type = "b", lty = 3, pch = 18,col = 3)
title(xlab = "",ylab = "", main = "g=1,h=1", mgp=c(1,1.1,0),cex.main=0.8,cex.lab = 0.8)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42),mgp=c(3,1,0))
```

The extra posterior density plots of $\boldsymbol{\Pi}$ conditional on the summarized clustering for the ZINB-SBM and ZIP-SBM $K=3$ cases can be checked by the code below, respectively.

``` r
# SS2 ZINB-SBM posterior density plots of Pi conditional on Z_s
SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredPi <- 
  array(unlist(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi),
        dim = c(nrow(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi[[1]]),
                ncol(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi[[1]]),
                length(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi)))[,,20001:40001]
plot(density(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredPi[1,],bw=0.01),col = 2, xlim=c(0,0.7),ylim=c(0,10), ylab = "",xlab="", main = TeX(r'(Posterior Density of $\pi_k$|\widetilde{\textbf{z}})'))
lines(density(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredPi[2,],bw=0.01),col = 3)
lines(density(SS2_ZINBSBM_N75_K3_Fixed_K3_T40000_1_InferredPi[3,],bw=0.01),col = 4)
abline(v=0.3,col = 2,ylim=c(0,9),lty=2)
abline(v=0.3,col = 3,ylim=c(0,9),lty=3)
abline(v=0.4,col = 4,ylim=c(0,9),lty=2)
legend("topright", legend=c(TeX(r'($\pi_1$)'),TeX(r'($\pi_2$)'),TeX(r'($\pi_3$)')),
       col=2:4, lty = 1, cex=0.6)

# SS2 ZIP-SBM posterior density plots of Pi conditional on Z_s
SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_InferredPi <- 
  array(unlist(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi),
        dim = c(nrow(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi[[1]]),
                ncol(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi[[1]]),
                length(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_Further40000InferCondZ_s$Pi)))[,,20001:40001]
plot(density(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_InferredPi[1,],bw=0.01),col = 2, xlim=c(0,0.7),ylim=c(0,10), ylab = "",xlab="", main = TeX(r'(Posterior Density of $\pi_k$|\widetilde{\textbf{z}})'))
lines(density(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_InferredPi[2,],bw=0.01),col = 3)
lines(density(SS2_ZIPSBM_N75_K3_Fixed_K3_T40000_1_InferredPi[3,],bw=0.01),col = 4)
abline(v=0.3,col = 2,ylim=c(0,9),lty=2)
abline(v=0.3,col = 3,ylim=c(0,9),lty=3)
abline(v=0.4,col = 4,ylim=c(0,9),lty=2)
legend("topright", legend=c(TeX(r'($\pi_1$)'),TeX(r'($\pi_2$)'),TeX(r'($\pi_3$)')),
       col=2:4, lty = 1, cex=0.6)
```
