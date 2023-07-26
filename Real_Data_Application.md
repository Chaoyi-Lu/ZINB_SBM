# Real Data Application

This markdown file provides the code and implementations for the real data application (RDA) illustrated in the paper *Zero-Inflated Negative-Binomial Stochastic Block Model*.
This RDA focuses on the `UKfaculty` real dataset which provides the personal friendship network of a faculty of a UK university.
The network can be loaded from the library `igraphdata` by the code:

``` r
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
# Real Data Application UKfaculty Dataset
library(igraphdata)
data("UKfaculty")
```

The adjacency matrix of such a real network can be extracted by the functions provided in the library `ìgraph`, that is,

``` r
library(igraph)
UKfaculty_adj <- as.matrix(get.adjacency(UKfaculty,attr = "weight")) # extract the observed adjacency matrix
```

Then we can check the proportion of zero interactions and the distribution of non-zero edge weights of the network:

``` r
# Check proportion of zeros in Y
(sum(UKfaculty_adj==0)-dim(UKfaculty_adj)[1])/(prod(dim(UKfaculty_adj))-dim(UKfaculty_adj)[1]) # sparse network
# Check edge weight distribution (without diagonal entries of the adjacency matrix)
table(c(UKfaculty_adj[upper.tri(UKfaculty_adj)],UKfaculty_adj[lower.tri(UKfaculty_adj)]))
```

The affiliation of each member within the network is treated as the reference clustering and can be extracted via:

``` r
# Check affiliation which is treated as the reference clustering
V(UKfaculty)$Group
```

Similar to the simulation studies, we first apply the label switching on the reference clustering in order to ensure the uniqueness of the clustering.

``` r
# Label switch the reference clustering
Real_data_application_UKfaculty_Z <- matrix(0,length(V(UKfaculty)$Group),length(table(V(UKfaculty)$Group)))
for (k in 1:length(table(V(UKfaculty)$Group))){
  Real_data_application_UKfaculty_Z[V(UKfaculty)$Group==k,k] <- 1 # transform to matrix form
}
Real_data_application_UKfaculty_Z <- LabelSwitching_SG2003_ZINBSBM(Z = list(Real_data_application_UKfaculty_Z))$Z[[1]]
```

We can also check the mean and variance of the interaction strengths:

``` r
# Check the mean and var of the interaction weights (diagonal removed)
mean(c(UKfaculty_adj[upper.tri(UKfaculty_adj)],
       c(UKfaculty_adj[lower.tri(UKfaculty_adj)])))
var(c(UKfaculty_adj[upper.tri(UKfaculty_adj)],
      c(UKfaculty_adj[lower.tri(UKfaculty_adj)])))
```

## The Implementations of the RDA

Recall here that we apply the partially collapsed Metropolis within Gibbs algorithm for the inference of the ZINB-SBM on the `UKfaculty` real network.
We consider fixed $K=2$ to $K=8$, and each fixed $K$ case is proposed to be implemented for $10$ times/rounds in order to reduce the effect of the mixing or the bad initial state problem.
Moreover, since multiple implementations are applied, we also propose to implement each round for $20,000$ iterations which are believed to be long enough to converge according to the performance of the simulation studies.
The prior settings are almost the same as those simulation study ZINB-SBM cases except the $p$ prior which is set to be more concentrated around $0.1$ in order for better calibration of the $p$ posterior chain, that is, $p \sim \text{Beta}(20,180)$.
Otherwise, the $p$ posterior samples may stuck around the high values which are believed to be around the local posterior mode, and it might take quite long time to move away from the local posterior mode due to the high dimensionality of the ZINB-SBM.
The implementations are applied one by one from $K=2$ round $1$ to $K=8$ round $10$ as shown below.
The reference implementation time is provided, and we don't suggest the readers to implement all of them again because it will take around $7$ days to finish.
Instead, we provide in this repository the outputs of the best pick we obtained in our experiments shown in the paper as an example later in this section.

``` r
# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 2, T = 20000, Round 1 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_1 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 2, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_1_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_1_time # Time difference of 1.148557 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_R1.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 2, T = 20000, Round 2 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_2 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 2, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_2_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_2_time # Time difference of 1.075736 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_R2.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 2, T = 20000, Round 3 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_3 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 2, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_3_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_3_time # Time difference of 1.071965 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_R3.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 2, T = 20000, Round 4 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_4 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 2, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_4_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_4_time # Time difference of 1.050778 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_R4.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 2, T = 20000, Round 5 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_5 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 2, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_5_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_5_time # Time difference of 1.046858 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_R5.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 2, T = 20000, Round 6 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_6 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 2, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_6_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_6_time # Time difference of 1.120932 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_R6.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 2, T = 20000, Round 7 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_7 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 2, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_7_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_7_time # Time difference of 1.08436 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_R7.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 2, T = 20000, Round 8 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_8 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 2, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_8_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_8_time # Time difference of 1.08621 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_R8.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 2, T = 20000, Round 9 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_9 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 2, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_9_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_9_time # Time difference of 1.08532 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_R9.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 2, T = 20000, Round 10 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_10 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 2, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_10_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_10_time # Time difference of 1.078473 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K2_Prior_p_Beta_20_180_T20000_R10.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 3, T = 20000, Round 1 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_1 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 3, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_1_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_1_time # Time difference of 1.571608 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_R1.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 3, T = 20000, Round 2 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_2 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 3, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_2_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_2_time # Time difference of 1.430506 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_R2.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 3, T = 20000, Round 3 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_3 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 3, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_3_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_3_time # Time difference of 1.436394 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_R3.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 3, T = 20000, Round 4 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_4 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 3, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_4_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_4_time # Time difference of 1.481325 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_R4.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 3, T = 20000, Round 5 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")

start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_5 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 3, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_5_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_5_time # Time difference of 1.421465 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_R5.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 3, T = 20000, Round 6 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_6 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 3, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_6_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_6_time # Time difference of 1.450154 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_R6.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 3, T = 20000, Round 7 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_7 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 3, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_7_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_7_time # Time difference of 1.460419 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_R7.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 3, T = 20000, Round 8 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_8 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 3, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_8_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_8_time # Time difference of 1.488658 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_R8.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 3, T = 20000, Round 9 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_9 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 3, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_9_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_9_time # Time difference of 1.455838 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_R9.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 3, T = 20000, Round 10 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_10 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 3, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_10_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_10_time # Time difference of 1.507368 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K3_Prior_p_Beta_20_180_T20000_R10.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 4, T = 20000, Round 1 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_1 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 4, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_1_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_1_time # Time difference of 1.915564 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_R1.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 4, T = 20000, Round 2 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_2 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 4, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_2_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_2_time # Time difference of 1.831008 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_R2.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 4, T = 20000, Round 3 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_3 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 4, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_3_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_3_time # Time difference of 1.816023 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_R3.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 4, T = 20000, Round 4 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_4 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 4, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_4_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_4_time # Time difference of 1.80309 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_R4.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 4, T = 20000, Round 5 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_5 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 4, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_5_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_5_time # Time difference of 1.892018 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_R5.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 4, T = 20000, Round 6 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_6 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 4, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_6_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_6_time # Time difference of 1.883947 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_R6.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 4, T = 20000, Round 6 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_7 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 4, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_7_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_7_time # Time difference of 1.830547 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_R7.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 4, T = 20000, Round 8 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_8 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 4, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_8_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_8_time # Time difference of 1.843263 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_R8.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 4, T = 20000, Round 9 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_9 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 4, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_9_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_9_time # Time difference of 1.871217 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_R9.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 4, T = 20000, Round 10 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_10 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 4, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_10_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_10_time # Time difference of 1.877589 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K4_Prior_p_Beta_20_180_T20000_R10.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 5, T = 20000, Round 1 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_1 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj,
                         K = 5, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_1_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_1_time # Time difference of 2.344764 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R1.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 5, T = 20000, Round 2 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_2 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj,
                         K = 5, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_2_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_2_time # Time difference of 2.326425 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R2.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 5, T = 20000, Round 3 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_3 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj,
                         K = 5, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_3_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_3_time # Time difference of 2.297538 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R3.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 5, T = 20000, Round 4 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj,
                         K = 5, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_time # Time difference of 2.271394 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R4.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 5, T = 20000, Round 5 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_5 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj,
                         K = 5, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_5_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_5_time # Time difference of 2.228516 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R5.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 5, T = 20000, Round 6 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_6 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 5, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_6_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_6_time # Time difference of 2.241017 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R6.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 5, T = 20000, Round 7 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_7 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 5, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_7_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_7_time # Time difference of 2.250606 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R7.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 5, T = 20000, Round 8 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_8 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 5, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_8_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_8_time # Time difference of 2.386486 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R8.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 5, T = 20000, Round 9 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_9 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 5, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_9_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_9_time # Time difference of 2.274145 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R9.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 5, T = 20000, Round 10 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_10 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 5, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_10_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_10_time # Time difference of 2.260947 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R10.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 6, T = 20000, Round 1 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_1 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj,
                         K = 6, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_1_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_1_time # Time difference of 2.794361 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_R1.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 6, T = 20000, Round 2 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_2 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj,
                         K = 6, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_2_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_2_time # Time difference of 2.739674 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_R2.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 6, T = 20000, Round 3 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_3 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj,
                         K = 6, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_3_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_3_time # Time difference of 2.785648 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_R3.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 6, T = 20000, Round 4 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_4 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj,
                         K = 6, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_4_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_4_time # Time difference of 2.767582 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_R4.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 6, T = 20000, Round 5 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_5 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj,
                         K = 6, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_5_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_5_time # Time difference of 2.75266 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_R5.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 6, T = 20000, Round 6 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_6 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 6, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_6_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_6_time # Time difference of 2.703092 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_R6.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 6, T = 20000, Round 7 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_7 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 6, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_7_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_7_time # Time difference of 2.697717 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_R7.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 6, T = 20000, Round 8 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_8 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 6, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_8_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_8_time # Time difference of 2.649882 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_R8.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 6, T = 20000, Round 9 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_9 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 6, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_9_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_9_time # Time difference of 2.772085 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_R9.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 6, T = 20000, Round 10 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_10 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 6, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_10_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_10_time # Time difference of 2.664274 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K6_Prior_p_Beta_20_180_T20000_R10.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 7, T = 20000, Round 1 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_1 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 7, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_1_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_1_time # Time difference of 3.148067 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_R1.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 7, T = 20000, Round 2 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_2 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 7, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_2_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_2_time # Time difference of 3.086314 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_R2.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 7, T = 20000, Round 3 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_3 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 7, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_3_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_3_time # Time difference of 3.020107 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_R3.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 7, T = 20000, Round 4 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_4 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 7, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_4_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_4_time # Time difference of 3.014228 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_R4.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 7, T = 20000, Round 5 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_5 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 7, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_5_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_5_time # Time difference of 3.069248 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_R5.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 7, T = 20000, Round 6 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_6 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 7, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_6_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_6_time # Time difference of 3.281996 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_R6.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 7, T = 20000, Round 7 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_7 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 7, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_7_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_7_time # Time difference of 3.205415 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_R7.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 7, T = 20000, Round 8 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_8 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 7, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_8_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_8_time # Time difference of 3.252727 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_R8.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 7, T = 20000, Round 9 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_9 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 7, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_9_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_9_time # Time difference of 3.213382 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_R9.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 7, T = 20000, Round 10 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_10 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 7, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_10_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_10_time # Time difference of 3.253955 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K7_Prior_p_Beta_20_180_T20000_R10.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 8, T = 20000, Round 1 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_1 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 8, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_1_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_1_time # Time difference of 3.768132 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_R1.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 8, T = 20000, Round 2 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_2 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 8, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_2_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_2_time # Time difference of 3.487546 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_R2.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 8, T = 20000, Round 3 with p prior Beta(20,180)

rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_3 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 8, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_3_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_3_time # Time difference of 3.505183 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_R3.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 8, T = 20000, Round 4 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_4 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 8, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_4_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_4_time # Time difference of 3.519141 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_R4.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 8, T = 20000, Round 5 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_5 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 8, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_5_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_5_time # Time difference of 3.495286 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_R5.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 8, T = 20000, Round 6 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_6 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 8, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_6_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_6_time # Time difference of 3.945458 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_R6.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 8, T = 20000, Round 7 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_7 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 8, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_7_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_7_time # Time difference of 3.911281 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_R7.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 8, T = 20000, Round 8 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_8 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 8, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_8_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_8_time # Time difference of 3.840511 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_R8.RData")

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 8, T = 20000, Round 9 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_9 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 8, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_9_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_9_time # Time difference of 3.871979 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_R9.RData")
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 8, T = 20000, Round 10 with p prior Beta(20,180)
rm(list=ls())
gc()
source("Functions_for_ZINB_SBM.R")
start.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_10 <- 
  Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj, 
                         K = 8, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
end.time <- Sys.time()
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_10_time <- end.time - start.time
RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_10_time # Time difference of 3.800605 hours
# save.image("RDA_UKfaculty_ZINBSBM_Fixed_K8_Prior_p_Beta_20_180_T20000_R10.RData")
```

In the case of that the memory of the laptop is enough to load all the data, the practitioners can consider not to store each implementation above and then remove the data in order to save the memory for the next implementation.
Otherwise, it's suggested to do so if practitioners prefer recovering all the work above.

Though the summarizing processes are almost the same as those of the simulation study ZINB-SBM cases, we briefly go through the processes again in this section.
Since the processes are also similar for all the implementations in RDA, we take the best round implementation shown in Section $5$ of the paper as an example here, that is, the $K=5$ round $4$ case.
Once we obatined the outputs from the PCMwG function, we can first apply the label switching on the clustering and those clustering dependent parameters.

``` r
# # # Real Data Application ZINB-SBM Metropolis within Gibbs algorithm for UKfaculty dataset, K = 5, T = 20000, Round 4 with p prior Beta(20,180)
# rm(list=ls())
# gc()
# source("Functions_for_ZINB_SBM.R")
# start.time <- Sys.time()
# RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4 <- 
#   Directed_ZINBSBM_PCMwG(Y = UKfaculty_adj,
#                          K = 5, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180)
# end.time <- Sys.time()
# RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_time <- end.time - start.time
# RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_time # Time difference of 2.271394 hours
# # save.image("RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R4.RData")

# # Summarize the outputs
## Apply label switching
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS <-
  LabelSwitching_SG2003_ZINBSBM(Z = RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$Z,
                                Pi = RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$Pi,
                                R = RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$R,
                                Q = RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$Q,
                                Acceptance_count_R = RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$Acceptance_count_R)
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$Z <- c()#in order to save memory
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$Pi <- c()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$R <- c()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$Q <- c()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$Acceptance_count_R <- c()
gc()
```

Here we provide the corresponding data in this repository and the data can be loaded by the file `RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R4.RData` shown below.
The data contains $4$ variables in the environment.
`RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS` stores the label-switched posterior samples of the clustering and those clustering dependent parameters as well as the acceptance count for each $r_{gh}$ M-H step, that is, posterior samples of $\boldsymbol{z},\boldsymbol{\Pi},\boldsymbol{R},\boldsymbol{Q}$ and the `Acceptance_count_R` defined in the function `Directed_ZINBSBM_PCMwG()`.
`RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4`stores the posterior samples of clustering indepdent parameters and variables, that is, posterior samples of $\boldsymbol{\nu}, p, \boldsymbol{X}$.
`RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_time` stores the implementation time and `RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s` stores the further $20,000$-iteration inference outputs we obtained in our experiments conditional on the summarized clustering $\tilde{\boldsymbol{z}}$.

``` r
# Load the data containing the posterior chains and further posterior chains conditional on the summarized clustering Z_s
load("RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_R4.RData")
```

After the label switching process, we can check the mixing of the clustering samples by the rand index plots and we can also check some specific clustering states shown below.

``` r
## check rand index for each iteration
require("fossil")
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_RI <- c()
for (t in 1:20001){
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_RI <-
    c(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_RI,
      rand.index(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[t]]%*%c(1:ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[t]])),
                 Real_data_application_UKfaculty_Z%*%c(1:4)))
}
plot(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_RI,ylim = c(0.5,1),type = "l",xlab = "",ylab = "", main = "Rand Index",cex.axis = 0.8)
# Check some specific clustering states
table(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[1]]%*%c(1:5),Real_data_application_UKfaculty_Z%*%c(1:4)) # initial state
table(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[2]]%*%c(1:5),Real_data_application_UKfaculty_Z%*%c(1:4)) # first state
table(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[3]]%*%c(1:5),Real_data_application_UKfaculty_Z%*%c(1:4)) # second state
table(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[20001]]%*%c(1:5),Real_data_application_UKfaculty_Z%*%c(1:4)) # end state
```

The marginal posterior mode of the clustering chain can be obtained following exactly the same processes as simulation studies, that is,

``` r
## Obtain the marginal posterior mode of the Z chain
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_States <- list()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesIteration <- list()
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_IterationLoop <- 10001:20001
StatesLabelIndicator = 0 # states order indicator
while (length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_IterationLoop)!=0){
  StatesLabelIndicator <- StatesLabelIndicator + 1
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_FirstState <- RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_IterationLoop[1]]]
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_States[[StatesLabelIndicator]] <- RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_FirstState
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesIteration_temp <- c()
  for (t in RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_IterationLoop){
    if (sum(c(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[t]]%*%1:ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[t]]))==
            c(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_FirstState%*%1:ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_FirstState)))==nrow(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_FirstState)){
      RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesIteration_temp <- c(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesIteration_temp,t)
    }
  }
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesIteration[[StatesLabelIndicator]] <- RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesIteration_temp
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_IterationLoop <- (10001:20001)[-(unlist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesIteration)-10000)]
}
length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_States) # check the number of different clustering states
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesFrequency <- c()
for (t in 1:length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_States)){
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesFrequency <- 
    c(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesFrequency,
      length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesIteration[[t]]))
}
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesFrequency # check the frequency for each state
which.max(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesFrequency) # find the marginal posterior mode, i.e. the most frequent clustering state
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ <- RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_States[[which.max(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_StatesFrequency)]]
table(c(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%1:ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ)),Real_data_application_UKfaculty_Z%*%c(1:4),dnn = c("",""))
```

Then we can apply the greedy algorithm proposed by [Rastelli and Friel (2018)](https://doi.org/10.1007/s11222-017-9786-y) starting form the marginal posterior mode to obtain the summarized clustering and check whether it's the same as the marginal posterior mode.

``` r
require(GreedyEPL) # obtain the summarized Z by the greedy algorithm proposed by Rastelli and Friel (2018)
Z_temp <- c()
for (t in 10001:20001){
  Z_temp <- rbind(Z_temp,c(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[t]]%*%1:ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[t]])))
}
output <- MinimiseEPL(Z_temp, list(Kup = 10, loss_type = "VI",
                                   decision_init = c(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%
                                                       1:ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ))))
table(output$decision,Real_data_application_UKfaculty_Z%*%c(1:4),dnn = c("","")) # the output agrees with the marginal posterior mode
output$EPL # VI posterior loss: 0.03585164
```

The summarized clustering can be visualized via the adjacency matrix plot based on such clustering, that is,

``` r
## Check the summarized clustering
image(t(UKfaculty_adj)[order(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)),rev(order(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)))])
group_counts <- (as.numeric(table(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5))))
abline(v = -1/(2*(nrow(UKfaculty_adj)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(UKfaculty_adj)-1))))
abline(h = 1-(-1/(2*(nrow(UKfaculty_adj)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(UKfaculty_adj)-1)))))
```

The acceptance rate of each $r_{gh}$ M-H step can be checked by:

``` r
## Check acceptance rate
Reduce(`+`, RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Acceptance_count_R[10001:20001])/10001
```

and we can also first summarize the posterior samples of $p$ and $\boldsymbol{\nu}$ as:

``` r
## Summarize p
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizedp <-
  mean(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$p[10001:20001])
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizedp
plot(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$p[10001:20001],type = "l")
hist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$p[10001:20001])

## Summarize nu
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizednu <- 
  matrix(0,nrow(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$nu[[1]]),
         ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$nu[[1]]))
for (t in 10001:20001){
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizednu <- 
    RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizednu +
    RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$nu[[t]]
}
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizednu <- 
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizednu/10001
```

The summarized $\tilde{\boldsymbol{\nu}}$ can be visualized based on the summarized clustering:

``` r
## Check the plot of summarized nu which corresponds to the proportion of the times the corresponding y_ij is classified as missing zero
image(t(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizednu)[order(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)),
                                                                                  rev(order(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)))])
group_counts <- (as.numeric(table(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5))))
abline(v = -1/(2*(nrow(UKfaculty_adj)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(UKfaculty_adj)-1))))
abline(h = 1-(-1/(2*(nrow(UKfaculty_adj)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(UKfaculty_adj)-1)))))
```

Once we obtained the summarized $\tilde{\boldsymbol{\nu}}$ and $\tilde{\boldsymbol{z}}$, we can approximate the summarized $\widetilde{\boldsymbol{P_{m0}}}$ as:

``` r
## Summarize the P_m0
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedProbObs0Missing0 <- matrix(0,5,5)
for (k1 in 1:5){
  for (k2 in 1:5){
    RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedProbObs0Missing0[k1,k2] <- sum(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizednu[RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)==k1,RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)==k2])/
      (length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizednu[RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)==k1,RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)==k2][UKfaculty_adj[RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)==k1,RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)==k2]==0])-sum((RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)==k1)*(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)==k2)))
  }
}
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedProbObs0Missing0
```

By obtaining the posterior mean of the $\boldsymbol{R}$ as the initial state of evaluating the ICPCL:

``` r
# ## Obtain the initial R to be optimized for the ICPCL by the posterior mean of R
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanR <-
  apply(array(unlist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$R),
              dim = c(nrow(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$R[[1]]),
                      ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$R[[1]]),
                      length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$R)))[,,10001:20001],1:2,mean)
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanR 
```

we can obtain the ICPCL value via:

``` r
# Obtain the approximate ICPCL
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_ICPCLandOptR <- 
  Directed_ZINBSBM_ApproximateICPCL_OptR(Y = UKfaculty_adj,
                                         ProbObs0Missing0 = RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedProbObs0Missing0,
                                         Z = RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ,
                                         R_0 = RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanR,
                                         alpha=1, beta1=20,beta2=180, betaq1=1,betaq2=1)
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_ICPCLandOptR$value # -6087.947
```

We apply all the processes above for each round of each fixed $K$ and finally obtain the ICPCL values for all the implementations as the table shown below.

\begin{table}[H]
\centering
\scalebox{0.9}{
\begin{tabular}[t]{|c|c c c c c c c|}
\hline
 ICPCL    & K2       & K3         & K4       & K5        & K6           & K7  	     & K8       \\
\hline
Round 1   & -6396.096 & \textcolor{blue}{-6453.414} & -6414.270 & -6294.696 & -6400.087 & -6585.189 & -7724.104 \\
\hline
Round 2   & -6454.013 & \textcolor{blue}{-6444.560} &\textcolor{blue}{-6217.934} & -6401.470 & -6889.146 & \textcolor{blue}{-7212.347} & -6853.329\\
\hline
Round 3   & -6450.135 & \textcolor{blue}{-6423.273} & \textcolor{blue}{-6505.771} & -6500.255 & -7031.210 & \textcolor{blue}{-6744.229} & \textcolor{red}{-6730.280}\\
\hline
Round 4   & -6446.579 & -6431.902 & -6304.791  &\textcolor{red}{-6087.947} & -6549.627 & -6583.900 & -7823.358\\
\hline
Round 5   & \textcolor{red}{-6373.261} & \textcolor{blue}{-6501.593} & \textcolor{blue}{-7022.227} & \textcolor{blue}{-6296.962} & -6295.793 & -6469.524 & -7085.496\\
\hline
Round 6   & -6386.550 & \textcolor{blue}{-6440.589} & -6459.649 & -6244.250 & -6798.479 & -6710.810 & -7340.897\\
\hline
Round 7   & -6437.208 & -6380.615 & -6291.410 & -6298.549 & -6695.855 & \textcolor{red}{-6290.709} & -7248.848\\
\hline
Round 8   & \textcolor{blue}{-6364.818} & \textcolor{red}{-6348.968} & \textcolor{blue}{-6238.089} & \textcolor{blue}{-6457.362} & -6701.610 & -6517.335 & -6996.099\\
\hline
Round 9   & -6428.155 & -6377.889 & \textcolor{blue}{-6276.172} & -6433.268 & \textcolor{red}{-6259.197} & -6499.246 & -7007.193\\
\hline
Round 10   & -6428.155 & -6401.935 & \textcolor{red}{-6238.031} & -6294.832 & \textcolor{blue}{-6510.474} & -6679.074 & -6570.362\\
\hline
\end{tabular}
}
\end{table}
