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

The affiliation of each member within the network is treated as the reference clustering and can be extract via:

``` r
# Check affiliation which is treated as the reference clustering
V(UKfaculty)$Group
```

Similar to the simulation studies, we first apply the label switching on the reference clustering in order to ensure the uniqueness.

``` r
# Label switch the reference clustering
Real_data_application_UKfaculty_Z <- matrix(0,length(V(UKfaculty)$Group),length(table(V(UKfaculty)$Group)))
for (k in 1:length(table(V(UKfaculty)$Group))){
  Real_data_application_UKfaculty_Z[V(UKfaculty)$Group==k,k] <- 1 # transform to matrix form
}
Real_data_application_UKfaculty_Z <- LabelSwitching_SG2003_ZINBSBM(Z = list(Real_data_application_UKfaculty_Z))$Z[[1]]
c(Real_data_application_UKfaculty_Z%*%1:4)
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

Recall here that we apply the partially collapsed Metropolis within Gibbs algorithm for the ZINB-SBM for the inference on the `UKfaculty` real network.
We consider fixed $K=2$ to $K=8$, and each fixed $K$ case is implemented for $10$ times/rounds in order to reduce the effect of the mixing or the bad initial state problem.
The prior settings are almost the same as simulation study ZINB-SBM cases except the $p$ prior which is set to be more concentrated around $0.1$ in order for the better calibration of the $p$ posterior chain.
Otherwise, the $p$ posterior samples may stuck around the high values which are believed to be local posterior mode, and it might take quite long time to move away from the local posterior mode due to the high dimensionality of the ZINB-SBM model.
The implementations are applied one by one from $K=2$ round $1$ to $K=8$ round $10$ as shown below.
The reference implementation time is provided, and we don't suggest the readers to implement all of them again because it will take around $7$ days to finish.
Instead, we provide in this repository the outputs of the best pick we obtained in our experiments as an example later in this section.

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

We know that it seems to be not very brilliant to apply the implementations and then store and delete the data to move to the next.
However, our laptop memory only allows us to load few cases once and this doesn't include the further inference data.

