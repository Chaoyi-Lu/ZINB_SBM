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

The adjacency matrix plots shown as Figure $10$ in Section $5$ of the paper can be recovered by the following code.

``` r
par(mfrow=c(1,3),mai = c(0.5, 0.05, 0.25, 0.05),mgp=c(0.1,0.1,0))
image(t(UKfaculty_adj),axes = FALSE,main = "Adjacency Matrix Y")

image(t(UKfaculty_adj)[order(Real_data_application_UKfaculty_Z%*%1:4),rev(order(Real_data_application_UKfaculty_Z%*%1:4))],axes = FALSE,main = TeX(r'($Y|$ Affiliation $z^*$)',bold=TRUE))
group_counts <- (as.numeric(table(Real_data_application_UKfaculty_Z%*%1:4)))
abline(v = -1/(2*(nrow(UKfaculty_adj)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(UKfaculty_adj)-1))))
abline(h = 1-(-1/(2*(nrow(UKfaculty_adj)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(UKfaculty_adj)-1)))))

image(t(1*(UKfaculty_adj!=0))[order(Real_data_application_UKfaculty_Z%*%1:4),rev(order(Real_data_application_UKfaculty_Z%*%1:4))],axes = FALSE,main = TeX(r'($Y\neq 0|$ Affiliation $z^*$)',bold=TRUE))
group_counts <- (as.numeric(table(Real_data_application_UKfaculty_Z%*%1:4)))
abline(v = -1/(2*(nrow(UKfaculty_adj)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(UKfaculty_adj)-1))))
abline(h = 1-(-1/(2*(nrow(UKfaculty_adj)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(UKfaculty_adj)-1)))))
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42),mgp = c(3,1,0))
```

## The Implementations of the RDA

Recall here that we apply the partially collapsed Metropolis within Gibbs algorithm for the inference of the ZINB-SBM fit to the `UKfaculty` real network.
We consider fixed $K=2$ to $K=8$, and each fixed $K$ case is proposed to be implemented for $10$ times/rounds, that is, multiple implementations, in order to reduce the effect of the mixing or the bad initial state problem.
Moreover, since multiple implementations are applied, we also propose to implement each round for $20,000$ iterations which are believed to be long enough for the posterior chains to converge according to the performance of the simulation studies.
The prior settings are almost the same as those simulation study ZINB-SBM cases except the $p$ prior which is set to be more concentrated around $0.1$ in order for better calibration of the $p$ posterior chain, that is, $p \sim \text{Beta}(20,180)$.
Otherwise, the $p$ posterior samples may stuck around high values which are believed to be around the local posterior mode, and it might take quite long time to move away from the local posterior mode due to the high dimensionality of the ZINB-SBM.
The implementations are applied one by one from $K=2$ round $1$ to $K=8$ round $10$ as shown below.
The reference implementation time is provided, and we don't suggest the readers to implement all of them again because it will take around $7$ days to finish.
Instead, we provide in this repository the outputs of the best case, the $K=5$ round $4$ case, we obtained in our experiments and showed in the paper as an example later in this section.

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

In the case of that the memory of the laptop is enough to load all the data, the practitioners can consider not to store each implementation above and then remove the data for the memory saving for the next implementation.
Otherwise, it's suggested to do so if practitioners prefer recovering all the work above.

## The Summarizing Processes of the RDA

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
`RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS` stores the label-switched posterior samples of the clustering and those clustering dependent parameters as well as the acceptance count for each $r_{gh}$ M-H step, that is, label-switched posterior samples of $\boldsymbol{Z},\boldsymbol{\Pi},\boldsymbol{R},\boldsymbol{Q}$ and the `Acceptance_count_R` defined in the function `Directed_ZINBSBM_PCMwG()`.
`RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4`stores the posterior samples of clustering indepdent parameters and variables, that is, posterior samples of $\boldsymbol{\nu}, p, \boldsymbol{X}$.
`RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_time` stores the implementation time and `RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s` stores the further $20,000$-iteration inference outputs we obtained in our experiments conditional on the summarized clustering $\tilde{\boldsymbol{Z}}$.

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

The summarized clustering and the summarized $\tilde{\boldsymbol{\nu}}$ can be visualized as shown in Figure $13$ of the paper and such a figure can be recovered by:

``` r
par(mfrow=c(1,2),mai = c(0.25, 0.05, 0.275, 0.05), mgp=c(1.25,0.35,0))
image(t(UKfaculty_adj)[order(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)),rev(order(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)))],axes = FALSE, main = TeX(r'($Y|$ Summarized $\widetilde{z}$)',bold=TRUE))
group_counts <- (as.numeric(table(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5))))
abline(v = -1/(2*(nrow(UKfaculty_adj)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(UKfaculty_adj)-1))))
abline(h = 1-(-1/(2*(nrow(UKfaculty_adj)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(UKfaculty_adj)-1)))))

image(t(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizednu)[order(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)),rev(order(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5)))],axes = FALSE, main = TeX(r'(Summarized $\widetilde{\nu}|$ Summarized $\widetilde{z}$)',bold=TRUE))
group_counts <- (as.numeric(table(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ%*%c(1:5))))
abline(v = -1/(2*(nrow(UKfaculty_adj)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(UKfaculty_adj)-1))))
abline(h = 1-(-1/(2*(nrow(UKfaculty_adj)-1)) + cumsum(group_counts/sum(group_counts))*(1+2/(2*(nrow(UKfaculty_adj)-1)))))
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42),mgp=c(3,1,0))
```

Once we obtained the summarized $\tilde{\boldsymbol{\nu}}$ and $\tilde{\boldsymbol{Z}}$, we can approximate the summarized $\widetilde{\boldsymbol{P_{m0}}}$ as:

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

| ICPCL | K2       | K3    | K4       | K5        | K6           | K7  	     | K8       |
| :---: |:---: |:---: |:---: |:---: |:---: |:---: |:---: |
|Round 1   | $-6396.096$ | $\textcolor{blue}{-6453.414}$ | $-6414.270$ | $-6294.696$ | $-6400.087$ | $-6585.189$ | $-7724.104$ |
|Round 2   | $-6454.013$ | $\textcolor{blue}{-6444.560}$ |$\textcolor{blue}{-6217.934}$ | $-6401.470$ | $-6889.146$ | $\textcolor{blue}{-7212.347}$ | $-6853.329$|
|Round 3   | $-6450.135$ | $\textcolor{blue}{-6423.273}$ | $\textcolor{blue}{-6505.771}$ | $-6500.255$ | $-7031.210$ | $\textcolor{blue}{-6744.229}$ | $\textcolor{red}{-6730.280}$|
|Round 4   | $-6446.579$ | $-6431.902$ | $-6304.791$  |$\textcolor{red}{-6087.947}$ | $-6549.627$ | $-6583.900$ | $-7823.358$|
|Round 5   | $\textcolor{red}{-6373.261}$ | $\textcolor{blue}{-6501.593}$ | $\textcolor{blue}{-7022.227}$ | $\textcolor{blue}{-6296.962}$ | $-6295.793$ | $-6469.524$ | $-7085.496$|
|Round 6   | $-6386.550$ | $\textcolor{blue}{-6440.589}$ | $-6459.649$ | $-6244.250$ | $-6798.479$ | $-6710.810$ | $-7340.897$|
|Round 7   | $-6437.208$ | $-6380.615$ | $-6291.410$ | $-6298.549$ | $-6695.855$ | $\textcolor{red}{-6290.709}$ | $-7248.848$|
|Round 8   | $\textcolor{blue}{-6364.818}$ | $\textcolor{red}{-6348.968}$ | $\textcolor{blue}{-6238.089}$ | $\textcolor{blue}{-6457.362}$ | $-6701.610$ | $-6517.335$ | $-6996.099$|
|Round 9   | $-6428.155$ | $-6377.889$ | $\textcolor{blue}{-6276.172}$ | $-6433.268$ | $\textcolor{red}{-6259.197}$ | $-6499.246$ | $-7007.193$|
|Round 10  | $-6428.155$ | $-6401.935$ | $\textcolor{red}{-6238.031}$ | $-6294.832$ | $\textcolor{blue}{-6510.474}$ | $-6679.074$ | $-6570.362$|

The corresponding minimized expected posterior loss table is shown as:

| minExPosVIloss | K2       | K3    | K4       | K5        | K6           | K7  	     | K8       |
| :---: |:---: |:---: |:---: |:---: |:---: |:---: |:---: |
| Round 1 |$0.0621$ |$\textcolor{blue}{0.3450}$|$0.0363$|$0.0687$|$0.0770$|$0.0473$|$0.0795$|
| Round 2 |$0.0582$|$\textcolor{blue}{0.3518}$|$\textcolor{blue}{0.3516}$|$0.1135	$|$0.0436$|$\textcolor{blue}{2.69e\text{-}13}$|$0.1353$|
| Round 3 |$0.0575$|$\textcolor{blue}{0.3490}$|$\textcolor{blue}{0.3934}$|$0.0404$|$0.0973$|$\textcolor{blue}{0.2482}$|$\textcolor{red}{0.1141}$|
| Round 4 |$0.0639$|$0.0900$|$0.0507$|$\textcolor{red}{0.0359}$|$0.1225$|$0.0515$|$0.1031$|
| Round 5 |$\textcolor{red}{0.0553}$|$\textcolor{blue}{0.3284}$|$\textcolor{blue}{0.0025}$|$\textcolor{blue}{0.4050}$|$0.0705$|$0.0729$|$0.1525$|
| Round 6 |$0.0615$|$\textcolor{blue}{0.3194}$|$0.0917$|$0.1034$|$0.1369$|$0.0523$|$0.1057$|
| Round 7 |$0.0623$|$0.0267$|$0.0786$|$0.0700$|$0.0417$|$\textcolor{red}{0.0765}$|$0.0666$|
| Round 8 |$\textcolor{blue}{0.0106}$|$\textcolor{red}{0.1112}$|$\textcolor{blue}{0.3330}$|$\textcolor{blue}{4.52e\text{-}05}$|$0.1053$|$0.0312$|$0.0442$|
| Round 9 |$0.0599$|$0.1335$|$\textcolor{blue}{0.0045}$|$0.1766$|$\textcolor{red}{0.0778}$|$0.1962$|$0.0422$|
| Round 10 |$0.0607$|$0.0381$|$\textcolor{red}{0.1531}$|$0.0308$|$\textcolor{blue}{0.3075}$|$0.1057$|$0.0449$|


Recall here that we propose in this real data application to post-process all the outputs by first discarding those implementation outputs whose expected posterior loss is significantly higher or lower than other implementations.
The blue colors shown in the two tables above are the outputs we propose to no longer consider.
We can observe that all those blue rounds provide the minimized expected posterior loss which are significantly different from other implementations.
Moreover, if we focus on the same fixed $K$ cases, those ICPCL values in blue are all comparable or worse than the highest values of other rounds which we propose to consider.
By discarding those blue rounds, we still have $\geq 5$ rounds left for each fixed $K$ case.
The red colors correspond to those rounds we consider which provide the highest ICPCL values within the same fixed $K$ rounds.
After the post-processing, the boxplots of the ICPCL values shown as Figure $11$ can be recovered by the code shown below.

``` r
par(mfrow=c(1,1),mai = c(0.25, 0.3, 0.25, 0.1),mgp=c(1,0.25,0))
boxplot(c(-6396.096, -6454.013, -6450.135, -6446.579, -6373.261, -6386.550, -6437.208, -6428.155, -6428.155),
        c(-6431.902, -6380.615, -6348.968, -6377.889, -6401.935),
        c(-6414.270, -6304.791, -6459.649, -6291.410, -6238.031),
        c(-6294.696, -6401.470, -6500.255, -6087.947, -6244.250, -6298.549, -6433.268, -6294.832),
        c(-6400.087, -6889.146, -7031.210, -6549.627, -6295.793, -6798.479, -6695.855, -6701.610, -6259.197),
        c(-6585.189, -6583.900, -6469.524, -6710.810, -6290.709, -6517.335, -6499.246, -6679.074),
        c(-7724.104, -6853.329, -6730.280, -7823.358, -7085.496, -7340.897, -7248.848, -6996.099, -7007.193, -6570.362),xlab = "",ylab = "", main = "",cex.axis = 0.6)
title(xlab = "",ylab = "ICPCL", main = "ICPCL Boxplot for each K = 2,3,...,8", mgp=c(0.9,1,0),cex.main=0.8,cex.lab = 0.6)
axis(1, at=c(1:7), labels = c("K2","K3","K4","K5","K6","K7","K8"),cex.axis=0.75)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42),mgp=c(3,1,0))
```

Compared to other rounds of $K=5$ case, our best pick, $K=5$ round $4$ case, seems to be a "lucky" round. 
However, based on the irreducibility property of the MCMC sampler, there is no doubt that any implementations would finally reach and to stay around the converged clustering states returned by the $K=5$ round $4$ case if we are able to implement the algorithm for long enough, because the $K=5$ round $4$ case returns the ICPCL value which is significantly better than all other cases.
While, in practice, we never know how many iterations are enough for the posterior chains to reach such stationary distribution or whether the posterior chains would stuck around some local posterior mode.
The mixing problem or bad initial state problem and so on might make the posterior chains take very long time to reach the stationary distribution which is believed to be around the global posterior mode.
These motivate us to apply the multiple implementations in this real data application.

The results show that the multiple implementations successfully help us find the $K=5$ round $4$ case whose summarized clustering almost perfectly recovers the affiliation of the real UKfaculty network and also reveals some surther underlying clustering structure which is not observed.
Moreover, our ZINB-SBM model also successfully provides good inference and interpretation of the missing zero structure of such a real network as discussed in Section $5$ of the paper.

By treating the $K=5$ round $4$ case as the best case in our experiments, we apply further inference conditional on the summarized clustering as:

``` r
# # Further inference of R,Q,Pi conditional on the summarized z
# RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s <-
#   Directed_ZINBSBM_PCMwG_FixedZ(Y = UKfaculty_adj,
#                                 K = 5, T = 20000, eps_R = 0.175,beta1 = 20, beta2 = 180,
#                                 Z_0 = RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ)
```

The further $20,000$-iteration inference outputs are already loaded at the beginning of this section and we can summarize those clustering dependent parameters as well the the mean and variance of the Negative-Binomial distribution assumed for the edge weights by the following code.

``` r
## Summarize Pi
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedPi <-
  apply(array(unlist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Pi),
              dim = c(nrow(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Pi[[1]]),
                      ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Pi[[1]]),
                      length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Pi)))[,,10001:20001],1,mean)
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedPi

# ## Summarize R
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedR <-
  apply(array(unlist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$R),
              dim = c(nrow(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$R[[1]]),
                      ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$R[[1]]),
                      length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$R)))[,,10001:20001],1:2,mean)
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedR
# RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanR # Compare with posterior mean R
# RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_ICPCLandOptR$R # Compare with the optimized R used for ICPCL
# # acceptance rate for R in the further inference conditional on Z_s
# Reduce(`+`, RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Acceptance_count_R[10001:20001])/10001

## Summarize Q
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedQ <-
  apply(array(unlist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Q),
              dim = c(nrow(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Q[[1]]),
                      ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Q[[1]]),
                      length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Q)))[,,10001:20001],1:2,mean)
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedQ

# Summarized mean 
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedR*
  (1-RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedQ)/
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedQ
# Summarized Var
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedR*
  (1-RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedQ)/
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedQ^2
```

The posterior clustering rand index plot, the posterior density plot of $\boldsymbol{\Pi}|\tilde{\boldsymbol{Z}}$ as well as the posterior trace plot and posterior density plot of $p$ shown as Figure $12$ can be recovered by:

``` r
par(mfrow=c(2,2),mai = c(0.25, 0.25, 0.25, 0.05), mgp=c(1,0.25,0),cex.main=0.8,cex.lab = 0.8)
# Posterior clustering rand index plot
# require("fossil")
# RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_RI <- c()
# for (t in 1:20001){
#   RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_RI <-
#     c(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_RI,
#       rand.index(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[t]]%*%c(1:ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Z[[t]])),
#                  Real_data_application_UKfaculty_Z%*%c(1:4)))
# }
plot(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LSZ_RI,ylim = c(0.5,1),type = "l",xlab = "",ylab = "", main = "Posterior Clustering Rand Index",cex.axis = 0.8)

# Posterior Pi density plots
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredPi <-
  array(unlist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Pi),
        dim = c(nrow(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Pi[[1]]),
                ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Pi[[1]]),
                length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Pi)))[,,10001:20001]
plot(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredPi[1,],bw=0.008),col = 2, xlim=c(0,0.55),ylim=c(0,17), ylab = "",xlab="", main = TeX(r'(Posterior Density of $\pi_k$|\widetilde{\textbf{z}})'))
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredPi[2,],bw=0.008),col = 3)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredPi[3,],bw=0.008),col = 4)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredPi[4,],bw=0.008),col = 5)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredPi[5,],bw=0.008),col = 6)
legend("topright", legend=c(TeX(r'($\pi_1$)'),TeX(r'($\pi_2$)'),TeX(r'($\pi_3$)'),TeX(r'($\pi_4$)'),TeX(r'($\pi_5$)')),
       col=2:6, lty = 1, cex=0.7)

# Posterior p trace plot
plot(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$p[10001:20001], type = "l", main = "p Mixing without Burn-in",xlab = "",ylab = "", font.main = 1)

# Posterior p density plots
hist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4$p[10001:20001],ylab = "",xlab="", main = "Posterior Density of p", font.main = 1)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42),mgp=c(3,1,0))
```

## Extra Summary Statistics and Plots

Similar as we illustrated in simulation studies, we can also obtain some further summary statistics which are not included in the paper and compare them with the ones we provide in the paper.
Based on the posterior samples, we can obtain the posterior mean of $\boldsymbol{\Pi}, \boldsymbol{Q}$ and compare them with summarized $\tilde{\boldsymbol{\Pi}}, \tilde{\boldsymbol{Q}}$:

``` r
## Posterior mean Pi
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanPi <-
  apply(array(unlist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Pi),
              dim = c(nrow(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Pi[[1]]),
                      ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Pi[[1]]),
                      length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Pi)))[,,10001:20001],1,mean)
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanPi
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedPi # Compare with summarized Pi
## Posterior mean Q
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanQ <-
  apply(array(unlist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Q),
              dim = c(nrow(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Q[[1]]),
                      ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Q[[1]]),
                      length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_LS$Q)))[,,10001:20001],1:2,mean)
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanQ
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedQ # Compare with summarized Q
```

We can also obtain the posterior mean of $p|\tilde{\boldsymbol{Z}}$ and compare with the summarized $\tilde{p}$:

``` r
# p Posterior mean conditional on Z_s
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanp_CondZ_s <-
  mean(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$p[10001:20001])
hist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$p[10001:20001])
plot(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$p[10001:20001]))
plot(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$p[10001:20001], type = "l")
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanp_CondZ_s
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizedp # Compare with summarized p
```

Based on posterior mean of $\boldsymbol{R}$ and $\boldsymbol{Q}$, we can alternatively approximate the mean and variance of the Negative-Binomial distribution assumed for the edge weights shown as below.
If we compare with the mean and variance approximated by the summarized parameters shown in Section $5$ of the paper, we can observe that they are close to each other.

``` r
# Check distribution mean based on PosteriorMeanR and PosteriorMeanQ
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanR*
  (1-RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanQ)/
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanQ
# Check distribution var based on PosteriorMeanR and PosteriorMeanQ
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanR*
  (1-RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanQ)/
  RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanQ^2
```

The summary statistic $\boldsymbol{P_{m0}}$ can also be approximated by several ways for comparison and they are all close to each other as expected:

``` r
# Compare P_m0 evaluated by posterior mean
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizedp/
  (RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizedp+
     (1-RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Summarizedp)*
     dnbinom(0,RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanR,
             RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanQ))
# with P_m0 evaluated by posterior mean conditional on Z_s
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanp_CondZ_s/
  (RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanp_CondZ_s+
     (1-RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_PosteriorMeanp_CondZ_s)*
     dnbinom(0,RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedR,
             RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedQ))
# and with summarized P_m0
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedProbObs0Missing0
```

The plots of the concaveness of the log $\boldsymbol{Y}_{gh}$ term for the $K=5$ round $4$ case ICPCL can be checked by the code below and the performance is as expected as we discussed in Section $4.1$ of the paper.

``` r
# check concaveness of ICPCL for each log Y_gh term w.r.t. r_gh
res <- Directed_ZINBSBM_ApproximateICPCL_logY_gh_term_CheckConcave(Y = UKfaculty_adj,
                                                             ProbObs0Missing0 = RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedProbObs0Missing0,
                                                             Z = RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_SummarizedZ,
                                                             alpha=1, beta1=20,beta2=180, betaq1=1,betaq2=1)
par(mfrow=c(5,5),mai = c(0.15, 0.15, 0.15, 0.05),mgp = c(2,0.25,0))
plot(seq(0.001,20,0.001),res[1,1,],type = "l", main = "Block 1,1",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[1,2,],type = "l", main = "Block 1,2",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[1,3,],type = "l", main = "Block 1,3",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[1,4,],type = "l", main = "Block 1,4",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[1,5,],type = "l", main = "Block 1,5",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[2,1,],type = "l", main = "Block 2,1",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[2,2,],type = "l", main = "Block 2,2",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[2,3,],type = "l", main = "Block 2,3",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[2,4,],type = "l", main = "Block 2,4",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[2,5,],type = "l", main = "Block 2,5",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[3,1,],type = "l", main = "Block 3,1",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[3,2,],type = "l", main = "Block 3,2",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[3,3,],type = "l", main = "Block 3,3",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[3,4,],type = "l", main = "Block 3,4",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[3,5,],type = "l", main = "Block 3,5",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[4,1,],type = "l", main = "Block 4,1",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[4,2,],type = "l", main = "Block 4,2",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[4,3,],type = "l", main = "Block 4,3",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[4,4,],type = "l", main = "Block 4,4",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[4,5,],type = "l", main = "Block 4,5",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[5,1,],type = "l", main = "Block 5,1",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[5,2,],type = "l", main = "Block 5,2",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[5,3,],type = "l", main = "Block 5,3",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[5,4,],type = "l", main = "Block 5,4",xlab = "",ylab = "")
plot(seq(0.001,20,0.001),res[5,5,],type = "l", main = "Block 5,5",xlab = "",ylab = "")
```

The posterior density plot of each element of $\boldsymbol{R}$ and $\boldsymbol{Q}$ conditional on the summarized clustering from the further inference is similar to those illustrated in simulation study $1$ ZINB-SBM cases and can be checked by the plots one by one provided below.

``` r
# Posterior density of r_gh|Z_s
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR <-
  array(unlist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$R),
        dim = c(nrow(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$R[[1]]),
                ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$R[[1]]),
                length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$R)))[,,10001:20001]

plot(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[1,1,],bw=0.005),col = 1, xlim=c(-0.1,0.5),ylim=c(0,65), ylab = "",xlab="", main = "r_gh outputs density")
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[1,2,],bw=0.005),col = 2)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[1,3,],bw=0.005),col = 3)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[1,4,],bw=0.005),col = 4)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[1,5,],bw=0.005),col = 5)
legend("topright", legend=c("1,1","1,2","1,3", "1,4","1,5"),
       col=c(1:5), lty = 1, cex=0.4)

plot(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[2,1,],bw=0.005),col = 1, xlim=c(-0.1,0.5),ylim=c(0,65), ylab = "",xlab="", main = "r_gh outputs density")
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[2,2,],bw=0.005),col = 2)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[2,3,],bw=0.005),col = 3)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[2,4,],bw=0.005),col = 4)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[2,5,],bw=0.005),col = 5)
legend("topright", legend=c("2,1","2,2","2,3", "2,4","2,5"),
       col=c(1:5), lty = 1, cex=0.4)

plot(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[3,1,],bw=0.005),col = 1, xlim=c(-0.1,0.5),ylim=c(0,65), ylab = "",xlab="", main = "r_gh outputs density")
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[3,2,],bw=0.005),col = 2)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[3,3,],bw=0.005),col = 3)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[3,4,],bw=0.005),col = 4)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[3,5,],bw=0.005),col = 5)
legend("topright", legend=c("3,1","3,2","3,3", "3,4","3,5"),
       col=c(1:5), lty = 1, cex=0.4)

plot(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[4,1,],bw=0.005),col = 1, xlim=c(-0.1,0.5),ylim=c(0,20), ylab = "",xlab="", main = "r_gh outputs density")
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[4,2,],bw=0.005),col = 2)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[4,3,],bw=0.005),col = 3)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[4,4,],bw=0.005),col = 4)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[4,5,],bw=0.005),col = 5)
legend("topright", legend=c("4,1","4,2","4,3", "4,4","4,5"),
       col=c(1:5), lty = 1, cex=0.4)

plot(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[5,1,],bw=0.03),col = 1, xlim=c(-0.1,3),ylim=c(0,10), ylab = "",xlab="", main = "r_gh outputs density")
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[5,2,],bw=0.03),col = 2)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[5,3,],bw=0.03),col = 3)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[5,4,],bw=0.03),col = 4)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredR[5,5,],bw=0.03),col = 5)
legend("topright", legend=c("5,1","5,2","5,3", "5,4","5,5"),
       col=c(1:5), lty = 1, cex=0.4)

# Posterior density of q_gh|Z_s
RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ <-
  array(unlist(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Q),
        dim = c(nrow(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Q[[1]]),
                ncol(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Q[[1]]),
                length(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_Further20000InferCondZ_s$Q)))[,,10001:20001]

plot(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[1,1,]),col = 1, xlim=c(-0.1,1.1),ylim=c(0,20), ylab = "",xlab="", main = "q_gh outputs density")
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[1,2,]),col = 2)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[1,3,]),col = 3)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[1,4,]),col = 4)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[1,5,]),col = 5)
legend("topright", legend=c("1,1","1,2","1,3", "1,4","1,5"),
       col=c(1:5), lty = 1, cex=0.4)

plot(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[2,1,]),col = 1, xlim=c(-0.1,1.1),ylim=c(0,30), ylab = "",xlab="", main = "q_gh outputs density")
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[2,2,]),col = 2)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[2,3,]),col = 3)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[2,4,]),col = 4)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[2,5,]),col = 5)
legend("topright", legend=c("2,1","2,2","2,3", "2,4","2,5"),
       col=c(1:5), lty = 1, cex=0.4)

plot(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[3,1,]),col = 1, xlim=c(-0.1,1),ylim=c(0,35), ylab = "",xlab="", main = "q_gh outputs density")
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[3,2,]),col = 2)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[3,3,]),col = 3)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[3,4,]),col = 4)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[3,5,]),col = 5)
legend("topright", legend=c("3,1","3,2","3,3", "3,4","3,5"),
       col=c(1:5), lty = 1, cex=0.4)

plot(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[4,1,]),col = 1, xlim=c(-0.1,1),ylim=c(0,10), ylab = "",xlab="", main = "q_gh outputs density")
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[4,2,]),col = 2)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[4,3,]),col = 3)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[4,4,]),col = 4)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[4,5,]),col = 5)
legend("topright", legend=c("4,1","4,2","4,3", "4,4","4,5"),
       col=c(1:5), lty = 1, cex=0.4)

plot(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[5,1,]),col = 1, xlim=c(-0.1,1),ylim=c(0,15), ylab = "",xlab="", main = "q_gh outputs density")
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[5,2,]),col = 2)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[5,3,]),col = 3)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[5,4,]),col = 4)
lines(density(RDA_UKfaculty_ZINBSBM_Fixed_K5_Prior_p_Beta_20_180_T20000_4_InferredQ[5,5,]),col = 5)
legend("topright", legend=c("5,1","5,2","5,3", "5,4","5,5"),
       col=c(1:5), lty = 1, cex=0.4)
```

We finally end this real data application with providing the summarized clustering of the best round $7$ of $K=7$ cases shown below.

|  Affiliation \\ $\tilde{\boldsymbol{z}}$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
| :---: |:---: |:---: |:---: |:---: |:---: |:---: |:---: |
| **1** | 18 | 0 | 0 | 0 | 0 | 0 | 1 |
| **2** | 1  | 18 | 0 | 0 | 6 | 5 | 3 |
| **3** | 0 | 0 | 12 | 15 | 0 | 0 | 0 |
| **4** | 0 | 2 | 0 | 0 | 0 | 0 | 0 |

Such a clustering also provides a very nice recovery of the true affiliations as well as further splitting of the "Central group" and the "Surrounding group" for both affiliation $2$ and affiliation $3$.
It also includes the two underlying clusters (cluster no. $5$ and no. $7$ above) which are provided by our best $K=5$ round $4$ case's summarized clustering shown in the paper.
However, as we discussed at the end of the RDA Section $5$ of the paper, though such a clustering is also interpretable and is shown to provide a good clustering, the model complexity brings the ICPCL value being $-6290.709$ which is significantly worse than that of the $K=5$ round $4$ case which gives $-6087.947$.
Thus it's clear that the $K=7$ round $7$ clustering shown above is not the best pick, but it might be worthwhile to also be considered when analysing the UKfaculty real network.
