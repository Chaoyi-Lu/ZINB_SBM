#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# The partially collapsed Metropolis-Hastings within Gibbs algorithm for the Zero-inflated Negative-Binomial Stochastic Block Model
# This algorithm updates r_gh by M-H step and update q_gh by Gibbs step
# Missing data imputation applied
Directed_ZINBSBM_PCMwG <- function(Y,K,T, alpha=1, beta1=1,beta2=1, betaq1=1,betaq2=1, eps_R=0.1, Z_0=NA){
  # Y: the observed adjacency matrix
  # K: the fixed number of clusters
  # T: the total number of iterations to implement
  # alpha: prior parameter for Pi
  # beta1, beta2: prior parameters for p
  # betaq1, betaq2: prior parameters for q_gh
  # eps_R: proposal epsilon for R
  # Z_0: the initial clustering state if provided
  require(DirichletReg) # for rdirichlet
  N <- nrow(Y)
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initial states
  R_0 = matrix(0.5,K,K) # initial state of r_gh for each block g,h
  Q_0 = matrix(betaq1/(betaq1+betaq2),K,K) # initial state of q for each block g,h
  Pi_0 <- t(matrix(rep(1/K,K))) # Membership probability
  p_0 <- beta1/(beta1+beta2) # set initial state of p as E(Prior(p))
  #--------------------------------------------------------------------------------------------------------------------------------------------
  if (is.na(sum(Z_0))){ # if no initial clustering state provided
    fitY <- kmeans(Y, centers = K) # use k-means output as initial state of Z
    Z_0 <- matrix(0,N,K)
    for (i in 1:K){
      Z_0[fitY$cluster==i,i] <- 1 # transform to matrix form
    }
  }else{
    Z_0 <- Z_0
  }
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Label switch Z_0
  Z_initial <- c(Z_0%*%c(1:K)) # store the initial state at Z_initial and transform from the form of latent clustering matrix to the form of latent clustering vector
  Z_0_new <- Z_initial # Z_0_new would be our final label-switched output for the initial clustering and is initialized here
  member_list <- list()
  for (k in 1:K){
    member_list[[k]] <- which(Z_initial==k) # store those nodes which are in the same cluster for each k
  }
  K_list <- 1:K # store the cluster labels 1,2,3,...,K
  flag <- 1 # denote the cluster label we currently assign
  count <- 1 # loop from the first node
  Z_0_new[member_list[[Z_initial[count]]]] <- flag # assign all the nodes which are in the same cluster as node 1 to the cluster 1
  K_list <- K_list[K_list!=Z_initial[count]] # delete the cluster label when it is already reassigned
  while(count < N){ # loop over all nodes
    count <- count + 1
    if (Z_initial[count] != Z_initial[count-1] & Z_initial[count]%in%K_list){ # if the the current node is not in the same cluster as the previous node and if the cluster of the current node has not relabeled
      flag <- flag + 1 # set the new label the cluster should be relabeled to
      Z_0_new[member_list[[Z_initial[count]]]] <- flag # relabel the cluster for all the nodes in such a cluster
      K_list <- K_list[K_list!=Z_initial[count]] # delete the cluster label in initial clustering which is already relabeled 
    }
  }
  new_membership <- matrix(0,N,K)
  for (k in 1:K){
    new_membership[Z_0_new==k,k] <- 1 # transform back to the matrix form of the latent clustering variable
  }
  Z_0 <- new_membership
  #--------------------------------------------------------------------------------------------------------------------------------------------
  nu_0 <- matrix(0,N,N) # generate the initial state of the augmented latent missing zero indicator variable \nu based on Z_0, R_0, Q_0, p_0
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (Y[i,j]==0){
        R_gh_0 <- Z_0[i,]%*%R_0%*%Z_0[j,]
        Q_gh_0 <- Z_0[i,]%*%Q_0%*%Z_0[j,]
        nu_0[i,j] <- rbinom(1,1,p_0/(p_0 + (1-p_0)*dnbinom(0,R_gh_0,Q_gh_0))) # Binomial(n=1) is equivalent to Bernoulli
      }
      if (Y[j,i]==0){
        R_hg_0 <- Z_0[j,]%*%R_0%*%Z_0[i,]
        Q_hg_0 <- Z_0[j,]%*%Q_0%*%Z_0[i,]
        nu_0[j,i] <- rbinom(1,1,p_0/(p_0 + (1-p_0)*dnbinom(0,R_hg_0,Q_hg_0)))
      }
    }
  }
  #--------------------------------------------------------------------------------------------------------------------------------------------
  X_0 <- Y # generate the initial state of the latent fully observed adjacency matrix X
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (nu_0[i,j] == 1){
        X_0[i,j] <- rnbinom(1,Z_0[i,]%*%R_0%*%Z_0[j,],Z_0[i,]%*%Q_0%*%Z_0[j,])
      }
      if (nu_0[j,i] == 1){
        X_0[j,i] <- rnbinom(1,Z_0[j,]%*%R_0%*%Z_0[i,],Z_0[j,]%*%Q_0%*%Z_0[i,])
      }
    }#end j
  }
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize the posterior chains
  Z_list <- list(Z_0)
  nu_list <- list(nu_0)
  R_list <- list(R_0)
  Q_list <- list(Q_0)
  p_list <- c(p_0)
  Pi_list <- list(Pi_0)
  X_list <- list(X_0)
  Acceptance_count_R_list <- list(matrix(0,K,K)) # store whether accept or not in each iteration for the Metropolis-Hastings step of each r_gh where g,h=1,2,..,K
  #--------------------------------------------------------------------------------------------------------------------------------------------
  for (t in 1:T){
    if ((t%%100) == 0){
      print(t) # monitor the implementation
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Sampling/full-conditional step to update nu
    nu_new <- nu_list[[t]]
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        if (Y[i,j]==0){ # \nu_{ij} is only updated when the corresponding y_{ij} = 0
          R_gh_t <- Z_list[[t]][i,]%*%R_list[[t]]%*%Z_list[[t]][j,] # extract R[z_i,z_j]
          Q_gh_t <- Z_list[[t]][i,]%*%Q_list[[t]]%*%Z_list[[t]][j,] # extract Q[z_i,z_j]
          P_nuij_0 <- dbinom(0,1,p_list[t]) * dnbinom(Y[i,j],R_gh_t,Q_gh_t) # calculate (1-p)*NB(0|r_gh,q_gh)
          P_nuij_1 <- dbinom(1,1,p_list[t]) # this is p
          nu_new[i,j] <- rbinom(1,1,P_nuij_1/(P_nuij_1 + P_nuij_0)) # simulate from Bernoulli(p/[p+(1-p)*NB(0|r_gh,q_gh)])
        }
        if (Y[j,i]==0){ # we focus on directed network so we need to implement the above steps for nu_ji
          R_hg_t <- Z_list[[t]][j,]%*%R_list[[t]]%*%Z_list[[t]][i,]
          Q_hg_t <- Z_list[[t]][j,]%*%Q_list[[t]]%*%Z_list[[t]][i,]
          P_nuji_0 <- dbinom(0,1,p_list[t]) * dnbinom(Y[j,i],R_hg_t,Q_hg_t)
          P_nuji_1 <- dbinom(1,1,p_list[t])
          nu_new[j,i] <- rbinom(1,1,P_nuji_1/(P_nuji_1 + P_nuji_0))
        }
      }
    }
    nu_list[[t+1]] <- nu_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update X
    X_new <- Y # recall here that X = Y + X_m
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        if (nu_new[i,j] == 1){ # impute missing zeros
          X_new[i,j] <- rnbinom(1,Z_list[[t]][i,]%*%R_list[[t]]%*%Z_list[[t]][j,],Z_list[[t]][i,]%*%Q_list[[t]]%*%Z_list[[t]][j,])
        }
        if (nu_new[j,i] == 1){
          X_new[j,i] <- rnbinom(1,Z_list[[t]][j,]%*%R_list[[t]]%*%Z_list[[t]][i,],Z_list[[t]][j,]%*%Q_list[[t]]%*%Z_list[[t]][i,])
        }
      }#end j
    }
    X_list[[t+1]] <- X_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Full-conditional step to update Z
    Z_new <- Z_list[[t]] # update from previous state
    for (i in sample.int(N,N)){ # update all nodes in random order; note here that the experiments show that there is no significant difference if we update one by one from node 1 to node N
      P_Zi <- c() # note that this will be in log form
      for (k in 1:K){ # loop for all z_i = k
        P_Zi_k <- log(Pi_list[[t]][k]) # set multinational part of the full-conditional probability of z_i = k
        Z_temp <- Z_new # set a temporary copy of the current clustering state
        Z_temp[i,] <- (1:K==k)*1 # assign node i in Z_temp to cluster k
        for (j in 1:N){ # loop over j = 1,..,N where j != i
          if (j != i){
            P_Zi_k <- P_Zi_k + # calculate the Negative-Binomial part of the full-conditional probability of z_i = k
              dnbinom(X_new[i,j],Z_temp[i,]%*%R_list[[t]]%*%Z_temp[j,],Z_temp[i,]%*%Q_list[[t]]%*%Z_temp[j,],log = TRUE)+
              dnbinom(X_new[j,i],Z_temp[j,]%*%R_list[[t]]%*%Z_temp[i,],Z_temp[j,]%*%Q_list[[t]]%*%Z_temp[i,],log = TRUE)
          }
        }
        P_Zi <- c(P_Zi,P_Zi_k) # store the full-conditional probability of all k and note that this is in log form
      }
      cum_Zi <- c() # aim to calculate the cumulative sum of the exp(P_Zi) where P_Zi is in log form
      for (k1 in 1:K){ # a mathematical trick to calculate cum_Zi
        P_Zi_denominator <- 0
        for (k2 in 1:K){
          P_Zi_denominator <- P_Zi_denominator + exp(P_Zi[k2]-P_Zi[k1])
        }
        cum_Zi <- c(cum_Zi,1/P_Zi_denominator)
      }
      cum_Zi <- cumsum(cum_Zi) # obtain the cumulative sum
      u <- runif(1) 
      Z_new[i,] <- (1:K==which(cum_Zi>u)[1])*1 # determine the new membership
    }
    Z_list[[t+1]] <- Z_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Gibbs step to update Q
    Q_new <- Q_list[[t]]
    for (g in 1:K){
      for (h in 1:K){
        betaq1_new <- R_list[[t]][g,h]*(length(X_new[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h]) - 
                                          sum((Z_new%*%c(1:K)==g)*(Z_new%*%c(1:K)==h))) + betaq1
        betaq2_new <- sum(X_new[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h]) + betaq2
        Q_new[g,h] <- rbeta(1,betaq1_new,betaq2_new)
      }
    }
    Q_list[[t+1]] <- Q_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # MH step updates R
    R_new <- R_list[[t]]
    Acceptance_count_R_list[[t+1]] <- matrix(0,K,K)
    for (g in 1:K){
      for (h in 1:k){
        # Propose the next state of r_gh
        R_new[g,h] <- runif(1, max(0,R_list[[t]][g,h]-eps_R),R_list[[t]][g,h]+eps_R) # bounded uniform proposal
        # Then we calculate the min(log(1),alpha_ratio_right_R) which is in log form
        alpha_ratio_right_R <- sum(dnbinom(X_new[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h][diag(1,N,N)[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h]==0],R_new[g,h],Q_new[g,h],log = TRUE))-
          sum(dnbinom(X_new[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h][diag(1,N,N)[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h]==0],R_list[[t]][g,h],Q_new[g,h],log = TRUE)) +
          log(R_list[[t]][g,h]+eps_R - max(0,R_list[[t]][g,h]-eps_R)) - log(R_new[g,h]+eps_R - max(0,R_new[g,h]-eps_R)) # log bounded uniform proposal ratio
        alpha_ratio_R <- min(log(1),alpha_ratio_right_R)
        if (log(runif(1)) < alpha_ratio_R){ # accept or not
          Acceptance_count_R_list[[t+1]][g,h] <- Acceptance_count_R_list[[t+1]][g,h] + 1
        }else{
          R_new[g,h] <- R_list[[t]][g,h] # if not accept, stay at current state
        }
      }
    }
    R_list[[t+1]] <- R_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Gibbs step to update p
    p_new <- p_list[t]
    beta1_new <- sum(nu_new) + beta1
    beta2_new <- N*(N-1)- sum(nu_new) + beta2
    p_new <- rbeta(1,beta1_new,beta2_new)
    p_list[t+1] <- p_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Gibbs step to update Pi
    n_k <- c()
    for (k in 1:K){
      n_k <- c(n_k,sum((Z_new%*%c(1:K)==k)))
    }
    Pi_list[[t+1]] <- rdirichlet(1,n_k+alpha)
  }# end T
  return(list(Z = Z_list, nu = nu_list, R = R_list, Q = Q_list, Acceptance_count_R = Acceptance_count_R_list, p = p_list, Pi = Pi_list,X = X_list))
}
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# The partially collapsed Metropolis-Hastings within Gibbs algorithm for the Zero-inflated Negative-Binomial Stochastic Block Model
# Note here that this implementation fixes the clustering as Z_0
Directed_ZINBSBM_PCMwG_FixedZ <- function(Y,K,T, alpha=1, beta1=1,beta2=1, betaq1=1,betaq2=1, eps_R = 0.1, Z_0=NA){
  require(DirichletReg)
  N <- nrow(Y)
  
  R_0 = matrix(0.5,K,K)
  Q_0 = matrix(betaq1/(betaq1+betaq2),K,K)
  Pi_0 <- t(matrix(rep(1/K,K)))
  p_0 <- beta1/(beta1+beta2)
  #--------------------------------------------------------------------------------------------------------------------------------------------
  Z_initial <- c(Z_0%*%c(1:K))
  Z_0_new <- Z_initial
  member_list <- list()
  for (k in 1:K){
    member_list[[k]] <- which(Z_initial==k)
  }
  K_list <- 1:K
  flag <- 1
  count <- 1
  Z_0_new[member_list[[Z_initial[count]]]] <- flag
  K_list <- K_list[K_list!=Z_initial[count]]
  while(count < N){
    count <- count + 1
    if (Z_initial[count] != Z_initial[count-1] & Z_initial[count]%in%K_list){
      flag <- flag + 1
      Z_0_new[member_list[[Z_initial[count]]]] <- flag
      K_list <- K_list[K_list!=Z_initial[count]]
    }
  }
  new_membership <- matrix(0,N,K)
  for (k in 1:K){
    new_membership[Z_0_new==k,k] <- 1
  }
  Z_0 <- new_membership
  #--------------------------------------------------------------------------------------------------------------------------------------------
  nu_0 <- matrix(0,N,N)
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (Y[i,j]==0){
        R_gh_0 <- Z_0[i,]%*%R_0%*%Z_0[j,]
        Q_gh_0 <- Z_0[i,]%*%Q_0%*%Z_0[j,]
        nu_0[i,j] <- rbinom(1,1,p_0/(p_0 + (1-p_0)*dnbinom(0,R_gh_0,Q_gh_0))) # Binomial(n=1) is equivalent to Bernoulli
      }
      if (Y[j,i]==0){
        R_hg_0 <- Z_0[j,]%*%R_0%*%Z_0[i,]
        Q_hg_0 <- Z_0[j,]%*%Q_0%*%Z_0[i,]
        nu_0[j,i] <- rbinom(1,1,p_0/(p_0 + (1-p_0)*dnbinom(0,R_hg_0,Q_hg_0)))
      }
    }
  }
  #--------------------------------------------------------------------------------------------------------------------------------------------
  X_0 <- Y
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (nu_0[i,j] == 1){
        X_0[i,j] <- rnbinom(1,Z_0[i,]%*%R_0%*%Z_0[j,],Z_0[i,]%*%Q_0%*%Z_0[j,])
      }
      if (nu_0[j,i] == 1){
        X_0[j,i] <- rnbinom(1,Z_0[j,]%*%R_0%*%Z_0[i,],Z_0[j,]%*%Q_0%*%Z_0[i,])
      }
    }
  }
  #--------------------------------------------------------------------------------------------------------------------------------------------
  nu_list <- list(nu_0)
  R_list <- list(R_0)
  Q_list <- list(Q_0)
  p_list <- c(p_0)
  Pi_list <- list(Pi_0)
  Acceptance_count_R_list <- list(matrix(0,K,K))
  X_list <- list(X_0)
  Z <- Z_0 # fix Z
  #--------------------------------------------------------------------------------------------------------------------------------------------
  for (t in 1:T){
    if ((t%%100) == 0){
      print(t)
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    nu_new <- nu_list[[t]]
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        if (Y[i,j]==0){
          R_gh_t <- Z[i,]%*%R_list[[t]]%*%Z[j,]
          Q_gh_t <- Z[i,]%*%Q_list[[t]]%*%Z[j,]
          P_nuij_0 <- dbinom(0,1,p_list[t]) * dnbinom(Y[i,j],R_gh_t,Q_gh_t)
          P_nuij_1 <- dbinom(1,1,p_list[t])
          nu_new[i,j] <- rbinom(1,1,P_nuij_1/(P_nuij_1 + P_nuij_0))
        }
        if (Y[j,i]==0){
          R_hg_t <- Z[j,]%*%R_list[[t]]%*%Z[i,]
          Q_hg_t <- Z[j,]%*%Q_list[[t]]%*%Z[i,]
          P_nuji_0 <- dbinom(0,1,p_list[t]) * dnbinom(Y[j,i],R_hg_t,Q_hg_t)
          P_nuji_1 <- dbinom(1,1,p_list[t])
          nu_new[j,i] <- rbinom(1,1,P_nuji_1/(P_nuji_1 + P_nuji_0))
        }
      }
    }
    nu_list[[t+1]] <- nu_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    X_new <- Y
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        if (nu_new[i,j] == 1){
          X_new[i,j] <- rnbinom(1,Z[i,]%*%R_list[[t]]%*%Z[j,],Z[i,]%*%Q_list[[t]]%*%Z[j,])
        }
        if (nu_new[j,i] == 1){
          X_new[j,i] <- rnbinom(1,Z[j,]%*%R_list[[t]]%*%Z[i,],Z[j,]%*%Q_list[[t]]%*%Z[i,])
        }
      }
    }
    X_list[[t+1]] <- X_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    Z_new <- Z # fix Z
    #--------------------------------------------------------------------------------------------------------------------------------------------
    Q_new <- Q_list[[t]]
    for (g in 1:K){
      for (h in 1:K){
        betaq1_new <- R_list[[t]][g,h]*(length(X_new[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h]) - 
                                          sum((Z_new%*%c(1:K)==g)*(Z_new%*%c(1:K)==h))) + betaq1
        betaq2_new <- sum(X_new[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h]) + betaq2
        Q_new[g,h] <- rbeta(1,betaq1_new,betaq2_new)
      }
    }
    Q_list[[t+1]] <- Q_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    R_new <- R_list[[t]]
    Acceptance_count_R_list[[t+1]] <- matrix(0,K,K)
    for (g in 1:K){
      for (h in 1:k){
        R_new[g,h] <- runif(1, max(0,R_list[[t]][g,h]-eps_R),R_list[[t]][g,h]+eps_R)
        alpha_ratio_right_R <- sum(dnbinom(X_new[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h][diag(1,N,N)[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h]==0],R_new[g,h],Q_new[g,h],log = TRUE))-
          sum(dnbinom(X_new[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h][diag(1,N,N)[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h]==0],R_list[[t]][g,h],Q_new[g,h],log = TRUE)) +
          log(R_list[[t]][g,h]+eps_R - max(0,R_list[[t]][g,h]-eps_R)) - log(R_new[g,h]+eps_R - max(0,R_new[g,h]-eps_R))
        alpha_ratio_R <- min(log(1),alpha_ratio_right_R)
        if (log(runif(1)) < alpha_ratio_R){
          Acceptance_count_R_list[[t+1]][g,h] <- Acceptance_count_R_list[[t+1]][g,h] + 1
        }else{
          R_new[g,h] <- R_list[[t]][g,h]
        }
      }
    }
    R_list[[t+1]] <- R_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    p_new <- p_list[t]
    beta1_new <- sum(nu_new) + beta1
    beta2_new <- N*(N-1)- sum(nu_new) + beta2
    p_new <- rbeta(1,beta1_new,beta2_new)
    p_list[t+1] <- p_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    n_k <- c()
    for (k in 1:K){
      n_k <- c(n_k,sum((Z_new%*%c(1:K)==k)))
    }
    Pi_list[[t+1]] <- rdirichlet(1,n_k+alpha)
  }
  return(list(nu = nu_list, R = R_list, Q = Q_list, Acceptance_count_R = Acceptance_count_R_list, p = p_list, Pi = Pi_list,X = X_list))
}
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# The Gibbs sampler for the Zero-inflated Poisson Stochastic Block Model
# Missing data imputation applied
Directed_ZIPSBM_Gibbs <- function(Y,K,T, alpha=1, beta1=1,beta2=1, alpha1=1,alpha2=1, Z_0=NA){
  # Y: the observed adjacency matrix
  # K: the fixed number of clusters
  # T: the total number of iterations to implement
  # alpha: prior parameter for Pi
  # beta1, beta2: prior parameters for p
  # alpha1, alpha2: prior parameters for Lambda
  # Z_0: the initial clustering state if provided
  require(DirichletReg) # for rdirichlet
  N <- nrow(Y)
  # Initial states
  Lambda_0 <- matrix(alpha1/alpha2,K,K) # set all initial state of lambda_gh = E(Prior(Lambda))
  Pi_0 <- t(matrix(rep(1/K,K))) # membership probability
  p_0 <- beta1/(beta1+beta2) # set initial state of p as E(Prior(p))
  #--------------------------------------------------------------------------------------------------------------------------------------------
  if (is.na(sum(Z_0))){
    fitY <- kmeans(Y, centers = K)
    Z_0 <- matrix(0,N,K)
    for (i in 1:K){
      Z_0[fitY$cluster==i,i] <- 1
    }
  }else{
    Z_0 <- Z_0
  }
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Label switch Z_0
  Z_initial <- c(Z_0%*%c(1:K))
  Z_0_new <- Z_initial
  member_list <- list()
  for (k in 1:K){
    member_list[[k]] <- which(Z_initial==k)
  }
  K_list <- 1:K
  flag <- 1
  count <- 1
  Z_0_new[member_list[[Z_initial[count]]]] <- flag
  K_list <- K_list[K_list!=Z_initial[count]]
  while(count < N){
    count <- count + 1
    if (Z_initial[count] != Z_initial[count-1] & Z_initial[count]%in%K_list){
      flag <- flag + 1
      Z_0_new[member_list[[Z_initial[count]]]] <- flag
      K_list <- K_list[K_list!=Z_initial[count]]
    }
  }
  new_membership <- matrix(0,N,K)
  for (k in 1:K){
    new_membership[Z_0_new==k,k] <- 1
  }
  Z_0 <- new_membership
  #--------------------------------------------------------------------------------------------------------------------------------------------
  nu_0 <- matrix(0,N,N) # generate the initial state of the augmented latent missing zero indicator variable \nu based on Z_0, Lambda_0, p_0
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (Y[i,j]==0){
        Lambda_gh_0 <- Z_0[i,]%*%Lambda_0%*%Z_0[j,]
        nu_0[i,j] <- rbinom(1,1,p_0/(p_0 + (1-p_0)*dpois(0,Lambda_gh_0))) # Pois(0|r_gh) here if z_i=g, z_j=h
      }
      if (Y[j,i]==0){
        Lambda_hg_0 <- Z_0[j,]%*%Lambda_0%*%Z_0[i,]
        nu_0[j,i] <- rbinom(1,1,p_0/(p_0 + (1-p_0)*dpois(0,Lambda_hg_0))) # Pois(0|r_hg) here
      }
    }
  }
  #--------------------------------------------------------------------------------------------------------------------------------------------
  X_0 <- Y # generate the initial state of the latent fully observed adjacency matrix X
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (nu_0[i,j] == 1){
        X_0[i,j] <- rpois(1,Z_0[i,]%*%Lambda_0%*%Z_0[j,]) # sample from Poisson distribution here
      }
      if (nu_0[j,i] == 1){
        X_0[j,i] <- rpois(1,Z_0[j,]%*%Lambda_0%*%Z_0[i,])
      }
      
    }#end j
  }
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Initialize the posterior chains
  Z_list <- list(Z_0)
  nu_list <- list(nu_0)
  Lambda_list <- list(Lambda_0)
  p_list <- c(p_0)
  Pi_list <- list(Pi_0)
  X_list <- list(X_0)
  #--------------------------------------------------------------------------------------------------------------------------------------------
  for (t in 1:T){
    if ((t%%100) == 0){
      print(t)
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Sampling/full-conditional step updates nu
    nu_new <- nu_list[[t]]
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        if (Y[i,j]==0){
          nu_new[i,j] <- rbinom(1,1,p_list[t]/(p_list[t] + (1-p_list[t])*dpois(0,Z_list[[t]][i,]%*%Lambda_list[[t]]%*%Z_list[[t]][j,])))
        }
        if (Y[j,i]==0){
          nu_new[j,i] <- rbinom(1,1,p_list[t]/(p_list[t] + (1-p_list[t])*dpois(0,Z_list[[t]][j,]%*%Lambda_list[[t]]%*%Z_list[[t]][i,])))
        }
      }
    }
    nu_list[[t+1]] <- nu_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update X
    X_new <- Y
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        if (nu_new[i,j] == 1){
          X_new[i,j] <- rpois(1,Z_list[[t]][i,]%*%Lambda_list[[t]]%*%Z_list[[t]][j,])
        }
        if (nu_new[j,i] == 1){
          X_new[j,i] <- rpois(1,Z_list[[t]][j,]%*%Lambda_list[[t]]%*%Z_list[[t]][i,])
        }
      }
    }
    X_list[[t+1]] <- X_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Full-conditional step updates Z
    Z_new <- Z_list[[t]]
    for (i in sample.int(N,N)){
      P_Zi <- c()
      for (k in 1:K){
        P_Zi_k <- log(Pi_list[[t]][k])
        Z_temp <- Z_new
        Z_temp[i,] <- (1:K==k)*1
        for (j in 1:N){
          if (j != i){
            P_Zi_k <- P_Zi_k + # calculate the Poisson part of the full-conditional probability of z_i = k
              dpois(X_new[i,j],Z_temp[i,]%*%Lambda_list[[t]]%*%Z_temp[j,],log = TRUE)+
              dpois(X_new[j,i],Z_temp[j,]%*%Lambda_list[[t]]%*%Z_temp[i,],log = TRUE)
          }
        }
        P_Zi <- c(P_Zi,P_Zi_k)
      }
      cum_Zi <- c()
      for (k1 in 1:K){
        P_Zi_denominator <- 0
        for (k2 in 1:K){
          P_Zi_denominator <- P_Zi_denominator + exp(P_Zi[k2]-P_Zi[k1])
        }
        cum_Zi <- c(cum_Zi,1/P_Zi_denominator)
      }
      cum_Zi <- cumsum(cum_Zi)
      u <- runif(1)
      Z_new[i,] <- (1:K==which(cum_Zi>u)[1])*1
    }
    Z_list[[t+1]] <- Z_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update Lambda
    Lambda_new <- Lambda_list[[t]]
    for (g in 1:K){
      for (h in 1:K){
        alpha1_new <- sum(X_new[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h]) + alpha1
        alpha2_new <- length(X_new[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h]) -
          sum((Z_new%*%c(1:K)==g)*(Z_new%*%c(1:K)==h)) + alpha2 
        Lambda_new[g,h] <- rgamma(1,alpha1_new,alpha2_new)
      }
    }
    Lambda_list[[t+1]] <- Lambda_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update p
    p_new <- p_list[t]
    beta1_new <- sum(nu_new) + beta1
    beta2_new <- N*(N-1)- sum(nu_new) + beta2
    p_new <- rbeta(1,beta1_new,beta2_new)
    p_list[t+1] <- p_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Update Pi
    n_k <- c()
    for (k in 1:K){
      n_k <- c(n_k,sum((Z_new%*%c(1:K)==k)))
    }
    Pi_list[[t+1]] <- rdirichlet(1,n_k+alpha)
  }
  return(list(Z = Z_list, nu = nu_list, Lambda = Lambda_list, p = p_list, Pi = Pi_list,X = X_list))
}
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# The Gibbs sampler for the Zero-inflated Poisson Stochastic Block Model
# Note here that this implementation fixes the clustering as Z_0
Directed_ZIPSBM_Gibbs_FixedZ <- function(Y,K,T, alpha=1, beta1=1,beta2=1, alpha1=1,alpha2=1, Z_0){
  require(DirichletReg)
  N <- nrow(Y)
  # Initial state
  Lambda_0 <- matrix(alpha1/alpha2,K,K)
  Pi_0 <- t(matrix(rep(1/K,K)))
  p_0 <- beta1/(beta1+beta2)
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Label switch Z_0
  Z_initial <- c(Z_0%*%c(1:K))
  Z_0_new <- Z_initial
  member_list <- list()
  for (k in 1:K){
    member_list[[k]] <- which(Z_initial==k)
  }
  K_list <- 1:K
  flag <- 1
  count <- 1
  Z_0_new[member_list[[Z_initial[count]]]] <- flag
  K_list <- K_list[K_list!=Z_initial[count]]
  while(count < N){
    count <- count + 1
    if (Z_initial[count] != Z_initial[count-1] & Z_initial[count]%in%K_list){
      flag <- flag + 1
      Z_0_new[member_list[[Z_initial[count]]]] <- flag
      K_list <- K_list[K_list!=Z_initial[count]]
    }
  }
  new_membership <- matrix(0,N,K)
  for (k in 1:K){
    new_membership[Z_0_new==k,k] <- 1
  }
  Z_0 <- new_membership
  #--------------------------------------------------------------------------------------------------------------------------------------------
  nu_0 <- matrix(0,N,N)
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (Y[i,j]==0){
        Lambda_gh_0 <- Z_0[i,]%*%Lambda_0%*%Z_0[j,]
        nu_0[i,j] <- rbinom(1,1,p_0/(p_0 + (1-p_0)*dpois(0,Lambda_gh_0))) # Binomial(n=1) is equivalent to Bernoulli
      }
      if (Y[j,i]==0){
        Lambda_hg_0 <- Z_0[j,]%*%Lambda_0%*%Z_0[i,]
        nu_0[j,i] <- rbinom(1,1,p_0/(p_0 + (1-p_0)*dpois(0,Lambda_hg_0)))
      }
    }
  }
  #--------------------------------------------------------------------------------------------------------------------------------------------
  X_0 <- Y
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      if (nu_0[i,j] == 1){
        X_0[i,j] <- rpois(1,Z_0[i,]%*%Lambda_0%*%Z_0[j,])
      }
      if (nu_0[j,i] == 1){
        X_0[j,i] <- rpois(1,Z_0[j,]%*%Lambda_0%*%Z_0[i,])
      }
    }
  }
  #--------------------------------------------------------------------------------------------------------------------------------------------
  # Posterior chains
  nu_list <- list(nu_0)
  Lambda_list <- list(Lambda_0)
  p_list <- c(p_0)
  Pi_list <- list(Pi_0)
  X_list <- list(X_0)
  #--------------------------------------------------------------------------------------------------------------------------------------------
  for (t in 1:T){
    if ((t%%100) == 0){
      print(t)
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    nu_new <- nu_list[[t]]
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        if (Y[i,j]==0){
          nu_new[i,j] <- rbinom(1,1,p_list[t]/(p_list[t] + (1-p_list[t])*dpois(0,Z_0[i,]%*%Lambda_list[[t]]%*%Z_0[j,])))
        }
        if (Y[j,i]==0){
          nu_new[j,i] <- rbinom(1,1,p_list[t]/(p_list[t] + (1-p_list[t])*dpois(0,Z_0[j,]%*%Lambda_list[[t]]%*%Z_0[i,])))
        }
      }
    }
    nu_list[[t+1]] <- nu_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    X_new <- Y
    for (i in 1:(N-1)){
      for (j in (i+1):N){
        if (nu_new[i,j] == 1){
          X_new[i,j] <- rpois(1,Z_0[i,]%*%Lambda_list[[t]]%*%Z_0[j,])
        }
        
        if (nu_new[j,i] == 1){
          X_new[j,i] <- rpois(1,Z_0[j,]%*%Lambda_list[[t]]%*%Z_0[i,])
        }
      }
    }
    X_list[[t+1]] <- X_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    Z_new <- Z_0 # fix Z
    #--------------------------------------------------------------------------------------------------------------------------------------------
    Lambda_new <- Lambda_list[[t]]
    for (g in 1:K){
      for (h in 1:K){
        alpha1_new <- sum(X_new[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h]) + alpha1
        alpha2_new <- length(X_new[Z_new%*%c(1:K)==g,Z_new%*%c(1:K)==h]) -
          sum((Z_new%*%c(1:K)==g)*(Z_new%*%c(1:K)==h)) + alpha2 
        Lambda_new[g,h] <- rgamma(1,alpha1_new,alpha2_new)
      }
    }
    Lambda_list[[t+1]] <- Lambda_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    p_new <- p_list[t]
    beta1_new <- sum(nu_new) + beta1
    beta2_new <- N*(N-1)- sum(nu_new) + beta2
    p_new <- rbeta(1,beta1_new,beta2_new)
    p_list[t+1] <- p_new
    #--------------------------------------------------------------------------------------------------------------------------------------------
    n_k <- c()
    for (k in 1:K){
      n_k <- c(n_k,sum((Z_new%*%c(1:K)==k)))
    }
    Pi_list[[t+1]] <- rdirichlet(1,n_k+alpha)
  }
  return(list(nu = nu_list, Lambda = Lambda_list, p = p_list, Pi = Pi_list,X = X_list))
}
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# Simulation from ZINB-SBM
Simulation_Directed_ZINBSBM <- function(N,K,Pi,p,R,Q){
  # N: the number of nodes
  # K: the number of clusters
  # Pi: an 1 X K matrix which corresponds to the cluster assignment probability
  # p: missing zero probability
  # R,Q: two K X K matrices which are Negative-Binomial parameters
  
  Z <- t(rmultinom(N,1,Pi)) # simulate the clustering
  nu <- matrix(0, N, N)# sample the latent missing zero indicator variable nu
  for (i in 1:N){
    for (j in 1:N){
      if (i != j){
        nu[i,j] <- rbinom(1,1,p) # Binomial(n=1) is equivalent to Bernoulli
      }
    }
  }
  X <- matrix(0,N,N) # sample the latent fully observed adjacency matrix X
  res <- X # initialize the output observed Y
  for (i in 1:N){
    for (j in 1:N){
      if(i != j){
        X[i,j] <- rnbinom(1, Z[i,]%*%R%*%Z[j,], Z[i,]%*%Q%*%Z[j,])
        if (nu[i,j] == 0){
          res[i,j] <- X[i,j] # not missing zero
        }else{
          res[i,j] <- 0 # missing 0
        }
      }
    }
  }
  return(list(Z = Z, nu = nu, Y = res, X = X))
}

# Simulation from ZIP-SBM
Simulation_Directed_ZIPSBM <- function(N,K,Pi,p,Lambda){
  # Lambda: a K X K matrix which corresponds to the Poisson parameter
  Z <- t(rmultinom(N,1,Pi))
  nu <- matrix(0, N, N)
  for (i in 1:N){
    for (j in 1:N){
      if (i != j){
        nu[i,j] <- rbinom(1,1,p)
      }
    }
  }
  X <- matrix(0,N,N)
  res <- X
  for (i in 1:N){
    for (j in 1:N){
      if(i != j){
        X[i,j] <- rpois(1, Z[i,]%*%Lambda%*%Z[j,])
        if (nu[i,j] == 0){
          res[i,j] <- X[i,j]
        }else{
          res[i,j] <- 0
        }
      }
    }
  }
  return(list(Z = Z, nu = nu, Y = res, X = X))
}
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Label switching function for ZINB-SBM
LabelSwitching_SG2003_ZINBSBM <- 
  function(Z = NA, Pi  = NA, R  = NA, Q  = NA, Acceptance_count_R  = NA){
    # Z, Pi, R, Q are four posterior chains of the latent variable or parameters
    # Acceptance_count_R stores the number of acceptance of the Metropolis-Hastings updates for each r_{gh} where g,h = 1,2,...,K
    T <- length(Z) # extract the length of the chain
    N <- nrow(Z[[T]]) # extract the number of nodes
    K <- ncol(Z[[T]]) # extract the number of clusters
    for (t in 1:T){
      if (t%%100==0){
        print(t) # monitor the process
      }
      Z_initial <- c(Z[[t]]%*%c(1:K)) # store the initial clustering and transform the clustering membership matrix to the clustering membership vector
      Z_new <- Z_initial # Z_new would be our final label-switched output for the t'th clustering and is initialized here
      if (dim(Z[[t]])[2]!=1){# Check if there is only 1 cluster
        # Label switch Z first
        member_list <- list()
        for (k in 1:K){
          member_list[[k]] <- which(Z_initial==k) # Store those nodes which are in the same cluster for each k
        }
        K_list <- 1:K # store the cluster labels 1,2,3,...,K
        flag <- 1 # denote the cluster label we currently assign
        count <- 1 # loop from the first node
        Z_new[member_list[[Z_initial[count]]]] <- flag # assign all the nodes which in the same cluster as node 1 to the cluster 1
        K_list <- K_list[K_list!=Z_initial[count]] # delete the cluster label when it is already reassigned
        reassignment <- matrix(c(Z_initial[count],flag),1,2) # store the assignment rule as the rows of a matrix with 2 columns, e.g. cluster 3 being relabeled as cluster 1 is stored by c(3, 1); this variable is used for the label switching of the clustering dependent parameters
        while(count < N){ # loop over all nodes
          count <- count + 1
          if (Z_initial[count] != Z_initial[count-1] & Z_initial[count]%in%K_list){ # if the the current node is not in the same cluster as the previous node and if the cluster of the current node has not relabeled
            flag <- flag + 1 # set the new label the cluster should be relabeled to
            Z_new[member_list[[Z_initial[count]]]] <- flag # relabel the cluster for all the nodes in such a cluster
            K_list <- K_list[K_list!=Z_initial[count]] # delete the cluster label in initial clustering which is already relabeled 
            reassignment <- rbind(reassignment,c(Z_initial[count],flag)) # store the assignment rule
          }
        }
        new_membership <- matrix(0,N,K) # transform back to the matrix form of the clustering variable
        for (k in 1:K){
          new_membership[Z_new==k,k] <- 1
        }
        Z[[t]] <- new_membership
        
        if ( (!is.na(Pi[1])) & (!is.na(R[1])) & (!is.na(Q[1])) & (!is.na(Acceptance_count_R[1]))){ # if all other clustering dependent parameters are not NA
          # Label switch clustering dependent parameters
          if (nrow(reassignment)>1){ # there may be a case where all nodes in the same group so that the reassignment only have 1 row
            reassignment_del <- c() # in order to prevent that the exchanging happens (e.g. cluster label 2 --> cluster label 1, and then cluster label 1 --> cluster label 2) leading to repeat reassignment
            for (k1 in 1:(nrow(reassignment)-1)){ # check each pair of assignment rule
              for (k2 in (k1+1):nrow(reassignment)){
                if (reassignment[k1,1]==reassignment[k2,2] & reassignment[k1,2] == reassignment[k2,1]){
                  reassignment_del <- c(reassignment_del,k2)
                }
              }
            }
            if (length(reassignment_del)>0){
              reassignment <- matrix(reassignment[-reassignment_del,],ncol = 2) # delete replicated assignment rule
            }
          }
          
          for (k in 1:nrow(reassignment)){ # apply the label switching for the clustering dependent parameters
            if (reassignment[k,1]!=reassignment[k,2]){
              Pi[[t]][c(reassignment[k,1],reassignment[k,2])] <- Pi[[t]][c(reassignment[k,2],reassignment[k,1])]
              R[[t]][c(reassignment[k,1],reassignment[k,2]),] <- R[[t]][c(reassignment[k,2],reassignment[k,1]),]
              R[[t]][,c(reassignment[k,1],reassignment[k,2])] <- R[[t]][,c(reassignment[k,2],reassignment[k,1])]
              Q[[t]][c(reassignment[k,1],reassignment[k,2]),] <- Q[[t]][c(reassignment[k,2],reassignment[k,1]),]
              Q[[t]][,c(reassignment[k,1],reassignment[k,2])] <- Q[[t]][,c(reassignment[k,2],reassignment[k,1])]
              Acceptance_count_R[[t]][c(reassignment[k,1],reassignment[k,2]),] <- Acceptance_count_R[[t]][c(reassignment[k,2],reassignment[k,1]),]
              Acceptance_count_R[[t]][,c(reassignment[k,1],reassignment[k,2])] <- Acceptance_count_R[[t]][,c(reassignment[k,2],reassignment[k,1])]
              reassignment[,1][reassignment[,1]==reassignment[k,2]] <- reassignment[k,1] # update the newest labeling for all the rest relabel rule
            }
          }
        }
      }
    }# end t
    return(list(Z = Z, Pi = Pi, R = R, Q = Q,Acceptance_count_R = Acceptance_count_R))
}

# Label switching function for ZIP-SBM
LabelSwitching_SG2003_ZIPSBM <- function(Z = NA, Pi  = NA, Lambda  = NA){
  T <- length(Z)
  N <- nrow(Z[[T]])
  K <- ncol(Z[[T]])
  for (t in 1:T){
    if (t%%100==0){
      print(t)
    }
    Z_initial <- c(Z[[t]]%*%c(1:K))
    Z_new <- Z_initial
    if (dim(Z[[t]])[2]!=1){
      member_list <- list()
      for (k in 1:K){
        member_list[[k]] <- which(Z_initial==k)
      }
      K_list <- 1:K
      flag <- 1
      count <- 1
      Z_new[member_list[[Z_initial[count]]]] <- flag
      K_list <- K_list[K_list!=Z_initial[count]]
      reassignment <- matrix(c(Z_initial[count],flag),1,2)
      while(count < N){
        count <- count + 1
        if (Z_initial[count] != Z_initial[count-1] & Z_initial[count]%in%K_list){
          flag <- flag + 1
          Z_new[member_list[[Z_initial[count]]]] <- flag
          K_list <- K_list[K_list!=Z_initial[count]]
          reassignment <- rbind(reassignment,c(Z_initial[count],flag))
        }
      }
      new_membership <- matrix(0,N,K)
      for (k in 1:K){
        new_membership[Z_new==k,k] <- 1
      }
      Z[[t]] <- new_membership
      
      if (!is.na(Pi[1]) & !is.na(Lambda[1])){
        if (nrow(reassignment)>1){
          reassignment_del <- c()
          for (k1 in 1:(nrow(reassignment)-1)){
            for (k2 in (k1+1):nrow(reassignment)){
              if (reassignment[k1,1]==reassignment[k2,2] & reassignment[k1,2] == reassignment[k2,1]){
                reassignment_del <- c(reassignment_del,k2)
              }
            }
          }
          if (length(reassignment_del)>0){
            reassignment <- matrix(reassignment[-reassignment_del,],ncol = 2)
          }
        }
        for (k in 1:nrow(reassignment)){
          if (reassignment[k,1]!=reassignment[k,2]){
            Pi[[t]][c(reassignment[k,1],reassignment[k,2])] <- Pi[[t]][c(reassignment[k,2],reassignment[k,1])]
            Lambda[[t]][c(reassignment[k,1],reassignment[k,2]),] <- Lambda[[t]][c(reassignment[k,2],reassignment[k,1]),]
            Lambda[[t]][,c(reassignment[k,1],reassignment[k,2])] <- Lambda[[t]][,c(reassignment[k,2],reassignment[k,1])]
            reassignment[,1][reassignment[,1]==reassignment[k,2]] <- reassignment[k,1]
          }
        }
      }
    }
  }
  return(list(Z = Z, Pi = Pi, Lambda = Lambda))
}

# #--------------------------------------------------------------------------------------------------------------------------------------------
# #--------------------------------------------------------------------------------------------------------------------------------------------
# Evaluate ICPCL for ZINB-SBM and optimize the approximate ICPCL w.r.t. R
Directed_ZINBSBM_ApproximateICPCL_OptR <- function(Y, ProbObs0Missing0, Z, R_0 ,alpha=1, beta1=1,beta2=1, betaq1=1,betaq2=1){ # in log form
  # OptR means optimizing w.r.t. R
  # Y is the N X N observed adjacency matrix
  # ProbObs0Missing0 is the summarized P_{m0}, K X K matrix
  # Z is summarized Z which is in the form of N X K matrix
  # R_0 is the initial/current state of R and is a K X K matrix
  # others are prior parameters
  K = ncol(Z) # extract the number of clusters
  N = nrow(Y) # extract the number of nodes
  ProbObs0Missing0[is.na(ProbObs0Missing0)] <- 0 # it's possible that some element of the P_{m0} can be NA, so we transform NA to 0; NA occurs when the summarized clustering has empty clusters
  
  # Evaluate the approximate number of ones in each block of nu
  Ones_in_nu <- matrix(0,K,K)
  for (g in 1:K){
    for (h in 1:K){
      Ones_in_nu[g,h] <- round((length(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]==0])-
                                  sum((Z%*%c(1:K)==g)*(Z%*%c(1:K)==h)))*ProbObs0Missing0[g,h])
    }
  }
  
  ## Optimize the log Y_gh term w.r.t. r_gh
  # Evaluate the log Y_gh term w.r.t. the current r_gh
  log_Y_part_0 <- matrix(0,K,K) # stores the log Y_gh term for each pair of g,h = 1,2,...,K
  for (g in 1:K){
    for (h in 1:K){
      betaq1_0 <- R_0[g,h]*(length(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h])-sum((Z%*%c(1:K)==g)*(Z%*%c(1:K)==h))-Ones_in_nu[g,h]) + betaq1
      betaq2_0 <- sum(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]) + betaq2
      log_Y_part_0[g,h] <- sum(lgamma(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]+R_0[g,h])) -
        sum(log(factorial(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]))) -
        length(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0])*lgamma(R_0[g,h]) +
        lbeta(betaq1_0,betaq2_0) - lbeta(betaq1,betaq2)
    }
  }
  flag <- matrix(1,K,K) # indicate whether any of the blocks increases or not; flag == 1 means increasing and flag == 0 means no increasing
  eps <- 0.01 # initialize the searching range, i.e. checking r_gh + eps and r_gh - eps
  while(sum(flag) != 0){ # if any block increases
    R_new_plus <- R_0 + eps # check two directions of each r_gh
    R_new_minus <- R_0 - eps
    log_Y_part_new_plus <- log_Y_part_0 # initialize the log Y_gh term for the proposed r_gh+eps, r_gh-eps and r_gh
    log_Y_part_new_minus <- log_Y_part_0
    log_Y_part_new <- log_Y_part_0
    R_new <- R_0 # initialize the updated R
    for (g in 1:K){ # loop over each pair of clusters
      for (h in 1:K){
        if (flag[g,h]==1){ # check if the corresponding log Y_gh term increased or not in the previous update
          # Calculate +eps blocks
          betaq1_new_plus <- R_new_plus[g,h]*(length(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h])-sum((Z%*%c(1:K)==g)*(Z%*%c(1:K)==h))-Ones_in_nu[g,h]) + betaq1
          betaq2_new_plus <- sum(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]) + betaq2
          log_Y_part_new_plus[g,h] <- sum(lgamma(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]+R_new_plus[g,h])) -
            sum(log(factorial(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]))) -
            length(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0])*lgamma(R_new_plus[g,h]) +
            lbeta(betaq1_new_plus,betaq2_new_plus) - lbeta(betaq1,betaq2)
          
          # Calculate -eps blocks
          if (R_new_minus[g,h] < 0){R_new_minus[g,h]<-0.000001} # check whether r_gh-eps < 0; if so, assign a small enough positive value to r_gh-eps
          betaq1_new_minus <- R_new_minus[g,h]*(length(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h])-sum((Z%*%c(1:K)==g)*(Z%*%c(1:K)==h))-Ones_in_nu[g,h]) + betaq1
          betaq2_new_minus <- sum(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]) + betaq2
          log_Y_part_new_minus[g,h] <- sum(lgamma(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]+R_new_minus[g,h])) -
            sum(log(factorial(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]))) -
            length(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0])*lgamma(R_new_minus[g,h]) +
            lbeta(betaq1_new_minus,betaq2_new_minus) - lbeta(betaq1,betaq2)
          
          # Compare the log Y_gh term evaluated at r_gh+eps, r_gh-eps and current r_gh; pick the one which provides highest log Y_gh term as the new state of r_gh 
          log_Y_part_new[g,h] <- max(log_Y_part_new_plus[g,h],log_Y_part_new_minus[g,h],log_Y_part_0[g,h])
          if (log_Y_part_new[g,h] == log_Y_part_0[g,h]){ # first check "if" the log Y_gh term doesn't increase
            flag[g,h] <- 0
            R_new[g,h] <- R_0[g,h]
          }else if (log_Y_part_new[g,h] == log_Y_part_new_plus[g,h]){ # "if" not then check "else if" the log Y_gh term increases
            flag[g,h] <- 1
            R_new[g,h] <- R_new_plus[g,h]
          }else{ # "else if" not then check "else"
            flag[g,h] <- 1
            R_new[g,h] <- R_new_minus[g,h]
          }
        }#end if
      }
    }# end g
    if (sum(flag) == 0 & eps > 0.00001){ # decease the searching width eps until eps <= 0.00001
      flag <- matrix(1,K,K)
      eps <- eps/2
    }
    R_0 <- R_new # update R
    log_Y_part_0 <- log_Y_part_new # update log Y_gh term
  }
  res <- sum(log_Y_part_new) - 0.5*(K^2)*log(N*(N-1)) # add the penalized term
  
  # Evaluate the log nu part
  beta1_new <- sum(Ones_in_nu) + beta1
  beta2_new <- N*(N-1) - sum(Ones_in_nu) + beta2
  res <- res + lbeta(beta1_new,beta2_new) - lbeta(beta1,beta2)
  
  # Evaluate the log Z part
  n_k_list <- c()
  for (k in 1:K){
    n_k_list <- c(n_k_list,sum((Z%*%c(1:K)==k)))
  }
  res <- res + sum(lgamma(n_k_list+alpha))+lgamma(K*alpha)-lgamma(N+K*alpha)-K*lgamma(alpha)
  
  return(list(value = res,R = R_new)) # return the optimaized ICPCL and the corresponding R
}

# The following function aims to check the concaveness of the ICPCL w.r.t. r_gh
Directed_ZINBSBM_ApproximateICPCL_logY_gh_term_CheckConcave <- 
  function(Y, ProbObs0Missing0, Z, alpha=1, beta1=1,beta2=1, betaq1=1,betaq2=1, rUB = 10){ # in log form
    # OptR means optimizing w.r.t. R
    # Y is the N X N observed adjacency matrix
    # ProbObs0Missing0 is the summarized P_{m0}, K X K matrix
    # Z is summarized Z which is in the form of N X K matrix
    # rUP is the upper bound of the range of each r_gh
    # others are prior parameters
    K = ncol(Z)
    N = nrow(Y)
    ProbObs0Missing0[is.na(ProbObs0Missing0)] <- 0
    Ones_in_nu <- matrix(0,K,K)
    for (g in 1:K){
      for (h in 1:K){
        Ones_in_nu[g,h] <- round((length(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]==0])-
                                    sum((Z%*%c(1:K)==g)*(Z%*%c(1:K)==h)))*ProbObs0Missing0[g,h])
      }
    }
    # Evaluate the log Y_gh parts at the sequences of r_gh for g,h=1,2,...,K
    res <- array(0,dim = c(K,K,length(seq(0.001,rUB,0.001)))) # e.g. res[3,3,1] corresponds to the log Y_33 term evaluated at r_11 = 0.001
    for (g in 1:K){
      for (h in 1:K){
        count <- 0
        for (r in seq(0.001,rUB,0.001)){
          count <- count + 1
          res[g,h,count] <- sum(lgamma(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]+r)) -
            sum(log(factorial(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]))) -
            length(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0])*lgamma(r) +
            lbeta(r*(length(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h])-sum((Z%*%c(1:K)==g)*(Z%*%c(1:K)==h))-Ones_in_nu[g,h]) + betaq1,
                  sum(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]) + betaq2) - lbeta(betaq1,betaq2)
        }
      }
    }
    return(res)
  }

# #--------------------------------------------------------------------------------------------------------------------------------------------
# #--------------------------------------------------------------------------------------------------------------------------------------------
# Evaluate the Exact ICL for ZIP-SBM
Directed_ZIPSBM_ExactICL <- function(Y, ProbObs0Missing0, Z, alpha1,alpha2, beta1,beta2, alpha){ # in log form
  # Y is a N X N matrix
  # ProbObs0Missing0 is a K X K matrix
  # Z is a N X K matrix
  # others are prior parameters
  K = ncol(Z)
  N = nrow(Y)
  ProbObs0Missing0[is.na(ProbObs0Missing0)] <- 0
  Ones_in_nu <- matrix(0,K,K)
  for (g in 1:K){
    for (h in 1:K){
      Ones_in_nu[g,h] <- round((length(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]==0])-
                                  sum((Z%*%c(1:K)==g)*(Z%*%c(1:K)==h)))*ProbObs0Missing0[g,h])
    }
  }
  
  # Evaluate the log Z part
  n_k_list <- c()
  for (k in 1:K){
    n_k_list <- c(n_k_list,sum((Z%*%c(1:K)==k)))
  }
  res <- sum(lgamma(n_k_list+alpha))+lgamma(K*alpha)-lgamma(N+K*alpha)-K*lgamma(alpha)
  
  # Evaluate the log Y part
  for (g in 1:K){
    for (h in 1:K){
      if (sum(Z%*%c(1:K)==g) != 0 & sum(Z%*%c(1:K)==h) != 0){# if cluster g and h are not empty
        alpha1_new <- sum(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]) + alpha1
        alpha2_new <- length(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h])-sum((Z%*%c(1:K)==g)*(Z%*%c(1:K)==h))-Ones_in_nu[g,h] + alpha2 
        res <- res - sum(log(factorial(Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h][Y[Z%*%c(1:K)==g,Z%*%c(1:K)==h]!=0]))) +
          alpha1*log(alpha2) - lgamma(alpha1) + lgamma(alpha1_new) - alpha1_new*log(alpha2_new)
      }
    }
  }
  
  # Evaluate the log nu part
  beta1_new <- sum(Ones_in_nu) + beta1
  beta2_new <- N*(N-1) - sum(Ones_in_nu) + beta2
  res <- res + lbeta(beta1_new,beta2_new) - lbeta(beta1,beta2)
  
  return(res)
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Load the Library required
library(igraphdata) # for UKfaculty data
library(igraph) # igraph class for the networks
library(latex2exp) # for writing LaTeX on the plots
library(fossil) # for evaluating rand index
library(GreedyEPL) # for summarizing the clustering chains so that the expected posterior loss is minimized
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Load simulation study 1 artificial dataset and the corresponding reference latent variables and parameters
SS1_ZINBSBM_N75_K3 <- 
  list(Y = as.matrix(read.csv("SS1_ZINBSBM_N75K3_obsY.csv",header = TRUE)),
       X = as.matrix(read.csv("SS1_ZINBSBM_N75K3_obsX.csv",header = TRUE)),
       nu = as.matrix(read.csv("SS1_ZINBSBM_N75K3_obsnu.csv",header = TRUE)),
       Z = as.matrix(read.csv("SS1_ZINBSBM_N75K3_obsZ.csv",header = TRUE)))
colnames(SS1_ZINBSBM_N75_K3$Y) <- c()
colnames(SS1_ZINBSBM_N75_K3$X) <- c()
colnames(SS1_ZINBSBM_N75_K3$nu) <- c()
colnames(SS1_ZINBSBM_N75_K3$Z) <- c()
# Label switch Z
SS1_ZINBSBM_N75_K3_LSZ <- LabelSwitching_SG2003_ZINBSBM(Z = list(SS1_ZINBSBM_N75_K3$Z))$Z[[1]]
# Label switch initial R,Q,Pi and evaluate the initial mean and variance which are treated as reference statistics
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
SS1_ZINBSBM_N75_K3_obs_Initialmean <- SS1_ZINBSBM_N75_K3_obs_InitialR*(1-SS1_ZINBSBM_N75_K3_obs_InitialQ)/SS1_ZINBSBM_N75_K3_obs_InitialQ
SS1_ZINBSBM_N75_K3_obs_Initialvar <- SS1_ZINBSBM_N75_K3_obs_InitialR*(1-SS1_ZINBSBM_N75_K3_obs_InitialQ)/SS1_ZINBSBM_N75_K3_obs_InitialQ^2
# Set initial p
SS1_ZINBSBM_N75_K3_obs_Initialp <- 0.15
# Evaluate the initial P_m0
SS1_ZINBSBM_N75_K3_obs_InitialProbObs0Missing0 <- matrix(0,3,3)
for (k1 in 1:3){
  for (k2 in 1:3){
    SS1_ZINBSBM_N75_K3_obs_InitialProbObs0Missing0[k1,k2] <- SS1_ZINBSBM_N75_K3_obs_Initialp/(SS1_ZINBSBM_N75_K3_obs_Initialp + (1-SS1_ZINBSBM_N75_K3_obs_Initialp)*dnbinom(0,SS1_ZINBSBM_N75_K3_obs_InitialR[k1,k2],SS1_ZINBSBM_N75_K3_obs_InitialQ[k1,k2]))
  }
}

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Load simulation study 2 artificial dataset and the corresponding reference latent variables and parameters
SS2_ZIPSBM_N75_K3 <- 
  list(Y = as.matrix(read.csv("SS2_ZIPSBM_N75K3_obsY.csv",header = TRUE)),
       X = as.matrix(read.csv("SS2_ZIPSBM_N75K3_obsX.csv",header = TRUE)),
       nu = as.matrix(read.csv("SS2_ZIPSBM_N75K3_obsnu.csv",header = TRUE)),
       Z = as.matrix(read.csv("SS2_ZIPSBM_N75K3_obsZ.csv",header = TRUE)))
colnames(SS2_ZIPSBM_N75_K3$Y) <- c()
colnames(SS2_ZIPSBM_N75_K3$X) <- c()
colnames(SS2_ZIPSBM_N75_K3$nu) <- c()
colnames(SS2_ZIPSBM_N75_K3$Z) <- c()
# Label switch Z
SS2_ZIPSBM_N75_K3_LSZ <- LabelSwitching_SG2003_ZIPSBM(Z = list(SS2_ZIPSBM_N75_K3$Z))$Z[[1]]
# Label switch initial Lambda, Pi
res <- LabelSwitching_SG2003_ZIPSBM(Z = list(SS2_ZIPSBM_N75_K3$Z),
                                    Pi = list(c(0.3,0.4,0.3)),
                                    Lambda = list(matrix(c(2,0.1,0.5,
                                                           0.7,2.0,1.1,
                                                           0.9,0.3,2.5),3,3)))
SS2_ZIPSBM_N75_K3_obs_InitialLambda <- res$Lambda[[1]]
SS2_ZIPSBM_N75_K3_obs_InitialPi <- res$Pi[[1]]
# Set initial p
SS2_ZIPSBM_N75_K3_obs_Initialp <- 0.15
# Evaluate the initial P_m0
SS2_ZIPSBM_N75_K3_obs_InitialProbObs0Missing0 <- matrix(0,3,3)
for (k1 in 1:3){
  for (k2 in 1:3){
    SS2_ZIPSBM_N75_K3_obs_InitialProbObs0Missing0[k1,k2] <- SS2_ZIPSBM_N75_K3_obs_Initialp/(SS2_ZIPSBM_N75_K3_obs_Initialp + (1-SS2_ZIPSBM_N75_K3_obs_Initialp)*dpois(0,SS2_ZIPSBM_N75_K3_obs_InitialLambda[k1,k2]))
  }
}
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
# Load real data application "UKfaculty" dataset and the corresponding reference clustering
data("UKfaculty")
UKfaculty_adj <- as.matrix(get.adjacency(UKfaculty,attr = "weight"))
colnames(UKfaculty_adj) <- c()
rownames(UKfaculty_adj) <- c()
# Label switch reference clustering
Real_data_application_UKfaculty_Z <- matrix(0,length(V(UKfaculty)$Group),length(table(V(UKfaculty)$Group)))
for (k in 1:length(table(V(UKfaculty)$Group))){
  Real_data_application_UKfaculty_Z[V(UKfaculty)$Group==k,k] <- 1 # transform to matrix form
}
Real_data_application_UKfaculty_Z <- LabelSwitching_SG2003_ZINBSBM(Z = list(Real_data_application_UKfaculty_Z))$Z[[1]]
