.hrpPortfolio <- function(Sigma, control = list()) {
  
  n <- dim(Sigma)[1]
  ctr <- .ctrPortfolio(n, control)
  
  # Create Hierarchy
  S <- S_mat(Sigma)
  P <- P_mat(S)
  
  # Fetch constraints
  UB <- ctr$UB
  LB <- ctr$LB
  
  w <- rep(NA, nrow(S))
  
  # Loop through nodes
  for (i in 1:nrow(S)) {
    
    if (is.na(w[i])) {
      
      # Identify adjacent nodes
      idx <- apply(P, 1, function(x) all(x == P[i,]))
      
      # Identify parent node
      parNode <- which(rowSums(as.matrix(S[, S[i,]==1])) == max(rowSums(as.matrix(S[, S[i,]==1]))))
      parNode <- parNode[rowSums(S)[parNode] > sum(S[i, ])]
      if (length(parNode)>0) parNode <- parNode[rowSums(S)[parNode] == min(rowSums(S)[parNode])]
      
      # Calculate node-specific constraints
      UBsub <- as.numeric(S %*% UB)[idx] / ifelse(length(parNode)==0, 1, x_rec(w, S)[parNode])
      LBsub <- as.numeric(S %*% LB)[idx] / ifelse(length(parNode)==0, 1, x_rec(w, S)[parNode])
      
      # Check if constraints are binding and set weights equal constraints to save an optimisation step
      if (round(sum(LBsub), 6) == 1) {
        
        w[idx] <- LBsub
        
      } else if (round(sum(UBsub), 6) == 1) {
        
        w[idx] <- UBsub
        
      } else {
        
        # Calculate cluster variance
        alpha <- sapply(which(idx), function(x) {
          ivp <- 1 / diag(as.matrix(Sigma[S[x,] == 1, S[x,] == 1]))
          ivp <- ivp / sum(ivp)
          ivp <- ivp %*% Sigma[S[x,] == 1, S[x,] == 1] %*% ivp
        })
        
        # Note: LdP does not take the sqrt() of alpha! This is for comparison with invvol
        alpha <- 1 / sqrt(alpha) / sum(1 / sqrt(alpha))
        
        # Resolve constraint violations
        # The idea is simple: if constraint is violated, set to binding limit (i.e. make equal to constraint)
        # Subsequently distribute any excess allocation to the remaining weights in proportion of their inverse variance
        delta <- pmax(0, LBsub - alpha) - pmax(0, alpha - UBsub)
        maxit <- 1000
        niter <- 0
        while (any(abs(delta) > 0) && niter < maxit) {
          
          alpha[abs(delta) > 0] <- alpha[abs(delta) > 0] + delta[abs(delta) > 0]
          alpha[delta == 0] <- alpha[delta == 0] + (1 - sum(alpha)) * alpha[delta == 0] / sum(alpha[delta == 0])
          alpha <- alpha / sum(alpha)
          delta <- pmax(0, LBsub - alpha) - pmax(0, alpha - UBsub)
          niter <- niter + 1
          
        }
        
        w[idx] <- alpha
        
      }
      
    }
    
  }
  
  # Collapse weights to the asset dimension
  w.out <- c()
  for (i in 1:ncol(S)) w.out[i] <- prod(w[S[,i]==1])
  
  return(w.out)
  
}

S_mat <- function(Sigma, clust.method = "AGNES") {
  
  corr <- cov2cor(Sigma)
  
  distmat <- ((1 - corr) / 2)^0.5
  if (clust.method == "AGNES") {
    clust <- cluster::agnes(dist(distmat), method = "ward")
  }
  if (clust.method == "DIANA") {
    clust <- cluster::diana(dist(distmat))
  }
  
  order <- clust$order
  
  bisect.inner <- function(x.vec, type, corr) {
    
    split.pos <- floor(length(x.vec) / 2)
    
    x.vec.1 <- 1:split.pos
    x.vec.2 <- c(1:length(x.vec))[-x.vec.1]
    
    vec.1 <- vec.2 <- rep(0, ncol(Sigma))
    vec.1[x.vec[x.vec.1]] <- 1
    vec.2[x.vec[x.vec.2]] <- 1
    
    S.mat <<- cbind(S.mat, vec.1, vec.2)
    
    if (length(x.vec.1) > 1) bisect.inner(x.vec[x.vec.1], type, corr[x.vec.1, x.vec.1])
    if (length(x.vec.2) > 1) bisect.inner(x.vec[x.vec.2], type, corr[x.vec.2, x.vec.2])
    
  }
  
  S.mat <- c()
  
  bisect.inner(order, type, corr)
  
  return(t(S.mat))
  
}

P_mat <- function(S) {
  
  Sfull <- rbind(1, S)
  
  SparentSums <- sapply(1:nrow(S), function(x) rowSums(as.matrix(Sfull[, S[x,]==1])))
  Sparent <- S
  
  for (i in 1:nrow(S)) {
    
    idx <- SparentSums[i + 1, i] == SparentSums[, i]
    idx[i + 1] <- F
    
    Ssize <- rowSums(Sfull)
    Ssize[!idx] <- max(Ssize) + 1
    
    Sparent[i,] <- Sfull[which.min(Ssize),]
    
  }
  
  P <- matrix(0, nrow(S), nrow(S))
  for (i in 1:nrow(S)) {
    
    P[i,] <- apply(Sparent, 1, function(x) all(x == Sparent[i,]))
    
  }
  
  return(P)
  
}

x_rec <- function(x, S) {
  
  x.out <- x
  for (i in 1:nrow(S)) {
    
    temp <- rowSums(as.matrix(S[, S[i, ]==1]))
    idx <- temp >= temp[i]
    x.out[i] <- prod(x[idx])
    
  }
  
  return(x.out)
  
}