.hminvolPortfolio <- function(Sigma, control = list()) {
  
  n <- dim(Sigma)[1]
  ctr <- .ctrPortfolio(n, control)
  
  # Create Hierarchy
  S <- S_mat(Sigma)
  P <- P_mat(S)
  
  Sigma <- S %*% Sigma %*% t(S)
  
  # Fetch constraints
  UB <- ctr$UB
  LB <- ctr$LB
  
  w <- rep(NA, nrow(S))
  
  # Loop through nodes
  for (i in 1:nrow(S)) {
    
    if (is.na(w[i])) {
      
      # Identify adjacent nodes
      idx <- apply(P, 1, function(x) all(x == P[i,]))
      
      cov <- Sigma[idx, idx]
      
      # Identify parent node
      parNode <- which(rowSums(as.matrix(S[, S[i,]==1])) == max(rowSums(as.matrix(S[, S[i,]==1]))))
      parNode <- parNode[rowSums(S)[parNode] > sum(S[i, ])]
      if (length(parNode)>0) parNode <- parNode[rowSums(S)[parNode] == min(rowSums(S)[parNode])]
      
      UBsub <- as.numeric(S %*% UB)[idx] / ifelse(length(parNode)==0, 1, x_rec(w, S)[parNode])
      LBsub <- as.numeric(S %*% LB)[idx] / ifelse(length(parNode)==0, 1, x_rec(w, S)[parNode])
      
      # Set up arguments for quadratic programmer
      Amat <- cbind(1, diag(sum(idx)), -diag(sum(idx)))
      bvec <- c(1, LBsub, -UBsub)
      
      # Check if constraints are binding
      if (round(sum(LBsub), 6) == 1) {
        
        w[idx] <- LBsub
        
      } else if (round(sum(UBsub), 6) == 1) {
        
        w[idx] <- UBsub
        
      } else {
        
        # Calculate optimal weights
        w[idx] <- quadprog::solve.QP(Dmat = as.matrix(cov), dvec = rep(0, sum(idx)), 
                                     Amat = Amat, bvec = bvec, meq = 1)$solution
        
      }
      
    }
    
  }
  
  # Collapse weights to asset mapping
  w.out <- c()
  for (i in 1:ncol(S)) w.out[i] <- prod(w[S[,i]==1])
  
  return(w.out)
  
}