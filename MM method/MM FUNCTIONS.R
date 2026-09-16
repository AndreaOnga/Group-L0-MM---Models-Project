# Version 2 of the CV functions use warm starts
library(glmnet)


#### GAUSSIAN ####
mm_gaussian <- function(X, y, lambda, rho, tau3 = 1e-6, 
                        max_iter_out = 100, max_iter_in = 100, tol = 1e-6,
                        warm_start=NULL,intercept=T) {
  # intercept=T corresponds to an unpenalized intercept (both by L0 and ridge)
  X <- as.matrix(X)
  
  if(intercept){
    X=cbind(1,X)
  }
  
  y <- as.vector(y)
  n <- nrow(X)
  p <- ncol(X)
  
  # 1. Precompute the squared norms of the columns (scale handling)
  v <- colSums(X^2)
  
  # 2. "Warm Start" initialization using Ridge Regression
  # Starting exactly from 0 would trap the MM weights at extremely large values
  if(is.null(warm_start) | length(warm_start)!=p){
    if (p < n) {
      beta <- as.vector(solve(t(X) %*% X + lambda * diag(p)) %*% t(X) %*% y)
    } else {
      beta <- as.vector(t(X) %*% solve(X %*% t(X) + ifelse(lambda==0,1e-4,lambda) * diag(n)) %*% y) # beta transposed
    }
  }
  else{
    beta <- warm_start # this for sure contains an intercept, if intercept is equal to 1 and warm_start is given
  }
  
  if(any(is.na(beta))){ # this happens if X^t X is singular even if p<n
    beta <- as.vector(solve(t(X) %*% X + ifelse(lambda==0,1e-4,lambda) * diag(p)) %*% t(X) %*%  y)
  }
  
  
  # Initial global residual
  r <- y - as.vector(X %*% beta)
  
  # Constant in the denominator of the log-sum penalty
  log_denom <- log(1 + 1 / tau3)
  
  # Outer loop: Majorization-Minimization (MM) step
  for (iter_out in 1:max_iter_out) {
    beta_old_out <- beta
    
    # ADAPTIVE WEIGHT UPDATE (Linearization of the log-sum penalty)
    # phi_j = rho / ( log(1 + 1/tau3) * (|beta_j| + tau3) )
    phi <- rho / (log_denom * (abs(beta) + tau3))
    
    if(intercept){
      phi[1]=0 # unpenalized intercept
    }
    
    # Inner loop: Coordinate Descent (minimization of the convex surrogate)
    for (iter_in in 1:max_iter_in) {
      beta_old_in <- beta
      
      for (j in 1:p) {
        # Compute the local correlation using the global residual:
        # x_j^T * r_{-j} = x_j^T * (r + x_j * beta_j) = x_j^T * r + v_j * beta_j
        a_j <- sum(X[, j] * r) + v[j] * beta[j]
        
        # Adaptive Soft-Thresholding threshold
        threshold <- 0.5 * phi[j]
        
        # Soft-Thresholding operator with Ridge penalty (v_j + lambda)
        beta_new <- 0
        if (abs(a_j) > threshold) {
          beta_new <- sign(a_j) * (abs(a_j) - threshold) / (v[j] + lambda*ifelse(intercept & j==1,0,1))
        }
        
        # If the coefficient changes, dynamically update the global residual r
        if (beta_new != beta[j]) {
          r <- r - X[, j] * (beta_new - beta[j])
          beta[j] <- beta_new
        }
      }
      
      # Check inner convergence (Coordinate Descent)
      if (max(abs(beta - beta_old_in)) < tol) {
        break
      }
    }
    
    # Check global convergence (MM step)
    if (max(abs(beta - beta_old_out)) < tol) {
      break
    }
  }
  
  return(list(beta = beta, iterations = iter_out,intercept=intercept))
}

cv_mm_gaussian <- function(X, y, lambda_seq=0, rho_seq, folds = NULL,nfolds=5,
                           tau3 = 1e-6,max_iter_out = 100,max_iter_in = 100,tol = 1e-6,
                           best_fit=F,intercept=T) {
  # BEST-fit indicates whether the best model must be estimated or not
  n <- nrow(X)
  p <- ncol(X)
  
  if(intercept){
    p=p+1
  }
  
  par_seq=expand.grid(lambda_seq, rho_seq)
  
  n_pars=nrow(par_seq)
  
  # Balanced fold creation
  if( is.null(folds)){
    folds <- sample(rep(1:nfolds, length.out = n))
  }
  
  nfolds=length(unique(folds))
  
  # Error matrix
  cv_error_matrix <- matrix(0, nrow = nfolds, ncol = n_pars)
  
  cat("Avvio Cross-Validation (", nfolds, "-folds) per", n_pars, "combinazioni di parametri...\n")
  
  for (f in 1:nfolds) {
    cat("  Elaborazione Fold", f, "/", nfolds, "\n")
    
    # Train / validation split
    train_idx <- which(folds != f)
    val_idx   <- which(folds == f)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_val   <- X[val_idx, , drop = FALSE]
    y_val   <- y[val_idx]
    
    if(intercept){
      X_val=cbind(1,X_val)
    }
    
    for (r_idx in 1:n_pars) {
      fit_train <- mm_gaussian(
        X = X_train, y = y_train, 
        lambda = par_seq[r_idx,1], rho = par_seq[r_idx,2],
        tau3 = tau3,max_iter_out = max_iter_out,max_iter_in = max_iter_in,tol=tol,
        intercept=intercept
      )
      
      # Prediction on the validation set
      y_pred <- drop(X_val %*% fit_train$beta)
      
      # ERROR
      cv_error_matrix[f, r_idx] <- mean((y_val - y_pred)^2)
    }
  }
  
  # Compute mean performance metrics and standard errors
  mean_mse <- colMeans(cv_error_matrix)
  se_mse   <- apply(cv_error_matrix, 2, sd) / sqrt(nfolds)
  
  # Identify the optimal rho (minimum MSE rule)
  best_idx <- which.min(mean_mse)
  best_par  <- par_seq[best_idx,]
  
  if(best_fit){
    best_fit=mm_gaussian(
      X = X, y = y, 
      lambda =as.numeric(best_par[1]), rho = as.numeric(best_par[2]),intercept=intercept
    )
  }
  return(list(
    best_par = best_par,
    par_seq = par_seq,
    mean_mse = mean_mse,
    se_mse = se_mse,
    best_fit=best_fit,intercept=intercept
  ))
}


cv_mm_gaussian2 <- function(X, y, lambda_seq=0, rho_seq, folds = NULL,nfolds=5,
                            tau3 = 1e-6,max_iter_out = 100,max_iter_in = 100,tol = 1e-6,
                            best_fit=F,intercept=T) {
  # BEST-fit indicates whether the best model must be estimated or not
  n <- nrow(X)
  p <- ncol(X)
  
  if(intercept){
    p=p+1
  }
  
  par_seq=expand.grid(lambda_seq, rho_seq)
  
  n_pars=nrow(par_seq)
  
  # Balanced fold creation
  if( is.null(folds)){
    folds <- sample(rep(1:nfolds, length.out = n))
  }
  
  nfolds=length(unique(folds))
  
  # Error matrix
  cv_error_matrix <- matrix(0, nrow = nfolds, ncol = n_pars)
  
  cat("Avvio Cross-Validation (", nfolds, "-folds) per", n_pars, "combinazioni di parametri...\n")
  
  for (f in 1:nfolds) {
    cat("  Elaborazione Fold", f, "/", nfolds, "\n")
    
    # Train / validation split
    train_idx <- which(folds != f)
    val_idx   <- which(folds == f)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_val   <- X[val_idx, , drop = FALSE]
    y_val   <- y[val_idx]
    
    if(intercept){
      X_val=cbind(1,X_val)
    }
    
    beta=NULL
    
    for (r_idx in 1:n_pars) {
      fit_train <- mm_gaussian(
        X = X_train, y = y_train, 
        lambda = par_seq[r_idx,1], rho = par_seq[r_idx,2],
        tau3 = tau3,max_iter_out = max_iter_out,max_iter_in = max_iter_in,tol=tol,
        warm_start = beta,intercept=intercept
      )
      
      beta=fit_train$beta
      
      # Prediction on the validation set
      y_pred <- drop(X_val %*% fit_train$beta)
      
      # ERROR
      cv_error_matrix[f, r_idx] <- mean((y_val - y_pred)^2)
    }
  }
  
  # Compute mean performance metrics and standard errors
  mean_mse <- colMeans(cv_error_matrix)
  se_mse   <- apply(cv_error_matrix, 2, sd) / sqrt(nfolds)
  
  # Identify the optimal rho (minimum MSE rule)
  best_idx <- which.min(mean_mse)
  best_par  <- par_seq[best_idx,]
  
  if(best_fit){
    best_fit=mm_gaussian(
      X = X, y = y, 
      lambda =as.numeric(best_par[1]), rho = as.numeric(best_par[2]),
      tau3 = tau3,intercept = intercept
    )
  }
  return(list(
    best_par = best_par,
    par_seq = par_seq,
    mean_mse = mean_mse,
    se_mse = se_mse,
    best_fit=best_fit,intercept=intercept
  ))
}

#### GROUPED GAUSSIAN ####
group_update <- function(z, penalty, lambda_ridge) {
  nz <- sqrt(sum(z^2))
  
  if (nz <= penalty) {
    return(rep(0, length(z)))
  }
  
  shrink <- 1 - penalty / nz
  beta <- (shrink / (1 + lambda_ridge)) * z
  
  return(beta)
}

mm_grouped_gaussian <- function(X, y, groups, lambda, rho, tau3 = 1e-6, 
                                max_iter_out = 100, max_iter_in = 100,
                                tol = 1e-6, group_weights = NULL,
                                intercept=T) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  group_ids <- sort(unique(groups))
  G <- length(group_ids)
  if(intercept){
    p=p+1
    G=G+1
    group_ids=c(max(groups)+1,group_ids) # intercept is placed before all other groups
    groups=c(max(groups)+1,groups)
    X=cbind(1,X)
  }# new group for the unpenalized intercept
  
  
  # Initialization of group weights (default: group size)
  if (is.null(group_weights)) {
    group_weights <- numeric(G)
    for (g in seq_len(G)) {
      group_weights[g] <-(sum(groups == group_ids[g])) # or 1 if desired
      if(intercept & g==1){
        group_weights[g]=0 
      }
    }
  }
  
  if(length(group_weights)<G){# in case of an intercept, with user-specified group_weights
    group_weights=c(0,group_weights)
  }
  
  #### This removes the L0 penalty from the intercept, while keeping the other penalties
  
  
  Xorth <- vector("list", G)
  Rinv  <- vector("list", G) # Store the inverse for computational efficiency
  theta <- vector("list", G)
  
  # --- PHASE 1: PREPROCESSING AND ORTHOGONALIZATION ---
  for (g in seq_len(G)) {
    idx <- which(groups == group_ids[g])
    Xg <- X[, idx, drop = FALSE]
    
    XtX <- crossprod(Xg)
    
    # Efficiency check: is the matrix diagonal? (Categorical variables / dummies)
    # If the sum of the absolute values equals the trace, all off-diagonal elements are zero.
    is_diag <- (sum(abs(XtX)) - sum(abs(diag(XtX)))) < 1e-10
    
    if (is_diag && length(idx) > 1) {
      # FAST PATH: Avoid QR decomposition and only rescale by the norm.
      d <- sqrt(diag(XtX))
      d[d == 0] <- 1e-10 # Avoid division by zero
      
      Xorth[[g]] <- sweep(Xg, 2, d, FUN = "/")
      # The inverse of a diagonal matrix is the diagonal matrix of reciprocals
      Rinv[[g]]  <- diag(1 / d, nrow = length(d)) 
      
    } else {
      # STANDARD PATH: QR decomposition for correlated variables
      qrg <- qr(Xg)
      Xorth[[g]] <- qr.Q(qrg)
      Rg <- qr.R(qrg)
      Rinv[[g]] <- backsolve(Rg, diag(length(idx)))
      
      #Rinv[[g]] <- solve(Rg + tol * diag(ncol(Rg))) # more stable in theory
    }
    
  }
  
  theta_tmp=rep(NA,ncol(X))
  if(ncol(X) <=nrow(X)){
    X_tmpss=do.call(cbind,Xorth)
    
    theta_tmp=lm.fit(X_tmpss,y)$coeff
    
  }
  if(any(is.na(theta_tmp))){ # some coefficients might be NA
    X_tmpss=do.call(cbind,Xorth)
    if(p>n){
      theta_tmp=as.vector(t(X_tmpss) %*% solve(X_tmpss %*% t(X_tmpss) + ifelse(lambda==0,1e-4,lambda) * diag(n)) %*% y)
    }
    else{
      theta_tmp=as.vector( solve(t(X_tmpss) %*% X_tmpss+ ifelse(lambda==0,1e-4,lambda) * diag(p)) %*%t(X_tmpss) %*% y)
    }
    
  }# some coefficients might be NA
  
  for (g in seq_len(G)) {
    idx <- which(groups == group_ids[g])
    
    
    theta[[g]] <- theta_tmp[idx]
  }
  
  # Initialization of the residual
  r <- y
  for (g in seq_len(G)) {
    r <- r - Xorth[[g]] %*% theta[[g]]
  }
  
  log_denom <- log(1 + 1 / tau3)
  
  # --- PHASE 2: OPTIMIZATION (MM + BCD) ---
  for (iter_out in seq_len(max_iter_out)) {
    
    old_theta_out <- unlist(theta)
    phi <- numeric(G)
    
    #print(max(abs(old_theta_out)))
    
    # 1. MM step: update the weights based on the current estimate
    for (g in seq_len(G)) {
      # Compute the norm of beta_g in the ORIGINAL, non-orthogonalized space
      beta_g_curr <- drop(Rinv[[g]] %*% theta[[g]])
      norm_beta <- sqrt(sum(beta_g_curr^2))
      
      # Adaptive weight for group g
      phi[g] <- rho / (log_denom * (norm_beta + tau3))
    }
    
    # 2. BCD step: Block Coordinate Descent
    for (iter_in in seq_len(max_iter_in)) {
      max_change <- 0
      
      for (g in seq_len(G)) {
        Qg <- Xorth[[g]]
        old_theta_g <- theta[[g]]
        
        # Add the contribution of the current group to the residual
        r <- r + Qg %*% old_theta_g
        
        # Correlation vector
        z <- drop(crossprod(Qg, r))
        
        # Dynamic threshold: 1/2 * phi_g * group weight
        penalty_g <- 0.5 * phi[g] * group_weights[g]
        
        # Group update
        new_theta_g <- group_update(z=z,penalty= penalty_g,lambda= lambda*ifelse(intercept & g==1,0,1)) # intercept
        
        theta[[g]] <- new_theta_g
        
        # Remove the new contribution from the residual
        r <- r - Qg %*% new_theta_g
        
        max_change <- max(max_change, max(abs(new_theta_g - old_theta_g)))
      }
      
      # Early stopping of BCD if inner convergence is reached
      if (max_change < tol) break
    }
    
    # Early stopping of MM if global convergence is reached
    if (max(abs(unlist(theta) - old_theta_out)) < tol) break
  }
  
  # --- PHASE 3: POST-PROCESSING (Inverse Rotation) ---
  beta <- numeric(p)
  for (g in seq_len(G)) {
    idx <- which(groups == group_ids[g])
    beta[idx] <- drop(Rinv[[g]] %*% theta[[g]])
    if(max(abs(beta[idx]))>1e4) print(g)
  }
  
  fitted <- drop(X %*% beta)
  
  list(
    beta = beta,
    fitted = fitted,
    residuals = y - fitted,
    iterations_MM = iter_out,intercept=intercept
  )
}

mm_grouped_gaussian____ <- function(X, y,Xorth,Rinv,theta, groups, lambda, rho, tau3 = 1e-6, 
                                    max_iter_out = 100, max_iter_in = 100,
                                    tol = 1e-6,group_weights) {# DO NOT USE
  if (any(!is.finite(unlist(theta)))) {
    cat(
      "iter =", c(lambda,rho),
      " max_theta =", max(abs(unlist(theta))),
      " max_beta =", max(sapply(seq_len(G), function(g)
        max(abs(Rinv[[g]] %*% theta[[g]]))
      )),
      "\n"
    )
    stop("theta exploded")
  }
  
  n <- nrow(X)
  p <- ncol(X)
  G=length(unique(groups))
  group_ids <- unique(groups) # this way the intercept will be in front (if present) (since it is the
  # first value of groups, and the rest are sorted)
  
  # Initialization of the residual
  r <- y
  for (g in seq_len(G)) {
    r <- r - Xorth[[g]] %*% theta[[g]]
  }
  
  log_denom <- log(1 + 1 / tau3)
  
  # --- PHASE 2: OPTIMIZATION (MM + BCD) ---
  for (iter_out in seq_len(max_iter_out)) {
    
    old_theta_out <- unlist(theta)
    phi <- numeric(G)
    
    # 1. MM step: update the weights based on the current estimate
    for (g in seq_len(G)) {
      # Compute the norm of beta_g in the ORIGINAL, non-orthogonalized space
      beta_g_curr <- drop(Rinv[[g]] %*% theta[[g]])
      norm_beta <- sqrt(sum(beta_g_curr^2))
      
      # Adaptive weight for group g
      phi[g] <- rho / (log_denom * (norm_beta + tau3))
    }
    
    # 2. BCD step: Block Coordinate Descent
    for (iter_in in seq_len(max_iter_in)) {
      max_change <- 0
      
      for (g in seq_len(G)) {
        Qg <- Xorth[[g]]
        old_theta_g <- theta[[g]]
        
        # Add the contribution of the current group to the residual
        r <- r + Qg %*% old_theta_g
        
        # Correlation vector
        z <- drop(crossprod(Qg, r))
        
        # Dynamic threshold: 1/2 * phi_g * group weight
        penalty_g <- 0.5 * phi[g] * group_weights[g]
        
        # Group update
        new_theta_g <- group_update(z=z, penalty=penalty_g, lambda=lambda*ifelse(group_weights[g]==0 & g==1,0,1)) # intercept
        
        theta[[g]] <- new_theta_g
        
        # Remove the new contribution from the residual
        r <- r - Qg %*% new_theta_g
        
        max_change <- max(max_change, max(abs(new_theta_g - old_theta_g)))
      }
      
      # Early stopping of BCD if inner convergence is reached
      if (max_change < tol) break
    }
    
    # Early stopping of MM if global convergence is reached
    if (max(abs(unlist(theta) - old_theta_out)) < tol) break
  }
  
  # --- PHASE 3: POST-PROCESSING (Inverse Rotation) ---
  beta <- numeric(p)
  for (g in seq_len(G)) {
    idx <- which(groups == group_ids[g])
    beta[idx] <- drop(Rinv[[g]] %*% theta[[g]])
  }
  
  fitted <- drop(X %*% beta)
  
  list(
    beta = beta,
    fitted = fitted,
    residuals = y - fitted,
    iterations_MM = iter_out,
    theta=theta
  )
}

cv_mm_grouped_gaussian <- function(X, y,groups,group_weights = NULL, lambda_seq=0, rho_seq, folds = NULL,nfolds=5,
                                   tau3 = 1e-6,max_iter_out = 100,max_iter_in = 100,tol = 1e-6,
                                   best_fit=F,intercept=T) {
  # intercept=T means that the intercept is unpenalized by L0 (but penalized by ridge)
  n <- nrow(X)
  p=ncol(X)
  
  group_ids <- sort(unique(groups))
  
  G <- length(group_ids)
  
  if(intercept){
    p=p+1
    G=G+1
    group_ids=c(max(groups)+1,group_ids) # intercept is placed before all other groups
    groups=c(max(groups)+1,groups)
  }# new group for the unpenalized intercept
  
  
  # Initialization of group weights (default: group size)
  if (is.null(group_weights)) {
    group_weights <- numeric(G)
    for (g in seq_len(G)) {
      group_weights[g] <-(sum(groups == group_ids[g])) # or 1 if desired
      if(intercept & g==1){
        group_weights[g]=0 
      }
    }
  }
  
  if(length(group_weights)<G){# in case of an intercept, with user-specified group_weights, this removes penalization for the intercept
    group_weights=c(0,group_weights)
  }
  
  
  par_seq=expand.grid(lambda_seq, rho_seq)
  
  n_pars=nrow(par_seq)
  
  # Balanced fold creation
  if( is.null(folds)){
    folds <- sample(rep(1:nfolds, length.out = n))
  }
  
  nfolds=length(unique(folds))
  
  # Error matrix
  cv_error_matrix <- matrix(0, nrow = nfolds, ncol = n_pars)
  
  cat("Avvio Cross-Validation (", nfolds, "-folds) per", n_pars, "combinazioni di parametri...\n")
  
  for (f in 1:nfolds) {
    cat("  Elaborazione Fold", f, "/", nfolds, "\n")
    
    # Train / validation split
    train_idx <- which(folds != f)
    val_idx   <- which(folds == f)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_val   <- X[val_idx, , drop = FALSE]
    y_val   <- y[val_idx]
    
    if(intercept){
      X_train=cbind(1,X_train)
      X_val=cbind(1,X_val)
    }
    
    # Orthogonalize beforehand
    
    Xorth <- vector("list", G)
    Rinv  <- vector("list", G)
    theta=vector("list", G)
    
    # --- PHASE 1: PREPROCESSING AND ORTHOGONALIZATION ---
    for (g in seq_len(G)) {
      idx <- which(groups == group_ids[g])
      Xg <- X_train[, idx, drop = FALSE]
      
      XtX <- crossprod(Xg)
      
      # Efficiency check: is the matrix diagonal? (Categorical variables / dummies)
      # If the sum of the absolute values equals the trace, all off-diagonal elements are zero.
      is_diag <- (sum(abs(XtX)) - sum(abs(diag(XtX)))) < 1e-10
      
      if (is_diag && length(idx) > 1) {
        # FAST PATH: Avoid QR decomposition and only rescale by the norm.
        d <- sqrt(diag(XtX))
        d[d == 0] <- 1e-10 # Avoid division by zero
        
        Xorth[[g]] <- sweep(Xg, 2, d, FUN = "/")
        # The inverse of a diagonal matrix is the diagonal matrix of reciprocals
        Rinv[[g]]  <- diag(1 / d, nrow = length(d)) 
        
      } else {
        # STANDARD PATH: QR decomposition for correlated variables
        qrg <- qr(Xg)
        Xorth[[g]] <- qr.Q(qrg)
        Rg <- qr.R(qrg)
        Rinv[[g]] <- backsolve(Rg, diag(length(idx)))
      }
    }
    theta_tmp=rep(NA,ncol(X))
    X_tmpss=do.call(cbind,Xorth)
    if(ncol(X_tmpss) <=nrow(X_tmpss)){
      
      theta_tmp=lm.fit(X_tmpss,y_train)$coeff
      
    }
    
    if(any(is.na(theta_tmp))){ # some coefficients might be NA
      
      if(ncol(X_tmpss) >nrow(X_tmpss)){
        theta_tmp=as.vector(t(X_tmpss) %*% solve(X_tmpss %*% t(X_tmpss) + ifelse(min(lambda_seq)==0,1e-4,min(lambda_seq)) * 
                                                   diag(nrow(X_tmpss))) %*% y_train)
      }
      else{
        theta_tmp=as.vector( solve(t(X_tmpss) %*% X_tmpss+ ifelse(min(lambda_seq)==0,1e-4,min(lambda_seq)) * diag(p)) %*%t(X_tmpss) %*% y_train)
      }
    }# some coefficients might be NA
    
    for (g in seq_len(G)) {
      idx <- which(groups == group_ids[g])
      
      
      theta[[g]] <- theta_tmp[idx]
    }
    
    
    for (r_idx in 1:n_pars) {
      fit_train <- mm_grouped_gaussian____(
        X = X_train, y = y_train, Xorth=Xorth,Rinv = Rinv,theta,groups = groups,
        lambda = par_seq[r_idx,1], rho = par_seq[r_idx,2],
        tau3 = tau3,max_iter_out = max_iter_out,max_iter_in = max_iter_in,tol=tol,
        group_weights = group_weights
      )
      
      # Prediction on the validation set
      y_pred<- X_val %*% fit_train$beta
      
      # Warm start
      theta<- fit_train$theta
      
      # ERROR
      cv_error_matrix[f, r_idx] <- mean((y_val - y_pred)^2)
    }
  }
  
  # Compute mean performance metrics and standard errors
  mean_mse <- colMeans(cv_error_matrix)
  se_mse   <- apply(cv_error_matrix, 2, sd) / sqrt(nfolds)
  
  # Identify the optimal rho (minimum MSE rule)
  best_idx <- which.min(mean_mse)
  best_par  <- par_seq[best_idx,]
  
  if(best_fit){
    if(intercept){
      group_weights=group_weights[-1] # remove the first 0
      groups=groups[-1]# remove the intercept group
    }
    best_fit=mm_grouped_gaussian(
      X = X, y = y, groups = groups,group_weights = group_weights,
      lambda =as.numeric(best_par[1]), rho = as.numeric(best_par[2]),tau3 = tau3,
      intercept = intercept
    )
  }
  return(list(
    best_par = best_par,
    par_seq = par_seq,
    mean_mse = mean_mse,
    se_mse = se_mse,
    best_fit=best_fit,intercept=intercept
  ))
}
#### LOGIT ####
MISCerror=function(true,predicted){
  mean(true!=predicted)
}

mm_logit <- function(X, y, lambda, rho, tau3 = 1e-6, 
                     max_iter_out = 100, max_iter_in = 100, tol = 1e-6,
                     warm_start=NULL,intercept=T) {
  # intercept=T means that the intercept is not penalized by L0 (but it is penalized by ridge)
  X <- as.matrix(X)
  y <- as.vector(y)
  n <- nrow(X)
  p <- ncol(X)
  
  if(intercept){
    p=p+1
    X=cbind(1,X)
  }
  
  v <- colSums(X^2)
  
  # Controlled initialization (avoids exact zeros due to the log-sum penalty)
  if(is.null(warm_start) | length(warm_start)!=p){
    if (p < n) {
      beta <- glm.fit(X,y,family = binomial(link="logit"))$coef
    } else {
      beta <- as.matrix(glmnet(X,y,family = "binomial",alpha = 0,lambda = ifelse(lambda==0,1e-4,lambda),intercept = F)$beta)
    }
  }
  else{
    beta <- warm_start
  }
  
  if(any(is.na(beta))){# happens if X is singular with n<0
    beta <- as.matrix(glmnet(X,y,family = "binomial",alpha = 0,lambda = ifelse(lambda==0,1e-4,lambda),intercept = F)$beta)
  }
  
  log_denom <- log(1 + 1 / tau3)
  
  # OUTER LOOP: MM (Coupled Bohning approximation + Spike-and-Slab weights update)
  for (iter_out in 1:max_iter_out) {
    beta_old_out <- beta
    
    # 1. Compute current probabilities (with safety clipping for numerical stability)
    eta <- as.vector(X %*% beta)
    prob <- 1 / (1 + exp(-eta))
    prob <- pmin(pmax(prob, 1e-15), 1 - 1e-15)
    
    # 2. Compute Bohning pseudo-response
    z <- eta - 4 * (prob - y)
    
    # 3. Update Spike-and-Slab weights
    phi <- rho / (log_denom * (abs(beta) + tau3))
    
    if(intercept){
      phi[1]=0
    }
    
    # 4. INNER LOOP: Coordinate Descent on the pseudo-response
    # Residual based on the pseudo-response z
    r <- z - as.vector(X %*% beta)
    
    for (iter_in in 1:max_iter_in) {
      beta_old_in <- beta
      
      for (j in 1:p) {
        a_j <- sum(X[, j] * r) + v[j] * beta[j]
        
        # Rescaled threshold due to the 1/4 Bohning factor
        threshold <- 2 * phi[j]
        
        beta_new <- 0
        if (abs(a_j) > threshold) {
          # Rescaled denominator: v_j + 4*lambda
          beta_new <- sign(a_j) * (abs(a_j) - threshold) / (v[j] + 4 * lambda*ifelse(intercept & j==1,0,1))# intercept
        }
        
        if (beta_new != beta[j]) {
          r <- r - X[, j] * (beta_new - beta[j])
          beta[j] <- beta_new
        }
      }
      if (max(abs(beta - beta_old_in)) < tol) break
    }
    
    if (max(abs(beta - beta_old_out)) < tol) break
  }
  
  return(list(beta = beta, iterations = iter_out,intercept=intercept))
}

cv_mm_logit <- function(X, y, lambda_seq=0, rho_seq, folds = NULL,nfolds=5,
                        tau3 = 1e-6,max_iter_out = 100,max_iter_in = 100,tol = 1e-6,
                        best_fit=F,loss=MISCerror,intercept=T) {
  # loss is a function of true,predicted; beware of the coding
  n <- nrow(X)
  
  par_seq=expand.grid(lambda_seq, rho_seq)
  
  n_pars=nrow(par_seq)
  
  # Balanced fold creation
  if( is.null(folds)){
    folds <- sample(rep(1:nfolds, length.out = n))
  }
  
  nfolds=length(unique(folds))
  
  # Error matrix
  cv_error_matrix <- matrix(0, nrow = nfolds, ncol = n_pars)
  
  cat("Avvio Cross-Validation (", nfolds, "-folds) per", n_pars, "combinazioni di parametri...\n")
  
  for (f in 1:nfolds) {
    cat("  Elaborazione Fold", f, "/", nfolds, "\n")
    
    # Train / validation split
    train_idx <- which(folds != f)
    val_idx   <- which(folds == f)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_val   <- X[val_idx, , drop = FALSE]
    y_val   <- y[val_idx]
    
    if(intercept){
      X_val=cbind(1,X_val)
    }
    
    for (r_idx in 1:n_pars) {
      fit_train <- mm_logit(
        X = X_train, y = y_train, 
        lambda = par_seq[r_idx,1], rho = par_seq[r_idx,2],
        tau3 = tau3,max_iter_out = max_iter_out,max_iter_in = max_iter_in,tol=tol,
        intercept = intercept
      )
      
      # Prediction on the validation set
      eta=drop(X_val %*% fit_train$beta)
      y_pred <- exp(eta)/(1+exp(eta))
      
      # ERROR
      cv_error_matrix[f, r_idx] <- loss(true=y_val==1,predicted=y_pred>0.5)
    }
  }
  
  # Compute mean performance metrics and standard errors
  mean_err <- colMeans(cv_error_matrix)
  se_err   <- apply(cv_error_matrix, 2, sd) / sqrt(nfolds)
  
  # Identify the optimal rho (minimum error rule)
  best_idx <- which.min(mean_err)
  best_par  <- par_seq[best_idx,]
  
  if(best_fit){
    best_fit=mm_logit(
      X = X, y = y, 
      lambda =as.numeric(best_par[1]), rho = as.numeric(best_par[2]),tau3 = tau3,
      intercept = intercept
    )
  }
  return(list(
    best_par = best_par,
    par_seq = par_seq,
    mean_err = mean_err,
    se_err = se_err,
    best_fit=best_fit,intercept=intercept
  ))
}


cv_mm_logit2 <- function(X, y, lambda_seq=0, rho_seq, folds = NULL,nfolds=5,
                         tau3 = 1e-6,max_iter_out = 100,max_iter_in = 100,tol = 1e-6,
                         best_fit=F,loss=MISCerror,intercept=T) {
  # loss is a function of true,predicted; beware of the coding
  n <- nrow(X)
  
  par_seq=expand.grid(lambda_seq, rho_seq)
  
  n_pars=nrow(par_seq)
  
  # Balanced fold creation
  if( is.null(folds)){
    folds <- sample(rep(1:nfolds, length.out = n))
  }
  
  nfolds=length(unique(folds))
  
  # Error matrix
  cv_error_matrix <- matrix(0, nrow = nfolds, ncol = n_pars)
  
  cat("Avvio Cross-Validation (", nfolds, "-folds) per", n_pars, "combinazioni di parametri...\n")
  
  for (f in 1:nfolds) {
    cat("  Elaborazione Fold", f, "/", nfolds, "\n")
    
    # Train / validation split
    train_idx <- which(folds != f)
    val_idx   <- which(folds == f)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_val   <- X[val_idx, , drop = FALSE]
    y_val   <- y[val_idx]
    
    if(intercept){
      X_val=cbind(1,X_val)
    }
    
    beta <- NULL
    
    for (r_idx in 1:n_pars) {
      fit_train <- mm_logit(
        X = X_train, y = y_train, 
        lambda = par_seq[r_idx,1], rho = par_seq[r_idx,2],
        tau3 = tau3,max_iter_out = max_iter_out,max_iter_in = max_iter_in,tol=tol,
        warm_start = beta,intercept = intercept
      )
      
      beta=fit_train$beta
      
      # Prediction on the validation set
      eta=drop(X_val %*% fit_train$beta)
      y_pred <- exp(eta)/(1+exp(eta))
      
      # ERROR
      cv_error_matrix[f, r_idx] <- loss(true=y_val==1,predicted=y_pred>0.5)
    }
  }
  
  # Compute mean performance metrics and standard errors
  mean_err <- colMeans(cv_error_matrix)
  se_err   <- apply(cv_error_matrix, 2, sd) / sqrt(nfolds)
  
  # Identify the optimal rho (minimum error rule)
  best_idx <- which.min(mean_err)
  best_par  <- par_seq[best_idx,]
  
  if(best_fit){
    best_fit=mm_logit(
      X = X, y = y, 
      lambda =as.numeric(best_par[1]), rho = as.numeric(best_par[2]),tau3 = tau3,
      intercept = intercept
    )
  }
  return(list(
    best_par = best_par,
    par_seq = par_seq,
    mean_err = mean_err,
    se_err = se_err,
    best_fit=best_fit,intercept=intercept
  ))
}
#### GROUP LOGIT ####

mm_grouped_logit <- function(X, y, groups, lambda, rho, tau3 = 1e-6, 
                             max_iter_out = 100, max_iter_in = 100,
                             tol = 1e-6, group_weights = NULL,
                             intercept=T) {
  n <- nrow(X)
  p <- ncol(X)
  
  group_ids <- sort(unique(groups))
  G <- length(group_ids)
  
  if(intercept){
    p=p+1
    G=G+1
    group_ids=c(max(groups)+1,group_ids) #intercept is before all the others
    groups=c(max(groups)+1,groups)
    X=cbind(1,X)
  }#new group for the unpenalized intercept
  
  # Initialization of group weights (default: group size)
  if (is.null(group_weights)) {
    group_weights <- numeric(G)
    for (g in seq_len(G)) {
      group_weights[g] <-(sum(groups == group_ids[g])) #or 1 if you want
      if(intercept & g==1){
        group_weights[g]=0 
      }
    }
  }
  
  if(length(group_weights)<G){#in case of an intercept, with user's group_weights
    group_weights=c(0,group_weights)
  }
  
  Xorth <- vector("list", G)
  Rinv  <- vector("list", G)
  theta <- vector("list", G)
  
  z0 <- 4 * (y - 0.5)
  
  # --- PHASE 1: LOGISTIC PREPROCESSING (Performed only once!) ---
  for (g in seq_len(G)) {
    idx <- which(groups == group_ids[g])
    Xg <- X[, idx, drop = FALSE]
    XtX <- crossprod(Xg)
    
    # Efficiency check: if it is a categorical dummy variable, XtX is diagonal
    is_diag <- (sum(abs(XtX)) - sum(abs(diag(XtX)))) < 1e-10
    
    if (is_diag && length(idx) > 1) {
      # FAST PATH for categorical variables (No QR decomposition)
      d <- sqrt(diag(XtX))
      d[d == 0] <- 1e-10
      Xorth[[g]] <- sweep(Xg, 2, d, FUN = "/")
      Rinv[[g]]  <- diag(1 / d, nrow = length(d))
    } else {
      # STANDARD PATH (QR decomposition for correlated quantitative variables)
      qrg <- qr(Xg)
      Xorth[[g]] <- qr.Q(qrg)
      Rinv[[g]] <- backsolve(qr.R(qrg), diag(length(idx)))
    }
  }
  
  theta_tmp=rep(NA,ncol(X))
  if(ncol(X) <=nrow(X)){
    X_tmpss=do.call(cbind,Xorth)
    
    theta_tmp=glm.fit(X_tmpss,y,family=binomial(link="logit"))$coeff
    
  }
  if(any(is.na(theta_tmp))){#some might be na
    X_tmpss=do.call(cbind,Xorth)
    
    theta_tmp=as.matrix(glmnet(X_tmpss,y,family = "binomial",alpha = 0,
                               lambda = ifelse(lambda==0,1e-4,lambda),intercept = F)$beta)
    
  }
  
  for (g in seq_len(G)) {
    idx <- which(groups == group_ids[g])
    
    
    theta[[g]] <- theta_tmp[idx]
  }
  
  log_denom <- log(1 + 1 / tau3)
  
  # --- PHASE 2: TWO-LEVEL MM + BCD ---
  for (iter_out in seq_len(max_iter_out)) {
    old_theta_out <- unlist(theta)
    
    # 1. Compute the linear predictor in the transformed space
    eta <- numeric(n)
    for (g in seq_len(G)) {
      eta <- eta + as.vector(Xorth[[g]] %*% theta[[g]])
    }
    
    # 2. Logistic probabilities and Bohning pseudo-response
    prob <- 1 / (1 + exp(-eta))
    prob <- pmin(pmax(prob, 1e-15), 1 - 1e-15)
    z <- eta - 4 * (prob - y)
    
    # 3. Compute weights based on the norm in the ORIGINAL space
    phi <- numeric(G)
    for (g in seq_len(G)) {
      beta_g_curr <- drop(Rinv[[g]] %*% theta[[g]])
      norm_beta <- sqrt(sum(beta_g_curr^2))
      phi[g] <- rho / (log_denom * (norm_beta + tau3))
    }
    
    # 4. Inner BCD on the pseudo-response
    r <- z - eta 
    for (iter_in in seq_len(max_iter_in)) {
      max_change <- 0
      for (g in seq_len(G)) {
        Qg <- Xorth[[g]]
        old_theta_g <- theta[[g]]
        
        r <- r + Qg %*% old_theta_g
        zg <- drop(crossprod(Qg, r))
        
        # Rescaled threshold for Bohning
        penalty_g <- 2 * phi[g] * group_weights[g]
        nz <- sqrt(sum(zg^2))
        
        # Rescaled shrinkage denominator: 1 + 4*lambda
        if (nz <= penalty_g) {
          new_theta_g <- rep(0, length(zg))
        } else {
          new_theta_g <- (zg / (1 + 4 * lambda*ifelse(intercept & g==1,0,1))) * (1 - penalty_g / nz)#intercept
        }
        
        theta[[g]] <- new_theta_g
        r <- r - Qg %*% new_theta_g
        max_change <- max(max_change, max(abs(new_theta_g - old_theta_g)))
      }
      if (max_change < tol) break
    }
    
    if (max(abs(unlist(theta) - old_theta_out)) < tol) break
  }
  
  # --- PHASE 3: POST-PROCESSING (Inverse Rotation) ---
  beta <- numeric(p)
  for (g in seq_len(G)) {
    idx <- which(groups == group_ids[g])
    beta[idx] <- drop(Rinv[[g]] %*% theta[[g]])
  }
  
  return(list(beta = beta, iterations_MM = iter_out,intercept=intercept))
}

mm_grouped_logit____ <- function(X, y,Xorth,Rinv,theta, groups, lambda, rho, tau3 = 1e-6, 
                                 max_iter_out = 100, max_iter_in = 100,
                                 tol = 1e-6, group_weights = NULL) {#DO NOT USE
  n <- nrow(X)
  p <- ncol(X)
  group_ids <- unique(groups) #this way the intercept will be in front(if present) (since is the
  #first value of groups, and the rest are sorted)
  
  G=length(unique(groups))
  log_denom <- log(1 + 1 / tau3)
  
  # --- PHASE 2: TWO-LEVEL MM + BCD ---
  for (iter_out in seq_len(max_iter_out)) {
    old_theta_out <- unlist(theta)
    
    # 1. Compute the linear predictor in the transformed space
    eta <- numeric(n)
    for (g in seq_len(G)) {
      eta <- eta + as.vector(Xorth[[g]] %*% theta[[g]])
    }
    
    # 2. Logistic probabilities 
    prob <- 1 / (1 + exp(-eta))
    prob <- pmin(pmax(prob, 1e-15), 1 - 1e-15)
    z <- eta - 4 * (prob - y)
    
    # 3. Compute weights based on the norm in the ORIGINAL space
    phi <- numeric(G)
    for (g in seq_len(G)) {
      beta_g_curr <- drop(Rinv[[g]] %*% theta[[g]])
      norm_beta <- sqrt(sum(beta_g_curr^2))
      phi[g] <- rho / (log_denom * (norm_beta + tau3))
    }
    
    # 4. Inner BCD on the pseudo-response
    r <- z - eta 
    for (iter_in in seq_len(max_iter_in)) {
      max_change <- 0
      for (g in seq_len(G)) {
        Qg <- Xorth[[g]]
        old_theta_g <- theta[[g]]
        
        r <- r + Qg %*% old_theta_g
        zg <- drop(crossprod(Qg, r))
        #rescaling
        penalty_g <- 2 * phi[g] * group_weights[g]
        nz <- sqrt(sum(zg^2))
        
        # Rescaled shrinkage denominator: 1 + 4*lambda
        if (nz <= penalty_g) {
          new_theta_g <- rep(0, length(zg))
        } else {
          new_theta_g <- (zg / (1 + 4 * lambda*ifelse(group_weights[g]==0 & g==1,0,1))) * (1 - penalty_g / nz) #intercept
        }
        
        theta[[g]] <- new_theta_g
        r <- r - Qg %*% new_theta_g
        max_change <- max(max_change, max(abs(new_theta_g - old_theta_g)))
      }
      if (max_change < tol) break
    }
    
    if (max(abs(unlist(theta) - old_theta_out)) < tol) break
  }
  
  # --- PHASE 3: POST-PROCESSING (Inverse Rotation) ---
  beta <- numeric(p)
  for (g in seq_len(G)) {
    idx <- which(groups == group_ids[g])
    beta[idx] <- drop(Rinv[[g]] %*% theta[[g]])
  }
  
  return(list(beta = beta, iterations_MM = iter_out,theta=theta))
}

cv_mm_grouped_logit <- function(X, y,groups,group_weights = NULL, lambda_seq=0, rho_seq, folds = NULL,nfolds=5,
                                tau3 = 1e-6,max_iter_out = 100,max_iter_in = 100,tol = 1e-6,
                                best_fit=F,loss=MISCerror,intercept=T) {
  #intercept=T means that the intercept is unpenalized by L0 (but it is penalized by ridge)
  n <- nrow(X)
  p=ncol(X)
  
  group_ids <- sort(unique(groups))
  
  G <- length(group_ids)
  
  if(intercept){
    p=p+1
    G=G+1
    group_ids=c(min(groups)-1,group_ids) #intercept is before all the others
    groups=c(min(groups)-1,groups)
  }#new group for the unpenalized intercept
  
  
  # Initialization of group weights (default: square root of the group size)
  if (is.null(group_weights)) {
    group_weights <- numeric(G)
    for (g in seq_len(G)) {
      group_weights[g] <-(sum(groups == group_ids[g])) #or 1 if you want
      if(intercept & g==1){
        group_weights[g]=0 
      }
    }
  }
  
  if(length(group_weights)<G){#in case of an intercept, with user's group_weights, this removes penalization for intercept
    group_weights=c(0,group_weights)
  }
  
  
  
  par_seq=expand.grid(lambda_seq, rho_seq)
  
  n_pars=nrow(par_seq)
  
  # Balanced fold creation
  if( is.null(folds)){
    folds <- sample(rep(1:nfolds, length.out = n))
  }
  
  nfolds=length(unique(folds))
  
  # Error matrix
  cv_error_matrix <- matrix(0, nrow = nfolds, ncol = n_pars)
  
  cat("Avvio Cross-Validation (", nfolds, "-folds) per", n_pars, "combinazioni di parametri...\n")
  
  for (f in 1:nfolds) {
    cat("  Elaborazione Fold", f, "/", nfolds, "\n")
    
    # Train / Validation split
    train_idx <- which(folds != f)
    val_idx   <- which(folds == f)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_val   <- X[val_idx, , drop = FALSE]
    y_val   <- y[val_idx]
    
    if(intercept){
      X_train=cbind(1,X_train)
      X_val=cbind(1,X_val)
    }
    
    #orthogonalize beforehand
    
    Xorth <- vector("list", G)
    Rinv  <- vector("list", G)
    theta=vector("list", G)
    
    z0 <- 4 * (y_train - 0.5)
    
    # --- PHASE 1: PREPROCESSING AND ORTHOGONALIZATION ---
    for (g in seq_len(G)) {
      idx <- which(groups == group_ids[g])
      Xg <- X_train[, idx, drop = FALSE]
      
      XtX <- crossprod(Xg)
      
      # Efficiency check: is the matrix diagonal? (Categorical variables / dummies)
      # If the sum of the absolute values equals the trace, all off-diagonal elements are zero.
      is_diag <- (sum(abs(XtX)) - sum(abs(diag(XtX)))) < 1e-10
      
      if (is_diag && length(idx) > 1) {
        # FAST PATH: Avoid QR decomposition and only rescale by the norm.
        d <- sqrt(diag(XtX))
        d[d == 0] <- 1e-10 # Avoid division by zero
        
        Xorth[[g]] <- sweep(Xg, 2, d, FUN = "/")
        # The inverse of a diagonal matrix is the diagonal matrix of reciprocals
        Rinv[[g]]  <- diag(1 / d, nrow = length(d)) 
        
      }
      else {
        # STANDARD PATH: QR decomposition for correlated variables
        qrg <- qr(Xg)
        Xorth[[g]] <- qr.Q(qrg)
        Rg <- qr.R(qrg)
        Rinv[[g]] <- backsolve(Rg, diag(length(idx)))
      }
    }
    
    theta_tmp=rep(NA,ncol(X))
    if(ncol(X_train) <=nrow(X_train)){
      X_tmpss=do.call(cbind,Xorth)
      
      theta_tmp=glm.fit(X_tmpss,y_train,family=binomial(link="logit"))$coeff
      
    }
    if(any(is.na(theta_tmp))){
      X_tmpss=do.call(cbind,Xorth)
      
      theta_tmp=as.matrix(glmnet(x=X_tmpss,y=y_train,family = "binomial",alpha = 0,
                                 lambda = ifelse(min(lambda_seq)==0,1e-4,min(lambda_seq)),intercept = F)$beta)
      
    }#some might be na
    
    for (g in seq_len(G)) {
      idx <- which(groups == group_ids[g])
      
      
      theta[[g]] <- theta_tmp[idx]
    }
    
    #probabilmente serve farlo anche con ncol(X)>n con una ridge penalizzata poco
    
    
    r_idx=1
    
    fit_train <- mm_grouped_logit____(
      X = X_train, y = y_train, Xorth=Xorth,Rinv = Rinv,theta=theta ,groups = groups,
      lambda = par_seq[r_idx,1], rho = par_seq[r_idx,2],
      tau3 = tau3,max_iter_out = max_iter_out,max_iter_in = max_iter_in,tol=tol,
      group_weights = group_weights
    )
    
    # Prediction on the validation set
    eta=drop(X_val %*% fit_train$beta)
    y_pred <- exp(eta)/1+exp(eta)
    
    #warm_start
    
    theta<- fit_train$theta
    
    for (r_idx in 2:n_pars) {
      if(r_idx>n_pars){ next }#this happens if n_pars was equal to 1
      
      fit_train <- mm_grouped_logit____(
        X = X_train, y = y_train, Xorth=Xorth,Rinv = Rinv,theta=theta ,groups = groups,
        lambda = par_seq[r_idx,1], rho = par_seq[r_idx,2],
        tau3 = tau3,max_iter_out = max_iter_out,max_iter_in = max_iter_in,tol=tol,
        group_weights = group_weights
      )
      
      # Prediction on the validation set
      eta=drop(X_val %*% fit_train$beta)
      y_pred <- exp(eta)/(1+exp(eta))
      
      #warm_start
      theta<- fit_train$theta
      
      
      # ERROR
      cv_error_matrix[f, r_idx] <- loss(true=y_val==1,predicted=y_pred>0.5)
    }
  }
  
  # Compute mean performance metrics and standard errors
  mean_err <- colMeans(cv_error_matrix)
  se_err   <- apply(cv_error_matrix, 2, sd) / sqrt(nfolds)
  
  # Identify the optimal rho (minimum error rule)
  best_idx <- which.min(mean_err)
  best_par  <- par_seq[best_idx,]
  
  if(best_fit){
    if(intercept){
      group_weights=group_weights[-1] #remove the first 0 (of the intercept)
      groups=groups[-1]#remove the intercept one
    }
    
    best_fit=mm_grouped_logit(
      X = X, y = y, groups = groups,group_weights = group_weights,
      lambda =as.numeric(best_par[1]), rho = as.numeric(best_par[2]),tau3=tau3,
      intercept = intercept
    )
  }
  return(list(
    best_par = best_par,
    par_seq = par_seq,
    mean_err = mean_err,
    se_err = se_err,
    best_fit=best_fit,intercept=intercept
  ))
}