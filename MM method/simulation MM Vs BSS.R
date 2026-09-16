# ============================================================================
# SIMULATION STUDY: Best Subset Selection (BSS) vs an MM-algorithm surrogate
#
# This script compares an exhaustive best-subset-selection procedure (which
# minimizes an L0-penalized criterion, e.g. AIC/BIC-style) against a proposed
# MM (Majorize-Minimize) algorithm that approximates the same L0-penalized
# solution, across several settings:
#   1) Gaussian linear regression (ungrouped variables)
#   2) Logistic regression (ungrouped variables)
#   3) Gaussian linear regression with GROUPED variables
#   4) Logistic regression with GROUPED variables
#
# For each setting, B=1000 Monte Carlo replications are run, coefficients are
# estimated by both methods, and the results are compared (bias, sparsity
# pattern recovery, and the value of the penalized loss) and visualized.
# ============================================================================

source("MM FUNCTIONS.R")   # loads mm_gaussian, mm_logit, mm_grouped_gaussian, mm_grouped_logit, etc.
library(mvtnorm)                          # for rmvnorm(), used to generate correlated design matrices

# ----------------------------------------------------------------------------
# best_subset(): exhaustive best subset selection
#
# Tries every possible subset of (groups of) variables, up to size max_size,
# fits the model on each subset with fit_fun, evaluates a user-supplied
# criterion_fun (e.g. penalized RSS or penalized negative log-likelihood) on
# the FULL coefficient vector (zero-padded outside the subset), and returns
# the subset achieving the minimum criterion.
#
# Arguments:
#   X            - design matrix (n x p)
#   y            - response vector
#   groups       - optional grouping vector of length p, mapping each column
#                  of X to a group id (defaults to 1:ncol(X), i.e. no grouping)
#   fit_fun      - function(X_S, y) that fits a model on the subset of
#                  columns X_S and returns an object with $coefficients
#   criterion_fun- function(beta, X, y) returning a scalar loss to minimize,
#                  evaluated at the full-length coefficient vector beta
#   max_size     - maximum number of groups to include in a subset (capped at G)
# ----------------------------------------------------------------------------
best_subset <- function(X, y,groups=1:ncol(X), fit_fun, criterion_fun,
                        max_size = ncol(X)) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  group_ids <- sort(unique(groups))
  G <- length(group_ids)
  
  # Massimo numero di gruppi di variabili da includere in un sottoinsieme
  # (Maximum number of variable groups to include in a subset)
  max_size <- min(max_size, G)
  
  # Oggetti per salvare i risultati
  # (Containers to store results for every subset tried)
  results <- vector("list", 0)
  
  counter <- 1
  
  # ============================================================
  # ESAUSTIVO: provo tutti i sottoinsiemi
  # (EXHAUSTIVE SEARCH: try every subset of groups, size 1 to max_size)
  # ============================================================
  
  for (k in 1:max_size) {
    
    # All combinations of k groups out of G
    subsets <- combn(G, k, simplify = FALSE)
    
    for (S in subsets) {
      
      # Restrict the design matrix to columns belonging to the
      # groups in the current subset S
      X_S <- X[, groups %in% S, drop = FALSE]
      
      # --------------------------------------------------------
      # STIMA (FIT the model on the restricted design matrix)
      # --------------------------------------------------------
      fit <- fit_fun(X_S, y)
      
      # --------------------------------------------------------
      # CRITERIO (Build the full-length coefficient vector, padding
      # zeros for excluded variables, then evaluate the penalized
      # criterion on it)
      # --------------------------------------------------------
      
      beta=rep(0,ncol(X))
      beta[groups %in% S]=fit$coefficients
      
      criterion <- criterion_fun(beta, X, y)
      
      # Salvo tutto (store subset, its size, the fit, and its criterion value)
      results[[counter]] <- list(
        subset = S,
        size = length(S),
        fit = fit,
        criterion = criterion
      )
      
      counter <- counter + 1
    }
  }
  
  # ------------------------------------------------------------
  # Trova il migliore
  # (Find the subset with the smallest criterion value -> the "best" model)
  # ------------------------------------------------------------
  
  criteria <- sapply(results, function(x) x$criterion)
  best <- which.min(criteria) #minimizziamo (we minimize the criterion)
  best=results[[best]]
  
  # Flatten all results into a matrix: each column is (subset indices, criterion)
  results=sapply(results, function(x) c(x$subset,x$criterion) )
  
  return(list(
    best = best,          # best model found (subset + fit + criterion)
    all_results = results,# summary of every subset tried
    groups=groups         # the grouping vector used
  ))
}

# Ordinary least squares fit via base R's lm.fit (used as fit_fun for Gaussian case)
fit_ols <- function(X, y) {
  lm.fit(X, y)
}

#LAMBDA (regularization / penalty weight per included variable)
##### gaussian #### -----------------------------------------------------
# SETTING 1: Gaussian linear regression, UNGROUPED variables
# Compares exhaustive best subset selection vs. mm_gaussian() over B=1000
# simulated datasets with p=12 correlated predictors, only 4 of which are
# truly non-zero.
# -------------------------------------------------------------------------
set.seed(125)
n <- 100
p <- 12
rho_seq=log(n)/2*10   # penalty weight (BIC-like, scaled up) for the L0 term

# Penalized RSS criterion: RSS + rho * (number of non-zero coefficients)
criterion_mse <- function(beta, X, y) {
  
  rss <- sum( (y-X%*%beta)^2)+rho_seq*sum(beta!=0)
  rss
}

#X <- matrix(rnorm(n * p), n, p)
# Build a design matrix X with equicorrelated columns (correlation 0.5, unit variance)
M=matrix(0.5,p,p)
diag(M)=1
X=rmvnorm(n=n,sigma = M)
X=scale(X)                     # standardize columns (mean 0, sd 1)

# True coefficient vector: 4 non-zero signals, rest are zero
beta <- c(2, 1.5, -2, -1.5)
beta=c(beta,rep(0,p-length(beta)))

B=1000   # number of Monte Carlo replications
MM_fit=matrix(NA,nrow = p,ncol = B)    # storage for MM algorithm estimates
BSS_fit=matrix(NA,nrow = p,ncol = B)   # storage for best-subset estimates
losses=matrix(NA,nrow=B,ncol = 2)      # storage for criterion values (BSS vs MM)
colnames(losses)=c("BSS","MM")

for(i in 1:B){
  # Simulate a new response for the same X: linear signal + Gaussian noise
  y <- X %*% beta + rnorm(n)
  
  # Exhaustive best subset selection (exact optimum of the penalized criterion)
  res <- best_subset(
    X = X,
    y = y,
    fit_fun = fit_ols,
    criterion_fun = criterion_mse)
  
  beta_tmp=rep(0,p)
  beta_tmp[res$best$subset]=res$best$fit$coefficients
  BSS_fit[,i]=beta_tmp
  
  # Proposed MM-algorithm approximate solution to the same penalized problem
  fit_tmp= mm_gaussian(X,y,lambda =0,rho = rho_seq, max_iter_out = 300,
                       max_iter_in = 500)
  
  MM_fit[,i]=fit_tmp$beta
  
  # Record the penalized-criterion value achieved by each method
  losses[i,]=c(criterion_mse(beta_tmp,X=X,y=y),criterion_mse(fit_tmp$beta,X=X,y=y))
  print(i)
}

# --- Quick numerical comparisons -----------------------------------------
cbind(beta,rowMeans(BSS_fit),rowMeans(MM_fit))   # true beta vs average estimates
rowMeans(BSS_fit!=0)   # empirical selection frequency per variable, BSS
rowMeans(MM_fit!=0)    # empirical selection frequency per variable, MM
hist(losses[,1]-losses[,2])   # distribution of (BSS loss - MM loss): how close MM gets to the exact optimum
differ=BSS_fit-MM_fit
rowMeans(differ)     # average difference in estimated coefficients (bias of MM relative to BSS)
rowMeans(differ^2)   # mean squared difference per coefficient

##### grafico ##### (PLOT: visualize estimates and MM-vs-BSS discrepancy)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) # Per unire i due grafici (to combine the two panels)

# --- 1. Preparazione dei dati ---
# (Data preparation)
# Assicuriamoci che BSS_fit e MM_fit siano matrici 1000 x 12
# (Make sure BSS_fit and MM_fit are 1000 x 12 matrices:
#  1000 simulations x 12 coefficients)
p <- length(beta)
sim_data <- data.frame()
for (j in 1:p) {
  bss_j <- BSS_fit[j, ]
  mm_j  <- MM_fit[j, ]
  
  # Build one summary row per coefficient: mean and 95% empirical interval
  # for both methods, plus the mean squared deviation between them
  sim_data <- rbind(sim_data, data.frame(
    Coef = factor(paste0("beta[", j, "]"), levels = paste0("beta[", 1:p, "]")),
    Coef_idx = j,
    True_Beta = beta[j],
    # Stime e Intervalli al 95% per BSS (estimates & 95% intervals for BSS)
    BSS_Mean = mean(bss_j),
    BSS_Low  = quantile(bss_j, 0.025),
    BSS_High = quantile(bss_j, 0.975),
    # Stime e Intervalli al 95% per MM (estimates & 95% intervals for MM)
    MM_Mean  = mean(mm_j),
    MM_Low   = quantile(mm_j, 0.025),
    MM_High  = quantile(mm_j, 0.975),
    # Mean Squared Deviation tra MM e BSS (MSD between MM and BSS)
    MSD      = mean((mm_j - bss_j)^2)
  ))
}

# Format lungo per il primo grafico (stima delle medie e intervalli)
# (Reshape to long format for the first plot: point estimates + intervals)
plot_estimates <- bind_rows(
  sim_data %>% transmute(Coef, Coef_idx, True_Beta, Method = "MM (Proposed)", Mean = MM_Mean, Low = MM_Low, High = MM_High),
  sim_data %>% transmute(Coef, Coef_idx, True_Beta, Method = "BSS (Exact)", Mean = BSS_Mean, Low = BSS_Low, High = BSS_High)
)

# --- 2. Grafico A: Stime e Intervalli vs Valori Veri ---
# (Panel A: point estimates + 95% intervals vs. true parameter values)
p1 <- ggplot(plot_estimates, aes(x = Coef, y = Mean, color = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  # Valore Vero (Rombo Nero) -- true value shown as a black diamond
  geom_point(aes(y = True_Beta), shape = 18, size = 4, color = "black",
             position = position_nudge(x = 0), show.legend = FALSE) +
  # Intervalli di confidenza 95% e medie dei metodi
  # (95% intervals and point estimates for each method)
  geom_pointrange(aes(ymin = Low, ymax = High),
                  position = position_dodge(width = 0.5), size = 0.5) +
  scale_color_manual(values = c("BSS (Exact)" = "#2B5C8F", "MM (Proposed)" = "#D95F02")) +
  scale_x_discrete(labels = parse(text = paste0("beta[", 1:p, "]"))) +
  labs(
    title = "Gaussian data",
    subtitle = "Black diamonds represent true parameter values",
    x = NULL,
    y = "Coefficient Value",
    color = "Method"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# --- 3. Grafico B: Scarto Medio Quadratico (MSD) MM vs BSS ---
# (Panel B: Mean Squared Deviation between MM and BSS, per coefficient)
p2 <- ggplot(sim_data, aes(x = Coef, y = MSD)) +
  geom_col(fill = "#7570B3", width = 0.55, alpha = 0.85) +
  geom_text(aes(label = round(MSD, 4)), vjust = -0.5, size = 3, color = "gray20") +
  scale_x_discrete(labels = parse(text = paste0("beta[", 1:p, "]"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    subtitle = "Mean Squared Deviation (MM vs BSS)",
    x = "Parameters",
    y = "MSD"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# --- 4. Unione e Salvataggio ---
# (Combine panels A and B into one figure and optionally save it)
poster_plot1 <- p1 / p2 + plot_layout(heights = c(2, 1))

# Visualizza il grafico (display the figure)
print(poster_plot1)

# Per salvarlo in alta risoluzione per il poster (PDF o PNG)
# (To save it in high resolution for the poster, PDF or PNG)
# ggsave("simulation_results.pdf", poster_plot, width = 10, height = 8)

##### logit ##### -------------------------------------------------------
# SETTING 2: Logistic regression, UNGROUPED variables
# Same comparison (BSS vs MM) but for a binary response modeled via logistic
# regression, using mm_logit() as the MM-algorithm counterpart.
# -------------------------------------------------------------------------
loglik_logit <- function(beta, X, y) {
  
  eta <- as.vector(X %*% beta)
  
  # Log-verosimiglianza (Bernoulli/logistic log-likelihood, computed in a
  # numerically stable way via log1p)
  loglik <- sum(y * eta - log1p(exp(eta)))
  
  return(loglik)
}

# Logistic regression fit via base R's glm.fit (used as fit_fun for logit case)
fit_glm <- function(X, y) {
  glm.fit(X,y,family = binomial(link = "logit"))
}

#LAMBDA
set.seed(125)
n <- 100
p <- 12
rho_seq=log(n)/8 #(BIC style)   -- penalty weight for the L0 term

# Penalized negative log-likelihood criterion: -loglik + rho * (# non-zero coefs)
criterion_loglik <- function(beta, X, y) {
  loss <- -loglik_logit(beta=beta,X=X,y=y)+rho_seq*sum(beta!=0)
  loss
}

# Same correlated design matrix structure as in the Gaussian case
M=matrix(0.5,p,p)
diag(M)=1
X=rmvnorm(n=n,sigma = M)
X=scale(X)
beta <- c(2, 1.5, -2, -1.5)
beta=c(beta,rep(0,p-length(beta)))

B=1000
MM_fit=matrix(NA,nrow = p,ncol = B)
BSS_fit=matrix(NA,nrow = p,ncol = B)
losses=matrix(NA,nrow=B,ncol = 2)
colnames(losses)=c("BSS","MM")

for(i in 1:B){
  eta=X %*% beta #+ rnorm(n)
  # Simulate a binary response from the logistic model with linear predictor eta
  y <- rbinom(n=n,size=1,prob=exp(eta)/(1+exp(eta)))
  
  res <- best_subset(
    X = X,
    y = y,
    fit_fun = fit_glm,
    criterion_fun = criterion_loglik
  )
  
  beta_tmp=rep(0,p)
  beta_tmp[res$best$subset]=res$best$fit$coefficients
  BSS_fit[,i]=beta_tmp
  
  fit_tmp= mm_logit(X,y,lambda =0,rho = rho_seq, max_iter_out = 300,
                    max_iter_in = 500)
  
  MM_fit[,i]=fit_tmp$beta
  
  losses[i,]=c(criterion_loglik(beta_tmp,X=X,y=y),criterion_loglik(fit_tmp$beta,X=X,y=y))
  
  print(i)
}

# --- Quick numerical comparisons -----------------------------------------
rowMeans(BSS_fit!=0)
rowMeans(MM_fit!=0)
cbind(beta,rowMeans(BSS_fit),rowMeans(MM_fit))
differ=BSS_fit-MM_fit
rowMeans(differ)
rowMeans(differ^2)
hist(losses[,1]-losses[,2])

##### grafico ##### (same plotting logic as before, for the logistic case)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) # Per unire i due grafici

# --- 1. Preparazione dei dati ---
# Assicuriamoci che BSS_fit e MM_fit siano matrici 1000 x 12
# (1000 simulazioni x 12 coefficienti)
p <- length(beta)
sim_data <- data.frame()
for (j in 1:p) {
  bss_j <- BSS_fit[j, ]
  mm_j  <- MM_fit[j, ]
  
  sim_data <- rbind(sim_data, data.frame(
    Coef = factor(paste0("beta[", j, "]"), levels = paste0("beta[", 1:p, "]")),
    Coef_idx = j,
    True_Beta = beta[j],
    # Stime e Intervalli al 95% per BSS
    BSS_Mean = mean(bss_j),
    BSS_Low  = quantile(bss_j, 0.025),
    BSS_High = quantile(bss_j, 0.975),
    # Stime e Intervalli al 95% per MM
    MM_Mean  = mean(mm_j),
    MM_Low   = quantile(mm_j, 0.025),
    MM_High  = quantile(mm_j, 0.975),
    # Mean Squared Deviation tra MM e BSS
    MSD      = mean((mm_j - bss_j)^2)
  ))
}

# Format lungo per il primo grafico (stima delle medie e intervalli)
plot_estimates <- bind_rows(
  sim_data %>% transmute(Coef, Coef_idx, True_Beta, Method = "MM (Proposed)", Mean = MM_Mean, Low = MM_Low, High = MM_High),
  sim_data %>% transmute(Coef, Coef_idx, True_Beta, Method = "BSS (Exact)", Mean = BSS_Mean, Low = BSS_Low, High = BSS_High)
)

# --- 2. Grafico A: Stime e Intervalli vs Valori Veri ---
p1 <- ggplot(plot_estimates, aes(x = Coef, y = Mean, color = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  # Valore Vero (Rombo Nero)
  geom_point(aes(y = True_Beta), shape = 18, size = 4, color = "black",
             position = position_nudge(x = 0), show.legend = FALSE) +
  # Intervalli di confidenza 95% e medie dei metodi
  geom_pointrange(aes(ymin = Low, ymax = High),
                  position = position_dodge(width = 0.5), size = 0.5) +
  scale_color_manual(values = c("BSS (Exact)" = "#2B5C8F", "MM (Proposed)" = "#D95F02")) +
  scale_x_discrete(labels = parse(text = paste0("beta[", 1:p, "]"))) +
  labs(
    title = "Binary data",
    subtitle = "Black diamonds represent true parameter values",
    x = NULL,
    y = "Coefficient Value",
    color = "Method"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# --- 3. Grafico B: Scarto Medio Quadratico (MSD) MM vs BSS ---
p2 <- ggplot(sim_data, aes(x = Coef, y = MSD)) +
  geom_col(fill = "#7570B3", width = 0.55, alpha = 0.85) +
  geom_text(aes(label = round(MSD, 4)), vjust = -0.5, size = 3, color = "gray20") +
  scale_x_discrete(labels = parse(text = paste0("beta[", 1:p, "]"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    subtitle = "Mean Squared Deviation (MM vs BSS)",
    x = "Parameters",
    y = "MSD"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# --- 4. Unione e Salvataggio ---
poster_plot2 <- p1 / p2 + plot_layout(heights = c(2, 1))

# Visualizza il grafico
print(poster_plot2)

# Per salvarlo in alta risoluzione per il poster (PDF o PNG)
# ggsave("simulation_results.pdf", poster_plot, width = 10, height = 8)

##### grouped gaussian ##### ---------------------------------------------
# SETTING 3: Gaussian linear regression, GROUPED variables
# Variables are organized into groups; the L0 penalty now applies at the
# GROUP level (a group is either entirely zero or entirely non-zero), and
# best_subset() searches over subsets of groups rather than single columns.
# mm_grouped_gaussian() is the MM-algorithm counterpart. An additional
# "MMad" (MM + adjustment) estimator is also computed: it refits OLS on the
# variables selected by MM_fit, to see if a post-selection refit closes the
# gap with BSS.
# -------------------------------------------------------------------------
fit_ols <- function(X, y) {
  lm.fit(X, y)
}

set.seed(125)
n <- 200
p <- 12 #number of groups
rho_seq=log(n)/2

# Grouping vector: groups 1..(p-5) each contain a single variable, and the
# last 5 groups each contain 3 variables (so total variables > p)
groups=c(1:(p-5) ,sapply((p-4):p,function(x) rep(x,3)))

M=matrix(0.5,length(groups),length(groups))
diag(M)=1
X=rmvnorm(n=n,sigma = M)
X=scale(X)

# True coefficients: several non-zero blocks matching some of the groups
beta <- c(2.0,  1.5 ,-2.0 ,-1.5,0,0,0,-2,1.5,-1.5,rep(0,3),
          2,-1.5,1.5)
beta=c(beta,rep(0,length(groups)-length(beta)))

# Group sizes, used to weight the penalty by how many variables a group contains
group_weights=table(groups)

# Penalized RSS criterion at the GROUP level: RSS + rho * (sum of sizes of
# groups that are active, i.e. have at least one non-zero coefficient)
criterion_mse <- function(beta, X, y) {
  
  AA= sapply(unique(groups),function(x) any(beta[groups==x]!=0) )
  
  rss <- sum( (y-X%*%beta)^2)+rho_seq*sum(group_weights[AA])
  rss
}

B=1000
MM_fit=matrix(NA,nrow = ncol(X),ncol = B)     # MM-algorithm estimates
MMad_fit=matrix(NA,nrow = ncol(X),ncol = B)   # MM + post-selection OLS refit
BSS_fit=matrix(NA,nrow = ncol(X),ncol = B)    # exact best-subset (group-wise) estimates
losses=matrix(NA,nrow=B,ncol = 3)
colnames(losses)=c("BSS","MM","MMadj")

for(i in 1:B){
  y <- X %*% beta + rnorm(n)
  
  # Exact group best-subset selection
  res <- best_subset(
    X = X,
    y = y,
    fit_fun = fit_ols,
    criterion_fun = criterion_mse,groups = groups
  )
  
  beta_tmp=rep(0,ncol(X))
  beta_tmp[groups %in% res$best$subset]=res$best$fit$coefficients
  BSS_fit[,i]=beta_tmp
  
  # MM-algorithm approximate group-sparse solution
  fit_tmp= mm_grouped_gaussian(X,y,groups=groups,lambda =0,rho = rho_seq, max_iter_out = 300,
                               max_iter_in = 500,group_weights = group_weights)
  
  MM_fit[,i]=fit_tmp$beta
  
  # Post-selection OLS refit: take the variables selected by MM and refit
  # OLS on just those columns (an "adjusted" MM estimator)
  Xs=X[,MM_fit[,i]!=0]
  MMad_fit[,i]=0
  MMad_fit[MM_fit[,i]!=0,i]=lm.fit(Xs,y)$coef
  
  losses[i,]=c(criterion_mse(beta_tmp,X=X,y=y),criterion_mse(fit_tmp$beta,X=X,y=y),
               criterion_mse(MMad_fit[,i],X=X,y=y))
  print(i)
}

# --- Quick numerical comparisons -----------------------------------------
rowMeans(BSS_fit!=0)
rowMeans(MM_fit!=0)
#save(BSS_fit,MM_fit,MMad_fit,losses,file = "20agoGAUSS.Rdata")
load("20agoGAUSS.Rdata")

##### grafico a 2 ##### (2-panel plot: coefficient-level and group-level comparison)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# --- 1. Parametri di simulazione e gruppi ---
# Assumiamo che 'groups' e 'beta' siano già definiti nel tuo ambiente
# Es: groups <- c(1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4)
p <- length(beta)
unique_groups <- unique(groups)
num_groups <- length(unique_groups)

# Calcolo dei punti di interruzione dinamici (indipendentemente dalla dimensione)
# (Compute where vertical separator lines should go between groups, on the x-axis)
group_counts <- table(factor(groups, levels = unique_groups))
group_separators <- cumsum(group_counts)[-num_groups] + 0.5

# --- 2. Dati al livello di singolo parametro (per Pannello A) ---
# (Per-coefficient summary data, for Panel A)
param_data <- data.frame()
for (j in 1:p) {
  bss_j <- BSS_fit[j, ]
  mm_j  <- MM_fit[j, ]
  
  param_data <- rbind(param_data, data.frame(
    Coef = factor(paste0("beta[", j, "]"), levels = paste0("beta[", 1:p, "]")),
    Group = factor(paste0("G[", groups[j], "]")), # <-- Assegna il gruppo dinamico (assign dynamic group label)
    Coef_idx = j,
    True_Beta = beta[j],
    BSS_Mean = mean(bss_j),
    BSS_Low  = quantile(bss_j, 0.025),
    BSS_High = quantile(bss_j, 0.975),
    MM_Mean  = mean(mm_j),
    MM_Low   = quantile(mm_j, 0.025),
    MM_High  = quantile(mm_j, 0.975)
  ))
}

plot_estimates <- bind_rows(
  param_data %>% transmute(Coef, Group, Coef_idx, True_Beta, Method = "MM (Proposed)", Mean = MM_Mean, Low = MM_Low, High = MM_High),
  param_data %>% transmute(Coef, Group, Coef_idx, True_Beta, Method = "BSS (Exact)", Mean = BSS_Mean, Low = BSS_Low, High = BSS_High)
)

# --- 3. Dati aggregati AL LIVELLO DI GRUPPO (per Pannello B) ---
# (Group-level aggregated data, for Panel B)
group_data <- data.frame()
for (g in unique_groups) {
  # Trova esattamente quali coefficienti appartengono al gruppo g
  # (Find which coefficients belong to group g)
  cols_in_group <- which(groups == g)
  
  # Matrici ridotte per il gruppo g (restrict to the columns in group g)
  bss_group <- BSS_fit[cols_in_group, ]
  mm_group  <- MM_fit[cols_in_group, ]
  
  # Errore medio all'interno del gruppo g (Mean Squared Deviation)
  # (Average within-group MSD between MM and BSS)
  group_msd <- mean((mm_group - bss_group)^2)
  
  group_data <- rbind(group_data, data.frame(
    Group_Label = factor(paste0("G[", g, "]")),
    Group_Index = g,
    Group_MSD   = group_msd
  ))
}

# --- 4. Pannello A: Stime dei Coefficienti e Separatori ---
# (Panel A: coefficient estimates with group separator lines)
p1 <- ggplot(plot_estimates, aes(x = Coef, y = Mean, color = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  # === LINEE SEPARATRICI DINAMICHE === (dynamic separator lines between groups)
  geom_vline(xintercept = group_separators, linetype = "dotted", color = "gray60", size = 0.9) +
  # Valore Vero (Rombo Nero)
  geom_point(aes(y = True_Beta), shape = 18, size = 4, color = "black",
             position = position_nudge(x = 0), show.legend = FALSE) +
  # Intervalli di confidenza 95% e medie dei metodi
  geom_pointrange(aes(ymin = Low, ymax = High),
                  position = position_dodge(width = 0.5), size = 0.5) +
  scale_color_manual(values = c("BSS (Exact)" = "#2B5C8F", "MM (Proposed)" = "#D95F02")) +
  scale_x_discrete(labels = parse(text = paste0("beta[", 1:p, "]"))) +
  labs(
    title = "Gaussian data",
    subtitle = "Black diamonds represent true parameter values. Dotted lines separate groups.",
    x = NULL,
    y = "Coefficient Value",
    color = "Method"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# --- 5. Pannello B: Errore Medio Entro-Gruppo (Within-Group MSD) ---
# (Panel B: within-group Mean Squared Deviation, MM vs BSS)
group_data$Group_MSD=round(group_data$Group_MSD,3)
# --- 5. Pannello B: Errore Medio Entro-Gruppo (Within-Group MSD) ---
p2 <- ggplot(group_data, aes(x = Group_Label, y = Group_MSD)) +
  geom_col(fill = "#5C6BC0", width = 0.5, alpha = 0.9) +
  # Etichette numeriche sopra ogni barra (numeric labels above each bar)
  geom_text(aes(label = sprintf("%.3f", Group_MSD)), vjust = -0.5, size = 2, fontface = "bold", color = "gray20") +
  scale_x_discrete(labels = parse(text = paste0("G[", unique_groups, "]"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    subtitle = "Within-Group Mean Squared Deviation (MM vs BSS)",
    x = "Variable Groups",
    y = "Group MSD"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# --- 6. Composizione finale --- (final composition of panels A and B)
poster_plot3 <- p1 / p2 + plot_layout(heights = c(2, 1))
print(poster_plot3)

##### grafico a 3 ##### (numeric comparisons of BSS vs MM vs MMad)
cbind(beta,rowMeans(BSS_fit),rowMeans(MM_fit),rowMeans(MMad_fit))
differ=BSS_fit-MM_fit
rowMeans(differ)
rowMeans(differ^2)
hist(losses[,1]-losses[,2])
hist(losses[,1]-losses[,3])

##### grouped logit ##### -------------------------------------------------
# SETTING 4: Logistic regression, GROUPED variables
# Same group-sparse comparison as Setting 3, but for a binary response,
# using mm_grouped_logit() as the MM-algorithm counterpart.
# ---------------------------------------------------------------------------
loglik_logit <- function(beta, X, y) {
  
  eta <- as.vector(X %*% beta)
  
  # Log-verosimiglianza (logistic log-likelihood)
  loglik <- sum(y * eta - log1p(exp(eta)))
  
  return(loglik)
}

fit_glm <- function(X, y) {
  glm.fit(X,y,family = binomial(link = "logit"))
}

set.seed(12)
n <- 200
p <- 12 #number of groups
rho_seq=log(n)/2
groups=c(1:(p-5) ,sapply((p-4):p,function(x) rep(x,3)))
M=matrix(0.5,length(groups),length(groups))
diag(M)=1
X=rmvnorm(n=n,sigma = M)
X=scale(X)
beta <- c(2.0,  1.5 ,-2.0 ,-1.5,0,0,0,-2,1.5,-1.5,rep(0,3),
          2,-1.5,1.5)
beta=c(beta,rep(0,length(groups)-length(beta)))
group_weights=table(groups)

# Penalized negative log-likelihood at the group level
criterion_loglik <- function(beta, X, y) {
  AA= sapply(unique(groups),function(x) any(beta[groups==x]!=0) )
  
  loss <- -loglik_logit(beta=beta,X=X,y=y)+rho_seq*sum(group_weights[AA])
  loss
}

B=1000
MM_fit=matrix(NA,nrow = ncol(X),ncol = B)
MMad_fit=matrix(NA,nrow = ncol(X),ncol = B)
BSS_fit=matrix(NA,nrow = ncol(X),ncol = B)
losses=matrix(NA,nrow=B,ncol = 3)
colnames(losses)=c("BSS","MM","MMadj")

for(i in 1:B){
  eta=X %*% beta #+ rnorm(n)
  y <- rbinom(n=n,size=1,prob=exp(eta)/(1+exp(eta)))
  
  res <- best_subset(
    X = X,
    y = y,
    fit_fun = fit_glm,
    criterion_fun = criterion_loglik,groups = groups
  )
  
  beta_tmp=rep(0,ncol(X))
  beta_tmp[groups %in% res$best$subset]=res$best$fit$coefficients
  BSS_fit[,i]=beta_tmp
  
  fit_tmp= mm_grouped_logit(X,y,groups=groups,lambda =0,rho = rho_seq, max_iter_out = 300,
                            max_iter_in = 500,group_weights = group_weights,tau3=1e-12)
  
  MM_fit[,i]=fit_tmp$beta
  
  # Post-selection refit via GLM on the variables selected by MM
  Xs=X[,MM_fit[,i]!=0]
  MMad_fit[,i]=0
  MMad_fit[MM_fit[,i]!=0,i]=glm.fit(Xs,y,family = binomial(link = "logit"))$coefficients
  
  print(i)
  
  losses[i,]=c(criterion_loglik(beta_tmp,X=X,y=y),criterion_loglik(fit_tmp$beta,X=X,y=y),
               criterion_loglik(MMad_fit[,i],X=X,y=y))
}

#save(BSS_fit,MM_fit,MMad_fit,losses,file = "20agologit.Rdata")
rowMeans(BSS_fit!=0)
rowMeans(MM_fit!=0)
cbind(beta,rowMeans(BSS_fit),rowMeans(MM_fit),rowMeans(MMad_fit))
differ=BSS_fit-MM_fit
rowMeans(differ)
rowMeans(differ^2)
hist(losses[,1]-losses[,2])
hist(losses[,1]-losses[,3])
load("20agologit.Rdata")

##### grafico a 2 ##### (same 2-panel plot logic, for grouped logistic case)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# --- 1. Parametri di simulazione e gruppi ---
# Assumiamo che 'groups' e 'beta' siano già definiti nel tuo ambiente
# Es: groups <- c(1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4)
p <- length(beta)
unique_groups <- unique(groups)
num_groups <- length(unique_groups)

# Calcolo dei punti di interruzione dinamici (indipendentemente dalla dimensione)
group_counts <- table(factor(groups, levels = unique_groups))
group_separators <- cumsum(group_counts)[-num_groups] + 0.5

# --- 2. Dati al livello di singolo parametro (per Pannello A) ---
param_data <- data.frame()
for (j in 1:p) {
  bss_j <- BSS_fit[j, ]
  mm_j  <- MM_fit[j, ]
  
  param_data <- rbind(param_data, data.frame(
    Coef = factor(paste0("beta[", j, "]"), levels = paste0("beta[", 1:p, "]")),
    Group = factor(paste0("G[", groups[j], "]")), # <-- Assegna il gruppo dinamico
    Coef_idx = j,
    True_Beta = beta[j],
    BSS_Mean = mean(bss_j),
    BSS_Low  = quantile(bss_j, 0.025),
    BSS_High = quantile(bss_j, 0.975),
    MM_Mean  = mean(mm_j),
    MM_Low   = quantile(mm_j, 0.025),
    MM_High  = quantile(mm_j, 0.975)
  ))
}

plot_estimates <- bind_rows(
  param_data %>% transmute(Coef, Group, Coef_idx, True_Beta, Method = "MM (Proposed)", Mean = MM_Mean, Low = MM_Low, High = MM_High),
  param_data %>% transmute(Coef, Group, Coef_idx, True_Beta, Method = "BSS (Exact)", Mean = BSS_Mean, Low = BSS_Low, High = BSS_High)
)

# --- 3. Dati aggregati AL LIVELLO DI GRUPPO (per Pannello B) ---
group_data <- data.frame()
for (g in unique_groups) {
  # Trova esattamente quali coefficienti appartengono al gruppo g
  cols_in_group <- which(groups == g)
  
  # Matrici ridotte per il gruppo g
  bss_group <- BSS_fit[cols_in_group, ]
  mm_group  <- MM_fit[cols_in_group, ]
  
  # Errore medio all'interno del gruppo g (Mean Squared Deviation)
  group_msd <- mean((mm_group - bss_group)^2)
  
  group_data <- rbind(group_data, data.frame(
    Group_Label = factor(paste0("G[", g, "]")),
    Group_Index = g,
    Group_MSD   = group_msd
  ))
}

# --- 4. Pannello A: Stime dei Coefficienti e Separatori ---
p1 <- ggplot(plot_estimates, aes(x = Coef, y = Mean, color = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  # === LINEE SEPARATRICI DINAMICHE ===
  geom_vline(xintercept = group_separators, linetype = "dotted", color = "gray60", size = 0.9) +
  # Valore Vero (Rombo Nero)
  geom_point(aes(y = True_Beta), shape = 18, size = 4, color = "black",
             position = position_nudge(x = 0), show.legend = FALSE) +
  # Intervalli di confidenza 95% e medie dei metodi
  geom_pointrange(aes(ymin = Low, ymax = High),
                  position = position_dodge(width = 0.5), size = 0.5) +
  scale_color_manual(values = c("BSS (Exact)" = "#2B5C8F", "MM (Proposed)" = "#D95F02")) +
  scale_x_discrete(labels = parse(text = paste0("beta[", 1:p, "]"))) +
  labs(
    title = "Binary data",
    subtitle = "Black diamonds represent true parameter values. Dotted lines separate groups.",
    x = NULL,
    y = "Coefficient Value",
    color = "Method"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

group_data$Group_MSD=round(group_data$Group_MSD,2)
# --- 5. Pannello B: Errore Medio Entro-Gruppo (Within-Group MSD) ---
p2 <- ggplot(group_data, aes(x = Group_Label, y = Group_MSD)) +
  geom_col(fill = "#5C6BC0", width = 0.5, alpha = 0.9) +
  # Etichette numeriche sopra ogni barra
  geom_text(aes(label = sprintf("%.2f", Group_MSD)), vjust = -0.5, size = 2, fontface = "bold", color = "gray20") +
  scale_x_discrete(labels = parse(text = paste0("G[", unique_groups, "]"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    subtitle = "Within-Group Mean Squared Deviation (MM vs BSS)",
    x = "Variable Groups",
    y = "Group MSD"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# --- 6. Composizione finale ---
poster_plot4 <- p1 / p2 + plot_layout(heights = c(2, 1))
print(poster_plot4)

# Combine the two group-level poster figures (Gaussian + Binary) side by side
TOT=(poster_plot3 | poster_plot4) + plot_layout(ncol = 2, nrow = 1)

ggsave(filename = "plot1.pdf",plot=poster_plot3,
       width = 30, height = 20, units = "cm")
ggsave(filename = "plot2.pdf",plot=poster_plot4,
       width = 30, height = 20, units = "cm")
ggsave(filename = "plot3.pdf",plot=TOT,
       width = 40, height = 20, units = "cm")

##### grafico a 3#### -----------------------------------------------------
# 3-METHOD COMPARISON PLOT (grouped logistic case): BSS vs MM vs MMad
# Adds the post-selection-refit estimator (MMad) to the comparison, both at
# the coefficient level (Panel A) and as a group-level distance-from-BSS
# trend (Panel B).
# ---------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# --- 1. Parametri di simulazione e gruppi ---
# Assumiamo che 'groups', 'beta', 'BSS_fit', 'MM_fit', e 'MMad_fit' siano già definiti
p <- length(beta)
unique_groups <- unique(groups)
num_groups <- length(unique_groups)

# Calcolo dei punti di interruzione dinamici (indipendentemente dalla dimensione)
group_counts <- table(factor(groups, levels = unique_groups))
group_separators <- cumsum(group_counts)[-num_groups] + 0.5

# --- 2. Dati al livello di singolo parametro (per Pannello A) ---
param_data <- data.frame()
for (j in 1:p) {
  bss_j  <- BSS_fit[j, ]
  mm_j   <- MM_fit[j, ]
  mmad_j <- MMad_fit[j, ] # <-- Aggiunto MMad (added MMad column)
  
  param_data <- rbind(param_data, data.frame(
    Coef = factor(paste0("beta[", j, "]"), levels = paste0("beta[", 1:p, "]")),
    Group = factor(paste0("G[", groups[j], "]")),
    Coef_idx = j,
    True_Beta = beta[j],
    BSS_Mean = mean(bss_j),
    BSS_Low  = quantile(bss_j, 0.025),
    BSS_High = quantile(bss_j, 0.975),
    MM_Mean  = mean(mm_j),
    MM_Low   = quantile(mm_j, 0.025),
    MM_High  = quantile(mm_j, 0.975),
    MMad_Mean = mean(mmad_j),
    MMad_Low  = quantile(mmad_j, 0.025),
    MMad_High = quantile(mmad_j, 0.975)
  ))
}

# Combiniamo per ggplot (ora con 3 metodi)
# (Combine into long format for ggplot, now with 3 methods)
plot_estimates <- bind_rows(
  param_data %>% transmute(Coef, Group, Coef_idx, True_Beta, Method = "MMad", Mean = MMad_Mean, Low = MMad_Low, High = MMad_High),
  param_data %>% transmute(Coef, Group, Coef_idx, True_Beta, Method = "MM", Mean = MM_Mean, Low = MM_Low, High = MM_High),
  param_data %>% transmute(Coef, Group, Coef_idx, True_Beta, Method = "BSS (Exact)", Mean = BSS_Mean, Low = BSS_Low, High = BSS_High)
)

# Impostiamo l'ordine dei livelli per la legenda
# (Set factor level order for the legend)
plot_estimates$Method <- factor(plot_estimates$Method, levels = c("BSS (Exact)", "MM", "MMad"))

# --- 3. Dati aggregati AL LIVELLO DI GRUPPO (per Pannello B) ---
group_data <- data.frame()
for (g in unique_groups) {
  cols_in_group <- which(groups == g)
  
  bss_group  <- BSS_fit[cols_in_group, ]
  mm_group   <- MM_fit[cols_in_group, ]
  mmad_group <- MMad_fit[cols_in_group, ] # <-- Aggiunto MMad
  
  # Calcolo errori quadratici medi rispetto a BSS
  # (Compute mean squared error of MM and MMad relative to BSS, per group)
  msd_mm   <- mean((mm_group - bss_group)^2)
  msd_mmad <- mean((mmad_group - bss_group)^2)
  
  group_data <- rbind(group_data,
                      data.frame(Group_Label = factor(paste0("G[", g, "]")), Group_Index = g, Method = "MM", MSD = msd_mm),
                      data.frame(Group_Label = factor(paste0("G[", g, "]")), Group_Index = g, Method = "MMad", MSD = msd_mmad)
  )
}

# --- 4. Pannello A: Stime dei Coefficienti e Separatori (3 Metodi) ---
# (Panel A: coefficient estimates and separators, 3 methods)
# Colori scelti: Blu (BSS), Arancione (MM), Verde (MMad)
# (Chosen colors: Blue = BSS, Orange = MM, Green = MMad)
color_palette <- c("BSS (Exact)" = "#2B5C8F", "MM" = "#D95F02", "MMad" = "#1B9E77")

p1 <- ggplot(plot_estimates, aes(x = Coef, y = Mean, color = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  # Linee separatrici (separator lines)
  geom_vline(xintercept = group_separators, linetype = "dotted", color = "gray60", size = 0.6) +
  # Valore Vero (true value)
  geom_point(aes(y = True_Beta), shape = 18, size = 5, color = "black",
             position = position_nudge(x = 0), show.legend = FALSE) +
  # Intervalli (allargato il dodge width per far spazio a 3 metodi)
  # (Intervals; dodge width widened to fit 3 methods side by side)
  geom_pointrange(aes(ymin = Low, ymax = High),
                  position = position_dodge(width = 0.7), size = 0.5) +
  scale_color_manual(values = color_palette) +
  scale_x_discrete(labels = parse(text = paste0("beta[", 1:p, "]"))) +
  labs(
    title = "Binary data",
    subtitle = "Black diamonds represent true parameter values. Dotted lines separate groups.",
    x = NULL,
    y = "Coefficient Value",
    color = "Method"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# --- 5. Pannello B: Errore Entro-Gruppo rispetto a BSS (Linee e Punti) ---
# (Panel B: within-group error relative to BSS, shown as lines + points)
p2 <- ggplot(group_data, aes(x = Group_Label, y = MSD, color = Method, group = Method)) +
  # Linee e punti collegati (connected lines and points)
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 3) +
  # Testo per far vedere il valore (aggiustato per non sovrapporsi troppo)
  # (Text labels showing the value, adjusted to reduce overlap)
  geom_text(aes(label = sprintf("%.4f", MSD)),
            vjust = -1.2, size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = color_palette[c("MM", "MMad")]) +
  scale_x_discrete(labels = parse(text = paste0("G[", unique_groups, "]"))) +
  # Espandiamo un po' l'asse Y verso l'alto per non tagliare le etichette di testo
  # (Expand the y-axis a bit upward so text labels aren't cut off)
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.3))) +
  labs(
    subtitle = "Distance from BSS (Mean Squared Deviation) by Group",
    x = "Variable Groups",
    y = "Group MSD"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none", # Rimuovo la legenda qui perché è già nel p1
    # (Legend removed here because it's already shown in p1)
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# --- 6. Composizione finale --- (final composition)
poster_plot4 <- p1 / p2 + plot_layout(heights = c(2, 1))
print(poster_plot3)
print(poster_plot4)

