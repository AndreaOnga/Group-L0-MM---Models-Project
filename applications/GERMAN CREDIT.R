# ============================================================
# 1. "GERMAN CREDIT" DATA PREPARATION
# Building the grouped design matrix:
# - numeric vars -> spline bases (one group per variable)
# - categorical vars -> dummy encoding (one group per variable)
# ============================================================

if (!requireNamespace("rchallenge", quietly = TRUE)) install.packages("rchallenge")
library(rchallenge)
library(splines)
library(gglasso)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

data("german", package = "rchallenge")

# Binary response (1 = "bad" credit risk, 0 = "good")
y_bin <- as.numeric(german$credit_risk == "bad")
predictors <- german[, setdiff(names(german), "credit_risk")]

# Function to create the group-structured design matrix
build_group_design <- function(data, df_spline = 4) {
  X_list <- list()
  group_vec <- integer(0)
  group_labels <- character(0)
  g <- 0L
  
  for (nm in names(data)) {
    v <- data[[nm]]
    g <- g + 1L
    
    if (is.numeric(v)) {
      # Spline bases for continuous variables (default df=4)
      cols <- matrix(splines::ns(v, df = df_spline), ncol = df_spline)
      colnames(cols) <- paste0(nm, "_bs", seq_len(df_spline))
    } else {
      # Full dummy encoding (without intercept) for categorical variables
      v <- droplevels(as.factor(v))
      cols <- model.matrix(~ v - 1)
      colnames(cols) <- paste0(nm, "_", levels(v))
    }
    
    X_list[[nm]] <- cols
    group_vec <- c(group_vec, rep(g, ncol(cols)))
    group_labels[g] <- nm
  }
  
  list(X = do.call(cbind, X_list), group = group_vec, group_labels = group_labels)
}

built <- build_group_design(predictors)
X <- built$X
groups <- built$group
group_labels <- built$group_labels

# ============================================================
# 2. ANALYSIS SETUP AND CROSS-VALIDATION
# ============================================================

# source("CROSS VALIDATION FUNCTIONS.R") # Make sure your functions are loaded

set.seed(1)
train <- 1:1000 # Using all 1000 obs. for the fit
y_bin <- y_bin[train]
X <- X[train,]
n <- length(y_bin)

# Creating balanced folds for Cross-Validation
set.seed(1)
K <- 5
folds <- sample(1:K, replace = TRUE, size = n)

# ============================================================
# 3. MODEL 1 ESTIMATION: GROUP L0 (MM LOGIT)
# ============================================================

lambda_seq <- c(0)
rho_seq <- seq(log(n)/5, log(n)/4, length.out = 100) # Must increase for warm_start

out_L0 <- cv_mm_grouped_logit(
  X = X, y = y_bin, groups = groups, 
  lambda_seq = lambda_seq, rho_seq = rho_seq, best_fit = TRUE, 
  tau3 = 1e-12, folds = folds, intercept = TRUE, 
  group_weights = NULL, max_iter_out = 100, max_iter_in = 100
)

# Estimated L0 coefficients (including intercept)
beta_L0 <- out_L0$best_fit$beta 

# ============================================================
# 4. SCALING FUNCTIONS FOR GGLASSO (L1)
# Standardize the entire group while preserving relative geometry
# ============================================================

scale_X_grouped <- function(X, group) {
  n <- nrow(X)
  center <- colMeans(X)
  Xc <- sweep(X, 2, center, "-")
  scale_vec <- numeric(ncol(X))
  
  for (g in unique(group)) {
    idx <- which(group == g)
    s_g <- sqrt(sum(Xc[, idx]^2) / (n * length(idx)))
    scale_vec[idx] <- ifelse(s_g == 0, 1, s_g) # Avoid division by zero
  }
  
  list(X = sweep(Xc, 2, scale_vec, "/"), center = center, scale = scale_vec, group = group)
}

unscale_coef_grouped <- function(beta_s, center, scale_vec, intercept_s = 0) {
  beta_orig <- beta_s / scale_vec
  intercept_orig <- intercept_s - sum(beta_s * center / scale_vec)
  list(intercept = intercept_orig, beta = beta_orig)
}

# ============================================================
# 5. MODEL 2 ESTIMATION: GROUP L1 (GGLASSO)
# ============================================================

X_sc <- scale_X_grouped(X, group = groups)

# gglasso with "logit" loss requires y in {-1, 1} format
y_gg <- ifelse(y_bin == 1, 1, -1) 

lss <- cv.gglasso(
  x = X_sc$X, y = y_gg, group = groups, foldid = folds, 
  loss = "logit", lambda = exp(seq(-5, -3, length.out = 100)), intercept = TRUE
)

# Extracting coefficients on the standardized scale
beta_s <- coef(lss$gglasso, s = lss$lambda.min)[-1]                  
intercept_s <- coef(lss$gglasso, s = lss$lambda.min)[1]

# Unscaling to the original dimension
back_L1 <- unscale_coef_grouped(beta_s, X_sc$center, X_sc$scale, intercept_s)
beta_L1 <- c(back_L1$intercept, back_L1$beta)

# The final coefficients to compare are `beta_L0` and `beta_L1`

# ============================================================
# 6. RESULTS VISUALIZATION 
# Combining Variable Selection and Functional Effect
# ============================================================

# ------------------------------------------------------------
# A. PREPARE FUNCTIONAL EFFECT DATA ("amount")
# ------------------------------------------------------------
focus_var <- "amount"
idx_amt <- which(group_labels[groups] == focus_var)

# Get ordered unique values and their corresponding spline basis
X_func_amt <- X[order(as.vector(predictors$amount)), idx_amt]
X_unq_amt <- unique(X_func_amt)

# Compute linear predictor components (excluding intercept)
y_L0_amt <- as.vector(X_unq_amt %*% beta_L0[-1][idx_amt])
y_L1_amt <- as.vector(X_unq_amt %*% beta_L1[-1][idx_amt])

df_amount <- data.frame(
  Amount = unique(sort(as.vector(predictors$amount))),
  L0 = y_L0_amt,
  L1 = y_L1_amt
)

# Plot for Functional Effect
p_amount <- ggplot(df_amount, aes(x = Amount)) +
  geom_smooth(aes(y = L0, color = "Group L0 (MM)", linetype = "Group L0 (MM)"), method = "loess", se = FALSE, span = 0.5, linewidth = 1.2) +
  geom_smooth(aes(y = L1, color = "Group Lasso (L1)", linetype = "Group Lasso (L1)"), method = "loess", se = FALSE, span = 0.5, linewidth = 1.2) +
  
  scale_color_manual(name = NULL, values = c("Group L0 (MM)" = "black", "Group Lasso (L1)" = "red")) +
  scale_linetype_manual(name = NULL, values = c("Group L0 (MM)" = "solid", "Group Lasso (L1)" = "dashed")) +
  
  coord_cartesian(ylim = c(-1, 2)) +
  labs(x = "Credit Amount", y = expression(f(x)), title = "Functional Effect of 'Amount'") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.25, 0.85),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
    legend.margin = margin(t = 2, r = 6, b = 2, l = 6),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ------------------------------------------------------------
# B. PREPARE VARIABLE SELECTION DATA
# ------------------------------------------------------------
df_coeff <- data.frame(
  Group = group_labels[groups],
  beta_L0 = beta_L0[-1],
  beta_L1 = beta_L1[-1]
)

# Aggregate by group to see which are selected
selection_df <- df_coeff %>%
  group_by(Group) %>%
  summarise(
    `Group L0 (MM)` = any(beta_L0 != 0),
    `Group Lasso (L1)` = any(beta_L1 != 0),
    .groups = "drop"
  ) %>%
  # Filter out variables NOT selected by either method
  filter(`Group L0 (MM)` | `Group Lasso (L1)`) %>%
  # Reshape data for ggplot
  pivot_longer(
    cols = c(`Group L0 (MM)`, `Group Lasso (L1)`), 
    names_to = "Method", 
    values_to = "Selected"
  )

# Plot for Variable Selection 
p_selection <- ggplot(selection_df, aes(x = Method, y = Group, color = Selected, shape = Selected)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("TRUE" = "darkblue", "FALSE" = "grey80")) +
  scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 4)) +
  labs(title = "Selected Predictors", x = "", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major.x = element_blank(), 
    legend.position = "none" 
  )

# ------------------------------------------------------------
# C. COMBINE AND EXPORT PLOTS
# ------------------------------------------------------------
# Side by side: Selection on the left, Functional Effect on the right
final_poster_plot <- p_selection + p_amount + 
  plot_layout(widths = c(1, 1.5)) # The right plot gets slightly more space

print(final_poster_plot)




