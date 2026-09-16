# ============================================================
# 1. "STUDENT PERFORMANCE" DATA PREPARATION
# ============================================================
library(splines)
library(gglasso)
library(dplyr)
library(tidyr)

# Download and load data directly from UCI repository
temp <- tempfile()
download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/00320/student.zip", temp)
student <- read.table(unz(temp, "student-mat.csv"), sep = ";", header = TRUE, stringsAsFactors = TRUE)
unlink(temp)

# Response variable: G3 (final math grade, 0 to 20)
y <- student$G3

# Exclude intermediate grades (G1, G2) to evaluate purely demographic/social predictors
drop_cols <- c("G1", "G2", "G3")
predictors <- student[, setdiff(names(student), drop_cols)]

# Remove degenerate columns (only one unique value)
is_degenerate <- sapply(predictors, function(v) length(unique(v)) < 2)
if (any(is_degenerate)) {
  predictors <- predictors[, !is_degenerate]
}

# ============================================================
# 2. DESIGN MATRIX BUILDING (Handling Zero-Inflated variables)
# ============================================================
build_group_design <- function(data, df_spline = 4, unique_threshold = 6, zero_inflated_threshold = 0.29) {
  X_list <- list()
  group_vec <- integer(0)
  group_labels <- character(0)
  g <- 0L
  
  for (nm in names(data)) {
    v <- data[[nm]]
    
    # Treat numeric variables with few unique values as categorical
    if (is.numeric(v) && length(unique(v)) < unique_threshold) v <- as.factor(v)
    
    if (is.numeric(v)) {
      prop_mode <- max(table(v)) / length(v)
      
      # Handle Zero-Inflated Variables: Split into binary indicator + spline for non-mass part
      if (prop_mode >= zero_inflated_threshold) {
        mode_val <- as.numeric(names(which.max(table(v))))
        is_nonmode <- v != mode_val
        
        # Group 1: Binary indicator "different from mass value"
        g <- g + 1L
        ind_col <- matrix(as.numeric(is_nonmode), ncol = 1)
        colnames(ind_col) <- paste0(nm, "_nonzero")
        X_list[[paste0(nm, "_ind")]] <- ind_col
        group_vec <- c(group_vec, g)
        group_labels[g] <- paste0(nm, "_ind")
        
        # Group 2: Spline on the non-mass part
        n_nonmode <- sum(is_nonmode)
        df_v <- min(df_spline, max(1, n_nonmode - 1))
        spline_part <- matrix(0, nrow = length(v), ncol = df_v)
        spline_part[is_nonmode, ] <- splines::ns(v[is_nonmode], df = df_v)
        colnames(spline_part) <- paste0(nm, "_bs", seq_len(df_v))
        
        g <- g + 1L
        X_list[[paste0(nm, "_spl")]] <- spline_part
        group_vec <- c(group_vec, rep(g, ncol(spline_part)))
        group_labels[g] <- paste0(nm, "_spl")
        next
      }
      
      # Standard case: Spline on the whole variable
      g <- g + 1L
      cols <- matrix(splines::ns(v, df = df_spline), ncol = df_spline)
      colnames(cols) <- paste0(nm, "_bs", seq_len(df_spline))
      
    } else {
      # Categorical: Full dummy encoding (drop intercept)
      g <- g + 1L
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

built <- build_group_design(predictors, unique_threshold = 6)
X <- built$X
groups <- built$group
group_labels <- built$group_labels

# ============================================================
# 3. SETUP & CROSS-VALIDATION
# ============================================================
# source("CROSS VALIDATION FUNCTIONS.R") # Ensure this is loaded

n <- length(y)
set.seed(1)
K <- 5
folds <- sample(1:K, replace = TRUE, size = n)

# ============================================================
# 4. MODEL 1: GROUP L0 (MM GAUSSIAN)
# ============================================================
lambda_seq <- c(0)
rho_seq <- seq(log(n)/0.1, log(n)/0.05, length.out = 50)

out_L0 <- cv_mm_grouped_gaussian(
  X = X, y = y, groups = groups,
  lambda_seq = lambda_seq, rho_seq = rho_seq, best_fit = TRUE,
  folds = folds, intercept = TRUE, group_weights = NULL
)

beta_L0 <- out_L0$best_fit$beta

# ============================================================
# 5. MODEL 2: GROUP L1 (GGLASSO) & SCALING
# ============================================================
scale_X_grouped <- function(X, group) {
  n <- nrow(X)
  center <- colMeans(X)
  Xc <- sweep(X, 2, center, "-")
  scale_vec <- numeric(ncol(X))
  for (g in unique(group)) {
    idx <- which(group == g)
    s_g <- sqrt(sum(Xc[, idx]^2) / (n * length(idx)))
    scale_vec[idx] <- ifelse(s_g == 0, 1, s_g)
  }
  list(X = sweep(Xc, 2, scale_vec, "/"), center = center, scale = scale_vec, group = group)
}

unscale_coef_grouped <- function(beta_s, center, scale_vec, intercept_s = 0) {
  beta_orig <- beta_s / scale_vec
  intercept_orig <- intercept_s - sum(beta_s * center / scale_vec)
  list(intercept = intercept_orig, beta = beta_orig)
}

X_sc <- scale_X_grouped(X, group = groups)

# gglasso defaults to Gaussian loss
lss <- cv.gglasso(
  X_sc$X, y = y, group = groups, foldid = folds, 
  lambda = exp(seq(-2, -1, length.out = 50)), intercept = TRUE
)

beta_s <- coef(lss$gglasso, s = lss$lambda.min)[-1]                  
intercept_s <- coef(lss$gglasso, s = lss$lambda.min)[1]

back_L1 <- unscale_coef_grouped(beta_s, X_sc$center, X_sc$scale, intercept_s)
beta_L1 <- c(back_L1$intercept, back_L1$beta)

# ============================================================
# 7. GENERATING TABLE METRICS (L2 NORMS & CV MSE)
# ============================================================

df_coeff <- data.frame(
  Group   = group_labels[groups],
  beta_L0 = beta_L0[-1],
  beta_L1 = beta_L1[-1]
)

#  Compute group-wise L2 norms and selection status
group_summary <- df_coeff %>%
  group_by(Group) %>%
  summarise(
    L2_Norm_L0  = sqrt(sum(beta_L0^2)),
    L2_Norm_L1  = sqrt(sum(beta_L1^2)),
    Selected_L0 = any(beta_L0 != 0),
    Selected_L1 = any(beta_L1 != 0),
    .groups     = "drop"
  )

group_summary_ordered <- group_summary %>%
  arrange(desc(Selected_L0), desc(Selected_L1), desc(L2_Norm_L0))

final_table_data <- data.frame(
  Group       = group_summary_ordered$Group,
  L2_Norm_L0  = round(group_summary_ordered$L2_Norm_L0, 2),
  L2_Norm_L1  = round(group_summary_ordered$L2_Norm_L1, 2)
)
print(final_table_data)

# 4. Extract Cross-Validation MSE for both models
cv_mse_L0 <- min(out_L0$mean_mse)
cv_mse_L1 <- min(lss$cvm)
