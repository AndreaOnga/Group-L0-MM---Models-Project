# ============================================================
## prepare_student.R
##
## Loads the UCI Student Performance Data Set (Math course) and 
## builds a group-structured design matrix ready for 
## group_l0_mm() / cv_group_l0_mm() (gaussian family).
##
## Group structure:
##   - continuous covariates -> natural-spline basis (df = 4)
##   - factor covariates -> full dummy encoding (drop intercept)
## ============================================================

library(splines)

DAT=gglasso::bardet

y=DAT$y

X <- DAT$x
groups <- sapply(1:20, function(x) rep(x,5))#20 geni 5 b_splines
group_labels <-1:20


scale_X_grouped <- function(X, group) {
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(length(group) == p)
  
  center <- colMeans(X)
  Xc <- sweep(X, 2, center, "-")
  
  group_ids <- unique(group)
  scale_vec <- numeric(p)
  
  for (g in group_ids) {
    idx <- which(group == g)
    Xg <- Xc[, idx, drop = FALSE]
    s_g <- sqrt(sum(Xg^2) / (n * length(idx)))
    if (s_g == 0) s_g <- 1  # evita divisione per zero (gruppo costante)
    scale_vec[idx] <- s_g
  }
  
  Xs <- sweep(Xc, 2, scale_vec, "/")
  
  list(X = Xs, center = center, scale = scale_vec, group = group)
}

## ============================================================
## unscale_coef_grouped() : riporta i coefficienti stimati sul
##                           modello con X standardizzata a gruppi
##                           alla scala originale di X.
##
## Modello fittato:   y = beta0 + Xs %*% beta_s + eps
## dove Xs = (X - center) / scale_vec  (scale_vec costante entro
## ogni gruppo, per costruzione di scale_X_grouped)
##
## Stessa algebra di prima: la formula è invariata rispetto al
## caso colonna-per-colonna, perché scale_vec è comunque un
## vettore lungo p (solo con valori ripetuti entro gruppo).
## ============================================================
unscale_coef_grouped <- function(beta_s, center, scale_vec, intercept_s = 0) {
  beta_orig <- beta_s / scale_vec
  intercept_orig <- intercept_s - sum(beta_s * center / scale_vec)
  
  list(intercept = intercept_orig, beta = beta_orig)
}

##### ANALYSIS #####
source("MM FUNCTIONS.R")

n <- length(y)

set.seed(1)
K <- 5
folds <- sample(1:K, replace = TRUE, size = length(y))

lambda_seq <- c(0)
rho_seq <- seq(log(n)/5, log(n)/2, l = 50)


X_sc <- scale_X_grouped(X,group = groups)

out <- cv_mm_grouped_gaussian(X, y = y, groups = groups,
                              lambda_seq = lambda_seq,
                              rho_seq = rho_seq, best_fit = TRUE,
                              folds = folds, intercept = TRUE, group_weights = NULL,
                              max_iter_out = 500,max_iter_in = 500)



out$best_par #### okk
lambda_seq
rho_seq
cbind(colnames(cbind(1, X)), out$best_fit$beta)

#beta_s <- out$best_fit$beta[-1]                   
#intercept_s <- out$best_fit$beta[1]

#backL0 <- unscale_coef_grouped(beta_s, X_sc$center, X_sc$scale, intercept_s)
#backL0=c(backL0$intercept,backL0$beta)

backL0=out$best_fit$beta

library(gglasso)

X_sc <- scale_X_grouped(X,group = groups)

lss <- cv.gglasso(X_sc$X, y = y, group = groups, foldid = folds, lambda = exp(seq(-7, -2, l = 50)), intercept = TRUE)

plot(lss)

beta_s <- coef(lss$gglasso,s=lss$lambda.min)[-1]                   
intercept_s <- coef(lss$gglasso,s=lss$lambda.min)[1]

backL1 <- unscale_coef_grouped(beta_s, X_sc$center, X_sc$scale, intercept_s)
backL1=c(backL1$intercept,backL1$beta)

plot(lss)
min(lss$cvm) ### ok (ho zoomato dentro)
min(out$mean_mse)

eta <- cbind(1, X) %*% backL0
sum((y - eta)^2)

plot(y,ylim=range(c(y,eta)))
points(eta)

eta <- drop(cbind(1, X) %*% backL1)
sum((y - eta)^2)

plot(y)
points(eta)


##### INTERPRETATION ####


A <- round(backL0, 30)
B <- round(backL1, 30)
View(cbind(colnames(cbind(1, X)), A, B))

Q <- cbind(c("int", group_labels[groups]), A, B)
unique(Q[Q[, 2] != 0, 1])[-1]
unique(Q[Q[, 3] != 0, 1])[-1]


##### SELECTION #####
View(cbind(colnames(cbind(1,X)),A,B) ) #FICO

Q=cbind(c("int",group_labels[groups]), A,B)

U=group_labels

unique(Q[Q[,2]!=0,1])[-1]

unique(Q[Q[,3]!=0,1])[-1]

library(dplyr)

# 1. Creazione dataframe escludendo l'intercetta
df_coeff <- data.frame(
  Group  = group_labels[groups],
  beta_A = A[-1],
  beta_B = B[-1]
)

# 2. Aggregazione per gruppo: Selezione (TRUE/FALSE) e Norma L2
tabella_gruppi <- df_coeff %>%
  group_by(Group) %>%
  summarise(
    Selected_A = any(beta_A != 0),
    L2_Norm_A  = sqrt(sum(beta_A^2)),
    Selected_B = any(beta_B != 0),
    L2_Norm_B  = sqrt(sum(beta_B^2)),
    .groups    = "drop"
  )

# Per visualizzarla
View(tabella_gruppi)

tabella_gruppi[tabella_gruppi$Selected_A | tabella_gruppi$Selected_B,]
dim(tabella_gruppi)

# Metodo con order()
sub_tabella <- tabella_gruppi[tabella_gruppi$Selected_A | tabella_gruppi$Selected_B, ]
risultato <- sub_tabella[order(-sub_tabella$Selected_A, -sub_tabella$Selected_B), ]

# Oppure unendo le due selezioni con rbind()
risultato <- rbind(
  tabella_gruppi[tabella_gruppi$Selected_A, ],
  tabella_gruppi[!tabella_gruppi$Selected_A & tabella_gruppi$Selected_B, ]
)

View(risultato)







# Variabile target per l'effetto funzionale (es. engine_size, molto correlata al prezzo)
focus <- 4
X_func <- X[, group_labels[groups] == focus]
plot(unique(X_func) %*% A[-1][group_labels[groups] == focus],
     type = "l", xlab = focus, ylab = "y")

lines(unique(X_func) %*% B[-1][group_labels[groups] == focus],
      type = "l", col = "red")
legend("topleft", lty = 1, col = c("black", "red"), legend = c("L0", "L1"))


# Variabile target per l'effetto funzionale (es. engine_size, molto correlata al prezzo)
focus <- 15
X_func <- X[, group_labels[groups] == focus]
plot(unique(X_func) %*% A[-1][group_labels[groups] == focus],
     type = "l", xlab = focus, ylab = "y")

lines(unique(X_func) %*% B[-1][group_labels[groups] == focus],
      type = "l", col = "red")
legend("topleft", lty = 1, col = c("black", "red"), legend = c("L0", "L1"))
