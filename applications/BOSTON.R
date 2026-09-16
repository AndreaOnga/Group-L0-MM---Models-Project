## ============================================================
## prepare_boston_housing.R
##
## Loads BostonHousing2 (mlbench) and builds a group-structured
## design matrix ready for group_l0_mm() / cv_group_l0_mm()
## (gaussian family).
##
## Group structure:
##   - each CONTINUOUS covariate -> expanded into a B-spline basis
##     (df = 4); all its basis columns form one group
##   - "chas" (Charles River dummy, binary factor) -> its own
##     single-column group
##   - "rad" (index of accessibility to radial highways) is
##     conventionally treated as categorical (9 discrete values in
##     the original data) -> encoded as dummies, one group
##   - "town", "tract", "lon", "lat" are identifiers / geographic
##     coordinates from the extended mlbench version; town/tract are
##     dropped here, lon/lat are kept as an optional spline group
##     (set INCLUDE_LONLAT <- FALSE to drop them and reproduce the
##     classic 1970s Boston housing variable set)
##
## Requires: install.packages("mlbench"); install.packages("splines")
## (splines ships with base R, no install needed)
## ============================================================


library(mlbench)
library(splines)

data("BostonHousing2", package = "mlbench")
str(BostonHousing2)

INCLUDE_LONLAT <- TRUE  # set TRUE to add lon/lat as a spline group

## response: cmedv = corrected median value of owner-occupied homes
y <- BostonHousing2$cmedv


?BostonHousing2
## drop identifiers / the two response columns (medv is the
## uncorrected version of cmedv - keep only cmedv as the response)
drop_cols <- c("town", "tract", "medv", "cmedv")
if (!INCLUDE_LONLAT) drop_cols <- c(drop_cols, "lon", "lat")

predictors <- BostonHousing2[, setdiff(names(BostonHousing2), drop_cols)]

build_group_design <- function(data, df_spline = 4) {
  X_list <- list()
  group_vec <- integer(0)
  group_labels <- character(0)
  g <- 0L
  
  for (nm in names(data)) {
    v <- data[[nm]]
    g <- g + 1L
    
    if (is.numeric(v)) {
      ## spline basis for this numeric variable, gdl specificabili
      ## df_spline può essere uno scalare (stesso per tutte le var. numeriche)
      ## o un vettore/lista nominata con gdl diversi per variabile
      df_v <- if (is.list(df_spline) || !is.null(names(df_spline))) {
        if (!is.null(df_spline[[nm]])) df_spline[[nm]] else df_spline[[1]]
      } else {
        df_spline
      }
      cols <- splines::ns(v, df = df_v)
      cols <- matrix(cols, ncol = ncol(cols))
      colnames(cols) <- paste0(nm, "_bs", seq_len(ncol(cols)))
    } else {
      v <- droplevels(as.factor(v))
      ## full dummy encoding for this factor alone (drop intercept)
      cols <- model.matrix(~ v - 1)
      colnames(cols) <- paste0(nm, "_", levels(v))
    }
    
    X_list[[nm]] <- cols
    group_vec <- c(group_vec, rep(g, ncol(cols)))
    group_labels[g] <- nm
  }
  
  X <- do.call(cbind, X_list)
  list(X = X, group = group_vec, group_labels = group_labels)
}

built <- build_group_design(predictors)

str(predictors)
X <- built$X
groups <- built$group
group_labels <- built$group_labels

cat(sprintf("\nDesign matrix: %d observations x %d columns, %d groups\n",
            nrow(X), ncol(X), length(group_labels)))

data.frame(group_id = seq_along(group_labels),
           label = group_labels,
           n_columns = as.integer(table(factor(groups, levels = seq_along(group_labels)))))

save.image("boston_fit.Rdata")
##### ANALYSIS #####

source("MM FUNCTIONS.R")



n=length(y)
set.seed(1)
K=5
folds=sample(1:K,replace = T,size = length(y))

#group_weights = sqrt(table(groups))#the same of gglasso
lambda_seq = c(0,exp(seq(-15,-12,l=1)))
rho_seq = seq(log(n)/17,log(n)/14,l=50)
#devono crescere sennoò cv spegne tutto col warm_start

out=cv_mm_grouped_gaussian(X,y=y,groups=groups,
                           lambda_seq = lambda_seq,
                           rho_seq = rho_seq,best_fit = T
                           ,folds = folds,intercept = T,group_weights = NULL,
)

out$best_par # 0 0.3726739
lambda_seq
rho_seq

cbind(colnames(cbind(1,X)),out$best_fit$beta)
library(gglasso)

lss=cv.gglasso(X,y=y,group=groups,foldid = folds,lambda = exp(seq(-7,-5,l=100)),intercept=T)


plot(lss) #okk


min(lss$cvm)
min(out$mean_mse)

eta=drop(cbind(1,X) %*% out$best_fit$beta)

sum((y-eta)^2)

eta=drop(cbind(1,X) %*% coef(lss$gglasso.fit,s=lss$lambda.min))

sum((y-eta)^2)


##### INTERPETRATION ####
load("boston_fit.Rdata")

A=round(out$best_fit$beta,30)
B=round(coef(lss$gglasso.fit,s=lss$lambda.min),30)

View(cbind(colnames(cbind(1,X)),A,B) ) #FICO

Q=cbind(c("int",group_labels[groups]), A,B)

unique(Q[Q[,2]!=0,1])[-1]

unique(Q[Q[,3]!=0,1])[-1]

focus="zn"

X_func=X[order(as.vector(predictors$zn)),group_labels[groups]==focus]




plot(unique(X_func)%*% A[-1][group_labels[groups]==focus],
     type="l")
lines(unique(X_func)%*% B[-1][group_labels[groups]==focus],
      type="l",col="red")

focus="lon"

X_func=X[order(as.vector(predictors$lon)),group_labels[groups]==focus]

plot(unique(X_func)%*% A[-1][group_labels[groups]==focus],
     type="l",xlab=focus,ylab="y")
lines(unique(X_func)%*% B[-1][group_labels[groups]==focus],
      type="l",col="red")
legend("topleft",lty=1,col=c("black","red"),legend = c("L0","L1"))



