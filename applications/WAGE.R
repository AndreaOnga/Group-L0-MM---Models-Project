## ============================================================
## prepare_wage.R
##
## Loads Wage (ISLR2) and builds a group-structured design
## matrix ready for group_l0_mm() / cv_group_l0_mm()
## (gaussian family).
##
## Group structure:
##   - each CONTINUOUS covariate (age, year) -> expanded into a
##     B-spline / natural-spline basis (df = 4); all its basis
##     columns form one group
##   - "maritl", "race", "education", "jobclass", "health",
##     "health_ins" (all categorical factors) -> encoded as
##     dummies, one group each
##   - "region" is dropped: in ISLR2::Wage it has a single level
##     (all observations are "2. Middle Atlantic"), so it carries
##     no information and would create a zero-column group
##   - "logwage" is dropped: it is just log(wage), i.e. a
##     transformed copy of the response, and must not appear
##     among the predictors
##
## Requires: install.packages("ISLR2"); install.packages("splines")
## (splines ships with base R, no install needed)
## ============================================================

library(ISLR2)
library(splines)

data("Wage", package = "ISLR2")
str(Wage)

## response: wage (raw wage, in $1000s)
y <- Wage$wage

## drop identifiers / transformed-response / degenerate columns
drop_cols <- c("wage", "logwage", "region")
predictors <- Wage[, setdiff(names(Wage), drop_cols)]

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

save.image("wage_fit.Rdata")

##### ANALYSIS #####
source("MM FUNCTIONS.R")

set.seed(10)
train=sample(1:length(y),length(y)*0.8,replace = F)
train=1:length(y)
y_verifica=y[-train]
X_verifica=X[-train,]
y=y[train]
X=X[train,]
n=length(y)

set.seed(1)
K=5
folds=sample(1:K,replace = T,size = length(y))

#group_weights = sqrt(table(groups))#the same of gglasso

lambda_seq = c(0,exp(seq(-15,-12,l=1)))
rho_seq = seq(log(n)/17,log(n)/14,l=50)#devono crescere sennoò cv spegne tutto col warm_start

out=cv_mm_grouped_gaussian(X,y=y,groups=groups,
                           lambda_seq = lambda_seq,
                           rho_seq = rho_seq,best_fit = T
                           ,folds = folds,intercept = T,group_weights = NULL,)

out$best_par # ...
lambda_seq
rho_seq
cbind(colnames(cbind(1,X)),out$best_fit$beta)

library(gglasso)
lss=cv.gglasso(X,y=y,group=groups,foldid = folds,lambda = exp(seq(-7,-5,l=100)),intercept=T)
plot(lss) #ok
kmin(lss$cvm)
min(out$mean_mse)

eta=drop(cbind(1,X) %*% out$best_fit$beta)
sum((y-eta)^2)
eta=drop(cbind(1,X) %*% coef(lss$gglasso.fit,s=lss$lambda.min))
sum((y-eta)^2)

##### INTERPRETATION ####
load("wage_fit.Rdata")

A=round(out$best_fit$beta,30)
B=round(coef(lss$gglasso.fit,s=lss$lambda.min),30)
View(cbind(colnames(cbind(1,X)),A,B) ) #FICO

Q=cbind(c("int",group_labels[groups]), A,B)
unique(Q[Q[,2]!=0,1])[-1]
unique(Q[Q[,3]!=0,1])[-1]

focus="age"
X_func=X[order(as.vector(predictors$age)),group_labels[groups]==focus]
plot(unique(X_func)%*% A[-1][group_labels[groups]==focus],
     type="l",xlab=focus,ylab="y")
lines(unique(X_func)%*% B[-1][group_labels[groups]==focus],
      type="l",col="red")
legend("topleft",lty=1,col=c("black","red"),legend = c("L0","L1"))

focus="year"
X_func=X[order(as.vector(predictors$year)),group_labels[groups]==focus]
plot(unique(X_func)%*% A[-1][group_labels[groups]==focus],
     type="l",xlab=focus,ylab="y")
lines(unique(X_func)%*% B[-1][group_labels[groups]==focus],
      type="l",col="red")
legend("topleft",lty=1,col=c("black","red"),legend = c("L0","L1"))