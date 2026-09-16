## ============================================================
## prepare_ames.R
##
## Loads Ames Housing (AmesHousing::make_ames()) and builds a
## group-structured design matrix ready for group_l0_mm() /
## cv_group_l0_mm() (gaussian family).
##
## Large-p setting: ~80 raw predictors -> after dummy expansion
## of the many categorical variables, the design matrix has
## several hundred columns organized into ~80 groups (one per
## original variable), good for stress-testing group selection
## at larger p than Boston/Wage.
##
## Group structure:
##   - each CONTINUOUS covariate -> B-spline / natural-spline
##     basis (df = 4); all its basis columns form one group
##   - each FACTOR covariate -> full dummy encoding (drop
##     intercept), one group per variable, group size = number
##     of levels
##   - a handful of identifier / near-constant / leakage columns
##     are dropped explicitly (see drop_cols below)
##
## Requires: install.packages("AmesHousing"); install.packages("splines")
## (splines ships with base R, no install needed)
## ============================================================

library(AmesHousing)
library(splines)

ames <- AmesHousing::make_ames()
str(ames)

## response: Sale_Price
y <- ames$Sale_Price

## drop identifiers / near-constant / response-adjacent columns
## - PID: parcel identifier, not a predictor
## - Sale_Price: the response itself
## - Latitude/Longitude: kept as an optional spline group (set
##   INCLUDE_LATLON <- FALSE to drop and match the "classic" Ames
##   variable set used in most textbook analyses)
INCLUDE_LATLON <- TRUE

drop_cols <- c("PID", "Sale_Price")
if (!INCLUDE_LATLON) drop_cols <- c(drop_cols, "Latitude", "Longitude")

predictors <- ames[, setdiff(names(ames), drop_cols)]

## drop any column that ends up with a single level (would create
## a zero-column / collinear group, as region did for Wage)
is_degenerate <- sapply(predictors, function(v) {
  if (is.factor(v) || is.character(v)) length(unique(v)) < 2 else FALSE
})
if (any(is_degenerate)) {
  cat("Dropping degenerate (single-level) columns:\n")
  print(names(predictors)[is_degenerate])
  predictors <- predictors[, !is_degenerate]
}

build_group_design <- function(data, df_spline = 4, unique_threshold = 5,
                               zero_inflated_threshold = 0.5) {
  X_list <- list()
  group_vec <- integer(0)
  group_labels <- character(0)
  g <- 0L
  
  for (nm in names(data)) {
    v <- data[[nm]]
    print(nm)
    
    ## se una variabile numerica ha meno di 'unique_threshold' valori
    ## unici, trattala come categoriale
    if (is.numeric(v) && length(unique(v)) < unique_threshold) {
      v <- as.factor(v)
    }
    
    if (is.numeric(v)) {
      
      ## variabile zero-inflated: una grossa massa di punti concentrata
      ## su un unico valore (tipicamente 0) più una coda continua sparsa
      ## -> split in (i) indicatore "diverso dalla massa" + (ii) spline
      ##    fittata solo sulla parte non-massa (0 altrove)
      prop_mode <- max(table(v)) / length(v)
      
      if (prop_mode >= zero_inflated_threshold) {
        mode_val <- as.numeric(names(which.max(table(v))))
        is_nonmode <- v != mode_val
        
        ## gruppo 1: indicatore binario "diverso dal valore di massa"
        g <- g + 1L
        ind_col <- matrix(as.numeric(is_nonmode), ncol = 1)
        colnames(ind_col) <- paste0(nm, "_nonzero")
        X_list[[paste0(nm, "_ind")]] <- ind_col
        group_vec <- c(group_vec, g)
        group_labels[g] <- paste0(nm, "_ind")
        
        ## gruppo 2: spline sulla parte non-massa, 0 dove v == mode_val
        df_v <- if (is.list(df_spline) || !is.null(names(df_spline))) {
          if (!is.null(df_spline[[nm]])) df_spline[[nm]] else df_spline[[1]]
        } else {
          df_spline
        }
        n_nonmode <- sum(is_nonmode)
        df_v <- min(df_v, max(1, n_nonmode - 1))  ## evita df_spline > punti disponibili
        
        spline_part <- matrix(0, nrow = length(v), ncol = df_v)
        spline_part[is_nonmode, ] <- splines::ns(v[is_nonmode], df = df_v)
        colnames(spline_part) <- paste0(nm, "_bs", seq_len(df_v))
        
        g <- g + 1L
        X_list[[paste0(nm, "_spl")]] <- spline_part
        group_vec <- c(group_vec, rep(g, ncol(spline_part)))
        group_labels[g] <- paste0(nm, "_spl")
        
        next
      }
      
      ## caso normale: spline su tutta la variabile
      g <- g + 1L
      df_v <- if (is.list(df_spline) || !is.null(names(df_spline))) {
        if (!is.null(df_spline[[nm]])) df_spline[[nm]] else df_spline[[1]]
      } else {
        df_spline
      }
      cols <- splines::ns(v, df = df_v)
      cols <- matrix(cols, ncol = ncol(cols))
      colnames(cols) <- paste0(nm, "_bs", seq_len(ncol(cols)))
      
    } else {
      g <- g + 1L
      v <- droplevels(as.factor(v))
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

apply(predictors, 2, function(x) length(unique(x)))

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

#save.image("ames_fit.Rdata")

##### ANALYSIS #####

source("MM FUNCTIONS.R")

set.seed(1)
K=5
folds=sample(1:K,replace = T,size = length(y))

#group_weights = sqrt(table(groups))#the same of gglasso

lambda_seq = c(0,exp(seq(-15,-12,l=1)))
rho_seq = seq(log(n)/13,log(n)/1,l=10)#devono crescere sennoò cv spegne tutto col warm_start

out=cv_mm_grouped_gaussian(X,y=y,groups=groups,
                           lambda_seq = lambda_seq,
                           rho_seq = rho_seq,best_fit = T
                           ,folds = folds,intercept = T,group_weights = NULL,)

out$best_par
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
load("ames_fit.Rdata")

A=round(out$best_fit$beta,30)
B=round(coef(lss$gglasso.fit,s=lss$lambda.min),30)
View(cbind(colnames(cbind(1,X)),A,B) )

Q=cbind(c("int",group_labels[groups]), A,B)
unique(Q[Q[,2]!=0,1])[-1]
unique(Q[Q[,3]!=0,1])[-1]

focus="Gr_Liv_Area"
X_func=X[order(as.vector(predictors$Gr_Liv_Area)),group_labels[groups]==focus]
plot(unique(X_func)%*% A[-1][group_labels[groups]==focus],
     type="l",xlab=focus,ylab="y")
lines(unique(X_func)%*% B[-1][group_labels[groups]==focus],
      type="l",col="red")
legend("topleft",lty=1,col=c("black","red"),legend = c("L0","L1"))