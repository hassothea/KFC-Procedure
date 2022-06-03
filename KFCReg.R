#### -------------------------------------------------------------------------- ####
#### -------------------- KFC procedure for regression ------------------------ ####
#### -------------------------------------------------------------------------- ####


# Check if package "fontawesome" is already installed 

lookup_packages <- installed.packages()[,1]
if(!("fontawesome" %in% lookup_packages))
  install.packages("fontawesome")

# ---------------------------------------------------------------------------------#

# Check if package "pacman" is already installed 

lookup_packages <- installed.packages()[,1]
if(!("pacman" %in% lookup_packages))
  install.packages("pacman")


# To be installed or loaded
pacman::p_load(magrittr)
pacman::p_load(tidyverse)

## package for "generateMachines"
pacman::p_load(tree)
pacman::p_load(glmnet)
pacman::p_load(randomForest)
pacman::p_load(FNN)
pacman::p_load(xgboost)
pacman::p_load(keras)
pacman::p_load(pracma)
pacman::p_load(latex2exp)
pacman::p_load(plotly)
pacman::p_load(parallel)
pacman::p_load(foreach)
pacman::p_load(doParallel)
rm(lookup_packages)


# ---------------------------------------------------------------------------------#

euclidDiv <- function(X., y., deg = NULL){
  res <- sweep(X., 2, y.)
  return(rowSums(res^2))
}
gklDiv <- function(X., y., deg = NULL){
  res <- c("/", "-") %>%
    map(.f = ~ sweep(X., 2, y., FUN = .x))
  return(rowSums(X.*log(res[[1]]) - res[[2]]))
}
logDiv <- function(X., y., deg = NULL){
  res <-  map2(.x = list(X., 1-X.),
               .y = list(y., 1-y.),
               .f = ~ sweep(.x, 2, .y, FUN = "/"))
  return(rowSums(X.*log(res[[1]])+(1-X.)*log(res[[2]])))
}
itaDiv <- function(X., y., deg = NULL){
  res <- sweep(X., 2, y., FUN = "/")
  return(rowSums(res-log(res) - 1))
}
expDiv <- function(X., y., deg = NULL){
  exp_y <- exp(y.)
  res <- sweep(1+X., 2, y.) %>%
    sweep(2, exp_y, FUN = "*")
  return(rowSums(exp(X.)-res))
}
polyDiv <- function(X., y., deg = 3){
  S <- map2(.x = list(X., X.^deg),
            .y = list(y., y.^deg),
            .f = ~ sweep(.x, 
                         MARGIN = 2, 
                         STATS = .y,
                         FUN = "-"))
  if(deg %% 2 == 0){
    Tem <- sweep(S[[1]], 2, y.^(deg-1), FUN = "*")
    res <- rowSums(S[[2]] - deg * Tem)
  }
  else{
    Tem <- sweep(S[[1]], 2, sign(y.) * y.^(deg-1), FUN = "*")
    res <- rowSums(S[[2]] - deg * Tem)
  }
  return(res)
}
lookup_div <- list(
  euclidean = euclidDiv,
  gkl = gklDiv,
  logistic = logDiv,
  itakura = itaDiv,
  exponential = expDiv,
  polynomial = polyDiv
)

# ---------------------------------------------------------------------------------#

BregmanDiv <- function(X., 
                       C., 
                       div = c("euclidean",
                               "gkl",
                               "logistic",
                               "itakura",
                               "exponential",
                               "polynomial"),
                       deg = 3){
  div <- match.arg(div)
  d_c <- dim(C.)
  if(is.null(d_c)){
    C <- matrix(C., nrow = 1, byrow = TRUE)
  } else{
    C <- as.matrix(C.)
  }
  if(is.null(dim(X.))){
    X <- matrix(X., nrow = 1, byrow = TRUE)
  } else{
    X <- as.matrix(X.)
  }
  dis <-  map_dfc(.x = 1:dim(C)[1],
                  .f = ~ tibble('{{.x}}' := lookup_div[[div]](X, C[.x,], deg = deg)))
  return(dis)
}

# ---------------------------------------------------------------------------------#

findClosestCentroid <- function(x., centroids., div, deg = 3){
  dist <- BregmanDiv(x., centroids., div, deg)
  clust <- 1:nrow(x.) %>%
    map_int(.f = ~ which.min(dist[.x,]))
  return(clust)
}
newCentroids <- function(x., clusters.){
  centroids <- unique(clusters.) %>%
    map_dfr(.f = ~ colMeans(x.[clusters. == .x, ]))
  return(centroids)
}

# ---------------------------------------------------------------------------------#

kmeansBD <- function(train_input,
                     K,
                     n_start = 5,
                     maxIter = 500,
                     deg = 3,
                     scale_input = FALSE,
                     div = "euclidean",
                     splits = 1,
                     epsilon = 1e-10,
                     center_ = NULL,
                     scale_ = NULL,
                     id_shuffle = NULL){
  start_time <- Sys.time()
  # Distortion function
  X <- as.matrix(train_input)
  N <- dim(X)
  if(scale_input){
    if(!(is.null(center_) & is.null(scale_))){
      if(length(center_) == 1){
        center_ <- rep(center_, N[2])
      }
      if(length(scale_) == 1){
        scale_ <- rep(scale_, N[2])
      }
    } else{
      min_ <- apply(X, 2, FUN = min)
      c_ <- abs(colMeans(X)/5)
      center_ <- min_ - c_
      scale_ <- apply(X, 2, FUN = max) - center_ + 1
    }
    X <- scale(X, center = center_, scale = scale_)
  }
  if(is.null(id_shuffle)){
    train_id <- rep(TRUE, N[1])
    if(splits < 1){
      train_id[sample(N[1], floor(N[1]*(1-splits)))] <- FALSE
    }
  } else{
    train_id <- id_shuffle
  }
  X_train1 <- X[train_id,]
  X_train2 <- X[!train_id,]
  mu <- as.matrix(colMeans(X_train1))
  distortion <- function(clus){
    cent <- newCentroids(X_train1, clus)
    var_within <- 1:K %>%
      map(.f = ~ BregmanDiv(X_train1[clus == .x,], 
                            cent[.x,], 
                            div, 
                            deg)) %>%
      map(.f = sum) %>%
      Reduce("+", .)
    return(var_within)
  }
  # Kmeans algorithm
  kmeansWithBD <- function(x., k., maxiter., eps.) {
    n. <- nrow(x.)
    # initialization
    init <- sample(n., k.)
    centroids_old <- x.[init,]
    i <- 0
    while(i < maxIter){
      # Assignment step
      clusters <- findClosestCentroid(x., centroids_old, div, deg)
      # Recompute centroids
      centroids_new <- newCentroids(x., clusters)
      if ((sum(is.na(centroids_new)) > 0) |
          (nrow(centroids_new) != k.)) {
        init <- sample(n., k.)
        centroids_old <- x.[init,]
        warning("NA produced -> reinitialize centroids...!")
      }
      else{
        if(sum(abs(centroids_new - centroids_old)) > eps.){
          centroids_old <- centroids_new
        } else{
          break
        }
      }
      i <- i + 1
    }
    return(clusters)
  }
  results <- 1:n_start %>% 
    map_dfc(.f = ~ tibble("{{.x}}" := kmeansWithBD(X_train1, 
                                                   K,
                                                   maxIter, 
                                                   epsilon)))
  opt_id <- 1:n_start %>%
    map_dbl(.f = ~ distortion(results[[.x]])) %>%
    which.min
  cluster <- clusters <- results[[opt_id]]
  j <- 1
  ID <- unique(cluster)
  for (i in ID) {
    clusters[cluster == i] = j
    j =  j + 1
  }
  centroids = newCentroids(X_train1, clusters)
  time_taken <- Sys.time() - start_time
  return(
    list(
      centroids = centroids,
      clusters = clusters,
      train_data = list(X_train = X_train1,
                        X_remain = X_train2,
                        id_remain = !train_id),
      parameters = list(div = div,
                        deg = deg,
                        center_ = center_,
                        scale_ = scale_),
      running_time = time_taken
    )
  )
}

# ---------------------------------------------------------------------------------#

fitLocalModels <- function(kmeans_BD,
                           train_response,
                           model = "lm",
                           formula = NULL){
  start_time <- Sys.time()
  X_train <- kmeans_BD$train_data$X_train
  y_train <- train_response[!(kmeans_BD$train_data$id_remain)]
  X_remain <- kmeans_BD$train_data$X_remain
  y_remain <- NULL
  if(!is.null(X_remain)){
    y_remain <- train_response[kmeans_BD$train_data$id_remain]
  }
  pacman::p_load(tree)
  pacman::p_load(randomForest)
  model_ <- ifelse(model == "tree", tree::tree, model)
  K <- nrow(kmeans_BD$centroids)
  if (is.null(formula)){
    form <- formula(target ~ .)
  }
  else{
    form <- update(formula, target ~ .)
  }
  data_ <- bind_cols(X_train, "target":= y_train)
  fit_lookup <- list(lm = "fitted.values",
                     rf = "predicted")
  if(is.character(model_)){
    model_lookup <- list(lm = lm,
                         rf = randomForest::randomForest)
    mod <- map(.x = 1:K, 
               .f = ~ model_lookup[[model_]](formula = form, 
                                             data = data_[kmeans_BD$clusters == .x, ]))
  } else{
    mod <- map(.x = 1:K, 
               .f = ~ model_(formula = form, 
                             data = data_[kmeans_BD$clusters == .x,]))
  }
  pred0 <- NULL
  if(!is.null(X_remain)){
    pred0 <- vector(mode = "numeric", 
                    length = length(y_remain))
    clus <- findClosestCentroid(x. = X_remain,
                                centroids. = kmeans_BD$centroids,
                                div = kmeans_BD$parameters$div,
                                deg = kmeans_BD$parameters$deg)
    for(i_ in 1:K){
      pred0[clus == i_] <- predict(mod[[i_]],
                                   as.data.frame(X_remain[clus == i_,]))
    }
  }
  time_taken <- Sys.time() - start_time
  return(list(
    local_models = mod,
    kmeans_BD = kmeans_BD,
    data_remain = list(fit = pred0,
                       response = y_remain),
    running_time = time_taken
  ))
}

# ---------------------------------------------------------------------------------#

localPredict <- function(localModels,
                         newData){
  kmean_BD <- localModels$kmeans_BD
  K <- nrow(kmean_BD$centroids)
  newData_ <- newData
  if(!(is.null(kmean_BD$parameters$center_))){
    newData_ <- scale(newData,
                      center = kmean_BD$parameters$center_,
                      scale = kmean_BD$parameters$scale_)
    id0 <- (newData_ <= 0)
    if(sum(id0) > 0){
      min_ <- min(newData_[id0])
      newData_[id0] <- runif(sum(id0), min(1e-3, min_/10), min_)
    }
  }
  clus <- findClosestCentroid(x. = newData_,
                              centroids. = kmean_BD$centroids,
                              div = kmean_BD$parameters$div,
                              deg = kmean_BD$parameters$deg)
  pred0 <- vector(mode = "numeric", length = nrow(newData_))
  for(i_ in 1:K){
    pred0[clus == i_] <- predict(localModels$local_models[[i_]],
                                 as.data.frame(newData_[clus == i_,]))
  }
  pred0 <- as_tibble(pred0)
  names(pred0) <- ifelse(kmean_BD$parameters$div == "polynomial",
                         paste0("polynomial", kmean_BD$parameters$deg),
                         kmean_BD$parameters$div)
  return(pred0)
}



# ---------------------------------------------------------------------------------#

pacman::p_load(devtools)
### Kernel based consensual aggregation
source_url("https://raw.githubusercontent.com/hassothea/AggregationMethods/main/KernelAggReg.R")
### MixCobra
source_url("https://raw.githubusercontent.com/hassothea/AggregationMethods/main/MixCobraReg.R")

# ---------------------------------------------------------------------------------#


stepK = function(K,
                 n_start = 5,
                 maxIter = 300,
                 deg = 3,
                 scale_input = FALSE,
                 div = NULL,
                 splits = 0.75,
                 epsilon = 1e-10,
                 center_ = NULL,
                 scale_ = NULL){
  return(list(K = K,
              n_start = n_start,
              maxIter = maxIter,
              deg = deg,
              scale_input = scale_input,
              div = div,
              splits = splits,
              epsilon = epsilon,
              center_ = center_ ,
              scale_ = scale_))
}

stepF = function(formula = NULL, 
                 model = "lm"){
  return(list(formula = formula, 
              model = model))
}

stepC = function(n_cv = 5,
                 method = c("cobra", "mixcobra"),
                 opt_methods = c("grad", "grid"),
                 kernels = "gaussian",
                 scale_features = FALSE){
  return(list(n_cv = n_cv,
              method = method,
              opt_methods = opt_methods,
              kernels = kernels,
              scale_features = scale_features))
}


# ---------------------------------------------------------------------------------#


KFCreg = function(train_input,
                  train_response,
                  test_input,
                  test_response = NULL,
                  n_cv = 5,
                  parallel = TRUE,
                  inv_sigma = sqrt(.5),
                  alp = 2,
                  K_step = stepK(splits = .5),
                  F_step = stepF(),
                  C_step = stepC(),
                  setGradParamAgg = setGradParameter(),
                  setGridParamAgg = setGridParameter(),
                  setGradParamMix = setGradParameter_Mix(),
                  setGridParamMix = setGridParameter_Mix()){
  start_time <- Sys.time()
  lookup_div_names <- c("euclidean",
                        "gkl",
                        "logistic",
                        "itakura",
                        "polynomial")
  div_ <- K_step$div
  ### K step: Kmeans clustering with BDs
  if (is.null(K_step$div)){
    divergences <- lookup_div_names
    warning("No divergence provided! All of them are used!")
  }
  else{
    divergences <- K_step$div %>% 
      map_chr(.f = ~ match.arg(arg = .x, 
                               choices = lookup_div_names))
  }
  div_list <- divergences %>% 
    map(.f = (\(x) if(x != "polynomial") return(x) else return(rep("polynomial", length(K_step$deg))))) %>%
    unlist
  deg_list <- rep(NA, length(div_))
  deg_list[div_list == "polynomial"] <- K_step$deg
  div_names <- map2_chr(.x = div_list,
                        .y = deg_list,
                        .f = (\(x, y) if(is.na(y)) return(x) else return(paste0(x,y))))
  ### Step K: Kmeans clustering with Bregman divergences
  dm <- dim(train_input)
  id_shuffle <- vector(length = dm[1])
  n_train <- floor(K_step$splits * dm[1])
  id_shuffle[sample(dm[1], n_train)] <- TRUE
  if(parallel){
    numCores <- parallel::detectCores()
    doParallel::registerDoParallel(numCores) # use multicore, set to the number of our cores
    kmean_ <- foreach(i=1:length(div_names)) %dopar% {
      devtools::source_url("https://raw.githubusercontent.com/hassothea/KFC-Procedure/master/kmeanBD.R")
      kmeansBD(train_input = train_input,
               K = K_step$K,
               div = div_list[i],
               n_start = K_step$n_start,
               maxIter = K_step$maxIter,
               deg = deg_list[i],
               scale_input = K_step$scale_input,
               splits = K_step$splits,
               epsilon = K_step$epsilon,
               center_ = K_step$center_,
               scale_ = K_step$scale_,
               id_shuffle = id_shuffle)
    }
    doParallel::stopImplicitCluster()
  } else{
    kmean_ <- map2(.x = div_list,
                   .y = deg_list,
                   .f = ~ kmeansBD(train_input = train_input,
                                   K = K_step$K,
                                   div = .x,
                                   n_start = K_step$n_start,
                                   maxIter = K_step$maxIter,
                                   deg = .y,
                                   scale_input = K_step$scale_input,
                                   splits = K_step$splits,
                                   epsilon = K_step$epsilon,
                                   center_ = K_step$center_,
                                   scale_ = K_step$scale_,
                                   id_shuffle = id_shuffle))
  }
  names(kmean_) <- div_names
  ### F step: Fitting the corresponding model on each observed cluster
  model_ <- div_names %>%
    map(.f = ~ fitLocalModels(kmean_[[.x]],
                              train_response = train_response,
                              model = F_step$model,
                              formula = F_step$formula))
  names(model_) <- div_names
  pred_combine <- model_ %>%
    map_dfc(.f = ~ .x$data_remain$fit)
  y_remain <- train_response[!id_shuffle]
  pred_test <- div_names %>%
    map_dfc(.f = ~ localPredict(model_[[.x]],
                                test_input))
  names(pred_test) <- names(pred_combine) <- div_names
  # C step: Consensual regression aggregation method with kernel-based COBRA
  list_method_agg <- list(mixcobra = function(pred){MixCobraReg(train_input = train_input[!id_shuffle,],
                                                                train_response = y_remain,
                                                                test_input = test_input,
                                                                train_predictions = pred,
                                                                test_predictions = pred_test,
                                                                test_response = test_response,
                                                                scale_input = K_step$scale_input,
                                                                scale_machine = C_step$scale_features,
                                                                n_cv = C_step$n_cv,
                                                                inv_sigma = inv_sigma,
                                                                alp = alp,
                                                                kernels = C_step$kernels,
                                                                optimizeMethod = C_step$opt_methods,
                                                                setGradParam = setGradParamMix,
                                                                setGridParam = setGridParamMix)},
                          cobra = function(pred){kernelAggReg(train_design = pred,
                                                              train_response = y_remain,
                                                              test_design = pred_test,
                                                              test_response = test_response,
                                                              scale_input = K_step$scale_input,
                                                              scale_machine = C_step$scale_features,
                                                              build_machine = FALSE,
                                                              machines = NULL,
                                                              n_cv = C_step$n_cv,
                                                              inv_sigma = sqrt(.5),
                                                              alp = 2,
                                                              kernels = C_step$kernels,
                                                              optimizeMethod = C_step$opt_methods,
                                                              setGradParam = setGradParamAgg,
                                                              setGridParam = setGridParamAgg)})
  res <- map(.x = C_step$method,
             .f = ~ list_method_agg[[.x]](pred_combine))
  list_agg_methods <- list(cobra = "cob",
                           mixcobra = "mix")
  names(res) <- C_step$method
  ext_fun <- function(L, nam){
    tab <- L$fitted_aggregate
    names(tab) <- paste0(names(tab), "_", nam)
    return(tab)
  }
  pred_fin <- C_step$method %>%
    map_dfc(.f = ~ ext_fun(res[[.x]], list_agg_methods[[.x]]))
  time.taken <- Sys.time() - start_time
  ### To finish
  if(is.null(test_response)){
    return(list(
      predict_final = pred_fin,
      predict_local = pred_test,
      agg_method = res,
      running_time = time.taken
    ))
  } else{
    error <- cbind(pred_test, pred_fin) %>%
      dplyr::mutate(y_test = test_response) %>%
      dplyr::summarise_all(.funs = ~ (. - y_test)) %>%
      dplyr::select(-y_test) %>%
      dplyr::summarise_all(.funs = ~ mean(.^2))
    return(list(
      predict_final = pred_fin,
      predict_local = pred_test,
      agg_method = res,
      mse = error,
      running_time = time.taken
    ))
  }
}



# ---------------------------------------------------------------------------------#

