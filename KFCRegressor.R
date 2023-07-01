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

devtools::source_url("https://raw.githubusercontent.com/hassothea/KFC-Procedure/master/kmeanBD.R")

# ---------------------------------------------------------------------------------#

# fitLocalModels <- function(kmeans_BD,
#                           train_response,
#                           model = "lm",
#                           formula = NULL){
#  start_time <- Sys.time()
#  X_train <- kmeans_BD$train_data$X_train
#  y_train <- train_response[!(kmeans_BD$train_data$id_remain)]
#  X_remain <- kmeans_BD$train_data$X_remain
#  nr_remain <- nrow(X_remain)
#  y_remain <- NULL
#  if(nr_remain != 0){
#    y_remain <- train_response[kmeans_BD$train_data$id_remain]
#  }
#  pacman::p_load(tree)
#  pacman::p_load(randomForest)
#  pacman::p_load(xgboost)
#  model_ <- ifelse(model == "tree", tree::tree, model)
#  K <- nrow(kmeans_BD$centroids)
#  if (is.null(formula)){
#    form <- formula(target ~ .)
#  }
#  else{
#    form <- update(formula, target ~ .)
#  }
#  data_ <- bind_cols(X_train, "target":= y_train)
#  fit_lookup <- list(lm = "fitted.values",
#                     rf = "predicted")
#  if(is.character(model_)){
#    model_lookup <- list(lm = lm,
#                         rf = randomForest::randomForest,
#                         xgboost = xgboost)
#    if(model_ == "xgboost"){
#      mat <- as.matrix(X_train)
#      mod <- map(.x = 1:K, 
#                 .f = ~ model_lookup[[model_]](data = mat[kmeans_BD$clusters == .x,],
#                                               label = y_train[kmeans_BD$clusters == .x],
#                                               objective = "reg:squarederror",
#                                               verbose = 0,
#                                               nrounds = 500))
#    } else{
#      mod <- map(.x = 1:K, 
#                 .f = ~ model_lookup[[model_]](formula = form, 
#                                               data = data_[kmeans_BD$clusters == .x, ]))
#    }
#  } else{
#    mod <- map(.x = 1:K, 
#               .f = ~ model_(formula = form, 
#                             data = data_[kmeans_BD$clusters == .x,]))
#  }
#  pred0 <- NULL
#  if(nr_remain != 0){
#    pred0 <- vector(mode = "numeric",
#                    length = length(y_remain))
#    clus <- findClosestCentroid(x. = X_remain,
#                                centroids. = kmeans_BD$centroids,
#                                div = kmeans_BD$parameters$div,
#                                deg = kmeans_BD$parameters$deg)
#    if(model_ == "xgboost"){
#      mat_ <- as.matrix(X_remain)
#      for(i_ in 1:K){
#        pred0[clus == i_] <- predict(mod[[i_]],
#                                     mat_[clus == i_,])
#      }
#    } else{
#      for(i_ in 1:K){
#        pred0[clus == i_] <- predict(mod[[i_]],
#                                     as.data.frame(X_remain[clus == i_,]))
#      }
#    }
#  }
#  time_taken <- Sys.time() - start_time
#  return(list(
#    local_models = mod,
#    kmeans_BD = kmeans_BD,
#    data_remain = list(fit = pred0,
#                       response = y_remain),
#    model_type = model_,
#    running_time = time_taken))
#}

fitLocalModels <- function(kmeans_BD,
                           train_response,
                           model = "lm",
                           formula = NULL){
  start_time <- Sys.time()
  X_train <- kmeans_BD$train_data$X_train
  y_train <- train_response[!(kmeans_BD$train_data$id_remain)]
  X_remain <- kmeans_BD$train_data$X_remain
  nr_remain <- nrow(X_remain)
  y_remain <- NULL
  if(nr_remain != 0){
    y_remain <- train_response[kmeans_BD$train_data$id_remain]
  }
  pacman::p_load(tree)
  pacman::p_load(randomForest)
  pacman::p_load(xgboost)
  model_ <- ifelse(model == "tree", tree::tree, model)
  K <- nrow(kmeans_BD$centroids)
  form <- formula(target ~ .)
  if (!is.null(formula)){
    form <- update(form, formula)
  }
  data_ <- bind_cols(X_train, "target":= y_train)
  fit_lookup <- list(lm = "fitted.values",
                     rf = "predicted")
  if(is.character(model_)){
    model_lookup <- list(lm = lm,
                         rf = randomForest::randomForest,
                         xgboost = xgboost)
    if(model_ == "xgboost"){
      mat <- as.matrix(X_train)
      mod <- map(.x = 1:K, 
                 .f = ~ model_lookup[[model_]](data = mat[kmeans_BD$clusters == .x,],
                                               label = y_train[kmeans_BD$clusters == .x],
                                               objective = "reg:squarederror",
                                               verbose = 0,
                                               nrounds = 300))
    } else{
      mod <- map(.x = 1:K, 
                 .f = ~ model_lookup[[model_]](formula = form, 
                                               data = data_[kmeans_BD$clusters == .x, ]))
    }
    
  } else{
    mod <- map(.x = 1:K, 
               .f = ~ model_(formula = form, 
                             data = data_[kmeans_BD$clusters == .x,]))
  }
  pred0 <- NULL
  if(nr_remain != 0){
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
    model_type = model_,
    running_time = time_taken))
}


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
  if(localModels$model_type == "xgboost"){
    mat <-  as.matrix(newData_)
    for(i_ in 1:K){
      pred0[clus == i_] <- predict(localModels$local_models[[i_]],
                                   mat[clus == i_,])
    }
  } else{
    for(i_ in 1:K){
      pred0[clus == i_] <- predict(localModels$local_models[[i_]],
                                   as.data.frame(newData_[clus == i_,]))
    }
  }
  pred0 <- as_tibble(pred0)
  names(pred0) <- ifelse(kmean_BD$parameters$div == "polynomial",
                         paste0("polynomial", kmean_BD$parameters$deg),
                         kmean_BD$parameters$div)
  return(pred0)
}

# ---------------------------------------------------------------------------------#

# localPredict <- function(localModels,
#                          newData){
#   kmean_BD <- localModels$kmeans_BD
#   K <- nrow(kmean_BD$centroids)
#   newData_ <- newData
#   if(!(is.null(kmean_BD$parameters$center_))){
#     newData_ <- scale(newData,
#                       center = kmean_BD$parameters$center_,
#                       scale = kmean_BD$parameters$scale_)
#     id0 <- (newData_ <= 0)
#     if(sum(id0) > 0){
#       min_ <- min(newData_[id0])
#       newData_[id0] <- runif(sum(id0), min(1e-3, min_/10), min_)
#     }
#   }
#   clus <- findClosestCentroid(x. = newData_,
#                               centroids. = kmean_BD$centroids,
#                               div = kmean_BD$parameters$div,
#                               deg = kmean_BD$parameters$deg)
#   pred0 <- vector(mode = "numeric", length = nrow(newData_))
#   for(i_ in 1:K){
#     pred0[clus == i_] <- predict(localModels$local_models[[i_]],
#                                  as.data.frame(newData_[clus == i_,]))
#   }
#   pred0 <- as_tibble(pred0)
#   names(pred0) <- ifelse(kmean_BD$parameters$div == "polynomial",
#                          paste0("polynomial", kmean_BD$parameters$deg),
#                          kmean_BD$parameters$div)
#   return(pred0)
# }

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
  if(localModels$model_type == "xgboost"){
    mat <-  as.matrix(newData_)
    for(i_ in 1:K){
      pred0[clus == i_] <- predict(localModels$local_models[[i_]],
                                   mat[clus == i_,])
    }
  } else{
    for(i_ in 1:K){
      pred0[clus == i_] <- predict(localModels$local_models[[i_]],
                                   as.data.frame(newData_[clus == i_,]))
    }
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
source_url("https://raw.githubusercontent.com/hassothea/AggregationMethods/main/GradientCOBRARegressor.R")
### MixCobra
source_url("https://raw.githubusercontent.com/hassothea/AggregationMethods/main/MixCobraRegressorressor.R")

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
                  setGridParamMix = setGridParameter_Mix(),
                  silent = FALSE){
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
  list_method_agg <- list(mixcobra = function(pred){MixCobraRegressor(train_input = train_input[!id_shuffle,],
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
                                                                setGridParam = setGridParamMix,
                                                                silent = silent)},
                          cobra = function(pred){GradientCOBRARegressor(train_design = pred,
                                                              train_response = y_remain,
                                                              test_design = pred_test,
                                                              test_response = test_response,
                                                              scale_input = K_step$scale_input,
                                                              scale_machine = C_step$scale_features,
                                                              build_machine = FALSE,
                                                              machines = NULL,
                                                              n_cv = C_step$n_cv,
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
