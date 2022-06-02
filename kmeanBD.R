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
#pacman::p_load(parallel)
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