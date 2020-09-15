

# ================
# From: 09/09/2020
# ================

# All Bregman divergences
# =======================
# euclidDist = function(x, y){
#   return((x-y)%*% (x-y))
# }
# 
# gklDist = function(x, y){
#   return(sum(x*log(abs(x/y))-(x-y)))
# }
# 
# logisticDist = function(x, y){
#   x0 = x/sum(abs(x))
#   y0 = y/sum(abs(y))
#   x1 = 1 - x0
#   y1 = 1 - y0
#   return(sum(x0*log(abs(x0/y0))+x1*log(abs(x1/y1))))
# }
# 
# itakuraDist = function(x, y){
#   ratio = x/y
#   return(sum(ratio-log(abs(ratio)) - 1))
# }
# 
# expDist = function(x, y){
#   y0 = exp(y)
#   return(sum(exp(x)-y0-(x-y)*y0))
# }
# 
# polyDist = function(x, y, deg){
#   if(deg %% 2 == 0){
#     res_dis = sum(x^deg - y^deg - deg * (x-y) * y^(deg - 1))
#   }
#   else{
#     res_dis = sum(x^deg - y^deg - deg * (x-y) * sign(y) * y^(deg - 1))
#   }
#   return(res_dis)
# }


require(Rcpp)
sourceCpp(".../BregmanDivergences.cpp")  # set directory for "BregmanDivergences.cpp" 

# ======================================================================================================================

# * bregmanDiverg function
# ========================
# Input: 
# ------
# - x      : matrix or data frame of the first set of data points with rows of individuals and columns of variables.
# - y      : matrix of data frame of the second set of data points with rows of individuals and columns of variables.
# - method : Bregman divergence to be used. It must be an element of the set {"euclidean", "gkl", "logistic", "itakura", "exponential", "polynomial"}.
# - deg    : degree of the polynomial divergence. It is NULL by default and it has to be non-null if the method = "polynomial".
# Output:
# -------
# - d      : matrix (d_ij) whose d_ij is the Bregman divergence between point (row) i of matrix x and point j (row) of matrix y.

bregmanDiverg = function(x,
                         y,
                         method = c("euclidean",
                                    "gkl",
                                    "logistic",
                                    "itakura",
                                    "exponential",
                                    "polynomial"),
                         deg = NULL) {
  method = match.arg(method)
  if (method == "polynomial")
    if (is.null(deg))
      stop("degree of polynomial is missing")
  if ((ncol(x) == 1) | (ncol(y) == 1)) {
    x = t(x)
    y = t(y)
  }
  d = switch(method,
             euclidean = apply(
               y,
               1,
               FUN = function(a) {
                 return(apply(
                   x,
                   1,
                   FUN = function(b) {
                     return(euclidDist(a, b))
                   }
                 ))
               }
             ),
             gkl = apply(
               y,
               1,
               FUN = function(a) {
                 return(apply(
                   x,
                   1,
                   FUN = function(b) {
                     return(gklDist(a, b))
                   }
                 ))
               }
             ),
             logistic = apply(
               y,
               1,
               FUN = function(a) {
                 return(apply(
                   x,
                   1,
                   FUN = function(b) {
                     return(logisticDist(a, b))
                   }
                 ))
               }
             ),
             itakura = apply(
               y,
               1,
               FUN = function(a) {
                 return(apply(
                   x,
                   1,
                   FUN = function(b) {
                     return(itakuraDist(a, b))
                   }
                 ))
               }
             ),
             exponential = apply(
               y,
               1,
               FUN = function(a) {
                 return(apply(
                   x,
                   1,
                   FUN = function(b) {
                     return(expDist(a, b))
                   }
                 ))
               }
             ),
             polynomial = apply(
               y,
               1,
               FUN = function(a) {
                 return(apply(
                   x,
                   1,
                   FUN = function(b) {
                     return(polyDist(a, b, deg))
                   }
                 ))
               }
             )
  )
  return(d)
}

# ======================================================================================================================

# * findClosestCentroid function
# ==============================
# Input:
# ------
# - x         : matrix or data frame containing data points to be assigned to the one of the given centroids.
# - centroids : matrix or data frame containing the centroids.
# - method    : Bregman divergence to be used. It must be an element of the set {"euclidean", "gkl", "logistic", "itakura", "exponential", "polynomial"}.
# - deg       : degree of the polynomial divergence. It is NULL by default and it has to be non-null if the method = "polynomial".
# Output:
# -------
# - clusters  : a vector of containing new clusters of the data points of x.

findClosestCentroid = function(x, centroids, method, deg = NULL) {
  distances = bregmanDiverg(x, centroids, method, deg = deg)
  clusters = apply(distances, 1, which.min)
  return(clusters)
}

# ======================================================================================================================

# * new Centroids computation : newCentroids
# ==========================================
# Inputs:
# -------
# - x           : matrix or data frame of row points to used recompute the centroids according to the "clusters".
# - clusters    : vector containing clusters.
# Output:
# -------
# - centroids : matrix containing the row points of the new centroids.

newCentroids = function(x, clusters) {
  centroids = apply(
    x,
    2,
    FUN = function(c)
      return(tapply(c, clusters, mean))
  )
  return(centroids)
}

# ======================================================================================================================

# * K-step: K-means algorithm with Bregman divergences
#=====================================================
# Input:
# ------
# - x            : matrix of row points to be clustered.
# - K            : the number of clusters.
# - n.strat      : the number of implementation of the algorithm to avoid being stuck in the local clustering structure.
# - maxIter      : the number of iteration in each implementation of the method (to avoid being trapped in the loop when the algorithm does not converge).
# - method       : Bregman divergence to be used.
# - epsilon      : stopping criterion.
# - figure       : make a 2 or 3-dimensional plot of the data, it is a boolean type controls whether or not to plot the result of the algorithm.
# Output: a list of the following outputs,
# -------
# - centroids    : matrix of centers of the K groups.
# - clusters     : vector of length equals to the number of points containing the cluster of each point.
# - method       : the Bregman divergence used.
# - Ranning.time : the

kmeansBregman = function(x,
                         K,
                         n.start = 10,
                         maxIter = 500,
                         deg = NULL,
                         method = "euclidean",
                         epsilon = 1e-10,
                         figure = F) {
  x = as.matrix(x)
  start.time <- Sys.time()
  # Distortion function
  distortion = function(l) {
    cent = apply(x, 2,
                  FUN = function(c)
                    return(tapply(c, l, mean)))
    if(ncol(cent) != ncol(x))
      cent = t(cent)
    mu = t(as.matrix(colMeans(x)))
    var.between = sum((bregmanDiverg(cent, mu, method = method, deg = deg)) ^ 2)
    return(var.between)
  }
  
  # Kmeans algorithm
  kmeansWithBD = function(x, k, n, maxiter, eps) {
    # initialization
    init = sample(n, k)
    centroids.old = x[init,]
    i = 0
    while (i < maxiter) {
      # Assignment step
      clusters = findClosestCentroid(x, centroids.old, method, deg = deg)
      # Recompute centroids
      centroids.new = newCentroids(x, clusters)
      if ((sum(is.na(centroids.new)) > 0) |
          (nrow(centroids.new) < k)) {
        init = sample(n, k)
        centroids.old = x[init,]
        warning("NA produced -> reinitialize centroids...!")
      }
      else{
        if (sum(abs(centroids.new - centroids.old)) > eps) {
          centroids.old = centroids.new
        }
        else{
          break
        }
      }
      i = i + 1
    }
    return(clusters)
  }
  size = dim(x)
  results = sapply(
    1:n.start,
    FUN = function(i)
      return(kmeansWithBD(x, K, size[1], maxIter, epsilon))
  )
  vars = apply(results, 2, FUN = distortion)
  Opt.ind = which.max(vars)
  cluster = clusters = results[, Opt.ind]
  k = 1
  ID = unique(cluster)
  for (i in ID) {
    clusters[which(cluster == i)] = k
    k =  k + 1
  }
  centroids = newCentroids(x, clusters)
  if (figure) {
    par(mfrow = c(1, 1))
    if (ncol(x) == 2) {
      plot(x,
           col = 1 + clusters,
           main = "K-means Algorithm",
           pch = 19)
      points(
        centroids,
        cex = 2,
        pch = 10,
        col = 1 + unique(clusters)
      )
    }
    if (ncol(x) == 3) {
      s3d = scatterplot3d::scatterplot3d(x,
                                         color = 1 + clusters,
                                         main = "K-means Algorithm",
                                         pch = 19)
      s3d$points3d(
        centroids,
        cex = 2,
        pch = 10,
        col = 1 + unique(clusters)
      )
    }
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  return(
    list(
      centroids = centroids,
      clusters = clusters,
      method = method,
      running.time = time.taken
    )
  )
}


# ======================================================================================================================

# * F-step: Fitting linear regression model on each cluster
# =========================================================
# Inputs:
# -------
# - data          : matrix or data frame containing the dataset.
# - target        : response varaible.
# - methods       : list of K-means algorithms givem by the function "kmeansBregman".
# - formula       : forumula for linear regression model fitting on each cluster. It is a class of formula and NULL by default.
# Output:
# -------
# - MSE           : mean square error of training data.
# - Coefficients  : list of coefficients of all constructed local models.
# - Fitted.values : predictions given by each candidate model corresponding to each Bregman divergence.
# - Running.Time  : the running time of the algortihm.

fitClusterRegressionModel = function(data, target, methods, formula = NULL) {
  start.time <- Sys.time()
  K = nrow(methods[[1]]$centroids)
  m = length(methods)
  N = length(target)
  if (is.null(formula))
    form = formula(temp.y ~ .)
  else
    form = update(formula, temp.y ~ .)
  mse = matrix(0, m, K)
  coef = list()
  fit = matrix(0, N, m)
  for (i in 1:m) {
    cluster = as.numeric(methods[[i]]$clusters)
    #    print(table(cluster))
    coef[[i]] = list()
    for (j in 1:K) {
      id = which(cluster == j)
      temp.y = target[id]
      model = lm(form, data = as.data.frame(data[id, ]))
      mse[i, j] = sum((model$residuals) ^ 2)
      fit[id, i] = model$fitted.values
      coef[[i]][[j]] = model$coefficients
    }
  }
  mse = rowSums(mse) / N
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  return(list(
    MSE = mse,
    Coefficients = coef,
    Fitted.values = fit,
    Running.Time = time.taken
  ))
}

# ======================================================================================================================

# * generateMachines: generate machines to be combined in consensual aggregation method
# =====================================================================================
# Input:
# ------
# - train.input    : matrix or data frame of training input data.
# - train.response : vector of training response.
# - test.input     : matrix or data frame of testing data.
# - machines       : choices of machines which must be a subset of {"lm","lars", "ridge", "knn", "tree", "rf"}.
# - splits         : the proportion taking value between (0,1) to split the data. By default it is NULL and all the five machiens are built.
# - k              : the paramter "k" of kNN method. 
# - ntree          : the paramter "ntree" of rf method.
# - mtry           : number of variable to split at each step of random forest.
# Output:
# -------
# - train.machine  : matrix of predictions of training data given by all the machines.
# - test.machine   : matrix of predictions of testing data given by all the machines.
# - id.val         : the random index of the data after shuffled.
# - machines       : the vector of machines used.

generateMachines = function(train.input, train.response, test.input, machines = NULL, splits = NULL, k = 10, ntree = 300, mtry = NULL){
  require(lars)
  require(tree)
  require(ridge)
  require(randomForest)
  require(FNN)
  d = dim(train.input)
  COBRA.mlr <- function(x){
    dff = cbind(train.response.1,train.input.1)
    colnames(dff) = c("y", paste0("x",1:ncol(train.input.1)))
    colnames(x) = paste0("x",1:ncol(train.input.1))
    a <- lm(y~., data = as.data.frame(dff))
    res <- predict(a, as.data.frame(x))
    return(res)
  }
  COBRA.lars <- function(x) {
    a <- lars(as.matrix(train.input.1), train.response.1, intercept = TRUE, use.Gram = G)
    res <- predict.lars(a, x, s = length(a$df))$fit
    return(res)
  }
  COBRA.tree <- function(x) {
    a <- tree(as.formula(paste("train.response.1~", 
                               paste(colnames(as.data.frame(x)), sep = "", 
                                     collapse = "+"), collapse = "", sep = "")), 
              data = as.data.frame(train.input.1))
    res <- as.vector(predict(a, as.data.frame(x)))
    return(res)
  }
  COBRA.ridge <- function(x) {
    a <- linearRidge(as.formula(paste("train.response.1~", 
                                      paste(colnames(as.data.frame(x)), sep = "", 
                                            collapse = "+"), collapse = "", sep = "")), 
                     data = as.data.frame(train.input.1))
    res <- as.vector(predict(a, as.data.frame(x)))
    return(res)
  }
  COBRA.randomForest <- function(x) {
    if(is.null(mtry))
      a <- randomForest(x = train.input.1, y = train.response.1, ntree = ntree)
    else
      a <- randomForest(x = train.input.1, y = train.response.1, ntree = ntree, mtry = mtry)
    res <- as.vector(predict(a, x))
    return(res)
  }
  COBRA.knn <- function(x) {
    a <- knn.reg(train = train.input.1, test = x, y = train.response.1, k = k)
    res = a$pred
    return(res)
  }
  all.machines = list(lm = COBRA.mlr ,lars = COBRA.lars, ridge = COBRA.ridge, knn = COBRA.knn, tree = COBRA.tree, rf = COBRA.randomForest)
  if(is.null(machines))
    mach = c("lm", "lars", "ridge", "knn", "tree", "rf")
  else
    mach = machines
  M = length(mach)
  n1 = d[1]
  n2 = nrow(test.input)
  id.rand = sample(n1)
  G = ifelse(d[1]<d[2], FALSE, TRUE)
  if(is.null(splits)){
    I = floor(n1/2)
    train.input.1 = train.input[id.rand[1:I],]
    train.input.2 = train.input[id.rand[(I+1):n1],]
    train.response.1 = train.response[id.rand[1:I]]
  }
  else{
    I = floor(n1*splits)
    train.input.1 = train.input[id.rand[1:I],]
    train.input.2 = train.input[id.rand[(I+1):n1],]
    train.response.1 = train.response[id.rand[1:I]]
  }
  res.tr = matrix(nrow = n1-I, ncol = M, data = 0)
  res.te = matrix(nrow = n2, ncol = M, data = 0)
  cat("Build machines\n==============\n")
  cat("\n* Machines...")
  if(M < 11){
    for(i in 1:M){
      if(i == 1)
        cat(" | --->",i)
      else
        cat(" --->",i)
    }
    cat("\n* Processing: | ")
  }
  else{
    cat("\n")
  }
  for(m in 1:M){
    df = as.matrix(rbind(train.input.2, test.input))
    pred = all.machines[[mach[m]]](df)
    res.tr[1:(n1-I),m] = pred[1:(n1-I)]
    res.te[1:n2,m] = pred[(n1-I+1):(n1+n2-I)]
    if(M < 11){
      if(m < M)
        cat("~~~>",m,"")
      else
        cat("~~~>",m,"done!\n\n")
    }
    else{
      if(m < M)
        cat("* Processing:",m,"out of",M,"machines have been built.\r")
      else
        cat("* Processing:",m,"out of",M,"machines have been built.\n\n")
    }
  }
  colnames(res.te) = mach
  colnames(res.tr) = mach
  return(list(train.machine = res.tr, test.machine = res.te, id.val = id.rand[(I+1):n1],
              machines = mach))
}

# ======================================================================================================================

# * setGD: setting of Gradient descent algorithm to estimate the window parameter of the consensual aggregation method
# ====================================================================================================================
# Input:
# ------
# - initial.ep : initial value of the window parameter. It is NULL by default.
# - n.tries    : the number of tries for the initial value of the parameter if the initial.ep = NULL. It is equal to 10 by default.
# - rate       : either the numerical learning rate or string corresponding to the speed of the learning rate. It must be an element of 
#                the set {"auto", "logarithm", "sqrtroot", "linear", "polynomial", "exponential"}. By default, it is NULL.
# - threshold  : threshold of relative change of the parameter. By defualt it is equal to 1e-10.
# - max.iter   : the maximum number of iteration of the method to stop if it does not converg. by default, it is equal to 500 iterations.
# - print.step : boolean type defining whether or not to print the values of the parameter, gradient and relative change in each iteration.
#                By default, it is TURE.
# - figure     : boolean type defining whether or not to plot the graph of parameter vs gradient. It is FALSE by default.
# - n.cv       : the number of folds in cross valisation method. It equals to 10 by default.
# - coef.auto  : constant learning rate of the algorithm.
# - coef.log   : coefficient of the logarithm rate (a) in where rate ~ a*log(i).
# - coef.sqrt  : coefficient of the square root rate (a) in where rate ~ a*sqrt(i).
# - coef.lm    : coefficient of the linear rate (a) in where rate ~ a*i.
# - deg.poly   : degree of the polynomial rate (a) in where rate ~ i^a.
# - base.exp   : base of the exponential rate (a) in where rate ~ a^i.
# Output: list of all the inputs.
# -------

setGD = function(initial.ep = NULL, n.tries = 10, rate = NULL, threshold = 1e-10, max.iter = 500, print.step = TRUE, figure = FALSE, n.cv = 10,
                 coef.auto = 0.1,
                 coef.log = 1,
                 coef.sqrt = 1,
                 coef.lm = 1,
                 deg.poly = 2,
                 base.exp = 1.5) {
  return(
    list(
      initial.ep = initial.ep,
      n.tries = n.tries,
      rate = rate,
      threshold = threshold,
      max.iter = max.iter,
      n.cv = n.cv,
      print.step = print.step,
      figure = figure,
      coef.auto = coef.auto,
      coef.log = coef.log,
      coef.sqrt = coef.sqrt,
      coef.lm = coef.lm,
      deg.poly = deg.poly,
      base.exp = base.exp
    )
  )
}

# * setGrid: setting of Grid search algorithm
# ===========================================
# Input:
# ------
# - min.ep    : minimum of the grid of paramters. It equals to 1e-100 by default.
# - max.ep    : maximum of the grid of paramters. It equals to NULL by default.
# - n.ep      : the number of discretization of the grid. It equals to 300 by default.
# - paramters : the grid in case of non-uniform grid.
# Output: list of all the inputs.
# -------

setGrid = function(min.ep = 1e-100, max.ep = NULL, n.ep = 300, parameters = NULL){
  return(list(min.ep = min.ep,
              max.ep = max.ep,
              n.ep = n.ep,
              parameters = parameters))
}

# * setKernels: setting of kernel function
# ========================================
# Input:
# ------
# - kernels : vector of string defining the kernels used. It must a subset of {"naive", "uniform", "epanechnikov", "biweight", "triweight", 
#             "triangular", "compact_gaussian", "gaussian", "exp4", "cauchy"}. 
# - sigma   : parameter sigma of Gaussian-type and exp4 kernels.
# - rho     : parameter defining the support of the compactly supported Gaussian kernel.
# Output: list of all the inputs.
# -------

setKernels = function(kernels = NULL, sigma = 1, rho = 1){
  return(list(kernels = kernels, sigma = sigma, rho = rho))
}

# ======================================================================================================================

# * optimizationSmoothParameter: optimization method to estimate the window parameter
# ===================================================================================
# Input:
# ------
# - train.machines  : matrix of predictions of training data.
# - train.responses : vector of response variable.
# - test.machines   : matrix of predictions of testing data.
# - opt.method      : optimization method. It is either "GD" or "grid" which corresponds to gradient descent and grid search resp.
# - set.GD          : "setGD" function to define the parameters of the gradient descent algorithm.
# - set.Grid        : "setGrid" function to define the parameters of the grid search algorithm.
# - set.Kernels     : "setKernels" function to define the parameters of the kernel method.
# Output: list of the following values,
# -------
# - epsilon         : the observed value of paramter.
# - n.machines      : the number of retaining machines in classical COBRA. It is NA for other cases.
# - cv.error        : the optimal cross-validation quadratic error.

optimizeSmoothParameter = function(train.machines,
                                   train.responses,
                                   test.machines,
                                   opt.method = c("grid", "GD"),
                                   set.GD = setGD(initial.ep = NULL, 
                                                  n.tries = 10,
                                                  rate = "auto", 
                                                  threshold = 1e-10,
                                                  max.iter = 500,
                                                  n.cv = 10,
                                                  print.step = FALSE,
                                                  figure = FALSE,
                                                  coef.auto = .001,
                                                  coef.log = 1,
                                                  coef.sqrt = 1,
                                                  coef.lm = 1,
                                                  deg.poly = 2,
                                                  base.exp = 2),
                                   set.Grid = setGrid(min.ep = 1e-300,
                                                      max.ep = NULL,
                                                      n.ep = 500,
                                                      parameters = NULL),
                                   set.Kernels = setKernels(kernels = "naive",
                                                            sigma = 1, rho = rho)
) {
  # parameters Grid search
  min.ep = set.Grid$min.ep
  max.ep = set.Grid$max.ep
  n.ep = set.Grid$n.ep
  parameters = set.Grid$parameters
  # Parameters GD
  initial.GD = set.GD$initial.ep
  n.tries = set.GD$n.tries
  rate.GD = set.GD$rate
  err.GD = set.GD$threshold
  max.iter = set.GD$max.iter
  n.cv = set.GD$n.cv
  print.step = set.GD$print.step
  figure = set.GD$figure
  coef.auto = set.GD$coef.auto
  coef.log = set.GD$coef.log
  coef.sqrt = set.GD$coef.sqrt
  coef.lm = set.GD$coef.lm
  deg.poly = set.GD$deg.poly
  base.exp = set.GD$base.exp
  # Kernels
  kernel = set.Kernels$kernels
  sigma = set.Kernels$sigma
  rho = set.Kernels$rho
  
  
  d = dim(train.machines)
  n.folds = floor(d[1] / n.cv)
  CV = sample(d[1])
  opt.method = match.arg(opt.method)
  kernel = match.arg(kernel, c("naive", "uniform", "epanechnikov", "biweight", "triweight", 
                               "triangular", "compact_gaussian", "gaussian", "exp4", "cauchy"))
  Dist = list()
  for (k in 1:n.cv) {
    q = ifelse(k < n.cv, (n.folds * k), d[1])
    id.test = CV[((k - 1) * n.folds + 1):q]
    id.remain = setdiff(CV, id.test)
    # Dist is of dim = id.remain * id.test
    Dist[[k]] = sapply(id.test, function(i) {
      temp = sweep(train.machines[id.remain, ], 2, train.machines[i, ])
      if (kernel %in% c("naive", "uniform"))
        return(abs(temp))
      else{
        if (kernel == "triangular")
          return(rowSums(abs(temp)))
        else{
          if (kernel == "exp4")
            return(rowSums((temp) ^ 4))
          else
            return(rowSums((temp) ^ 2))
        }
      }
    })
    colnames(Dist[[k]]) = id.test
    if (kernel %in% c("naive", "uniform"))
      rownames(Dist[[k]]) = rep(id.remain, d[2])
    else
      rownames(Dist[[k]]) = id.remain
  }
  
  # The kernel functions
  naive.weight = function(L, eps) {
    id0 = apply(L, 2, function(x)
      return(rowSums(matrix(x, ncol = d[2], byrow = FALSE) <= eps)))
    return(id0)
  }
  epan.weight = function(L, eps) {
    val =  L / eps ^ 2
    id0 = val <= 1
    val[!id0] = 1
    res = (1 - val)
    res[is.na(res)] = 0
    return(res)
  }
  bi.weight = function(L, eps) {
    val =  L / eps ^ 2
    id0 = val <= 1
    val[!id0] = 1
    res = ((1 - val) ^ 2)
    res[is.na(res)] = 0
    return(res)
  }
  tri.weight = function(L, eps) {
    val =  L / eps ^ 2
    id0 = val <= 1
    val[!id0] = 1
    res = ((1 - val) ^ 3)
    res[is.na(res)] = 0
    return(res)
  }
  triang.weight = function(L, eps) {
    val =  L / eps
    id0 = val <= 1
    val[!id0] = 1
    res = (1 - val)
    res[is.na(res)] = 0
    return(res)
  }
  cGaus.weight = function(L, eps) {
    val =  L / (2 * (sigma * eps) ^ 2)
    id0 = val <= rho
    res = exp(-val)
    res[!id0 | is.na(res)] = 0
    return(res)
  }
  gaus.weight = function(L, eps) {
    val =  L / (2 * (sigma * eps) ^ 2)
    res = exp(-val)
    res[is.na(res)] = 0
    return(res)
  }
  exp4.weight = function(L, eps) {
    val =  L / (2 * (sigma * eps) ^ 4)
    res = exp(-val)
    res[is.na(res)] = 0
    return(res)
  }
  cauchy.weight = function(L, eps) {
    val =  1 / (1 + L / eps ^ 2)
    val[is.na(val)] = 0
    return(val)
  }
  
  # define the right weight for the chosen kernel
  weight = switch(
    kernel,
    naive = naive.weight,
    uniform = naive.weight,
    epanechnikov = epan.weight,
    biweight = bi.weight,
    triweight = tri.weight,
    triangular = triang.weight,
    compact_gaussian = cGaus.weight,
    gaussian = gaus.weight,
    exp4 = exp4.weight,
    cauchy = cauchy.weight
  )
  
  # Cross-validation error evaluated at a given epsilon
  costFunction = function(D, eps) {
    res = sapply(D,
                 FUN = function(L) {
                   wgt = weight(L, eps)
                   id.test = as.numeric(colnames(L))
                   if (kernel %in% c("naive", "uniform")) {
                     id.remain = as.numeric(names(L[1:(d[1] - length(id.test)), 1]))
                     temp = sapply(
                       1:length(id.test),
                       FUN = function(i) {
                         pre = sapply(
                           1:d[2],
                           FUN = function(id)
                             return(wgt[, i] == id)
                         )
                         temp0 = train.responses[id.remain] %*% pre / colSums(pre)
                         temp0[is.na(temp0)] = 0
                         return(temp0 - train.responses[id.test[i]])
                       }
                     )
                     r = rowSums(temp ^ 2)
                     return(r)
                   }
                   else{
                     id.remain = as.numeric(rownames(L))
                     temp = train.responses[id.remain] %*% wgt / colSums(wgt)
                     temp[is.na(temp)] = 0
                     r = sum((temp - train.responses[id.test])^2)
                     return(r)
                   }
                 }
    )
    if (kernel %in% c("naive", "uniform"))
      return(rowMeans(res))
    else
      return(mean(res))
  }
  
  # Smoothing parameters
  if (is.null(parameters))
    epsilon = seq(min.ep, max.ep, length.out = n.ep)
  else
    epsilon = parameters
  
  # Optimization step:
  # ==================
  if ((opt.method == "grid") | (kernel %in% c("naive", "uniform"))){
    if(print.step)
      cat("\n   ~> Grid search with",kernel,"kernel...")
    cv.risks = sapply(epsilon,
                      FUN = function(ep)
                        return(costFunction(Dist, ep))
    )
    if (kernel %in% c("naive", "uniform")){
      if (opt.method == "GD") {
        warning("Naive kernel is not suitable to be optimized using GD! Grid search method is used!")
      }
      id.opt = which(cv.risks == min(cv.risks), arr.ind = TRUE)[1,]
      op.eps = epsilon[id.opt[2]]
      n.M = id.opt[1]
      op.risk = cv.risks[id.opt[1], id.opt[2]]
    }
    else{
      id.opt = which.min(cv.risks)
      op.eps = epsilon[id.opt]
      n.M = NULL
      op.risk = cv.risks[id.opt]
    }
    if(print.step)
      cat("done!")
  }
  else{
    collect.ep = c()
    gradients = c()
    if (is.null(initial.GD)){
      epsilons = seq(min.ep, max.ep, length.out = n.tries)
      gr = sapply(epsilons,
                  FUN = function(x)
                    costFunction(Dist, x)
      )
      ep_ = ep0 = epsilons[which.min(gr)]
      grad = pracma::grad(
        f = function(x)
          costFunction(Dist, x),x0 = ep0, heps = .Machine$double.eps ^ (1 / 3))
    }
    else
      ep_ = ep0 = initial.GD
    if(print.step){
      cat("\n   ~> Gradient descent with",kernel,"kernel...")
      cat("\n  Initialization:   \t~ epsilon:", format(round(ep0, 5) , nsmall = 5, digits = 10), "\t~ gradient:", format(round(grad, 10) , digits = 10), 
          "  \t~        threshold:", err.GD, "\n")
    }
    if (is.numeric(rate.GD))
      lambda0 = rate.GD
    else{
      r0 = coef.auto / abs(grad)
      rate.GD = match.arg(rate.GD, c("auto", "logarithm", "sqrtroot", "linear", "polynomial", "exponential"))
      lambda0 = switch(rate.GD,
                       auto = r0,
                       logarithm = function(i)
                         coef.log * log(2 + i) * r0,
                       sqrtroot = function(i)
                         coef.sqrt * sqrt(i) * r0,
                       linear = function(i)
                         coef.lm * (i) * r0,
                       polynomial = function(i)
                         i ^ deg.poly * r0,
                       exponential = function(i)
                         base.exp ^ i * r0
      )
    }
    i = 0
    if (is.numeric(rate.GD) | rate.GD == "auto") {
      while (i < max.iter) {
        ep = ep0 - lambda0 * grad
        if (ep < 0){
          ep = ep0 - runif(1,ep0/5,ep0/2)
        }
        if(is.na(ep))
          ep = ep0 + runif(1, - ep0/3, ep0/3)
        if(i > 5){
          if(sign(grad)*sign(grad0) < 0){
            lambda0 = lambda0 / 1.5
          }
        }
        relative = abs((ep - ep0) / ep0)
        test.threshold = max(relative, abs(grad))
        if (test.threshold > err.GD){
          ep0 = ep
          grad0 = grad
        }
        else{
          break
        }
        grad = pracma::grad(
          f = function(x)
            costFunction(Dist, x), x0 = ep0, heps = .Machine$double.eps ^ (1 / 3)
        )
        i = i + 1
        if(print.step){
          cat("       Iteration:", i, "\t~ epsilon:", format(round(ep, 5) , nsmall = 5), "\t~ gradient:", grad, 
              "     \t~ max(change,grad):", test.threshold, "\r")
        }
        # cat("\n*", i, "\t\t|\t", grad, "\t|\t", ep)
        collect.ep = c(collect.ep, ep)
        gradients = c(gradients, grad)
      }
    }
    else{
      while (i < max.iter) {
        ep = ep0 - lambda0(i) * grad
        if (ep < 0){
          ep = ep0 - runif(1,ep0/5,ep0/2)
        }
        if(is.na(ep))
          ep = ep0 + runif(1, - ep0/3, ep0/3)
        if(i > 5){
          if(sign(grad)*sign(grad0) < 0)
            r0 = r0 / 1.5
        }
        relative = abs((ep - ep0) / ep0)
        test.threshold = max(relative, abs(grad))
        if (test.threshold > err.GD){
          ep0 = ep
          grad0 = grad
        }
        else{
          break
        }
        grad = pracma::grad(
          f = function(x)
            costFunction(Dist, x), x0 = ep0, heps = .Machine$double.eps ^ (1 / 3)
        )
        i = i + 1
        if(print.step){
          cat("       Iteration:", i, "\t~ epsilon:", format(round(ep, 5) , nsmall = 5), "\t~ gradient:", grad, 
              "     \t~ max(change,grad):",test.threshold, "\r")
        }
        collect.ep = c(collect.ep, ep)
        gradients = c(gradients, grad)
      }
    }
    op.eps = ep
    n.M = NULL
    op.risk = costFunction(Dist, op.eps)
    if(print.step){
      cat("                                                                                                                  \r")
      if(grad == 0){
        cat("          Stopped:", i, "\t~ epsilon:", format(round(ep, 5) , nsmall = 5), "\t~ gradient:", grad, 
            "\t\t\t~ max(change,grad):",test.threshold)
      }
      else{
        cat("         Stopped:", i, "\t~ epsilon:", format(round(ep, 5) , nsmall = 5), "\t~ gradient:", grad, 
            "    \t~ max(change,grad):",test.threshold)
      }
    }
  }
  if (figure) {
    if (opt.method == "grid") {
      if(kernel %in% c("naive", "uniform")){
        library(scatterplot3d)
        persp(
          x = epsilon,
          y = 1:d[2],
          z = t(cv.risks),
          main = paste("Grid seearch with", kernel, "kernel"),
          theta = -10,
          xlab = "Epsilon",
          ylab = "No. of machines",
          zlab = "MSE",
          col = "lightblue3",
          zlim = c(0, max(cv.risks) + 1)
        ) -> img3d
        points(
          trans3d(
            x = op.eps,
            y = n.M,
            z =  0,
            pmat = img3d
          ),
          col = "red3",
          pch = "*",
          cex = 1.5
        )
      }
      else{
        plot(
          epsilon,
          cv.risks,
          main = paste("Grid seearch with", kernel, "kernel"),
          xlab = "parameter",
          ylab = "CV risk",
          type = "l",
          col = "red",
          cex = 2
        )
        points(
          x = op.eps,
          y = op.risk,
          pch = 19,
          col = "blue",
          cex = 1.5
        )
      }
    }
    else{
      I = length(collect.ep)
      plot(collect.ep, gradients, main = paste("Gradient descent with", kernel, "kernel"),
           xlab = "Parameter",
           ylab = "Gradient",
           type = "l",
           col = "red",
           cex = 2)
      points(collect.ep[I], gradients[I], pch = 19, col = "blue", cex = 1.5)
    }
  }
  return(list(
    epsilon = op.eps,
    n.machines = n.M,
    cv.error = op.risk
  ))
}




# * setMachines: setting machines in COBRA method
# ===============================================
# Input:
# ------
# - build.machines : boolean type define whether or not to build the machines.
# - machines       : choices of machines which must be a subset of {"lars", "ridge", "knn", "tree", "rf"}.
# - splits         : the proportion taking value between (0,1) to split the data. By default it is NULL and all the five machiens are built.
# - scale.machines : boolean type defining whether or not to scale the machines. It must be a subset of {"normalization","standardization"}.
# - k              : the parameter k of kNN method.
# - ntree          : the parameter ntree of rf method.

setMachines = function(build.machines = TRUE,
                       machines = NULL,
                       splits = NULL,
                       scale.machines = NULL,
                       k = 10, 
                       ntree = 300,
                       mtry = NULL){
  return(list(build.machines = build.machines,
              machines = machines,
              splits = splits,
              scale.machines = scale.machines,
              k = k, 
              ntree = ntree,
              mtry = mtry))
}


# ======================================================================================================================

# *** kernelCOBRA: the kernel-based COBRA method
# ==============================================
# Input:
# ------
# - train.design   : matrix or data frame of training data or machines.
# - train.response : vector of training response.
# - test.design    : matrix or data frame of testing data or machines.
# - opt.methods    : optimation method.
# - scale.input    : boolean type defining whether or not to scale the input. It must be a subset of {"normalization","standardization"}.
# - opt.method     : optimization method. It is either "GD" or "grid" which corresponds to gradient descent and grid search resp.
# - set.GD         : "setGD" function to define the parameters of the gradient descent algorithm.
# - set.Grid       : "setGrid" function to define the parameters of the grid search algorithm.
# - set.Kernels    : "setKernels" function to define the parameters of the kernel method.
# - set.Machines   : "setMachines" function to define the parameters of the machines to be built.
# Output: list of predictions and optimal values of parameters.
# -------

kernelCOBRA = function(train.design,
                       train.response,
                       test.design,
                       opt.methods = c("grid", "GD"),
                       scale.input = NULL,
                       set.Machines = setMachines(build.machines = TRUE,
                                                  machines = NULL,
                                                  splits = NULL,
                                                  scale.machines = NULL,
                                                  k = 10, ntree = 300,
                                                  mtry = NULL),
                       set.GD = setGD(initial.ep = NULL, 
                                      n.tries = 10,
                                      rate = "auto",
                                      threshold = 1e-10,
                                      max.iter = 500,
                                      n.cv = 10,
                                      print.step = FALSE,
                                      figure = FALSE,
                                      coef.auto = .001,
                                      coef.log = 1,
                                      coef.sqrt = 1,
                                      coef.lm = 1,
                                      deg.poly = 2,
                                      base.exp = 2),
                       set.Grid = setGrid(min.ep = 1e-300,
                                          max.ep = NULL,
                                          n.ep = 300,
                                          parameters = NULL),
                       set.Kernels = setKernels(kernels = c("naive", "uniform", "epanechnikov", "biweight", "triweight", 
                                                            "triangular", "compact_gaussian", "gaussian", "exp4", "cauchy"),
                                                sigma = 1, rho = 1),
                       print.info = TRUE, 
                       print.train.results = TRUE) {
  train.design = as.matrix(train.design)
  test.design = as.matrix(test.design)
  
  # parameters Grid search
  min.eps = set.Grid$min.ep
  max.eps = set.Grid$max.ep
  n.eps = set.Grid$n.ep
  param.list = set.Grid$parameters
  
  # Parameters GD
  initial.GD = set.GD$initial.ep
  n.tries = set.GD$n.tries
  rate.GD = set.GD$rate
  err.GD = set.GD$threshold
  max.iter = set.GD$max.iter
  n.cv = set.GD$n.cv
  print.GD.step = set.GD$print.step
  figure = set.GD$figure
  coef.auto = set.GD$coef.auto
  coef.log = set.GD$coef.log
  coef.sqrt = set.GD$coef.sqrt
  coef.lm = set.GD$coef.lm
  deg.poly = set.GD$deg.poly
  base.exp = set.GD$base.exp
  
  # Kernels
  kernels = set.Kernels$kernels
  sigma = set.Kernels$sigma
  rho = set.Kernels$rho
  
  cat("\n*** Kernel-based COBRA Method ***\n=================================\n")
  
  # Sale input
  train.responses = train.response
  train.machine0 = train.design
  test.machine0 = test.design
  if(is.null(scale.input)){
    train.machines = train.design
    test.machines = test.design
  }
  else{
    scale.input = match.arg(scale.input, c("standardization", "normalization"))
    if(scale.input == "standardization"){
      SD = apply(train.design, 2, sd)
      SD.inv = diag(1/SD)
      colnames(SD.inv) = colnames(train.design)
      train.machines = train.design %*% SD.inv
      test.machines = test.design %*% SD.inv
    }
    if(scale.input == "normalization"){
      max.inp = max(train.design)
      min.inp = min(train.design)
      train.machines = (train.design - min.inp)/(max.inp - min.inp)
      test.machines = (test.design - min.inp)/(max.inp - min.inp)
    }
  }
  x.train = train.design
  y.train = train.response
  # Building machines
  if(set.Machines$build.machines){
    build = generateMachines(train.input = train.machines, train.response = train.responses, 
                             test.input = test.machines, machines = set.Machines$machines, 
                             splits = set.Machines$splits,
                             k = set.Machines$k, ntree = set.Machines$ntree, mtry = set.Machines$mtry)
    train.machines = train.machine0 = build$train.machine
    train.responses = train.response[build$id.val]
    test.machines = test.machine0 = build$test.machine
    max.sd = max(train.machines)
    x.train = train.design[build$id.val,]
    y.train = train.response[build$id.val]
    # Scalling 
    if(!is.null(set.Machines$scale.machines)){
      scale_ = match.arg(set.Machines$scale.machines, c("normalization", "standardization"))
      if(scale_ == "standardization"){
        SD = apply(train.machines , 2, sd)
        SD.inv = diag(1/SD)
        colnames(SD.inv) = colnames(train.machines)
        train.machines = train.machines %*% SD.inv
        test.machines = test.machines %*% SD.inv
        max.sd = max(SD)
      }
      if(scale_ == "normalization"){
        max.sd = max(train.machines)
        train.machines = train.machines / max.sd
        test.machines = test.machines / max.sd
      }
    }
  }
  
  # Maximum epsilon
  if(is.null(max.eps))
    max.eps = 1.5*max.sd*ncol(train.machines)
  
  # All kernels are used
  #=====================
  MU = 0
  d = ifelse(is.null(set.Machines$machines), 5, length(set.Machines$machines))
  all.kernels = c(
    "naive",
    "uniform",
    "epanechnikov",
    "biweight",
    "triweight",
    "triangular",
    "compact_gaussian",
    "gaussian",
    "exp4",
    "cauchy"
  )
  ALL = FALSE
  kers = kernels
  ker = 0
  L = length(kers)
  if (("all" %in% kernels) | L == 9) {
    ALL = TRUE
    ker = "all"
    kers = setdiff(all.kernels, "uniform")
    L = 9
    opt.methods = rep("grid", 9)
  }
  if (print.info) {
    if(length(opt.methods) == 1){
      if(opt.methods == "grid"){
        cat("* Optimization method: grid search.\n* Parameters:")
        if (is.null(param.list)) {
          eps = seq(min.eps, max.eps, length.out = n.eps)
          cat(
            "\n\t-",
            n.eps,
            "values of epsilon in [",
            min.eps,
            ",",
            max.eps,
            "]\n\t-",
            n.cv,
            "fold cross-validation for each epsilon\n"
          )
        }
        else{
          eps = param.list
          n.eps = length(eps)
          cat(
            "\n\t-",
            length(param.list),
            "values of epsilon in [",
            min(param.list),
            ",",
            max(param.list),
            "]\n\t-",
            n.cv,
            "fold cross-validation for each epsilon\n")
        }
      }
      else{
        cat("* Optimization method: gradient descent.\n* Parameters:")
        if (is.numeric(rate.GD)) {
          cat("\n\t- learning rate:", rate.GD)
        }
        else{
          cat("\n\t- learning rate is of", rate.GD, "speed")
        }
        cat(
          "\n\t- initial try:",
          n.tries,
          "\n\t- stopping threshold:",
          err.GD,
          "\n\t- maximum iteration:",
          max.iter,
          "\n"
        )
      }
    }
    else{
      cat("* Optimization method: mixed.\n  ~ Parameters for grid search:")
      if (is.null(param.list)) {
        eps = seq(min.eps, max.eps, length.out = n.eps)
        cat(
          "\n\t-",
          n.eps,
          "values of epsilon in [",
          min.eps,
          ",",
          max.eps,
          "]\n\t-",
          n.cv,
          "fold cross-validation for each epsilon\n"
        )
      }
      else{
        eps = param.list
        n.eps = length(eps)
        cat(
          "\n\t-",
          length(param.list),
          "values of epsilon in [",
          min(param.list),
          ",",
          max(param.list),
          "]\n\t-",
          n.cv,
          "fold cross-validation for each epsilon\n")
      }
      cat("  ~ Parameters for gradient descent:")
      if (is.numeric(rate.GD)) {
        cat("\n\t- learning rate:", rate.GD)
      }
      else{
        cat("\n\t- learning rate is of", rate.GD, "speed")
      }
      cat(
        "\n\t- initial try:",
        n.tries,
        "\n\t- stopping threshold:",
        err.GD,
        "\n\t- maximum iteration:",
        max.iter,
        "\n"
      )
    }
    if ((L > 1) | ALL) {
      if (!ALL) {
        if (L == 2)
          cat("\t-", kers[1], "and", kers[2], "kernels are used\n")
        else{
          cat("\t- ")
          cat(kers[1:(L - 1)], sep = ", ")
          cat(" and", kers[L], "kernels are used\n")
        }
      }
      else
        cat("\t-", ker, "kernels are used\n")
    }
    else
      cat("\t-", kers, "kernel is used\n")
  }
  
  # Predictions on testing data
  predict.test = function(x,
                          ep,
                          n.M,
                          kerns = c(
                            "naive",
                            "uniform",
                            "epanechnikov",
                            "biweight",
                            "triweight",
                            "triangular",
                            "compact_gaussian",
                            "gaussian",
                            "exp4",
                            "cauchy")) {
    dif = sweep(train.machines, 2, x)
    pred.naive = function(eps, M) {
      id = rowSums((abs(dif) <= eps))
      val = sum(train.responses[id == M]) / sum(id == M)
      return(ifelse(is.na(val), MU, val))
    }
    
    pred.epan = function(eps) {
      val = rowSums(dif ^ 2) / eps ^ 2
      id0 = val <= 1
      wgt = (1 - val)[id0]
      res = train.responses[id0] %*% wgt / sum(wgt)
      return(ifelse(is.na(res), MU, res))
    }
    
    pred.bi = function(eps) {
      val = rowSums(dif ^ 2) / eps ^ 2
      id0 = val <= 1
      wgt = ((1 - val) ^ 2)[id0]
      res = train.responses[id0] %*% wgt / sum(wgt)
      return(ifelse(is.na(res), MU, res))
    }
    
    pred.tri = function(eps) {
      val = rowSums(dif ^ 2) / eps ^ 2
      id0 = val <= 1
      wgt = ((1 - val) ^ 3)[id0]
      res = train.responses[id0] %*% wgt / sum(wgt)
      return(ifelse(is.na(res), MU, res))
    }
    
    pred.triang = function(eps) {
      val = rowSums(abs(dif)) / eps ^ 2
      id0 = val <= 1
      wgt = (1 - val)[id0]
      res = train.responses[id0] %*% wgt / sum(wgt)
      return(ifelse(is.na(res), MU, res))
    }
    
    pred.CGaus = function(eps) {
      val = rowSums(dif ^ 2) / (2 * (sigma * eps) ^ 2)
      id0 = val <= rho
      wgt = exp(-val[id0])
      res = train.responses[id0] %*% wgt / sum(wgt)
      return(ifelse(is.na(res), MU, res))
    }
    
    pred.gaus = function(eps) {
      wgt = exp(-rowSums(dif ^ 2) / (2 * (sigma * eps) ^ 2))
      res = wgt %*% train.responses / sum(wgt)
      return(ifelse(is.na(res), MU, res))
    }
    
    pred.exp4 = function(eps) {
      val =  rowSums(dif ^ 4) / (2 * (sigma * eps) ^ 4)
      wgt = exp(-val)
      wgt[is.na(wgt)] = 0
      res = train.responses %*% wgt / sum(wgt)
      return(ifelse(is.na(res), MU, res))
    }
    
    pred.cauchy = function(eps) {
      wgt =  1 / (1 + rowSums(dif ^ 2) / eps ^ 2)
      wgt[is.na(wgt)] = 0
      res = train.responses %*% wgt / sum(wgt)
      return(ifelse(is.na(res), MU, res))
    }
    kerns = match.arg(kerns)
    pred = switch(
      kerns,
      naive = pred.naive(ep, n.M),
      uniform = pred.naive(ep, n.M),
      epanechnikov = pred.epan(ep),
      biweight = pred.bi(ep),
      triweight = pred.tri(ep),
      triangular = pred.triang(ep),
      compact_gaussian = pred.CGaus(ep),
      gaussian = pred.gaus(ep),
      exp4 = pred.exp4(ep),
      cauchy = pred.cauchy(ep)
    )
    return(pred)
  }
  if(length(opt.methods) == 1){
    parameter = sapply(kers,
                       FUN = function(kern)
                         return(
                           optimizeSmoothParameter(
                             train.machines = train.machines,
                             train.responses = train.responses,
                             test.machines = test.machines,
                             opt.method = opt.methods,
                             set.GD = setGD(initial.ep = initial.GD, 
                                            n.tries = n.tries,
                                            rate = rate.GD,
                                            threshold = err.GD,
                                            max.iter = max.iter,
                                            n.cv = n.cv,
                                            print.step = print.GD.step,
                                            figure = figure,
                                            coef.auto = coef.auto,
                                            coef.log = coef.log,
                                            coef.sqrt = coef.sqrt,
                                            coef.lm = coef.lm,
                                            deg.poly = deg.poly,
                                            base.exp = base.exp),
                             set.Grid = setGrid(min.ep = min.eps, 
                                                max.ep = max.eps, 
                                                n.ep = n.eps, 
                                                parameters = param.list),
                             set.Kernels = setKernels(kernels = kern, sigma = sigma, rho = rho)
                           )
                         )
    )
  }
  else{
    if(length(opt.methods) == L){
      parameter = mapply(
        FUN = function(kern, opt.meth)
          return(
            optimizeSmoothParameter(
              train.machines = train.machines,
              train.responses = train.responses,
              test.machines = test.machines,
              opt.method = opt.meth,
              set.GD = setGD(
                initial.ep = initial.GD,
                n.tries = n.tries,
                rate = rate.GD,
                threshold = err.GD,
                max.iter = max.iter,
                n.cv = n.cv,
                print.step = print.GD.step,
                figure = figure,
                coef.auto = coef.auto,
                coef.log = coef.log,
                coef.sqrt = coef.sqrt,
                coef.lm = coef.lm,
                deg.poly = deg.poly,
                base.exp = base.exp
              ),
              set.Grid = setGrid(
                min.ep = min.eps,
                max.ep = max.eps,
                n.ep = n.eps,
                parameters = param.list
              ),
              set.Kernels = setKernels(
                kernels = kern,
                sigma = sigma,
                rho = rho
              )
            )
          ),
        kers,
        opt.methods
      )
    }
    else{
      opt.methods = rep("grid", L)
      parameter = mapply(
        FUN = function(kern, opt.meth)
          return(
            optimizeSmoothParameter(
              train.machines = train.machines,
              train.responses = train.responses,
              test.machines = test.machines,
              opt.method = opt.meth,
              set.GD = setGD(
                initial.ep = initial.GD,
                n.tries = n.tries,
                rate = rate.GD,
                threshold = err.GD,
                max.iter = max.iter,
                n.cv = n.cv,
                print.step = print.GD.step,
                figure = figure,
                coef.auto = coef.auto,
                coef.log = coef.log,
                coef.sqrt = coef.sqrt,
                coef.lm = coef.lm,
                deg.poly = deg.poly,
                base.exp = base.exp
              ),
              set.Grid = setGrid(
                min.ep = min.eps,
                max.ep = max.eps,
                n.ep = n.eps,
                parameters = param.list
              ),
              set.Kernels = setKernels(
                kernels = kern,
                sigma = sigma,
                rho = rho
              )
            )
          ),
        kers,
        opt.methods
      )
      warning("opt.methods and kernels have different lengths. 
              Grid search method will be implemented.")
    }
  }
  
  if (print.train.results){
    cat("\n\n* Training summary:\n")
    print(t(parameter))
  }
  predictions = t(as.matrix(apply(
    test.machines,
    1,
    FUN = function(x) {
      res = mapply(function(ker, ep, n.M)
        return(predict.test(x, ep, n.M, ker)),
        kers, parameter["epsilon", ], parameter["n.machines", ])
    }
  )))
  if ((L == 1) & !ALL) {
    predictions = t(predictions)
    colnames(predictions) = kers
  }
  if(ALL)
    kers = "all"
  return(list(
    prediction = predictions,
    machines = test.machine0,
    epsilon = unlist(parameter["epsilon", ]),
    opt.train.risk = unlist(parameter["cv.error", ]),
    full.pred = rbind(train.machines, test.machines),
    train.split.input = x.train,
    train.split.output = y.train,
    kernels = kers
  ))
}


# ======================================================================================================================

# *** New version of regression MixCOBRA 
# ======================================
# Input:
# ------
# - train.input         : matrix or data frame of training data.
# - test.input          : matrix or data frame of testing data.
# - train.response      : vector of training response.
# - train.machines      : matrix or data frame of training machines.
# - test.machines       : matrix or data frame of testing machines.
# - test.response       : vector of testing response. It is null by default.
# - scale.input         : boolean type defining whether or not to scale the input. It must be a subset of {"normalization","standardization"}.
# - set.Machines        : "setMachines" function to define the parameters of the machines to be built.
# - param.list          : matrix or list containing parameters (alpha, beta) of the method.
# - min.alpha           : minimum of alpha in the uniform grid.
# - max.alpha           : maximum of alpha in the uniform grid.
# - n.alpha             : size of the grid of alpha.
# - min.beta            : minimum of beta in the uniform grid.
# - max.beta            : maximum of beta in the uniform grid.
# - n.beta              : size of the grid of beta.
# - figure              : boolean type defining whether or not to plot the graph of error vs parameters.
# - print.info          : boolean type defining whether or not to print the progress of the method.
# - print.train.results : boolean type defining whether or not to print the observed values of the method.
# - n.cv                : number of the folds in cross-validation method.
# - kernels             : vector of string defining the kernels used. It must a subset of {"naive", "uniform", "epanechnikov", "biweight", "triweight", 
#                         "triangular", "compact_gaussian", "gaussian", "exp4", "cauchy"}. 
# - sigma               : parameter of Gaussian-type kernels.
# - rho                 : parameter of compact support Gaussian kernel.
# Output: list of predictions and optimal values of parameters.
# -------

mixCOBRAreg = function(train.input,
                       test.input,
                       train.responses,
                       train.machines = NULL,
                       test.machines = NULL,
                       test.response = NULL,
                       scale.input = NULL,
                       set.Machines = setMachines(build.machines = TRUE,
                                                  machines = NULL,
                                                  splits = NULL,
                                                  scale.machines = NULL,
                                                  k = 10, 
                                                  ntree = 300,
                                                  mtry = NULL),
                       param.list = NULL,
                       min.alpha = 1e-300,
                       max.alpha = 10,
                       n.alpha = 100,
                       min.beta = 1e-300,
                       max.beta = 10,
                       n.beta = 100,
                       figure = FALSE,
                       print.info = TRUE,
                       print.train.results = TRUE,
                       n.cv = 5,
                       kernels = "uniform",
                       sigma = 1,
                       rho = 1,
                       angle = 10) {
  train.input = as.matrix(train.input)
  test.input = as.matrix(test.input)
  ker.names = c(
    "naive",
    "uniform",
    "epanechnikov",
    "biweight",
    "triweight",
    "triangular",
    "compact_gaussian",
    "gaussian",
    "exp4",
    "cauchy"
  )
  pred.train0 = train.machines
  pred.test0 = test.machines
  if(!is.null(scale.input)){
    scale.input = match.arg(scale.input, c("standardization", "normalization"))
    if(scale.input == "standardization"){
      SD = apply(train.input, 2, sd)
      SD.inv = diag(1/SD)
      colnames(SD.inv) = colnames(train.input)
      train.input = train.input %*% SD.inv
      test.input = test.input %*% SD.inv
    }
    if(scale.input == "normalization"){
      max.inp = max(train.input)
      min.inp = min(train.input)
      dom = (max.inp - min.inp)
      train.input = (train.input - min.inp + .1)/dom
      test.input = (test.input - min.inp + .1)/dom
    }
  }
  if(is.null(train.machines) | is.null(test.machines)){
    build = generateMachines(train.input = train.input, 
                             train.response = train.responses, 
                             test.input = test.input, 
                             machines = set.Machines$machines, 
                             splits = set.Machines$splits,
                             k = set.Machines$k, 
                             ntree = set.Machines$ntree,
                             mtry = set.Machines$mtry)
    train.machines = pred.train0 = build$train.machine
    train.responses = train.responses[build$id.val]
    train.input = train.input[build$id.val,]
    test.machines = pred.test0 = build$test.machine
  }
  
  if(!is.null(set.Machines$scale.machines)){
    scale_ = match.arg(set.Machines$scale.machines, c("normalization", "standardization"))
    if(scale_ == "standardization"){
      S2 = apply(train.machines, 2, sd)
      inv.S2 = diag(1/S2)
      colnames(inv.S2) = colnames(train.machines)
      train.machines = train.machines %*% inv.S2
      test.machines = test.machines %*% inv.S2
    }
    if(scale_ == "normalization"){
      max.sd = max(train.machines)
      train.machines = train.machines / max.sd
      test.machines = test.machines / max.sd
    }
  }
  
  if (is.null(param.list)) {
    alpha = seq(min.alpha, max.alpha, length.out = n.alpha)
    beta = seq(min.beta, max.beta, length.out = n.beta)
  }
  else{
    if(is.matrix(param.list)){
      if (dim(param.list)[2] > 2) {
        alpha = param.list[1,]
        beta = param.list[2,]
      }
      else{
        alpha = param.list[, 1]
        beta = param.list[, 2]
      }
    }
    if(is.list(param.list)){
      alpha = param.list[[1]]
      beta = param.list[[2]]
    }
  } 
  ALL = FALSE
  if (("all" %in% kernels) | ("All" %in% kernels) | ("ALL" %in% kernels)){
    M = 9
    kern = ker.names[-2]
    ALL = TRUE
  }
  else{
    kern = kernels
    M = length(kernels)
  }
  cat("\n*** MixCOBRA Method ***\n=======================\n")
  if(print.info){
    if (ALL) {
      cat(
        "* Parameters:\n\t-",
        n.alpha,
        "values of alpha in [",
        min(alpha),
        ",",
        max(alpha),
        "]\n\t-",
        n.beta,
        "values of beta in [",
        min(beta),
        ",",
        max(beta),
        "]\n\t-",
        n.cv,
        "folds cross-validation for each (alpha, beta)")
    }
    else{
      cat(
        "* Parameters:\n\t-",
        n.alpha,
        "values of alpha in [",
        min(alpha),
        ",",
        max(alpha),
        "]\n\t-",
        n.beta,
        "values of beta in [",
        min(beta),
        ",",
        max(beta),
        "]\n\t-",
        n.cv,
        "folds cross-validation for each pair of (alpha, beta)\n")
    }
    if ((M > 1) | ALL) {
      if (!ALL) {
        if (M == 2)
          cat("\t-", kern[1], "and", kern[2], "kernels are used.\n")
        else{
          cat("\t- ")
          cat(kern[1:(M - 1)], sep = ", ")
          cat(" and", kern[M], "kernels are used.\n")
        }
      }
      else
        cat("\t- All kernels are used.\n")
    }
    else
      cat("\t-", kern, "kernel is used.\n")
  }
  d1 = dim(train.input)
  d2 = dim(train.machines)
  
  # ================================ to be modified ==================================
  
  # The kernel functions
  naive.weight = function(L, al, be, ob) {
    ob0 = abs(L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))
    w = rowSums(ob0 <= 1)
    resp = train.responses[id.remain]
    res = sapply(1:(d1[2]+d2[2]), FUN = function(t){
      y.hat = mean(resp[w == t])
      tem = (ifelse(is.na(y.hat), 0, y.hat) - ob)^2
      return(tem)
    }, simplify = "vector")
    return(res)
  }
  epan.weight = function(L, al, be, ob) {
    val = rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^2)
    id0 = val < 1
    res = (1 - val[id0])
    y.hat = ifelse(sum(res) == 0, 0, res %*% train.responses[id.remain][id0] / sum(res))
    return((y.hat - ob)^2)
  }
  bi.weight = function(L, al, be, ob) {
    val = rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^2)
    id0 = val < 1
    res = (1 - val[id0])^2
    y.hat = ifelse(sum(res) == 0, 0, res %*% train.responses[id.remain][id0] / sum(res))
    return((y.hat - ob)^2)
  }
  tri.weight = function(L, al, be, ob) {
    val = rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^2)
    id0 = val < 1
    res = (1 - val[id0])^3
    y.hat = ifelse(sum(res) == 0, 0, res %*% train.responses[id.remain][id0] / sum(res))
    return((y.hat - ob)^2)
  }
  triang.weight = function(L, al, be, ob) {
    val = rowSums(abs(L %*% diag(1/rep(c(al,be),c(d1[2],d2[2])))))
    id0 = val < 1
    res = (1 - val[id0])
    y.hat = ifelse(sum(res) == 0, 0, res %*% train.responses[id.remain][id0] / sum(res))
    return((y.hat - ob)^2)
  }
  cGaus.weight = function(L, al, be, ob) {
    val = rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^2)/(2*sigma^2)
    id0 = val < rho
    res = exp(-val[id0])
    y.hat = ifelse(sum(res) == 0, 0, res %*% train.responses[id.remain][id0] / sum(res))
    return((y.hat - ob)^2)
  }
  gaus.weight = function(L, al, be, ob) {
    val = rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^2)/(2*sigma^2)
    res = exp(-val)
    y.hat = ifelse(sum(res) == 0, 0, res %*% train.responses[id.remain] / sum(res))
    return((y.hat - ob)^2)
  }
  exp4.weight = function(L, al, be, ob) {
    val =  rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^4)/(2*sigma^4)
    res = exp(-val)
    y.hat = ifelse(sum(res) == 0, 0, res %*% train.responses[id.remain] / sum(res))
    return((y.hat - ob)^2)
  }
  cauchy.weight = function(L, al, be, ob) {
    res =  1 / (1 + rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^2))
    y.hat = ifelse(sum(res) == 0, 0, res %*% train.responses[id.remain] / sum(res))
    return((y.hat - ob)^2)
  }
  
  all.kernels = list(naive.weight, naive.weight, epan.weight, bi.weight, tri.weight, triang.weight, 
                     cGaus.weight, gaus.weight, exp4.weight, cauchy.weight)
  names(all.kernels) = ker.names
  n.folds = floor(d2[1] / n.cv)
  CV = sample(d2[1])
  results = list()
  error = list()
  for(p in 1:M){
    ker = match.arg(kern[p], ker.names)
    weight = all.kernels[[ker]]
    kern[p] = ker
    error[[ker]] = 0
    if(print.info){
      cat("\n* ",ker,"kernel:\n   ~ processing")
      Per = 100L/n.cv
      t = 1
    } 
    for(k in 1:n.cv) {
      q = ifelse(k < n.cv, (n.folds * k), d2[1])
      id.test = CV[((k - 1) * n.folds + 1):q]
      id.remain = setdiff(CV, id.test)
      er = sapply(beta, FUN = function(bet){
        er0 = sapply(alpha, FUN = function(alp) {
          er1 = sapply(id.test, function(idd) {
            temp = cbind(sweep(train.input[id.remain,], 2, train.input[idd,]), 
                         sweep(train.machines[id.remain, ], 2, train.machines[idd,]))
            res0 = weight(temp, alp, bet, train.responses[idd])
            return(res0)
          }, simplify = "matrix")
          if (ker %in% c("naive", "uniform")){
            return(rowSums(er1))
          }
          else{
            return(sum(er1))
          }
        }, simplify = "array")
        return(er0)
      }, simplify = "array")
      error[[ker]] = error[[ker]] + er
      if(print.info){
        if(n.cv < 11){
          if(k < n.cv)
            cat(paste0("...",Per*k,"%"))
          else
            cat(paste0("...100%."))
        }
        else{
          gap = round(t*n.cv/10)
          if((k == gap) & (k < n.cv)){
            cat(paste0("...",t*10,"%"))
            t = t + 1
          }
          if(k == n.cv)
            cat(paste0("...100%."))
        }
      }
    }
    error[[ker]] = error[[ker]]/n.cv
    id.op = which(error[[ker]] == min(error[[ker]]), arr.ind = TRUE)
    if (ker %in% c("naive", "uniform")){
      results[[ker]]$alpha = alpha[id.op[2]]
      results[[ker]]$beta = beta[id.op[3]]
      results[[ker]]$opt.error = min(error[[ker]])
      results[[ker]]$n.machines = id.op[1]
    }
    else{
      results[[ker]]$alpha = alpha[id.op[1]]
      results[[ker]]$beta = beta[id.op[2]]
      results[[ker]]$opt.error = min(error[[ker]])
    }
  }
  # =============================== prediction part ============================= #
  
  # For prediction:
  naive.weight.pred = function(L, al, be, m = 0) {
    ob = abs(L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))
    w = rowSums(ob  <= 1)
    res = mean(train.responses[w == m])
    return(ifelse(is.na(res), 0, res))
  }
  epan.weight.pred = function(L, al, be, m = 0) {
    val = rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^2)
    id0 = val < 1
    res = (1 - val[id0])
    return(ifelse(sum(res) == 0, 0, res %*% train.responses[id0] / sum(res)))
  }
  bi.weight.pred = function(L, al, be, m = 0) {
    val = rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^2)
    id0 = val < 1
    res = (1 - val[id0])^2
    return(ifelse(sum(res) == 0, 0, res %*% train.responses[id0] / sum(res)))
  }
  tri.weight.pred = function(L, al, be, m = 0) {
    val = rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^2)
    id0 = val < 1
    res = (1 - val[id0])^3
    return(ifelse(sum(res) == 0, 0, res %*% train.responses[id0] / sum(res)))
  }
  triang.weight.pred = function(L, al, be, m = 0) {
    val = rowSums(abs(L %*% diag(1/rep(c(al,be),c(d1[2],d2[2])))))
    id0 = val < 1
    res = (1 - val[id0])
    return(ifelse(sum(res) == 0, 0, res %*% train.responses[id0] / sum(res)))
  }
  cGaus.weight.pred = function(L, al, be, m = 0) {
    val = rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^2)/(2*sigma^2)
    id0 = val < rho
    res = exp(-val[id0])
    return(ifelse(sum(res) == 0, 0, res %*% train.responses[id0] / sum(res)))
  }
  gaus.weight.pred = function(L, al, be, m = 0) {
    val = rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^2)/(2*sigma^2)
    res = exp(-val)
    return(ifelse(sum(res) == 0, 0, res %*% train.responses / sum(res)))
  }
  exp4.weight.pred = function(L, al, be, m = 0) {
    val =  rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^4)/(2*sigma^4)
    res = exp(-val)
    return(ifelse(sum(res) == 0, 0, res %*% train.responses / sum(res)))
  }
  cauchy.weight.pred = function(L, al, be, m = 0) {
    res =  1 / (1 + rowSums((L %*% diag(1/rep(c(al,be),c(d1[2],d2[2]))))^2))
    return(ifelse(sum(res) == 0, 0, res %*% train.responses / sum(res)))
  }
  
  all.kernels.pred = list(naive.weight.pred, naive.weight.pred, epan.weight.pred, bi.weight.pred, tri.weight.pred, triang.weight.pred, 
                          cGaus.weight.pred, gaus.weight.pred, exp4.weight.pred, cauchy.weight.pred)
  names(all.kernels.pred) = ker.names
  prep = cbind(train.input, train.machines)
  pred = sapply(kern,FUN = function(ker){
    weight = all.kernels.pred[[ker]]
    L = sapply(1:nrow(test.input), function(id) {
      obj = c(test.input[id,],test.machines[id,])
      temp = sweep(prep, 2, obj)
      return(weight(temp, results[[ker]]$alpha, results[[ker]]$beta, results[[ker]]$n.machines))
    }, simplify = "vector")
  }, simplify = "matrix")
  colnames(pred) = kern
  if(!is.null(test.response)){
    if(ncol(pred) > 1)
      err = colMeans((sweep(pred, 1, test.response))^2)
    else
      err = mean((pred-test.response)^2)
  }
  param = data.frame(alpha = sapply(kern, FUN = function(ker) results[[ker]]$alpha, simplify = "vector"),
                     beta =sapply(kern, FUN = function(ker) results[[ker]]$beta, simplify = "vector"),
                     n.machine = sapply(kern, FUN = function(ker) ifelse(is.null(results[[ker]]$n.machines),"NULL",results[[ker]]$n.machines), simplify = "vector"),
                     opt.cv.error = sapply(kern, FUN = function(ker) results[[ker]]$opt.error, simplify = "vector"))
  if (print.train.results) {
    cat("\n\n* Training summary:\n")
    print(param)
  }
  # ===============================
  
  if (figure) {
    if(M < 4)
      par(mfrow=c(1,M))
    else{
      if(M == 4)
        par(mfrow=c(2,2))
      else{
        if(M < 7)
          par(mfrow=c(2,floor(M/2)+1))
        else{
          par(mfrow=c(3,3))
        }
      }
    }
    for(i in 1:M){
      ker = kern[i]
      if(ker %in% c("uniform", "naive")){
        persp(
          x = alpha,
          y = beta,
          z = error[[ker]][results[[ker]]$n.machines,,],
          main = paste("MSE of", ker, "kernel"),
          theta = angle,
          xlab = "Alpha",
          ylab = "Beta",
          zlab = "MSE",
          col = "lightblue3"
        ) -> img3d
        points(
          trans3d(
            x = results[[ker]]$alpha,
            y = results[[ker]]$beta,
            z = results[[ker]]$opt.error,
            pmat = img3d
          ),
          col = "red3",
          pch = "*",
          cex = 2
        )
      }
      else{
        persp(
          x = alpha,
          y = beta,
          z = error[[ker]],
          main = paste("MSE of", ker, "kernel"),
          theta = angle,
          xlab = "Alpha",
          ylab = "Beta",
          zlab = "MSE",
          col = "lightblue3"
        ) -> img3d
        points(
          trans3d(
            x = results[[ker]]$alpha,
            y = results[[ker]]$beta,
            z = results[[ker]]$opt.error,
            pmat = img3d
          ),
          col = "red3",
          pch = "*",
          cex = 2
        )
      }
    }
  }
  if(ALL)
    kern = "all"
  if(is.null(test.response))
    return(
      list(
        prediction = pred,
        pred.machines = list(pred.tain = pred.train0,
                             pred.test0),
        optimal.alpha = param$alpha,
        optimal.beta = param$beta,
        kernels = kern
      )
    )
  else
    return(
      list(
        test.mse = err,
        prediction = pred,
        pred.machines = list(pred.tain = pred.train0,
                             pred.test0),
        optimal.alpha = param$alpha,
        optimal.beta = param$beta,
        kernels = kern
      )
    )
}



# ================================================================================
# ======================== Crieteria for each step of KFC ========================
# ================================================================================

# * Kstep: setting parameters for K-means algorithm in the K-step of the procedure
# ================================================================================
# Input:
# ------
# - K         : the number of clusters in K-means algorithm.
# - max.iter  : maximum iteration of K-means.
# - n.start   : the number of independent runs of the algorithm to avoid the local clustering structure.
# - method    : the divergence used in the algorithm.
# - deg       : degree of polynomial divergence.
# - threshold : the threshold to stop the algorithm when the change of clusters is less than this value.
# - figure    : boolean type define whether or not to plot the result of the algorithm in 2D or 3D cases.
# Output: list of all the inputs.
# -------

Kstep = function(K, max.iter = 300, n.start = 5, method = NULL, deg = 3, 
                 threshold = 1e-10, figure = FALSE){
  return(list(K = K, max.iter = max.iter, n.start = n.start, method = method,
              deg = deg, threshold = threshold, figure = FALSE))
}

# ===============================================================================================================================
# Fstep: setting the parameters for the F-step of the procedure
# =============================================================
# Input:
# ------
# - formula      : the formula of linear regression on each cluster.
# - single.model : boolean type whether or not to peform multiple linear regression with all the predictors without clustering just to compare with the procedure.
# Output: list of the two inputs.
# -------
Fstep = function(formula = NULL, single.model = FALSE){
  return(list(formula = formula, single.model = single.model))
}

# ===============================================================================================================================
# * Cstep: setting the parameters for the C-step of the procedure
# ===============================================================
# Input:
# ------
# - method      : the consensual aggregation method to perform in the C-step.
# - opt.method  : the optimization method in estimating the window parameter. It is either "grid" or "GD" or a vetor of the same
#                 length as the "kernels" which indicates the optimization algo for each kernel.
# - kernels     : vector of kernels used in the combination. It must be a subset of {"naive","uniform","epanechnikov","biweight",
#                 "triweight", "triangular", "exp4", "compact_gaussian", "gaussian", "cauchy"} or "all" which means that all the kernels
#                 will be used. The case of "uniform" or "naive" is the clasical COBRA method.
# - paramters   : parameters to be used in grid search algorithm.
# ....
# others are parameters for kernel-based COBRA or MixCOBRA...
# ...
# figure        : boolean type defining whether or not to plot the graph of error vs the parameter.
# show.info     : boolean type defining whether or not to print the observed information and the progress of the algorithm.

Cstep = function(method = c("cobra", "mixcobra"),
                 opt.methods = c("grid", "GD"),
                 kernels = "uniform",
                 parameters = NULL,
                 min.eps = 1e-10,
                 max.eps = 50,
                 n.eps = 100,
                 min.alpha = 1e-10,
                 max.alpha = 10,
                 n.alpha = 100,
                 min.beta = 1e-10,
                 max.beta = 10,
                 n.beta = 100,
                 n.cv = 5,
                 sigma = 1,
                 rho = 1,
                 figure = FALSE,
                 show.info = FALSE){
  return(list(method = method,
              opt.methods = opt.methods,
              kernels = kernels,
              parameters = parameters,
              min.eps = min.eps,
              max.eps = max.eps,
              n.eps = n.eps,
              min.alpha = min.alpha,
              max.alpha = max.alpha,
              n.alpha = n.alpha,
              min.beta = min.beta,
              max.beta = max.beta,
              n.beta = n.beta,
              n.cv = n.cv,
              sigma = sigma,
              rho = rho,
              figure = figure,
              show.info = show.info))
}


# ===============================================================
# =============== KFC procedure for regression ==================
# ===============================================================
# Input:
# ------
# - train.input         : matrix or data frame of the training data.
# - train.responses     : vector of response variable.
# - test.input          : matrix or data frame of the testing data.
# - test.responses      : vector of testing data. If it is NULL, the output will the predictions of the this vector.
# - scale.input.cluster : boolean type defining whether or not to scale the input before clsutering. It is either "normalizatoin" or "standardization".
# - scale.input.combine : boolean type defining whether or not to scale the input before combine. It is either "normalizatoin" or "standardization".
# - scale.machines      : boolean type defining whether or not to scale the machines. It is either "normalizatoin" or "standardization".
# - K_step              : "Kstep" function to set parameters for the K-step.
# - F_step              : "Fstep" function to set parameters for the F-step.
# - C_step              : "Cstep" function to set parameters for the C-step.
# - set_GD              : "setGD" fucntion to set parameters for gradient descent algorithm.
# Output: list of predictions, error, or/and values of observed parameters and the running time.
# -------

KFCreg = function(train.input,
                  train.responses,
                  test.input,
                  test.responses = NULL,
                  scale.input.cluster = NULL,
                  scale.input.combine = NULL,
                  scale.machines = NULL,
                  K_step = Kstep(K),
                  F_step = Fstep(formula = NULL),
                  C_step = Cstep(kernels = "uniform"),
                  set_GD = setGD(n.tries = 10, max.iter = 500),
                  angle = 10) {
  start.time <- Sys.time()
  x.train = as.matrix(train.input)
  P = ncol(x.train)
  y.train = train.responses
  n.train = length(y.train)
  
  x.test = as.matrix(test.input)
  y.test = test.responses
  x1.test = cbind(1, x.test)
  n.test = nrow(x.test)
  all.kmeans.methods = c("euclidean",
                         "gkl",
                         "logistic",
                         "itakura",
                         "polynomial")
  if(!is.null(scale.input.cluster)){
    scale.input.cluster = match.arg(scale.input.cluster, c("standardization", "normalization"))
    if(scale.input.cluster == "standardization"){
      SD = apply(x.train, 2, sd)
      SD.inv = diag(1/SD)
      colnames(SD.inv) = colnames(x.train)
      x.train = x.train %*% SD.inv
      x.test = x.test %*% SD.inv
    }
    if(scale.input.cluster == "normalization"){
      max.inp = max(x.train)
      min.inp = min(x.train)
      dom = (max.inp - min.inp)
      x.train = (x.train - min.inp + .1)/dom
      x.test = (x.test - min.inp + .1)/dom
    }
  }
  # K step: Kmeans clustering with BDs
  if (is.null(K_step$method)) {
    cluster.methods = all.kmeans.methods
    kms = lapply(cluster.methods, FUN = function(meth){
      km = kmeansBregman(x.train, K = K_step$K,
                         n.start = K_step$n.start,
                         maxIter = K_step$max.iter,
                         deg = K_step$deg,
                         method = meth,
                         epsilon = K_step$threshold,
                         figure = K_step$figure)
      return(list(clusters = km$clusters, centroids = km$centroids, method = km$method, deg = NULL))
    })
    METH = 5
    warning("No divergence provided! All of them are used!")
  }
  else{
    cluster.methods = K_step$method
    if(("polynomial" %in% cluster.methods) | ("poly" %in% cluster.methods)){
      kms1 = list()
      for(de in 1:length(K_step$deg)){
        km = kmeansBregman(x.train, K = K_step$K,
                           n.start = K_step$n.start,
                           maxIter = K_step$max.iter,
                           deg = K_step$deg[de],
                           method = "polynomial",
                           epsilon = K_step$threshold,
                           figure = K_step$figure)
        kms1[[de]] = list(clusters = km$clusters, centroids = km$centroids, method = km$method, 
                          deg = K_step$deg[de])
      }
      kms2 = lapply(setdiff(cluster.methods,"polynomial"), FUN = function(meth){
        km = kmeansBregman(x.train, K = K_step$K,
                           n.start = K_step$n.start,
                           maxIter = K_step$max.iter,
                           deg = K_step$deg,
                           method = meth,
                           epsilon = K_step$threshold,
                           figure = K_step$figure)
        return(list(clusters = km$clusters, centroids = km$centroids, method = km$method, deg = NULL))
      })
      kms = union(kms2, kms1)
    }
    else{
      kms = lapply(cluster.methods, FUN = function(meth){
        km = kmeansBregman(x.train, K = K_step$K,
                           n.start = K_step$n.start,
                           maxIter = K_step$max.iter,
                           deg = K_step$deg,
                           method = meth,
                           epsilon = K_step$threshold,
                           figure = K_step$figure)
        return(list(clusters = km$clusters, centroids = km$centroids, method = km$method, deg = NULL))
      })
    }
    METH = length(kms)
  }
  clusters = lapply(kms, function(km){
    return(list(clusters = findClosestCentroid(x.test, km$centroids, km$method, deg = km$deg)))
  })
  
  if (is.null(K_step$method)){
    k.meth.name = all.kmeans.methods
    names(kms) = k.meth.name
    names(clusters) = k.meth.name
  }
  else{
    if(("polynomial" %in% cluster.methods) | ("poly" %in% cluster.methods)){
      poly.id = which(c(cluster.methods == "polynomial",cluster.methods == "poly"))
      k.meth.name = c(setdiff(K_step$method, cluster.methods[poly.id]), paste0("poly",K_step$deg))
      names(kms) = k.meth.name
      names(clusters) = k.meth.name
    }
    else{
      k.meth.name = K_step$method
      names(kms) = k.meth.name
      names(clusters) = k.meth.name
    }
  }
  # F step: Fitting the corresponding model on each observed cluster
  model = fitClusterRegressionModel(x.train, y.train, kms, formula = F_step$formula)
  pred.train = y.hat.train = matrix(model$Fitted.values, nrow = n.train, byrow = FALSE)
  colnames(y.hat.train) = colnames(pred.train) = k.meth.name
  # Single model without clustering
  if(is.null(F_step$formula)){
    single.model = lm(y.train ~ ., data = as.data.frame(x.train))
    pred = predict(single.model, as.data.frame(x.test))
    single.mse = mean((pred - y.test) ^ 2)
  }
  else{
    single.model = lm(F_step$formula, data = as.data.frame(x.train))
    pred = predict(single.model, as.data.frame(x.test))
    single.mse = mean((pred - y.test) ^ 2)
  }
  
  mse = c()
  y.hat.test = matrix(0, n.test, METH)
  for (j in 1:METH) {
    clust = clusters[[j]]$clusters
    for (l in 1:K_step$K) {
      test.id = which(clust == l)
      x.temp = as.matrix(x1.test[test.id,])
      if (ncol(x.temp) != length(model$Coefficients[[j]][[l]]))
        x.temp = t(x.temp)
      y.hat.test[test.id, j] = x.temp %*% as.numeric(model$Coefficients[[j]][[l]])
    }
  }
  pred.test = y.hat.test
  BOTH = FALSE
  # C step: Consensual regression aggregation method with kernel-based COBRA
  if (METH > 1) {
    if (C_step$method %in% c("cobra", "cob", "COB", "COBRA")){
      pred = kernelCOBRA(train.design = y.hat.train,
                         train.response = y.train,
                         test.design = y.hat.test,
                         opt.methods = C_step$opt.methods,
                         scale.input = scale.machines,
                         set.Machines = setMachines(build.machines = FALSE,
                                                    scale.machines = NULL),
                         set.GD = setGD(initial.ep = set_GD$initial.ep, 
                                        n.tries = set_GD$n.tries,
                                        rate = set_GD$rate,
                                        threshold = set_GD$threshold,
                                        max.iter = set_GD$max.iter,
                                        n.cv = C_step$n.cv,
                                        print.step = set_GD$print.step,
                                        figure = C_step$figure,
                                        coef.auto = set_GD$coef.auto,
                                        coef.log = set_GD$coef.log,
                                        coef.sqrt = set_GD$coef.sqrt,
                                        coef.lm = set_GD$coef.lm,
                                        deg.poly = set_GD$coef.poly,
                                        base.exp = set_GD$coef.exp),
                         set.Grid = setGrid(min.ep = C_step$min.eps,
                                            max.ep = C_step$max.eps,
                                            n.ep = C_step$n.eps,
                                            parameters = C_step$parameters),
                         set.Kernels = setKernels(kernels = C_step$kernels,
                                                  sigma = C_step$sigma, rho = C_step$rho),
                         print.info = C_step$show.info,
                         print.train.results = C_step$show.info
      )
      ALL = ifelse("all" %in% pred$kernels, TRUE, FALSE)
    }
    else {
      if (C_step$method %in% c("mix", "mixcobra","MixCOBRA", "MIXCOBRA")){
        pred = mixCOBRAreg(
          train.input = x.train,
          train.machines = y.hat.train,
          test.input = x.test,
          test.machines = y.hat.test,
          train.responses = y.train,
          kernels = C_step$kernels,
          scale.input = scale.input.combine,
          set.Machines = setMachines(build.machines = FALSE,
                                     scale.machines = scale.machines),
          param.list = C_step$parameters,
          min.alpha = C_step$min.alpha,
          max.alpha = C_step$max.alpha,
          n.alpha = C_step$n.alpha,
          min.beta = C_step$min.beta,
          max.beta = C_step$max.beta,
          n.beta = C_step$n.beta,
          rho = C_step$rho,
          figure = C_step$figure,
          n.cv = C_step$n.cv,
          sigma = C_step$sigma,
          print.info = C_step$show.info,
          print.train.results = C_step$show.info,
          angle = angle
        )
        ALL = ifelse("all" %in% pred$kernels, TRUE, FALSE)
      }
      else{
        if(C_step$method %in% c("both", "all", "mixed")){
          pred.cob = kernelCOBRA(train.design = y.hat.train,
                                 train.response = y.train,
                                 test.design = y.hat.test,
                                 opt.methods = C_step$opt.methods,
                                 scale.input = scale.machines,
                                 set.Machines = setMachines(build.machines = FALSE,
                                                            scale.machines = NULL),
                                 set.GD = setGD(initial.ep = set_GD$initial.ep, 
                                                n.tries = set_GD$n.tries,
                                                rate = set_GD$rate,
                                                threshold = set_GD$threshold,
                                                max.iter = set_GD$max.iter,
                                                n.cv = C_step$n.cv,
                                                print.step = set_GD$print.step,
                                                figure = C_step$figure,
                                                coef.auto = set_GD$coef.auto,
                                                coef.log = set_GD$coef.log,
                                                coef.sqrt = set_GD$coef.sqrt,
                                                coef.lm = set_GD$coef.lm,
                                                deg.poly = set_GD$coef.poly,
                                                base.exp = set_GD$coef.exp),
                                 set.Grid = setGrid(min.ep = C_step$min.eps,
                                                    max.ep = C_step$max.eps,
                                                    n.ep = C_step$n.eps,
                                                    parameters = C_step$parameters),
                                 set.Kernels = setKernels(kernels = C_step$kernels,
                                                          sigma = C_step$sigma, rho = C_step$rho),
                                 print.info = C_step$show.info,
                                 print.train.results = C_step$show.info
          )
          pred.mix = mixCOBRAreg(
            train.input = x.train,
            train.machines = y.hat.train,
            test.input = x.test,
            test.machines = y.hat.test,
            train.responses = y.train,
            kernels = C_step$kernels,
            scale.input = scale.input.combine,
            set.Machines = setMachines(build.machines = FALSE,
                                       scale.machines = scale.machines),
            param.list = C_step$parameters,
            min.alpha = C_step$min.alpha,
            max.alpha = C_step$max.alpha,
            n.alpha = C_step$n.alpha,
            min.beta = C_step$min.beta,
            max.beta = C_step$max.beta,
            n.beta = C_step$n.beta,
            rho = C_step$rho,
            figure = C_step$figure,
            n.cv = C_step$n.cv,
            sigma = C_step$sigma,
            print.info = C_step$show.info,
            print.train.results = C_step$show.info,
            angle = angle
          )
          BOTH = TRUE
          ALL = ifelse("all" %in% pred.mix$kernels, TRUE, FALSE)
        }
        else 
          warning("Combining method is not found!")
      }
    }
  }
  else{
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    if (is.null(y.test)) {
      colnames(pred.test) = k.meth.name
      return(list(prediction = pred.test, Run.time = time.taken))
    }
    else{
      return(list(
        Table = data.frame(MSE = mean((pred.test - y.test) ^ 2), row.names = k.meth.name),
        Run.time = time.taken
      ))
    }
  }
  name.all = c( "naive",
                "uniform",
                "epanechnikov",
                "biweight",
                "triweight",
                "triangular",
                "exp4",
                "compact_gaussian",
                "gaussian",
                "cauchy")
  colnames(y.hat.test) = colnames(pred.test) = k.meth.name
  if (ALL) {
    if (is.null(y.test)) {
      if(BOTH){
        y.hat.test = list(pred.candidate = pred.test, 
                          pred.cobra = pred.cob$prediction, 
                          pred.mixcobra = pred.mix$prediction)
      }
      else{
        colnames(pred$prediction) = name.all[-2]
        y.hat.test = cbind(pred.test, pred$prediction)
      }
    }
    else{
      if(BOTH){
        name = c(k.meth.name, paste0(name.all[-2], ".cob"), paste0(name.all[-2], ".mix"))
        mse = colMeans((sweep(as.matrix(pred.test), 1, y.test)) ^ 2)
        names(mse) = k.meth.name
        cob.mse = colMeans((sweep(pred.cob$prediction, 1, y.test)) ^ 2)
        names(cob.mse) = colnames(pred.cob$prediction)
        mix.mse = colMeans((sweep(pred.mix$prediction, 1, y.test)) ^ 2)
        names(mix.mse) = colnames(pred.mix$prediction)
        y.hat.test = list(pred.candidate = pred.test, 
                          pred.cobra = pred.cob$prediction, 
                          pred.mixcobra = pred.mix$prediction)
        tab = data.frame(
          mse.test = c(mse, cob.mse, mix.mse),
          mse.single = single.mse,
          row.names = name
        )
      }
      else{
        name = c(k.meth.name, name.all[-2])
        mse = colMeans((sweep(as.matrix(pred.test), 1, y.test)) ^ 2)
        names(mse) = k.meth.name
        cob.mse = colMeans((sweep(pred$prediction, 1, y.test)) ^ 2)
        names(cob.mse) = colnames(pred$prediction)
        y.hat.test = cbind(pred.test, pred$prediction)
        names(y.hat.test) = name
        tab = data.frame(
          mse.test = c(mse, cob.mse),
          mse.single = single.mse,
          row.names = name
        )
      }
    }
  }
  else{
    if (is.null(y.test)) {
      if(BOTH){
        colnames(pred.cob$prediction) = colnames(pred.mix$prediction) = C_step$kernels
        y.hat.test = list(pred.candidate = pred.test, 
                          pred.cobra = pred.cob$prediction, 
                          pred.mixcobra = pred.mix$prediction)
      }
      else{
        y.hat.test = cbind(pred.test, pred$prediction)
      }
    }
    else{
      if(BOTH){
        name = c(k.meth.name, paste0(C_step$kernels, ".cob"), paste0(C_step$kernels, ".mix"))
        colnames(pred.cob$prediction) = colnames(pred.mix$prediction) = C_step$kernels
        mse = colMeans((sweep(as.matrix(pred.test), 1, y.test)) ^ 2)
        names(mse) = k.meth.name
        cob.mse = colMeans((sweep(pred.cob$prediction, 1, y.test)) ^ 2)
        names(cob.mse) = paste0(C_step$kernels, ".cob")
        mix.mse = colMeans((sweep(pred.mix$prediction, 1, y.test)) ^ 2)
        names(mix.mse) = paste0(C_step$kernels, ".mix")
        y.hat.test = list(pred.candidate = pred.test,  
                          pred.cobra = pred.cob$prediction, 
                          pred.mixcobra = pred.mix$prediction)
        tab = data.frame(
          mse.test = c(mse, cob.mse, mix.mse),
          mse.single = single.mse,
          row.names = name
        )
      }
      else{
        mse = colMeans((sweep(as.matrix(pred.test), 1, y.test)) ^ 2)
        names(mse) = k.meth.name
        com.mse = colMeans((sweep(pred$prediction, 1, y.test)) ^ 2)
        names(com.mse) = colnames(pred$prediction)
        name = c(k.meth.name, C_step$kernels)
        y.hat.test = cbind(pred.test, pred$prediction)
        names(y.hat.test) = name
        test.mse = c(mse, com.mse)
        tab = data.frame(
          mse.test = test.mse,
          mse.single = single.mse,
          row.names = name
        )
      }
    }
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  if (is.null(y.test)) {
    return(list(train.predict = pred.train,
                test.predict = y.hat.test, 
                run.time = time.taken))
  }
  else{
    return(list(
      table = tab,
      train.predict = pred.train,
      test.predict = y.hat.test,
      run.time = time.taken
    ))
  }
}


# ============================ Application ========================

x1 = c(rnorm(500,1,3),rnorm(500,5,2), rnorm(500,-5,5))
x2 = c(rnorm(500,3,3),rnorm(500,-5,3), rnorm(500,5,2))
x3 = c(rnorm(500,-3,2),rnorm(500,1,3), rnorm(500,0,2))
tab = data.frame(x1,x2,x3)
y = (-2*exp(x1)+3*exp(-x2)-5*sin(x3)+3*log(abs(x1+x2+x3)))[1:500]
y = c(y, (log(abs(rowSums(tab)))+5*sin(x1)-2*x2^2+x3^2)[501:1000])
y = c(y, (exp(-x2)-cos(x1^2)/cos(x3)+2*exp(x2+x3))[1001:1500])


res = c()
for(i in 1:10){
train = sample(1500,0.8*1500)
kfc = KFCreg(train.input = tab[train,], train.responses = y[train], test.input = tab[-train,],
             test.responses = y[-train],
             #  scale.input = "standard",
             scale.machines = "standard",
             K_step = Kstep(K = 3, max.iter = 500, method = c("exp" ,"euclidean", "polynomial"), deg = c(3,4)),
             F_step = Fstep(single.model = TRUE),
             C_step = Cstep(method = "mixed", kernels = c("gaus"), 
                            opt.methods = c("GD"),
                            max.eps = 3, n.eps = 300, figure = TRUE, show.info = TRUE,
                            min.alpha = 0.001, min.beta = 0.001, min.eps = 0.001, max.alpha = 3,
                            max.beta = 5, n.alpha = 30, n.beta = 30,
                            n.cv = 10),
             set_GD = setGD(print.step = T, rate = "log", coef.log = 0.3, figure = T, max.iter = 1000, threshold = 1e-15))

library(randomForest)
rf = randomForest(tab[train,], y = y[train], xtest = tab[-train,], ntree = 500)

library(gbm)
boost = gbm(y[train]~., data = tab[train,], n.trees = 500, distribution = "gaussian")
pred = predict.gbm(boost, tab[-train,])

# MSE
mse = c(kfc$table$mse.test, mean((rf$test$predicted - y[-train])^2), mean((pred - y[-train])^2))


names(mse) = c(rownames(kfc$table), "randomForest", "boost")
res = rbind(res, mse)
}