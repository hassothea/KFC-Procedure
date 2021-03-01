
# ================
# From: 27/01/2021
# ================

# All Bregman divergences
# =======================
euclidDist = function(x, y){
  return((x-y)%*% (x-y))
}

gklDist = function(x, y){
  return(sum(x * log(abs(x/y)) - (x-y)))
}

logisticDist = function(x, y){
  x0 = x/sum(abs(x))
  y0 = y/sum(abs(y))
  x1 = 1 - x0
  y1 = 1 - y0
  return(sum(x0*log(abs(x0/y0))+x1*log(abs(x1/y1))))
}

itakuraDist = function(x, y){
  ratio = x/y
  return(sum(ratio-log(abs(ratio)) - 1))
}

expDist = function(x, y){
  y0 = exp(y)
  return(sum(exp(x)-y0-(x-y)*y0))
}

polyDist = function(x, y, deg){
  if(deg %% 2 == 0){
    res_dis = sum(x^deg - y^deg - deg * (x-y) * y^(deg - 1))
  }
  else{
    res_dis = sum(x^deg - y^deg - deg * (x-y) * sign(y) * y^(deg - 1))
  }
  return(res_dis)
}

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

# * F-step: Fitting linear logistic regression model on each cluster
# ==================================================================
# Inputs:
# -------
# - data          : matrix or data frame containing the dataset.
# - target        : response varaible.
# - methods       : list of K-means algorithms givem by the function "kmeansBregman".
# - formula       : forumula for linear regression model fitting on each cluster. It is a class of formula and NULL by default.
# Output:
# -------
# - table         : table containing all the measures of fit such as accuracy, balanced accuracy, precision, recall and f1 score.
# - coefficients  : list of coefficients of all constructed local models.
# - fitted.prob   : probabilities of being in class 0 and 1 given by each candidate model corresponding to each Bregman divergence.
# - fitted.y      : predicted classes given by each candidate model corresponding to each Bregman divergence.
# - running.Time  : the running time of the algortihm.

fitClusterLogitModel = function(data, target, methods, formula = NULL){
  start.time <- Sys.time()
  require(stats)
  K = nrow(methods[[1]]$centroids)
  m = length(methods)
  N = length(target)
  if(!is.factor(target))
    target = as.factor(target)
  if(is.null(formula))
    form = formula(temp.y ~ .)
  else
    form = formula
  err = c()
  coef = list()
  p.fit = matrix(0, N, m)
  y.fit = matrix(0, N, m)
  for(i in 1:m){
    cluster = as.numeric(methods[[i]]$clusters)
    coef[[i]] = list()
    for(j in 1:K){
      id = which(cluster == j)
      n1 = length(id)
      temp.x = data[id,]
      temp.y = target[id]
      logistic.Model = glm(formula = form, binomial(link='logit'), data=as.data.frame(temp.x))
      p.fit[id,i] = logistic.Model$fitted.values
      y.fit[id,i] = ifelse(p.fit[id,i] >= 0.5,1,0)
      coef[[i]][[j]] = logistic.Model$coefficients
    }
    err = c(err, mean(target != y.fit[,i]))
  }
  
  name = c("euclidean", "gkl", "logistic", "itakura-saito")
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  return(list(error = err, 
              coefficients = coef,
              fitted.prob = p.fit, 
              fitted.y = y.fit,
              running.time = time.taken))
}

# ======================================================================================================================

# * generateMachines: generate machines to be combined in consensual aggregation method
# =====================================================================================
# Input:
# ------
# - train.input    : matrix or data frame of training input data.
# - train.response : vector of training response.
# - test.input     : matrix or data frame of testing data.
# - machines       : choices of machines which must be a subset of {"lda", "qda", "knn", "tree", "rf", "logit"}.
# - splits         : the proportion taking value between (0,1) to split the data. By default it is NULL and all the five machiens are built.
# - k              : the paramter "k" of kNN method.
# - ntree          : the paramter "ntree" of rf method.
# Output:
# -------
# - train.machine  : matrix of predictions of training data given by all the machines.
# - test.machine   : matrix of predictions of testing data given by all the machines.
# - id.val         : the random index of the data after shuffled.
# - machines       : the vector of machines used.

generateMachines.class = function(train.input, train.response, test.input, machines = NULL, 
                                  splits = NULL, k = 10, ntree = 300, mtry = NULL){
  require(tree)
  require(randomForest)
  require(FNN)
  require(MASS)
  require(stats)
  d = dim(train.input)
  COBRA.lda <- function(x) {
    a <- lda(as.factor(train.response.1)~., as.data.frame(train.input.1))
    res <- predict(a, as.data.frame(x))$class
    return(res)
  }
  COBRA.qda <- function(x) {
    a <- qda(as.factor(train.response.1)~., as.data.frame(train.input.1))
    res <- predict(a, as.data.frame(x))$class
    return(res)
  }
  COBRA.tree <- function(x) {
    a <- tree(as.factor(train.response.1)~.,data = as.data.frame(train.input.1))
    res <- predict(a, as.data.frame(x), type = "class")
    return(res)
  }
  COBRA.randomForest <- function(x) {
    if(is.null(mtry))
      a <- randomForest(x = as.data.frame(train.input.1), y = as.factor(train.response.1), ntree = ntree)
    else
      a <- randomForest(x = as.data.frame(train.input.1), y = as.factor(train.response.1), ntree = ntree, mtry = mrty)
    res <- predict(a, as.data.frame(x), type = "class")
    return(res)
  }
  COBRA.knn <- function(x) {
    a <- knn(train = train.input.1, test = x, cl = as.factor(train.response.1), k = k, prob = T)
    res = a[1:length(a)]
    return(res)
  }
  COBRA.logit = function(x){
    a <- glm(as.factor(train.response.1)~., family = binomial(link='logit'), data=as.data.frame(train.input.1))
    res = predict.glm(a, as.data.frame(x), type = "response")
    return(as.factor(as.numeric(res > 0.5)))
  }
  
  all.machines = list(lda = COBRA.lda, qda = COBRA.qda, knn = COBRA.knn, tree = COBRA.tree, rf = COBRA.randomForest, logit = COBRA.logit)
  
  if(is.null(machines))
    mach = c("lda", "qda", "knn", "tree", "rf", "logit")
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
  return(list(train.machine = res.tr - 1, test.machine = res.te - 1, id.val = id.rand[(I+1):n1],
              machines = mach))
}


# =============================================================================================================================================================


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

setGridClass = function(min.h = 1e-100, max.h = NULL, n.h = 300, parameters = NULL, n.cv = 10){
  return(list(min.h = min.h,
              max.h = max.h,
              n.h = n.h,
              n.cv = n.cv,
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

setKernels = function(kernels = NULL, sigma = 1, rho = 1, full.naive = TRUE, expo = 2){
  return(list(kernels = kernels, sigma = sigma, rho = rho, full.naive = full.naive, expo = expo))
}


# ======================================================================================================================

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

kernelAggregateClassification = function(train.design,
                                         train.response,
                                         test.design,
                                         test.response = NULL,
                                         scale.input = NULL,
                                         set.Machines = setMachines(build.machines = TRUE,
                                                                    machines = NULL,
                                                                    splits = NULL,
                                                                    scale.machines = NULL,
                                                                    k = 10, ntree = 300,
                                                                    mtry = NULL),
                                         set.Grid = setGridClass(min.h = 1e-300,
                                                                 max.h = NULL,
                                                                 n.h = 300,
                                                                 n.cv = 10,
                                                                 parameters = NULL),
                                         set.Kernels = setKernels(
                                           kernels = c(
                                             "naive",
                                             "uniform",
                                             "epanechnikov",
                                             "biweight",
                                             "triweight",
                                             "triangular",
                                             "compact_exp",
                                             "exponential",
                                             "cauchy"),
                                           sigma = 1,
                                           rho = 1,
                                           expo = 2,
                                           full.naive = FALSE), 
                                         print.info = TRUE,
                                         print.train.results = TRUE,
                                         figure = TRUE) {
  d = dim(train.design)
  # parameters Grid search
  min.h = set.Grid$min.h
  max.h = set.Grid$max.h
  n.h = set.Grid$n.h
  n.cv = set.Grid$n.cv
  param.list = set.Grid$parameters
  
  # Kernels
  kernels = set.Kernels$kernels
  sigma = set.Kernels$sigma
  rho = set.Kernels$rho
  expo = set.Kernels$expo
  
  cat("\n*** Kernel-based Consensual Classification Aggregation ***\n==========================================================\n")
  
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
  
  # Building machines
  if(set.Machines$build.machines){
    M = ifelse(is.null(set.Machines$machines), 5, length(set.Machines$machines))
    build = generateMachines.class(train.input = train.machines, train.response = train.responses,
                                   test.input = test.machines, machines = set.Machines$machines,
                                   splits = set.Machines$splits,
                                   k = set.Machines$k, 
                                   ntree = set.Machines$ntree,
                                   mtry = set.Machines$mtry)
    train.machines = train.machine0 = build$train.machine
    train.responses = train.response[build$id.val]
    test.machines = test.machine0 = build$test.machine
    max.sd = max(train.machines)
    d = dim(train.machines)
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
  if(is.null(max.h))
    max.h = 1.5*max.sd*ncol(train.machines)
  
  # All kernels are used
  #=====================
  MU = 0
  ker.names = c(
    "naive",
    "uniform",
    "epanechnikov",
    "biweight",
    "triweight",
    "triangular",
    "compact_exp",
    "exponential",
    "cauchy"
  )
  if(is.null(param.list)){
    param.list = seq(min.h, max.h, length.out = n.h)
  }
  ALL = FALSE
  kers = kernels
  ker = 0
  L = length(kers)
  if (("all" %in% kernels) | L == 8) {
    ALL = TRUE
    ker = "all"
    kers = setdiff(ker.names, "uniform")
    L = 8
  }
  if (print.info) {
    cat("* Optimization method: grid search.\n* Parameters:")
    cat(
      "\n\t-",
      n.h,
      "values of h in [",
      min.h,
      ",",
      max.h,
      "]\n\t-",
      n.cv,
      "fold cross-validation for each h\n"
    )
    if ((length(intersect(kers, c("naive", "uniform"))) > 0) & set.Kernels$full.naive)
      cat("\t- full-machine is used for naive/unifrom kernel\n")
    else{
      if (length(intersect(kers, c("naive", "uniform"))) > 0)
        cat("\t- alpha-machine is used for naive/unifrom kernel\n")
    }
    if ((L > 1) | ALL) {
      if (!ALL) {
        if (L == 2)
          cat("\t-", kers[1], "and", kers[2], "kernels are used.\n")
        else{
          cat("\t- ")
          cat(kers[1:(L - 1)], sep = ", ")
          cat(" and", kers[L], "kernels are used.\n")
        }
      }
      else
        cat("\t-", ker, "kernels are used.\n")
    }
    else
      cat("\t-", kers, "kernel is used.\n")
  }
  major = as.numeric(names(table(train.responses)[1]))
  M = ncol(train.machines)
  # Predictions on testing data
  # ===========================
  f.naive = function(y){
    w = apply(train.k, 1, function(x) sum(x == y))
    res = sapply(1:M, FUN = function(t){
      y.hat = ifelse(mean(y.k[w == t]) > 0.5, 1, 0)
      return(ifelse(is.na(y.hat), 0, y.hat))
    }, simplify = "vector")
    return(res)
  }
  
  f.epan = function(y,h){
    val1 = apply(train.k, 1, function(x) sum(x != y))/h
    idx = which(val1 <= 1)
    w = (1 - (val1[idx])^2)
    y.hat = ifelse(y.k[idx] %*% w/sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat),  0,y.hat))
  }
  
  f.biw = function(y,h){
    val1 = apply(train.k, 1, function(x) sum(x != y))/h
    idx = which(val1 <= 1)
    w = (1 - (val1[idx])^2)^2
    y.hat = ifelse(y.k[idx] %*% w/sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat),  0,y.hat))
  }
  
  f.triw = function(y,h){
    val1 = apply(train.k, 1, function(x) sum(x != y))/h
    idx = which(val1 <= 1)
    w = (1 - (val1[idx])^2)^3
    y.hat = ifelse(y.k[idx] %*% w/sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat),  0, y.hat))
  }
  
  f.triang = function(y,h){
    val1 = apply(train.k, 1, function(x) sum(x != y))/h
    idx = which(val1 <= 1)
    w = (1 - val1[idx])
    y.hat = ifelse(y.k[idx] %*% w/sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat),  0, y.hat))
  }
  
  f.exp.compact = function(y,h){
    val1 = apply(train.k, 1, function(x) sum(x != y))
    val2 = (val1/(sigma*h))^expo
    indx = (val2 <= rho)
    w = exp(-val2[indx])
    y.hat = ifelse(y.k[indx] %*% w/sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat),  0, y.hat))
  }
  
  f.exp = function(y,h){
    val = apply(train.k, 1, function(x) sum(x != y))
    w = exp(-(val/(sigma*h))^expo)
    y.hat = ifelse(y.k %*% w/sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat),  0, y.hat))
  }
  
  f.cauchy = function(y,h){
    val1 = apply(train.k, 1, function(x) sum(x != y))
    w = 1/(1+(val1/h)^2)
    y.hat = ifelse(y.k %*% w/sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat),  0,y.hat))
  }
  
  f.naive.full.pred = function(y){
    w = apply(train.machines, 1, function(x) sum(x != y))
    y.hat = ifelse(mean(train.responses[w == M]) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat),  0, y.hat))
  }
  
  # ==================================================================================================
  
  all.kernels = list(f.naive, f.naive, f.epan, f.biw, f.triw, f.triang, f.exp.compact, f.exp, f.cauchy)
  names(all.kernels) = ker.names
  n.folds = floor(d[1] / n.cv)
  CV = sample(d[1])
  results = list()
  error = list()
  Remain = vector("logical", length = L)
  for(p in 1:L){
    ker = match.arg(kers[p], ker.names)
    agg = all.kernels[[ker]]
    kers[p] = ker
    if(ker %in% c("naive", "uniform") & set.Kernels$full.naive){
      error[[ker]] = NULL
      results[[ker]]$h = NULL
      results[[ker]]$n.machines = M
      results[[ker]]$opt.error = NULL
      if(print.info){
        cat("\n* ",ker,"kernel:\n   ~ processing")
        Per = 100L/n.cv
        t = 1
      }
      pred.full = apply(test.machines, 1, function(y){
        return(f.naive.full.pred(y))
      })
      Remain[p] = TRUE
      if(print.info)
        cat("...done!")
    }
    else{
      error[[ker]] = 0
      if(print.info){
        cat("\n* ",ker,"kernel:\n   ~ processing")
        Per = 100L/n.cv
        t = 1
      }
      for (k in 1:n.cv) {
        q = ifelse(k < n.cv, (n.folds * k), d[1])
        id.test = CV[((k - 1) * n.folds + 1):q]
        id.remain = setdiff(CV, id.test)
        train.k = train.machines[id.remain,]
        y.k = train.responses[id.remain]
        if(ker %in% c("naive", "uniform")){
          err0 = sapply(id.test, FUN = function(idd){
            err1 = agg(train.machines[idd,])
            return(as.numeric(err1 != train.responses[idd]))
          })
          er = rowMeans(err0)
        }
        else{
          er = sapply(param.list, FUN = function(h){
            er1 = sapply(id.test, function(idd) {
              res0 = agg(train.machines[idd,],h)
              return(as.numeric(res0 != train.responses[idd]))
            }, simplify = "matrix")
            return(mean(er1))
          }, simplify = "array")
        }
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
        error[[ker]] = error[[ker]] + er
      }
      error[[ker]] = error[[ker]]/n.cv
      id.op = which.min(error[[ker]])
      if (ker %in% c("naive", "uniform")){
        results[[ker]]$n.machines = id.op
        results[[ker]]$opt.error = min(error[[ker]])
      }
      else{
        results[[ker]]$h = param.list[id.op]
        results[[ker]]$opt.error = min(error[[ker]])
      }
    }
  }
  # =============================== prediction part ============================= #
  
  # For prediction:
  f.naive.pred = function(y,h){
    w = apply(train.machines, 1, function(x) sum(x == y))
    y.hat = ifelse(mean(train.responses[w == h]) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat),  0, y.hat))
  }
  
  f.epan.pred = function(y,h){
    val1 = apply(train.machines, 1, function(x) sum(x != y))/h
    idx = which(val1 <= 1)
    w = (1 - (val1[idx])^2)
    y.hat = ifelse(train.responses[idx] %*% w / sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat), major, y.hat))
  }
  
  f.biw.pred = function(y,h){
    val1 = apply(train.machines, 1, function(x) sum(x != y))/h
    idx = which(val1 <= 1)
    w = (1 - (val1[idx])^2)^2
    y.hat = ifelse(train.responses[idx] %*% w / sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat), major, y.hat))
  }
  
  f.triw.pred = function(y,h){
    val1 = apply(train.machines, 1, function(x) sum(x != y))/h
    idx = which(val1 <= 1)
    w = (1 - (val1[idx])^2)^3
    y.hat = ifelse(train.responses[idx] %*% w / sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat), major, y.hat))
  }
  
  f.triang.pred = function(y,h){
    val1 = apply(train.machines, 1, function(x) sum(x != y))/h
    idx = which(val1 <= 1)
    w = (1 - val1[idx])
    y.hat = ifelse(train.responses[idx] %*% w / sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat), major, y.hat))
  }
  
  f.exp.compact.pred = function(y,h){
    val1 = apply(train.machines, 1, function(x) sum(x != y))
    val2 = (val1/(sigma*h))^expo
    indx = (val2 <= rho)
    w = exp(-val2[indx])
    y.hat = ifelse(train.responses[indx] %*% w / sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat), major, y.hat))
  }
  
  f.exp.pred = function(y,h){
    val = apply(train.machines, 1, function(x) sum(x != y))
    w = exp(-(val/(sigma*h))^expo)
    y.hat = ifelse(train.responses %*% w / sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat), major, y.hat))
  }
  
  f.cauchy.pred = function(y,h){
    val1 = apply(train.machines, 1, function(x) sum(x != y))
    w = 1/(1+(val1/h)^2)
    y.hat = ifelse(train.responses %*% w / sum(w) > 0.5, 1, 0)
    return(ifelse(is.na(y.hat), major, y.hat))
  }
  
  all.kernels.pred = list(f.naive.pred, f.naive.pred, f.epan.pred, f.biw.pred, f.triw.pred, 
                          f.triang.pred, f.exp.compact.pred, f.exp.pred, f.cauchy.pred)
  names(all.kernels.pred) = ker.names
  if((L == 1) & (sum(Remain) > 0)){
    pred = matrix(pred.full, nrow = nrow(test.design), byrow = TRUE)
    colnames(pred) = kers
  }
  else{
    pred0 = sapply(kers[!Remain],FUN = function(ker){
      final.agg = all.kernels.pred[[ker]]
      if(ker %in% c("naive", "uniform"))
        op.h = results[[ker]]$n.machines
      else
        op.h = results[[ker]]$h
      pre0 = apply(test.machines, 1, function(y) {
        return(final.agg(y, op.h))
      })
    }, simplify = "matrix")
    pred = matrix(nrow = nrow(test.design), ncol = L, data = NA)
    if(sum(Remain) > 0){
      pred[,Remain] = pred.full
      j = 1
      for(i in 1:L){
        if(!Remain[i]){
          pred[,i] = pred0[,j]
          j = j + 1
        }
      }
    }
    else
      pred = pred0
    colnames(pred) = kers
  }
  if(!is.null(test.response)){
    err0 = colMeans(sweep(test.machines, 1, test.response, FUN = "!="))
    if(L > 1){
      err1 = colMeans(sweep(pred, 1, test.response, FUN = "!="))
    }
    else{
      err1 = mean(pred != test.response)
    }
    err = c(err0, err1)
  }
  param = data.frame(parameter = sapply(kers, FUN = function(ker) ifelse(is.null(results[[ker]]$h),"NULL",results[[ker]]$h), simplify = "vector"),
                     n.machine = sapply(kers, FUN = function(ker) ifelse(is.null(results[[ker]]$n.machines),"NULL",results[[ker]]$n.machines), simplify = "vector"),
                     opt.cv.error = sapply(kers, FUN = function(ker) ifelse(is.null(results[[ker]]$opt.error),"NULL", results[[ker]]$opt.error), simplify = "vector"))
  if (print.train.results) {
    cat("\n\n* Training summary:\n")
    print(param)
  }
  # ===============================
  ex = ifelse(sum(Remain) > 0, 1, 0)
  pic = L - ex
  if (figure & (pic > 0)){
    if(pic < 4)
      par(mfrow=c(1,pic))
    else{
      if(pic == 4)
        par(mfrow=c(2,2))
      else{
        if(pic < 7)
          par(mfrow=c(2,floor(pic/2)+1))
        else{
          par(mfrow=c(3,3))
        }
      }
    }
    for(i in 1:L){
      ker = kers[i]
      if(ker %in% c("uniform", "naive") & ex == 0){
        plot(1:M, error[[ker]], type = "l", col = "blue", cex = 1.5, main = ker, xlab = "No. machines", ylab = "Misclassification error")
        points(results[[ker]]$n.machines, results[[ker]]$opt.error, col = "red", cex = 1.5, pch = 19)
      }
      else{
        plot(param.list, error[[ker]], type = "l", col = "blue", cex = 1.5, main = ker, xlab = "Parameter", ylab = "Misclassification error")
        points(results[[ker]]$h, results[[ker]]$opt.error, col = "red", cex = 1.5, pch = 19)
      }
    }
    par(mfrow=c(1,1))
  }
  if(ALL)
    kers = "all"
  if(is.null(test.response))
    return(
      list(
        prediction = pred,
        pred.machines = list(pred.train = train.machine0,
                             pred.test = test.machine0),
        optimal.param = param$parameter,
        kernels = kers
      )
    )
  else
    return(
      list(
        test.error = err,
        prediction = pred,
        pred.machines = list(pred.train = train.machine0,
                             pred.test = test.machine0),
        optimal.param = param$parameter,
        kernels = kers
      )
    )
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

mixCOBRAclass = function(train.input,
                         test.input,
                         train.responses,
                         train.machines = NULL,
                         test.machines = NULL,
                         test.response = NULL,
                         scale.input = NULL,
                         set.Machines = setMachines(
                           build.machines = TRUE,
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
                         n.cv = 5,
                         set.Kernels = setKernels(
                           kernels = "naive",
                           sigma = 1,
                           rho = 1,
                           expo = 2,
                           full.naive = TRUE),
                         figure = FALSE,
                         print.info = FALSE,
                         print.train.results = FALSE,
                         angle = 10) {
  sigma = set.Kernels$sigma
  rho = set.Kernels$rho
  expo = set.Kernels$expo
  kernels = set.Kernels$kernels
  train.input = as.matrix(train.input)
  test.input = as.matrix(test.input)
  ker.names = c(
    "naive",
    "uniform",
    "epanechnikov",
    "biweight",
    "triweight",
    "triangular",
    "compact_exp",
    "exponential",
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
    build = generateMachines.class(train.input = train.input, train.response = train.responses,
                                   test.input = test.input, machines = set.Machines$machines,
                                   splits = set.Machines$splits,
                                   k = set.Machines$k, mtry = set.Machines$mtry,
                                   ntree = set.Machines$ntree)
    train.machines = pred.train0 = build$train.machine
    train.responses = train.responses[build$id.val]
    test.machines = pred.test0 = build$test.machine
    train.input = train.input[build$id.val,]
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
  train.machines = as.matrix(train.machines)
  test.machines = as.matrix(test.machines)
  d1 = dim(train.input)
  d2 = dim(train.machines)
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
    L = 8
    kers = ker.names[-2]
    ALL = TRUE
  }
  else{
    kers = kernels
    L = length(kernels)
  }
  cat("\n*** Classification MixCOBRA ***\n===============================\n")
  if(print.info){
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
      "folds cross-validation for each (alpha, beta)\n")
    if ((length(intersect(kers, c("naive", "uniform"))) > 0) & set.Kernels$full.naive)
      cat("\t- full-machine is used for naive/unifrom kernel\n")
    else{
      if (length(intersect(kers, c("naive", "uniform"))) > 0)
        cat("\t- alpha-machine is used for naive/unifrom kernel\n")
    }
    if ((L > 1) | ALL) {
      if (!ALL) {
        if (L == 2)
          cat("\t-", kers[1], "and", kers[2], "kernels are used.\n")
        else{
          cat("\t- ")
          cat(kers[1:(L - 1)], sep = ", ")
          cat(" and", kers[L], "kernels are used.\n")
        }
      }
      else
        cat("\t-", kers, "kernels are used.\n")
    }
    else
      cat("\t-", kers, "kernel is used.\n")
  }
  
  # ================================ to be modified ==================================
  M = d1[2]+d2[2]
  
  # The kernel functions
  f.naive = function(L,ob) {
    w = rowSums(abs(L) <= 1)
    resp = train.responses[id.remain]
    res = sapply(1:M, FUN = function(t){
      y.hat = mean(resp[w == t])
      tem = ifelse(y.hat>0.5, 1, 0)
      return(ifelse(is.na(tem),0,tem))
    }, simplify = "vector")
    return(as.numeric(res != ob))
  }
  f.naive.full = function(L, ob) {
    w = rowSums(abs(L) <= 1)
    resp = train.responses[id.remain]
    y.hat = mean(resp[w == M])
    tem = ifelse(y.hat>0.5, 1, 0)
    tem = ifelse(is.na(tem),0,tem)
    return(as.numeric(tem != ob))
  }
  f.epan = function(L, ob) {
    val = rowSums(L^2)
    idx = (val < 1)
    res = (1-val[idx])
    y.hat = ifelse(res %*% train.responses[id.remain][idx] / sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), 0, y.hat)
    return(as.numeric(y.hat != ob))
  }
  f.biw = function(L,ob) {
    val = rowSums(L^2)
    idx = (val < 1)
    res = (1-val[idx])^2
    y.hat = ifelse(res %*% train.responses[id.remain][idx] / sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), 0, y.hat)
    return(as.numeric(y.hat != ob))
  }
  f.triw = function(L, ob) {
    val = rowSums(L^2)
    idx = (val < 1)
    res = (1-val[idx])^3
    y.hat = ifelse(res %*% train.responses[id.remain][idx] / sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), 0, y.hat)
    return(as.numeric(y.hat != ob))
  }
  f.triang = function(L,ob) {
    val = rowSums(abs(L))
    idx = (val < 1)
    res = (1-val[idx])
    y.hat = ifelse(res %*% train.responses[id.remain][idx] / sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), 0, y.hat)
    return(as.numeric(y.hat != ob))
  }
  f.exp.compact = function(L,ob) {
    val = rowSums(L^2)/(2*sigma^2)
    idx = val < rho
    res = exp(-val[idx])
    y.hat = ifelse(res %*% train.responses[id.remain][idx] / sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), 0, y.hat)
    return(as.numeric(y.hat != ob))
  }
  f.exp = function(L, ob) {
    val = rowSums(L^2)/(2*sigma^2)
    res = exp(-val)
    y.hat = ifelse(res %*% train.responses[id.remain]/sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), 0, y.hat)
    return(as.numeric(y.hat != ob))
  }
  f.cauchy = function(L,ob) {
    res =  1 / (1 + rowSums(L^2))
    y.hat = ifelse(res %*% train.responses[id.remain]/sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), 0, y.hat)
    return(as.numeric(y.hat != ob))
  }
  
  # ==================================================================================================
  
  all.kernels = list(f.naive, f.naive, f.epan, f.biw, f.triw, f.triang, f.exp.compact, f.exp, f.cauchy)
  names(all.kernels) = ker.names
  n.folds = floor(d2[1] / n.cv)
  CV = sample(d2[1])
  results = list()
  error = list()
  FULL_ = FALSE
  Naive__ = 0
  for(p in 1:L){
    ker = match.arg(kers[p], ker.names)
    agg = all.kernels[[ker]]
    if(ker %in% c("naive", "uniform") & !set.Kernels$full.naive){
      Naive_ = TRUE
      Naive__ = TRUE
    }
    else
      Naive_ = FALSE
    if(ker %in% c("naive", "uniform") & set.Kernels$full.naive){
      FULL_ = TRUE
      agg = f.naive.full
    }
    kers[p] = ker
    error[[ker]] = 0
    if(print.info){
      cat("\n* ",ker,"kernel:\n   ~ processing")
      Per = 100L/n.cv
      count = 1
    } 
    for (k in 1:n.cv) {
      q = ifelse(k < n.cv, (n.folds * k), d2[1])
      id.test = CV[((k - 1) * n.folds + 1):q]
      id.remain = setdiff(CV, id.test)
      err = sapply(beta, FUN = function(bet){
        er0 = sapply(alpha, FUN = function(alp) {
          er1 = sapply(id.test, function(idd) {
            temp0 = cbind(sweep(train.input[id.remain,], 2, train.input[idd,]), 
                          sweep(train.machines[id.remain, ], 2, train.machines[idd,]))
            temp = temp0 %*% diag(1/rep(c(alp,bet),c(d1[2],d2[2])))
            res0 = agg(temp, train.responses[idd])
            return(res0)
          }, simplify = "matrix")
          if(Naive_){
            return(rowMeans(er1))
          }
          else{
            return(mean(er1))
          }
        }, simplify = "array")
        return(er0)
      }, simplify = "array")
      error[[ker]] = error[[ker]] + err
      if(print.info){
        if(n.cv < 11){
          if(k < n.cv)
            cat(paste0("...",round(Per*k,2),"%"))
          else
            cat(paste0("...100%."))
        }
        else{
          gap = round(t*n.cv/10)
          if((k == gap) & (k < n.cv)){
            cat(paste0("...",t*10,"%"))
            count = count + 1
          }
          if(k == n.cv)
            cat(paste0("...100%."))
        }
      }
    }
    error[[ker]] = error[[ker]]/n.cv
    id.op = which(error[[ker]] == min(error[[ker]]), arr.ind = TRUE)[1,]
    if (Naive_){
      results[[ker]]$alpha = alpha[id.op[2]]
      results[[ker]]$beta = beta[id.op[3]]
      results[[ker]]$opt.error = min(error[[ker]])
      results[[ker]]$n.machines = id.op[1]
    }
    else{
      if(FULL_)
        results[[ker]]$n.machines = M
      results[[ker]]$alpha = alpha[id.op[1]]
      results[[ker]]$beta = beta[id.op[2]]
      results[[ker]]$opt.error = min(error[[ker]])
    }
  }
  
  # =============================== prediction part ============================= #
  major = as.numeric(names(which.max(table(train.responses))))
  # For prediction:
  f.naive.pred = function(L, m = 0) {
    w = rowSums(abs(L) <= 1)
    y.hat = mean(train.responses[w == m])
    tem = ifelse(is.na(y.hat), major, ifelse(y.hat>0.5, 1, 0))
    return(tem)
  }
  f.epan.pred = function(L, m = 0) {
    val = rowSums(L^2)
    idx = (val < 1)
    res = (1-val[idx])
    y.hat = ifelse(res %*% train.responses[idx]/sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), major, y.hat)
    return(y.hat)
  }
  f.biw.pred = function(L, m = 0) {
    val = rowSums(L^2)
    idx = (val < 1)
    res = (1-val[idx])^2
    y.hat = ifelse(res %*% train.responses[idx]/sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), major, y.hat)
    return(y.hat)
  }
  f.triw.pred = function(L, m = 0) {
    val = rowSums(L^2)
    idx = (val < 1)
    res = (1-val[idx])^3
    y.hat = ifelse(res %*% train.responses[idx]/sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), major, y.hat)
    return(y.hat)
  }
  f.triang.pred = function(L, m = 0) {
    val = rowSums(abs(L))
    idx = (val < 1)
    res = (1-val[idx])
    y.hat = ifelse(res %*% train.responses[idx]/sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), major, y.hat)
    return(y.hat)
  }
  f.exp.compact.pred = function(L, m = 0) {
    val = rowSums(L^2)/(2*sigma^2)
    idx = val < rho
    res = exp(-val[idx])
    y.hat = ifelse(res %*% train.responses[idx]/sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), major, y.hat)
    return(y.hat)
  }
  f.exp.pred = function(L, m = 0) {
    val = rowSums(L^2)/(2*sigma^2)
    res = exp(-val)
    y.hat = ifelse(res %*% train.responses/sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), major, y.hat)
    return(y.hat)
  }
  f.cauchy.pred = function(L, m = 0) {
    res =  1 / (1 + rowSums(L^2))
    y.hat = ifelse(res %*% train.responses/sum(res) > 0.5, 1, 0)
    y.hat = ifelse(is.na(y.hat), major, y.hat)
    return(y.hat)
  }
  
  all.kernels.pred = list(f.naive.pred, f.naive.pred, f.epan.pred, f.biw.pred, 
                          f.triw.pred, f.triang.pred, f.exp.compact.pred, f.exp.pred, f.cauchy.pred)
  names(all.kernels.pred) = ker.names
  prep = cbind(train.input, train.machines)
  pred = sapply(kers,FUN = function(ker){
    agg.final = all.kernels.pred[[ker]]
    L = sapply(1:nrow(test.input), function(id){
      obj = c(test.input[id,],test.machines[id,])
      temp = sweep(prep, 2, obj) %*% diag(1/rep(c(results[[ker]]$alpha,results[[ker]]$beta),c(d1[2],d2[2])))
      return(agg.final(temp, results[[ker]]$n.machines))
    }, simplify = "vector")
  }, simplify = "matrix")
  colnames(pred) = kers
  err0 = colMeans(sweep(test.machines, 1, test.response, FUN = "!="))
  if(!is.null(test.response)){
    if(L > 1){
      err1 = colMeans(sweep(pred, 1, test.response, FUN = "!="))
    }
    else{
      err1 = mean(pred != test.response)
    }
    names(err1) = kers
  }
  err = c(err0, err1)
  param = data.frame(alpha = sapply(kers, FUN = function(ker) results[[ker]]$alpha, simplify = "vector"),
                     beta =sapply(kers, FUN = function(ker) results[[ker]]$beta, simplify = "vector"),
                     n.machine = sapply(kers, FUN = function(ker) ifelse(is.null(results[[ker]]$n.machines),"NULL",results[[ker]]$n.machines), simplify = "vector"),
                     opt.cv.error = sapply(kers, FUN = function(ker) results[[ker]]$opt.error, simplify = "vector"))
  if (print.train.results) {
    cat("\n\n* Training summary:\n")
    print(param)
  }
  # ===============================
  
  if (figure) {
    if(L < 4)
      par(mfrow=c(1,L))
    else{
      if(L == 4)
        par(mfrow=c(2,2))
      else{
        if(L < 7)
          par(mfrow=c(2,floor(L/2)+1))
        else{
          par(mfrow=c(3,3))
        }
      }
    }
    for(i in 1:L){
      ker = kers[i]
      if(Naive__){
        persp(
          x = alpha,
          y = beta,
          z = error[[ker]][results[[ker]]$n.machines,,],
          main = paste("Error of", ker, "kernel"),
          theta = angle,
          xlab = "Alpha",
          ylab = "Beta",
          zlab = "Misclassification error",
          col = "lightblue3"
        ) -> img3d
        points(
          trans3d(
            x = results[[ker]]$alpha,
            y = results[[ker]]$beta,
            z =  results[[ker]]$opt.error,
            pmat = img3d
          ),
          col = "red3",
          pch = "*",
          cex = 2
        )
        Naive__ = FALSE
      }
      else{
        persp(
          x = alpha,
          y = beta,
          z = error[[ker]],
          main = paste("Error of", ker, "kernel"),
          theta = angle,
          xlab = "Alpha",
          ylab = "Beta",
          zlab = "Misclassification error",
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
    par(mfrow=c(1,1))
  }
  if(ALL)
    kers = "all"
  if(is.null(test.response))
    return(
      list(
        prediction = pred,
        pred.machines = list(pred.tain = pred.train0,
                             pred.test = pred.test0),
        optimal.alpha = param$alpha,
        optimal.beta = param$beta,
        kernels = kers
      )
    )
  else
    return(
      list(
        test.error = err,
        prediction = pred,
        pred.machines = list(pred.tain = pred.train0,
                             pred.test = pred.test0),
        optimal.alpha = param$alpha,
        optimal.beta = param$beta,
        kernels = kers
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
                 kernels = "uniform",
                 parameters = NULL,
                 min.h = 1e-10,
                 max.h = 50,
                 n.h = 100,
                 min.alpha = 1e-10,
                 max.alpha = 10,
                 n.alpha = 100,
                 min.beta = 1e-10,
                 max.beta = 10,
                 n.beta = 100,
                 n.cv = 5,
                 sigma = 1,
                 rho = 1,
                 expo = 2,
                 full.naive = FALSE,
                 figure = FALSE,
                 show.info = FALSE){
  return(list(method = method,
              kernels = kernels,
              parameters = parameters,
              min.h = min.h,
              max.h = max.h,
              n.h = n.h,
              min.alpha = min.alpha,
              max.alpha = max.alpha,
              n.alpha = n.alpha,
              min.beta = min.beta,
              max.beta = max.beta,
              n.beta = n.beta,
              n.cv = n.cv,
              sigma = sigma,
              rho = rho,
              expo = expo,
              full.naive = full.naive,
              figure = figure,
              show.info = show.info))
}


# =======================================================================================
# Sigmoid function: logistic regression
# =====================================
sigmoid = function(x){
  return(1/(1+exp(-x)))
}

# ===============================================================
# =============== KFC procedure for regression ==================
# ===============================================================
# Input:
# ------
# - train.input     : matrix or data frame of the training data.
# - train.responses : vector of response variable.
# - test.input      : matrix or data frame of the testing data.
# - test.responses  : vector of testing data. If it is NULL, the output will the predictions of the this vector.
# - scale.input     : boolean type defining whether or not to scale the input. It is either "normalizatoin" or "standardization".
# - split_data      : the percentage of training data randomly splitted for clustering in K step. The remaining part will be used in the C step.
#                   : It is NULL by defalt and all the data will be used in both steps K and C step. If the given value is not in [0,1], the value of 0.5 will be used.
# - scale.machines  : boolean type defining whether or not to scale the machines. It is either "normalizatoin" or "standardization".
# - K_step          : "Kstep" function to set parameters for the K-step.
# - F_step          : "Fstep" function to set parameters for the F-step.
# - C_step          : "Cstep" function to set parameters for the C-step.
# - set_GD          : "setGD" fucntion to set parameters for gradient descent algorithm.
# Output: list of predictions, error, or/and values of observed parameters and the running time.
# -------

KFCclass = function(train.input,
                    train.responses,
                    test.input,
                    test.responses = NULL,
                    scale.input = NULL,
                    split_data = NULL,
                    K_step = Kstep(K = 3),
                    F_step = Fstep(formula = NULL),
                    C_step = Cstep(kernels = "uniform"),
                    set_GD = setGD(n.tries = 10, max.iter = 500),
                    angle = 10) {
  start.time <- Sys.time()
  dim_input = dim(train.input)
  x.train0 = x.train1 = x.train2 = as.matrix(train.input)
  y.train0 = y.train1 = y.train2 = train.responses
  n.train0 = n.train1 = length(y.train0)
  
  P = ncol(x.train0)
  x.test = as.matrix(test.input)
  y.test = test.responses
  x1.test = cbind(1, x.test)
  n.test = nrow(x.test)
  
  all.kmeans.methods = c("euclidean",
                         "gkl",
                         "logistic",
                         "itakura",
                         "polynomial")
  
  if(!is.null(scale.input)){
    scale.input = match.arg(scale.input, c("standardization", "normalization"))
    if(scale.input == "standardization"){
      SD = apply(x.train0, 2, sd)
      SD.inv = diag(1/SD)
      colnames(SD.inv) = colnames(x.train0)
      x.train1 = x.train0 %*% SD.inv
      x.test = x.test %*% SD.inv
    }
    if(scale.input == "normalization"){
      max.inp = max(x.train0)
      min.inp = min(x.train0)
      dom = (max.inp - min.inp)
      x.train1 = (x.train0 - min.inp + .1)/dom
      x.test = (x.test - min.inp + .1)/dom
    }
  }
  
  SPD = FALSE
  if(!(is.null(split_data))){
    SPD = TRUE
    if((split_data < 0) | (split_data > 1)){
      warning("The split is not between 0 and 1! 50% split is done.")
      split_data = 0.5
    } 
    train_id0 = sample(n.train0, floor(split_data * n.train0))
    x.train1 = x.train0[train_id0,]
    y.train1 = train.responses[train_id0]
    n.train1 = length(y.train1)
    x.train2 = x.train0[-train_id0,]
    y.train2 = train.responses[-train_id0]
  }
  
  # K step: Kmeans clustering with BDs
  if (is.null(K_step$method)) {
    cluster.methods = all.kmeans.methods
    kms = lapply(cluster.methods, FUN = function(meth){
      km = kmeansBregman(x.train1, K = K_step$K,
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
        km = kmeansBregman(x.train1, K = K_step$K,
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
        km = kmeansBregman(x.train1, K = K_step$K,
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
        km = kmeansBregman(x.train1, K = K_step$K,
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
  clusters2 = lapply(kms, function(km){
    return(list(clusters = findClosestCentroid(x.train2, km$centroids, km$method, deg = km$deg)))
  })
  
  clusters_test = lapply(kms, function(km){
    return(list(clusters = findClosestCentroid(x.test, km$centroids, km$method, deg = km$deg)))
  })
  
  if (is.null(K_step$method)){
    k.meth.name = all.kmeans.methods
    names(kms) = k.meth.name
    names(clusters2) = names(clusters_test) = k.meth.name
  }
  else{
    if(("polynomial" %in% cluster.methods) | ("poly" %in% cluster.methods)){
      poly.id = which(c(cluster.methods == "polynomial",cluster.methods == "poly"))
      k.meth.name = c(setdiff(K_step$method, cluster.methods[poly.id]), paste0("poly",K_step$deg))
      names(kms) = k.meth.name
      names(clusters2) = names(clusters_test) = k.meth.name
    }
    else{
      k.meth.name = K_step$method
      names(kms) = k.meth.name
      names(clusters2) = names(clusters_test) = k.meth.name
    }
  }
  
  # Single model without clustering
  if(is.null(F_step$formula)){
    single.model = glm(y.train0 ~ ., data = as.data.frame(x.train0), family = binomial(link='logit'))
    pred = predict.glm(single.model, as.data.frame(x.test), type = "response")
    single.error = mean(as.numeric(pred > 0.5) != y.test)
  }
  else{
    single.model = glm(F_step$formula, data = as.data.frame(x.train0), family = binomial(link='logit'))
    pred = predict.glm(single.model, as.data.frame(x.test), type = "response")
    single.error = mean(as.numeric(pred > 0.5) != y.test)
  }
  
  # F step: Fitting the corresponding model on each observed cluster
  model = fitClusterLogitModel(x.train1, y.train1, kms, formula = F_step$formula)
  pred.prob0 = y.hat.train0 = matrix(model$fitted.prob, nrow = length(y.train1), byrow = FALSE)
  pred.train.y0 = ifelse(y.hat.train0 > 0.5, 1, 0)
  colnames(y.hat.train0) = colnames(pred.train.y0) = k.meth.name
  
  mse = c()
  y.hat.test = matrix(0, n.test, METH)
  y.hat.train = matrix(0, length(y.train2), METH)
  colnames(y.hat.train) = colnames(y.hat.test) = k.meth.name
  x1.train2 = cbind(1, x.train2)
  for (j in 1:METH) {
    clust2 = clusters2[[j]]$clusters
    clu_test = clusters_test[[j]]$clusters
    for (l in 1:K_step$K) {
      test.id2 = which(clust2 == l)
      test.id_test = which(clu_test == l)
      x.temp2 = as.matrix(x1.train2[test.id2,])
      x.temp_test = as.matrix(x1.test[test.id_test,])
      # if (ncol(x.temp) != length(model$coefficients[[j]][[l]]))
      #   x.temp = t(x.temp)
      y.hat.train[test.id2, j] = sigmoid(x.temp2 %*% as.numeric(model$coefficients[[j]][[l]]))
      y.hat.test[test.id_test, j] = sigmoid(x.temp_test %*% as.numeric(model$coefficients[[j]][[l]]))
    }
  }
  pred.train.prob = y.hat.train
  pred.train.y = ifelse(y.hat.train > 0.5,1,0)
  pred.test.prob = y.hat.test
  pred.test.y = ifelse(y.hat.test > 0.5,1,0)
  BOTH = FALSE
  # C step: Consensual regression aggregation method with kernel-based COBRA
  if (METH > 1) {
    if (C_step$method %in% c("cobra", "cob", "COB", "COBRA")){
      pred = kernelAggregateClassification(train.design = pred.train.y,
                                           train.response = y.train2,
                                           test.design = pred.test.y,
                                           test.response = y.test,
                                           set.Machines = setMachines(build.machines = FALSE),
                                           set.Grid = setGridClass(min.h = C_step$min.h,
                                                                   max.h = C_step$max.h,
                                                                   n.h = C_step$n.h,
                                                                   parameters = C_step$parameters,
                                                                   n.cv = C_step$n.cv),
                                           set.Kernels = setKernels(kernels = C_step$kernels,
                                                                    sigma = C_step$sigma, 
                                                                    rho = C_step$rho,
                                                                    expo = C_step$expo,
                                                                    full.naive = C_step$full.naive),
                                           print.info = C_step$show.info,
                                           print.train.results = C_step$show.info,
                                           figure = C_step$figure
      )
      ALL = ifelse("all" %in% pred$kernels, TRUE, FALSE)
    }
    else {
      if (C_step$method %in% c("mix", "mixcobra","MixCOBRA", "MIXCOBRA")){
        pred = mixCOBRAclass(
          train.input = x.train2,
          train.machines = pred.train.y,
          test.input = x.test,
          test.machines = pred.test.y,
          train.responses = y.train2,
          test.response = y.test,
          scale.input = NULL,
          set.Machines = setMachines(build.machines = FALSE),
          set.Kernels = setKernels(kernels = C_step$kernels,
                                   sigma = C_step$sigma,
                                   rho = C_step$rho,
                                   expo = C_step$expo,
                                   full.naive = C_step$full.naive),
          param.list = C_step$parameters,
          min.alpha = C_step$min.alpha,
          max.alpha = C_step$max.alpha,
          n.alpha = C_step$n.alpha,
          min.beta = C_step$min.beta,
          max.beta = C_step$max.beta,
          n.beta = C_step$n.beta,
          figure = C_step$figure,
          n.cv = C_step$n.cv,
          print.info = C_step$show.info,
          print.train.results = C_step$show.info,
          angle = angle
        )
        ALL = ifelse("all" %in% pred$kernels, TRUE, FALSE)
      }
      else{
        if(C_step$method %in% c("both", "all", "mixed")){
          pred.cob = kernelAggregateClassification(train.design = pred.train.y,
                                                   train.response = y.train2,
                                                   test.design = pred.test.y,
                                                   test.response = y.test,
                                                   set.Machines = setMachines(build.machines = FALSE),
                                                   set.Grid = setGridClass(min.h = C_step$min.h,
                                                                           max.h = C_step$max.h,
                                                                           n.h = C_step$n.h,
                                                                           parameters = C_step$parameters,
                                                                           n.cv = C_step$n.cv),
                                                   set.Kernels = setKernels(kernels = C_step$kernels,
                                                                            sigma = C_step$sigma, 
                                                                            rho = C_step$rho,
                                                                            expo = C_step$expo,
                                                                            full.naive = C_step$full.naive),
                                                   print.info = C_step$show.info,
                                                   print.train.results = C_step$show.info,
                                                   figure = C_step$figure
          )
          pred.mix = mixCOBRAclass(
            train.input = x.train2,
            train.machines = pred.train.y,
            test.input = x.test,
            test.machines = pred.test.y,
            train.responses = y.train2,
            test.response = y.test,
            scale.input = NULL,
            set.Machines = setMachines(build.machines = FALSE),
            set.Kernels = setKernels(kernels = C_step$kernels,
                                     sigma = C_step$sigma,
                                     rho = C_step$rho,
                                     expo = C_step$expo,
                                     full.naive = C_step$full.naive),
            param.list = C_step$parameters,
            min.alpha = C_step$min.alpha,
            max.alpha = C_step$max.alpha,
            n.alpha = C_step$n.alpha,
            min.beta = C_step$min.beta,
            max.beta = C_step$max.beta,
            n.beta = C_step$n.beta,
            figure = C_step$figure,
            n.cv = C_step$n.cv,
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
      colnames(pred.test.y) = k.meth.name
      return(list(prediction = pred.test.y, 
                  run.time = time.taken))
    }
    else{
      return(list(
        table = data.frame(error = mean(pred.test.y != y.test), 
                           row.names = k.meth.name),
        run.time = time.taken
      ))
    }
  }
  name.all = c(
    "naive",
    "uniform",
    "epanechnikov",
    "biweight",
    "triweight",
    "triangular",
    "compact_exp",
    "exponential",
    "cauchy"
  )
  # colnames(pred.test.y) = k.meth.name
  if (ALL) {
    if (is.null(y.test)) {
      if(BOTH){
        y.hat.test = list(pred.candidate = pred.test.y,
                          pred.cobra = pred.cob$prediction,
                          pred.mixcobra = pred.mix$prediction)
      }
      else{
        colnames(pred$prediction) = name.all[-2]
        y.hat.test = cbind(pred.test.y, pred$prediction)
      }
    }
    else{
      if(BOTH){
        name = c(k.meth.name, paste0(name.all[-2], ".cob"), paste0(name.all[-2], ".mix"))
        err = colMeans(sweep(as.matrix(pred.test.y), 1, y.test, FUN = "!="))
        cob.err = colMeans(sweep(pred.cob$prediction, 1, y.test, FUN = "!="))
        mix.err = colMeans(sweep(pred.mix$prediction, 1, y.test, FUN = "!="))
        y.hat.test = list(pred.candidate = pred.test.y,
                          pred.cobra = pred.cob$prediction,
                          pred.mixcobra = pred.mix$prediction)
        tab = data.frame(
          test.error = c(err, cob.err, mix.err),
          single.error = single.error,
          row.names = name
        )
      }
      else{
        name = c(k.meth.name, name.all[-2])
        err = colMeans(sweep(as.matrix(pred.test.y), 1, y.test, FUN = "!="))
        cob.err = colMeans(sweep(pred$prediction, 1, y.test, FUN = "!="))
        y.hat.test = cbind(pred.test.y, pred$prediction)
        colnames(y.hat.test) = name
        tab = data.frame(
          test.error = c(err, cob.err),
          single.error = single.error,
          row.names = name
        )
      }
    }
  }
  else{
    if (is.null(y.test)) {
      if(BOTH){
        colnames(pred.cob$prediction) = colnames(pred.mix$prediction) = C_step$kernels
        y.hat.test = list(pred.candidate = pred.test.y,
                          pred.cobra = pred.cob$prediction,
                          pred.mixcobra = pred.mix$prediction)
      }
      else{
        y.hat.test = cbind(pred.test.y, pred$prediction)
      }
    }
    else{
      if(BOTH){
        name = c(k.meth.name, paste0(C_step$kernels, ".cob"), paste0(C_step$kernels, ".mix"))
        colnames(pred.cob$prediction) = colnames(pred.mix$prediction) = C_step$kernels
        err = colMeans(sweep(as.matrix(pred.test.y), 1, y.test, FUN = "!="))
        cob.err = colMeans(sweep(pred.cob$prediction, 1, y.test, FUN = "!="))
        mix.err = colMeans(sweep(pred.mix$prediction, 1, y.test, FUN = "!="))
        y.hat.test = list(pred.candidate = pred.test.y,
                          pred.cobra = pred.cob$prediction,
                          pred.mixcobra = pred.mix$prediction)
        tab = data.frame(
          test.error = c(err, cob.err, mix.err),
          single.error = single.error,
          row.names = name
        )
      }
      else{
        err = colMeans(sweep(as.matrix(pred.test.y), 1, y.test, FUN = "!="))
        com.err = colMeans(sweep(pred$prediction, 1, y.test, FUN = "!="))
        name = c(k.meth.name, C_step$kernels)
        y.hat.test = cbind(pred.test.y, pred$prediction)
        colnames(y.hat.test) = name
        test.mse = c(err, com.err)
        tab = data.frame(
          test.error = test.mse,
          single.error = single.error,
          row.names = name
        )
      }
    }
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  if (is.null(y.test)) {
    return(list(train.predict = pred.train.y,
                test.predict = y.hat.test,
                run.time = time.taken))
  }
  else{
    return(list(
      table = tab,
      train.predict = pred.train.y,
      test.predict = y.hat.test,
      run.time = time.taken
    ))
  }
}


# ====================================================================================


dataSimulation = function(m, sigma, prop.test = 0.2, distrib, model, figure = F){
  start.time <- Sys.time()
  DIM3 = F
  LOG = F
  M = round((1 + prop.test) * m)
  sigmoid = function(x){
    return(1/(1+exp(-x)))
  }
  if((distrib == "Exp") || (distrib == "exp") ||(distrib == "exponential") || 
     (distrib == "Exponential")){
    x1 = c(rexp(M, 0.05) + 9, rexp(M, 0.5), rexp(M,0.1) + 8)
    x2 = c(rexp(M, 0.5), rexp(M, 0.05), rexp(M,0.1) + 8)
    x = cbind(x1, x2) + 1
  }
  if((distrib == "Pois") || (distrib == "pois") ||(distrib == "poisson") || 
     (distrib == "poisson")){
    x1 = c(rpois(M, 3), rpois(M,10), rpois(M,13) + 2)
    x2 = c(rpois(M, 11), rpois(M,2), rpois(M,12) + 2)
    x = cbind(x1, x2) + 1
  }
  if((distrib == "geom") || (distrib == "Geom") ||(distrib == "Geometric") || 
     (distrib == "geometric")){
    x1 = c(rgeom(M, 0.07)+9, rgeom(M, 0.55),rgeom(M, 0.15)+7)
    x2 = c(rgeom(M, 0.35), rgeom(M, 0.07)+17,rgeom(M, 0.15)+12)
    x = cbind(x1, x2) + 1
  }
  if((distrib == "gaussian2") || (distrib == "Gaussian2") ||(distrib == "Normal2") || 
     (distrib == "normal2")){
    x1 = c(rnorm(M, 4, 1), rnorm(M, 22, 2),rnorm(M, 10, 2))
    x2 = c(rnorm(M, 12, 1), rnorm(M, 9, 1),rnorm(M, 5, 2))
    x = cbind(x1, x2)
    x = x + abs(min(x)) + 1
  }
  if((distrib == "gaussian3") || (distrib == "Gaussian3") ||(distrib == "Normal3") || 
     (distrib == "normal3")){
    x1 = c(rnorm(M, 6, 1), rnorm(M, 5, 2),rnorm(M, 8, 1))
    x2 = c(rnorm(M, 14, 2), rnorm(M, 10, 1),rnorm(M, 6, 1))
    x3 = c(rnorm(M, 6, 1), rnorm(M, 15, 2),rnorm(M, 14, 2))
    x = cbind(x1,x2,x3) 
    x = x + abs(min(x)) + 5
    DIM3 = T
  }
  if((model == "logit") || (model == "logistic") || (model == "Logit") || 
     (model == "Logistic")){
    LOG = T
  }
  id.train = c(1:m,(M+1):(M+m),(2*M+1):(2*M+m))
  x.train = x[id.train,]
  x.test = x[-id.train,]
  df = as.matrix(cbind(x0 = 1, x))
  p = ncol(x)
  if(DIM3){
    b1 = c(-10,3,7)
    b2 = c(7,5,-12)              # We fixed the coeffiencets of the 3D case
    b3 = c(6,-11,10)
    beta = c(b1,b2,b3)
  }
  else{
    b1 = c(-8,3)
    b2 = c(-6,-5)                 # We fixed the coeffiencets of the 2D case
    b3 = c(5,-7)
    beta = c(b1,b2,b3)
  }
  
  if(LOG){
    beta1 = c(-colMeans(x.train[1:m,]) %*% beta[1:p], beta[1:p])
    beta2 = c(-colMeans(x.train[(m+1):(2*m),]) %*% beta[(p+1):(2*p)], 
              beta[(p+1):(2*p)])
    beta3 = c(-colMeans(x.train[(2*m+1):(3*m),]) %*% beta[(2*p+1):(3*p)], 
              beta[(2*p+1):(3*p)])
    y.1 = df[1:M,] %*% beta1
    y.2 = df[(M+1):(2*M),] %*% beta2 
    y.3 = df[(2*M+1):(3*M),] %*% beta3
    z = sigmoid(c(y.1,y.2,y.3)+ rnorm(3*M,0,sigma))
    target = factor(as.numeric(z > 0.5), levels = c(1,0))
  }
  else{
    beta1 = c(-15,beta[1:p])
    beta2 = c(25,beta[(p+1):(2*p)])
    beta3 = c(-10,beta[(2*p+1):(3*p)])
    y.1 = df[1:M,] %*% beta1
    y.2 = df[(M+1):(2*M),] %*% beta2 
    y.3 = df[(2*M+1):(3*M),] %*% beta3
    target = c(y.1, y.2,y.3) + rnorm(3*M,0,sigma)
  }
  cluster = c(rep(1,M), rep(2,M), rep(3,M))
  if(LOG){
    y.train = target[id.train]
    y.test = target[-id.train]
    col.train = as.numeric(y.train) + 2
  }
  else{
    y.train = target[id.train]
    y.test = target[-id.train]
    col.train = c(rep(2,m), rep(3,m), rep(4,m))
  }
  if(figure){
    if(DIM3){
      par(mfrow=c(1,1))
      scatterplot3d::scatterplot3d(x.train, color = cluster[id.train] + 1, 
                                   main = "Simulated Predictor", xlab = "x1", 
                                   ylab = "x2", zlab = "x3")
    }
    else{
      par(mfrow=c(1,2))
      if(LOG){
        plot(x.train, col=cluster[id.train] + 1, main = "Predictors", xlab = "x1", ylab = "x2")
        scatterplot3d::scatterplot3d(x=x.train[,1], y=x.train[,2],z=y.train, color = as.numeric(y.train) + 5, main = "Completed Data", xlab = "x1", ylab = "x2", zlab = "y")
      }
      else{
        plot(x.train[1:m,], col=2, main = "Predictors", xlim = c(min(x.train[,1])-5, max(x.train[,1])+5), ylim = c(min(x.train[,2])-5, max(x.train[,2])+5), pch = 1, cex.main = 2)
        points(x.train[(m+1):(2*m),], col=3, pch=17)
        points(x.train[(2*m+1):(3*m),], col=4, pch=8)
        scatterplot3d::scatterplot3d(x=c(x.train[,1],x.train[,1]), y=c(x.train[,2],x.train[,2]),z=c(y.train + max(y.train), rep(min(y.train),length(y.train))), color = rep(cluster[id.train] + 1,2), main = "Complete Data", xlab = "x1", ylab = "x2", zlab = "y",  cex.main = 2, angle = 20)
      }
    }
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  return(list(x.train = as.data.frame(x.train), 
              y.train = y.train,
              x.test = as.data.frame(x.test),
              y.test = y.test,
              BETA = cbind(beta1, beta2, beta3), 
              CLUSTERS = cluster, Running.Time = time.taken))
}

# ============================ Application ==================================



  df = dataSimulation(500, 5, 0.3, distrib = "normal2", model = "logit" ,figure = TRUE)
  dftrain = df$x.train
  ytrain = as.integer(df$y.train) - 1
  dftest = df$x.test
  ytest = as.integer(df$y.test) - 1
  
  kfc = KFCclass(train.input = dftrain, 
                 train.responses = ytrain, 
                 test.input = dftest,
                 test.responses = ytest, 
                 split_data = 0.5,
                 scale.input = NULL,
                 K_step = Kstep(K = 3, max.iter = 500, method = c("euclidean", "gkl" , "itakura", "logistic"), figure = TRUE),
                 F_step = Fstep(single.model = TRUE),
                 C_step = Cstep(method = "mixed", kernels = c("expo"), min.h = 0.0001, full.naive = FALSE,
                                max.h = 5, n.h = 100, figure = TRUE, show.info = TRUE, n.cv = 10,
                                min.alpha = 0.01, min.beta = 0.01, max.alpha = 5, max.beta = 5,
                                n.alpha = 20, n.beta = 20))

library(randomForest)
rf = randomForest(x = as.data.frame(dftrain), y = as.factor(ytrain), ntree = 100)
pred1 = predict(rf, as.data.frame(dftest), type = "class")
library(adabag)
dftrain$y = as.factor(ytrain)
boost <- boosting(y ~., data = as.data.frame(dftrain), mfinal = 100)
pred2 <- predict(boost, as.data.frame(dftest), type = "class")

dftest$y = as.factor(ytest)
err = c(comb$test.error, mean(pred1 != Y[-train]), mean(pred2$class != Y[-train]))
names(err) = c(names(comb$test.error), "randomForest", "boost")
err


