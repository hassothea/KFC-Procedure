# 1. Install all the necessary libraries.
# 2. Run all the codes in the "KFC_regression.R" file
# 3. Run the following codes.

# Example for regression problem
# ==============================
df.norm2.reg = simulateData(1000, 10, 0.3, "normal2", "reg")


data.train = as.matrix(df.norm2.reg$Training.Data[,1:2])
y.train = df.norm2.reg$Training.Data[,3]
data.test = as.matrix(df.norm2.reg$Testing.Data[,1:2])
y.test = df.norm2.reg$Testing.Data[,3]

kfc = KFCreg(train.input = data.train, 
             train.responses = y.train, 
             test.input = data.test,
             test.responses = y.test,
             split_data = 0.5,
             scale.input.cluster = NULL,
             scale.input.combine = "standard",
             scale.machines = 'standard',
             K_step = Kstep(K = 3, max.iter = 500, method = c("euclidean", "gkl" , "itakura", "logistic"), deg = c(3,4), figure = TRUE),
             F_step = Fstep(single.model = FALSE),
             C_step = Cstep(method = "mixed", kernels = c("gaussian"), 
                            opt.methods = c("GD"),
                            max.eps = 1, n.eps = 300, figure = TRUE, show.info = TRUE,
                            min.alpha = 0.05, min.beta = 0.05, min.eps = 0.001, max.alpha = 1,
                            max.beta = 1, n.alpha = 30, n.beta = 30,
                            n.cv = 10),
             set_GD = setGD(print.step = TRUE, rate = "log", coef.auto = 0.005, figure = TRUE, max.iter = 300, threshold = 1e-15))


library(randomForest)
rf = randomForest(data.train, y = y.train, xtest = data.test, ntree = 500)

library(gbm)
boost = gbm(y.train~., data = as.data.frame(data.train), n.trees = 500, distribution = "gaussian")
pred = predict.gbm(boost, as.data.frame(data.test))

# MSE
mse = c(kfc$table$mse.test, mean((rf$test$predicted - y.test)^2), mean((pred - y.test)^2))
names(mse) = c(rownames(kfc$table), "randomForest", "boost")
mse


# Example for classification problem.
# ===================================


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

