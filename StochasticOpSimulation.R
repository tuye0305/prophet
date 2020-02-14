# We are trying to solve the following optimization problem
# max f(x) = E(F(x,U)) s.t. E(-g(x,Vk)) <= 0, 0.999 =< sum(x in each cohort) <= 1.001
# X is not stochastic here.
# U,Vk are random variable.

rm(list = ls())      # Clear all variables
graphics.off()    # Close graphics windows

# Loading required packages
library(rosqp)
library(MASS)
library(Matrix)
library(data.table)
library(jsonlite)
library(ggplot2)
library(foreach)
library(doParallel)

########################### Helper Functions ##########################
evaluateConstraints <- function(VList, xList, etaList, k, d, J, N){
  # evaluate k constriants E(-g(x,Vk)) <= etak and return the violated ones
  #
  # Args:
  #   VList: a list of values for the constraints.
  #   xList: a list represents the probabilistic assignments of treatments.
  #   etaList: a list of eta represent the thresholds of constraints.
  #   k: the number of constriants.
  #   d: dimension of x (number of treatment * number of cohorts)
  #   J: the number of samples to be evaluated
  #   N: the number of total iterations
  #
  # Returns:
  #   A list of index and a list of values for the constraints that got violated.
  VListSingle <- lapply(VList, function(x) { x[k,]})
  etaListSingle <- lapply(etaList, function(x) { x[k,]})
  
  randIdxs <- sample(1 : N, J, replace = TRUE)
  violateList <- list()
  idxList <- list()
  cnt <- 1
  sumV <- rep(0) * d
  for (i in 1 : length(etaList)) {
    muV <- colMeans(VList[[i]][randIdxs,])
    if (- muV %*% xList[[k]] > etaListSingle[[i]]) {
      violateList[[cnt]] <- as.matrix(VListSingle[[i]])
      idxList[[cnt]] <- i
      cnt <- cnt + 1
    }
  }
  list(idxList, violateList)
}

addCohortLevelConstraints <- function(A.matrix, l, m){
  # per each cohort add a column in matrix A
  #
  # Args:
  #   A.metrix: a matrix represent the constraints
  #   l: the number of treatments
  #   m: the numberof cohorts
  # Returns:
  #   a revised A matrix.
  baseCol <- rep(0, l * m)
  for (i in 0 : (m - 1)) {
    col <- baseCol
    col[c((l * i + 1) : (l * (i + 1)))] <- 1
    A.matrix <- rbind(A.matrix, t(col))
  }
  A.matrix
}


csa.solve <- function(muU, sigmaU, muVList, sigmaVList, N, l, m, gamma, etaList, alpha, s, J, type = "SGD") {
  # The code below solves Max x^tU s.t -x^Vk <= 0, 0 <= x <= 1, sum(x per cohort) = 1
  #
  # Args:
  #   muU: the mean of objective utility for each treatment and cohort.
  #   sigmaU: the variance of objective utility for each treatment and cohort.
  #   muVList: a list of the mean of constraint utilities for each treatment and cohort.
  #   sigmaVList: a list of the variance of constraint utilities for each treatment and cohort.
  #   N: number of iterations
  #   l: dimension of treatments
  #   m: dimension of cohorts
  #   gamma: initializations of the optimization.
  #   etaList: a list of eta represent the thresholds of constraints.
  #   alpha: hyperparameter of the learning rate.
  #   s: the starting iteration for evaluate whether the solution satisfy the constraints.
  #   J: the number of samples to be evaluated.
  #   type: the choice of optimization method "SGD" or "Adagrad" (default method is "SGD")
  # Returns:
  #   Results of the optimization in terms of
  # 1)final treatment assignment;
  # 2) treatment assignments for each iteration;
  # 3) the number of solutions that satisfy all the constraints.
  U <- mvrnorm(n = N, mu = muU, Sigma = sigmaU)
  VList <- list()
  for (i in 1 : length(sigmaVList)) {
    VList[[i]] <- mvrnorm(n = N, mu = muVList[[i]], Sigma = sigmaVList[[i]])
  }
  
  d = length(muU)
  xList <- list()
  xList[[1]] <- rep(0, d)
  
  if (type == "SGD") {
    # First we try regular Hk = I/alpha
    Hk <- 2 * sparseMatrix(i = 1 : d, j = 1 : d, x = rep(1, d)) / alpha
  }
  
  HkMat <- mat.or.vec(d, d)
  direction <- "0"
  
  for (k in 1 : N) {
    violations <- evaluateConstraints(VList, xList, etaList, k, d, J, N)
    violatedIdx <- violations[[1]]
    violatedList <- violations[[2]]
    if (length(violatedList) == 0) {
      gk <- - U[k,]
      direction <- "0"
    } else {
      # randomly choose a violated constriant as the gradient (gk)
      randIdx <- unlist(sample(1 : length(violatedList), 1, replace = FALSE))
      gk <- - as.matrix(violatedList[[randIdx]])
      direction <- as.character(violatedIdx[[randIdx]])
    }
    if (type == "Adagrad") {
      HkMat <- HkMat + gk %*% t(gk)
      Hk <- 2 * sparseMatrix(i = 1 : d, j = 1 : d, x = Re(sqrt(diag(HkMat))) / alpha)
    }
    results <- solve_osqp(P = Hk,
                          q = (gamma[k] * gk - Hk %*% xList[[k]]),
                          A = addCohortLevelConstraints(sparseMatrix(i = 1 : d, j = 1 : d, x = rep(1, d)), l, m),
                          l = c(rep(0, d), rep(0.999, m)),
                          u = c(rep(1, d), rep(1.001, m)),
                          osqpSettings(eps_abs = 1e-6, eps_rel = 1e-6, verbose = TRUE))
    xList[[k + 1]] <- results$x
  }
  
  numerator <- 0
  denominator <- 0
  count <- 0
  for (k in s : N) {
    violatedList <- evaluateConstraints(VList, xList, etaList, k, d, J, N)[[2]]
    if (length(violatedList) == 0) {
      numerator <- numerator + gamma[k] * xList[[k]]
      denominator <- denominator + gamma[k]
      count <- count + 1;
    }
  }
  
  xFinal <- numerator / denominator
  return(list(xFinal, xList, count))
}

csa.solveFull <- function(muU, sigmaU, muVList, sigmaVList, N, l, m, tau, alpha, J){
  # solve for the parameter assignments (Max x^tU s.t x^Vk >= 0, 0 <= x <= 1, sum(x) = 1)
  #
  # Args:
  #   muU: the mean of objective utility for each treatment and cohort.
  #   sigmaU: the variance of objective utility for each treatment and cohort.
  #   muVList: a list of the mean of constraint utilities for each treatment and cohort.
  #   sigmaVList: a list of the variance of constraint utilities for each treatment and cohort.
  #   N: number of iterations.
  #   l: dimension of treatments.
  #   m: dimension of cohorts.
  #   tau: hyperparameter for constraints violation level.
  #   alpha: hyperparameter of the learning rate.
  #   J: the number of samples to be evaluated.
  # Returns:
  #   the optimization solutions with SGD and Adagrad options.
  alpha <- 1
  s <- N / 2
  d <- l * m
  
  Dx <- sqrt(d) / alpha
  MF <- sqrt(sum(diag(sigmaU)) + sum(muU ^ 2))
  MGList <- list()
  for (i in 1 : length(muVList)) {
    MGList[[i]] <- sqrt(sum(diag(sigmaVList[[i]])) + sum(muVList[[i]] ^ 2))
  }
  
  gamma <- Dx / ((MF + sum(unlist(MGList))) * sqrt(1 : N))
  etaList <- list()
  for (i in 1 : length(MGList)) {
    etaList[[i]] <- as.matrix(tau * Dx * (MF + MGList[[i]]) / sqrt(1 : N))
  }
  
  csaSGD <- csa.solve(muU, sigmaU, muVList, sigmaVList, N, l, m, gamma, etaList, alpha, s, J, "SGD")
  csaAdagrad <- csa.solve(muU, sigmaU, muVList, sigmaVList, N, l, m, gamma, etaList, alpha, s, J, "Adagrad")
  list(csaSGD, csaAdagrad)
}

trim <- function (x){
  # trim a string to remove white spaces
  # Args: x: the raw string.
  #
  # Returns:
  #   a trimmed string.
  gsub("^\\s+|\\s+$", "", x)
}

processUtility <- function(colIdx, data, idx){
  # process the utility input data (mean pm error margin) and extract the mean & error margin information
  # Args:
  #   colIdx: column index of the utility in the dataframe.
  #   data: the dataframe with all utility information.
  #   idx: the index of the item within the column.
  # Returns:
  #   A matrix form of the processed utilities.
  utility <- c()
  for (i in 1 : length(data[, colIdx])) {
    utility[[i]] <- as.numeric(trim(strsplit(data[, colIdx][[i]], ";")[[1]][idx]))
  }
  as.matrix(utility)
}

writeDetails <- function(outputPath, result, muU, muVList, tau, alpha, metricDesiredDirectionsList, fileName) {
  # write details of the QP solutions in json format
  # Args:
  #   outputPath: the path of output directory.
  #   result: the result of solved treatment assignments from the optimization.
  #   muU: the mean of objective utility for each treatment and cohort.
  #   muVList: a list of the mean of constraint utilities for each treatment and cohort.
  #   tau: hyperparameter for constraints violation level.
  #   alpha: hyperparameter of the learning rate.
  #   metricDesiredDirectionsList: a list of the desired direction of metrics of interets
  #   fileName: the file name for the output.
  # Returns:
  # NULL
  detailResults <- list()
  detailResults[[1]] <- paste0("Objective Value: ", result[[1]] %*% muU)
  cnt <- 2
  for (i in 1 : length(muVList)) {
    metricDirection <- metricDesiredDirectionsList[[1]][i + 1]
    if (metricDirection == "up")
      detailResults[[cnt]] <- paste0("Constraint: ", result[[1]] %*% muVList[[i]])
    else detailResults[[cnt]] <- paste0("Constraint: ", - result[[1]] %*% muVList[[i]])
    cnt <- cnt + 1
  }
  detailResults[[cnt]] <- paste0("Chosen Output: ")
  detailResults[[cnt + 1]] <- result[[1]]
  detailResults[[cnt + 2]] <- paste0(" count:", result[[3]])
  jsonResults <- toJSON(detailResults, auto_unbox = TRUE)
  write(prettify(jsonResults), file = paste0(outputPath, '/', fileName, tau, "-", alpha))
}

plotResults <- function(output, constraintDirection, fxSGDMean, fxAdagradMean, fxSGDStd, fxAdagradStd, fxSGDConstraintMean, fxAdagradConstraintMean, fxSGDConstraintStd, fxAdagradConstraintStd){
  # plotting the Expectation of effect on objective metric along with the increasing iterations
  # Args:
  #   output: the path of the output.
  #   constraintDirection: the desired direction of the constraint.
  #   fxSGDMean: the mean of F(x,U) on the objective utility given the SGD solution.
  #   fxAdagradMean: the mean of F(x,U) on the objective utility given the Adagrad solution.
  #   fxSGDStd: the standard deviation of F(x,U) on the objective utility given the SGD solution.
  #   fxAdagradStd: the standard deviation of F(x,U) on the objective utility  given the Adagrad solution.
  #   fxSGDConstraintMean: the mean of F(x,V) on the constriant utility given the SGD solution.
  #   fxAdagradConstraintMean: the mean of F(x,V) on the constriant utility given the Adagrad solution.
  #   fxSGDConstraintStd: the standard deviation of F(x,V) on the constriant utility given the SGD solution.
  #   fxAdagradConstraintStd: the standard deviation of F(x,V) on the constriant utility given the Adagrad solution.
  # Returns:
  #   NULL
  x <- c(1 : length(fxSGDMean))
  if (constraintDirection == "down") {
    fxAdagradConstraintMean = - fxAdagradConstraintMean
  }
  
  plt <- ggplot() +
    geom_line(aes(x, fxSGDMean, colour = "SGD Objective"), size = 1.5) +
    geom_line(aes(x, fxAdagradMean, colour = "Adagrad Objective"), size = 1.5) +
    geom_ribbon(aes(x, ymin = fxSGDMean - fxSGDStd, ymax = fxSGDMean + fxSGDStd), fill = "#C77CFF", alpha = 0.3) +
    geom_ribbon(aes(x, ymin = fxAdagradMean - fxAdagradStd, ymax = fxAdagradMean + fxAdagradStd), fill = "#7CAE00", alpha = 0.3) +
    geom_line(aes(x, fxSGDConstraintMean, colour = "SGD Constraint"), size = 1.5) +
    geom_line(aes(x, fxAdagradConstraintMean, colour = "Adagrad Constraint"), size = 1.5) +
    geom_ribbon(aes(x, ymin = fxSGDConstraintMean - fxSGDConstraintStd, ymax = fxSGDConstraintMean + fxSGDConstraintStd), fill = "#00BFC4", alpha = 0.3) +
    geom_ribbon(aes(x, ymin = fxAdagradConstraintMean - fxAdagradConstraintStd, ymax = fxAdagradConstraintMean + fxAdagradConstraintStd), fill = "#F8766D", alpha = 0.3) +
    xlab("Iterations") +
    ylab("Estimations of Effects") +
    theme(text = element_text(size = 40),
          axis.text.x = element_text(size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size = 20),
          legend.position = c(0.7, 0.4),
          legend.key.size = unit(4, 'lines'),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  ggsave(output, plot = plt, width = 9, height = 7.5)
}


randPercent <- function(N, M) {
  # generate a list of N percentages that will sum up to a fixed values M
  # Args:
  #   N: the number of values in the output
  #   M: the number of the sum of the output values
  # Returns:
  #   A list of numeric values
  vec <- abs(rnorm(N))
  vec / sum(vec) * M
}

summarizeResults <- function(resultList) {
  # process the results (list format) and calculate the mean and standard deviations and output the summary stats
  # Args:
  #   resultList: A list of results
  #
  # Returns:
  #   A list with the mean and standard deviation of the results
  colSize = length(resultList)
  rowSize = length(resultList[[1]])
  resultMatrix <- matrix(unlist(resultList), ncol = colSize, nrow = rowSize)
  resultMean <- rowMeans(resultMatrix)
  resultStd <- unlist(apply(resultMatrix, 1, function(x) { sd(x)}))
  list(resultMean, resultStd)
}

########################### Preparing Simulation Input Data ##########################
l <- 10 # dimension of treatments
N <- 200 #iters
J <- 50 #number of samples used for evaluate the constraints
m <- 3 # dimension of cohorts

set.seed(5)
cohorts <- rep(unlist(lapply(1 : m, function(x) {paste("Cohort: ", x)})), l)
treatments <- rep(unlist(lapply(1 : l, function(x) {paste("Treatment: ", x)})), m)

#Assume the input data has the mean, standard deviation of the effects, followed by the sample size % for each cohort
metric.1 <- unlist(mapply(function(x, y, z) {paste(x, ";", y, ";", z)},
                          rnorm(l * m), abs(rnorm(l * m)), rep(randPercent(m , 1), l)))
metric.2 <- unlist(mapply(function(x, y, z) {paste(x, ";", y, ";", z)},
                          rnorm(l * m), abs(rnorm(l * m)), rep(randPercent(m , 1), l)))
metric.3 <- unlist(mapply(function(x, y, z) {paste(x, ";", y, ";", z)},
                          rnorm(l * m), abs(rnorm(l * m)), rep(randPercent(m , 1), l)))
utilityData <- data.frame(
  "cohort" = cohorts,
  "treatment" = treatments,
  "metric.1" = metric.1,
  "metric.2" = metric.2,
  "metric.3" = metric.3
)
utilityData$metric.1 <- as.character(utilityData$metric.1)
utilityData$metric.2 <- as.character(utilityData$metric.2)
utilityData$metric.3 <- as.character(utilityData$metric.3)

utilityData <- utilityData[order(utilityData$cohort, utilityData$treatment),]
d <- nrow(utilityData)

metricDesiredDirections <- "up-up-down" #preferred metric directions for metric 1, 2, 3

########################### Running the Optimization per each cohort ##########################
outputPath <- Sys.getenv("HOME")
constraintIdx <- 1 # one of the important constraints that will be visualized in the result plot
constraintDirection <- "up"
numCores <- detectCores()
# Extract the utility data (both point estimates and variances)
metricDesiredDirectionsList <- strsplit(metricDesiredDirections, "-")
muU <- processUtility(3, utilityData, 1) * processUtility(3, utilityData, 3)
sigmaUv <- processUtility(3, utilityData, 2) * processUtility(3, utilityData, 3) ^ 2
sigmaU <- diag(as.list(sigmaUv))

muVList <- list()
sigmaVList <- list()

for (i in 4 : ncol(utilityData)) {
  cnt <- i - 3
  metricDirection <- metricDesiredDirectionsList[[1]][cnt + 1]
  if (metricDirection == "up")
    muVList[[cnt]] <- processUtility(i, utilityData, 1) * processUtility(i, utilityData, 3)
  else if (metricDirection == "down")
    muVList[[cnt]] <- - processUtility(i, utilityData, 1) * processUtility(i, utilityData, 3)
  sigmaV <- processUtility(i, utilityData, 2) * processUtility(i, utilityData, 3) ^ 2
  sigmaVList[[cnt]] <- diag(as.list(sigmaV))
}

taus <- c(0.001) # hyperparameter for constraints violation level
alphas <- c(0.01) # hyperparameter of the learning rate
for (tau in taus) {
  for (alpha in alphas) {
    results <- list()
    registerDoParallel(numCores)
    # Parallelly run the stochastic optimizations
    results <- foreach (j = 1 : 10, .combine = "c") %dopar% {
      csa.solveFull(muU, sigmaU, muVList, sigmaVList, N, l, m, tau, alpha, J)
    }
    #stop parallel runs
    registerDoSEQ()
    
    csaSGDList <- list()
    csaAdagradList <- list()
    fxSGDList <- list()
    fxAdagradList <- list()
    fxSGDConstraintList <- list()
    fxAdagradConstraintList <- list()
    cnt <- 1
    for (i in seq(1, 20, by = 2)) {
      csaSGD <- results[[i]]
      csaAdagrad <- results[[i + 1]]
      csaSGDList[[cnt]] <- csaSGD
      csaAdagradList[[cnt]] <- csaAdagrad
      
      fxSGDList[[cnt]] <- unlist(lapply(csaSGD[[2]], function(x) { x %*% muU}))
      fxSGDConstraintList[[cnt]] <- unlist(lapply(csaSGD[[2]], function(x) { x %*% muVList[[constraintIdx]]}))
      fxAdagradList[[cnt]] <- unlist(lapply(csaAdagrad[[2]], function(x) { x %*% muU}))
      fxAdagradConstraintList[[cnt]] <- unlist(lapply(csaAdagrad[[2]], function(x) { x %*% muVList[[constraintIdx]]}))
      cnt <- cnt + 1
    }
    #extract summary stats for the results
    fxSGDMean <- summarizeResults(fxSGDList)[[1]]
    fxSGDStd <- summarizeResults(fxSGDList)[[2]]
    fxAdagradMean <- summarizeResults(fxAdagradList)[[1]]
    fxAdagradStd <- summarizeResults(fxAdagradList)[[2]]
    
    fxSGDConstraintMean <- summarizeResults(fxSGDConstraintList)[[1]]
    fxSGDConstraintStd <- summarizeResults(fxSGDConstraintList)[[2]]
    fxAdagradConstraintMean <- summarizeResults(fxAdagradConstraintList)[[1]]
    fxAdagradConstraintStd <- summarizeResults(fxAdagradConstraintList)[[2]]
    
    #save all the outputs into files
    resultDF <- data.frame("cohort" = utilityData$cohort, "treatment" = utilityData$treatment, "SGD" = csaSGD[[1]], "Adagrad" = csaAdagrad[[1]])
    writeDetails(outputPath, csaSGD, muU, muVList, tau, alpha, metricDesiredDirectionsList, "SGDDetails")
    writeDetails(outputPath, csaAdagrad, muU, muVList, tau, alpha, metricDesiredDirectionsList, "AdagradDetails")
    write.csv(resultDF, file = paste0(outputPath, '/FinalAssignments', tau, "-", alpha, '.csv'), row.names = FALSE, quote = FALSE)
    
    output <- paste0(outputPath, "/ObjVsOneConstraintPlot", tau, "-", alpha, ".pdf")
    plotResults(output, constraintDirection, fxSGDMean, fxAdagradMean, fxSGDStd, fxAdagradStd, fxSGDConstraintMean, fxAdagradConstraintMean, fxSGDConstraintStd, fxAdagradConstraintStd)
  }
}