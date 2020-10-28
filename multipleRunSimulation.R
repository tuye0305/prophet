rm(list = ls())      # Clear all variables
graphics.off()    # Close graphics windows
gc()
# Loading required packages
library(simcausal)
library(causalTree)
library(rosqp)
library(MASS)
library(Matrix)
library(data.table)
library(jsonlite)
library(ggplot2)
library(foreach)
library(doParallel)
library(rattle)
library(lpSolve)
library(osqp)
library(olpsR)

checkResultData <- function(data){
  (is.data.frame(data) && nrow(data)!=0) || (!is.data.frame(data) && !is.na(data))
}

homePath <- "/Users/yetu/prophet/" #TODO: Please specify the home directory for the project
outcomeList <- c("Y1", "Y2", "Y3")
treatmentList <- c("A1", "A2", "A3")

source(paste0(homePath, "generateSimulationDataUtils.R"))
source(paste0(homePath, "memberLevelProphetUtils.R"))
source(paste0(homePath, "cohortLevelProphetUtils.R"))
source(paste0(homePath, "mergeTreeUtils.R"))

seed <- 12345

set.seed(seed)
size <- 500
uncertaintyW <- 1
w1List <- c(rnorm(1, mean = 1, sd = 1), rnorm(1, mean = 1, sd = 1), rnorm(1, mean = 1, sd = 1))
w2List <- c(rnorm(1, mean = 1, sd = 1), rnorm(1, mean = 1, sd = 1), rnorm(1, mean = 1, sd = 1))
w3List <-  c(rnorm(1, mean = 0, sd = 1), rnorm(1, mean = 0, sd = 1), rnorm(1, mean = 0, sd = 1))
w4List <-  c(rnorm(1, mean = 0, sd = 1), rnorm(1, mean = 0, sd = 1), rnorm(1, mean = 0, sd = 1))
w5List <-  c(rnorm(1, mean = -0.5, sd = 1), rnorm(1, mean = -0.5, sd = 1), rnorm(1, mean = -0.5, sd = 1))
w6List <-   c(rnorm(1, mean = -0.5, sd = 1), rnorm(1, mean = -0.5, sd = 1), rnorm(1, mean = -0.5, sd = 1))

minLeafSize<- 150
numTreatment <- length(treatmentList)
objective.outcome <- "Y1"

uncertaintyWList <- c(5,10,15,20)
numEpoc <- 10

for (uncertaintyW in uncertaintyWList){
  tau <- if (uncertaintyW <= 1) 0.005 else 0.001/(uncertaintyW^2) # hyperparameter for constraints violation level
  alpha <- 0.01 # hyperparameter of the learning rate
    
  causaltreeMergedStochasticPolicyResultDF <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
  causaltreeStochasticPolicyResultDF <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
  causalTreeDeterministicPolicyResultDF <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
  causalForestDeterministicPolicyResultDF <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
  deltaModelPolicyResultDF <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
  fixedParameterPolicyResultDF <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
  heuristicStochasticPolicyResultDF <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
  for (e in 1:numEpoc){
    # one epoc
    seedList <- lapply( 1:numTreatment, function(x) c(x * e * 5, x * e * 1999))
    generateSimulationData(homePath, treatmentList, w1List, w2List, w3List, w4List, w5List, w6List, seedList, uncertaintyW)
    trainCausalTree(homePath, outcomeList, treatmentList, minLeafSize)
    
    #merging trees
    Y1.ATE <- -Inf
    causaltreeMergedStochasticPolicyResults <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
    
    dt.list <-list()
    cnt <- 1
    for (treatment in treatmentList){
      for (outcome in outcomeList){
        decisionTable <- fread(paste0(homePath, treatment, outcome, 'decisionTable.csv') , header = TRUE)
        dt.list[[cnt]] <- decisionTable
        cnt <- cnt + 1
      }
    }
    #impute data
    dt.out.list <- imputeCols(dt.list)
    dt.table.final <- mergeTrees(dt.out.list)
    decisionTable<- dt.table.final
    write.csv(dt.table.final, file = paste(homePath, 'mergedDecisionTable.csv', sep = ""), row.names = FALSE, quote = FALSE)
    
    utilityDF <- calculateUtilitiesMergedTree(dt.table.final, homePath, outcomeList[[1]], treatmentList)
    utilityDF <- subset(utilityDF, select = c("cohort", "treatment"))
    for (outcome in outcomeList){
      utilities <- calculateUtilitiesMergedTree(dt.table.final, homePath, outcome, treatmentList)
      utilityDF[[outcome]] <- as.character(utilities[["utilities"]])  
    }
    utilityDF <- utilityDF[order(utilityDF$cohort),]
    write.csv(utilityDF, file = paste(homePath, 'utilityDF.csv', sep = ""), row.names = FALSE, quote = FALSE)
    
    stochasticOptimization(homePath, utilityDF, numCohorts = nrow(decisionTable), numTreatments = length(treatmentList), "CT.M", c(tau), c(alpha))
    finalAssignment <- readData(paste0(homePath, '/FinalAssignments',"CT.M", tau, "-", alpha, '.csv'))
    policy.adagrad <- scoreTestDataGeneratePolicyMerged(decisionTable, homePath, treatmentList,finalAssignment, method = "Adagrad")
    policy.sgd<- scoreTestDataGeneratePolicyMerged(decisionTable, homePath, treatmentList,finalAssignment, method = "SGD")
    
    newPolicyResult.ada <- evaluatePolicy(policy.adagrad, homePath, treatmentList, outcomeList)
    newPolicyResult.sgd <- evaluatePolicy(policy.sgd, homePath, treatmentList, outcomeList)
    newPolicyResult <- if(is.na(newPolicyResult.sgd[1])) newPolicyResult.ada else newPolicyResult.sgd
    causaltreeMergedStochasticPolicyResults <- newPolicyResult
    
    #member level causal tree 
    generateCausalTreeEstimations(homePath, outcomeList, treatmentList, minLeafSize)
    causalTreeDeterministicpolicy <- generatePolicy(homePath, 'causalTreeEstimationsTraining.csv', 'causalTreeEstimationsTesting.csv', treatmentList, outcomeList, objective.outcome, 0)
    causalTreeDeterministicpolicyResult <- evaluatePolicy(causalTreeDeterministicpolicy, homePath, treatmentList, outcomeList)
    
    #member level causal forest 
    generateCausalForestEstimations(homePath, outcomeList, treatmentList, minLeafSize)
    causalForestDeterministicpolicy <- generatePolicy(homePath, 'causalForestEstimationsTraining.csv', 'causalForestEstimationsTesting.csv', treatmentList, outcomeList, objective.outcome, 0)
    causalForestDeterministicpolicyResult <- evaluatePolicy(causalForestDeterministicpolicy, homePath, treatmentList, outcomeList)
    
    #member level delta model
    generateTwoDeltaModelEstimations(homePath, outcomeList, treatmentList, minLeafSize)
    deltaModelpolicy <- generatePolicy(homePath, 'deltaModelEstimationsTraining.csv', 'deltaModelEstimationsTesting.csv', treatmentList, outcomeList, objective.outcome, 0)
    deltaModelpolicyResult <- evaluatePolicy(deltaModelpolicy, homePath, treatmentList, outcomeList)
    
    #heuristic cohorts
    heuristic.cohorts <- data.frame(createHeuristCohorts(c("H1", "H2", "H3", "H4"), homePath, treatmentList))
    utilityDF <- calculateUtilities(heuristic.cohorts, homePath, outcomeList[[1]], treatmentList)
    utilityDF <- subset(utilityDF, select = c("cohort", "treatment"))
    for (outcome in outcomeList){
      utilities <- calculateUtilities(heuristic.cohorts, homePath, outcome, treatmentList)
      utilityDF[[outcome]] <- as.character(utilities[["utilities"]])  
    }
    utilityDF <- utilityDF[order(utilityDF$cohort),]
    write.csv(utilityDF, file = paste(homePath, 'heuristicutilityDF.csv', sep = ""), row.names = FALSE, quote = FALSE)
    
    stochasticOptimization(homePath, utilityDF, numCohorts = nrow(heuristic.cohorts), numTreatments = length(treatmentList), "heuristic", c(tau), c(alpha))
    finalAssignment <- readData(paste0(homePath, '/FinalAssignments',"heuristic", tau, "-", alpha, '.csv'))
    policy.adagrad <- scoreTestDataGeneratePolicyMerged(decisionTable, homePath, treatmentList,finalAssignment, method = "Adagrad")
    policy.sgd<- scoreTestDataGeneratePolicyMerged(decisionTable, homePath, treatmentList,finalAssignment, method = "SGD")
    
    newPolicyResult.ada <- evaluatePolicy(policy.adagrad, homePath, treatmentList, outcomeList)
    newPolicyResult.sgd <- evaluatePolicy(policy.sgd, homePath, treatmentList, outcomeList)
    newPolicyResult <- if(is.na(newPolicyResult.sgd[1])) newPolicyResult.ada else newPolicyResult.sgd
    heuristicStochasticPolicyResults <- newPolicyResult
    
    #fixed parameter
    policyA1 <- data.frame("A1" = rep(1, size * 2), "A2" = rep(0, size * 2), "A3" = rep(0, size * 2))
    policyA2 <- data.frame("A1" = rep(0, size * 2), "A2" = rep(1, size * 2), "A3" = rep(0, size * 2))
    policyA3 <- data.frame("A1" = rep(0, size * 2), "A2" = rep(0, size * 2), "A3" = rep(1, size * 2))

    Y1.ATE <- -Inf
    fixedParameterPolicyResults <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
    for (newPolicy in list(policyA1, policyA2, policyA3)){
      newPolicyResult <- evaluatePolicy(newPolicy, homePath, treatmentList, outcomeList)
      
      if (newPolicyResult$Y1[[1]] > Y1.ATE &
          newPolicyResult$Y2[[1]] >= -0.1 & 
          newPolicyResult$Y3[[1]] >= -0.1){
        Y1.ATE <- newPolicyResult$Y1[[1]]
        fixedParameterPolicyResults <- newPolicyResult
      }
    }
    heuristicStochasticPolicyResultDF <- rbind(heuristicStochasticPolicyResultDF, heuristicStochasticPolicyResults)
    fixedParameterPolicyResultDF <- rbind(fixedParameterPolicyResultDF, fixedParameterPolicyResults)
    causaltreeMergedStochasticPolicyResultDF <- rbind(causaltreeMergedStochasticPolicyResultDF, causaltreeMergedStochasticPolicyResults)
    causalTreeDeterministicPolicyResultDF <- rbind(causalTreeDeterministicPolicyResultDF, causalTreeDeterministicpolicyResult)
    causalForestDeterministicPolicyResultDF <- rbind(causalForestDeterministicPolicyResultDF, causalForestDeterministicpolicyResult)
    
    deltaModelPolicyResultDF <- rbind(deltaModelPolicyResultDF, deltaModelpolicyResult)
    write.csv(heuristicStochasticPolicyResultDF, file = paste(homePath, uncertaintyW, 'heuristicStochasticPolicyResults.csv', sep = ""), row.names = FALSE, quote = FALSE)
    write.csv(causaltreeMergedStochasticPolicyResultDF, file = paste(homePath, uncertaintyW, 'causaltreeMergedStochasticPolicyResults.csv', sep = ""), row.names = FALSE, quote = FALSE)
    write.csv(fixedParameterPolicyResultDF, file = paste(homePath, uncertaintyW, 'fixedParameterResults.csv', sep = ""), row.names = FALSE, quote = FALSE)
    write.csv(causalTreeDeterministicPolicyResultDF, file = paste(homePath, uncertaintyW, 'causaltreeDeterministicPolicyResults.csv', sep = ""), row.names = FALSE, quote = FALSE)
    write.csv(causalForestDeterministicPolicyResultDF, file = paste(homePath, uncertaintyW, 'causalForestDeterministicPolicyResults.csv', sep = ""), row.names = FALSE, quote = FALSE)
    write.csv(causaltreeStochasticPolicyResultDF, file = paste(homePath, uncertaintyW, 'causaltreeStochasticPolicyResults.csv', sep = ""), row.names = FALSE, quote = FALSE)
    write.csv(deltaModelPolicyResultDF, file = paste(homePath, uncertaintyW, 'deltaModelPolicyResults.csv', sep = ""), row.names = FALSE, quote = FALSE)
  }
}


fullResultSummary <- data.frame("method" = as.character(), "uncertaintyW" = as.numeric(),"outcome" = as.character(),
                                "ATE" = as.numeric(), "SD" = as.numeric())


for (uncertaintyW in uncertaintyWList){
  heuristicStochasticPolicyResultDF <- readData(paste0(homePath, uncertaintyW, 'heuristicStochasticPolicyResults.csv', sep = ""))
  causaltreeMergedStochasticPolicyResultDF <- readData(paste0(homePath, uncertaintyW, 'causaltreeMergedStochasticPolicyResults.csv', sep = ""))
  causalTreeDeterministicPolicyResultDF <- readData(paste0(homePath, uncertaintyW, 'causaltreeDeterministicPolicyResults.csv', sep = ""))
  causalForestDeterministicPolicyResultDF <- readData(paste0(homePath, uncertaintyW, 'causalForestDeterministicPolicyResults.csv', sep = ""))
  deltaModelPolicyResultDF <- readData(paste0(homePath, uncertaintyW, 'deltaModelPolicyResults.csv', sep = ""))
  fixedParameterPolicyResultDF <- readData(paste0(homePath, uncertaintyW, 'fixedParameterResults.csv', sep = ""))
  

  if (checkResultData(fixedParameterPolicyResultDF)){
    fullResultSummary <- rbind(fullResultSummary, getSummaryStats("Global",uncertaintyW, fixedParameterPolicyResultDF, size))
  }
  
  if (checkResultData(heuristicStochasticPolicyResultDF)){
    fullResultSummary <- rbind(fullResultSummary, getSummaryStats("HT.ST", uncertaintyW, heuristicStochasticPolicyResultDF, size))
  }
  
  if (checkResultData(causaltreeMergedStochasticPolicyResultDF)){
    fullResultSummary <- rbind(fullResultSummary, getSummaryStats("CT.ST", uncertaintyW, causaltreeMergedStochasticPolicyResultDF, size))
  }
  
  if (checkResultData(causalForestDeterministicPolicyResultDF)){
    fullResultSummary <- rbind(fullResultSummary, getSummaryStats("CF.DT", uncertaintyW, causalForestDeterministicPolicyResultDF, size))
  }
  
  if (checkResultData(deltaModelPolicyResultDF)){
    fullResultSummary <- rbind(fullResultSummary, getSummaryStats("TM.DT", uncertaintyW, deltaModelPolicyResultDF, size))
  }
}

plotResultSummary(subset(fullResultSummary, outcome == "Y1"), paste0("Evaluation on the objective metric Y0"))

plotResultSummary(subset(fullResultSummary, outcome == "Y2"), paste0("Evaluation on the constraint metric Y1"))

plotResultSummary(subset(fullResultSummary, outcome == "Y3"), paste0("Evaluation on the constraint metric Y2"))

