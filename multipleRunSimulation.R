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

homePath <- "" #TODO: Please specify the home directory for the project

source(paste0(homePath, "generateSimulationDataUtils.R"))
source(paste0(homePath, "memberLevelProphetUtils.R"))
source(paste0(homePath, "cohortLevelProphetUtils.R"))

outcomeList <- c("Y1", "Y2", "Y3")
treatmentList <- c("A1", "A2", "A3", "A4")


seed <- 313134
set.seed(seed)
size <- 500
uncertaintyW <- 1
w1List <- c(rnorm(1), rnorm(1), rnorm(1), rnorm(1))
w2List <-  c(rnorm(1), rnorm(1), rnorm(1), rnorm(1))

minLeafSize<- 50
numTreatment <- length(treatmentList)
objective.outcome <- "Y1"

numEpoc <- 2
for (uncertaintyW in c(1, 2, 3, 4, 5, 6)){
  causaltreeStochasticPolicyResultDF <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
  causalTreeDeterministicPolicyResultDF <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
  deltaModelPolicyResultDF <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
  fixedParameterPolicyResultDF <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
  for (e in 1:numEpoc){
    # one epoc
    seedList <- lapply( 1:numTreatment, function(x) c(x * e * 5, x * e * 1999))
    generateSimulationData(homePath, treatmentList, w1List, w2List, seedList, uncertaintyW)
    trainCausalTree(homePath, outcomeList, treatmentList, minLeafSize)
    
    Y1.ATE <- -Inf
    causaltreeStochasticPolicyResults <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
    for (treatment in treatmentList){
      decisionTable <- readData(paste0(homePath, treatment, objective.outcome, 'decisionTable.csv'))
      if (!is.na(decisionTable)){
        utilityDF <- calculateUtilities(decisionTable, homePath, outcomeList[[1]], treatmentList)
        utilityDF <- subset(utilityDF, select = c("cohort", "treatment"))
        for (outcome in outcomeList){
          utilities <- calculateUtilities(decisionTable, homePath, outcome, treatmentList)
          utilityDF[[outcome]] <- as.character(utilities[["utilities"]])  
        }
        utilityDF <- utilityDF[order(utilityDF$cohort),]
        tau <- 0.001 # hyperparameter for constraints violation level
        alpha <- 0.01 # hyperparameter of the learning rate
        
        stochasticOptimization(homePath, utilityDF, numCohorts = nrow(decisionTable), numTreatments = length(treatmentList))
        finalAssignment <- readData(paste0(homePath, '/FinalAssignments', tau, "-", alpha, '.csv'))
        policy <- scoreTestDataGeneratePolicy(decisionTable, homePath, treatmentList,finalAssignment, method = "SGD")
        newPolicyResult <- evaluatePolicy(policy, homePath, treatmentList, outcomeList)
        
        if (!is.na(newPolicyResult$Y1[[1]]) & newPolicyResult$Y1[[1]] > Y1.ATE){
          Y1.ATE <- newPolicyResult$Y1[[1]]
          causaltreeStochasticPolicyResults <- newPolicyResult
        }
      }
    }
    
    #member level causal tree 
    generateCausalTreeEstimations(homePath, outcomeList, treatmentList, minLeafSize)
    causalTreeDeterministicpolicy <- generatePolicy(homePath, 'causalTreeEstimations.csv', treatmentList, outcomeList, objective.outcome)
    
    causalTreeDeterministicpolicyResult <- evaluatePolicy(causalTreeDeterministicpolicy, homePath, treatmentList, outcomeList)
    #member level delta model
    generateTwoDeltaModelEstimations(homePath, outcomeList, treatmentList, minLeafSize)
    deltaModelpolicy <- generatePolicy(homePath, 'deltaModelEstimations.csv', treatmentList, outcomeList, objective.outcome)
    deltaModelpolicyResult <- evaluatePolicy(deltaModelpolicy, homePath, treatmentList, outcomeList)
    
    #fixed parameter
    policyA1 <- data.frame("A1" = rep(1, size * 2), "A2" = rep(0, size * 2), "A3" = rep(0, size * 2), "A4" = rep(0, size * 2))
    policyA2 <- data.frame("A1" = rep(0, size * 2), "A2" = rep(1, size * 2), "A3" = rep(0, size * 2), "A4" = rep(0, size * 2))
    policyA3 <- data.frame("A1" = rep(0, size * 2), "A2" = rep(0, size * 2), "A3" = rep(1, size * 2), "A4" = rep(0, size * 2))
    policyA4 <- data.frame("A1" = rep(0, size * 2), "A2" = rep(0, size * 2), "A3" = rep(0, size * 2), "A4" = rep(1, size * 2))
    
    Y1.ATE <- -Inf
    fixedParameterPolicyResults <- data.frame("Y1" = as.numeric(), "Y2" = as.numeric(), "Y3" = as.numeric())
    for (newPolicy in list(policyA1, policyA2, policyA3, policyA4)){
      newPolicyResult <- evaluatePolicy(newPolicy, homePath, treatmentList, outcomeList)
      
      if (newPolicyResult$Y1[[1]] > Y1.ATE &
          newPolicyResult$Y2[[1]] >= -0.05 & 
          newPolicyResult$Y3[[1]] >= -0.05){
        Y1.ATE <- newPolicyResult$Y1[[1]]
        fixedParameterPolicyResults <- newPolicyResult
      }
    }
    fixedParameterPolicyResultDF <- rbind(fixedParameterPolicyResultDF, fixedParameterPolicyResults)
    causalTreeDeterministicPolicyResultDF <- rbind(causalTreeDeterministicPolicyResultDF, causalTreeDeterministicpolicyResult)
    causaltreeStochasticPolicyResultDF <- rbind(causaltreeStochasticPolicyResultDF, causaltreeStochasticPolicyResults)
    deltaModelPolicyResultDF <- rbind(deltaModelPolicyResultDF, deltaModelpolicyResult)
    write.csv(fixedParameterPolicyResultDF, file = paste(homePath, uncertaintyW, 'fixedParameterResults.csv', sep = ""), row.names = FALSE, quote = FALSE)
    write.csv(causalTreeDeterministicPolicyResultDF, file = paste(homePath, uncertaintyW, 'causaltreeDeterministicPolicyResults.csv', sep = ""), row.names = FALSE, quote = FALSE)
    write.csv(causaltreeStochasticPolicyResultDF, file = paste(homePath, uncertaintyW, 'causaltreeStochasticPolicyResults.csv', sep = ""), row.names = FALSE, quote = FALSE)
    write.csv(deltaModelPolicyResultDF, file = paste(homePath, uncertaintyW, 'deltaModelPolicyResults.csv', sep = ""), row.names = FALSE, quote = FALSE)
  }
}


fullResultSummary <- data.frame("method" = as.character(), "uncertaintyW" = as.numeric(),"outcome" = as.character(),
                                "ATE" = as.numeric(), "SD" = as.numeric())

checkResultData <- function(data){
  (is.data.frame(data) && nrow(data)!=0) || (!is.data.frame(data) && !is.na(data))
}
for (uncertaintyW in c(1, 2, 3, 4, 5, 6)){
  causaltreeStochasticPolicyResultDF <- readData(paste0(homePath, uncertaintyW, 'causaltreeStochasticPolicyResults.csv', sep = ""))
  causalTreeDeterministicPolicyResultDF <- readData(paste0(homePath, uncertaintyW, 'causaltreeDeterministicPolicyResults.csv', sep = ""))
  deltaModelPolicyResultDF <- readData(paste0(homePath, uncertaintyW, 'deltaModelPolicyResults.csv', sep = ""))
  fixedParameterPolicyResultDF <- readData(paste0(homePath, uncertaintyW, 'fixedParameterResults.csv', sep = ""))
  
  if (checkResultData(causaltreeStochasticPolicyResultDF)){
    fullResultSummary <- rbind(fullResultSummary, getSummaryStats("CT.ST", uncertaintyW, causaltreeStochasticPolicyResultDF, size))
  }
  if (checkResultData(causalTreeDeterministicPolicyResultDF)){
    fullResultSummary <- rbind(fullResultSummary, getSummaryStats("CT.DT",uncertaintyW, causalTreeDeterministicPolicyResultDF, size))
  }
  if (checkResultData(deltaModelPolicyResultDF)){
    fullResultSummary <- rbind(fullResultSummary, getSummaryStats("TM.DT",uncertaintyW, deltaModelPolicyResultDF, size))
  }
  if (checkResultData(fixedParameterPolicyResultDF)){
    fullResultSummary <- rbind(fullResultSummary, getSummaryStats("Global",uncertaintyW, fixedParameterPolicyResultDF, size))
  }
}

plotResultSummary(subset(fullResultSummary, outcome == "Y1"), paste0("Evaluation on the objective metric Y0"))

plotResultSummary(subset(fullResultSummary, outcome == "Y2"), paste0("Evaluation on the constraint metric Y1"))

plotResultSummary(subset(fullResultSummary, outcome == "Y3"), paste0("Evaluation on the constraint metric Y2"))

