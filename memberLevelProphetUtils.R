generateCausalForestEstimations <- function(homePath, outcomeList, treatmentList, minLeafSize){
  splitAlpha <- 0.5
  cvHonest <- T
  cvAlpha <- 0.5
  for (treatment in treatmentList){
    # Read in input data
    dataNew <- fread(paste0(homePath, treatment, 'SimulationDataTraining.csv'), sep = ",", stringsAsFactors = FALSE, header = TRUE)
    testData <- fread(paste0(homePath, treatment, 'SimulationDataTest.csv'), sep = ",", stringsAsFactors = FALSE, header = TRUE)
    
    names(dataNew)[names(dataNew) == "A"] <- "lixTreatment"
    
    for (outcome in outcomeList){
      names(dataNew)[names(dataNew) == outcome] <- "label"
      
      # Create a formula for the tree model
      cm <- colnames(dataNew)
      rm.list <- c("ID", "label", "lixTreatment", outcomeList[!outcomeList %in% c(outcome)])
      rm.cm <- which(cm %in% rm.list)
      # Create a simple linear formula
      ff <- as.formula(paste0("label ~ ", paste(cm[- rm.cm], collapse = " + ")))
      
      truePropen = (nrow(subset(dataNew, lixTreatment == 1)) * 1.0) / nrow(dataNew)
      
      # Train the forest
      c.forest <- causal_forest(data.matrix(dataNew[,c("H1", "H2", "H3", "H4")]), dataNew$label, dataNew$lixTreatment)
      
      estimations.testing <- predict(c.forest, data.matrix(testData[,c("H1", "H2", "H3", "H4")]))
      estimations.training <- predict(c.forest, data.matrix(dataNew[,c("H1", "H2", "H3", "H4")]))
      write.csv(estimations.training, file = paste(homePath, treatment, outcome, 'causalForestEstimationsTraining.csv', sep = ""), row.names = FALSE, quote = FALSE)
      write.csv(estimations.testing, file = paste(homePath, treatment, outcome, 'causalForestEstimationsTesting.csv', sep = ""), row.names = FALSE, quote = FALSE)
    }
  }
}

generateCausalTreeEstimations <- function(homePath, outcomeList, treatmentList, minLeafSize){
  splitAlpha <- 0.5
  cvHonest <- T
  cvAlpha <- 0.5
  for (treatment in treatmentList){
    # Read in input data
    dataNew <- fread(paste0(homePath, treatment, 'SimulationDataTraining.csv'), sep = ",", stringsAsFactors = FALSE, header = TRUE)
    testData <- fread(paste0(homePath, treatment, 'SimulationDataTest.csv'), sep = ",", stringsAsFactors = FALSE, header = TRUE)
    
    names(dataNew)[names(dataNew) == "A"] <- "lixTreatment"
    
    for (outcome in outcomeList){
      names(dataNew)[names(dataNew) == outcome] <- "label"
      
      # Create a formula for the tree model
      cm <- colnames(dataNew)
      rm.list <- c("ID", "label", "lixTreatment", outcomeList[!outcomeList %in% c(outcome)])
      rm.cm <- which(cm %in% rm.list)
      # Create a simple linear formula
      ff <- as.formula(paste0("label ~ ", paste(cm[- rm.cm], collapse = " + ")))
      
      truePropen = (nrow(subset(dataNew, lixTreatment == 1)) * 1.0) / nrow(dataNew)
      
      # Train the tree
      tree <- causalTree(ff, data = dataNew, treatment = dataNew$lixTreatment, split.alpha = splitAlpha, split.Rule = "CT", cv.option = "TOT", split.Honest = T, cv.Honest = cvHonest, split.Bucket = T, xval = 3, cp = 0, minsize = minLeafSize, propensity = truePropen)
      
      # Get the minimal xerror
      opcp <- tree$cptable[, 1][which.min(tree$cptable[, 4])]
      
      # Prune the tree based on the minimal xerror
      opfit <- prune(tree, opcp)
      estimations.testing <- predict(opfit, testData)
      estimations.training <- predict(opfit, dataNew)
      
      write.csv(estimations.training, file = paste(homePath, treatment, outcome, 'causalTreeEstimationsTraining.csv', sep = ""), row.names = FALSE, quote = FALSE)
      write.csv(estimations.testing, file = paste(homePath, treatment, outcome, 'causalTreeEstimationsTesting.csv', sep = ""), row.names = FALSE, quote = FALSE)
    }
  }
}

generateTwoDeltaModelEstimations <- function(homePath, outcomeList, treatmentList, minLeafSize){
  for (treatment in treatmentList){
    # Read in input data
    dataNew <- fread(paste0(homePath, treatment, 'SimulationDataTraining.csv'), sep = ",", stringsAsFactors = FALSE, header = TRUE)
    testData <- fread(paste0(homePath, treatment, 'SimulationDataTest.csv'), sep = ",", stringsAsFactors = FALSE, header = TRUE)
    
    names(dataNew)[names(dataNew) == "A"] <- "lixTreatment"
    
    for (outcome in outcomeList){
      names(dataNew)[names(dataNew) == outcome] <- "label"
      
      # Create a formula for the tree model
      cm <- colnames(dataNew)
      rm.list <- c("ID", "label", "lixTreatment", outcomeList[!outcomeList %in% c(outcome)])
      rm.cm <- which(cm %in% rm.list)
      # Create a simple linear formula
      ff <- as.formula(paste0("label ~ ", paste(cm[- rm.cm], collapse = " + ")))
      
      # Train the tree
      controlSet <- rpart.control(minsplit = minLeafSize,  xval = 3, cp = 0)
      #fitTreatment <- rpart(ff, data = subset(dataNew, lixTreatment == 1.0), method = "anova", control =	controlSet)
      #fitControl <- rpart(ff, data = subset(dataNew, lixTreatment == 0.0), method = "anova", control = controlSet)
      
      fitTreatment <- randomForest(ff, data = subset(dataNew, lixTreatment == 1.0))
      fitControl <- randomForest(ff, data = subset(dataNew, lixTreatment == 0.0))
      
      # Get the minimal xerror
      #opcpTreatment <- fitTreatment$cptable[, 1][which.min(fitTreatment$cptable[, 4])]
      #opcpControl <- fitControl$cptable[, 1][which.min(fitControl$cptable[, 4])]
      
      # Prune the tree based on the minimal xerror
      #opfitTreatment <- prune(fitTreatment, opcpTreatment)
      #opfitControl <- prune(fitControl, opcpControl)
      
      estimations.testing <- predict(fitTreatment, testData) - predict(fitControl, testData)
      estimations.training <- predict(fitTreatment, dataNew) - predict(fitControl, dataNew)
      
      #estimations.testing <- predict(opfitTreatment, testData) - predict(opfitControl, testData)
      #estimations.training <- predict(opfitTreatment, dataNew) - predict(opfitControl, dataNew)
      
      write.csv(estimations.training, file = paste(homePath, treatment, outcome, 'deltaModelEstimationsTraining.csv', sep = ""), row.names = FALSE, quote = FALSE)
      write.csv(estimations.testing, file = paste(homePath, treatment, outcome, 'deltaModelEstimationsTesting.csv', sep = ""), row.names = FALSE, quote = FALSE)
    }
  }
}

getEstimations <- function(homePath, treatmentList, outcomeList, fileName){
  
  f.obj <- c()
  f.conFull <- c()  
  constraintList <- outcomeList[!outcomeList %in% c(objective.outcome)]
  numTreatments <- length(treatmentList)
  
  for (treatment in treatmentList){
    estimations <- fread(paste0(homePath, treatment, objective.outcome, fileName), sep = ",", stringsAsFactors = FALSE, header = TRUE)
    f.obj <- rbind(f.obj, -estimations)
  }
  numMembers <- nrow(estimations)
  
  for (outcome in constraintList){
    f.con <- c()
    for (treatment in treatmentList){
      estimations <- fread(paste0(homePath, treatment, outcome, fileName), sep = ",", stringsAsFactors = FALSE, header = TRUE)
      f.con <- rbind(f.con, -estimations)
    }
    f.conFull <- cbind(f.conFull, f.con)
  }
  
  return(list(f.obj, f.conFull, numMembers))
}

generatePolicy <- function(homePath, trainingfileName, testingfileName, treatmentList, outcomeList, objective.outcome, tau){
  f.list <- getEstimations(homePath, treatmentList, outcomeList, trainingfileName)
  f.obj <- f.list[[1]]
  f.conFull <- f.list[[2]]
  numMembers <- f.list[[3]]
  
  memberLevelConstraint <- c()
  for (treatment in treatmentList){
    
    memberLevelConstraint <- rbind(memberLevelConstraint, diag(x = 1, numMembers, numMembers))
  }
  f.conFull <-cbind(f.conFull, memberLevelConstraint)

  f.lhs <- c( tau, tau, rep(0.9999, numMembers))
  f.rhs <- c( Inf, Inf, rep(1.0001, numMembers))
  
  settings <- osqpSettings(verbose = TRUE)
  # Solve with OSQP
  gamma <- 0.01
  P <- diag(x = gamma*2, numMembers*3, numMembers*3)
  A <- t(as.matrix(f.conFull))
  l <-  f.lhs
  u <-  f.rhs

  q <- unlist(f.obj)
  res <- solve_osqp(P, q, A, l, u, settings)
  
  lambda <- c(res$y[[1]], res$y[[2]])
  
  f.list.testing <- getEstimations(homePath, treatmentList, outcomeList, testingfileName)
  f.obj.testing <- f.list.testing[[1]]
  f.conFull.testing <- f.list.testing[[2]]
  #dual to primal
  raw.assignment <- (-f.obj.testing - lambda*f.conFull.testing[,1] - lambda*f.conFull.testing[,2])/gamma
  
  res <- projsplx(unlist(raw.assignment), b= 3)
  
  results <- matrix(unlist(raw.assignment), nrow = numMembers)
  policy <- data.frame("A1" = c(), "A2" = c(), "A3" = c())
  for (i in 1:numMembers){
    project.results <- projsplx_2(results[i,])
    policy <- rbind(policy, data.frame("A1" = c(project.results[1]), "A2" = c(project.results[2]), "A3" = c(project.results[3])))
  }
  return(policy)
}

