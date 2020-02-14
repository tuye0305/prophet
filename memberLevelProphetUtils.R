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
      estimations <- predict(opfit, testData)
      write.csv(estimations, file = paste(homePath, treatment, outcome, 'causalTreeEstimations.csv', sep = ""), row.names = FALSE, quote = FALSE)
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
      fitTreatment <- rpart(ff, data = subset(dataNew, lixTreatment == 1.0), method = "anova", control =	controlSet)
      fitControl <- rpart(ff, data = subset(dataNew, lixTreatment == 0.0), method = "anova", control = controlSet)
                                                  
      # Get the minimal xerror
      opcpTreatment <- fitTreatment$cptable[, 1][which.min(fitTreatment$cptable[, 4])]
      opcpControl <- fitControl$cptable[, 1][which.min(fitControl$cptable[, 4])]
      
      # Prune the tree based on the minimal xerror
      opfitTreatment <- prune(fitTreatment, opcpTreatment)
      opfitControl <- prune(fitControl, opcpControl)
      
      
      estimations <- predict(opfitTreatment, testData) - predict(opfitControl, testData)
      write.csv(estimations, file = paste(homePath, treatment, outcome, 'deltaModelEstimations.csv', sep = ""), row.names = FALSE, quote = FALSE)
    }
  }
}

generatePolicy <- function(homePath, fileName, treatmentList, outcomeList, objective.outcome){
  f.obj <- c()
  f.conFull <- c()  
  constraintList <- outcomeList[!outcomeList %in% c(objective.outcome)]
  numTreatments <- length(treatmentList)

  for (treatment in treatmentList){
    estimations <- fread(paste0(homePath, treatment, objective.outcome, fileName), sep = ",", stringsAsFactors = FALSE, header = TRUE)
    f.obj <- rbind(f.obj, - estimations)
  }
  numMembers <- nrow(estimations)
  
  for (outcome in constraintList){
    f.con <- c()
    for (treatment in treatmentList){
      estimations <- fread(paste0(homePath, treatment, outcome, fileName), sep = ",", stringsAsFactors = FALSE, header = TRUE)
      f.con <- rbind(f.con, - estimations)
    }
    f.conFull <- cbind(f.conFull, f.con)
  }
  memberLevelConstraint <- c()
  for (treatment in treatmentList){
    
    memberLevelConstraint <- rbind(memberLevelConstraint, diag(x = 1, numMembers, numMembers))
  }
  f.conFull <-cbind(f.conFull, cbind( memberLevelConstraint, - memberLevelConstraint) )

  f.dir <- rep("<=", length(constraintList) + 2 * numMembers)
  f.rhs <- c( 0, 0, rep(1.001, numMembers), rep(- 0.999, numMembers))

  sol <- lp("min", unlist(f.obj), t(f.conFull), f.dir, f.rhs)

  results <- matrix(sol$solution, nrow = numMembers)

  policy <- data.frame("A1" = results[, 1], "A2" = results[, 2], "A3" = results[, 3],  "A4" = results[, 4])
  
  return(policy)
}

