########################### Generate Simulation Data Helper Function ##########################

createCasualDAG <- function(seed, size, uncertaintyW, w1, w2, w3, w4, w5, w6, pow){
  D <- DAG.empty()
  D <- D +
    node("H1",
         distr = "rnorm",
         mean = 0.5,
         sd = 1) +      
    node("H2",
         distr = "rnorm",
         mean = (H1^2 * 2),
         sd = 2.5) +
    node("H3",
         distr = "rnorm",
         mean = ( 2 * H1 + 1.4 * H2),
         sd = 2) +
    node("H4",
         distr = "rnorm",
         mean = (-0.5 + 0.7 * H3 + 0.3 * H2),
         sd = 1
         ) +
    node("A",
         distr = "rbern",
         prob = plogis(3.5 + 0.3 * H1 + 1.3 * H2 + 0.2 * H3)) +
    node("U.Y",
         distr = "rnorm",
         mean = 3,
         sd = 1.5 * uncertaintyW) +
    node("Y1",
         distr = "rnorm",
         mean = (U.Y + w1 * A * H2^pow + w2 * A * H3^pow + w3 * A * H4^pow), 
         sd = 0.5 * uncertaintyW) +
    node("Y2",
         distr = "rnorm",
         mean = (w4 * A * H3^pow + U.Y + w5 * A *H4),
         sd = 1.5 * uncertaintyW) +
    node("Y3",
         distr = "rnorm",
         mean = (w6 * A * H1^pow + U.Y * H3),
         sd = 2.5 * uncertaintyW )
  Dset <- set.DAG(D, latent.v = c("U.Y"))
  
  A1 <- node("A", distr = "rbern", prob = 1)
  Dset <- Dset + action("A1", nodes = A1)
  A0 <- node("A", distr = "rbern", prob = 0)
  Dset <- Dset + action("A0", nodes = A0)
  plotDAG(Dset, xjitter = 0.3, yjitter = 0.04,
          edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
          vertex_attrs = list(size = 12, label.cex = 0.8))
  

  Xdat1 <- sim(DAG = Dset, actions = c("A1", "A0"), n = size, rndseed = seed)
  
  Xdat1.A1 <- Xdat1[["A1"]]
  Xdat1.A2 <- Xdat1[["A0"]]
  
  Xdat1.A1["Y1.0"] <- Xdat1.A2["Y1"]
  Xdat1.A1["Y2.0"] <- Xdat1.A2["Y2"]
  Xdat1.A1["Y3.0"] <- Xdat1.A2["Y3"]
  
  
  Xdat1.A2["Y1.0"] <- Xdat1.A2["Y1"]
  Xdat1.A2["Y2.0"] <- Xdat1.A2["Y2"]
  Xdat1.A2["Y3.0"] <- Xdat1.A2["Y3"]
  
  return(list(rbind(Xdat1[["A1"]], Xdat1[["A0"]]), rbind(Xdat1.A1, Xdat1.A2))) 
}

generateSimulationData <- function(homePath, treatmentList, w1List, w2List, w3List, w4List, w5List, w6List, seedList, uncertaintyW){
  typeList <- c("Training", "Test")
  for (t in 1: length(treatmentList)){
    treatment <- treatmentList[[t]]
    for (i in 1: length(typeList)){
      seed <- seedList[[t]][[i]] 
      type <- typeList[[i]]
      simulationData <- createCasualDAG(seed, size, uncertaintyW, w1List[[t]], w2List[[t]], w3List[[t]], w4List[[t]], w5List[[t]], w6List[[t]], pow = 2)
      write.csv(simulationData[[1]], file = paste0(homePath, treatment, 'SimulationData',type, '.csv'), row.names = FALSE, quote = FALSE)
      write.csv(simulationData[[2]], file = paste0(homePath, treatment, 'SimulationDataAll',type, '.csv'), row.names = FALSE, quote = FALSE)
    }
  }
}

readData <- function(inputPath){
  tryCatch( fread(inputPath, sep = ",", stringsAsFactors = FALSE, header = TRUE),
            error = function(e){NA})
}

########################### Process result data Helper Function ##########################


getSummaryStats <- function(method, uncertaintyW, inputDf, size){
  #"method" = as.character(), "uncertaintyW" = as.numeric(),"outcome" = as.character(),"ATE" = as.numeric(), "SD" = as.numeric()
  resultDF <- data.frame()
  for (outcome in c("Y2", "Y3")){
    inputDf <- inputDf[ which(inputDf[[outcome]] >= -0.1),]
  }
  for (outcome in outcomeList){
    resultDF <- rbind(resultDF, data.frame("method" = method, "uncertaintyW" =uncertaintyW, "outcome" = outcome,
                                           "ATE" =  mean(inputDf[[outcome]], na.rm = T),  "SD" = sd(inputDf[[outcome]][!is.na(inputDf[[outcome]])])/(size^(0.5)), na.rm = T)) 
  }
  return(resultDF)
}

plotResultSummary <- function(dfSummary, titleText){
  # Default bar plot
  p<- ggplot(dfSummary, aes(x=uncertaintyW, y=ATE, fill=method, group = method)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(x = uncertaintyW, ymin=ATE-SD, ymax=ATE+SD, group = method), width= 2,
                  position=position_dodge(4.5)) 
  # Finished bar plot
  p <- p+labs(x="uncertainty weight", y = "ITE") + theme_bw() + theme(text = element_text(size=16))
  return(p)
}

