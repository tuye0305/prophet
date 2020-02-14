########################### Generate Simulation Data Helper Function ##########################

createCasualDAG <- function(seed, size, uncertaintyW, w1, w2){
  D <- DAG.empty()
  D <- D +
    node("H1",
         distr = "rcat.b1",
         probs = c(0.5, 0.25, 0.25)) +
    node("H2",
         distr = "rnorm",
         mean = ifelse(H1 == 1, 0, ifelse(H1 == 2, 3, 10)),
         sd = 0.5 * uncertaintyW) +
    node("H3",
         distr = "rnorm",
         mean = (- 2 * H1 + 1.4 * H2),
         sd = 1 * uncertaintyW) +
    node("H4",
         distr = "rbern",
         prob = plogis(-0.5 - 0.7 * H3 + 0.3 * H2)) +
    node("A",
         distr = "rbern",
         prob = plogis(-0.5 - 0.3 * H1 - 0.3 * H2 - 0.2 * H3)) +
    node("U.Y",
         distr = "rnorm",
         mean = 3,
         sd = 0.5 * uncertaintyW) +
    node("Y1",
         distr = "rnorm",
         mean = (-0.1 + 2 * U.Y - w1 * 1.2 * A + w1 * A * H2 + 2 * w2 * A * H3 + 3 * A * H4 * H2 + 0.3 * H2 + 0.2 * H4), 
         sd = 1 * uncertaintyW) +
    node("Y2",
        distr = "rbern",
        prob = plogis(-3 + -w1 * 0.6 * A + 0.1 * H3 + 2* U.Y + 0.3 * H4)) +
    node("Y3",
         distr = "rnorm",
         mean = ( - w2 * 0.5 * A + H1 + U.Y * H3),
         sd = 0.2 * uncertaintyW )
  Dset <- set.DAG(D, latent.v = c("U.Y"))
  
  A1 <- node("A", distr = "rbern", prob = 1)
  Dset <- Dset + action("A1", nodes = A1)
  A0 <- node("A", distr = "rbern", prob = 0)
  Dset <- Dset + action("A0", nodes = A0)
  plotDAG(Dset, xjitter = 0.3, yjitter = 0.04,
          edge_attrs = list(width = 0.5, arrow.width = 0.4, arrow.size = 0.8),
          vertex_attrs = list(size = 12, label.cex = 0.8))
  

  Xdat1 <- sim(DAG = Dset, actions = c("A1", "A0"), n = size, rndseed = seed)
  
  return(rbind(Xdat1[["A1"]], Xdat1[["A0"]]))
}

generateSimulationData <- function(homePath, treatmentList, w1List, w2List, seedList, uncertaintyW){
  typeList <- c("Training", "Test")
  for (t in 1: length(treatmentList)){
    treatment <- treatmentList[[t]]
    for (i in 1: length(typeList)){
      seed <- seedList[[t]][[i]] 
      type <- typeList[[i]]
      simulationData <- createCasualDAG(seed, size, uncertaintyW, w1List[[t]], w2List[[t]])
      write.csv(simulationData, file = paste0(homePath, treatment, 'SimulationData',type, '.csv'), row.names = FALSE, quote = FALSE)
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
  for (outcome in outcomeList){
    resultDF <- rbind(resultDF, data.frame("method" = method, "uncertaintyW" =uncertaintyW, "outcome" = outcome,
                                           "ATE" =  mean(inputDf[[outcome]], na.rm = T),  "SD" = sd(inputDf[[outcome]])/(size^(0.5)), na.rm = T)) 
  }
  return(resultDF)
}

plotResultSummary <- function(dfSummary, titleText){
  # Default bar plot
  p<- ggplot(dfSummary, aes(x=uncertaintyW, y=ATE, fill=method, group = method)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(x = uncertaintyW, ymin=ATE-SD, ymax=ATE+SD, group = method), width= 0.4,
                  position=position_dodge(0.9)) 
  # Finished bar plot
  p <- p+labs(x="uncertainty weight", y = "ATE") + theme_bw() + theme(text = element_text(size=16))
  return(p)
}

