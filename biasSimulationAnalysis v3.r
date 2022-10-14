#Load packages
library(dplyr)
library(pracma)
library(chutils)


#set your simulated k value here
kValues <- c(0.07)

mainDir <- getwd()

for(kVal in kValues) {
  subDir <- paste("k=", kVal)
  ch.newDir (mainDir, subDir)

  source(paste(mainDir, "biasSimulationRTxDelay.r", sep = "/"))
  source(paste(mainDir, "biasSimulationExtractEqualValuesFromData.r", sep = "/"))

  setwd(mainDir)
}
