#Load packages
library(dplyr)
library(pracma)
library(RRW)
library(chutils)
library(foreach)
library(doParallel)

setUpParallel <- function () {

		cores <- parallel::detectCores()
		cl <- parallel::makeCluster(cores, outfile='log.txt')
	  doParallel::registerDoParallel(cl)

  print(paste("Parallel Setup: Number of Cores = ", cores, sep=": "))
	return (cl)
}

endParallel <- function (cl) {
	parallel::stopCluster(cl)
	rm(cl)
	print("cluster stopped")
}

getDelay <- function (val, k, cv) {
  delay = ((100/val) - 1)/k
  delay = delay + rnorm(length(delay),0,cv*delay)
  return(delay)
}

kFit <- function(stand, delay, k) {
  out <-stand/(1+(k*delay))
  return(out)
}

# Create rounding function
  round_any <-  function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
  round.val <- 0.01

#set simulation constants
N.ObsInDist <- 20000
N.Samples <- 20000
#Set the coefficient of variation
CV1 <- .2
#set your simulated k value here
kValues <- c(0.07)
#kVal <- 0.035
#set startValue biases here
startValue <- c(0, 0.25, 0.5, 0.75, -0.25, -0.5, -0.75)

rrwLoops <- 1000
nSD <- 2
boundary <- 15
bootSimLoops <- 25

RTScaleFactor <- 8
#Set switch for test runs versus full runs
testing <- FALSE

  ### Create simulated data for money
  if(testing) {
    DollarVec <- c(35,55)
    startValue <- c(0, 0.5, -0.5)
  } else {
    DollarVec <- seq(5,95, 10)
    startValue <- c(0, 0.25, 0.5, 0.75, -0.25, -0.5, -0.75)
  }

  DollarVec <- trunc(DollarVec)
  MeanNames <- DollarVec
  MeansVec <- DollarVec
  Dollars100 <- 100

  ProbeMeansRound <- round_any(MeansVec,round.val)
  #Produce a dataframe with names and rounded values of simulated data
  Probes <- data.frame(MeanNames, ProbeMeansRound)

  ### Create data for delays
  #extend Dollar vector
  DollarVecExt <- c(0.5, 1, 2.5,DollarVec)
  DelayMeansVec <- Dollars100*2 - DollarVecExt
  DelayMeanRound <- round_any(DelayMeansVec,round.val)


  cl<- setUpParallel ()

  mainDir <- getwd()

  for(kVal in kValues) {
      subDir <- paste("k=", kVal)
      ch.newDir (mainDir, subDir)

      #Calculate the value for each delay assuming the "k function" accurately describes delay discounting
      DelayNames <- getDelay(DollarVecExt,kVal,CV1)
      Delays <- data.frame(DelayNames, DelayMeanRound)


     ### Set variables to null
      Probe <- NULL
      AHigher <- NULL
      Overlap <- NULL
      ProbeMean <- NULL
      ProbeVal <- NULL
      direction <- NULL
      Delay <- NULL
      DelayDat.tmp <- NULL

      for(x in Probes$MeanNames) {
        for(y in Delays$DelayNames) {

            Probe <- x
            ProbeVal <- Probes[Probes$MeanNames == x,"ProbeMeansRound"]

            DelayName <- y
            DelayVal <- Delays[Delays$DelayNames == y, "DelayMeanRound"]
            DelayVal <- as.numeric(DelayVal)

            ### Delayed Probe Distribution- Create Probe and Cost Distrubutions then combine them
            #Create distributions

            df.boot.out <- foreach::foreach(j=1:bootSimLoops, .combine="rbind") %dopar% {
              B.tmp <- rnorm(N.ObsInDist,mean=ProbeVal,sd=(ProbeVal*CV1))
              ADist <- rnorm(N.ObsInDist,mean=Dollars100,sd=(Dollars100*CV1))
              BCost <- rnorm(N.ObsInDist,mean=DelayVal,sd=(DelayVal*CV1))
              BDist <- rowMeans(cbind(B.tmp,BCost))

              A <- sample(ADist,N.Samples,replace=T)
              #Select value from large Distribution
              B <- sample(BDist,N.Samples,replace=T)

              AHigher <- ifelse(A > B,1,ifelse(A < B, 0, .5))

              # Proportion of draws for which A > B
              data.frame(MeanP = mean(AHigher), mADist = mean(ADist), mBDist = mean(BDist))
            }

            MeanP <- mean(df.boot.out$MeanP)
            StandardVal <- mean(df.boot.out$mADist)
            ProbeVal <- mean(df.boot.out$mBDist)

            # Calculate Overlap
            Overlap <- 1-(abs(MeanP-0.5)/0.5)
            ProbeMean <- x
            # Directional Overlap
            direction <- ifelse(MeanP > .5, 1, -1)

            tmp.df <- data.frame(DelayName=DelayName, DelayVal=DelayVal, Standard = "100", StandardVal = StandardVal, Probe=Probe,ProbeVal=ProbeVal,Overlap=Overlap,direction = direction)
          if (is.null(DelayDat.tmp)) {
            DelayDat.tmp <- tmp.df
          } else {
            DelayDat.tmp <- rbind(DelayDat.tmp, tmp.df)
          }
        }
      }


      ### Split DelayDat.tmp by direction to introduce the correct bias, etc.
      ### When direction ==1, then the standard is the HVO
      ### When direction == -1, then the delay is the HVO
      ### If you want the p(HVO) to actually equal the p(Standard), then the simulation has to be run separately for
      ### direction == 1 and direction == -1.

      DelayDat.tmp.Delay <- DelayDat.tmp[DelayDat.tmp$direction == 1,]
      DelayDat.tmp.Now <- DelayDat.tmp[DelayDat.tmp$direction == -1,]


    ##### Run the RRW for each sv x delay for now vs delay.
    HVOtimes <- c("now", "delay")
    delayDat <- NULL
    for(sv in startValue) {
      for(ht in HVOtimes) {
        if(ht == "now") {
          dat.tmp.1 <- DelayDat.tmp.Now
          svIn <- sv
        } else {
          dat.tmp.1 <- DelayDat.tmp.Delay
          svIn <- sv * -1
        }
        cat("starting ", ht, " ", sv, "\n")
        probes.tmp <- unique(dat.tmp.1$Probe)
        df.s.tmp <- NULL
        for(x in probes.tmp) {
          print(x)
          dat.tmp <- dat.tmp.1[dat.tmp.1$Probe == x,]

          #HVO is probe.
          ### do for loop to get stable estimates
          tmp.out.1 <- foreach::foreach(tmL=1:bootSimLoops, .combine="rbind", .export=c('getMomentsOfRRWoverlap'), .packages="RRW") %dopar% {


            getMomentsOfRRWoverlap(dat.tmp$Overlap,boundary=boundary, startValue= svIn, loops = rrwLoops,noiseSD = nSD)
          }

          tmp.out <- data.frame(tmp.out.1 %>% group_by ( overlap,boundary,startValue,noiseSD,decayAsymptote,decayBeta, correct) %>% summarize(mean = mean(mean, na.rm = T), sd = sd(sd, na.rm = T), Q25 = mean(Q25, na.rm = T), Q50 = mean(Q50, na.rm = T), Q75 = mean(Q75, na.rm = T), pCross = mean(pCross, na.rm = T) ))

          df.s.tmp <- merge(dat.tmp, tmp.out, ,by.x="Overlap",by.y="overlap")
          #convert pCross p(HVO) to pStandard (which is the LVO) so an "error" is pStandard
#          df.s.tmp$pStandard <- ifelse(df.s.tmp$correct == T, 1-df.s.tmp$pCross, df.s.tmp$pCross)
          df.s.tmp$startVal <- sv
          df.s.tmp$bias <- ifelse(sv==0, "unbiased", ifelse(sv > 0, "immediacy", "delayed"))
          df.s.tmp$Delay <- df.s.tmp$DelayName
          if (is.null(delayDat)) {
            delayDat <- df.s.tmp
          } else {
            delayDat <- rbind(delayDat, df.s.tmp)
          }
        }
      }
    }

    ### recode the "correct" column to indicate the stimulus the participant chose (now vs delay)
    delayDat$choice <- ifelse((delayDat$direction < 0 & delayDat$correct == TRUE) | (delayDat$direction > 0 & delayDat$correct == FALSE), "now", "delay")
    delayDat$pStandard <- ifelse(delayDat$choice  == "now", pCross, 1-pCross)
    write.table(delayDat, file = "DelayDat.txt", quote = F, sep = "\t", row.names = F)

    setwd(mainDir)
  }

endParallel (cl)
