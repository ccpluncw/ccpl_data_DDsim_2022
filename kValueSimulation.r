#Load packages
library(dplyr)
library(pracma)
library(RRW)

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

#Set the coefficient of variation
CV1 <- .2
#set your simulated k value here
kVal <- 0.07

rrwLoops <- 500
#Set switch for test runs versus full runs
testing <- FALSE

  ### Create simulated data for money
  DollarVec <- seq(5,95, 10)
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
  #Calculate the value for each delay assuming the "k function" accurately describes delay discounting
  DelayNames <- getDelay(DollarVecExt,kVal,CV1)
  DelayMeansVec <- Dollars100*2 - DollarVecExt

  DelayMeanRound <- round_any(DelayMeansVec,round.val)
  Delays <- data.frame(DelayNames, DelayMeanRound)

  #set simulation constants
  N.ObsInDist <- 20000
  N.Samples <- 20000

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
      j <- 1

        Probe <- x
        ProbeVal <- Probes[Probes$MeanNames == x,"ProbeMeansRound"]

        DelayName <- y
        DelayVal <- Delays[Delays$DelayNames == y, "DelayMeanRound"]
        DelayVal <- as.numeric(DelayVal)

        ### Delayed Probe Distribution- Create Probe and Cost Distrubutions then combine them
        #Create distributions

        B.tmp <- rnorm(N.ObsInDist,mean=ProbeVal,sd=(ProbeVal*CV1))
        ADist <- rnorm(N.ObsInDist,mean=Dollars100,sd=(Dollars100*CV1))
        BCost <- rnorm(N.ObsInDist,mean=DelayVal,sd=(DelayVal*CV1))
        BDist <- rowMeans(cbind(B.tmp,BCost))

        A <- sample(ADist,N.Samples,replace=T)
        #Select value from large Distribution
        B <- sample(BDist,N.Samples,replace=T)

        AHigher <- ifelse(A > B,1,ifelse(A < B, 0, .5))

        # Proportion of draws for which A > B
        MeanP <- mean(AHigher)
        # Calculate Overlap
        Overlap <- 1-(abs(MeanP-0.5)/0.5)
        ProbeMean <- x
        # Directional Overlap
        direction <- ifelse(MeanP > .5, 1, -1)
        j <- j + 1

        tmp.df <- data.frame(DelayName=DelayName, DelayVal=DelayVal, Standard = "100", StandardVal = mean(ADist), Probe=Probe,ProbeVal=mean(BDist),Overlap=Overlap,direction = direction)
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

###Unbiased

  print("starting Unbiased Standard")
  probes.tmp <- unique(DelayDat.tmp.Now$Probe)
  df.s.tmp <- NULL
  unbiased.Now <- NULL
  for(x in probes.tmp) {
    print(x)
    dat.tmp <- DelayDat.tmp.Now[DelayDat.tmp.Now$Probe == x,]
    #HVO is probe.
    tmp.out <- getMomentsOfRRWoverlap(dat.tmp$Overlap,boundary=15, startValue= 0, loops = rrwLoops,noiseSD = 2)
    df.s.tmp <- merge(dat.tmp, tmp.out, ,by.x="Overlap",by.y="overlap")
    #convert pCross p(HVO) to pStandard (which is the LVO)
    df.s.tmp$pStandard <- 1-df.s.tmp$pCross
    if (is.null(df.s.tmp)) {
      unbiased.Now <- df.s.tmp
    } else {
      unbiased.Now <- rbind(unbiased.Now, df.s.tmp)
    }
  }


  print("starting Unbiased Delay")
  probes.tmp <- unique(DelayDat.tmp.Delay$Probe)
  df.s.tmp <- NULL
  unbiased.Delay <- NULL
  for(x in probes.tmp) {
    print(x)
    dat.tmp <- DelayDat.tmp.Delay[DelayDat.tmp.Delay$Probe == x,]
    #HVO is standard.
    tmp.out <- getMomentsOfRRWoverlap(dat.tmp$Overlap,boundary=15, startValue= 0, loops = rrwLoops,noiseSD = 2)
    df.s.tmp <- merge(dat.tmp, tmp.out, ,by.x="Overlap",by.y="overlap")
    #set pCross p(HVO) to pStandard (which is the HVO)
    df.s.tmp$pStandard <- df.s.tmp$pCross
    if (is.null(df.s.tmp)) {
      unbiased.Delay <- df.s.tmp
    } else {
      unbiased.Delay <- rbind(unbiased.Delay, df.s.tmp)
    }
  }

    unbiased <- rbind(unbiased.Now,unbiased.Delay)
    DelayDat.unbiased <- subset(unbiased, select = -c(sd,Q25,Q50,Q75,boundary,startValue,noiseSD,decayAsymptote,decayBeta))
### Unbiased Plot
    DelayDat.unbiased$Delay <- DelayDat.unbiased$DelayName
    DelayDat.unbiased <-DelayDat.unbiased[order(DelayDat.unbiased$Delay),]

    numMeans <- length(MeanNames)
    c1 <- vector(length=numMeans)
    for(ct in 1:numMeans) {
      c1[ct] <- (.9/numMeans) * ct
    }

   if(!testing) {
### Immediacy Bias
    print("starting Immediacy Bias Standard")
    probes.tmp <- unique(DelayDat.tmp.Now$Probe)
    df.s.tmp <- NULL
    Immediacy.Now <- NULL
    for(x in probes.tmp) {
      print(x)
      dat.tmp <- DelayDat.tmp.Now[DelayDat.tmp.Now$Probe == x,]
      tmp.out <- getMomentsOfRRWoverlap(dat.tmp$Overlap,boundary=15, startValue=0.5, loops = rrwLoops,noiseSD = 2)
      df.s.tmp <- merge(dat.tmp, tmp.out, ,by.x="Overlap",by.y="overlap")
      df.s.tmp$pStandard <- 1-df.s.tmp$pCross
      if (is.null(df.s.tmp)) {
        Immediacy.Now <- df.s.tmp
      } else {
        Immediacy.Now <- rbind(Immediacy.Now, df.s.tmp)
      }
    }


    print("starting Immediacy Bias Delay")
    probes.tmp <- unique(DelayDat.tmp.Delay$Probe)
    df.s.tmp <- NULL
    Immediacy.Delay <- NULL
    for(x in probes.tmp) {
      print(x)
      dat.tmp <- DelayDat.tmp.Delay[DelayDat.tmp.Delay$Probe == x,]
      tmp.out <- getMomentsOfRRWoverlap(dat.tmp$Overlap,boundary=15, startValue= -0.5, loops = rrwLoops,noiseSD = 2)
      df.s.tmp <- merge(dat.tmp, tmp.out, ,by.x="Overlap",by.y="overlap")
      df.s.tmp$pStandard <- df.s.tmp$pCross
      if (is.null(df.s.tmp)) {
        Immediacy.Delay <- df.s.tmp
      } else {
        Immediacy.Delay <- rbind(Immediacy.Delay, df.s.tmp)
      }
    }

      Immediacy <- rbind(Immediacy.Now,Immediacy.Delay)
      DelayDat.immediacy <- subset(Immediacy, select = -c(sd,Q25,Q50,Q75,boundary,startValue,noiseSD,decayAsymptote,decayBeta))
### Immediacy Plot
      DelayDat.immediacy$Delay <- DelayDat.immediacy$DelayName
      DelayDat.immediacy <- DelayDat.immediacy[order(DelayDat.immediacy$Delay),]


### Bias Towards Delay
      print("starting Delay Bias Standard")
      probes.tmp <- unique(DelayDat.tmp.Now$Probe)
      df.s.tmp <- NULL
      delayed.Now <- NULL
      for(x in probes.tmp) {
        print(x)
        dat.tmp <- DelayDat.tmp.Now[DelayDat.tmp.Now$Probe == x,]
        tmp.out <- getMomentsOfRRWoverlap(dat.tmp$Overlap,boundary=15, startValue=-0.5, loops = rrwLoops,noiseSD = 2)
        df.s.tmp <- merge(dat.tmp, tmp.out, ,by.x="Overlap",by.y="overlap")
        df.s.tmp$pStandard <- 1-df.s.tmp$pCross
        if (is.null(df.s.tmp)) {
          delayed.Now <- df.s.tmp
        } else {
          delayed.Now <- rbind(delayed.Now, df.s.tmp)
        }
      }


      print("starting Delay Bias Delay")
      probes.tmp <- unique(DelayDat.tmp.Delay$Probe)
      df.s.tmp <- NULL
      delayed.Delay <- NULL
      for(x in probes.tmp) {
        print(x)
        dat.tmp <- DelayDat.tmp.Delay[DelayDat.tmp.Delay$Probe == x,]
        tmp.out <- getMomentsOfRRWoverlap(dat.tmp$Overlap,boundary=15, startValue=0.5, loops = rrwLoops,noiseSD = 2)
        df.s.tmp <- merge(dat.tmp, tmp.out, ,by.x="Overlap",by.y="overlap")
        df.s.tmp$pStandard <- df.s.tmp$pCross
        if (is.null(df.s.tmp)) {
          delayed.Delay <- df.s.tmp
        } else {
          delayed.Delay <- rbind(delayed.Delay, df.s.tmp)
        }
      }

        Delayed <- rbind(delayed.Now,delayed.Delay)
        DelayDat.delayed <- subset(Delayed, select = -c(sd,Q25,Q50,Q75,boundary,startValue,noiseSD,decayAsymptote,decayBeta))
      ### Bias Towards Delay Plot

        DelayDat.delayed$Delay <- DelayDat.delayed$DelayName
        DelayDat.delayed <-DelayDat.delayed[order(DelayDat.delayed$Delay),]
}

### Get x y pair for each probe and plot by delay (value on y and delay on x), fit k with nonlinear least squares
### Use approx to find k values
        ## Unbiased
        ApproxVal.1 <- NULL
        for(a in unique(DelayDat.unbiased$Probe)){
          tmp.dat <- DelayDat.unbiased[DelayDat.unbiased$Probe == a & DelayDat.unbiased$correct == T,]
          x <- tmp.dat$Delay
          y <- tmp.dat$pStandard
          approx_val <- approx(x,y,  n=1000)
          indif.df <- data.frame(approx_val$x,approx_val$y)
          val.delay <- indif.df[which.min(abs(approx_val$y - 0.50)),1]
          val.pStandard <- indif.df[which.min(abs(approx_val$y - 0.50)),2]
          df.tmp <- data.frame(DolVal = a, delay = val.delay, pStand = val.pStandard)

          if (is.null(ApproxVal.1)) {
            ApproxVal.1 <- df.tmp
          } else {
            ApproxVal.1 <- rbind(ApproxVal.1, df.tmp)
          }
        }

if(!testing) {

		## Bias Towards immediacy
        ApproxVal.2 <- NULL
        for(a in unique(DelayDat.immediacy$Probe)){
          tmp.dat <- DelayDat.immediacy[DelayDat.immediacy$Probe == a & DelayDat.immediacy$correct == T,]
          x <- tmp.dat$Delay
          y <- tmp.dat$pStandard
          approx_val <- approx(x,y, n=10000)
          indif.df <- data.frame(approx_val$x,approx_val$y)
          val.delay <- indif.df[which.min(abs(approx_val$y - 0.50)),1]
          val.pStandard <- indif.df[which.min(abs(approx_val$y - 0.50)),2]
          df.tmp <- data.frame(DolVal = a, delay = val.delay, pStand = val.pStandard)

          if (is.null(ApproxVal.2)) {
            ApproxVal.2 <- df.tmp
          } else {
            ApproxVal.2 <- rbind(ApproxVal.2, df.tmp)
          }
        }

        ## Bias Towards Delay
        ApproxVal.3 <- NULL
        for(a in unique(DelayDat.delayed$Probe)){
          tmp.dat <- DelayDat.delayed[DelayDat.delayed$Probe == a & DelayDat.delayed$correct == T,]
          x <- tmp.dat$Delay
          y <- tmp.dat$pStandard
          approx_val <- approx(x,y, n=1000)
          indif.df <- data.frame(approx_val$x,approx_val$y)
          val.delay <- indif.df[which.min(abs(approx_val$y - 0.50)),1]
          val.pStandard <- indif.df[which.min(abs(approx_val$y - 0.50)),2]
          df.tmp <- data.frame(DolVal = a, delay = val.delay, pStand = val.pStandard)

          if (is.null(ApproxVal.3)) {
            ApproxVal.3 <- df.tmp
          } else {
            ApproxVal.3 <- rbind(ApproxVal.3, df.tmp)
          }
        }
      }
      ## Calculate k values
      ### Unbiased
        ApproxVal.1 <- ApproxVal.1[ApproxVal.1$pStand < .53 & ApproxVal.1$pStand > 0.47,]
        nls.sum.tmp1 <- nls(formula=DolVal ~ 100/(1+(k*delay)),data=ApproxVal.1,start=list(k=0.1))
        nls.sum.unbiased <- summary(nls.sum.tmp1)
        #k value
        ApproxVal.1$k <- nls.sum.unbiased$coef[1]
        ApproxVal.1$Bias <- "unbiased"
        ApproxVal.1$Fit <- fitted(nls.sum.tmp1)

if(!testing) {
    ### Bias towards Immediacy
        ApproxVal.2 <- ApproxVal.2[ApproxVal.2$pStand < .53 & ApproxVal.1$pStand > 0.47,]
        nls.sum.tmp2 <- nls(formula=DolVal ~ 100/(1+(k*delay)),data=ApproxVal.2,start=list(k=0.00001))
        nls.sum.immediacy <- summary(nls.sum.tmp2)
        ApproxVal.2$k <- nls.sum.immediacy$coef[1]
        ApproxVal.2$Bias <- "immediacy"
        ApproxVal.2$Fit <- fitted(nls.sum.tmp2)

      ### Bias Towards Delay
        ApproxVal.3 <- ApproxVal.3[ApproxVal.3$pStand < .53 & ApproxVal.1$pStand > 0.47,]
        nls.sum.tmp3 <- nls(formula=DolVal ~ 100/(1+(k*delay)),data=ApproxVal.3,start=list(k=0.00001))
        nls.sum.delayed <- summary(nls.sum.tmp3)
        ApproxVal.3$k <- nls.sum.delayed$coef[1]
        ApproxVal.3$Bias <- "delayed"
        ApproxVal.3$Fit <- fitted(nls.sum.tmp3)

        k.df <- rbind(ApproxVal.1,ApproxVal.2,ApproxVal.3)
}
if(testing) {
  k.df <- ApproxVal.1
}
###points

    pdf(file="kValues x Bias.pdf")
        xHigh <- 1.1*max(k.df$delay)
        k.tmp <- k.df[k.df$Bias == "unbiased",]
        with(k.tmp,plot(DolVal ~ delay ,pch=16,frame.plot=F, xlim=c(0,300),ylim=c(0,100),xlab="Delay",ylab="Dollars",las=1))
        with(k.tmp,lines(c(1:xHigh), kFit(100,c(1:xHigh), k.tmp$k[1]) ))

if(!testing) {
        k.tmp <- k.df[k.df$Bias == "immediacy",]
        with(k.tmp,points(DolVal ~ delay, pch=16,col="red"))
        with(k.tmp,lines(c(1:xHigh), kFit(100,c(1:xHigh), k.tmp$k[1]), col="red") )

        k.tmp <- k.df[k.df$Bias == "delayed",]
        with(k.tmp,points(DolVal ~ delay,pch=16,col="blue"))
        with(k.tmp,lines(c(1:xHigh), kFit(100,c(1:xHigh), k.tmp$k[1]),col="blue" ))
        legend(200,80, legend=c("Bias Towards Immediacy","Unbiased","Bias Towards Delay"),title="Start Point Bias",col=c("red","black","blue"),pch=16,cex=0.75)
}
    dev.off()

      ## Write dataframes to a txt file
          sink("kvalues.txt",append=F)
          print(k.df)
          sink(NULL)
