library(dplyr)
library(pracma)
library(chutils)

# in this code pStandard is the probablity of responding now: p(now)
#return estimated delay from fit function
getEqualValueDelay <- function(data, delayColumn, fitColumn, pStandardColumn, unbiased = F, intervals = 5000) {
  data <- data[order(data[[delayColumn]]),]
  x.n <- data[[delayColumn]]
  fit <- data[[fitColumn]]
  pS <- df.tmp[[pStandardColumn]]

  approx_val.fit <- approx(x.n,fit,  n=intervals)
  approx_val.pS <- approx(x.n,pS,  approx_val.fit$x)

  indif.df.pS <- data.frame(Delay = approx_val.pS$x,pStandard = approx_val.pS$y)
  indif.df.fit <- data.frame(Delay = approx_val.fit$x,Fit = approx_val.fit$y)

  indif.df <- merge (indif.df.pS, indif.df.fit)

  if(unbiased == F) {
    df.tmp.fit <- indif.df[order(indif.df$Fit),]
    df.tmp.fit <- df.tmp.fit[complete.cases(df.tmp.fit),]
    df.tmp.fit.max <- df.tmp.fit[nrow(df.tmp.fit), ]
  } else {
    df.tmp.fit.max <- indif.df[which.min(abs(indif.df$pStandard - 0.50)),]
  }

  return(df.tmp.fit.max)
}

# Create rounding function
  round_any <-  function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
  round.val <- 0.01

  source("../biasSimulationExtractKvalues.r")
  actualKs <-read.table("kvalues Actual from Simulation.txt", sep = "\t", header = T)
  trueK <- mean(actualKs$k, na.rm=T)

  delayDat <- read.table("delayDat.txt", sep = "\t", header = T)
  startValue <- sort(unique(delayDat$startVal))
  delayDat$lDelay <- round(log(delayDat$Delay),2)
  delayDat$pStandard <- ifelse(delayDat$choice  == "now", delayDat$pCross, 1-delayDat$pCross)
  delayDat$RT <- scale(delayDat$mean)


  ############# Here, for each sv, we plot the p(now) x RT and p(delay) x RT. then we fit a polynomial to
  ######  each set of data.  This shows the relation of p(now) to p(delay) and how that relates to RT.

  sink("RT x pChoose fits.text", append =F)
    cat("\n ********** Loess fits of the RT by p(choose stimulus) **********\n\n")
  sink(NULL)

  biases <- c("unbiased", "delayed", "immediacy")
  pdf(file="RT x Probe x Bias.pdf")
  df.fit <- NULL
  maxY <- 1.1*max(delayDat$RT, na.rm=T)
  minY <- .9*min(delayDat$RT, na.rm=T)
  maxY <- round_any(maxY, 1, ceiling)
  minY <- round_any(minY, 1, floor)

  maxX.ld <- 1.1*max(delayDat$lDelay, na.rm=T)
  minX.ld <- .9*min(delayDat$lDelay, na.rm=T)
  maxX.ld <- round_any(maxX.ld, 1, ceiling)
  minX.ld <- round_any(minX.ld, 1, floor)

  df.FitFromData <- NULL
  for(sv in startValue) {
      tmp.dat <- delayDat[delayDat$startVal == sv, ]
      bias <- unique(tmp.dat$bias)
      title1 <- paste(bias,sv[1],"All Probes")
      tmp.df.d <- tmp.dat[tmp.dat$choice == "delay",]
      tmp.df.d <- tmp.df.d[order(tmp.df.d$pStandard),]
      fit.delay.all <- loess(RT~pStandard, span = .5,data=tmp.df.d)
      tmp.df.d$Fit <- predict(fit.delay.all,tmp.df.d)

      tmp.df.n <- tmp.dat[tmp.dat$choice == "now",]
      tmp.df.n <- tmp.df.n[order(tmp.df.n$pStandard),]
      fit.now.all <- loess(RT~pStandard, span = .5, data=tmp.df.n)
      tmp.df.n$Fit <- predict(fit.now.all,tmp.df.d)


      #test for difference in RT at pStandard between .45 nd .55 with predicted fits
      now50 <- predict(fit.now.all,seq(.45,.55,.005))
      delay50 <- predict(fit.delay.all,seq(.45,.55,.005))
      t.out <- t.test(now50,delay50)

      #if the pStandard at 50% is significantly different, then there is a bias, so do the following
      if(t.out$p.value < 0.05) {
        #select the now dataset if there is a bias toward delay (these are the higher RT datasets)
        if(t.out$statistic > 1) {
          df.eq.dat <- tmp.df.n
        } else {
          #select the delay dataset if there is a bias toward now (these are the higher RT datasets)
          df.eq.dat <- tmp.df.d
        }

        #get equal valued data from predicted RT in the Fit function
        #Basically, for each probe, identify the highest fit RT, and that delay is the equal value item
        title <- paste(bias,sv[1])
        with(df.tmp, plot(NULL,NULL, xlim=c(minX.ld,maxX.ld),ylim=c(minY,maxY), pch=16, col="black", frame.plot=F,xlab="log(Delay)",ylab="RT",las=1, type="p", main = title))
        i <- 10
        for(prb in sort(unique(df.eq.dat$Probe))) {
          df.tmp <- df.eq.dat[df.eq.dat$Probe == prb,]
          df.tmp <- df.tmp[order(df.tmp$lDelay),]
          col1 = paste("grey", i, sep = "")
          with(df.tmp, lines(Fit ~ lDelay, col = col1))

          eqDelay <- getEqualValueDelay(df.tmp, "Delay", "Fit", "pStandard", unbiased = F)
          df.tmp.fit.max <- data.frame(Probe = prb, startVal = sv, pStandard = eqDelay$pStandard[1], Delay = eqDelay$Delay[1], Fit = eqDelay$Fit[1])

          if (is.null(df.fit)) {
            df.fit <- df.tmp.fit.max
          } else {
            df.fit <- rbind(df.fit, df.tmp.fit.max)
          }
          i <- i + 10
        }

      } else {
        df.eq.dat <- tmp.df.n
        #if there is no bias, then get equal value data from pStandard
        title <- paste(bias,sv[1])
        with(df.tmp, plot(NULL, NULL, xlim=c(minX.ld,maxX.ld),ylim=c(0,1), pch=16, col="black", frame.plot=F,xlab="log(Delay)",ylab="p(Now)",las=1, type="p", main = title))
        abline(h=0.5, lty = 3, col="grey")
        i <- 10
        for(prb in sort(unique(df.eq.dat$Probe))) {
          df.tmp <- df.eq.dat[df.eq.dat$Probe == prb,]
          df.tmp <- df.tmp[order(df.tmp$lDelay),]
          col1 = paste("grey", i, sep = "")
          with(df.tmp, lines(pStandard ~ lDelay, col=col1))

          eqDelay <- getEqualValueDelay(df.tmp, "Delay", "Fit", "pStandard", unbiased = T)
          df.tmp.fit <- data.frame(Probe = prb, startVal = sv, pStandard = eqDelay$pStandard[1], Delay = eqDelay$Delay[1], Fit = eqDelay$Fit[1])

          i <- i + 10
          if (is.null(df.fit)) {
            df.fit <- df.tmp.fit
          } else {
            df.fit <- rbind(df.fit, df.tmp.fit)
          }
        }
      }

      ### Plot the all probes for start value
        with(tmp.df.d, plot(RT ~ pStandard,xlim=c(0,1),ylim=c(minY,maxY), pch=16, col="black", frame.plot=F,xlab="p(Now)",ylab="RT",las=1, type="p", main = title1))

        with(tmp.df.d, lines(pStandard, predict(fit.delay.all,tmp.df.d), col='black'))
        with(tmp.df.n, points(RT ~ pStandard, pch=16, col="grey"))
        with(tmp.df.n, lines(pStandard, predict(fit.now.all,tmp.df.n), col='grey'))

        eq.val.1 <- tmp.df.n[tmp.df.n$Overlap > .99,]
        eq.val.t <- mean(eq.val.1[eq.val.1$correct == T,"pStandard"])
        abline(v=eq.val.t, lty=3, col="grey25")

        equalValProbes <- df.fit[df.fit$startVal == sv,]
        eq.val.dat <- mean(equalValProbes$pStandard, na.rm = T)
        abline(v=eq.val.dat, lty=3, col="grey75")

      #Fit the hyperbolic function to start value
      nls.sum.tmp1 <- nls(formula=Probe ~ 100/(1+(k*Delay)),data=equalValProbes,start=list(k=0.1))
      sink("RT x pChoose fits.text", append = T)
        cat("\n ******** k Value Fit for  ******** \n\n")
        cat("\n\n **** Bias = ", bias, ";Start Value = ", sv, " ****\n")
        print(summary(nls.sum.tmp1))
        cat("\n k = ", coef(nls.sum.tmp1), "\n\n")
      sink(NULL)
      equalValProbes$Fit <- fitted(nls.sum.tmp1)
      equalValProbes <- equalValProbes[order(equalValProbes$Delay),]
      with(equalValProbes, plot(Probe ~ Delay, xlim=c(0,300),ylim=c(0,100), pch=16, col="black", frame.plot=F,xlab="Delay",ylab="Probe",las=1, type="p", main = title1))
      with(equalValProbes,lines(Delay,Fit))

      df.FitFromData.1 <- data.frame(startValue = sv, k = coef(nls.sum.tmp1), pStandardEqualValue = eq.val.dat, trueK = trueK, pStandardEqualValueTrue = eq.val.t)
      df.FitFromData <- ch.rbind(df.FitFromData, df.FitFromData.1)
  }
  plot(NULL,NULL, xlim=c(0,300),ylim=c(0,100), pch=16, col="black", frame.plot=F,xlab="Delay",ylab="Probe",las=1, type="p", main = "All Fits")
  i <- 12
  for(sv in sort(startValue)) {
    equalValProbes <- df.fit[df.fit$startVal == sv,]
    nls.sum.tmp1 <- nls(formula=Probe ~ 100/(1+(k*Delay)),data=equalValProbes,start=list(k=0.1))
    equalValProbes$Fit <- fitted(nls.sum.tmp1)
    equalValProbes <- equalValProbes[order(equalValProbes$Delay),]
    col1 = paste("grey", i, sep = "")
    with(equalValProbes,lines(Delay,Fit, col = col1))
    i <- i + 12
  }

  with(df.FitFromData, plot(pStandardEqualValueTrue ~ startValue, xlim=c(-1,1),ylim=c(0,1), pch=16, col="black", frame.plot=F,xlab="StartValue",ylab="Equal Valued p(Now) Points",las=1, type="p", main = "Response Bias"))

  dev.off()
  write.table(df.FitFromData, file = "FitsFromData.txt", quote = F, sep = "\t", row.names = F)
  write.table(df.fit, file = "equalValueProbesFromData.txt", quote = F, sep = "\t", row.names = F)
