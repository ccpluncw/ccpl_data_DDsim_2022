library(chutils)
library(chMorals)
library(pracma)
library(dplyr)

### can't do by subject because each subject only received one instance of each trial combination
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

findPlateau <- function (data, yCol, threshold = 0.05, lag = 2) {
  absoluteThreshold <- max(abs(data[[yCol]]))*threshold
  data$diff <- ifelse(as.numeric(rownames(data)) > lag, diff(data[[yCol]], lag), NA)
  data$belowThresh <- ifelse(abs(data$diff) < absoluteThreshold & !is.na(data$diff), TRUE, FALSE)
  data$cumSumTrue <- cumsum(data$belowThresh)

  if(max(data$cumSumTrue) > 0) {
    inflectIndx <- min(which(data$cumSumTrue == 1))
    i <- 2
    done <- FALSE
    while(done == FALSE) {
      #check to see if the inflection fails on next row
      if(inflectIndx + 1 > length(data$belowThresh)) {
        done <- TRUE
        inflectIndx <- NA
      } else {

        if (data$belowThresh[inflectIndx + 1] == FALSE) {
          #if so, check to see if there is a flattening later
          if(max(data$cumSumTrue) >= i) {
            #if so, make that the new inflection point
            inflectIndx <- min(which(data$cumSumTrue == i))
            i <- i + 1
          } else {
            #if there is no flattening later, then plateau fails
            inflectIndx <- NA
            done <- TRUE
          }
        #if the inflection doesn't fail later, then we are done
        } else {
          done <- TRUE
        }

      }
    }
  } else {
    inflectIndx <- NA
  }
  return(inflectIndx)
}

threshold <- 0.02
lag <- 2
intervals <- 1000
prettyRoundPrecision <- 0.25

originalDataFile <- "amasinoEtAl_behavior.csv"
replicationDataFile <- "amasinoEtAl_behavior_rep.csv"

sink("Amasino.out.txt", append = F)
  cat("\n ******** Amasino Data Analysis ******** \n\n")
sink(NULL)

pdfFileName <- paste("Amasino.pdf")
pdf(file=pdfFileName)

df.EVP.comp <- NULL
for (dataSet in 1:2) {
  if(dataSet == 1) {
    fileName <- originalDataFile
    tag <- "original"
  } else {
    fileName <- replicationDataFile
    tag <- "replication"
  }
  delayDat <- read.csv(fileName, header = F)
  colnames(delayDat) <- c("sn", "SSamount", "LLamount", "SStime", "LLtime", "choice", "RT", "side", "condition")
  delayDat$lDelay <- round(log(delayDat$LLtime),2)

  delayDat <- delayDat[complete.cases(delayDat),]
  #remove the single equal value trial
  delayDat <- delayDat[delayDat$SSamount != 10, ]
  delayDat$RT.res <- log(delayDat$RT)
  delayDat$chooseNow <- ifelse(delayDat$choice == 0, 1, 0)

  sink("Amasino.out.txt", append = T)
    cat("\n ******** ", tag, " ******** \n\n")
  sink(NULL)

    df.pnow <- data.frame(delayDat %>% group_by(SSamount, LLtime) %>% summarize (pNow = mean(chooseNow, na.rm=T)))
    df.pnow$choice <- 0
    df.pnow.1 <- df.pnow
    df.pnow.1$choice <- 1
    df.pnow.1$pNow <- 1 - df.pnow.1$pNow
    df.pnow <- rbind(df.pnow, df.pnow.1)

    df.rt <- data.frame(delayDat %>% group_by(SSamount, LLtime, choice) %>% summarize (mRT = mean(RT.res, na.rm=T), n = length(RT.res)))
 
    df.all <- merge(df.rt, df.pnow)

    tmp.df.d <- df.all[df.all$choice == 1, ]
    fit.delay.all <- loess(mRT~pNow, span = 2,data=tmp.df.d)
    tmp.df.d$Fit <- predict(fit.delay.all,tmp.df.d)

    tmp.df.n <- df.all[df.all$choice == 0, ]
    fit.now.all <- loess(mRT~pNow, span = 2,data=tmp.df.n)
    tmp.df.n$Fit <- predict(fit.now.all,tmp.df.n)

    #test for difference in RT at pStandard between .45 nd .55 with predicted fits
    now50 <- predict(fit.now.all,seq(.45,.55,.005))
    delay50 <- predict(fit.delay.all,seq(.45,.55,.005))
    t.out <- t.test(now50,delay50, paired = T)

    #if the pStandard at 50% is significantly different, then there is a bias, so do the following
    if(t.out$p.value < 0.1) {
      #select the now dataset if there is a bias toward delay (these are the higher RT datasets)
      if(t.out$statistic > 1) {
        df.eq.dat <- tmp.df.n
      } else {
        #select the delay dataset if there is a bias toward now (these are the higher RT datasets)
        df.eq.dat <- tmp.df.d
      }
    }

    #get equal valued data from predicted RT in the Fit function
    #Basically, for each probe, identify the highest fit RT, and that delay is the equal value item
    title <- paste("Amasino RT x Probe Dataset =", tag)
    maxY <- 1.1*max(df.eq.dat$Fit, na.rm=T)
    minY <- .9 * abs(min(df.eq.dat$Fit, na.rm=T)) * sign(min(df.eq.dat$Fit, na.rm=T))
    maxY <- ch.round_any(maxY, prettyRoundPrecision, ceiling, center = F)
    minY <- ch.round_any(minY, prettyRoundPrecision, floor, center = F)

    maxX.ld <- 1.1*max(df.all$LLtime, na.rm=T)
    minX.ld <- .9*min(df.all$LLtime, na.rm=T)
    maxX.ld <- ch.round_any(maxX.ld, prettyRoundPrecision, ceiling)
    minX.ld <- ch.round_any(minX.ld, prettyRoundPrecision, floor)

    with(df.all, plot(NULL,NULL, xlim=c(minX.ld,maxX.ld),ylim=c(minY,maxY), pch=16, col="black", frame.plot=F,xlab="log(Delay)",ylab="RT",las=1, type="p", main = title))
    i <- 100
    df.fit <- NULL

    for(prb in sort(unique(df.eq.dat$SSamount))) {
      df.tmp <- df.eq.dat[df.eq.dat$SSamount == prb,]
      df.tmp <- df.tmp[order(df.tmp$LLtime),]
			
      col1 = paste("grey", i, sep = "")
      with(df.tmp, points(Fit ~ LLtime, col = col1, type="l", main = prb))

      eqDelay <- getEqualValueDelay(df.tmp, "LLtime", "Fit", "pNow", unbiased = F, intervals= intervals)
      df.tmp.fit.max <- data.frame(Probe = prb, pNow = eqDelay$pStandard[1], Delay = eqDelay$Delay[1], Fit = eqDelay$Fit[1])

      if (is.null(df.fit)) {
        df.fit <- df.tmp.fit.max
      } else {
        df.fit <- rbind(df.fit, df.tmp.fit.max)
      }
      i <- i - 5
    }

    title <- paste("Amasino Dataset =", tag)
    maxY <- 1.1*max(c(tmp.df.d$mRT, tmp.df.n$mRT), na.rm=T)
    minY <- .9 * abs(min(c(tmp.df.d$mRT, tmp.df.n$mRT), na.rm=T)) * sign(min(c(tmp.df.d$mRT, tmp.df.n$mRT), na.rm=T))
    maxY <- ch.round_any(maxY, prettyRoundPrecision, ceiling, center = F)
    minY <- ch.round_any(minY, prettyRoundPrecision, floor, center = F)
		
    with(tmp.df.d, plot(mRT ~ pNow,xlim=c(0,1),ylim=c(minY,maxY), col="black", frame.plot=F,xlab="p(Now)",ylab="RT",las=1, type="p", main = title))
      with(tmp.df.d, points(pNow, predict(fit.delay.all,tmp.df.d), pch=16, col='black'))
      with(tmp.df.n, points(mRT ~ pNow, col="grey"))
      with(tmp.df.n, points(pNow, predict(fit.now.all,tmp.df.n), pch=16, col='grey'))

    df.fit <- df.fit[complete.cases(df.fit),]
  print(df.fit)
    inflectionIndx <- findPlateau(df.fit, "Fit", threshold = threshold, lag = lag)
    if(!is.na(inflectionIndx)) {
      equalValProbes <- df.fit[inflectionIndx:nrow(df.fit),]
      eq.val.dat <- mean(equalValProbes$pNow, na.rm = T)
      abline(v=eq.val.dat, lty=3, col="black")

      #Fit the hyperbolic function to start value
      nls.sum.tmp1 <- nls(formula=Probe ~ 10/(1+(k*Delay)),data=equalValProbes,start=list(k=0.1))
      sink("Amasino.out.txt", append = T)
        cat("\n ******** df.fit ******** \n\n")
        print(df.fit)
        cat("\n\n ******** Equal Value p(Now) point ******** \n\n")
        print(eq.val.dat)
        cat("\n\n ******** k Value Fit for  ******** \n\n")
        print(summary(nls.sum.tmp1))
        cat("\n k = ", coef(nls.sum.tmp1), "\n\n")
        cat("\n log(k) = ", log(coef(nls.sum.tmp1)), "\n\n")
      sink(NULL)
      equalValProbes$Fit <- fitted(nls.sum.tmp1)
      equalValProbes <- equalValProbes[order(equalValProbes$Delay),]

      title <- paste("Amasino k Dataset =", tag)
      with(equalValProbes, plot(Probe ~ Delay, xlim=c(0,400),ylim=c(0,10), pch=16, col="black", frame.plot=F,xlab="Delay",ylab="Probe",las=1, type="p", main = title))
      with(equalValProbes,lines(Delay,Fit))
    } else {
      print("inflection Index is NA")
      print(df.fit)
    }
		equalValProbes$dataset <- tag
		df.EVP.comp <- ch.rbind(df.EVP.comp, equalValProbes)
}

  title <- paste("Amasino k dataset comparison")
  with(df.EVP.comp[df.EVP.comp$dataset == "replication",], plot(Probe ~ Delay, xlim=c(0,400),ylim=c(0,10), pch=4, col="grey", frame.plot=F,xlab="Delay",ylab="Probe",las=1, type="p", main = title))
  with(df.EVP.comp[df.EVP.comp$dataset == "replication",],lines(Delay,Fit, col="grey"))
  with(df.EVP.comp[df.EVP.comp$dataset == "original",], points(Probe ~ Delay,  col="black"))
  with(df.EVP.comp[df.EVP.comp$dataset == "original",],lines(Delay,Fit, col="black"))


dev.off()
