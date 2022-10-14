library(dplyr)
library(pracma)
library(chutils)


# Create rounding function
  round_any <-  function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
  round.val <- 0.01

  delayDat <- read.table("delayDat.txt", sep = "\t", header = T)
  startValue <- sort(unique(delayDat$startVal))
  delayDat$lDelay <- round(log(delayDat$Delay),2)
  delayDat$pStandard <- ifelse(delayDat$choice  == "now", delayDat$pCross, 1-delayDat$pCross)
  delayDat$RT <- scale(delayDat$mean)

  biases <- c("unbiased", "delayed", "immediacy")

  maxY <- 1.1*max(delayDat$RT, na.rm=T)
  minY <- .9*min(delayDat$RT, na.rm=T)
  maxY <- round_any(maxY, 1, ceiling)
  minY <- round_any(minY, 1, floor)

  maxX.ld <- 1.1*max(delayDat$lDelay, na.rm=T)
  minX.ld <- .9*min(delayDat$lDelay, na.rm=T)
  maxX.ld <- round_any(maxX.ld, 1, ceiling)
  minX.ld <- round_any(minX.ld, 1, floor)


  ################### RT by delay #############
  biases <- c("unbiased", "delayed", "immediacy")
  pdf(file="RT x Delay.pdf")
  df.eqval <- NULL
  df.comb.diff.n <- NULL
  df.comb.diff.d <- NULL
  for(sv in startValue) {
      tmp.dat <- delayDat[delayDat$startVal == sv, ]
      for(prb in sort(unique(tmp.dat$Probe))) {
        bias <- unique(tmp.dat$bias)
        title <- paste(bias,sv[1],prb[1])
        df.tmp <- tmp.dat[tmp.dat$Probe == prb,]
        tmp.df.1 <- df.tmp[df.tmp$choice == "delay",]
        tmp.df.1 <- tmp.df.1[order(tmp.df.1$lDelay),]
        x.tmp <- diff(tmp.df.1$RT)
        tmp.df.1$RTdiff <- append(x.tmp, NA, after = 0)
        with(tmp.df.1, plot(RT ~ lDelay,ylim=c(minY, maxY), pch=16, col="black", frame.plot=F,xlab="Delay",ylab="RT",las=1, type="l", main = title))

        tmp.df.2 <- df.tmp[df.tmp$choice == "now",]
        tmp.df.2 <- tmp.df.2[order(tmp.df.2$lDelay),]
        x.tmp <- diff(tmp.df.2$RT)
        tmp.df.2$RTdiff <- append(x.tmp, NA, after = 0)
        with(tmp.df.2, lines(RT ~ lDelay, pch=16, col="grey"))

        eq.val <- df.tmp[which.min(abs(df.tmp$Overlap - 1)),]
        abline(v=eq.val$lDelay[1], lty=3, col="black")


        maxY.t <- 1.1*max(tmp.df.2$RTdiff,tmp.df.1$RTdiff, na.rm=T)
        minY.t <- .9*min(tmp.df.2$RTdiff,tmp.df.1$RTdiff, na.rm=T)
        maxY.t <- round_any(maxY.t, 1, ceiling)
        minY.t <- round_any(minY.t, 1, floor)
        fit.now <- NULL
        if(nrow(tmp.df.1[!is.na(tmp.df.1$RTdiff),]) > 2) {
          tryCatch(
              expr = {
                  fit.now <- loess(RTdiff~lDelay, data=tmp.df.1)
              },
              error = function(e){
                }
          )
        }
        with(tmp.df.1, plot(RTdiff ~ lDelay,ylim=c(minY.t,maxY.t), pch=16, col="black", frame.plot=F,xlab="Delay",ylab="RTdiff",las=1, type="p", main = paste("DELAY RTdiff", title)))
        if(!is.null(fit.now)) {
          with(tmp.df.1, lines(lDelay, predict(fit.now,tmp.df.1), col='black'))
          tmp.df.1$fit <- predict(fit.now,tmp.df.1)
        } else {
          tmp.df.1$fit <- NA
        }
        abline(v=eq.val$lDelay[1], lty=3, col="black")
        abline(h=0, lty=3, col="black")

        fit.now <- NULL
        if(nrow(tmp.df.2[!is.na(tmp.df.2$RTdiff),]) > 2) {
          tryCatch(
              expr = {
                fit.now <- loess(RTdiff~lDelay, data=tmp.df.2)
              },
              error = function(e){
                }
          )
        }
        with(tmp.df.2, plot(RTdiff ~ lDelay,ylim=c(minY.t,maxY.t), pch=16, col="grey", frame.plot=F,xlab="Delay",ylab="RTdiff",las=1, type="p", main = paste("NOW RTdiff", title)))
        if(!is.null(fit.now)) {
          with(tmp.df.2, lines(lDelay, predict(fit.now,tmp.df.2), col='grey'))
          tmp.df.2$fit <- predict(fit.now,tmp.df.2)
        } else {
          tmp.df.2$fit <- NA
        }
        abline(v=eq.val$lDelay[1], lty=3, col="black")
        abline(h=0, lty=3, col="black")

        tmp.df.2$eqVal <- eq.val$lDelay[1]
        tmp.df.2$delayDist <- tmp.df.2$lDelay - tmp.df.2$eqVal

        tmp.df.1$eqVal <- eq.val$lDelay[1]
        tmp.df.1$delayDist <- tmp.df.1$lDelay - tmp.df.1$eqVal

        if (is.null(df.eqval)) {
          df.eqval <- eq.val
        } else {
          df.eqval <- rbind(df.eqval, eq.val)
        }

        if (is.null(df.comb.diff.n)) {
          df.comb.diff.n <- tmp.df.2
        } else {
          df.comb.diff.n <- rbind(df.comb.diff.n, tmp.df.2)
        }

        if (is.null(df.comb.diff.d)) {
          df.comb.diff.d <- tmp.df.1
        } else {
          df.comb.diff.d <- rbind(df.comb.diff.d, tmp.df.1)
        }

      }
  }

  time1 <- c("Now", "Delay")
  for(tm in time1) {
    if(tm == "Now") {
      df.comb.diff <- df.comb.diff.n
    } else {
      df.comb.diff <- df.comb.diff.d
    }
    maxY <- 1.1*max(df.comb.diff$RTdiff, na.rm=T)
    minY <- .9*min(df.comb.diff$RTdiff, na.rm=T)
    maxY <- round_any(maxY, 1, ceiling)
    minY <- round_any(minY, 1, floor)
    maxX <- 1.1*max(df.comb.diff$delayDist, na.rm=T)
    minX <- .9*min(df.comb.diff$delayDist, na.rm=T)
    maxX <- round_any(maxX, 1, ceiling)
    minX <- round_any(minX, 1, floor)
    with(df.comb.diff, plot(NULL, NULL,ylim=c(minY,maxY), xlim=c(minX,maxX),pch=16, col="grey", frame.plot=F,xlab="log(Delay) - log(Delay Equal Value Point) ",ylab="Delay_n_RT - Delay_n-1_RT",las=1, type="p", main = paste(tm)))
    for(sv in startValue) {
        for(prb in sort(unique(tmp.dat$Probe))) {
          df.tmp <- df.comb.diff[df.comb.diff$Probe == prb & df.comb.diff$startVal == sv,]
          with(df.tmp, lines(delayDist, fit, col='grey'))
        }
    }
    abline(v=0, lty=3, col="black")
    abline(h=0, lty=3, col="black")
  }

  dev.off()
  write.table(df.eqval, file = "eqValues.txt", quote = F, sep = "\t", row.names = F)
