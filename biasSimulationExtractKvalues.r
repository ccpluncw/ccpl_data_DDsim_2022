library(dplyr)
library(pracma)
library(chutils)

kFit <- function(stand, delay, k) {
  out <-stand/(1+(k*delay))
  return(out)
}

# Create rounding function
  round_any <-  function(x, accuracy, f=round){f(x/ accuracy) * accuracy}
  round.val <- 0.01

#set your simulated k value here
  delayDat <- read.table("delayDat.txt", sep = "\t", header = T)
  startValue <- sort(unique(delayDat$startVal))
  delayDat$lDelay <- round(log(delayDat$Delay),2)
  delayDat$pStandard <- ifelse(delayDat$choice  == "now", delayDat$pCross, 1-delayDat$pCross)
  delayDat$RT <- scale(delayDat$mean)


  ############# Here, for each sv, we plot the p(now) x RT and p(delay) x RT. then we fit a polynomial to
  ######  each set of data.  This shows the relation of p(now) to p(delay) and how that relates to RT.


  maxY <- 1.1*max(delayDat$RT, na.rm=T)
  minY <- .9*min(delayDat$RT, na.rm=T)
  maxY <- round_any(maxY, 1, ceiling)
  minY <- round_any(minY, 1, floor)

  maxX.ld <- 1.1*max(delayDat$lDelay, na.rm=T)
  minX.ld <- .9*min(delayDat$lDelay, na.rm=T)
  maxX.ld <- round_any(maxX.ld, 1, ceiling)
  minX.ld <- round_any(minX.ld, 1, floor)


  ### Get x y pair for each probe and plot by delay (value on y and delay on x), fit k with nonlinear least squares
  ### Use approx to find k values
          ## Unbiased
    #these are true k values from simulated data
    ApproxVal.ov <- NULL
    #these are estimated k values using traditional method
    ApproxVal <- NULL
    biases <- c("unbiased", "delayed", "immediacy")
      for(sv in startValue) {
          tmp.dat <- delayDat[delayDat$startVal == sv, ]

          ApproxVal.ov.1 <- NULL
          ApproxVal.1 <- NULL
          for(a in unique(tmp.dat$Probe)){
            #here, we are interpolating points to find the indifferent point
            tmp.dat.n <- tmp.dat[tmp.dat$Probe == a & tmp.dat$choice == "now",]
            tmp.dat.d <- tmp.dat[tmp.dat$Probe == a & tmp.dat$choice == "delay",]
            tmp.dat.n <- tmp.dat.n[order(tmp.dat.n$Delay),]
            tmp.dat.d <- tmp.dat.d[order(tmp.dat.d$Delay),]
            x.n <- tmp.dat.n$Delay
            x.d <- tmp.dat.d$Delay
            y <- tmp.dat.n$pCross
            pVal <- tmp.dat.n$ProbeVal
            overlp <- tmp.dat.n$Overlap
            rt.n <- tmp.dat.n$RT
            rt.d <- tmp.dat.d$RT

            approx_val.pc <- approx(x.n,y,  n=5000)
            approx_val.pv <- approx(x.n,pVal,  approx_val.pc$x)
            approx_val.ovlp <- approx(x.n,overlp,  approx_val.pc$x)
            approx_val.rt.n <- approx(x.n,rt.n,  approx_val.pc$x)
            approx_val.rt.d <- approx(x.d,rt.d,  approx_val.pc$x)

            indif.df.pc <- data.frame(Delay = approx_val.pc$x,pCross = approx_val.pc$y)
            indif.df.pv <- data.frame(Delay = approx_val.pv$x,ProbeVal = approx_val.pv$y)
            indif.df.ovlp <- data.frame(Delay = approx_val.ovlp$x,overlap = approx_val.ovlp$y)
            indif.df.rt.n <- data.frame(Delay = approx_val.rt.n$x,RT.now = approx_val.rt.n$y)
            indif.df.rt.d <- data.frame(Delay = approx_val.rt.d$x,RT.delay = approx_val.rt.d$y)
            indif.df.1 <- merge (indif.df.pc, indif.df.rt.n)
            indif.df.2 <- merge (indif.df.1, indif.df.rt.d)
            indif.df.3 <- merge (indif.df.2, indif.df.ovlp)
            indif.df <- merge (indif.df.3, indif.df.pv)
            indif.df$DolVal <- a
            indif.df$bias <- unique(tmp.dat$bias)
            indif.df$startVal <- sv

            #indifference points
            df.tmp <- indif.df[which.min(abs(indif.df$pCross - 0.50)),]
            if (is.null(ApproxVal.1)) {
              ApproxVal.1 <- df.tmp
            } else {
              ApproxVal.1 <- rbind(ApproxVal.1, df.tmp)
            }

            #true equal value points
            df.tmp.ov <- indif.df[which.min(abs(indif.df$overlap - 1)),]
            if (is.null(ApproxVal.ov.1)) {
              ApproxVal.ov.1 <- df.tmp.ov
            } else {
              ApproxVal.ov.1 <- rbind(ApproxVal.ov.1, df.tmp.ov)
            }
          }

          #indifference points
          ApproxVal.2 <- ApproxVal.1[ApproxVal.1$pCross < .53 & ApproxVal.1$pCross > 0.47,]
          nls.sum.tmp1 <- nls(formula=DolVal ~ 100/(1+(k*Delay)),data=ApproxVal.2,start=list(k=0.1))
          nls.sum.unbiased <- summary(nls.sum.tmp1)
          #k value
          ApproxVal.2$k <- nls.sum.unbiased$coef[1]
          ApproxVal.2$Fit <- fitted(nls.sum.tmp1)

          if (is.null(ApproxVal)) {
            ApproxVal <- ApproxVal.2
          } else {
            ApproxVal <- rbind(ApproxVal, ApproxVal.2)
          }

          #true equal value points
          ApproxVal.ov.2 <- ApproxVal.ov.1[ApproxVal.ov.1$overlap < 1.03 & ApproxVal.ov.1$overlap > 0.97,]
          nls.sum.tmp1 <- nls(formula=DolVal ~ 100/(1+(k*Delay)),data=ApproxVal.ov.2,start=list(k=0.1))
          nls.sum.unbiased <- summary(nls.sum.tmp1)
          #k value
          ApproxVal.ov.2$k <- nls.sum.unbiased$coef[1]
          ApproxVal.ov.2$Fit <- fitted(nls.sum.tmp1)
          if (is.null(ApproxVal.ov)) {
            ApproxVal.ov <- ApproxVal.ov.2
          } else {
            ApproxVal.ov <- rbind(ApproxVal.ov, ApproxVal.ov.2)
          }
      }

    ApproxVal$rt.diff <- ApproxVal$RT.now - ApproxVal$RT.delay
    ApproxVal.ov$rt.diff <- ApproxVal.ov$RT.now - ApproxVal.ov$RT.delay

  #plot how RTnow vs RTdelay for the indifferent points
    maxY <- 1.1*max(ApproxVal$RT.delay,ApproxVal$RT.now, na.rm=T)
    minY <- .9*min(ApproxVal$RT.delay,ApproxVal$RT.now, na.rm=T)
    maxY <- round_any(maxY, 1, ceiling)
    minY <- round_any(minY, 1, floor)
    maxY.d <- 1.1*max(ApproxVal$rt.diff,ApproxVal$rt.diff, na.rm=T)
    minY.d <- .9*min(ApproxVal$rt.diff,ApproxVal$rt.diff, na.rm=T)
    maxY.d <- round_any(maxY.d, 1, ceiling)
    minY.d <- round_any(minY.d, 1, floor)

    fileName <- "indiff RT x Bias x StartValue.pdf"
    pdf(file=fileName)
    for(sv in startValue) {
        df.tmp <- ApproxVal[ApproxVal$startVal == sv,]
        bias <- unique(ApproxVal$bias)
        tmp.lm <- with(df.tmp, lm(RT.delay ~ Delay))
        with(df.tmp, plot(RT.delay ~ Delay, ylim = c(minY,maxY),xlim=c(0,300), pch=16, col="black", frame.plot=F,xlab="Delay",ylab="RT",las=1, main = paste(bias[1], sv[1])))
        abline(tmp.lm, col="black")
        tmp.lm <- with(df.tmp, lm(RT.now ~ Delay))
        with(df.tmp, points(RT.now ~ Delay, pch=16, col="grey"))
        abline(tmp.lm, col="grey")
        with(df.tmp, boxplot(RT.now, RT.delay, main = paste(bias[1], sv[1]),frame.plot=F,names = c("Now","Delayed") , ylim = c(minY,maxY)))
        with(df.tmp, boxplot(rt.diff, main = paste(bias[1], sv[1]),frame.plot=F,names = c("Now -Delayed") , ylim = c(minY.d,maxY.d)))
        abline(h=0, col="grey")

    }
    with(ApproxVal, boxplot(rt.diff~startVal,frame.plot=F, ylab="rt.now - rt.delay", ylim = c(minY.d,maxY.d)))
    abline(h=0, col="grey")
    dev.off()

    #plot how RTnow vs RTdelay for the true equal values points
    maxY <- 1.1*max(ApproxVal.ov$RT.delay,ApproxVal.ov$RT.now, na.rm=T)
    minY <- .9*min(ApproxVal.ov$RT.delay,ApproxVal.ov$RT.now, na.rm=T)
    maxY <- round_any(maxY, 1, ceiling)
    minY <- round_any(minY, 1, floor)
    maxY.d <- 1.1*max(ApproxVal.ov$rt.diff,ApproxVal.ov$rt.diff, na.rm=T)
    minY.d <- .9*min(ApproxVal.ov$rt.diff,ApproxVal.ov$rt.diff, na.rm=T)
    maxY.d <- round_any(maxY.d, 1, ceiling)
    minY.d <- round_any(minY.d, 1, floor)

    fileName <- "true equal RT x Bias x StartValue.pdf"
    pdf(file=fileName)
    for(sv in startValue) {
        df.tmp <- ApproxVal.ov[ApproxVal.ov$startVal == sv,]
        bias <- unique(ApproxVal.ov$bias)
        tmp.lm <- with(df.tmp, lm(RT.delay ~ Delay))
        with(df.tmp, plot(RT.delay ~ Delay, ylim = c(minY, maxY),xlim=c(0,300), pch=16, col="black", frame.plot=F,xlab="Delay",ylab="RT",las=1, main = paste(bias[1], sv[1])))
        abline(tmp.lm, col="black")
        tmp.lm <- with(df.tmp, lm(RT.now ~ Delay))
        with(df.tmp, points(RT.now ~ Delay, pch=16, col="grey"))
        abline(tmp.lm, col="grey")
        with(df.tmp, boxplot(RT.now, RT.delay, main = paste(bias[1], sv[1]),frame.plot=F,names = c("Now","Delayed") , ylim = c(minY, maxY)))
        with(df.tmp, boxplot(rt.diff, main = paste(bias[1], sv[1]),frame.plot=F,names = c("Now -Delayed") , ylim = c(minY.d, maxY.d)))
        abline(h=0, col="grey")

    }
    with(ApproxVal.ov, boxplot(rt.diff~startVal,frame.plot=F, ylab="rt.now - rt.delay", ylim = c(minY.d, maxY.d)))
    abline(h=0, col="grey")
    dev.off()
###points

    pdf(file="kValues x Bias Indifference points.pdf")
    startVs <- unique(abs(startValue[startValue>0]))
    for(sv in startVs) {
        k.df <- ApproxVal[ApproxVal$startVal == 0,]
        xHigh <- 1.1*max(k.df$Delay)
        k.tmp <- k.df[k.df$bias == "unbiased",]
        with(k.tmp,plot(DolVal ~ Delay ,pch=16,frame.plot=F, xlim=c(0,300),ylim=c(0,100),xlab="Delay",ylab="Dollars",las=1, main = sv[1]))
        with(k.tmp,lines(c(1:xHigh), kFit(100,c(1:xHigh), k.tmp$k[1]) ))

        k.tmp <- ApproxVal[ApproxVal$startVal == sv,]
        with(k.tmp,points(DolVal ~ Delay, pch=16,col="red"))
        with(k.tmp,lines(c(1:xHigh), kFit(100,c(1:xHigh), k.tmp$k[1]), col="red") )

        k.tmp <- ApproxVal[ApproxVal$startVal == (sv*-1),]
        with(k.tmp,points(DolVal ~ Delay,pch=16,col="blue"))
        with(k.tmp,lines(c(1:xHigh), kFit(100,c(1:xHigh), k.tmp$k[1]),col="blue" ))
        legend(200,80, legend=c("Bias Towards Immediacy","Unbiased","Bias Towards Delay"),title="Start Point Bias",col=c("red","black","blue"),pch=16,cex=0.75)

    }
    dev.off()


    pdf(file="kValues x Bias True Equal Value.pdf")
    startVs <- unique(abs(startValue[startValue>0]))
    for(sv in startVs) {
        k.df <- ApproxVal.ov[ApproxVal.ov$startVal == 0,]
        xHigh <- 1.1*max(k.df$Delay)
        k.tmp <- k.df[k.df$bias == "unbiased",]
        with(k.tmp,plot(DolVal ~ Delay ,pch=16,frame.plot=F, xlim=c(0,300),ylim=c(0,100),xlab="Delay",ylab="Dollars",las=1, main = sv[1]))
        with(k.tmp,lines(c(1:xHigh), kFit(100,c(1:xHigh), k.tmp$k[1]) ))

        k.tmp <- ApproxVal.ov[ApproxVal.ov$startVal == sv,]
        with(k.tmp,points(DolVal ~ Delay, pch=16,col="red"))
        with(k.tmp,lines(c(1:xHigh), kFit(100,c(1:xHigh), k.tmp$k[1]), col="red") )

        k.tmp <- ApproxVal.ov[ApproxVal.ov$startVal == (sv*-1),]
        with(k.tmp,points(DolVal ~ Delay,pch=16,col="blue"))
        with(k.tmp,lines(c(1:xHigh), kFit(100,c(1:xHigh), k.tmp$k[1]),col="blue" ))
        legend(200,80, legend=c("Bias Towards Immediacy","Unbiased","Bias Towards Delay"),title="Start Point Bias",col=c("red","black","blue"),pch=16,cex=0.75)

    }
    dev.off()

    write.table(ApproxVal, file = "kvalues Traditional Estimated.txt", quote = F, sep = "\t", row.names = F)
    write.table(ApproxVal.ov, file = "kvalues Actual from Simulation.txt", quote = F, sep = "\t", row.names = F)
