# Start point bias plots
    # Set overlap values
    overlaps <- seq(0,1.5,0.025)
    # Run RRW based on overlaps and bias
    # Unbiased start point
    unbiased <- getMomentsOfRRWoverlap(overlaps,boundary=10, startValue=0, loops = 500,noiseSD = 2)
    #Positive start point bias
    posbias1 <- getMomentsOfRRWoverlap(overlaps,boundary=10, startValue=0.25, loops = 500,noiseSD = 2)
    posbias2 <- getMomentsOfRRWoverlap(overlaps,boundary=10, startValue=0.50, loops = 500,noiseSD = 2)
    # Negative start point bias
    negbias1 <- getMomentsOfRRWoverlap(overlaps,boundary=10, startValue= -0.25, loops = 500,noiseSD = 2)
    negbias2 <- getMomentsOfRRWoverlap(overlaps,boundary=10, startValue= -0.50, loops = 500,noiseSD = 2)

    # Generate plot found in paper
     with(unbiased[unbiased$correct==T,],plot(pCross ~ overlap, pch=16, col="black",frame.plot=F,xlab="Overlap",ylab="pHO",ylim=c(0,1),las=1,cex=0.75,xaxt='n'))
       abline(0.5,0,lty=2)
       abline(v=1,lty=2)
       axis(1, at = seq(0, 1.5, by = .1), las=2)
       with(posbias1[posbias1$correct==T,],points(pCross ~ overlap, pch=16, col="sky blue",cex=0.75))
       with(posbias2[posbias2$correct==T,],points(pCross ~ overlap, pch=16, col="blue",cex=0.75))
       with(negbias1[negbias1$correct==T,],points(pCross ~ overlap, pch=16, col="light pink",cex=0.75))
       with(negbias2[negbias2$correct==T,],points(pCross ~ overlap, pch=16, col="red",cex=0.75))
       legend(0.2,0.4, legend=c("0.50","0.25","0 (unbiased)","-0.25","-0.50"),title="Start Point",col=c("blue", "skyblue","black","light pink","red"),pch=16,cex=0.75)
       # Write to a pdf
       dev.copy(pdf,"pHVOxOverlapVal.pdf")
       dev.off()
