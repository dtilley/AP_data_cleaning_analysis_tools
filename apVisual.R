# Generic plot fucntion with figure publication standards from the JGP
# EPS vector graphics work best for cross platform compatibility with AI and Inkscape.
jgpplotter <- function(x,y,mar=c(5.1,5.1,5.1,2.1),lwdpt=0.75,...){
    # convert pt to lwd
    lwd <- lwdpt*(1/0.75)
    par(mar=mar,lwd=lwd)
    plot(x,y,bty="n",...)
}


visualizeAPDs <- function(t,mV,apVoltages,apDurations,...){
    jgpplotter(t,mV,type="l",...)
    abline(v=t[apVoltages$peakV.ndx],col="blue")
    abline(v=t[apDurations$first.ndx],lty=1,col="seagreen")
    abline(v=t[apDurations$last.ndx],lty=2,col="magenta")
    colors <- rainbow(n=9)
    for (i in 1:nrow(apDurations)){
        segments(x0=apDurations$apd10.t1[i],y0=apDurations$ap10[i],x1=apDurations$apd10.t2[i],y1=apDurations$ap10[i],col=colors[1])
        segments(x0=apDurations$apd20.t1[i],y0=apDurations$ap20[i],x1=apDurations$apd20.t2[i],y1=apDurations$ap20[i],col=colors[2])
        segments(x0=apDurations$apd30.t1[i],y0=apDurations$ap30[i],x1=apDurations$apd30.t2[i],y1=apDurations$ap30[i],col=colors[3])
        segments(x0=apDurations$apd40.t1[i],y0=apDurations$ap40[i],x1=apDurations$apd40.t2[i],y1=apDurations$ap40[i],col=colors[4])
        segments(x0=apDurations$apd50.t1[i],y0=apDurations$ap50[i],x1=apDurations$apd50.t2[i],y1=apDurations$ap50[i],col=colors[5])
        segments(x0=apDurations$apd60.t1[i],y0=apDurations$ap60[i],x1=apDurations$apd60.t2[i],y1=apDurations$ap60[i],col=colors[6])
        segments(x0=apDurations$apd70.t1[i],y0=apDurations$ap70[i],x1=apDurations$apd70.t2[i],y1=apDurations$ap70[i],col=colors[7])
        segments(x0=apDurations$apd80.t1[i],y0=apDurations$ap80[i],x1=apDurations$apd80.t2[i],y1=apDurations$ap80[i],col=colors[8])
        segments(x0=apDurations$apd90.t1[i],y0=apDurations$ap90[i],x1=apDurations$apd90.t2[i],y1=apDurations$ap90[i],col=colors[9])
    }

}


visualizeAPDs.SAP <- function(t, mV, APVsAPDs, ...) {
    jgpplotter(t, mV, type="l", ...)
    abline(v=t[APVsAPDs$peakV.ndx],col="blue")
    abline(v=t[APVsAPDs$first.ndx],lty=1,col="seagreen")
    abline(v=t[APVsAPDs$last.ndx],lty=2,col="magenta")
    colors <- rainbow(n=9)

    segments(x0=APVsAPDs$apd10.t1, y0=APVsAPDs$ap10, x1=APVsAPDs$apd10.t2,
             y1=APVsAPDs$ap10, col=colors[1])
    segments(x0=APVsAPDs$apd20.t1, y0=APVsAPDs$ap20, x1=APVsAPDs$apd20.t2,
             y1=APVsAPDs$ap20, col=colors[2])
    segments(x0=APVsAPDs$apd30.t1, y0=APVsAPDs$ap30, x1=APVsAPDs$apd30.t2,
             y1=APVsAPDs$ap30, col=colors[3])
    segments(x0=APVsAPDs$apd40.t1, y0=APVsAPDs$ap40, x1=APVsAPDs$apd40.t2,
             y1=APVsAPDs$ap40, col=colors[4])
    segments(x0=APVsAPDs$apd50.t1, y0=APVsAPDs$ap50, x1=APVsAPDs$apd50.t2,
             y1=APVsAPDs$ap50, col=colors[5])
    segments(x0=APVsAPDs$apd60.t1, y0=APVsAPDs$ap60, x1=APVsAPDs$apd60.t2,
             y1=APVsAPDs$ap60, col=colors[6])
    segments(x0=APVsAPDs$apd70.t1, y0=APVsAPDs$ap70, x1=APVsAPDs$apd70.t2,
             y1=APVsAPDs$ap70, col=colors[7])
    segments(x0=APVsAPDs$apd80.t1, y0=APVsAPDs$ap80, x1=APVsAPDs$apd80.t2,
             y1=APVsAPDs$ap80, col=colors[8])
    segments(x0=APVsAPDs$apd90.t1, y0=APVsAPDs$ap90, x1=APVsAPDs$apd90.t2,
             y1=APVsAPDs$ap90, col=colors[9])

}
