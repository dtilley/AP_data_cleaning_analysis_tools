
apvs <- function(t, mV) {
    ## Calculates the voltages in an AP train 
    ## dt: period between samples and defaults to 0.1 ms
    ## ndxDF: a dataframe contraining indices from apfind 
    
    ## Default values for AP properties
    bri <- NA
    mdp <- NA 
    mdp.ndx <- NA 
    peakV <- NA 
    peakV.ndx <- NA

    ## Calculates peakV
    peakV <- max(mV)
    peakV.ndx <- which(mV == peakV)
    ## Calculates the minium between the first AP and the following AP
    mdp <- min(mV)
    mdp.ndx <- which(mV == mdp)

    bri <- t[length(t)] - t[1]
    rtndf <- as.data.frame(cbind(bri,mdp,mdp.ndx,peakV,peakV.ndx))
    return(rtndf)
}


apds <- function(t, mV, preThreshold.dt=200) {
    ## t: time vector 
    ## mV: membrane voltage vector
    ## preThreshold.dt: the time prior to the peak voltage to calculate diastolic membrane voltage 
    ## apVoltages: a dataframe of indices from apVoltages
    apVoltages <- apvs(t, mV)

    preThreshold.V <- NA
    
    first.ndx <- NA
    last.ndx <- NA
    amp <- NA

    ap10 <- NA
    ap20 <- NA
    ap30 <- NA
    ap40 <- NA
    ap50 <- NA
    ap60 <- NA
    ap70 <- NA
    ap80 <- NA
    ap90 <- NA

    apd10.t1 <- NA
    apd20.t1 <- NA
    apd30.t1 <- NA
    apd40.t1 <- NA
    apd50.t1 <- NA
    apd60.t1 <- NA
    apd70.t1 <- NA
    apd80.t1 <- NA
    apd90.t1 <- NA

    apd10.t2 <- NA
    apd20.t2 <- NA
    apd30.t2 <- NA
    apd40.t2 <- NA
    apd50.t2 <- NA
    apd60.t2 <- NA
    apd70.t2 <- NA
    apd80.t2 <- NA
    apd90.t2 <- NA

    apd10 <- NA
    apd20 <- NA    
    apd30 <- NA
    apd40 <- NA
    apd50 <- NA
    apd60 <- NA
    apd70 <- NA
    apd80 <- NA    
    apd90 <- NA
    
    ## Calculates first.ndx
    time.firstndx <- t[apVoltages$peakV.ndx]-preThreshold.dt
    first.ndx <- which.min(abs(t-time.firstndx))

    ## Calculates last.ndx
    last.ndx <- length(t)

    ## Calculates APDs
    ## Calculates preThreshold.V over 20 indices forward from first.ndx        
    preThreshold.V <- mV[first.ndx]
    amp <- apVoltages$peakV-preThreshold.V
    ap10 <- preThreshold.V+0.9*amp
    ap20 <- preThreshold.V+0.8*amp
    ap30 <- preThreshold.V+0.7*amp
    ap40 <- preThreshold.V+0.6*amp
    ap50 <- preThreshold.V+0.5*amp
    ap60 <- preThreshold.V+0.4*amp
    ap70 <- preThreshold.V+0.3*amp
    ap80 <- preThreshold.V+0.2*amp
    ap90 <- preThreshold.V+0.1*amp
    
    ## Calculate apd10
    ## pre-peak points to interpolate
    p1.ndx <- min(which(mV[first.ndx:apVoltages$peakV.ndx]>=ap10))+first.ndx-1
    p2.ndx <- p1.ndx-1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t1 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap10)
    apd10.t1 <- t1$y
    ## post-peak points to interpolate
    p1.ndx <- max(which(mV[apVoltages$peakV.ndx:last.ndx]>=ap10))+apVoltages$peakV.ndx-1
    p2.ndx <- p1.ndx+1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t2 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap10)
    apd10.t2 <- t2$y
    apd10 <- apd10.t2-apd10.t1
    
    ## Calculate apd20
    ## pre-peak points to interpolate
    p1.ndx <- min(which(mV[first.ndx:apVoltages$peakV.ndx]>=ap20))+first.ndx-1
    p2.ndx <- p1.ndx-1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t1 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap20)
    apd20.t1 <- t1$y
    ## post-peak points to interpolate
    p1.ndx <- max(which(mV[apVoltages$peakV.ndx:last.ndx]>=ap20))+apVoltages$peakV.ndx-1
    p2.ndx <- p1.ndx+1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t2 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap20)
    apd20.t2 <- t2$y
    apd20 <- apd20.t2-apd20.t1
    
    ## Calculate apd30
    ## pre-peak points to interpolate
    p1.ndx <- min(which(mV[first.ndx:apVoltages$peakV.ndx]>=ap30))+first.ndx-1
    p2.ndx <- p1.ndx-1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t1 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap30)
    apd30.t1 <- t1$y
    ## post-peak points to interpolate
    p1.ndx <- max(which(mV[apVoltages$peakV.ndx:last.ndx]>=ap30))+apVoltages$peakV.ndx-1
    p2.ndx <- p1.ndx+1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t2 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap30)
    apd30.t2 <- t2$y
    apd30 <- apd30.t2-apd30.t1
    
    ## Calculate apd40
    ## pre-peak points to interpolate
    p1.ndx <- min(which(mV[first.ndx:apVoltages$peakV.ndx]>=ap40))+first.ndx-1
    p2.ndx <- p1.ndx-1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t1 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap40)
    apd40.t1 <- t1$y
    ## post-peak points to interpolate
    p1.ndx <- max(which(mV[apVoltages$peakV.ndx:last.ndx]>=ap40))+apVoltages$peakV.ndx-1
    p2.ndx <- p1.ndx+1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t2 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap40)
    apd40.t2 <- t2$y
    apd40 <- apd40.t2-apd40.t1
    
    ## Calculate apd50
    ## pre-peak points to interpolate
    p1.ndx <- min(which(mV[first.ndx:apVoltages$peakV.ndx]>=ap50))+first.ndx-1
    p2.ndx <- p1.ndx-1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t1 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap50)
    apd50.t1 <- t1$y
    ## post-peak points to interpolate
    p1.ndx <- max(which(mV[apVoltages$peakV.ndx:last.ndx]>=ap50))+apVoltages$peakV.ndx-1
    p2.ndx <- p1.ndx+1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t2 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap50)
    apd50.t2 <- t2$y
    apd50 <- apd50.t2-apd50.t1
    
    ## Calculate apd60
    ## pre-peak points to interpolate
    p1.ndx <- min(which(mV[first.ndx:apVoltages$peakV.ndx]>=ap60))+first.ndx-1
    p2.ndx <- p1.ndx-1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t1 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap60)
    apd60.t1 <- t1$y
    ## post-peak points to interpolate
    p1.ndx <- max(which(mV[apVoltages$peakV.ndx:last.ndx]>=ap60))+apVoltages$peakV.ndx-1
    p2.ndx <- p1.ndx+1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t2 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap60)
    apd60.t2 <- t2$y
    apd60 <- apd60.t2-apd60.t1

    ## Calculate apd70
    ## pre-peak points to interpolate
    p1.ndx <- min(which(mV[first.ndx:apVoltages$peakV.ndx]>=ap70))+first.ndx-1
    p2.ndx <- p1.ndx-1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t1 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap70)
    apd70.t1 <- t1$y
    ## post-peak points to interpolate
    p1.ndx <- max(which(mV[apVoltages$peakV.ndx:last.ndx]>=ap70))+apVoltages$peakV.ndx-1
    p2.ndx <- p1.ndx+1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t2 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap70)
    apd70.t2 <- t2$y
    apd70 <- apd70.t2-apd70.t1
    
    ## Calculate apd80
    ## pre-peak points to interpolate
    p1.ndx <- min(which(mV[first.ndx:apVoltages$peakV.ndx]>=ap80))+first.ndx-1
    p2.ndx <- p1.ndx-1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t1 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap80)
    apd80.t1 <- t1$y
    ## post-peak points to interpolate
    p1.ndx <- max(which(mV[apVoltages$peakV.ndx:last.ndx]>=ap80))+apVoltages$peakV.ndx-1
    p2.ndx <- p1.ndx+1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t2 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap80)
    apd80.t2 <- t2$y
    apd80 <- apd80.t2-apd80.t1
    
    ## Calculate apd90
    ## pre-peak points to interpolate
    p1.ndx <- min(which(mV[first.ndx:apVoltages$peakV.ndx]>=ap90))+first.ndx-1
    p2.ndx <- p1.ndx-1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t1 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap90)
    apd90.t1 <- t1$y
    ## post-peak points to interpolate
    p1.ndx <- max(which(mV[apVoltages$peakV.ndx:last.ndx]>=ap90))+apVoltages$peakV.ndx-1
    p2.ndx <- p1.ndx+1
    
    t.interpolate <- t[c(p1.ndx,p2.ndx)]
    mV.interpolate <- mV[c(p1.ndx,p2.ndx)]
    t2 <- approx(x=mV.interpolate,y=t.interpolate,xout=ap90)
    apd90.t2 <- t2$y
    apd90 <- apd90.t2-apd90.t1
    
    rtndf <- as.data.frame(cbind(preThreshold.V, amp, first.ndx, last.ndx, apd10,
                                 apd20, apd30, apd40, apd50, apd60, apd70, apd80,
                                 apd90, ap10, ap20, ap30, ap40, ap50, ap60, ap70,
                                 ap80, ap90, apd10, apd20, apd30, apd40,
                                 apd50,apd60,apd70,apd80,apd90,
                                 apd10.t1,apd20.t1,apd30.t1,
                                 apd40.t1,apd50.t1,apd60.t1,
                                 apd70.t1,apd80.t1,apd90.t1,
                                 apd10.t2,apd20.t2,apd30.t2,
                                 apd40.t2,apd50.t2,apd60.t2,
                                 apd70.t2,apd80.t2,apd90.t2))
    return(rtndf)
    
}
