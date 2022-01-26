source("~/projects/code-share-112018/jgpplotter.R")

plotIPSC.params <- function(params, cols=NULL, ...) {
    param.labels <- names(params)[1:14]
    y.labels  <- c(0.01, 0.1, 1, 10)
    y.mat <- as.matrix(params[,1:14])
    n <- nrow(y.mat)
    x.mat <- as.matrix(cbind(rep(1,n),rep(2,n),rep(3,n),rep(4,n),rep(5,n),
                             rep(6,n),rep(7,n),rep(8,n),rep(9,n),rep(10,n),
                             rep(11,n),rep(12,n),rep(13,n),rep(14,n)))
    if ( length(cols)==0 ) {
        jgpplotter(x=jitter(x.mat,factor=0.5),y=log(y.mat), axes=FALSE, xlab="", ylab="Scale Factor",
                   ...)
        axis(side=1, at=1:14, labels=param.labels, las=2)
        axis(side=2, at=log(y.labels), labels=y.labels)
    } else if ( length(cols) == nrow(params) ) {
        cols.mat <- NULL
        for (i in seq(length(cols))) {
            cols.mat <- rbind(cols.mat, rep(cols[i], 14))
        }
        cols.mat <- as.matrix(cols.mat)
        jgpplotter(x=jitter(x.mat,factor=0.5),y=log(y.mat), axes=FALSE, xlab="",
                   ylab="Scale Factor", col=cols.mat,  ... )
        axis(side=1, at=1:14, labels=param.labels, las=2)
        axis(side=2, at=log(y.labels), labels=y.labels)
    } else {
        jgpplotter(x=jitter(x.mat,factor=0.5),y=log(y.mat), axes=FALSE, xlab="", ylab="Scale Factor",
                   ...)
        axis(side=1, at=1:14, labels=param.labels, las=2)
        axis(side=2, at=log(y.labels), labels=y.labels)
    }
}

plotIPSC.paramsMSD <- function(means, sds, cols=NULL, ...) {
    param.labels <- names(means)[1:14]
    y.labels  <- c(0.01, 0.1, 1, 10)
    y.mat <- as.matrix(means[,1:14])
    y.lb <- y.mat - as.matrix(sds[,1:14])
    y.ub <- y.mat + as.matrix(sds[,1:14])
    n <- nrow(y.mat)
    x.mat <- as.matrix(cbind(rep(1,n),rep(2,n),rep(3,n),rep(4,n),rep(5,n),
                             rep(6,n),rep(7,n),rep(8,n),rep(9,n),rep(10,n),
                             rep(11,n),rep(12,n),rep(13,n),rep(14,n)))
    x.mat <- jitter(x.mat,factor=0.5)
    if ( length(cols)==0 ) {
        jgpplotter(x=x.mat,y=y.mat, axes=FALSE, xlab="", ylab="Scale Factor",
                   pch=20, ...)
        arrows(x0=x.mat, y0=y.lb, x1=x.mat, y1=y.ub, code=0, length=0)
        axis(side=1, at=1:14, labels=param.labels, las=2)
        axis(side=2, at=log(y.labels), labels=y.labels)
    } else if ( length(cols) == nrow(means) ) {
        cols.mat <- NULL
        for (i in seq(length(cols))) {
            cols.mat <- rbind(cols.mat, rep(cols[i], 14))
        }
        cols.mat <- as.matrix(cols.mat)
        jgpplotter(x=x.mat,y=y.mat, axes=FALSE, xlab="",
                   pch=20, ylab="Scale Factor", col=cols.mat,  ... )
        arrows(x0=x.mat, y0=y.lb, x1=x.mat, y1=y.ub, code=0, length=0, col=cols.mat)
        axis(side=1, at=1:14, labels=param.labels, las=2)
        axis(side=2, at=log(y.labels), labels=y.labels)
    } else {
        jgpplotter(x=x.mat,y=y.mat, axes=FALSE, xlab="", ylab="Scale Factor",
                   pch=20, ...)
        arrows(x0=x.mat, y0=y.lb, x1=x.mat, y1=y.ub, code=0, length=0)
        axis(side=1, at=1:14, labels=param.labels, las=2)
        axis(side=2, at=log(y.labels), labels=y.labels)
    }
}


num2col <- function(clstrs, cols) {
    tmp <- unique(clstrs)
    if (length(tmp)==length(cols)) {
        for (i in tmp) {
            clstrs[which(clstrs==i)] <- cols[i]
        }
    }
    return(clstrs)
}

apply.paramsFUN <- function(paramsLIST, FUN) {
    tmp <- NULL
    for (i in seq(length(paramsLIST))) {
        tmp <- rbind(tmp, sapply(paramsLIST[[i]], FUN=FUN))
    }
    return(tmp)
}
