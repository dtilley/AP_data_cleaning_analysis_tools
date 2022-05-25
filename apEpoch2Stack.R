## description 
source("~/projects/code-share-112018/jgpplotter.R")

APepoch2stack <- function(epoch.list, APVsAPDs.list) {
    
    epochs.aps.list <- list()

    for (i in 1:length(epoch.list)) {
        epoch <- epoch.list[[i]]
        apvs.apds <- APVsAPDs.list[[i]]
        
        aps.list <- list()
        for (j in 1:nrow(apvs.apds)) {
            t <- epoch[,1][apvs.apds$first.ndx[j]:apvs.apds$last.ndx[j]]
            mV <- epoch[,2][apvs.apds$first.ndx[j]:apvs.apds$last.ndx[j]]
            dVdt.max <- which.max(diff(mV)/diff(t))
            t <- t-t[dVdt.max]
            aps.list[[j]] <- as.data.frame(cbind(t,mV))
        }

        epochs.aps.list[[i]] <- aps.list
    }

    return(epochs.aps.list)
}

plot.APepoch2stack <- function(stack.list, cols, ltys, ...) {
    tmp.ap <- stack.list[[1]][[1]]
    jgpplotter(tmp.ap[,1], tmp.ap[,2], type="n", ...)
    for (i in 1:length(stack.list)) {
        ap_stack <- stack.list[[i]]
        for (j in 1:length(ap_stack)) {
            lines(ap_stack[[j]][,1], ap_stack[[j]][,2], lty=ltys[i], col=cols[i])
        }
    }
}
