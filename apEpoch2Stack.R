source("~/projects/code-share-112018/apVisual.R")

cntrl_tmV <- read.table("./cell_1_APVsAPDs_tmV/cntrl_ishi_122420_3_drew.txt", header=T)
nICaL_tmV <- read.table("./cell_1_APVsAPDs_tmV/last15s_ical_decrease_122420_3_drew.txt", header=T)
pICaL_tmV <- read.table("./cell_1_APVsAPDs_tmV/t30-45s_ical_increase_122420_3_drew.txt", header=T)
nIto_tmV <- read.table("./cell_1_APVsAPDs_tmV/last15s_ito_decrease_122420_3_drew.txt", header=T)
pIto_tmV <- read.table("./cell_1_APVsAPDs_tmV/last15s_ito_increase_122420_3_drew.txt", header=T)

cntrl_APVsAPDs <- read.table("./cell_1_APVsAPDs_tmV/APDsAPVs_cntrl_ishi_122420_3_drew.txt", header=T)
nICaL_APVsAPDs <- read.table("./cell_1_APVsAPDs_tmV/APDsAPVs_ical_decrease_122420_3_drew.txt", header=T)
pICaL_APVsAPDs <- read.table("./cell_1_APVsAPDs_tmV/APDsAPVs_ical_increase_122420_3_drew.txt", header=T)
nIto_APVsAPDs <- read.table("./cell_1_APVsAPDs_tmV/APDsAPVs_ito_decrease_122420_3_drew.txt", header=T)
pIto_APVsAPDs <- read.table("./cell_1_APVsAPDs_tmV/APDsAPVs_ito_increase_122420_3_drew.txt", header=T)

tmV_list <- list(cntrl_tmV, nICaL_tmV, pICaL_tmV, nIto_tmV, pIto_tmV)
APVsAPDs_list <- list(cntrl_APVsAPDs, nICaL_APVsAPDs, pICaL_APVsAPDs,
                      nIto_APVsAPDs, pIto_APVsAPDs)
                 

cols <- c("black", "blue", "blue", "orangered3", "orangered3")
ltys <- c("solid", "dotdash", "solid", "solid", "dotdash")


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


