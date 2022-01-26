sort.ipscParam.clusters <- function(params, clusters, returnDF=FALSE,
                                    file.prefix="cell1_clstr_", write.files=FALSE) {
    if (nrow(params) != nrow(clusters)) {
        print("Number of individuals in params and clusters does not match.")
        return()
    }
    tmp <- as.data.frame(cbind(params, clusters))
    tmp.sort <- tmp[order(tmp$clstrs,tmp$fitness),]
    if (write.files) {
        nclusters <- unique(clusters)$clstrs
        for (i in nclusters) {
            filename <- paste(file.prefix, i, ".txt", sep="")
            tmp.sub <- subset(tmp.sort, clstrs == i, select=phi:clstrs)
            write.table(x=tmp.sub, file=filename, row.names=FALSE, quote=FALSE)
        }
    }
    if (returnDF) {
        return(tmp.sort)
    }
}


sep.ipscParam.clusters <- function(params, clusters) {
    if (nrow(params) != nrow(clusters)) {
        print("Number of individuals in params and clusters does not match.")
        return()
    }
    tmp <- as.data.frame(cbind(params, clusters))
    tmp.sort <- tmp[order(tmp$clstrs,tmp$fitness),]
    tmp.list <- list()
    nclusters <- unique(clusters)$clstrs
    for (i in nclusters) {
        tmp.list[[i]] <- subset(tmp.sort, clstrs == i, select=phi:clstrs)
    }
    return(tmp.list)
}

get.top.clusters <- function(params, clusters) {
    if (nrow(params) != nrow(clusters)) {
        print("Number of individuals in params and clusters does not match.")
        return()
    }
    tmp <- as.data.frame(cbind(params, clusters))
    tmp.sort <- tmp[order(tmp$clstrs,tmp$fitness),]
    tmp.list <- list()
    nclusters <- unique(clusters)$clstrs
    topCLSTRS <- NULL
    for (i in nclusters) {
        tmp <- subset(tmp.sort, clstrs == i, select=phi:clstrs)
        topCLSTRS <- rbind(topCLSTRS, tmp[1,])
    }
    topCLSTRS <- topCLSTRS[order(topCLSTRS$fitness),]
    return(topCLSTRS)
}
