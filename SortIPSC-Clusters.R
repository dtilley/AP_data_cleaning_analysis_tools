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
