source("~/projects/code-share-112018/jgpplotter.R")
## Reads in action potentials from evaluated models of iPSC EA fits

read_apset  <- function(file_seq=seq(0,4), file_prefix="cell_2_") {
    if (file_seq[1] != 0) {print("Warning: Sequencing did not begin at 0.")}
    ap.files <- list.files(pattern=file_prefix)
    models  <- list()
    for (i in file_seq) {
        ndx.pattern  <- paste("_",as.character(i),".txt",sep="")
        filenames  <- ap.files[grep(pattern=ndx.pattern, x=ap.files)]
        tmp  <- list()
        ## Assign APs to index
        tmp[[1]] <- read.table(file=filenames[grep(pattern="_cntrl_", x=filenames)],header=T)
        tmp[[2]] <- read.table(file=filenames[grep(pattern="_-0.15_ical_", x=filenames)],header=T)
        tmp[[3]] <- read.table(file=filenames[grep(pattern="_0.7_ical_", x=filenames)],header=T)
        tmp[[4]] <- read.table(file=filenames[grep(pattern="_-0.25_ikr_", x=filenames)],header=T)
        tmp[[5]] <- read.table(file=filenames[grep(pattern="_0.9_ikr_", x=filenames)],header=T)
        tmp[[6]] <- read.table(file=filenames[grep(pattern="_-0.9_ito_", x=filenames)],header=T)
        tmp[[7]] <- read.table(file=filenames[grep(pattern="_1.5_ito_", x=filenames)],header=T)
        tmp[[8]] <- read.table(file=filenames[grep(pattern="_10_iks_", x=filenames)],header=T)
        tmp[[9]] <- read.table(file=filenames[grep(pattern="_4_iks_", x=filenames)],header=T)

        models[[(i+1)]]  <- tmp
    }
    return(models)
}

get_apset_order  <- function() {
    return(c("cntrl", "-0.15_ical", "0.7_ical",
             "-0.25_ikr", "0.9_ikr", "-0.9_ito",
             "1.5_ito", "10_iks", "4_iks"))
}

get_apset_g_scale <- function() {
    return(c("+100% IK1", "-15% ICaL", "+70% ICaL",
             "-25% IKr", "+90% IKr", "-90% Ito",
             "+150% Ito", "+1000% IKs", "+400% IKs"))
}

get_apset_pApF <- function(experssion_format=FALSE) {
    if (experssion_format) {
        return(c(expression(paste(+1.875," pA/pF ", I[K1])),
                 expression(paste(-0.0462," pA/pF ", I[CaL])),
                 expression(paste(+0.2156," pA/pF ", I[CaL])),
                 expression(paste(-0.0545," pA/pF ", I[Kr])),
                 expression(paste(+0.1962," pA/pF ", I[Kr])),
                 expression(paste(-0.1061," pA/pF ", I[to])),
                 expression(paste(+0.1768," pA/pF ", I[to])),
                 expression(paste(+0.077," pA/pF ", I[Ks])),
                 expression(paste(+0.0308," pA/pF ", I[Ks]))))
    }
    else {
        return(c("+2.5 pA/pF IK1", "-0.0462 pA/pF ICaL", "+0.2156 pA/pF ICaL",
             "-0.0545 pA/pF IKr", "+0.1962 pA/pF IKr", "-0.1061 pA/pF Ito",
             "+0.1768 pA/pF Ito", "+0.077 pA/pF IKs", "+0.0308 pA/pF IKs"))
    }
}

get_apset_cols  <- function() {
    return(c("dimgrey", "blue", "blue4",
             "darkcyan", "cyan2", "orangered4",
             "orangered", "magenta2", "darkmagenta"))
}

## Color by error
err2col <- function(err, na.col="dimgrey") {
    colramp_fun <- colorRampPalette(c("#1b98e0", "red"))
    ## 100 colors for % normalized error
    err_cols <- colramp_fun(101)
    norm_err <- (err/max(err, na.rm = T))*100 ## in %
    cols <- list()
    for (i in seq_len(length(norm_err))) {
        if (is.na(norm_err[i])) {
            cols[i] <- na.col
        } else {
            cols[i] <- err_cols[(as.integer(round(norm_err[i], digits = 0)))+1]
        }
    }
    return(unlist(cols))
}

get_apset_rgb  <- function(scaled=FALSE) {
    cols <- c("dimgrey", "blue", "blue4",
              "darkcyan", "cyan2", "orangered4",
              "orangered", "magenta2", "darkmagenta")
    cols.rgb <- as.data.frame(sapply(X=cols, FUN = col2rgb))
    if (scaled) {
        return(cols.rgb/255)
    } else {
        return(cols.rgb)
    }
}

get_apset_ltys  <- function() {
    return(c("solid", "dotdash", "solid",
             "solid", "dotdash", "solid",
             "dotdash", "dotdash", "solid"))
}

get_apset_borders <- function() {
    border <- c(NA, "blue", NA,
                NA, "cyan2", NA,
                "orangered", "magenta2", NA)
    return(border)
}

get_apset_shading <- function() {
    
    alpha <- c(0.4, 1.0, 0.4,
               0.4, 1.0, 0.4,
               1.0, 1.0, 0.4)
    
    density <- c(NA, 35, NA,
                 NA, 35, NA,
                 35, 35, NA)

    angle <- c(NA, 45, NA,
               NA, -45, NA,
               -45, 45, NA)

    shading <- as.data.frame(cbind(alpha, density, angle))
    return(shading)
}

get_apset_pchs  <- function() {
    return(c(19, 6, 17, 17, 6, 15, 12, 12, 15))
}

plot_model_set  <- function(modelSet, ap.index=1, cols="black", ref.col="red", ...) {
    ap_label <- get_apset_order()[ap.index]
    print(paste("AP: ", ap_label, sep=""))
    ap  <- as.data.frame(modelSet[[1]][[ap.index]])
    jgpplotter(ap$t, ap$mV_cell, type="n", axes=FALSE, ...)
    for (i in seq(length(modelSet))) {
        m <- as.data.frame(modelSet[[i]][[ap.index]])
        if (length(cols) == length(modelSet)) {
            lines(m$t, m$mV_simu, col=cols[i])
        } else {
            lines(m$t, m$mV_simu, col=cols)
        }
    }
    lines(ap$t, ap$mV_cell, col=ref.col, lty=2)
}

##  Euclidean Distance from origin
get.euclidean <- function(X, ref=NULL) {
    if (is.null(ref)) { 
        d  <- as.data.frame((apply((X*X), MARGIN=1, FUN=sum))^0.5)
    }
    else if (ncol(X) == length(ref)){
        d <- as.data.frame((apply(((X-ref)^2), MARGIN=1, FUN=sum))^0.5)
    } else {
        stop("Dimensions of X and ref did not match.")
    }
    names(d)  <- c("euclidean")
    return(d)
}

##  Get fitness values from scrs
get.fitness  <- function(scrs) {
    fitness <- as.data.frame(apply(X=scrs, MARGIN=1, FUN=sum))
    names(fitness)  <- c("fitness")
    return(fitness)
}
