library(inline)

deltaSVD_slow = function(x, w=rep(1, nrow(x)),
        dw=NULL, dx=NULL, xsvd=svd(x)) {
    require(Matrix)
    if (is.vector(w)) {w = Diagonal(x=w)}
    if (length(dw)>0 && is.vector(dw)) {dw = Diagonal(x=dw)}
    u = xsvd$u
    v = xsvd$v
    d = xsvd$d
    if (nrow(x) >= ncol(x)) {
        n = if (length(dx) > 0) {
            t(dx) %*% x + t(x) %*% dx
        } else if (length(dw) > 0) {
            2 * t(x) %*% w %*% dw %*% x
        }
        vtnv = t(v) %*% n %*% v
        dd = diag(vtnv) / (2*d)
        invD = diag(1 / d)
        dv = matrix(0, nrow(v), ncol(v))
         for (f in 1:nrow(v)) {
            for (r in 1:ncol(v)) {
                for (q in setdiff(1:ncol(v), r)) {
                    dv[f, r] = dv[f, r] +
                            v[f, q] * vtnv[q, r] / (d[r]^2 - d[q]^2)
                }
            }
        }
        du = if (length(dx) > 0) {
            dx %*% v %*% invD -
            u %*% diag(dd) %*% invD -
            u %*% diag(d) %*% t(dv) %*% v %*% invD
        } else if (length(dw) > 0) {
            dw %*% x %*% v %*% invD -
            u %*% diag(dd) %*% invD -
            u %*% diag(d) %*% t(dv) %*% v %*% invD
        }
    } else if (ncol(x) > nrow(x)) {
        xxt = x %*% t(x)
        m = if (length(dx) > 0) {
            dx %*% t(x) + x %*% t(dx)
        } else if (length(dw) > 0) {
            dw %*% xxt %*% w + w %*% xxt %*% dw
        }
        utmu = t(u) %*% m %*% u
        dd = diag(utmu) / (2*d)
        invD = diag(1 / d)
        du = matrix(0, nrow(u), ncol(u))
        for (k in 1:nrow(u)) {
            for (r in 1:ncol(u)) {
                for (q in setdiff(1:ncol(u), r)) {
                    du[k, r] = du[k, r] +
                            u[k, q] * utmu[q, r] / (d[r]^2 - d[q]^2)
                }
            }
        }
        ut = t(u)
        dv = if (length(dx) > 0) {
            t(
                invD %*% ut %*% dx -
                invD %*% ut %*% du %*% diag(d) %*% t(v) -
                invD %*% diag(dd) %*% t(v)
            )
        } else if (length(dw) > 0) {
            t(
                invD %*% ut %*% dw %*% x -
                invD %*% ut %*% du %*% diag(d) %*% t(v) -
                invD %*% diag(dd) %*% t(v)
            )
        }
    }
    return(list(d=dd, u=du, v=dv))
}


deltaSVD_cvCode = '
int vnr0 = (*vnr);
int vnc0 = (*vnc);
for (int f=0; f<vnr0; f++) {
    int fOffset = f * vnr0;
    for (int r=0; r<vnc0; r++) {
        for (int q=0; q<vnc0; q++) {
            if (q != r) {
                dvFlat[fOffset + r] += (vFlat[fOffset + q] *
                        vtnvFlat[q*vnc0 + r] / (d[r]*d[r] - d[q]*d[q]));
            }
        }
    }
}
'
deltaSVD_cv = cfunction(
    signature(
        vnr = "numeric",
        vnc = "numeric",
        dvFlat = "numeric",
        vFlat = "numeric",
        vtnvFlat = "numeric",
        d = "numeric"
    ),
    body = deltaSVD_cvCode,
    language = "C",
    convention = ".C"
)

deltaSVD_cuCode = '
int unr0 = (*unr);
int unc0 = (*unc);
for (int k=0; k<unr0; k++) {
    int kOffset = k * unr0;
    for (int r=0; r<unc0; r++) {
        for (int q=0; q<unc0; q++) {
            if (q != r) {
                duFlat[kOffset + r] += (uFlat[kOffset + q] *
                        utmuFlat[q*unc0 + r] / (d[r]*d[r] - d[q]*d[q]));
            }
        }
    }
}
'
deltaSVD_cu = cfunction(
    signature(
        unr = "numeric",
        unc = "numeric",
        duFlat = "numeric",
        uFlat = "numeric",
        utmuFlat = "numeric",
        d = "numeric"
    ),
    body = deltaSVD_cuCode,
    language = "C",
    convention = ".C"
)



deltaSVD = function(x, w=rep(1, nrow(x)),
        dw=NULL, dx=NULL, xsvd=svd(x)) {
    require(Matrix)
    if (is.vector(w)) {w = Diagonal(x=w)}
    if (length(dw)>0 && is.vector(dw)) {dw = Diagonal(x=dw)}
    u = xsvd$u
    v = xsvd$v
    d = xsvd$d
    if (nrow(x) >= ncol(x)) {
        n = if (length(dx) > 0) {
            t(dx) %*% x + t(x) %*% dx
        } else if (length(dw) > 0) {
            2 * t(x) %*% w %*% dw %*% x
        }
        vtnv = t(v) %*% n %*% v
        dd = diag(vtnv) / (2*d)
        invD = diag(1 / d)
        dv = t(matrix(deltaSVD_cv(
            vnr = nrow(v),
            vnc = ncol(v),
            dvFlat = numeric(nrow(v) * ncol(v)),
            vFlat = as.vector(t(v)),
            vtnvFlat = as.vector(t(vtnv)),
            d = d
        )$dvFlat, nrow=nrow(v), ncol=ncol(v)))
        du = if (length(dx) > 0) {
            dx %*% v %*% invD -
            u %*% diag(dd) %*% invD -
            u %*% diag(d) %*% t(dv) %*% v %*% invD
        } else if (length(dw) > 0) {
            dw %*% x %*% v %*% invD -
            u %*% diag(dd) %*% invD -
            u %*% diag(d) %*% t(dv) %*% v %*% invD
        }
    } else if (ncol(x) > nrow(x)) {
        xxt = x %*% t(x)
        m = if (length(dx) > 0) {
            dx %*% t(x) + x %*% t(dx)
        } else if (length(dw) > 0) {
            dw %*% xxt %*% w + w %*% xxt %*% dw
        }
        utmu = t(u) %*% m %*% u
        dd = diag(utmu) / (2*d)
        invD = diag(1 / d)
        du = t(matrix(deltaSVD_cu(
            unr = nrow(u),
            unc = ncol(u),
            duFlat = numeric(nrow(u) * ncol(u)),
            uFlat = as.vector(t(u)),
            utmuFlat = as.vector(t(utmu)),
            d = d
        )$duFlat, nrow=nrow(u), ncol=ncol(u)))
        ut = t(u)
        dv = if (length(dx) > 0) {
            t(
                invD %*% ut %*% dx -
                invD %*% ut %*% du %*% diag(d) %*% t(v) -
                invD %*% diag(dd) %*% t(v)
            )
        } else if (length(dw) > 0) {
            t(
                invD %*% ut %*% dw %*% x -
                invD %*% ut %*% du %*% diag(d) %*% t(v) -
                invD %*% diag(dd) %*% t(v)
            )
        }
    }
    return(list(d=dd, u=du, v=dv))
}


dwScale = function(dw, w=rep(1, length(dw)), norm=1e-3) {
    ## project w out of dw
    dw = dw - (sum(dw*w) / sum(w*w)) * w
    ## now scale to have specified norm
    dw = norm * dw / sqrt(sum(dw*dw))
    return(dw)
}


froeb = function(m) {sqrt(sum(m^2))}


pseudoFMediods = function(x, cl) {
    x = as.matrix(x)
    clunique = unique(cl)
    names(clunique) = as.character(clunique)
    overallMediod = apply(x, 2, median)
    xSwept = sweep(x, 2, overallMediod, `-`)
    xDistSum = sum(abs(xSwept))
    xInGroupDist = sum(sapply(
        clunique,
        function(clIndex) {
            xcl = x[cl==clIndex, , drop=FALSE]
            mediod = apply(xcl, 2, median)
            swept = sweep(xcl, 2, mediod, `-`)
            return(sum(abs(swept)))
        }
    ))
    return(xDistSum / xInGroupDist)
}


silhouetteIndex = function(x, cl, ...) {
    require(clusterSim)
    return(clusterSim::index.S(d=dist(x, ...), cl=cl))
}


ldaLlhdRatio = function(x, cl, diag=FALSE) {
    x = as.matrix(x)
    ## don't need to apply add mlDenAdj to scores b/c cancels
    ## mlDenAdj = (log(nrow(x)-1) - log(nrow(x))) / 2
    noGroupScore = if (diag) {
        -sum(log(apply(x, 2, var)))
    } else {
        -log(abs(det(cov(x)))) / 2
    }
    for (cluster in unique(cl)) {
        clin = (cl == cluster)
        x[clin, ] = sweep(x[clin, ], 2, colMeans(x[clin, ]), `-`)
    }
    groupedScore = if (diag) {
        -sum(log(apply(x, 2, var)))
    } else {
        -log(abs(det(cov(x)))) / 2
    }
    return(groupedScore - noGroupScore)
}


centerScale = function(x, center=c("none", "column", "row", "both"),
        scale=c("none", "column", "row")) {
    center = match.arg(center)
    scale = match.arg(scale)
    if (center %in% c("row", "both")) {
        x = sweep(x, 1, STATS=rowMeans(x))
    }
    if (center %in% c("column", "both")) {
        x = sweep(x, 2, STATS=colMeans(x))
    }
    if (scale == "row") {
        x = sweep(x, 1, STATS=apply(X=x, MARGIN=1, FUN=sd), FUN=`/`)
    } else if (scale == "column") {
        x = sweep(x, 2, STATS=apply(X=x, MARGIN=2, FUN=sd), FUN=`/`)
    }
    return(x)
}


dPcaPseudoF = function(x, y, dw, k=3,
        w=rep(1, ncol(x)), norm=1e-3,
        center=c("both", "column", "row", "none"),
        scale=c("none", "column", "row"),
        xsvd=NULL, cq=NULL) {
    require(clusterSim)
    require(Matrix)
    if (length(cq) == 0) {cq = clusterSim::index.G1}
    if (is.character(dw)) {
        dwInds = dw
        dw = rep(0, ncol(x))
        dw[colnames(x) %in% dwInds] = 1
    }
    dw = dwScale(dw=dw, w=w, norm=norm)
    if (center!="none" || scale!="none") {
        x = centerScale(x, center=center, scale=scale)
    }
    if (length(xsvd) == 0) {
        xsvd = svd(t(x))
    }
    dxsvd = deltaSVD(x=t(x), w=w, dw=dw, xsvd=xsvd)
    v0 = xsvd$v[ , 1:k, drop=FALSE] %*% Diagonal(x=xsvd$d[1:k])
    v1 = v0
    v1 = v1 + dxsvd$v[ , 1:k, drop=FALSE] %*% Diagonal(x=xsvd$d[1:k])
    v1 = v1 + xsvd$v[ , 1:k, drop=FALSE] %*% Diagonal(x=dxsvd$d[1:k])
    if (is.character(y)) {y = factor(y)}
    if (is.factor(y)) {y = as.integer(y)}
    cq0 = cq(v0, y)
    cq1 = cq(v1, y)
    return((cq1 - cq0) / norm)
}


ggDeltaSVD = function(
        x,
        y,
        dw,
        small = 1e-4,
        mag = 2.5,
        center = c("both", "column", "row", "none"),
        scale = c("none", "column", "row"),
        rlab = FALSE,
        clab = FALSE,
        cshow = ncol(x),
        rsize = 4,
        csize = 2,
        lsize = 3,
        ralpha0 = 0.75,
        ralpha1 = 0.4,
        calpha0 = 0.75,
        calpha1 = 0.4,
        rname = "Sample",
        cname = "Variable",
        lname = "",
        grid = FALSE,
        print = TRUE,
        colscale,
        ...) {
    require(ggplot2)
    center = match.arg(center)
    scale = match.arg(scale)
    if (length(rlab)==1 && is.logical(rlab)) {
        rlab = if (rlab) {rownames(x)} else {rep("", nrow(x))}
    }
    if (length(clab)==1 && is.logical(clab)) {
        clab = if (clab) {colnames(x)} else {rep("", ncol(x))}
    }
    if (!missing(y)) {
        if (is.character(y)) {y = factor(y, levels=unique(y))}
        if (length(names(y)) == 0) {names(y) = rownames(x)}
        classLevels = c(cname, levels(y))
        y = structure(as.character(y), names=names(y))
    } else {
        classLevels = c(cname, rname)
    }
    x = x[ , sapply(x, function(z) {!any(is.na(z))}), drop=FALSE]
    if (is.character(dw)) {
        dwn = dw
        dw = structure(rep(0, ncol(x)), names=colnames(x))
        dw[names(dw) %in% dwn] = 1
        dw = dwScale(dw, norm=small)
    }
    x = centerScale(x, center=center, scale=scale)
    xsvd0 = svd(t(x))
    dxsvd = deltaSVD(t(x), dw=dw, xsvd=xsvd0)
    xsvd1 = list(
        d = xsvd0$d + (mag / small) * dxsvd$d,
        u = xsvd0$v + (mag / small) * dxsvd$v,
        v = xsvd0$u + (mag / small) * dxsvd$u
    )
    xsvd0 = list(d=xsvd0$d, u=xsvd0$v, v=xsvd0$u)
    rsf = max(xsvd0$u[ , 1]) - min(xsvd0$u[ , 1])
    csf = max(xsvd0$v[ , 1]) - min(xsvd0$v[ , 1])
    ggdata = data.frame(
        PC1 = c(xsvd0$u[ , 1], xsvd1$u[ , 1]) / rsf,
        PC2 = c(xsvd0$u[ , 2], xsvd1$u[ , 2]) / rsf,
        label = c(rlab, rep("", nrow(x))),
        size = rsize,
        alpha = c(rep(ralpha0, nrow(x)), rep(ralpha1, nrow(x))),
        stringsAsFactors = FALSE
    )
    if (cshow > 0) {
        cdata = data.frame(
            PC1 = c(xsvd0$v[ , 1], xsvd1$v[ , 1]) / csf,
            PC2 = c(xsvd0$v[ , 2], xsvd1$v[ , 2]) / csf,
            label = c(clab, rep("", ncol(x))),
            size = csize,
            alpha = c(rep(calpha0, ncol(x)), rep(calpha1, ncol(x))),
            stringsAsFactors = FALSE
        )
        if (cshow < ncol(x)) {
            cscores = cdata$PC1^2 + cdata$PC2^2
            cscores[ncol(x):length(cscores)] = 0
            names(cscores) = rep(colnames(x), 2)
            cscores = sort(cscores, decreasing=TRUE)
            ctoshow = which(colnames(x) %in% names(cscores[1:cshow]))
            ctoshow = c(ctoshow, ctoshow + ncol(x))
            cdata = cdata[ctoshow, ]
        }
        ggdata = rbind(cdata, ggdata)
    }
    cclass = rep(cname, times=2*cshow)
    if (!missing(y)) {
        ggdata$class = factor(
            c(cclass, rep(y, times=2)),
            levels = classLevels
        )
    } else {
        ggdata$class = factor(
            c(cclass, rep(rname, times=2*nrow(x))),
            levels = classLevels
        )
    }
    ggobj = ggplot(
        aes(
            x = PC1,
            y = PC2,
            color = class,
            size = size,
            alpha = alpha,
            label = label
        ), 
        data = ggdata
    ) + geom_point() + theme_bw()
    ggobj = ggobj + geom_text(vjust=-1.1, show_guide=FALSE, size=lsize)
    if (missing(colscale) && (length(unique(ggdata$class)) < 8)) {
        colscale = c("black", "firebrick", "dodgerblue",
                "goldenrod", "seagreen", "darkorchid4", "darkgray")[
                1:length(unique(ggdata$class))]
        if (length(colscale) <= 3) {colscale[2] = "red"}
    }
    if (all(classLevels %in% names(colscale))) {
        colscale = colscale[classLevels]
    }
    ggobj = ggobj + scale_color_manual(values=colscale, name=lname)
    ggobj = ggobj + scale_size_continuous(
            guide=FALSE, range=sort(c(rsize, csize)))
    ggobj = ggobj + scale_alpha_continuous(guide=FALSE, range=c(
        min(ralpha0, ralpha1, calpha0, calpha1),
        max(ralpha0, ralpha1, calpha0, calpha1)
    ))
    ggobj = ggobj + xlab(paste0(
        "PC1 (",
        round(100 * xsvd0$d[1]^2 / sum(xsvd0$d^2), 1),
        "% explained var.)"
    ))
    ggobj = ggobj + ylab(paste0(
        "PC 2 (",
        round(100 * xsvd0$d[2]^2 / sum(xsvd0$d^2), 1),
        "% explained var.)"
    ))
    if (!grid) {
        ggobj = ggobj + theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.background = element_blank()
        )
    }
    if (print) {
        print(ggobj)
    }
    invisible(ggobj)
}
