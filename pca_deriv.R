source("deltaSVD.R")

nr = 10
nc = 5
x = matrix(rnorm(nr*nc), nr, nc)
xsvd = svd(x)

xtx = t(x) %*% x
xxt = x %*% t(x)

(xxt %*% xsvd$u) / (xsvd$u %*% diag(xsvd$d) %*% diag(xsvd$d))
(xtx %*% xsvd$v) / (xsvd$v %*% diag(xsvd$d) %*% diag(xsvd$d))


## -----------------------------------------------------------------
## try weighting...
## -----------------------------------------------------------------
w = diag(c(rep(1, nr-1), sqrt(2)))
wx = w %*% x
wxsvd = svd(wx)

x2 = rbind(
    x,
    x[5, ]
)
x2svd = svd(x2)

## sum(abs(wxsvd$d - x2svd$d)) < 1e-10
## sum(abs(wxsvd$v - x2svd$v)) < 1e-10

w = rep(1, nr)
set.seed(1234)
dw = 0.0001 * rnorm(nr)

xdw = diag(w + dw) %*% x
xdwsvd = mapply(`-`, svd(xdw), svd(x))
xdsvd = deltaSVD(x, w, dw, xsvd=xsvd)

mapply(`/`, xdwsvd, xdsvd)

dx = matrix(0, nrow(x), ncol(x))
dx[1, 2] = 0.001
xdx = x + dx
xdxsvd = mapply(`-`, svd(xdx), svd(x))
xdxsvd2 = deltaSVD(x, dx=dx, xsvd=xsvd)

mapply(`/`, xdxsvd, xdxsvd2)



## -----------------------------------------------------------------
## test eigenvalue equations
## -----------------------------------------------------------------
froeb = function(m) {sqrt(sum(m^2))}

froeb(xxt %*% xsvd$u - xsvd$u %*% diag(xsvd$d^2))
froeb(
    (xdw %*% t(xdw)) %*% (xsvd$u + xdwsvd$u) -
    (xsvd$u + xdwsvd$u) %*% diag((xsvd$d + xdwsvd$d)^2)
)

froeb(
    (xdw %*% t(xdw)) %*% (xsvd$u + xdsvd$u) -
    (xsvd$u + xdsvd$u) %*% diag((xsvd$d + xdsvd$d)^2)
)
froeb(
    (t(xdw) %*% xdw) %*% (xsvd$v + xdsvd$v) -
    (xsvd$v + xdsvd$v) %*% diag((xsvd$d + xdsvd$d)^2)
)

froeb((xdw %*% t(xdw)) %*% xdsvd$u)
froeb(xdsvd$u %*% diag((xsvd$d + xdsvd$d)^2))
