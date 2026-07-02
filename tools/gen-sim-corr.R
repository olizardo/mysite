gen.sim.corr.abs <- function(x, sigma = 0.001) {
        library(expm) #package required for matrix multiplication in R
        r <- nrow(x) #Number of rows in the input matrix
        c <- ncol(x) #Number of columns in the input matrix
        r.c <- diag(r) #Identity matrix of dimensions r X r
        c.c <- diag(c) #Identity matrix of dimensions c X c
        r.m <- rowMeans(x) #row means
        c.m <- colMeans(x) #column means
        d.r.c <- 1 #initializing difference values
        d.c.c <- 1 #initializing difference values
        k <- 1 #initializing iteration counter
        while (d.r.c > sigma | d.c.c > sigma) { #iterate until difference is less than sigma
                p.r.c <- r.c
                p.c.c <- c.c
                for (i in 1: r) { #compute pairwise generalized similarities for rows
                        for (j in 1:r) {
                                if (i != j) {
                                        r.x <- x[i, ] - r.m[i]
                                        r.y <- x[j, ] - r.m[j]
                                        r.xy <- r.x %*% c.c * t(r.y)
                                        r.xx <- r.x %*% c.c * t(r.x)
                                        r.yy <- r.y %*% c.c * t(r.y) 
                                        r.num <- sum(r.xy)
                                        r.den <- sqrt(sum(r.xx)) * sqrt(sum(r.yy))
                                        r.c[i, j] <- r.num / r.den
                                }
                        }
                }
                for (i in 1: c) {
                        for (j in 1:c) { #compute pairwise generalized similarities for columns
                                if (i != j) {
                                        c.x <- x[, i] - c.m[i]
                                        c.y <- x[, j] - c.m[j]
                                        c.xy <- c.x %*% r.c * t(c.y) 
                                        c.xx <- c.x %*% r.c * t(c.x) 
                                        c.yy <- c.y %*% r.c * t(c.y) 
                                        c.num <- sum(c.xy)
                                        c.den <- sqrt(sum(c.xx)) * sqrt(sum(c.yy))
                                        c.c[i, j] <- c.num / c.den
                                }
                        }
                }
                d.r.c <- sum(rowSums(abs(p.r.c - r.c)))
                d.c.c <- sum(rowSums(abs(p.c.c - c.c)))
                k <- k + 1
        }
        rownames(r.c) <- rownames(x)
        colnames(r.c) <- rownames(x)
        rownames(c.c) <- colnames(x)
        colnames(c.c) <- colnames(x)
        return(list(row.sims = r.c, col.sims = c.c, n.iter = k)) #return rowwise and columnwise generalized similarities
}

