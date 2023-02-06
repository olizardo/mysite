reflections <- function(x, iter = 20, ref.type = "hidalgo") { #x is a matrix with people as rows and group.s as columns iters is number of reflections
   library(expm) #package for matrix operations
   #Number of persons and groups
   p <- nrow(x)
   g <- ncol(x)
   
   #Labels for persons, groups, and reflections
   pnames <- rownames(x)
   gnames <- colnames(x)
   rnames <- paste("R", c(1:iter), sep = "")
   
   #Person and group degree centrality vectors
   p.c <- rowSums(x) 
   g.c <- colSums(x) 
   
   #Person and group degree matrices
   p.dm <- diag(p.c, p, p)
   g.dm <- diag(g.c, g, g)
   
   #Initializing trajectory matrices
   p.ref0 <- matrix(0, p, iter) 
   p.ref1 <- matrix(0, p, iter) 
   p.ref2 <- matrix(0, p, iter) 
   g.ref0 <- matrix(0, g, iter)
   g.ref1 <- matrix(0, g, iter) 
   g.ref2 <- matrix(0, g, iter) 
   
   #Assigning person and group row names to trajectory matrices
   rownames(p.ref0) <- pnames
   rownames(p.ref1) <- pnames
   rownames(p.ref2) <- pnames
   rownames(g.ref0) <- gnames
   rownames(g.ref1) <- gnames
   rownames(g.ref2) <- gnames
   
   #Assigning reflection number column names to trajectory matrices
   colnames(p.ref0) <- rnames
   colnames(p.ref1) <- rnames
   colnames(p.ref2) <- rnames
   colnames(g.ref0) <- rnames
   colnames(g.ref1) <- rnames
   colnames(g.ref2) <- rnames
   
   #Initializing first reflection of trajectory matrices with person and group degrees
   p.ref0[, 1] <- p.c
   p.ref1[, 1] <- p.c
   p.ref2[, 1] <- p.c
   g.ref0[, 1] <- g.c
   g.ref1[, 1] <- g.c
   g.ref2[, 1] <- g.c
   
   #Creating similarity matrices
   p.s0 <- x %*% t(x) #projection matrix
   p.s1 <- x %*% solve(g.dm) %*% t(x) #degree-weighted similarity
   p.s2 <- solve(p.dm) %*% x %*% solve(g.dm) %*% t(x) #reflective similarity
   g.s0 <- t(x) %*% x #projection matrix
   g.s1 <- t(x) %*% solve(p.dm) %*% x #degree-weighted similarity
   g.s2 <- solve(g.dm) %*% t(x) %*% solve(p.dm) %*% x #reflective similarity
   
   # Computing centralities
   p.eig0 <- as.numeric(eigen(p.s0)$vectors[, 1] * -1)
   p.eig1 <- as.numeric(eigen(p.s1)$vectors[, 2])
   p.eig2 <- as.numeric(eigen(p.s2)$vectors[, 2])
   g.eig0 <- as.numeric(eigen(g.s0)$vectors[, 1] * -1)
   g.eig1 <- as.numeric(eigen(g.s1)$vectors[, 2])
   g.eig2 <- as.numeric(eigen(g.s2)$vectors[, 2])
   
   # Storing centralities in data frame
   p.eig.dat <- round(data.frame(cbind(p.eig0, p.eig1, p.eig2)), 3)
   rownames(p.eig.dat) <- pnames
   colnames(p.eig.dat) <- c("PP", "PP(S)", "PP(R)")
   g.eig.dat <- round(data.frame(cbind(g.eig0, g.eig1, g.eig2)), 3)
   rownames(g.eig.dat) <- gnames
   colnames(g.eig.dat) <- c("GG", "GG(S)", "GG(R)")
   
   #While loop to populate reflection matrices
   k <- 1 #initializing counter
   while (k < iter) {
      m <- k + 1
      for(i in 1:p) {
         p.ref0[i, m] <- sum(x[i, ] * g.ref0[, k])
         p.ref2[i, m] <- sum(x[i, ] * g.ref2[, k]) * 1/p.c[i] 
      } #end person loop
      for(j in 1:g) {
         g.ref0[j, m] <- sum(x[, j] * p.ref0[, k])
         g.ref2[j, m] <- sum(x[, j] * p.ref2[, k]) * 1/g.c[j] 
      } #end group loop
      k <- k + 1 #increase counter
   } #end while loop
   
   persons <- list(ref = list(ref0 = p.ref0, ref1 = p.ref1, ref2 = p.ref2),
                   sim = list(sim0 = p.s0, sim1 = p.s1, sim2 = p.s2),
                   cen = p.eig.dat)
   
   groups <- list(ref = list(ref0 = g.ref0, ref1 = g.ref1, ref2 = g.ref2),
                   sim = list(sim0 = g.s0, sim1 = g.s1, sim2 = g.s2),
                   cen = g.eig.dat)
   return(list(persons = persons, groups = groups))
} #end function