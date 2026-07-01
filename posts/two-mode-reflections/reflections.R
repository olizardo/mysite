reflections <- function(x, iter = 20) { #x is a matrix with people as rows and groups as columns iters is number of reflections
   p <- nrow(x)
   g <- ncol(x)
   p.c <- matrix(0, p, iter) #initialize person centralities
   g.c <- matrix(0, g, iter) #initialize group centralities trajectory matrix
   
   p.s <- matrix(0, p, iter) #initialize person centralities
   g.s <- matrix(0, g, iter) #initialize group centralities trajectory matrix
   
   p.r <- matrix(0, p, iter) #initialize person centralities
   g.r <- matrix(0, g, iter) #initialize group centralities trajectory matrix
   
   rownames(p.c) <- rownames(x)
   rownames(g.c) <- colnames(x)
   
   rownames(p.s) <- rownames(x)
   rownames(g.s) <- colnames(x)
   
   rownames(p.r) <- rownames(x)
   rownames(g.r) <- colnames(x)
   
   colnames(p.c) <- paste("Cr", c(1:iter), sep = "")
   colnames(g.c) <- paste("Cr_", c(1:iter), sep = "")
   
   colnames(p.s) <- paste("Cr_", c(1:iter), sep = "")
   colnames(g.s) <- paste("Cr_", c(1:iter), sep = "")
   
   colnames(p.r) <- paste("Cr_", c(1:iter), sep = "")
   colnames(g.r) <- paste("Cr_", c(1:iter), sep = "")
   
   p.c[, 1] <- rowSums(x) #person degree centrality 
   g.c[, 1] <- colSums(x) #group degree centrality 
   
   k <- 1 #initializing counter
   while (k < iter) {
      m <- k + 1
      for(i in 1:p) {
         p.c[i, m] <- sum(x[i, ] * g.c[, k]) * (1/p.c[i, 1]) #assign person avg. centrality groups they belong to
      } #end person loop
      for(j in 1:g) {
         g.c[j, m] <- sum(x[, j] * p.c[, k]) * (1/g.c[j, 1]) #assign genre avg. person centrality of people in group
      } #end group loop
      k <- k + 1 #increase counter
   } #end while loop
   for (j in 1:iter) {
      p.s[, j] <- scale(p.c[, j]) #rescaling person reflections
      g.s[, j] <- scale(g.c[, j]) #rescaling group reflections
   }
   for (j in 1:iter) {
      p.r[, j] <- rank(p.s[, j]) #ranking rescaled person reflections
      g.r[, j] <- rank(g.s[, j]) #ranking rescaled group reflections
   }
   p.r <- (max(p.r) + 1) - p.r #reversing rank values so that smaller is highest
   return(list(p.c = p.c, g.c = g.c, p.s = p.s, g.s = g.s, p.r = p.r, g.r = g.r))
} #end function