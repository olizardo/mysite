ref.long.dat <- function(x, max.iter) {
     library(tidyverse)
     evens <- function(x) subset(x, x %% 2 == 0)
     odds <- function(x) subset(x, x %% 2 != 0)
     pr.long <- data.frame(x$person$rr)
     pr.long <- pr.long %>% 
          mutate(person = rownames(pr.long)) %>% 
          pivot_longer(
               cols = 1:ncol(x$person$rr),
               names_to = "iter",
               values_to = "ref"
          ) %>% 
          mutate(iter = factor(iter, ordered = TRUE, levels = names(pr.long))) %>%  
          mutate(n.iter = as.integer(iter))
     pr.long.e <- filter(pr.long, n.iter %in% c(evens(2:max.iter)))
     pr.long.o <- filter(pr.long, n.iter %in% c(odds(3:max.iter)))
     
     gr.long <- data.frame(x$group$rr)
     gr.long <- gr.long %>% 
          mutate(group = rownames(gr.long)) %>% 
          pivot_longer(
               cols = 1:ncol(x$group$rr),
               names_to = "iter",
               values_to = "ref"
          ) %>% 
          mutate(iter = factor(iter, ordered = TRUE, levels = names(gr.long))) %>%  
          mutate(n.iter = as.integer(iter))
     gr.long.e <- filter(gr.long, n.iter %in% c(evens(2:max.iter)))
     gr.long.o <- filter(gr.long, n.iter %in% c(odds(3:max.iter)))
     return(list(person = list(all = pr.long, even = pr.long.e, odd = pr.long.o), 
                 group = list(all = gr.long, even = gr.long.e, odd = gr.long.o)))
}