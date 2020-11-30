
library(markovchain)
library(microbenchmark)


probs <- runif(2, 0, 0.01)

transition_matrix <- matrix(c(1 - probs[1], probs[1], probs[2], 1 - probs[2]), 2, 2, byrow = TRUE)

chain <- new('markovchain', states = c("0", "1"), 
             transitionMatrix = transition_matrix, byrow = TRUE)

n = 100000000


microbenchmark({transition_matrix^n}, {chain^n})
