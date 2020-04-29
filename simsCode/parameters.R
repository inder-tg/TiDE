
library(stepR)

MSE <- function(estimate, true){
  mean(( estimate - true )^2)
}

createSignalPeak <- function(leftValue, rightValue, cp1, cp2, value) {
  function(t) ifelse(t < cp1, rep(leftValue, length(t)),
                     ifelse(t < cp2, rep(value, length(t)), rep(rightValue, length(t))))
}
