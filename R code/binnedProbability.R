#' binnedProbability computes sample and ideal statistics for a vector of
#' discrete probabilities summed over a specified number of samples, over a
#' large number of trials
#'
#' @param probabilities is a vector of discrete probabilities between 0 and 1
#'   whose sum is less than or equal to one
#' @param samples summed over per trial
#' @param trials is the total number of trials
#' @param weights are optional additive weights to create a weighted summed
#'   value of each class (should be the same as the length of the probability
#'   vector)
#' @param added_constant to be added to the weighted sum to create the desired
#'   parameter
#'
#' @returns a list containing: count (matrix of counts for each probability),
#'   the mean value averaged over all the trials, the ideal covariance, the
#'   sample covariance, the weighted value over the probability cases, and the
#'   weighted sample and ideal covariance matrices
require(stats)

binnedProbability <- function(probabilities,samples,trials=1000,weights=0,
                              added_constant=0) {
  cumPct <- c(0,cumsum(probabilities),1)
  sampleMat <- array(as.numeric(apply(t(runif(samples*trials)),c(1),
                cut,breaks=cumPct,labels=1:(length(probabilities)+1))),
                dim=c(samples,trials))
  eventCt <- array(0,dim=c(length(probabilities),trials))
  for (i in (1:length(probabilities))) {
    eventCt[i,] <- colSums(sampleMat==i)
  }
  meanProb <-colMeans(t(eventCt))/samples
  C <- ((eventCt/samples-meanProb) %*% t(eventCt/samples-meanProb))/trials
  idealC <- -probabilities %*% t(probabilities)
  idealC <- (idealC+diag(probabilities*(1-probabilities)-diag(idealC)))/samples
  if (length(weights)==1) {
    result <- list(count=eventCt,sample_mean=colMeans(t(eventCt))/samples,
                   sample_cov=C,ideal_cov=idealC)
  }
  else {
    result <- list(count=eventCt,mean=colMeans(t(eventCt))/samples,
                   sample_cov=C,
                   values=t((weights/samples) %*% eventCt)+added_constant,
                   ideal_cov=idealC,
                   weighted_sample_cov=(weights %*% t(weights)) * C,
                   weighted_ideal_cov=(weights %*% t(weights)) * idealC)
  }
  result
}
