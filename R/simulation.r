#' Simulate a dataset which can be analyzed by SLTCA
#'
#' @description Simulate a dataset with longitudinal observations.
#' @param n Sample size.
#' @author Teng Fei. Email: tfei@emory.edu
#' @return Returns a data frame with 6 longitudinal features y.1 - y.6, including count (y.1 and y.2), binary (y.3 and y.4) and continuous (y.5 and y.6) type. Variable baselinecov is the baseline risk factor of latent classes. Variable latent is the true latent class labels.
#' @references Hart, K.R., Fei, T. and Hanfelt, J.J. (2020), Scalable and robust latent trajectory class analysis using artificial likelihood. Biometrics. Accepted Author Manuscript <doi:10.1111/biom.13366>.
#' @examples
#'
#' dat <- simulation(500)
#'
#' @importFrom stats as.formula binomial coef dist fitted gaussian poisson rbinom rmultinom rpois runif weights
#' @export


simulation <- function(n){
  #require(mvtnorm)
  requireNamespace("mvtnorm")
  output <- simulate(n)
  return(output)
}
