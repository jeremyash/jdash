#' Correlation coefficient for vectors
#'
#' Calculates a correlation coefficient for two independent sets of vectors and performs a significance test 
#'
#' @param W1 the first set of vectors expressed by their scalar components in a two-column matrix or data frame (u and v, the change along the x and y axes, respectively) 
#'
#' @param W2 the second set of vectors expressed by their scalar components in a two-column matrix or data frame (u and v, the change along the x and y axes, respectively)
#'
#' @details The correlation coefficient is a generalization of the square of the one-dimensional correaltion coefficient.  The scale of the coefficient is from 0 (no correlation) to 2 (perfect correlation) and the probability value is derived from a chi-square distributoin with 4 degrees of freedom for sample size of 64 or greater.  Smaller sample sizes require modfified sampling distributions (see reference for additional detail).
#'
#' @references Crosby, D.S., Breaker, L.C., & W.H. Gemmill. (1993). A proposed definition for vector correlation in geophysics: theory and application.  Journal of Atmospheric and Oceanic Technology, 10(3), 355-367.
#'
#' @return The correlation coefficient and the probability value (for sample sizes >= 64)
#'
#' @examples 
#' x <- matrix(rnorm(200, 2, 1), ncol=2)
#' y <- matrix(rnorm(200, 3, 1), ncol=2)
#'
#' vector_corr(x, y) 
#' @export
vector_corr <- 
function(W1, W2) {
  if(dim(W1)[2] > 2 | dim(W2)[2] > 2) stop("input data has more than 2 columns")
  if(dim(W1)[1] < 64) warning("no p-value reported due to insufficient sample size")

  u1 <- W1[,1] 
  v1 <- W1[,2] 

  u2 <- W2[,1]  
  v2 <- W2[,2]


  f <- cov(u1, u1) * (cov(u2, u2) * cov(v1, v2)**2 + cov(v2, v2) * cov(v1,u2)**2) + 
      cov(v1, v1) * (cov(u2, u2) * cov(u1, v2)**2 + 
      cov(v2, v2) * cov(u1, u2)**2) + 
      2 * (cov(u1, v1) * cov(u1, v2) * cov(v1, u2) * cov(u2, v2)) + 
      2 * (cov(u1, v1) * cov(u1, u2) * cov(v1, v2) * cov(u2, v2)) - 
      2 * (cov(u1, u1) * cov(v1, u2) * cov(v1, v2) * cov(u2, v2)) - 
      2 * (cov(v1, v1) * cov(u1, u2) * cov(u1, v2) * cov(u2, v2)) - 
      2 * (cov(u2, u2) * cov(u1, v1) * cov(u1, v2) * cov(v1, v2)) - 
      2 * (cov(v2, v2) * cov(u1, v1) * cov(u1, u2) * cov(v1, u2))

  g <- (cov(u1, u1) * cov(v1, v1) - cov(u1, v1)**2) * (cov(u2, u2) * cov(v2, v2) - cov(u2, v2)**2)

  corr_stat <- f/g 
  p_value <- dchisq(corr_stat, 4)

  print(paste(expression("Correlation coefficient = "), corr_stat)) 
  if(dim(W1)[1] >= 64) print(paste(expression("Probability value = "), p_value))
}
