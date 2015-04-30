#' Correlation coefficient for vectors
#'
#' Calculates a correlation coefficient and p-value for two independent sets of vectors.
#'
#' @param W1 The first set of vectors expressed by their scalar components in a two-column matrix (u and v, the change along the x and y axes, respectively) 
#'
#' @param W2 The second set of vectors expressed by their scalar components in a two-column matrix (u and v, the change along the x and y axes, respectively)
#'
#'
#' @export
#' @return The correlation coefficient and a probability value 
#'
#' @examples 
#' x <- matrix(rnorm(60, 2, 1), ncol=2)
#' y <- matrix(rnorm(60, 2, 1), ncol=2)
#'
#' vector_corr(x, y) 
vector_corr <- 
function(W1, W2) {
  #function based on Crosby et al. 1993 J Atm Oce Tech to calculate the correlation between two-dimensional vectors
  #WARNING ABOUT SMALL SAMPLE SIZES BASED ON LENGTH OF VECTOR

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

  print(corr_stat) 
  print(p_value)
}