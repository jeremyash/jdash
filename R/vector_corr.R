#' Correlation coefficient for vectors
#'
#' Calculates a correlation coefficient and p-value for two independent sets of vectors.
#'
#' @param W1 The first set of vectors expressed by their scalar components in a two-column matrix or data frame (u and v, the change along the x and y axes, respectively) 
#'
#' @param W2 The second set of vectors expressed by their scalar components in a two-column matrix or data frame (u and v, the change along the x and y axes, respectively)
#'
#' @export
#' @return The correlation coefficient and a probability value 
#'
#' @examples 
#' x <- matrix(rnorm(200, 2, 1), ncol=2)
#' y <- matrix(rnorm(200, 3, 1), ncol=2)
#'
#' vector_corr(x, y) 
vector_corr <- 
function(W1, W2) {
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


