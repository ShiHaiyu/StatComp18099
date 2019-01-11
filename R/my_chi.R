#' @title computes the chi-square test statistic for two vectors
#' @description a faster version of chisq.test() that only computes the chi-square test statistic when the input is two numeric vectors with no missing values
#' @param x,y two vectors
#' @return chi-square test statistic
#' @examples
#' \dontrun{
#' x<-runif(10)
#'y<-runif(10)
#'my_chi(x,y)
#'}
#' @export
my_chi<-function(x,y)
{
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  OK <- complete.cases(x, y)
  x <- factor(x[OK])
  y <- factor(y[OK])
  if ((nlevels(x) < 2L) || (nlevels(y) < 2L))
    stop("'x' and 'y' must have at least 2 levels")
  x <- table(x, y)
  p = rep(1/length(x), length(x))
  p <- p/sum(p)
  n<-sum(x)
  E <- n * p
  V <- n * p * (1 - p)
  STATISTIC <- sum((x - E)^2/E)
  STATISTIC
}

