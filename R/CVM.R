#' @title  p value of two-sample Cramer-von Mises test for equal distributions
#' @description  two-sample Cramer-von Mises test for equal distributions
#' @param x random samples from F
#' @param y random samples from G
#' @return p value
#' @examples
#' \dontrun{
#' attach(chickwts)
#'x <- sort(as.vector(weight[feed == "soybean"]))
#'y <- sort(as.vector(weight[feed == "linseed"]))
#'CVM(x,y)
#' }
#' @export
CVM<-function(x,y){
  N=999#permutation numbers
  n<-length(x);
  m<-length(y);
  f<-ecdf(x);
  g<-ecdf(y);
  z<-c(x,y);
  t<-numeric(N);
  t0<-m*n/(m+n)^2*(sum((f(x)-g(x))^2)+sum((f(y)-g(y))^2));
  for(i in 1:N){
    k<-sample(length(z),size=n,replace = FALSE);
    x1<-z[k];
    y1<-z[-k];
    f1<-ecdf(x1);
    g1<-ecdf(y1);
    t[i]<-m*n/(m+n)^2*(sum((f1(x)-g1(x))^2)+sum((f1(y)-g1(y))^2))
  }
  p <- mean(c(t0,t) >= t0)
  return(p)
}
