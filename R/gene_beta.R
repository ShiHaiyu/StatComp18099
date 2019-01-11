#' @title A beta(a,b) sampler using R
#' @description A function to produce n samples of Beta(a,b)
#' @param n  :the number of samples
#' @param a,b  :the param of ditribution Beta(a,b)
#' @return a random sample of size \code{n}
#' @examples
#'
#' x=gene_beta(3,2,1000)
#' hist(x,prob=TRUE,main="1/beta(2,3)*x^2*(1-x)")
#' y=seq(0,1,0.01)
#' lines(y,1/beta(3,2)*y^2*(1-y))
#'
#' @export
gene_beta<-function(a,b,n)
{k=0
y=c()
while(k<=n) #stop until n datas are produced
{u<-runif(1)
x<-runif(1) #random variate from g(x)
if(u<=x^(a-1)*(1-x)^(b-1))  #accpet condition
{k=k+1
y=c(y,x)
}
}
return(y) #output a vector
}

