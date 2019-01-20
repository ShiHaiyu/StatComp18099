## ------------------------------------------------------------------------
a<-matrix(1:9,3,3)
b<-matrix(2:4,3,1)
a%*%b

## ------------------------------------------------------------------------
x=rnorm(10)
y=rnorm(10)
plot(x,y,main="example2")

## ------------------------------------------------------------------------
layout(matrix(1:9,3,3))
layout.show(9)


## ------------------------------------------------------------------------
library(ggplot2)
set.seed(123)#guarantee that the program can be repeated with the same result
x <- 0:4#define r.v. x
p <- c(0.1,0.2,0.2,0.2,0.3)#the relative probability of x
cp <- cumsum(p)#cumulative sum of probability
m <- 1000#size of samples
r <- numeric(m)#construct a zero vector with length of 1000
r <- x[findInterval(runif(m),cp)+1]
# Generate 1000 sample from a uniform distribution and find where does every sample fall in the interval of cp.Every interval corresponds a value of random variable x.
table_r <- table(r)
print(table_r)

#theoretical probability
t_p <- p*1000
names(t_p) <- x
print(t_p)#show the theoretical 1000 samples
a <- data.frame(x,freq = c(105,193,205,200,297,100,200,200,200,300),type = rep(c('random sample','theoretical sample'),each = 5))
ggplot(a,aes(x = x ,y =freq ,fill = type))+geom_col(position = 'dodge')

## ------------------------------------------------------------------------
set.seed(300)
beta_arm <- function(n,a,b){
  beta_pdf <- function(x){
    (1/beta(a,b))*(x)^(a-1)*(1-x)^(b-1)#define the pdf of beta distribution
  }
  j <- 0#the time of iteration
  k <- 0#the available sample
  y <- numeric(n)#a zero vector of length n
  while (k < n){
    u <- runif(1)
    j <- j + 1
    x <- runif(1) #random variate from g
    if ((a-1)*x^(a-2) * (b-1)*(1-x)^(b-2) > u) {
      #accept x
      k <- k + 1
      y[k] <- x
    }
  }
  return(y)
}
sample_beta <- beta_arm(1000,3,2)#record the sample
head(sample_beta,20)
sample_beta_theo <- rbeta(1000,3,2)#generate the sample by using the function "rbeta"
data_beta <- data.frame(sample = c(sample_beta,sample_beta_theo),class = rep(c('empirical','theoretical'),each = 1000))#construst a data frame in order to generate a figure
ggplot(data_beta,aes(x = sample,fill = class))+geom_histogram(position="identity", alpha=0.5)

## ------------------------------------------------------------------------
set.seed(100)
n <- 1000
r <- 4
beta <- 2#define the parameters
lambda <- rgamma(n, r, beta)#generate the random sample from the gamma distribution with r=4, beta=2
x <- rexp(n, lambda)
head(x,30)#list the first 30 samples
x_data <- data.frame(x1 = x)
ggplot(x_data,aes(x = x1))+geom_histogram(fill = 'skyblue',alpha =0.7,binwidth = 0.4)

## ------------------------------------------------------------------------
MC_Beta<-function(x)
{
  m<-1e5
  y<-runif(m,min=0,max=x)
  F_beta_hat<-mean(x*y^2*(1-y)^2/beta(3,3))
  return(signif(F_beta_hat,4))
}

## ------------------------------------------------------------------------
pbeta_1<-function(x)
{
  return(pbeta(x,3,3))
}
x=seq(0.1,0.9,0.1)
cdf=sapply(x,MC_Beta)
est_cdf=sapply(x,pbeta_1)
A=signif(rbind(x,cdf,est_cdf),3)
print(A)

## ------------------------------------------------------------------------
sample_Ray<-function(m,sigma)  ##X/2+X'/2
{
  x=runif(m)
  g=sqrt(-2*sigma^2*log(1-x))/2+sqrt(-2*sigma^2*log(x))/2
  return(g)
}
Asample_Ray<-function(m,sigma)  ##X1/2+X2/2
{
  x1=runif(m)
  x2=rinif(m)
  g=sqrt(-2*sigma^2*log(1-x1))/2+sqrt(-2*sigma^2*log(1-x2))/2
  return(g)
}


## ------------------------------------------------------------------------
m=1e5
g<-function(x)
{
  x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)
}

u=runif(m)   ##f_1
x=sqrt(1-2*log(1-u))
fg=g(x)/(exp(1/2)*x*exp(-x^2/2))
MC_1=sd(fg)

x=rnorm(m)  ##f_2
fg=sqrt(2*pi)*g(x)/exp(-x^2/2)
MC_2=sd(fg)
se=c(MC_1,MC_2)
print(se)  ##output the sd of f_1 and f_2

## ------------------------------------------------------------------------
m=1e5
g<-function(x)
{
  x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)
}

u=runif(m)   ##f_1
x=sqrt(1-2*log(1-u))
fg=g(x)/(exp(1/2)*x*exp(-x^2/2))
theta_hat=mean(fg)
print(theta_hat)

## ------------------------------------------------------------------------
n=20
m=1000
Gi_sam1<-numeric(m)
Gi_sam2<-numeric(m)
Gi_sam3<-numeric(m)
for(i in 1:m)
{
  x<-sort(rlnorm(n,0,1))
  y<-sort(runif(n,0,1))
  z<-sort(rbinom(n,1,0.5))
  J=2*seq(1:n)-n-1
  Gi_sam1[i]=(J%*%x)/(n^2*mean(x))
  Gi_sam2[i]=(J%*%y)/(n^2*mean(y))
  Gi_sam3[i]=(J%*%z)/(n^2*mean(z))
}
log_norm<-c(mean(Gi_sam1),median(Gi_sam1),quantile(Gi_sam1,seq(0,1,0.1)))
unif_norm<-c(mean(Gi_sam2),median(Gi_sam1),quantile(Gi_sam2,seq(0,1,0.1)))
Bern_norm<-c(mean(Gi_sam2),median(Gi_sam1),quantile(Gi_sam3,seq(0,1,0.1)))
A=signif(rbind(log_norm,unif_norm,Bern_norm),3)
colnames(A)=c("mean","median",names(quantile(Gi_sam1,seq(0,1,0.1))))
print(A)


## ----echo=FALSE----------------------------------------------------------
hist(Gi_sam1,prob=TRUE,main="replicates for lognorm")
hist(Gi_sam2,prob=TRUE,main="replicates for uniform")
hist(Gi_sam3,prob=TRUE,main="replicates for Bernoulli")

## ------------------------------------------------------------------------
n=30
m=100
alpha=0.025
Gi_sam1=numeric(m)
UCL_low<-numeric(1000)
UCL_upp<-numeric(1000)
Mean<-numeric(1000)
for(k in 1:1000)
{ for(i in 1:m) #generate m replitcates
{
  x<-sort(rlnorm(n,0,1))
  J=2*seq(1:n)-n-1
  Gi_sam1[i]=(J%*%x)/(n^2*mean(x))
}

  UCL_low[k]<-  mean(Gi_sam1)-qnorm(1-alpha)*sd(Gi_sam1)/sqrt(m)
  UCL_upp[k]<-  mean(Gi_sam1)+qnorm(1-alpha)*sd(Gi_sam1)/sqrt(m)
  Mean[k]<-mean(Gi_sam1)
}
fin_mean<-mean(Mean)
k=sum(UCL_low<fin_mean & UCL_upp>fin_mean)
mean1<-k/1000

## ----echo=FALSE----------------------------------------------------------
mean1

## ------------------------------------------------------------------------
library(MASS)
m=1000
n=200
ratio<-matrix(0,3,10)
for(j in 1:10)
{
  mean<-c(0,1)
  sigma=matrix(c(1,j/20,j/20,1),2,2)
  p_pearson<-numeric(m)
  p_Spearman<-numeric(m)
  p_Kendall<-numeric(m)
  for(i in 1:m)
  {
    sample=mvrnorm(n,mean,sigma)
    x=sample[,1]
    y=sample[,2]
    test1<-cor.test(x,y,method="pearson")
    test2<-cor.test(x,y,method="kendall")
    test3<-cor.test(x,y,method="spearman")
    p_pearson[i]<-test1$p.value
    p_Spearman[i]<-test2$p.value
    p_Kendall[i]<-test3$p.value
  }

  ratio[,j]=c(mean(p_pearson<0.05),mean(p_Spearman<0.05),mean(p_Kendall<0.05))
}
rownames(ratio)<-c("pearson","speqrman","Kendall")
print(ratio)
plot((1:10)/20,ratio[1,],type="l",col="RED",main="power of three method")
lines((1:10)/20,ratio[2,],type="l",col="BLACK")
lines((1:10)/20,ratio[3,],type="l",col="YELLOW")

## ------------------------------------------------------------------------
library(MASS)
m=1000
n=200
mean<-c(0,1)
sigma=matrix(c(1,0,0,1),2,2)
p_pearson<-numeric(m)
p_Spearman<-numeric(m)
p_Kendall<-numeric(m)
for(i in 1:m)
{
  sample=mvrnorm(n,mean,sigma)
  x=sample[,1]
  y=sample[,2]
  test1<-cor.test(x,y,method="pearson")
  test2<-cor.test(x,y,method="kendall")
  test3<-cor.test(x,y,method="spearman")
  p_pearson[i]<-test1$p.value
  p_Spearman[i]<-test2$p.value
  p_Kendall[i]<-test3$p.value
}

ratio=c(mean(p_pearson<0.05),mean(p_Spearman<0.05),mean(p_Kendall<0.05))
names(ratio)<-c("peason","Speqrman","Kendall")
print(ratio)

## ------------------------------------------------------------------------
library(bootstrap)
n=length(law$GPA)
theta.hat<-cor(law$LSAT,law$GPA)
theta.jack<-numeric(n)
for(i in 1:n)      #leave one out
{
  theta.jack[i]<-cor(law$LSAT[-i],law$GPA[-i])
}
bias.jack<-(mean(theta.jack)-theta.hat)*(n-1)
se.jack<-sqrt((n-1)/n*sum((theta.jack-mean(theta.jack))^2))
round(c(original=theta.hat,bias=bias.jack,se=se.jack),3)

## ------------------------------------------------------------------------
library(boot)
hour<-aircondit
B=200     #number of repicates
data<-hour
theta.boot<-function(data,i) #compute statisic
{
  mean(data[i,])
}
boot.obj<-boot(data,statistic = theta.boot,R=2000)

## ------------------------------------------------------------------------
print(boot.obj)

## ------------------------------------------------------------------------
print(boot.ci(boot.obj,type=c("basic","norm","perc","bca")))

## ------------------------------------------------------------------------
n=nrow(scor)
theta.an<-function(x) #compute theta
{
  eigen(x)$values[1]/sum(eigen(x)$values)
}
theta.hat1<-theta.an(var(scor))
theta.jack1<-numeric(n)
for(i in 1:n)
{
  theta.jack1[i]<-theta.an(var(scor[-i,]))
}
bias.jack1<-(mean(theta.jack1)-theta.hat1)*(n-1)
se.jack1<-sqrt((n-1)/n*sum((theta.jack1-mean(theta.jack1))^2))
round(c(original=theta.hat1,bias=bias.jack1,se=se.jack1),5)

## ------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n=length(magnetic)
e1<-e2<-e3<-e4<-numeric()
#leave-two-out

for(k in 1:(n-1))
{ l=k+1
y=magnetic[c(-k,-l)]
x=chemical[c(-k,-l)]

J1<-lm(y ~ x)
yhat_a1<-J1$coef[1] + J1$coef[2]*chemical[k]

yhat_b1<-J1$coef[1] + J1$coef[2]*chemical[l]
e1<-c(e1,magnetic[k]-yhat_a1)
e1<-c(e1,magnetic[l]-yhat_b1)

J2<-lm(y ~ x+I(x^2))
yhat_a2<-J2$coef[1] + J2$coef[2]*chemical[k]+J2$coef[3]*chemical[k]^2
yhat_b2<-J2$coef[1] + J2$coef[2]*chemical[l]+J2$coef[3]*chemical[l]^2
e2<-c(e2,magnetic[k]-yhat_a2)
e2<-c(e2,magnetic[l]-yhat_b2)


J3<-lm(log(y) ~ x)
logyhat_a3<-J3$coef[1] + J3$coef[2]*chemical[k]
logyhat_b3<-J3$coef[1] + J3$coef[2]*chemical[l]
e3<-c(e3,magnetic[k]-exp(logyhat_a3))
e3<-c(e3,magnetic[l]-exp(logyhat_b3))

J4<-lm(log(y) ~ log(x))
logyhat_a4<-J4$coef[1] + J4$coef[2]*log(chemical[k])
logyhat_b4<-J4$coef[1] + J4$coef[2]*log(chemical[l])
e4<-c(e4,magnetic[k]-exp(logyhat_a4))
e4<-c(e4,magnetic[l]-exp(logyhat_b4))
}
c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2))



## ------------------------------------------------------------------------
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


## ------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
fin <- CVM(x,y)
fin
detach(chickwts)

## ------------------------------------------------------------------------
set.seed(1)
##install.packages("RANN")
##install.packages("energy")
##install.packages("Ball")
library(RANN)
library(boot)
library(energy)
library(Ball)

m <- 50
k<-3
p<-2
n1 <- n2 <- 50
R<-999
n <- n1+n2
N = c(n1,n2)

Tn <- function(z, ix, sizes,k)  {
  n1 <- sizes[1]
  n2 <- sizes[2]
  n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0)
  z <- z[ix, ]
  NN <- nn2(data=z, k=k+1)
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,mean = 1,sd = 0.6),ncol=p);
  y <- matrix(rnorm(n2*p,mean = 1,sd = 0.85),ncol=p);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed = i*12)$p.value
}

alpha <- 0.1;
pow1 <- colMeans(p.values<alpha)
print(pow1)

## ------------------------------------------------------------------------
set.seed(1)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,mean = 0.4,sd = 0.6),ncol=p);
  y <- matrix(rnorm(n2*p,mean = 0.5,sd = 0.85),ncol=p);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed = i*12)$p.value
}
alpha <- 0.1;
pow2<- colMeans(p.values<alpha)
print(pow2)

## ------------------------------------------------------------------------
set.seed(1)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  # t distribution with 1 df (heavy-tailed distribution)
  x <- matrix(rt(n1*p,df = 1),ncol=p);
  #bimodel distribution (mixture of two normal distributions)
  y <- cbind(rnorm(n2,mean = 0.4),rnorm(n2,mean = 0.5));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed = i*12)$p.value
}
alpha <- 0.1
pow3 <- colMeans(p.values<alpha)
print(pow3)

## ------------------------------------------------------------------------
set.seed(1)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- c(rnorm(n1,mean = 1,sd = 1))
  y <- c(rnorm(n2,mean = 2,sd = 2))
  z <- c(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed = i*12)$p.value
}
alpha <- 0.1
pow4 <- colMeans(p.values<alpha)
print(pow4)

## ------------------------------------------------------------------------
f <- function(x,theta,eta)
{
  stopifnot(theta > 0)
  return(1/(theta*pi*(1 + ((x-eta)/theta)^2)))
}
m=5000
theta=1
eta=0
sigma <- 1
x <- numeric(m)
x[1] <- rnorm(1, 0, sigma)
k=0
u <- runif(m)

for(i in 2:m)
{
  xt <- x[i-1]
  y<- rnorm(1, mean = xt, sd = sigma)
  num <-f(y,theta,eta) *dnorm(xt, mean = y, sd = sigma)
  den <-f(xt, theta, eta) * dnorm(y, mean = xt, sd = sigma)
  if(u[i] <= num/den) x[i] <- y
  else{
    x[i]<-xt
    k=k+1
  }
}
print(k/m)


## ------------------------------------------------------------------------
b0=1001
index=b0:m
y=x[index]
plot(y~index , type="l", main="",ylab="x")

## ------------------------------------------------------------------------
p=seq(0.1,0.9,0.1)
round( rbind(quantile(y,p), qcauchy(p)),3)


## ------------------------------------------------------------------------
hist(y,prob=T,main='',xlab='x',breaks=50)
xarg=seq(min(y),max(y),0.1)
lines(xarg,f(xarg,theta,eta))

## ------------------------------------------------------------------------
set.seed(906)

f <- function( theta,x )  {
  if (theta<0 || theta>1 ) return (0)
  (2+theta)^x[1] * (1-theta)^(x[2]+x[3]) * theta^x[4]
}
xdata = c( 125, 18, 20, 34 )
m = 10000
th = numeric(m)
th[1] = runif(1)
k = 0
u = runif(m)

for (t in 2:m) {
  xt = th[t-1]
  alph = xt/(1-xt)
  y <- rbeta(1, shape1=alph, shape2=1 )
  numer = f( y,xdata ) * dbeta( xt, y/(1-y), 1)
  denom = f( xt,xdata ) * dbeta( y, alph, 1)
  if ( u[t] <= numer/denom )
    th[t] = y else {
      th[t] = th[t-1]
      k = k + 1
    }
}

## ------------------------------------------------------------------------
theta.hat = mean( th[2001:m] )
print(theta.hat)

## ------------------------------------------------------------------------
ft<-function(a,df){
s1=pt(q=sqrt(a^2*(df-1)/(df-a^2)),df=(df-1),lower =   FALSE)
s2=pt(q=sqrt(a^2*df/(df+1-a^2)),df=df,lower= FALSE)
return(s1-s2)
 }
 kt=c(4:25,100,500,1000)
out=numeric(length(kt))
for(i in 1:length(kt))
{
out[i]=uniroot(ft,lower = 0.001,upper =2,df=kt[i])$root
}
rbind(kt,out)

## ------------------------------------------------------------------------
attach(mtcars)
xt<-list(disp,I(1/disp),disp+wt,I(1/disp)+wt)
##for loop
lm1<-list()
for(i in seq_along(xt)){
  lm1[[i]]<-lm(mpg~xt[[i]])
}
lm1
##lapply
lm2<-lapply(xt,function(x) lm(mpg~x))
lm2
detach(mtcars)

## ------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
##lappply
lm3<-lapply(bootstraps, function(x) lm(x$mpg~x$disp))
lm3
##for loop
lm4<-list()
for(i in seq_along(bootstraps)){
  lm4[[i]]<-lm(bootstraps[[i]]$mpg~bootstraps[[i]]$disp)
}
lm4

## ------------------------------------------------------------------------
##Exercise 3
sapply(lm1,function(mod) summary(mod)$r.squared)
sapply(lm2,function(mod) summary(mod)$r.squared)
##Exercise 4
sapply(lm3,function(mod) summary(mod)$r.squared)
sapply(lm4,function(mod) summary(mod)$r.squared)

## ------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

sapply(trials,function(x) x$p.value)

## ------------------------------------------------------------------------
options(warn=-1)
set.seed(123)
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

x<-runif(10)
y<-runif(10)

my_chi(x,y)
chisq.test(x,y)$statistic

library(microbenchmark)
microbenchmark(
  my_chi(x,y),
  chisq.test(x,y)$statistic
)

## ------------------------------------------------------------------------
options(warn=-1)
my_table<-function(x,y){
  args<-list(x,y)
  bin <- 0L
  lens <- NULL
  dims <- integer()
  pd <- 1L
  dn <- NULL
  for (a in args) {

    fact.a <- is.factor(a)
    if (!fact.a) {
      a0 <- a
      a <- factor(a)
    }
    ll <- levels(a)
    a <- as.integer(a)
    nl <- length(ll)
    dims <- c(dims, nl)
    if (prod(dims) > .Machine$integer.max)
      stop("attempt to make a table with >= 2^31 elements")
    dn <- c(dn, list(ll))
    bin <- bin + pd * (a - 1L)
    pd <- pd * nl
  }
  bin <- bin[!is.na(bin)]
  if (length(bin))
    bin <- bin + 1L
  y <- array(tabulate(bin, pd), dims, dimnames = dn)
  class(y) <- "table"
  y
}


library(microbenchmark)
x<-1:4
y<-2:5
my_table(x,y)
table(x,y)
microbenchmark(
  my_table(x,y),
  table(x,y)
)


## ------------------------------------------------------------------------
options(warn=-1)
myy_chi<-function(x,y)
{
  if (length(x) != length(y))
  stop("'x' and 'y' must have the same length")
  OK <- complete.cases(x, y)
  x <- factor(x[OK])
  y <- factor(y[OK])
  if ((nlevels(x) < 2L) || (nlevels(y) < 2L))
  stop("'x' and 'y' must have at least 2 levels")
  x <- my_table(x,y)
  p = rep(1/length(x), length(x))
  p <- p/sum(p)
  n<-sum(x)
  E <- n * p
  V <- n * p * (1 - p)
  STATISTIC <- sum((x - E)^2/E)
  STATISTIC
}

x<-runif(10)
y<-runif(10)

microbenchmark(
myy_chi(x,y),
my_chi(x,y),
chisq.test(x,y)$statistic
)



