library(coda)


y<-c(16,9,10,13,19,20,18,17,35,55)
c.oth<-c(58,90,48,57,103,57,86,112,273,64)
n<-y+c.oth
m <- length(n)
par(mfrow=c(1,1))
## plot the data
plot(n, y/n, bty="n", main="", xlim=c(60,330),ylim=c(0.1,0.6),
     ylab="proportion of bicycles", xlab="total number of vehicles")
text(n,y/n,labels=as.character(y),pos=4)
points(n[m], y[m], col=2, pch=19)
abline(0, y[m]/n[m], col=2)

## ab.marg.llik:
##
## function for calculating the marginal log-likelihood
## of the alpha and beta parameters to the hierarchical
## beta-binomial model

ab.marg.llik <- function(alpha, beta, y, n)
{
  m <- length(y)
  const <- m*(lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta))
  lprod <- sum(lgamma(alpha + y) + lgamma(beta + n - y) - lgamma(alpha + beta + n))
  return(const + lprod)
}


## ab.lprior:
##
## function for evaluating the log-prior of alpha and
## beta for the hierarchical beta-binomial model

ab.lprior <- function(alpha, beta)
{
  return((alpha + beta)^(-5/2))
}


## set up the MCMC, allocate space, etc
S <- 100000
alpha <- beta <- rep(NA, S)
alpha[1] <- 1
beta[1] <- 1 
lpost <- ab.marg.llik(alpha[1], beta[1], y, n) + ab.lprior(alpha[1], beta[1])
L <- 1
U <- 2

## MCMC rounds for alpha and beta
for(s in 2:S) {
  
  ## propose new alpha and calculate new log posteiror
  aprime <- runif(1, L*alpha[s-1]/U, L*alpha[s-1]*U)
  lprime <- ab.marg.llik(aprime, beta[s-1], y, n) + ab.lprior(aprime, beta[s-1])
  
  ## accept or reject the proposal
  if(runif(1) < exp(lprime - lpost) * alpha[s-1]/aprime) {  ## accept
    alpha[s] <- aprime
    lpost <- lprime
  } else alpha[s] <- alpha[s-1]  ## reject
  
  
  ## propose new BETA and calculate new log posteiror
  bprime <- runif(1, L*beta[s-1]/U, U*beta[s-1]*L)
  lprime <- ab.marg.llik(alpha[s], bprime, y, n) + ab.lprior(alpha[s], bprime)
  
  ## accept or reject the proposal
  if(runif(1) < exp(lprime - lpost) * beta[s-1]/bprime) {  ## accept
    beta[s] <- bprime
    lpost <- lprime
  } else beta[s] <- beta[s-1]  ## reject
  
}


##
## trace plots of the sampled parameters, and ESS
##

## alpha
plot(alpha[seq(1,S,100)], type="l", bty="n", xlab="s", ylab="alpha")
effectiveSize(alpha[-(1:10)])

## beta
plot(beta[seq(1,S,100)], type="l", bty="n", xlab="s", ylab="beta")
effectiveSize(beta[-(1:10)])

## plot the joint samples
plot(alpha[seq(100,S,100)], beta[seq(100,S,100)],
     bty="n", xlab="alpha", ylab="beta")

## plot a useful transformation of the joint samples
plot(log(alpha/beta)[seq(100,S,100)], log(alpha + beta)[seq(100,S,100)], bty="n",
     xlab="log(alpha/beta)", ylab="log(alpha + beta)")


##
## converting samples of alpha and beta into samples of theta
##

## allocate space for the samples
theta <- matrix(NA, nrow=S, ncol=m)

## gather the samples by direct MC
for(j in 1:m) {
  theta[,j] <- rbeta(S, alpha + y[j], beta + n[j] - y[j])
}

## calculate the posterior means of the theta samples


## plot the posterior means to assess shrinkage
plot(theta.mean, y/n, bty="n", ylab="theta-bar", xlab="y/n")
points(theta.mean[m], y[m]/n[m], col="red", pch=19)
abline(0,1)

## plotting the difference between theta-bar
## and y/n as a function of n
plot(n, theta.mean - y/n, bty="n", xlab="n", ylab="theta-bar - y/n")
points(n[m], theta.mean[m] - y[m]/n[m], col="red", pch=19)
abline(h=0)

## comparing the ms trial to the previous 70
theta.pooled <- apply(theta[-(1:100),], 1, mean)
dm <- density(theta[-(1:100),m])
dab <- density((alpha/(alpha+beta))[-(1:100)])
ylim <- range(c(dm$y, dab$y))
plot(dm, lwd=2, bty="n", ylim=ylim, main="",
     xlab="theta")
lines(dab, col=2, lty=2, lwd=2)
legend("topleft", c("rat F344", "a/(a+b)"), col=1:2,
       lwd=2, lty=1:2, bty="n")


## posterior expectations
mean(theta[-(1:100),m])
mean((alpha/beta)[-(1:100)])
mean(theta[-(1:100),m] > (alpha/(alpha+beta))[-(1:100)])