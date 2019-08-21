# Bayesian linear regression


# p number of parameters
# X: n by p matrix
# y: n by 1 matrix
# beta: p by 1 matrix

# Data generation

# construct data by y ~ N(Xb, sigma^2)

n=500
p=5

# generate X
ev <- 1:p # choose some positive eigenvalues
U <- svd(matrix(runif(p*p), nc=p))$u  # an orthogonal matrix
pdM <- U %*% diag(ev) %*% t(U) 
#x = mvrnorm(n, mu = c(1:p)/2, pdM)

x = matrix(rnorm(p*n),ncol=p)
#x = cbind(rep(1,n),x)
beta = -2:2
y = x%*%matrix(beta,ncol=1) + rnorm(n,0,1)

library(MASS)

##### MCMC #####

# @input mu0 = p dim vector, Lambda0 = p by p vector
# @output sig2 = scalar, mu = p by 1 matrix
conjugate.prior = function(Lambda0, a0=1, b0=1) {
  sigma2.draw = b0/rgamma(1,shape=a0,scale=1)
  mu0.draw = mvrnorm(n=1,rep(0,dim(Lambda0)[1]),Sigma=solve(Lambda0)*sigma2.draw)
  return(list(sig2 = sigma2.draw, mu = matrix(mu0.draw,ncol=1)))
}


# @posterior calculation
# @output: list of beta and sigma2
posterior.mcmc = function(ndraw,x,y, Lambda0,a0=0.1, b0=0.1){
  n = length(y)
  an = a0 + n/2
  p = dim(x)[2]
  Lambdan = t(x)%*%x + Lambda0
  mun = solve((t(x)%*%x+Lambda0))%*%(t(x)%*%y)
  bn = as.numeric(b0+0.5*(sum(y^2)-t(mun)%*%Lambdan%*%mun))
  sigma2.pos = bn/rgamma(ndraw,an,scale = 1)
  beta.pos = unlist(lapply(1:ndraw,function(x) mvrnorm(1,mun,sigma2.pos[x]*solve(Lambdan))))
  beta.pos = matrix(beta.pos,ncol=p,byrow = TRUE)
  return(list(beta = beta.pos, sigma2 = sigma2.pos))
}

mcmc.linear = function(niter=2000, y, x ){
  p = dim(x)[2]
  posterior.mcmc(niter, x,y,Lambda0 = diag(1/100,p))
}


##### VB ##### 

# @ theta {beta, log sig2} (p+1) dim
# @ lambda {mu, C^{-T}C^{-1}}, C is a lower triangular matrix with positive
# @ diagnal entry

# @ input: p, dim of beta
# @ output: a vectorised diagonal lower triangular matrix. The diagonal is log diag(C)
# @ C is T' in the paper
MattoVec.C = function(C){
  C[lower.tri(C,TRUE)]
}

# @input: p, number of dimension of beta, dim(C) = p+1
# @output: a lower triangular matrix 
VectoMat.C = function(p,C){
  ind = matrix(1:(p+1)^2,ncol=p+1,byrow=FALSE)
  ind[seq(1,(p+1)^2,1)[lower.tri(ind,TRUE)]] = C
  ind[upper.tri(ind,FALSE)] = 0
  return(ind)
}

gradient.h = function(y,x,theta,Lambda0 = diag(1/100,dim(x)[2]),a0=0.1,b0=0.1){
  n = length(y)
  p = length(theta)-1
  beta = matrix(theta[1:p],ncol=1)
  sig2 = exp(theta[p+1])
  grad.beta =  matrix(apply(sweep(x,1,(y - x%*%beta),"*"),2,sum),ncol=1)/sig2 - beta*diag(Lambda0)/sig2
  grad.lsig2 = -n/2 + sum((y - x%*%beta)^2)/(2*sig2) -p/2+sum(beta^2*diag(Lambda0))/sig2/2 - (a0+1) + b0/sig2
  c(grad.beta,grad.lsig2)
}


ELBO = function(y,x,theta,mu,C,Lambda0 = diag(1/100,dim(x)[2]),a0=0.1,b0=0.1){
  p = length(theta) - 1
  beta = matrix(theta[1:p],ncol=1)

  sig2 =exp(theta[p+1])
  invSIGMA = C%*%t(C)
  loglik = sum(dnorm(y,mean=x%*%beta,sd = sqrt(sig2),log=TRUE))
  logprior = -(1+a0)*log(sig2) - b0/sig2 + sum(dnorm(beta,mean = rep(0,p),sd = sqrt(1/diag(Lambda0)[0])*sig2),log=TRUE)
  logq = 0.5*log(det(invSIGMA))- 0.5*( t(theta-mu)%*%invSIGMA%*%(theta-mu))
  return(loglik+logprior - logq)
}


VB.linear = function(x,y,p = dim(x)[2],rho = 0.9, ep = 1e-6){
  # dimension of lambda
  dim.lambda = p+1 + (1+p+1)*(p+1)/2
  
  diff.lambda = rep(NA,dim.lambda)
  
  mu.trace = matrix(0,ncol=p+1,nrow=1)
  
  C = diag(1,p+1)
  Cp = diag(0,p+1)
  C.trace= matrix(MattoVec.C(C),nrow=1)
  
  Eg.mu2 = Edelta.mu2 = Eg_Cp2 = Edelta.Cp2 = 0  
  LB = 0
  
  t = 1
  
  while (any(diff.lambda>0.001) || t<5000){
    # read lastest mu and C
    mu = matrix(t(mu.trace[t,]),ncol=1)
    C.vec = C.trace[t,]
    C = VectoMat.C(p,C.vec)
    #Cp = C
    #diag(Cp) = log(diag(C))
    # step 1
    s = matrix(rnorm(p+1),ncol=1)
    
    # step 2
    theta = mu + t(solve(C))%*%s
    
    # step 3
    g_mu = gradient.h(y,x,theta) + C%*%s
    
    # step 4
    Eg.mu2 = rho * Eg.mu2 + (1-rho)*g_mu^2
    
    # step 5
    delta.mu = sqrt(Edelta.mu2+ep)/sqrt(Eg.mu2+ep)*g_mu
    
    # step 6
    Edelta.mu2 = rho * Edelta.mu2 + (1-rho)*delta.mu^2
    
    # step 7
    mu = mu + delta.mu
    
     
    # step 8
    g_Cp = -solve(t(C))%*%s%*%t(solve(C)%*%g_mu)
    
    # step 9
    diag(g_Cp) = diag(g_Cp)*diag(C)
    
    # step 10
    Eg_Cp2 = rho * Eg_Cp2 + (1-rho)*MattoVec.C(g_Cp)^2
    
    # step 11
    delta.Cp = sqrt(Edelta.Cp2 + ep)/sqrt(Eg_Cp2 + ep) * MattoVec.C(g_Cp)
    
    # step 12
    Edelta.Cp2 = rho*Edelta.Cp2 + (1-rho)*delta.Cp^2
    
    # step 13
    Cp = MattoVec.C(Cp) + delta.Cp
    Cp = VectoMat.C(p,Cp)
    
    # step 14
    C = Cp
    
    
    # step 15
    diag(C) = exp(diag(Cp))
    
    # combine results
    mu.trace = rbind(mu.trace,c(mu))
    C.trace = rbind(C.trace,MattoVec.C(C))
    
    t = t+1
  
    diff.lambda = c(diff(mu.trace[(t-1):t,]),diff(C.trace[(t-1):t,]))
    LB = c(LB,ELBO(y,x,theta,mu,C))
  }
  return(list(mu=mu.trace,C = C.trace,LB=LB))
}


#### VB result
VB.out = VB.linear(x=x,y=y)


#### MCMC result
mcmc.output = mcmc.linear(y=y,x=x)

#### result transformation

mu.VB = tail(VB.out$mu,1)
C.VB = VectoMat.C(p,tail(VB.out$C,1))
SIGMA.VB = solve(C.VB%*%t(C.VB))

s = t(mvrnorm(2000,rep(0,p+1),diag(1,p+1)))

# theta ~ N (mu, inv(CC'))
theta.VB = solve(t(C.VB))%*%s 
theta.VB = sweep(theta.VB,1,mu.VB,"+")

# mu of beta
mubeta.VB = mu.VB[1:p]

# covariance matrix of beta
covbeta.VB = cov(t(theta.VB[1:p,]))

# sigma (not need that)
sigma2.VB = exp(mu.VB[p+1]+diag(SIGMA.VB)[p+1]/2)
var.sigma2.VB = sigma2.VB^2 * (exp(diag(SIGMA.VB)[p+1])-1)

### plot function 
pdf("VBMCMC.linear_regression.pdf")

for (i in 1:(p+1)){
plot(VB.out$mu[,i], main = paste0("trace plot of mu",i),ylab="",xlab="iter",type="l")
}

plot(VB.out$LB,main="Lower Bound of VB (based on 1 sample",type="l")

for (i in 1:p){
  plot(density(mcmc.output$beta[,i]),col=2,main=paste0("density plot of beta",i))
  lines(density(theta.VB[i,]))
  legend("topright",legend=c("MCMC","VB"),col=c(2,1),lty=c(1,1))
}


plot(density(mcmc.output$sigma2),col=2,main=expression(paste("density plot of ", sigma^2)))
lines(density(exp(theta.VB[p+1,])))

dev.off()
