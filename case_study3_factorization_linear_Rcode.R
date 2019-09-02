# Bayesian linear regression


# p number of parameters
# X: n by p matrix
# y: n by 1 matrix
# beta: p by 1 matrix


# construct data by y ~ N(Xb, sigma^2)

n=500
p=10

# generate X
ev <- 1:p # choose some positive eigenvalues
U <- svd(matrix(runif(p*p), nc=p))$u  # an orthogonal matrix
# pdM <- U %*% diag(ev) %*% t(U) 
# x = mvrnorm(n, mu = c(1:p)/2, pdM)

x = matrix(rnorm(p*n),ncol=p)
#x = cbind(rep(1,n),x)
beta = seq(from=-p/2,to=p/2-1, length.out = p)
y = x%*%matrix(beta,ncol=1) + rnorm(n,0,1)


##### helper function #####
get.index.linear = function(x, p){
  m  = dim(x)[2]+1
  total.para  = m + sum(m:(m-p+1))+ m
  mu.index = 1:m
  B.index = (m+1):(m+sum(m:(m-p+1)))
  d.index = (m + sum(m:(m-p+1))+1): total.para
  list(mu.index = mu.index, B.index = B.index, d.index = d.index)
}


generate.theta = function(x,VB.linear,nsample = 2000 ){
  m = dim(x)[2]+1
  nf = VB.linear$nf
  index = get.index.linear(x,nf)
  mu.index = index$mu.index
  B.index = index$B.index
  d.index = index$d.index
  
  mu.VB.linear = tail(VB.linear$result[,mu.index],1)
  B.VB.linear = VectoMat.C(nr =m, nc = nf, tail(VB.linear$result[,B.index],1))
  d.VB.linear = tail(VB.linear$result[,d.index],1)
  
  SIGMA.linear = B.VB.linear%*%t(B.VB.linear) 
  diag(SIGMA.linear) = diag(SIGMA.linear) + d.VB.linear^2
  
  theta.VB.linear = mvrnorm(n=2000, mu = mu.VB.linear, SIGMA.linear)
  return(theta.VB.linear) 
}



library(MASS)

##### MCMC #####

# @input mu0 = p dim vector, Lambda0 = p by p vector
# @output sig2 = scalar, mu = p by 1 matrix
conjugate.prior.linear = function(Lambda0, a0=1, b0=1) {
  sigma2.draw = b0/rgamma(1,shape=a0,scale=1)
  mu0.draw = mvrnorm(n=1,rep(0,dim(Lambda0)[1]),Sigma=solve(Lambda0)*sigma2.draw)
  return(list(sig2 = sigma2.draw, mu = matrix(mu0.draw,ncol=1)))
}


# @posterior calculation
# @output: list of beta and sigma2
posterior.mcmc.linear = function(ndraw,x,y, Lambda0,a0=0.1, b0=0.1){
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
  posterior.mcmc.linear(niter, x,y,Lambda0 = diag(1/100,p))
}




##### VB Linear regression ##### 

# @ theta {beta, log sig2} (p+1) dim
# @ lambda {mu, B, d}, B is a lower triangular matrix 
# @ d diagnol entry of D

# @ input: p, dim of beta
# @ output: a vectorised diagonal lower triangular matrix. The diagonal is log diag(C)
# @ C is T' in the paper
MattoVec.C = function(C){
  C[lower.tri(C,TRUE)]
}

# @input: nr number of rows, nc: number of cols, C a vector
# @output: a lower triangular matrix with nr rows, nc cols
VectoMat.C = function(nr,nc,C){
  ind = matrix(1:(nr*nc),ncol=nc,byrow=FALSE)
  ind[seq(1,nr*nc,1)[lower.tri(ind,TRUE)]] = C
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


ELBO.linear = function(y,x,theta,mu,B, D ,Lambda0 = diag(1/100,dim(x)[2]),a0=0.1,b0=0.1){
  p = length(theta) - 1
  beta = matrix(theta[1:p],ncol=1)

  sig2 =exp(theta[p+1])
  invSIGMA = solve(B%*%t(B) + D%*%D )
  loglik = sum(dnorm(y,mean=x%*%beta,sd = sqrt(sig2),log=TRUE))
  logprior = -(1+a0)*log(sig2) - b0/sig2 + sum(dnorm(beta,mean = rep(0,p),sd = sqrt(1/diag(Lambda0)[0])*sig2),log=TRUE)
  logq = 0.5*log(det(invSIGMA))- 0.5*( t(theta-mu)%*%invSIGMA%*%(theta-mu))
  return(loglik+logprior - logq)
}


VAFC.linear = function(x,y,m = dim(x)[2]+1,zeta = 0.95, ep0 = 1e-6, p = 2){
  # dimension of lambda
  dim.lambda = m + sum(m:(m-p+1))+ m
  
  diff.lambda = rep(NA,dim.lambda)
  
  result.trace = matrix(0,dim.lambda, nrow = 1)
  
  mu.index = 1:m
  
  B.index = (m+1):(m+sum(m:(m-p+1)))
  result.trace[1,B.index] = 0
  
  d.index = (m + sum(m:(m-p+1))+1):dim.lambda
  result.trace[1,d.index] = 1
  
  Eg2 = Edelta2 = rep(0,dim.lambda)
  LB = 0
  
  t = 1
  
  while (any(diff.lambda>0.001) || t<5000){
    # read lastest mu and C
    mu = matrix(t(result.trace[t,mu.index]),ncol=1)
    B.vec = matrix(t(result.trace[t,B.index]),ncol=1)
    B = VectoMat.C(nr = m,nc = p,B.vec)
    d = matrix(t(result.trace[t,d.index]),ncol=1)
    D = diag(c(d),length(d))

    
    # step 1
    ep = matrix(rnorm(m),ncol=1)
    z = matrix(rnorm(p),ncol=1)
    
    theta = mu + B%*%z + d*ep
    
    inv.D2 = diag((1/c(d))^2,m)
    inv.sigma = inv.D2 - inv.D2%*%B%*%solve(diag(p)+t(B)%*%inv.D2%*%B)%*%t(B)%*%inv.D2
    
    # step 2
    g_mu = matrix(gradient.h(y,x,theta),ncol=1)
    
    
    #g_B = g_mu %*% t(z)+ solve(B%*%t(B)+D%*%D)%*%(B%*%z + d*ep)%*%t(z)
    g_B = g_mu %*% t(z)+ inv.sigma%*%(B%*%z + d*ep)%*%t(z)
    g_B = MattoVec.C(g_B)
    
    g_d = diag(g_mu%*%t(ep)+  inv.sigma%*%(B%*%z + d*ep)%*%t(ep))
    #g_d =  diag(g_mu%*%t(ep)+  solve(B%*%t(B)+D%*%D)%*%(B%*%z + d*ep)%*%t(ep))
    
    # step 3
    
    
    ## part for mu
    Eg2[mu.index] = zeta*Eg2[mu.index] + (1-zeta)*t(g_mu^2)
    
    rho_mu = sqrt(Edelta2[mu.index]+ep0)/sqrt(Eg2[mu.index]+ep0)
    
    delta.lamba2.mu = rho_mu * g_mu
    
    mu = mu + delta.lamba2.mu
    
    Edelta2[mu.index] = zeta * Edelta2[mu.index] + (1-zeta)*delta.lamba2.mu^2
    
    ## part for B
    Eg2[B.index] = zeta*Eg2[B.index] + (1-zeta)*g_B^2
    
    rho_B = sqrt(Edelta2[B.index]+ep0)/sqrt(Eg2[B.index]+ep0)
    
    delta.lamba2.B = rho_B * g_B
    
    B.vec = B.vec + delta.lamba2.B
    
    B = VectoMat.C(nr = m, nc =p, B.vec)
    
    Edelta2[B.index] = zeta * Edelta2[B.index] + (1-zeta)*delta.lamba2.B^2

    
    ## part for d
    Eg2[d.index] = zeta*Eg2[d.index] + (1-zeta)*g_d^2
    
    rho_d = sqrt(Edelta2[d.index]+ep0)/sqrt(Eg2[d.index]+ep0)
    
    delta.lamba2.d = rho_d * g_d
    
    diag(D) = d + delta.lamba2.d
    
    d = diag(D)
    
    Edelta2[d.index] = zeta * Edelta2[d.index] + (1-zeta)*delta.lamba2.d^2
    
  
    # step 8: combine results
    result.trace = do.call(rbind, list(result.trace, c(mu,B.vec,d)))
 
    t = t+1
    diff.lambda = diff(result.trace[(t-1):t,])
    LB = c(LB,ELBO.linear(y,x,theta,mu,B,D))
  }
  return(list(result = result.trace,LB=LB,nf = p))
}


#### VB result
VB.linear_f5 = VAFC.linear(x=x,y=y, p = 5)

VB.linear_f2 = VAFC.linear(x=x,y=y, p = 2)
VB.linear_f1 = VAFC.linear(x=x,y=y, p = 1)


tail(VB.linear_f1$LB)

tail(VB.linear_f2$LB)
tail(VB.linear_f5$LB)

#ELBO is almost at the same level, chooose nfactor =2 for demonstration purpose


#### MCMC result #####
mcmc.output = mcmc.linear(y=y,x=x)

##### Analysis of Linear regression ##### 

# transformation of VB result 

theta.f1 = generate.theta(x,VB.linear_f1)
theta.f2 = generate.theta(x,VB.linear_f2)

theta.f5 = generate.theta(x,VB.linear_f5)


### plot function 
pdf("Case3_VBMCMC.factorization_linear.pdf")


### for traceplot, take num factor = 2
index = get.index.linear(x,2)
mu.index = index$mu.index
B.index = index$B.index
d.index = index$d.index
for (i in mu.index){
  plot(VB.linear_f1$result[,i], main = paste("trace plot of mu",i),ylab="",xlab="iter",type="l")
  lines(VB.linear_f2$result[,i],col=2)
  lines(VB.linear_f5$result[,i],col=4)
  legend("topright",col=c(1,2,4), legend = c("1 factor", "2 factor","5 factor"),lty=rep(1,4))
}



plot(VB.linear_f1$LB,main="Lower Bound of VB (based on 1 sample) ",type="l",col=1)
lines(VB.linear_f2$LB,type="l",col=2)
lines(VB.linear_f5$LB,type="l",col=4)
legend("topright",col=c(1,2,4), legend = c("1 factor", "2 factor","5 factor"),lty=rep(1,4))



# variance of VB is too narrow
for (i in 1:(length(mu.index)-1)){
  plot(density(mcmc.output$beta[,i]),col=2,main=paste0("density plot of beta",i),lty=2)
  lines(density(theta.f1[,i]),col=1)
  lines(density(theta.f2[,i]),col=2)
  lines(density(theta.f5[,i]),col=4)
  legend("topright",col=c(2,1,2,4,6), legend = c("MCMC","1 factor", "2 factor","5 factor"),lty=c(2,rep(1,4)))
}


plot(density(mcmc.output$sigma2),col=2,main=expression(paste("density plot of ", sigma^2)),lty=2)
lines(density(exp(theta.f1[,length(mu.index)])),col=1)
lines(density(exp(theta.f2[,length(mu.index)])),col=2)
lines(density(exp(theta.f5[,length(mu.index)])),col=4)
legend("topright",col=c(2,1,2,4), legend = c("MCMC","1 factor", "2 factor","5 factor"),lty=c(2,rep(1,4)))

dev.off()
