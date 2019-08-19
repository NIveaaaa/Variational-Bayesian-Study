
# y ~ bernoulli (theta)

# MCMC to get inference for theta

# VB (reparameterization trick)

# VB standard operation

# VB natural gradient to be filled

# VB ADADELA/ADAM with reparameterization trick

# compute log likelihood of a bernoulli distribution
loglik = function(y,eta){
  -sum(log(1+exp((1-2*y)*eta)))
}

# sigmoid function
sigmoid = function(z){
  1/(1+exp(-z))
}

# proposal for eta (MCMC)
prop.eta= function(eta,sigma = 0.1){
  eta.prop = rnorm(1,eta,sigma)
  return(eta.prop)
}

# prior for eta
prior.eta = function(eta){
  exp(-eta)/(1+exp(-eta))^2
}


# update lambda - fixed learning rate
update.lambda = function(lambda, grad, rho){
  lambda.new = lambda + rho*grad
  return (lambda.new)
}



# MCMC
mcmc=function(y,iter=2000,burn_in=2000){
  eta.mcmc = rep(NA,burn_in+iter)
  eta.init = 0.5
  
  for (i in 1:(burn_in+iter)){
    if (i ==1 ) eta.cur = eta.init
    #propose an eta
    eta.prop = prop.eta(eta.cur)
    
    log_accept = min(0, loglik(y,eta.prop)+log(prior.eta(eta.prop))
                     -loglik(y,eta.cur)-log(prior.eta(eta.prop)))
    if (log_accept>log(runif(1))){
      #accept it
      eta.cur = eta.prop
    }
    eta.mcmc[i] = eta.cur
    
  }
  return(eta.mcmc)
}

# gradient of reparameterization trick
gradient.repar=function(mu,logsig,y){
  ep = rnorm(1)
  sig = exp(logsig)
  theta.hat = mu+sig*ep
  grad = (sum((2*y-1)/(1+exp((2*y-1)*theta.hat)))) +
    (exp(-theta.hat)-1)/(1+exp(-theta.hat)) + (theta.hat-mu)/sig/sig
  grad = grad*c(1,ep*sig)
  
  return(grad)
}


# gradient of standard VB
gradient.q = function(theta.hat,mu,sigma,y){
  sigma2 = sigma*sigma
  a = (theta.hat - mu)/sigma2
  b = -0.5 + 0.5*(theta.hat-mu)^2/sigma2
  prod = loglik(y,theta.hat) + log(prior.eta(theta.hat))-dnorm(theta.hat,mu,sigma,log=TRUE)
  return(c(a,b)*prod)
}


# natural gradient

natural.gradient.q = function(theta.hat,mu,sigma,y){
  grad = gradient.q(theta.hat,mu,sigma,y)
  Fisher = matrix(c(1/sigma/sigma,0,0,0.5),ncol=2,nrow=2)
  natural.gradient = solve(Fisher)%*%grad
  return(natural.gradient)
}

ELBO = function(theta,mu,sd,y){
  loglik(y,theta) + log(prior.eta(theta)) - dnorm(theta,mu,sd,log=TRUE)
}

# Adam
adam.theta = function(theta,grad, t, mt, vt, alpha = 0.001, beta1 = 0.99, beta2 = 0.999,ep = 1e-8){
  grad = -grad
  mt.new = beta1*mt + (1-beta1)*grad
  vt.new = beta2*vt + (1-beta2)*(grad^2)
  mt.hat = mt.new/(1-beta1^(t+1))
  vt.hat = vt.new/(1-beta2^(t+1))
  # print(alpha*mt.hat/(sqrt(vt.hat)+ep))
  return(list(theta=theta-alpha*mt.hat/(sqrt(vt.hat)+ep),mt = mt.new, vt=vt.new))
}


# Adadelta
adadelta.theta = function(x,grad,Egrad2,Edelta.x2,epsilon = 1e-4, rho = 0.9){
  Egrad2 = (1-rho)*grad*grad + rho*Egrad2
  delta.x = -sqrt(Edelta.x2+epsilon)/sqrt(Egrad2+epsilon)*grad
  Edelta.x2 = rho*Edelta.x2 + (1-rho)*delta.x^2
  return(list(x=x+delta.x, Egrad2 = Egrad2,Edelta.x2 = Edelta.x2))
}


# standard VB fails for small S value,
standard.VB = function(y,S=50){
  diff = c(0,0)
  t = 0
  rho = 1/(t+5)/5
  mu.trace = mu = 0
  logsig2.trace = logsig2 = 0
  elbo.trace = mean(ELBO(rnorm(S,mu.trace,exp(logsig2/2)),mu.trace,exp(logsig2/2),y))
  
  while(any(diff>0.01)||t<100){
    mu = tail(mu.trace,1)
    lsigma2 = tail(logsig2.trace,1)
    sigma = exp(lsigma2/2)
    eta = rnorm(S,mu,sigma)
    
    
    out = lapply(1:S,function(x) gradient.q(eta[x],mu,sigma,y))
    out.reorder = lapply(1:2,function (x) unlist(lapply(out, `[[`, x)))
    
    
    #return mean and log sigma2 gradient
    grad.L = unlist(lapply(out.reorder,mean))
    lambda.new = update.lambda(c(mu,lsigma2),grad.L,rho=rho)
    
    eta.new = rnorm(1000,lambda.new[1],exp(lambda.new[2]/2))
    elbo.trace= c(elbo.trace,mean(ELBO(eta.new,lambda.new[1],exp(lambda.new/2),y)))
    diff[1] = abs(lambda.new[1] - tail(mu.trace,1))
    diff[2] = abs(lambda.new[2]- tail(logsig2.trace,1))
    
    mu.trace=c(mu.trace,lambda.new[1])
    logsig2.trace=c(logsig2.trace,lambda.new[2])
    
    t = t+1
    rho = 1/(t+5)/5
  }
  return (list(mu.trace = mu.trace,logsig2.trace = logsig2.trace,elbo.trace = elbo.trace))
}

# natural gradient VB
natural.VB = function(y,S=50){
  diff = c(0,0)
  t = 0
  rho = 1/(t+5)
  mu.trace = mu = 0
  logsig2.trace = logsig2 = 0
  elbo.trace = mean(ELBO(rnorm(S,mu.trace,exp(logsig2/2)),mu.trace,exp(logsig2/2),y))
  
  while(any(diff>0.01)||t<100){
    mu = tail(mu.trace,1)
    lsigma2 = tail(logsig2.trace,1)
    sigma = exp(lsigma2/2)
    eta = rnorm(S,mu,sigma)
    
    out = lapply(1:S,function(x) natural.gradient.q(eta[x],mu,sigma,y))
    out.reorder = lapply(1:2,function (x) unlist(lapply(out, `[[`, x)))
    
    #return mean and log sigma2 gradient
    grad.L = unlist(lapply(out.reorder,mean))
    lambda.new = update.lambda(c(mu,lsigma2),grad.L,rho=rho)
    
    eta.new = rnorm(1000,lambda.new[1],exp(lambda.new[2]/2))
    elbo.trace= c(elbo.trace,mean(ELBO(eta.new,lambda.new[1],exp(lambda.new/2),y)))
    diff[1] = abs(lambda.new[1] - tail(mu.trace,1))
    diff[2] = abs(lambda.new[2]- tail(logsig2.trace,1))
    
    mu.trace=c(mu.trace,lambda.new[1])
    logsig2.trace=c(logsig2.trace,lambda.new[2])
    
    t = t+1
    rho = 1/(t+5)
  }
  return (list(mu.trace = mu.trace,logsig2.trace = logsig2.trace,elbo.trace = elbo.trace))
}


# reparameterization VB
repara.VB = function (y,S=10){
  diff = c(0,0)
  t = 1
  rho = 1/(t+5)
  mu.trace = 0
  logsig.trace = 0
  
  while (any(diff>0.01) || t<100){
    out = lapply(1:S,function(x) gradient.repar(tail(mu.trace,1),tail(logsig.trace,1),y))
    out.reorder = lapply(1:2,function (x) unlist(lapply(out, `[[`, x)))
    gradient = unlist(lapply(out.reorder,mean))
    lambda.new = update.lambda(c(tail(mu.trace,1),tail(logsig.trace,1)),gradient, rho)
    
    diff[1] = abs(lambda.new[1] - tail(mu.trace,1))
    diff[2] = abs(lambda.new[2] - tail(logsig.trace,1))
    
    mu.trace = c(mu.trace,lambda.new[1])
    logsig.trace = c(logsig.trace,lambda.new[2])
    
    t = t+1
    rho = 1/(t+5)
  }
  return (list(mu.trace = mu.trace,logsig.trace = logsig.trace))
}



# reparameterization VB ADAM
repara.VB.adam = function (y,S=10){
  diff = c(0,0)
  t = 0
  mu.trace = 0
  logsig.trace = 0
  mt.old =  0
  vt.old =  0
  elbo.trace = mean(ELBO(rnorm(S,mu.trace,exp(logsig.trace)),mu.trace,exp(logsig.trace),y))
  
  while (any(diff>0.01) || t<4000){
    mu = tail(mu.trace,1)
    logsig = tail(logsig.trace,1)
    out = lapply(1:S,function(x) gradient.repar(tail(mu.trace,1),tail(logsig.trace,1),y))
    out.reorder = lapply(1:2,function (x) unlist(lapply(out, `[[`, x)))
    gradient = unlist(lapply(out.reorder,mean))
    
    adam.result = adam.theta(c(mu,logsig),gradient,t,mt.old,vt.old)
    mt.old = adam.result$mt
    vt.old = adam.result$vt
    
    lambda.new = adam.result$theta
    
    diff[1] = abs(lambda.new[1] - tail(mu.trace,1))
    diff[2] = abs(lambda.new[2] - tail(logsig.trace,1))
    
    eta.new = rnorm(1000,lambda.new[1],exp(lambda.new[2]))
    elbo.trace= c(elbo.trace,mean(ELBO(eta.new,lambda.new[1],exp(lambda.new[2]),y)))
    
    mu.trace = c(mu.trace,lambda.new[1])
    logsig.trace = c(logsig.trace,lambda.new[2])
    
    t = t+1
  }
  return (list(mu.trace = mu.trace,logsig.trace = logsig.trace,elbo.trace = elbo.trace))
}


# reparemeterization VB ADADELTA

repara.VB.adadelta = function (y,S=10){
  diff = c(0,0)
  t = 0
  mu.trace = 0
  logsig.trace = 0
  elbo.trace = mean(ELBO(rnorm(S,mu.trace,exp(logsig.trace)),mu.trace,exp(logsig.trace),y))
  Egrad2 =  0
  Edelta.x2 = c(0,0)
  while (any(diff>0.01) || t<100){
    mu = tail(mu.trace,1)
    logsig = tail(logsig.trace,1)
    out = lapply(1:S,function(x) gradient.repar(tail(mu.trace,1),tail(logsig.trace,1),y))
    out.reorder = lapply(1:2,function (x) unlist(lapply(out, `[[`, x)))
    gradient = unlist(lapply(out.reorder,mean))
    
    adadelta.result = adadelta.theta(c(mu,logsig),-gradient,Egrad2,Edelta.x2)
    Egrad2 = adadelta.result$Egrad2
    Edelta.x2 = adadelta.result$Edelta.x2
    
    lambda.new = adadelta.result$x
    
    # print(lambda.new)
    diff[1] = abs(lambda.new[1] - tail(mu.trace,1))
    diff[2] = abs(lambda.new[2] - tail(logsig.trace,1))
    
    eta.new = rnorm(1000,lambda.new[1],exp(lambda.new[2]))
    elbo.trace= c(elbo.trace,mean(ELBO(eta.new,lambda.new[1],exp(lambda.new[2]),y)))
    
    mu.trace = c(mu.trace,lambda.new[1])
    logsig.trace = c(logsig.trace,lambda.new[2])
    
    t = t+1
  }
  return (list(mu.trace = mu.trace,logsig.trace = logsig.trace,elbo.trace = elbo.trace))
}


# density generation
generate.density = function(mutrace,lsigtrace,n=2000){
  muhat = tail(mutrace,1)
  sigmahat = exp(tail(lsigtrace,1))
  y = rnorm(n,muhat,sigmahat)
  theta = sigmoid(y)
  return(theta)
}


#data
set.seed(12)
y50 = rbinom(50,size=1,prob=0.3)
y200 = rbinom(200,size=1,prob=0.3)

mcmc_n_50 = mcmc(y50)
mcmc_n_200=mcmc(y200)


theta.n50.mcmc = sigmoid(mcmc_n_50)
theta.n200.mcmc = sigmoid(mcmc_n_200)

# reparameterization
repar.n50 = repara.VB(y50,S=10)
repar.n200 = repara.VB(y200,S=10)

theta.n50.repara.vb=generate.density(mutrace=repar.n50$mu.trace,lsigtrace = repar.n50$logsig.trace)
theta.n200.repara.vb=generate.density(mutrace=repar.n200$mu.trace,lsigtrace = repar.n200$logsig.trace)

# standardize VB
std.VB.n50 = standard.VB(y50,S=100)
std.VB.n200 = standard.VB(y200,S=100)

theta.n50.std.VB = generate.density(std.VB.n50$mu.trace,std.VB.n50$logsig2.trace/2)
theta.n200.std.VB = generate.density(std.VB.n200$mu.trace,std.VB.n200$logsig2.trace/2)


#### natural gradient VB
natural.VB.n50 = natural.VB(y50,S=50)
theta.n50.natural.VB = generate.density(natural.VB.n50$mu.trace,natural.VB.n50$logsig2.trace/2)


natural.VB.n200 = natural.VB(y50,S=50)
theta.n200.natural.VB = generate.density(natural.VB.n200$mu.trace,natural.VB.n200$logsig2.trace/2)



#### repara VB + ADAM
repar.n50.adam = repara.VB.adam(y50,S=1)
theta.n50.repara.vb.adam = generate.density(repar.n50.adam$mu.trace,repar.n50.adam$logsig.trace)


repar.n200.adam = repara.VB.adam(y200,S=1)
theta.n200.repara.vb.adam = generate.density(repar.n200.adam$mu.trace,repar.n200.adam$logsig.trace)


#### repara VB + ADADELTA
repar.n50.adadelta = repara.VB.adadelta(y50,S=1)
theta.n50.repara.vb.adadelta = generate.density(repar.n50.adadelta$mu.trace,repar.n50.adadelta$logsig.trace)


repar.n200.adadelta = repara.VB.adadelta(y200,S=1)
theta.n200.repara.vb.adadelta = generate.density(repar.n200.adadelta$mu.trace,repar.n200.adadelta$logsig.trace)


### plot function ---

pdf("VBMCMC.bernoulli.pdf")

par(mfrow=c(2,1))
plot(sigmoid(std.VB.n50[[1]]),ylab="theta",type='l',
     main="standardize VB (S=100, n=50)")
plot(std.VB.n50[[3]],type='l',ylab = "ELBO")

plot(sigmoid(std.VB.n200[[1]]),ylab="theta",type='l',
     main="standardize VB (S=100, n=200)" )
plot(std.VB.n200[[3]],type='l',ylab = "ELBO")

plot(sigmoid(natural.VB.n50[[1]]),ylab="theta",type='l',
     main="natural gradient VB (S=10,n=50)")
plot(natural.VB.n50[[3]],type='l',ylab = "ELBO")

plot(sigmoid(natural.VB.n200[[1]]),ylab="theta",type='l',
     main="natural gradient VB (S=10,n=200)")
plot(natural.VB.n200[[3]],type='l',ylab = "ELBO")


plot(sigmoid(repar.n50.adam$mu.trace),type='l',
     main="reparameterizaiton trick+ ADAM (n=50, S=1)")
plot(repar.n50.adam$elbo.trace,type='l')


plot(sigmoid(repar.n200.adam$mu.trace),type='l',
     main="reparameterizaiton trick+ ADAM (n=50, S=1)")
plot(repar.n200.adam$elbo.trace,type='l')


plot(sigmoid(repar.n50.adadelta$mu.trace),type='l',
     main="reparameterization trick + ADADELTA (n=50, S=1)")
abline(h=sum(y50)/50)
plot(repar.n50.adadelta$elbo.trace,type='l')


plot(sigmoid(repar.n200.adadelta$mu.trace),type='l',
     main="reparameterization trick + ADADELTA (n=200, S=1)")
abline(h=sum(y200)/200)
plot(repar.n200.adadelta$elbo.trace,type='l')


 par(mfrow=c(1,2))
 plot(mcmc_n_50,type='l')
 plot(mcmc_n_200,type='l')

 par(mfrow=c(2,2))
 plot(repar.n50$mu.trace,type='l',ylab=expression(mu),main="reparemeterization VB, trace plot of mu (n=50)")
 plot(exp(repar.n50$logsig.trace)^2,type="l",ylab=expression(sigma^2),main="trace plot of sigma^2 (n=50)")

 plot(repar.n200$mu.trace,type='l',ylab=expression(mu),main="trace plot of mu (n=200)")
 plot(exp(repar.n200$logsig.trace)^2,type="l",ylab=expression(sigma^2),main="trace plot of sigma^2 (n=200)")

 par(mfrow=c(1,1))
#plot the density of theta
 plot(density(theta.n50.repara.vb),xlab=expression(theta),main = "n=50")
 lines(density(theta.n50.mcmc),col=2)
 lines(density(theta.n50.std.VB),col=3)
 lines(density(theta.n50.natural.VB),col=4)
 lines(density(theta.n50.repara.vb.adam),col=5)
 lines(density(theta.n50.repara.vb.adadelta),col=6)
 legend("topright",
        legend = c("VB.repara (S=10)","MCMC","VB.standard (S=100)","VB.natural.gradient (S=50)",
                   "VB.repara.adam (S=1)","VB,repara.adadelta (S=1"),lty=rep(1,6),col=c(1,2,3,4,5,6))
 plot(density(theta.n200.repara.vb),xlab=expression(theta),main = "n=200",xlim=c(0.2,0.4))
 lines(density(theta.n200.mcmc),col=2)
 lines(density(theta.n200.std.VB),col=3)
 lines(density(theta.n200.natural.VB),col=4)
 lines(density(theta.n200.repara.vb.adam),col=5)
 lines(density(theta.n200.repara.vb.adadelta),col=6)
 legend("topright",
        legend = c("VB.repara (S=10)","MCMC","VB.standard (S=100)","VB.natural.gradient (S=50)",
                   "VB.repara.adam (S=1)","VB,repara.adadelta (S=1"),lty=rep(1,6),col=c(1,2,3,4,5,6))
 dev.off()
