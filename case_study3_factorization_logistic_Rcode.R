# Bayesian logistic regression - factorization




##### helper function #####
expit = function(eta){
  1/(1+exp(-eta))
}

get.index.logistic = function(x, p){
  m  = dim(x)[2]
  total.para  = m + sum(m:(m-p+1))+ m
  mu.index = 1:m
  B.index = (m+1):(m+sum(m:(m-p+1)))
  d.index = (m + sum(m:(m-p+1))+1): total.para
  list(mu.index = mu.index, B.index = B.index, d.index = d.index)
}



##### data generation #####

# p number of parameters
# X: n by p matrix
# y: n by 1 matrix
# beta: p by 1 matrix

# construct data by logit(y) = Xb

n=500
p=10

# generate X
# ev <- 1:p # choose some positive eigenvalues
# U <- svd(matrix(runif(p*p), nc=p))$u  # an orthogonal matrix
# pdM <- U %*% diag(ev) %*% t(U) 
# x = mvrnorm(n, mu = c(1:p)/2, pdM)

x = matrix(rnorm(p*n),ncol=p)

beta.true = seq(from=-p/2,to=p/2-1, length.out = p)/4
Pry  = expit(x%*%matrix(beta.true,ncol=1))

y = matrix(rbinom(n,1, Pry),ncol=1)




library(MASS)

##### MCMC : independent MH #####


loglik = function(beta, y ,x){
  t(y)%*%x%*%beta - sum(log(1+ exp(x%*%beta)))
}


get.hyperpara = function(y, x){
  glm.fit = glm(y ~ x+0, family = binomial)
  beta.mle =  matrix(glm.fit$coefficients,ncol=1)
  k = nrow(beta.mle)
  sigmainv.hat = diag(k)
  scalor = exp(x%*%beta.mle)/(1+exp(x%*%beta.mle))^2
  for ( i in 1:dim(x)[1]){
    sigmainv.hat = sigmainv.hat + x[i,]%*%t(x[i,])*scalor[i]
  }
  sigma.hat = solve(sigmainv.hat)
  return (list(mu.hat = beta.mle, sigma.hat = sigma.hat, sigmainv.hat=sigmainv.hat))
}

logprior = function(beta, mu.hat, sigma.hat, sigmainv.hat){
  k = nrow(beta)
  -0.5*k*log(2*pi)- 0.5* log(det(sigma.hat))- 0.5*t(beta-mu.hat)%*%sigmainv.hat%*%(beta-mu.hat)
}





# @posterior calculation
# @output: list of beta 
# @ assume prior of beta ~ N(0, Ip)

mcmc.logistic = function(niter=2000,burn_in = 2000, y, x ){
  i = 1
  p = dim(x)[2]
  n = dim(x)[1]
  beta.matrix = matrix(0,nrow = niter+burn_in, ncol = p)
  beta.curr = matrix(rnorm(p),ncol=1)
  
  hyper.para = get.hyperpara(y,x)
  mu.hat = hyper.para$mu.hat
  sigma.hat = hyper.para$sigma.hat
  sigmainv.hat = hyper.para$sigmainv.hat
   
  while (i < niter+burn_in){
  
    beta.prop = mvrnorm(n=1,mu=mu.hat, Sigma =  sigma.hat)
    
    log.accept = loglik(beta.prop,y,x)+
      logprior(beta.curr,mu.hat, sigma.hat,sigmainv.hat)-
      loglik(beta.curr,y,x) - logprior(beta.prop,mu.hat,sigma.hat,sigmainv.hat )
    
    if (min(1,exp(log.accept))>runif(1)){
      beta.curr = beta.prop
    }
    beta.matrix[i,] = beta.curr
    i = i+1
  }
  return(beta.matrix)
  
}



##### VB logistic regression ##### 

# @ theta {beta} (p) dim
# @ lambda {mu, B, d}, B is a lower triangular matrix 
# @ d diagnol entry of D

# @ input: p, dim of beta
# @ output: a vectorised diagonal lower triangular matrix. 
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

# output: return a p by 1 matrix
gradient.h = function(y,x,theta,mu.hat,sigma.hat,sigmainv.hat){
  p = length(theta)
  beta.est = matrix(theta[1:p],ncol=1)
  # grad.theta = rep(0,p)
  # for (i in 1:nrow(y)){
  #   grad.theta = grad.theta+ c(y[i,])*x[i,] - c(exp(x[i,]%*%beta.est)/(1+exp(x[i,]%*%beta.est))) * x[i,]
  # }
  # grad.theta = grad.theta  -sigmainv.hat%*%(beta.est-mu.hat)

  grad.theta = matrix(apply(sweep(x,1,-exp(x%*%beta.est)/(1+exp(x%*%beta.est)),"*") ,2,sum),ncol=1) + matrix(t(y)%*%x,ncol=1)-sigmainv.hat%*%(beta.est-mu.hat)
  return(grad.theta)
}


ELBO.logistic = function(y,x,theta,mu,B, D ,mu.hat,sigma.hat, sigmainv.hat){
  p = length(theta)
  beta = matrix(theta[1:p],ncol=1)
  
  invSIGMA = solve(B%*%t(B) + D%*%D )
  loglik = loglik(beta,y,x)
  logprior = logprior(beta,mu.hat,sigma.hat,sigmainv.hat)
  logq = 0.5*log(det(invSIGMA))- 0.5*( t(theta-mu)%*%invSIGMA%*%(theta-mu))
  return(loglik+logprior - logq)
}


VAFC.logistic = function(x,y,m = dim(x)[2],zeta = 0.95, ep0 = 1e-6, p = 2){
  # dimension of lambda
  dim.lambda = m + sum(m:(m-p+1))+ m
  
  diff.lambda = rep(NA,dim.lambda)
  
  result.trace = matrix(0,dim.lambda, nrow = 1)
  
  mu.index = 1:m

  B.index = (m+1):(m+sum(m:(m-p+1)))
  
  d.index = (m + sum(m:(m-p+1))+1):dim.lambda
  result.trace[d.index] = 1
  Eg2 = Edelta2 = rep(0,dim.lambda)
  LB = 0
  
  t = 1
  
  hyper.para = get.hyperpara(y,x)
  mu.hat = hyper.para$mu.hat
  sigma.hat = hyper.para$sigma.hat
  sigmainv.hat = hyper.para$sigmainv.hat
  
  while (any(diff.lambda>0.001) || t<5000){
    # read lastest mu and C
    mu = matrix(t(result.trace[t,mu.index]),ncol=1)
    B.vec = matrix(t(result.trace[t,B.index]),ncol=1)
    B = VectoMat.C(nr = m,nc = p,B.vec)
    d = matrix(t(result.trace[t,d.index]),ncol=1)
    D = diag(c(d))
    
    # step 1
    ep = matrix(rnorm(m),ncol=1)
    z = matrix(rnorm(p),ncol=1)
    
    theta = mu + B%*%z + d*ep
    
    # step 2
    g_mu = matrix(gradient.h(y,x,theta,mu.hat,sigma.hat,sigmainv.hat),ncol=1)
    
    g_B = g_mu %*% t(z)+ solve(B%*%t(B)+D%*%D)%*%(B%*%z + d*ep)%*%t(z)
    
    
    g_d =  diag(g_mu%*%t(ep)+  solve(B%*%t(B)+D%*%D)%*%(B%*%z + d*ep)%*%t(ep))
    
    # step 3
    
    
    ## part for mu
    Eg2[mu.index] = zeta*Eg2[mu.index] + (1-zeta)*t(g_mu^2)
    
    rho_mu = sqrt(Edelta2[mu.index]+ep0)/sqrt(Eg2[mu.index]+ep0)
    
    delta.lamba2.mu = rho_mu * g_mu
    
    mu = mu + delta.lamba2.mu
    
    Edelta2[mu.index] = zeta * Edelta2[mu.index] + (1-zeta)*delta.lamba2.mu^2
    
    ## part for B
    Eg2[B.index] = zeta*Eg2[B.index] + (1-zeta)*MattoVec.C(g_B)^2
    
    rho_B = sqrt(Edelta2[B.index]+ep0)/sqrt(Eg2[B.index]+ep0)
    
    delta.lamba2.B = rho_B * MattoVec.C(g_B)
    
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
    LB = c(LB,ELBO.logistic(y,x,theta,mu,B,D,mu.hat,sigma.hat,sigmainv.hat))
  }
  return(list(result = result.trace,LB=LB))
}


#### VB result
VB.logistic_f1 = VAFC.logistic(x=x,y=y, p = 1)
VB.logistic_f2 = VAFC.logistic(x=x,y=y, p = 2)
VB.logistic_f3 = VAFC.logistic(x=x,y=y, p = 3)
VB.logistic_f4 = VAFC.logistic(x=x,y=y, p = 4)

tail(VB.logistic_f1$LB)
tail(VB.logistic_f2$LB)
tail(VB.logistic_f3$LB)
tail(VB.logistic_f4$LB)



VB.logsitic = VB.logistic_f1

#ELBO is best in terms of 1 factor
numfactor = 1

#### MCMC result
mcmc.output = mcmc.logistic(y=y,x=x)

##### Analysis of Logistic regression ##### 

# transformation of VB result 
transform.result = function(){}

m = dim(x)[2]
index = get.index.logistic(x=x,p=numfactor)
mu.index = index$mu.index
B.index = index$B.index
d.index = index$d.index

mu.VB.logistic = tail(VB.logsitic$result[,mu.index],1)
B.VB.logistic = VectoMat.C(nr =m, nc = numfactor, tail(VB.logsitic$result[,B.index],1))
d.VB.logistic = tail(VB.logsitic$result[,d.index],1)


SIGMA.logistic = B.VB.logistic%*%t(B.VB.logistic) 
diag(SIGMA.logistic) = diag(SIGMA.logistic) + d.VB.logistic^2

theta.VB.logistic = mvrnorm(n=2000, mu = mu.VB.logistic, SIGMA.logistic)

### plot function 
pdf("Case3_VBMCMC.factorization_logistic.pdf")

for (i in mu.index){
plot(VB.logsitic$result[,i], main = paste("trace plot of mu",i),ylab="",xlab="iter",type="l")
}

plot(VB.logsitic$LB,main="Lower Bound of VB (based on 1 sample",type="l")

# variance of VB is too narrow
for (i in 1:m){
  plot(density(mcmc.output[,i]),col=2,main=paste0("density plot of beta",i))
  lines(density(theta.VB.logistic[,i]))
  legend("topright",legend=c("MCMC","VB"),col=c(2,1),lty=c(1,1))
}


dev.off()
