# Bayesian logistic regression - Gaussian Copula

##### helper function #####

generate.theta = function(x,VB,nsample = 2000 ){
  m = dim(x)[2]
  nf = VB$nf
  index = get.index.logistic(x,nf)
  mu.index = index$mu.index
  B.index = index$B.index
  d.index = index$d.index
  gamma.index = index$gamma.index
  
  mu.VB = tail(VB$result[,mu.index],1)
  B.VB= VectoMat.C(nr =m, nc = nf, tail(VB$result[,B.index],1))
  d.VB= tail(VB$result[,d.index],1)
  gamma.VB = tail(VB$result[,gamma.index],1)
  
  SIGMA = B.VB%*%t(B.VB) 
  diag(SIGMA) = diag(SIGMA) + d.VB^2
  
  psi.VB = mvrnorm(n=2000, mu = mu.VB, SIGMA)
  res=apply(psi.VB,1,function(x) sapply(1:m, function(i) t_inv_psi(x[i],gamma.VB[i])  ))
  return(t(res)) 
}


expit = function(eta){
  1/(1+exp(-eta))
}

get.index.logistic = function(x, p){
  # m: number of beta
  m  = dim(x)[2]
  total.para  = m + sum(m:(m-p+1))+ m + m
  mu.index = 1:m
  B.index = (m+1):(m+sum(m:(m-p+1)))
  d.index = (m + sum(m:(m-p+1))+1): (2*m + sum(m:(m-p+1)))
  gamma.index = (2*m + sum(m:(m-p+1))+1) : total.para
  list(mu.index = mu.index, B.index = B.index, d.index = d.index,
       gamma.index = gamma.index, total.para = total.para)
}

# compute d gamma / d eta, where eta = log (gamma / (2-gamma))
dev_logitgamma = function(x){
  2*exp(x)/(1+exp(x))^2
}

# functions of Y-J transformation (table 1 of Smith & Nott 2019)


t_gamma = function(theta, gamma){
  if (theta < 0){
    theta.bar = 1- theta
    output = (theta.bar^(2-gamma)-1)/(gamma-2)
  }else{
    output = ((theta+1)^gamma - 1 )/gamma
  }
  return(output)
}

t_inv_psi = function(psi, gamma){
  if (psi < 0){
    output = 1- (1-psi*(2-gamma))^(1/(2-gamma))
  }else{
    output = (1 + psi*gamma)^(1/gamma) -1
  }
  return (output)
}

dev_theta_t= function( theta, gamma){
  if (theta < 0){
    output = (1-theta)^(1-gamma)
  }else{
    output = (theta+1)^(gamma-1)
  }
  return (output)
}

dev2_theta2_t = function(theta, gamma){
  if (theta < 0){
    theta.bar = 1 - theta
    output = (gamma - 1)*(theta.bar)^(-gamma)
  }else{
    output = (gamma - 1)*(theta + 1)^(gamma - 2)
  }
  return (output)
}

dev_psi_tinv= function(psi, gamma){
  if (psi < 0){
    output = (1 - psi *(2-gamma))^((gamma-1)/(2-gamma))  
  }else{
    output = (1 + psi * gamma)^((1-gamma)/gamma)
  }
  return (output)
}


dev_gamma_t = function(theta, gamma){
  if (theta < 0){
    theta.bar = 1 - theta
    output = ((2-gamma)*theta.bar^(2-gamma)*log(theta.bar)- theta.bar^(2-gamma)+1)/(2-gamma)/(2-gamma)
  }else{
    output = (gamma * (1+theta)^gamma * log(1+theta)- (1+theta)^gamma +1)/gamma/gamma
  }
  return (output)
}

dev_gamma_tinv = function(psi, theta, gamma){
    output = - dev_psi_tinv(psi,gamma) * dev_gamma_t(theta,gamma)
    return(output)
}

dev_gamma_t = function(theta, gamma){
  if (theta < 0){
    theta.bar = 1 - theta
    output = -(theta.bar^(1-gamma))*log(theta.bar)
  }else{
    output = (theta + 1)^(gamma-1)*log(theta+1)
  }
  return(output)
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
    
    log.accept = loglik(beta.prop,y,x)-
      loglik(beta.curr,y,x)
    
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

gradient.q = function(theta, psi,mu,B,D,gamma){
  m = nrow(theta)
  numer = sapply(1:m, function(x) dev2_theta2_t(theta[x],gamma[x]))
  denom = sapply(1:m, function(x) dev_theta_t(theta[x],gamma[x]))
  # replace by woodbury formula later
  invSIGMA = solve(B%*%t(B) + D%*%D )
  dpsi_dtheta = sapply(1:m, function(x) dev_theta_t(theta[x],gamma[x]))
  matrix(-c(invSIGMA%*%(psi-mu))*dpsi_dtheta + numer/denom,ncol=1)
  }

ELBO.logistic = function(y,x,psi,theta,mu,B, D, gamma, mu.hat,sigma.hat, sigmainv.hat){
  p = length(theta)
  
  psi = matrix(psi,ncol=1)
  beta = matrix(theta[1:p],ncol=1)
  
  # in terms of psi
  invSIGMA = solve(B%*%t(B) + D%*%t(D) )

  loglik = loglik(beta,y,x)
  
  logprior = logprior(beta,mu.hat,sigma.hat,sigmainv.hat)
  
  logq = 0.5*log(det(invSIGMA))- 0.5*( t(psi-mu)%*%invSIGMA%*%(psi-mu))+ sum(log(sapply(1:p, function(x) dev_theta_t(theta[x],gamma[x]))))
  return(loglik+logprior - logq)
}


GC.logistic = function(x,y,m = dim(x)[2],zeta = 0.95, ep0 = 1e-6, p = 2){
  # dimension of lambda
  index = get.index.logistic(x, p)
  
  diff.lambda = rep(NA,index$total.para)
  result.trace = matrix(0,ncol=index$total.para, nrow = 1)
  
  mu.index = index$mu.index
  B.index = index$B.index
  d.index = index$d.index
  gamma.index = index$gamma.index
  
  #initialize gamma and d to 1
  result.trace[gamma.index] = 1 
  result.trace[d.index] = 1
  
  Eg2 = Edelta2 = rep(0,index$total.para)
  
  LB = 0
  t = 1
  
  # about prior hyper parameter
  hyper.para = get.hyperpara(y,x)
  mu.hat = matrix(hyper.para$mu.hat,ncol=1)
  sigma.hat = hyper.para$sigma.hat
  sigmainv.hat = hyper.para$sigmainv.hat
  
  while (any(diff.lambda>0.001) || t<5000){
    # read lastest mu and C
    mu = matrix(result.trace[t,mu.index],ncol=1)
    B.vec = matrix(result.trace[t,B.index],ncol=1)
    B = VectoMat.C(nr = m,nc = p,B.vec)
    d = matrix(result.trace[t,d.index],ncol=1)
    D = diag(c(d))
    gamma = matrix(result.trace[t,gamma.index],ncol=1)
    
    # step 1
    ep = matrix(rnorm(m),ncol=1)
    z = matrix(rnorm(p),ncol=1)
    
    psi = mu + B%*%z + d*ep
    theta = sapply(1:m, function(x) t_inv_psi(psi[x],gamma[x]))
    
    # put in matrix form
    theta = matrix(theta,ncol=1)
    
    #dev_psi_tinv_value = sapply(1:m, function(x) dev_psi_tinv(psi[x],gamma[x]))
    
    dev_gamma_tinv_value = dev_psi_tinv_value = rep(NA,m)
    for (i in 1:m){
      dev_psi_tinv_value[i] = dev_psi_tinv(psi[i],gamma[i])
      dev_gamma_tinv_value[i] = dev_gamma_tinv(psi[i],theta[i],gamma[i])
    }
    
    #rewrite this function, sapply can not be used here
    #dev_gamma_tinv_value = sapply(1:m, function(x) dev_gamma_tinv(psi[x],theta[x],gamma[x]))
    
    # step 2
    
    grad.h  = gradient.h(y,x,theta,mu.hat,sigma.hat,sigmainv.hat)
    grad.q = gradient.q(theta,psi,mu,B,D,gamma)
    diff.grad = grad.h - grad.q
    # gradient wrt mu
    g_mu = diff.grad*dev_psi_tinv_value
    
    # gradient wrt B
    g_B = (dev_psi_tinv_value*diff.grad) %*% t(z) 
    g_B = MattoVec.C(g_B)
    # gradient wrt d
    
    g_d = dev_psi_tinv_value*diff.grad *ep
    
    # gradient wrt gamma
    eta = log(gamma/(2-gamma))
    g_eta = diff.grad * dev_gamma_tinv_value * dev_logitgamma(eta)
    
    # step 3 : let's update all para at one time
    # rewrite as follows:
    g_vec = c(g_mu,g_B,g_d,g_eta)
    
    ## 
    Eg2 = zeta*Eg2 + (1-zeta)*g_vec^2
    
    rho = sqrt(Edelta2+ep0)/sqrt(Eg2+ep0)
    
    delta.lambda2 = rho * g_vec
    
    update_para = result.trace[t,] + delta.lambda2
    
    Edelta2 = zeta * Edelta2 + (1-zeta)*delta.lambda2^2
    
    mu = update_para[mu.index]
    B = VectoMat.C(nr = m,nc = p,update_para[B.index])
    D = diag(c(update_para[d.index]))
    gamma = expit(update_para[gamma.index])*2
    update_para[gamma.index] = gamma
    # step 8: combine results
    result.trace = do.call(rbind, list(result.trace, update_para))
 
    t = t+1
    diff.lambda = diff(result.trace[(t-1):t,])
    LB = c(LB,ELBO.logistic(y,x,psi,theta,mu,B,D,gamma,mu.hat,sigma.hat,sigmainv.hat))
  }
  return(list(result = result.trace,LB=LB, nf = p))
}


#### VB result
VB.logistic_f4 = GC.logistic(x=x,y=y, p = 4)
theta.f4 = generate.theta(x,VB.logistic_f4)
mcmc.output = mcmc.logistic(y=y,x=x)


pdf("Case4_VBMCMC.gaussian_copula.pdf")

plot(VB.logistic_f4$LB,type="l", main="ELBO")

par(mfrow=c(2,5),xpd=FALSE)

for (i in 1:dim(x)[2]){
  plot(density(theta.f4[,i]),col=1,main=paste0("beta",i))
  lines(density(mcmc.output[,i]),col=2)
  #legend("topright",legend=c("VB","MCMC"),col=c(1,2),lty=c(1,1))
}

for (i in 1:dim(VB.logistic_f4$result)[2]){
  plot(VB.logistic_f4$result[,i],col=1,main=paste0("para",i),type='l')
  #legend("topright",legend=c("VB","MCMC"),col=c(1,2),lty=c(1,1))
}



dev.off()
