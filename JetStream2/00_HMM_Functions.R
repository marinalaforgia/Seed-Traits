### HMM FUNCTIONS ###

# matrix of Y transitions (Y = seedbank)
# Zt (combined state) is a function of the seedbank (Y) at time t-1 and flora (x) at time t)
calculateQY = function(param) # calculate QY is a function of param; QY transition matrix 
{
  result = matrix(0,2,2) # result is a matrix of zeros with 2 columns and 2 rows
  result[1,1] = 1 - param$c # probability of going from A to A: that plants don't colonize
  result[2,1] = (1 - param$c)*(1 - param$s)*(1 - param$g*param$r) # probability of going from P to A: the second row of the first column is (1-c) * (1-s) * (1-gF) = probability that seeds that don't colonize, dont survive, and don't germinate
  result[,2] = 1 - result[,1] # column 2 is the opposite of column one so that each row sums to 1; probability of going from A to P is the probability that seeds colonize
  # probability of going from P to P is 1 - probability that seeds don't colonize, don't survive and don't germinate 
  return(result) 
}


calculateQZ = function(param) # Qz transition matrix
{
  result = matrix(0,3,3) # 3 by 3 matrix of 0s
  result[1,1] = 1 - param$c # AA to AA: probability seeds don't colonize
  result[2,1] = (1 - param$c)*(1 - param$s) # PA to AA: probability seeds don't colonize, and seeds that are there don't survive
  result[3,1] = (1 - param$c)*(1 - param$s)*(1 - param$r) # PP to AA; If r = 1, then this will always be zero. 
  result[,2] = (1 - result[,1])*(1 - param$g) #[1,2] = AA to PA: (1-g)*c; [2,2] = PA to PA = (1-g)(1-(1-c)(1-s)); [3,2] = PP to PA = (1-g)
  result[,3] = (1 - result[,1])*param$g #[1,3] = AA to PP = gc; [2,3] = PA to PP = g(1-(1-c)(1-s)); [3,3] = PP to PP = g
  return(result)
}


calculatePeq = function(Q)
{
  Vmat = eigen(t(Q))$vectors # extract the eigenvectors from QZ
  #We look for an eigenvector whose components have the same sign and we renormalize it so that the sum makes 1 because all probabilities must sum to 1
  n = dim(Q)[1] # n = number of rows in QZ = number of possible states
  for (i in 1:n)
  {
    v = Vmat[,i] 
    if ( all(sign(v)>=0) | all(sign(v)<=0) ) return(v/sum(v))
  }
}


makeParametersCalculations = function(param)
{
  param$QY = calculateQY(param) # calculate seed bank transitions (QY)
  param$QZ = calculateQZ(param) # calculate state transition probability matrix (QZ)
  param$PZeq = calculatePeq(param$QZ) # calculate the invariant of QZ transition matrix; equilibrium state?
  param$PYeq = param$PZeq[2] + param$PZeq[3] 
  param$pZ0 = c(1-param$p0, param$p0*(1-param$g), param$p0*param$g) # initial distribution vector
  return(param) 
}


pZXprespast = function(X,param)  # forward algorithm for one time series
{
  tf = length(X) # time frame is length of flora
  result = matrix(0,3,tf) # result is a matrix of 0s with three rows (each state: AA, PA, PP) and tf columns
  for (t in 1:tf) # for each column (time step)
  {
    if (t == 1) pZXpast = param$pZ0 else pZXpast = t(param$QZ) %*% result[,t-1] # if it's the first time step, use the initial distribution vector, otherwise transpose the state transition probability matrix and do matrix multiplication on the result matrix from the previous time step
    result[,t] = pZXpast # paste it in for the current time step
    if (X[t] == 1) result[1:2,t] = 0 else result[3,t] = 0 # if flora observed in the current time step, the result for both rows one and two are 0 (because it's observed so it cant be absent belowground) otherwise the third row is 0)
  }
  return(result)
}



pXfuturegivenZ = function(X,param) #backward algorithm for one time series
{
  tf = length(X) # number of time steps
  result = matrix(0,3,tf) # same matrix as before
  result[,tf] = 1 # starting at 1
  for (t in tf:2) 
  {
    PZXgivenOldZ = param$QZ # take the state transition probability matrix
    if (X[t] == 1) PZXgivenOldZ[,1:2] = 0 else PZXgivenOldZ[,3] = 0 # if flora is observed in that time step, the first two columns are 0, otherwise the 3rd column is 0
    result[,t-1] = PZXgivenOldZ %*% result[,t]
  }
  return(result)
}

vectForwardBackward = function(X,param) #log-likelihood and lambda coefficients for one time series; 
{
  tf = length(X) # length of the time series
  pZXpp = pZXprespast(X,param)
  pXfgZ = pXfuturegivenZ(X,param)
  pZX = pZXprespast(X,param)*pXfuturegivenZ(X,param)   #P(Zt,whole vector X)
  pZZX = array(rep(param$QZ,tf-1),c(3,3,tf-1))           #P(Zt,Zt+1,whole vector X)
  for (j in 1:3) pZZX[,j,] = pZZX[,j,]*pZXpp[,1:(tf-1)]
  for (i in 1:3) pZZX[i,,] = pZZX[i,,]*pXfgZ[,2:tf]
  pZZX[,1:2,X[2:tf] == 1] = 0
  pZZX[,3,X[2:tf] == 0] = 0
  likelihood = mean(apply(pZX,2,sum))
  result = list()
  result$ll = log(likelihood)
  result$coeffZ = pZX[,1]/likelihood
  result$coeffQ = apply(pZZX,1:2,sum)/likelihood
  return(result)
}

ForwardBackward = function(X,param) #same function extended for a N*T matrix of time series
{
  d = dim(X)
  if (length(d) == 0) return(vectForwardBackward(X,param)) else
  {
    N = d[1]
    result = list()
    llvect = apply(X,1,function(Xvect) vectForwardBackward(Xvect,param)$ll)  #N-vector,i_th component is the log-likelihood of patch i
    result$ll = sum(llvect)                                                 #sum on the N patches
    coeffZmat = apply(X,1,function(Xvect) vectForwardBackward(Xvect,param)$coeffZ) #3*N-matrix,i-th column is the lambda(z) vector of patch i
    result$coeffZ = apply(coeffZmat,1,sum)    #sum on the N patches
    coeffQarray = array(apply(X,1,function(Xvect) vectForwardBackward(Xvect,param)$coeffQ),c(3,3,N)) #3*3*N-matrix,[,,i] is the lambda(z,z') matrix of patch i
    result$coeffQ = apply(coeffQarray,1:2,sum)       #sum on the N patches
    return(result)
  }
}


logLikelihood = function(X,param) ForwardBackward(X,param)$ll # extract log likelihood from final model


Maximization = function(coeffQ, coeffZ, p0 = 'free', g = 'free', c = 'free', s = 'free', r = 'free'){
  # takes into account cases of non-identifiability
  if (p0 == 'free') p0 = sum(coeffZ[2:3])/sum(coeffZ)
  if (g == 'free') g = (coeffZ[3] + sum(coeffQ[,3]))/(sum(coeffZ[2:3]) + sum(coeffQ[,2:3]))
  if (r == 'free') c2 = sum(coeffQ[3,2:3])/sum(coeffQ[3,])      #calculation of c double prime (first step)
  
  if (s == 'free'){                                             #calculation of c prime (first step)
    c1 = sum(coeffQ[2,2:3])/sum(coeffQ[2,])
    if (r == 'free') if(c2 < c1) r = 0
    if (r == 0) c1 = sum(coeffQ[2:3,2:3])/sum(coeffQ[2:3,])
    }
  
  if (c == 'free'){                                              #estimation of c
    c0 = sum(coeffQ[1,2:3])/sum(coeffQ[1,])                     #calculation of c (first step)
    if(s == 'free') if(c1 < c0) s = 0
    if(s == 0){
      if(r == 0) c0 = sum(coeffQ[,2:3])/sum(coeffQ) else c0 = sum(coeffQ[1:2,2:3])/sum(coeffQ[1:2,])
      c1 = c0
      }
    c = c0
    }
  
  if (s == 'free') if(c1 < c) s = 0 else s = 1 - (1-c1)/(1-c)   #estimation of s
  c1 = 1- (1-c)*(1-s)
  if (r == 'free') if(c2 < c1) r = 0 else r = 1 - (1-c2)/(1-c1) #estimation of r
  
  result = list()
  result$p0 = p0
  result$g = g
  result$c = c
  result$s = s
  result$r = r
  result$c1 = c1
  return(result)
}


EMiteration = function(X, param, p0 = 'free', g = 'free', c = 'free', s = 'free', r = 'free')
{
  FBresult = ForwardBackward(X,param)
  coeffZ = FBresult$coeffZ
  coeffQ = FBresult$coeffQ
  newparam = Maximization(coeffQ,coeffZ,p0,g,c,s,r)
  newparam = makeParametersCalculations(newparam)
  EMresult = list()
  EMresult$ll = FBresult$ll
  EMresult$newparam = newparam
  return(EMresult)
}


EMestimation = function(X, nIterations = 150, precision = 10^(-5), p0 = 'free', g = 'free', c = 'free', s = 'free', r = 'free')
{
  alwaysSeedSurvival = FALSE
  if (s == 1 | c == 1 | (g == 1 & r == 1))
  {
    alwaysSeedSurvival = TRUE
    if (r == 'free') print('r non identifiable')
    r = 1
    if (s == 'free') print('s non identifiable')
    s = 1
  }
  if (p0 == 1 & alwaysSeedSurvival)
  {
    if (c == 'free') print('c non identifiable')
    c = 1
  }
  if (g == 1)
  {
    if (s == 'free' & r == 'free') print('r non identifiable from s')
    r = 0
  }
  lllist = rep(0,nIterations)
  p0list = rep(0,nIterations)
  glist = rep(0,nIterations)
  clist = rep(0,nIterations)
  slist = rep(0,nIterations)
  rlist = rep(0,nIterations)
  param = list()
  if (p0 == 'free') param$p0 = runif(1) else param$p0 = p0
  if (g == 'free') param$g = runif(1) else param$g = g
  if (c == 'free') param$c = runif(1) else param$c = c
  if (s == 'free') param$s = runif(1) else param$s = s
  if (r == 'free') param$r = runif(1) else param$r = r
  param = makeParametersCalculations(param)
  oldLogLikelihood = -Inf
  for (k in 1:nIterations)
  {
    print(t(param[c('p0','g','c','s','r')]))
    EMresult = EMiteration(X,param,p0,g,c,s,r)
    logLikelihood = EMresult$ll
    lllist[k] = logLikelihood
    p0list[k] = param$p0
    glist[k] = param$g
    clist[k] = param$c
    slist[k] = param$s
    rlist[k] = param$r
    print(logLikelihood)
    if(logLikelihood < oldLogLikelihood + precision) break
    print(k)
    oldLogLikelihood = logLikelihood
    param = EMresult$newparam
  }
  EMresult = list()
  EMresult$lllist = lllist
  EMresult$p0list = p0list
  EMresult$glist = glist
  EMresult$clist = clist
  EMresult$slist = slist
  EMresult$rlist = rlist
  EMresult$ll = logLikelihood
  EMresult$param = param
  return(EMresult)
  
}
