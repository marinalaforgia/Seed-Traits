
simulation = function(param, N = 50, tf = 19) # tf = time steps, N = number of patches/plots, this is simulates a single species
{
  Y = matrix(0,N,tf) # matrix of zeroes with sites in rows and time steps in columns
  X = matrix(0,N,tf)  # matrix of zeroes with sites in rows and time steps in columns (same as above)
  Y[,1] = rbinom(N, 1, param$p0) 
  X[,1] = rbinom(N, 1, Y[,1]*param$g) 
  for(t in 2:tf) 
  {
    C = rbinom(N, 1, param$c)           #has this species colonized into the site from outside
    R = rbinom(N, 1, X[,t-1]*param$r)   #Given that a species was present aboveground in the previous time step, did it reproduce
    S = rbinom(N, 1, Y[,t-1]*param$s)   #given that a seed of a species was present belowground in the previous timestep, did it survive
    Y[,t] = 1 - (1 - C)*(1 - R)*(1 - S) # tying it all together, Y (is the species present belowground) is 1 minus the joint probability of total failure i.e. the species failed to colonize from outside the plot, seed survival failed, and there was no existing flora aboveground to provide new seeds; 1 minus this gives whether or not the plant was belowground; if colonization OR reproduction OR seed surival occurred in the previous timestep, then there would be seeds in the seedbank in year t 
    X[,t] = rbinom(N,1,Y[,t]*param$g) 
  }
  return(X)
}

