
simulation = function(param, N = 50, tf = 19) # tf = time steps, N = sites, this is simulations of a single species
{
  #Add the matrices of transitions and invariant probabilitiess (see other script)
  #param$p0 = param$PYeq       #to start from the state of equilibrium
  Y = matrix(0,N,tf) # matrix of zeroes with sites in rows and time steps in columns
  X = matrix(0,N,tf)  # matrix of zeroes with sites in rows and time steps in columns (same as above)
  Y[,1] = rbinom(N, 1, param$p0) # return 300 observations with 1 trial each at a probability = p0, this corresponds to the first time step (initial probability of seeds) = Starting state of whether or not that species is present as a seed in that site
  X[,1] = rbinom(N, 1, Y[,1]*param$g) # now return 300 observations with 1 trial each at a probability = whether a seed was present initially, times the probability of that seed germinating? # given that a species is present belowground at the start of the time series, does it germinate/is it present aboveground? starting state of whether that species is present aboveground at that site; note that a species must be present belowground (Y[,1] = 1) for it to germinate; this messes with me a little because these are in the same time step, cant a species be present aboveground without being present belowground? or does this include colonization? are there seeds there and do they germinate?
  for(t in 2:tf) # then for every subsequent time step calculate whether or not it's present  given rates in teh previous timestep
  {
    C = rbinom(N, 1, param$c)           #has this species colonized into the site from outside
    R = rbinom(N, 1, X[,t-1]*param$r)   #Given that a species was present aboveground in the previous time step, did it reproduce
    S = rbinom(N, 1, Y[,t-1]*param$s)   #given that a seed of a species was present belowground in the previous timestep, did it survive
    Y[,t] = 1 - (1 - C)*(1 - R)*(1 - S) # tying it all together, Y (is the species present belowground) is 1 minus the joint probability of total failure i.e. the species failed to colonize from outside the plot, seed survival failed, and there was no existing flora abveground to provide new seeds; 1 minus this gives whether or not the plant was belowground; if colonization OR reproduction OR seed surival occurred in the previous timestep, then there would be seeds in the seedbank in year t 
    X[,t] = rbinom(N,1,Y[,t]*param$g) #germination plants are aboveground in year t if plants are belowground in year t and they germinate; plants can only germinate if tehy are present; presence of existing species aboveground depends only on the seed bank status in the previous year
  }
  return(X)
}

