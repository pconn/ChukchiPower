#simulate density for Bearded and Ribbon seals for Chukchi power analysis
sim_density_Chukchi<-function(){
  load("Knot_cell_distances.Rdata")  
  load("Grid_chukchi.Rda")
  load("Grid_chukchi_ExpDensities.Rda")
  
  Adj=Matrix(Grid.chukchi$Adj)
  diag(Adj)=2
  Adj=Adj/rowSums(Adj)
  D.bearded=as.vector(Grid_chukchi[["D.bearded"]]%*%Adj)
  D.ringed=as.vector(Grid_chukchi[["D.ringed"]]%*%Adj)
 
  S=length(Grid_chukchi)
  Area.adjust=1-Grid_chukchi[["land_cover"]]
  Log.area.adjust=log(Area.adjust)
   
  n.species=2
  n.knots=ncol(K.data$K)
  #generate true abundance by group
  G.true=matrix(0,n.species,S)
  Lambda.true=G.true
  Eta.true=G.true
  Sim.pars=list(tau.epsilon=c(100,100),tau.eta=c(1,1))
  X=matrix(0,n.species,S)
  X[1,]=log(625*(D.bearded+0.0001))
  X[2,]=log(625*(D.ringed+0.0001))
  for(isp in 1:n.species){
    Q=Matrix(Sim.pars$tau.eta[isp]*K.data$Q.knot)
    Alpha=rrw(Q)
    Eta.true[isp,]=K.data$K%*%Alpha
    Lambda.true[isp,]=Area.adjust*exp(X[isp,]+Eta.true[isp,]+rnorm(S,0,1/sqrt(Sim.pars$tau.epsilon[isp])))
  }
  Out=list(Lambda.true=Lambda.true,Eta.true=Eta.true)
}
