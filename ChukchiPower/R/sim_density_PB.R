#' simulate density for seals and polar bears for Chukchi power analysis
#' @param Grid_chukchi A list object consisting of "Grid" - a SpatialPolygonsDataFrame with housing
#'        different sampling units (grid cells), and "Adj" - an adjacency matrix describing neighborhood structure
#' @return a list giving expected abundance (Lambda.true) and simulated random effects (Eta.true)
#' @export
#' @import Matrix
#' @keywords simulation abundance
#' @author Paul B. Conn
sim_density_PB<-function(Grid_chukchi){
  Adj=Matrix(Grid_chukchi$Adj)
  diag(Adj)=2
  Adj=Adj/apply(Adj,1,'sum')
  D.PB=as.vector(Grid_chukchi$Grid[["D.pbear"]]%*%Adj)
  D.ringed=as.vector(Grid_chukchi$Grid[["D.ringed"]]%*%Adj)
  D.bearded=as.vector(Grid_chukchi$Grid[["D.bearded"]]%*%Adj)
  
  S=length(Grid_chukchi$Grid)
  Area.adjust=1-Grid_chukchi$Grid[["land_cover"]]
  Log.area.adjust=log(Area.adjust)
   
  n.species=3
  n.knots=ncol(Grid_chukchi$K.data$K)
  #generate true abundance by group
  G.true=matrix(0,n.species,S)
  Lambda.true=G.true
  Eta.true=G.true
  Sim.pars=list(tau.epsilon=c(10,10,10),tau.eta=c(5,5,5))
  X=matrix(0,n.species,S)
  X[1,]=log(625*(D.bearded+0.000001))
  X[2,]=log(625*(D.ringed+0.000001))
  X[3,]=log(625*(D.PB+0.000001))
  for(isp in 1:n.species){
    Q=Matrix(Sim.pars$tau.eta[isp]*Grid_chukchi$K.data$Q.knot)
    Alpha=rrw(Q)
    Eta.true[isp,]=Grid_chukchi$K.data$K%*%Alpha
    Lambda.true[isp,]=Area.adjust*exp(X[isp,]+Eta.true[isp,]+rnorm(S,0,1/sqrt(Sim.pars$tau.epsilon[isp])))
  }
  Out=list(Lambda.true=Lambda.true,Eta.true=Eta.true)
}
