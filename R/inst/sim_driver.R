#master simulation driver for Chukchi power analysis
library(Matrix)
source('../hierarchicalDS/hierarchicalDS/R/spat_funcs.R')
source('./R/R/simulate_density_Chukchi.R')
source('./R/R/simulate_counts_Chukchi.R')
source('../BOSSst/BOSSst/R/util_funcs.R')
source('../BOSS/BOSS/R/hierarchical_boss.R')
source('../BOSS/BOSS/R/mcmc_boss.R')
source('../BOSS/BOSS/R/util_funcs.R')

load("Grid_chukchi.Rda")
load("Grid_chukchi_ExpDensities.Rda")
#load("Knot_cell_distances.Rdata")  
load("Haulout_samples.Rdat")

Adj=Matrix(Grid.chukchi$Adj)

n.sims=10
n.species=3  #the 3rd is for anomalies/other
SIM=TRUE  #if true, simulate each set of true abundance data; otherwise load from memory
set.seed(20935)
S=length(Grid_chukchi)

if(SIM){
  for(isim in 1:n.sims){
    Sim.abund=sim_density_Chukchi()
    cur.file=paste("./sim_data/Lambda_chukchi_",isim,'.Rda',sep='')
    save(Sim.abund,file=cur.file)
  }
}

#take a look of some plots of 'noisy' and expected seal densities
Temp=vector("list",1)
Temp[[1]]=Grid_chukchi
plot_N_map(1,matrix(Sim.abund$Lambda.true[1,],S,1),Grid=Temp)

Temp=vector("list",1)
Temp[[1]]=Grid_chukchi
plot_N_map(1,matrix(Grid_chukchi[["D.ringed"]],S,1),Grid=Temp)

##generate posterior thinning samples
Thin=array(0,dim=c(3,1,1,1000))
P=rbeta(1000,67,5)
Thin[1,1,1,]=Haulout.samples[[1]][29,23,]
Thin[2,1,1,]=Thin[1,1,1,]*(0.65/mean(Thin[1,1,1,])) #rescale ringed to have a mean of 0.65
Thin=Thin*P
Thin[3,1,1,]=1
Thin.sim=apply(Thin,1,'mean')
rm(Haulout.samples)

#formulate Psi matrix for estimation (assume not species misID)
Psi=array(0,dim=c(3,3,1))
Psi[,,1]=diag(3)

#simulation independent MCMC options
ZIP=FALSE
misID=FALSE
spat.ind=FALSE
n.obs.cov=0
Obs.cov=NULL
Hab.cov=Grid_chukchi@data[,3:5]
Hab.cov$sqrt.mainland=sqrt(Grid_chukchi[["dist_mainland"]])
Hab.cov$EastNorth=Hab.cov$Easting*Hab.cov$Northing
Hab.pois.formula=vector("list",3)
Hab.bern.formula=NULL
for(isp in 1:2)Hab.pois.formula[[isp]]=~dist_mainland+sqrt.mainland+Easting+Northing+EastNorth
Hab.pois.formula[[3]]=~1
Cov.prior.parms=array(0,dim=c(n.species,2,1))
Cov.prior.parms[,1,1]=0.1  
Cov.prior.parms[,2,1]=0.1
Cov.prior.fixed=matrix(0,n.species,dim(Cov.prior.parms)[3])
Cov.prior.pdf=matrix(0,n.species,1)  
Cov.prior.pdf[,1]="pois1"  #model group size as a zero truncated poisson
Cov.prior.n=matrix(2,n.species,1)
fix.tau.nu=FALSE
srr=TRUE
srr.tol=0.5
grps=TRUE
post.loss=FALSE
Control=list(iter=100100,burnin=100,thin=50,MH.nu=matrix(.2,n.species,S),adapt=2000)
Inits=list(tau.nu=rep(100,n.species))
Prior.pars=list(a.eta=1,b.eta=.01,a.nu=1,b.nu=.01,beta.tau=0.01) #(1,.01) prior makes it closer to a uniform distribution near the origin
adapt=TRUE
Area.hab=1-Grid_chukchi[["land_cover"]]

for(isim in 1:n.sims){
  cur.file=paste("./sim_data/Lambda_chukchi_",isim,'.Rda',sep='')
  load(cur.file)
  load("./sim_data/Sim_Effort.Rda")
  for(ifl in 1:9){
     #simulate count data conditional on true abundance, flight profile
    Dat=simulate_counts_Chukchi(Lambda=Sim.abund$Lambda.true,Effort=Effort[[ifl]],Thin=Thin.sim)
    #now set up MCMC call
    n.transects=nrow(Effort[[ifl]])
    Prop.photo=rep(0.2,n.transects)
    Mapping=Effort[[ifl]][,"poly"]
    Area.trans=Effort[[ifl]][,"area"]
    DayHour=matrix(1,n.transects,2)
    colnames(DayHour)=c("day","hour")
    
    Out=hierarchical_boss(Dat=Dat,Adj=Adj,Area.hab=Area.hab,Mapping=Mapping,Area.trans=Area.trans,DayHour=DayHour,Thin=Thin,Prop.photo=Prop.photo,Hab.cov=Hab.cov,Obs.cov=Obs.cov,n.obs.cov=n.obs.cov,Hab.pois.formula=Hab.pois.formula,Hab.bern.formula=Hab.bern.formula,Cov.prior.pdf=Cov.prior.pdf,Cov.prior.parms=Cov.prior.parms,Cov.prior.fixed=Cov.prior.fixed,Cov.prior.n=Cov.prior.n,ZIP=ZIP,spat.ind=spat.ind,fix.tau.nu=fix.tau.nu,srr=srr,srr.tol=srr.tol,Inits=Inits,grps=grps,n.species=n.species,Control=Control,adapt=adapt,Prior.pars=Prior.pars,Psi=Psi,post.loss=post.loss)
    save(Out,file=paste('./sim_data/Out_Fl',ifl,'_sim',isim,'.Rda',sep=''))  
    #estimate abundance and save MCMC results
  }
}


#Cur.abund=matrix(apply(Out$Post$N[1,200:2000,],2,'mean'),S,1)
#plot_N_map(1,Cur.abund,Grid=Temp)
#plot_N_map(1,matrix(Sim.abund$Lambda.true[2,],S,1),Grid=Temp)

