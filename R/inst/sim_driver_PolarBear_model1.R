#master simulation driver for Chukchi power analysis
library(Matrix)
library(GenKern)  #to model sampling density
source('../hierarchicalDS/hierarchicalDS/R/spat_funcs.R')
source('./R/R/simulate_density_PB.R')
source('./R/R/simulate_counts_PB.R')
source('../BOSSst/BOSSst/R/util_funcs.R')
source('./R/R/Chukchi_mcmc.R')
source('../BOSS/BOSS/R/util_funcs.R')

load("Grid_chukchi.Rda")
load("Grid_chukchi_ExpDensities_PB.Rdata")
#load("Knot_cell_distances.Rdata")  
load("Haulout_samples.Rdat")
load("Chukchi_centroids.Rda")


Adj=Matrix(Grid.chukchi$Adj)

n.sims=1
n.species=3  #the 3rd is for anomalies/other
n.models=12
SIM=FALSE  #if true, simulate each set of true abundance data; otherwise load from memory
set.seed(20935)
S=length(Grid_chukchi)

Cur.abund=matrix(Grid_chukchi[["D.pbear"]]*625,S,1)
Temp=vector("list",1)
Temp[[1]]=Grid_chukchi
#plot_N_map(1,Cur.abund,Grid=Temp)


if(SIM){
  for(isim in 1:n.sims){
    Sim.abund=sim_density_PB()
    cur.file=paste("./sim_data/Lambda_PB_",isim,'.Rda',sep='')
    save(Sim.abund,file=cur.file)
  }
}

#take a look of some plots of 'noisy' and expected seal densities
Temp=vector("list",1)
Temp[[1]]=Grid_chukchi
#plot_N_map(1,matrix(Sim.abund$Lambda.true[3,],S,1),Grid=Temp)

Temp=vector("list",1)
Temp[[1]]=Grid_chukchi
#plot_N_map(1,matrix(Grid_chukchi[["D.ringed"]],S,1),Grid=Temp)
#plot_N_map(1,matrix(Hab.cov$Seals,S,1),Grid=Temp,leg.title="Rel seal biomass")

##generate posterior thinning samples
Thin=array(0,dim=c(4,1,1000))
P=rbeta(1000,67,5)
Thin[1,1,]=Haulout.samples[[1]][29,23,]
Thin[2,1,]=Thin[1,1,]*(0.65/mean(Thin[1,1,])) #rescale ringed to have a mean of 0.65
Thin[3,1,]=1
Thin=Thin*P
Thin[4,1,]=1
Thin.sim=apply(Thin,1,'mean')
rm(Haulout.samples)


#simulation independent MCMC options
Strata=as.factor(Grid_chukchi@data[,"D.bearded"])
Hab.cov=Grid_chukchi
Hab.cov@data$Strata=Strata
Hab.cov@data$Seals=Grid_chukchi[["D.ringed"]]+5*Grid_chukchi[["D.bearded"]]
Hab.cov@data$Seals=Hab.cov$Seals/mean(Hab.cov$Seals)
Hab.cov@data$sqrt.mainland=sqrt(Grid_chukchi[["dist_mainland"]])
Hab.cov@data$EastNorth=Hab.cov$Easting*Hab.cov$Northing
Formula=vector("list",n.models)
Spat.mod=matrix(0,n.models,3)
Spat.mod[(n.models/2+1):n.models,]=2
for(imod in 1:n.models)Formula[[imod]]=vector("list",3)
for(ispace in 1:2){
  for(isp in 1:2){
    Formula[[(ispace-1)*n.models/2+1]][[isp]]=~dist_mainland+sqrt.mainland+Easting+Northing+EastNorth
    Formula[[(ispace-1)*n.models/2+2]][[isp]]=~dist_mainland+sqrt.mainland+Easting+Northing+EastNorth+Samp.dens
    Formula[[(ispace-1)*n.models/2+3]][[isp]]=~Strata
    Formula[[(ispace-1)*n.models/2+4]][[isp]]=~dist_mainland+sqrt.mainland+Easting+Northing+EastNorth+Strata
    Formula[[(ispace-1)*n.models/2+5]][[isp]]=~Strata+Samp.dens
    Formula[[(ispace-1)*n.models/2+6]][[isp]]=~Strata+dist_mainland+sqrt.mainland+Easting+Northing+EastNorth+Samp.dens
  }
  Formula[[(ispace-1)*n.models/2+1]][[3]]=~Seals
  Formula[[(ispace-1)*n.models/2+2]][[3]]=~Seals+Samp.dens
  Formula[[(ispace-1)*n.models/2+3]][[3]]=~Strata
  Formula[[(ispace-1)*n.models/2+4]][[3]]=~Seals+Strata
  Formula[[(ispace-1)*n.models/2+5]][[3]]=~Strata+Samp.dens
  Formula[[(ispace-1)*n.models/2+6]][[3]]=~Strata+Seals+Samp.dens 
}  
Models.with.strata=c(3,4,5,6,9,10,11,12)
#dist_mainland+sqrt.mainland+Easting+Northing+EastNorth
Inits=list(tau.nu=rep(100,n.species))
Prior.pars=list(a.eta=1,b.eta=.01,a.eps=1,b.eps=.01,beta.tau=0.01) #(1,.01) prior makes it closer to a uniform distribution near the origin
adapt=TRUE
Area.hab=1-Grid_chukchi[["land_cover"]]

#simulate count data conditional on true abundance, flight profile

load("./sim_data/Sim_Effort.Rda")
if(SIM){
  for(isim in 1:n.sims){
    cur.file=paste("./sim_data/Lambda_PB_",isim,'.Rda',sep='')
    load(cur.file)
    for(ifl in 1:9){
      Dat=simulate_counts_PB(Lambda=Sim.abund$Lambda.true,Effort=Effort[[ifl]],Thin=Thin.sim)
      save(Dat,file=paste("./sim_data/Counts_fl_",ifl,"_sim",isim,".rda",sep=''))
    }
  }
}

for(imod in 8:8){
  min.fl=1
  if(imod %in% Models.with.strata)min.fl=2
  for(ifl in 2:2){
    Out.table=vector("list",3)
    for(isp in 1:3){
      Out.table[[isp]]=matrix(0,n.sims,6)
      colnames(Out.table[[isp]])=c("Mean","Median","Quant025","Quant05","Quant95","Quant975")
    }
    for(isim in 1:n.sims){      
      cat(paste("Simulation #",isim,", Flight",ifl,"\n"))
      load(paste("./sim_data/Counts_fl_",ifl,"_sim",isim,".rda",sep=''))
      
      #now set up MCMC call
      n.transects=nrow(Effort[[ifl]])
      Mapping=Effort[[ifl]][,"poly"]
      Area.trans=Effort[[ifl]][,"area"]*0.84  # reduce area to that where photographs are available
      DayHour=matrix(1,n.transects,2)
      colnames(DayHour)=c("day","hour")
      
      XY=Chukchi_centroids[Mapping,]
      X=XY$x-min(XY$x)
      Y=XY$y-min(XY$y)
      X=X/max(X)
      Y=Y/max(Y)
      KDE=KernSur(x=X,y=Y) #,xbandwidth=0.03,ybandwidth=0.03)
      #extract values for each X,Y
      Samp.dens=rep(0,S)
      XY=Chukchi_centroids
      X=XY$x-min(XY$x)
      Y=XY$y-min(XY$y)
      X=X/max(X)
      Y=Y/max(Y)
      for(i in 1:S)Samp.dens[i]=KDE$zden[nearest(KDE$xords,X[i]),nearest(KDE$yords,Y[i])]
      Hab.cov@data$Samp.dens=Samp.dens
      Hab.cov@data$Log.samp.dens=log(Samp.dens)
      Hab.cov@data$Samp.dens=Hab.cov$Samp.dens/mean(Hab.cov$Samp.dens)
      Hab.cov@data$Log.samp.dens=Hab.cov$Log.samp.dens/mean(Hab.cov$Log.samp.dens)
      
      for(isp in 1:3){  
        Dat.sp=Dat[which(Dat[,"Obs"]==isp),]
        Count=tabulate(Dat.sp[,"Transect"])
        #add in other zero counts not included w/ tabulate call
        if(length(Mapping)>length(Count))Count=c(Count,rep(0,length(Mapping)-length(Count)))#add in other zero counts
        Effort.mcmc=data.frame(Counts=Count,Mapping=Mapping)
        Thin.mc=array(Thin[isp,,],dim=c(1,1,dim(Thin)[3]))
        Control=list(iter=3100,burnin=100,thin=100,MH.nu=rep(.6,S),adapt=TRUE,srr.tol=0.5,fix.tau.epsilon=FALSE)
        Out<-Chukchi_mcmc(formula=Formula[[imod]][[isp]],Data=Hab.cov,Effort=Effort.mcmc,spat.mod=Spat.mod[imod,isp],Offset=Area.trans,Area.adjust=Area.hab,Control=Control,Assoc=Adj,K=NULL,DayHour=DayHour,Thin=Thin.mc,Names.gam=NULL,Knots.gam=NULL,Prior.pars=Prior.pars,Precision.pars=NULL)
        Control=Out$Control
        Control$adapt=FALSE
        if(isp<3){
          Control$iter=60100
          Control$thin=25
          start.iter=400
          end.iter=2399          
        }
        if(isp==3){
          Control$iter=450000
          start.iter=250
          Control$thin=200 
          end.iter=2249
        }
        Out<-Chukchi_mcmc(formula=Formula[[imod]][[isp]],Data=Hab.cov,Effort=Effort.mcmc,spat.mod=Spat.mod[imod,isp],Offset=Area.trans,Area.adjust=Area.hab,Control=Control,Assoc=Adj,K=NULL,DayHour=DayHour,Thin=Thin.mc,Names.gam=NULL,Knots.gam=NULL,Prior.pars=Prior.pars,Precision.pars=NULL)
        save(Out,file=paste('MCMC_mod9_fl',ifl,"_sp",isp,".Rda",sep=''))
        N.hat=apply(Out$MCMC$Pred[,start.iter:end.iter],2,'sum')
        Out.table[[isp]][isim,]=c(mean(N.hat),quantile(N.hat,c(0.5,0.025,0.05,0.95,0.975)))        
      }
    }
    for(isp in 1:3)write.csv(Out.table[[isp]],file=paste('./sim_data/OutPB_mod_',imod,'_isp_',isp,'Fl',ifl,'.csv',sep='')) 
  }
}

#N.hat=apply(Out$MCMC$Pred,2,'sum')
#plot(N.hat)
#Cur.abund=matrix(apply(Out$MCMC$Pred,1,'median'),S,1)
#plot_N_map(1,Cur.abund,Grid=Temp)
#plot_N_map(1,matrix(Sim.abund$Lambda.true[2,],S,1),Grid=Temp)
#plot_N_map(1,matrix(Hab.cov$Samp.dens,S,1),Grid=Temp)

#N.tot=apply(Out$Post$N[3,200:2000,],1,'sum')
#hist(N.tot)

