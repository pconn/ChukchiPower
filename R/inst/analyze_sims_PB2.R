#analyze power simulation results
library(xtable)
set.seed(23456)
load("./sim_data/Sim_Effort.Rda")

Models=c(1,2,3,4,5,6,7,9,11)
n.sims=100
n.flights=9
n.models=length(Models)

Bias=vector("list",3)
Cov=MSE=CV=Lam.true=Bias
for(isp in 1:3){
  if(isp<3){
    Bias[[isp]]=data.frame(matrix(NA,n.models,n.flights))
    colnames(Bias[[isp]])=c("A1","A2","B1","B2","B3","C1","C2","C3","C4")
    Cov[[isp]]=MSE[[isp]]=CV[[isp]]=Bias[[isp]]                       
    Lam.true[[isp]]=rep(0,n.sims)
  }
  else{
    Bias[[isp]]=data.frame(matrix(NA,6,n.flights))
    colnames(Bias[[isp]])=c("A1","A2","B1","B2","B3","C1","C2","C3","C4")
    Cov[[isp]]=MSE[[isp]]=CV[[isp]]=Bias[[isp]]                       
    Lam.true[[isp]]=rep(0,n.sims)    
  }
}

#load('./sim_data/Out_Fl3_sim1.Rda')
#cur.pl=1
for(isim in 1:n.sims){
  load(paste("./sim_data/Lambda_PB_",isim,".Rda",sep=''))
  for(isp in 1:3){
    lam.true=sum(Sim.abund$Lambda.true[isp,])
    Lam.true[[isp]][isim]=lam.true
  }
}

estimator="Mean"
for(isp in 1:3){
  for(imodel in 1:n.models){   
    imod=Models[imodel]
    for(ifl in 1:9){
      if((isp<3 & (ifl>1 | imod%in%c(1,2,7)))| (isp==3 & imod%in%c(1,2,3,5) & ifl>2)){
        Tmp=read.csv(paste("./sim_data/OutPB_mod_",imod,"_isp_",isp,"Fl",ifl,".csv",sep=''),header=TRUE)
        Bias[[isp]][imodel,ifl]=median((Tmp[,estimator]-Lam.true[[isp]])/Lam.true[[isp]])
        MSE[[isp]][imodel,ifl]=sqrt(mean((Tmp[,estimator]-Lam.true[[isp]])^2))
        CV[[isp]][imodel,ifl]=median(Tmp[,"SE"]/Tmp[,estimator])
        Cov[[isp]][imodel,ifl]=mean((Lam.true[[isp]]>Tmp[,"Quant05"] & Lam.true[[isp]]<Tmp[,"Quant95"]))
      }
    }
  }
}
Bias[[3]]=Bias[[3]][c(1,2,3,5),]
MSE[[3]]=MSE[[3]][c(1,2,3,5),]
CV[[3]]=CV[[3]][c(1,2,3,5),]
Cov[[3]]=Cov[[3]][c(1,2,3,5),]


#Thought I could do models as rownames but you can't have duplicated row names so just adding model w/ cbind
Model=c("covs","covs + samp.dens","stratum","covs + stratum","stratum + samp.dens","covs + stratum + samp.dens","covs + RE","stratum + RE","stratum + samp.dens + RE")
for(isp in 1:2){
  Bias[[isp]]=cbind(Model,Bias[[isp]])
  MSE[[isp]]=cbind(Model,MSE[[isp]])
  CV[[isp]]=cbind(Model,CV[[isp]])
  Cov[[isp]]=cbind(Model,Cov[[isp]])
}
Model=c("seals","seals + samp.dens","stratum","stratum + samp.dens")
isp=3
Bias[[isp]]=cbind(Model,Bias[[isp]])
MSE[[isp]]=cbind(Model,MSE[[isp]])
CV[[isp]]=cbind(Model,CV[[isp]])
Cov[[isp]]=cbind(Model,Cov[[isp]])



Bias.table=xtable(rbind(Bias[[1]],Bias[[2]],Bias[[3]]))
MSE.table=xtable(rbind(MSE[[1]],MSE[[2]],MSE[[3]]),digits=0)
CV.table=xtable(rbind(CV[[1]],CV[[2]],CV[[3]]))
Cov.table=xtable(rbind(Cov[[1]],Cov[[2]],Cov[[3]]))

print(Bias.table, include.rownames=FALSE)
print(CV.table, include.rownames=FALSE)
print(MSE.table, include.rownames=FALSE)
print(Cov.table, include.rownames=FALSE)

