### plot preferential sampling plots
source('../SpatPred/SpatPred/R/util_funcs.R')
library(sp)
library(maptools)
library(rgeos)

isim=3
ifl=7
S=505
load(paste("./sim_data/Counts_fl_",ifl,"_sim",isim,'.Rda',sep=''))
load(paste("./sim_data/Lambda_PB_",isim,'.Rda',sep=''))
load(paste("./Grid_chukchi.Rda"))
load("./sim_data/Sim_Effort.Rda")

load(paste("./sim_data/MCMC_mod9_fl",ifl,"_sp1.Rda",sep='')) 
Bearded.MCMC=Out$MCMC
load(paste("./sim_data/MCMC_mod9_fl",ifl,"_sp2.Rda",sep='')) 
Ringed.MCMC=Out$MCMC
load(paste("./sim_data/MCMC_mod9_fl",ifl,"_sp3.Rda",sep='')) 
PB.MCMC=Out$MCMC

Mapping=Effort[[ifl]]$poly

center=apply(bbox(Grid.chukchi$Grid), 1, mean) #gCentroid(Grid.chukchi$Grid)
Grid.rotate=elide(Grid.chukchi$Grid,rotate=-90,center=center)


Lam.1=matrix(Sim.abund$Lambda.true[1,],S,1)
Lam.2=matrix(Sim.abund$Lambda.true[2,],S,1)
Lam.3=matrix(Sim.abund$Lambda.true[3,],S,1)

N.1=matrix(apply(Bearded.MCMC$Pred[,401:2400],1,'median'),S,1)
N.2=matrix(apply(Ringed.MCMC$Pred[,401:2400],1,'median'),S,1)
N.3=matrix(apply(PB.MCMC$Pred[,250:2249],1,'median'),S,1)



library(ggplot2)
library(RColorBrewer)
library(gridExtra)
greenPalette <- colorRampPalette(brewer.pal(9, "Greens"))
YlOrBr<- colorRampPalette(brewer.pal(9, "YlOrBr"))
txt.size=8

Grid.list=vector("list",1)
Grid.list[[1]]=Grid.rotate
#max.N=max(c(N.sys,N.propX,N.propN))
max.N=400
Lam.plot.1=plot_N_map(1,Lam.1,Grid=Grid.list)+ggtitle("N, bearded")+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   
N.plot.1=plot_N_map(1,N.1,Grid=Grid.list)+ggtitle(expression(paste(hat(N),", bearded")))+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   
max.N=10000
Lam.plot.2=plot_N_map(1,Lam.2,Grid=Grid.list)+ggtitle("N, ringed")+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   
N.plot.2=plot_N_map(1,N.2,Grid=Grid.list)+ggtitle(expression(paste(hat(N),", ringed")))+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   
max.N=12
Lam.plot.3=plot_N_map(1,Lam.3,Grid=Grid.list)+ggtitle("N, polar bear")+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   
N.plot.3=plot_N_map(1,N.3,Grid=Grid.list)+ggtitle(expression(paste(hat(N),", polar bear")))+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=YlOrBr(100),limits=c(0,max.N))   

Counts=matrix(0,3,length(Mapping))
for(i in 1:3){
  Cur.dat=Dat[which(Dat[,"Obs"]==i),]
  Count=tabulate(Cur.dat[,"Transect"])
  #add in other zero counts not included w/ tabulate call
  if(length(Mapping)>length(Count))Count=c(Count,rep(0,length(Mapping)-length(Count)))#add in other zero counts
  Counts[i,]=Count
}

Count=matrix(rep(NA,S),ncol=1)
Count[Mapping]=Counts[1,]
Count.plot.1=plot_N_map(1,matrix(Count,ncol=1),Grid=Grid.list)+ggtitle("Counts, bearded")+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=greenPalette(100),limits=c(0,10))   
Count[Mapping]=Counts[2,]
Count.plot.2=plot_N_map(1,matrix(Count,ncol=1),Grid=Grid.list)+ggtitle("Counts, ringed")+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=greenPalette(100),limits=c(0,500))   
Count[Mapping]=Counts[3,]
Count.plot.3=plot_N_map(1,matrix(Count,ncol=1),Grid=Grid.list)+ggtitle("Counts, polar bear")+theme(plot.title = element_text(hjust = 0),legend.title = element_blank(),text=element_text(size=txt.size),plot.margin=unit(c(0,0,0,0),"lines"))+scale_fill_gradientn(colours=greenPalette(100),limits=c(0,2))   

  

pdf(file="PB_sim_maps.pdf")
grid.arrange(arrangeGrob(Lam.plot.1,Lam.plot.2,Lam.plot.3,Count.plot.1,Count.plot.2,Count.plot.3,N.plot.1,N.plot.2,N.plot.3,widths=unit(0.33,"npc"),heights=unit(0.33,"npc"),nrow=3))
dev.off()




    