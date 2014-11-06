#given underlying density, effort, and observation process information, simulate animal counts

simulate_counts_Chukchi<-function(Lambda,Effort,Thin){
  n.sampled=nrow(Effort)
  Covered=matrix(0,2,n.sampled)
  Avail=Covered
  Covered[1,]=rpois(n.sampled,Lambda[1,Effort[,"poly"]]*Effort[,"area"])
  Covered[2,]=rpois(n.sampled,Lambda[2,Effort[,"poly"]]*Effort[,"area"]) 
  Avail[1,]=rbinom(n.sampled,Covered[1,],Thin[1])
  Avail[2,]=rbinom(n.sampled,Covered[2,],Thin[2])
  lam.anomaly=mean(Lambda)*0.2
  Anomaly=rpois(n.sampled,lam.anomaly*Effort[,"area"])
  Dat=data.frame(matrix(0,10000,4))
  colnames(Dat)=c("Transect","Photo","Obs","Group")
  Dat[,"Obs"]=NA
  Dat[,"Group"]=NA
  cur.pl=1
  for(icell in 1:n.sampled){
    for(isp in 1:2){
      if(Avail[isp,icell]>0){
        #determine how many photographed
        n.photo=rbinom(1,Avail[isp,icell],0.8)
        if(n.photo>0){
          Dat[cur.pl:(cur.pl+n.photo-1),1]=icell
          Dat[cur.pl:(cur.pl+n.photo-1),2]=1
          Dat[cur.pl:(cur.pl+n.photo-1),3]=isp
          Dat[cur.pl:(cur.pl+n.photo-1),4]=1   
          cur.pl=cur.pl+n.photo
        }
        n.no.photo=Avail[isp,icell]-n.photo
        if(n.no.photo>0){
          Dat[cur.pl:(cur.pl+n.no.photo-1),1]=icell
          cur.pl=cur.pl+n.no.photo #advance for # not photographed
        }
      }
    }
    #now for anomalies
    if(Anomaly[icell]>0){
      n.photo=rbinom(1,Anomaly[icell],0.8)
      if(n.photo>0){
        Dat[cur.pl:(cur.pl+n.photo-1),1]=icell
        Dat[cur.pl:(cur.pl+n.photo-1),2]=1
        Dat[cur.pl:(cur.pl+n.photo-1),3]=3
        Dat[cur.pl:(cur.pl+n.photo-1),4]=1   
        cur.pl=cur.pl+n.photo
      }
      n.no.photo=Anomaly[icell]-n.photo
      if(n.no.photo>0){
        Dat[cur.pl:(cur.pl+n.no.photo-1),1]=icell
        cur.pl=cur.pl+n.no.photo #advance for # not photographed
      }
    }
  }
  Dat=Dat[1:(cur.pl-1),]
}