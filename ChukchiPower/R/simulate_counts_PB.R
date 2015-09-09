#' given underlying density, effort, and observation process information, simulate animal counts
#' @param Lambda A (3 x n.cells) matrix holding expected abundance for bearded seals, ringed seals, and polar bear in each grid cell
#' @param Effort A matrix or data.frame giving effort in sampled cells.  Columns are "poly", giving which cells are sampled, and "area", the proportion of each sampled grid cell that is surveyed
#' @param Thin A vector giving an additional thinning factor for each species (corresponding to detection and availability)
#' @return A data frame holding animal observation data.  Columns include 
#'     -"Transect" (which cell the observation occurred in)
#'     -"Photo" (whether a photo was obtained - here these are all set to 1.0)
#'     -"Obs" (what species was observed), and
#'     -"Group" (group size for each hotspot) - here assumed to by 1.0
#' @export
#' @keywords simulation, count data
#' @author Paul B. Conn
simulate_counts_PB<-function(Lambda,Effort,Thin){
  n.sampled=nrow(Effort)
  Covered=matrix(0,3,n.sampled)
  Avail=Covered
  Covered[1,]=rpois(n.sampled,Lambda[1,Effort[,"poly"]]*Effort[,"area"])
  Covered[2,]=rpois(n.sampled,Lambda[2,Effort[,"poly"]]*Effort[,"area"]) 
  Covered[3,]=rpois(n.sampled,Lambda[3,Effort[,"poly"]]*Effort[,"area"])   
  Avail[1,]=rbinom(n.sampled,Covered[1,],Thin[1])
  Avail[2,]=rbinom(n.sampled,Covered[2,],Thin[2])
  Avail[3,]=rbinom(n.sampled,Covered[3,],Thin[3])  
  lam.anomaly=mean(Lambda)*0.2
  Anomaly=rpois(n.sampled,lam.anomaly*Effort[,"area"])
  Dat=data.frame(matrix(0,10000,4))
  colnames(Dat)=c("Transect","Photo","Obs","Group")
  Dat[,"Obs"]=NA
  Dat[,"Group"]=NA
  cur.pl=1
  for(icell in 1:n.sampled){
    for(isp in 1:3){
      if(Avail[isp,icell]>0){
        #determine how many photographed
        n.photo=rbinom(1,Avail[isp,icell],0.84)
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
        Dat[cur.pl:(cur.pl+n.photo-1),3]=4
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