#' generate initial values for misID model if not already specified by user
#' @param t.steps number of time steps
#' @param Surveyed  If provided, which entries in "Mapping" to use in extrapolating density (i.e. not including 'extra' zeros)
#' @param DM.hab   a list vector of design matrices for the fixed effects model (one for each species)
#' @param N.hab.par  vector giving number of parameters for the fixed effects model for each species
#' @param G.transect a matrix of the number of groups of animals in area covered by each transect; each row gives a separate species		
#' @param thin.mean Vector giving average thinning proportion for each species
#' @param Area.trans	a vector giving the proportion of a strata covered by each transect
#' @param Area.hab	a vector of the relative areas of each strata
#' @param Mapping	a vector mapping each transect to it's associated strata
#' @param spat.ind  is spatial independence assumed? (TRUE/FALSE)
#' @param grp.mean  a vector giving the pois1 parameter for group size (one entry for each species)
#' @return a list of initial parameter values
#' @export
#' @keywords initial values, mcmc
#' @author Paul B. Conn
generate_inits_BOSSst<-function(t.steps,Surveyed=NULL,DM.hab,N.hab.par,G.transect,thin.mean=thin.mean,Area.trans,Area.hab,Mapping,spat.ind,grp.mean){		
  n.species=nrow(G.transect)
  n.cells=length(Area.hab)
  S=n.cells/t.steps
  if(is.null(Surveyed))G.tot=ceiling(S*apply(G.transect,1,'sum')/(sum(Area.trans)*thin.mean))
  else G.tot=ceiling(S*apply(G.transect[,Surveyed],1,'sum')/(sum(Area.trans[Surveyed])*thin.mean))
  hab=matrix(0,n.species,max(N.hab.par))
  for(isp in 1:n.species){
    hab[isp,1:N.hab.par[isp]]=nlminb(start=rep(0,N.hab.par[isp]),multinom_logL,X=matrix(DM.hab[[isp]][Mapping,],length(Mapping),ncol(DM.hab[[isp]])),C=G.transect[isp,],p=Area.trans*Area.hab[Mapping])$par
      #solve(crossprod(DM.hab[[isp]][Mapping,]),t(DM.hab[[isp]][Mapping,]))%*%log(G.transect[isp,]/Area.trans+.1)
  }
  Par=list(G.tot=G.tot,hab=hab,Eta=matrix(rnorm(n.species*n.cells),n.species,n.cells),
           tau.eta=runif(n.species,0.5,2),tau.eps=runif(n.species,0.5,2))
  Par$G=matrix(0,n.species,n.cells)
  for(isp in 1:n.species){
    Pi=exp(DM.hab[[isp]]%*%hab[isp,1:N.hab.par[isp]])*Area.hab
    for(it in 1:t.steps)Par$G[isp,]=rmultinom(1,G.tot[isp],Pi[((it-1)*S+1):((it-1)*S+S)])
  }
  I.error=(G.transect>Par$G[Mapping])
  if(sum(I.error)>0){
    Which.error=which(G.transect>Par$G[Mapping])
    for(i in 1:sum(I.error))Par$G[Mapping[Which.error[i]]]=Par$G[Mapping[Which.error[i]]]+G.transect[Which.error[i]]
  }
  Par$N=Par$G
  for(isp in 1:n.species)Par$N[isp,]=Par$G[isp,]+rpois(n.cells,grp.mean[isp]*Par$G[isp,])
  if(spat.ind==1)Par$Eta=0*Par$Eta
  Par$Omega=matrix(0,n.species,n.cells)
  Par
}

#function to compute log likelihood for habitat covariate values
multinom_logL<-function(Par,X,C,p){
  Pi=exp(X%*%Par)*p
  -dmultinom(C,size=sum(C),prob=Pi,log=TRUE)
}

sample_nophoto_sp<-function(itrans,Lam,n.sp)sample(n.sp,1,prob=Lam[,itrans])

#function to match Cur.vec with a row of the matrix Pointer
get_place<-function(Cur.vec,Pointer){  #requires that columns of Cur.vec match up with Pointer!!
  n.compare=ncol(Pointer)
  I.match=(Pointer[,1]==Cur.vec[1])
  if(n.compare>1){
    for(itmp in 2:n.compare){
      I.match=I.match*(Pointer[,itmp]==Cur.vec[itmp])
    }  
  }
  which(I.match==1)
}

#functions to extract matrix and list entries given index vectors for rows and column
get_mat_entries<-function(tmp,Mat,Row,Col)Mat[Row[tmp],Col[tmp]]
get_conf_entries<-function(tmp,Conf,List.num,Row,Col)Conf[[List.num[tmp]]][Row[tmp],Col[tmp]]


#' function to convert BOSSst MCMC list vector (used in estimation) into an mcmc object (cf. coda package) 
#' @param MCMC list vector providing MCMC samples for each parameter type 
#' @param N.hab.par see help for mcmc_ds.R
#' @param Cov.par.n see help for mcmc_ds.R
#' @param Hab.names see help for mcmc_ds.R
#' @param Cov.names see help for mcmc_ds.R
#' @param fix.tau.eps see help for mcmc_ds.R
#' @param spat.ind see help for mcmc_ds.R
#' @export
#' @keywords MCMC, coda
#' @author Paul B. Conn
convert.BOSSst.to.mcmc<-function(MCMC,N.hab.par,Cov.par.n,Hab.names,Cov.names,fix.tau.eps=FALSE,spat.ind=TRUE){
  require(coda)
  n.species=nrow(MCMC$Hab)
  n.iter=length(MCMC$Hab[1,,1])
  n.col=n.species*3+sum(N.hab.par)+(1-spat.ind)*n.species+sum(Cov.par.n)*n.species
  n.cells=dim(MCMC$G)[3]
  Mat=matrix(0,n.iter,n.col)
  Mat[,1:n.species]=t(MCMC$N.tot)
  counter=n.species
  col.names=paste("Abund.sp",c(1:n.species),sep='')
  for(isp in 1:n.species){
    Mat[,counter+isp]=MCMC$G.tot[isp,] #total abundance of groups
    col.names=c(col.names,paste("Groups.sp",isp,sep=''))
  }
  counter=counter+n.species
  for(isp in 1:n.species){  #habitat parameters
    Mat[,(counter+1):(counter+N.hab.par[isp])]=MCMC$Hab[isp,,1:N.hab.par[isp]]
    col.names=c(col.names,paste("Hab.sp",isp,Hab.names[[isp]],sep=''))
    counter=counter+sum(N.hab.par[isp])
  }
  if(spat.ind==FALSE){
    Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.eta)
    col.names=c(col.names,paste("tau.eta.sp",c(1:n.species),sep=''))
    counter=counter+n.species
  }
  #if(fix.tau.nu==FALSE){
    Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.eps)
    col.names=c(col.names,paste("tau.eps.sp",c(1:n.species),sep=''))
    counter=counter+n.species
  #}
  if(is.null(Cov.par.n)==FALSE){
    max.par=max(Cov.par.n)
    for(isp in 1:n.species){
      for(ipar in 1:length(Cov.par.n)){
        Mat[,(counter+1):(counter+Cov.par.n[ipar])]=MCMC$Cov.par[isp,,((ipar-1)*max.par+1):((ipar-1)*max.par+Cov.par.n[ipar])]
        counter=counter+Cov.par.n[ipar]
        col.names=c(col.names,paste("Cov.sp",isp,".",Cov.names[[ipar]],sep=''))
      }
    }
  }
  colnames(Mat)=col.names
  Mat=mcmc(Mat)
  Mat
}


#' function to calculate observation probability matrix
#' @param Cur.mat holds observation probability matrix (nrows=number of species, cols=number of observation types)
#' @param Alpha Confusion matrix - rows are for true species - cols are for observed species
#' @param C  Array holding c_{skm} = the probability that an individual truly of species s but observed to be species k will have certainty covariate m
#' @param n.species  number of species
#' @param n.conf number of confidence categories
#' @export
#' @author Paul B. Conn
get_misID_mat<-function(Cur.mat,Alpha,C,n.species,n.conf){
  for(isp in 1:n.species){
    for(iobs in 1:(n.species+1))
      for(iconf in 1:n.conf){
        Cur.mat[isp,n.conf*(iobs-1)+iconf]=Alpha[isp,iobs]*C[isp,iobs,iconf]
     }
  }
  Cur.mat[,ncol(Cur.mat)-n.conf+1]=0
  Cur.mat[,ncol(Cur.mat)-n.conf+1]=1-apply(Cur.mat,1,'sum')
  Cur.mat
}


#' function to calculate posterior predictive loss given the output object from BOSS analysis
#' @param Out Output object from running hierarchicalDS
#' @param burnin Any additional #'s of values from beginning of chain to discard before calculating PPL statistic (default is 0)
#' @return A matrix with posterior variance (P), sums of squares (G) for the posterior mean and median predictions (compared to Observations), and total posterior loss (D)
#' @export
#' @keywords Posterior predictive loss
#' @author Paul B. Conn
post_loss_boss<-function(Out,burnin=0){
  dims.Pred=dim(Out$Pred.det)
  median.Pred=array(0,dim=dims.Pred[2:3])
  mean.Pred=median.Pred
  var.Pred=mean.Pred
  for(itrans in 1:dims.Pred[2]){
    for(iobs in 1:dims.Pred[3]){
        median.Pred[itrans,iobs]=median(Out$Pred.det[(burnin+1):dims.Pred[1],itrans,iobs])
        mean.Pred[itrans,iobs]=mean(Out$Pred.det[(burnin+1):dims.Pred[1],itrans,iobs])
        var.Pred[itrans,iobs]=var(Out$Pred.det[(burnin+1):dims.Pred[1],itrans,iobs])
    }  
  }
  sum.sq.mean=sum((Out$Obs.det-mean.Pred)^2)
  sum.sq.median=sum((Out$Obs.det-median.Pred)^2)
  Loss=matrix(0,2,3)
  colnames(Loss)=c("P","G","D")
  rownames(Loss)=c("mean","median")
  Loss[,1]=sum(var.Pred)
  Loss[1,2]=sum.sq.mean
  Loss[2,2]=sum.sq.median
  Loss[,3]=rowSums(Loss[1:2,1:2])
  Loss
}

plot_covar<-function(DM,MCMC,Vars,n.species,n.points=20,Sp.names,const.tau=NULL,bern=FALSE){
  #assumes additive covariates (no interactions) and that polynomial terms have names like "ice2"
  #determine number of plots that will be required
  #determine which covariates have the keyword "var" in them
  require(ggplot2)
  my_grep<-function(cur.eff,cur.var)(length(grep(cur.var,cur.eff))>0)
  n.var=length(Vars)
  mcmc.str="Hab.pois.sp"
  if(bern==TRUE)mcmc.str="Hab.bern.sp"
  plot.df=data.frame(matrix(0,n.species*n.points*n.var,4))
  colnames(plot.df)=c("Species","Cov.name","Value","Abundance")
  if(is.null(const.tau)){
    Tau=rep(0,n.species)
    for(isp in 1:n.species)Tau[i]=mean(eval(parse(text=paste("MCMC[,'tau.nu.sp",isp,"']",sep=''))))
  }
  else Tau=rep(const.tau,n.species)
  cur.pl=1
  for(isp in 1:n.species){
    #pull out relevant columns from MCMC object
    I.col=sapply(colnames(MCMC),'my_grep',cur.var=paste(mcmc.str,isp,sep=''))
    Cur.MCMC=MCMC[,which(I.col==TRUE)]
    for(ivar in 1:n.var){
      cur.nchar=nchar(Vars[ivar])
      cur.min=min(DM[[isp]][,Vars[ivar]])
      cur.max=max(DM[[isp]][,Vars[ivar]])
      Cur.x=seq(cur.min,cur.max,length.out=n.points)
      Cur.DM=matrix(0,n.points,ncol(DM[[isp]]))
      Cur.DM[,which(colnames(DM[[isp]])==Vars[ivar])]=Cur.x
      I.effect=sapply(colnames(DM[[isp]]),"my_grep",cur.var=Vars[ivar])  
      if(sum(I.effect)>1){ # fill in polynomial effects if >1
        Cur.which=which(I.effect==TRUE)
        for(i in 1:sum(I.effect)){
          if(colnames(DM[[isp]])[Cur.which[i]]!=Vars[ivar]){
            pol.eff=as.numeric(substr(colnames(DM[[isp]])[Cur.which[i]],cur.nchar+1,nchar(colnames(DM[[isp]])[Cur.which[i]])))
            Cur.DM[,Cur.which[i]]=Cur.x^pol.eff
          }
        }
      }
      Cur.which=which(I.effect==FALSE)
      for(i in 1:length(Cur.which)){
        Cur.DM[,Cur.which[i]]=mean(DM[[isp]][,Cur.which[i]])
      }
      Cur.resp=exp(Cur.DM%*%t(Cur.MCMC)+0.5/Tau[isp])
      plot.df[cur.pl:(cur.pl+n.points-1),"Species"]=Sp.names[isp]
      plot.df[cur.pl:(cur.pl+n.points-1),"Cov.name"]=Vars[ivar]
      plot.df[cur.pl:(cur.pl+n.points-1),"Value"]=Cur.x
      plot.df[cur.pl:(cur.pl+n.points-1),"Abundance"]=apply(Cur.resp,1,'mean') #mean posterior prediction
      cur.pl=cur.pl+n.points
    }
  }
  plot.df[,1]=as.factor(plot.df[,1])
  plot.df[,2]=as.factor(plot.df[,2])
  myPlot=ggplot(plot.df)+geom_line(aes(x=Value,y=Abundance,color=Species),size=1.5)+facet_wrap(~Cov.name,scales='free')
  myTheme=theme(text=element_text(size=20),axis.text=element_text(size=11))
  myPlot=myPlot+myTheme
  myPlot
}

#' Produce an RW1 adjacency matrix for a rectangular grid for use with areal spatial models (queens move)
#' @param x number of cells on horizontal side of grid
#' @param y number of cells on vertical side of grid
#' @param byrow If TRUE, cell indices are filled along rows (default is FALSE)
#' @return adjacency matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn \email{paul.conn@@noaa.gov}
rect_adj <- function(x,y,byrow=FALSE){
  Ind=matrix(c(1:(x*y)),y,x,byrow=byrow)
  if(byrow==TRUE)Ind=t(Ind)
  n.row=nrow(Ind)
  n.col=ncol(Ind)
  Adj=matrix(0,x*y,x*y)
  for(i in 1:n.row){
    for(j in 1:n.col){
      if(i==1 & j==1){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
      }
      if(i==1 & j>1 & j<n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
      }
      if(i==1 & j==n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
      }
      if(i>1 & i<n.row & j==1){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
      }
      if(i>1 & i<n.row & j>1 & j<n.col){
        cur.nums=c(Ind[i,j]-n.row-1,Ind[i,j]-n.row,Ind[i,j]-n.row+1,Ind[i,j]-1,Ind[i,j]+1,Ind[i,j]+n.row-1,Ind[i,j]+n.row,Ind[i,j]+n.row+1)
        Adj[Ind[i,j],cur.nums]=1
      }
      if(i>1 & i<n.row & j==n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
        
      }
      if(i==n.row & j==1){
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
      }
      if(i==n.row & j>1 & j<n.col){
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
      }
      if(i==n.row & j==n.col){
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
      }
    }
  }
  if(byrow==TRUE)Adj=t(Adj)
  return(Adj)
}

#' Produce weights for bivariate normal kernel for upper left triangular quadrant
#' To save computing time, some building blocks are passed into
#' the function
#' @param Tmp.vec A length n vector for holding integrals of bivariate normal
#' @param XY 2xn matrix giving x and y distances
#' @param Sigma vector of length 2 giving sd of bivariate normal in the x and y directions 
#' @return a filled vector that holds unnormalized redistribution kernel values
#' @export 
#' @keywords bivariate normal, kernel weight
#' @author Paul Conn \email{paul.conn@@noaa.gov}
d_biv_normal<-function(Tmp.vec,XY,Sigma){
  return(dnorm(XY[1,],0,Sigma[1])*dnorm(XY[2,],0,Sigma[2]))
}
  

plot_N_map<-function(cur.t,N,Grid,highlight=NULL,leg.title="Abundance"){
  require(rgeos)
  require(ggplot2)
  library(RColorBrewer)
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  Tmp=Grid[[1]]
  if(is.null(highlight)==FALSE){
    midpoints=data.frame(gCentroid(Tmp[highlight,],byid=TRUE))
    colnames(midpoints)=c("Easting","Northing")
  }
  Abundance=N[,cur.t]
  Cur.df=cbind(data.frame(gCentroid(Tmp,byid=TRUE)),Abundance)
  new.colnames=colnames(Cur.df)
  new.colnames[1:2]=c("Easting","Northing")
  colnames(Cur.df)=new.colnames
  tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  p1=ggplot(Cur.df)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme+scale_fill_gradientn(colours=myPalette(100),name=leg.title)
  if(is.null(highlight)==FALSE){
    #p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067,xmax=Easting,ymin=Northing,ymax=Northing+25067))
    p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
  }
  p1
}
  
#' Calculate partial derivative vector for for Langevin omega MH updates (for CPIF model)
#' (Assuming function is called separately for each multinomial update)
#' @param Omega Current values of Omega
#' @param Counts Animal counts for current time step being updated
#' @param Mu Expected values for Omega
#' @param tau precision for exchangable errors
#' @param G.sum total group abundance for sampled cells at the current time step
#' @param Cur.thin A vector of values giving availability * proportion of cell surveyed for cells being evaluated
#' @return Vector of partial derivates evaluated at Omega
#' @export 
#' @keywords logit, expit
#' @author Paul Conn \email{paul.conn@@noaa.gov}
d_logP_omega<-function(Omega,Counts,Mu,tau,G.sum,Cur.thin){
  Omega.exp=exp(Omega)
  cur.sum=sum(Omega.exp)
  #cur.sum2=cur.sum-exp(Omega)
  #return(Counts - (Omega-Mu)*tau - exp(Omega)/cur.sum + (G.sum-sum(Counts))*Cur.thin*cur.sum2*exp(Omega)/(cur.sum*((Cur.thin-1)*exp(Omega)-cur.sum2)))
  return(Counts+(Mu-Omega)*tau-Counts*Omega.exp/cur.sum)
} 

#' function to sample species using apply
#' @param probability weights for selecting each species
#' @param n.species number of species to sample from
#' @return A sampled species value
#' @export
#' @keywords species classification
#' @author Paul B. Conn
sample_species<-function(probs,n.species)sample(c(1:n.species),1,prob=probs)

#' function to concatenate posterior predictions from multiple files
#' @param fname Char string providing base file name (no extension)
#' @param n.files Number of files (naming convention is fname#.RData)
#' @return Merged list object with posterior predictions
#' @export
#' @author Paul B. Conn
cat_preds<-function(fname,n.files){
  load(paste0(fname,"1",".RData"))
  Dim=dim(Post$N)
  Tmp=array(0,dim=c(Dim[1],Dim[2]*n.files,Dim[3]))
  Tmp[,1:Dim[2],]=Post$N
  if(n.files>1){
    for(ifile in 2:n.files){
      load(paste0(fname,ifile,".RData"))
      Tmp[,(Dim[2]*(ifile-1)+1):(Dim[2]*ifile),]=Post$N
    }
  }
  Tmp
}



