#' generate initial values for misID model if not already specified by user
#' @param DM.hab.pois   a list of design matrices for the Poisson habitat model (elements are named sp1,sp2, etc.)
#' @param DM.hab.bern   If a hurdle model, a list of design matrices for the Bernoulli habitat model (elements are named sp1,sp2, etc.) (NULL if not hurdle)
#' @param N.hab.pois.par  vector giving number of parameters in the Poisson habitat model for each species
#' @param N.hab.bern.par  vector giving number of parameters in the Bernoulli habitat model for each species (NULL if not hurdle)
#' @param G.transect a matrix of the number of groups of animals in area covered by each transect; each row gives a separate species		
#' @param Area.trans	a vector giving the proportion of a strata covered by each transect
#' @param Area.hab	a vector of the relative areas of each strata
#' @param Mapping	a vector mapping each transect to it's associated strata
#' @param spat.ind  is spatial independence assumed? (TRUE/FALSE)
#' @param grp.mean  a vector giving the pois1 parameter for group size (one entry for each species)
#' @return a list of initial parameter values
#' @export
#' @keywords initial values, mcmc
#' @author Paul B. Conn
generate_inits_BOSS<-function(DM.hab.pois,DM.hab.bern,N.hab.pois.par,N.hab.bern.par,G.transect,Area.trans,Area.hab,Mapping,spat.ind,grp.mean){		
  i.hurdle=1-is.null(DM.hab.bern)
  n.species=nrow(G.transect)
  n.cells=length(Area.hab)
  hab.pois=matrix(0,n.species,max(N.hab.pois.par))
  hab.bern=NULL
  tau.eta.bern=NULL
  Eta.bern=NULL
  if(i.hurdle==1){
    hab.bern=matrix(0,n.species,max(N.hab.bern.par))
    tau.eta.bern=runif(n.species,0.5,2)
    Eta.bern=matrix(rnorm(n.species*n.cells),n.species,n.cells)
  }
  Nu=matrix(0,n.species,n.cells)
  for(isp in 1:n.species){
    Nu[isp,]=log(max(G.transect[isp,])/mean(Area.trans)*exp(rnorm(length(Area.hab),0,0.1)))
  }
  Par=list(hab.pois=hab.pois,hab.bern=hab.bern,
           Nu=Nu,Eta.pois=matrix(rnorm(n.species*n.cells),n.species,n.cells),Eta.bern=Eta.bern,
           tau.eta.pois=runif(n.species,0.5,2),tau.eta.bern=tau.eta.bern,tau.nu=runif(n.species,0.5,2))
  Par$hab.pois[,1]=log(apply(G.transect,1,'mean')/(mean(Area.trans)*mean(Area.hab))*exp(rnorm(n.species,0,1)))
  Par$G=round(exp(Par$Nu)*Area.hab*exp(rnorm(length(Par$Nu))))
  I.error=(G.transect>Par$G[Mapping])
  if(sum(I.error)>0){
    Which.error=which(G.transect>Par$G[Mapping])
    for(i in 1:sum(I.error))Par$G[Mapping[Which.error[i]]]=Par$G[Mapping[Which.error[i]]]+G.transect[Which.error[i]]
  }
  for(isp in 1:n.species)Par$N[isp,]=Par$G[isp,]+rpois(n.cells,grp.mean[isp]*Par$G[isp,])
  if(length(spat.ind)==1){
    if(spat.ind==1){
      Par$Eta.bern=0*Par$Eta.bern
      Par$Eta.pois=0*Par$Eta.pois
    }
  }
  else{
    for(isp in 1:n.species){
      if(spat.ind[isp]==1){
        if(i.hurdle==1)Par$Eta.bern[isp,]=0
        Par$Eta.pois[isp,]=0
      }
    }
  }
  Par
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




#' function to convert BOSS MCMC list vector (used in estimation) into an mcmc object (cf. coda package) 
#' @param MCMC list vector providing MCMC samples for each parameter type 
#' @param N.hab.pois.par see help for mcmc_ds.R
#' @param N.hab.bern.par see help for mcmc_ds.R
#' @param Cov.par.n see help for mcmc_ds.R
#' @param Hab.pois.names see help for mcmc_ds.R
#' @param Hab.bern.names see help for mcmc_ds.R
#' @param Cov.names see help for mcmc_ds.R
#' @param fix.tau.nu see help for mcmc_ds.R
#' @param spat.ind see help for mcmc_ds.R
#' @export
#' @keywords MCMC, coda
#' @author Paul B. Conn
convert.BOSS.to.mcmc<-function(MCMC,N.hab.pois.par,N.hab.bern.par,Cov.par.n,Hab.pois.names,Hab.bern.names,Det.names,Cov.names,fix.tau.nu=FALSE,spat.ind=TRUE,point.ind=TRUE){
  require(coda)
  if(min(spat.ind)==0)s.i=0
  i.ZIP=!is.na(N.hab.bern.par)[1]
  n.species=nrow(MCMC$Hab.pois)
  n.iter=length(MCMC$Hab.pois[1,,1])
  n.col=n.species*2+sum(N.hab.pois.par)+(1-s.i)*n.species+(1-fix.tau.nu)*n.species+sum(Cov.par.n)*n.species+1
  if(i.ZIP)n.col=n.col+sum(N.hab.bern.par)+(1-s.i)*n.species #for ZIP model
  n.cells=dim(MCMC$G)[3]
  Mat=matrix(0,n.iter,n.col)
  Mat[,1:n.species]=t(MCMC$N.tot)
  counter=n.species
  col.names=paste("Abund.sp",c(1:n.species),sep='')
  for(isp in 1:n.species){
    Mat[,counter+isp]=rowSums(as.matrix(MCMC$G[isp,,],nrow=n.iter,ncol=n.cells)) #total abundance of groups
    col.names=c(col.names,paste("Groups.sp",isp,sep=''))
  }
  counter=counter+n.species
  for(isp in 1:n.species){  #habitat parameters
    Mat[,(counter+1):(counter+N.hab.pois.par[isp])]=MCMC$Hab.pois[isp,,1:N.hab.pois.par[isp]]
    col.names=c(col.names,paste("Hab.pois.sp",isp,Hab.pois.names[[isp]],sep=''))
    counter=counter+sum(N.hab.pois.par[isp])
  }
  if(i.ZIP){
    for(isp in 1:n.species){  #habitat parameters
      Mat[,(counter+1):(counter+N.hab.bern.par[isp])]=MCMC$Hab.bern[isp,,1:N.hab.bern.par[isp]]
      col.names=c(col.names,paste("Hab.bern.sp",isp,Hab.bern.names[[isp]],sep=''))
      counter=counter+sum(N.hab.bern.par[isp])
    }   
  }
  if(s.i==FALSE){
    Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.eta.pois)
    col.names=c(col.names,paste("tau.eta.pois.sp",c(1:n.species),sep=''))
    counter=counter+n.species
  }
  if(s.i==FALSE & i.ZIP){
    Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.eta.bern)
    col.names=c(col.names,paste("tau.eta.bern.sp",c(1:n.species),sep=''))
    counter=counter+n.species
  }
  if(fix.tau.nu==FALSE){
    Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.nu)
    col.names=c(col.names,paste("tau.nu.sp",c(1:n.species),sep=''))
    counter=counter+n.species
  }
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
  Mat[,counter+1]=MCMC$Thin.pl
  col.names=c(col.names,"Thin.index")
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

#' function to sample species using apply
#' @param probability weights for selecting each species
#' @param n.species number of species to sample from
#' @return A sampled species value
#' @export
#' @keywords species classification
#' @author Paul B. Conn
sample_species<-function(probs,n.species)sample(c(1:n.species),1,prob=probs)

