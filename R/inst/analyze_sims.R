#analyze power simulation results
set.seed(23456)
n.sims=10
start.iter=500
end.iter=2000

Result.table=vector("list",2)
Est.table=Result.table
Lam.true=Result.table
for(isp in 1:2){
  Result.table[[isp]]=matrix(0,9,4)
  colnames(Result.table[[isp]])=c("Bias","MSE","CV","Coverage90")
  Est.table[[isp]]=data.frame(matrix(0,90,3))
  colnames(Est.table[[isp]])=c("Flight","Simulation","Estimate")
  Lam.true[[isp]]=rep(0,n.sims)
}

#load('./sim_data/Out_Fl3_sim1.Rda')
cur.pl=1
for(isim in 1:n.sims){
  load(paste("./sim_data/Lambda_chukchi_",isim,".Rda",sep=''))
  for(isp in 1:2){
    lam.true=rpois(1,sum(Sim.abund$Lambda.true[isp,]))
    Lam.true[[isp]][isim]=lam.true
    for(ifl in 1:9){
      load(paste("./sim_data/Out_Fl",ifl,"_sim",isim,".Rda",sep=''))
      Result.table[[isp]][ifl,1]=Result.table[[isp]][ifl,1]+(median(Out$MCMC[start.iter:end.iter,isp])-lam.true)/lam.true
      Result.table[[isp]][ifl,2]=Result.table[[isp]][ifl,2]+(median(Out$MCMC[start.iter:end.iter,isp])-lam.true)^2
      Result.table[[isp]][ifl,3]=Result.table[[isp]][ifl,3]+sqrt(var(Out$MCMC[start.iter:end.iter,isp]))/median(Out$MCMC[,isp])
      Quant=quantile(Out$MCMC[start.iter:end.iter,isp],c(0.05,0.95))
      Result.table[[isp]][ifl,4]=Result.table[[isp]][ifl,4]+(Quant[1]<lam.true & lam.true<Quant[2])
      Est.table[[isp]][cur.pl,]=c(ifl,isim,median(Out$MCMC[start.iter:end.iter,isp]))
      cur.pl=cur.pl+1
    }
  }
}

#average over number of sims
for(isp in 1:2)Result.table[[isp]]=Result.table[[isp]]/n.sims

