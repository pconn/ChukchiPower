library(maptools)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(gridExtra)

laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
AK_albers_proj <-paste("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0"," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
Polar_stereo<-paste("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
#Polar_stereo<-paste("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-155 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")

load('Grid_chukchi.Rda')
#Grid=Grid.chukchi$Grid
Grid=spTransform(Grid.chukchi$Grid,CRS(Polar_stereo))

Area.hab=1-Grid[["land_cover"]]

long1<-readShapeSpatial("./shapefiles/to_offshore_1",proj=CRS(AK_albers_proj))
names(long1)=tolower(names(long1))
n=length(long1)
long2<-readShapeSpatial("./shapefiles/to_offshore_2",proj=CRS(AK_albers_proj))
names(long2)=tolower(names(long2))
long2<- spChFIDs(long2,as.character((n+1):(n+length(long2))))   #for all datasets after the first, need to come up with new line identifiers
n=n+length(long2)
short1<-readShapeSpatial("./shapefiles/to_coastal_1",proj=CRS(AK_albers_proj))
names(short1)=tolower(names(short1))
short1<- spChFIDs(short1,as.character((n+1):(n+length(short1))))   #for all datasets after the first, need to come up with new line identifiers
n=n+length(short1)
short2<-readShapeSpatial("./shapefiles/to_coastal_2",proj=CRS(AK_albers_proj))
names(short2)=tolower(names(short2))
short2<- spChFIDs(short2,as.character((n+1):(n+length(short2))))   #for all datasets after the first, need to come up with new line identifiers
n=n+length(short2)
short3<-readShapeSpatial("./shapefiles/to_coastal_3",proj=CRS(AK_albers_proj))
names(short3)=tolower(names(short3))
short3<- spChFIDs(short3,as.character((n+1):(n+length(short3))))   #for all datasets after the first, need to come up with new line identifiers
n=n+length(short3)

long1<-spTransform(long1,CRS(Polar_stereo))
long2<-spTransform(long2,CRS(Polar_stereo))
short1<-spTransform(short1,CRS(Polar_stereo))
short2<-spTransform(short2,CRS(Polar_stereo))
short3<-spTransform(short3,CRS(Polar_stereo))


plot(Grid.chukchi$Grid)
Grid=spChFIDs(Grid, as.character(c(1:length(Grid))))

#formulate Mapping and area surveyed for each configuration
Effort=vector("list",9)
F.list=Effort
F.list[[1]]=long1[c(1,2,3,6),]
F.list[[2]]=spRbind(long1[c(1,7),],short1[c(1,2),])
F.list[[3]]=long1
F.list[[4]]=spRbind(long1[c(1,2,3,4,5,7),],short1[c(1,2),])
F.list[[5]]=spRbind(spRbind(long1[c(1,2,3,6),],short1[1:2,]),short3[2:3,])
F.list[[6]]=spRbind(long1,long2[c(1,2,4,7),])
F.list[[7]]=spRbind(spRbind(long1,long2[1:2,]),short1[1:2,])
F.list[[8]]=spRbind(spRbind(long1,short1[1:2,]),short3[2:3,])
F.list[[9]]=spRbind(spRbind(long1[c(1,2,3,4,5,7),],short1),short3)



make.plot=TRUE
if(make.plot){
  center=apply(bbox(Grid.chukchi$Grid), 1, mean) #gCentroid(Grid.chukchi$Grid)
  Grid.rotate=elide(Grid.chukchi$Grid,rotate=-90,center=center)
  F.plot=F.list
  for(ifl in 1:9)F.plot[[ifl]]=elide(F.list[[ifl]],rotate=-90,center=center)  
  
  #Cur.df=data.frame(gCentroid(Grid.rotate,byid=TRUE))
  #colnames(Cur.df)=c("Easting","Northing")
  #tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  #p1=ggplot(Cur.df)+aes(Easting,Northing)+geom_raster(fill="gray")+tmp.theme
  #p1
  
  
  pdf(file="sim_flights_Chukchi.pdf") 
  par(mfrow=c(4,3),mar=c(0.1,0.1,0.1,0.1),mgp=c(0,0,0))
  plot(Grid.rotate,col='darkgray')
  lines(F.plot[[1]],col='white',lwd=1.5)
  text(-2300000,1400000,"A1",cex=1.3)
  plot(Grid.rotate,col='darkgray')
  lines(F.plot[[3]],col='yellow',lwd=1.5)
  text(-2300000,1400000,"B1",cex=1.3)
  plot(Grid.rotate,col='darkgray')
  lines(F.plot[[6]],col='green',lwd=1.5)
  text(-2300000,1400000,"C1",cex=1.3)  
  plot(Grid.rotate,col='darkgray')
  lines(F.plot[[2]],col='white',lwd=1.5)
  text(-2300000,1400000,"A2",cex=1.3) 
  plot(Grid.rotate,col='darkgray')
  lines(F.plot[[4]],col='yellow',lwd=1.5)
  text(-2300000,1400000,"B2",cex=1.3)  
  plot(Grid.rotate,col='darkgray')
  lines(F.plot[[7]],col='green',lwd=1.5)
  text(-2300000,1400000,"C2",cex=1.3)   
  frame()
  plot(Grid.rotate,col='darkgray')
  lines(F.plot[[5]],col='yellow',lwd=1.5)
  text(-2300000,1400000,"B3",cex=1.3)   
  plot(Grid.rotate,col='darkgray')
  lines(F.plot[[8]],col='green',lwd=1.5)
  text(-2300000,1400000,"C3",cex=1.3)  
  frame()
  frame()
  plot(Grid.rotate,col='darkgray')
  lines(F.plot[[9]],col='green',lwd=1.5)  
  text(-2300000,1400000,"C4",cex=1.3) 
  dev.off()
}



for(ifl in 1:9){
  int=gIntersects(F.list[[ifl]],Grid,byid=TRUE)
  vec <- vector(mode="list", length=nrow(F.list[[ifl]]))
  for (i in seq(along=vec)) vec[[i]] <- try(gIntersection(F.list[[ifl]][i,],Grid[int[,i],], byid=TRUE))
  out <- do.call("rbind", vec)
  rn <- row.names(out)
  nrn <- do.call("rbind", strsplit(rn, " "))
  Length.df <- data.frame(Fl=nrn[,1], poly=as.numeric(as.character(nrn[,2])), len=gLength(out,byid=TRUE))
  #combine data when different flights go through the same cell
  Which.dup=which(duplicated(Length.df[,"poly"]))
  if(length(Which.dup)>0){
    for(i in 1:length(Which.dup)){
      which.first=which(Length.df[,"poly"]==Length.df[Which.dup[i],"poly"])[1]
      Length.df[which.first,"len"]=Length.df[which.first,"len"]+Length.df[Which.dup[i],"len"]
    }
    Length.df=Length.df[-Which.dup,]
  }
  
  #calculate area surveyed for each of these flight_segment * grid cell combos
  #Row.index=rep(0,nrow(Length.df))
  #for(i in 1:length(Row.index))Row.index[i]=which(F.list[[ifl]][["seg_id"]]==as.character(Length.df[i,"Fl"]))

  
  Diam=280  #280 meter swath width for aero commander
  Length.df["area"]=Length.df[,"len"]*0.001*(Diam*.001)/625  #proportional area surveyed for each cell
  
  Effort[[ifl]]=Length.df[,c(2,4)]
  Order=rank(Effort[[ifl]][,"poly"])
  Effort[[ifl]][Order,]=Effort[[ifl]]
  rownames(Effort[[ifl]])=c(1:nrow(Effort[[ifl]]))
  
}

save(Effort,file="Sim_effort.Rda")
  


