library(maptools)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
source('c:/users/paul.conn/git/hierarchicalDS/hierarchicalDS/R/spat_funcs.R') #for rect_adj fcn

laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
AK_albers_proj <-paste("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0"," +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
Polar_stereo<-paste("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")

load('Grid_chukchi.Rda')
Grid=Grid.chukchi$Grid
Area.hab=1-Grid[["land_cover"]]
Grid_chukchi=spChFIDs(Grid, as.character(c(1:length(Grid))))


#read in strata
Strata<-readShapeSpatial("./shapefiles/chukchi_strata",proj=CRS(AK_albers_proj))
names(Strata)=tolower(names(Strata))
n=length(Strata)
strata<-spTransform(Strata,CRS(Polar_stereo))

Strata_lookup=c(12,11,4,10,9,7,6,5,1,8,2,3,4)

#calculate density per strata (calculated as area weighted values for pack and non-pack ice)
D.r=read.csv('BengtsonRinged.csv',header=FALSE,colClasses=c("factor",rep("numeric",8)))
#make an entry for pack and non pack for each strata even if one doesn't exist
Temp=data.frame(matrix(0,24,9))
Temp[1:7,]=D.r[1:7,]
Temp[8,]=D.r[7,]
Temp[9:21,]=D.r[8:20,]
Temp[22,]=D.r[20,]
Temp[23,]=D.r[21,]
Temp[24,]=D.r[21,]
Temp[6,6:9]=Temp[5,6:9]
Temp[16,6:9]=Temp[15,6:9]
Temp[18,6:9]=Temp[17,6:9]
D.r=Temp
D.ringed=rep(0,12)
for(i in 1:12)D.ringed[i]=0.5*(D.r[(i-1)*2+1,4]*D.r[(i-1)*2+1,2]+D.r[(i-1)*2+2,4]*D.r[(i-1)*2+2,2])/(D.r[(i-1)*2+1,2]+D.r[(i-1)*2+2,2])+0.5*(D.r[(i-1)*2+1,8]*D.r[(i-1)*2+1,6]+D.r[(i-1)*2+2,8]*D.r[(i-1)*2+2,6])/(D.r[(i-1)*2+1,6]+D.r[(i-1)*2+2,6])

D.r=read.csv('BengtsonBearded.csv',header=FALSE,colClasses=c("factor",rep("numeric",8)))
#make an entry for pack and non pack for each strata even if one doesn't exist
Temp=data.frame(matrix(0,24,9))
Temp[1:7,]=D.r[1:7,]
Temp[8,]=D.r[7,]
Temp[9:21,]=D.r[8:20,]
Temp[22,]=D.r[20,]
Temp[23,]=D.r[21,]
Temp[24,]=D.r[21,]
Temp[6,6:9]=Temp[5,6:9]
Temp[16,6:9]=Temp[15,6:9]
Temp[18,6:9]=Temp[17,6:9]
D.r=Temp
D.bearded=rep(0,12)
for(i in 1:12)D.bearded[i]=0.5*(D.r[(i-1)*2+1,4]*D.r[(i-1)*2+1,2]+D.r[(i-1)*2+2,4]*D.r[(i-1)*2+2,2])/(D.r[(i-1)*2+1,2]+D.r[(i-1)*2+2,2])+0.5*(D.r[(i-1)*2+1,8]*D.r[(i-1)*2+1,6]+D.r[(i-1)*2+2,8]*D.r[(i-1)*2+2,6])/(D.r[(i-1)*2+1,6]+D.r[(i-1)*2+2,6])


#intersect grid with strata; those that don't intersect will be associated with SpDF #1 (Bengtson strata 12)
I.intersect=gIntersects(strata,Grid_chukchi,byid=TRUE)
Cum.intersect=apply(I.intersect,1,'sum')
Grid.strata=rep(0,length(Grid_chukchi))
Grid.strata[which(Cum.intersect==0)]=1
Which.intersect=which(Cum.intersect>0)
for(i in 1:length(Which.intersect))Grid.strata[Which.intersect[i]]=which(rmultinom(1,1,I.intersect[Which.intersect[i],])==1)

plot(Grid_chukchi,col=Grid.strata)

#attach to SpDF and export
Grid_chukchi@data$D.bearded=D.bearded[Strata_lookup[Grid.strata]]
Grid_chukchi@data$D.ringed=D.ringed[Strata_lookup[Grid.strata]]
plot(Grid_chukchi,col=gray(1-Grid_chukchi$D.bearded/max(Grid_chukchi$D.bearded)))
plot(Grid_chukchi,col=gray(1-Grid_chukchi$D.ringed/max(Grid_chukchi$D.ringed)))

save(Grid_chukchi,file="Grid_chukchi_ExpDensities.Rda")

#establish knots for process convolution
Coords=coordinates(Grid_chukchi)
x.min=min(Coords[,1])-50000
x.max=max(Coords[,1])+50000
y.min=min(Coords[,2])-50000
y.max=max(Coords[,2])+50000

X=x.min+(x.max-x.min)/5*c(0:5)
Y=y.min+(y.max-y.min)/5*c(5:0)
XY=expand.grid(x=X,y=Y)
Knots=SpatialPoints(coords=XY,proj4string=CRS(Polar_stereo))
Distances=gDistance(Knots,Grid_chukchi,byid=TRUE)
Distances=apply(Distances,2,'min')
my.buffer=100000
Which.include=which(Distances<my.buffer)

Knot.cell.distances=gDistance(Knots[Which.include,],Grid_chukchi,byid=TRUE)
diff.x=(x.max-x.min)/5
diff.y=(y.max-y.min)/5
sigma=(diff.x+diff.y)/2

Knot.Adj=rect_adj(6,6)
Knot.Adj=Knot.Adj[Which.include,Which.include]
Q.knot=-Knot.Adj
diag(Q.knot)=apply(Knot.Adj,2,'sum')
Q.knot=Matrix(Q.knot)

K=dnorm(Knot.cell.distances,0,sigma)
K=K/apply(K,1,'sum')
K.data=list(K=K,Q.knot=Q.knot)
save(K.data,file="Knot_cell_distances.Rdata")

pdf('Chukchi_Grid_wKnots.pdf')
plot(Knots[Which.include,],col='red',pch=20)
plot(Grid_chukchi,add=TRUE)
dev.off()

Knots=Knots[Which.include,]
save(Knots,file="Chukchi_Knots_SP.Rda")


#plot(strata,col='red',add=TRUE)
#plot(Grid_chukchi,add=TRUE)
#plot(short_lines,add=TRUE,col="blue")

#plot(Grid_chukchi)
#plot(strata,col='red',add=TRUE)
#plot(Grid_chukchi,add=TRUE)
#plot(long_lines,add=TRUE,col="blue")

