### points of ngsAdmix output of SNPs
rm(list=ls())
library(maps);library(fields);library(reshape);library(mapdata);library(colorRamps);library(plotrix)
dat <- read.delim('data/NGSadmix.runs/k07run1.qopt',sep=" ",header=F)
meta <- read.csv('data/dipsmill.351inds.txt',header=F)[,1]
pop <- substr(meta,1,3)
pop <- factor(pop)
siteorder.to.use <- c('fut','aka','usu','hay','waj','shk','mou','sou','mng','cnt','hik',
                      'pmo','ptw','eld','bob','tmb', # WUSA
                      'grb','ris','lhp','nyc','uho','wwf','osh','cfj','ebp','sav', # EUSA
                      'dns','dhh','dhd','fpl','fdm','fpa','gal','ave','far') #EUR
siteorder <- match(pop,siteorder.to.use)

kcol <- c("black",
          "red",
          "darkgreen",
          "gainsboro",
          "yellow",
          "deepskyblue",
          "brown")#,
colorder <- c(5,1,2,3,6,4,8,7) #ks=8
colorder <- c(2,1,7,4,6,5,3) #ks=7
dat <- dat[,-(dim(dat)[2])]
dat <- dat[,order(colorder)]
colnames(dat) <- paste("cluster",1:7,sep="")
dat$ID <- pop

#### summarize the mean cluster contribution per site
md <- melt(dat,id="ID")
xbar <- cast(md,ID~variable,mean)

### add lat and long to dataset and do plots
metapop <- read.csv("data/meta-SNP.pops_edited.V2.csv")
xbar$lat <- metapop$Latitude[match(xbar$ID,metapop$Site.Abb)]
xbar$lon <- metapop$Longitude[match(xbar$ID,metapop$Site.Abb)]
xbar$latpie <- metapop$lat.name[match(xbar$ID,metapop$Site.Abb)]
xbar$lonpie <- metapop$long.name[match(xbar$ID,metapop$Site.Abb)]
#### move pies
xbar[xbar$ID=="mou",c("latpie","lonpie")] <- c(39.9535,144.2987)
xbar[xbar$ID=="pmo",c("latpie","lonpie")] <- c(56.7,-124.1)
xbar[xbar$ID=="ptw",c("latpie","lonpie")] <- c(52.5,-131.2)
xbar[xbar$ID=="eld",c("latpie","lonpie")] <- c(47.7,-133.9)
xbar[xbar$ID=="bob",c("latpie","lonpie")] <- c(37.7,-132.8)
xbar[xbar$ID=="tmb",c("latpie","lonpie")] <- c(31.0,-130.1)
xbar[xbar$ID=='sav',c('latpie','lonpie')] <- c(23.733,-82)
xbar[xbar$ID=='ebp',c('latpie','lonpie')] <- c(25.733,-75.478)
xbar[xbar$ID=='cfj',c('latpie','lonpie')] <- c(30.614,-73.513)
xbar[xbar$ID=='osh',c('latpie','lonpie')] <- c(34.474,-69.293)
xbar[xbar$ID=='wwf',c('latpie','lonpie')] <- c(35.75,-61.313)
xbar[xbar$ID=='uho',c('latpie','lonpie')] <- c(39.853,-65.193)
xbar[xbar$ID=='nyc',c('latpie','lonpie')] <- c(43.70,-78.577)
xbar[xbar$ID=='lhp',c('latpie','lonpie')] <- c(48.691,-75.201)
xbar[xbar$ID=='ris',c('latpie','lonpie')] <- c(51.477,-69.533)
xbar[xbar$ID=='grb',c('latpie','lonpie')] <- c(44.753,-65.916)
xbar[xbar$ID=='far',c('latpie','lonpie')] <- c(34.033,-15.911)
xbar[xbar$ID=='ave',c('latpie','lonpie')] <- c(39.437,-17.5)
xbar[xbar$ID=='gal',c('latpie','lonpie')] <- c(43.930,-15.223)
xbar[xbar$ID=='fpa',c('latpie','lonpie')] <- c(47.457,-9.720)
xbar[xbar$ID=='fdm',c('latpie','lonpie')] <- c(55.780,-4.905)
xbar[xbar$ID=='fpl',c('latpie','lonpie')] <- c(52.290,-11.096)
xbar[xbar$ID=='dns',c('latpie','lonpie')] <- c(58.465,3.005)
xbar[xbar$ID=='dhd',c('latpie','lonpie')] <- c(59.002,10.571)
xbar[xbar$ID=='dhh',c('latpie','lonpie')] <- c(57.391,17)
###### convert colors #########
xbar2 <- round(xbar[,substr(colnames(xbar),1,7)=="cluster"],2)*1000
xbar2[xbar2==0] <- 1

### Native range
png('output/ngsAdmix-points.in.color.japan.png',width=10,height=8,units="in",res=500)
map("worldHires","China",bg="lightsteelblue1",fill=T,col="white",xlim=c(125,150),ylim=c(29,47),xaxt="n",yaxt="n",ylab="",xlab="")
map("worldHires","South Korea",add=TRUE,col="white",fill=TRUE)
map("worldHires","North Korea",add=TRUE,col="white",fill=TRUE)
map("worldHires","Japan",add=TRUE,col="white",fill=TRUE)
map("worldHires","USSR",col="white",fill=TRUE,add=TRUE)
#box()
segments(x0 = xbar$lon,y0=xbar$lat,
         x1 = xbar$lonpie,y1=xbar$latpie,lwd=3)

for (i in 1:dim(xbar)[1])
{
floating.pie(xpos = xbar$lonpie[i],ypos = xbar$latpie[i],
             x = as.numeric(xbar2[i,]),col=kcol)
}
dev.off()

### Non-native range
png('output/ngsAdmix-points.in.color.non-native.png',width=10,height=5,units="in",res=500)
map('worldHires',xlim=c(-150,20),ylim=c(19,65),bg="lightsteelblue1",fill=T,col="white")
segments(x0 = xbar$lon,y0=xbar$lat,
         x1 = xbar$lonpie,y1=xbar$latpie,lwd=3)

for (i in 1:dim(xbar)[1])
{
  floating.pie(xpos = xbar$lonpie[i],ypos = xbar$latpie[i],
               x = as.numeric(xbar2[i,]),col=kcol,radius = 3)
}
dev.off()
