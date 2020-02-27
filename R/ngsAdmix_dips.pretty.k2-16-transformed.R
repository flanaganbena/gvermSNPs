rm(list=ls())
pdf('output/ngsAdmix_dips.pretty.k2-16-transformed.pdf',width=12,height=10)
par(mfrow=c(17,1),mar=c(0,0,0,0),xpd = TRUE)
ks <- c("k02run1","k03run1","k04run1","k05run1","k06run1","k07run1","k08run1","k09run1","k10run1","k11run1","k12run1","k13run1","k14run1","k15run1","k16run1")
indmeta <- read.csv('data/dipsmill.351inds.txt',header=F)[,1]
pop <- substr(indmeta,1,3)
pop[pop=="fjs"] <- "cfj"
pop <- factor(pop)
siteorder.to.use <- c('fut','aka','usu','hay','waj','shk','mou','sou','mng','cnt','hik',
                      'pmo','ptw','eld','bob','tmb', # WUSA
                      'grb','ris','lhp','nyc','uho','wwf','osh','cfj','ebp','sav', # EUSA
                      'dns','dhh','dhd','fpl','fdm','fpa','gal','ave','far') #EUR
siteorder <- match(pop,siteorder.to.use)
meta <- read.csv('data/meta-SNP.pops_edited.V2.csv')
pop.correctednames <- meta$Site.Abb.2[match(siteorder.to.use,meta$Site.Abb)]

kcol <- c("black", #1
          "red",#2
          "darkgreen",#3
          "gainsboro",#4
          "yellow",#5
          "deepskyblue",#6
          "brown",#7
          "dodgerblue4",#8
          "darkorchid",#9
          "burlywood",#10
          "aquamarine",#11
          "chocolate",#12
          "khaki",#13
          "darksalmon",#14
          "floralwhite",#15
          "springgreen",#16
          "goldenrod",
          "firebrick",
          "deeppink",
          "forestgreen")
colorder <- list(
  c(1,2), #ks=2
  c(1,2,3), #ks=3
  c(2,1,4,3), #ks=4
  c(5,3,1,4,2), #ks=5
  c(4,3,1,5,2,6), #ks=6
  c(2,1,7,4,6,5,3), #ks=7
  c(5,1,2,3,6,4,8,7), #ks=8
  c(7,4,5,8,9,3,1,2,6), #ks=9
  c(3,1,5,6,9,4,2,10,7,8), #ks=10
  c(3,7,11,6,4,10,2,9,5,1,8), #ks=11
  c(10,2,6,3,5,1,12,8,7,11,9,4), #ks=12
  c(12,1,11,6,3,4,13,8,7,10,9,2,5), #ks=13
  c(3,7,8,9,2,5,1,6,12,10,11,4,14,13), #ks=14
  c(8,4,10,11,9,5,6,7,3,12,1,14,2,15,13), #ks=15
  c(5,16,11,1,15,6,2,4,9,14,10,3,12,7,8,13), #ks=16
  1:17, #ks=17
  1:18, #ks=18
  1:19, #ks=19
  1:20 #ks=20
)



for (k in 1:length(ks))
{
  dat <- read.delim(paste("data/NGSadmix.runs/",ks[k],".qopt",sep=""),sep=" ",header = F)
  dat <- dat[order(siteorder),]
  dat <- dat[,-(dim(dat)[2])]
  dat <- dat[,order(colorder[[k]])]
  fig <- barplot(t(dat),col=kcol[1:length(colorder[[k]])],space=0,border=NA,xlab="",ylab="",
          names.arg = rep("",nrow(dat)),horiz=F)#,main=substr(ks[k],1,3))
  mtext(substr(ks[k],1,3),side=2,line=-1.5,at=.5)
  for(i in 1:length(siteorder.to.use))
  {
    pop2 <- pop[order(siteorder)]
    x <- 1:dim(dat)[1]
    xtmp <- x[pop2==siteorder.to.use[i]]
    segments(max(xtmp),1,max(xtmp),-0.05,col="white")
  }
  if(k==length(ks)){
    fig <- barplot(t(dat),col="white",space=0,border=NA,xlab="",ylab="",
                   names.arg = rep("",nrow(dat)),horiz=F)
    for(i in 1:length(siteorder.to.use))
    {
      pop2 <- pop[order(siteorder)]
      x <- 1:dim(dat)[1]
      xtmp <- x[pop2==siteorder.to.use[i]]
      text(x=mean(xtmp),y=.5,pop.correctednames[i],cex=1)
  }}
}

dev.off()