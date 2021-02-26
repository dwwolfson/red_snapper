##
## Example of using a CTMC model for movement
##
## Steps:
##  1. Fit Quasi-Continuous Path Model to telemetry data (done using Buderman et al 2015)
##  2. Create covariate raster objects (the CTMC will be on the raster
##     grid cells)
##  3. Impute a quasi-continuous path
##  4. Turn quasi-continuous path into a CTMC discrete-space path
##  5. Turn discrete-space path into latent Poisson GLM format
##  6. Fit a Poisson GLM model to the data
##


library(ctmcmove)
load("lion.pair.Rdata")
attach(lion.pair)
lion1$ID
lion2$ID



AF79=lion.pair[[1]]
AM80=lion.pair[[2]]

## subsetting time (three-week observation window)

mintime=4100.085
maxtime=4120.085

time.idx.1=which(AF79$locs.orig[,3] >= mintime & AF79$locs.orig[,3]<maxtime)
time.idx.2=which(AM80$locs.orig[,3] >= mintime & AM80$locs.orig[,3]<maxtime)


#########################################################
##
## 1.  fitting functional movement model to telemetry data
##
#########################################################

AM80$locs=AM80$locs.orig[time.idx.2,]
xy=AM80$locs[,1:2]
t=AM80$locs[,3]

library(fda)

## define time points where the quasi-continuous path will be sampled
knots = seq(mintime,maxtime,by=1/24/12)
## create B-spline basis vectors used to approximate the path
b=create.bspline.basis(c(mintime,maxtime),breaks=knots,norder=3)
## define the sequence of times on which to sample the imputed path
tpred=seq(mintime,maxtime,by=1/24/60)


## Fit latent Gaussian model using MCMC
out=mcmc.fmove(xy,t,b,tpred=knots,QQ="CAR2",n.mcmc=400,a=1,r=1,num.paths.save=10,sigma.fixed=NA)
str(out)

## plot 3 imputed paths
plot(xy,type="p",xlim=c(455000,462000),pch=20,cex=1)
points(out$pathlist[[1]]$xy,col="red",type="l",lwd=3)
points(out$pathlist[[2]]$xy,col="blue",type="l",lwd=3)
points(out$pathlist[[3]]$xy,col="green",type="l",lwd=3)
points(xy,type="p",pch=20,cex=2,lwd=3)


## pdf("telem.pdf",width=2,height=2)
## par(mar=c(1,1,1,1))
## plot(xy,type="b",pch=20,axes=F,lwd=2,xlab="",ylab="",cex=3)

## savePlot("telem.jpg")
## dev.off()

## pdf("path1.pdf",width=2,height=2)
## par(mar=c(1,1,1,1))
## plot(out$pathlist[[1]]$xy,type="l",axes=F,col="red",lwd=2,xlab="",ylab="")
## points(xy,type="p",pch=20,cex=1)
## dev.off()

## pdf("path2.pdf",width=2,height=2)
## par(mar=c(1,1,1,1))
## plot(out$pathlist[[2]]$xy,type="l",axes=F,col="blue",lwd=2,xlab="",ylab="")
## points(xy,type="p",pch=20,cex=1)
## dev.off()

## pdf("path3.pdf",width=2,height=2)
## par(mar=c(1,1,1,1))
## plot(out$pathlist[[1]]$xy,type="l",axes=F,col="green",lwd=2,xlab="",ylab="")
## points(xy,type="p",pch=20,cex=1)
## dev.off()




#########################################################
##
## 2. Get Rasters for static covariates
##
#########################################################

intercept=stack.cropped[[1]]
values(intercept) <- 1

NotForest=stack.cropped[[1]]+stack.cropped[[2]]+stack.cropped[[3]]+stack.cropped[[4]]+stack.cropped[[6]]
names(NotForest) <- "NotForest"



##
##
## Get Rasters for dynamic covariate (distance to nearest potential kill site)
##
##

source("FindCenters.r")
cl=FindCenters(AM80,mintime=mintime,maxtime=maxtime,max.dist=200,max.time.dist=27,nightonly=T)
cl
kill.sites.rast=intercept
values(kill.sites.rast) <- NA
cent.cells=cellFromXY(kill.sites.rast,cl$clusters)
values(kill.sites.rast)[cent.cells] <- 0
d2kill=distance(kill.sites.rast)
plot(d2kill)

path=out$pathlist[[1]]

plot(d2kill,main="Background: Distance to Potential Kill Site")
points(path$xy,type="l",col="blue")
points(AM80$locs,type="p",pch=20)
##savePlot("d2k.jpg")

dynamic.stack=stack(d2kill,stack.cropped[[7]])
names(dynamic.stack) <- c("grad.Dist2Kill","grad.Elevation")

static.cover=stack(intercept,d2kill,NotForest)
names(static.cover) <- c("Int","Dist2Kill","NotForest")

plot(static.cover)
plot(dynamic.stack)

##
##
##
##
##  Example #1: Constant response (in time) to covariates
##
##
##
##


#########################################################
##
## 4. Turn continuous space path into a CTMC discrete space path
##
#########################################################

path=out$pathlist[[1]]
ctmc=path2ctmc(path$xy,path$t,dynamic.stack)


#########################################################
##
## 5. Turn CTMC discrete path into latent Poisson GLM data
##
#########################################################

glm.data=ctmc2glm(ctmc,static.cover,dynamic.stack)

str(glm.data)
summary(glm.data)

## now repeat 3-5 for 10 imputations
for(i in 2:10){
        path=out$pathlist[[i]]
        ctmc=path2ctmc(path$xy,path$t,intercept)
        glm.data=rbind(glm.data,ctmc2glm(ctmc,static.cover,dynamic.stack))
    }

str(glm.data)


## remove transitions that are nearly instantaneous
idx.0=which(glm.data$tau<10^-5)
idx.0
if(length(idx.0)>0){
    glm.data=glm.data[-idx.0,]
}

str(glm.data)

#########################################################
##
## 6. Fit models
##
#########################################################


fit=glm(z~NotForest+Dist2Kill+grad.Dist2Kill+grad.Elevation+crw,family="poisson",offset=log(tau),weights=rep(1/10,nrow(glm.data)),data=glm.data)
summary(fit)

library(mgcv)

fit=gam(z~s(Dist2Kill)+s(grad.Dist2Kill)+s(I(t%%1),bs="cc")+crw,weights=rep(1/10,nrow(glm.data)),family="poisson",offset=log(tau),data=glm.data)
summary(fit)

par(mfrow=c(1,3))
plot(fit)

