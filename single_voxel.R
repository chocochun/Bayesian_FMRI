#### correlated voxel

v=20

roi1.y <- fmridata$ROI.timeseriesN[[1]][[1]]
roi2.y <- fmridata$ROI.timeseriesN[[2]][[1]]
roi3.y <- fmridata$ROI.timeseriesN[[3]][[1]]


data.run1.roi1 <- list(y=t(roi1.y),
                       x1 = fmridata$X[[1]][[1]][,1], 
                       x2 = fmridata$X[[1]][[1]][,2], 
                       x3 = fmridata$X[[1]][[1]][,3], 
                       Omega=diag(rep(0.1,v)),
                       T=128 ,v=20, mean = c(0,0,0), prec=diag(rep(1.0E-6,3)))

#mu.zero <- matrix(1, nrow=20, ncol=128)
#mu.zero <- rep(0,128)
mu.zero <- matrix(1, nrow=v, ncol=128)

inits <- list(list(mu.beta=c(0,0,0), tauY=1.0, R = diag(rep(1,v)),
                          alpha=1.0, beta=0.95, tau2=1.0, w0=20.0),
              list(mu.beta=c(10,10,10), tauY=2.0, R = diag(rep(0.1,v)),
                   alpha=-1.0, beta=-0.95, tau2=2.0,w0=10.0),
              list(mu.beta=c(-10,-10,-10), tauY=3.0, R = diag(rep(10,v)),
                   alpha=0, beta=0, tau2=3.0,w0=5.0))

parameters <- c("mu.beta") 

setwd("/Users/Michaelzmc/Dropbox/Study/course/Bios372/Project/")

fmrimodel <- function(){
  for(i in 1:v){
    mu[1,i] <- alpha + beta*w0
    w[1,i]  ~  dnorm(mu[1,i],tau2)
    muy[1,i] <-  mu.beta[1] * x1[1] + mu.beta[2] * x2[1] +  mu.beta[3] * x3[1] + w[1,i]
  }  
   
    y[1,1:v] ~ dmnorm(muy[1:v,1], R[ , ])
  
    for(t in 2 : T) {
      for(i in 1:v){
        mu[t,i] <- alpha + beta*w[t-1,i]
        w[t,i]  ~ dnorm(mu[t,i],tau2)
        muy[t,i] <-  mu.beta[1] * x1[t] + mu.beta[2] * x2[t] +  mu.beta[3] * x3[t] + w[t,i]      
      }
      y[t,1:v] ~ dmnorm(muy[1:v,t], R[ , ])
  }
    
  R[1 : v , 1 : v] ~ dwish(Omega[ , ], v)  
  alpha ~ dnorm(0.0, 1.0E-6)
  beta  ~ dnorm(0.0, 1.0E-6)
  tau2  ~ dgamma(0.01, 0.01)
  w0    ~ dnorm(0.0,1.0E-6)
  
  mu.beta[1 : 3] ~ dmnorm(mean[], prec[ , ])
  tauY  ~ dgamma(0.01, 0.01)
}


# write the model code out to a file
write.model(fmrimodel, "fmrimodel.txt")
model.file <- paste(getwd(),"fmrimodel.txt", sep="/")

#set the WINE working directory and the directory to OpenBUGS - change the OpenBUGS.exe location as necessary
WINE="/opt/local/bin/wine"
WINEPATH="/opt/local/bin/winepath"
OpenBUGS.pgm="/Users/Michaelzmc/.wine/drive_c/Program Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe"

fmri.sim1 <- bugs(data.run1.roi1, inits, parameters, n.chains=3, n.iter=1000,
                  model.file = model.file, OpenBUGS.pgm=OpenBUGS.pgm, 
                  WINE=WINE, WINEPATH=WINEPATH, useWINE=T, debug=TRUE)

#beta <- fmri.sim1$sims.list$mu.beta[,1]-fmri.sim1$sims.list$mu.beta[,2]



#sum(fmri.sim1$sims.list$mu.beta[,1]-fmri.sim1$sims.list$mu.beta[,2] > 0) /length(fmri.sim1$sims.list$mu.beta[,1]-fmri.sim1$sims.list$mu.beta[,2])
