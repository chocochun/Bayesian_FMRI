### uncorrelated voxel
v=2
data.run1.roi1 <- list(y1=fmridata$ROI.timeseriesN[[1]][[1]][1,], 
                       y2=fmridata$ROI.timeseriesN[[1]][[1]][2,], 
                       x1 = fmridata$X[[1]][[1]][,1], 
                       x2 = fmridata$X[[1]][[1]][,2], 
                       x3 = fmridata$X[[1]][[1]][,3], 
                       T=128 , v=v,  mean = c(0,0,0), prec=diag(rep(1.0E-6,3)))

mu.zero <- matrix(1, nrow=v, ncol=128)

inits <- function() {list(mu.beta=c(0,0,0), tauY=1, 
                          alpha=rep(1.0,v), beta=rep(0.95v), tau2=1.0,
                          w=t(mu.zero),
                          mu=t(mu.zero), muy=t(mu.zero), w0=rep(20,v)}
parameters <- c("mu.beta") 

setwd("/Users/Michaelzmc/Dropbox/Study/course/Bios372/Project/")

fmrimodel <- function(){
    # ar(1) for voxel i
    mu[1,1] <- alpha + beta*w0
    w[1,1]  ~  dnorm(mu[1,i],tau2[i])
    muy[1,1] <-  mu.beta[1] * x1[1] + mu.beta[2] * x2[1] +  mu.beta[3] * x3[1] + w[1,i]
    y1[1] ~ dnorm(muy[1,i], tauY[i])
    
    mu[1,i] <- alpha[i] + beta[i]*w0[i]
    w[1,i]  ~  dnorm(mu[1,i],tau2[i])
    muy[1,i] <-  mu.beta[1] * x1[1] + mu.beta[2] * x2[1] +  mu.beta[3] * x3[1] + w[1,i]
    y1[1] ~ dnorm(muy[1,i], tauY[i])
    
    for(t in 2 : T) {
      mu[t,i] <- alpha[i] + beta[i]*w[t-1,i]
      w[t,i]  ~ dnorm(mu[t,i],tau2[i])
      muy[t,i] <-  mu.beta[1] * x1[t] + mu.beta[2] * x2[t] +  mu.beta[3] * x3[t] + w[t,i]
      y[t,i] ~ dnorm(muy[t,i], tauY[i]) 
    }
    

  
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


fmri.sim1 <- bugs(data.run1.roi1, inits, parameters, n.chains=1, n.iter=10,
                  model.file = model.file, OpenBUGS.pgm=OpenBUGS.pgm, 
                  WINE=WINE, WINEPATH=WINEPATH, useWINE=T, debug=TRUE)
