
data.run1.roi1 <- list(y=t(fmridata$ROI.timeseriesN[[1]][[1]]), x = fmridata$X[[1]][[1]][,1:2], 
                       Omega=diag(rep(0.1,20)), T=128 , v=20,  mean = c(0,0,0), prec=diag(rep(1.0E-6,3)))

mu.zero <- matrix(1, nrow=20, ncol=128)

inits <- function() {list(beta=c(0,0,0), tauC = t(mu.zero), b0=1, b1=1, w=t(mu.zero),
                          R = diag(rep(1,20)), mu=t(mu.zero), muy=t(mu.zero), w0=rep(1,20))}

parameters <- c("beta") 

setwd("/Users/Michaelzmc/Dropbox/Study/course/Bios372/Project/")

fmrimodel <- function(){
  for(i in 1:v){
    y[1,1:v] ~ dmnorm(muy[1,1:v], R[ , ])
    muy[1,i] <- mu.beta[1] + mu.beta[2] * x[1,1] + mu.beta[3] * x[1,2] + w[1,i]
    w[1,i] ~ dnorm(mu[1,i], tauC[1,i])
    mu[1,i] <- b0 + b1 * w0[i]
    tauC[1,i] ~ dgamma(0.001, 0.001)
    w0[i] ~ dbeta(1,1)

    for(t in 2 : T) {
      y[t,1:v] ~ dmnorm(muy[t,1:v], R[ , ])
      muy[t,i] <- mu.beta[1] + mu.beta[2] * x[t,1] + mu.beta[3] * x[t,2] + w[t,i]
      w[t,i] ~ dnorm(mu[t,i], tauC[t,i])
      mu[t,i] <- b0 + b1 * w[t-1,i]
      tauC[t,i] ~ dgamma(0.001, 0.001)
    }
  }
  
  mu.beta[1 : 3] ~ dmnorm(mean[], prec[ , ])
  R[1 : v , 1 : v] ~ dwish(Omega[ , ], v)
  b0 ~ dnorm(0,0.001)
  b1 ~ dnorm(0,0.001)
}


# write the model code out to a file
write.model(fmrimodel, "fmrimodel.txt")
model.file <- paste(getwd(),"fmrimodel.txt", sep="/")

#set the WINE working directory and the directory to OpenBUGS - change the OpenBUGS.exe location as necessary
WINE="/opt/local/bin/wine"
WINEPATH="/opt/local/bin/winepath"
OpenBUGS.pgm="/Users/Michaelzmc/.wine/drive_c/Program Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe"

fmri.sim1 <- bugs(data.run1.roi1, inits, parameters, n.chains=1, n.iter=10,
                  model.file = model.file, OpenBUGS.pgm=OpenBUGS.pgm, 
                  WINE=WINE, WINEPATH=WINEPATH, useWINE=T, debug=TRUE)
