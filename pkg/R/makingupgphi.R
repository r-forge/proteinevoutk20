#true.phi<-runif(106, min=min(rokasPhi$SEMPPR, na.rm=TRUE), max=max(rokasPhi$SEMPPR, na.rm=TRUE))
true.phi <- rlnorm(100,mean(cphi$SEMPPR),sd(cphi$SEMPPR)) # random samples from log normal distribution
gphi<-rnorm(length(true.phi), 10*true.phi, 0.1)
print(range(true.phi))
print(range(gphi))
#plot(gphi, true.phi)
bad.phi<-abs(rnorm(length(true.phi), true.phi, 1))
awful.phi<-runif(100,min=exp(min(cphi$SEMPPR)),max=exp(max(cphi$SEMPPR)))
#plot(true.phi, bad.phi)
par(mfrow=c(2,3))
phi.range<-range(c(true.phi, bad.phi, awful.phi))
g.range<-range(c(gphi/true.phi, gphi/bad.phi, gphi/awful.phi))
plot(true.phi, gphi, xlim=phi.range)
plot(bad.phi, gphi, xlim=phi.range)
plot(awful.phi, gphi, xlim=phi.range)
plot(true.phi, gphi/true.phi, xlim=phi.range, ylim=g.range, log="y")
plot(bad.phi, gphi/bad.phi, xlim=phi.range, ylim=g.range, log="y")
plot(awful.phi, gphi/awful.phi, xlim=phi.range, ylim=g.range, log="y")

plot(true.phi, gphi)
plot(bad.phi, gphi)
plot(awful.phi, gphi)
plot(true.phi, gphi/true.phi, log="y")
plot(bad.phi, gphi/bad.phi, log="y")
plot(awful.phi, gphi/awful.phi,log="y")
#plot(true.phi, gphi)

phi <- rokasPhi
phi$Ingolia <- NULL
phi$gene <- NULL
phi.log <- log(phi)
phi.log.center <- phi.log
for(i in 1:3){
  phi.log.center[,i] <- phi.log[,i] - mean(phi.log[,i],na.rm=TRUE)
}
phi.log.mean <- colMeans(phi.log.center,na.rm=TRUE)

nas <- which(is.na(phi.log.center$SEMPPR))
gphi <- s_max[-nas]
phi.log.center <- phi.log.center[-nas,]
cphi <- phi.log.center #rename the centered log phi values as cphi -- shorter name