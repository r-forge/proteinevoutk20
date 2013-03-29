# #true.phi<-runif(106, min=min(rokasPhi$SEMPPR, na.rm=TRUE), max=max(rokasPhi$SEMPPR, na.rm=TRUE))
# phi.sd <- 0.1
# #true.phi <- phidata$SEMPPR # random samples from log normal distribution
# true.phi <- rnorm(500,1,0.2)
# gs <- 0.5+2*true.phi #Grantham sensitivity
# lgphi<-rnorm(length(true.phi), gs+true.phi, 0.1)
# 
# phi1<-rnorm(length(true.phi), true.phi, phi.sd)
# phi2<-rnorm(length(true.phi), true.phi, phi.sd)
# phi3<-rnorm(length(true.phi), true.phi, phi.sd)
# allphi <- c(phi1,phi2,phi3)
# mkp.data <- data.frame(Beyer=phi1,SEMPPR=phi2,Ing=phi3,lgphi=lgphi)
# pairs(mkp.data)
# phidata.lm(mkp.data)
# lm.mkp <- lm(lgphi~Beyer+SEMPPR+Ing,data=mkp.data)
# print(summary(lm.mkp))
# print(sum(lm.mkp$coefficients[-1])-1)
# 
# 
# 
# lm1 <- lm(lgphi~phi1)
# summary(lm1)
# 
# predict.gphi <- lm1$fitted.values # est(gphi)
# #weighted.phi <- data.matrix(phidata[,1:3]) %*% lm.phi$coefficients[-1] #weighted phi
# est.gval <- lm1$coefficients[1]+(lm1$coefficients[2]-1)*lm1$model[[2]] #est(g): a + b*log(phi)
# plot((predict.gphi-phi1)~phi1,main="est(gphi)/obs(phi) vs. obs(phi) on log scale",xlab="obs(phi)",ylab="est(gphi)/obs(phi)")
# plot(est.gval~phi1,main="est(g) vs. obs(phi) on log scale",xlab="obs(phi)",ylab="est(g)")
# plot(est.gval~gs,main="est(g) vs true(g) on log scale",xlab="true g", ylab="estimated g")
# plot((lgphi-phi1)~phi1,main="obs(gphi)/obs(phi) vs. obs(phi) on log scale",xlab="obs(phi)",ylab="obs(gphi)/obs(phi)")
# 


# mkp.data$Beyer <- log(phi1)
# mkp.data$SEMPPR <- log(phi2)
# mkp.data$Ing <- log(phi3)
# mkp.data$lgphi <- log(gphi)

########These are not simulated under the same assumption as the likelihood approach
true.phi<-runif(106, min=min(rokasPhi$SEMPPR, na.rm=TRUE), max=max(rokasPhi$SEMPPR, na.rm=TRUE))
gphi<-abs(rnorm(length(true.phi), 10*true.phi, 0.1))
print(range(true.phi))
print(range(gphi))
#plot(gphi, true.phi)
bad.phi<-rnorm(length(true.phi), true.phi, 0.1)
awful.phi<-runif(106)
est.phi <- cbind(true.phi,bad.phi,awful.phi)
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
#plot(true.phi, gphi)