setwd("~/proteinevoutk20/pkg/Plot/aadist_vs_phi")
data<-read.table("rokas_data_frame.txt", stringsAsFactors=F)
data<-data[-which(is.na(data[,2])),]
data<-data[-which(is.na(data[,4])),]
data$g<-data$gPhi/data$Phi

PlotVsExpression<-function(data) {
  x.val<-data$Phi
  axis.offset<-0.7
  plot(c(0.5*min(x.val),max(x.val)),c(0,1), type="n", bty="n",xlab="expression", ylab="proportion different from observed aa", yaxt="n", log="x")
  axis(2,pos=axis.offset*min(x.val))
  points(x.val, data$pdiff_us, col=rgb(1,0,0,0.7),pch=16)
  points(rep(axis.offset*min(x.val), length(x.val)), data$pdiff_us, col=rgb(1,0,0,0.7),pch="-")
  points(x.val, data$pdiff_wag, col=rgb(0,0,1,0.7),pch=15)
  points(rep(axis.offset*min(x.val), length(x.val)), data$pdiff_wag, col=rgb(0,0,1,0.7),pch="-")
  us<-density(data$pdiff_us,from=0, to=1)
  wag<-density(data$pdiff_wag, from=0, to=1)
  #lines(x=3*us$y/max(us$y)+axis.offset*min(x.val), y=us$x, col=rgb(1,0,0))
  #lines(x=axis.offset*min(x.val) - 3*us$y/max(us$y), y=us$x, col=rgb(1,0,0))
  
  polygon(x=c(3*us$y/max(us$y, wag$y)+axis.offset*min(x.val), (axis.offset*min(x.val)-3*us$y/max(us$y, wag$y))), y=c(us$x, (us$x)), col=rgb(1,0,0,0.4), border=NA)
  polygon(x=c(3*wag$y/max(us$y, wag$y)+axis.offset*min(x.val), (axis.offset*min(x.val)-3*wag$y/max(us$y, wag$y))), y=c(wag$x, (wag$x)), col=rgb(0,0,1,0.4), border=NA)
}

PlotUsVsWag<-function(data) {
  #note that color is inferred magnitude of selection
  logGPhi<-log(data$gPhi)
  rescaled.selection<-(logGPhi-min(logGPhi))/max(logGPhi-min(logGPhi))
  
  
 plot(c(0,1),c(0,1),type="n", bty="n", xlab="WAG", ylab="New model", yaxt="n", xaxt="n", main="Proportion simulated aa\nmatching observed aa")
# points(1-data$pdiff_wag, 1-data$pdiff_us, pch=16)
 #text(1-data$pdiff_wag, 1-data$pdiff_us, labels=round(10*log(data$gPhi),0),col=hsv(1,.8,sqrt(rescaled.selection)))
 points(1-data$pdiff_wag, 1-data$pdiff_us, pch=16,col=hsv(1,.8,sqrt(rescaled.selection)))
  points(rep(0, length(data$pdiff_us)), 1-data$pdiff_us, pch="+",col=hsv(1,.8,sqrt(rescaled.selection)))
  points(1-data$pdiff_wag, rep(0, length(data$pdiff_wag)), pch="+",col=hsv(1,.8,sqrt(rescaled.selection)))
  
 lines(x=c(0,1), y=c(0,1), lty="dashed")
  axis(1,pos=0)
  axis(2,pos=0)
}

PlotUsVsWag(data)
#PlotVsExpression(data)