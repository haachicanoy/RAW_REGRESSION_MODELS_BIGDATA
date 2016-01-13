
multiProfile <- function(data, profiles, variable,xlim0=NULL) {
  var <- data[,variable]
  if(is.factor(var)) {
    p <- as.data.frame(profiles[variable])
    se <- apply(p, 1, function(x) 1.96*sd(x, na.rm=TRUE)/sqrt(ncol(p)))
    mean <- rowMeans(p)
    mean <- as.data.frame(mean - mean(mean))
    mean <- cbind(mean, table(var))
    mean <- cbind(mean, se)
    names(mean) <- c("Mean", "Class", "Freq", "SE")
    
    limits <- transform(mean,lower = Mean - SE, upper=Mean + SE)
    
    m <- ggplot(mean, aes(x=Class, y=Mean, fill=Freq))
    m <- m + geom_bar(stat="identity") + ylab("Mean effect on output")+
      coord_flip() + geom_errorbar(aes(ymax = lower, ymin=upper), width=0.25,data=limits) +
      theme_bw() + ggtitle(paste("Individual influence of", variable, "(with", ncol(p), "profiles)")) +
      scale_fill_gradient2("Count", low = "red", high = "green", midpoint=0)
    suppressWarnings(print(m))
    
  } else {
    p <- as.data.frame(profiles[variable])
    plot(y=data[,ncol(data)], x=var, xlab=variable,
         ylab=names(data)[ncol(data)], pch=19,
         main=paste("Individual influence of", variable, "(with", ncol(p), "profiles)"),
         ylim=c(min(p), max(p)), col=0, xlim=xlim0)
    rug(var,ticksize = 0.06,  lwd = 0.8)
    #sapply(p, function(x) lines(y=x, x=row.names(p), lwd=0.5, col=rgb(0,0,0,0.2)))
    lines(x=row.names(p), y=rowMeans(p), type="l", col="green", lwd=4)
  }
}