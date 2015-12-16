
multiProfile <- function(data, profiles, variable, xlim0=NULL, xla=variable, yla=names(data)[ncol(data)]) {
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
    m <- m + geom_bar(stat="identity") + ylab("Mean effect on output")
    m <- m + coord_flip() + geom_errorbar(aes(ymax = lower, ymin=upper), width=0.25, data=limits)
    m <- m + theme_bw() + ggtitle(paste("Individual influence of", variable))
    m <- m + theme(axis.text.x = element_text(size=7),
                   axis.text.y = element_text(size=7),
                   axis.title.x = element_text(face="bold",size=7),
                   axis.title.y = element_text(face="bold",size=7),
                   legend.text = element_text(size=6),
                   legend.title = element_text(face="bold",size=7),
                   plot.title = element_text(face="bold", size=10))
    m <- m + scale_fill_gradient2("Count", low="red", high="green", midpoint=0)
    suppressWarnings(print(m))
    
  } else {
    p <- as.data.frame(profiles[variable])
    par(mar=c(5, 4, 2, 2) + 0.1)
    plot(x=var, y=data[,ncol(data)],
         xlab=xla, ylab=yla, pch=19,
         main='', # paste("Individual influence of", variable),
         ylim=c(min(p), max(p)), col=0, xlim=xlim0)
    rug(var, ticksize=0.06,  lwd=0.8, col=rgb(red=0, green=0, blue=0, alpha=0.4))
    sapply(p, function(x) lines(y=x, x=row.names(p), lwd=0.5, col=rgb(red=0, green=0, blue=0, alpha=0.3)))
    lines(x=row.names(p), y=rowMeans(p), type="l", col="green", lwd=4)
  }
}