
multiProfile <- function(data, profiles, variable,xlim0=NULL,pp.szmain=15,pp.sztxtx=15,
                pp.sztxty=15,pp.szlbty=20,pp.szlbtx=20,pp.lgndtx=15) {
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
      scale_fill_gradient2("Count", low = "red", high = "green", midpoint=0)+
      theme(plot.title = element_text(vjust=3,size=pp.szmain),
            axis.text.x = element_text(size = pp.sztxtx),
            axis.text.y = element_text(size = pp.sztxty),
            axis.title.x = element_text(size = pp.szlbty),
            axis.title.y = element_text(size = pp.szlbtx),
            legend.text=element_text(size=pp.lgndtx))
    
    suppressWarnings(print(m))
    
  } else {
    p <- as.data.frame(profiles[variable])
    
    datf0 <- data.frame(x=var, y=data[,ncol(data)])
    datf  <- data.frame(x=as.numeric(row.names(p)), y=as.numeric(rowMeans(p)))
    
    partialDep <-  ggplot(datf0,aes(x=x,y=y))+geom_line(col="green",lwd=1.5,aes(x=x,y=y),data=datf)+
        ylim(c(min(p), max(p)))+theme_bw()+ggtitle(paste("Individual influence of",
                                         variable, "(with", ncol(p), "profiles)"))+
         geom_rug(sides="b")+
            theme(plot.title = element_text(vjust=3,size=pp.szmain),
                axis.text.x = element_text(size = pp.sztxtx),
                axis.text.y = element_text(angle = 90,hjust = 0.5,size = pp.sztxty),
                axis.title.x = element_text(size = pp.szlbty),
                axis.title.y = element_text(size = pp.szlbtx))+
                ylab(names(data)[ncol(data)])+xlab(variable)
        
    
      #sapply(p, function(x) lines(y=x, x=row.names(p), lwd=0.5, col=rgb(0,0,0,0.2)))
  print(partialDep)
  }
}