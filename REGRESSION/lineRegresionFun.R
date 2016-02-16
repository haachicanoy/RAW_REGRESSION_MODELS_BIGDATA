


lineaRegresionFun <- function(variety,dirLocation=paste0(getwd(),"/"),ylabs="Yield (Kg/HA)" )
{
  ngw <- nchar(dirLocation)
  if( substring(dirLocation,ngw-16,ngw)=="VARIETY_ANALYSIS/" ){}else{return(cat("Aun no se encuentra en la carpeta VARIETY_ANALYSIS\nUtilize la funcion setwd para dirigirse a este carpeta"))}
  
  library(relaimpo)
  require(caret)
  
  dirDataSet <- paste0(dirLocation,variety,"/DATA_SETS/",variety,"_reduced.csv")
  dirDataSetCom <- paste0(dirLocation,variety,"/DATA_SETS/",variety,"_complet.csv")
  
  dirSave    <- paste0(dirLocation,variety,"/LINEAR_REGRESSION/")  
  
  dataSets   <- lapply(dirDataSet,function(x){read.csv(x,row.names=1)})
  dataSets2   <- lapply(dirDataSetCom,function(x){read.csv(x,row.names=1)})
  
  
  for(i in 1:length(variety))
  { 
    
    data <- dataSets[[i]]
    dataComp <- dataSets2[[i]]
    
    if(sum(
      unlist(
        lapply(
          1:ncol(dataComp),function(x){is.factor(dataComp[,x])}
        )
      )
    )>0){return(
      cat("lineRegresionFun solo funciona cuando no hay variables cualitativas en los datos")
    )}else{}
    
    
    formLm <- formula(paste(names(data)[ncol(data)], "~ ."))
    
    mode <- lm(formLm, data = data)
    
    sink(paste0(dirSave[i],variety[i],"_summaryLm.txt"))
    
    print(summary(mode))
    
    sink()
    
    crli <- calc.relimp(mode)
    
    relev <-  round(sort(crli$lmg,decreasing = T)*100,2)
    
    library(ggplot2)
    
    #leer las etiquetas y mantener el orden
    namVar = factor(names (relev),levels=names (relev))
    
    relvar = as.numeric(relev) #lee los valores
    
    #armar el ser de datos para graficar
    metriData = data.frame(namVar,relvar)
    
    tiff(paste0(dirSave[i],variety[i],"_InputRelvance.tiff"), wid=800,hei=600, pointsize = 20)
    
    k <-  ggplot(metriData,aes(namVar, relvar))+geom_bar( stat="identity",fill="royalblue1",width=0.6,
                                                          color="black")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(axis.title.y=element_text(
                                                            size=15),axis.text= element_text(colour = "gray26",size = 13))+xlab("")+ylab("% variance explained")+theme(panel.background =element_rect(colour = "gray38")) +theme(panel.background = element_rect(colour = "gray32"))
    kf <- k + theme(panel.background = element_rect(fill = 0,colour = "gray"))+ geom_hline(h=0)
    
    kf <- kf + ggtitle(paste("Metricas"))
    
    print(kf) #ejecutar el gracico en Plots y en la ventana
    
    dev.off()
    
    
    ###CON BARCKWARD
    
    mode1 <- step(lm(formLm , data = data),direction = "backward")
    
    sink(paste0(dirSave[i] ,variety[i],"_summaryLmBackward.txt"))
    
    
    print(summary(mode1))
    
    sink()
    
    crli <- calc.relimp(mode1)
    
    relev <-  round(sort(crli$lmg,decreasing = T)*100,2)
    
    
    #leer las etiquetas y mantener el orden
    namVar = factor(names (relev),levels=names (relev))
    
    relvar = as.numeric(relev) #lee los valores
    
    #armar el ser de datos para graficar
    metriData = data.frame(namVar,relvar)
    
    tiff(paste0(dirSave[i],variety[i],"_InputRelvanceBackward.tiff"), wid=800,hei=600, pointsize = 20)
    
    k <-  ggplot(metriData,aes(namVar, relvar))+geom_bar( stat="identity",fill="royalblue1",width=0.6,
                                                          color="black")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(axis.title.y=element_text(
                                                            size=15),axis.text= element_text(colour = "gray26",size = 13))+xlab("")+ylab("% variance explained")+theme(panel.background =element_rect(colour = "gray38")) +theme(panel.background = element_rect(colour = "gray32"))
    kf <- k + theme(panel.background = element_rect(fill = 0,colour = "gray"))+ geom_hline(h=0)
    
    kf <- kf + ggtitle(paste("Metricas"))
    
    print(kf) #ejecutar el gracico en Plots y en la ventana
    
    dev.off()
    
    #EXTRACCION DE PERFILES
    
    namPredic <- names(mode1$coefficients)[-1]
    dataSetStep <- data[namPredic]
    nColums <- ncol(dataSetStep)
    
    mindss <- apply(dataSetStep,2,min)
    maxdss <- apply(dataSetStep,2,max)
    
    
    matMin <- matrix(mindss,100,nColums,byrow = T)
    matQ1  <- matrix(apply(dataSetStep,2,function(x){quantile(x,0.25)}),100,nColums,byrow = T)
    matMed <- matrix(apply(dataSetStep,2,median),100,nColums,byrow = T)
    matQ3  <- matrix(apply(dataSetStep,2,function(x){quantile(x,0.75)}),100,nColums,byrow = T)
    
    
    for(j in 1:nColums)
    {  
      matMax <- matrix(maxdss,100,nColums,byrow = T)
      
      xProf <- seq(mindss[j],maxdss[j],length.out = 100)  
      listMats <- list(matMin,matQ1,matMed,matQ3,matMax)
      listFitted <- apply(sapply(listMats,function(x){z <- as.data.frame(x);z[,j] <- xProf ;colnames(z) <- namPredic;y <- predict.lm(mode1,z);return(y)}),1,median)
      
      tiff(paste(dirSave[i],"pefil_",namPredic[j],".tiff",sep=""),width = 700, height = 350)
      layout(matrix(c(1,2),ncol=2,nrow=1))
      plot(xProf,listFitted,pch=21,cex=0.8,bg="blue",col="blue",ylab=ylabs,ylim=c(min(data[,ncol(data)]),max(data[,ncol(data)])),xlab=namPredic[j],main=namPredic[j])
      points(dataSetStep[,j],data[,ncol(data)],bg="azure3",pch=21)
      plot(xProf,listFitted,type="l",col=0,ylab=ylabs,xlab=namPredic[j],main="Profile Zoom")
      lines(supsmu(xProf,listFitted),lwd=2,col="green")
      dev.off()
    }
  }
  
}
