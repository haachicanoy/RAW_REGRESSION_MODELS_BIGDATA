

multilayerPerceptronFun <- function(variety, dirLocation=paste0(getwd(),"/"), nb.it=100, HQplots=FALSE,
                                    ylabs="Yield (Kg/HA)", pertuRelevance=TRUE, ncores=20,
                                    saveWS=FALSE, uncorrset=TRUE)
{
  
  ngw <- nchar(dirLocation)
  if( substring(dirLocation,ngw-16,ngw)=="VARIETY_ANALYSIS/" ){}else{return(cat("Aun no se encuentra en la carpeta VARIETY_ANALYSIS\nUtilize la funcion setwd para dirigirse a este carpeta"))}
  
  library(snowfall)
  require(caret)
  require(nnet)
  
  
  #CARET PART 1
  sfInit(parallel=T,cpus=ncores)
  sfLibrary(caret)
  sfLibrary(nnet)
  
  if(uncorrset){
    dirDataSet <- paste0(dirLocation,variety,"/DATA_SETS/",variety,"_reduced.csv")
  }else if(uncorrset==FALSE){
    dirDataSet <- paste0(dirLocation,variety,"/DATA_SETS/",variety,"_complet.csv")
  }else{return(print("ERROR: uncorrset should be a logical value"))}
  
  #  dirDataSet <- paste0(dirLocation,variety,"/DATA_SETS/andresMatrizDePrueba.csv")#----------Activar
  
  dirSave    <- paste0(dirLocation,variety,"/ARTIFICIAL_NEURAL_NETWORK/")  
  
  dataSets   <- lapply(dirDataSet,function(x){read.csv(x,row.names=1)})
  
  
  #FUNCION TRAIN EN CARET
  
  mlpModel <- function(x)
  {
    
    #CREATING PARTITION
    
    output <- ncol(normMat)
    
    inTrain  <- createDataPartition(y=normMat[,output], p=0.7, list=F)
    training <- normMat[inTrain,]
    testing  <- normMat[-inTrain,]
    
    ctrl <- expand.grid(size=1:15, decay=(1:10/100) )
    
    
    model <- train(training[,-ncol(training)] ,training[,output], method="nnet",
                   tuneGrid=ctrl, trControl=trainControl(method="repeatedcv", number=5), maxit=1000, linOut=TRUE)
    rmseVals <- RMSE(predict(model, testing[,-ncol(training)]), testing[,output])      
    
    rsquare <- R2(predict(model, testing[,-ncol(training)]), testing[,output]) * 100
    
    return(list(model,rmseVals,rsquare,training,testing))
    
  } 
  
  #CARET LOAD FUNCTION PART 2
  
  for(i in 1:length(variety))
  { 
    cat(paste("Variety",variety[i]),"\n")
    data <- dataSets[[i]]
    
    dimData <- dim(data)
    
    mind   <- as.vector(apply(data,2,min))
    rangd  <- as.vector(apply(data,2,max)-mind)
    
    
    normMat <- (data-matrix(mind,nrow = dimData[1],ncol= dimData[2],byrow=T))*1/matrix(rangd,nrow = dimData[1],ncol= dimData[2],byrow=T)-0
    
    rmseVals <- 0
    
    #EJECUCION EN PARALELO DE LOS MODELOS EN CARET
    
    sfExport("normMat")
    
    cat("Starting model process: ",paste(variety[i]),"\n")
    Sys.time()->start
    allModelsAndRMSE <- sfLapply(1:nb.it,mlpModel)
    cat("Finished process: ",paste(variety[i]),"\n")
    print(Sys.time()-start)
    sfRemove("normMat")  
    
    allRMSE   <- unlist(lapply(allModelsAndRMSE,function(x){x[[2]]}))    
    bestModels <- order(allRMSE)[1:nb.it]
    
    allModels <- lapply(allModelsAndRMSE,function(x){x[[1]]})[bestModels]
    allR2     <- unlist(lapply(allModelsAndRMSE,function(x){x[[3]]}))[bestModels]
    allTraining <- lapply(allModelsAndRMSE,function(x){x[[4]]})[bestModels]
    allTesting  <- lapply(allModelsAndRMSE,function(x){x[[5]]})[bestModels]
    
    topModel    <- allModels
    
    cat(paste0("Computing Metrics ",variety[i]),"\n")
    
    #t0 <- proc.time()
    
    #proc.time() - t0
    
    #-------------------------------------------------------------------------------------------------
    
    if(pertuRelevance==T)
    {
      
      Sys.time()->start
      pertuImport       <- lapply(topModel,function(x){varImportance(x)})  
      print(Sys.time()-start)
      
      currentVarImp <-   do.call(cbind,pertuImport)
      
      pertuImportMedian <- sort(apply(do.call(cbind,pertuImport),1,mean),decreasing = F)
      
      #SCALING
      
      scale <- allR2 / as.numeric(apply(currentVarImp,2,sum))
      
      scaledVarImp <-  t(t(currentVarImp) * scale)
      
      v <- as.data.frame(scaledVarImp) 
      
      write.csv(v,paste0(dirSave[i],"weighMatrix.csv"))
      
      ordered <- sort(apply(v,1, median), decreasing=F)
      
      perf1 <- signif(mean(allR2), 5)
      
      se <- apply(v, 1, function(x){ 1.96*sd(x, na.rm=TRUE)/sqrt(ncol(v))})
      se <- data.frame(se,names(se))
      names(se) <- c("se","Variable")
      
      mean <- as.data.frame(ordered)
      mean <- cbind(mean, names(ordered))
      names(mean) <- c("Mean", "Variable")
      mean$Variable <- factor(names(ordered), levels= names(ordered))
      
      stadistc <- merge(se,mean,by.x="Variable",by.y="Variable")
      
      stadistc <- stadistc[order(stadistc$Mean,decreasing=F),]
      
      mean$se <- array(0.3,nrow(mean))
      
      errBars <- transform(stadistc, lower=Mean-se,upper=Mean+se )
      
      if(HQplots==FALSE){
        
        png(paste0(dirSave[i], variety[i], "_InputRelvancePerturbatuion.png"), wid=800, hei=500, pointsize=20)
        
        m <- ggplot(mean, aes(x=Variable, y=Mean))
        m <- m + geom_bar(stat="identity", width=0.5, fill="blue") + ylab("Mean importance")+
          geom_errorbar(aes(ymax = lower, ymin=upper), width=0.25,data=errBars) + coord_flip() +
          theme_bw() + 
          ggtitle(paste("Importance of variables (with a mean R2 of", perf1, "%)")) +
          theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))
        
        suppressWarnings(print(m))
        dev.off()
        
      } else {
        
        m <- ggplot(mean, aes(x=Variable, y=Mean))
        m <- m + geom_bar(stat="identity", width=0.5, fill="blue") + ylab("Mean importance")+
          geom_errorbar(aes(ymax = lower, ymin=upper), width=0.25,data=errBars) + coord_flip() +
          theme_bw() + 
          ggtitle(paste("Importance of variables (with a mean R2 of", perf1, "%)")) +
          theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0))
        
        wid = 6.67; hei = 10.67
        ggsave(filename=paste0(dirSave[i], variety[i], "_InputRelvancePerturbatuion.png"), plot=m, width=wid, height=hei, units='in')
        system(paste("convert -verbose -density 300 ", dirSave[j], variety[i], "_InputRelvancePerturbatuion.pdf -quality 100 -sharpen 0x1.0 -alpha off ", dirSave[j], variety[i], "_InputRelvancePerturbatuion.png", sep=""), wait=TRUE)
        
      }
      
    }else{}
    
    sink(paste0(dirSave[i] ,variety[i],"_RSQUARE.txt"))
    cat("R SQUARTE MODELS: \n\n")
    print(summary(allR2))
    cat("\n RMSE SQUARTE MODELS: \n")
    print(summary(allRMSE))
    sink()
    
    #EXTRACCION DE PERFILES
    
    namPredic   <- names(normMat)[-ncol(normMat)]
    dataSetStep <- normMat[namPredic]
    nColums     <- ncol(dataSetStep)
    
    mindss <- apply(dataSetStep,2,min)
    maxdss <- apply(dataSetStep,2,max)
    
    matMin <- matrix(mindss,100,nColums,byrow = T)
    matQ1  <- matrix(apply(dataSetStep,2,function(x){quantile(x,0.25)}),100,nColums,byrow = T)
    matMed <- matrix(apply(dataSetStep,2,median),100,nColums,byrow = T)
    matQ3  <- matrix(apply(dataSetStep,2,function(x){quantile(x,0.75)}),100,nColums,byrow = T)
    
    for(j in 1:nColums)
    {
      matMax     <- matrix(maxdss,100,nColums,byrow = T)
      
      xProf      <- seq(mindss[j],maxdss[j],length.out = 100)  
      xProfDesn  <- (xProf-0)*(rangd[j])/1+mind[j]
      listMats <- list(matMin,matQ1,matMed,matQ3,matMax)
      listFitted <- apply(sapply(listMats,function(x){z <- as.data.frame(x);z[,j] <- xProf ;colnames(z) <- namPredic;y <- predict(topModel[[1]],z);return(y)}),1,median)
      fitteDesn <-(listFitted+0)*(rangd[ncol(normMat)])/1+mind[ncol(normMat)]
      
      png(paste(dirSave[i], "PROFILES/pefil_", namPredic[j], ".png", sep=""), width=700, height=350)
      
      layout(matrix(c(1, 2), ncol=2, nrow=1))
      plot(data[,j], data[,ncol(data)], pch=21, cex=0.8, bg="azure3", col="azure3", ylab=ylabs, ylim=c(min(data[,ncol(data)]), max(data[,ncol(data)])), xlab=namPredic[j], main=namPredic[j])
      points(xProfDesn, fitteDesn, bg="blue", col="blue", pch=21)
      plot(xProfDesn, fitteDesn, type="l", col=0, ylab=ylabs, xlab=namPredic[j], main="Profile Zoom")
      lines(supsmu(xProfDesn, fitteDesn), lwd=2, col="green")
      
      dev.off()
      
    }
    
    rm(allModelsAndRMSE)
    
    if(HQplots==TRUE){
      # Remove pdf files
      setwd(dirSave[j])
      pdfFiles <- list.files(path=getwd(), pattern='*.pdf$', full.names=TRUE)
      file.remove(pdfFiles) # system(paste('find . ! -name "*.png" ! -name "weighMatrix.csv" -type f -delete',sep=''))
    }
    
  }
  
  sfStop()
  
}
