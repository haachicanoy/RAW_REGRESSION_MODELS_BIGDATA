#--------------------------------------------------------------
#--------------------------FENNIX FUNCTION---------------------

completeSummary <-function(table)
{
  
  ncTable <- ncol(table)
  
  factCol <- unlist(lapply(1:ncTable,function(x){is.factor(table[,x])}))
  
  if(sum(factCol)>0)
  {
    
    factCol <- which(factCol)
    
    cuantVar <- table[,-factCol]
    
    cualVar  <- as.data.frame(table[,factCol])
    
    if(ncol(cualVar)==1){names(cualVar) <- names(table)[factCol];summCual <- list(table(cualVar));names(summCual) <- names(table)[factCol]}
    else{summCual <- apply(cualVar,2,table)}
    
    
  }else{cuantVar  <- table}  
  
  min <- apply(cuantVar,2,min)
  Q1  <- apply(cuantVar,2,function(x){quantile(x,0.25)})
  mad <- apply(cuantVar,2,median)
  med <- apply(cuantVar,2,mean)
  Q3  <- apply(cuantVar,2,function(x){quantile(x,0.75)})
  max <- apply(cuantVar,2,max)
  sd  <- apply(cuantVar,2,sd)
  cv <-  apply(cuantVar,2,function(x){if(mean(x)!=0){sd(x)/mean(x)*100}else{cv=NA}})
  #if(mean(x!=0))
  #{  
  # cv  <- apply(table,2,function(x){sd(x)/mean(x)*100})
  #}else{cv=NA}
  
  summCuant <- as.data.frame(round(as.matrix( data.frame(max,min,Q1,Q3,mad,med,sd,cv)),3) )
  
  var <- row.names(summCuant)
  
  summCuant <- data.frame(var,summCuant)
  
  if(sum(factCol)>0){summaryDataSet <- list(summCuant,summCual)}else{summaryDataSet <- summCuant}
  
  
  
  return(summaryDataSet)
}


corScheme <- function(compMatrix,unCorMatrix,dirSave)
{
  namUnMatrix <- names(unCorMatrix)
  namMatrix   <- names(compMatrix)
  
  #sink(paste0(dirSave,"corScheme.txt"))
  fileNam     <- paste0(dirSave,"corScheme2.txt")
  fileCon  <- file(fileNam)
  writeLines("correlation Matrix",fileCon)
  close(fileCon)
  
  for(i in 1:(ncol(unCorMatrix)-1))
  {  
    target <- unCorMatrix[namUnMatrix[i]][,1]
    inputs <- compMatrix[-which(namUnMatrix[i] == namMatrix)]
    
    
    corOth <- cor(x=inputs,y=target)[,1]
    namVar <- names(corOth)
    
    corColum <- factor("")
    
    if(sum(corOth >0.69)>0 & sum(corOth < -0.69)>0)
    {
      nvN <- namVar[corOth < -0.69]
      nvP <- namVar[corOth >0.69]
      
      corColun <- c(nvN,paste0(namUnMatrix[i],"*"),nvP)
      
    }else if(sum(corOth < -0.69)>0)
    {
      nvN <- namVar[corOth < -0.69]
      corColun <- c(nvN,paste0(namUnMatrix[i],"*"))
    }else if(sum(corOth > 0.69)>0){   
      nvP <- namVar[corOth >0.69]
      corColun <- c(paste0(namUnMatrix[i],"*"),nvP)
    }else{  corColun <- c(paste0(namUnMatrix[i],"*")) }
    cat(paste(corColun,collapse='\t'), file=fileNam, append = T,   sep = "\n")
  }
  #sink()
}


unCorrMatrix <- function(data,dirSav,cor.reduce="caret")
{ 
  
  if(sum(round(apply(data,2,sd),2)==0)>0)
  {  
    data <- data[,-which(round(apply(data,2,sd),2)==0)]
  }else{data}
  corMat <- abs(cor(data[,-ncol(data)]))
  
  
  if(cor.reduce=="clasic")
  { 
    sink(paste0(dirSav,"saveTest.txt")) 
    while(sum(corMat>=0.7 & corMat<1)!=0 )
    {
      i=1
      while(ncol(data)>i)
      { 
        
        print(paste("Step",i))
        outPut    <- data[,ncol(data)]
        colCor    <- corMat[,i]
        
        #ABRE
        nam <- names(data)
        #CIERRA      
        
        if(sum(colCor>=0.69 & colCor<1)==0)
        {   
          data <- data
        }else{
          corCol    <- c(i,which(colCor>=0.69 & colCor<1))
          #ABRE
          print("Correlated variables:")
          print(data.frame(nam[corCol],colCor[corCol]))
          #CIERRA
          corOutput <- as.vector(abs(cor(data[,corCol],outPut)))
          corMax    <- corCol[which(corOutput==max(corOutput))[1]]
          
          if(i == corMax)
          {
            quitCol <-  corCol[-1]
          }else{
            quitCol <-  corCol[1]
          }
          #quitCol   <- corCol[which(corOutput!=max(corOutput))] #HE QUITADO [1]
          
          data      <- data[,-quitCol]
          #ABRE
          print("variables removed:")
          
          print(nam[quitCol])
          print("---------------*--------------")
          #CIERRA      
        }
        corMat <- abs(cor(data))    
        i=i+1
        #print(c(i,ncol(data)))
      }
    }
    sink()
  }else if(cor.reduce=="caret"){
    require(caret)
    
    highlyCorDescr <- findCorrelation(corMat, cutoff = 0.69)
    if(length(highlyCorDescr)>0)
    {
      namHCD <- colnames(data)[highlyCorDescr]
      data  <- data[, -highlyCorDescr]
      cat("\n\nBY CORRELATION:\n\n",file=paste0(dirSav,"RemovedVariables.txt"),append=T)
      cat(paste(namHCD,collapse="\n"),file=paste0(dirSav,"RemovedVariables.txt"),append=T)
    }else{data  <- data}
    
  }else{}
  return(data)
}


createFolders <- function(dirFol,variety)
{
  if(!file.exists(paste0(dirFol,"VARIETY_ANALYSIS"))){dir.create("VARIETY_ANALYSIS");setwd("VARIETY_ANALYSIS")}else if( substring(getwd(),nchar(getwd())-15,nchar(getwd()))=="VARIETY_ANALYSIS" ){}else{setwd("VARIETY_ANALYSIS")}
  
  
  for(i in 1:length(variety))
  {  
    if(!file.exists(as.character(variety[i]))){dir.create(as.character(variety[i]))}else{}
    if(!file.exists(paste0(variety[i],"/LINEAR_REGRESSION"))){dir.create(paste0(variety[i],"/LINEAR_REGRESSION"))}else{}
    if(!file.exists(paste0(variety[i],"/DESCRIPTIVE_ANALYSIS"))){dir.create(paste0(variety[i],"/DESCRIPTIVE_ANALYSIS"))}else{}
    if(!file.exists(paste0(variety[i],"/ARTIFICIAL_NEURAL_NETWORK"))){dir.create(paste0(variety[i],"/ARTIFICIAL_NEURAL_NETWORK"))}else{}
    if(!file.exists(paste0(variety[i],"/RANDOM_FOREST"))){dir.create(paste0(variety[i],"/RANDOM_FOREST"))}else{}
    if(!file.exists(paste0(variety[i],"/STOC_GRAD_BOOS"))){dir.create(paste0(variety[i],"/STOC_GRAD_BOOS"))}else{}
    if(!file.exists(paste0(variety[i],"/C_FOREST"))){dir.create(paste0(variety[i],"/C_FOREST"))}else{}
    if(!file.exists(paste0(variety[i],"/DATA_SETS"))){dir.create(paste0(variety[i],"/DATA_SETS"))}else{}
    
    dirFennixVar <- paste0(variety[i],"/ARTIFICIAL_NEURAL_NETWORK/")
    
    if(!file.exists(paste0(dirFennixVar,"/PROFILES"))){dir.create(paste0(dirFennixVar,"/PROFILES"))}else{}
  }
  
}



descriptiveGraphics <- function(variety, dataSet, inputs, segme, output,
                                ylabel="Rendimiento (kg/ha)", smooth=FALSE,
                                smoothInd=NULL, ghrp="box", res=NA, sztxt=15, szlbls=15,
                                colbox="skyblue", colpts="greenyellow", colsmot="red",
                                szpts=4, szdts=1.5)
{
  namsDataSet <- names(dataSet)
  
  if(variety!="All")
    {
    dataSetV    <- subset(dataSet,dataSet[,segme]==variety)[,-segme]
    } else {
      dataSetV <- dataSet[,-segme]
    }
  
  ncTable <- ncol(dataSetV)
  outputVar <- dataSetV[,(output-1)]
  factCol0 <- unlist(lapply(1:ncTable,function(x){is.factor(dataSetV[,x])}))
  
  if(sum(factCol0)>0)
  {
    
    factCol <- which(factCol0)
    cuantVar <- as.data.frame(dataSetV[,-factCol])
    cualVar  <- as.data.frame(dataSetV[,factCol])
    
    if(ncol(cualVar)==1){
      names(cualVar) <- names(dataSetV)[factCol]
      }
  }else {
    cuantVar  <- dataSetV
  }
  
  nInputs <- ncol(cuantVar)
  namCuanInputs <- names(cuantVar)[-nInputs]
  
  tme <- theme_bw() + theme(legend.position="none")
  tme <- tme + theme(axis.title.x = element_text(size=szlbls),
                     axis.title.y = element_text(size=szlbls),
                     axis.text.x  = element_text(size=sztxt),
                     axis.text.y  = element_text(size=sztxt, hjust=0.5)) # angle=90
  
  for(namC in namCuanInputs)
  {
    
    plotData <- data.frame(x=cuantVar[,namC], y=outputVar)
    
    plo <- ggplot()+ geom_point(aes(x=x, y=y), data=plotData, colour="grey4", fill=colpts, size=szpts, shape=21) + ylab(ylabel) + xlab(namC) + tme # , alpha=0.2
    box <- qplot(y=x ,x=namC, alpha=I(0.00001), data=plotData) + geom_boxplot(width=0.3, fill=colbox) + xlab("") + ylab(namC) + tme
    
    if(!isTRUE(smooth))
    {
      
      # png(paste0(variety, "/DESCRIPTIVE_ANALYSIS/", variety, "_", namC, ".png"), width=775, height=280, res=res)
      # grid.arrange(plo, box, ncol=2, top=textGrob(paste0(variety, " - ", namC), gp=gpar(cex=1.1, fontface='bold')))
      g = arrangeGrob(plo, box, ncol=2, top=textGrob(paste0(variety, " - ", namC), gp=gpar(cex=1.1, fontface='bold')))
      g2 = gTree(children=gList(g), cl=c("arrange", "ggplot"))
      ggsave(filename=paste0(variety, "/DESCRIPTIVE_ANALYSIS/", variety, "_", namC, ".pdf"), plot=g2, width=10.33, height=5.33, units='in', dpi=300)
      # dev.off()
      
    } else{
      
      # png(paste0(variety, "/DESCRIPTIVE_ANALYSIS/", variety, "_", namC, ".png"), width=775, height=280, res=res)
      
      lo <- suppressWarnings(loess(outputVar~cuantVar[,namC]))
      xl <- seq(min(cuantVar[,namC]), max(cuantVar[,namC]), (max(cuantVar[,namC]) - min(cuantVar[,namC]))/1000)
      predicLoess <- try(predict(lo, xl), silent=T)
      
      if(is(predicLoess)[1]!="try-error")
      {
        
        loesMat <- data.frame(xl, predicLoess)
        plo <- plo + geom_line(aes(x=xl, y=predicLoess), lwd=1.3, col=colsmot, data=loesMat)
        
      } else{
        
        lows <- lowess(outputVar~cuantVar[,namC])
        lowsMat <- data.frame(x=lows$x, y=lows$y)
        plo <- plo + geom_line(aes(x=x, y=y), lwd=1.3, col=colsmot, data=lowsMat)
        
      }
      
    }
    
    if(ghrp=="box")
    {
      # grid.arrange(plo, box, ncol=2, top=textGrob(paste0(variety, " - ", namC), gp=gpar(cex=1.1, fontface='bold')))
      g = arrangeGrob(plo, box, ncol=2, top=textGrob(paste0(variety, " - ", namC), gp=gpar(cex=1.1, fontface='bold')))
      g2 = gTree(children=gList(g), cl=c("arrange", "ggplot"))
      
    }else if(ghrp=="varPoints"){
      
      datf  <- data.frame(x=dataSet[,namC], y=dataSet[,output])
      datf1 <- data.frame(x=cuantVar[,namC], y=outputVar)
      
      plo1 <- ggplot() + geom_point(aes(x=x, y=y), colour="grey4", fill="grey", data=datf, size=szpts, shape=21) # ,alpha =0.6
      plo1 <- plo1 + geom_point(aes(x=x, y=y), data=datf1, colour="grey4", fill="red", size=szpts, shape=21) # ,alpha =0.6
      plo1 <- plo1 + theme_bw()+ ylab(ylabel) + xlab(namC) + tme
      
      # grid.arrange(plo1, box, ncol=2, top=textGrob(paste0(variety, " - ", namC), gp=gpar(cex=1.1, fontface='bold')))
      g = arrangeGrob(plo1, box, ncol=2, top=textGrob(paste0(variety, " - ", namC), gp=gpar(cex=1.1, fontface='bold')))
      g2 = gTree(children=gList(g), cl=c("arrange", "ggplot"))
      
    } else{
      stop("ghrp Invlid")
      }
    # dev.off()
    ggsave(filename=paste0(variety, "/DESCRIPTIVE_ANALYSIS/", variety, "_", namC, ".pdf"), plot=g2, width=10.33, height=5.33, units='in', dpi=300)
  }
  
  # png(paste0(variety, "/DESCRIPTIVE_ANALYSIS/", variety, "_Yield.png"), width=775, height=280, res=res)
  otpvar <- data.frame(x=outputVar)
  
  ydplot <- ggplot()+geom_histogram(aes(x=x), data=otpvar, breaks=seq(min(outputVar), max(outputVar), length.out=nclass.Sturges(outputVar)), colour="gray1", fill=colbox)
  ydplot <- ydplot + tme + xlab(ylabel) + ylab("Frecuencia")
  
  boxdf <- data.frame(x=ylabel,y=outputVar)
  box <- qplot(y=y, x=x, alpha=I(0.00001), data=boxdf) + geom_boxplot(width=0.3, fill=colbox) + xlab("") + ylab(ylabel) + tme
  
  # grid.arrange(ydplot, box, ncol=2, top=textGrob(paste(variety, "_Yield"), gp=gpar(cex=1.1, fontface='bold')))
  g = arrangeGrob(ydplot, box, ncol=2, top=textGrob(paste(variety, "_Yield"), gp=gpar(cex=1.1, fontface='bold')))
  g2 = gTree(children=gList(g), cl=c("arrange", "ggplot"))
  
  # dev.off()
  
  summ <- completeSummary(dataSetV)
  
  if(sum(factCol0)>0)
  { 
    namCualVar <- names(cualVar)
    
    write.table(summ[[1]], paste0(variety, "/DESCRIPTIVE_ANALYSIS/", variety, "summaryVariables.csv"), append=F, sep=",", row.names=F)
    
    for(namC in namCualVar)
    {
      dp <- data.frame(x=cualVar[,namC],y=outputVar)
      dotplot <- ggplot() + geom_dotplot(aes(x=x, y=y), data=dp, binaxis="y", stackdir="center", colour="grey4", fill=colpts, dotsize=szdts, na.rm=T, binwidth=diff(range(dp$y)/30))
      dotplot < dotplot + ylab(ylabel) + xlab(namC) + tme + theme(axis.text.x = element_text(angle=45, hjust=1))
      
      boxplot <- ggplot() + geom_boxplot(width=0.3, aes(y=y, x=x), data=dp, fill=colbox) + tme + xlab(namC) + ylab(ylabel) + theme(axis.text.x = element_text(angle=45, hjust=1))
      
      # png(paste0(variety, "/DESCRIPTIVE_ANALYSIS/", variety, "_", namC, ".png"), width=775, height=280, res=res)
      # grid.arrange(dotplot, boxplot, ncol=2, top=textGrob(paste(variety, "_", namC), gp=gpar(cex=1.1, fontface="bold")))
      g = arrangeGrob(dotplot, boxplot, ncol=2, top=textGrob(paste(variety, "_", namC), gp=gpar(cex=1.1, fontface="bold")))
      g2 = gTree(children=gList(g), cl=c("arrange", "ggplot"))
      # dev.off()
      
      cat("\n", namC, "\n", file=paste0(variety, "/DESCRIPTIVE_ANALYSIS/", variety, "summaryVariables.csv"), append=T)
      write.table(summ[[2]][[namC]],paste0(variety,"/DESCRIPTIVE_ANALYSIS/",variety,"summaryVariables.csv"),append=T, row.names=F, col.names=F, sep=",")
      
    }
    
  } else{
    write.table(summ, paste0(variety, "/DESCRIPTIVE_ANALYSIS/", variety, "summaryVariables.csv"), append=F, sep=",", row.names=F)
    }
  cat("The exploratory analysis is available in:", paste0("\t\t", "VARIETY_ANALYSIS/", variety, "/DESCRIPTIVE_ANALYSIS/"))
}



dataSetProces <- function(variety,dataSet,segme,corRed)
{
  
  nInputList <- 0
  
  for(i in 1:length(variety))
  { 
    dirSav <- paste0(variety[i],"/DATA_SETS/")
    cat(paste("List of variables removed:\n"),file=paste0(dirSav,"RemovedVariables.txt"),append=F)
    if(variety[i]=="All")
    {
      listVari <- dataSet[,-segme]
    }else{listVari   <- subset(dataSet[,-segme],as.character(dataSet[,segme])==variety[i])}
    
    nvz <- nearZeroVar(listVari)
    
    if(length(nvz)>0){
      
      allData <- listVari[,-nvz]
      
      
      
      cat(paste("\nVARIABLES WITH VARIANCE NEAR TO ZERO:\n\n"),file=paste0(dirSav,"RemovedVariables.txt"),append=T)
      cat(paste(names(allData)[nvz],collapse = "\n"),file=paste0(dirSav,"RemovedVariables.txt"),append=T)
    }else{allData <- listVari}
    
    #     if(sum(round(apply(listVari,2,sd),2)==0)>0)
    #     {  
    #       allData <- listVari[,-which(round(apply(listVari,2,sd),2)==0)]
    #     }else{allData <- listVari }
    
    listVari2 <- listVari
    
    listVari2$ID <- row.names(listVari)
    
    listVari2 <- listVari2[,c(length(listVari2),1:(length(listVari2)-1))]
    
    write.csv(listVari2,paste0(dirSav,variety[i],"_complet.csv"), row.names = F)
    
    
    cualVar <- unlist(lapply(1:ncol(allData),function(x){!is.factor(allData[,x])}))
    
    if(sum(cualVar)!=length(cualVar))
    { 
      library(earth)
      
      cuantAllData <- allData[,cualVar]
      
      form <- formula(paste(names(allData)[ncol(allData)],"~ ."))
      
      namImp  <- names(listVari)
      
      listVari1 <- do.call(data.frame,lapply(1:(length(listVari)-1),function(x){z<-listVari[,x];if(!is.factor(z)){z}else{droplevels(z)} }))
      
      listVari[,1:(length(listVari)-1)] <- listVari1
      
      names(listVari)  <- namImp
      
      dummies <- dummyVars(form, data = listVari)
      
      tranfVars <- predict(dummies, newdata = listVari)
      
      unlistvari <- as.data.frame(unCorrMatrix(tranfVars,dirSav,cor.reduce = corRed))
      
      unlistvari$output <-  allData[,ncol(allData)]
      
      corScheme(as.data.frame(tranfVars),unlistvari,dirSav)
      
    }else{cuantAllData <- allData ; unlistvari <- as.data.frame(unCorrMatrix(listVari,dirSav,cor.reduce = corRed)) ; corScheme(as.data.frame(cuantAllData ),unlistvari,dirSav)}
    
    
    
    
    corMatrix <- cor(cuantAllData)
    
    write.csv(corMatrix,paste0(dirSav,variety[i],"_corMatrix.csv"))
    
    
    
    #saveFennixFormat(unlistvari,variety[i],paste0(variety[i],"/ARTIFICIAL_NEURAL_NETWORK/"))
    #nInputList[i] <- ncol(unlistvari)-1 
    #summVar <- completeSummary(unlistvari)
    
    write.csv(unlistvari,paste0(dirSav,variety[i],"_reduced.csv"))
    
  }
  #return(nInputList)
}




#VSURF ALTERNATIVE

vSurFun <- function(variety,dirLocation=paste0(getwd(),"/"),nCor=1)
{
  require(VSURF)  
  
  dirDataSet <- paste0(dirLocation,variety,"/DATA_SETS/",variety,"_complet.csv")
  dirSave    <- paste0(dirLocation,variety,"/VSURF/")
  
  dataSets   <- lapply(dirDataSet,function(x){read.csv(x,row.names=1)})
  
  for(i in 1:length(variety))
  {
    listVari  <- dataSets[[i]]
    inpVarMat <- listVari[,1:(ncol(listVari)-1)]
    outVarMat <- listVari[,ncol(listVari)]    
    
    vSurfResoults <- VSURF.parallel(x=inpVarMat,y=outVarMat,ncores = nCor,ntree = 2000,mtry=dim(inpVarMat)/3)
    
    relvar <- round( vSurfResoults$ord.imp$x/sum(vSurfResoults$ord.imp$x),3)*100  
    namVar <- names(inpVarMat)[vSurfResoults$ord.imp$ix]  
    
    namVar <- factor(namVar,levels=namVar)
    threes <-  length(vSurfResoults$varselect.thres)+0.5
    
    metriData  <- data.frame(namVar,relvar,threes)     
    
    
    cols <- array("royalblue1",ncol(inpVarMat))
    cols[namVar %in% names(inpVarMat)[vSurfResoults$varselect.interp ]] <- "red"
    
    tiff(paste0(dirSave[i],variety[i],"_relevance_vsurf.tiff"),
         width = 800, height = 600)
    
    k <-  ggplot(metriData,aes(namVar, relvar))+
      geom_bar( stat="identity",fill=cols,width=0.6,color="black")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme(axis.title.y=element_text(size=15),
            axis.text= element_text(colour = "gray26",size = 13))+xlab("")+
      ylab("% Sensitivity (Relevance)")+
      theme(panel.background = element_rect(colour = "gray38")) +
      theme(panel.background = element_rect(colour = "gray32"))
    
    kf <- k + theme(panel.background = element_rect(fill = 0,colour = "gray"))+
      geom_line(aes(x = threes, y = relvar),linetype=2) + geom_hline(h=0)
    
    kf <- kf + ggtitle(paste("Vsurf - Variety ", variety[i]))
    
    print(kf)
    
    dev.off()
    print(variety[i])
  }  
  
}

varImportance <- function(model, pred.data=model$trainingData, ..., scale=T) 
{
  library(caret)
  inputs <- pred.data[,-ncol(pred.data)]
  response <- pred.data[,ncol(pred.data)]
  n <- nrow(pred.data)
  
  errors <- rep(0, ncol(inputs))
  errors <- as.data.frame(errors)
  row.names(errors) <- names(inputs)
  
  for(input in names(inputs)) {
    
    xname <- ifelse(is.character(input),
                    input,
                    ifelse(is.name(input),
                           deparse(input),
                           eval(input)))
    
    n.pt = min(length(unique(pred.data[, xname])), 50)
    
    err <- 0
    
    xv <- pred.data[, xname]
    if (is.factor(xv) && !is.ordered(xv)) {
      x.pt <- levels(xv)
      mode <- xv[which.max(tabulate(match(xv, x.pt)))]
      x.data <- pred.data
      x.data[, xname] <- rep(mode, n)
      mode.err <- RMSE(predict(model, x.data), response)
      for (i in seq(along = x.pt)) {
        x.data <- pred.data
        x.data[, xname] <- factor(rep(x.pt[i], n), levels = x.pt)
        input.cancellation <- (RMSE(predict(model, x.data), response) - mode.err)^2
        if(!is.na(input.cancellation)) err <- err + input.cancellation
      }
    } else {
      if (is.ordered(xv)) 
        xv <- as.numeric(xv)
      x.pt <- seq(min(xv), max(xv), length = n.pt)
      mean <- mean(xv)
      x.data <- pred.data
      x.data[, xname] <- rep(mean, n)
      mean.err <- RMSE(predict(model, x.data), response)
      for (i in seq(along = x.pt)) {
        x.data <- pred.data
        x.data[, xname] <- rep(x.pt[i], n)
        input.cancellation <- (RMSE(predict(model, x.data), response) - mean.err)^2
        if(!is.na(input.cancellation)) err <- err + input.cancellation
      }
    }
    errors[input,] <- err
  }
  
  if(scale) {
    errors <- errors - min(errors, na.rm = TRUE)
    errors <- errors/max(errors, na.rm = TRUE) * 100
  }
  
  invisible(errors)
  return(errors)
}
