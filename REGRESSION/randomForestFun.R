## All subsequent models are then run in parallel

randomForestFun <- function(variety,dirLocation=paste0(getwd(),"/"),saveWS=F,
                            barplot=FALSE,col.grap="lightskyblue",nb.it = 100,
                            wid=500,hei=800,ab=7,iz=4.1,ar=4.1,de=2.1,ncores=21,
                            sztxty=15,sztxtx=15,szlbty=15,szlbtx=15,szmain=15,
                            pp.szmain=15,pp.sztxtx=15,pp.sztxty=18,pp.szlbty=18,
                            pp.szlbtx=15,pp.lgndtx=15)
{
    ngw <- nchar(dirLocation)
    if( substring(dirLocation,ngw-16,ngw)=="VARIETY_ANALYSIS/" ){}else{return(cat("Aun no se encuentra en la carpeta VARIETY_ANALYSIS\nUtilize la funcion setwd para dirigirse a este carpeta"))}
    
  sfInit(parallel=T,cpus=ncores)
  sfLibrary(caret)
  sfLibrary(party)
  
  dirDataSet <- paste0(dirLocation,variety,"/DATA_SETS/",variety,"_complet.csv")
  dirSave    <- paste0(dirLocation,variety,"/RANDOM_FOREST/")  
  
  dataSets   <- lapply(dirDataSet,function(x){read.csv(x,row.names=1)})
  
  cat(paste("random Forest with Conditional Importance:\n"))
  
  cForestCaret <- function(x)
  {  
    setseed <- .Random.seed[1:nb.it]
    nOutPut <- ncol(data)
    
    #  pb <- winProgressBar(title="Progress bar", label="0% done", min=0, max=100, initial=0)    
    
    
    
    #  info <- sprintf("%d%% done", round((i/(nb.it)*100)))
    # setWinProgressBar(pb, paste(i/(nb.it)*100), label=info)
    
    inTrain <- createDataPartition(y=data[,nOutPut], p=0.7, list=F)
    training <- data[inTrain,]
    testing <- data[-inTrain,]
    
    
    
    grid <- expand.grid(mtry=round((ncol(training)-1)/3))
    
    
    model <- train(training[,-ncol(training)], training[,nOutPut],
                   method="rf", tuneGrid=grid,importance = TRUE,ntree = 2000)
    
    performance <- R2(predict(model, testing), testing[,nOutPut]) * 100
    
    vaRelevance <- varImp(model, scale=T)$importance
    
    return(list(model,performance,vaRelevance))
    
  }
  # close(pb) 
  sfExport("cForestCaret")
  
  sfExport("nb.it")
  
  for(j in 1:length(variety))
  { 
    

    cat(paste(j,"- Variety:",variety[j],"\n"))
    
    data0 <- dataSets[[j]]
    nvz <- nearZeroVar(data0)
    if(length(nvz)==0){data <- data0}else{data <- data0[,-nvz] }
    v <- integer()
    
    profiles <- list()
    length(profiles) <- length(names(data)[-ncol(data)])
    names(profiles) <- names(data)[-ncol(data)]

   
    sfExport("data")
  

    
    cat(paste("Running ", nb.it, "models in cross validation\n"))
    
    Sys.time()->start
    cForestModels <- sfLapply(1:nb.it,cForestCaret)
    print(Sys.time()-start)
    
    allModels <- lapply(cForestModels,function(x){x[[1]]})
    #allRMSE   <- unlist(lapply(allModelsAndRMSE,function(x){x[[2]]}))
    performance     <- unlist(lapply(cForestModels,function(x){x[[2]]}))
    relevances    <- lapply(cForestModels,function(x){x[[3]]})
    
    
    currentVarImp <- do.call(cbind,relevances)
    
    #sort(apply(do.call(cbind,currentVarImp),1,mean),decreasing = T)
    
    scale <- performance / as.numeric(apply(currentVarImp,2,sum))
    
    
    scaledVarImp <-  t(t(currentVarImp) * scale)
    
    ord <- list(0) ; for(k in 1:ncol(scaledVarImp)){ord[[k]] <- scaledVarImp[,k]}
    ordered <- lapply(ord,function(x){sort(x,decreasing = T)})
    
    
    
    princVar <- lapply(ordered,function(x){names(x)[1:3]})
    
    cat(paste("Computing profiles\n"))

    #profilesList <- sfLapply(1:100,function(x){ profLis <- list(0,0,0);names(profLis) <- princVar[[x]] ;for(n in princVar[[x]]){ profil <- profilePlot(allModels[[x]], n, data, F) ; profLis[[n]] <- data.frame(profil$y) ; row.names(profLis[[n]]) <- profil$x };return(profLis)})    
    
    Sys.time()->star
    profilesList <- lapply(1:nb.it,function(x){ profLis <- list(0,0,0);names(profLis) <- princVar[[x]] ;for(n in princVar[[x]]){ profil <- profilePlot(allModels[[x]], n, data, F) ; profLis[[n]] <- data.frame(profil$y) ; row.names(profLis[[n]]) <- profil$x };return(profLis)})    
    print(Sys.time()-start)

    
    
    for(z in  1:length(ordered))
    {  
      
      toProfile <- profilesList[[z]]
      
      for(n in names(toProfile)) {
        
        profile <- toProfile[n]
        
        if(length(profiles[[n]]) == 0) {
          
          profiles[[n]] <- as.data.frame(profile)
          
          #names(profiles)[n] <- n
          #row.names(profiles[[n]]) <- profile$x
          
        } else {
          profiles[[n]] <- cbind(profiles[[n]],as.data.frame(profile))
        }
      }
    }
    
    
 
    
    
    
    v <- as.data.frame(scaledVarImp)
    write.csv(v,paste0(dirSave[j],"weighMatrix.csv"))
    
    perf1 <- signif(sum(performance) / nb.it, 5)
    
    if(barplot){
        #Comienzo de barPlot
        
        se <- apply(v, 1, function(x){ 1.96*sd(x, na.rm=TRUE)/sqrt(ncol(v))})
        se <- data.frame(se,names(se))
        names(se) <- c("se","Variable")
        
        ordered <- sort(apply(v,1, median), decreasing=F)
        
        mean <- as.data.frame(ordered)
        mean <- cbind(mean, names(ordered))
        names(mean) <- c("Mean", "Variable")
        mean$Variable <- factor(names(ordered), levels= names(ordered))
        
        stadistc <- merge(se,mean,by.x="Variable",by.y="Variable")
        
        stadistc <- stadistc[order(stadistc$Mean,decreasing=F),]
        
        errBars <- transform(stadistc, lower=Mean-se,upper=Mean+se )
        
        
        png(paste0(dirSave[j],"InputRelvance.png"),width = wid, hei = hei,
            pointsize = 20,res=80)
        
        m <- ggplot(mean, aes(x=Variable, y=Mean))
        m <- m + geom_bar(stat="identity", width=0.5, fill="slategray1") +
             ylab("Mean importance")+ geom_errorbar(aes(ymax = lower, 
             ymin=upper), width=0.25,data=errBars) + coord_flip() + theme_bw() +
             ggtitle(paste("Importance of variables \n(with a mean R2 of",
             perf1, "%)")) +theme(plot.title = element_text(size = szmain, 
             face = "bold", colour = "black", vjust = 1.5),
             axis.text.y =element_text(size = sztxty),
             axis.text.x =element_text(size = sztxtx),
             axis.title.x = element_text(size = szlbty),
             axis.title.y = element_text(size = szlbtx))
         suppressWarnings(print(m))
        
        dev.off()
    }else{
        require(cowplot)
        #Comienzo boxplot
        
        newV <-  melt(t(v))[,-1]
        
        names(newV) <-  c("variable","value")
        
        medOrdenada <- names(with(newV,sort(tapply(value,variable,median))))
        
        newV$variable <- factor(newV$variable,levels=medOrdenada)
        
        noParameOut <- kruskal(newV$value,newV$variable,group = T)
        
        groupsData <- data.frame(noParameOut$groups$trt,noParameOut$groups$M)
        
        groupsData$noParameOut.groups.trt <- str_replace_all(groupsData$noParameOut.groups.trt, pattern=" ", repl="")
        
        maxDist <- {maxDis <- tapply(newV$value,newV$variable,max)+4.5;data.frame(nam=names(maxDis),max=maxDis)}
        
        groupsData <- merge(groupsData,maxDist,by.x="noParameOut.groups.trt",by.y="nam",all=T,sort=F)
        
        newV1 <- merge(newV,groupsData,by.x="variable",by.y="noParameOut.groups.trt",all.x=T,all.y=F,sort = F)
        
        
        png(paste0(dirSave[j],"InputRelvance.png"),width = wid, hei = hei, pointsize = 20)
        m <- ggplot(newV1, aes(x=variable, y=value))
        m <- m + geom_boxplot(fill=col.grap) + ylab("Importance")+
            xlab("Input variable")+theme_bw() + 
            ggtitle(paste("Importance of variables\n(with a mean R2 of", perf1,
            "%)")) + theme(axis.text.x = element_text(angle=0, hjust=0.5, 
            vjust=0,size=sztxtx),plot.title = element_text(vjust=3,size=szmain),
            axis.text.y =element_text(size = sztxty),
            axis.title.x = element_text(size = szlbty),
            axis.title.y = element_text(size = szlbtx))+ coord_flip()+
            geom_text(aes(y = max,label = noParameOut.groups.M))
        
        print(ggdraw(switch_axis_position(m, 'x')))
        dev.off()
    }
    #Fin del grafico boxplot
    
    
    
    namSort <- names(sort(apply(v,1, median), decreasing=T))
    
    profData     <- unlist(lapply(profiles,function(x){!is.null(x)}))
    profRealData <- names(profData)[profData]
    
    limProf <- if(length(profRealData) < 5){ length(profRealData)}else{5}
    
    
    for(i in 1:limProf)
    {
        if(!is.null(unlist(profiles[namSort[i]])))
        { 
            png(paste0(dirSave[j],"MultiProfile_",namSort[i],".png"),width =,650, hei =410 , pointsize = 40)
            multiProfile(data,profiles,namSort[i],pp.szmain=pp.szmain,
                         pp.sztxtx=pp.sztxtx,pp.sztxty=pp.sztxty,
                         pp.szlbty=pp.szlbty,pp.szlbtx=pp.szlbtx,
                         pp.lgndtx=pp.lgndtx)
            dev.off()
        } else{print(paste("Few profiles references for:",namSort[i]))}
    }
    
    if(saveWS==T){save(list = ls(all = TRUE), file = paste0(dirSave[j],"workSpace.RData"))}else{}
  

  }
  sfStop()
}
