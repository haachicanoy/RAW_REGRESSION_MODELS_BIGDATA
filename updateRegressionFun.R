
dirServ <- "C:/Users/haachicanoy/Documents/GitHub/RAW_REGRESSION_MODELS_BIGDATA/"

lisFil <- list.files(dirServ, full.names=T, pattern=".R", recursive=TRUE)
lisFil <- lisFil[-grep("updateRegressionFun",lisFil)]
lisFil <- lisFil[-grep("multiProfileFix",lisFil)]
lisFil <- lisFil[-grep("All-Functions-AEPS_BD_updated",lisFil)]

lapply(lisFil, source)

save(list=ls(all.names=TRUE), file="C:/Users/haachicanoy/Documents/GitHub/RAW_REGRESSION_MODELS_BIGDATA/All-Functions-AEPS_BD_updated.RData")
