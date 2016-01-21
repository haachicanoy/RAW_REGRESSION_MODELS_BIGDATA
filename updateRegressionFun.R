
dirServ <- "C:/Users/haachicanoy/Documents/GitHub/RAW_REGRESSION_MODELS_BIGDATA"

lisFil <- list.files(dirServ,full.names = T,pattern = "*.R$",recursive = T)
lisFil <- lisFil[-grep("updateRegressionFun",lisFil)]
lisFil <- lisFil[-grep("multiProfileFix",lisFil)]

lapply(lisFil,source)

# save(list = ls(all.names = TRUE), file = "//dapadfs/workspace_cluster_6/TRANSVERSAL_PROJECTS/MADR/COMPONENTE_2/OPEN_BIGDATA_AEPS/REGRESSION_MODELS/All-Functions-AEPS_BD.RData")
save(list = ls(all.names = TRUE), file = "C:/Users/haachicanoy/Documents/GitHub/RAW_REGRESSION_MODELS_BIGDATA/All-Functions-AEPS_BD.RData")
