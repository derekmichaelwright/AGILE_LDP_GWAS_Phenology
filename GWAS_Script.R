##################################################
root <- "D:/gitfolder/AGILE_LDP_GWAS_Phenology"
setwd(root)
library(agData)
##################################################
# Phenotypes
myY <- read.csv("myY.csv")
# Genotypes
myG <- read.csv("myG_LDP.csv", header = F)
#myG$V3 <- gsub("Lcu.2RBY.Chr", "", myG$V3)
#c(head(myG$V3), tail(myG$V3))
#myG <- myG %>% filter(V3 %in% c("chrom", "1", "2", "3", "4", "5", "6", "7"))
unique(myG$V3)
setdiff(myY$Name, t(myG[1,]))
setdiff(t(myG[1,]), myY$Name)
#
myCV <- myY %>% select(Name, PTModel_b, PTModel_c, Su17_Tb, Su17_Pc)
##################################################
colnames(myY)
#dtf_traits <- colnames(myY)[grepl("Name|DTF",colnames(myY))]
# Load library
library(GAPIT3) # devtools::install_github("jiabowang/GAPIT3",force=TRUE)
setwd(paste0(root,"/Results"))
myGAPIT <- GAPIT(
  Y = myY,
  G = myG,
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 4
)
setwd(paste0(root,"/Results_b"))
myGAPIT <- GAPIT(
  Y = myY%>%select(-PTModel_b),
  G = myG,
  CV = myCV[,c("Name","PTModel_b")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 0
)
setwd(paste0(root,"/Results_c"))
myGAPIT <- GAPIT(
  Y = myY%>%select(-PTModel_c),
  G = myG,
  CV = myCV[,c("Name","PTModel_c")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 0
)
setwd(paste0(root,"/Results_c"))
myGAPIT <- GAPIT(
  Y = myY%>%select(-PTModel_c),
  G = myG,
  CV = myCV[,c("Name","PTModel_c")],
  model = "MLMM",
  PCA.total = 0
)

setwd(paste0(root,"/Results_Tb"))
myGAPIT <- GAPIT(
  Y = myY%>%select(-Su17_Tb),
  G = myG,
  CV = myCV[,c("Name","Su17_Tb")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 0
)
setwd(paste0(root,"/Results_Pc"))
myGAPIT <- GAPIT(
  Y = myY%>%select(-Su17_Pc),
  G = myG,
  CV = myCV[,c("Name","Su17_Pc")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 0
)










setwd(paste0(root,"/Results"))
myGAPIT <- GAPIT(
  Y = myY[c("Name","PC2")],
  G = myG,
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 4
)
setwd(paste0(root,"/Results_b"))
myGAPIT <- GAPIT(
  Y = myY[,c("Name","PCA_PC2")],
  G = myG,
  CV = myCV[,c("Name","PTModel_b")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 0
)
setwd(paste0(root,"/Results_c"))
myGAPIT <- GAPIT(
  Y = myY%>%select(-PTModel_c),
  G = myG,
  CV = myCV[,c("Name","PTModel_c")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 0
)

setwd(paste0(root,"/Results"))
myGAPIT <- GAPIT(
  Y = myY[c("Name","Ne17_DTF")],
  G = myG,
  model = c("FarmCPU"),
  PCA.total = 4
)












setwd(paste0(root,"/Results_PC1"))
myGAPIT <- GAPIT(
  Y = myY[,dtf_traits],
  G = myG,
  CV = myCV_clusters[,c("Name","PC1")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 4
)
setwd(paste0(root,"/Results_PC2"))
myGAPIT <- GAPIT(
  Y = myY[,dtf_traits],
  G = myG,
  CV = myCV_clusters[,c("Name","PC2")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 4
)
setwd(paste0(root,"/Results_PC"))
myGAPIT <- GAPIT(
  Y = myY[,dtf_traits],
  G = myG,
  CV = myCV_clusters[,c("Name","PC1","PC2")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 4
)
#
setwd(paste0(root,"/Results_bPC"))
myGAPIT <- GAPIT(
  Y = myY[,dtf_traits],
  G = myG,
  CV = myCV[,c("Name","b")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 4
)
setwd(paste0(root,"/Results_cPC"))
myGAPIT <- GAPIT(
  Y = myY[,dtf_traits],
  G = myG,
  CV = myCV[,c("Name","c")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 4
)
setwd(paste0(root,"/Results_bc"))
myGAPIT <- GAPIT(
  Y = myY[,dtf_traits],
  G = myG,
  CV = myCV[,c("Name","b","c")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 0
)






traits <- c("Name","Ne17_DTF")
setwd(paste0(root,"/Results_new"))
myGAPIT <- GAPIT(
  Y = myY[,traits],
  G = myG,
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 4
)





# Run GWAS (Results moved to folder `Results/`)
traits <- c("Name","Model_a", "Model_b", "Model_c", 
            "PCA_PC1", "PCA_PC2", "PCA_PC3", "PCA_PC4", "PCA_PC5", "PCA_PC6")
traits <- c("Name","Ne17_DTF")
setwd(paste0(root,"/Results"))
myGAPIT <- GAPIT(
  Y = myY[,traits],
  G = myG,
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 4
)
myGAPIT <- GAPIT(
  Y = myY,
  G = myG,
  model = "Blink",
  PCA.total = 4
)
myGAPIT <- GAPIT(
  Y = myY[,traits],
  G = myG,
  model = "FarmCPU",
  PCA.total = 4
)
phen_traits <- c("Name",
  "Ro16_DTF", "Ro17_DTF", "Su16_DTF", "Su17_DTF", "Su18_DTF", "Us18_DTF",
  "In16_DTF", "In17_DTF", "Ba16_DTF", "Ba17_DTF", "Ne16_DTF", "Ne17_DTF",
  "Mo16_DTF", "Mo17_DTF", "Sp16_DTF", "Sp17_DTF", "It16_DTF", "It17_DTF",
  
  "Ro16_DTS", "Ro17_DTS", "Su16_DTS", "Su17_DTS", "Su18_DTS", "Us18_DTS",
  "In16_DTS", "In17_DTS", "Ba16_DTS", "Ba17_DTS", "Ne16_DTS", "Ne17_DTS",
  "Mo16_DTS", "Mo17_DTS", "Sp16_DTS", "Sp17_DTS", "It16_DTS", "It17_DTS",
  
  "Ro16_DTM", "Ro17_DTM", "Su16_DTM", "Su17_DTM", "Su18_DTM", "Us18_DTM",
  "In16_DTM", "In17_DTM", "Ba16_DTM", "Ba17_DTM", "Ne16_DTM", "Ne17_DTM",
  "Mo16_DTM", "Mo17_DTM", "Sp16_DTM", "Sp17_DTM", "It16_DTM", "It17_DTM",
  
  "Ro16_REP", "Ro17_REP", "Su16_REP", "Su17_REP", "Su18_REP", "Us18_REP",
  "In16_REP", "In17_REP", "Ba16_REP", "Ba17_REP", "Ne16_REP", "Ne17_REP",
  "Mo16_REP", "Mo17_REP", "Sp16_REP", "Sp17_REP", "It16_REP", "It17_REP"
)
# GWAS with b covariate (Results moved to folder `Results_b/`)
setwd(paste0(root,"/Results_b"))
myGAPIT <- GAPIT(
  Y = myY[,phen_traits],
  G = myG,
  CV = myCV[,c("Name","b")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 0
)
myGAPIT <- GAPIT(
  Y = myY[,phen_traits],
  G = myG,
  CV = myCV[,c("Name","b","Extra1","Extra2")],
  model = "Blink",
  PCA.total = 0
)
# GWAS with c covariate (Results moved to folder `Results_c/`)
setwd(paste0(root,"/Results_c"))
myGAPIT <- GAPIT(
  Y = myY[,phen_traits],
  G = myG,
  CV = myCV[,c("Name","c")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 0
)
myGAPIT <- GAPIT(
  Y = myY[,phen_traits],
  G = myG,
  CV = myCV[,c("Name","c")],
  model = "Blink",
  PCA.total = 0
)
myGAPIT <- GAPIT(
  Y = myY[,phen_traits[c(1,42:73)]],
  G = myG,
  CV = myCV[,c("Name","c")],
  model = "FarmCPU",
  PCA.total = 0
)
#
setwd(paste0(root,"/Results_ab"))
myGAPIT <- GAPIT(
  Y = myY[,phen_traits],
  G = myG,
  CV = myCV[,c("Name","a","b")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 0
)
#
phen_traits_dtf <- c("Name",
                 "Ro16_DTF", "Ro17_DTF", "Su16_DTF", "Su17_DTF", "Su18_DTF", "Us18_DTF",
                 "In16_DTF", "In17_DTF", "Ba16_DTF", "Ba17_DTF", "Ne16_DTF", "Ne17_DTF",
                 "Mo16_DTF", "Mo17_DTF", "Sp16_DTF", "Sp17_DTF", "It16_DTF", "It17_DTF"
)
setwd(paste0(root,"/Results_PCb"))
myGAPIT <- GAPIT(
  Y = myY[,phen_traits_dtf],
  G = myG,
  CV = myCV[,c("Name","b")],
  model = "Blink",
  PCA.total = 4
)
myGAPIT <- GAPIT(
  Y = myY[,phen_traits_dtf],
  G = myG,
  CV = myCV[,c("Name","b")],
  model = c("MLM","MLMM","FarmCPU"),#,"Blink"
  PCA.total = 4
)
#
setwd(paste0(root,"/Results_PCc"))
myGAPIT <- GAPIT(
  Y = myY[,phen_traits_dtf],
  G = myG,
  CV = myCV[,c("Name","b")],
  model = "Blink",
  PCA.total = 4
)
myGAPIT <- GAPIT(
  Y = myY[,phen_traits_dtf],
  G = myG,
  CV = myCV[,c("Name","b")],
  model = c("MLM","MLMM","FarmCPU","Blink"),
  PCA.total = 4
)
