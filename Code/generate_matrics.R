# 
# setwd('C:/Users/yojia/Box/Data/simulated_data/spot_pattern/result_file')
# 
# file_tot = list.files('.', pattern = ".RData")
# 
# for(i in 1:length(file_tot)){
#   load(file_tot[i])
#   if(method == "spark") {method = "SPARK";rank_raw_pvalue = rank(raw_pvalue,ties.method = 'average')
#   rank_adjusted_pvalue = rank(adjusted_pvalue,ties.method = 'average')}
#   if(method == 'binSpect_kmeans') {method = "BinSpect_kmeans";rank_raw_pvalue = rank(raw_pvalue,ties.method = 'average')
#   rank_adjusted_pvalue = rank(adjusted_pvalue,ties.method = 'average')}
#   if(method == 'binSpect_rank') {method = "BinSpect_rank";rank_raw_pvalue = rank(raw_pvalue,ties.method = 'average')
#   rank_adjusted_pvalue = rank(adjusted_pvalue,ties.method = 'average')}
#   if(method == 'spatialDE') {method = "SpatialDE"}
#   if(method == 'trendsceek') {method = "Trendsceek"}
#   if(is.list(version)) {version = as.character(version)}
#   
#   save(raw_pvalue,rank_raw_pvalue,adjusted_pvalue,rank_adjusted_pvalue,gamma,method, version,file = file_tot[i])
# }
auc_bm = function(ap,gamma){
  TPR = rep(NA,101)
  FPR = rep(NA,101)
  treshd = seq(0,1,0.01)
  for(i in 1:101){
    adjusted = ifelse(ap <= treshd[i], TRUE,FALSE)
    # print(adjusted)
    confuse <- table(adjusted,gamma)
    # print(confuse)
    if(nrow(confuse)==1 &rownames(confuse)[1]=='TRUE'){
      confuse = rbind( c(0,0),confuse )
      row.names(confuse) = c("FALSE","TRUE")
    }
    if(nrow(confuse)==1& rownames(confuse)[1]=='FALSE'){
      confuse = rbind( confuse,c(0,0) )
      row.names(confuse) = c("FALSE","TRUE")
    }
    #print(confuse)
    TN = confuse[1]
    FP = confuse[2]
    FN = confuse[3]
    TP = confuse[4]
    TPR[i] = TP/(TP+FN)
    FPR[i] = FP/(FP+TN)
  }
  auc= 0
  for(j in 1:100){
    auc = auc + (FPR[j+1]-FPR[j])*(TPR[j+1]+TPR[j])/2
    #print(auc)
  }
  return(auc)
}

auc_bm_bs = function(ap,gamma){
  aa = range(ap,na.rm = T)
  treshd = seq(aa[1],aa[2],length.out = 1000)
  TPR = rep(NA,1000)
  FPR = rep(NA,1000)
  for(i in 1:1000){
    if(g==3){
      adjusted = ifelse(ap >=treshd[i]|ap <= 1/treshd[i], TRUE,FALSE)
    }
    if(g!=3){
      adjusted = ifelse(ap >=treshd[i], TRUE,FALSE)
    }
    
    # print(adjusted)
    confuse <- table(adjusted,gamma)
    # print(confuse)
    if(nrow(confuse)==1 &rownames(confuse)[1]=='TRUE'){
      confuse = rbind( c(0,0),confuse )
      row.names(confuse) = c("FALSE","TRUE")
    }
    if(nrow(confuse)==1& rownames(confuse)[1]=='FALSE'){
      confuse = rbind( confuse,c(0,0) )
      row.names(confuse) = c("FALSE","TRUE")
    }
    #print(confuse)
    TN = confuse[1]
    FP = confuse[2]
    FN = confuse[3]
    TP = confuse[4]
    TPR[i] = TP/(TP+FN)
    FPR[i] = FP/(FP+TN)
  }
  auc= 0
  for(j in 2:1000){
    auc = auc + (FPR[j-1]-FPR[j])*(TPR[j-1]+TPR[j])/2
    #print(auc)
  }
  
  return(auc)
}


perf_metric <- function(confusion_matrix){
  TN = confusion_matrix[1]
  FP = confusion_matrix[2]
  FN = confusion_matrix[3]
  TP = confusion_matrix[4]
  # colnames(result_sub) = c("Sensitivity","Specificity" ,"F1 score" ,"FDR", "AUC" ,"MCC")
  Sensitivity = ifelse(TP==0,0,TP/(TP+FN))
  Specificity = ifelse(TN==0,0,TN/(TN+FP))
  F1_score = ifelse(TP==0,0,2*TP/(2*TP+FP+FN))
  FDR = ifelse(FP==0,0,FP/(FP+TP))
  AUC = auc_bm(adjusted_pvalue,gamma)
  MCC = ifelse((TP+FP)==0|(TP+FN)==0|(TN+FP)==0|(TN+FN)==0,0,
               (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  return(list(Sensitivity =Sensitivity,Specificity = Specificity,F1_score = F1_score,FDR = FDR, AUC = AUC,MCC = MCC ))
}

perf_metric_bs <- function(confusion_matrix){
  TN = confusion_matrix[1]
  FP = confusion_matrix[2]
  FN = confusion_matrix[3]
  TP = confusion_matrix[4]
  # colnames(result_sub) = c("Sensitivity","Specificity" ,"F1 score" ,"FDR", "AUC" ,"MCC")
  Sensitivity = ifelse(TP==0,0,TP/(TP+FN))
  Specificity = ifelse(TN==0,0,TN/(TN+FP))
  F1_score = ifelse(TP==0,0,2*TP/(2*TP+FP+FN))
  FDR = ifelse(FP==0,0,FP/(FP+TP))
  AUC = auc_bm_bs(bf,gamma)
  MCC = ifelse((TP+FP)==0|(TP+FN)==0|(TN+FP)==0|(TN+FN)==0,0,
               (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  return(list(Sensitivity =Sensitivity,Specificity = Specificity,F1_score = F1_score,FDR = FDR, AUC = AUC,MCC = MCC ))
}
perf_metric_tk <- function(confusion_matrix){
  TN = confusion_matrix[1]
  FP = confusion_matrix[2]
  FN = confusion_matrix[3]
  TP = confusion_matrix[4]
  # colnames(result_sub) = c("Sensitivity","Specificity" ,"F1 score" ,"FDR", "AUC" ,"MCC")
  Sensitivity = ifelse(TP==0,0,TP/(TP+FN))
  Specificity = ifelse(TN==0,0,TN/(TN+FP))
  F1_score = ifelse(TP==0,0,2*TP/(2*TP+FP+FN))
  FDR = ifelse(FP==0,0,FP/(FP+TP))
  AUC = auc_bm(adjusted_pvalue[,s],gamma)
  MCC = ifelse((TP+FP)==0|(TP+FN)==0|(TN+FP)==0|(TN+FN)==0,0,
               (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  return(list(Sensitivity =Sensitivity,Specificity = Specificity,F1_score = F1_score,FDR = FDR, AUC = AUC,MCC = MCC ))
}

#-----------------------------------------------------------
# geneate table 
#-----------------------------------------------------------
pattern = c("bc_pattern","linear_pattern","mob_i_pattern","mob_ii_pattern","spot_pattern")
tt <- data.frame(pattern = NA,zero_setting = NA, replicate=NA,method=NA,metrics=NA,value=NA)
for(g in 1:5){
  setwd(paste0("C:\\Users\\yojia\\Box\\Data\\simulated_data\\",pattern[g],"\\result_file"))
  file_tot = list.files(".",pattern = ".RData")
  library(mltools)
  library(stringr)
  # bin_km <- file_tot[startsWith(file_tot,"binSpect_kmeans" )]
  # bin_rk <- file_tot[startsWith(file_tot,"binSpect_rank" )]
  # spark <- file_tot[startsWith(file_tot, "spark")]
  # spatialde <-file_tot[startsWith(file_tot, "spatialde")]
  trendsceek <- file_tot[startsWith(file_tot, "trendsceek")] # handle independently
  non_trendsceek <- file_tot[!startsWith(file_tot, "trendsceek")&!startsWith(file_tot, "boost_ising")]
  #method <- c("BinSpect_kmeans","BinSpect_rank","SPARK","Trendsceek","SpatialDE")
  boost_ising <- file_tot[startsWith(file_tot, "boost_ising")] # handle independently
  #colnames(tt) <- c("pattern","zero_setting", "replicate","method","metrics","value")
  for(i in 1:length(non_trendsceek)){
    load(non_trendsceek[i])
    #print(non_trendsceek[i])
    ww = ifelse(adjusted_pvalue <0.05, TRUE,FALSE)
    confuse_matrix <- table(ww,gamma)
    if(nrow(confuse_matrix)==1 & rownames(confuse_matrix)[1]=='TRUE'){
      confuse_matrix = rbind( c(0,0),confuse_matrix )
      row.names(confuse_matrix) = c("FALSE","TRUE")
    }
    if(nrow(confuse_matrix)==1& rownames(confuse_matrix)[1]=='FALSE'){
      confuse_matrix = rbind( confuse_matrix,c(0,0) )
      row.names(confuse_matrix) = c("FALSE","TRUE")
    }
    
    # sensitivity 
    result <- t(as.matrix(perf_metric(confuse_matrix)))
    colnames(result) = c("Sensitivity","Specificity" ,"F1 score" ,"FDR", "AUC" ,"MCC")
    for(j in 1:6){
      replicate = str_sub(file_tot[i],start =str_locate(file_tot[i],"replicate")[2]+2,end =str_locate(file_tot[i],".RData")[1]-1)
      zero_set <- str_sub(file_tot[i],start =str_locate(file_tot[i],"zero_")[2]+1,end =str_locate(file_tot[i],"_replicate")[1]-1)
      tt[(nrow(tt)+1),] = c(pattern[g],zero_set,replicate,method,colnames(result)[j],result[j])
    }
    rm(raw_pvalue)
    rm(adjusted_pvalue)
    rm(method)
    rm(version)
    rm(confuse_matrix)
    rm(ww)
  }
  tk_method <- c("Treendsceek_Emark","Treendsceek_Vmark","Treendsceek_markcorr","Treendsceek_markvario")
  for(i in 1:length(trendsceek)){
    load(trendsceek[i])
    for( s in 1:4){
      #print(non_trendsceek[i])
      
      ts = ifelse(adjusted_pvalue[,s] <0.05, TRUE,FALSE)
      confuse_matrix <- table(ts,gamma)
      if(nrow(confuse_matrix)==1 & rownames(confuse_matrix)[1]=='TRUE'){
        confuse_matrix = rbind( c(0,0),confuse_matrix )
        row.names(confuse_matrix) = c("FALSE","TRUE")
      }
      if(nrow(confuse_matrix)==1& rownames(confuse_matrix)[1]=='FALSE'){
        confuse_matrix = rbind( confuse_matrix,c(0,0) )
        row.names(confuse_matrix) = c("FALSE","TRUE")
      }
      
      # sensitivity 
      result <- t(as.matrix(perf_metric_tk(confuse_matrix)))
      colnames(result) = c("Sensitivity","Specificity" ,"F1 score" ,"FDR", "AUC" ,"MCC")
      for(j in 1:6){
        replicate = str_sub(file_tot[i],start =str_locate(file_tot[i],"replicate")[2]+2,end =str_locate(file_tot[i],".RData")[1]-1)
        zero_set <- str_sub(file_tot[i],start =str_locate(file_tot[i],"zero_")[2]+1,end =str_locate(file_tot[i],"_replicate")[1]-1)
        tt[(nrow(tt)+1),] = c(pattern[g],zero_set,replicate,tk_method[s],colnames(result)[j],result[j])
      }
      
    }
    rm(raw_pvalue)
    rm(adjusted_pvalue)
    rm(method)
    rm(version)
  }
  for(i in 1:length(boost_ising)){
    load(boost_ising[i])
    if(g==3){
      ww = ifelse(bf>150 | bf<1/150, TRUE,FALSE)
    }
    if(g!=3){
      ww = ifelse(bf >150, TRUE,FALSE)
    }
    confuse_matrix <- table(ww,gamma)
    if(nrow(confuse_matrix)==1 & rownames(confuse_matrix)[1]=='TRUE'){
      confuse_matrix = rbind( c(0,0),confuse_matrix )
      row.names(confuse_matrix) = c("FALSE","TRUE")
    }
    if(nrow(confuse_matrix)==1& rownames(confuse_matrix)[1]=='FALSE'){
      confuse_matrix = rbind( confuse_matrix,c(0,0) )
      row.names(confuse_matrix) = c("FALSE","TRUE")
    }
    
    # sensitivity 
    result <- t(as.matrix(perf_metric_bs(confuse_matrix)))
    colnames(result) = c("Sensitivity","Specificity" ,"F1 score" ,"FDR", "AUC" ,"MCC")
    for(j in 1:6){
      replicate = str_sub(file_tot[i],start =str_locate(file_tot[i],"replicate")[2]+2,end =str_locate(file_tot[i],".RData")[1]-1)
      zero_set <- str_sub(file_tot[i],start =str_locate(file_tot[i],"zero_")[2]+1,end =str_locate(file_tot[i],"_replicate")[1]-1)
      tt[(nrow(tt)+1),] = c(pattern[g],zero_set,replicate,method,colnames(result)[j],result[j])
    }
    rm(raw_pvalue)
    rm(adjusted_pvalue)
    rm(method)
    rm(version)
    rm(confuse_matrix)
    rm(ww)
  }
}

tt<- tt[-1,]
tt_total <- aggregate(tt,by = list(tt$pattern ,tt$zero_setting,tt$method,tt$metrics), FUN = mean)
tt_total <- tt_total[,-c(5:9)]
colnames(tt_total) <- c("pattern","zero_setting","method","metrics","value")
tt_total$value <- round(tt_total$value,3)
Path <- "C:/Users/yojia/Box/Data/result.csv"
write.csv(tt_total,file = Path)  

