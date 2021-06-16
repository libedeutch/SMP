

#-----------------------------------------------------------
# geneate P value table 
#-----------------------------------------------------------
library(stringr)
library(tidyr)
pattern = c("bc_pattern","linear_pattern","mob_i_pattern","mob_ii_pattern","spot_pattern")
tt <- data.frame(gene_name = NA,pattern = NA,zero_setting = NA, replicate=NA,method=NA,P_value = NA)
for(g in 1:5){
  #g=5
  setwd(paste0("C:\\Users\\yojia\\Box\\Data\\simulated_data\\",pattern[g],"\\result_file"))
  #setwd(paste0("/home/jxy190005/simulated_data/",pattern[g],"\\result_file"))
  file_tot = list.files(".",pattern = ".RData")
  trendsceek <- file_tot[startsWith(file_tot, "trendsceek")] # handle independently
  #non_trendsceek <- file_tot[!startsWith(file_tot, "trendsceek")&!startsWith(file_tot, "spark")]
  non_trendsceek <- file_tot[!startsWith(file_tot, "trendsceek")&!startsWith(file_tot,"boost_ising")]
  print(non_trendsceek)
  for(i in 1:length(non_trendsceek)){
      load(non_trendsceek[i])
      replicate = str_sub(file_tot[i],start =str_locate(file_tot[i],"replicate")[2]+2,end =str_locate(file_tot[i],".RData")[1]-1)
      zero_set <- str_sub(file_tot[i],start =str_locate(file_tot[i],"zero_")[2]+1,end =str_locate(file_tot[i],"_replicate")[1]-1)
      for (j in 1:100){
        tt[(nrow(tt)+1),] = c(paste0("gene",j),pattern[g],zero_set,replicate,method,round(adjusted_pvalue[j],3))
      }
      
    rm(raw_pvalue)
    rm(adjusted_pvalue)
    rm(method)
    rm(version)
  }

  tk_method <- c("Treendsceek_Emark","Treendsceek_Vmark","Treendsceek_markcorr","Treendsceek_markvario")
  for(i in 1:length(trendsceek)){
    print(i)
    load(trendsceek[i])
    for( s in 1:4){
      for(j in 1:100){
        replicate = str_sub(file_tot[i],start =str_locate(file_tot[i],"replicate")[2]+2,end =str_locate(file_tot[i],".RData")[1]-1)
        zero_set <- str_sub(file_tot[i],start =str_locate(file_tot[i],"zero_")[2]+1,end =str_locate(file_tot[i],"_replicate")[1]-1)
        tt[(nrow(tt)+1),] = c(paste0("gene",j),pattern[g],zero_set,replicate,tk_method[s],round(adjusted_pvalue[j],3))
      }
      
    }
    rm(raw_pvalue)
    rm(adjusted_pvalue)
    rm(method)
    rm(version)
  }
}


tt<- tt[-1,]
attach(tt)
tt_total <- tt %>% spread(method,P_value)

Path <- "C:/Users/yojia/Box/Data/pvalue_result.csv"

write.csv(tt_total,file = Path, row.names = FALSE)

# with boost_ising
library(stringr)
library(tidyr)
library(dplyr)
pattern = c("bc_pattern","linear_pattern","mob_i_pattern","mob_ii_pattern","spot_pattern")
tt <- data.frame(gene_name = NA,pattern = NA,zero_setting = NA, replicate=NA,method=NA,P_value = NA)
for(g in 1:5){
  
  setwd(paste0("C:\\Users\\yojia\\Box\\Data\\simulated_data\\",pattern[g],"\\result_file"))
  
  file_tot = list.files(".",pattern = ".RData")
  boost_ising <- file_tot[startsWith(file_tot, "boost_ising")] # handle independently
  
  for(i in 1:length(boost_ising)){
    load(boost_ising[i])
    replicate = str_sub(file_tot[i],start =str_locate(file_tot[i],"replicate")[2]+2,end =str_locate(file_tot[i],".RData")[1]-1)
    zero_set <- str_sub(file_tot[i],start =str_locate(file_tot[i],"zero_")[2]+1,end =str_locate(file_tot[i],"_replicate")[1]-1)
    for (j in 1:100){
      tt[(nrow(tt)+1),] = c(paste0("gene",j),pattern[g],zero_set,replicate,method,round(bf[j],3))
    }
    
  }
  
}

tt<- tt[-1,]
tt <- tt[,-5]
pp = read.csv("C:/Users/yojia/Box/Data/pvalue_result.csv")
colnames(tt) <- c(colnames(tt)[1:4],"BOOST-Ising(beyas factor)")
mode(tt$zero_setting) <- "integer"
mode(tt$replicate) <- "integer"
tt_total <- left_join(pp,tt)
Path <- "C:/Users/yojia/Box/Data/pvalue_result.csv"

write.csv(tt_total,file = Path, row.names = FALSE)