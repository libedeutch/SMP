# normalization 
# count matrix pxn
# library(SPARK)
Path = "C:\\Users\\yojia\\Box\\Data\\simulated_data"
patterns = c("bc_pattern","linear_pattern","mob_i_pattern","mob_ii_pattern","spot_pattern")
for(w in 1:5){
  file_tot <- list.files(paste0(Path,'/',patterns[w],'/data'),pattern = '.RData')
  print(file_tot)
  for(i in 1:length(file_tot)){

  varx = apply(count, 2, var)
  meanx = apply(count, 2, mean)
  phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
  ## regress out log total count
  norm_counts <- log(count + 1/(2 * phi))
  total_counts <- apply(count, 1, sum)
  res_norm_counts <- t(apply(norm_counts, 1, function(x){resid(lm(x ~ log(total_counts)))} ))
  
  # method 2 
  load(paste0(Path,'/',patterns[w],'/data','/',file_tot[i]))
 #normalized_count = NormalizeVST(t(count))
  #normalized_count = t(normalized_count)
  save(normalized_count,loc,file = paste0(Path,'/',patterns[w],'/normalized_data/',"normalized_",file_tot[i]))                    
 }
}




# ST
Path = "C:\\Users\\yojia\\Box\\Data\\st"
patterns = c("human_breast_cancer","human_heart","human_melanoma","human_prostate","mouse_olfactory_bulb")
for(w in 2:2){
  file_tot <- list.files(paste0(Path,'/',patterns[w+3],'/processed_data'),pattern = '.Rdata')
  for(i in 1:length(file_tot)){
    # # method 1 
    varx = apply(count, 2, var)
    meanx = apply(count, 2, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
    ## regress out log total count
    norm_counts <- log(count + 1/(2 * phi))
    total_counts <- apply(count, 1, sum)
    res_norm_counts <- t(apply(norm_counts, 1, function(x){resid(lm(x ~ log(total_counts)))} ))
    # method 2 
    load(paste0(Path,'/',patterns[w+3],'/processed_data','/',file_tot[i]))
    normalized_count = NormalizeVST(t(count))
    normalized_count = t(as.matrix(normalized_count))
    save(normalized_count,loc,gene_name,gene_id,file = paste0(Path,'/',patterns[w+3],'/normalized_data/',"normalized_",file_tot[i]))                    
  }
}

#fish
Path = "C:\\Users\\yojia\\Box\\Data\\fish"
patterns = c("mouse_hippocampus","mouse_subventricular_zone_and_olfactory_bulb")
for(w in 1:2){
  file_tot <- list.files(paste0(Path,'/',patterns[w],'/processed_data'),pattern = '.Rdata')
  print(file_tot)
  for(i in 1:length(file_tot)){
    # # method 1 
    varx = apply(count, 2, var)
    meanx = apply(count, 2, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
    ## regress out log total count
    norm_counts <- log(count + 1/(2 * phi))
    total_counts <- apply(count, 1, sum)
    res_norm_counts <- t(apply(norm_counts, 1, function(x){resid(lm(x ~ log(total_counts)))} ))
    # method 2 
    load(paste0(Path,'/',patterns[w],'/processed_data','/',file_tot[i]))
    normalized_count = NormalizeVST(t(count))
    normalized_count = t(normalized_count)
    save(normalized_count,loc,gene_name,gene_id,file = paste0(Path,'/',patterns[w],'/normalized_data/',"normalized_",file_tot[i]))                    
  }
}