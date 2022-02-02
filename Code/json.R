

if(!require(rjson)){
    install.packages("rjson")
    library(rjson)
}

if(!require(RJSONIO)){
    install.packages("RJSONIO")
    library(RJSONIO)
}

if(!require(stringr)){
    install.packages("stringr")
    library(stringr)
}


#--------------------------------------------------------
# Author: Jie Yang
# Date: 07/22/2021
# Revision Date: 02/01/2022
#--------------------------------------------------------

setwd("~/Box/Data")
first_dir <- list.files(".")

image<- function(path){#paste(getwd(),"/","fish/",second_dir_fish[i],"/pathology_image",sep = '')
  if(file.exists(path)){
    return(paste0("st/",second_dir_st[j],"/pathology_image/",list.files(path)))
  }
  else{
  return(list.files(path))
}
}
# cohorts 
second_dir_fish <- list.files(paste(getwd(),"/","fish",sep = ''))
second_dir_st <-  list.files(paste(getwd(),"/","st",sep = ''))
second_dir_hdst <-  list.files(paste(getwd(),"/","hdst",sep = ''))
second_dir_tenx <-  list.files(paste(getwd(),"/","10x",sep = ''))
second_dir_other <- list.files(paste(getwd(),"/","microarray_singlecell",sep = ''))
second_dir_slideseqv2 <-  list.files(paste(getwd(),"/","slide_seqv2",sep = ''))
second_dir_starmap <- list.files(paste(getwd(),"/","starmap",sep = ''))

# organisms
fish_organism <- c("mouse","mouse","mouse")
st_organism <-c("Human","Human","Human","Human","Mouse")
hdst_organism <- c("Human", "Mouse")
tenx_organism <- c(rep("Human",15), rep("Mouse",4))
other_organism<- "Human"
slideseqv2_organism<- rep("Mouse",5)
starmap_organism <-  rep("Mouse",3)

#tissue
fish_tissue <-c("Brain-hippocampus","Brain-hypothalamus","Brain-olfactory bulb","Brain-subventricular zone")
st_tissue   <-c("Breast","Heart","Skin","Prostate","Brain-olfactory bulb")
hdst_tissue<- c("Breast", "Brain-olfactory bulb")
#tenx_tissue <-c(str_to_title(str_remove_all(str_remove_all(str_remove(second_dir_tenx,"human_|mouse_|_cancer|_ffpe|_serial"),"_cancer"),'_ffpe|2|normal')))
tenx_tissue <-c("Brain-cerebral cortex", "Breast","Breast","Breast","Breast-cerebellum", "Colon/rectum",
                "Brain-prefrontal cortex","Brain","Heart", "Breast", "Lymph node", "Prostate","Ovarian", "Prostate",
                "Spinal cord", "Brain", "Brain", "Brain", "Kidney")
other_tissue <- "Pancreas"
slideseqv2_tissue <- c("Stem cell-embryo", "Stem cell-brain","Brain-hippocampus","Brain-olfactory bulb","Brain-somatosensory cortex")
starmap_tissue <- c("Brain-prefrontal cortex", rep("Brain-primary visual cortex",2))

# function to calculate zero proportion 
zero_propor_and_gene_location <- function(count_matrix_path){
  load(count_matrix_path)
  
  if(exists("count")) { return( list(number_of_gene = dim(count)[2],
                                     number_of_location = dim(count)[1],
                                     percentage_of_zero = sum(count ==0,na.rm = TRUE ) /(dim(count)[1]*dim(count)[2]) ) )}
  if(exists("count_norm")) { return(list(number_of_gene = dim(count_norm)[2],
                                         number_of_location = dim(count_norm)[1],
                                         percentage_of_zero = sum(count_norm ==0,na.rm = TRUE ) /(dim(count_norm)[1]*dim(count_norm)[2]) )) }
}

# calculate location, gene, zero proportion for fish data sets
zero_gene_location_fish_tt = list(list())
for(i in 1:length(second_dir_fish)){
  processed_data = list.files(paste(getwd(),"/","fish/",second_dir_fish[i],"/processed_data",sep = ''),pattern = ".Rdata")
  #print(processed_data)
  para_fish = matrix(nrow = 3, ncol = length(processed_data))
  row.names(para_fish) = c("number_of_genes","number_of_location","percentage_of_zero")
  for(j in 1:length(processed_data)){
  Path = paste(getwd(),"/","fish/",second_dir_fish[i],"/processed_data","/",processed_data[j],sep = '')
  #print(Path)
   tt = zero_propor_and_gene_location(Path)
   para_fish[1,j] = tt$number_of_gene
   para_fish[2,j] = tt$number_of_location
   para_fish[3,j] = tt$percentage_of_zero
  }
  #print(para_fish)
  zero_gene_location_fish_tt[[i]] = list(number_of_gene = para_fish[1,],number_of_location = para_fish[2,],percentage_of_zero = (para_fish[3,]))
}

# calculate location, gene, zero proportion for st data sets
zero_gene_location_st_tt = list(list())

for(i in 1:length(second_dir_st)){
  processed_data = list.files(paste(getwd(),"/","st/",second_dir_st[i],"/processed_data",sep = ''),pattern = ".Rdata")
  para_st = matrix(nrow = 3, ncol = length(processed_data))
  row.names(para_st) = c("number_of_genes","number_of_location","percentage_of_zero")
  for(j in 1:length(processed_data)){
    Path = paste(getwd(),"/","st/",second_dir_st[i],"/processed_data","/",processed_data[j],sep = '')
    #print(Path)
    tt = zero_propor_and_gene_location(Path)
    para_st[1,j] = tt$number_of_gene
    para_st[2,j] = tt$number_of_location
    para_st[3,j] = tt$percentage_of_zero
  }
  #print(para_st)
  zero_gene_location_st_tt[[i]] = list(number_of_gene = para_st[1,],number_of_location = para_st[2,],percentage_of_zero = para_st[3,])
}

# calculate location, gene, zero proportion for 10x data sets
zero_gene_location_tenx_tt = list(list())
for(i in 1:length(second_dir_tenx)){
  processed_data = list.files(paste(getwd(),"/","10x/",second_dir_tenx[i],"/processed_data",sep = ''),pattern = ".RData")
  para_tenx = matrix(nrow = 3, ncol = length(processed_data))
  row.names(para_tenx) = c("number_of_genes","number_of_location","percentage_of_zero")
  for(j in 1:length(processed_data)){
    Path = paste(getwd(),"/","10x/",second_dir_tenx[i],"/processed_data","/",processed_data[j],sep = '')
    tt = zero_propor_and_gene_location(Path)
    para_tenx[1,j] = tt$number_of_gene
    para_tenx[2,j] = tt$number_of_location
    para_tenx[3,j] = tt$percentage_of_zero
  }
  zero_gene_location_tenx_tt[[i]] = list(number_of_gene = para_tenx[1,],number_of_location = para_tenx[2,],percentage_of_zero = para_tenx[3,])
}

#----------------------------------------------------------------------------------------
# JSON UPDATE- 
# Add link, technology 
#----------------------------------------------------------------------------------------
link_fish <- c("http://dx.doi.org/10.1016/j.neuron.2016.10.001","http://dx.doi.org/10.1126/science.aau5324",
               "https://doi.org/10.1038/s41586-019-1049-y","https://doi.org/10.1038/s41586-019-1049-y")
tech_fish <- c("seqFISH", "MERFISH", "seqFISH","seqFISH")

# JSON for fish 
pre_json <- list(list())
for(i in 1:length(second_dir_fish)){
 pre_json[[i]]<- list(
    dataset_name =second_dir_fish[i], 
    display_name = gsub("_",' ',str_to_title(second_dir_fish[i])),
    dataset_type = "Real Data",
    disease  = NA,
    technology = tech_fish[i] ,
    organism = "Mouse" ,
    tissue= fish_tissue[i] ,
   # dataset_description = "Author: Jie Yang (jxy190005@utdallas.edu)",
    raw_data = paste0("fish/",second_dir_fish[i],"/raw_data/",list.files(paste(getwd(),"/","fish/",second_dir_fish[i],"/raw_data",sep = ''))),
    processed_data = paste0("fish/",second_dir_fish[i],"/processed_data/",list.files(paste(getwd(),"/","fish/",second_dir_fish[i],"/processed_data",sep = ''))),
    normalized_expression = paste0("fish/",second_dir_fish[i],"/normalized_data/",list.files(paste(getwd(),"/","fish/",second_dir_fish[i],"/normalized_data",sep = ''))),
    
    number_of_genes = (zero_gene_location_fish_tt[[i]][1])$number_of_gene ,
    number_of_location = (zero_gene_location_fish_tt[[i]][2])$number_of_location ,
    zero_percentage = round(100*(zero_gene_location_fish_tt[[i]][3]$percentage_of_zero),2),
    source_name = str_to_title(str_remove(list.files(paste(getwd(),"/","fish/",second_dir_fish[i],"/paper",sep = ''))[1],'.pdf')),
    source_link =link_fish[i] ,
    image_file = image(paste(getwd(),"/","fish/",second_dir_fish[i],"/pathology_image",sep = ''))
    #dataset_code = paste0("fish/",second_dir_fish[i],"/code/",list.files(paste(getwd(),"/","fish/",second_dir_fish[i],"/code",sep = '')))
  )
}
le = 0 # track the number of processed datasets 
# JSON for st datasets
link_st = c("http://dx.doi.org/10.1126/science.aaf2403","http://dx.doi.org/10.1038/s41598-017-13462-5", "http://dx.doi.org/10.1158/0008-5472.CAN-18-0747","http://dx.doi.org/10.1038/s41467-018-04724-5","10.1126/science.aaf2403")
disease_st <-c("breast cancer", "heart failure","melanoma","prostate cancer",NA)
for(j in (1):(length(second_dir_st))){
  le = le + length(second_dir_fish)
  pre_json[[le +j]]<- list(
    dataset_name =second_dir_st[j], 
    display_name = gsub("_",' ',str_to_title(second_dir_st[j])),
    dataset_type = "Real Data",
    disease = str_to_title(disease_st[j]),
    technology = "ST",#"Spatial Transcriptomic", 
    organism = st_organism[j] ,
    tissue = st_tissue[j],
  #  dataset_description = "Author: Jie Yang (jxy190005@utdallas.edu)",
    raw_data = paste0("st/",second_dir_st[j],"/raw_data/",list.files(paste(getwd(),"/","st/",second_dir_st[j],"/raw_data",sep = ''))),
    
    processed_data = paste0("st/",second_dir_st[j],"/processed_data/",list.files(paste(getwd(),"/","st/",second_dir_st[j],"/processed_data",sep = ''))),
    normalized_expression = paste0("st/",second_dir_st[j],"/normalized_data/",list.files(paste(getwd(),"/","st/",second_dir_st[j],"/normalized_data",sep = ''))),
   
    number_of_genes = (zero_gene_location_st_tt[[j]][1])$number_of_gene ,
    number_of_location = (zero_gene_location_st_tt[[j]][2])$number_of_location,
    zero_percentage =round(100*(zero_gene_location_st_tt[[j]][3]$percentage_of_zero),2),
    
    source_name =str_to_title(str_remove(list.files(paste(getwd(),"/","st/",second_dir_st[j],"/paper",sep = ''))[1],'.pdf')) ,
    source_link =link_st[j],
    image_file = image(paste(getwd(),"/","st/",second_dir_st[j],"/pathology_image",sep = ''))
    #dataset_code = paste0("st/",second_dir_st[j],"/code/",list.files(paste(getwd(),"/","st/",second_dir_st[j],"/code",sep = '')))
  )
}

# JSON for hdst datasets
hdst_number_of_genes <- matrix(c(22133,19068,22732,16174,19950,19407), ncol = 2,byrow = FALSE)
hdst_number_of_location <-  matrix(c(178800,130996,177707,118551,181367,176840), ncol = 2,byrow = FALSE)
hdst_zero_percent <- matrix(c(99.987,99.985,99.983,99.981,99.960,99.969), ncol = 2, byrow = FALSE)
disease_hdst = c("breast cancer", NA)
for(j in (1):(length(second_dir_hdst))){
  le = le + length(second_dir_st) 
  pre_json[[le+j]]<- list(
    dataset_name =second_dir_hdst[j], 
    display_name = gsub("_",' ',str_to_title(second_dir_hdst[j])),
    dataset_type = "Real Data", 
    technology = "HDST",#"Hign Definition Spatial Transcriptomic", 
    disease = str_to_title(disease_hdst[j]), 
    organism = hdst_organism[j] ,
    tissue = hdst_tissue[j],
  #  dataset_description = "Author: Jie Yang (jxy190005@utdallas.edu)",
    raw_data = paste0("hdst/",second_dir_hdst[j],"/raw_data/",list.files(paste(getwd(),"/","hdst/",second_dir_hdst[j],"/raw_data",sep = ''), pattern = '.tsv')),
   
    processed_data = paste0("hdst/",second_dir_hdst[j],"/processed_data/",list.files(paste(getwd(),"/","hdst/",second_dir_hdst[j],"/processed_data",sep = ''))),
   # normalized_expression = paste0("hdst/",second_dir_hdst[i],"/normalized_data/",list.files(paste(getwd(),"/","hdst/",second_dir_hdst[i],"/normalized_data",sep = ''))),
 
    number_of_genes = hdst_number_of_genes[,j] ,
    number_of_location =hdst_number_of_location[,j] ,
    zero_percentage = hdst_zero_percent[,j], 
    source_name =str_to_title("vickovic_2019") ,
    source_link = 'https://doi.org/10.1038/s41592-019-0548-y' ,
    image_file = paste0("hdst/",second_dir_hdst[j],"/pathology_image/",list.files(paste(getwd(),"/","hdst/",second_dir_hdst[j],"/pathology_image",sep = '')))
   # dataset_code = paste0("hdst/",second_dir_hdst[j],"/code/",list.files(paste(getwd(),"/","hdst/",second_dir_hdst[j],"/code",sep = '')))
  )
}

# JSON for 10x datasets
link_10x <- list()
link_10x[[1]] =  c("https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Human_Brain_Section_1",
                   "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Human_Brain_Section_2")
link_10x[[2]] = "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.3.0/Visium_FFPE_Human_Breast_Cancer" 
link_10x[[3]] =    c("https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Breast_Cancer_Block_A_Section_1",
                     "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Breast_Cancer_Block_A_Section_2") 
link_10x[[4]] = "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.2.0/Parent_Visium_Human_BreastCancer"
link_10x[[5]] = "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.2.0/Parent_Visium_Human_Cerebellum"
link_10x[[6]] = "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.2.0/Parent_Visium_Human_ColorectalCancer"
link_10x[[7]] = "https://doi.org/10.1038/s41593-020-00787-0"
link_10x[[8]] ="https://support.10xgenomics.com/spatial-gene-expression/datasets/1.2.0/Parent_Visium_Human_Glioblastoma"
link_10x[[9]] = "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Human_Heart"
link_10x[[10]]=  "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.2.0/V1_Human_Invasive_Ductal_Carcinoma"
link_10x[[11]]="https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Human_Lymph_Node"
link_10x[[12]]="https://support.10xgenomics.com/spatial-gene-expression/datasets/1.3.0/Visium_FFPE_Human_Normal_Prostate"
link_10x[[13]]= "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.2.0/Parent_Visium_Human_OvarianCancer"
link_10x[[14]]= "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.3.0/Visium_FFPE_Human_Prostate_Cancer"
link_10x[[15]]= "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.2.0/Parent_Visium_Human_SpinalCord"
link_10x[[16]]= "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain"
link_10x[[17]]= c("https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain_Coronal_Section_1",
                  "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain_Coronal_Section_2")
link_10x[[18]]= c("https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Mouse_Brain_Sagittal_Anterior",
                 "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Mouse_Brain_Sagittal_Anterior_Section_2",
                 "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Mouse_Brain_Sagittal_Posterior",
                 "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Mouse_Brain_Sagittal_Posterior_Section_2")
link_10x[[19]]= "https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Mouse_Kidney"

disease_tenx = c(NA,"breast cancer", "brease cancer","brease cancer",
                 NA,"colorectal cancer","glioblastoma",NA,
                 "breast cancer-invasive ductal carcinoma",NA,"normal","ovarian cancer",
                 "prostate cancer",
                 NA,NA,NA,NA,NA)
display_name_tenx <-c("Adult Human Brain (Cerebral Cortex)", "Human Breast Cancer: Ductal Carcinoma In Situ, Invasive Carcinoma (FFPE)","Human Breast Cancer (Block A)",
                      "Human Breast Cancer: Whole Transcriptome Analysis","Human Cerebellum: Whole Transcriptome Analysis","Human Colorectal Cancer: Whole Transcriptome Analysis",
                      "Human Dorsolateral Prefrontal Cortex",
                      "Human Glioblastoma: Whole Transcriptome Analysis", "Human Heart", "Invasive Ductal Carcinoma Stained With Fluorescent CD3 Antibody",
                      "Human Lymph Node","Normal Human Prostate (FFPE)", "Human Ovarian Cancer: Whole Transcriptome Analysis. Stains: DAPI, Anti-PanCK, Anti-CD45",
                      "Human Prostate Cancer, Adenocarcinoma with Invasive Carcinoma (FFPE)","Human Spinal Cord: Whole Transcriptome Analysis. Stains: DAPI, Anti-SNAP25, Anti-GFAP, Anti-Myelin CNPase",
                      "Mouse Brain Section (Coronal)", "Adult Mouse Brain (Coronal). Stains: DAPI, Anti-NeuN", "Mouse Brain Serial Section (Sagitaal)","Mouse Kidney Section (Coronal)")
for(j in (1):(length(second_dir_tenx))){
  le = le + length(second_dir_hdst)
  pre_json[[le+j]]<- list(
    dataset_name =second_dir_tenx[j], 
    display_name = paste(display_name_tenx[j],"(10x)"),
    dataset_type = "Real Data",
    technology = "ST",#"Spatial Transcriptomic", 
    disease = str_to_title(disease_tenx[j]),
    organism = tenx_organism[j] ,
    tissue = tenx_tissue[j],
  #  dataset_description = "Author: Jie Yang (jxy190005@utdallas.edu)",
    raw_data = paste0("10x/",second_dir_tenx[j],"/raw_data/",list.files(paste(getwd(),"/","10x/",second_dir_tenx[j],"/raw_data",sep = ''))),
    
    processed_data = paste0("10x/",second_dir_tenx[j],"/processed_data/",list.files(paste(getwd(),"/","10x/",second_dir_tenx[j],"/processed_data",sep = ''))),
    # normalized_expression = ifelse(dir.exists(paste(getwd(),"/","10x/",second_dir_tenx[i],"/normalized_data",sep = '')),
    #                                paste0("10x/",second_dir_tenx[i],"/normalized_data/",list.files(paste(getwd(),"/","10x/",second_dir_tenx[i],"/normalized_data",sep = ''))),
    #                                NULL),
    
    number_of_genes = as.vector((zero_gene_location_tenx_tt[[j]][1])$number_of_gene),
    number_of_location = as.vector((zero_gene_location_tenx_tt[[j]][2])$number_of_location),
    zero_percentage =as.vector(round(100*(zero_gene_location_tenx_tt[[j]][3]$percentage_of_zero),2)),
    source_name = "10x Genomics Spatial Gene Expression Datasets",
    source_link = link_10x[[j]],  
    image_file = paste0("10x/",second_dir_tenx[j],"/pathology_image/",list.files(paste(getwd(),"/","10x/",second_dir_tenx[j],"/pathology_image",sep = '')))
    #dataset_code = paste0("10x/",second_dir_tenx[j],"/code/",list.files(paste(getwd(),"/","10x/",second_dir_tenx[j],"/code",sep = '')))
  )
}

#--------------------------------------------------------------------------------------------
# JSON UPDATE - MORE DATASETS
#--------------------------------------------------------------------------------------------

# JSON for slide-SeqV2
slide_number_of_genes <-list(list())
slide_number_of_genes[[1]] = c(23124,	46473	,23124,	23117,19600)
slide_number_of_genes[[2]] = c(22683,19653)
slide_number_of_genes[[3]] = c(23264,22457)
slide_number_of_genes[[4]] = 21220
slide_number_of_genes[[5]] = 22542
#slide_number_of_location <-  c(51649,	23055,	51649,	56047,	33611,	20143,	37844,53208	,34199,	21724	,42550)
slide_number_of_location <-list(list())
slide_number_of_location[[1]] <- c(51649,	23055,	51649,	56047,37844)
slide_number_of_location[[2]] <- c(33611,	20143)
slide_number_of_location[[3]] <- c(53208	,34199)
slide_number_of_location[[4]] <- 21724
slide_number_of_location[[5]] <- 42550
slide_zero_percent = list(list())
slide_zero_percent[[1]] = c(98.05	,98.22,	98.05,	99.23,99.69)
slide_zero_percent[[2]] = c(	98.55,	99.05)
slide_zero_percent[[3]] = c(98.19,	98.49	)
slide_zero_percent[[4]] =98.38
slide_zero_percent[[5]] = 98.47
# generate temporary  file and then put together to avoid long running time
#pre_json = list(list())
for(j in (1):(length(second_dir_slideseqv2))){
  le = le + length(second_dir_tenx)
  pre_json[[le+j]]<- list(
    dataset_name =second_dir_slideseqv2[j], 
    display_name = gsub("_",' ',str_to_title(second_dir_slideseqv2[j])),
    dataset_type = "Real Data", 
    technology = "Slide-seqV2",
    disease = NA, 
    organism =slideseqv2_organism,
    tissue = slideseqv2_tissue[j],
    #  dataset_description = "Author: Jie Yang (jxy190005@utdallas.edu)",
    raw_data = paste0("slide_seqv2/",second_dir_slideseqv2[j],"/raw_data/",list.files(paste(getwd(),"/","slide_seqv2/",second_dir_slideseqv2[j],"/raw_data",sep = ''))),
    
    processed_data = paste0("slide_seqv2/",second_dir_slideseqv2[j],"/processed_data/",list.files(paste(getwd(),"/","slide_seqv2/",second_dir_slideseqv2[j],"/processed_data",sep = ''))),
    # normalized_expression = paste0("slide/",second_dir_slide[i],"/normalized_data/",list.files(paste(getwd(),"/","slide/",second_dir_slide[i],"/normalized_data",sep = ''))),
    
    number_of_genes = slide_number_of_genes[[j]] ,
    number_of_location =slide_number_of_location[[j]] ,
    zero_percentage = slide_zero_percent[[j]], 
    source_name =str_to_title("stickels_2021") ,
    source_link = 'https://doi.org/10.1038/s41587-020-0739-1' ,
    image_file = NULL
    # dataset_code = paste0("slide/",second_dir_slide[j],"/code/",list.files(paste(getwd(),"/","slide/",second_dir_slide[j],"/code",sep = '')))
  )
}

# JSON for other 
other_number_of_genes <-list(list())
other_number_of_genes[[1]] = c( 19738,    19738,    19738,    19738,    19738,    19738,    19738,    19738,    19738,     5758)
other_number_of_location <-list(list())
other_number_of_location[[1]] = c(428,      288,      316,      288,      224,      321,      242,      460,      359,      359)
other_zero_percent = list(list())
other_zero_percent[[1]] = c( 95.16,    94.61,    95.46,    97.84,    94.09,    97.31,    96.36,    91.25,    96.10,    90.96)
for(j in (1):(length(second_dir_other))){
  le = le + length(second_dir_slideseqv2)
  pre_json[[le+j]]<- list(
    dataset_name =second_dir_other[j], 
    display_name = gsub("_",' ',str_to_title(second_dir_other[j])),
    dataset_type = "Real Data", 
    technology = "Other",
    disease = NA, 
    organism =other_organism[j],
    tissue = other_tissue[j],
    #  dataset_description = "Author: Jie Yang (jxy190005@utdallas.edu)",
    raw_data = paste0("microarray_singlecell/",second_dir_other[j],"/raw_data/",list.files(paste(getwd(),"/","microarray_singlecell/",second_dir_other[j],"/raw_data",sep = ''))),
    
    processed_data = paste0("microarray_singlecell/",second_dir_other[j],"/processed_data/",list.files(paste(getwd(),"/","microarray_singlecell/",second_dir_other[j],"/processed_data",sep = ''))),
    # normalized_expression = paste0("slide/",second_dir_slide[i],"/normalized_data/",list.files(paste(getwd(),"/","slide/",second_dir_slide[i],"/normalized_data",sep = ''))),
    
    number_of_genes = other_number_of_genes[[j]] ,
    number_of_location =other_number_of_location[[j]] ,
    zero_percentage = other_zero_percent[[j]], 
    source_name =str_to_title("Moncada_2020") ,
    source_link = 'https://doi.org/10.1038/s41587-019-0392-8' ,
    image_file = NULL
    # dataset_code = paste0("slide/",second_dir_slide[j],"/code/",list.files(paste(getwd(),"/","slide/",second_dir_slide[j],"/code",sep = '')))
  )
}

# JSON for starmap

zero_gene_location_starmap = list(list())
for(k in 1:length(second_dir_starmap)){
  processed_data = list.files(paste(getwd(),"/","starmap/",second_dir_starmap[k],"/processed_data",sep = ''),pattern = ".Rdata")
  #print(processed_data)
  para_starmap = matrix(nrow = 3, ncol = length(processed_data))
  row.names(para_starmap) = c("number_of_genes","number_of_location","percentage_of_zero")
  for(j in 1:length(processed_data)){
    Path = paste(getwd(),"/","starmap/",second_dir_starmap[k],"/processed_data","/",processed_data[j],sep = '')
    tt = zero_propor_and_gene_location(Path)
    para_starmap[1,j] = tt$number_of_gene
    para_starmap[2,j] = tt$number_of_location
    para_starmap[3,j] = tt$percentage_of_zero
  }
  zero_gene_location_starmap[[k]] = list(number_of_gene = para_starmap[1,],number_of_location = para_starmap[2,],percentage_of_zero = (para_starmap[3,]))
}
for(j in (1):(length(second_dir_starmap))){
  le = le + length(second_dir_other)
  pre_json[[le+j]]<- list(
    dataset_name =second_dir_starmap[j], 
    display_name = gsub("_",' ',str_to_title(second_dir_starmap[j])),
    dataset_type = "Real Data", 
    technology = "STARmap",
    disease = NA, 
    organism =starmap_organism[j],
    tissue = starmap_tissue[j], 
    #  dataset_description = "Author: Jie Yang (jxy190005@utdallas.edu)",
    raw_data = paste0("starmap/",second_dir_starmap[j],"/raw_data/",list.files(paste(getwd(),"/","starmap/",second_dir_starmap[j],"/raw_data",sep = ''))),
    
    processed_data = paste0("starmap/",second_dir_starmap[j],"/processed_data/",list.files(paste(getwd(),"/","starmap/",second_dir_starmap[j],"/processed_data",sep = ''))),
    # normalized_expression = paste0("slide/",second_dir_slide[i],"/normalized_data/",list.files(paste(getwd(),"/","slide/",second_dir_slide[i],"/normalized_data",sep = ''))),
    
    number_of_genes = (zero_gene_location_starmap[[j]][1])$number_of_gene ,
    number_of_location = (zero_gene_location_starmap[[j]][2])$number_of_location ,
    zero_percentage = round(100*(zero_gene_location_starmap[[j]][3]$percentage_of_zero),2),
    source_name =str_to_title("Wang_2018") ,
    source_link = 'https://doi.org/10.1126/science.aat5691' ,
    image_file = NULL
    # dataset_code = paste0("slide/",second_dir_slide[j],"/code/",list.files(paste(getwd(),"/","slide/",second_dir_slide[j],"/code",sep = '')))
  )
}

tt<- toJSON(pre_json)

write(tt, file = "real_data.json")




#--------------------------------------------------------------------------------
# simulation data
#--------------------------------------------------------------------------------

# zero_setting

zero_real <- read.csv("C:/Users/yojia/Box/Data/zero_setting.csv")[,-1]

simulation_dataset_name = c("bc_pattern","linear_pattern","mob_i_pattern","mob_ii_pattern","spot_pattern")
patterns = c("Breast Cancer Pattern", "Linear Pattern" ,
             "Mouse Olfactory Bulb i Pattern","Mouse Olfactory Bulb ii Pattern","Spot Pattern")
zero_setting<- c("0","10","30","50")
display_name = c("Breast cancer pattern","Linear pattern","Mouse olfactory bulb i pattern","Mouse olfactory bulb ii pattern","Spot pattern")

data_discription =c("Sample_size:n=250, Distribution: ZeroInflatedNegativeBinomial",
                    "Sample_size:n=256, Distribution: ZeroInflatedNegativeBinomial",
                    "Sample_size:n=260, Distribution: ZeroInflatedNegativeBinomial",
                    "Sample_size:n=260, Distribution: ZeroInflatedNegativeBinomial",
                    "Sample_size:n=256, Distribution: ZeroInflatedNegativeBinomial" ) 

number_of_genes = 100
number_of_location = c(250,256,260,260,256)




pre_json <- list(list())
for(s in 1:5){
  for(j in 1:4){
    pre_json[[(s-1)*4+j]] <- list(
    pattern = patterns[s],
    dataset_name =paste0(simulation_dataset_name[s],"_",zero_setting[j]), 
    display_name = paste(display_name[s],zero_setting[j]),
    dataset_type = "Simulation Data",
  #  dataset_description = data_discription[s],
    processed_data = paste0("simulated_data/",simulation_dataset_name[s],"/data/",
                      list.files(paste("C:/Users/yojia/Box/Data/simulated_data/",simulation_dataset_name[s],"/data",sep = ""))[(30*(j-1)+1):(30*j)]),
    image_file =paste0("simulated_data/",simulation_dataset_name[s],"/images/", 
                       list.files(paste("C:/Users/yojia/Box/Data/simulated_data/",simulation_dataset_name[s],"/images",sep = ""))),
    # dataset_code = paste0("simulated_data/",simulation_dataset_name[s],"/code/",
    #                       list.files(paste("C:/Users/yojia/Box/Data/simulated_data/",simulation_dataset_name[s],"/code",sep = ""),
    #                                  pattern = ".R|.Rdata")),
    source_name = 'Li_2021',
    source_link = "https://doi.org/10.1093/bioinformatics/btab455",
    number_of_genes = rep(100,30) ,
    number_of_location = rep(number_of_location[s],30) ,
    zero_percentage = 100*zero_real[((30*(j-1)+1):(30*j)),s]
    
  )
  }
  
}

tt<- toJSON(pre_json)

write(tt, file = "simulated_data.json")
