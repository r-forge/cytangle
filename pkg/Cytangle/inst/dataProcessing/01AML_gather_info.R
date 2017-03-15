# AML Marker Check Script
# Created 03/11/16 by R. McGee
# Last update: 10/4/16 by R. McGee
# Takes all .fcs in specified folder and converts them to .csv files.
# Inspired by: https://support.bioconductor.org/p/32759/

# Set working directory (and location of .fcs data files), load libraries

.libPaths("~/Software/RLib")

#setwd("./Desktop/CyTOF_Analysis/code/experiments") 
if (!grepl(paste0(getwd(),"$"),"/home/mcgee.278/Dropbox/Research/CyTOF_Analysis/code/experiments",fixed = TRUE)){
  setwd("~/Dropbox/Research/CyTOF_Analysis/code/experiments")
}

###############################################################################################################
###############################################################################################################
run_date <- 091216
root <- '~/Desktop'  #'/scratch/mcgee.278' #  '../..'
dataset <-'AML_experiment_51184' # Nature_viSNE_data
###############################################################################################################
###############################################################################################################

folder <- paste(root,"/CyTOF_Analysis/data/",dataset,"/FCS_edit", sep="")
outFolder <- paste0("../results/")  #signaling, cell_cycle

library("tools")
library("base")
library("flowCore")

# Collect data files
allfiles <- list.files(path = folder, pattern = ".fcs")

# Non-MDS and Stimulated files
MDS_indices <- grep('MDS|(GCS|SC)F|Tm',allfiles)
if (length(MDS_indices) == 0){
  files <- allfiles
}else{
  files <- allfiles[-MDS_indices]
}

# Conversion loop

#AML_indices <- c(1,11,14,20,35,56,139:141,152)
all_markers <- character(0) 
#common_markers <- character(0)
pot_surf <- character(0)
tubeA_markers <- character(0)
tubeB_markers <- character(0)

barcodes <- c("Pd102","Pd104","Pd105","Pd106","Pd108","Pd110")
others <- c("Uridine","Ce140","CD56",'DNA-1')

groups <- array(0,c(length(files),4))
colnames(groups) <- c("Identifier","Group Number","AML Subgroup","Tube")

num_clusts <- 483
num_grps <- 10

for (idx in seq_along(files)){

  myFrame <- read.FCS(file.path(folder, files[idx]), transformation=FALSE)
  base <- file_path_sans_ext(files[idx])
  short <- gsub('(.GCSF|.SCF|)_M[A-z0-9]*.',"",gsub('T(A|B)-AMLx(0|1|2|3|4)_(Tm|)',"",base))
  
  # initialize ion lookup array
  if (idx == 1){
    ions <- colnames(myFrame)[3:length(colnames(myFrame))]
    ion_lookup <- array(0,c(length(ions),3))
    colnames(ion_lookup) <- c("Ion","Tube A","Tube B")
    ion_lookup[,1] <- ions
  }
  
  # save tube information
  groups[idx,4] <- substr(base,2,2)
  
  # key for naming convention
  if (grepl(short,"AML34")){
    groups[idx,1] <- "APL5"
  }else if(grepl(short,"AML24")){
    groups[idx,1] <- "APL4"
  }else if(grepl(short,"APL2")){
    groups[idx,1] <- "AML39"
  }else if(grepl(short,"Rec1")){
    groups[idx,1] <- "AML40"
  }else if(grepl(short,"Rec4")){
    groups[idx,1] <- "AML41"
  }else if(grepl(short,"AML16")){
    groups[idx,1] <- "CR6"
  }else if(grepl(short,"Rec3")){
    groups[idx,1] <- "CR3"
  }else if(grepl(short,"AML11")){
    groups[idx,1] <- "CR5"
  }else if(grepl(short,"AML17")){
    groups[idx,1] <- paste0('CR7','_x',substr(base,8,8)) #more than one file with this short name
  }else if(grepl(short,"Rec2")){
    groups[idx,1] <- "CR2"
  }else if(grepl(short,'Nl-3|Nl-4|Nl-5|NL-6|CR7|AML25|AML8|Nl-3|Nl-4|Nl-5|NL-6')) 
    groups[idx,1] <- paste0(short,'_x',substr(base,8,8)) #more than one file with this short name
  else{
    groups[idx,1] <- short
  } 
  
  if (grepl('AML26|AML35|AML5|AML10|AML32',groups[idx,1],ignore.case = TRUE)){ #(groups[idx,1] %in% c('AML26','AML35','AML5','AML10','AML32')){
    groups[idx,2] <- 1
    groups[idx,3] <- 'CBF-AML'
  }else if (grepl('APL',groups[idx,1],ignore.case = TRUE)){
    groups[idx,2] <- 2
    groups[idx,3] <- 'APL'
  }else if (grepl('AML27|AML7|AML9|AML20|AML21|AML23|AML30|AML39|AML37|AML38|AML40',groups[idx,1],ignore.case = TRUE)){ #(groups[idx,1] %in% c('AML27','AML7','AML9','AML20','AML21','AML23','AML30','AML39','AML37','AML38','AML40')){
    groups[idx,2] <- 3
    groups[idx,3] <- 'NK-AML FLT3-ITD'
  }else if (grepl('AML13|AML36|AML33|AML42|AML14|AML15',groups[idx,1],ignore.case = TRUE)){ # (groups[idx,1] %in% c('AML13','AML36','AML33','AML42','AML14','AML15')){
    groups[idx,2] <- 4
    groups[idx,3] <- 'NK-AML FLT3wt'
  }else if(grepl('AML18|AML4|AML19|AML29|AML41|AML25',groups[idx,1],ignore.case = TRUE)){ #(groups[idx,1] %in% c('AML18','AML4','AML19','AML29','AML41','AML25')){
    groups[idx,2] <- 5
    groups[idx,3] <- 'Adverse-risk'
  }else if(grepl('AML8|AML6|AML31',groups[idx,1],ignore.case = TRUE)){ # (groups[idx,1] %in% c('AML8','AML6','AML31')){
    groups[idx,2] <- 6
    groups[idx,3] <- 'Normal Karotype AML FLT3 tyrosine kinase domain mutation'
  }else if(grepl('AML12',groups[idx,1],ignore.case = TRUE)){ # %in% c('AML12')){
    groups[idx,2] <- 7
    groups[idx,3] <- 'Granulylocytic sarcoma'
  }else if(grepl('AML22',groups[idx,1],ignore.case = TRUE)){ #(groups[idx,1] %in% c('AML22')){
    groups[idx,2] <- 8
    groups[idx,3] <- 'MLL rearrangement'
  }else if(grepl('CR',groups[idx,1],ignore.case = TRUE)){
    groups[idx,2] <- 9
    groups[idx,3] <- 'Complete Recovery'
  }else if(grepl('Nl',groups[idx,1],ignore.case = TRUE)){
    groups[idx,2] <- 10
    groups[idx,3] <- 'Normal'
  }
  
  # grab markers
  temp <- markernames(myFrame)
  
  sharedInd <- grepl('-',temp) & !(temp %in% c('DNA-1','DNA-2','Ki-67','HLA-DR')) # Logical
  
  # strsplit('-',temp)
  if (grepl('TA',base)){  #If it's Tube A
    temp[sharedInd] <- gsub('-[A-z0-9]*',"",temp[sharedInd]) 
    tubeA_markers <- unique(c(tubeA_markers,temp))
  }else{ # If it's Tube B
    temp[sharedInd] <- gsub(".*-","",temp[sharedInd])
    tubeB_markers <- unique(c(tubeB_markers,temp))
  }
  # take barcodes and others out here?
  
  if(idx == 1){
    tubeA_panel <- temp
  }else if(idx == length(files)){
    tubeB_panel <- temp
  }
  
  # Marker check
  if (length(all_markers)==0){
    all_markers <- temp
    #common_markers <- temp
    pot_surf <- temp[!sharedInd]
  }else{
    tempInd <- temp %in% all_markers 
    all_markers <- c(all_markers,temp[!tempInd])
    #commInd <- temp %in% common_markers
    #common_markers <- temp[commInd]
    pot_surf <- unique(c(pot_surf,temp[!sharedInd]))
  }
  
  #write.table(all_markers, markerFile, sep="\t", quote=FALSE, append=TRUE)
  #rm(temp,tempInd,markerFile)
}
tubeA_markers <- tubeA_markers[!tubeA_markers %in% c(barcodes,others)]
tubeB_markers <- tubeB_markers[!tubeB_markers %in% c(barcodes,others)]

# Remove extraneous 
extraKeys1 <- grep('empty|File|bead|barcode',all_markers)
aml_markers <- all_markers[-extraKeys1]

# Surface markers
surfaceInd <- grep('CD([0-4,6-8])|HLA-DR',pot_surf)  # Change this so that it does non-hypen
surface_markers <- pot_surf[surfaceInd]

# Common Markers
#extraKeys2 <- grep('empty|File|bead|barcode',common_markers)
#common_markers <- common_markers[-extraKeys2]
common_markers2 <- intersect(tubeA_markers,tubeB_markers) # common_markers and common_markers2 SHOULD be equal
extraKeys3 <- grep('empty|File|bead|barcode',common_markers2)
common_markers2 <- common_markers2[-extraKeys3]
#extraKeys4 <- grep(c(others,barcodes),common_markers2)
#common_markers2 <- common_markers2[-extraKeys4]
commonsurfInd <- common_markers2 %in% surface_markers
common_nonsurface_markers <- common_markers2[!commonsurfInd]
#common_markers <- common_markers2
common_markers <- common_markers2

# Markers unique to each tube
extraKeysAB <- grep('empty|File|bead|barcode',tubeA_markers)
tubeA_markers <- tubeA_markers[-extraKeysAB]
tubeB_markers <- tubeB_markers[-extraKeysAB]
tubeA_featInd <- !(tubeA_markers %in% c(common_markers,barcodes,others)) #| (grepl("pS6",tubeA_markers))
tubeA_features <- tubeA_markers[tubeA_featInd] 
tubeB_featInd <- !(tubeB_markers %in% c(common_markers,barcodes,others)) #| (grepl("pS6",tubeB_markers))
tubeB_features <- tubeB_markers[tubeB_featInd]

#find surfInd and featInd
tempInd <- temp %in% surface_markers
#change #edits have an added entry
surfInd <- c(1:length(temp))[tempInd]+2 # Cell_length is included in markernames, makes sure this isn't the case with all files
tempInd2 <- temp %in% c(tubeB_features,common_nonsurface_markers) #temp is a tube B file at this point
featInd <- c(1:length(temp))[tempInd2]+2

# finish ion lookup array

ion_lookup[,2:3] <- c(tubeA_panel,tubeB_panel)

# save 
marker_info <- paste(outFolder,"/AML_info.rda",sep="")
save(run_date,root,dataset,groups,num_clusts,num_grps,surfInd,featInd,aml_markers,common_markers,common_nonsurface_markers,surface_markers,tubeA_panel,tubeB_panel,tubeA_features,tubeB_features,barcodes,others,ion_lookup,file=marker_info)
