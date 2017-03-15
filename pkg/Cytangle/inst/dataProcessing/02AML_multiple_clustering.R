# AML Clustering Script

# Created before 5/26/16
# Updated on 9/16/16

# Set working directory (and location of .fcs data files), load libraries

.libPaths("~/Software/RLib")

#setwd("./Desktop/CyTOF_Analysis/code/experiments") 
if (!grepl(paste0(getwd(),"$"),"/home/mcgee.278/Dropbox/Research/CyTOF_Analysis/code/experiments",fixed = TRUE)){
  setwd("~/Dropbox/Research/CyTOF_Analysis/code/experiments")
}

# Load AML Info
load("../results/AML_info.rda") 

folder <- paste(root,"/CyTOF_Analysis/data/",dataset,"/FCS_edit", sep="") 
outFolder <- paste(root,"/CyTOF_Analysis/data/",dataset,"/spade",run_date, sep="")

#library("tools") # check that these actually aren't necessary anymore
#library("base")
library("flowCore") #colnames(data_obj) will be null on line 47 if this isn't loaded and you don't use
#library("class")
library("spade")

# Collect data files
allfiles <- list.files(path = folder, pattern = ".fcs")

# Non-MDS and Stimulated files
MDS_indices <- grep('MDS|(GCS|SC)F|Tm',allfiles)
if (length(MDS_indices) == 0){
  files <- allfiles
}else{
  files <- allfiles[-MDS_indices]
}

# gather file paths

src <- file.path(folder, files)
src1 <- file.path(folder, files[111])
data_obj <- SPADE.read.FCS(src1, transformation=FALSE) #08/17 check that this works, or load flowCore

# identify surface markers for clustering

context_vars <- colnames(data_obj)[surfInd]
rm(data_obj)

# Spade

#num_clusts <- 483 #floor(sqrt(dim(data_obj)[[1]]))

data_driver <- SPADE.driver(src, out_dir = outFolder, cluster_cols = context_vars,
                            transforms=flowCore::arcsinhTransform(a=0,b=0.2), k = num_clusts)

# Plotting

plot_path <- paste0(outFolder,"/plots")
if (!dir.exists(plot_path)){
  dir.create(plot_path)
} 

mst_graph <- igraph:::read.graph(paste(outFolder,"mst.gml",sep=.Platform$file.sep),format="gml")

layout <- read.table(paste(outFolder,"layout.table",sep=.Platform$file.sep))

SPADE.plot.trees(mst_graph, outFolder, file_pattern = "*anno.Rsave", out_dir = plot_path, layout = igraph::layout.kamada.kawai,
                 attr_pattern = "percent|medians|fold|cvs", scale = NULL, pctile_color=c(0.02,0.98), normalize="global", size_scale_factor=1,
                 edge.color="grey", bare=FALSE, palette="bluered")


