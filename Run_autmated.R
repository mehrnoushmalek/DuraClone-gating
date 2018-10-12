#Author: Mehrnoush Malek
#Last updated: October 2018
 

#Panel, dataset, and type needs to be changed manually to the dataset you wish to analyze.
#MAke sure all the path corresponds to your local path, where the files and codes are.
rm(list=ls())
library(flowCore)
library(flowDensity)
library(tcltk)
#library(pryr)
###This needs to be modified
source('~/PathToDuraClone-gatingFiles/HelperFunctions.R', echo=F)

#Choose which panel, type and centre should be analyzed
panel <- c("Basic","B cell")[1]
dataset <-c("whole blood","Fortessa vs. Navios")[2]
#For whole blood the type is either healthy or HSCT, and for Fortessa vs. Navios, either Navios or Fortessa
#Make sure the typpe corresponds to the propert dataset, and the infromation above
type <-c("healthy volunteer","post HSCT","Navios","Fortessa")[3]

#change this to directory where all the files are
###This needs to be modified
dir.path <- "~/PathToFlowRepositoryDownloadedFolder"
#Change this to a directory where you want the files to be saved. This would be the parent folder.
#Te analysis for each data set will be saved under their own name in this folder.
###This needs to be modified
output.path <- "~/PathToSaveTheResults"
#Files are in Bioinformatics folder
path <- paste(dir.path,dataset,"sample",type,sep="/")

dir.create(paste(output.path,dataset,"Results",type,sep="/"),recursive=T)
res.path <- paste(output.path,dataset,"Results",type,sep="/")
#path to FCS files
filenames<- list.files(path,pattern=".LMD|.fcs",recursive = T)
#This is the metadata in flowRepository Experiment
meta.data<- read.csv("~/MetaData for FlowRepository uploads.csv")
sub.inds<- which(meta.data$Panel[match(filenames,meta.data$Filename)]==panel)
fcs.path <- list.files(path,full.names=T,pattern=".LMD|.fcs",recursive = T)[sub.inds]
#Make sure all codes are within the same folder, otheriwise proper calling is needed.
source(paste0("~/PathToDuraClone-gating/",panel, "-Gating.R"))

