rm(list=ls())
library(flowCore)
library(flowDensity)
library(tcltk)
#library(pryr)
source('~/Project-SVN/OneStudy/HelperFunctions.R', echo=F)

#Choose which panel, type and centre should be analyzed
panel <- c("ONE_1","ONE_5")[2]
#type is "" for all datas except the CNTRP that has whole blood samples and PBMCs
type <-c("","WBC")[2]
centre <-c("CNTRP","CNTRP-PhaseII","Navios","Fortessa")[2]

#change this to directory where all the files are
dir.path <- "/mnt/f/FCS data/OneStudy/"
#Change this to a directory where you want the files to be saved
res.path <- "/mnt/data/One-Study/"
#Files are in Bioinformatics folder
path <- paste(dir.path,centre,panel,sep="/")
if (centre=="CNTRP")
  path <- paste(path,type,sep="/")

dir.create(paste(res.path,centre,"Results",type,panel,sep="/"),recursive=T)
res <- paste("/mnt/data/One-Study",centre,"Results",type,panel,sep="/")
#path to FCS files
fcs.path <- list.files(path,full.names=T,pattern=".LMD|.fcs",recursive = T)
source(paste0("~/Project-SVN/OneStudy/",panel, "-Gating.R"))

