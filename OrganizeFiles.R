#Author: Mehrnoush Malek
#Last updated: October 2018

#This code needs to be run before the automated analysis, and after file download,
# to move the files in proper folders and formatting

library(readr)
flowrepos.path<-"!/flowRepositoryExperimentDownloaded/"
desired.path<-"~/PathToSaveFilesintoFolders"
#This csv files is insdie the flowRepository Experiment
meta.data<- read.csv("~PathTo/MetaData for FlowRepository uploads.csv")
unique(meta.data$Folder)
grbg<- sapply(unique(meta.data$Folder),function(x){
  tmp<- meta.data[meta.data$Folder==x,]
  types<-unique(meta.data[which(meta.data$Folder==x),"Filetype"])
  grbg<- sapply(types[types!="sample"], function(y) {
    samples<-meta.data[which(meta.data$Folder==x & meta.data$Filetype==y),"Filename"]
    suppressWarnings(dir.create(paste(desired.path,x,y,sep="/"),recursive=T))
    grbg<-file.copy(from=paste0(flowrepos.path,samples),to = paste(desired.path,x,y,samples,sep="/"))
  })
  sub.type <- unique(meta.data[which(meta.data$Folder==x & meta.data$Filetype=="sample"),"Sampletype"])
  grbg<- sapply(sub.type, function(z) {
    samples<-meta.data[which(meta.data$Folder==x & meta.data$Filetype=="sample" & meta.data$Sampletype==z),"Filename"]
    
  suppressWarnings(dir.create(paste(desired.path,x,"sample",z,sep="/"),recursive=T))
    grbg<-file.copy(from=paste0(flowrepos.path,samples),to = paste(desired.path,x,"sample",z,samples,sep="/"))
  })
})
       
                
