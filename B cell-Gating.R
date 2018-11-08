#Panel 5
#Author: Mehrnoush Malek
#Date: Last revised in October 2018



#Read the second dataset in the LMD file as it's not transformed or compensated
#In here we only keep the markers names from the first data set in order to use it for the second dataset within each FCS file
#Reading files
#------------------------------------------------------------------------------------

f <- read.FCS(fcs.path[1],dataset = 1)  
markers <- unlist(lapply(as.vector(f@parameters@data[,2]), function(x) unlist(strsplit(x," "))[1]))
markers[grep(x=markers,pattern = "CD*")]<-toupper(markers[grep(x=markers,pattern = "CD*")])
if (type!="Fortessa")
{
  markers <- markers[-c(which(markers=="FS"),which(markers=="SS"))]
  
  fs <- tryCatch(read.flowSet(fcs.path,dataset = 2),error=function(x) {return(1)})
  if (mode(fs)=="numeric")
  {
    fs <- as(lapply(fcs.path,function(x) {
      f<-read.FCS(x,dataset=2)
      f <- f[,c("FS-A","FS-W","SS-A","FL1-A","FL2-A","FL3-A","FL5-A","FL6-A","FL7-A","FL8-A","FL9-A","FL10-A")]
    }), Class="flowSet")
  }
  fs<- fsApply(fs, function(x) {
    x@parameters@data[,2] <- c("FSC-A","FSC-W","SSC-A",markers);
    return(x)
  })
  
}else{
  fs <- tryCatch(read.flowSet(fcs.path,dataset = 2),error=function(x) {return(1)})
  markers <- markers[-which(is.na(markers))]
  if (mode(fs)=="numeric")
  {
     fs <- as(lapply(fcs.path,function(x) {
       f<-read.FCS(x,dataset=2)
     markers <- markers[-c(which(is.na(markers)) , which(markers=="BLANK"))]
    stop("This should be rewritten.")
  }))
  }
  fs<- fsApply(fs, function(x) {
    x@parameters@data[,2] <- c("FS-A","FSC-H","FS-W","SS-A","SSC-H","SSC-W",markers,"Time");
    return(x)
  })
}


sampleNames(fs)<-unlist(lapply(fcs.path, function(x) return(unlist(strsplit(x,split=paste("/",type,"/",sep="",collapse="")))[2])))
#Reading compensation files
#-----------------------------------------------------------------------------------

if(type=="healthy volunteer")
{
  
  comp.path <- list.files(paste(dir.path,dataset,"compensation",sep = "/"),full.names=T,pattern="phase 1.csv")
  names(comp.path) <- unlist(lapply(comp.path, function(x) unlist(strsplit(unlist(strsplit(x,"gating "))[2]," phase 1.csv"))[1]))
  sub.centre <- unlist(lapply(sampleNames(fs), function(x) return(unlist(strsplit(unlist(strsplit(x,"ONE_05 "))[2],"_"))[1])))
  names(sub.centre) <-sampleNames(fs)
  
}

#Marker names
#---------------------------------------------------------------------------------------
#****************************************************************************************
#                    Finding markers
#****************************************************************************************
channels.ind <- c(1:ncol(fs[[1]]))
if(type!="Fortessa"){
  names(channels.ind) <- c("FS-A","FS-W","SS-A",markers)
}else{
  names(channels.ind) <- c("FS-A","FSC-H","FS-W","SS-A","SSC-H","SSC-W",markers)
}
#Compensating and transforming files
#----------------------------------------------------------------------------------------
#*****************************************************************************************
#                        /Compensation/
#*****************************************************************************************
#Using the compensation matrix provided

if(type=="healthy volunteer")
{
  fs <- fsApply(fs, function(x) {
    comp <-  read.csv(comp.path[sub.centre[identifier(x)]],check.names = F)
    comp <- comp[,-1]
    return(compensate(x, comp))})
}
if(type=="post HSCT")
{
  comp<- read.csv(paste(dir.path,dataset,"compensation","Compensation for automated gating all sites phase 2.csv",sep = "/"),check.names = F)
  comp <- comp[,-1]
  fs <- fsApply(fs, function(x) return(compensate(x, comp)))
}

if(type=="Navios")
{
  comp<- read.csv(paste(dir.path,dataset,"compensation","2018Jun13_Navios comp adj. for Basic.csv",sep="/"),check.names = F)
  comp <- comp[,-1]
  fs <- fsApply(fs, function(x) return(compensate(x, comp)))
}
if(type=="Fortessa")
{
  comp<- read.csv(paste(dir.path,dataset,"compensation","2018Jun13 Fortessa adj. Basic.csv",sep="/"),check.names = F)
  comp <- comp[,-1]
  fs <- fsApply(fs, function(x) return(compensate(x, comp)))
}

#Doublet removal
#----------------------------------------------------------------------------------------
#*****************************************************************************************
#                      Removing doublets
#*****************************************************************************************
if(type!="Fortessa")
{
  n<-fsApply(fs, function(x) ifelse(sd(exprs(x)[,channels.ind["FS-W"]])>35,yes = 1.8,no = 3.8)) 
  
  fsw.gate <-fsApply(fs,function(x) return(tail(getPeaks(x,channels.ind["FS-W"],tinypeak.removal=.4)$Peaks,1)+n[identifier(x),1]*sd(exprs(x)[,channels.ind["FS-W"]])))
  names(fsw.gate) <- sampleNames(fs)
  Singlet <- fsApply_pb(fs, function(x) {
    message("Gating singlets")
    print(identifier(x))
    gate <- max(exprs(x)[,channels.ind["FS-A"]])-2
    exprs(x)[which(exprs(x)[,channels.ind["SS-A"]]> max(exprs(x)[,channels.ind["SS-A"]])-2),] <-NA
    singlets <- flowDensity(x,channels=channels.ind[c("FS-W","FS-A")],position=c(F,NA),gates=c(fsw.gate[identifier(x)],NA))
    return(singlets)
  })
}else{
  Singlet <- fsApply_pb(fs, function(x) {
    message("Gating singlets")
    channels <- c("FSC-A","FSC-H")
    rot <- rotate.data(x, channels)
    gate.1 <- max(deGate(rot$data, channels[2],percentile=NA,tinypeak.removal = .98,twin.factor = .2,magnitude = .1, upper=T, alpha=.005,verbose = F),
                  deGate(rot$data, channels[2],percentile=NA,tinypeak.removal = .98,twin.factor = .2))
    gate.2 <- min(deGate(rot$data, channels[2],percentile=NA,tinypeak.removal = .98,twin.factor = .2,magnitude = .1, upper=F, alpha=.0005,verbose = F),
                  deGate(rot$data, channels[2],tinypeak.removal = .98,percentile=.003))
    temp<- flowDensity(rot$data,channels,position = c(NA,F),gates=c(NA,gate.1),verbose = F)
    singlets<- rotate.fd(flowDensity(temp,channels,position = c(NA,T),gates=c(NA,gate.2),verbose = F),angle = rot$theta)
    
    return(singlets)
  })
}

fs.sngl <-as(lapply(Singlet, function(x) x@flow.frame), Class="flowSet")

fs.temp <-as(lapply(Singlet, function(x) getflowFrame(x)), Class="flowSet")
#*****************************************************************************************
#                      /Transformation/
#*****************************************************************************************


if (type=="healthy volunteer")
{
  lgl <- lapply(unique(sub.centre), function(x){
    return(estimate.logicle(fs.temp[which(sub.centre==x)],med=T,trans.chans=grep(colnames(fs[[1]]),pattern="FL"),estimate=T,return.set =F))
  })
  names(lgl) <-  unique(sub.centre)
  fs.sngl<-fsApply(fs.temp, function(x){
    x<- transform(x,  lgl[[as.vector(sub.centre[identifier(x)])]])
    y<-fs.sngl[[identifier(x)]]
    na.inds <- which(!is.na(exprs(y)[,1]))
    exprs(y)[na.inds,] <- exprs(x)
    return(y)
  })
  
}else if (type=="Fortessa"){
  lgl <- estimate.logicle(fs.temp,med=T,estimate=T,return.set =F)
  fs.sngl<-fsApply(fs.temp, function(x){
    x<- transform(x,  lgl)
    y<-fs.sngl[[identifier(x)]]
    na.inds <- which(!is.na(exprs(y)[,1]))
    exprs(y)[na.inds,] <- exprs(x)
    return(y)
  })
}else{
  lgl <- estimate.logicle(fs.temp,med=T,trans.chans=grep(colnames(fs[[1]]),pattern="FL"),estimate=T,return.set =F)
  fs.sngl<-fsApply(fs.temp, function(x){
    x<- transform(x,  lgl)
    y<-fs.sngl[[identifier(x)]]
    na.inds <- which(!is.na(exprs(y)[,1]))
    exprs(y)[na.inds,] <- exprs(x)
    return(y)
  })
}

#Gating CD45, and Granulocytes
#----------------------------------------------------------------------------------------
#*****************************************************************************************
#                      Gating CD45
#*****************************************************************************************

cd45.gate <- fsApply(fs.sngl,function(x) {return(deGate(x,channels.ind["CD45"],percentile=NA, tinypeak.removal=1/80,upper=F,alpha=.005,verbose = F))})

##############Attention###########################
#Changed tinypeak.removal from 1/50 to 1/8 and added all.cuts=T, should be checked for the PhaseI
##################################################

names(cd45.gate) <- sampleNames(fs.sngl)

CD45 <- fsApply_pb(fs.sngl, function(x) {
  message("Gating CD45")
  if (cd45.gate[identifier(x)]>2.6)
    cd45.gate[identifier(x)]<- deGate(x,channels.ind["CD45"],percentile=NA, tinypeak.removal=.99,upper=F,alpha=.01,verbose = F)
  cd45 <- flowDensity(x, channels.ind[c("CD45","SS-A")],position=c(T,NA),gates=c(cd45.gate[identifier(x)],NA))
  cd45.lo <- deGate(getflowFrame(cd45),channels.ind["CD45"],use.upper = T,upper=F,tinypeak.removal = .1,verbose = F)
  cd45 <- flowDensity(x, channels.ind[c("CD45","SS-A")],position=c(T,NA),gates=c(cd45.lo*.95,NA))
  return(cd45)
})

fs.45 <-lapply(CD45, function(x) x@flow.frame)
fs.45 <- as(object=fs.45, Class="flowSet")

PBMC<- fsApply_pb(fs.45, function(x){
  message("Gating PBMC")
  x2<- rotate.data(x,chans = channels.ind[c("FS-A","SS-A")],theta =-pi/4 )$data
  fs.gate <- min(deGate(x2,channels.ind["FS-A"],upper=F,tinypeak.removal = .1,verbose = F),-170000)
  temp <- flowDensity(x2,channels.ind[c("FS-A","SS-A")],position = c(T,NA),gates=c(fs.gate, NA))
  ss.gate <- tail(deGate(temp@flow.frame,channels.ind["SS-A"],tinypeak.removal = 1/20,all.cuts = T,verbose = F),1)
  if(CD45[[identifier(x)]]@cell.count<2){
    ss.gate <-NA
  }else{
    if(ss.gate < deGate(temp@flow.frame,channels.ind["SS-A"],tinypeak.removal = 1/20,percentile = .1,use.percentile = T,verbose = F))
      ss.gate <- deGate(temp@flow.frame,channels.ind["SS-A"],tinypeak.removal = 1/20,use.upper=T,upper=T,alpha=.05,verbose = F)
  }
  filt <- rotate.fd(flowDensity(x2,channels.ind[c("FS-A","SS-A")],position = c(T,F),gates=c(fs.gate, ss.gate)),angle =-pi/4 )@filter
  not.gran <-flowDensity(x,channels.ind[c("FS-A","SS-A")],position = c(T,F),filter=filt)
  
  peaks<- getPeaks(not.gran,channels.ind["FS-A"])
  main.p <- peaks$Peaks[which.max(peaks$P.h)]*.95
  all.gates <- deGate(not.gran,channels.ind["FS-A"],all.cuts = T,tinypeak.removal = 1/50,twin.factor = .9,upper=F,alpha=.05,verbose = F)
  if (length(which(all.gates<main.p))>0 & type %in% c("Navios","Fortessa"))
  {
    deb.gate <-tail(all.gates[which(all.gates<main.p)],1)
    filt <- flowDensity(not.gran,channels.ind[c("FS-A","SS-A")],position = c(T,NA),gates=c(deb.gate,NA))@filter
    not.debris <- flowDensity(not.gran,channels.ind[c("FS-A","SS-A")],position = c(T,NA),filter=filt)
  }else{
    not.debris <- not.gran
  }
  return(not.debris)
})
fs.pbmc<-lapply(PBMC, function(x) x@flow.frame)
fs.pbmc <- as(object=fs.pbmc, Class="flowSet")
#*****************************************************************************************
#                      Gating Lymph
#*****************************************************************************************
multi <-2

fs.hi <- fsApply(fs.pbmc,deGate,channels.ind["FS-A"],upper=T,percentile=NA,tinypeak.removal=.9,verbose = F)
Lymph <- fsApply_pb(fs.pbmc, function(x) {
  message("Gating Lymphocytes")
  x2 <- rotate.data(x,chans = channels.ind[c("FS-A","SS-A")],theta =pi/6)$data
  peaks<- getPeaks(x2,channels.ind["FS-A"],tinypeak.removal = .05)$Peaks
  peaks <- peaks[which(peaks>45000)]
  fs.gate <-ifelse(test = length(peaks)>2,yes = deGate(x2,channels.ind["FS-A"],tinypeak.removal = .05,all.cuts =T,verbose = F)[2],
                   no=ifelse((length(peaks)>1 &length(peaks)==2 & peaks[2]>50000),no = peaks[1]+multi*sd(exprs(x2)[,channels.ind["FS-A"]],na.rm = T),
                             yes=deGate(x2,channels.ind["FS-A"],verbose = F,tinypeak.removal = .05)))
  fs.hi <- deGate(x2,channels.ind["FS-A"],upper=T,percentile=NA,tinypeak.removal=.9,verbose = F)*1.05
  fs.hi <-ifelse(fs.hi<40000,yes = deGate(x2,channels.ind["FS-A"],upper=T,percentile=NA,tinypeak.removal=.2,verbose = F)*1.05,no = fs.hi)
  fs.gate <- min(fs.gate, fs.hi,na.rm=T)
  ss.gate <- deGate(x,channels.ind["SS-A"], tinypeak.removal=.9,percentile=NA,upper=T,alpha=.01,verbose = F)
  temp <- flowDensity(x2,channels.ind[c("FS-A","SS-A")],position=c(F,NA), gates=c(fs.gate,NA))
  no.mono <- rotate.data(getflowFrame(temp),chans = channels.ind[c("FS-A","SS-A")],theta =- pi/6)$data
  filt <- flowDensity(no.mono, channels.ind[c("FS-A","SS-A")],position=c(NA,F), gates=c(NA,ss.gate*1.2))@filter
  lymph <- flowDensity(x,channels.ind[c("FS-A","SS-A")],position=c(NA,F),filter=filt)
  return(lymph)
  
})
fs.lymph <-lapply(Lymph, function(x) x@flow.frame)
fs.lymph <- as(object=fs.lymph, Class="flowSet")
#*****************************************************************************************
#                      Gating Bcells
#*****************************************************************************************
cd19.gate <- averageGates(fsApply(fs.lymph, deGate,channels.ind["CD19"],tinypeak.removal=.05,upper=T,percentile=NA,alpha=.005,verbose = F),sd.coeff=2)
names(cd19.gate) <- sampleNames(fs.lymph)
Bcells<- fsApply_pb(fs.lymph, function(x) {
  message("Gating Bcells")
  temp<- flowDensity(x,channels.ind[c("CD19","SS-A")],position=c(T,NA),gates=c(cd19.gate[identifier(x)],NA))
  cd19.hi <- max(deGate(temp@flow.frame,channels.ind["CD19"], percentile = NA,upper=F,alpha=.05,tinypeak.removal = 1/10,bimodal=T,verbose = F)
                 ,cd19.gate[identifier(x)]*1.1)
  bcell <-  flowDensity(x, channels.ind[c("CD19","SS-A")],position=c(T,F),gates=c(cd19.hi,NA),percentile=c(F,.99), use.percentile=c(F,T))
  return(bcell)
})
fs.bcell <-lapply(Bcells, function(x) x@flow.frame)
fs.bcell <- as(object=fs.bcell, Class="flowSet")
which.low <- which(unlist(lapply(Bcells,function(x) x@cell.count))<500)
if (length(which.low)==length(fs.bcell)){
  stop("All samples have less than 500 Bcells. Pipeline stops here!")
  
}else if (length(which.low)>0)
{ print(paste(paste(sampleNames(fs.bcell)[which.low],sep=","),"less than 1000 Bcells --> will be removed"),sep=" ")
  for ( i in which.low)
    fs.bcell[[i]] <- getflowFrame(Bcells[[i]])[1,]
}

#*****************************************************************************************
#                      Gating Bcell subsets
#*****************************************************************************************
cd27.gate <- averageGates(fsApply(fs.pbmc,deGate,channels.ind["CD27"],verbose = F),sd.coeff=2)
names(cd27.gate)  <- sampleNames(fs.bcell)

CD27.g <- fsApply_pb(fs.pbmc, function(x) {
  temp<-getflowFrame(flowDensity(x, channels.ind[c("CD27","SS-A")],position=c(T,NA),gates=c(cd27.gate[identifier(x)],NA)))
  gate <- tail(getPeaks(temp,channels.ind["CD27"],tinypeak.removal = .1)$Peaks,1)-1*sd(exprs(temp[,channels.ind["CD27"]]))
 
  return(gate)
})

igd.gate <- averageGates(as.numeric(fsApply(fs.bcell,function(x) return(deGate(x,channels.ind["IgD"],after.peak = F,percentile=NA,upper=F,verbose = F)[1]))),sd.coeff = 2)
cd27.gate <- averageGates(fsApply(fs.pbmc,deGate,channels.ind["CD27"],verbose = F),sd.coeff=2)
names(cd27.gate)  <- names(igd.gate) <-sampleNames(fs.bcell)
cd27.pos<- CD27.g[,1]
names(cd27.pos) <- sampleNames(fs.bcell)
naive <- fsApply_pb(fs.bcell, function(x) {
  message("Gating Naive population")
  quad.1 <- flowDensity(x, channels.ind[c("CD27","IgD")],position=c(F,F),gates=c(cd27.gate[identifier(x)]*1.05,igd.gate[identifier(x)]))
  quad.2 <- flowDensity(x,channels.ind[c("CD27","IgD")],position=c(T,F), gates=quad.1@gates)
  quad.3 <- flowDensity(x, channels.ind[c("CD27","IgD")],position=c(F,T), gates=quad.1@gates)
  quad.4 <- flowDensity(x, channels.ind[c("CD27","IgD")],position=c(T,T), gates=quad.1@gates)
  cd27.neg <-flowDensity(fs.pbmc[[identifier(x)]], channels.ind[c("CD27","IgD")],position=c(F,NA),gates=c(cd27.gate[identifier(x)],igd.gate[identifier(x)]))
  return(list(quad1=quad.1,quad2=quad.2,quad3=quad.3, quad4=quad.4,cd27neg=cd27.neg))
})
#*****************************************************************************************
#                        Gating IgD/IgM
#*****************************************************************************************
igm.gate <- averageGates(as.numeric(fsApply(fs.bcell,deGate,channels.ind["IgM"],percentile=NA,upper=F,verbose = F)),sd.coeff = 2)
names(igm.gate)  <- sampleNames(fs.bcell)
IG <- fsApply_pb(fs.bcell, function(x) {
  message("Gating IgM/IgD")
  temp <- flowDensity(x, channels.ind[c("IgD","IgM")],position=c(T,NA),gates=c(igd.gate[identifier(x)],NA))
  igm.gate <- deGate ( temp@flow.frame, channels.ind["IgM"],percentile=NA,tinypeak.removal = .99,use.upper=T,upper=F,alpha=.8 ,verbose = F)     
  quad.1 <- flowDensity(x, channels.ind[c("IgD","IgM")],position=c(F,F),gates=c(igd.gate[identifier(x)],igm.gate))
  quad.2 <- flowDensity(x,channels.ind[c("IgD","IgM")],position=c(T,F), gates=quad.1@gates)
  quad.3 <- flowDensity(x, channels.ind[c("IgD","IgM")],position=c(F,T), gates=quad.1@gates)
  quad.4 <- flowDensity(x, channels.ind[c("IgD","IgM")],position=c(T,T), gates=quad.1@gates)
  igm <- flowDensity(x, channels.ind[c("IgD","IgM")],position=c(NA,T), gates=quad.1@gates)
  return(list(quad1=quad.1,quad2=quad.2,quad3=quad.3, quad4=quad.4, igm=igm))
})
fs.ig1 <- lapply(IG, function(x) x$quad1@flow.frame)
fs.ig1 <- as(object=fs.ig1, Class="flowSet")
fs.ig4 <- lapply(IG, function(x) x$quad4@flow.frame)
fs.ig4 <- as(object=fs.ig4, Class="flowSet")
fs.igm <- lapply(IG, function(x)x$igm@flow.frame)
fs.igm <- as(object=fs.igm, Class="flowSet")
which.low <- which(unlist(lapply(IG,function(x) x$quad1@cell.count))<100)
if (length(which.low)==length(fs.ig1)){
  print("All samples have less than 100 IgD-IgM-. Pipeline stops here!")
}else if (length(which.low)>0)
{ print(paste(paste(sampleNames(fs.ig1)[which.low],sep=","),"less than 100 IgD-IgM- --> will be removed"),sep=" ")
  for ( i in which.low)
  {
    fs.ig1[[i]] <- getflowFrame(IG[[i]]$quad1)[1,]
    IG[[i]]$quad1@index<-0
  }
}
#############
cd24.gate <- averageGates(fsApply(fs.sngl,deGate,channels.ind["CD24"],verbose = F),sd.coeff=1)
names(cd24.gate) <- sampleNames(fs.45)

cd38.gate <- fsApply(fs.pbmc, function(x){
  temp <- getflowFrame(flowDensity(x,channels = channels.ind[c("CD27","CD24")],position = c(F,F),gates=c(cd27.gate[identifier(x)],cd24.gate[identifier(x)])))
  cd38.gate <- deGate(temp,channels.ind["CD38"],upper=F,percentile=NA,alpha=.2,verbose = F)
  cd38.hi <- max(getPeaks(temp,channels.ind["CD38"],tinypeak.removal = .5)$Peaks)+.3*sd(exprs(temp)[,channels.ind["CD38"]])
  cd38.hi2 <- max(getPeaks(temp,channels.ind["CD38"],tinypeak.removal = .5)$Peaks)+.45*sd(exprs(temp)[,channels.ind["CD38"]])
  cd38.up <- deGate(temp, channels.ind["CD38"],use.upper = T, upper=T, magnitude = .7, alpha=.5,verbose = F)
  return(c(cd38.gate,cd38.hi,cd38.hi2,cd38.up))
})
#*****************************************************************************************
#                        Gating plasmablast on IgD-IgM-
#*****************************************************************************************

Plasma.switch <- fsApply_pb(fs.ig1, function(x) {
  message("Gating plasma")
  quad.1 <- flowDensity(x,channels.ind[c("CD38","CD27")],position=c(F,F),gates=c(cd38.gate[identifier(x),3],cd27.gate[identifier(x)]))
  quad.2 <- flowDensity(x, channels.ind[c("CD38","CD27")],position=c(T,F),gates=quad.1@gates)
  quad.3 <- flowDensity(x,channels.ind[c("CD38","CD27")],position=c(F,T),gates=quad.1@gates)
  quad.4 <- flowDensity(x,channels.ind[c("CD38","CD27")],position=c(T,T),gates=quad.1@gates)
  plasma <- flowDensity(x, channels.ind[c("CD38","CD27")],position=c(T,T),gates=c(quad.1@gates[1],cd27.pos[identifier(x)]))
  return(list(quad1=quad.1,quad2=quad.2,switched=quad.3, quad4=quad.4,plasma=plasma))
})
#*****************************************************************************************
#                        Gating IgM Memory
#*****************************************************************************************
Memory <- fsApply_pb(fs.igm, function(x) {
  message("Gating IgM memory")
  quad.1 <- flowDensity(x, channels.ind[c("CD38","CD27")],position=c(F,F),gates=c(cd38.gate[identifier(x),2],cd27.gate[identifier(x)]))
  quad.2 <- flowDensity(x,channels.ind[c("CD38","CD27")],position=c(T,F),gates=quad.1@gates)
  quad.3 <- flowDensity(x, channels.ind[c("CD38","CD27")],position=c(F,T),gates=quad.1@gates)
  quad.4 <- flowDensity(x,channels.ind[c("CD38","CD27")],position=c(T,T),gates=quad.1@gates)
  cd27.neg  <- flowDensity(x, channels.ind[c("CD38","CD27")],position=c(NA,F),gates=quad.1@gates)
  return(list(quad1=quad.1,quad2=quad.2,quad3=quad.3, quad4=quad.4,cd27=cd27.neg))
})
fs.27neg <- lapply(Memory, function(x) x$cd27@flow.frame)
fs.27neg <- as(object=fs.27neg, Class="flowSet")
#*****************************************************************************************
#                        Gating Transitional on IgM+CD27-
#*****************************************************************************************
#Setting CD24hi on CD45+ cells

CD24.g <- fsApply(fs.45, function(x) {
  
  temp<- flowDensity(x, channels.ind[c("CD24","SS-A")],position=c(T,NA),gates=c(cd24.gate[identifier(x)],NA))
  gate <- getPeaks(temp@flow.frame,channels.ind["CD24"],tinypeak.removal = .9)$Peaks-1.6*sd(exprs(getflowFrame(temp)[,channels.ind["CD24"]]))
  return(gate)
})
cd24.hi <- averageGates(CD24.g[,1],sd.coeff=1)
names(cd24.hi) <-  sampleNames(fs.27neg)

Trans <- fsApply_pb(fs.27neg, function(x) {
  message("Gating Transitional")
  trans <- flowDensity(x, channels.ind[c("CD38","CD24")],position=c(T,T),gates=c(cd38.gate[identifier(x),2],cd24.hi[identifier(x)]))
  return(trans)
  
})
#*****************************************************************************************
#                        Gating CD21Lo
#*****************************************************************************************

cd21.gate <- as.numeric(fsApply(fs.pbmc,deGate,channels.ind["CD21"],percentile=NA, upper=T,twin.factor=.2,alpha=.008,verbose = F))


cd38.lo<- averageGates(fsApply(fs.pbmc,deGate,channels.ind["CD38"],percentile=NA, upper=F, alpha=.99,verbose = F),sd.coeff=1)
names(cd21.gate)<- names(cd38.lo)<-sampleNames(fs.ig4)

CD21Lo <- fsApply_pb(fs.bcell, function(x) {
  message("Gating CD21")
  cd21 <- flowDensity(x, channels.ind[c("CD38","CD21")],position=c(F,F),gates=c(cd38.gate[identifier(x),1],cd21.gate[identifier(x)]))
  return(cd21)
  
})

cell.prop <-sapply(1:length(Bcells), function(fd){
return( c((Bcells[[fd]]@cell.count/Lymph[[fd]]@cell.count*100),
                                (c(naive[[fd]]$quad3@cell.count,naive[[fd]]$quad4@cell.count,
                                   Trans[[fd]]@cell.count,naive[[fd]]$quad2@cell.count,
                                   naive[[fd]]$quad1@cell.count,Memory[[fd]]$cd27@cell.count,Memory[[fd]]$quad3@cell.count, 
                                   Plasma.switch[[fd]]$switched@cell.count, Plasma.switch[[fd]]$plasma@cell.count,
                                   CD21Lo[[fd]]@cell.count)/Bcells[[fd]]@cell.count*100 ),Bcells[[fd]]@cell.count,IG[[fd]]$quad1@cell.count))})

rownames(cell.prop) <- c("CD19+","Naive B cells","Marginal zone B cells","Transitional B", "IgD-CD27+","IgD-CD27-", 
                         "IgM+CD27-","IgM memory","Class switched memory",
                         "Plasmablasts","CD21Low","CD19+ counts","IgM-IgD- counts")
colnames(cell.prop) <- names(Bcells)
write.csv(cell.prop,paste0(res.path,"/",panel,"_Proportions_",format(Sys.time(), "%Y-%m-%d"), ".csv"))
