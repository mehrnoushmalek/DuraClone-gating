#Panel 1
#Author: Mehrnoush Malek
#Date: Last revised in October 2018

#Read the second dataset in the LMD file as it's not transformed or compensated
#In here we only keep the markers names from the first data set in order to use it for the second dataset within each FCS file
#Reading files
#------------------------------------------------------------------------------------
f <- read.FCS(fcs.path[1],dataset = 1)  
markers <- unlist(lapply(as.vector(f@parameters@data[,2]), function(x) unlist(strsplit(x," "))[1]))
if (type!="Fortessa")
{
  markers <- markers[-c(which(markers=="FS"),which(markers=="SS"))]
  fs <- tryCatch(read.flowSet(fcs.path,dataset = 2),error=function(x) {return(1)})
  if (mode(fs)=="numeric")
  {
    stop("All your files need to be consistent. Order or number of channels mismatch.")
    
  }else{
  fs<- fsApply(fs, function(x) {
    x@parameters@data[,2] <- c("FS-A","FS-W","SS-A",markers);
    return(x)
  })
  }
}else{
  markers <- markers[-which(is.na(markers))]
  fs <-  tryCatch(read.flowSet(fcs.path,dataset = 2),error=function(x) {return(1)})
  ##Temporary for new fortessa files, order of channles are messed up
  if (mode(fs)=="numeric")
  {
   warning("All your files need to be consistent. Order or number of channels mismatch.")
    fs <- as(lapply(fcs.path,function(x) {
        f<-read.FCS(x,dataset=2)
         f <- f[,c("FSC-A","FSC-H","FSC-W","SSC-A","SSC-H","SSC-W","FITC-A","PE-A","ECD-A","PE-Cy5.5-A","PE-Cy7-A","APC-A","Alexa Fluor 700-A","APC-Alexa 750-A","Pacific Blue-A","V500-A","Time")]
       }), Class="flowSet")

  }
  fs<- fsApply(fs, function(x) {
    x@parameters@data[,2] <- c("FS-A","FSC-H","FS-W","SS-A","SSC-H","SSC-W",markers,"Time");
    return(x)
  })
}

sampleNames(fs)<-unlist(lapply(fcs.path, function(x) return(unlist(strsplit(x,split=paste("/",type,"/",sep="",collapse="")))[2])))
#Reading compensation files
print(sampleNames(fs))
#-----------------------------------------------------------------------------------

if(type=="healthy volunteer")
{

  comp.path <- list.files(paste(dir.path,dataset,"compensation",sep = "/"),full.names=T,pattern="phase 1.csv")
  names(comp.path) <- unlist(lapply(comp.path, function(x) unlist(strsplit(unlist(strsplit(x,"gating "))[2]," phase 1.csv"))[1]))
  sub.centre <- unlist(lapply(sampleNames(fs), function(x) return(unlist(strsplit(unlist(strsplit(x,"ONE_01 "))[2],"_"))[1])))
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
    gate <- max(exprs(x)[,1])-2
    exprs(x)[which(exprs(x)[,3]> max(exprs(x)[,3])-2),] <-NA
    singlets <- flowDensity(x,channels=channels.ind[c(2,1)],position=c(F,F),gates=c(fsw.gate[identifier(x)],gate))
    return(singlets)
  })
}else{
  Singlet <- fsApply_pb(fs, function(x) {
    message("Gating singlets")
    channels <- c("FSC-A","FSC-H")
    y <- nmRemove(x, c("FSC-A","SSC-A"))
    rot <- rotate.data(y, channels)
    gate.1 <- max(deGate(rot$data, channels[2],percentile=NA,tinypeak.removal = .98,twin.factor = .2,magnitude = .1, upper=T, alpha=.005,verbose = F),
                  deGate(rot$data, channels[2],percentile=NA,tinypeak.removal = .98,twin.factor = .2))
    gate.2 <- min(deGate(rot$data, channels[2],percentile=NA,tinypeak.removal = .98,twin.factor = .2,magnitude = .1, upper=F, alpha=.0005,verbose = F),
                  deGate(rot$data, channels[2],tinypeak.removal = .98,percentile=.001))
    temp<- flowDensity(rot$data,channels,position = c(NA,F),gates=c(NA,gate.1),verbose = F)
    filt<- rotate.fd(flowDensity(temp,channels,position = c(NA,T),gates=c(NA,gate.2),verbose = F),angle = rot$theta)@filter
    singlets <- flowDensity(x,channels,position = c(NA,T),filter=filt)
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
  
}else if (type=="Navios"){
  lgl <- estimate.logicle(fs.temp,med=F,trans.chans=grep(colnames(fs[[1]]),pattern="FL"),estimate=T,m=5,return.set =F)
  fs.sngl<-fsApply(fs.temp, function(x){
    x<- transform(x,  lgl)
    y<-fs.sngl[[identifier(x)]]
    na.inds <- which(!is.na(exprs(y)[,1]))
    exprs(y)[na.inds,] <- exprs(x)
    return(y)
  })
}else if (type=="Fortessa"){
  trans.chans<- channels.ind[-c(grep("fs", tolower(names(channels.ind))),grep("ss", tolower(names(channels.ind))),which(is.na(names(channels.ind))))]
  lgl <- estimate.logicle(fs.temp,trans.chans = trans.chans,med=T,estimate=T, return.set =F)
  fs.sngl<-fsApply(fs.temp, function(x){
    x<- transform(x,  lgl)
    y<-fs.sngl[[identifier(x)]]
    na.inds <- which(!is.na(exprs(y)[,1]))
    exprs(y)[na.inds,] <- exprs(x)
    return(y)
  })
}else {
  lgl <- estimate.logicle(fs.temp,med=T,trans.chans=grep(colnames(fs[[1]]),pattern="FL"),estimate=T,return.set =F)
  fs.sngl<-fsApply(fs.temp, function(x){
    x<- transform(x,  lgl)
    y<-fs.sngl[[identifier(x)]]
    na.inds <- which(!is.na(exprs(y)[,1]))
    exprs(y)[na.inds,] <- exprs(x)
    return(y)
  })
}

#----------------------------------------------------------------------------------------
#*****************************************************************************************
#                      Gating CD45
#*****************************************************************************************
if (type %in% c("Navios","Fortessa"))
{
lo.ss <- fsApply(fs.sngl, function(x){
  return(flowDensity(x, channels.ind[c("CD45","SS-A")],position = c(NA,F),gates=c(NA, 200000))@flow.frame)
})
cd45.gate <- averageGates(fsApply(lo.ss,function(x) deGate(x,channels.ind["CD45"],percentile=NA,
                                                    tinypeak.removal=1/60,upper=F,alpha=.01,verbose=F)),sd.coeff = 2)
}else{
  cd45.gate <- averageGates(fsApply(fs.sngl,function(x) deGate(x,channels.ind["CD45"],percentile=NA, 
                                                        tinypeak.removal=1/60,upper=F,all.cuts=T,alpha=.01)[1]),sd.coeff = 2.5)
  
}
names(cd45.gate) <- sampleNames(fs.sngl)

CD45 <- fsApply_pb(fs.sngl, function(x) {
  message("Gating CD45")
  y <- nmRemove(x, channels.ind["SS-A"])
  cd45.nomargin <- flowDensity(y, channels.ind[c("CD45","SS-A")],position=c(T,NA),gates=c(cd45.gate[identifier(x)],NA))
  cd45 <- flowDensity(x, channels.ind[c("CD45","SS-A")],position=c(T,NA),filter=cd45.nomargin@filter)
  return(cd45)
})
fs.45 <-lapply(CD45, function(x) x@flow.frame)
fs.45 <- as(object=fs.45, Class="flowSet")

ss.gate <- averageGates(fsApply(fs.45,deGate,channels.ind["SS-A"],percentile=.999,tinypeak.removal=1/5,verbose=F))
gran.gate <- averageGates(fsApply(fs.45,deGate,channels.ind["CD14"],verbose=F) )
names(ss.gate) <- names(gran.gate) <- sampleNames(fs.45)


PBMC <- fsApply_pb(fs.45, function(x){
  message("Gating PBMC")
  temp.1 <- flowDensity(x,channels.ind[c("CD14","SS-A")],position = c(F,F),gates=c(gran.gate[identifier(x)], ss.gate[identifier(x)]))
  ss.lo<- deGate(getflowFrame(temp.1),channels.ind["SS-A"],percentile = NA,upper=T,alpha=.01,verbose=F)
  ss.lo <-ifelse(ss.lo < deGate(getflowFrame(temp.1),channels.ind["SS-A"],percentile = .5,use.percentile=T),no = ss.lo,
                 yes = deGate(getflowFrame(temp.1),channels.ind["SS-A"],upper=T,use.upper =T,tinypeak.removal=.9,alpha=.05) )
  s.gate <- min(ss.lo, ss.gate[identifier(x)]*.93)
  temp.2 <- notSubFrame(x,channels.ind[c("CD14","SS-A")],position = c(F,F),gates=c(gran.gate[identifier(x)], s.gate))
  rot<- rotate.data(getflowFrame(temp.2),channels.ind[c("CD14","SS-A")],min.max = T)
  cd14.up <- deGate(rot$data, channels.ind["CD14"],percentile = NA,tinypeak.removal = 1/15,upper=F,alpha=.1,verbose=F)
  gran <- rotate.fd(flowDensity(rot$data,channels.ind[c("CD14","SS-A")],position = c(F,NA),gates=c(cd14.up, NA)),angle = rot$theta)
  not.gran <- notSubFrame(x,channels.ind[c("CD14","SS-A")],position = c(F,T),filter=gran@filter)
  not.gran@filter <- gran@filter
  plotDens(x,channels.ind[c("CD14","SS-A")],main=paste("CD45+",identifier(x)))
  lines(gran@filter,lty=3,lwd=3)
  return(not.gran)
})

fs.pbmc<-lapply(PBMC, function(x) x@flow.frame)
fs.pbmc <- as(object=fs.pbmc, Class="flowSet")

#----------------------------------------------------------------------------------------
#*****************************************************************************************
#                      Gating Lymph
#*****************************************************************************************
multi <- 1.2
fs.hi <- fsApply_pb(fs.pbmc,deGate,channels.ind["FS-A"],upper=T,percentile=NA,tinypeak.removal=.9,verbose=F)
Lymph <- fsApply(fs.pbmc, function(x) {
  message("Gating Lymphocytes")
  peaks<- getPeaks(x,channels.ind["FS-A"])
  main.p <- peaks$Peaks[which(peaks$P.h>max(peaks$P.h)*.95)[1]]
  all.gates <- c(deGate(x,channels.ind["FS-A"],use.upper=T,upper=F,alpha=.05,verbose=F),
                 deGate(x,channels.ind["FS-A"],all.cuts = T,tinypeak.removal = 1/70,twin.factor = .9,upper=F,alpha=.05,verbose=F))
  if (length(which(all.gates<main.p))>0)
  {
    deb.gate <-tail(all.gates[which(all.gates<main.p)],1)
    filt <- flowDensity(x,channels.ind[c("FS-A","SS-A")],position = c(T,NA),gates=c(deb.gate,NA))@filter
    not.debris <- flowDensity(x,channels.ind[c("FS-A","SS-A")],position = c(T,NA),filter=filt)@flow.frame
  }else{
    not.debris <- x
  }
  x2 <- rotate.data(not.debris,chans = channels.ind[c("FS-A","SS-A")],theta =pi/6)$data
  peaks<- getPeaks(x2,channels.ind["FS-A"],tinypeak.removal = .05)$Peaks
  peaks <- peaks[which(peaks>48000)]
  fs.gate <-ifelse(test = length(peaks)>2,yes = deGate(x2,channels.ind["FS-A"],tinypeak.removal = .05,all.cuts =T)[2],no=peaks[1]+multi*sd(exprs(x2)[,channels.ind["FS-A"]],na.rm = T))
  fs.hi <- deGate(x2,channels.ind["FS-A"],upper=T,percentile=NA,tinypeak.removal=.9,verbose=F)*1.05
  fs.hi <-ifelse(fs.hi<400000,yes = deGate(x2,channels.ind["FS-A"],upper=T,percentile=NA,tinypeak.removal=.2,verbose=F)*1.05,no = fs.hi)
  fs.gate <- min(fs.gate, fs.hi)
  ss.gate <- deGate(x,channels.ind["SS-A"], tinypeak.removal=.9,percentile=NA,upper=T,alpha=.01,verbose=F)
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
Bcells <- fsApply_pb(fs.lymph, function(x) {
  message("Gating Bcell/Tcell")
  
  temp <-  getflowFrame(flowDensity(x, channels.ind[c("CD19","CD3")],position=c(NA,F)))
  peaks <-getPeaks(x,channels.ind["CD19"],tinypeak.removal = 1/45)$Peaks
  flag <-0
  if (length(peaks)<2)
    flag<-1
  cd19.gate<- deGate(temp,channels.ind["CD19"],percentile=NA,after.peak = T,tinypeak.removal = 1/20, upper=T,alpha=0.01,verbose=F)
  temp <- removeMargins(temp, channels.ind["CD19"],neg = .5)
  cd19.hi<- deGate(temp,channels.ind["CD19"],percentile=NA, upper=T,tinypeak.removal=1/50,alpha=0.001,verbose=F)
  quad.2 <- flowDensity(x, channels.ind[c("CD19","CD3")],position=c(T,F),gates= c(cd19.hi,NA))
  quad.3 <- flowDensity(x, channels.ind[c("CD19","CD3")],position=c(F,T),gates= c(cd19.gate,NA))
  return(list(quad2=quad.2,quad3=quad.3,flag=flag))
})

fs.tcell <- lapply(Bcells, function(x) x$quad3@flow.frame)
fs.tcell <- as(object=fs.tcell, Class="flowSet")
#*****************************************************************************************
#                      Gating CD4
#*****************************************************************************************

CD4.8 <- fsApply_pb(fs.tcell, function(x) {
  message("Gating CD4/8")
  cd4.gate <- deGate(x,channels.ind["CD4"],verbose=F)
  temp <- flowDensity(x, channels.ind[c("CD4","CD8")],position=c(T,NA),gates=c(cd4.gate,NA))
  cd8.gate <- deGate(temp@flow.frame,channels.ind["CD8"],use.upper=T,upper=T,percentile=NA,tinypeak.removal = .9,verbose=F)*1.05
  temp.2 <- flowDensity(x, channels.ind[c("CD4","CD8")],position=c(F,NA),gates=c(cd4.gate,NA))
  if(length(getPeaks(temp.2@flow.frame,channels.ind["CD8"],tinypeak.removal = 1/40)$Peaks)>1)
    cd8.gate.lo<-deGate(temp.2@flow.frame,channels.ind["CD8"],tinypeak.removal = 1/40,percentile=.07,verbose=F)
  else
    cd8.gate.lo <- min(deGate(temp.2@flow.frame,channels.ind["CD8"],tinypeak.removal = 1/40,percentile=.07,verbose=F),
                       deGate(temp.2@flow.frame,channels.ind["CD8"],percentile=NA,tinypeak.removal =.5, upper=F,alpha=.005,verbose=F))
  cd4<-flowDensity(x, channels.ind[c("CD4","CD8")],position=c(T,F),gate=c(cd4.gate,cd8.gate))
  dn<-flowDensity(x, channels.ind[c("CD4","CD8")],position=c(F,F),gate=c(cd4.gate*.95,cd8.gate.lo*.98))
  cd8 <- flowDensity(x, channels.ind[c("CD4","CD8")],position=c(F,T),gates=c(cd4.gate,cd8.gate.lo))
  return(list(cd4=cd4,cd8=cd8, dn=dn))
})
#*****************************************************************************************
#                      Gating CD14
#*****************************************************************************************
cd14.gate <- averageGates(fsApply(fs.pbmc,deGate, channels.ind["CD14"],percentile=NA,tinypeak.removal=.9, upper=T),sd.coeff = 1.5)
names(cd14.gate) <- sampleNames(fs.lymph)

CD14 <- fsApply_pb(fs.pbmc, function(x) {
  message("Gating CD14/64")
  x2 <- rotate.data(x,chans = channels.ind[c("CD16","CD14")],theta = -pi/18)$data
  cd14.hi <- deGate(x2,channels.ind["CD14"],percentile=NA,tinypeak.removal = .99,upper=T,alpha=.05,verbose=F)
  cd14.gate <- deGate(x2,channels.ind["CD14"],percentile=NA,verbose=F)
  cd14.gate <- min(cd14.gate,cd14.hi)
  filt<- flowDensity(x2, channels.ind[c("CD16","CD14")],position=c(NA,T),gates=c(NA,cd14.gate))@filter
  filt <- rotate.data(filt,theta=pi/18)$data
  cd14<- flowDensity(x, channels.ind[c("CD16","CD14")],position=c(NA,T),filter=filt)
  return(cd14)
})
fs.14 <-lapply(CD14, function(x) x@flow.frame)
fs.14 <- as(object=fs.14, Class="flowSet")
#*****************************************************************************************
#                      Gating CD14
#*****************************************************************************************
cd16.hi <- fsApply(fs.14,deGate,channels.ind["CD16"],upper=T,percentile=NA,tinypeak.removal=.8,alpha=.09,verbose=F)
names(cd16.hi) <- sampleNames(fs.14)
CD14.hi <- fsApply_pb(fs.14, function(x) {
  message("Gating CD14++")
  
  temp <- flowDensity(x, channels.ind[c("CD16","CD14")],position=c(F,NA),gates=c(cd16.hi[identifier(x)],NA))
  cd14.gate <- deGate(temp@flow.frame, channels.ind["CD14"],upper=F,percentile=NA,tinypeak.removal = .99,alpha=.1,verbose=F)
  quad.1 <- flowDensity(x, channels.ind[c("CD16","CD14")],position=c(F,F),gates=c(cd16.hi[identifier(x)],cd14.gate*1.02))
  quad.2 <- flowDensity(x, channels.ind[c("CD16","CD14")],position=c(T,F),gates=quad.1@gates)
  quad.3 <- flowDensity(x, channels.ind[c("CD16","CD14")],position=c(F,T),gates=quad.1@gates)
  quad.4 <- flowDensity(x, channels.ind[c("CD16","CD14")],position=c(T,T),gates=quad.1@gates)
  flag <-0
  if (quad.4@proportion>20)
    flag<-1
  return(list(quad1=quad.1,quad2=quad.2,quad3=quad.3, quad4=quad.4,flag=flag))
})

#*****************************************************************************************
#                      Gating CD56
#*****************************************************************************************
cd3.gate <- fsApply(fs.lymph,deGate,channels.ind["CD3"],verbose=F)
names(cd3.gate) <- sampleNames(fs.lymph)

CD56 <- fsApply_pb(fs.lymph, function(x) {
  message("Gating NK/NKT")
  
  temp <- flowDensity(x,channels=channels.ind[c('CD56','CD3')],position=c(NA,F),gates=c(NA, cd3.gate[identifier(x)]))
  cd56.gate <- deGate(temp@flow.frame,channels.ind['CD56'],percentile=NA,bimodal = T,sd.threshold = T,verbose=F)
  if (cd56.gate<deGate(temp@flow.frame,channels.ind['CD56'],percentile=.4,use.percentile = T,verbose=F))
    cd56.gate <- deGate(temp@flow.frame,channels.ind['CD56'],percentile=NA,after.peak = T,sd.threshold = T,verbose=F)
  temp <- flowDensity(x,channels=channels.ind[c('CD56','CD3')],position=c(F,T), gates=c(cd56.gate, cd3.gate[identifier(x)]))
  cd56.gate.lo <- deGate(temp@flow.frame,channels.ind['CD56'],upper=T,percentile=NA,tinypeak.removal = .99,alpha=.01,verbose=F)
  if (cd56.gate.lo>cd56.gate)
    cd56.gate.lo<-cd56.gate
  quad.2 <- flowDensity(x,channels=channels.ind[c('CD56','CD3')],position=c(T,F), gates=c(cd56.gate, cd3.gate[identifier(x)]))
  quad.4 <- flowDensity(x,channels=channels.ind[c('CD56','CD3')],position=c(T,T), gates=c(cd56.gate.lo, cd3.gate[identifier(x)]))
  cd56.gate.hi <- deGate(quad.4@flow.frame,channels.ind['CD56'],percentile =.99995,use.percentile=T)
  return(list(nk=quad.2, nkt=quad.4,gate=cd56.gate.hi))
})
fs.nk <-lapply(CD56, function(x) x$nk@flow.frame)
fs.nk <- as(object=fs.nk, Class="flowSet")
#*****************************************************************************************
#                      Gating CD16/56
#*****************************************************************************************

CD56.nk <- fsApply_pb(fs.nk, function(x) {
  message("Gating CD56")
  cd56.gate <- deGate(x,channels.ind["CD56"],twin.factor = .7,percentile = NA,magnitude = .5,upper=T,alpha=.1,verbose=F)
  if(cd56.gate<getPeaks(x,channels.ind["CD56"],tinypeak.removal = .9)$Peaks)
    cd56.gate <-NA
  cd56.upper<- deGate(x,channels.ind["CD56"],percentile=NA,upper=T,use.upper=T,tinypeak.removal = .99,alpha=.99,magnitude=.3,verbose=F)
  cd56.gate <- min(cd56.gate,cd56.upper*1.03,na.rm = T)
  if (cd56.gate>3.9)
    cd56.gate <- deGate(x,channels.ind["CD56"],percentile=NA,upper=T,use.upper=T,tinypeak.removal = .99,alpha=.99,magnitude=.5,verbose=F)

  cd56.hi <-  flowDensity(x,channels=channels.ind[c("CD56","CD16")],position=c(T,NA),gates=c(cd56.gate,NA))

  return(cd56.hi)
})
fs.56 <- lapply(CD56.nk, function(x) x@flow.frame)
fs.56 <- as(object=fs.56, Class="flowSet")
#*****************************************************************************************
#                      Gating CD16/64
#*****************************************************************************************
if (type %in% c("healthy volunteer","post HSCT"))
{
  cd64.gate <- fsApply(fs.14,deGate,channels.ind["CD64"],tinypeak.removal=.5,percentile=NA,upper=F)
  names(cd64.gate)  <- sampleNames(fs.14)
 
  CD64.14 <- fsApply_pb(fs.14, function(x) {
    message("Gating CD64/CD16")
    
    cd64 <- flowDensity(x,channels=channels.ind[c("CD16","CD64")],position=c(T,T),gates=c(cd16.hi[identifier(x)],cd64.gate[identifier(x)]))
    return(cd64)
  })
}
#****************************************************************************************
#****************************************************************************************

CD56.16 <- fsApply_pb(fs.56, function(x) {
  message("Gating 56+16+")
  cd16.neg <- flowDensity(x,channels=channels.ind[c("CD56","CD16")],position=c(NA,F),gates=c(NA,cd16.hi[identifier(x)]))
  return(cd16.neg)
})

  if (type %in% c("healthy volunteer","post HSCT"))
  {
    cell.prop <-sapply(1:length(Lymph), function(fd){
      return(c(Lymph[[fd]]@proportion,
                                    Bcells[[fd]]$quad2@proportion,Bcells[[fd]]$quad3@proportion,
                                    CD4.8[[fd]]$dn@proportion,CD4.8[[fd]]$cd4@proportion,CD4.8[[fd]]$cd8@proportion,
                                    CD14[[fd]]@proportion, CD14.hi[[fd]]$quad2@proportion,CD14.hi[[fd]]$quad3@proportion,
                                    CD14.hi[[fd]]$quad4@proportion,CD56[[fd]]$nk@proportion,CD56[[fd]]$nkt@proportion,
                                    CD56.nk[[fd]]@proportion, CD64.14[[fd]]@proportion,CD56.16[[fd]]@proportion,Bcells[[fd]]$flag))})
  }else{
    cell.prop <-sapply(1:length(Lymph), function(fd){
      return(c(Lymph[[fd]]@proportion,
               Bcells[[fd]]$quad2@proportion,Bcells[[fd]]$quad3@proportion,
               CD4.8[[fd]]$dn@proportion,CD4.8[[fd]]$cd4@proportion,CD4.8[[fd]]$cd8@proportion,
               CD14[[fd]]@proportion, CD14.hi[[fd]]$quad2@proportion,CD14.hi[[fd]]$quad3@proportion,
               CD14.hi[[fd]]$quad4@proportion,CD56[[fd]]$nk@proportion,CD56[[fd]]$nkt@proportion,
               CD56.nk[[fd]]@proportion,CD56.16[[fd]]@proportion,Bcells[[fd]]$flag))})
  }
  
if (type %in% c("healthy volunteer","post HSCT"))
{
  rownames(cell.prop) <- c("Lymphocytes","Bcells","Tcells","DN T cells","CD4+ T cells","CD8+ T cells","CD14+",
                           "CD14+CD16+","CD14++CD16-","CD14++CD16+","NK cells","NKT cells","CD56bright NK cells","CD64++CD16+","CD56bright CD16-","Bcell count Flag")
}else{
  rownames(cell.prop) <- c("Lymphocytes","Bcells","Tcells","DN T cells","CD4+ T cells","CD8+ T cells","CD14+",
                           "CD14+CD16+","CD14++CD16-","CD14++CD16+","NK cells","NKT cells","CD56bright NK cells","CD56bright CD16-","Bcell count Flag")
}
colnames(cell.prop) <- sampleNames(fs)
write.csv(cell.prop,paste0(res.path,"/",panel,"_Proportions_",format(Sys.time(), "%Y-%m-%d"), ".csv"))
