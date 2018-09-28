########################################################
########################################################
#Outlier detection and replacements after usingdeGate for all FCS files in the flowSet

averageGates <- function(vec, med =T, global.var= NA,sd.coeff = 2, name_of_gate="-"){
  # Args:
  #   vec: a numeric vector of values
  #   med: to use median
  #   global.var: a number used instead of median
  # Value:
  #   vec with outlying values replaced by the median
  # (http://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation)
  #s = sd(x)/c4(N) is estimator of SD in vector
  vec[which(is.infinite(vec))]<-NA
  vec <- as.vector(vec,mode = "numeric")
  m <- ifelse(med,median(vec,na.rm=T),global.var)
  if ( any(is.na(vec)) ){
    cat("Changing NAs", which(is.na(vec)), "from", name_of_gate, "\n", sep=" ")
    vec[is.na(vec)] <- m
  }
  
  N <- length(vec)
  c4 <- sqrt(2/(N-1)) * gamma(N/2)/gamma((N-1)/2)
  sdev <- sd(vec,na.rm=T)/c4
  outliers <- which(abs(vec - m)/sdev > sd.coeff)
  if (any(outliers)){
    cat("Changing outliers", outliers, "from", name_of_gate, "\n", sep=" ")
    vec[outliers] <- median(vec[-outliers])
  }
  return (vec)
}
########################################################
########################################################
#Rotating a 2-D data
rotate.data <- function(data, chans=NULL, theta=NULL,min.max=F)
{
  if(nrow(data)<3)
  {
    print("Cannot rotate a matrix with less than 3 rows, returning back the original matrix")
    return(list(data=data,theta=theta))
  }else{
    if (class(data)== "flowFrame" & !is.null(chans))
    {
      if(all(is.na(exprs(data)[,1])))
        return("Cannot rotate a flowFrame with all cells NA")
      no.nas <- which(!is.na(exprs(data)[,chans[1]]))
      data.new <- exprs(data)[no.nas,chans]
      if (is.null(theta))
      {
        reg.slope <- atan(lm(data.new[,1] ~ data.new[,2])$coefficients[2])
        slope <-  atan((max(data.new[,2])-min(data.new[,2]))/(max(data.new[,1])-min(data.new[,1])))
        theta <- ifelse(min.max,no = pi/2-reg.slope,yes = slope-pi/2)
      }
      data.new <- data.new %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
      exprs(data)[no.nas,chans] <- data.new
    }else{
      data <- data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
    }
    return(list(data=data,theta=theta))
  }
}

########################################################
########################################################
#Rotating back a flowDensity object
rotate.fd <- function(fd.object, angle)
{
  new.f <-new(Class = "CellPopulation")
  no.na<-which(!is.na(exprs(fd.object@flow.frame)[,1]))
  dat <- fd.object@flow.frame
  temp <-rotate.data(dat[no.na,], fd.object@channels,theta=-angle)$data
  exprs(dat)[no.na,] <-exprs(temp)
  new.f@flow.frame<-dat
  new.f@filter <-rotate.data(fd.object@filter ,theta =-angle)$data
  colnames(new.f@filter)<-colnames(fd.object@filter )
  new.f@channels<-fd.object@channels
  new.f@proportion<-fd.object@proportion
  new.f@cell.count <-fd.object@cell.count
  return(new.f)
}
########################################################
########################################################
#' Estimates a common logicle transformation for a flowSet.
#'
#' Of the negative values for each channel specified, the median of the specified
#' quantiles are used.
#'
#' @param flow_set object of class 'flowSet'
#' @param channels character vector of channels to transform
#' @param m TODO -- default value from .lgclTrans
#' @param q quantile
#' @return TODO
estimateMedianLogicle <- function(flow_set, channels, m = 4.5, q = 0.05) {
  if (!is(flow_set, "flowSet")) {
    stop("flow_set has to be an object of class 'flowSet'")
  }
  if (missing(channels)) {
    stop("Please specify the channels to be logicle transformed")
  }
  indx <- channels %in% unname(colnames(exprs(flow_set[[1]])))
  if (!all(indx)) {
    stop(paste("Channels", channels[!indx], "were not found in flow_set "))
  }
  
  neg_marker_quantiles <- fsApply(flow_set, function(sample) {
    apply(exprs(sample), 2, function(markers) {
      quantile(markers[markers < 0], probs = q)
    })
  })
  # Replaces 'r' in flowCore:::.lgclTrans
  neg_marker_quantiles <- apply(neg_marker_quantiles, 2,
                                median, na.rm = TRUE)[channels]
  
  # In the case that no negative markers are present, we set this quantile to the
  # default value of 1/2.
  neg_marker_quantiles <- replace(neg_marker_quantiles,
                                  is.na(neg_marker_quantiles), 0.5)
  
  # Replaces 't' in flowCore:::.lgclTrans
  max_range <- do.call(rbind, lapply(fsApply(flow_set, range), function(x) {
    x[2, channels]
  }))
  max_range <- apply(max_range, 2, max)
  
  # Replaces 'w' in flowCore:::.lgclTrans
  w <- (m - log10(max_range / abs(neg_marker_quantiles))) / 2
  if (any(w<0))
     w[which(w<0)]<-0.5
  transformation <- lapply(channels, function(channel) {
    transId <- paste(channel, "medianLogicleTransform", sep = "_")
    
    logicleTransform(transformationId = transId, w = w[channel],
                     t = max_range[channel], m = m, a = 0)
  })
  
  transformList(channels, transformation,
                transformationId = "medianLogicleTransform")
}
########################################################
########################################################
#Transformation for a flowSet
estimate.logicle<-function(fs.raw,talk= TRUE,return.set=TRUE,med=TRUE, m=NA,trans.chans=NULL,estimate=T,return.raw=F)
  
{
  temp <- fsApply(fs.raw, function(x){
    inds <- which(is.na(exprs(x)[,1]))
    y<-x
    if (length(inds)>0)
      y <- x[-inds,]
    return(list(frame=y, na.inds=inds))
  })
  names(temp) <- sampleNames(fs.raw)
  fs <- as(lapply(temp, function(x)
    return(x$frame)),"flowSet")
  sampleNames(fs ) <- sampleNames(fs.raw)
  f<-fs[[1]]
  if (is.null(trans.chans))
  {
    log.channels <- paste("P",1:ncol(f),"DISPLAY",sep="")
    trans.chans <- which(f@description[log.channels]=="LOG")
    if ( length(trans.chans) ==0){
      trans.chans <- setdiff (1:length(f@parameters@data$desc),
                              c(grep("fsc", tolower(f@parameters@data$desc)),
                                grep("ssc", tolower(f@parameters@data$desc)),
                                grep("time", tolower(f@parameters@data$desc)) ) )
    }
    if(any(is.null(trans.chans)) | any(is.na(trans.chans)) | (length(trans.chans)==0))
    {
      return(cat("You have to find the channels manually, there's no information on FCS file","\n"))
    }
  }
  if( talk)
    cat("Channels to be transformed ", trans.chans, "\n")
  if(estimate)
  {
    if (med==TRUE){
      lgl <-estimateMedianLogicle(flow_set=fs,channels=colnames(f)[trans.chans])
      
    }else if(!is.na(m)){
      
      lgl <- estimateLogicle(getGlobalFrame(fs),channels=colnames(f)[trans.chans],m=m)
    }else{
      
      lgl<-tryCatch(estimateLogicle(getGlobalFrame(fs),channels=colnames(f)[trans.chans]), error=function(x) {return(1)})
      if (mode(lgl)=="numeric")
        lgl <- estimateLogicle(getGlobalFrame(fs),channels=colnames(f)[trans.chans],type="data")
    }
    if(return.set)
    {
      if(!return.raw)
      {
        return(fsApply(fs,function(x) transform(x,lgl)))
      }else{
        fs.temp <-fsApply(fs.raw,function(x) {
          tr <- transform(fs[[identifier(x)]],lgl);
          exprs(x)[-temp[[identifier(x)]]$na.inds,]<-exprs(tr)
          w<-Make.FCS(markers = as.vector(x@parameters@data[,2]),data = exprs(x))
          return(w)
        })
      }
    }
    else
      
      return(lgl)
  }else{
    print("Logicle should of been used but it is commented out")
    library(Logicle)
    trans.fs <-fsApply(fs,function(frame) {
      lg<-Logicle::create(T=262144, W=0.5)
      x<-exprs(frame)[,trans.chans]
      y<-Logicle::scale(lg, x)
      exprs(frame)[,trans.chans]<-y
      # detach(package:Logicle, unload=TRUE)
      return (frame)
    })
  }
}

########################################################
########################################################

fsApply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   
  
  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- fsApply(X, wrapper, ...)
  close(pb)
  res
}

########################################################
########################################################
removeMargins<- function(f,chans,sens=1, debris=FALSE,neg=500, verbose = T,return.gate=F)
{
  neg <-rep(neg,length(chans))
  data <- exprs(f)
  margins <- c()
  marg.gates<-c()
  if(!debris)
  {
    for(chan in chans)
    {
      if (is.character(chan))
        chan <-which(colnames(f)==chan)
      stain.max <-max(data[,chan],na.rm = T)
      margins <- which ( data[, chan] >= stain.max*sens)
      data <- data[ -margins, ]
      if(verbose == T){print(paste(length(margins), "margin events in",colnames(f)[chan], "will be removed.",sep =" "))}
      marg.gates <- append(marg.gates,stain.max*sens-1 )
    }
    if(return.gate)
      return(marg.gates)
  }else
  {
    for(chan in chans)
    {
      if (is.character(chan))
        chan <-which(colnames(f)==chan)
      stain.min <-min(data[,chan],na.rm=T)
      margins <- which ( data[, chan] <= stain.min*sens)
      data <- data[ -margins, ]
      if(verbose == T){print(paste(length(margins), "debris events in",colnames(f)[chan], "will be removed.",sep =" "))}
    }
  }
  for(i in 1:length(chans))
  {
    if (neg[i]<500)
    {
      ch<-ifelse (is.character(chans[i]),yes = which(colnames(f)==chans[i]),no = chans[i])
      negs <- which ( data[, ch] < neg[i])
      margins <- negs
      if (length(margins)!=0){
        data <- data[ -margins, ]
      }
      if(verbose == T){print(paste(length(margins), "negative events in",colnames(f)[ch], "will be removed.",sep =" "))}
    }
  } 
  exprs(f) <- data
  
  return(f)
}
#########################################
getGlobalFrame <- function(fs, sample.length=NA, all.cells=F){
  if (is(fs, 'flowFrame')){
    return (frame)
  }
  n <- length(fs)
  sample.n <- ifelse(is.na(sample.length),n,sample.length)
  global.frame <- fsApply(fs[sample(n, sample.n)],
                          function(frame){
                            m <- nrow(frame)
                            sample.size <- ifelse (all.cells,yes = m,no =ceiling(m/sample.n ))
                            
                            exprs(frame )<- exprs(frame)[sample(m, sample.size),]
                            return (frame)
                          })
  global.frame <- as(global.frame, 'flowFrame')
  return (global.frame)
}