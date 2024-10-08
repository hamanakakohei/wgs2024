#! /usr/bin/env Rscript

#CNView: a visualization and annotation tool for copy number variation from whole-genome sequencing

# Copyright (c) 2017 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

#Main plotting function, called later by Rscript (see bottom of script)
CNView <- function(chr,start,end,            #region to be plotted
                   sampleID,                 #Character vector of IDs of samples to plot
                   covmatrix,                #Absolute path of coverage matrix. Header with sample IDs required
                   compression="optimize",   #compression factor for rebinning, if desired
                   highlight=NA,             #list of coordinate pairs; intervals to highlight; defaults to query interval; NULL disables
                   highlightcol="gold",      #vector of colors to shade each highlighted interval
                   window=NA,                #distance to append to both sides of input interval for viewing; NA = 61.8% on either side
                   yscale="optimize",        #vector of values to be represented on y axis
                   normDist=5000000,         #distance outside region to normalize (both sides). Must either be int or "genome"
                   subsample=200,            #will only load this many samples into memory; useful to reduce runtime & memory reqs for very large cohorts
                   UCSCtracks=c("Gene",      #append UCSC sequence context information; choose either NULL
                                "SegDup",    #or up to three among "Gene", "Gap", "RepMask", "blacklist", and "SegDup"
                                "Gap"),
                   genesymbols=TRUE,         #print gene symbols below UCSC gene body annotations
                   probs=TRUE,               #option to add CNV probabilities below each highlighted interval
                   gridlines=TRUE,           #option to draw horizontal gridlines behind plot
                   gcex=1,                   #global scaling for all fonts
                   title=NULL,               #option to add custom title. Overrides default
                   panelnames=NA,            #optional vector for custom names printed above each plot; only works for multi-sample plots
                   legend=T,                 #logical option to plot legend
                   output=NULL,              #path to output as pdf. If NULL, will plot to active device
                   plot=TRUE,                #logical option to disable plot step; mandatory for output!=NULL
                   tabix=FALSE,              #logical option to use tabix to index into coverage matrix
                   noUnix=FALSE,             #logical option to specify a non-unix OS (i.e. no awk, needed to read data)
                   returnData=FALSE,         #logical option to return df of all normalized coverage values
                   quiet=FALSE){             #logical option to disable verbose output
  
  ##Sanity check input##
  if(!(is.numeric(c(start,end)) & end > start)){
    stop("INPUT ERROR: Improper input coordinates")}
  if(!(is.character(sampleID))){
    stop("INPUT ERROR: Improper sampleID")}
  if(compression!="optimize"){
    if(!(compression >= 1 & all.equal(compression,as.integer(compression)))){
      stop('INPUT ERROR: compression parameter must either be "optomize" or positive integer > 1')}}
  suppressWarnings(if(!(is.na(highlight)) & !(is.null(highlight))){
    if(length(highlightcol)!=length(highlight)){
      stop("INPUT ERROR: highlightcol must be same length as intervals to highlight")}})
  if(is.na(window)){
    window <- round(0.618*(end-start),0)
  }
  if(normDist!="genome"){
    if(window>normDist | !(all.equal(window,as.integer(window)))){
      stop('INPUT ERROR: window must be a positive, whole number less than normDist')}
  }
  if(!(is.null(UCSCtracks)) & length(as.vector(UCSCtracks)) > 3){
    stop('INPUT ERROR: UCSCtracks must be a vector of no more than three track names (see documentation for options)')}
  if(plot==F && !(is.null(output))){
    stop('INPUT ERROR: plot must be TRUE for output other than NULL')}
  if(file.exists(covmatrix)==F){
    stop('INPUT ERROR: coverage matrix file not found')
  }
  if(noUnix==T){
    warning('noUnix parameter specified as TRUE; operation speed will be substantially slower')
  }
  if(length(sampleID)>subsample){
    subsample <- length(sampleID)
  }
  
  ##Replace All Hyphens with Periods in Sample IDs##
  sampleID <- gsub("-",".",sampleID)
  
  ##Set preferences##
  options(scipen=1000, #disables scientific notation
          warn=-1) #disables warnings
  
  ##Counts number of samples to plot##
  nsamp=length(sampleID)
  
  ##Interleaver Helper FX##
  interleave <- function(v1,v2){
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
  }
  
  ##Loads required packages##
  require(RMySQL)
  require(plyr)
  require(MASS)
  require(metap)
  require(Rsamtools)
  
  ##Parameter cleanup##
  if(!(is.null(highlight))){
    if(is.na(highlight)){
      highlight <- list(c(start,end))
      highlightcol <- "gold"
    }
  }
  if(!(is.null(UCSCtracks))){
    UCSCtracks <- rev(UCSCtracks)
  }
  
  ##Subset & Load Plotting Values##
  if(quiet==F){cat("Filtering & loading coverage matrix...")}
  if(noUnix==TRUE){
    cov <- read.table(covmatrix,header=T,sep="\t",check.names=F,comment.char="")
    cov <- cov[which(cov[,1]==chr & cov[,2]<=end & cov[,3]>=start),]
  }else{
    if(normDist!="genome"){
      subcovmatrix <- tempfile()
      if(tabix==T){
        tab.gr <- GRanges(chr,IRanges(start=max(start-normDist,0),end=end+normDist))
        tab.file <- open(TabixFile(file=covmatrix))
        header <- unlist(strsplit(headerTabix(tab.file)$header,split="\t"))
        header[1:3] <- c("chr","start","end")
        close(tab.file)
        tab.file <- open(TabixFile(file=covmatrix))
        cov.dat <- lapply(unlist(scanTabix(tab.file,param=tab.gr)),strsplit,split="\t")
        cov.dat <- lapply(cov.dat,function(vals){
          coords <- unlist(vals)[1:3]
          cov.vals <- as.numeric(unlist(vals)[-c(1:3)])
          return(c(coords,cov.vals))
        })
        cov <- as.data.frame(t(matrix(as.vector(unlist(cov.dat)),ncol=length(cov.dat))))
        colnames(cov) <- header
        rownames(cov) <- NULL
        cov[,-1] <- apply(cov[,-1],2,as.numeric)
        close(tab.file)
      }else{
        if(summary(file(covmatrix))$class != "gzfile"){
          system(paste("head -n1 ",covmatrix," > ",subcovmatrix,sep=""))
          system(paste("awk -v OFS=\"\t\" '{ if ($1==\"",chr,"\" && $2<=",end+normDist," && $3>=",start-normDist,") print $0 }' ",covmatrix," >> ",
                       subcovmatrix,sep=""))
        }else{
          system(paste("zcat ",covmatrix," | head -n1 > ",subcovmatrix,sep=""))
          system(paste("zcat ",covmatrix," | awk -v OFS=\"\t\" '{ if ($1==\"",chr,"\" && $2<=",end+normDist," && $3>=",start-normDist,") print $0 }' >> ",
                       subcovmatrix,sep=""))          
        }
      }
    }else{
      subcovmatrix <- covmatrix
    }
    if(tabix==F){
      cov <- read.table(subcovmatrix,header=T,sep="\t",check.names=F,comment.char="")
    }
  }
  if(quiet==F){cat(" Complete\n")}
  
  ##Drop Columns to Specified Sample Size##
  cov <- cov[,unique(c(1:3,as.vector(sapply(head(unique(c(sampleID,sample(names(cov[,-c(1:3)])))),n=subsample),function(val){which(val==colnames(cov))}))))]
  
  ##Rebin Helper FX##
  rebin <- function(df,compression){
    Chr <- df[1,1]
    Start <- df[1,2]
    End <- df[compression,3]
    for(i in 2:(floor(nrow(df)/compression))) {
      Chr <- c(Chr,as.character(df[((i-1)*compression)+1,1]))
      Start <- c(Start,as.integer(df[((i-1)*compression)+1,2]))
      End <- c(End,as.integer(df[i*compression,3]))
    }
    newvals <- apply(df[,4:ncol(df)],2,
                     function(vals,compression){
                       newcol <- sum(vals[1:compression])
                       for(i in 2:(floor(length(vals)/compression))) {
                         newcol <- c(newcol,as.integer(sum(vals[(((i-1)*compression)+1):(i*compression)])))
                       }
                       return(newcol)
                     },compression)
    return(as.data.frame(cbind(Chr,Start,End,newvals)))
  }
  
  ##Rebins values##
  if(quiet==F){cat("Compressing coverage matrix [")}
  obinsize <- cov[1,3]-cov[1,2]
  if(compression=="optimize"){
    compression <- round(((end-start+(2*window))/160)/obinsize,0)
    #compression <- round(((end-start+(2*window))/80)/obinsize,0) # hamanaka
  }
  if(compression<1){
    compression <- 1
  }
  #binsize <- 100000
  binsize <- compression*obinsize # hamanaka
  if(quiet==F){cat(paste(prettyNum(binsize,big.mark=",")," bp bins]... ",sep=""))}
  if(compression>1){
    res <- rebin(cov,compression)    
  } else {
    res <- cov
  }
  if(quiet==F){cat(" Complete\n")}
  
  ##Scale each col within that sample by median##
  if(quiet==F){cat("Performing intra-sample normalization...")}
  res[,4:ncol(res)] <- apply(res[,4:ncol(res)],2,
                             function(vals){
                               nvals <- as.numeric(vals)/median(as.numeric(vals[which(vals>0)]))
                               nvals[is.infinite(nvals)] <- NA
                               return(nvals)
                             })
  if(quiet==F){cat(" Complete\n")}
  
  ##Normalize each row across all samples##
  if(quiet==F){cat("Performing inter-sample normalization...")}
  colnames(res)[1:3] <- c("Chr","Start","End")
  names <- colnames(res)
  oncol <- ncol(res)
  res[,4:oncol] <- data.frame(t(apply(res[,4:oncol],1,function(vals){
    return( as.numeric(vals) / mean(as.numeric(vals)) * 2 ) # hamanaka
    #return(scale(as.numeric(vals)))
  })))
  res[is.na(res)] <- 0
  res$mean <- apply(res[,4:oncol],1,mean)
  res$sd <- apply(res[,4:oncol],1,sd)
  res$median <- apply(res[,4:oncol],1,median)
  res$mad <- apply(res[,4:oncol],1,mad)
  if(quiet==F){cat(" Complete\n")}
  
  ##Subset View Window##
  plotSet <- as.data.frame(apply(res[which(as.integer(as.character(res$Start)) <= end+window & 
                                             as.integer(as.character(res$End)) >= start-window & 
                                             as.character(res$Chr) == as.character(chr)),],
                                 2,function(col){
                                   return(as.numeric(as.character(col)))
                                 }))
  plotSet[is.na(plotSet)] <- 0
  
  ##Get Sample Indexes##
  sampIdx <- as.vector(sapply(as.vector(sampleID),function(val){grep(paste("\\b",as.character(val),"\\b",sep=""),
                                                                     colnames(plotSet),ignore.case=F)}))
  

  ##Output Options##  ##Output Options##
  if(plot==T){
    if(!(is.null(output))){
      if(quiet==F){cat(paste("Plotting samples to ",output,"...",sep=""))}
      if(nsamp==1){
        pdf(output,width=5.1,height=2.5)
      }else{
        pdf(output,width=17,height=(2.7+(1.4*nsamp)))
      }
    } else {
      if(quiet==F){cat("Plotting samples to screen...")}
    }
    par(mar=c(0,0,0,0),
        oma=c(5,5,4,3))
    if(!(is.null(UCSCtracks))){
      layout(matrix(c(1:(nsamp+1)),byrow=T),heights=c(rep(4,nsamp),1))
    } else {
      par(mfrow=c(nsamp,1))
    }
    
    ##Plot per sample##
    for(k in 1:nsamp){
      
      ##Generate Colors##
      if(pt(plotSet[1,sampIdx[k]],df=(ncol(plotSet)-8)) >= 1-(0.05/(ncol(plotSet)-7))){
        colval <- "blue"
      }else if(pt(plotSet[1,sampIdx[k]],df=(ncol(plotSet)-8)) <= 0.05/(ncol(plotSet)-7)){
        colval <- "red"
      }else{
        colval <- "gray40"
      }
      for(i in 2:nrow(plotSet)){
        if(p.adjust(pt(plotSet[i,sampIdx[k]],df=(ncol(plotSet)-8),
                       lower.tail=F),
                    method="fdr",n=ncol(plotSet)-7) <= 0.95){
          colval <- c(colval,"blue")
        } else if(p.adjust(pt(plotSet[i,sampIdx[k]],df=(ncol(plotSet)-8),
                              lower.tail=T),
                           method="fdr",n=ncol(plotSet)-7) <= 0.05){
          colval <- c(colval,"red")
        } else {
          colval <- c(colval,"gray40")
        }
      }
      colval <- interleave(colval,colval)
      for(i in 2:length(colval)){
        if(colval[i]=="blue"){
          colval[i-1] <- "blue"
        }else if(colval[i]=="red"){
          colval[i-1] <- "red"
        }
      }
      
      ####Plot main####
      plot(as.numeric(plotSet$Start),
           as.numeric(plotSet[,sampIdx[k]]),
           type="n",xlim=c(max(start-window,0),
                           max(plotSet$Start)),
           if(yscale=="optimize"){
             aaa = max(abs(min(plotSet[,sampIdx])), abs(max(plotSet[,sampIdx]))) # hamanaka
             ylim=c( 0, aaa + 1 ) # hamanaka
             #ylim=c(min(0,min(plotSet[,sampIdx]))-2, max(0,max(plotSet[,sampIdx]))+2)
           }else{
             ylim=yscale
           },
           xlab="",xaxt="n",ylab="",xaxs="i",yaxt="n",
           panel.first=c(if(gridlines==T){abline(h=seq(round_any(par("usr")[3],2),
                                      round_any(par("usr")[4],2),by=2),
                                col="gray80")},
                         #polygon(y=c(2*plotSet$mad, rev(-2*plotSet$mad)),
                         polygon(y=c(2 + 2*plotSet$sd, 2 + rev(-2*plotSet$sd)), # hamanaka
                                 x=c(as.numeric(as.character(plotSet$Start)),
                                     rev(as.numeric(as.character(plotSet$Start)))),
                                 col="gray85",border="gray60"),
                         #polygon(y=c(plotSet$mad, rev(-plotSet$mad)),
                         polygon(y=c(2 + plotSet$mad, 2 + rev(-plotSet$mad)), # hamanaka
                                 x=c(as.numeric(as.character(plotSet$Start)),
                                     rev(as.numeric(as.character(plotSet$Start)))),
                                 col="gray75",border="gray55"),
                         points(as.integer(as.character(plotSet$Start)),
                                as.numeric(as.character(plotSet$median)),
                                lty=2,type="l"),
                         if(!(is.null(highlight))){
                           for(i in 1:length(highlight)){
                             rect(highlight[[i]][1],
                                  par("usr")[3],
                                  highlight[[i]][2],
                                  par("usr")[4],
                                  col=adjustcolor(highlightcol[i],alpha=0.2),
                                  border=NA)
                           }
                         },
                         segments(x0=interleave(plotSet[seq(1,(nrow(plotSet)-1)),2],
                                                plotSet[seq(2,(nrow(plotSet))),2]),
                                  y0=interleave(plotSet[seq(1,(nrow(plotSet)-1)),sampIdx[k]],
                                                plotSet[seq(1,(nrow(plotSet)-1)),sampIdx[k]]),
                                  x1=interleave(plotSet[seq(2,(nrow(plotSet))),2],
                                                plotSet[seq(2,(nrow(plotSet))),2]),
                                  y1=interleave(plotSet[seq(1,(nrow(plotSet)-1)),sampIdx[k]],
                                                plotSet[seq(2,(nrow(plotSet))),sampIdx[k]]),
                                  lwd=3,
                                  col='gray40'),
                                  #col=colval), # hamanaka
                         if(!(is.null(highlight))){
                           for(i in 1:length(highlight)){
                             abline(v=highlight[[i]][1],
                                    lty=3)
                             abline(v=highlight[[i]][2],
                                    lty=3)
                           }
                         }))
      #Add probabilites, if optioned. t-test for each bin then Fisher's method to combine p-values
      if(probs==T){
        if(!(is.null(highlight))){
          for(i in 1:length(highlight)){
            pDEL <- p.adjust(sumlog(pt(plotSet[which(plotSet$End>=highlight[[i]][1] & 
                                                       plotSet$Start<=highlight[[i]][2]),
                                               sampIdx[k]],
                                       df=ncol(plotSet)-7,lower.tail=T))$p,
                             method="fdr",
                             n=length(plotSet[which(plotSet$End>=highlight[[i]][1] & 
                                                      plotSet$Start<=highlight[[i]][2]),
                                              sampIdx[k]]))
            pDEL <- paste(round(as.numeric(strsplit(format(pDEL,scientific=T),split="e")[[1]][1]),3),
                          "E",strsplit(format(pDEL,scientific=T),split="e")[[1]][2],sep="")
            pDUP <- p.adjust(sumlog(pt(plotSet[which(plotSet$End>=highlight[[i]][1] & 
                                                       plotSet$Start<=highlight[[i]][2]),
                                               sampIdx[k]],
                                       df=ncol(plotSet)-7,lower.tail=F))$p,
                             method="fdr",
                             n=length(plotSet[which(plotSet$End>=highlight[[i]][1] & 
                                                      plotSet$Start<=highlight[[i]][2]),
                                              sampIdx[k]]))
            pDUP <- paste(round(as.numeric(strsplit(format(pDUP,scientific=T),split="e")[[1]][1]),3),
                          "E",strsplit(format(pDUP,scientific=T),split="e")[[1]][2],sep="")
            text(x=mean(highlight[[i]]),y=par("usr")[3],pos=3,cex=gcex,
                 labels=paste("q(Del) = ",pDEL,"\nq(Dup) = ",pDUP,sep=""))
          }
        }
      }
      #Y Axis
      axis(2,at=seq(round_any(par("usr")[3],2),
                    round_any(par("usr")[4],2),by=2),
           las=2,cex.axis=gcex)
      mtext(paste("Norm. Depth t Score",sep=""),
            side=2,line=2.3,cex=gcex)
      #Print Sample ID if >1 sample
      if(nsamp>1){
        if(is.na(panelnames)){
          text(x=mean(par("usr")[1:2]),y=par("usr")[4],
               labels=names(plotSet)[sampIdx[k]],
               cex=1.2*gcex,pos=1,font=2)
        }else{
          text(x=mean(par("usr")[1:2]),y=par("usr")[4],
               labels=panelnames[k],
               cex=1.2*gcex,pos=1,font=2)
        }
      }
      
      ##Title if first sample
      if(k==1){
        #Title
        if(!(is.null(title))){
          mtext(text=title,outer=T,side=3,line=1,font=2,cex=1.1*gcex)
        }else{
          if(nsamp==1){
            mtext(text=paste("Normalized Sequencing Depth of ",sampleID,sep=""),
                  outer=T,side=3,line=1,font=2,cex=1.1*gcex)
          }else{
            mtext(text=paste("Normalized Sequencing Depth of ",nsamp," Samples",sep=""),
                  outer=T,side=3,line=1,font=2,cex=1.1*gcex)
          }
        }
        mtext(text=paste("chr",chr," : ",prettyNum(max((start-window),0),big.mark=",")," - ",
                         prettyNum(end+window,big.mark=","),sep=""),
              outer=T,side=3,line=0,cex=0.7*gcex)
        #Legend & resolution
        if(nsamp==1){
          rcex=1
        }else{
          rcex=1.3
        }
        if(max(plotSet[,sampIdx])+min(plotSet[,sampIdx]) >= 0){
          text(x=par("usr")[1],
               y=0.9*par("usr")[4],
               labels=paste(prettyNum(binsize,big.mark=",")," bp Bins",sep=""),
               font=4,pos=4,cex=rcex*gcex*0.1)
               #font=4,pos=4,cex=rcex*gcex) #hamanaka
        }else{
          text(x=par("usr")[1],
               y=0.9*par("usr")[3],
               labels=paste(prettyNum(binsize,big.mark=",")," bp Bins",sep=""),
               font=4,pos=4,cex=rcex*gcex*0.1)
               #font=4,pos=4,cex=rcex*gcex) #hamanaka
        }
        if(legend==T){
          if(nsamp==1){
            lcex=0.6
          }else{
            lcex=0.8
          }
          if(max(plotSet[,sampIdx])+min(plotSet[,sampIdx]) >= 0){
            legend("topright",
                   legend=c(paste("q(Dup) < 0.05 (df=",ncol(plotSet)-8,")",sep=""),
                            paste("q(Del) < 0.05 (df=",ncol(plotSet)-8,")",sep=""),
                            "Median t Score",
                            "1 * MAD",
                            "2 * MAD"),
                   pch=c(NA,NA,NA,15,15),pt.cex=c(1,1,1,1.5,1.5)*gcex,
                   lty=c(1,1,2,NA,NA),lwd=c(4,4,1,NA,NA),
                   col=c("blue","red","black","gray54","lightgray"),
                   bg="white",cex=lcex*gcex)
          } else if(max(plotSet[,sampIdx])+min(plotSet[,sampIdx]) < 0){
            legend("bottomright",
                   legend=c(paste("q(Dup) < 0.05 (df=",ncol(plotSet)-8,")",sep=""),
                            paste("q(Del) < 0.05 (df=",ncol(plotSet)-8,")",sep=""),
                            "Median t Score",
                            "1 * MAD",
                            "2 * MAD"),
                   pch=c(NA,NA,NA,15,15),pt.cex=c(1,1,1,1.5,1.5)*gcex,
                   lty=c(1,1,2,NA,NA),lwd=c(4,4,1,NA,NA),
                   col=c("blue","red","black","gray54","lightgray"),
                   bg="white",cex=lcex*gcex)
          }
        }
      }
    }
    if(quiet==F){cat(" Complete\n")}
    
    ##UCSC plot##
    if(nsamp==1){
      UCSClabcex=0.7
    }else{
      UCSClabcex=1
    }
    if(!(is.null(UCSCtracks))){
      ##Connects to UCSC##
      if(quiet==F){cat("Attempting to connect to UCSC Genome Browser...")}
      UCSC <- dbConnect(MySQL(),
                        user='genome',
                        dbname='hg19',
                        host='genome-mysql.cse.ucsc.edu')
      if(quiet==F){cat(" Complete\n")}
      
      if(quiet==F){cat("Appending UCSC tracks...")}
      plot(plotSet$Start,
           c(1,2,rep(3,nrow(plotSet)-2)),
           ylim=c(0,3),
           type="n",xaxt="n",yaxt="n",
           ylab="",xlab="",xaxs="i")
      if("Gap" %in% UCSCtracks){
        gaps <- dbGetQuery(UCSC,paste("SELECT chromStart, chromEnd, type FROM gap WHERE `chrom` = 'chr",chr,"' ",
                                      "AND `chromStart` <= ",end+window," ",
                                      "AND `chromEnd` >= ",start-window,sep=""))
        if(nrow(gaps) > 0){
          for(i in 1:nrow(gaps)){
            rect(xleft=gaps$chromStart[i],
                 ybottom=grep("Gap",UCSCtracks)-.9,
                 xright=gaps$chromEnd[i],
                 ytop=grep("Gap",UCSCtracks)-.1,
                 col="cadetblue3",border=NA)
          }
        }
      }
      if("SegDup" %in% UCSCtracks){
        segDups <- dbGetQuery(UCSC,paste("SELECT chromStart, chromEnd FROM genomicSuperDups WHERE `chrom` = 'chr",chr,"' ",
                                         "AND `chromStart` <= ",end+window," ",
                                         "AND `chromEnd` >= ",start-window,sep=""))
        if(nrow(segDups) > 0){
          for(i in 1:nrow(segDups)){
            rect(xleft=segDups$chromStart[i],
                 ybottom=grep("SegDup",UCSCtracks)-.9,
                 xright=segDups$chromEnd[i],
                 ytop=grep("SegDup",UCSCtracks)-.1,
                 border=NA,col="darkorange")
          }
        }
      }
      if("RepMask" %in% UCSCtracks){
        repeatMasker <- dbGetQuery(UCSC,paste("SELECT genoStart, genoEnd, repClass FROM rmsk WHERE `genoName` = 'chr",chr,"' ",
                                              "AND `genoStart` <= ",end+window," ",
                                              "AND `genoEnd` >= ",start-window,sep=""))
        if(nrow(repeatMasker) > 0){
          for(i in 1:nrow(repeatMasker)){
            rect(xleft=repeatMasker$genoStart[i],
                 ybottom=grep("RepMask",UCSCtracks)-.9,
                 xright=repeatMasker$genoEnd[i],
                 ytop=grep("RepMask",UCSCtracks)-.1,
                 border=NA,col="lightsteelblue4")
          }
        }
      }
      if("Gene" %in% UCSCtracks){
        genes <- unique(dbGetQuery(UCSC,paste("SELECT txStart, txEnd, name2, strand, exonStarts, exonEnds FROM refGene WHERE `chrom` = 'chr",chr,"' ",
                                              "AND `txStart` <= ",end+window," ",
                                              "AND `txEnd` >= ",start-window,sep="")))
        if(genesymbols==T){
          genesqueeze=0.6
        }else{
          genesqueeze=0.9
        }
        if(nrow(genes) > 0){
          for(i in 1:nrow(genes)){
            rect(xleft=genes$txStart[i],
                 ybottom=grep("Gene",UCSCtracks)-genesqueeze,
                 xright=genes$txEnd[i],
                 ytop=grep("Gene",UCSCtracks)-.1,
                 border=NA,col="lightgreen")
            segments(x0=genes$txStart[i],
                     x1=genes$txEnd[i],
                     y0=mean(c(grep("Gene",UCSCtracks)-.1,
                               grep("Gene",UCSCtracks)-genesqueeze)),
                     y1=mean(c(grep("Gene",UCSCtracks)-.1,
                               grep("Gene",UCSCtracks)-genesqueeze)),
                     col="darkgreen")
          }
          exons <- unique(data.frame("start"=as.numeric(as.character(unlist(sapply(genes$exonStarts,
                                                                                   function(string){return(strsplit(string,split=","))})))),
                                     "end"=as.numeric(as.character(unlist(sapply(genes$exonEnds,
                                                                                 function(string){return(strsplit(string,split=","))}))))))
          
          for(i in 1:nrow(exons)){
            rect(xleft=exons[i,1],
                 xright=exons[i,2],
                 ybottom=grep("Gene",UCSCtracks)-(genesqueeze-0.05),
                 ytop=grep("Gene",UCSCtracks)-0.15,
                 border=NA,col="darkgreen")
          }
          if(genesymbols==T){
            for(i in unique(genes$name2)){
              text(x=(max(c(min(genes[which(genes$name2==i),1]),par("usr")[1]))+min(c(max(genes[which(genes$name2==i),2]),par("usr")[2])))/2,
                   y=grep("Gene",UCSCtracks)-.8,
                   labels=i,
                   cex=0.75*gcex,font=4)
            }
          }
        }
      }
      axis(2,at=c(seq(0.5,length(UCSCtracks)-0.5)),
           labels=UCSCtracks,cex.axis=UCSClabcex,las=1,tick=F)
      if(quiet==F){cat(" Complete\n")}
    }
    
    ##Adds X axis##
    axis(1,at=seq(min(plotSet$Start),
                  max(plotSet$Start),
                  by=(max(plotSet$Start)-min(plotSet$Start))/8),
         cex.axis=UCSClabcex,
         labels=prettyNum(seq(min(plotSet$Start),
                              max(plotSet$Start),
                              by=(max(plotSet$Start)-min(plotSet$Start))/8),
                          big.mark=","))
    mtext(paste("chr",chr," Coordinate (bp)",sep=""),
          side=1,outer=T,line=2)
    
    ##Close Output##
    if(!(is.null(output))){
      dev.off()
    }
  }
  
  ##Disconnect from UCSC, if necessary##
  if(exists("UCSC")){
    dbDisconnect(UCSC)
  }
  
  ##Return Values as df if specified#
  if(returnData==T){
    if(quiet==F){cat("Returning normalized values\n")}
    if(quiet==F){cat(paste("\n** FINISHED ON ",date()," **\n\n",sep=""))}
    return(plotSet)
  }
  
  #Remove temporary file##
  if(normDist!="genome"){
    system(paste("rm ",subcovmatrix,sep=""))
  }
  
  ##Finish up##
  if(quiet==F){cat(paste("\n** FINISHED ON ",date()," **\n\n",sep=""))}
}

####################################################
####Rscript functionality for command line usage####
####################################################

#Disables factor default
options(stringsAsFactors=F)

#Loads optparse
require(optparse)

#list of Rscript options
option_list <- list(
  make_option(c("-c", "--compression"), type="integer", default=NULL,
              help="compression scalar for rebinning, if desired [default '%default']", 
              metavar="integer"),
  make_option(c("-i","--highlight"), type="character", default=NA,
              help="tab-delimited list of coordinate pairs for intervals to highlight and color as third column; NULL disables highlighting [default %default]",
              metavar="character"),
  make_option(c("-w","--window"), type="integer", default=NULL,
              help="distance to append to both sides of input interval for viewing [default 61.8% of plot interval]",
              metavar="integer"),
  make_option("--ymin", type="integer", default=NULL,
              help="minimum value for y axis [default %default]", 
              metavar="integer"),
  make_option("--ymax", type="integer", default=NULL,
              help="maximum value for y axis [default %default]", 
              metavar="integer"),
  make_option(c("-n","--normDist"), type="integer", default=5000000,
              help="distance outside region to use for normalization (both sides) [default %default]",
              metavar="integer"),
  make_option(c("-s","--subsample"), type="integer", default=200,
              help="truncate coverage matrix to [s] samples; useful for very large cohorts [default %default]",
              metavar="integer"),
  make_option("--gcex", type="integer", default=1,
              help="scalar applied to all fonts and legend [default %default]",
              metavar="integer"),
  make_option("--names", type="character", default=NA,
              help="list of custom names to be applied to each plot panel (e.g. 'mother', 'father', 'child', rather than actual sample IDs) [default %default]",
              metavar="character"),
  make_option(c("-t","--title"), type="character", default=NULL,
              help="custom title for plot [default %default]",
              metavar="character"),
  make_option(c("-p","--probs"), action="store_true", default=FALSE,
              help="add CNV probabilities below each higlighted interval [default %default]"),
  make_option(c("-g","--gridlines"), action="store_true", default=FALSE,
              help="add horizontal gridlines to plot background [default %default]"),
  make_option(c("-u","--noUCSC"), action="store_true", default=FALSE,
              help="disable UCSC track plotting [default %default]"),
  make_option(c("-G","--nogenesymbols"), action="store_true", default=FALSE,
              help="disable gene symbol printing below gene bodies in UCSC tracks [default %default]"),
  make_option(c("--tabix"), action="store_true", default=FALSE,
              help="use tabix to index into coverage matrix [default %default]"),
  make_option(c("--noUnix"), action="store_true", default=FALSE,
              help="disable use of unix coreutils [default %default]"),
  make_option(c("-q","--quiet"), action="store_true", default=FALSE,
              help="disable verbose output [default %default]"),
  make_option(c("-l", "--nolegend"), action="store_false", default=TRUE,
              help="disable legend on plot [default %default]")
)

#Get command-line arguments & options
parser <- OptionParser(usage="%prog [options] chr start end samples.list covmatrix.bed outfile",
                       option_list=option_list,add_help_option=T)
args <- parse_args(parser,positional_arguments=TRUE)
opts <- args$options

#checks for appropriate positional arguments
if(length(args$args) != 6) {
  print_help(parser)
  stop("Incorrect number of required positional arguments")
}

#Prints startup message
if(opts$quiet==F){cat("+-------------------+\n| CNView Visualizer |\n|     (c) 2017      |\n+-------------------+\n")}

#Processes arguments & cleans options
if(file.exists(args$args[4])==T){
  if(opts$quiet==F){cat(paste("Reading sample IDs from ",args$args[4],"\n",sep=""))}
  samps <- read.table(args$args[4])[,1]
  #Add "X" before IDs beginning with number (R doesn't like data frame colnames beginning with numbers)
  samps <- as.vector(unlist(sapply(samps,function(ID){
    if(substr(as.character(ID),0,1) %in% 0:9){
      return(as.character(paste("X",ID,sep="")))
    }else{
      return(ID)
    }
  })))
}else{
  if(opts$quiet==F){cat(paste("Sample ID file '",args$args[4],"' not found, assuming single sample ID provided\n",sep=""))}
  samps <- args$args[4]
  #Add "X" before IDs beginning with number (R doesn't like data frame colnames beginning with numbers)
  if(substr(as.character(samps),0,1) %in% 0:9){
    samps <- paste("X",samps,sep="")
  }
}
#suppressWarnings(if(!(is.na(opts$highlight)) & !(is.null(opts$highlight))){
#  hightable <- read.table(opts$highlight,sep="\t",header=F,comment.char="")[,1:3]
#  highopt <- as.list(as.data.frame(t(hightable[,1:2])))
#  highcolopt <- as.vector(hightable[,3])
#}else{
  highopt <- opts$highlight
  highcolopt <- "gold"
#})
if(!(is.null(opts$ymin)) & !(is.null(opts$ymax))){
  yscale=c(opts$ymin,opts$ymax)
}else{
  yscale="optimize"
}
if(is.null(opts$compression)){
  compression <- "optimize"
}else{
  compression <- opts$compression
}
if(is.null(opts$window)){
  window=round(0.618*(as.numeric(args$args[3])-as.numeric(args$args[2])),0)
}else{
  window=round(opts$window,0)
}
if(opts$noUCSC==TRUE){
  UCSCtracks=NULL
}else{
  UCSCtracks=c("Gene","SegDup","Gap")
}
if(!(is.na(opts$names)) & file.exists(opts$names)==T){
  names <- as.vector(read.table(opts$names,header=F)[,1])
}else{
  names <- opts$names
}

#Runs CNView function
if(opts$quiet==T){
  sink("/dev/null")
}
CNView(chr=args$args[1],
       start=as.numeric(args$args[2]),
       end=as.numeric(args$args[3]),
       sampleID=as.vector(samps),
       covmatrix=args$args[5],
       compression=compression,
       highlight=highopt,
       highlightcol=highcolopt,
       window=window,
       yscale=yscale,
       normDist=opts$normDist,
       subsample=opts$subsample,
       UCSCtracks=UCSCtracks,
       genesymbols=!(opts$nogenesymbols),
       probs=opts$probs,
       gridlines=opts$gridlines,
       gcex=opts$gcex,
       title=opts$title,
       panelnames=names,
       legend=opts$nolegend,
       output=args$args[6],
       plot=T,
       tabix=opts$tabix,
       noUnix=opts$noUnix,
       returnData=F,
       quiet=opts$quiet)
