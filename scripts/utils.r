
## PDB to dynamic network function
pdb2mat2dynnet <- function(seqq, xcoord, ycoord, zcoord, cutoff, naa, outnetDirectory){

  uniqueids <- unique(seqq)
  lastPositions <- length(seqq)-match(unique(seqq),rev(seqq))+1
  corData <- data.frame(X=xcoord,Y=ycoord,Z=zcoord)

  # get pairwise 3d euclidean distance
  allDist <- as.matrix(dist(corData))

  # select unique pairs for a given sequence position
  # keeping only the minimum position-pair-wise distances
  amat <- matrix(0,nrow=length(uniqueids),ncol=length(uniqueids))
  start <- 1
  loop <- length(uniqueids)-1


  for(k in 1:loop){
    
    i <- k+1
    
    for(j in i:length(uniqueids)){
      
      startc <- lastPositions[j-1]+1
      endc <- lastPositions[j]
      
      startr <- start
      endr <- lastPositions[k]
      
      tempmat <- allDist[startr:endr,startc:endc]
      amat[k,j] <- min(tempmat)
      amat[j,k] <- min(tempmat)
    }
    
    start <- lastPositions[k]+1
    
  }

  colnames(amat) <- uniqueids
  rownames(amat) <- uniqueids


  #### dynamic networks
  diag(amat) <- diag(amat)+100

  amat[amat <= cutoff] <- 1
  amat[amat > cutoff] <- 0

  seqq <- rownames(amat)
  subs <- seq(0, length(seqq), naa)
  if(length(seqq) %in% subs){loop <- subs}else{
    loop <- subs
    loop <- c(loop, length(seqq))
  }

  count <- 1

  for(k in 2:length(loop)){

    wh <- which(rownames(amat) %in% seqq[seq(1,as.numeric(loop[k]))])
    ttadm <- amat[wh,wh]
    g  <- graph.adjacency(t(ttadm), mode="undirected")
    subg <- get.data.frame(g)
    fwrite(subg, paste0(outnetDirectory,'/',count,'.txt'),quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)
    count  <- count + 1
  }

}

pdb2mat2dynnetpct <- function(seqq, xcoord, ycoord, zcoord, cutoff, pct, outnetDirectory){

  uniqueids <- unique(seqq)
  lastPositions <- length(seqq)-match(unique(seqq),rev(seqq))+1
  corData <- data.frame(X=xcoord,Y=ycoord,Z=zcoord)

  # get pairwise 3d euclidean distance
  allDist <- as.matrix(dist(corData))

  # select unique pairs for a given sequence position
  # keeping only the minimum position-pair-wise distances
  amat <- matrix(0,nrow=length(uniqueids),ncol=length(uniqueids))
  start <- 1
  loop <- length(uniqueids)-1


  for(k in 1:loop){
    
    i <- k+1
    
    for(j in i:length(uniqueids)){
      
      startc <- lastPositions[j-1]+1
      endc <- lastPositions[j]
      
      startr <- start
      endr <- lastPositions[k]
      
      tempmat <- allDist[startr:endr,startc:endc]
      amat[k,j] <- min(tempmat)
      amat[j,k] <- min(tempmat)
    }
    
    start <- lastPositions[k]+1
    
  }

  colnames(amat) <- uniqueids
  rownames(amat) <- uniqueids


  #### dynamic networks
  diag(amat) <- diag(amat)+100

  amat[amat <= cutoff] <- 1
  amat[amat > cutoff] <- 0

  seqq <- rownames(amat)
  naa <- floor(length(seqq)/pct)

  count <- 1
  start <- 1
  end <- naa
  loop <- pct-1

  for(k in 1:pct){

      temps <- seqq[start:end]
      wh <- which(rownames(amat) %in% temps)
      ttadm <- amat[wh,wh]
      g  <- graph.adjacency(t(ttadm), mode="undirected")
      subg <- get.data.frame(g)
      fwrite(subg, paste0(outnetDirectory,'/',count,'.txt'),quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)
      count  <- count + 1
      if(k == loop){end <- length(seqq)}else{end <- end+naa}
 
  }

}

## PDB to unweighted network function
pdb2mat2uwnet <- function(seqq, xcoord, ycoord, zcoord, cutoff, naa, outnetDirectory, fname){

  uniqueids <- unique(seqq)
  lastPositions <- length(seqq)-match(unique(seqq),rev(seqq))+1
  corData <- data.frame(X=xcoord,Y=ycoord,Z=zcoord)

  # get pairwise 3d euclidean distance
  allDist <- as.matrix(dist(corData))

  # select unique pairs for a given sequence position
  # keeping only the minimum position-pair-wise distances
  amat <- matrix(0,nrow=length(uniqueids),ncol=length(uniqueids))
  start <- 1
  loop <- length(uniqueids)-1


  for(k in 1:loop){
    
    i <- k+1
    
    for(j in i:length(uniqueids)){
      
      startc <- lastPositions[j-1]+1
      endc <- lastPositions[j]
      
      startr <- start
      endr <- lastPositions[k]
      
      tempmat <- allDist[startr:endr,startc:endc]
      amat[k,j] <- min(tempmat)
      amat[j,k] <- min(tempmat)
    }
    
    start <- lastPositions[k]+1
    
  }

  colnames(amat) <- uniqueids
  rownames(amat) <- uniqueids


  #### unweighted networks
  loop1 <- length(rownames(amat))-1
  loop2 <- loop1+1
  p1 <- c()
  p2 <- c()

  for(k in 1:loop1){
    i <- k+1
    for(j in i:loop2){ 
      if(amat[k,j] <= cutoff){
        p1 <- c(p1, rownames(amat)[k])
        p2 <- c(p2, rownames(amat)[j])
      }
    }
  }

  Data <- data.frame(x=p1, y=p2)
  fwrite(Data,paste0(outnetDirectory,'/',fname,'.txt'), quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)

  command <- paste("./scripts/list2leda.sh",paste0(outnetDirectory,'/',fname,'.txt'),">",paste0(outnetDirectory,'/',fname,'.gw'))
  system(command)

}


## PDB to weighted network function
pdb2mat2wnet <- function(seqq, xcoord, ycoord, zcoord, cutoff, naa, outnetDirectory, fname){

  uniqueids <- unique(seqq)
  lastPositions <- length(seqq)-match(unique(seqq),rev(seqq))+1
  corData <- data.frame(X=xcoord,Y=ycoord,Z=zcoord)

  # get pairwise 3d euclidean distance
  allDist <- as.matrix(dist(corData))

  # select unique pairs for a given sequence position
  # keeping only the minimum position-pair-wise distances
  amat <- matrix(0,nrow=length(uniqueids),ncol=length(uniqueids))
  start <- 1
  loop <- length(uniqueids)-1


  for(k in 1:loop){
    
    i <- k+1
    
    for(j in i:length(uniqueids)){
      
      startc <- lastPositions[j-1]+1
      endc <- lastPositions[j]
      
      startr <- start
      endr <- lastPositions[k]
      
      tempmat <- allDist[startr:endr,startc:endc]
      amat[k,j] <- min(tempmat)
      amat[j,k] <- min(tempmat)
    }
    
    start <- lastPositions[k]+1
    
  }

  colnames(amat) <- uniqueids
  rownames(amat) <- uniqueids


  #### unweighted networks
  loop1 <- length(rownames(amat))-1
  loop2 <- loop1+1
  p1 <- c()
  p2 <- c()
  seq1 <- c()
  seq2 <- c()
  wt <- c()

  for(k in 1:loop1){
    i <- k+1
    for(j in i:loop2){ 
      if(amat[k,j] <= cutoff){

        p1 <- c(p1, rownames(amat)[k])
        p2 <- c(p2, rownames(amat)[j])

        b <- as.numeric(gsub('.*?([0-9]+).*$','\\1',rownames(amat)[k]))
        a <- as.numeric(gsub('.*?([0-9]+).*$','\\1',rownames(amat)[j]))
      
        seq1 <- c(seq1, as.character(a))
        seq2 <- c(seq2, as.character(b))
        wt <- c(wt, as.character(amat[k,j]))

      }
    }
  }

  diff <- sqrt(abs(as.numeric(seq2)-as.numeric(seq1))/as.numeric(wt))

  Data <- data.frame(x=p1, y=p2, z=diff)
  fwrite(Data,paste0(outnetDirectory,'/',fname,'.txt'), quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)

  command <- paste("./scripts/list2leda.sh",paste0(outnetDirectory,'/',fname,'.txt'),">",paste0(outnetDirectory,'/',fname,'.gw'))
  system(command)

}


removeZero <- function(df){
  df[, sapply(df, function(x) sum(x) == 0)]
}

removeZeroVar3 <- function(df){
  df[, !sapply(df, function(x) min(x) == max(x))]
}

# get correlation
getcor = function(m){
    corm = cov(m)
    dcorm= diag(corm)
    for(i in 2:nrow(corm)){
        for(j in 1:(i-1)){
            if(corm[i,j]!=0){
                corm[i,j]=corm[i,j]/sqrt(dcorm[i]*dcorm[j])
            }
        }
    }
    vct = corm[lower.tri(corm)]
    return(vct)
}


# dynamic GDVM to vector
dgdvm2vec <- function(outdir, indir, outfile){

  files = list.files(indir, pattern='*dcgdv*', full.names=TRUE)

  temp <- fread(files[1], sep=' ', header=FALSE)
  loop <- length(temp) - 1

  remove <- seq(1,loop)

  for(k in 1:length(files)){

      f1 <- fread(files[k], sep=' ', header=FALSE)
      f2 <- f1[,V1:=NULL]
      f3 <- removeZero(f2)
      wh <- which(f3 == TRUE)
      remove <- intersect(remove, wh)

  }

  remove <- remove+1


  ofile <- paste0(outdir,'/vec-',basename(outfile))
  if(file.exists(ofile)){file.remove(ofile)}

  xx <- file.create(ofile)
  fileConn <- file(ofile, open='wt')

  for(k in 1:length(files)){

      m = as.matrix(fread(files[k], drop=remove))
      m1= m[,1] #Labels
      m = m[,2:ncol(m)]
      m = apply(m,2,as.numeric)
      cm = getcor(m)
      cm1 <- c(basename(files[k]), cm)
      cm2 <- paste(cm1, collapse='\t')
      writeLines(cm2, fileConn)

  }

  close(fileConn)

}



## eliminate zero variation columns
elimcol <- function(indir, outdir){

  f1 <- fread(indir, sep='\t', header=FALSE)
  label <- f1[[1]]
  f2 <- f1[,V1:=NULL]
  f3 <- removeZeroVar3(f2)
  wh <- which(f3 == TRUE)
  f4 <- f2[, ..wh]
  f5 <- cbind(data.frame(label), f4)
  fwrite(f5, outdir, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)

}



## PCA function
topca <- function(inputfile, variationThreshold, outputpcfile){

  dat <- fread(inputfile, stringsAsFactors=FALSE, header=FALSE, sep='\t')
  x <- as.matrix(dat[, -c(1)])
  noRows <- nrow(x)
  # cat("No of objects:",noRows,"\n")
  noColumns <- ncol(x)
  # cat("No of original dimensions:",noColumns,"\n")
  x.pca <- prcomp(x, scale.=TRUE)
  sumEigen <- sum(x.pca$sdev^2)
  sumEigenCurrent <- 0
  for(i in 1:noColumns)
  {
    sumEigenCurrent <- sumEigenCurrent + x.pca$sdev[i]^2
    # cat(round(sumEigenCurrent*100/sumEigen,digits=0),"% variation ( EV: ",x.pca$sdev[i]^2,") with",i,"dimension\n")
    if(sumEigenCurrent/sumEigen >= variationThreshold)
    {
      break
    }
  }
  noColumnsReduced <- i
  # cat("No of reduced dimensions:",noColumnsReduced,"\n")
  if(noColumnsReduced < 2) {
    noColumnsReduced <- 2
    # cat("No of reduced dimensions is corrected to 2 to ensure cosine similarity-based comparison\n")
  }
  
  ndata <- cbind(dat[,c(1)], x.pca$x[,1:noColumnsReduced])
  fwrite(ndata,file=outputpcfile,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

}


## create feature matrix
tofeaturematrix <- function(indir, outdir, idens){

  pca <- fread(paste0(indir, '/pca-', basename(idens)), sep='\t',colClasses='character', header=FALSE)
  pca$V1 <- substr(pca$V1,1,nchar(pca$V1)-16)

  temp <- fread(idens, sep='\t', header=FALSE, colClasses='character')
  class <- temp[[1]]
  # print(class)
  ord <- unlist(lapply(strsplit(temp[[2]], '[.]'), '[[',1))

  pcaj <- pca[pca$V1 %in% as.factor(ord), ]
  pcaj1 <- pcaj[match(ord, pcaj$V1),]
  psan <- pcaj1$V1
  pcaj1 <- pcaj1[,V1:=NULL]
  pcaj2 <- cbind(data.frame(class), pcaj1)
  pcaj3 <- cbind(psan, pcaj2)

  fwrite(pcaj3, paste0(outdir,'/matrix-pca-', basename(idens)), quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

}



get_sampled_idfile = function(filename){
  df.ids = read.table(filename,stringsAsFactors=F)
  names(df.ids)=c('label','prot')
  df.ids = df.ids[sample(nrow(df.ids)),]
  rownames(df.ids) = NULL
  return(df.ids)
}

get_sampled_idfile_with_labels = function(protfile,labelfile){
  df.prot   = read.table(protfile,stringsAsFactors=F)
  names(df.prot)=c('prot')
  df.labels = read.table(labelfile,stringsAsFactors=F)
  names(df.labels)=c('label','prot')
  df.ids = merge(df.prot,df.labels,by='prot')
  df.ids = df.ids[sample(nrow(df.ids)),]
  rownames(df.ids) = NULL
  return(df.ids)
}

make_partition_list = function(num.folds){
  partitions = list()
  for(i in 1:num.folds){
    partitions[paste0('prot_idx',i)] = list(c())
    partitions[paste0('tr_prot_idx',i)] = list(c())
  }
}

split_indices = function(df.l,partitions, num.folds){
  for(i in 1:num.folds){
    partitions[[paste0('prot_idx',i)]]    = c(partitions[[paste0('prot_idx',i)]], df.l$prot[df.l$idx==i])
    partitions[[paste0('tr_prot_idx',i)]] = c(partitions[[paste0('tr_prot_idx',i)]], df.l$prot[df.l$idx!=i])
  }
  return(partitions)
}

divide_by_labels = function(df.ids,partitions,num.folds){
  labels = unique(df.ids$label)
  for(l in labels){
    df.l = df.ids[df.ids$label==l,]
    df.l$idx = rep(1:num.folds,length.out=nrow(df.l))
    partitions = split_indices(df.l,partitions,num.folds)
  }
  return(partitions)
}

get_partitions = function(df.ids, num.folds){
  partitions = make_partition_list(num.folds)
  partitions = divide_by_labels(df.ids,partitions,num.folds)
  return(partitions)
}


write_partitions = function(name.outputdir, partitions){
  dir.create(name.outputdir,showWarnings=F,recursive=T)
  for(i in 1:10){
    name.outfile = paste0(name.outputdir,names(partitions[i]),'.txt')
    vc.outfile = partitions[[i]]
    write.table(vc.outfile,name.outfile,quote=F,row.names=F,col.names=F)
  }
}


createPartition <- function(file, odir, num.folds){

  df.ids <- get_sampled_idfile(file)
  partitions <- get_partitions(df.ids, num.folds)
  name.outputdir <- paste0(odir,'/')
  dir.create(name.outputdir)
  write_partitions(name.outputdir,partitions)

  tr.files <- list.files(name.outputdir,pattern='tr_*',full.names=T)

  for(tr.file in tr.files){
    df.ids <- get_sampled_idfile_with_labels(tr.file,file)
    partitions <- get_partitions(df.ids, num.folds)
    name.outputdir <- paste0(odir,'/training_partitions/',strsplit(basename(tr.file),'[.]')[[1]][1],'/')
    dir.create(name.outputdir, recursive=T)
    write_partitions(name.outputdir,partitions)
  }

}