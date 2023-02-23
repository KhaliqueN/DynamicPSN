
runmaincmd <- function(cutoff, naa, choice, annotationFile, pdbDirectory, partitionFlag, partitionFolder){

	outputDirectory <- 'output'#strsplit(basename(annotationFile),'[.]')[[1]][1]

	#check if there is an output folder
	if(dir.exists(outputDirectory)){unlink(outputDirectory, recursive=TRUE)}

	ann <- fread(annotationFile, sep='\t', header=FALSE)

	pdbID <- unlist(lapply(strsplit(ann[[2]],'[_]'),'[[',1))
	chainID <- unlist(lapply(strsplit(ann[[2]],'[_]'),'[[',2))

	st <- strsplit(ann[[2]],'[+]')
	startID <- list()
	endID <- list()
	for(k in 1:length(st)){
		temp <- st[[k]]
		sid <- c()
		eid <- c()
		for(j in 1:length(temp)){
			temp1 <- strsplit(temp[j],'[_]')
			if(j == 1){
				s1 <- 3
				e1 <- 4
			}else{
				s1 <- 1
				e1 <- 2
			}
			sid <- c(sid, temp1[[1]][s1])
			eid <- c(eid, temp1[[1]][e1])
		}

		startID[[k]] <- sid
		endID[[k]] <- eid
	}


	for(k in 1:length(pdbID)){

		tfile <- readLines(paste0(pdbDirectory,"/",pdbID[k],".cif"))

		# Extract start and end positions of the coordinates entries
		wh <- which(tfile == "loop_")+1
		tfile0 <- trimws(tfile[wh])
		whh1 <- which(tfile0 == "_atom_site.group_PDB")
		start <- wh[whh1]
		end <- wh[whh1+1]-1-2

		# Extract the coordinates part of the PDB file
		tfile <- tfile[start:end]
		lineID <- word(tfile, 1)
		wh <- which(lineID == "ATOM" | lineID == "HETATM")

		# Extract the field entries
		whf <- setdiff(seq(1,length(tfile)), wh)
		fields <- trimws(tfile[whf])

		tfile1 <- trimws(tfile[wh])
		tfile2 <- read.table(textConnection(tfile1), sep='', colClasses = "character")

		#filter using chain
		chainPosition <- which(fields == "_atom_site.auth_asym_id")
		chain <- tfile2[[chainPosition]]
		wh <- which(chain == chainID[k])
		tfile3 <- tfile2[wh, ]

		#filter using model number. keep only the first model
		modelPosition <- which(fields == "_atom_site.pdbx_PDB_model_num")
		mdl <- unique(tfile3[[modelPosition]])
		wh <- which(tfile3[[modelPosition]] == mdl[1])
		tfile4 <- tfile3[wh, ]

		tfile5 <- list()

		seqPosition <- which(fields == "_atom_site.auth_seq_id")
		seq <- tfile4[[seqPosition]]

		#check if the whole chain is used
		if(as.numeric(gsub('[^0-9.-]','',as.character(endID[[k]][1]))) != 100000){

		# extract using start and end domain id
		for(j in 1:length(endID[[k]])){

			st <- startID[[k]][j]
			ed <- endID[[k]][j]

			whh3 <- which(seq == st)
			whh4 <- which(seq == ed)

			# if the start or end ID not present in the resolved part of the structure
			loop <- length(whh3)
			while(loop == 0){
				st <- as.numeric(gsub('[^0-9.-]','',st))
				whh3 <- which(seq == st)
				if(length(whh3) == 0){
					st <- as.numeric(gsub('[^0-9.-]','',st))+1
		    		whh3 <- which(seq == st)
				}
				loop <- length(whh3)
			}

			loop <- length(whh4)
			while(loop == 0){
				ed <- as.numeric(gsub('[^0-9.-]','',ed))
				whh4 <- which(seq == ed)
				if(length(whh4) == 0){
					ed <- as.numeric(gsub('[^0-9.-]','',ed))-1
			    	whh4 <- which(seq == ed)
				}
		    	loop <- length(whh4)		
			}

			wh1 <- min(whh3)
			wh2 <- max(whh4)
			tfile5[[j]] <- tfile4[wh1:wh2, ]

		}

		}else{

		tfile5[[1]] <- tfile4

		}

		# append the dataframes
		tfile6 <- do.call("rbind",tfile5)

		# extract coordinates for only the heavy atoms (S, O, C, N)
		atomPosition <- which(fields == "_atom_site.type_symbol")
		lineID2 <- tfile6[[atomPosition]]
		wh <- which(lineID2 == "C" | lineID2 == "S" | lineID2 == "O" | lineID2 == "N")
		tfile6 <- tfile6[wh, ]

		# keep only ATOM coordinates
		wh <- which(tfile6[[1]] == "ATOM")
		tfile6 <- tfile6[wh,]

		#call to create networks
		seqq <- tfile6[[seqPosition]]
		wh <- which(fields == "_atom_site.Cartn_x")
		xcoord <- tfile6[[wh]]
		wh <- which(fields == "_atom_site.Cartn_y")
		ycoord <- tfile6[[wh]]
		wh <- which(fields == "_atom_site.Cartn_z")
		zcoord <- tfile6[[wh]]

		outnetDirectory1 <- paste0(outputDirectory,'/dynamic-networks/',ann[[2]][k])
		dir.create(outnetDirectory1, recursive=TRUE, showWarnings=FALSE)
		pdb2mat2dynnet(seqq, xcoord, ycoord, zcoord, cutoff, naa, outnetDirectory1) 

	}

	cat("\nNetwork creation done.\n")


	if(choice == "Mode 2" | choice == "Mode 3"){

		## run dynamic graphlets
		outgraphletDirectory <- paste0(outputDirectory,'/dynamic-graphlets')
		dir.create(outgraphletDirectory)
		alldirs <- list.dirs(paste0(outputDirectory,"/dynamic-networks"), full.names=TRUE, recursive=FALSE)

		for(k in 1:length(alldirs)){

			tfiles <- gtools::mixedsort(list.files(alldirs[k], full.names=TRUE))
			tname <- paste(tfiles, collapse=",")

			command=paste("./bin/dcount",tname,6,4,1,paste0(outgraphletDirectory,"/",basename(alldirs[k])),"-v vectorfiles/graphlets_6_4.txt -g vectorfiles/orbits_6_4.txt -d '\t' -t")
			system(command)

		}

		## dynamic graphlet to vector form
		outgraphletvecDirectory <- paste0(outputDirectory,'/dyn-graphlets-vec')
		dir.create(outgraphletvecDirectory)
		dgdvm2vec(outgraphletvecDirectory, outgraphletDirectory, annotationFile)

		## eliminate non-changing columns
		elimcol(paste0(outgraphletvecDirectory,"/vec-",basename(annotationFile)), paste0(outgraphletvecDirectory,"/elim-",basename(annotationFile)))

		## PCA
		topca(paste0(outgraphletvecDirectory,"/elim-",basename(annotationFile)), 0.90, paste0(outgraphletvecDirectory,"/pca-",basename(annotationFile)))

		cat("\nFeature extraction done.\n")

	}


	if(choice == "Mode 3"){
		## convert to feature matrix
		outfeatureDirectory <- paste0(outputDirectory,'/feature-matrix')
		dir.create(outfeatureDirectory)
		tofeaturematrix(outgraphletvecDirectory, outfeatureDirectory, annotationFile)

		# check number of proteins
		temp <- fread(annotationFile, sep='\t', header=FALSE)
		numofcls <- unique(temp[[1]])
		cnt <- min(table(temp[[1]]))

		if(cnt >= 10 & length(numofcls) >= 2){

		if(partitionFlag == 0){
			## create partitions
			createPartition(annotationFile, paste0(outputDirectory,'/partitions'), 5)
		}else{
			cat("\nUsing the partitions provided by the user. Please make sure the format of the partitions are consistent with the program (see docmentation of the program for details).\n")
			system(paste('cp -r',partitionFolder, outputDirectory))
		}
		
		## run classification
		command=paste("python ./scripts/RunLR.py", outfeatureDirectory, 5, outputDirectory)
		system(command)

		cat("\nClassification done.\n")


		}else{

      cat("\nThe program could not perform classification.\n")
      cat("\nTo perform classification, the number of structural classes should be 2 or more and the number of proteins in each class should be 10 or more.\n")
    }

	}


}

