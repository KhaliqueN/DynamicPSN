args=commandArgs(TRUE)

if (length(args) < 6) {
  stop("#### Code should be run as 'Rscript maincmd.r [Distance cutoff] [Number of amino acids] 
    [Choice of task] [Path to the annotation file] [Path to the directory with all .cif files] [Partition flag] [Partition folder]'\n #### Check the program manual for more details ", call.=FALSE)
}


suppressMessages(library(tools))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(igraph))


source("./scripts/utils.r")
source("./scripts/maincmd.r")

cutoff <- as.numeric(args[1])
naa <- as.numeric(args[2])
choice <- args[3]
annotationFile <- args[4]
pdbDirectory <- args[5]
partitionFlag <- args[6]
partitionFolder <- args[7]

runmaincmd(cutoff, naa, choice, annotationFile, pdbDirectory, partitionFlag, partitionFolder) 

