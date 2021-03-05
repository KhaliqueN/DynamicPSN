args=commandArgs(TRUE)

if (length(args)< 5) {
  stop("#### Code should be run as 'Rscript maincmd.r [Distance cutoff] [Number of amino acids] 
    [Choice of task] [Path to the annotation file] [Path to the directory with all .cif files]'\n #### Check the manual for more details ", call.=FALSE)
}

# if (!require("shiny")) {
#   install.packages("shiny")
#   library(shiny)
# }

# if (!require("shinyjs")) {
#   install.packages("shinyjs")
#   library(shinyjs)
# }

if (!require("tools")) {
  install.packages("tools")
  library(tools)
}

if (!require("data.table")) {
  install.packages("data.table")
  library(data.table)
}

if (!require("stringr")) {
  install.packages("stringr")
  library(stringr)
}

if (!require("igraph")) {
  install.packages("igraph")
  library(igraph)
}

source("./scripts/utils.r")
source("./scripts/maincmd.r")

cutoff <- as.numeric(args[1])
naa <- as.numeric(args[2])
choice <- args[3]
annotationFile <- args[4]
pdbDirectory <- args[5]


runmaincmd(cutoff, naa, choice, annotationFile, pdbDirectory) 

