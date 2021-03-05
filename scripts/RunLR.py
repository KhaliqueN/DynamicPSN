import sys
from ClassifyingLR import classify

dataset = str(sys.argv[1])
nfolds = int(sys.argv[2])
outputDirectory = str(sys.argv[3])
classify(dataset, nfolds, outputDirectory)
