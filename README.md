# DynamicPSN

Please see documentation in docs/index.html


Folder structure
================

bin/ contains the dcount binary file that is used to compute dynamic graphlet degree vector matrix

docs/ contains documentation for DynaP (docs/index.html contains the main documentation file)

examples/ contains example annotation files

cif-directory/ contains all of the Crystallographic Information Files (CIFs) that you use

scripts/ contains R or python scripts 

src/ contains the source file, dcount.cpp, which computes dynamic graphlet degree vector matrix 

vectorfiles/ contains helper text files for the dynamic graphlets program

output/ if present, contains the output files of the latest program run

output/dynamic-networks contains dyanmic the constructed dynamic PSNs

output/feature-matrix contains features of dynamic PSNs based on dynamic graphlets

output/Classification-results contains prediction files and performance accuracy for each of the 5 folds

