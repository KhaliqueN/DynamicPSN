## IMPORTANT!
****If you downloaded the software before 23rd Feb. 2023, note that there was a bug in that code (briefly, the code was not taking the correct order of dynamic network snapshots for the dynamic graphlets computations). We have fixed the bug and uploaded the updated software. We suggest downloading the updated software and rerunning your analyses. Importantly, the bug did not affect the results reported in the paper, as the bug was introduced only after those results were produced, while aiming to create user-friendly version of the correct code that was used to produce the results in the paper.****

# Multilayer sequential PSNs

Please see documentation in docs/index.html


Folder structure
================

bin/ contains the dcount binary file that is used to compute dynamic graphlet degree vector matrix

docs/ contains documentation for our software (docs/index.html contains the main documentation file)

examples/ contains example annotation files

datasets/ contains all annotation files related to all datasets used in this study

cif-directory/ contains all of the Crystallographic Information Files (CIFs) that you use

scripts/ contains R or python scripts 

src/ contains the source file, dcount.cpp, which computes dynamic graphlet degree vector matrix 

vectorfiles/ contains helper text files for the dynamic graphlets program

output/ if present, contains the output files of the latest program run

output/dynamic-networks contains dyanmic the constructed dynamic PSNs

output/feature-matrix contains features of dynamic PSNs based on dynamic graphlets

output/Classification-results contains prediction files and performance accuracy for each of the 5 folds

====

This folder contains the original implementation of our software by Khalique Newaz, 2021. This folder is part of the publication "Multi-layer sequential protein structure networks improve protein structural classification", K. Newaz, J. Piland, P.L. Clark, S.J. Emrich, Jun Li, and T. Milenković (2021), Proteins: Structure, Function, and Bioinformatics.


