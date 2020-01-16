DYSE algorithm

SOFTWARE AUTHOR: Basma Abdelkarim
CONTACT: basmataher94@gmail.com

CITATIONS:
For ROSE:
Warren A. Whyte, David A. Orlando, Denes Hnisz, Brian J. Abraham, Charles Y. Lin, Michael H. Kagey, Peter B. Rahl, Tong Ihn Lee and Richard A. Young Cell 153, 307-319, April 11, 2013

Jakob Lovén, Heather A. Hoke, Charles Y. Lin, Ashley Lau, David A. Orlando, Christopher R. Vakoc, James E. Bradner, Tong Ihn Lee, and Richard A. Young
Cell 153, 320-334, April 11, 2013

For DYSE:
Abdelkarim, B. 2020, 'Dynamics of Super-Enhancers Throughout Myogenesis', MSc thesis, University of Ottawa, Ottawa
Developed using Python 2.7.13, SAMtools 0.1.18, BEDtools 2.26.0, R 3.5.0

Purpose: To compare super-enhancers between different stages in a given cell lineage

Extra Python packages used:
-matplotlib
-matplotlib_venn
-pandas
-joblib

1) REQUIREMENTS

DYSE must be run from the directory in which it is stored
ROSE algorithm package must be installed in the same directory with DYSE
(https://bitbucket.org/young_computation/rose/src)
.bam files of sequence reads for ranking factor of interest and control factor (if used) -> .bam files must be sorted and indexed
.gff files of constituent enhancers identified, for each factor
.gff files must have the same format described for the ROSE algorithm

2)RUN_FILES

DYSE_main.py: main program - distinguishes stage-specific enhancers from those that occur in multiple stages
ROSE_diffexp.py: uses differential expression data to identify SE-associated genes that are differentially expressed between stages

3) USAGE

DYSE is run using the following command:
python DYSE_main.py -i INPUT_FILES -g GENOME -t ANALYSIS_TYPE -o OUTPUT_FOLDER [OPTIONAL: -r ROSE_OUT_FOLDER -s STITCHING_DISTANCE]

NOTE: If the parameter -r ROSE_OUT_FOLDER is used, DYSE will call the ROSE algorithm first to identify SEs for each factor used then group together SEs from the same stage
otherwise, DYSE will consifder the input files to be a list of SE enhancer files and proceed without calling the ROSE algorithm, in that case, a .bam sequence reads file for the ranking factor is NOT necessary.

Required parameters:
-i  INPUT_FILES: .txt files (one for each stage) separated by commas. Each file consists of a list in the following format:
constituent_enhancers_file.gff,ranking_file.bam,control.bam (optional)
-g  GENOME: available USCS genome builds used for read mapping are hg18, hg19, mm8, mm9 or mm10
-t  ANALYSIS_TYPE: choices are either 1 (to compare SEs between adjacent stages only) or 2 (to compare SEs between any pairs of stages)
-o  OUTPUT_FOLDER: FOLDER to store output files

Optional parameters:
-r  ROSE_OUTPUT_FOLDER: Name of folder to store output from ROSE algorithm. Using this argument prompts DYSE to call the ROSE algorithm to identify SEs for each factor of interest.
-s  STITCHING_DISTANCE: maximum distance to be used to stitch two adjacent enhancers, the dfault is 12.5 kb

NOTE: TSS exclusion zone is set to 2500 bp as recommended by the ROSE algorithm

4) CODE RUN


