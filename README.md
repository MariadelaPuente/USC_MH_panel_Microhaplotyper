# USC_MHs_Microhaplotyper

Pipeline based on publicly avaliable software for microhaplotype calling for MPS data based on the VISAGE Enhanced Tool




Before you start, install the following programs:

1.- Burrows-Wheeler Aligner (BWA), avaliable at https://sourceforge.net/projects/bio-bwa/files/ 

2.- Samtools, avaliable at https://sourceforge.net/projects/samtools/files/samtools/

3.- R (and Rstudio), avaliable at https://www.r-project.org/ (and https://www.rstudio.com/products/RStudio/)

4.- Open R (in Rstudio) install the following packages from CRAN repository by running:

install.packages("tidyverse")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("data.table")

5.- Install in R the microhaplot package avaliable at https://github.com/ngthomas/microhaplot



Instructions for running MICROHAPLOTYPER

1.- Place a folder (preferably on your Desktop) containing all the .fastq from the samples to be analised. 
	Avoid negative controls or samples containing almost no reads as they might hinder the pipeline.
	Name the folder according to a RUN_NAME (preferably with no spaces).

2.- Copy and paste all the files from the MICROHAPLOTYPER folder into the folder with the .fastq files

3.- Open a terminal and change directory to the fastq folder by typing: 
	cd Desktop/RUN_NAME

4.- Check you entered the folder and run the script for obtaining the input files for R by typing:
	sh S5_microhaplotyper_USC_MHs_pipeline.sh
	
5.- Once it's ready, open the S5_microhaplotyper_USC_MHs_pipeline.R	and set the working directory to the fastq folder.
	Customize RUN_NAME (and folders path if your folder is not in Desktop).
	Customize filtering thresholds.
	Select and run the whole script.
	
6.- Folder should now contain, for each sample:
	A .sam file that can be used to inspect alignments.
	A .pdf MH profile.
	A genotypes.txt table (* in the genotype column means that more than 2 alleles were found).
	
