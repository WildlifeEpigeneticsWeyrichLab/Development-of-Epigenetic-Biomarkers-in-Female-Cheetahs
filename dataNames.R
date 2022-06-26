## Created 6/29/2015 by Daniel Beck
## Last modified 5/31/2016
#Modified by Johana P. 20/03/2019
#Modified by Tania G. 25/03/2019

## This is the configuration script for the MeDIP-seq analysis pipeline. It holds
## most analysis parameters and options. It also holds sample/filename information
## and defines which analyses are performed. It also performs some preliminary 
## checks to ensure the analysis can procede. 

################################
## General project attributes ##
################################
## This section sets the folder structure. The default is to have a single project
## folder with three main subfolders for code,f data, and results. An additional 
## genome directory holds the reference genome and is located within the data folder.
## It should be fine to modify this structure, however, extensive testing of any
## alternatives has not been done.

projectName <- "cheetah_methylation"

# Project directory (main folder)
#absolute path to directory
projectDirectory <- "/data/fg2/yasar/projects/Cheetah_T1/"
# Data directory (holds data files)
dataDirectory <- paste(projectDirectory, "data/new_mapped/SE_read1/", sep="")
# Code directory (holds all code used for the MeDIP-seq analysis)
codeDirectory <- paste(projectDirectory, "code/", sep="")
# Results directory (folder for all results)
resultsDirectory <- paste(projectDirectory, "results/", sep="") 
# Genome directory (holds reference genome and BSgenome package)
genomeDirectory <- paste(projectDirectory, "data/genome/", sep="")


################################
## BSgenome package variables ##
################################
## These options are used for the BSgenome package creation. The seqnames parameter
## is often the most tedious (especially for large numbers of scaffolds). This will
## need to be improved in the future. 
#Tania G. It can be easier to create the BSgenome using the forgeBSGenome2_Cheetah.R script
# in which case the output will be and work as a more 
#"strict" BSgenome package (look for notes along scripts for differences on its usage)
## BSgenome package. 08/07/2019

# The seed file name for BSgenome package creation
seedName <- "BSgenome.AciJubGCF.NCBI.Acijub2-seed"
# The package name (NO UNDERSCORES OR DASHES!!!)
bsgenomePackageName <- "BSgenome.AciJub.NCBI.Acijub2"
# Name for the reference genome
referenceName <- "AciJub"
# Species name
species <- "Acinonyx jubatus"
# Source of reference genome
provider <- "NCBI_RefSeq"
# Reference genome version
version <- "2_Rico"
# Release date
release_date <- "2018/11/13"
# Names of all scaffolds or chromosomes
seqnames <-  "paste('NW_0',seq(from=020834726,to=020837944,by=1),'.1', sep='')"
# Any prefix needed to regenerate filename from seqnames
seqfiles_prefix <- "AciJub_2_"
# Any suffix needed to regenerate filename from seqnames
seqfiles_suffix <- ".fasta"
# Source directory of scaffolds
seqs_srcdir <- genomeDirectory

###############################
## Raw file / sample pairing ##
###############################
## This section associates files with samples. The files are assumed to be located in the
## dataDirectory folder. The full path will be based on that directory.

# The seqFiles data frame holds all sample information. The ctFlag is an identifier for the 
# treatment group, but any attribute of the sample can be used. 
seqFiles <- data.frame(
  sampleName = c("FRA135", "FRA138", "FRA171", "FRP096", "FRP114", "A2462", "A2627", "A3657", "Z192", "Z194", "Z196", "Z197", "Z3037", "Z3044", "Z3052"), 
 # p1FileName = c("FRA135_R1.fastq.gz", "FRA138_R1.fastq.gz", "FRA171_R1.fastq.gz", "FRP096_R1.fastq.gz", "FRP114_R1.fastq.gz", "Z192_R1.fastq.gz", "Z194_R1.fastq.gz", "Z195_R1.fastq.gz", "Z196_R1.fastq.gz", "Z197_R1.fastq.gz"),
  #p2FileName = c("FRA135_R2.fastq.gz", "FRA138_R2.fastq.gz", "FRA171_R2.fastq.gz", "FRP096_R2.fastq.gz", "FRP114_R2.fastq.gz", "Z192_R2.fastq.gz", "Z194_R2.fastq.gz", "Z195_R2.fastq.gz", "Z196_R2.fastq.gz", "Z197_R2.fastq.gz"),
  ctFlag = c("FR","FR","FR","FR","FR","FR","FR","FR", "Zoo", "Zoo", "Zoo", "Zoo","Zoo","Zoo","Zoo"),
  stringsAsFactors = F
)


bamFileName <- paste(seqFiles$sampleName, ".bam", sep="")


#######################
## MEDIPS parameters ##
#######################
## These parameters modify the MEDIPS analysis. The defaults shown here have not been fully
## explored. It isn't clear what the ideal values might be. Most of these are inherited from
## Haque's code. See MEDIPS documentation for details.


# uniq = max allowed number of stacked reads per genomic position
uniq <- 1
# extend = all reads will be extended to a length of xx nt according to strand information
extend <- 350
# shift is alterntaive to extend: reads are not extended but shifted by the specified number of nucleotides
shift <- 0
# window size
ws = 200
# subselection of chromosomes -> usefull for quality control analysis
chr.select <- NULL 

#NW_020834836.1

# correct p.values for multiple testing. Available are: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none
p.adj <- "fdr"
# methods for calculating differential coverage (options: ttest and edgeR. ttest will be calculated only in case there are at least three replicates per group. edgeR can be applied when there are less than three replicates per group
diff.method <- "edgeR"


# This parameter determines, if a CpG density dependent relative methylation scores (rms) will be calculated for the MEDIPS SETs given at the slots MSet1 and MSet2.
MeDIP <- FALSE

#The coupling set must be created based on the same reference genome, the same set of chromosomes, and with the same window size used for the MEDIPS SETs. 
#For this, we specify the first MEDIPS SET in the hESCs object as reference, but any of the other MEDIPS SETs would be fine as well, 
#because all of them consist of the same set of chromosomes and have been generated with the same window size.
CScalc <- FALSE 

#  copy number variation will be tested by applying the package DNAcopy to the window-wise log-ratios calculated based on the the means per group. 
CNV <- FALSE

# threshold for a minimum sum of counts across all samples per window (default=10). Windows with lower coverage will not be tested for differential coverage.
minRowSum <- 10

# This vector holds all raw p-value thresholds to use for the analyses. 
pValues <- c(1e-03, 1e-04, 1e-05, 1e-06, 1e-07)
# This vector holds all multiple testing adjusted p-value thresholds to use for the analyses. 
MTCpValues <- c(0.3, 0.2, 0.1, 0.05)

####################
## DMR parameters ##  = Differentially methylated regions
####################
## These parameters define how DMR edges are defined. These values are completely arbitrarily
## chosen. It isn't clear how these values should be chosen to accuratly reflect a biologically
## significant feature.

# This p-value threshold defines DMR boundaries TANIA: AGAIN these settings need to be thorughly studied and chosen! 
dmrBoundPvalue <- 0.1
# Adjacency distance (=1 when windows must be exactly adjacent). This determines how far apart
# significant windows can be and remain in the same DMR.
adjDist <- 1000

# The maxDMRnum variable gives a maximum number of DMRs on which to calculate CpG density and other
# information. This speeds up the pipeline. However, this will need to be increased if the p-value 
# of interest has more DMR than this number.
maxDMRnum <- 5000

###########################
## Annotation parameters ##
###########################
## These variables hold information on the annotation. 

annotationType <- "gff"

annotationGFF <- "Acijub_2_GCF_003709585.1_genomic.gff"

# Description of annotation source
ann_source <- "asemmbly from Rico_2"

# Organism name
org_name <- "Aci jub"

# Taxonomy ID
taxID <- 32536

############################
## Comparisons / analyses ##
############################
## This section defines which comparisons should be made. All analyses will iterate through
## each of these.

# Names for the comparisons ### TANIA: Not sure how this affected the analysis, this was the "default" but I got the feeling that it should be smth else
comparisonNames <- c("all_R1_onlyF")

# Which samples are being compared. The pairs flag was included to allow for a pairwise
# type analysis. It isn't currently functional, but some code requires it.
#comparison <- list()
#comparison[[1]] <- data.frame(mset1=c(1:8), mset2=c(9:15,NA))
# Description added directly to report
comparisonDescription <- list()
comparisonDescription[[1]] <- "All free-ranging vs zoo cheetahs"

#################
## Data checks ##
#################
## These are very basic checks to detect problems in this configuration file. This section
## should be extended each time an unexpected error occurs in this script.

if (length(comparison)!=length(unique(comparisonNames))) { 
  warning("The number of comparison names does not match the number of comparisons.\n",
          "This could be caused by duplicate comparison names.")
}
if (CScalc & !MeDIP) {
  warning("The CSset coupling set is being calculated (CScalc = TRUE) 
          but not used (MeDIP = FALSE).")
}
if (!CScalc & MeDIP) {
  warning("The CSset coupling set is not calculated (CScalc = FALSE) 
          but is needed (MeDIP = TRUE).")
}

