########################################################################################
### MixOmics sPLS-DA - Endobiotics Cohort Meta-Omic Profiles Cluster by Inhibition
########################################################################################

### Usage
# metagenomics
# Rscript mixomics_sPLS-DA_iter_export_v4.r WGS_genus_mixomics.csv ./mixomics_MGX/

# metaproteomics
# Rscript mixomics_sPLS-DA_iter_export_v4.r metaproteomics_mixomics.csv ./mixomics_MPX/

########################################################################################
### main
########################################################################################

# Clean workspace 
rm(list=ls())

# Define required packages
required_packages <- c("mixOmics", "dplyr")

# Identify packages that are not installed
installed <- installed.packages()[,"Package"]
not_installed <- setdiff(required_packages, installed)

# Install missing packages
if(length(not_installed) > 0) {
  install.packages(not_installed, dependencies = TRUE)
}

# Load packages invisibly
invisible(suppressPackageStartupMessages(lapply(required_packages, library, character.only = TRUE)))

# Return an error if required packages are not loaded, invisible above hides any such errors
loaded_packages <- character(length(required_packages))
for (i in seq_along(required_packages)) {
  loaded_packages[i] <- tryCatch({
    library(required_packages[i], character.only = TRUE)
    required_packages[i]
  }, error = function(e) {
    message(paste("Error loading package:", required_packages[i], "-", conditionMessage(e)))
    NA
  })
}


# Load input data and output directory

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Insufficient arguments. Please provide the input file and output path")
}

in_csv <- as.character(args[1])

outdir <- as.character(args[2])


if (!dir.exists(outdir)){
  dir.create(outdir)
}

inhibition_csv <- "./inhibition_clustered.csv"

feature_df <- read.csv(in_csv)
rownames(feature_df) <- feature_df$sample_ID
feature_df$sample_ID <- NULL

inhibition_df <- read.csv(inhibition_csv)
rownames(inhibition_df) <- inhibition_df$sample_ID

# main mixOmics call
for(inhibition_col in colnames(inhibition_df)) {
    if (inhibition_col != 'sample_ID'){
        X <- feature_df # use sample feature data as X matrix
        Y <- inhibition_df[[inhibition_col]] # use inhibition capacities as Y matrix
        pdf(paste0(outdir,inhibition_col,"_mixomics_out.pdf"))
        result.splsda <- splsda(X, Y) # run the sPLS-DA method
        groupColors <- c('#666666','#8E6580') # specify group coloring
        plotIndiv(result.splsda, # plot
                comp=1:2,
                legend=TRUE,
                title=inhibition_col,
                legend.title = inhibition_col,
                style = 'ggplot2',
                col.per.group = groupColors, 
                ellipse=TRUE) # include 95% confidence ellipse for each class
            
    }
}

########################################################################################
### Session Log
########################################################################################
writeSessionInfo <- function(outputFilePath) {

  # Retrieve R, package versions
  rVersion <- R.version.string
  packageVersions <- sapply(required_packages, function(pkg) paste0(pkg, " version ", packageVersion(pkg)))

  # Combine all information
  sessionInfo <- c(rVersion, packageVersions)

  # Write to log
  writeLines(sessionInfo, outputFilePath)
}

# Specify output path
outputFilePath <- "./session_log.txt"

# Write session log 
writeSessionInfo(outputFilePath)