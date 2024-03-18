########################################################################################
### Metabolomics Heatmap Generation
########################################################################################

# Clean workspace 
rm(list=ls())

# Define required packages
required_packages <- c("pheatmap","dendextend", "RColorBrewer")

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


# Z Score calculation function
cal_z_score <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}


# Load input data 
fec_conj <- read.csv("glucuronide_aglycone_fecal_metabolites.csv")
fec_conj_num <- read.csv("name_drop_gluc_aglycone_fecal_metabolites.csv")

# Define Metabolite Class and Microbiota Vectors
Metabolite_Class <- rep(c("Inactive Glucuronide", "Active Aglycone"), times = c(12, 13))
Microbiota <- rep(c("Conventional", "GF"), each = 12)

# Associate row names
fecPhase <- fec_conj_num
rownames(fecPhase) <- fec_conj$Name

# Normalize the data by row
fecPhase_norm <- t(apply(fecPhase, 1, cal_z_score))

# Prepare annotation data frames with row names
fecPhase_Microbiota <- data.frame(Microbiota = factor(Microbiota))
row.names(fecPhase_Microbiota) <- colnames(fecPhase)
fecPhase_Metabolite_Class <- data.frame(Metabolite_Class = factor(Metabolite_Class))
row.names(fecPhase_Metabolite_Class) <- rownames(fecPhase)

# Define group annotation colors
ann_colors <- list(Microbiota = c(Conventional = "#5E7CE2",GF = "#FCD581"), 
                   Metabolite_Class = c("Inactive Glucuronide" = "#CCD1D1", "Active Aglycone" = "#7E5A6D"))

########################################################################################
### Heatmap Generation
########################################################################################
fecPhase_pheatmap <- pheatmap(fecPhase_norm,
                              color = colorRampPalette(c("#080357", "#EAEBE7", "#50C878"))(500),
                              annotation_col = fecPhase_Microbiotarow,
                              annotation_row = fecPhase_Metabolite_Classcol,
                              cutree_cols = 2,
                              cutree_rows = 2,
                              border_color = "black",
                              fontsize_row = 8,
                              fontsize_col = 7,
                              annotation_colors = ann_colors
                              )

# Function to save heatmap 
save_pheatmap_png <- function(x, filename, width=5000, height=3600, res = 600) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Save heatmap
save_pheatmap_png(fecPhase_pheatmap, "gf_v_conv_fecal_fc_heatmap.png")

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

