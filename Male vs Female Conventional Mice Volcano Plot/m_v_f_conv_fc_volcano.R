########################################################################################
### Metabolomics Volcano Plot Generation
########################################################################################

# Clean workspace 
rm(list=ls())

# Define required packages
required_packages <- c("ggplot2", "RColorBrewer")

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
data <- read.csv("mvf_fc.csv")

# Calculate the logarithm of Fold and -log10(P)
data$log2Fold <- log2(data$Fold)
data$minuslog10P <- -log10(data$P)

# Define the hexadecimal colors for each class
class_colors <- c("#7E5A6D", "#CCD1D1", "#0000FF")  

# Manually specify major x-axis tick positions in the volcano plot
custom_breaks <- c(-10,-5, 0, 5, 10)  

# Specify the desired width
width <- 35



# Create the plot
plot <- ggplot(data, aes(x = log2Fold, y = minuslog10P, color = Class)) +
  geom_point(size = 5) +
  geom_point(shape = 1, size = 5, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1.25) +
  labs(x = expression(log[2]("Fold Change")), y = expression(-log[10]("P"))) +
  theme_minimal() +
  scale_color_manual(values = class_colors) +
  scale_x_continuous(breaks = custom_breaks, limits = c(-max(abs(data$log2Fold)), max(abs(data$log2Fold)))) +  # Symmetrical limits
  scale_y_continuous(limits = c(0, NA)) +
  geom_hline(yintercept = 1.301029996, linetype = "dashed", color = "red", linewidth = 1.25) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 3),
        panel.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 30),
        axis.title = element_text(size = 30, face = "bold"),
        legend.position = "none")  # Remove legend
plot

# Save the plot
ggsave("gf_v_conv_fc_volcano.png", plot, width = 10, height = 10, units = "in")

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


