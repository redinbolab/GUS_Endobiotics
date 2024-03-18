########################################################################################
### Cladogram - Endobiotics Cohort Metagenomic GUSome
########################################################################################

### Usage
# Rscript ggtree_Endo_14_MGX_v3.r

########################################################################################
### main
########################################################################################

cat("loading packages...\n")

# Clean workspace 
rm(list=ls())

# Define required packages
required_packages <- c("ape", "Biostrings", "ggplot2", "ggtree", 
                       "tidyverse", "dplyr", "geiger")

# Identify packages that are not installed
installed <- installed.packages()[,"Package"]
not_installed <- setdiff(required_packages, installed)

# Install missing packages
if(length(not_installed) > 0) {
  install.packages(not_installed, dependencies = TRUE)
}

# Load packages invisibly
invisible(suppressPackageStartupMessages(lapply(required_packages, library, character.only = TRUE)))

# Mute output other than progress updates
# sink("/dev/null")


# Input data and output directories
cat("loading input data...\n")
data.file <- "./Endo_14_MGX_GUS_100_tax.afa.treefile"
tree <- read.tree(data.file)

outname <- 'ggtree_Endo_14_MGX.png'

loopDF_path <- "./endo_14_mgx_ggtree.csv"
loopDF <- read.csv(loopDF_path)

# Specify cladogram colors
values <- c(
  "Loop_1" = "#F08080",
  "Loop_2" = "#7B68EE",
  "Mini-Loop_1" = "#8FBC8F",
  "Mini-Loop_2" = "#FF7F50",
  "Mini-Loop_1_2" = "#DEB887",
  "No_Loop" = "#666666",
  "NTL" = "#00CCCC",
  "CTD" = "#FFC533")

cat('building plot...\n')

# Plot cladogram
treePlot <- ggtree(tree, aes(color=Loop_Class_Final,label=label), branch.length="none", layout="circular", size=1) %<+%
    loopDF +
    geom_tiplab(size=3, color="black",hjust=-.025,fontface='bold') + # label tips by name
    scale_fill_manual(na.translate=FALSE,values=values)+
    theme(
      plot.margin = margin(4,4,4,4, "cm"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA), legend.background=element_rect(fill = alpha("white", 0)), legend.key=element_rect(fill = alpha("white", 0)),
      legend.position="none",
      )+    

    xlim(-3.5, NA) + # modifies diam of center, changing shape of lines @ center. increase to increase diam of center
    
    scale_color_manual(
      values = c(
        "Loop_1" = "#F08080",
        "Loop_2" = "#7B68EE",
        "Mini-Loop_1" = "#8FBC8F",
        "Mini-Loop_2" = "#FF7F50",
        "Mini-Loop_1_2" = "#DEB887",
        "No_Loop" = "#666666",
        "NTL" = "#00CCCC",
        "CTD" = "#FFC533")
      )

# Build plot for further edits
pg <- ggplot_build(treePlot)

# Change default branch coloring
grey_to_smt <- function(x) { 
  if(x=='grey50') return('#B3B3B3')
  else return(x)
}

pg$data[[2]]$colour <- sapply(pg$data[[2]]$colour, grey_to_smt)
pg$data[[1]]$colour <- sapply(pg$data[[2]]$colour, grey_to_smt)

cat('coloring branches...\n')
for (x in pg$data[[2]]$node) {
  nodeLeafList1 <- tips(tree,x)
  str_x <- as.character(x)

  node_row <- pg$data[[2]] %>% filter(node == str_x)
  colorList <- c()
  if (length(nodeLeafList1)!=0){
    for (y in nodeLeafList1){
      str_y <- as.character(y)
      leaf_row <- pg$data[[2]] %>% filter(label == str_y)
      leafColor <- leaf_row$colour
      colorList <- c(colorList, leafColor)}
  uniqueColorList <- unique(colorList)
  colorCount <- length(uniqueColorList)
  if (colorCount == 1){
    color <- uniqueColorList[1]
    pg$data[[1]]$colour[x]<-color
    pg$data[[2]]$colour[x]<-color
}}}

q <- ggplot_gtable(pg)


# Plot and save
ggsave(outname, plot=q, device=png, height=10, width=14.5, bg = "transparent",dpi=1200)

########################################################################################
### Session Log
########################################################################################

writeSessionInfo <- function(outputFilePath) {
  # Retrieve R, package versions
  rVersion <- R.version.string
  packageVersions <- sapply(required_packages, function(pkg) paste0(pkg, " version ", packageVersion(pkg)))
  # Combine log information
  sessionInfo <- c(rVersion, packageVersions)
  # Write to log
  writeLines(sessionInfo, outputFilePath)
}

# Specify output path
outputFilePath <- "./session_log.txt"
# Write session log 
writeSessionInfo(outputFilePath)