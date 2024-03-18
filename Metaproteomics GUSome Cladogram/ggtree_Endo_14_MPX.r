########################################################################################
### Cladogram - Endobiotics Cohort Proteomic GUSome
########################################################################################

### Usage
# Rscript ggtree_Endo_14_MPX.r

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

########################################################################################
# Prepare data
########################################################################################

cat("loading input data...\n")

# REPLACE FILENAME
data.file <- "./Endo_14_MGX_100_GUS_MPX_tax.afa.treefile"
tree <- read.tree(data.file)

outname <- 'ggtree_Endo_14_MPX_labels.png'

loopDF_path <- "./endo_14_mpx_ggtree_correlationLabels.csv"
loopDF <- read.csv(loopDF_path)

loopDF$Correlation[loopDF$Correlation==""] <- NA
loopDF$Correlation[loopDF$Correlation==" "] <- NA
# Now exclude NA values from factor levels
loopDF$Correlation <- factor(loopDF$Correlation, exclude = NA)

shapeValuesMapped <- c(
  "CTD1" = 22,
  "CTD2" = 24
  )

colorMapping <- c(
        "CTD1" = "#FFC533",
        "CTD2" = "#FFC533"
        )

values <- c(
  "Loop_1" = "#F08080",
  "Loop_2" = "#7B68EE",
  "Mini-Loop_1" = "#8FBC8F",
  "Mini-Loop_2" = "#FF7F50",
  "Mini-Loop_1_2" = "#DEB887",
  "No_Loop" = "#666666",
  "NTL" = "#00CCCC",
  "CTD" = "#FFC533")

loops_uniq <- unique(loopDF$Loop_Class_Final)
N <- length(loops_uniq)

########################################################################################
# build plot
########################################################################################

treePlot <- ggtree(tree, aes(color=Loop_Class_Final,fill=Correlation,label=label), branch.length="none", layout="circular", size=1) %<+%
    loopDF +
    geom_tiplab(size=3, color="black",hjust=-.045,fontface='bold') + # label tips by name
    geom_tippoint(aes(shape=Correlation), size=3.8 ,na.rm=TRUE,color='black')+
    scale_shape_manual(
    values=shapeValuesMapped)+
    scale_fill_manual(
      na.translate=FALSE,
      values=colorMapping)+

    theme(
      plot.margin = margin(4,4,4,4, "cm"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      legend.position="none",)+    
      
    xlim(-3.5, NA) +

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

# rotate for loop class
treePlot <- ggtree::rotate(treePlot,66)
treePlot <- ggtree::rotate(treePlot,55)
treePlot <- ggtree::rotate(treePlot,65)

# # rotate for corr groups
treePlot <- ggtree::rotate(treePlot,74)
treePlot <- ggtree::rotate(treePlot,71)
treePlot <- ggtree::rotate(treePlot,69)
treePlot <- ggtree::rotate(treePlot,57)
treePlot <- ggtree::rotate(treePlot,96)
treePlot <- ggtree::rotate(treePlot,95)
treePlot <- ggtree::rotate(treePlot,67)
treePlot <- ggtree::rotate(treePlot,77)

treePlot <- ggtree::rotate(treePlot,71)
treePlot <- ggtree::rotate(treePlot,73)
treePlot <- ggtree::rotate(treePlot,60)
treePlot <- ggtree::rotate(treePlot,63)

treePlot <- ggtree::rotate(treePlot,68)

# big shift
treePlot <- ggtree::rotate(treePlot,53)

treePlot <- ggtree::rotate(treePlot,65)
treePlot <- ggtree::rotate(treePlot,66)

treePlot <- ggtree::rotate(treePlot,64)
treePlot <- ggtree::rotate(treePlot,92)

treePlot <- ggtree::rotate(treePlot,65)
treePlot <- ggtree::rotate(treePlot,93)

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

# plot & save
# plot(q)
ggsave(outname, plot=q, device=png, height=10, width=12.5, bg = "transparent",dpi=1200)


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