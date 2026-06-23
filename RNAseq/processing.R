# Processing of nfcore output for QC before DESEQ

# Load libraries:
library(limma)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(RColorBrewer)
library(viridis)
library(PCAtools)
library(patchwork)
library(ggpubr)
library(rstatix)


# Set working directory:
setwd("~/Documents/My_R")

setwd("/Users/velandn/Documents/My_R")


# Create output_dir:
output_dir <- "/Users/velandn/Documents/My_R/RNAplots_Dec2024/"

ifelse(!dir.exists(file.path(output_dir)),
       dir.create(file.path(output_dir)),
       "Directory Exists")

