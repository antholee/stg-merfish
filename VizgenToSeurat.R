# load libs ----
#remotes::install_github(repo = 'alikhuseynov/seurat', ref = 'feat/vizgen')
library(ggplot2)
library(Seurat)
library(dplyr)
library(magrittr)
library(BiocParallel)
library(progressr)
library(sfarrow)
library(spatstat)    
library(sf)
library(future)
plan("multisession", workers = 8)
library(data.table)
library(arrow)
library(scales) 
library(parallel)
library(tidyverse)
library(irlba)
library(Matrix)
library("leidenAlg")




# helper function to return Seurat object metadata
callmeta <- function (object = NULL) { 
  return(object@meta.data) 
}


extractSampleID <- function(pathname){
  # Check if pathname contains "aMTG"
  if (grepl("aMTG", pathname)) {
    # Extract the substring between "aMTG" and the next hyphen as the sample name
    sample_name <- gsub(".*aMTG(\\d+)-.*", "\\1", pathname)
    sampleID <- paste("aMTG", sample_name, sep = "")
  }
  # Check if pathname contains "pSTG"
  if (grepl("pSTG", pathname)) {
    # Extract the substring between "pSTG" and the next hyphen as the sample name
    sample_name <- gsub(".*pSTG(\\d+)-.*", "\\1", pathname)
    sampleID <- paste("pSTG", sample_name, sep = "")
  }
  # Check if pathname contains "aSTG"
  if (grepl("aSTG", pathname)) {
    # Extract the substring between "aSTG" and the next hyphen as the sample name
    sample_name <- gsub(".*aSTG(\\d+)-.*", "\\1", pathname)
    sampleID <- paste("aSTG", sample_name, sep = "")
  }
  # Check if pathname contains "mac"
  if (grepl("mac", pathname)) {
    # Extract the substring between "mac" and the next hyphen as the sample name
    sample_name <- gsub(".*mac(\\d+)-.*", "\\1", pathname)
    sampleID <- paste("mac", sample_name, sep = "")
  }
  
  return(sampleID)
}

## comment out the old samples and run just the newest one, then after SeuratObject created, uncomment all pathnames to include all samples/ filepaths to data:
pathnames <- c(
  # 'XXX'
)

##create SeuratObject for each merfish experiment ----
for (pathname in pathnames) {
  sample <- extractSampleID(pathname)
  
  start.time <- Sys.time()
  
  vizgen_seurat <- LoadVizgen(data.dir = pathname,  
               fov = sample, 
               assay = "Vizgen",
               metadata = c("volume", "fov"), # add cell volume info
               type = c("segmentations", "centroids"), # type of cell spatial coord matrices
               z = 3L,
               add.zIndex = TRUE, # add z slice section to a cell
               update.object = TRUE,
               use.BiocParallel = TRUE,
               workers.MulticoreParam = 10, # for `BiocParallel` processing
               verbose = TRUE
                )
  
  end.time <- Sys.time()
  message("Time taken to load object = ", 
          round(end.time - start.time, digits = 2), " ", 
          attr(c(end.time - start.time), "units"))
  vizgen_seurat
  vizgen_seurat %>% callmeta %>% str
  file_path <- file.path("~/Dropbox/Lab/MERFISH", paste0(sample, "raw.obj.rds"))
  saveRDS(vizgen_seurat, file=file_path)
  rm(vizgen_seurat)
}


## create list of SeuratObjects, do this after uncommenting pathnames above so that it includes all samples----

allrawMERFISH <- list()

# Loop over each pathname
for (pathname in pathnames) {
  # Dynamic variable name
  curr_dataset <-extractSampleID(pathname)
  file_path <- file.path("XXX"))
  
  # Load the object from the RDS file
  loaded_object <- readRDS(file=file_path)

  allrawMERFISH[curr_dataset] <- loaded_object

  rm(loaded_object)
}


##plot pre-filtering baseline data ----
# Extract nFeatures from each Seurat object in the list
nFeatures_list <- lapply(allrawMERFISH, function(seurat_object) {
  return(seurat_object$nFeature_Vizgen)
})

volumes_list <- lapply(allrawMERFISH, function(seurat_object) {
  return(seurat_object$volume)
})

nCounts_list <- lapply(allrawMERFISH, function(seurat_object) {
  return(seurat_object$nCount_Vizgen)
})

prefilter.prop_df <- do.call(rbind, lapply(seq_along(nFeatures_list), function(i) {
  data.frame(Sample = names(allrawMERFISH)[i], nFeatures = nFeatures_list[[i]], nCounts = nCounts_list[[i]], volumes = volumes_list[[i]])
}))

# Create violin plots
ggplot(prefilter.prop_df, aes(x = Sample, y = nFeatures)) +
  geom_violin() +
  labs(x = "Sample", y = "nFeatures") +
  scale_y_continuous(breaks = seq(0, max(prefilter.prop_df$nFeatures), by = 25)) +
  theme(axis.text.x = element_text(angle = 0)) +
  scale_y_continuous(trans = "log10", 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = function(x) format(x, scientific = FALSE),
                     limits = c(1, 10^3)) +
  geom_hline(yintercept = c(10), linetype = "dashed", color = "red")

ggplot(prefilter.prop_df, aes(x = Sample, y = nCounts)) +
  geom_violin() +
  labs(x = "Sample", y = "nCounts") +
  scale_y_continuous(breaks = seq(0, max(prefilter.prop_df$nCounts), by = 25)) +
  theme(axis.text.x = element_text(angle = 0)) +
  scale_y_continuous(trans = "log10", 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = function(x) format(x, scientific = FALSE),
                     limits = c(1, 10^3.5)) +
  geom_hline(yintercept = c(10), linetype = "dashed", color = "red")

ggplot(prefilter.prop_df, aes(x = Sample, y = volumes)) +
  geom_violin() +
  labs(x = "Sample", y = "volumes") +
  scale_y_continuous(breaks = seq(0, max(prefilter.prop_df$volumes), by = 25)) +
  theme(axis.text.x = element_text(angle = 0))  +
  scale_y_continuous(trans = "log10", 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = function(x) format(x, scientific = FALSE),
                     limits = c(10, 10^4)) +
  geom_hline(yintercept = c(300, 5000), linetype = "dashed", color = "red")



##QC ----

# Define a function to calculate the 98th quantile of nCount for a single Seurat object
calculate_quantile <- function(seurat_object) {
  return(quantile(seurat_object@meta.data$nCount, probs = 0.98))
}

# Use lapply to apply the function to each Seurat object in the list
quantiles <- lapply(allrawMERFISH, calculate_quantile)

# Print the quantiles for each Seurat object
names(quantiles) <- names(allrawMERFISH)  # Assign names to quantiles corresponding to Seurat objects
quantiles <- unlist(quantiles)
mean.quantile <- mean(quantiles)

# Define a function to subset Seurat objects based on nCount, volume, and nFeatures
subset_data <- function(seurat_object) {
   subset_condition <- seurat_object@meta.data$nCount_Vizgen >= 10 & 
    seurat_object@meta.data$nCount_Vizgen <= mean.quantile &
    seurat_object@meta.data$volume >= 300 &
    seurat_object@meta.data$volume <= 5000 &
    seurat_object@meta.data$nFeature_Vizgen >= 10 
  # Subset Seurat object based on the condition
  return(seurat_object[, subset_condition])
}








