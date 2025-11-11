# BiocManager::install("phyloseq")
library(curatedMetagenomicData)
library(ExperimentHub)
library(dplyr)
library(phyloseq)
library(vegan)


# Download the metadata for all studies
meta_data <- sampleMetadata


target_metadata <- meta_data %>%
  filter(
    menopausal_status %in% c("pre", "post"),
    gender == "female",
    body_site == "stool",
    !is.na(age) # Make sure we have age for our covariate
  )

# See how many samples are in each group:
#pre, post
cat("Found", nrow(target_metadata), "total samples.\n")
table(target_metadata$menopausal_status)

write.csv(target_metadata, file = "menopausal_dataset.csv", row.names = FALSE)
write.csv(meta_data, file = "all_data.csv", row.names = FALSE)

# Download the actual microbiome data for ONLY those samples
# This returns a list of SummarizedExperiment objects, one for each study
# We specify dataType = "relative_abundance"
merged_eset_list <- returnSamples(
  sampleMetadata = target_metadata,
  dataType = "relative_abundance",
  counts = FALSE
)

#Create bacteria matrix (Rows=Bacteria, Columns=Samples)
bacteria_matrix <- assay(merged_eset_list)
write.csv(bacteria_matrix, file = "bacteria_matrix.csv", row.names = TRUE)

final_metadata <- colData(merged_eset_list)
write.csv(final_metadata, file = "metadata.csv", row.names = TRUE)

#Create our metadata with menopausal status (Rows=Samples, Column=Menopausal_status)
filtered_metadata <- final_metadata[, c("menopausal_status")]
write.csv(filtered_metadata, file = "filtered_metadata.csv", row.names = TRUE)

#Check if didn't miss something and number of samples is the same for bacteria_matrix and metadata
dim(bacteria_matrix)
dim(filtered_metadata)

#Convert matrix to the right format
# phyloseq needs a "taxonomy table" and an "OTU table"
# For this data, we can just use the bacteria_matrix as the "OTU table"
OTU <- otu_table(bacteria_matrix, taxa_are_rows = TRUE)

#Convert filtered_metadata_df to dataframe type
filtered_metadata_df <- as.data.frame(filtered_metadata)

#Convert your metadata to the right format
META <- sample_data(filtered_metadata_df)

#Create the phyloseq object!
#This bundles everything together.
physeq <- phyloseq(OTU, META)

#Calculate Alpha-Diversity
plot_richness(physeq, x = "menopausal_status", measures = c("Shannon", "Simpson"))

alpha_data <- estimate_richness(physeq, measures = c("Shannon", "Simpson"))
head(alpha_data)

#Calculate Bray-Curtis distance (phyloseq has a function)
bray_dist <- phyloseq::distance(physeq, method = "bray")

#Prepare data for PERMANOVA (Permutational multivariate analysis of variance)
adonis_meta <- sample_data(physeq)

#Manually build a *new* data.frame from its components
adonis_meta_df <- data.frame(
  menopausal_status = adonis_meta$menopausal_status,
  row.names = rownames(adonis_meta) # Critically, we copy the row names
)

#Now, run adonis2 using the new data.frame
adonis_result <- adonis2(
  bray_dist ~ menopausal_status,
  data = adonis_meta_df  # <-- Use the .df version here
)

#Print the results
print(adonis_result)