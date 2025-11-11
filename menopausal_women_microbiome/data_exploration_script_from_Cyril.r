# install.packages("BiocManager")
# BiocManager::install("curatedMetagenomicData")

library(curatedMetagenomicData)

str(sampleMetadata) # 22588 obs. of 141 variables


colnames(sampleMetadata) # names of the 141 variables


table(sampleMetadata$menopausal_status) # samples with menopausal status
sample_meno <- sampleMetadata[
  !is.na(sampleMetadata$menopausal_status),
]

# quick overview of the 431 menopausal samples

table(sample_meno$study_name) # they come from three studies
table(sample_meno$PMID) # should look into these studies
length(unique(sample_meno$subject_id)) # they come from 426 subjects/people
table(sample_meno$body_site) # all stool samples
table(sample_meno$disease) # 342 healthy samples, some diseased samples
table(sample_meno$age) # from 19 to 80 years old with non-uniform coverage
table(sample_meno$gender) #330 female, 86 male??
table(sample_meno$country) #mostly from the UK
table(sample_meno$treatment) # a few under treatment
table(sample_meno$pregnant) # some are pregnant
table(sample_meno$smoker) # some are smoker
table(sample_meno$BMI) # probably something to correct for
table(sample_meno$alcohol)
