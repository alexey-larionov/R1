
# ==== Done till here ===== #

```{r todo}
samples.df %>% 
  filter(grepl("not_annotated", filter)) %>% 
  arrange(filter)

not_annotated <-
  samples.df %>% 
  filter(grepl("not_annotated", filter)) %>% 
  select(gwas_id)

not_annotated <- as.vector(not_annotated$gwas_id)

samples.df %>% 
  filter(gwas_id %in% not_annotated)

demographics.df %>% 
  select(labid.x,labid.y) %>% 
  filter(labid.x %in% not_annotated | labid.y %in% not_annotated)

demographics.df %>% 
  select(starts_with("labid")) %>% 
  filter(is.na(labid.x) | is.na(labid.y))

demographics.df %>% 
  select(starts_with("labid")) %>% 
  filter(labid.x %in% wes.done.not.annotated | labid.y %in% wes.done.not.annotated)



sum(chemo != chemo_self_mra)
sum(hormone != hormone_self_mra, na.rm=TRUE)
sum(status2 != status2)
sum(race)
sum(dsmiss)
sum(good_location)
sum(Deleterious)
sum(family_history)
sum(XRTBrCHAR != XRTBreast)
sum(sub_dx_age_cat)

detach(pheno.df)

```

# reshape_demographics



# select_relevant_phenotypes_annotations

```{r select_relevant_phenotype_annotations}

phenotypes.df <- 
  covar.df %>% 
  select(labid, cc, 
         sub_dx_age, offset, family_history, 
         chemo, hormone, XRTBreast, dose, 
         Eigen_1, Eigen_2, Eigen_3, Eigen_4, Eigen_5)

colnames(phenotypes.df)[1] <- "gwas_id" # for merging with samples.df

rm(covar.df)

```



# add_WES_cases_ids_to_phenotypes_table

```{r add_WES_cases_ids_to_phenotypes_table}

# Explore cases submitted for WES
sort(table(samples.df$gwas_id))

# Explore cases annotated in GWAS
sort(table(phenotypes.df$gwas_id))

# Merge samples and phenotypes tables
merged.df <- merge(samples.df, phenotypes.df, by="gwas_id")
dim(merged.df)

# Explore cases in the merged set cases missed and measured twice in wes
sort(table(merged.df$gwas_id))

# Cases without wes results (failes sequencing)
sort(table(merged.df$gwas_id))[ sort(table(merged.df$gwas_id)) == 0 ]

# Cases sequenced twice
sort(table(merged.df$gwas_id))[ sort(table(merged.df$gwas_id)) == 2 ]

# Look closer into two samples that had been measured twice
merged.df[merged.df$gwas_id %in% c("id259643", "id272715"),]

# For case sequenced twice: exclude measurement with lower coverage (as in lab records)
cases_to_exclude <- merged.df$merged_id %in% c("P4_F01_id259643", "P5_D07_id272715")
phenotypes.df <- merged.df[!cases_to_exclude, ]
dim(phenotypes.df)

# Sort phenotypes table
phenotypes.df <- phenotypes.df %>% arrange(wes_id)

# Clean-up
rm(merged.df, samples.df, cases_to_exclude)

```

At this point phenotypes table contains information about 498 samples.  
These are all annotated samples placed on WES plates.  
12 of 512 WES samples were removed because they had no phenotype annotations.
2 samples were removed because they were placed to the plates in duplicate.

# reshape_genotypes_table

```{r reshape_genotypes_table}

# Rename genotype table
genotypes.mx <- gt.pf.msk.lof.mx
rm(gt.pf.msk.lof.mx)

# Update column names in genotype table
colnames(genotypes.mx) <- sub(".GT$","",colnames(genotypes.mx))
colnames(genotypes.mx)

# Get IDs of annotated cases (from GWAS dataset)
annotated.cases <- as.vector(phenotypes.df$wes_id)
length(annotated.cases)

# Get IDs of successfully sequenced cases (from WES genotypes table)
sequenced.cases <- colnames(genotypes.mx)
length(sequenced.cases)

# Cases to analyse
selected.cases <- sort(intersect(annotated.cases, sequenced.cases))
length(selected.cases)
rm(annotated.cases, sequenced.cases)

# Select cases with annotated phenotypes
genotypes.mx <- genotypes.mx[,selected.cases]
dim(genotypes.mx)

# Use merged IDs for column names in genotype table
rownames(phenotypes.df) <- phenotypes.df$wes_id
colnames(genotypes.mx) <- phenotypes.df[selected.cases,]$merged_id

```

# update_phenotypes_table

```{r update_phenotypes_table}

# Remove ceses that failed sequencing
phenotypes.df <- phenotypes.df[selected.cases,]

# Reorder columns
phenotypes.df <- phenotypes.df[,c("merged_id", "wes_id", "gwas_id", "cc", "offset", "family_history", "sub_dx_age", "chemo", "hormone", "XRTBreast", "Eigen_1", "Eigen_2", "Eigen_3", "Eigen_4", "Eigen_5")]

# Place merged id to rownames
rownames(phenotypes.df) <- phenotypes.df[,1]

# Clean-up
rm(selected.cases)

```

At this point phenotypes and genotypes tables contains information about 483 samples.  
These are all annotated samples which were successfully sequenced:  
  - 2 samples were removed because they were placed to the plates in duplicate
- 12 of 512 WES samples were removed because they had no phenotype annotations
- 15 annotated cases failed sequencing

Variants and genotypes tables contain information about 5,439 high impact variants 
that have a potential to cause loss of the protein functuon.
