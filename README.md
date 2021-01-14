
---
title: "Parsing total and Phospho data"
author:
- name: Oliver Crook
date: '`r format(Sys.Date(), "%B %d, %Y")`'
abstract: Some abstract here
output:
  html_notebook:
    fig_caption: yes
    fig_retina: null
    code_folding: hide
    toc: yes
    toc_float: TRUE
    df_print: paged
---

Needed packages
```{r,}
require(tidyverse)
require(camprotR)
require(MSnbase)
```



Here we read in the (peptide-spectrum-match) PSM level data into R.
```{r,}
psm <- read.delim("Data/Anja_benchmarking_total_proteome_method1_PSMs.txt")

```

check psm data is what we expect to be
```{r,}
head(psm)
```



Make the cRAP list for filtering (this is a list of common contaminates)

```{r,}
get_fasta_ids <- function(fasta){
  # Load the FASTA 
  bs.fasta <- Biostrings::fasta.index(fasta, seqtype = "AA")
  
  # Extract the UniProt accessions 
  accessions <- bs.fasta %>% 
    pull(desc) %>% 
    stringr::str_extract_all("(?<=\\|).*?(?=\\|)") %>% 
    unlist()

  accessions
}

crap.accessions <- get_fasta_ids('shared_files/cRAP_FullIdentifiers.fasta')

```


We are working with yeast (sc) and human (hs) so we need to match to correct database
```{r,}
hs.accessions <- get_fasta_ids(
  'shared_files/CCP_SwissProt_homo_sapiens_proteome_without_isoforms_20180409.fasta')
sc.accessions <- get_fasta_ids(
  'shared_files/swissprot_yeast_database_june2020.fasta')

uniprot_2_species <- data.frame('id'=c(hs.accessions, sc.accessions),
                                'species'=c(rep('H.sapiens', length(hs.accessions)),
                                            rep('S.cerevisiae', length(sc.accessions))))
head(uniprot_2_species)

```


Now we need to parse data and filter PSMs that share PSMs with cRAP proteins

```{r,}

psm_parsed <- psm %>% parse_features(TMT=TRUE, level='PSM',
                 crap_proteins = crap.accessions, unique_master = FALSE)


```
Annotated tbe dataset with the correct species
```{r,}


psm_parsed_annt <- list(psm_parsed) %>% lapply(function(x){
  
 species_matches <- x %>% select(Protein.Accessions) %>%
  mutate(Protein.Accessions_sep=Protein.Accessions) %>%
  separate_rows(Protein.Accessions_sep) %>%
  merge(uniprot_2_species, by.x='Protein.Accessions_sep', by.y='id', all.x=TRUE) %>%
  group_by(Protein.Accessions) %>%
  summarise(all_species=paste0(unique(species), collapse='; ')) %>%
  mutate(species=ifelse(grepl(';', all_species), 'mixed', all_species))
 
 x %>% merge(species_matches, by='Protein.Accessions')
})

```
```{r,}
psm_parsed_annt %>% lapply(nrow) #check everything is still good 100119PSMs

```
Now check PSMs and proteins per species

```{r,}

  p1 <- psm_parsed_annt[[1]] %>%
    group_by(species) %>%
    tally() %>%
    ggplot(aes(species, n)) +
    geom_bar(stat='identity') +
    theme_camprot() +
    xlab('') +
    ylab('PSMs') + ggtitle("PSMs")
  
  p2 <- psm_parsed_annt[[1]] %>%
    select(Master.Protein.Accessions, species) %>%
    unique() %>%
    group_by(species) %>% tally() %>%
    ggplot(aes(species, n)) +
    geom_bar(stat='identity') +
    theme_camprot() +
    xlab('') +
    ylab('Proteins') +
    ggtitle("Proteins")
  
  print(p1)
  print(p2)
  
```
Put data into MSnSets

```{r,}
sample_info <- read.delim('Data/sample_info.tsv') %>% tibble::column_to_rownames('tag')

psm_res <- psm_parsed_annt %>% lapply(function(x){
  # Abundance columns for TMT PD-output start with Abundance 
  abundance_cols <- colnames(x)[grepl('Abundance.', colnames(x))]
  
  .e <- as.matrix(x[,abundance_cols])
  .f <- x[,setdiff(colnames(x), abundance_cols)]
  
  # update the column names to remove the 'Abundance.` prefix
  colnames(.e) <- gsub('Abundance.', '', colnames(.e))
  
  res <- MSnbase::MSnSet(exprs=.e, fData=.f, pData=sample_info)
  
  res
})

```

Plotting the distribution of tag intensities in each dataset and the single species subsets. Note that the tag intensities for yeast fall into the 3 groups we expect given the experimental design. The results are very similar suggesting succesful experiment

```{r,}

  all <- psm_res[[1]]
  hs <- all[fData(all)$species=='H.sapiens']
  sc <- all[fData(all)$species=='S.cerevisiae']
  
  slices <- list('All'=all, 'H.sapiens'=hs, 'S.cerevisiae'=sc)
  for(slice in names(slices)){
    p <- slices[[slice]] %>% log(base=2) %>% plot_quant() +
      ggtitle(sprintf('%s - %s', 1, slice)) +
      ylab('PSM intensity (log2)')
    print(p)
    
    p <- slices[[slice]] %>% log(base=2) %>% plot_quant(method='density') +
      xlab('PSM intensity (log2)') +
      ggtitle(sprintf('%s - %s', 1, slice))
    print(p)
  }



```
save for later
```{r,}
saveRDS(psm_res, 'results/psm_res.rds')
```

Now, we need to add in experimental design information
```{r,}
exp_design <- pData(psm_res[[1]]) %>%
  select(condition, S.cerevisiae=yeast, H.sapiens=human) %>%
  unique()
  
sc_spikes <- exp_design$S.cerevisiae
hs_spikes <- exp_design$H.sapiens

get_ground_truth <- function(sc_spikes, hs_spikes, ix_1, ix_2){
  comparison <- sprintf('%s vs %s', sc_spikes[ix_2], sc_spikes[ix_1])
  hs_ground_truth <- hs_spikes[ix_2]/hs_spikes[ix_1]
  sc_ground_truth <- sc_spikes[ix_2]/sc_spikes[ix_1]
  return(c(comparison, hs_ground_truth, sc_ground_truth))
}
require(gtools)

expected <- apply(permutations(n=3,r=2), 1, function(x){
  get_ground_truth(sc_spikes, hs_spikes, x[1], x[2])
}) %>% t() %>% data.frame() %>%
  setNames(c('comparison', 'H.sapiens', 'S.cerevisiae')) %>%
  mutate_at(vars(S.cerevisiae, 
                 H.sapiens), 
            funs(as.numeric)) %>%
  pivot_longer(-comparison, names_to='species', values_to='expected')
```
This gives the expected fold changes that that were induced by the experimental design
```{r,}
print(expected)
```
The crib for the conditions replicates and changes are given in the pData
```{r,}
pData(psm_res[[1]])

```
Hence in condition 2 vs 1 for human we expect to see 0.95 change for humans and 2 change for yeast.

# Aggregation to peptide or protein level

```{r,}
aggregate_psms <- function(psm, level='peptide'){
  
  if(level=='peptide'){
    feature_name <- 'Sequence'
    min_features <- 1
  } else{
    feature_name <- 'Master.Protein.Accessions'
    min_features <- 2
  }
  
  psm_filt <- psm %>% 
    MSnbase::filterNA(pNA = 0) %>% # remove PSMs with missing values 
    restrict_features_per_protein(min_features = min_features, # remove proteins with <2 PSMs
                                  master_protein_col=feature_name, plot=FALSE) 
  
  group <- MSnbase::fData(psm_filt)[[feature_name]]

  agg <- psm_filt %>% MSnbase::combineFeatures(groupBy = group, method = 'sum') # Note we are summing intesnisty other methods are available

  feature_counts <- count_features_per_protein(psm_filt)
  
  return(list('agg'=agg, 'feature_counts'=feature_counts))
}

count_features_per_protein <- function(obj, master_prot_col='Master.Protein.Accessions'){
  obj %>%
    exprs() %>%
    data.frame() %>%
    tibble::rownames_to_column('feature') %>%
    tidyr::pivot_longer(-.data$feature, names_to='sample') %>%
    dplyr::mutate(sample=remove_x(.data$sample)) %>%
    filter(is.finite(.data$value)) %>%
    merge(fData(obj)[,master_prot_col,drop=FALSE], by.x='feature', by.y='row.names') %>%
    group_by(.data$sample, !!sym(master_prot_col)) %>%
    tally()
}



```

aggregate to protien levle abundances

```{r,}

    # log center-median normalise PSMs before aggregation (this controls for the PSM level effect that would be a random effect in a linear model)
    x <- psm_res[[1]] %>% log(base=2) %>%
      MSnbase::normalise(method='diff.median')
    exprs(x) <- 2^exprs(x) # then exponentiate to get untransformed intensities
    
prot_res <- aggregate_psms(x, level='protein')



```

aggregate to peptide level

```{r,}

    # log center-median normalise PSMs before aggregation (this controls for the PSM level effect that would be a random effect in a linear model)
    x <- psm_res[[1]] %>% log(base=2) %>%
      MSnbase::normalise(method='diff.median')
    exprs(x) <- 2^exprs(x) # then exponentiate to get untransformed intensities
    
pep_res <- aggregate_psms(x, level='peptide')



```


save data for later use
```{r,}
saveRDS(prot_res, "results/prot_res.rds")
saveRDS(pep_res, "results/pep_res.rds")

```

We can see the induced fold changes using a boxplot at protein level or peptide level
```{r,}
boxplot(log2(exprs(prot_res$agg)[fData(prot_res$agg)$species == "S.cerevisiae",]))

```
```{r,}
boxplot(log2(exprs(pep_res$agg)[fData(pep_res$agg)$species == "S.cerevisiae",]))

```
There are small changes for human proteins by they are VERY small.
```{r,}
boxplot(log2(exprs(pep_res$agg)[fData(pep_res$agg)$species == "H.sapiens",]))
```
Now let's look at the phosphoPeptide data again reading in at PSM level

```{r,}
psm_phos <- read.delim("Data/Anja_benchmarking_phospho_proteome_PSMs.txt")

```

double check
```{r,}
head(psm_phos)
```
Parse and filter PSMs to remove cRAP protiens

```{r,}
psm_phos_parsed <- psm_phos %>% parse_features(TMT = TRUE, level = "PSM",
                                               crap_proteins = crap.accessions,
                                               unique_master = TRUE)

```
We need to annotated the data again with the correct species

```{r,}
 species_matches <- psm_phos_parsed %>% select(Protein.Accessions) %>%
  mutate(Protein.Accessions_sep=Protein.Accessions) %>%
  separate_rows(Protein.Accessions_sep) %>%
  merge(uniprot_2_species, by.x='Protein.Accessions_sep', by.y='id', all.x=TRUE) %>%
  group_by(Protein.Accessions) %>%
  summarise(all_species=paste0(unique(species), collapse='; ')) %>%
  mutate(species=ifelse(grepl(';', all_species), 'mixed', all_species))
 
psm_phos_parsed_annt <-  psm_phos_parsed %>% merge(species_matches, by='Protein.Accessions')  
  
  
  
```

check nothing weird happened and that the annotation is indeed there
```{r,}
nrow(psm_phos_parsed_annt)

#barplot with phospho PSMs per species
barplot(table(psm_phos_parsed_annt$species))

```
Time to make the MSnSets
```{r,}
sample_info_phos <- read.delim('Data/sample_info_phos.txt') %>% tibble::column_to_rownames('tag')

psm_phos_res <- list(psm_phos_parsed_annt) %>% lapply(function(x){
  # Abundance columns for TMT PD-output start with Abundance 
  abundance_cols <- colnames(x)[grepl('Abundance.', colnames(x))]
  
  .e <- as.matrix(x[,abundance_cols])
  .f <- x[,setdiff(colnames(x), abundance_cols)]
  
  # update the column names to remove the 'Abundance.` prefix
  colnames(.e) <- gsub('Abundance.', '', colnames(.e))
  
  res <- MSnbase::MSnSet(exprs=.e, fData=.f, pData=sample_info_phos)
  
  res
})

```

Plot itensites, note that humans have very small change in itensity between the two scenarios. Yeast has clear change between the two conditions

```{r,}
  all <- psm_phos_res[[1]]
  hs <- all[fData(all)$species=='H.sapiens']
  sc <- all[fData(all)$species=='S.cerevisiae']
  
  slices <- list('All'=all, 'H.sapiens'=hs, 'S.cerevisiae'=sc)
  for(slice in names(slices)){
    p <- slices[[slice]] %>% log(base=2) %>% plot_quant() +
      ggtitle(sprintf('%s - %s', 1, slice)) +
      ylab('PSM intensity (log2)')
    print(p)
    
    p <- slices[[slice]] %>% log(base=2) %>% plot_quant(method='density') +
      xlab('PSM intensity (log2)') +
      ggtitle(sprintf('%s - %s', 1, slice))
    print(p)
  }



```
Save for downstrema results
```{r,}
saveRDS(psm_phos_res, "results/psm_res")

```

Want to filter out proteins with poor phospho scores
```{r,}
phospho_psm <- parse_PTM_scores(psm_phos_parsed_annt, ptm_col="PhosphoRS.Best.Site.Probabilities", threshold=7)

```
```{r,}
phospho_psm_filtered <- phospho_psm %>% filter(filtered_score!="") # The filtering doesn't actually take place until this step 

```

Add useful phosphorylation information about the sequence

```{r,}

protein_fasta_inf <- 'shared_files/combined_species.fasta'


phospho_psm_filtered <- phospho_psm_filtered %>% add_PTM_positions(protein_fasta_inf)
phospho_psm_filtered <- add_site_sequence(phospho_psm_filtered, protein_fasta_inf)
phospho_psm_filtered <- phospho_psm_filtered %>% add_peptide_positions(protein_fasta_inf)
```


Time to make the MSnSets again
```{r,}
sample_info_phos <- read.delim('Data/sample_info_phos.txt') %>% tibble::column_to_rownames('tag')

psm_phos_res <- list(phospho_psm_filtered) %>% lapply(function(x){
  # Abundance columns for TMT PD-output start with Abundance 
  abundance_cols <- colnames(x)[grepl('Abundance.', colnames(x))]
  
  .e <- as.matrix(x[,abundance_cols])
  .f <- x[,setdiff(colnames(x), abundance_cols)]
  
  # update the column names to remove the 'Abundance.` prefix
  colnames(.e) <- gsub('Abundance.', '', colnames(.e))
  
  res <- MSnbase::MSnSet(exprs=.e, fData=.f, pData=sample_info_phos)
  
  res
})

```

Note experimatal design is here:
```{r,}
pData(psm_phos_res[[1]])
```

Get expected results etc
```{r,}
exp_design <- pData(psm_phos_res[[1]]) %>%
  select(condition, S.cerevisiae=yeast, H.sapiens=human) %>%
  unique()
  
sc_spikes <- exp_design$S.cerevisiae
hs_spikes <- exp_design$H.sapiens

get_ground_truth <- function(sc_spikes, hs_spikes, ix_1, ix_2){
  comparison <- sprintf('%s vs %s', sc_spikes[ix_2], sc_spikes[ix_1])
  hs_ground_truth <- hs_spikes[ix_2]/hs_spikes[ix_1]
  sc_ground_truth <- sc_spikes[ix_2]/sc_spikes[ix_1]
  return(c(comparison, hs_ground_truth, sc_ground_truth))
}
require(gtools)

expected_phos <- apply(permutations(n=3,r=2), 1, function(x){
  get_ground_truth(sc_spikes, hs_spikes, x[1], x[2])
}) %>% t() %>% data.frame() %>%
  setNames(c('comparison', 'H.sapiens', 'S.cerevisiae')) %>%
  mutate_at(vars(S.cerevisiae, 
                 H.sapiens), 
            funs(as.numeric)) %>%
  pivot_longer(-comparison, names_to='species', values_to='expected')
```

Match phos camprisons with total comparisons
```{r,}
expected_phos$comparison <- expected$comparison
```

Relative change in phospho compared to total, double check these number by hand!
```{r,}
expected_rel_phospho <- cbind(expected_phos[,c(1,2)], expected_phos$expected/expected$expected)
expected_rel_phospho
```
so let's talk through it. The comparison 2 vs 1, the yeast total doubles, and the yeast phospho increases 6 times. Hence the relative increase in phospho is 3. Now comparison 6 vs 1, the yeast total goes up by 6, and the yeast phospho increases 6 times, hence the relative increase is 1 (i.e. no change).

```{r,}
saveRDS(psm_phos_res, file = "results/psm_phos_res.rds")

```

Now wish wish to aggregate to peptide level and protein level
```{r,}
    # log center-median normalise PSMs before aggregation
    psm_phos_res <- psm_phos_res[[1]] %>% log(base=2) %>%
      MSnbase::normalise(method='diff.median')
    exprs(psm_phos_res) <- 2^exprs(psm_phos_res) # then exponentiate to get untransformed intensities
    
pep_phos_res <- aggregate_psms(psm_phos_res, level='peptide')

pep_phos_res

```
Let's us have a look at the phospho fold changes, note that log2(6) = 2.58 
```{r,}
boxplot(log2(exprs(pep_phos_res$agg)[fData(pep_phos_res$agg)$species == unique(fData(pep_phos_res$agg)$species)[1],]))
boxplot(log2(exprs(pep_phos_res$agg)[fData(pep_phos_res$agg)$species == unique(fData(pep_phos_res$agg)$species)[2],]))
```

Now wish wish to aggregate to peptide level and protein level
```{r,}
    # log center-median normalise PSMs before aggregation
    psm_phos_res <- psm_phos_res %>% log(base=2) %>%
      MSnbase::normalise(method='diff.median')
    exprs(psm_phos_res) <- 2^exprs(psm_phos_res) # then exponentiate to get untransformed intensities
    
prot_phos_res <- aggregate_psms(psm_phos_res, level='protein')

prot_phos_res

```
fold changes
```{r,}
boxplot(log2(exprs(prot_phos_res$agg)[fData(prot_phos_res$agg)$species == unique(fData(prot_phos_res$agg)$species)[1],]))
boxplot(log2(exprs(prot_phos_res$agg)[fData(prot_phos_res$agg)$species == unique(fData(prot_phos_res$agg)$species)[2],]))
```
```{r,}
Let look at protein level first but the analysis should probably be thought at the peptide level
```

```{r,}
yeastOI <- rownames(exprs(prot_phos_res$agg[fData(prot_phos_res$agg)$species == unique(fData(prot_phos_res$agg)$species)[2],]))

```

```{r,}
plot(log2(exprs(prot_res$agg)[rownames(prot_res$agg) %in% yeastOI,][1,]), col = "blue", pch = 19, ylim = c(5,14))
points(log2(exprs(prot_phos_res$agg)["P00359",]), col = "purple", pch = 19)

```
plot on fold change scale
```{r,}
plot(log2(exprs(prot_res$agg)[rownames(prot_res$agg) %in% yeastOI,][1,]) - mean(log2(exprs(prot_res$agg)[rownames(prot_res$agg) %in% yeastOI,][1,1:4])), col = "blue", cex = 2, pch = 19)
points(log2(exprs(prot_phos_res$agg)["P00359",]) - mean(log2(exprs(prot_phos_res$agg)["P00359",1:4])), col = "purple", pch = 19)


```
```{r,}
saveRDS(prot_phos_res, file = "prot_phos_res.rds")
saveRDS(pep_phos_res, file = "pep_phos_res.rds")
```
