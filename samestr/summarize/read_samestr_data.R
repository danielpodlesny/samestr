# libraries
suppressPackageStartupMessages(require(tidyverse))

# functions
read_distmat = function(directory, species, suffix) {
  fn = paste0(directory, '/', species, suffix)
  distmat = suppressWarnings(read.delim(fn, '\t', header=T, row.names=1, stringsAsFactors=F, check.names = F))
  return(distmat)
}

distmat_nonfinite_to_NA = function(distmat) {
  "
  Function takes distance matrix as data frame as input
  transforms non-finite values to NA, 
  returns distance matrix as data frame.
  "
  col_names <- colnames(distmat)
  row_names <- rownames(distmat)
  distmat <- do.call(data.frame, lapply(distmat, function(x) replace(x, !is.finite(x), NA)))
  colnames(distmat) <- col_names
  rownames(distmat) <- row_names
  return(distmat)
}

distmat_rm_upper_tri = function(distmat, rm_diag=FALSE) {
  "
  Function takes distance matrix as data frame as input
  and returns distance matrix as data frame without upper triangle.
  "  
  if (rm_diag == TRUE) {
    distmat[upper.tri(distmat, diag = TRUE)] <- NA
  } else {
    distmat[upper.tri(distmat, diag = FALSE)] <- NA
  }
  return(distmat)
}

distmat_to_long = function(distmat, value_name, rm_diag = TRUE) {
  "
  Function takes distance matrix as data frame as input
  and returns a data frame in long format,
  containing all values from the distance matrix:
  
  Requires: 
  At least lower triangle of distance matrix,
  modify_distmat script for distance matrix formatting.
  
  -- Output Format
  row col distance
  A B 2.5
  "
  distmat.long <- 
    distmat_rm_upper_tri(distmat, rm_diag = rm_diag) %>% 
    rownames_to_column('row') %>% 
    pivot_longer(names_to = 'col', 
                 values_to = value_name, 
                 cols = -matches('row'), 
                 values_drop_na = T)
  return(distmat.long)
}

get_merged_distmats = function(directory, species) {
  "
  Takes directory and species as input
  and returns a data.frame where file pairs for each species are merged.
  
  Requires:
  Input dir must contain samestr compare output files with
  extensions: .fraction.txt and .overlap.txt
  "
  distmat = read_distmat(directory, species, '.fraction.txt')
  overlap = read_distmat(directory, species, '.overlap.txt')
  
  # tag nan, inf as NA
  distmat = distmat_nonfinite_to_NA(distmat)
  overlap = distmat_nonfinite_to_NA(overlap)
  
  # wide to long
  distmat.long = distmat_to_long(distmat, value_name='similarity', rm_diag=T)
  overlap.long = distmat_to_long(overlap, value_name='overlap', rm_diag=T)
  
  if (nrow(distmat.long) & nrow(overlap.long)) {
    merged_distmat_overlap <- 
      
      full_join(distmat.long, overlap.long, 
                by = c('row','col')) %>% 
      mutate(species = species) %>% 
      filter(row != col)
  }
  return(merged_distmat_overlap)
}

# main
read_samestr_data <- function(sstr_dir) {
  "
  Takes samestr compare output directory as input and runs get_merged_distmats
  "
  distance_species_list <- gsub(list.files(path = sstr_dir, 
                                           pattern = '.fraction.txt', full.names = F), 
                                pattern = '.fraction.txt', replacement = '')
  overlap_species_list <- gsub(list.files(path = sstr_dir, 
                                          pattern = '.overlap.txt', full.names = F), 
                               pattern = '.overlap.txt', replacement = '')
  
  species_intersect = intersect(overlap_species_list, distance_species_list)
  print(paste0('Merging matrix for ', length(species_intersect), ' species.'))
  sstr_data = data.frame()
  for (species in species_intersect) {
    sstr_data = bind_rows(sstr_data, get_merged_distmats(directory=sstr_dir, species=species))
  }
  return(sstr_data)
}

