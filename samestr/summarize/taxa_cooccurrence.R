suppressPackageStartupMessages(require(tidyverse))

mp_long_to_matrix <- function(mp.long, taxonomic_level=NA) {
  mp.matrix <- 
    mp.long %>% 
    select(all_of(taxonomic_level), Name, rel_abund) %>% 
    filter(across(all_of(taxonomic_level), ~ !grepl(.x, pattern = 'unclassified'))) %>% 
    group_by_at(vars(Name, all_of(taxonomic_level))) %>% 
    summarize(rel_abund = as.numeric(sum(rel_abund, na.rm = T) > 0), .groups = 'drop') %>% 
    pivot_wider(id_cols = 'Name', names_from = all_of(taxonomic_level), values_from ='rel_abund') %>% 
    mutate_at(.vars = vars(everything(), -Name), 
              .funs = ~ replace_na(as.numeric(. > 0), 0)) %>%
    column_to_rownames('Name') %>% 
    as.matrix()
  return(mp.matrix)
}

mp_cooccurrence <- function(mp_species.long) {
  tax_counts <- NULL
  cooc_longs <- NULL
  print(paste0('Evaluating taxonomic co-occurrence in ', length(unique(mp_species.long$Name)), ' Samples'))
  for (tax_lev in c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')) {
    mp.matrix <- mp_long_to_matrix(mp_species.long, taxonomic_level = tax_lev)
    mp.cooc <- mp.matrix %*% t(mp.matrix)
    mp.cooc.long <- 
      distmat_to_long(distmat = as.data.frame(mp.cooc), 
                      value_name = paste0('n.shared'),
                      rm_diag = T) %>% 
      mutate(taxonomy = tax_lev)
      
    cooc_longs <- 
      bind_rows(cooc_longs, mp.cooc.long)
    
    mp.tax_count <- 
      data.frame(n = diag(mp.cooc)) %>% 
      rownames_to_column('Name') %>% 
      mutate(taxonomy = tax_lev)
    tax_counts = bind_rows(tax_counts, mp.tax_count)
  }
  
  return(list('tax_counts' = tax_counts, 'tax_cooc' = cooc_longs))
}
