# libraries 
suppressPackageStartupMessages(require(tidyverse))

# merge individual metaphlan tables
read_metaphlan_tables <- function(file_list) {
  print(paste('Loading', length(file_list), 'MetaPhlAn2 profiles.'))
  profiles  = data.frame()
  for (file_fp in file_list) {
    file_fn = basename(file_fp)
    profiles = suppressMessages(read_tsv(file_fp, col_names = c('taxonomy', 'rel_abund'), skip = 1)) %>% 
      mutate(file_fn = file_fn) %>% 
      bind_rows(profiles, .)
  }
  return(profiles)
}

# filter metaphlan table for taxonomy-level
filter_metaphlan_level <- function(metaphlan_table, taxonomy_level) {
  l = switch(taxonomy_level,
             kingdom = c('k__', 'p__'),
             phylum = c('p__', 'c__'), 
             class = c('c__', 'o__'),
             order = c('o__', 'f__'),
             family = c('f__', 'g__'),
             genus = c('g__', 's__'),
             species = c('s__', 't__'),
             strain = c('t__', ''))
  keep = l[1]
  drop = l[2]
  
  taxonomy_level.count = sum(grepl(metaphlan_table[, 'taxonomy'], pattern = keep))
  print(paste('Found', taxonomy_level.count, 'entries for', taxonomy_level))
  output_table = metaphlan_table %>% 
    filter(grepl(taxonomy, pattern = keep),
           !grepl(taxonomy, pattern = drop))
  return(output_table)
}

# split metaphlan taxonomy to columns
split_metaphlan_taxonomy <- function(metaphlan_table) {
  taxonomy_levels = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')
  n_split = str_count(metaphlan_table[1, 'taxonomy'], pattern = '\\|') + 1
  
  metaphlan_table.taxonomy = str_split_fixed(metaphlan_table[,'taxonomy'], pattern = '\\|', n = n_split) %>% 
    as.data.frame() %>% 
    set_names(taxonomy_levels[1:n_split]) %>% 
    mutate_all(.funs = ~ gsub(pattern = '^[k,p,c,o,f,g,s,t]__', replacement = '', .)) 
  
  output_table = cbind(metaphlan_table.taxonomy, 
                       select(metaphlan_table, setdiff(names(metaphlan_table), 'taxonomy')))
  return(output_table)
}

# wrapper function for loading metaphlan data on specific taxonomic level
get_metaphlan_level <- function(file_list, taxonomy_level, output_fn=F) {
  table = read_metaphlan_tables(file_list)
  table = filter_metaphlan_level(table, taxonomy_level)
  table = split_metaphlan_taxonomy(table)
  
  if(output_fn != F) write_tsv(table, output_fn)
  return(table)
}

