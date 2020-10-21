#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(tidyverse))

option_list = list(
  make_option(c("--input-dir"), action="store", default=NA, type='character',
              help = 'Path to `samestr compare` output directory'),
  make_option(c("--mp-profiles-dir"), action="store", default=NA, type='character',
              help = 'Path to directory with metaphlan profiles'),
  make_option(c("--mp-profiles-extension"), action="store", default='.profile.txt', type='character',
              help = 'File extension of metaphlan profiles'),
  make_option(c("--output-dir"), action="store", default=NA, type='character',
              help = 'Path to output directory'),
  make_option(c("--aln-pair-min-overlap"), action="store", default=5000, type='numeric',
              help = 'Minimum number of overlapping positions to evaluate alignment similarity'),
  make_option(c("--aln-pair-min-similarity"), action="store", default=0.999, type='numeric',
              help = 'Minimum pairwise alignment similarity to define a shared strain')
)
opt = parse_args(OptionParser(option_list=option_list))
if (any(unlist(lapply(opt, is.na)))) {
  stop(paste0('Missing arguments: ', 
              paste0(names(opt[unlist(lapply(opt, is.na))]), collapse = ', ')))
}

sstr_dir = opt$`input-dir`
mp_dir = opt$`mp-profiles-dir`
mp_ext = opt$`mp-profiles-extension`
min_overlap = opt$`aln-pair-min-overlap`
min_similarity = opt$`aln-pair-min-similarity`
output_dir = opt$`output-dir`

# mkdir -p output_dir
dir.create(file.path(output_dir), showWarnings = FALSE)

# load scripts
args = commandArgs(trailingOnly=F)
script.dir <- dirname(gsub(args[grepl(args, pattern = '^--file')],
                pattern = '--file=', replacement = ''))

source(paste0(script.dir, '/read_samestr_data.R'))
source(paste0(script.dir, '/read_metaphlan_data.R'))
source(paste0(script.dir, '/taxa_cooccurrence.R'))

# read metaphlan data
mp_file_list <- list.files(mp_dir, pattern = mp_ext, full.names = T)
if (!length(mp_file_list)) {
  stop(paste0('Could not find MetaPhlAn input files at: ', paste0(mp_dir, '*', mp_ext)))
}
mp_species.long <- get_metaphlan_level(file_list = mp_file_list, taxonomy_level = 'species') %>% 
  mutate(Name = str_remove(file_fn, pattern = mp_ext)) %>% 
  select(-file_fn)
mp_taxonomy <- mp_species.long %>% 
  select(kingdom:species) %>% 
  distinct()

# read samestr data
sstr_data = read_samestr_data(sstr_dir) %>% 
  mutate(analysis_level = ifelse(overlap > min_overlap, 'strain-level', 'species-level'), 
         shared = overlap > min_overlap & similarity > min_similarity, 
         event = case_when(
           shared ~ 'shared_strain',
           !shared & overlap > min_overlap ~ 'other_strain',
           T ~ 'same_species')) %>% 
  select(row, col, species, 
         similarity, overlap, 
         analysis_level, shared, event)

# taxa data
## taxa cooccurrence
### kingdom to species-level
cooccurrence <- mp_cooccurrence(mp_species.long)

### strain-level
strain_cooc <- sstr_data %>% 
  group_by(row, col) %>% 
  summarize(n.shared = sum(shared == T, na.rm = T), 
            n.analyzed_strains = sum(overlap > min_overlap, na.rm = T), 
            .groups = 'drop') %>% 
  mutate(taxonomy = 'strain') %>%
  ungroup()


## taxa counts
### kingdom to species-level
mp_counts <- 
  cooccurrence$tax_counts %>% 
  pivot_wider(names_from = 'taxonomy', 
              names_prefix = 'n.',
              values_from = 'n', 
              values_fill = 0)
  
## merge
cooc <- 
  bind_rows(cooccurrence$tax_cooc, strain_cooc) %>% 
  mutate(n.analyzed_strains = replace_na(n.analyzed_strains, 0)) %>% 
  pivot_wider(names_from = 'taxonomy', 
              names_prefix = 'n.shared_',
              values_from = 'n.shared', 
              values_fill = list(n.shared = 0),
	      values_fn = sum) %>% 
  group_by(row, col) %>%
  summarize_all(.funs = sum) %>%
  relocate(n.analyzed_strains, .after = last_col())

# save output
print('Saving output tables.')
write_tsv(cooc, path = paste0(output_dir, 'sstr_cooc.tsv'))
write_tsv(mp_counts, path = paste0(output_dir, 'mp_counts.tsv'))
write_tsv(mp_species.long %>% select(Name, species, rel_abund), paste0(output_dir, path = 'mp_species.long.tsv'))
write_tsv(mp_taxonomy, path = paste0(output_dir, 'mp_taxonomy.tsv'))
write_tsv(sstr_data, path = paste0(output_dir, 'sstr_data.tsv'))
print(paste0('Files saved to ', output_dir))
