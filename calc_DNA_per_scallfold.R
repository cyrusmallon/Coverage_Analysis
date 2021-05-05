library(tidyverse)

#load coverage data
#loaded data is named <df>
load("total_bp_sequenced.RData")

#examind df
head(df)

#remove first and empty row
df <- df[-1,]

#remove unneeded text from ID
df$ID
df$ID <- sub(".coverage.*", "", df$ID)

#remove double underscore (for single underscore)
#for easy separation
df$ID <- sub("__", "_", df$ID)

#fix D13 to D14
df$ID <- sub("D13", "D14", df$ID)

#examine data
df %>%
  as_tibble() %>%
  print(n = 100)

#calc approx # of cells (i.e., genomes) sequenced
data_seq_DNA <-
  df %>%
  #remove ancestral species
  filter(., !grepl("*.Ancestral.*",ID)) %>%
  separate(., ID, into = c("community", "replicate", "Day", "species", "scaffold"), sep = "_") %>%
  group_by(community,replicate,Day,species) %>%
  summarize(n= n(), length_genome = sum(length_scaffold),total_bp_seq = sum(Sum_bp_Total)) %>%
  mutate(Approx_Cells = total_bp_seq/length_genome) %>%
  mutate(species = sub("-Prokka","",species)) %>%
  mutate(mut_sp_origin = species) %>%
  unite("Label",1:3,sep = "_", remove=FALSE) %>%
  mutate(Label = str_replace(Label, "r", "R")) %>%
  mutate(replicate = str_replace(replicate, "r", "R")) %>%
  #filter contaminated samples
  filter(., !grepl("^ABC_R2.*|^ABC_R3.*|^AB_R1.*|^AB_R3.*|^C_R1.*", Label))

save(data_seq_DNA, file = "data_seq_DNA.RData")





