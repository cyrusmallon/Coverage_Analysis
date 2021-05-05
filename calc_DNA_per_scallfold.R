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

###NEW DF with added columns
data_cov <-
df %>%
  #remove ancestral species
  filter(., !grepl("*.Ancestral.*",ID)) %>%
  separate(., ID, into = c("community", "replicate", "Day", "species", "scaffold"), sep = "_") %>%
  #via the group, different scaffolds and their various levels of coverage are taken into account
  group_by(community,replicate,Day,species) %>%
  summarise(n = n(), 
            #Avg Coverage Across All Scaffolds;thus per genome, accounting for differences in coverage between scaffolds
            #Avg_Cov is another measure for average number of cells
            Avg_COV = mean(Average_Cov_Total, na.rm = TRUE), 
            #Avg_Cells lacks weighing each scaffolds coverage
            Avg_Cells = sum(Sum_bp_Total)/sum(length_scaffold), 
            #Total_DNA per species
            Total_DNA = sum(Sum_bp_Total), 
            length_genome = sum(length_scaffold)) %>%
  #print(n = 200) %>%
  mutate(richness = 
           ifelse(nchar(community) == 1, "1",
                  ifelse(nchar(community) == 2, "2",
                         ifelse(nchar(community) == 3, "3",
                                ifelse(nchar(community) == 4, "4"))))) %>%
  mutate(species = sub("-Prokka","",species)) %>%
  mutate(mut_sp_origin = species) %>%
  unite("Label",1:3,sep = "_", remove=FALSE) %>%
  mutate(Label = str_replace(Label, "r", "R")) %>%
  mutate(replicate = str_replace(replicate, "r", "R")) %>%
  #filter contaminated samples
  filter(., !grepl("^ABC_R2.*|^ABC_R3.*|^AB_R1.*|^AB_R3.*|^C_R1.*", Label)) %>%
  print(n = 100)


save(data_cov, file = "data_cov.RData")


