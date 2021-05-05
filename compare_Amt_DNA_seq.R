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

#how does average coverage change for species B
df %>%
  #remove ancestral species
  filter(., !grepl("*.Ancestral.*",ID)) %>%
  separate(., ID, into = c("community", "replicate", "Day", "species", "scaffold"), sep = "_") %>%
  group_by(community,replicate,Day,species) %>%
  filter(grepl("*B.*",community),grepl("*B-.*",species)) %>%
  group_by(community,replicate) %>%
  ggplot(.,aes(x = community,y=Average_Cov_Total,color = scaffold))+
  geom_point()+
  facet_grid(rows = vars(scaffold), cols = vars(replicate))

# Average Coverage per species per day
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
  mutate(Day = fct_relevel(Day, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(.,aes(x = richness,y = Avg_COV, color = species, fill = species)) +
  geom_point() +
  geom_boxplot() +
  facet_grid(cols = vars(Day))



# Average Coverage per species per day
df %>%
  #remove ancestral species
  filter(., !grepl("*.Ancestral.*",ID)) %>%
  filter(., !grepl("^ABC_R2.*|^ABC_R3.*|^AB_R1.*|^AB_R3.*|^C_R1.*",ID)) %>%
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
  mutate(Day = fct_relevel(Day, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(.,aes(x = richness,y = Avg_COV, color = species, fill = species)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2) +
  geom_boxplot(alpha = 0.1) +
  ggtitle("Average Number of Cells Sequenced")+
  theme(plot.title = element_text(hjust = 0.5))

# Total DNA sequenced per species
df %>%
  #remove ancestral species
  filter(., !grepl("*.Ancestral.*",ID)) %>%
  filter(., !grepl("^ABC_R2.*|^ABC_R3.*|^AB_R1.*|^AB_R3.*|^C_R1.*",ID)) %>%
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
  mutate(Day = fct_relevel(Day, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(.,aes(x = richness,y = Total_DNA, color = species, fill = species)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2) +
  geom_boxplot(alpha = 0.1) +
  ggtitle("Total Amount of DNA sequenced per species")+
  theme(plot.title = element_text(hjust = 0.5))



