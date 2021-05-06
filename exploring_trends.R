

load(file = "finaldata_5May2021.RData")

finaldata[,c(1,2,3,4, 10:ncol(finaldata))]

finaldata %>%
  as_tibble() %>%
  glimpse()

dat <- finaldata %>%
  as_tibble() %>%
  mutate(
    Mutations_Nrmlz = Mutations / Avg_COV,
    Mutations_S_Nrmlz = Mutations_S / Avg_COV,
    Mutations_cumulative_Nrmlz = Mutations_cumulative / Avg_COV,
    dN_cumulative_Nrmlz = dN_cumulative / Avg_COV,
    dS_cumulative_Nrmlz = dS_cumulative / Avg_COV,
    Mutations_Pos_Selection_Nrmlz = ((Mutations_cumulative/Avg_COV) - (dS_cumulative/Avg_COV))
  ) %>%
  #dN are all mutations, not just SNPs
  mutate(
    dNS_Nrmlz = Mutations_Pos_Selection_Nrmlz /dS_cumulative_Nrmlz
  )

#save data
save(dat,file = "finaldata_6May2021")

dat %>%
  glimpse()


#All cumulative Mutations
dat %>%
  filter(Mutations_cumulative_Nrmlz < 20) %>%
  mutate(Day.x = fct_relevel(Day.x, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(., aes(x = richness, y = Mutations_cumulative_Nrmlz,color = mut_sp_origin, fill = mut_sp_origin))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2)+
  geom_boxplot(alpha = 0.2)+
  ylab("Normalized Cumulative Mutations")+
  ggtitle("All Cumulative and Normalized Mutations")+
  facet_grid(cols = vars(Day.x), rows = vars(mut_sp_origin),scales = "free_y")+
  theme(plot.title = element_text(hjust = 0.5))

#All mutations normalized
dat %>%
  filter(Mutations_Nrmlz < 10) %>%
  mutate(Day.x = fct_relevel(Day.x, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(., aes(x = richness, y = Mutations_Nrmlz,color = mut_sp_origin,fill = mut_sp_origin))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2)+
  geom_boxplot(alpha = 0.2)+
  facet_grid(cols = vars(Day.x), rows = vars(mut_sp_origin))

#All nonsynon mutations normalized
dat %>%
  filter(Mutations_N_Nrmlz < 1.5) %>%
  mutate(Day.x = fct_relevel(Day.x, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(., aes(x = richness, y = Mutations_N_Nrmlz,color = mut_sp_origin,fill = mut_sp_origin))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2)+
  geom_boxplot(alpha = 0.2)+
  facet_grid(cols = vars(Day.x), rows = vars(mut_sp_origin))

#all synon mutations normalized
dat %>%
  filter(Mutations_S_Nrmlz < 1.5) %>%
  mutate(Day.x = fct_relevel(Day.x, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(., aes(x = richness, y = Mutations_S_Nrmlz,color = mut_sp_origin,fill = mut_sp_origin))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2)+
  geom_boxplot(alpha = 0.2)+
  facet_grid(cols = vars(Day.x), rows = vars(mut_sp_origin))

#all cumulative nonsynon mutations normalized
dat %>%
  filter(dN_cumulative_Nrmlz    < 2) %>%
  mutate(Day.x = fct_relevel(Day.x, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(., aes(x = richness, y = dN_cumulative_Nrmlz   ,color = mut_sp_origin,fill = mut_sp_origin))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2)+
  geom_boxplot(alpha = 0.2)+
  facet_grid(cols = vars(Day.x), rows = vars(mut_sp_origin))

#all cumulative synon mutations normalized
dat %>%
  filter(dS_cumulative_Nrmlz    < 1.5) %>%
  mutate(Day.x = fct_relevel(Day.x, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(., aes(x = richness, y = dS_cumulative_Nrmlz   ,color = mut_sp_origin,fill = mut_sp_origin))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2)+
  geom_boxplot(alpha = 0.2)+
  facet_grid(cols = vars(Day.x), rows = vars(mut_sp_origin))

#dN/dS ratio, with normalized mutations
dat %>%
  filter(dNS_Nrmlz < 3) %>%
  mutate(Day.x = fct_relevel(Day.x, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(., aes(x = richness, y = dNS_Nrmlz  ,color = mut_sp_origin,fill = mut_sp_origin))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2)+
  geom_boxplot(alpha = 0.2)+
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  ggtitle("dN/dS (with normalized # mutatations)")+
  ylab("dN/dS ratio")+
  facet_grid(cols = vars(Day.x), rows = vars(mut_sp_origin),scales = "free_y")+
  theme(plot.title = element_text(hjust = 0.5))

#Positive Selection Mutations by day and species
#Thus not just nonynon mutations but also Indels and other structural variants
dat %>%
  filter(Mutations_Pos_Selection_Nrmlz  < 15) %>%
  mutate(Day.x = fct_relevel(Day.x, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(., aes(x = richness, y = Mutations_Pos_Selection_Nrmlz  ,color = mut_sp_origin,fill = mut_sp_origin))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2)+
  geom_boxplot(alpha = 0.2)+
  ggtitle("All Non-Neutral Mutations")+
  ylab("Normalized and Non-Neutral Mutations")+
  facet_grid(cols = vars(Day.x), rows = vars(mut_sp_origin),scales = "free_y")+
  theme(plot.title = element_text(hjust = 0.5))

#All cumulative and normalized mutations
dat %>%
  mutate(Day.x = fct_relevel(Day.x, c("D7", "D14", "D21","D28","D31"))) %>%
  filter(Mutations_cumulative_Nrmlz < 8) %>%
  ggplot(., aes(x = richness, y = Mutations_cumulative_Nrmlz,color = mut_sp_origin, fill = mut_sp_origin))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2)+
  geom_boxplot(alpha = 0.2)+
  facet_grid(cols = vars(Day.x), rows = vars(mut_sp_origin), scales = "free_y")+  
  
#All Positive via richness
dat %>%
  filter(Mutations_Pos_Selection_Nrmlz  < 10) %>%
  mutate(Day.x = fct_relevel(Day.x, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(., aes(x = richness, y = Mutations_Pos_Selection_Nrmlz, color = richness ))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2)+
  geom_boxplot(alpha = 0.2)
  #facet_grid(cols = vars(Day.x), rows = vars(mut_sp_origin),scales = "free_y")

#All Cumulative Mutations via Richness
dat %>%
  filter(Mutations_cumulative_Nrmlz < 6) %>%
  mutate(Day.x = fct_relevel(Day.x, c("D7", "D14", "D21","D28","D31"))) %>%
  ggplot(., aes(x = richness, y = Mutations_cumulative_Nrmlz, color = richness ))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2)+
  geom_boxplot(alpha = 0.2)
  #facet_grid(cols = vars(Day.x), rows = vars(mut_sp_origin),scales = "free_y")


