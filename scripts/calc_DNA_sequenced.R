#Script by: Cyrus A. Mallon
#Purpose: To calculate the total number of basepairs sequenced for each scaffold,
#later this info will be used to normalize the number of mutations per species
#by the total amount of DNA that was sequenced for that particular species

#load libraries
library(here)
library(tidyverse)

#view all files
files <- list.files(here("data"))

#filter files
#this is needed because all four species were used as reference genomes in breseq
#thus there are coverage files for a species even if it doesn't appear in community
#the goal here is to thus filter all those uneneeded coverage files
filtered_files <- c()

for (i in 1:length(files)) {
  if (grepl(
    gsub(".*__(A|B|C|D)-.*","\\1",files[i]),
    gsub("^([ABCD]{1,4})[-_].*","\\1",files[i])
  )) {
    print(files[i])
    filtered_files[i] <- files[i]
  }
}

#remove NAs from filtered files
filtered_files <- filtered_files[!is.na(filtered_files)]

#set up empty dataframe as container
df <- data.frame(ID = NA, 
                 length_scaffold = NA, 
                 q0 = NA, 
                 q25 = NA, 
                 q50 = NA, 
                 q75 = NA, 
                 Average_Cov_TOP = NA, 
                 Average_Cov_BOT = NA, 
                 Average_Cov_Total = NA,
                 Sum_bp_TOP = NA,
                 Sum_bp_BOT = NA,
                 Sum_bp_Total = NA
                 ) 


#fill df file by file in loop
#load new coverage file, make calculations, remove, then load next file
for(i in 1:length(filtered_files)){
print(i)  
x <- read.csv(here(paste0("data/",filtered_files[i])), header = TRUE, sep = "")
x <- x[,1:2]
x$cov_per_bp <- rowSums(x[,1:2])


y <- data.frame(ID = filtered_files[i],
                 length_scaffold = nrow(x),
                 q0 = quantile(x$cov_per_bp)[[1]][1],
                 q25 = quantile(x$cov_per_bp)[[2]][1],
                 q50 = quantile(x$cov_per_bp)[[3]][1],
                 q75 = quantile(x$cov_per_bp)[[4]][1],
                 q100 = quantile(x$cov_per_bp)[[5]][1],
                 Average_Cov_TOP = mean(x$unique_top_cov),
                 Average_Cov_BOT = mean(x$unique_bot_cov),
                 Average_Cov_Total = mean(x$cov_per_bp),
                 Sum_bp_TOP = sum(x$unique_top_cov),
                 Sum_bp_BOT = sum(x$unique_bot_cov),
                 Sum_bp_Total = sum(x$cov_per_bp) 
                 )

df <- dplyr::bind_rows(df,y)

rm(x)

} 

#save_data
save(df, file = "total_bp_sequenced.RData")
















