library(tidyverse)

if(Sys.info()["nodename"]=="DDJ9Y7V9GC"){
  setwd("/Volumes/TracerX/")
}else{
  setwd("/camp/project/tracerX/")
}


df <- read_tsv("working/VCAM1_GnT/DATA/plate9/bam/overlap_metrics.tsv", show_col_types = FALSE)

# Proportion of wells that show stacking of 5+ reads anywhere
prop_wells_with_gt3 <- mean(df$prop_gt3 > 0)
prop_wells_with_gt3

# If you want a stricter definition, for example more than 1 percent of reads in 5+ stacks
prop_wells_with_gt3_1pct <- mean(df$prop_gt3 > 0.01)
prop_wells_with_gt3_1pct

ggplot(df, aes(x = prop_gt3)) +
  geom_histogram(bins = 50) +
  labs(x = "Fraction of alignments in stacks of 5 or more (overlap > 3)",
       y = "Number of wells")

ggplot(df %>% arrange(desc(prop_gt3)),
       aes(x = reorder(well, prop_gt3), y = prop_gt3)) +
  geom_col() +
  coord_flip() +
  theme_bw() + 
  labs(x = "Well", y = "Fraction in stacks of 5 or more")
