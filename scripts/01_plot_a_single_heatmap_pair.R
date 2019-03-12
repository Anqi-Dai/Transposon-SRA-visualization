# This script will plot a single clustered heatmap pair(positive and negative) using two input files, each has rows as bins and columns as samples.
setwd('scripts/')
library(magick)
library(tidyverse)
library(pheatmap)

# load the fake data
pos <- read_csv('../data/test_data_pos.csv') %>%
  column_to_rownames(var = 'bin') %>%
  as.matrix()

neg <- read_csv('../data/test_data_neg.csv') %>%
  column_to_rownames(var = 'bin') %>%
  as.matrix()


# get the max number of the two matrix that will be the max value for the two heatmaps
max_cnt <- max(max(pos), max(neg)) 

colGradient <- function( cols, length, cmax=255 )
{ ramp <- colorRamp(cols)
rgb( ramp(seq(0,1,length=length)), max=cmax )
}

yellow <- colGradient(c("#EBEBF3","#4C3C90"),length=15) 

blue <- colGradient(c("#DAE8F5","#0B559F"),length=15) 

# Create a heatmap of your genes using the entire training set.
positive <- t(pos)
negative <- t(neg)

sample_names <- rownames(positive)

annot <- data_frame(
  Sample =sample_names, 
  Group = c(rep('Disease', length(sample_names[grepl(pattern = '^D', sample_names)])),
            rep('Control', length(sample_names[grepl(pattern = '^C', sample_names)])) )
) %>%
  mutate(Group = factor(Group))%>%
  column_to_rownames('Sample')

# the palette for the Group variable
annotCol <- list(Group = c('#42B440', '#EC0000'))
names(annotCol$Group) = c('Control','Disease')

# breaks : which set the limit of the count and color mapping relationship
breaks <- seq(0, max_cnt, length.out = 1+length(blue))


# for the positive strand
pheatmap(
  positive,
  color= blue, 
  breaks = breaks,
  annotation_row = annot,
  annotation_colors = annotCol,
  border_color = NA,
  treeheight_row = 100,
  angle_col = 45,
  show_rownames = TRUE,
  show_colnames = TRUE,
  height = 5,
  width = 25,
  cluster_rows = T,
  cluster_cols = F,
  filename = '../figs/pos_heatmap.png',
  main = 'Positive strand'
)


# for the negative strand
pheatmap(
  negative,
  color= yellow, 
  breaks = breaks,
  annotation_row = annot,
  annotation_colors = annotCol,
  border_color = NA,
  treeheight_row = 100,
  angle_col = 45,
  show_rownames = TRUE,
  show_colnames = TRUE,
  height = 5,
  width = 25,
  cluster_rows = T,
  cluster_cols = F,
  filename = '../figs/neg_heatmap.png',
  main = 'Negative strand'
)


# combine the above two and make them one figure

positive_hm <- image_read('../figs/pos_heatmap.png')
negative_hm <- image_read('../figs/neg_heatmap.png')
imgs = c(positive_hm, negative_hm)
# concatenate them left-to-right (use 'stack=T' to do it top-to-bottom)
heatmaps = image_append(imgs, stack=T)
image_write(heatmaps, path = "../figs/heatmaps.jpg", format = "jpg")
```