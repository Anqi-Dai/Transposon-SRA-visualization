# This script will plot a single clustered heatmap pair(positive and negative) using two input files, each has rows as bins and columns as samples.

library(magick)
library(tidyverse)
library(pheatmap)
library(gridExtra)

# for the color palette 
colGradient <- function( cols, length, cmax=255 )
{ ramp <- colorRamp(cols)
rgb( ramp(seq(0,1,length=length)), max=cmax )
}

# prepare the count table to draw
prepare_table_for_heatmap <- function(TE_name, StrandNum){
  DF <- Test
  
  all_table <- DF %>%
    filter(TE == TE_name) %>%
    dplyr::select(Sample, Status, bin, cnt = StrandNum)
  
  cnt_wide <-  all_table %>%
    dplyr::select(Sample, bin, cnt) %>%
    spread(key = Sample, value =  cnt ) %>%
    column_to_rownames('bin') %>%
    t
  return(cnt_wide)
}

draw_clustered_heatmap <- function(TE_name, Num){
  

  
  POS_table <- prepare_table_for_heatmap(TE_name, str_glue('Pos{Num}'))
  NEG_table <- prepare_table_for_heatmap(TE_name, str_glue('Neg{Num}'))
  
  max_cnt <- max(max(POS_table), max(NEG_table)) 
  
  # the color scale for the heatmap
  yellow <- colGradient(c("#EBEBF3","#4C3C90"),length=9) 
  blue <- colGradient(c("#DAE8F5","#0B559F"),length=9) 
  
  # the palette for the Status variable in the legend
  annotCol <- list(Status = c('#42B440', '#EC0000'))
  names(annotCol$Status) = c('Control','Huntingtons')
  
  # breaks : which set the limit of the count and color mapping relationship
  breaks <- seq(0, max_cnt, length.out = 1+length(blue))
  
  # the df linking the sample and its corresponding group
  annot <- data_frame(
    Sample = rownames(POS_table)
  ) %>%
    full_join(all_table %>%
                distinct(Sample, .keep_all = T) %>%
                dplyr::select(Sample, Status), by = 'Sample') %>%
    mutate(Status = factor(Status))%>%
    column_to_rownames('Sample')
  
  # draw heatmap for the positive strand
  pos_hm <- pheatmap(
    POS_table,
    color= blue, 
    breaks = breaks,annotation_row = annot,annotation_colors = annotCol,
    border_color = NA,treeheight_row = 100,angle_col = 45,
    show_rownames = TRUE,show_colnames = TRUE,
    height = 15,width = 20,
    cluster_rows = T,cluster_cols = F,
    fontsize = 20, fontsize_row=15,fontsize_col=10,
    main = if_else(Num == 1, str_glue('{TE_name} length 18-23nt POSITIVE STRAND clustered heatmap'), if_else(Num == 2, str_glue('{TE_name} length 24-35nt POSITIVE STRAND clustered heatmap'), str_glue('{TE_name} any length POSITIVE STRAND clustered heatmap')))
  )
  # for the negative strand
  neg_hm <- pheatmap(
    NEG_table,
    color= yellow, 
    breaks = breaks,annotation_row = annot,annotation_colors = annotCol,
    border_color = NA,treeheight_row = 100,angle_col = 45,
    show_rownames = TRUE,show_colnames = TRUE,
    height = 15,width = 20,
    cluster_rows = T,cluster_cols = F,
    fontsize = 20, fontsize_row=15,fontsize_col=10,
    main = if_else(Num == 1, str_glue('{TE_name} length 18-23nt NEGATIVE STRAND clustered heatmap'), if_else(Num == 2, str_glue('{TE_name} length 24-35nt NEGATIVE STRAND clustered heatmap'), str_glue('{TE_name} any length NEGATIVE STRAND clustered heatmap')))
  )
  
  ret <- ggarrange(pos_hm[[4]], neg_hm[[4]], 
                   ncol = 1, nrow = 2) 
  return(ret)
}

#*********************
#TE_name = 'ALU:1-312'
#Num = 1
#*********************
#ret <- draw_clustered_heatmap(TE_name, Num)



