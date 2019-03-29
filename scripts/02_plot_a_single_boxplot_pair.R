## Per-bin boxplots
# make this a plotting a function
library(tidyverse)
library(ggpubr)

draw_per_bin_boxpot <- function(TE_name, Num) {
  # the input should be 
  box_pos <- DF %>%
    filter(TE == TE_name) %>%
    dplyr::select(bin, cnt = str_glue('Pos{Num}'))
  
  # change all the negative numbers on the neg strand to its absolute value.
  box_neg <- DF %>%
    filter(TE == TE_name) %>%
    dplyr::select(bin, cnt = str_glue('Neg{Num}')) %>%
    mutate(cnt = abs(cnt))
  
  # get the max number of the two matrix that will be the max value for the two heatmaps
  max_cnt <- max(max(box_pos$cnt), max(box_neg$cnt)) 
  
  
  # plot the TPM of each bin as a box plot. the y limit should be 0 to max count too
  # the postitive strand
  pos_box <- ggboxplot(box_pos, x = 'bin', y = 'cnt',
                       color = '#00468B',
                       xlab = 'Bin', ylab = 'Positive', 
                       title = if_else(Num == 1, str_glue('{TE_name} length 18-23nt per-bin boxplots'), if_else(Num == 2, str_glue('{TE_name} length 24-35nt per-bin boxplots'), str_glue('{TE_name} any length per-bin boxplots')))) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          text = element_text(size=25),
          plot.margin = margin(20,20,20,20),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
  
  # the negative strand
  neg_box <- ggboxplot(box_neg, x = 'bin', y = 'cnt',
                       color = '#00468B',
                       xlab = 'Bin', ylab = 'Negative' )  +
    scale_y_reverse()  + 
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          plot.margin = margin(20,20,20,20),
          text = element_text(size=25),
          axis.ticks.x=element_blank()) 
  
  # now assemble the two plots together in one plot
  ret <- ggarrange(pos_box, neg_box, 
            ncol = 1, nrow = 2) 
  return(ret)
}



