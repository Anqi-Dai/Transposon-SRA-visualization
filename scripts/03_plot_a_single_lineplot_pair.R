## Differential line plots
# make this a plotting a function
library(tidyverse)
library(ggpubr)
library(MASS)

# the function to generate the table for plotting since the pos and neg have many things in common
prepare_table_for_plotting <- function(TE_name, StrandNum) {
  DF <- Test
  
  all_table <- DF %>%
    filter(TE == TE_name) %>%
    dplyr::select(Sample, Status, bin, cnt = StrandNum)
  
  # the shaded area is the interquantile range
  ribbon_part <- all_table %>%
    dplyr::select(-Sample)  %>%
    group_by(Status, bin) %>%
    summarise(Median = median(cnt),
              quantile1 = quantile(cnt, 0.25),
              quantile3 = quantile(cnt, 0.75))
  
  # the position of the stars should be the max of the two medians and then + 10% of the median of each medians.
  defined_star_height <- ribbon_part %>%
    group_by(bin) %>%
    summarise(star_height = max(Median) ) %>%
    summarise(star_height = median(star_height)*0.1) %>%
    pull(star_height)
  
  star_position_part <-  ribbon_part %>%
    group_by(bin) %>%
    summarise(star_height = max(Median) + defined_star_height )
  
  # the regression part
  # it's like counts ~ group kind of glm regression
  pheno <- all_table %>%
    arrange(Sample) %>%
    distinct(Sample, .keep_all = T) %>%
    pull(Status)
  
  cnt_wide <-  all_table %>%
    dplyr::select(Sample, bin, cnt) %>%
    spread(key = Sample, value =  cnt ) %>%
    column_to_rownames('bin')
  
  regression_part <- apply(cnt_wide, 1, function(Row){
    fit <- try(suppressWarnings(glm.nb(as.numeric(Row) ~ pheno)))
    if(is(fit, "try-error")) return(1) else return(summary(fit)[[12]][2,4]) 
  }) %>% 
    as.data.frame %>% 
    mutate(bin = as.numeric(rownames(cnt_wide))) %>%
    rename(pval = names(.)[1]) %>%
    mutate(padj = p.adjust(pval, method = 'BH')) %>%
    mutate(significance = if_else(padj < 0.01, '***', ifelse(padj < 0.05, '**', if_else(padj < 0.1, '*', '')))) %>%
    full_join(star_position_part , by = 'bin')
  
  # add the regression part to the ribbon part
  final <- ribbon_part %>% 
    left_join(regression_part %>%
                dplyr::select(significance, star_height, bin), by = 'bin') 
  return(final)
}


# the actual function to draw the plot
draw_differential_line_plot <- function(TE_name, Num) {
  
  POS_table <- prepare_table_for_plotting(TE_name, str_glue('Pos{Num}'))
  NEG_table <- prepare_table_for_plotting(TE_name, str_glue('Neg{Num}'))

  # the max_cnt, the scale of the plot should be the max of the two quantile 3 times 1.2
  max_q3 <-max(max(POS_table$quantile3), max(NEG_table$quantile3))
  max_cnt <-  max_q3 * 1.2

  # make the differential line plot
  # for positive
  dl_pos <- ggplot(POS_table) + 
    geom_ribbon(aes(x=bin,
                    ymin=quantile1,
                    ymax=quantile3,
                    group=Status, 
                    fill= Status ),alpha=0.3) + 
    geom_line(aes(x=bin,
                  y=Median,
                  group=Status,
                  color=Status),lwd=1) +
    geom_point(aes(x=bin,
                   y=Median,
                   group=Status,
                   color=Status)) +
    geom_text(aes(x=bin, y=star_height, label = significance),angle = 90, nudge_x = 0.2) +
    ylim(0, max_cnt) + 
    labs(x = 'Bin',
         y = 'Positive',
         title = if_else(Num == 1, str_glue('{TE_name} length 18-23nt differential lineplot'), if_else(Num == 2, str_glue('{TE_name} length 24-35nt differential lineplot'), str_glue('{TE_name} any length differential lineplots')))) +
    theme_bw() +
    scale_x_continuous( breaks = POS_table$bin) +
    scale_color_manual(values = c('#00468B', '#EC0000'))  +
    scale_fill_manual(values = c('#00468B', '#EC0000'))  +
    theme(legend.position='top', 
          legend.justification='right',
          legend.direction='horizontal') +
    theme(axis.title.x=element_blank()) +
    theme(text = element_text(size=25),
          plot.margin = margin(20,20,20,20)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
  # negative strand
  dl_neg <- ggplot(NEG_table) + 
    geom_ribbon(aes(x=bin,
                    ymin=quantile1,
                    ymax=quantile3,
                    group=Status, 
                    fill= Status ),alpha=0.3) + 
    geom_line(aes(x=bin,
                  y=Median,
                  group=Status,
                  color=Status),lwd=1) +
    geom_point(aes(x=bin,
                   y=Median,
                   group=Status,
                   color=Status)) +
    geom_text(aes(x=bin, y=star_height, label = significance),angle = 90, nudge_x = 0.2) +
    labs(x = 'Bin',
         y = 'Negative') +
    theme_bw() +
    scale_color_manual(values = c('#00468B', '#EC0000'))  +
    scale_fill_manual(values = c('#00468B', '#EC0000'))  +
    #scale_y_reverse()  + 
    ylim(max_cnt, 0) +
    scale_x_continuous( breaks = NEG_table$bin)  +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          text = element_text(size=25),
          legend.position='none', 
          plot.margin = margin(20,20,20,20),
          axis.text.x=element_blank())  
   
  
  
  # now assemble the two plots together in one plot
  ret <- ggarrange(dl_pos, dl_neg, 
            ncol = 1, nrow = 2) 
  return(ret)
  }

#prepare_table_for_plotting('HERV49I:1-6331', 'Pos1')

#draw_differential_line_plot('HERV49I:1-6331', 1) +ggsave('figs/test.line.jpg',  width = 20, height = 12, dpi = 300)
