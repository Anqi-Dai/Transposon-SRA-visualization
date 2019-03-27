## Differential line plots
# make this a plotting a function
library(tidyverse)
library(ggpubr)
library(MASS)
# HERV49I:1-6331

draw_differential_line_plot <- function(TE_name, Num) {

  # the annot table with the group info of the samples
  # the shaded area is the interquantile range
  
  DF <- Test %>%
    filter(TE == 'ALU:1-312' ) %>%
    dplyr::select( Status, bin, Pos1, Neg1)
  

  line_pos <- DF %>%
    dplyr::select(-Neg1,  cnt = Pos1) %>%
    group_by(Status, bin) %>%
    summarise(Median = median(cnt),
              quantile1 = quantile(cnt, 0.25),
              quantile3 = quantile(cnt, 0.75))

  
  line_neg <- negative %>%
    as.data.frame() %>%
    rename_all(
      funs(paste0('b',.))
    ) %>%
    rownames_to_column() %>%
    rename(Sample =rowname) %>%
    gather(key = bin , value = cnt, b0:b99) %>%
    left_join(sample_annot, by = 'Sample')  %>%
    group_by(Group, bin) %>%
    summarise(Median = median(cnt),
              quantile1 = quantile(cnt, 0.25),
              quantile3 = quantile(cnt, 0.75))
  
  # the negative binomial distribution part
  
  # it's like counts ~ group kind of glm regression
  # the pheno label
  #all.equal(rownames(annot) , colnames(pos))
  pheno <- annot$Group
  pos <- as.data.frame(pos)
  # for positive
  pos_reg <- apply(pos, 1, function(Row){
    ret = summary(glm.nb(as.numeric(Row) ~ pheno))[[12]][2,4]  
    return(ret)
  }) %>% 
    as.data.frame %>% 
    rename(pval = names(.)[1]) %>%
    mutate(padj = p.adjust(pval, method = 'BH'))


  # the position of the stars should be the max of the two medians and then + 20.
  star_position <- line_pos %>%
    group_by(bin) %>%
    summarise(star_height = max(Median) + 50)
  
  star_position_neg <- line_neg %>%
    group_by(bin) %>%
    summarise(star_height = max(Median) + 50)
  
  # the bins that should annotate with the star signs for the significance level
  pos_reg <- pos_reg %>%
    mutate(significance = if_else(padj < 0.01, '***', ifelse(padj < 0.05, '**', if_else(padj < 0.1, '*', '')))) %>%
    mutate(bin = paste0('b', seq(0,99)))%>%
    full_join(star_position , by = 'bin')
  
  neg_reg <- neg_reg %>%
    mutate(significance = if_else(padj < 0.01, '***', ifelse(padj < 0.05, '**', if_else(padj < 0.1, '*', '')))) %>%
    mutate(bin = paste0('b', seq(0,99))) %>%
    full_join(star_position_neg , by = 'bin')

  # add the significance to the line_pos table
  Dline_pos <- line_pos %>% 
    left_join(pos_reg %>%
                dplyr::select(significance, star_height, bin), by = 'bin') %>%
    mutate(bin = str_replace_all(bin, '^b','')) %>%
    mutate(bin = as.numeric(bin)) %>%
    arrange(bin)
  
  Dline_neg <- line_neg %>% 
    left_join(neg_reg %>%
                dplyr::select(significance, star_height, bin), by = 'bin')%>%
    mutate(bin = str_replace_all(bin, '^b','')) %>%
    mutate(bin = as.numeric(bin)) %>%
    arrange(bin)
  

  # make the differential line plot
  # for positive
  dl_pos <- ggplot(Dline_pos) + 
    geom_ribbon(aes(x=bin,
                    ymin=quantile1,
                    ymax=quantile3,
                    group=Group),alpha=0.40,fill="#3985ff") + 
    geom_line(aes(x=bin,
                  y=Median,
                  group=Group,
                  color=Group),lwd=1) +
    geom_point(aes(x=bin,
                   y=Median,
                   group=Group,
                   color=Group)) +
    geom_text(aes(x=bin, y=star_height, label = significance),angle = 90, nudge_x = 0.2) +
    ylim(0, max_cnt) + 
    labs(x = 'Bin',
         y = 'Positive',
         title = 'Differential lineplot') +
    theme_bw() +
    scale_x_continuous( breaks = seq(0,99)) +
    scale_color_manual(values = c('#00468B', '#EC0000'))  +
    theme(legend.position="none") +
    theme(axis.title.x=element_blank()) +
    theme(text = element_text(size=25)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
  
  # negative strand
  dl_neg <- ggplot(Dline_neg) + 
    geom_ribbon(aes(x=bin,
                    ymin=quantile1,
                    ymax=quantile3,
                    group=Group),alpha=0.40,fill="#3985ff") + 
    geom_line(aes(x=bin,
                  y=Median,
                  group=Group,
                  color=Group),lwd=1) +
    geom_point(aes(x=bin,
                   y=Median,
                   group=Group,
                   color=Group)) +
    geom_text(aes(x=bin, y=star_height, label = significance),angle = 90, nudge_x = 0.2) +
    ylim(0, max_cnt) + 
    labs(x = 'Bin',
         y = 'Negative') +
    theme_bw() +
    scale_color_manual(values = c('#00468B', '#EC0000'))  +
    scale_y_reverse()  + 
    ylim(max_cnt, 0) +
    scale_x_continuous( breaks = seq(0,99)) +
    theme(legend.position="bottom") +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          text = element_text(size=25),
          axis.text.x=element_blank()) 
  
  
  # now assemble the two plots together in one plot
  ggarrange(dl_pos, dl_neg, 
            ncol = 1, nrow = 2) +
    ggsave('../figs/assembled_lineplot.jpg', width = 20, height = 12, dpi = 300)
  
}