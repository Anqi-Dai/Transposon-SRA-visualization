#! /usr/bin/env Rscript
'Plot TPM distribution in the bins from the TE alignment data in different groups of samples.
Usage:
  TE_viz_line_plot.R <input_fn> <output_dir>

Options:
  -h --help  Show this screen.
  -v --version  Show version.

Arguments:
  input_fn  The absolute or relative path of the input file. Please see the READ.md for more details of the format of the file.
  output_dir  Output directory path either absolute or relative, for example "../figs"
' -> doc

# install if the package is missing
list.of.packages <- c("tidyverse", "docopt")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = "http://cran.us.r-project.org")

# load the library
suppressMessages(library(tidyverse))
suppressMessages(library(docopt))

arguments <- docopt(doc, version = 'E_viz_line_plot v1.0\n\n')

fn <- arguments$input_fn
outDir <- arguments$output_dir
if(!dir.exists(outDir)){dir.create(outDir)}

## Read the input data and format it
dat <- read_csv(fn, col_types = 'cccc') %>%
  mutate(fileContent = map(Path,  ~ suppressWarnings(read_tsv(., col_names = F, skip = 1, col_types = 'cddddddd')))) %>%
  unnest %>%
  dplyr::select(-X1, -Path) %>%
  rename(bin = X2,  
         Pos1 = X3,
         Neg1 = X4,
         Pos2 = X5,
         Neg2 = X6,
         Pos3 = X7,
         Neg3 = X8) 

## Function to draw
draw_line_plot <- function(TE_name, Num, statusGroup, outDir){
  # Input: 
  #   TE_name: the name of the TE 
  #   Num: 1 for length 18-23nt; 2 for length 24-35nt; 3 for any length.
  #   statusGroup: string of either 'Huntingtons' or 'Control'
  #   outDir: the output directory for saving the pdf result
  # Output:
  #   A line plot of ONE status with both pos(red) and neg(blue) lines on the same plot
  
  all_table <- dat %>%
    filter(Status == statusGroup) %>%
    filter(TE == TE_name) %>%
    dplyr::select(Sample, 
                  Status, 
                  bin, 
                  pos_cnt = str_glue('Pos{Num}'), 
                  neg_cnt = str_glue('Neg{Num}'))
  
  # the shaded area is the interquantile range
  # the line(point) is showing the median
  ribbon_part_pos <- all_table %>%
    group_by(bin) %>%
    summarise(Median = median(pos_cnt),
              quantile1 = quantile(pos_cnt, 0.25),
              quantile3 = quantile(pos_cnt, 0.75)) %>%
    mutate(Strand = 'Positive')
  
  ribbon_part_neg <- all_table %>%
    group_by(bin) %>%
    summarise(Median = median(neg_cnt),
              quantile1 = quantile(neg_cnt, 0.25),
              quantile3 = quantile(neg_cnt, 0.75)) %>%
    mutate(Strand = 'Negative')
  
  # the table for plotting
  final <- bind_rows(ribbon_part_pos, ribbon_part_neg) %>%
    mutate(Strand = factor(Strand))
  
  # plot
  g <- final %>%
    ggplot() + 
    geom_line(aes(x=bin,
                  y=Median,
                  group=Strand,
                  color=Strand),lwd=1) +
    geom_ribbon(aes(x=bin,
                    ymin=quantile1,
                    ymax=quantile3,
                    group=Strand, 
                    fill= Strand ),alpha=0.3) + 
    labs(title = if_else(Num == 1, str_glue('{TE_name} length 18-23nt line plot in {statusGroup}'), if_else(Num == 2, str_glue('{TE_name} length 24-35nt line plot in {statusGroup}'), str_glue('{TE_name} any length line plot in {statusGroup}')))) +
    theme_classic() +
    #scale_x_continuous(breaks = final$bin) +
    scale_color_manual(values = c('#00468B', '#EC0000'))  +
    scale_fill_manual(values = c('#00468B', '#EC0000'))  +
    theme(legend.position='top', 
          legend.justification='right',
          legend.direction='horizontal',
          text = element_text(size=25),
          axis.title.x=element_blank(),
          plot.margin = margin(30,20,30,20),
          axis.text.x = element_text(size = 20)) +
    ggsave(if_else(Num == 1, str_glue('{outDir}/{TE_name} length 18-23nt line plot in {statusGroup}.pdf'), if_else(Num == 2, str_glue('{outDir}/{TE_name} length 24-35nt line plot in {statusGroup}.pdf'), str_glue('{outDir}/{TE_name} any length line plot in {statusGroup}.pdf'))), 
           device = 'pdf', 
           width = 20, 
           height = 12, 
           dpi = 300)

}

## generate a dataframe with columns as inputs for the plotting function
TEs <- names(table(dat$TE))
Nums <- seq(1,3)
statusGroups <- names(table(dat$Status))

input_params <- expand.grid(TEs, Nums, statusGroups, outDir) %>%
  rename(TE_name = names(.)[1],
         Num = names(.)[2],
         statusGroup = names(.)[3],
         outDir = names(.)[4])

## main function to run
main <- function(){
  invisible(apply(input_params,
        1,
        function(row){
          TE_name<-row[["TE_name"]]
          Num<-row[["Num"]]
          statusGroup<-row[["statusGroup"]]
          outDir<-row[["outDir"]]
          draw_line_plot(TE_name,Num,statusGroup,outDir)
        }))
}


main()
