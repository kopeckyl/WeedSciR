### RCBD map plot #######


library(tidyverse)


RCBD_design <- function(blocks_n, treat_n, plot_label = c(treatments, plots), print_plan = TRUE){
  # # set number of replications and treatments
  #blocks_n <- 6
  # treat_n <- 10
  #generate treatments
  treatments <- c()
  for (t in 1:treat_n) {
    block <- sample(rep(1:treat_n))
    treatments <- append(treatments, block)
  }
  treatments <- factor(treatments)

  # generate plot data frame treatments
  x1 <- rep(1:treat_n, each = treat_n)
  x2 <- x1 + 1

  y1 <- rep(1:treat_n, treat_n)
  y2 <- y1 + 1


  df <- tibble(x1,x2,y1,y2,treatments)

  # generate plot number sequences
  seq_test <- seq(from = 100, to = 100*blocks_n, by = 100)
  plots <- NULL
  for (i in 1:treat_n){
    for (block in seq_test){
      plots <- append(plots, block+i)
    }
    plots <- sort(plots)
  }

  # create experimental tables
  experiment_table <- df %>%
    top_n(n= treat_n*blocks_n) %>%
    mutate(plots = plots,
           ID = seq(from = 1, to = treat_n*blocks_n, by = 1)) %>%
    select(ID,blocks = x1, plots, treatments)

  #save table to workspace
  experiment_plan <<- experiment_table

  # generate block sequences
  blocks <- seq(1:blocks_n)
  block_regions_x <- NULL
  block_regions_y <- NULL
  block_names <- NULL
  for (i in blocks) {
    block_regions_x <- append(block_regions_x,i+.5)
    block_regions_y <- append(block_regions_y, treat_n+1.5)
    block_names <- append(block_names, paste("Block ",i, sep = ""))
  }

  #get block annotation
  blocks_ann <- tibble(x1 = block_regions_x, y1 =block_regions_y, block_names)

  # generate plot for design
  plot_design <- df %>%
    mutate(plots = c(experiment_table$plots, rep(NA,nrow(df) - nrow(experiment_table)))) %>%
      full_join(blocks_ann) %>%
    ggplot(aes(x1,y1)) +
    xlim(1,blocks_n + 1) +
    ylim(1,treat_n+2) +
    geom_rect(mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = treatments), alpha = 0.2, color = "black") +
    geom_text(aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label= {{ plot_label }} ), size=4) +
    theme_nothing() +
    theme(legend.position = "none",
          axis.title = element_blank()) +
    geom_text(aes(x = x1, y = y1, label = block_names))

  # save plot to workspeace
  plot_design <<- plot_design
  # return table with treatments
  if (print_plan == TRUE) {
    return(plot_design)
  }
}

# run function
RCBD_design(blocks_n = 6, treat_n = 9, plot_label = plots, print_plan = TRUE)


