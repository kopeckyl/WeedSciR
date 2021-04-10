### function to generate a RCBD table with randomized plots


library(tidyverse)

# function examples

test <- RCBD_table(n_reps = 4, n_treats = 8,  outtable = TRUE)

colnames(test)

101 %% 100

plot <- test %>%
  gather(key = "Block", value = "Plot", -treatments) %>%
  mutate(x = 1,
         y = 100,
         Rep = Plot %% 100) %>%
  ggplot(aes(x = x, y = y, fill = treatments)) +
  geom_col() +
  facet_grid(Rep ~ Block) +
  theme(strip.text.y = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  geom_text(aes(label = Plot),
            position = position_stack(vjust = .5))

library(plotly)

# plot of 24m2
plot_h <- 12
plot_w <- 12

ggplotly(test %>%
           gather(key = "Block", value = "Plot", -treatments) %>%
           mutate(x = plot_w/2,
                  y = plot_h,
                  Rep = Plot %% 10) %>%
           ggplot(aes(x = x, y = y, fill = treatments)) +
           geom_col() +
           facet_grid(Rep ~ Block) +
           geom_text(aes(label = Plot),
                     position = position_stack(vjust = .5)) +
           labs(x = "", y = "") +
           theme_classic() +
           theme(legend.position = "none"),
         tooltip="treatments")

# function: specify the number of blocks and treatments. Also allow to save a tidy version of the table
RCBD_table <- function(n_reps, n_treats, outtable = TRUE, save_spread = FALSE, save_tidy = FALSE, blk_matrix = FALSE) {
  # variables
  #n_reps <- 6
  #n_treats <- 8
  # n_pop <- 14

  # get 100 levels blocks names
  blk_initial <- seq(from = 100,to = (n_reps*100), by = 100)

  # create a list of block names
  blk_names <- NULL
  for (i in seq(1:n_reps)){
    blk_names <- append(blk_names,rep(i,n_treats))
  }

  # create plot numbers
  blocks <- NULL
  block <- NULL
  for(i in blk_initial){
    for(j in seq(1:n_treats)){
      plot = i + j
      block <- append(block,plot)
    }
    block <- sample(block)
    blocks <- append(blocks, block)
    block <- NULL
  }

  # create tidy and spread tables
  blocks_df <- tibble(blocks, blk_names) %>%
    rename(plots = blocks) %>%
    mutate(blk_names = factor(paste0("Block_",blk_names)),
           treatments = rep(seq(1:n_treats),n_reps))
  blocks_spread <- blocks_df %>%
    spread(key = blk_names, value = plots)

  if (save_tidy == TRUE) {
    blocks_df <<- blocks_df
  }
  # create spread
  if (save_spread == TRUE){
    blocks_spread <<- blocks_df %>%
      spread(key = blk_names, value = plots)
  }
  # Spread table not randomized
  if (blk_matrix == TRUE){
    block_matrix <<- blocks_spread %>%
      select(-treatments) %>%
      gather(key = "block", value = "plot") %>%
      arrange(plot) %>%
      mutate(Row = rep(seq(1:n_treats),n_reps)) %>%
      spread(block, plot) %>%
      select(-Row)
  }

  #print table?
  if (outtable == TRUE){
    return(blocks_spread)
  }

}




