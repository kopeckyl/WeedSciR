### Experimental design map plot #######


library(tidyverse)
library(ggmap)
# RCBD design function ------------------------
RCBD_design <- function(blocks_n, treat_n, plot_label = c(treatments, plots), seed = NA){
  # # set number of replications and treatments
  #blocks_n <- 6
  # treat_n <- 10
  #generate treatments

  # Save the old random seed and use the new one, if present
  if (!is.na(seed)) {
    if (exists(".Random.seed"))  { saved.seed <- .Random.seed }
    else                         { saved.seed <- NA }
    set.seed(seed)
  }


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
    dplyr::select(ID,blocks = x1, plots, treatments)

  #save table
  experiment_plan <- experiment_table

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
      full_join(blocks_ann, by = c("x1", "y1")) %>%
    ggplot(aes(x1,y1)) +
    xlim(1,blocks_n + 1) +
    ylim(1,treat_n+2) +
    geom_rect(mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = treatments), alpha = 0.2, color = "black",
              na.rm = T) +
    geom_text(aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label= {{ plot_label }} ), size=4,
              na.rm = T) +
    theme_nothing() +
    theme(legend.position = "none",
          axis.title = element_blank()) +
    geom_text(aes(x = x1, y = y1, label = block_names), na.rm = TRUE)

  # save plot
  plot_design <- plot_design
  # return table with treatments
  # if (print_plan == TRUE) {
  #   return(plot_design)
  # }

  # get final results
  experimental_design <- list(experiment_book = experiment_plan, block_ann = blocks_ann,
                              plot = plot_design)

  return(experimental_design)
  # Restore the old random seed, if present
  if (!is.na(seed) && !is.na(saved.seed)) { .Random.seed <- saved.seed }

  if (returnstrings) { return(squareid) }
  else               { return(allsq) }
}

# run function
test <- RCBD_design(blocks_n = 4, seed = 1234, treat_n = 9, plot_label = plots)
test




# Generate Latin square function ------------
## - len is the size of the latin square
## - reps is the number of repetitions - how many Latin squares to generate
## - seed is a random seed that can be used to generate repeatable sequences
## - returnstrings tells it to return a vector of char strings for each square,
##    instead of a big matrix. This option is only really used for checking the
##    randomness of the squares.

treat <- c("A","B","C","D")
length(treatments)

latinsquare <- function(treatments, reps=1, seed=NA, returnstrings=FALSE) {

  # Save the old random seed and use the new one, if present
  if (!is.na(seed)) {
    if (exists(".Random.seed"))  { saved.seed <- .Random.seed }
    else                         { saved.seed <- NA }
    set.seed(seed)
  }

  # get the length of treatments
  len <- length(treatments)

  # This matrix will contain all the individual squares
  allsq <- matrix(nrow=reps*len, ncol=len)

  # Store a string id of each square if requested
  if (returnstrings) {  squareid <- vector(mode = "character", length = reps) }

  # Get a random element from a vector (the built-in sample function annoyingly
  #   has different behavior if there's only one element in x)
  sample1 <- function(x) {
    if (length(x)==1) { return(x) }
    else              { return(sample(x,1)) }
  }

  # Generate each of n individual squares
  for (n in 1:reps) {

    # Generate an empty square
    sq <- matrix(nrow=len, ncol=len)

    # If we fill the square sequentially from top left, some latin squares
    # are more probable than others.  So we have to do it random order,
    # all over the square.
    # The rough procedure is:
    # - randomly select a cell that is currently NA (call it the target cell)
    # - find all the NA cells sharing the same row or column as the target
    # - fill the target cell
    # - fill the other cells sharing the row/col
    # - If it ever is impossible to fill a cell because all the numbers
    #    are already used, then quit and start over with a new square.
    # In short, it picks a random empty cell, fills it, then fills in the
    # other empty cells in the "cross" in random order. If we went totally randomly
    # (without the cross), the failure rate is much higher.
    while (any(is.na(sq))) {

      # Pick a random cell which is currently NA
      k <- sample1(which(is.na(sq)))

      i <- (k-1) %% len +1       # Get the row num
      j <- floor((k-1) / len) +1 # Get the col num

      # Find the other NA cells in the "cross" centered at i,j
      sqrow <- sq[i,]
      sqcol <- sq[,j]

      # A matrix of coordinates of all the NA cells in the cross
      openCell <-rbind( cbind(which(is.na(sqcol)), j),
                        cbind(i, which(is.na(sqrow))))
      # Randomize fill order
      openCell <- openCell[sample(nrow(openCell)),]

      # Put center cell at top of list, so that it gets filled first
      openCell <- rbind(c(i,j), openCell)
      # There will now be three entries for the center cell, so remove duplicated entries
      # Need to make sure it's a matrix -- otherwise, if there's just
      # one row, it turns into a vector, which causes problems
      openCell <- matrix(openCell[!duplicated(openCell),], ncol=2)

      # Fill in the center of the cross, then the other open spaces in the cross
      for (c in 1:nrow(openCell)) {
        # The current cell to fill
        ci <- openCell[c,1]
        cj <- openCell[c,2]
        # Get the numbers that are unused in the "cross" centered on i,j
        freeNum <- which(!(1:len %in% c(sq[ci,], sq[,cj])))

        # Fill in this location on the square
        if (length(freeNum)>0) { sq[ci,cj] <- sample1(freeNum) }
        else  {
          # Failed attempt - no available numbers
          # Re-generate empty square
          sq <- matrix(nrow=len, ncol=len)

          # Break out of loop
          break;
        }
      }
    }

    # Store the individual square into the matrix containing all squares
    allsqrows <- ((n-1)*len) + 1:len
    allsq[allsqrows,] <- sq

    # Store a string representation of the square if requested. Each unique
    # square has a unique string.
    if (returnstrings) { squareid[n] <- paste(sq, collapse="") }

  }

  # Restore the old random seed, if present
  if (!is.na(seed) && !is.na(saved.seed)) { .Random.seed <- saved.seed }

  if (returnstrings) { return(squareid) }
  else               { return(allsq) }}

latinsquare(treatments = treat, reps = 2)[1,2]

latinSquare_design <- function(treatments, reps=1, seed=NA, returnstrings=F) {
  square <- latinsquare(treatments, reps, seed, returnstrings)

  if (returnstrings == TRUE) {
    return(square)
  } else {
    new_square <- matrix(nrow = nrow(square), ncol = ncol(square))
    for (row in 1:nrow(square)) {
      for (col in 1:ncol(square)) {
        value <- square[row,col]
        new_square[row,col] <- treatments[value]
      }
    }
    new_square <- as_tibble(new_square)
    names(new_square) <- gsub(x = names(new_square), pattern = "V", replacement = "Col_")
    new_square <- new_square %>% mutate(Row = rep(1:4, reps)) %>% dplyr::select(Row, everything())
    if (reps > 1) {
      new_square <- new_square %>% mutate(reps = rep(1:reps, each = length(treatments)))
    }
    return(new_square)
  }
}


test <- latinSquare_design(treatments = treat, reps = 2, seed = 1234)

reps = 2
test[,1:length(treat)+1]

# generate plot data frame treatments
x1 <- rep(1:length(treat), each = length(treat))
x2 <- x1 + 1

y1 <- rep(1:nrow(test), nrow(test))
y2 <- y1 + 1
df <- tibble(x1,x2,y1,y2,treatments)


# generate block sequences
blocks <- seq(1:length(treat))
block_regions_x <- NULL
block_regions_y <- NULL
block_names <- NULL
for (i in blocks) {
  block_regions_x <- append(block_regions_x,i+.5)
  block_regions_y <- append(block_regions_y, length(treat)+1.5)
  block_names <- append(block_names, paste("Column ",i, sep = ""))
}

#get block annotation
blocks_ann <- tibble(x1 = block_regions_x, y1 =block_regions_y, block_names)

# generate plot for design
plot_design <- df %>%
  mutate(plots = c(experiment_table$plots, rep(NA,nrow(df) - nrow(experiment_table)))) %>%
  full_join(blocks_ann, by = c("x1", "y1")) %>%
  ggplot(aes(x1,y1)) +
  xlim(1,blocks_n + 1) +
  ylim(1,treat_n+2) +
  geom_rect(mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill = treatments), alpha = 0.2, color = "black",
            na.rm = T) +
  geom_text(aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label= {{ plot_label }} ), size=4,
            na.rm = T) +
  theme_nothing() +
  theme(legend.position = "none",
        axis.title = element_blank()) +
  geom_text(aes(x = x1, y = y1, label = block_names), na.rm = TRUE)

