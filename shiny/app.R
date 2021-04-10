
#packages
library(shiny)
library(tidyverse)
library(DT)
library(kableExtra)
library(plotly)

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

  # create tidy tibble
  blocks_df <- tibble(blocks, blk_names) %>%
    rename(plots = blocks) %>%
    mutate(blk_names = factor(paste0("Block ",blk_names)),
           treatments = rep(seq(1:n_treats),n_reps))

  # create spread tibble
  blocks_spread <- blocks_df %>%
    spread(key = blk_names, value = plots) %>%
    mutate_at(vars(matches("Block")), as.integer) %>%
    mutate(treatments = factor(treatments))

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


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("WeedSci - RCBD design"),

    # Sidebar panel for inputs ----
    sidebarPanel(
    # set inputs - RCBD function
      sliderInput(
        inputId = "n_reps",
        label = "Number of blocks: ",
        value = 4,
        min = 2,
        max = 9
          ),
      sliderInput(
        inputId = "n_treats",
        label = "Number of treatments: ",
        value = 8,
        min = 2,
        max = 15
      ),
      sliderInput(
        inputId = "plot_dim",
        label = "Plot dimension (m2): ",
        value = 24,
        min = 2,
        max = 100
      )
    ),
    mainPanel(
      # add block table
      dataTableOutput(outputId = "RCBD_table"),
      # add block plot
      plotlyOutput(outputId = "block_plot")
    )
)


# Define server logic required to draw a histogram
server <- function(input, output) {

  # set experimental ggplot variable
  block_plot <- reactiveVal()

  # set data for table and plot
  exp_data <- reactive({
    data <- RCBD_table(n_reps = input$n_reps, n_treats = input$n_treats,  outtable = TRUE)
    block_plot(ggplotly(data %>%
                       gather(key = "Block", value = "Plot", -treatments) %>%
                       mutate(x = input$plot_dim/2,
                              y = input$plot_dim,
                              Rep = Plot %% 100) %>%
                       ggplot(aes(x = x, y = y, fill = treatments)) +
                       geom_col() +
                       facet_grid(Rep ~ Block) +
                       geom_text(aes(label = Plot),
                                 position = position_stack(vjust = .5)) +
                       labs(x = "", y = "") +
                       theme_classic() +
                       theme(legend.position = "none") +
                       theme(axis.title=element_blank(), axis.text=element_blank(),
                             axis.ticks=element_blank()),
                     tooltip="treatments"))
    data
  })
  # set output of table
  output$RCBD_table <- renderDataTable({
    exp_data()
  })
  # set output of plotly
  output$block_plot <- renderPlotly({
    block_plot()
  })
}

# Run the application
shinyApp(ui = ui, server = server)

