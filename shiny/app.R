# SHINY DASHBOARD  - WEEDsciR

# Packages library ------------------------------
library(shiny)
library(shinydashboard)
library(tidyverse)
library(DT)
library(kableExtra)
library(plotly)
library(rhandsontable) # user filled table
library(drc)
# Functions -----------------------------------


# RCBD function =========================
# Specify the number of blocks and treatments.
#Also allow to save a tidy version of the table
RCBD_table <- function(n_reps, n_treats, outtable = TRUE,
                       save_spread = FALSE, save_tidy = FALSE,
                       blk_matrix = FALSE) {
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


# Dose response analysis function ===================

# packages used
library(drc)
library(grid)
library(tidyverse)
library(broom)
library(patchwork)

# load data
df <- S.alba
# function
FitDoseResponse <- function(dose, var_name, group, df,
                            fit_weibull = FALSE, get_AIC_table = FALSE,
                            get_summary = TRUE) {

  # check if grouping variable is a factor
  if (is.factor(df[,group]) == FALSE) {
    df$group <- as.factor(df[,group])
  }

  # extract variables
  dose <- df[,as.character(dose)]
  var_name <- df[,as.character(var_name)]
  group <- df[,as.character(group)]

  # fit models
  if (fit_weibull == FALSE) {
    # fit log logistics models

    #builds a model with three parameters
    ll3_model <- drm(var_name ~ dose, group, data= df,
                     fct=LL.3())
    #builds a model with four parameters
    ll4_model <- drm(var_name ~ dose, group, data= df,
                     fct=LL.4())
    # builds a model with five parameters
    ll5_model <- drm(var_name ~ dose, group, data= df,
                     fct=LL.5())
    # join models into list
    models_list <- list(ll3_model = ll3_model,
                        ll4_model = ll4_model,
                        ll5_model = ll5_model)
  } else{
    # fit log logistics models
    # three parameters
    ll3_model <- drm(var_name ~ dose, group, data= df,
                     fct=LL.3())
    # four parameters
    ll4_model <- drm(var_name ~ dose, group, data= df,
                     fct=LL.4())
    # five parameters
    ll5_model <- drm(var_name ~ dose, group, data= df,
                     fct=LL.5())

    models_logistic <- list(ll3_model = ll3_model,
                            ll4_model = ll4_model,
                            ll5_model = ll5_model)
    # build model with weillbul curves
    #W1
    W13_model <- drm(var_name ~ dose, group, data = df,
                     fct = W1.3())
    W14_model <- drm(var_name ~ dose, group, data = df,
                     fct = W1.4())
    #W2
    W23_model <- drm(var_name ~ dose, group, data = df,
                     fct = W2.3())
    W24_model <- drm(var_name ~ dose, group, data = df,
                     fct = W2.4())
    # join weibull models into a list
    models_weibull <- list(W13_model = W13_model,
                           W14_model = W14_model,
                           W23_model = W23_model,
                           W24_model = W24_model)

    # create a list with all modules

    models_list <- c(models_logistic, models_weibull)
  }

  # get better model using AIC and BIC
  model_table <- models_list %>%
    lapply(FUN=glance) %>%
    lapply(FUN=function(x) x[(names(x) %in% c("AIC", "BIC"))]) %>%
    bind_rows(.id = "id") %>%
    mutate(model = models_list,
           summary = map(model, tidy)) %>%
    arrange(AIC, BIC)

  # select best model
  selected_model <<- get(model_table[1,]$id)



  # print message with the selected model
  print(paste("According to AIC and BIC parameters, the model", model_table[1,1]$id, "is the best fit for this data",
              "(AIC =", round(model_table[1,]$AIC,2), "BIC =", paste0(round(model_table[1,]$BIC,2), ")")))

  # if user wants to print the AIC table
  if (get_AIC_table == TRUE) {
    # create a AIC table
    AIC_table <<- model_table %>% dplyr::select(id:BIC)
    print("Model fit results saved into global environment.")
  }

  if (get_summary == TRUE) {
    summ_table <<- model_table[1,] %>%
      unnest(summary) %>%
      dplyr::select(paramater = term, group = curve, estimate:p.value) %>%
      dplyr::mutate(paramater = dplyr::case_when(paramater == "b" ~ "slope",
                                                 paramater == "c" ~ "lower",
                                                 paramater == "d" ~ "upper",
                                                 paramater == "e" ~ "ED50",
                                                 TRUE ~ as.character(paramater)))
    print("Summary table:")
    return(summ_table)

  }


}

# fit models and get tables
FitDoseResponse("Dose", "DryMatter", "Herbicide",
                df, fit_weibull = TRUE, get_AIC_table = TRUE)


# Dose response assumptions function ======================

# function to check for normality
normalQQ_plot <- function (model) # argument: vector of numbers
{
  # following four lines from base R's qqline()
  vec <- residuals(model)
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]

  # data frame for residuals
  d <- data.frame(resids = vec)

  # calculate shapiro test
  Shap_test <- shapiro.test(residuals(model))
  Shap_pval <- Shap_test$p.value

  # create annotation with the shapiro test
  grob <- grobTree(textGrob(paste("Shapiro-wilk test p-value:", round(Shap_pval,4)), x=0.1,  y=0.95, hjust=0,
                            gp=gpar(col="red", fontsize=13, fontface="italic")))

  # create plot
  ggplot(d, aes(sample = resids)) +
    geom_qq_line(size = 1.2) +
    stat_qq(color = "#FF6666", size = 2) +
    theme_light() +
    labs(title = "QQ-plot for normality assumption",
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    annotation_custom(grob)

}

# function to check for homogeineity
homogTest_plot <- function(model){
  # run fligner test - homogeinety
  df <- tibble(model$data[c(1,2,4)])
  Group <-  rep("Lower",nrow(df)) #Creates a vector that repeats "Lower" n times
  Group[selected_model$data$var_name > median(selected_model$data$var_name)] <-  "Upper" #Changing the appropriate values to "Upper"
  Group <- as.factor(Group) #Changes it to a factor, which R recognizes as a grouping variable.
  df$Group <- Group
  the.FKtest <- fligner.test(residuals(model), df$Group)
  FK_pval <- the.FKtest$p.value

  # create annotation with the fligner test
  grob <- grobTree(textGrob(paste("Fligner test p-value:", round(FK_pval,4)), x=0.1,  y=0.95, hjust=0,
                            gp=gpar(col="red", fontsize=13, fontface="italic")))

  #  create plot
  ggplot(augment(selected_model, data = df),
         aes(.fitted, .resid)) + geom_point(color = "#FF6666", size = 2) +
    stat_smooth(method="loess", formula = 'y ~ x') +
    geom_hline(yintercept=0, col="red", linetype="dashed") +
    labs(x = "Fitted values", y = "Residuals", title = "Residual vs Fitted") +
    theme_light() +
    annotation_custom(grob)
}

# function to apply correction if needed
applyBoxCox <- function(model)
{
  # calculate shapiro-wilk
  Shap_test <- shapiro.test(residuals(model))
  Shap_pval <- Shap_test$p.value

  # run fligner test - homogeneity
  df <- tibble(model$data[c(1,2,4)])
  Group <-  rep("Lower",nrow(df))
  Group[df$var_name > median(df$var_name)] <-  "Upper"
  df$Group <- as.factor(Group)
  the.FKtest <- fligner.test(residuals(model), df$Group)
  FK_pval <- the.FKtest$p.value

  # apply correction
  if (FK_pval <= 0.05 | Shap_pval <= 0.05) {
    print("Box-Cox correction applied using the Anova method. New model saved into global environment.")
    corrected_model <<- boxcox(selected_model,method="anova", plotit = F)
  } else {
    print("No correction needed.")
  }
}

# apply correction
applyBoxCox(selected_model)

# function for outlier removal
ouliersRemoval <- function(model) {

  # check for outliers
  df <- tibble(model$data[c(1,2,4)])
  ei.s <- residuals(model)/sqrt(sum(residuals(model)^2)/(nrow(df) - length(model$coefficients)))
  alpha <- 0.1 ; n = nrow(df); p = length(model$coefficients)
  cutoff <- qt(1-alpha/(2*n), n -p )
  cutoff.deleted <- qt(1-alpha/(2*n), n -p -1 )
  outliers <- which(abs(ei.s) > cutoff)

  # create new data without the outliers

  # return new data
  if (length(outliers) == 0) {
    print("No outliers detected.")
  } else{
    new.data <- df[-outliers,]
    return(new.data)
  }

}


# apply functions
ouliersRemoval(selected_model)
homogTest_plot(selected_model)
normalQQ_plot(selected_model)

# apply functions to corrected model
homogTest_plot(corrected_model)
normalQQ_plot(corrected_model)


# get summary from corrected plot

get_tidy_summary <- function(model){
  tidy(model) %>%
    dplyr::select(paramater = term, group = curve, estimate:p.value) %>%
    dplyr::mutate(paramater = dplyr::case_when(paramater == "b" ~ "slope",
                                               paramater == "c" ~ "lower",
                                               paramater == "d" ~ "upper",
                                               paramater == "e" ~ "ED50",
                                               TRUE ~ as.character(paramater)))
}


get_tidy_summary(corrected_model)

# Tabs content -------------------------------------------

# Experimental design tab =============================================
design_tab <- tabItem(tabName = "Design",
                      h2("Experimental design"),
                      fluidRow(
                        box(dataTableOutput(outputId = "RCBD_table")),
                        box(plotlyOutput(outputId = "block_plot")),

                        box(
                          sliderInput(
                            inputId = "n_reps",
                            label = "Number of blocks: ",
                            value = 4,
                            min = 2,
                            max = 9),
                          sliderInput(
                            inputId = "n_treats",
                            label = "Number of treatments: ",
                            value = 8,
                            min = 2,
                            max = 15),
                          sliderInput(
                            inputId = "plot_dim",
                            label = "Plot dimension (m2): ",
                            value = 24,
                            min = 2,
                            max = 100)
                        )
                      )
)



# Report tab =======================================
report_tab <- tabItem(
  tabName = "Report",
  h2("This tab will have report options")
)



# DR analysis tab ======================================
DR_analysis_tab <- tabItem(
  tabName = "DR_analysis",
  h2("This tab will have Dose response analysis options"),
  #fluidRow(column(4,rHandsontableOutput("test_table"), actionButton("saveBtn", "Save"))) # test table need to work on it

)

# linear analysis tab ======================================
linear_reg_tab <- tabItem(
  tabName = "linear_analysis",
  h2("This tab will have linear response analysis options"),
  #fluidRow(column(4,rHandsontableOutput("test_table"), actionButton("saveBtn", "Save"))) # test table need to work on it

)

# mixed effect analysis tab ======================================
mixed_reg_tab <- tabItem(
  tabName = "mixed_analysis",
  h2("This tab will have mixed effect response analysis options"),
  #fluidRow(column(4,rHandsontableOutput("test_table"), actionButton("saveBtn", "Save"))) # test table need to work on it
)


# mixed effect analysis tab ======================================
genetic_seg_tab <- tabItem(
  tabName = "segregation_analysis",
  h2("This tab will have genetic segregation analysis options"),
  #fluidRow(column(4,rHandsontableOutput("test_table"), actionButton("saveBtn", "Save"))) # test table need to work on it
)

# SideBar content ---------------------------------------------------------------------------

sideBar_content <- dashboardSidebar(
  sidebarMenu(
    menuItem("Experimental Design", tabName = "Design", icon = icon("calendar")),
    menuItem("Dose response analysis", tabName = "DR_analysis", icon = icon("signal", lib = "glyphicon")),
    menuItem("Linear regression analysis", tabName = "linear_analysis", icon = icon("signal", lib = "glyphicon")),
    menuItem("Mixed-effect model analysis", tabName = "mixed_analysis", icon = icon("signal", lib = "glyphicon")),
    menuItem("Gene segregation analysis", tabName = "segregation_analysis", icon = icon("signal", lib = "glyphicon")),
    menuItem("Generate report", tabName = "Report", icon = icon("save"))
  )
)

# BODY content ------------------------------------------------------------------------------

body_content <- dashboardBody(
  tabItems(
    # to make changes do it at the tab section above
    design_tab,
    report_tab,
    DR_analysis_tab,
    linear_reg_tab,
    mixed_reg_tab,
    genetic_seg_tab
  )
)

# UI ----------------------------------------------------------------------------------------

ui <- dashboardPage(
  dashboardHeader(title = "WeedSciR"),
  ## Sidebar content
  sideBar_content,
  ## Body content
  body_content,
  ## Aesthetic
  skin = "green"
)

# Server ------------------------------------------------------------------------------------

server <- function(input, output) {
  # EXPERIMENTAL DESIGN OUTPUT ########################
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
  output$RCBD_table <- renderDT({
    datatable(exp_data(),editable = "cell")
  })
  # set output of plotly
  output$block_plot <- renderPlotly({
    block_plot()
  })

  # DR OUTPUT ############

}

# Run shiny app ---------------------------------------------------------------------------
shinyApp(ui, server)

