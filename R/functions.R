##### FUNCTIONS FOR PACKAGE #####

### LIBRARY ---------
library(devtools)
library(drc)
library(grid)
library(tidyverse)
library(broom)
library(patchwork)


# FUNCTIONS --------------------------------

# RCBD FUNCTION ======
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

# FIT DOSE RESPONSE FUNCTION =====

FitDoseResponse <- function(dose, var_name, variable, df,
                            fit_weibull = FALSE, get_AIC_table = FALSE,
                            get_summary = FALSE, applyBC = FALSE, remove_outliers = FALSE,
                            plot_curve = FALSE) {

  # check if grouping variable is a factor
  if (is.factor(df[,variable]) == FALSE) {
    df$variable <- as.factor(df[,variable])
  }

  # extract variables
  dose <- df[,as.character(dose)]
  var_name <- df[,as.character(var_name)]
  variable <- df[,as.character(variable)]

  # fit models
  if (fit_weibull == FALSE) {
    # fit log logistics models

    #builds a model with three parameters
    ll3_model <- do.call("drm", list(var_name ~ dose, variable, data= df,
                                     fct=LL.3()))
    #builds a model with four parameters
    ll4_model <- do.call("drm", list(var_name ~ dose, variable, data= df,
                                     fct=LL.4()))
    # builds a model with five parameters
    ll5_model <- do.call("drm", list(var_name ~ dose, variable, data= df,
                                     fct=LL.5()))
    # join models into list
    models_list <- list(ll3_model = ll3_model,
                        ll4_model = ll4_model,
                        ll5_model = ll5_model)
  } else{
    # fit log logistics models
    #builds a model with three parameters
    ll3_model <- do.call("drm", list(var_name ~ dose, variable, data= df,
                                     fct=LL.3()))
    #builds a model with four parameters
    ll4_model <- do.call("drm", list(var_name ~ dose, variable, data= df,
                                     fct=LL.4()))
    # builds a model with five parameters
    ll5_model <- do.call("drm", list(var_name ~ dose, variable, data= df,
                                     fct=LL.5()))

    models_logistic <- list(ll3_model = ll3_model,
                            ll4_model = ll4_model,
                            ll5_model = ll5_model)
    # build model with weillbul curves
    #W1
    W13_model <- do.call("drm", list(var_name ~ dose, variable, data= df,
                                     fct = W1.3()))
    W14_model <- do.call("drm", list(var_name ~ dose, variable, data= df,
                                     fct = W1.4()))
    #W2
    W23_model <- do.call("drm", list(var_name ~ dose, variable, data= df,
                                     fct = W2.3()))
    W24_model <- do.call("drm", list(var_name ~ dose, variable, data= df,
                                     fct = W2.4()))
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

  # remove outliers

  if (remove_outliers == TRUE) {
    print("Check for outliers...")
    ouliersRemoval(selected_model)
  }

  # boxCox transformation

  if (applyBC == TRUE) {
    print("Applying Box-Cox correction...")
    applyBoxCox(selected_model)
  }

  # if user wants to print the AIC table
  if (get_AIC_table == TRUE) {
    # create a AIC table
    AIC_table <<- model_table %>% dplyr::select(id:BIC)
    print("Model fit results saved as AIC_table .")
  }

  if (get_summary == TRUE) {
    summ_table <<- model_table[1,] %>%
      unnest(summary) %>%
      dplyr::select(paramater = term, variable = curve, estimate:p.value) %>%
      dplyr::mutate(paramater = dplyr::case_when(paramater == "b" ~ "slope",
                                                 paramater == "c" ~ "lower",
                                                 paramater == "d" ~ "upper",
                                                 paramater == "e" ~ "ED50",
                                                 TRUE ~ as.character(paramater)))
    print("Summary table saved as summ_table.")

  }

  # plot curve
  if (plot_curve == TRUE) {
    plot_DR(selected_model)
  }

}

# function to check for normality =====
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

# function to check for homogeinety =====
homogTest_plot <- function(dr_model){
  # run fligner test - homogeneity
  data_model <- data.frame(dr_model$data[c(1,2,4)])
  direction <- rep("Lower",nrow(data_model))
  direction[data_model$var_name > median(data_model$var_name)] <- "Upper"
  data_model$direction <- as.factor(direction)
  the.FKtest <- fligner.test(residuals(dr_model), data_model$direction)
  FK_pval <- the.FKtest$p.value

  # create annotation with the fligner test
  grob <- grobTree(textGrob(paste("Fligner test p-value:", round(FK_pval,4)), x=0.1,  y=0.95, hjust=0,
                            gp=gpar(col="red", fontsize=13, fontface="italic")))

  #  create plot
  ggplot(augment(selected_model, data = data_model),
         aes(.fitted, .resid)) + geom_point(color = "#FF6666", size = 2) +
    stat_smooth(method="loess", formula = 'y ~ x') +
    geom_hline(yintercept=0, col="red", linetype="dashed") +
    labs(x = "Fitted values", y = "Residuals", title = "Residual vs Fitted") +
    theme_light() +
    annotation_custom(grob)
}

# function to look for outliers ======
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

dr_model <- selected_model
# function to apply correction (not working properly) =====
applyBoxCox <- function(dr_model,data){
  # calculate shapiro-wilk
  Shap_test <- shapiro.test(residuals(dr_model))
  Shap_pval <- Shap_test$p.value

  # run fligner test - homogeneity
  data_model <- data.frame(dr_model$data[c(1,2,4)])
  direction <- rep("Lower",nrow(data_model))
  direction[data_model$var_name > median(data_model$var_name)] <- "Upper"
  data_model$direction <- as.factor(direction)
  the.FKtest <- fligner.test(residuals(dr_model), data_model$direction)
  FK_pval <- the.FKtest$p.value

  # apply correction
  if (FK_pval <= 0.05 | Shap_pval <= 0.05) {
    print("Box-Cox correction applied using the Anova method. New model saved into global environment.")
    corrected_model <<- boxcox(dr_model,method="anova", plotit = F)
  } else {
    print("No correction needed.")
  }
}

# function to call the first type of plot =====

plot_DR <- function(DR_model) {

  # general sd function
  std_mean <- function(x) sd(x,na.rm=TRUE)/sqrt(length(x))

  # correct model dataframe and create a dataframe to plot
  DR_data <- DR_model$data
  cols <- which(names(DR_data) == 'variable')
  names(DR_data)[cols] <- paste0('variable', seq_along(cols))
  df_plot <- tibble(DR_data) %>%
    group_by(dose, variable2) %>%
    summarize(var_value = mean(var_name, na.rm=TRUE),
              sd = std_mean(var_name)) %>%
    mutate(dose = dose + 0.1)

  # generate plot
  p1 <- ggplot(data = df_plot, aes(x = dose, y = var_value)) +
    geom_point(aes(color = variable2,
                   text = paste("Dose:", dose,
                                "\nHerbicide:", variable2,
                                "\nBiomass:", var_value))) +
    scale_x_log10() +
    geom_errorbar(mapping=aes(ymin=var_value-sd, ymax=var_value+sd,color = variable2), width=0.2, alpha = .4) +
    geom_smooth(aes(color = variable2),
                method = drm,
                method.args = list(fct = L.4()), se = F) +
    theme_light() +
    labs(title= "", x = "Dose (g a.i /ha)",  y = "Biomass") +
    theme(legend.position = "bottom") + guides(color=guide_legend(title="Group"))

  p1
}

