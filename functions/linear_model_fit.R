### Linear model fit #####


# libraries ---------

library(agricolae)
library(modelr)
library(lmtest)


# linear model function to create a data frame with all models ----------

fit_linear_model <- function(response_variable, not_to_include = NULL, data) {

  # specify predictor for the models
  pred <- response_variable

  # save query dataframe
  query_df <- data %>% dplyr::select(-not_to_include)

  # specify response variables
  targets <- names(query_df)[names(query_df) != pred]

  # get all possible combinations - if user wants to test which combinations of variables is the best

  combinations <- NULL # generate empty vector
  get_combinations <- function(x) { # function to get all possible combinations without interaction
    c <- combn(targets, m =x , paste0, collapse = '+')
    combinations <- c(combinations, c)
  }

  combinations <- unlist(sapply(1:length(targets), get_combinations)) # create vector with all formulas

  # create a nested dataframe with all models and their information ---------

  model_summary <- expand.grid(pred, combinations) %>%    # create combinations
    mutate(model_id = row_number(),                       # create model id
           frml = paste0(Var1, "~", Var2)) %>%            #create model formula
    group_by(model_id, Var1, Var2) %>%                    # group by the above
    nest() %>%                                            # nest data
    mutate(model = map(data, ~lm(as.formula(.$frml), data = query_df)), # fit models
           summary = map(model, ~tidy(.)), # get summary for all models
           model_fit = map(model, ~glance(.)), # get information about the model fit
           normality_shapiro = map(model, ~shapiro.test(residuals(.))$p.value), # calculate normality p-value
           heter_breush = map(model, ~ as.numeric(bptest(.)$p.value)), # calculate heteroskedasticity p-value
           residuals = map(model, ~ augment(.)) # save residuals
    ) %>%
    unnest(normality_shapiro) %>% # unnest normality information
    mutate(model_fit = map(model_fit, ~ .x %>%
                             mutate(normal_pval = normality_shapiro, # add normality info to model_fit
                                    hetersk_pval = as.numeric(heter_breush)))) %>% # add homogeneity info to model_fit
    select(-normality_shapiro, -heter_breush) %>%
    ungroup()

  return(model_summary)
}

# Test using the iris data set

model_summary <- fit_linear_model(response_variable = "Sepal.Width", data = iris)

# Function to choose model ------

get_best_linear_model <- function(model_df) {
  best_model <- model_summary %>%
    unnest(model_fit) %>%
    arrange(desc(adj.r.squared), AIC, BIC, logLik) %>%
    slice_head()

  return(best_model)
}

# test function
best_model <- get_best_linear_model(model_summary)

# function to get model summary ---------

get_model_summary <- function(best_model) {
  best_model %>%
  dplyr::select(model_id, model, summary) %>%
  unnest(summary)
}

# Create diagnostic plots

LMdiagnostic_plots <- function(best_model) {
  # choose model
  model_choice <- best_model$model_id[1]

  # create plot data frame
  resid_df <- best_model %>%
    ungroup() %>%
    filter(model_id == model_choice) %>%
    select(residuals) %>%
    unnest(residuals)

  # create plot linearity assumption ---------

  p1 <- resid_df %>%
    ggplot(aes(.fitted, .resid)) +
    geom_point(color = "#FF6666", size = 2) +
    stat_smooth(method="loess", formula = 'y ~ x') +
    geom_hline(yintercept=0, col="red", linetype="dashed") +
    labs(x = "Fitted values", y = "Residuals", title = "Residual vs Fitted") +
    theme_light()

  # create plot for normality ---------

  grob <- grobTree(textGrob(paste("Shapiro-wilk test p-value:", round(model_fit_choice$normal_pval,4)), x=0.1,  y=0.95, hjust=0,
                            gp=gpar(col="red", fontsize=13, fontface="italic")))

  p2 <- resid_df %>%
    ggplot(aes(sample = .resid)) +
    geom_qq_line(size = 1.2) +
    stat_qq(color = "#FF6666", size = 2) +
    theme_light() +
    labs(title = "QQ-plot for normality assumption",
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    annotation_custom(grob)

  # create plot for Scale location (homogeinety) ------

  # create annotation with the Breusch-pagan test
  grob <- grobTree(textGrob(paste("Breusch-pagan test p-value:", round(model_fit_choice$hetersk_pval,4)),
                            x=0.1,  y=0.95, hjust=0,
                            gp=gpar(col="red", fontsize=13, fontface="italic")))


  p3 <- ggplot(resid_df, aes(.fitted, sqrt(abs(.std.resid)))) +
    geom_point(na.rm=TRUE,color = "#FF6666", size = 2) +
    stat_smooth(method="loess", na.rm = TRUE, formula = 'y ~ x') +
    labs(x = "Fitted Value" , y = expression(sqrt("Standardized residuals")),
         title = "Scale-Location Plot") +
    theme_light()  +
    annotation_custom(grob)

  # Cook's distance plot -------

  influence_tr <- function(p,n) {
    return(4/(n - p - 1))
  }

  p4 <- ggplot(resid_df, aes(seq_along(.cooksd), .cooksd)) +
    geom_bar(stat="identity", position="identity",fill = "#FF6666") +
    labs(x = "Observations" , y = "Cook's distance",
         title = "Cook's Distance Plot") +
    theme_light() +
    geom_hline(aes(yintercept = influence_tr(n_predictors, nrow(resid_df)), colour = "Influential treshold"),
               linetype = "dashed") +
    scale_color_manual(name = "Statistics",
                       values = c(`Influential threshold` = "purple")) +
    theme(legend.position = "bottom")


  # Residual vs leverage -----------

  # high leverage points will have larger values than the threshold
  lev_treshold <- function(p, n) {
    val <- (2*(p+1))/n
    return(val)
  }

  n_predictors <- ncol(resid_df %>% select(-c(.fitted:.std.resid)))


  p5 <- ggplot(resid_df, aes(x = .hat, y = .std.resid)) +
    geom_point(aes(size=.cooksd), na.rm=TRUE, colour = "#FF6666") +
    stat_smooth(method="loess", na.rm=TRUE, formula = 'y ~ x') +
    labs(x = "Leverage", y = "Standardized Residuals",
         title = "Residual vs Leverage Plot") +
    scale_size_continuous("Cook's Distance", range=c(1,5)) +
    theme_light() +
    theme(legend.position="bottom") +
    geom_vline(aes(xintercept = lev_treshold(n_predictors, nrow(resid_df)),
                   colour = "Leverage threshold"), linetype = "dashed") +
    geom_hline(aes(yintercept = 3,
                   colour = "Outlier threshold"), linetype = "dashed") +
    geom_hline(aes(yintercept = -3,
                   colour = "Outlier threshold"), linetype = "dashed") +
    scale_color_manual(name = "Statistics",
                       values = c(`Leverage threshold` = "purple", `Outlier threshold` = "orange"))



  # Cook's vs leverage -----------

  p6 <- ggplot(resid_df, aes(.hat, .cooksd)) +
    geom_point(na.rm=TRUE, colour = "#FF6666") +
    stat_smooth(method="loess", na.rm=TRUE, formula = 'y ~ x') +
    labs(x = "Leverage hii", y = "Cook's Distance",
         title = "Cook's distance vs Leverage Plot") +
    geom_abline(slope=seq(0,3,0.5), color="gray", linetype="dashed") +
    theme_light()

  # combine plots ----------

  p7 <- (p1+p2+p3) / (p4+p5+p6)

  ## save plots into list

  diag_plots <- list(fitted_resid_plot = p1, qq_plot = p2, scale_location_plot = p3,
                     cooks_plot = p4, resid_leverage_plot = p5, cooks_leverage_plot = p6,
                     combined_plot = p7)

}

# test function
diag_plot <- LMdiagnostic_plots(best_model)


# NEXT:

## FOR NUMERIC VALUES

# CREATE A SUMMARY PLOT FOR GETTING THE EFFECT OF EACH VARIABLE

## FOR GROUP VARIABLES

# CREATE A FUNCTION FOR LSD AND TUKEY GROUPS
# CREATE A FUNCTION FOR BAR PLOT WITH GROUPS (GET FROM MASTERS PAPER)


## INCLUDE INTERACTION TERMS INTO SEARCH FOR MODEL AND OPTION FOR USER TO INCLUDE HIS OWN FORMULA(S)

## CREATE A CORRECTION FORMULA. GIVE THE OPTION FOR CHOOSING TRANSFORMATION (LOG, BOX-COX)
library(moments)
skewness(iris$Sepal.Length, na.rm = TRUE) # check skewness to choose best transformation
#square-root for moderate skew:
# sqrt(x) for positively skewed data,
#sqrt(max(x+1) - x) for negatively skewed data

#log for greater skew:
#  log10(x) for positively skewed data,
#log10(max(x+1) - x) for negatively skewed data

#inverse for severe skew:
#  1/x for positively skewed data
#1/(max(x+1) - x) for negatively skewed data

#Linearity and heteroscedasticity:
#  first try log transformation in a situation where the dependent variable starts to increase more rapidly with increasing independent variable values
#If your data does the opposite – dependent variable values decrease more rapidly with increasing independent variable values – you can first consider a square transformation.

# auto transformation
# Use Lambert W x Gaussian transform.
# The R package LambertW has an implementation for automatically transforming heavy or light tailed data with Gaussianize().
#Tukey’s Ladder of Powers.
# For skewed data, the implementation transformTukey()from the R package rcompanion uses Shapiro-Wilk tests iteratively to find at which lambda value the data is closest to normality and transforms it. Left skewed data should be reflected to right skew and there should be no negative values.
# Box-Cox
# Yeo-Johnson Transformation.



## CREATE A ONE-STEP FUNCTION
