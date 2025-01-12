# load libraries
library(drc)
library(grid)
library(tidyverse)
library(broom)
library(patchwork)

# load data
df <- S.alba
# choose variables
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

FitDoseResponse("Dose", "DryMatter", "Herbicide",
                df, fit_weibull = TRUE, get_AIC_table = TRUE)

# plot model assumption plots

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

# function for variance assumption
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
applyBoxCox <- function(model){
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

applyBoxCox()
selected_model$call
# final diagnostic plot
homogTest_plot(selected_model)
normalQQ_plot(selected_model)
#apply correction
applyBoxCox(selected_model)


# function to look for outliers
model <- selected_model
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


ouliersRemoval(selected_model)
