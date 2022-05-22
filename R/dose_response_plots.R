###### Dose response plots #####

library(drc)
library(tidyverse)
library(plotly)
library(broom)

# load data
df <- S.alba

# fit dose response

FitDoseResponse(dose = "Dose", "DryMatter", variable = "Herbicide",
                df = df, fit_weibull = TRUE, get_AIC_table = TRUE, get_summary = TRUE,
                remove_outliers = TRUE, applyBC = TRUE, plot_curve = TRUE)





applyBoxCox(selected_model)
t <- selected_model$call
# plot 1 - DRC version

plot(Curve_fit$model, bp=.2, bty="l",
     ylab="Biomass reduction (%)",
     xlab="Dicamba (g a.e /ha)",
     main="Biomass reduction dose response",
     xlim=c(0,100000),
     col = T,
     ylim = c(0,5),
     broken = T,
     pch = 1,
     lwd = 2.5)
arrows(.1, 50, 200, 50, code=0, lty=1, col="red")
arrows(200, 50, 200, 0, code=0, lty=1, col="red")


#plot_2
options(scipen=10000)
std_mean <- function(x) sd(x,na.rm=TRUE)/sqrt(length(x))

df_1 <- tibble(df) %>%
  group_by(Dose, Herbicide) %>%
  summarize(var_value = mean(DryMatter, na.rm=TRUE),
            sd = std_mean(DryMatter),
            ) %>%
  mutate(Dose = Dose + 0.1)

field_rate <- 320

p1 <- ggplot(data = df_1, aes(x = Dose, y = var_value)) +
  geom_point(aes(color = Herbicide,
                 text = paste("Dose:", Dose,
                              "\nHerbicide:", Herbicide,
                              "\nBiomass:", var_value))) +
  scale_x_log10() +
  geom_errorbar(mapping=aes(ymin=var_value-sd, ymax=var_value+sd,color = Herbicide), width=0.2, alpha = .4) +
  geom_smooth(aes(color = Herbicide),method = drm, method.args = list(fct = L.4()), se = F, fullrange=TRUE) +
  theme_light() +
  labs(title= "", x = "Dose (g a.e /ha)",  y = "Biomass") +
  theme(legend.position = "bottom")

t <- Curve_fit$data
cols <- which(names(selected_model$data) == 'variable')
names(t)[cols] <- paste0('variable', seq_along(cols))
names(t)

## create function for the first plot ---------
selected_model$data
plot_DR <- function(DR_model, legend_title = "Group", y_label = "Response", x_label = "Rate") {

  # general sd function
  std_mean <- function(x) sd(x,na.rm=TRUE)/sqrt(length(x))

  # correct model dataframe and create a dataframe to plot
  DR_data <- DR_model$model$data
  colnames(DR_data) <- c("dose", "var_name", "variable1", "variable2", "weights")
  #cols <- which(names(DR_data) == 'variable')
  #names(DR_data)[cols] <- paste0('variable', seq_along(cols))
  df_plot <- tibble(DR_data) %>%
    group_by(dose, variable2) %>%
    summarise(var_value = mean(var_name, na.rm=TRUE),
              sd = std_mean(var_name)) %>%
    mutate(dose = dose + 0.01)

  # get

  # generate plot
  p1 <- ggplot(data = df_plot, aes(x = dose, y = var_value)) +
    geom_point(aes(color = variable2)) +
    scale_x_log10() +
    geom_errorbar(mapping=aes(ymin= var_value - sd,
                              ymax= var_value + sd,
                              color = variable2),
                  width=0.2, alpha = .4) +
    geom_smooth(aes(color = variable2),
                method = "drm",
                method.args = list(fct = L.4()), se = F,
                formula = 'y ~ x') +
    theme_light() +
    labs(title= "", x = x_label,  y = y_label) +
    theme(legend.position = "bottom") +
    guides(color=guide_legend(title=legend_title))

  return(p1)
}

plotly_build(plot_DR(Curve_fit, x_label = "Herbicide (g ai/ha)"))  %>%
  layout(legend = list(orientation = "h",   # show entries horizontally
                       xanchor = "center",  # use center of legend as anchor
                       x = 0.5))

Curve_fit$model$data
