###### Dose response plots #####

library(drc)
library(tidyverse)
library(plotly)
library(broom)

# load data
df <- S.alba

# fit dose response

FitDoseResponse(dose = "Dose", "DryMatter", variable = "Herbicide",
                df = df, fit_weibull = TRUE, get_AIC_table = TRUE, get_summary = TRUE)

applyBoxCox(selected_model, data = df)
t <- selected_model$call
# plot 1 - DRC version

plot(selected_model, bp=.2, bty="l",
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
  geom_smooth(aes(color = Herbicide),method = drm, method.args = list(fct = L.4()), se = F) +
  theme_light() +
  labs(title= "", x = "Dose (g a.e /ha)",  y = "Biomass") +
  theme(legend.position = "bottom")


## create function for the first plot ---------
selected_model$data
plot_DR <- function(DR_model) {

  # general sd function
  std_mean <- function(x) sd(x,na.rm=TRUE)/sqrt(length(x))

  # create a dataframe to plot
  df_plot <- tibble(selected_model$data) %>%
    group_by(Dose, Herbicide) %>%
    summarize(var_value = mean(DryMatter, na.rm=TRUE),
              sd = std_mean(DryMatter)) %>%
    mutate(Dose = Dose + 0.1)

}



ggplotly(p1, tooltip = "text")


