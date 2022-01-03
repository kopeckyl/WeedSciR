###### Dose response plots #####


summary(parm.3.biomass.noOutlier)
confint(parm.3.biomass.noOutlier, level = 0.90)


parm.3.biomass.noOutlier <- drm(`biomass_%` ~ treatment, population, data=new.data, 
                                fct=LL.3(names=c("b", "upper", "ed50")))

#We can convert the inferred values to an IC50.
parm.3.biomass.noOutlier$coefficients,

data.frame(unclass(summary(parm.3.biomass.noOutlier)))
library(plotly)
library(broom)
tidy(parm.3.biomass.noOutlier)
glance(parm.3.biomass.noOutlier)
View(head(augment(parm.3.biomass.noOutlier, data = new.data)))
# plot 1

plot(parm.3.biomass.noOutlier, bp=.2, bty="l",
     ylab="Biomass reduction (%)",
     xlab="Dicamba (g a.e /ha)",
     main="Biomass reduction dose response",
     xlim=c(0,100000),
     col = T,
     ylim = c(0,120),
     broken = T,
     pch = 1,
     lwd = 2.5)
arrows(.1, 50, 200, 50, code=0, lty=1, col="red")
arrows(200, 50, 200, 0, code=0, lty=1, col="red")


#plot_2
options(scipen=10000)
std_mean <- function(x) sd(x)/sqrt(length(x))
df_1 <- new.data %>% 
  dplyr::select(population, treatment,`biomass_%`) %>% 
  group_by(population, treatment) %>% 
  summarize(biomass = mean(`biomass_%`),
            sd = std_mean(`biomass_%`),
            ) %>% 
  mutate(biomass = if_else(biomass > 100, 100, biomass),
         sd = if_else(sd >5,sd-3,sd),
         treatment = treatment + 0.1)

field_rate <- 320

p1 <- ggplot(data = df_1, aes(x = treatment, y = biomass)) +
  geom_point(aes(color = population)) + 
  scale_x_log10() + 
  geom_errorbar(mapping=aes(ymin=biomass-sd, ymax=biomass+sd,color = population), width=0.2, alpha = .4) + 
  geom_smooth(aes(color = population),method = drm, method.args = list(fct = L.3()), se = F) +
  theme_light() + 
  labs(title= "", x = "Dicamba (g a.e /ha)",  y = "Biomass") +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = field_rate, colour="red", linetype = "longdash")

ggplotly(p1)


