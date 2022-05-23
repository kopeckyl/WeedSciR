# inheritance functions


library(tidyverse) # data wrangling
library(readxl) # read excel files
library(ggpubr) # plot package
library(ggpmisc) # plot aesthetic package
library(cowplot)

# generate random data ------------
set.seed(1234)
f2 <- rnorm(2000, mean=10, sd=2)
sample(x = 1:5, size = 10, replace = TRUE)
hist(f2)

f1 <- rnorm(n = 500, mean = 10, sd = 0.2)

R <- rnorm(n = 500, mean = 15, sd = 0.2)
hist(R)
S <- rnorm(n = 500, mean = 5, sd = 0.5)
hist(S)

# Trial A
f2 <- data.frame(biomass = rnorm(2000, mean=10, sd=1), genotype = rep("F2", 2000))
f1 <- data.frame(biomass = rnorm(n = 500, mean = 10, sd = 0.2), genotype = rep("F1", 500))
Dom <- data.frame(biomass = rnorm(n = 500, mean = 15, sd = 0.2), genotype = rep("Dom", 500))
Res <- data.frame(biomass =rnorm(n = 500, mean = 5, sd = 0.5), genotype = rep("Res", 500))
df_A <- tibble(f2 %>% full_join(f1) %>% full_join(Dom) %>% full_join(Res))
# Trial B
f2 <- data.frame(biomass = rnorm(2000, mean=10, sd=1), genotype = rep("F2", 2000))
f1 <- data.frame(biomass = rnorm(n = 500, mean = 10, sd = 0.2), genotype = rep("F1", 500))
Dom <- data.frame(biomass = rnorm(n = 500, mean = 15, sd = 0.2), genotype = rep("Dom", 500))
Res <- data.frame(biomass =rnorm(n = 500, mean = 5, sd = 0.5), genotype = rep("Res", 500))
df_B <- tibble(f2 %>% full_join(f1) %>% full_join(Dom) %>% full_join(Res))

# all data
df <- bind_rows(A = df_A, B = df_B, .id = "group")
df$group <- factor(df$group)


# real data
df <- read_csv("R/df_polled.csv")
df

# Function to pull data -------------
# Description - Check if a data can be pulled via Levene's test

check_pull <- function(df,group, variable, log_var = FALSE){

    if (!is.factor(df[group])) {
     group <- factor(df %>% pull(group))
    }

  group <- df %>% pull(group)
  variable <- df %>% pull(variable)

  if (log_var == TRUE) {
    log_var <- log10(variable)
    t <- car::leveneTest(log_var, group)
  } else {
    t <- car::leveneTest(variable, group)
  }

  if (t$`Pr(>F)`[1] <= 0.05) {
    print("Homogeinety across groups not uniform. Analysis must to be done separately.")
    print("Spitting dataset:")
    df_split <- split(df, df %>% pull(group))
    return(df_split)
  } else {
    print("Data from different groups can be analyzed together")
  }

}

check_pull(df, group = 'population', variable = 'area_cm', log_var = T)




# data summary
df %>%
  group_by(population) %>%
  summarise(n_individuals = n(),
            mean = mean(area_cm),
            sd = sd(area_cm))


# create plot
library(ggridges)
df %>%
  ggplot(aes(x = log(Biomass), y = population, fill = population)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none", legend.title = element_blank(),
        axis.line = element_line(color="black", size = .5),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = expression ("Variable"), y = "Genotype")
