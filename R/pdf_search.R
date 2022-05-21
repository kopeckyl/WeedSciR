install.packages(c("dbplyr", "RSQLite"))
Product <- "Roundup PowerMax II"
AI <- "glyphosate"
company <- "Bayer"
conc_g_L <- 540
conc_lb_A <- 4.5
spray_vol_gal_A <- 20
spray_vol_L_ha <- spray_vol_gal_A*9.35 
ID <- 1
type <- "herbicide"
tank_volume_L <- 2
tank_volume_gal <- tank_volume_L/3.785

entry <- data.frame(ID, type, Product,company, AI,conc_g_L, conc_lb_A)

herbicides <- NULL

library(pdfsearch)
file <- "~/Downloads/roundup_label.pdf"
AI <- "glyphosate"
result <- keyword_search(file, 
                         keyword = c(AI, 'active ingredient'),
                         path = TRUE)

head(result$line_text)[[1]]


library(dplyr)
library(dbplyr)
herbicides <- DBI::dbConnect(RSQLite::SQLite(), "~/Desktop/WeedSci_package/weedSci.db")
src_dbi(herbicides)


