### Herbicide rate calculation #####


# dose response
type <- c("solid", "liquid") # type of product
concentration <- 2.9 # g ai/ha or lb ai/ha
unit <- c("metric", "US")
Volume <- 20 # l/ha or gal/A
tank <- 0.0105 # in Liters or gallons
rate_1X <- 22



# if dose response
DR <- TRUE
n_rates <- 9
tritation <- c("fold_log","fold_2","fold_10")


