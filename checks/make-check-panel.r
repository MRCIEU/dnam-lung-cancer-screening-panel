library(data.table)
library(tidyverse)
sites.filename <- "panel-reduced.csv"
output.filename <- "checks/panel-reduced-py.csv"

sites <- fread(file=sites.filename) %>%
			as.data.frame()
sites$spot.checked <- 0

set.seed(27435)
idx <- sample(1:nrow(sites))
my.sites <- sites[idx,]

write.csv(my.sites, file=output.filename, row.names=F)
