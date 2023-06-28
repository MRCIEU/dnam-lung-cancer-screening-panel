#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(tidyverse)
sites.filename <- args[1]
output.filename <- args[2]

sites <- fread(file=sites.filename) %>%
			as.data.frame()
sites$spot.checked <- 0

set.seed(27435)
idx <- sample(1:nrow(sites))
my.sites <- sites[idx,]

write.csv(my.sites, 
	file=output.filename,
	quote = F, 
	row.names=F)
