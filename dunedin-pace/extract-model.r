library(remotes)
remotes::install_github("danbelsky/DunedinPACE")

library(DunedinPACE)

model <- mPACE_Models$model_weights$DunedinPACE
write.csv(
  data.frame(
    cpg=names(model),
    coef=model),
  file="pace-model.csv",
  row.names=F
)
