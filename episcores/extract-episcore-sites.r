table1e <- read.csv("table1e.csv")
proteins <- c(
  "OSM", "S100A12", "CXCL9", "VEGFA", 
  "TGFA", "IGFBP1", "MMP12", "HGF") 
table1e <- table1e[table1e$gene %in% proteins,]
table(table1e$gene)
#CXCL9     HGF  IGFBP1   MMP12     OSM S100A12    TGFA   VEGFA 
#   56      72     120     257     125      82     100     159 
write.csv(table1e, file="episcore-sites.csv")