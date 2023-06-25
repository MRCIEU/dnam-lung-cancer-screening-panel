#European vs Mexican ancestry 

Rscript src/check-sites-GSE40279.r \
  output/sites.csv \
  src/check-sites-GSE40279.rmd \
  output/check-sites-GSE40279.html

#European vs African American ancestry

Rscript src/check-sites-GSE117861.r \
  output/sites.csv \
  src/check-sites-GSE117861.rmd \
  output/check-sites-GSE117861.html

#European vs African American ancestry (in cord blood)

Rscript src/check-sites-GSE64940.r \
  output/sites.csv \
  src/check-sites-GSE64940.rmd \
  output/check-sites-GSE64940.html

#Mexican vs Puerto Rican ancestry

Rscript src/check-sites-GSE77716.r \
  output/sites.csv \
  src/check-sites-GSE77716.rmd \
  output/check-sites-GSE77716.html
