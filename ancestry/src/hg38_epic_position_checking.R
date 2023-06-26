library(data.table)

## Read in the panel
panel <- fread("panel.csv")

## Download EPIC manifest
online_file <- "https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip"
download.file(url = online_file, destfile = "epic_manifest.zip")

## Unzip EPIC manifest
unzip("epic_manifest.zip")

## Read in EPIC manifest
manifest_epic <- fread("infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7)

## Keep relevant annotation information
manifest_epic <- manifest_epic[,c("Name", "CHR_hg38", "Start_hg38", "End_hg38")]
colnames(manifest_epic) <- c("name_manifest", "chr_hg38_manifest", "start_hg38_manifest", "end_hg38_manifest")
## Merge annotation information with panel (EPIC)
panel_epic <- merge(panel, manifest_epic, by.x = "details", by.y = "name_manifest")

## The panel start and end sit exactly in-between the manifest start and end
## Probably because it's chr:pos rather than a 3-base start and end
head(panel_epic[,c("start", "end", "start_hg38_manifest", "end_hg38_manifest")])

## Take the middle value of the manifest start and end
panel_epic$pos_hg38_manifest <- (panel_epic$start_hg38_manifest + panel_epic$end_hg38_manifest)/2

## Create a chr:pos variable for the panel and the manifest
panel_epic$panel_loc <- paste0(panel_epic$chr, ":", panel_epic$start)
panel_epic$manifest_loc <- paste0(panel_epic$chr_hg38_manifest, ":", panel_epic$pos_hg38_manifest)

table(panel_epic$panel_loc == panel_epic$manifest_loc)
# FALSE  TRUE 
# 36     3533 

## Investigate where false
panel_epic[which(panel_epic$panel_loc != panel_epic$manifest_loc),]

# All positions on the panel that aren't exactly in the middle of the start and 
# end of the manifest are either AT the start or end position - zero outliers