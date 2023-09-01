library(data.table)

panel <- fread("panel-reduced.csv")
panel$chr <- paste0("chr",panel$chr)
panel$start <- panel$start - 1

uncovered <- fread("twist-clarification/all_target_segments_not_covered_by_probes_withshifting.bed")

colnames(uncovered) <- c("chr","start","end")

uncovered$idx <- sapply(1:nrow(uncovered), function(i) 
    with(uncovered, 
         which(
             panel$chr==uncovered$chr[i]
             & panel$start <= uncovered$start[i]
             & panel$end >= uncovered$end[i])))

stopifnot(!is.list(uncovered$idx))
stopifnot(all(is.integer(uncovered$idx)))

uncovered$pct <- with(uncovered, (end-start+1)/(panel$end[idx]-panel$start[idx]+1))

panel$uncovered <- 1:nrow(panel) %in% uncovered$idx
panel$is.total <- 1:nrow(panel) %in% uncovered$idx[uncovered$pct >= 1]

## all uncovered regions that are more than a single CpG site were included to estimate blood cell counts
stopifnot(all(!(panel$uncovered & panel$end-panel$start > 1) | panel$source == "blood-cell-types"))

## for each blood cell type, > 50% of the regions specific to that cell type can be targeted 
stopifnot(all(
    sapply(
        unique(panel$details[panel$uncovered & panel$source == "blood-cell-types"]),
        function(type) {
            mean(panel$uncovered[panel$details==type])
        })
    < 0.5))

## for each ancestry, > 90% of the CpG sites specific to that ancestry can be targeted
ancestry <- fread("ancestry/output/sites.csv")
ancestry$uncovered <- ancestry$cpg %in% panel$details[panel$uncovered]
stopifnot(all(
    sapply(unique(ancestry$ancestry),function(name)
        mean(ancestry$uncovered[ancestry$ancestry == name]))
    < 0.1))

## for each episcore, > 97% of the CpG sites in the model can be targeted
episcores <- fread("episcores/episcore-sites.csv")
episcores$uncovered <- episcores$cpg %in% panel$details[panel$uncovered]
stopifnot(all(
    sapply(unique(episcores$gene), function(name)
        mean(episcores$uncovered[episcores$gene == name]))
    < 0.03))

## for all remaining sources (mainly predictive models), >93% of contributing CpG sites can be targeted
remaining <- unique(panel$source[panel$uncovered & ! panel$source %in% c("blood-cell-types","ancestry","episcores")])
pct.uncovered <- sapply(remaining, function(name) 
    mean(panel$uncovered[panel$source==name]))
stopifnot(all(pct.uncovered < 0.07))

library(knitr)
kable(sort(pct.uncovered))
## |                    |         x|
## |:-------------------|---------:|
## |bmi                 | 0.0029940|
## |alcohol-consumption | 0.0037879|
## |crp                 | 0.0087719|
## |smoking-cessation   | 0.0123457|
## |hdl                 | 0.0128205|
## |dunedin-pace        | 0.0194805|
## |breast-cancer       | 0.0400000|
## |dunedin-poam38      | 0.0540541|
## |cotinine            | 0.0625000|



