library(fgsea)
library(tidyverse)
library(dplyr)
library(here)
library(grid)
library(gridExtra)
library(fs)

outdir <- here("output", "TCGA-LUAD",
               "genset_results.R")
dir_create(outdir, recurse=TRUE)

comps <- here("code", c(
  'TCGA-LUAD/KRAS_codon_12-NTRK3_co_kras',
  'TCGA-LUAD/KRAS_G12C-KEAP1_co_kras',
  'TCGA-LUAD/KRAS_G12C-PTPRD_co_kras',
  'TCGA-LUAD/KRAS_G12C-RBM10_co_kras',
  'TCGA-LUAD/KRAS_G12C-TP53_co_kras',
  'TCGA-LUAD/KRAS_G12V-TP53_co_tar'
)) %>%
    paste0("_de.csv")
sets <- here("code", "sets", "sets.gmt") %>%
    gmtPathways()
##gl <- vector("list", length(comps))
read_data <- function(fname){
  de <- read.table(fname, sep=',', header=T, row.names=1)
  ranks <- -log10(de$pvalue) * sign(de$log2FoldChange)
  ranks[ranks == Inf] <- 300
  ranks[ranks == -Inf] <- -300
  names(ranks) <- rownames(de)
  ## Read in GSEA data
  ## Note, not subsetting by top and bottom. Subset by padj later on.
  fname <- str_replace(fname, "_de.csv", "_gsea.RDS")
  results <- readRDS(fname)
  result.list <- list(results=results, ranks=ranks)
  ##results$ranks <- ranks
  result.list
}
result.list <- map(comps, read_data)
kras <- strsplit(basename(comps), "-") %>%
    sapply("[", 1) %>%
    str_replace_all("_", " ")
comutant <- strsplit(basename(comps), "-") %>%
    sapply("[", 2) %>%
    strsplit("_") %>%
    sapply("[", 1)
plot_titles <- paste0(kras, " / ", comutant, " co-mutant") %>%
    paste0(" vs ", kras, " mutant LUAD")

set_names <- names(sets)
new_names <- set_names

# First convert gene sets too long by hand (perserving original format).
# This may need to be added to if tweaking p-values below.
conv <- here("code", "sets",
             "manual_set_conversion.txt") %>%
    read.table(sep="\t")
new_names[match(conv[,1], new_names)] <- conv[,2]

# Step through various search and replace to rename all sets
new_names <- gsub('_', ' ',  new_names, fixed=T)
new_names <- gsub('.', ' ',  new_names, fixed=T)
new_names <- gsub('HALLMARK', 'HM',  new_names, fixed=T)
new_names <- gsub('REACTOME', 'RT',  new_names, fixed=T)
new_names <- gsub('KEGG', 'KG',  new_names, fixed=T)
new_names <- gsub('BIOCARTA', 'BC',  new_names, fixed=T)
new_names <- str_to_title(new_names)
set_conv <- cbind(set_names, new_names) %>%
    as_tibble()

## Save results to table
##outfile <- here("public", "table", "fig6.Rmd",
##                "set_name_conversion.csv")
##fs::dir_create(dirname(outfile))
##write.table(set_conv, outfile, sep=',',
##            col.names=T,
##            row.names=F)
abbreviated_pathways <- function(x, abbrv){
    abbrv2 <- abbrv$new_names %>%
        setNames(abbrv$set_names)
    results <- x$results %>%
        mutate(pathway_old=pathway,
               pathway=abbrv2[pathway])
    x$results <- results
    x
}
result.list2 <- result.list %>%
    map(abbreviated_pathways, set_conv)
nms <- set_conv$new_names %>% setNames(set_conv$set_names)
names(sets) <- nms[names(sets)]

fig6.list <- list(plot_titles=plot_titles,
                  sets=sets,
                  result.list=result.list2)
saveRDS(fig6.list, file.path(outdir, "fig6_data.rds"))
