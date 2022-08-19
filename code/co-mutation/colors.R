library(magrittr)
library(tidyverse)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(ggpubr)
library(here)
results <- here("output",
            "summarize-com",
            "tables.R",
            "combined.rds")  %>%
    readRDS()

# Co-occurrence of mutations in non-RAS genes

For volcano plots, highlight genes with adjusted p-value less than 0.01 and color code the different cancers.

```{r significance}

saveRDS(highlight, here("public", "table", "colors.rds"))
```
