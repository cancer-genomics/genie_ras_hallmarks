library(genie.6.1)
library(tidyverse)
data(patient_universe)
library(here)
# Add cancer labels for manuscript figures
x <- tibble(cancer=unique(patient_universe$cancer),
            cancer_label="")%>%
    mutate(cancer_label = gsub(" cancer", "", cancer))%>%
    mutate(cancer_label = gsub(" tumor", "", cancer_label))%>%
  
    mutate(cancer_label = Hmisc::capitalize(cancer_label))%>%
    mutate(cancer_label = recode(cancer_label, 
                                 "Colorectal" = "CRC", 
                  "Non-small cell lung" = "NSCLC", 
                  "Cancer of unknown primary" = "CUP", 
                  "B-lymphoblastic leukemia/lymphoma" = "BLL/Lymphoma",
                  "Neuroendocrines" = "Neuroendocrine",
                  "Small cell lung" = "SCLC", 
                  "Myelodysplastic/myeloproliferative syndromes" = "MDS/MPN",
                  "Non-melanoma skins" = "Non-melanoma skin", 
                  "Mature t and nk neoplasms" = "T / NK neoplasms", 
                  "Adrenocortical carcinoma" = "Adrenocortical",
                  "Cns" = "CNS" ))

write_csv(x, here("data", "cancer_labels.csv"))
