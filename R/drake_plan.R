library("drake")
library("tidyverse")
library("assertr")
#library("readxl")
library("vegan")
library("rioja")
library("palaeoSig")
#library("Hmisc") # assumes you have mdb-tools installed

plan <- drake_plan(
  #load calibration data

  #load fossil data
  fos0 = readxl::read_xlsx("data/Alla kustkärnor med koder_20190903.xlsx", sheet = "Sheet1", skip = 1) %>%
    rename(site = ...1, sample = ...2, depth = ...3, date = ...4, countsum = ...5),
  fos_meta = fos0 %>% select(site, sample, depth, date, countsum),

  #calculate percent
  fos_percent = fos0 %>%
    #check countsum matches actual countsums
    verify(all.equal(
      countsum,
      fos0 %>%
        select(-site,-sample,-depth,-date,-countsum) %>%
        rowSums()
    )) %>%
    select(-countsum) %>%
    #merge taxa
    pivot_longer(cols = -(site:date), names_to = "taxa", values_to = "count") %>%
    mutate(taxa = str_remove(taxa, "\\.{3}\\d*$")) %>% #dups end with three dots and then digits
    group_by(site, sample, depth, date, taxa) %>%
    #merge duplicate taxa
    summarise(count = sum(count)) %>%
    #calculate percent
    mutate(percent = count/sum(count) * 100) %>%
    select(-count) %>%
    pivot_wider(names_from = taxa, values_from = percent) %>%
    #remove metadata
    ungroup() %>%
    select(-site,-sample,-depth,-date),

  #diagnostics
  #goodness of fit
  #analogue
  #reconstruction significance

  #reconstructions
  #model
  #predictions

  #rmarkdown
  output = rmarkdown::render(knitr_in("kustkarnor.Rmd"))
)

config <- drake_config(plan = plan)

config
