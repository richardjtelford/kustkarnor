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
    gather(key = taxa, value = count, -site, -sample,-depth,-date) %>%
    mutate(taxa = gsub("\\.\\.\\.\\d*$", "", taxa)) %>% #dups end with three dots and then digits
    group_by(site, sample, depth, date, taxa) %>%
    summarise(count = sum(count)) %>%
    #calculate percent
    mutate(percent = count/sum(count) * 100) %>%
    select(-count) %>%
    spread(key = taxa, value = percent) %>%
    #remove metadata
    ungroup() %>%
    select(-site,-sample,-depth,-date),

  #diagnostics

  #reconstructions

  #rmarkdown
  output = rmarkdown::render(knitr_in("kustkarnor.Rmd"))
)

config <- drake_config(plan = plan)

config
