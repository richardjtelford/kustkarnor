library("drake")
library("tidyverse")
library("assertr")
library("readxl")
library("vegan")
library("rioja")
library("palaeoSig")
library("Hmisc") # assumes you have mdb-tools installed

plan <- drake_plan(
  #load calibration data

  #load fossil data
  fos0 = read_xlsx("data/Alla kustkärnor med koder_20190903.xlsx", sheet = "Sheet1", skip = 1) %>%
    rename(site = ...1, sample = ...2, depth = ...3, date = ...4, countsum = ...5),
  fos_meta = fos0 %>% select(site, sample, depth, date, countsum),

  #calculate percent
  fos_percent = fos0 %>%
    verify(all.equal(#check countsum matches actual countsums
      countsum,
      fos0 %>%
        select(-site,-sample,-depth,-date,-countsum) %>%
        rowSums()
    )) %>%
    select(-site,-sample,-depth,-date) %>%
    mutate_at(vars(-countsum), ~ {
      . / countsum * 100
    }) %>%
    select(-countsum)

  #diagnostics

  #reconstructions

  #rmarkdown
)

config <- drake_config(plan = plan)

config
