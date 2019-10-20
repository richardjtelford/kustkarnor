library("drake")
library("tidyverse")
library("assertr")
#library("readxl")
library("vegan")
library("rioja")
library("palaeoSig")
library("here")
#library("Hmisc") # assumes you have mdb-tools installed

plan <- drake_plan(
  ####load calibration data####

  #read env data
  env = Hmisc::mdb.get(file_in("data/define.mdb"), tables = "fchem") %>%
    #remove annoying labels
    mutate_if(is.numeric, as.vector) %>%
    as_tibble() %>%
    #siteId to lower case
    mutate(siteId = tolower(siteId)) %>%
    #TN = TDN * 1.5
    mutate(TN = if_else(!is.na(TDN), true = TDN * 1.5, false = TN)) %>%
    #missing Norwegian salinities
    mutate(salinity = case_when(
      siteId == "mo-2" ~ 20, # wide range, estimated from range
      siteId == "s-9" ~ 16, # wide range, estimated from range
      TRUE ~ salinity)) %>%
    #zap missing data
    filter(!is.na(TN)),


  #transform env
  envT = env %>%
    select(siteId, salinity, depth, TN, TP, exposed, countryId, latitude, longitude) %>%
    mutate(
      salinity = sqrt(salinity),
      depth = log(depth),
      TP = log(TP),
      TN = log(TN),
      exposed = exposed == "1") %>% #no exposed sites?
    #sites with diatom data
    semi_join(spp, by = c(siteId = "sampleId")) %>%
    arrange(siteId),
  # all.env<-!is.na(rowSums(envT[,2:5]))
  # env$siteId[!all.env]

  site_map = {
    mp <- map_data("world", xlim = c(-10, 50), ylim = c(40, 75))
  ggplot(envT, aes(x = longitude, y = latitude, colour = countryId)) +
    geom_map(map = mp, data = mp, aes(map_id = region), inherit.aes = FALSE, fill = "grey70") +
    geom_point() +
    coord_quickmap()},

  #load species data
  spp = Hmisc::mdb.get(file_in("data/processCounts.mdb"), tables = "FinalPercent") %>%
    mutate(
      sampleId = tolower(sampleId),
      perc = as.vector(perc)
    ) %>%
    pivot_wider(
      names_from = "taxonCode",
      values_from = "perc",
      values_fill = list(perc = 0)) %>%
    #sites with chemistry
    semi_join(env, spp, by = c("sampleId" = "siteId" )) %>%
    arrange(sampleId),



  ####load fossil data ####
  fos0 = readxl::read_xlsx(
    file_in(here("data", "Alla kustkärnor med koder_20190903.xlsx")),
    sheet = "Sheet1", skip = 1) %>%
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



