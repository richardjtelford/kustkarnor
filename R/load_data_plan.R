## plan to load data

data_plan <- drake_plan(
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
    semi_join(spp0, by = c(siteId = "sampleId")) %>%
    arrange(siteId),
  # all.env<-!is.na(rowSums(envT[,2:5]))
  # env$siteId[!all.env]

  #load species data
  spp0 = Hmisc::mdb.get(file_in("data/processCounts.mdb"), tables = "FinalPercent", stringsAsFactors = FALSE) %>%
    mutate(
      sampleId = tolower(sampleId),
      perc = as.vector(perc)
    ) %>%
    #sites with chemistry
    semi_join(env, by = c("sampleId" = "siteId" )) %>%
    #merge taxa merged in fossil data
    mutate(taxonCode = case_when(
      taxonCode %in% c("DiaMon", "DiaTen") ~ "DiaCom",
      taxonCode %in% c("EpiAdn", "EpiSor", "EpiTur") ~ "EpiCom",
      TRUE ~ taxonCode
    )) %>%
    group_by(countryId, sampleId, taxonCode) %>%
    summarise(perc = sum(perc)) %>%
    ungroup() %>%
    pivot_wider(
      names_from = "taxonCode",
      values_from = "perc",
      values_fill = list(perc = 0)) %>%
    arrange(sampleId),

  #check then remove meta data
  spp = spp0 %>%
    verify(identical(sampleId, envT$siteId)) %>%
    select(-countryId, -sampleId),

  ####load fossil data ####
  fos0 = readxl::read_xlsx(
    file_in(!!here("data", "Alla kustkärnor med koder_20190903.xlsx")),
    sheet = "Sheet1", skip = 1) %>%
    rename(site = ...1, sample = ...2, depth = ...3, date = ...4, countsum = ...5) %>%
    #check countsum matches actual countsums
    verify(all.equal(
      countsum,
      select(., -(site:countsum)) %>%
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
    #delete taxa with all zero abundances
    group_by(taxa) %>%
    filter(sum(percent) > 0),

  #calculate percent
  fos_percent = fos0 %>%
    select(-count) %>%
    pivot_wider(names_from = taxa, values_from = percent),

  ## New fossil taxa
  fos_new_taxa = readxl::read_xlsx(file_in(!!here("data", "Alla kustkärnor med koder_20190903.xlsx")), sheet = "New taxa")
)
