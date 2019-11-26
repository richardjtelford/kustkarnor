## plan to load data

data_plan <- drake_plan(
  tidy_eval = FALSE,
  ####load calibration data####
  #database
  db = file_in(!!here("data", "define.sqlite")),

  #salinity limits
  salinity_lims = c(1, 10),

  #read env data
  env = tbl(src = conn(db), "fchem") %>%
    #TN = TDN * 1.5
    mutate(TN = if_else(!is.na(TDN), true = TDN * 1.5, false = TN)) %>%
    #missing Norwegian salinities
    mutate(salinity = case_when(
      siteId == "mo-2" ~ 20, # wide range, estimated from range
      siteId == "s-9" ~ 16, # wide range, estimated from range
      TRUE ~ salinity)) %>%
    filter(
      #zap missing data
      !is.na(TN),
      #limit salinity range
      between(salinity, !!salinity_lims[1], !!salinity_lims[2])
    ) %>%
    collect(),


  #transform env
  envT = env %>%
    select(siteId, salinity, depth, TN, TP, exposed, countryId, latitude, longitude) %>%
    mutate(
      depth = log(depth),
      TP = log(TP),
      TN = log(TN),
      exposed = exposed == "1") %>% #no exposed sites?
    #sites with diatom data
    semi_join(spp0, by = c("siteId")) %>%
    arrange(siteId) %>% collect(),
  # all.env<-!is.na(rowSums(envT[,2:5]))
  # env$siteId[!all.env]

  #load species data
  spp0 = {
    con <- conn(db)
    tbl(con, "counts") %>%
    #merge synonyms
    left_join(tbl(con, "Synonyms"), by = c("taxonCode" = "synonymCode")) %>%
    mutate(taxonCode = coalesce(correctCode, taxonCode)) %>%
    #get country code for merges & drop samples with no chem - "RIB16102b"
    inner_join(tbl(con, "fchem") %>% select(siteId, countryId), by  = "siteId") %>%
    left_join(tbl(con, "Merges"), by = c("taxonCode" = "oldtaxonCode", "countryId")) %>%
    mutate(
      taxonCode = coalesce(mergedTaxonCode, taxonCode)
    ) %>%
    #remove excluded taxa
    anti_join(tbl(con, "ExcludedTaxa"), by = "taxonCode") %>%
    #select required columns
    select(siteId, taxonCode, count) %>%
    #sum merged taxa
    group_by(siteId, taxonCode) %>%
    collect() %>%
    assert(not_na, everything()) %>%
    summarise(count = sum(count, na.rm = TRUE)) %>% # , na.rm = TRUE to avoid warning
    mutate(percent = count / sum(count, na.rm = TRUE) * 100) %>%
    #drop samples that don't meet salinity criteria
    semi_join(env, by  = "siteId")
    },

  #clean, check then remove meta data
  spp = spp0 %>%
    #merge taxa merged in fossil data
    mutate(taxonCode = case_when(
      taxonCode %in% c("DiaMon", "DiaTen") ~ "DiaCom",
      taxonCode %in% c("EpiAdn", "EpiSor", "EpiTur") ~ "EpiCom",
      TRUE ~ taxonCode
    )) %>%
    group_by(siteId, taxonCode) %>%
    summarise(percent = sum(percent)) %>%
    #remove rare taxa
    group_by(taxonCode) %>%
    filter(
      max(percent) >= 2, #max percent >= 2
      n() >= 2) %>% # at least 2 occurances
    ungroup() %>%
    pivot_wider(
      names_from = "taxonCode",
      values_from = "percent",
      values_fill = list(percent = 0)) %>%
    arrange(siteId) %>%
    #check identical sampleId to envT
    verify(identical(siteId, envT$siteId)) %>%
    select(-siteId),

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
