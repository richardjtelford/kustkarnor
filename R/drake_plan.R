library("drake")
library("tidyverse")
library("assertr")
#library("readxl")
library("vegan")
library("rioja")
library("palaeoSig")
library("here")
library("ggpalaeo")
library("magrittr")
#library("Hmisc") # assumes you have mdb-tools installed

source("R/coverage_plot.R")

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
    semi_join(spp0, by = c(siteId = "sampleId")) %>%
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
    rename(site = ...1, sample = ...2, depth = ...3, date = ...4, countsum = ...5),

  #calculate percent
  fos_percent = fos0 %>%
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
    select(-count) %>%
    #delete taxa with all zero abundances
    group_by(taxa) %>%
    filter(sum(percent) > 0) %>%
    pivot_wider(names_from = taxa, values_from = percent) ,

  #####diagnostics####
  #coverage
  coverage_plots = fos_percent %>%
    #nest
    group_by(site) %>%
    select(-(sample:date)) %>%
    #coverage plot
    coverage_plot(spp = spp, fos = .) +
    facet_wrap(~site),

  #goodness of fit
  reslen_plot = analogue::residLen(
    X = sqrt(spp),
    env = envT$TN,
    passive = fos_percent %>% select(-(site:date)) %>%
      sqrt()
  ) %>%
    autoplot(df = fos_percent %>% select(site:date),
                 x_axis = "date") +
    facet_wrap( ~ site) +
    labs(x = "Date CE", y = "Squared residual distance", fill = "Goodness of fit"),

  #analogue
  ana_qual = analogue_distances(
    spp = spp/100,
    fos = fos_percent %>%
      select(-(site:date)) %>%
      divide_by(100)
  ) %>%
    autoplot(df = fos_percent %>% select(site:date),
              x_axis = "date") +
    facet_wrap( ~ site) +
    labs(x = "Date CE", y = expression(Chord^2~distance)),

  #####Transfer function
  #WA
  mod_wa = WA(sqrt(spp), envT$TN, mono = TRUE) %>%
    crossval(),


  #reconstructions
  pred = predict(
    object = mod_wa,
    newdata = select(fos_percent, -(site:date))
  ) %>%
    pluck("fit") %>%
    as_tibble() %>%
    bind_cols(select(fos_percent, site:date)),

  pred_plot = ggplot(pred, aes(x = date, y = WA.m, colour = site)) +
    geom_line() +
    labs(x = "Date CE", y = "Log(TN?units)", colour = "Site"),

  #reconstruction significance
  recon_sig = fos_percent %>%
    group_by(site) %>%
    select(-(sample:date)) %>%
    nest() %>%
    mutate(data = map(data, sqrt)) %>%
    mutate(
      recon_sig = map(
        .x = data,
        ~randomTF(spp = sqrt(spp), fos = .x, env = envT$TN, fun = WA, mono = TRUE, col = 3, n = 999)),
      recon_sig_plot = map(recon_sig, autoplot, variable_names = "log(TN)"),
      recon_sig_plot = map2(.x = recon_sig_plot, .y = site, ~{.x + ggtitle(.y)})
      ),

  #rmarkdown
  output = rmarkdown::render(knitr_in("kustkarnor.Rmd"))
)

config <- drake_config(plan = plan)

if(FALSE){#testbed

.x <-   recon_sig %>% ungroup() %>% slice(1) %>% pull(recon_sig)
.x <- .x[[1]]
palaeoSig:::fortify_palaeosig(sim = .x$sim, variable_names = "log(TN)",
                  p_val = 0.05, nbins = 20, top = 0.7, PC1 = .x$MAX,
                  EX = .x$EX)
autoplot_sig(x_fort, xlab = "Proportion variance explained",
             xmin = 0)

  sort(setdiff(names(fos_percent), names(spp)))
  sort(setdiff(names(spp), names(fos_percent)))
}

config
