library("drake")
library("tidyverse")
library("assertr")
#library("readxl")
library("vegan")
#devtools::install_github("gavinsimpson/ggvegan", upgrade = FALSE)
library("ggvegan")
#devtools::install_github("nsj3/rioja", upgrade = FALSE)#need latest version
library("rioja")
library("palaeoSig")
library("here")
#devtools::install_github("richardjtelford/ggpalaeo", upgrade = FALSE)
library("ggpalaeo")
library("magrittr")
library("ggnewscale")
library("conflicted")
conflict_prefer("filter", "dplyr")

#database connection
conn <- function(db){
  DBI::dbConnect(
    RSQLite::SQLite(),
    dbname = db)
}

source("R/coverage_plot.R")#replacement for palaeoSig version allowing for grouped data

#import plans
source("R/load_data_plan.R")#data loading

#analyses
analysis_plan <- drake_plan(
  #site map
  site_map = {
    mp <- map_data("world", xlim = c(-10, 50), ylim = c(40, 75))
  ggplot(envT, aes(x = longitude, y = latitude, size = TN)) +
    geom_map(map = mp, data = mp, aes(map_id = region), inherit.aes = FALSE, fill = "grey70", colour = "grey50") +
    geom_point(alpha = 0.5, colour = "red") +
    coord_quickmap() +
    scale_size_area() +
    labs(x = "°E", y = "°N", size = "log(TN)")
    },

  ## Pairs plot
  pairs_plot = envT %>%
    select(-siteId, -countryId, -longitude, -exposed) %>%
    GGally::ggpairs(),


  ## modern ordination
  mod_decorana = decorana(sqrt(spp)),

  mod_cca = cca(sqrt(spp) ~ TN, data = envT),

  mod_cca1 = cca(sqrt(spp) ~ salinity + depth + TN + TP, data = envT),
  mod_cca1_plot = autoplot(mod_cca1),


  #####diagnostics####
  #abundance of new taxa
  fos_new_taxa_abun = fos0 %>%
    filter(!taxa %in% names(spp), count > 0) %>%
    left_join(fos_new_taxa, by = c("taxa" = "Code")) %>%
    group_by(taxa, Taxon) %>%
    summarise(n = n(), max= max(percent)) %>%
    arrange(desc(max)),

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

  #time-track
  time_track_plot = {
    tt <- analogue::timetrack(
      X = spp,
      passive = select(fos_percent,-(site:date)),
      env = envT, transform = "sqrt", formula = ~ TN)

  autoplot(tt$ordination, layers = c("sites", "biplot")) +
    new_scale_color() +
    geom_path(
      data = as_tibble(tt$fitted.values) %>% bind_cols(fos_percent),
      aes(x = CCA1, y = CA1, colour = site)) +
    scale_colour_viridis_d()
  },

  #rmarkdown
  output = rmarkdown::render(knitr_in("kustkarnor.Rmd"))
)

#### combine plans
plan <- bind_rows(data_plan, analysis_plan)


#### configure plan
config <- drake_config(plan = plan)

config
