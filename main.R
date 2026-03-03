# ============================================================
# LTAN and Climatology Tools
# ------------------------------------------------------------
# Description: Simple statistics and plotting tools for climatology data. 
#              Development ongoing, all comments and suggestions are appreciated!
# Author: Rachel Dietrich - rachel.dietrich@nrcan-rncan.gc.ca
# Organization: Canadian Forest Service
# Repository: https://github.com/CWFIS/LTAN_tools
# ============================================================


# ---- Packages ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  data.table, dplyr, tidyr, lubridate,
  ggplot2, purrr, rlang, RColorBrewer, openair
)

# ---- Load functions ----
source("functions.R")

# ---- User Inputs ----

fire_name <- "TESTFIRE2026"           # name of fire, 
data<-fread("insert_your_file_here.csv") # path to your CSV,you can put this csv
                                         # in the same folder as this R project!

#data<-fread("Sample_data.csv")  # you can uncomment and use the sample data

year_to_highlight <- 2025
start_year<-NULL # NULL= all data, change to a start year otherwise
stations_to_analyze <- NULL # NULL = all stations or use c("Station ID1"...)
variables <- c("FFMC","ISI","DMC","DC","BUI","FWI","ACCUM_PRECIP")
percentiles <- c(
  p10 = 0.10,
  p25 = 0.25,
  p50 = 0.50,
  p75 = 0.75,
  p90 = 0.90,
  p95 = 0.95
)

# Define spread thresholds as a list of expressions, you can add as many as 
# you'd like. Can choose from FFMC, ISI, DMC, DC, BUI, FWI. Does not support 
# windspeed yet. 
spread_thresholds <- list(
  quote(ISI >= 6 & BUI >= 50),
  quote(ISI >= 7.5 & BUI >= 50),
  quote(ISI >= 10 & BUI >= 50)
)

# Match your data columns EXACTLY, ensure all these variables exist in your data
# your all of your date / time columns should be numeric 
names(data) #will help with this
var_map <- list(
  STATION     = "STATION_NAME",
  FFMC        = "FINE_FUEL_MOISTURE_CODE",
  ISI         = "INITIAL_SPREAD_INDEX",
  DMC         = "DUFF_MOISTURE_CODE",
  DC          = "DROUGHT_CODE",
  BUI         = "BUILDUP_INDEX",
  FWI         = "FIRE_WEATHER_INDEX",
  h_TEMP      = "HOURLY_TEMPERATURE",
  h_RH        = "HOURLY_RELATIVE_HUMIDITY",
  h_WS        = "HOURLY_WIND_SPEED",
  h_WD        = "HOURLY_WIND_DIRECTION",
  PRECIP_24h  = "PRECIPITATION",
  DAY         = "DAY",
  MONTH       = "MONTH",
  YEAR        = "YEAR",
  HOUR        = "HOUR"
)

### No user input past this point ###
# ---- Reshape Data ----
wx_std <- import_weather_data(data, var_map, start_year = start_year)

if (!is.null(stations_to_analyze)) {
  wx_std <- wx_std %>% filter(STATION %in% stations_to_analyze)
}


wx_std <- calc_accum_precip(wx_std)
wx_daily <- aggregate_daily(wx_std)

# ---- Create folders ----
fire_folder <- fire_name
plot_folder <- file.path(fire_folder, "plots")
table_folder <- file.path(fire_folder, "tables")
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(table_folder, recursive = TRUE, showWarnings = FALSE)

# ---- Calculate monthly percentiles ----
wx_pct <- calc_monthly_percentiles(
  wx_std,
  variables   = variables,
  percentiles = percentiles
)
write.csv(wx_pct, file.path(table_folder, "monthly_percentiles.csv"), row.names = FALSE)

# ---- Calculate days over threshold for defined thresholds ----
thres_days_list <- purrr::map(spread_thresholds, function(cond) {
  stats <- calc_exceedance_stats(
    wx_daily,
    !!cond,
    percentiles = percentiles   # ← REQUIRED NOW
  ) %>%
    dplyr::mutate(Condition = rlang::as_label(cond))
  
  stats
})
thres_days <- bind_rows(thres_days_list)
write.csv(thres_days, file.path(table_folder, "threshold_days.csv"), row.names = FALSE)

# ---- Average wx for days over thresholds ----
wx_at_thresholds<- calc_monthly_weather_for_thresholds(wx_std, spread_thresholds)
write.csv(wx_at_thresholds, file.path(table_folder, "wx_at_thresholds.csv"), row.names = FALSE)

# ---- Generate percentile ribbon plots ----
plots_list <- wx_daily %>% 
  split(.$STATION) %>% 
  map(function(station_data) { 
    station_name <- unique(station_data$STATION)
    map(variables, function(var) { 
      plot_percentile_ribbons_station( 
        station_data, 
        !!sym(var), 
        year_line = year_to_highlight, 
        station_name = station_name ) 
      }) %>% set_names(variables) 
    })
## ---- View percentile plots ----
for(station in names(plots_list)) {
  for(var in names(plots_list[[station]])) {
    print(plots_list[[station]][[var]])
  }
}
## ---- Save percentile plots ----
for(station in names(plots_list)) {
  station_dir <- file.path(plot_folder, station)
  dir.create(station_dir, recursive = TRUE, showWarnings = FALSE)
  for(var in names(plots_list[[station]])) {
    ggsave(
      filename = file.path(station_dir, paste0(station, "_", var, ".png")),
      plot = plots_list[[station]][[var]],
      width = 10, height = 4,
      bg = "white"
    )
  }
}

# ---- Generate days over threshold plots ----
## you can change what months to display here! ##
threshold_plots <- wx_daily %>%
  split(.$STATION) %>%
  map(~ plot_weather_events(.x, !!!spread_thresholds,start_month = 4, end_month = 10))
## ---- View threshold plots ----
threshold_plots
## ---- Save threshold plots ----
for(station in names(threshold_plots)) {
  station_dir <- file.path(plot_folder, station)
  dir.create(station_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = file.path(station_dir, paste0(station, "_threshold_events.png")),
    plot = threshold_plots[[station]],
    width = 10, height = 4,
    bg = "white"
  )
}

# ---- Generate wind roses ----
plots_windroses <- generate_windroses(wx_std)
## ---- Plot wind roses ----
for (station in names(plots_windroses)) {
  grid::grid.newpage()
  
  # Draw the windrose plot
  grid::grid.draw(plots_windroses[[station]])
  grid::grid.text(
    label = station,
    x = 0.5,       
    y = 0.95,        
    gp = grid::gpar(fontsize = 16, fontface = "bold")
  )
}

## ---- Save wind roses ----
for (station in names(plots_windroses)) {
  station_dir <- file.path(plot_folder, station)
  dir.create(station_dir, recursive = TRUE, showWarnings = FALSE)
  png(
    filename = file.path(station_dir, paste0(station, "_WindSpeed_MONTH.png")),
    width = 4000, 
    height = 3000,
    res = 300,
    bg = "white"
  )
  grid::grid.newpage()
  grid::grid.draw(plots_windroses[[station]])
  grid::grid.text(
    label = station,
    x = 0.5,        
    y = 0.95,         
    gp = grid::gpar(fontsize = 24, fontface = "bold")  
  )
  dev.off()
}

