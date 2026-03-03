# ---- Packages ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, lubridate, ggplot2, tidyr, rlang,
  purrr, RColorBrewer, openair, grid, gridExtra
)

# ------------------------
# 1. Validate variable map
# ------------------------
validate_var_map <- function(data, var_map, required_vars) {
  missing_map <- setdiff(required_vars, names(var_map))
  if(length(missing_map) > 0) stop("Variable map is missing: ", paste(missing_map, collapse = ", "))
  
  missing_cols <- setdiff(unlist(var_map), names(data))
  if(length(missing_cols) > 0) stop("Columns in var_map not found in data: ", paste(missing_cols, collapse = ", "))
  invisible(TRUE)
}

# ------------------------
# 2. Import & standardize weather data
# ------------------------
import_weather_data <- function(data, var_map, start_year = NULL) {
  
  required_vars <- c(
    "STATION", "HOUR", "DAY", "MONTH", "YEAR",
    "FFMC", "ISI", "DMC", "DC", "BUI", "FWI",
    "h_TEMP", "h_RH", "h_WS", "h_WD",
    "PRECIP_24h"
  )
  
  validate_var_map(data, var_map, required_vars)
  
  # map columns to standard names
  transmute_exprs <- lapply(required_vars, function(v) {
    set_names(list(sym(var_map[[v]])), v)
  }) %>% unlist(recursive = FALSE)
  
  out <- transmute(data, !!!transmute_exprs) %>%
    mutate(
      YEAR = as.integer(YEAR),
      DATE = as.Date(sprintf("%04d-%02d-%02d", YEAR, MONTH, DAY))
    )
  
  # optional start year filter
  if (!is.null(start_year)) {
    out <- out %>% filter(YEAR >= start_year)
  }
  
  return(out)
}

# ------------------------
# 3. Accumulate precipitation
# ------------------------
calc_accum_precip <- function(wx) {
  wx %>%
    arrange(STATION, DATE) %>%
    group_by(STATION, YEAR) %>%
    mutate(ACCUM_PRECIP = cumsum(coalesce(PRECIP_24h, 0))) %>%
    ungroup()
}

# ------------------------
# 4. Safe first helper
# ------------------------
safe_first <- function(x) {
  x <- na.omit(x)
  if(length(x) == 0) NA_real_ else x[1]
}

# ------------------------
# 5. Daily aggregation
# ------------------------
aggregate_daily <- function(wx) {
  wx %>%
    group_by(STATION, DATE) %>%
    summarise(
      YEAR  = first(YEAR),
      MONTH = first(MONTH),
      DAY   = first(DAY),
      across(c(FFMC, ISI, DMC, DC, BUI, FWI, PRECIP_24h), safe_first),
      ACCUM_PRECIP = max(ACCUM_PRECIP, na.rm = TRUE),
      .groups = "drop"
    )
}

# ------------------------
# 6. Monthly percentiles
# ------------------------
calc_monthly_percentiles <- function(wx, variables, percentiles) {
  
  missing_vars <- setdiff(variables, names(wx))
  if (length(missing_vars) > 0) {
    stop(
      paste("Missing variables:", paste(missing_vars, collapse = ", "))
    )
  }
  
  probs <- unname(percentiles)
  pct_names <- names(percentiles)
  
  wx %>%
    dplyr::group_by(STATION, MONTH) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(variables),
        ~ list(stats::quantile(.x, probs = probs, na.rm = TRUE)),
        .names = "{.col}"
      ),
      .groups = "drop"
    ) %>%
    tidyr::unnest(cols = dplyr::all_of(variables)) %>%
    dplyr::mutate(
      PCTL = factor(
        rep(pct_names, length.out = dplyr::n()),
        levels = pct_names
      )
    )
}

# ------------------------
# 7. Exceedance statistics
# ------------------------
calc_exceedance_stats <- function(data, ..., percentiles) {
  
  conditions <- rlang::enquos(...)
  
  if (length(conditions) == 0)
    stop("Supply at least one condition, e.g. ISI >= 7.6")
  
  if (is.null(names(percentiles)))
    stop("percentiles must be a NAMED numeric vector")
  
  monthly_counts <- data %>%
    dplyr::mutate(
      EXCEED = Reduce(`&`, lapply(conditions, rlang::eval_tidy, data = data))
    ) %>%
    dplyr::group_by(STATION, YEAR, MONTH) %>%
    dplyr::summarise(
      N_EXCEED = sum(EXCEED, na.rm = TRUE),
      .groups = "drop"
    )
  
  monthly_counts %>%
    dplyr::group_by(STATION, MONTH) %>%
    dplyr::summarise(
      MEAN_EXCEED = mean(N_EXCEED, na.rm = TRUE),
      
      !!!purrr::imap(
        percentiles,
        ~ rlang::expr(quantile(N_EXCEED, !!.x, na.rm = TRUE))
      ),
      
      .groups = "drop"
    )
}

# ------------------------
# 8. Average wx on over threshold days
# ------------------------
calc_monthly_weather_for_thresholds <- function(wx_hourly, thresholds) {
  all_results <- tibble()
  
  for (threshold_expr in thresholds) {
    thresh_name <- as_label(threshold_expr)
    
    daily_flags <- wx_hourly %>%
      group_by(STATION, DATE, YEAR, MONTH) %>%
      summarise(
        across(c(ISI, BUI, FFMC, DMC, DC, FWI), ~ if(all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      mutate(EXCEED = eval_tidy(threshold_expr))
    
    exceed_days <- daily_flags %>% filter(EXCEED) %>% select(STATION, DATE, YEAR, MONTH)
    
    if(nrow(exceed_days)==0){
      placeholder <- wx_hourly %>%
        distinct(STATION, MONTH) %>%
        mutate(
          mean_max_temp = NA_real_, mean_min_temp = NA_real_, mean_min_rh = NA_real_,
          mean_precip = NA_real_, mean_min_ws = NA_real_, mean_max_ws = NA_real_,
          n_days = 0, threshold = thresh_name
        )
      all_results <- bind_rows(all_results, placeholder)
      next
    }
    
    wx_exceed <- wx_hourly %>% inner_join(exceed_days, by = c("STATION","DATE","YEAR","MONTH"))
    
    daily_weather <- wx_exceed %>%
      group_by(STATION, DATE, YEAR, MONTH) %>%
      summarise(
        max_temp = if(all(is.na(h_TEMP))) NA_real_ else max(h_TEMP, na.rm = TRUE),
        min_temp = if(all(is.na(h_TEMP))) NA_real_ else min(h_TEMP, na.rm = TRUE),
        min_rh   = if(all(is.na(h_RH))) NA_real_ else min(h_RH, na.rm = TRUE),
        precip   = if(all(is.na(PRECIP_24h))) NA_real_ else max(PRECIP_24h, na.rm = TRUE),
        min_ws   = if(all(is.na(h_WS))) NA_real_ else min(h_WS, na.rm = TRUE),
        max_ws   = if(all(is.na(h_WS))) NA_real_ else max(h_WS, na.rm = TRUE),
        .groups = "drop"
      )
    
    monthly_weather <- daily_weather %>%
      group_by(STATION, MONTH) %>%
      summarise(
        mean_max_temp = mean(max_temp, na.rm = TRUE),
        mean_min_temp = mean(min_temp, na.rm = TRUE),
        mean_min_rh   = mean(min_rh, na.rm = TRUE),
        mean_precip   = mean(precip, na.rm = TRUE),
        mean_min_ws   = mean(min_ws, na.rm = TRUE),
        mean_max_ws   = mean(max_ws, na.rm = TRUE),
        n_days = n(),
        .groups = "drop"
      ) %>%
      mutate(threshold = thresh_name)
    
    all_results <- bind_rows(all_results, monthly_weather)
  }
  
  return(all_results)
}

# ------------------------
# 9. Percentile ribbon plot
# ------------------------
plot_percentile_ribbons_station <- function(data, variable, year_line = NULL, station_name = NULL) {
  var <- enquo(variable)
  var_name <- as_label(var)
  
  df <- data %>% mutate(DOY = yday(DATE))
  percentiles_df <- df %>%
    group_by(DOY) %>%
    summarize(
      p10 = quantile(!!var, 0.10, na.rm=TRUE),
      p25 = quantile(!!var, 0.25, na.rm=TRUE),
      p50 = quantile(!!var, 0.50, na.rm=TRUE),
      p75 = quantile(!!var, 0.75, na.rm=TRUE),
      p90 = quantile(!!var, 0.90, na.rm=TRUE),
      p95 = quantile(!!var, 0.95, na.rm=TRUE),
      .groups="drop"
    ) %>% mutate(DATE = as.Date(DOY-1, origin="2000-01-01"))
  
  year_data <- df %>%
    filter(YEAR == year_line) %>%
    select(DOY, value = !!var) %>%
    mutate(DATE = as.Date(DOY-1, origin="2000-01-01"))
  
  fill_colors <- c("Min-Q10"="#d9ef8b","Q25-Q50"="#ffffbf","Q50-Q75"="#fee08b","Q75-Q90"="#fc8d59","Q90-Max"="#d73027")
  line_color <- setNames("black", paste(year_line,"Current"))
  
  ggplot() +
    geom_ribbon(data=percentiles_df, aes(x=DATE, ymin=p10, ymax=p25, fill="Min-Q10"), alpha=0.5) +
    geom_ribbon(data=percentiles_df, aes(x=DATE, ymin=p25, ymax=p50, fill="Q25-Q50"), alpha=0.5) +
    geom_ribbon(data=percentiles_df, aes(x=DATE, ymin=p50, ymax=p75, fill="Q50-Q75"), alpha=0.5) +
    geom_ribbon(data=percentiles_df, aes(x=DATE, ymin=p75, ymax=p90, fill="Q75-Q90"), alpha=0.5) +
    geom_ribbon(data=percentiles_df, aes(x=DATE, ymin=p90, ymax=p95, fill="Q90-Max"), alpha=0.5) +
    geom_line(data=year_data, aes(x=DATE, y=value, color=paste(year_line,"Current")), size=1) +
    scale_fill_manual(name="Percentiles", values=fill_colors, guide=guide_legend(reverse=TRUE)) +
    scale_color_manual(name=NULL, values=line_color) +
    scale_x_date(date_labels="%b", date_breaks="1 month") +
    labs(x="Date", y=var_name,
         title = ifelse(is.null(station_name),
                        paste("Percentile ribbons for", var_name),
                        paste("Percentile ribbons for", var_name, "at", station_name))) +
    theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1))
}

# ------------------------
# 10. Plot weather events
# ------------------------
plot_weather_events <- function(data, ..., start_month = 1, end_month = 12) {
  conditions <- enquos(...)
  if(length(conditions) == 0) stop("Provide at least one condition, e.g. ISI >= 8")
  
  df <- data %>%
    mutate(
      DATE_ONLY = as.Date(DATE),
      DAY_MONTH = as.Date(format(DATE, paste0("2000-%m-%d"))),
      MONTH = as.integer(format(DATE, "%m"))
    ) %>%
    filter(MONTH >= start_month & MONTH <= end_month)
  
  for(i in seq_along(conditions)) df[[paste0("cond_", i)]] <- eval_tidy(conditions[[i]], data = df)
  cond_cols <- paste0("cond_", seq_along(conditions))
  
  df <- df %>%
    rowwise() %>%
    mutate(
      EVENT_LEVEL = { vals <- c_across(all_of(cond_cols)); if(!any(vals, na.rm=TRUE)) NA_integer_ else max(which(vals)) }
    ) %>%
    ungroup()
  
  labels <- sapply(conditions, function(x) as_label(get_expr(x)))
  df <- df %>% mutate(EVENT = factor(EVENT_LEVEL, levels=seq_along(labels), labels=labels))
  
  n <- length(labels)
  cols <- colorRampPalette(c("#f4a300","#ff69b4"))(n)
  
  station_name <- unique(df$STATION)
  
  ggplot(df, aes(x = DAY_MONTH, y = factor(YEAR), fill = EVENT)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_x_date(date_labels = "%b-%d", date_breaks = "1 month") +
    scale_fill_manual(values = cols, na.value = "grey90") +
    labs(x="Date", y="Year", fill="Condition",
         title = paste("Days over Threshold for", station_name,
                       paste0("(Months ", start_month, "-", end_month, ")"))) +
    theme_minimal() +
    theme(panel.grid=element_blank(), axis.text.x = element_text(angle=45,hjust=1))
}

# ------------------------
# 11. Generate wind roses
# ------------------------
generate_windroses <- function(wx_data) {
  stations <- unique(wx_data$STATION)
  plots <- list()
  
  for(station in stations) {
    station_data <- filter(wx_data, STATION == station)
    
    monthly_ws_list <- lapply(1:12, function(m) {
      month_data <- filter(station_data, MONTH == m)
      if(nrow(month_data)==0) return(grid::nullGrob())
      
      p <- windRose(
        month_data,
        ws = "h_WS", wd = "h_WD",
        breaks = c(0,5,10,15,20,30,40),
        key.header = "Wind Speed", key.footer="(km/h)",
        paddle=FALSE, plot=FALSE,
        par.settings=list(fontsize=list(text=6, axis=6))
      )$plot
      p$main <- paste("Month", m)
      p
    })
    
    monthly_ws_grid <- gridExtra::arrangeGrob(grobs=monthly_ws_list, ncol=4)
    plots[[station]] <- monthly_ws_grid
  }
  
  return(plots)
}

# ------------------------
# 12. Update from Git
# ------------------------
auto_update <- function() {
  
  if (!dir.exists(".git")) return(invisible(FALSE))
  
  message("Checking for updates from GitHub...")
  
  # Check for local changes
  status <- system("git status --porcelain", intern = TRUE)
  
  if (length(status) > 0) {
    
    message("Local changes detected.")
    response <- readline(
      prompt = "Updating may overwrite local edits. Continue? (y/n): "
    )
    
    if (tolower(response) != "y") {
      message("Update cancelled.")
      return(invisible(FALSE))
    }
  }
  
  # Perform update
  system("git fetch origin main", intern = TRUE)
  system("git reset --hard origin/main", intern = TRUE)
  
  message("Auto-update complete.")
  invisible(TRUE)
}

# Run when sourced
auto_update()
