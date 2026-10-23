# Chargement des librairies
library(nasapower)
library(dplyr)
library(readr)
library(lubridate)
library(purrr)
library(janitor)
library(daymetr)
library(ggplot2)
library(slider)

# 1. Import des sites
sites <- read.csv("Data/Meteo/sites_meteo.csv", sep=";") %>%
  mutate(
    start_year = year(ymd(start)),
    end_year   = year(ymd(end))
  )
str(sites)

# 2. Données Daymet pour chaque site
get_daymet <- function(site, lat, lon, start, end) {
  download_daymet(
    site = site,
    lat = lat,
    lon = lon,
    start = start,
    end = end,
    internal = TRUE
  )$data %>%
    mutate(Site = site)
}

daymet <- purrr::pmap_dfr(
  list(
    sites$Site,
    sites$LAT,
    sites$LON,
    sites$start_year,
    sites$end_year
  ),
  get_daymet
) %>%
  rename(
    YEAR = year,
    DOY = yday,
    Day_length = dayl..s.,
    T_Max  = tmax..deg.c.,
    T_Min  = tmin..deg.c.,
    Precipitation   = prcp..mm.day.,
    Solar_radiation   = srad..W.m.2.,
    SWE    = swe..kg.m.2.,
    Vapor_pressure = vp..Pa.
  ) %>%
  mutate(DOY = DOY + 1)

# 3. Données NASA POWER pour chaque site
vars <- c("WS2M", "WD2M")

get_power_safe <- function(site, lon, lat, start, end) {
  tryCatch({
    pw <- nasapower::get_power(
      community = "AG",
      lonlat = c(lon, lat),
      pars = vars,
      dates = c(as.character(start), as.character(end)),
      temporal_api = "DAILY"
    )
    if (!"YYYYMMDD" %in% names(pw)) stop("Données manquantes")
    out <- pw %>%
      mutate(Site = site, Date = as.Date(YYYYMMDD)) %>%
      select(Site, Date, any_of(vars))
    missing_vars <- setdiff(vars, names(out))
    if(length(missing_vars)>0) out[missing_vars] <- NA
    out
  }, error = function(e){
    reslist <- lapply(vars, function(param) {
      pw <- tryCatch(
        nasapower::get_power("AG", lonlat = c(lon, lat), pars = param, dates = c(as.character(start), as.character(end)), temporal_api = "DAILY"),
        error = function(e) NULL
      )
      if (!is.null(pw) && "YYYYMMDD" %in% names(pw) && param %in% names(pw)) {
        data.frame(
          Site = site,
          Date = as.Date(pw$YYYYMMDD),
          setNames(list(pw[[param]]), param)
        )
      } else NULL
    })
    valid <- reslist[!sapply(reslist, is.null)]
    if (length(valid) > 0) {
      comb <- purrr::reduce(valid, full_join, by = c("Site", "Date"))
      missing_vars <- setdiff(vars, names(comb))
      if(length(missing_vars)>0) comb[missing_vars] <- NA
      comb
    } else {
      data.frame(Site = site, Date = as.Date(character()), 
                 WS2M = NA_real_, WD2M = NA_real_)
    }
  })
}

nasa_power <- purrr::pmap_dfr(
  list(sites$Site, sites$LON, sites$LAT, sites$start, sites$end),
  ~{
    message(sprintf("Site: %s", .x))
    get_power_safe(.x, .y, ..3, ..4, ..5)
  }
) %>%
  mutate(DOY = yday(Date-1)) %>%
  as.data.frame()

# 4. Fusion des deux jeux de données
Meteo <- daymet %>%
  left_join(nasa_power, by = c("Site", "DOY"))

# 5. Création des variables additionnelles
Meteo <- Meteo %>%
  arrange(Site, DOY) %>%
  group_by(Site) %>%
  mutate(
    Mean_temp = (T_Max + T_Min) / 2,                     
    Degree_days = pmax(T_Min, 0),                     
    Mean_temp_3 = slide_dbl(Mean_temp, mean, .before = 3, .after = 0, .complete = TRUE),
    Mean_temp_5 = slide_dbl(Mean_temp, mean, .before = 5, .after = 0, .complete = TRUE),
    Mean_temp_7 = slide_dbl(Mean_temp, mean, .before = 7, .after = 0, .complete = TRUE),
    Cum_rain_3 = slide_dbl(Precipitation, sum, .before = 3, .after = -1, .complete = TRUE),
    Cum_rain_5 = slide_dbl(Precipitation, sum, .before = 5, .after = -1, .complete = TRUE),
    Cum_rain_7 = slide_dbl(Precipitation, sum, .before = 7, .after = -1, .complete = TRUE),
    Relative_humidity = 100 * (Vapor_pressure/100) / (6.1078 * exp((17.27 * Mean_temp) / (Mean_temp + 237.3)))
  ) %>%
  rename(
    Wind_direction = WD2M,
    Wind_speed = WS2M
  ) %>%
  as.data.frame() %>%
  ungroup()

# 6. Visualisation et export
glimpse(Meteo)
head(Meteo)
Meteo %>% 
  ggplot(aes(x = DOY, y = Cum_rain_7, color = Site)) +
  geom_line() +
  facet_wrap(. ~ YEAR) +
  scale_color_manual(values = c(
    "black", "red", "blue", "darkgrey", "#6E8E84FF", '#CAC27EFF', "#1A477FFF",
    '#55752FFF','#7B92A8FF', '#2D6D66FF',  '#E37E00FF', 
    '#FFD200FF', '#A0527DFF'
  ))

write.csv(Meteo, "Data/Meteo/Meteo.csv")

