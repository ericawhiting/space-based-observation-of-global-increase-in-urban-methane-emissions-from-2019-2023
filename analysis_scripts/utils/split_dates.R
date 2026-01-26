split_dates <- function(date) {
  # passes in date that has form "yyyymmdd"
  year <- substr(date, 1, 4)
  month <- substr(date, 5, 6)
  day <- substr(date, 7, 8)
  date_string <- paste(year, month, day, sep = "-")
  date_type <- as.Date(date_string)
  return(date_type)
}

weekend_weekday_split_single <- function(date) {
  dayofweek <- weekdays(as.Date(df$date))
  if ((dayofweek > 0) & (dayofweek < 6)) {
    # 0 is sunday, 6 is saturday (1-5 are M-F)
    weekend <- FALSE
  } else {
    weekend <- TRUE
  }
  return(weekend)
}

weekend_weekday_split_df <- function(df) {
  # takes in df and runs each row through split_dates.weekend_weekday_split_single
  df$weekend <- FALSE
  for (i in 1:dim(df)[1]){
    df$weekend[i] <- weekend_weekday_split_single(df$datetimes[i])
  }
  return(df)
}

assign_season_year <- function(year, month) {
  if (month == 1 || month == 2) {
    season <- "DJF"
    season_year <- paste(year, season, sep = "_")
  } else if (month == 12) { # next years winter
    season <- "DJF"
    season_year <- paste(year + 1, season, sep = "_")
  } else if (month %in% c(3, 4, 5)) {
    season <- "MAM"
    season_year <- paste(year, season, sep = "_")
  } else if (month %in% c(6, 7, 8)) {
    season <- "JJA"
    season_year <- paste(year, season, sep = "_")
  } else if (month %in% c(9, 10, 11)) {
    season <- "SON"
    season_year <- paste(year, season, sep = "_")
  }
  return(season_year)
}
