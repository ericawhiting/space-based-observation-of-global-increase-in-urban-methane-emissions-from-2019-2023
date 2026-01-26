calc_dry_column_mixing_ratio <- function(surf_pres, gas_total_column, gas_total_column_prec, h2o_total_column) {
  # calculate dry column mixing ratio for (used in TROPOMI CH4:CO project for CO data)
  mair <- 28.96  #molar mass of air, g/mol
  mh2o <- 18.02  #molar mass of water vapor, g/mol
  g <- 9.80665   # gravitational acceleration, m/s^2

  # calculate scaling factor from surface pressure and h2o column density
  # using input surface pressure, alternatively could input and index (e.g pressure_levels[50,,])
  dry_scaling <- (surf_pres / (mair * g) * 1000  - h2o_total_column * mh2o / mair)^-1


  gas_total_column_dry <- gas_total_column * dry_scaling * 1e9 # X_gas (ppbv)
  gas_total_column_prec_dry <- gas_total_column_prec * dry_scaling * 1e9 # (ppbv)
  outputs <- list(gas_total_column_dry, gas_total_column_prec_dry)
  return(outputs)
}