.datatable.aware <- TRUE

#' Inverse Gaussian Random Variable (ADMB style)
#'
#' @param n number of deviates to generate
#' @param arithmetic_mean the standard mean of historical recruitment
#' @param harmonic_mean the harmonic mean of historical recruitment
#' @return vector of recruitment deviates
#' @export
inverse_gaussian <- function(n, arithmetic_mean, harmonic_mean) {

  # calculate shape pParameters (from ADMB LOCAL_CALCS)
  # gamma = AMean / HMean
  # delta = 1 / (gamma - 1)
  gamma_param = arithmetic_mean / harmonic_mean
  delta = 1 / (gamma_param - 1)
  beta = arithmetic_mean
  z = rnorm(n)
  u = runif(n)
  # apply Michaelies-Schwartz-Tweedie algorithm
  psi = z^2
  term = sqrt(4 * delta * psi + psi^2)
  omega = beta * (1 + (psi - term) / (2 * delta))
  zeta = beta * (1 + (psi + term) / (2 * delta))
  gtheta = beta / (beta + omega)
  rec_out = ifelse(u <= gtheta, omega, zeta)
  return(rec_out)
}

#' Calculate recruitment (supports SR models or mean recruitment)
#'
#' @param ssb spawning stock biomass (ignored if type = "Mean")
#' @param b0 unfished equilibrium biomass (ignored if type = "Mean")
#' @param h steepness (ignored if type = "Mean")
#' @param r0 unfished recruitment (used as parameter for SR curves)
#' @param r_mean average recruitment (used ONLY if type = "Mean")
#' @param h_mean harmonic mean recruitment (used ONLY if type = "IG" (Inverse Gaussian))
#' @param type "BH", "Ricker", "Mean", "or "IG"
#' @param rho autocorrelation coefficient
#' @param sigma_r log-scale standard deviation of recruitment
#' @param epsilon_prev the total log-deviation from the previous year
#' @param rec_dev current year random normal deviate
#'
#' @return a list containing calculated recruitment and the error term (chi)
#' @export
get_recruitment <- function(ssb, b0, h, r0, r_mean = NULL, h_mean = NULL,
                            type = "Mean",
                            rho = 0, sigma_r = 0.5,
                            epsilon_prev = 0, rec_dev = rnorm(1)) {

  # inverse gaussian to be like ADMB
  if (type == "IG") {
    if(is.null(r_mean) || is.null(h_mean)) {
      stop("Must provide r_mean and h_mean for Inverse-Gaussian.")
    }

    # generate deviate using the exact ADMB algorithm
    # note: do not use rho or sigma_r here, as ADMB Rec_Gen=1 doesn't.
    rec_val = inverse_gaussian(1, r_mean, h_mean)

    return(list(rec = rec_val, chi = 0)) # No chi/epsilon tracking needed
  }

  # lognormal
  rec_exp = 0
  if (type == "Mean") {
    rec_exp = r_mean
  } else if (type == "BH") {
    top = 4 * h * r0 * ssb
    bot = b0 * (1 - h) + ssb * (5 * h - 1)
    rec_exp = top / bot
  } else if (type == "Ricker") {
    rec_exp = (ssb / b0) * r0 * exp(h * (1 - ssb / b0))
  }

  # correlated log-normal error
  chi = rho * epsilon_prev + sqrt(1 - rho^2) * rec_dev * sigma_r
  rec_out = rec_exp * exp(chi - (sigma_r^2)/2)

  return(list(rec = rec_out, chi = chi))
}


#' Calculate F based on NPFMC Tier System (Amendment 56)
#'
#' @param ssb current spawning biomass
#' @param f_ref reference F (F40, F35, or Fabc)
#' @param b_ref reference biomass (B40, B35, or B100 depending on tier)
#' @param alpha slope cut-off (usually 0.05)
#' @param ssl_protection boolean. if TRUE, F cuts to 0 at higher B/Bref ratio (0.2)
#'
#' @return the target Fishing Mortality (F)
get_tier_f <- function(ssb, f_ref, b_ref, alpha = 0.05, ssl_protection = FALSE) {

  # ratio of current biomass to reference (B40 or Bmsy)
  ratio = ssb / b_ref

  if (ssl_protection) {
    # Steller Sea Lion specific rule (harder cutoff)
    if (ssb < 0.2 * b_ref) {
      f_out = 0
    } else if (ssb >= 0.2 * b_ref & ssb < b_ref) {
      # linear ramp
      f_out = f_ref * (ssb/b_ref - alpha) / (1 - alpha)
    } else {
      f_out = f_ref
    }
  } else {
    # standard Amendment 56
    if (ssb < alpha * b_ref) {
      f_out = 0
    } else if (ssb >= alpha * b_ref & ssb < b_ref) {
      # linear ramp
      f_out = f_ref * (ssb/b_ref - alpha) / (1 - alpha)
    } else {
      # above target
      f_out = f_ref
    }
  }

  return(pmax(0, f_out)) # ensure non-negative
}

#' Solve for F given a Target Catch
get_f_from_catch <- function(target_catch, n_at_age, m_at_age, sel_at_age, wt_at_age, scalar=1) {
  calc_catch_diff <- function(f_trial) {
    Z = m_at_age + f_trial * sel_at_age
    C_num = n_at_age * (f_trial * sel_at_age / Z) * (1 - exp(-Z))
    # apply scalar to calculated catch to match target units
    C_bio = sum(C_num * wt_at_age) * scalar
    return(C_bio - target_catch)
  }

  tryCatch({
    res <- uniroot(calc_catch_diff, interval = c(0, 5), tol = 1e-5)
    return(res$root)
  }, error = function(e) return(3.0))
}

#' Single Year Population Projection Step
#'
#' @param n_at_age vector of abundance at age at start of year t
#' @param m_at_age vector of M at age
#' @param sel_at_age vector of Selectivity at age
#' @param wt_at_age vector of Weight at age
#' @param mat_at_age vector of Maturity at age
#' @param f_current scalar Fishing mortality for year t
#' @param rec_params list of SR parameters (b0, h, r0, rho, sigma_r)
#' @param prev_chi the error term from year t-1 (for Rho calc)
#' @param spawn_frac fraction of year elapsed before spawning
#'
#' @return A list containing N_at_age (t+1), Catch (t), SSB (t), and new chi
project_step <- function(n_at_age, m_at_age, sel_at_age,
                         wt_at_age, mat_at_age,
                         f_current, rec_params, prev_chi, spawn_frac) {

  ages = length(n_at_age)

  # calculate SSB for current year
  Z_t = m_at_age + f_current * sel_at_age
  ssb_t = sum(n_at_age * wt_at_age * mat_at_age * exp(-spawn_frac * Z_t))

  # catch
  F_at_age = f_current * sel_at_age
  catch_numbers = n_at_age * (F_at_age / Z_t) * (1 - exp(-Z_t))
  catch_biomass = sum(catch_numbers * wt_at_age)

  # survival
  S_t = exp(-Z_t)
  n_next = numeric(ages)

  # ages 2 to max-1
  n_next[2:(ages-1)] <- n_at_age[1:(ages-2)] * S_t[1:(ages-2)]

  # plus group (accumulates survivors from last age and previous age)
  n_next[ages] <- (n_at_age[ages-1] * S_t[ages-1]) + (n_at_age[ages] * S_t[ages])

  # recruitment (age 1)
  rec_out <- get_recruitment(
    ssb = ssb_t,
    # structural parameters (Pass NULL if not using a curve)
    b0 = rec_params$b0,
    h  = rec_params$h,
    r0 = rec_params$r0,
    # arguments for mean recruitment
    r_mean = rec_params$r_mean,   # historical average (required if type="Mean")
    h_mean = rec_params$h_mean,
    type  = rec_params$rec_type, # "Mean", "BH", or "Ricker"
    # error structure
    rho = rec_params$rho,
    sigma_r = rec_params$sigma_r,
    epsilon_prev = prev_chi
  )

  n_next[1] = rec_out$rec

  return(list(
    n_at_age = n_next,
    catch_biomass = catch_biomass,
    ssb = ssb_t,
    chi = rec_out$chi
  ))
}


#' Scenario Logic
get_scenario_f <- function(scenario, ssb, year_idx, report,
                           current_F_sys = NULL) {

  F40 = report$F40
  F35 = report$F35
  B40 = report$B40
  B35 = report$B35

  # 5-year avg F
  term_idx <- length(report$Ft)
  hist_Fs <- report$Ft[(term_idx-4):(term_idx-1)]

  # get F for Current Year
  #    if we calculated a specific F for the fixed catch use it.
  #    otherwise fall back to the report's terminal F.
  last_F = if(!is.null(current_F_sys)) current_F_sys else report$Ft[term_idx]

  # 3. calculate average (includes the new F for last year)
  F_avg = mean(c(hist_Fs, last_F))

  f_out = 0

  if (scenario == 1) {
    # Alt 1: Max ABC (F40)
    f_out = get_tier_f(ssb, f_ref = F40, b_ref = B40)

  } else if (scenario == 2) {
    # Alt 2: Author ABC (Defaults to F40)
    f_out = get_tier_f(ssb, f_ref = F40, b_ref = B40)

  } else if (scenario == 3) {
    # Alt 3: 50% Max F
    f_out = get_tier_f(ssb, f_ref = F40, b_ref = B40) * 0.5

  } else if (scenario == 4) {
    # Alt 4: 5-Year Average (Updated)
    f_out = F_avg

  } else if (scenario == 5) {
    # Alt 5: No Fishing
    f_out = 0

  } else if (scenario == 6) {
    # Alt 6: OFL (F35)
    f_out = get_tier_f(ssb, f_ref = F35, b_ref = B35)

  } else if (scenario == 7) {
    # Alt 7: Approaching Overfished
    if (year_idx <= 3) {
      f_out = get_tier_f(ssb, f_ref = F40, b_ref = B40)
    } else {
      f_out = get_tier_f(ssb, f_ref = F35, b_ref = B35)
    }
  }
  return(f_out)
}

#' Run NPFMC Scenarios 1-7
#'
#' @param report The RTMB report object
#' @param future_catch vector of catches (e.g., c("2024" = 1170, "2025" = 1968))
#' @param yield_ratio A numeric value (e.g., 0.3598).
#'        If provided, Scenario 2 (Author F) catches for years 2 & 3 will be calculated as:
#'        Scenario 1 Catch * yield_ratio.
#' @param n_sims number of simulations
#' @param n_years number of projection years
#' @param unit_conversion Scaling factor to align Biomass units with Catch units.
#'        - If Nat=Millions, Waa=Grams, Catch=Tons -> Use 1.0
#'        - If Nat=Thousands, Waa=Grams, Catch=Tons -> Use 0.001
#' @param sex_ratio proportion of females for SSB calculation (default 0.5)
#' @param rec_model the recruitment type: "Mean", "IG" (mean using Inverse Gaussian),
#'                  "Ricker", or "BH" (Beverton-Holt).
#' @param sigma_r_override optional numeric to override report$sigmaR.
#' @export
run_projections <- function(report, future_catch = NULL,
                            yield_ratio = NULL, n_sims = 1000,
                            n_years = 14, unit_conversion = 1, sex_ratio = 0.5,
                            rec_model = "IG", sigma_r_override = NULL) {

  # setup
  spawn_frac = report$spawn_fract

  n_ages =length(report$waa)
  last_yr_idx = length(report$years)
  m_vec = rep(report$M, n_ages)
  mat_vec = report$maa
  wt_vec = report$waa
  n_start = report$Nat[, last_yr_idx]

  sel_vec <- report$slx[,1] # flag: this won't work with POP!!

  # recruitment parameters (Mean)
  hist_rec_val = report$recruits[report$years %in% (1977 + report$ages[[1]]):(tail(report$years, 1) - report$ages[[1]])]
  a_mean = mean(hist_rec_val, na.rm = TRUE) # arithmetic mean
  h_mean = 1 / mean(1 / hist_rec_val, na.rm = TRUE) # harmonic mean

  b0_val = if(!is.null(report$B0)) report$B0 else NA
  r0_val = exp(report$log_mean_R)
  h_val  = if(!is.null(report$steepness)){
    report$steepness
  } else{
      1.0
    }
  sig_r_val = if(!is.null(sigma_r_override)) {sigma_r_override
  } else{
      report$sigmaR
    }
  rho_val   = if(!is.null(report$rho)){ report$rho
  } else{
      0
    }

  rec_params <- list(
    rec_type = rec_model,
    # Stats for Mean/IG
    r_mean = a_mean,
    h_mean = h_mean,
    # Curve Params for Ricker/BH
    b0 = b0_val,
    r0 = r0_val,
    h = h_val,
    # Error
    sigma_r = sig_r_val,
    rho = rho_val
  )

  if (rec_model %in% c("Ricker", "BH") && (is.na(b0_val) || is.na(r0_val))) {
    warning("Recruitment model is Ricker/BH but B0 or R0 could not be found in report.")
  }

  # catch
  # Scenario 2 gets the full custom vector
  # others get only the first year (current yYear)
  catch_vec_author <- rep(NA, n_years)
  if (!is.null(future_catch)) {
    len = min(length(future_catch), n_years)
    catch_vec_author[1:len] = future_catch[1:len]
  }

  catch_vec_standard = rep(NA, n_years)
  if (!is.null(future_catch)) catch_vec_standard[1] = future_catch[1]

  # scenarios

  if (!is.null(yield_ratio)) {
    # get mean catch from Scenario 1 for years 2 and 3
    s1_catch = scen1_res %>%
      summarize(mean_c = mean(catch), .by = year) %>%
      pull(mean_c)

    # apply ratio: Author Catch = Max Catch * Ratio
    # Year 1 is fixed, Years 2 & 3 are calculated
    catch_vec_author[2] <- s1_catch[2] * yield_ratio
    catch_vec_author[3] <- s1_catch[3] * yield_ratio
    # (Optional: extend further ...)
  }
  tidytable(scenario = 1:7) %>%
    tidytable::mutate(res =  tidytable::map(scenario, function(scen) {

      fixed_catch_vec = if(scen == 2) catch_vec_author else catch_vec_standard

      tidytable::tidytable(sim_id = 1:n_sims) %>%
        mutate(projection_output = tidytable::map(sim_id, function(id) {

          n_curr = n_start
          chi_curr = 0
          output_list = list()
          f_year_1 = NULL
          for(y in 1:n_years) {

            # SSB
            # apply unit_conversion and sex_ratio
            ssb_curr = sum(n_curr * sex_ratio * wt_vec * mat_vec * exp(-m_vec * spawn_frac)) * unit_conversion
            if(is.na(ssb_curr)) ssb_curr = 0

            # determine F
            target_c = fixed_catch_vec[y]
            f_val = 0

            if (!is.na(target_c)) {
              # fixed catch mode
              expl_bio = sum(n_curr * wt_vec * sel_vec) * unit_conversion
              if(target_c > expl_bio) {
                f_val = 5.0
              } else {
                f_val = get_f_from_catch(target_c, n_curr, m_vec, sel_vec, wt_vec, unit_conversion)
              }
            } else {
              # HCR
              f_val = get_scenario_f(scen, ssb_curr, y, report, current_F_sys = f_year_1)
            }

            # project
            step_res <- project_step(
              n_at_age = n_curr, m_at_age = m_vec, sel_at_age = sel_vec,
              wt_at_age = wt_vec, mat_at_age = mat_vec, f_current = f_val,
              rec_params = rec_params, prev_chi = chi_curr, spawn_frac = spawn_frac
            )

            n_curr = step_res$n_at_age
            chi_curr = step_res$chi

            output_list[[y]] <- data.table::data.table(
              year = report$years[last_yr_idx] + y - 1,
              ssb = step_res$ssb * unit_conversion * sex_ratio,
              catch = step_res$catch_biomass * unit_conversion,
              recruitment = n_curr[1],
              f = f_val
            )
          }
          return(data.table::rbindlist(output_list))
        })) %>%
        tidytable::unnest(projection_output)
    })) %>%
    tidytable::unnest(res)
}


#' Format Projection Results (Dynamic Variable)
#'
#' @param projection_data The raw output from run_projections()
#' @param var The variable to summarize: "ssb", "catch", or "f"
#' @return A tidytable matching the ADMB csv format
#' @export
format_output <- function(projection_data, var = "ssb") {
  scen_map = data.table::data.table(scenario = 1:7,
    name = c("maxf", "authf", "half_maxf", "avg5f", "nof", "overf", "appoverf")
  )
  projection_data %>%
    tidytable::left_join(scen_map, by = "scenario") %>%
    tidytable::summarise(mean_val = mean(.data[[var]]), .by = c(year, name)) %>%
    tidytable::pivot_wider(names_from = name, values_from = mean_val) %>%
    tidytable::select(year, maxf, authf, half_maxf, avg5f, nof, overf, appoverf) %>%
    tidytable::mutate(tidytable::across(-year, ~ if(var == "f") round(.x, 4) else round(.x, 1)))
}


