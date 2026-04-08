.datatable.aware <- TRUE

#' Inverse Gaussian random variable (ADMB style)
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
    type  = rec_params$rec_type, # "mean", "BH", or "Ricker"
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

  # calculate average (includes the new F for last year)
  F_avg = mean(c(hist_Fs, last_F))

  f_out = 0

  if (scenario == 1) {
    # alt 1: Max ABC (F40)
    f_out = get_tier_f(ssb, f_ref = F40, b_ref = B40)

  } else if (scenario == 2) {
    # alt 2: Author ABC (Defaults to F40)
    f_out = get_tier_f(ssb, f_ref = F40, b_ref = B40)

  } else if (scenario == 3) {
    # alt 3: 50% Max F
    f_out = get_tier_f(ssb, f_ref = F40, b_ref = B40) * 0.5

  } else if (scenario == 4) {
    # alt 4: 5-Year Average
    f_out = F_avg

  } else if (scenario == 5) {
    # alt 5: No Fishing
    f_out = 0

  } else if (scenario == 6) {
    # alt 6: OFL (F35)
    f_out = get_tier_f(ssb, f_ref = F35, b_ref = B35)

  } else if (scenario == 7) {
    # alt 7: Approaching Overfished
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
#' @param scenarios the NPFMC scenarios to run (1:7)
#' @export
run_projections <- function(report, future_catch = NULL,
                            yield_ratio = NULL, n_sims = 1000,
                            n_years = 14, unit_conversion = 1, sex_ratio = 0.5,
                            rec_model = "IG", sigma_r_override = NULL,
                            scenarios = 1:7) {

  # setup
  spawn_frac = report$spawn_fract
  n_ages = length(report$waa)
  last_yr = length(report$years)
  m_vec = rep(report$M, n_ages)
  mat_vec = report$maa
  wt_vec = report$waa
  n_start = report$Nat[, last_yr]
  sel_vec = report$slx[, 1] # flag!!! won't work with many assessments...
  mat_wt_spawn_vec = sex_ratio * wt_vec * mat_vec * exp(-m_vec * spawn_frac)

  year_out_template = report$years[last_yr] + (1:n_years) - 1

  # recruitment parameters (mean)
  hist_rec_val = report$recruits[report$years %in% (1977 + report$ages[[1]]):(tail(report$years, 1) - report$ages[[1]])]
  a_mean = mean(hist_rec_val, na.rm = TRUE)
  h_mean = 1 / mean(1 / hist_rec_val, na.rm = TRUE)

  rec_params <- list(
    rec_type = rec_model,
    r_mean = a_mean,
    h_mean = h_mean,
    b0 = if(!is.null(report$B0)) report$B0 else NA,
    r0 = exp(report$log_mean_R),
    h = if(!is.null(report$steepness)) report$steepness else 1.0, # Removed double comma
    sigma_r = if(!is.null(sigma_r_override)) sigma_r_override else report$sigmaR,
    rho = if(!is.null(report$rho)) report$rho else 0
  )

  if (rec_model %in% c("Ricker", "BH") && (is.na(rec_params$b0) || is.na(rec_params$r0))) {
    warning("Recruitment model is Ricker/BH but B0 or R0 could not be found in report.")
  }

  # catch vector setup
  catch_vec_std = rep(NA, n_years)
  if (!is.null(future_catch)) catch_vec_std[1] = future_catch[1]

  catch_vec_author = rep(NA, n_years)
  if (!is.null(future_catch)) {
    len = min(length(future_catch), n_years)
    catch_vec_author[1:len] = future_catch[1:len]
  }

  # main loop ----
  results_list = list()
  s1_means = NULL
  run_order = sort(unique(scenarios))

  for(scen in run_order) {

    if (scen == 2 && !is.null(yield_ratio)) {
      if (is.null(s1_means)) {
        stop("Scenario 2 with yield_ratio requires Scenario 1 to be run simultaneously.")
      }
      catch_vec_author[2] = s1_means[2] * yield_ratio
      catch_vec_author[3] = s1_means[3] * yield_ratio
    }

    fixed_catch_vec = if(scen == 2) catch_vec_author else catch_vec_std

    # run sims
    sims <- lapply(1:n_sims, function(id) {
      n_curr = n_start
      chi_curr = 0
      f_year_1 = NULL

      ssb_out = catch_out = rec_out = f_out = numeric(n_years)

      for(y in 1:n_years) {
        ssb_curr = sum(n_curr * mat_wt_spawn_vec) * unit_conversion

        target_c = fixed_catch_vec[y]
        if (!is.na(target_c)) {
          expl_bio = sum(n_curr * wt_vec * sel_vec) * unit_conversion
          f_val = if(target_c > expl_bio) 5.0 else
            get_f_from_catch(target_c, n_curr, m_vec, sel_vec, wt_vec, unit_conversion)
        } else {
          f_val = get_scenario_f(scen, ssb_curr, y, report, current_F_sys = f_year_1)
        }

        if (y == 1) f_year_1 = f_val

        step_res = project_step(n_curr, m_vec, sel_vec, wt_vec, mat_vec, f_val, rec_params, chi_curr, spawn_frac)

        n_curr = step_res$n_at_age
        chi_curr = step_res$chi
        ssb_out[y] = step_res$ssb * unit_conversion * sex_ratio
        catch_out[y] = step_res$catch_biomass * unit_conversion
        rec_out[y] = n_curr[1]
        f_out[y] = f_val
      }
      # closes the 'y' for loop and returns the data.table for one simulation
      return(data.table::data.table(scenario = scen, sim_id = id, year = year_out_template,
                                    ssb = ssb_out, catch = catch_out,
                                    recruitment = rec_out, f = f_out))
    }) # close the lapply

    scen_dt <- data.table::rbindlist(sims)
    results_list[[as.character(scen)]] <- scen_dt

    # store means if this was scenario 1 to feed scenario 2
    if (scen == 1) {
      s1_means = scen_dt[, .(m = mean(catch)), by = year]$m
    }
  } # close the for loop

  return(data.table::rbindlist(results_list))
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

#' Generate Executive Summary and Projections
#' @param report The RTMB report object
#' @param year current assessment year
#' @param species species name string
#' @param region region name string
#' @param future_catch numeric vector of catches for the next 2 years
#' @param yield_ratio ratio to scale maxABC to Author ABC (e.g. 0.8)
#' @param output_dir where to save the CSVs
proj_rtmb <- function(report, year, species, region,
                              future_catch, yield_ratio,
                              output_dir = "processed") {

  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # run projections (runs scenarios 1-7)
  # use unit_conversion=1 assuming report and catch units match (e.g., Tons)
  raw_proj <- run_projections(
    report = report,
    future_catch = future_catch,
    yield_ratio = yield_ratio,
    n_sims = 1000,
    n_years = 14,
    rec_model = "IG"
  )

  # format and save standard scenario files
  ssb_tab = format_output(raw_proj, "ssb")
  f_tab   = format_output(raw_proj, "f")
  yld_tab = format_output(raw_proj, "catch")

  vroom::vroom_write(ssb_tab, file.path(output_dir, "mscen_ssb.csv"), ",")
  vroom::vroom_write(f_tab,   file.path(output_dir, "mscen_f.csv"), ",")
  vroom::vroom_write(yld_tab, file.path(output_dir, "mscen_yld.csv"), ",")

  # build executive summary table
  # pulling reference points directly from the report object
  y1 <- year + 1
  y2 <- year + 2

  # get year 1 and year 2 means for Scenario 1 (Max ABC) and Scenario 6 (OFL)
  summary_stats = raw_proj %>%
    tidytable::filter(year %in% c(y1, y2)) %>%
    tidytable::summarise(
      ssb = mean(ssb),
      abc = mean(catch),
      f   = mean(f),
      .by = c(year, scenario)
    )

  # xtract specific rows for the table
  # scenario 1 = Max ABC/F40, scenario 6 = OFL/F35
  exec_df = data.frame(
    item = c(
      "M (natural mortality)",
      "Tier",
      "Projected total biomass (t)",
      "Projected female spawning biomass (t)",
      "B100%", "B40%", "B35%",
      "FOFL", "maxFABC", "FABC",
      "OFL (t)", "maxABC (t)", "ABC (t)"
    ),
    y1_val = c(
      report$M,
      ifelse(summary_stats$ssb[summary_stats$year==y1 & summary_stats$scenario==1] > report$B40, "3a", "3b"),
      NA, # Total biomass calculation depends on your report structure
      summary_stats$ssb[summary_stats$year==y1 & summary_stats$scenario==1],
      report$B0, report$B40, report$B35,
      summary_stats$f[summary_stats$year==y1 & summary_stats$scenario==6],
      summary_stats$f[summary_stats$year==y1 & summary_stats$scenario==1],
      summary_stats$f[summary_stats$year==y1 & summary_stats$scenario==2],
      summary_stats$abc[summary_stats$year==y1 & summary_stats$scenario==6],
      summary_stats$abc[summary_stats$year==y1 & summary_stats$scenario==1],
      summary_stats$abc[summary_stats$year==y1 & summary_stats$scenario==2]
    ),
    y2_val = c(
      report$M,
      ifelse(summary_stats$ssb[summary_stats$year==y2 & summary_stats$scenario==1] > report$B40, "3a", "3b"),
      NA,
      summary_stats$ssb[summary_stats$year==y2 & summary_stats$scenario==1],
      report$B0, report$B40, report$B35,
      summary_stats$f[summary_stats$year==y2 & summary_stats$scenario==6],
      summary_stats$f[summary_stats$year==y2 & summary_stats$scenario==1],
      summary_stats$f[summary_stats$year==y2 & summary_stats$scenario==2],
      summary_stats$abc[summary_stats$year==y2 & summary_stats$scenario==6],
      summary_stats$abc[summary_stats$year==y2 & summary_stats$scenario==1],
      summary_stats$abc[summary_stats$year==y2 & summary_stats$scenario==2]
    )
  )

  # formatting
  colnames(exec_df) = c("item", as.character(y1), as.character(y2))
  vroom::vroom_write(exec_df, file.path(output_dir, "exec_summ.csv"), ",")

  return(exec_df)
}

# murky waters...
#' Apply TAC-ABC fitting logic
#' @param abc vector of calculated ABCs for the species in the complex
#' @param tacpar the list object from read_tacpar()
#' @param total_oy_cap The OY cap (e.g., 2000 for 2 million t)
apply_tac <- function(abc, tacpar, total_oy_cap = 1945) {
  # scale ABCs
  #in ADMB: abctmp = agg_abc(itacspp) / maxabc(itacspp)
  scaled_abc = abc / tacpar$maxabc

  # find nodes (indexing)
  # in ADMB: ijunk = min(nnodes, int(abctmp * nnodes))
  # note: R uses 1-based indexing, ADMB used 0-based
  node_idx = pmin(tacpar$nnodes, floor(scaled_abc * tacpar$nnodes)) + 1

  # compute aggregate TACs
  # in ADMB: agg_tac(itacspp) = abctmp * mfexp(theta(ijunk,itacspp))
  # extract the specific theta value for each species' current node
  weights = sapply(seq_along(abc), function(i) tacpar$theta[node_idx[i], i])
  agg_tac = scaled_abc * exp(weights)

  # constrain to OY Cap
  # in ADMB: agg_tac /= sum(agg_tac); agg_tac *= 1945.;
  final_tac = (agg_tac / sum(agg_tac)) * total_oy_cap

  return(final_tac)
}


# possible setup for adjusting
run_projections <- function(report,
                            future_catch = NULL,
                            yield_ratio = NULL,
                            tacpar = NULL,  # Add this here
                            ...) {

 #...

  for(scen in run_order) {

    if (scen == 2 && !is.null(tac_params)) {

      # 1. get current Max ABC catch for all species in the complex
      # assumes you are projecting the whole complex at once
      # abc_current <- [Calculated Catch at F40]

      # apply the fitted TAC logic
      # catch_vec_author <- apply_tac(abc_current, tac_params)
    }

    # more code...
  }
}


