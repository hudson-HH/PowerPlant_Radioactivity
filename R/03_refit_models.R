#!/usr/bin/env Rscript

required_packages(c("data.table", "fixest", "splines"))
suppressPackageStartupMessages({ library(data.table); library(fixest) })
fixest::setFixest_nthreads(max(1L, floor(parallel::detectCores(logical = TRUE) / 2)))
options(fixest_notes = FALSE)

model_file <- file.path(MODEL_DATA_DIR, "filter_period_primary_2001_2024.csv.gz")
if (!file.exists(model_file)) stop("Prepared model matrix is missing: ", model_file, call. = FALSE)
d <- fread(model_file, showProgress = TRUE)
assert_true(nrow(d) == 107407L, "Prepared matrix must contain 107407 formal filter periods.")
assert_true(uniqueN(d$loc_num) == 133L, "Prepared matrix must contain 133 monitors.")

fuels <- c("coal", "gas", "oil", "biomass")
metrics <- c("count", "heat")
radii <- seq(20L, 50L, 5L)

exposure_name <- function(fuel, metric, radius, direction = "upwind") {
  suffix <- if (metric == "count") "n_plants" else "heat_observed"
  paste0("fpmean_campd_", fuel, "_", direction, "_", suffix, "_", radius, "km_90deg")
}

rhs_main <- paste(c("exposure", "exposure:vel", "c_prop", "sunspot_monthly", "pm_model_centered",
  "inv_hpbl", "rhum", "soilm", "vel", "temp_centered", "temp2", "temp3",
  "sunspot_monthly:lat_diff_35", "surface_u_ppm_article:soilm"), collapse = " + ")

# H20B constructs spline bases on the full formal matrix BEFORE fuel/metric
# missingness is applied. Recomputing knots within each model changes results.
ns_vars <- c("c_prop","sunspot_monthly","pm_model_centered","inv_hpbl","rhum","soilm","vel","temp_centered")
ns_terms <- character()
for (vv in ns_vars) {
  basis <- splines::ns(d[[vv]], df=3)
  cols <- paste0("ns_", vv, "_", seq_len(ncol(basis)))
  for (jj in seq_along(cols)) set(d, j=cols[jj], value=basis[,jj])
  ns_terms <- c(ns_terms,cols)
}

fit_one <- function(fuel, metric, radius, natural_spline = FALSE) {
  ex <- exposure_name(fuel, metric, radius)
  x <- d[is.finite(get(ex))]
  x[, exposure := get(ex) / if (metric == "heat") 10000 else 1]
  wind <- median(x$vel, na.rm = TRUE)
  if (!natural_spline) {
    f <- as.formula(paste0("gross_beta_mBq_m3 ~ ", rhs_main, " | monitor_year_fe"))
  } else {
    f <- as.formula(paste0("gross_beta_mBq_m3 ~ exposure + exposure:vel + ",
      paste(c(ns_terms,"sunspot_monthly:lat_diff_35","surface_u_ppm_article:soilm"), collapse=" + "),
      " | monitor_year_fe"))
  }
  fit <- feols(f, data = x, vcov = ~loc_num, notes = FALSE, warn = TRUE)
  b <- coef(fit); v <- vcov(fit)
  beta <- unname(b["exposure"] + wind * b["exposure:vel"])
  var <- v["exposure","exposure"] + wind^2*v["exposure:vel","exposure:vel"] + 2*wind*v["exposure","exposure:vel"]
  se <- sqrt(var)
  data.table(fuel, metric, radius_km = radius, model = if (natural_spline) "natural_spline_df3" else "current",
             N_model = nobs(fit), N_monitors = uniqueN(x$loc_num), wind_reference = "sample_median",
             wind_ms = wind, estimate = beta, se = se, ci95_low = beta - 1.96*se, ci95_high = beta + 1.96*se)
}

# Full Figure 4 and H20B refit: four fuels, two metrics, seven radii,
# and both current and natural-spline covariate forms (112 models).
grid <- CJ(fuel = fuels, metric = metrics, radius = radii)
current <- rbindlist(lapply(seq_len(nrow(grid)), function(i) fit_one(grid$fuel[i], grid$metric[i], grid$radius[i], FALSE)))
bio_grid <- grid
bio_ns <- rbindlist(lapply(seq_len(nrow(bio_grid)), function(i) fit_one(bio_grid$fuel[i], bio_grid$metric[i], bio_grid$radius[i], TRUE)))
refit <- rbind(current, bio_ns, fill = TRUE)
fwrite(refit, file.path(OUTPUT_DIR, "tables", "model_refit_effects.csv"))

frozen <- read_csv(file.path(SOURCE_DIR, "sensitivity", "03_current_and_natural_spline_effects.csv"))[model == "current"]
cmp <- merge(current, frozen[, .(fuel, metric, radius_km, frozen_wind = wind_ms, frozen_estimate = estimate,
                                 frozen_low = ci95_low, frozen_high = ci95_high, frozen_N = N_model)],
             by = c("fuel", "metric", "radius_km"), all.x = TRUE)
cmp[, `:=`(diff_estimate = estimate - frozen_estimate, diff_low = ci95_low - frozen_low,
           diff_high = ci95_high - frozen_high, diff_wind = wind_ms - frozen_wind,
           N_equal = N_model == frozen_N)]
tol <- 1e-7
assert_true(nrow(cmp) == 56L && all(is.finite(cmp$frozen_estimate)) && all(is.finite(cmp$frozen_N)), "Missing or duplicated Figure 4 reference rows.")
cmp[, status := ifelse(pmax(abs(diff_estimate), abs(diff_low), abs(diff_high), abs(diff_wind), na.rm = TRUE) <= tol & N_equal, "PASS", "FAIL")]
fwrite(cmp, file.path(VALIDATION_DIR, "Figure4_refit_vs_frozen.csv"))
assert_true(all(cmp$status == "PASS"), paste0("Figure 4 numerical refit failed; max difference = ",
  max(abs(c(cmp$diff_estimate, cmp$diff_low, cmp$diff_high, cmp$diff_wind)), na.rm = TRUE)))

ns_frozen <- read_csv(file.path(SOURCE_DIR, "sensitivity", "03_current_and_natural_spline_effects.csv"))[
  model == "natural_spline"]
ns_cmp <- merge(bio_ns, ns_frozen[, .(fuel, metric, radius_km, frozen_N = N_model, frozen_wind = wind_ms,
                                      frozen_estimate = estimate, frozen_low = ci95_low, frozen_high = ci95_high)],
                by = c("fuel", "metric", "radius_km"), all.x = TRUE)
ns_cmp[, `:=`(diff_estimate = estimate - frozen_estimate, diff_low = ci95_low - frozen_low,
              diff_high = ci95_high - frozen_high, diff_wind = wind_ms - frozen_wind,
              N_equal = N_model == frozen_N)]
assert_true(nrow(ns_cmp) == 56L && all(is.finite(ns_cmp$frozen_estimate)) && all(is.finite(ns_cmp$frozen_N)), "Missing or duplicated H20B reference rows.")
ns_cmp[, status := ifelse(pmax(abs(diff_estimate), abs(diff_low), abs(diff_high), abs(diff_wind)) <= tol & N_equal, "PASS", "FAIL")]
fwrite(ns_cmp, file.path(VALIDATION_DIR, "H20B_all_fuels_refit_vs_frozen.csv"))
assert_true(all(ns_cmp$status == "PASS"), "H20B natural-spline refit failed numerical equality.")

# Figure 5: joint upwind/downwind models for every fuel/metric/radius.
linear_effect <- function(fit, terms, weights) {
  b <- coef(fit); v <- vcov(fit)
  est <- sum(weights * b[terms]); se <- sqrt(as.numeric(t(weights) %*% v[terms, terms, drop=FALSE] %*% weights))
  data.table(estimate = est, se = se, ci95_low = est - 1.96*se, ci95_high = est + 1.96*se)
}
fit_direction_20 <- function(fuel, metric, radius) {
  up <- exposure_name(fuel, metric, radius, "upwind"); down <- exposure_name(fuel, metric, radius, "downwind")
  x <- d[is.finite(get(up)) & is.finite(get(down))]
  scale <- if (metric == "heat") 10000 else 1
  x[, `:=`(up_exposure = get(up)/scale, down_exposure = get(down)/scale)]
  wind <- median(x$vel, na.rm=TRUE)
  form <- as.formula(paste0("gross_beta_mBq_m3 ~ up_exposure + down_exposure + up_exposure:vel + down_exposure:vel + ",
                            paste(c("c_prop","sunspot_monthly","pm_model_centered","inv_hpbl","rhum","soilm","vel",
                                    "temp_centered","temp2","temp3","sunspot_monthly:lat_diff_35","surface_u_ppm_article:soilm"),collapse=" + "),
                            " | monitor_year_fe"))
  fit <- feols(form, data=x, vcov=~loc_num, notes=FALSE, warn=TRUE)
  up_e <- linear_effect(fit, c("up_exposure","up_exposure:vel"), c(1,wind))
  down_e <- linear_effect(fit, c("down_exposure","down_exposure:vel"), c(1,wind))
  diff_e <- linear_effect(fit, c("up_exposure","up_exposure:vel","down_exposure","down_exposure:vel"), c(1,wind,-1,-wind))
  rbindlist(list(up_e[,component:="upwind"], down_e[,component:="downwind"], diff_e[,component:="upwind_minus_downwind"]))[
    , `:=`(fuel=fuel, metric=metric, radius_km=radius, wind_ms=wind, N_model=nobs(fit))]
}
f5_grid <- grid
f5_refit <- rbindlist(lapply(seq_len(nrow(f5_grid)), function(i) fit_direction_20(f5_grid$fuel[i], f5_grid$metric[i], f5_grid$radius[i])))
f5_frozen <- read_csv(file.path(SOURCE_DIR,"main_figures","Figure5.csv"))
f5_cmp <- merge(f5_refit, f5_frozen[,.(fuel,metric,radius_km,component,frozen_wind=wind_ms,frozen_estimate=estimate,
  frozen_low=ci95_low,frozen_high=ci95_high,n_pre=n_rows)], by=c("fuel","metric","radius_km","component"), all.x=TRUE)
f5_cmp[, `:=`(diff_estimate=estimate-frozen_estimate,diff_low=ci95_low-frozen_low,diff_high=ci95_high-frozen_high,
              diff_wind=wind_ms-frozen_wind)]
assert_true(nrow(f5_cmp)==168L && all(is.finite(f5_cmp$frozen_estimate)), "Missing Figure 5 reference rows.")
f5_cmp[, status:=ifelse(pmax(abs(diff_estimate),abs(diff_low),abs(diff_high),abs(diff_wind),na.rm=TRUE)<=tol & N_model<=n_pre,"PASS","FAIL")]
fwrite(f5_cmp,file.path(VALIDATION_DIR,"Figure5_all_radii_refit_vs_frozen.csv"))
assert_true(all(f5_cmp$status=="PASS"),"Figure 5 joint-direction refit failed numerical equality.")

# H21 complete refit from the distributed compact 2001–2023 GFED matrix.
h21_file <- file.path(MODEL_DATA_DIR,"filter_period_biomass_GFED_2001_2023.csv.gz")
h21d <- fread(h21_file,showProgress=FALSE)
h21d[, `:=`(pm_centered_10=pm_filter_period_official-10,temp_centered=temp-273.15,
            monitor_year_fe=interaction(loc_num,calendar_year,drop=TRUE,lex.order=TRUE))]
h21_rhs <- c("c_prop","sunspot_monthly","pm_centered_10","inv_hpbl","rhum","soilm","vel","temp_centered",
             "I(temp_centered^2)","I(temp_centered^3)","sunspot_monthly:lat_diff_35","surface_u_ppm_article:soilm")
h21_form <- function(fire=FALSE) as.formula(paste0("gross_beta_mBq_m3 ~ exposure_model + exposure_model:vel + ",
  paste(c(h21_rhs,if(fire)"fire_pm25_gfed_centered"),collapse=" + ")," | monitor_year_fe"))
fit_h21 <- function(metric,radius){
  ex<-exposure_name("biomass",metric,radius,"upwind"); vars<-c("gross_beta_mBq_m3","c_prop","sunspot_monthly","pm_centered_10","inv_hpbl","rhum","soilm","vel","temp_centered","lat_diff_35","surface_u_ppm_article","loc_num","monitor_year_fe","fire_pm25_gfed_mean",ex)
  x<-h21d[fire_complete & complete.cases(h21d[,..vars])]
  x[,exposure_model:=get(ex)/if(metric=="heat")10000 else 1]
  center<-mean(x$fire_pm25_gfed_mean); x[,fire_pm25_gfed_centered:=fire_pm25_gfed_mean-center]; wind<-median(x$vel)
  m0<-feols(h21_form(FALSE),data=x,vcov=~loc_num,notes=FALSE,warn=FALSE); m1<-feols(h21_form(TRUE),data=x,vcov=~loc_num,notes=FALSE,warn=FALSE)
  e0<-linear_effect(m0,c("exposure_model","exposure_model:vel"),c(1,wind)); e1<-linear_effect(m1,c("exposure_model","exposure_model:vel"),c(1,wind))
  firese<-sqrt(vcov(m1)["fire_pm25_gfed_centered","fire_pm25_gfed_centered"]); fireb<-coef(m1)["fire_pm25_gfed_centered"]
  data.table(metric,radius_km=radius,n_model_M0=nobs(m0),n_model_M1=nobs(m1),wind_ms=wind,fire_pm25_center=center,
    original_estimate=e0$estimate,original_ci95_low=e0$ci95_low,original_ci95_high=e0$ci95_high,
    fire_adjusted_estimate=e1$estimate,fire_adjusted_ci95_low=e1$ci95_low,fire_adjusted_ci95_high=e1$ci95_high,
    fire_pm25_coef=fireb,fire_pm25_ci95_low=fireb-1.96*firese,fire_pm25_ci95_high=fireb+1.96*firese)
}
h21_grid<-CJ(metric=metrics,radius=radii); h21_refit<-rbindlist(lapply(seq_len(nrow(h21_grid)),function(i)fit_h21(h21_grid$metric[i],h21_grid$radius[i])))
h21_frozen<-read_csv(file.path(SOURCE_DIR,"sensitivity","H21_biomass_GFED_full.csv"))
h21_cmp<-merge(h21_refit,h21_frozen,by=c("metric","radius_km"),suffixes=c("","_frozen"))
num_check<-c("wind_ms","fire_pm25_center","original_estimate","original_ci95_low","original_ci95_high","fire_adjusted_estimate","fire_adjusted_ci95_low","fire_adjusted_ci95_high","fire_pm25_coef","fire_pm25_ci95_low","fire_pm25_ci95_high")
assert_true(nrow(h21_cmp)==14L && all(is.finite(h21_cmp$original_estimate_frozen)), "Missing H21 reference rows.")
for(cc in num_check) h21_cmp[,(paste0("diff_",cc)):=abs(get(cc)-get(paste0(cc,"_frozen")))]
h21_cmp[,max_abs_diff:=do.call(pmax,c(.SD,list(na.rm=TRUE))),.SDcols=paste0("diff_",num_check)]
h21_cmp[,status:=ifelse(max_abs_diff<=tol & n_model_M0==n_model_M0_frozen & n_model_M1==n_model_M1_frozen,"PASS","FAIL")]
fwrite(h21_cmp,file.path(VALIDATION_DIR,"H21_refit_vs_frozen.csv")); assert_true(all(h21_cmp$status=="PASS"),"H21 GFED refit failed numerical equality.")

# H18 environment-only model: record actual estimation N (nobs), not merely the
# 107407 complete-row input count, and reproduce the published coefficients.
env_form<-as.formula(paste0("gross_beta_mBq_m3 ~ ",paste(c("c_prop","sunspot_monthly","pm_model_centered","inv_hpbl","rhum","soilm","vel","temp_centered","temp2","temp3","sunspot_monthly:lat_diff_35","surface_u_ppm_article:soilm"),collapse=" + ")," | monitor_year_fe"))
env_fit<-feols(env_form,data=d,vcov=~loc_num,notes=FALSE,warn=TRUE)
env_ct<-as.data.table(coeftable(env_fit),keep.rownames="term"); setnames(env_ct,c("Estimate","Std. Error","t value","Pr(>|t|)"),c("estimate","se","t","p_wald"))
env_ct[, p_wald := 2*pnorm(-abs(estimate/se))]
env_ct[,`:=`(ci95_low=estimate-1.96*se,ci95_high=estimate+1.96*se,Model_N=nobs(env_fit),Monitors=uniqueN(d$loc_num))]
env_frozen<-read_csv(file.path(ROOT,"provenance","H18_primary_raw_coefficients.csv"))
env_frozen[term=="pm_centered_10",term:="pm_model_centered"]
env_frozen[term=="temp_centered_sq",term:="temp2"]
env_frozen[term=="temp_centered_cu",term:="temp3"]
env_frozen[term=="soilm:surface_u_ppm_article",term:="surface_u_ppm_article:soilm"]
env_ct[term=="soilm:surface_u_ppm_article",term:="surface_u_ppm_article:soilm"]
env_cmp<-merge(env_ct,env_frozen[,.(term,frozen_estimate=estimate,frozen_low=ci95_low,frozen_high=ci95_high)],by="term",all.x=TRUE)
assert_true(nrow(env_cmp)==12L && all(is.finite(env_cmp$frozen_estimate)), "Missing H18 reference rows.")
env_cmp[,max_abs_diff:=pmax(abs(estimate-frozen_estimate),abs(ci95_low-frozen_low),abs(ci95_high-frozen_high),na.rm=TRUE)]
env_cmp[,status:=ifelse(max_abs_diff<=tol,"PASS","FAIL")]
fwrite(env_cmp,file.path(VALIDATION_DIR,"H18_refit_vs_frozen.csv")); assert_true(all(env_cmp$status=="PASS"),"H18 environment model refit failed numerical equality.")
units <- c(c_prop="per unit proportion",sunspot_monthly="per sunspot at latitude 35 degrees",pm_model_centered="per microgram/m3",inv_hpbl="per 1/m",rhum="per percentage point",soilm="per native NARR kg/m2 at zero surface uranium",vel="per m/s",temp_centered="per degree C (linear polynomial term)",temp2="per degree C squared",temp3="per degree C cubed",`sunspot_monthly:lat_diff_35`="per sunspot x degree latitude",`surface_u_ppm_article:soilm`="per ppm x native NARR kg/m2")
env_pub<-copy(env_ct)[, `:=`(Model_unit=unname(units[term]),Outcome_unit="mBq/m3",N_pre=nrow(d),Note="Conditional coefficients; temperature terms interpreted jointly; two-sided normal-Wald p")]
fwrite(env_pub,file.path(OUTPUT_DIR,"tables","Table_environmental_covariates_corrected.csv"))

summary <- data.table(check = c("formal_rows", "formal_monitors", "Figure4_models", "Figure4_max_abs_diff",
                                "Figure4_all_N_equal", "H20B_spline_models", "H20B_max_abs_diff", "H20B_all_N_equal"),
  observed = c(nrow(d), uniqueN(d$loc_num), nrow(current),
               max(abs(c(cmp$diff_estimate, cmp$diff_low, cmp$diff_high, cmp$diff_wind)), na.rm=TRUE),
               all(cmp$N_equal), nrow(bio_ns),
               max(abs(c(ns_cmp$diff_estimate, ns_cmp$diff_low, ns_cmp$diff_high, ns_cmp$diff_wind)), na.rm=TRUE),
               all(ns_cmp$N_equal)),
  status = "PASS")
fwrite(summary, file.path(VALIDATION_DIR, "model_refit_summary.csv"))
extra_summary<-data.table(check=c("Figure5_contrast_rows","Figure5_max_abs_diff","H21_models","H21_max_abs_diff","H18_actual_nobs","H18_max_abs_diff"),
 observed=c(nrow(f5_cmp),max(abs(c(f5_cmp$diff_estimate,f5_cmp$diff_low,f5_cmp$diff_high,f5_cmp$diff_wind)),na.rm=TRUE),nrow(h21_cmp),max(h21_cmp$max_abs_diff),nobs(env_fit),max(env_cmp$max_abs_diff)),status="PASS")
fwrite(rbind(summary,extra_summary,fill=TRUE),file.path(VALIDATION_DIR,"model_refit_summary.csv"))
cat("Figure 4, Figure 5, H20B all fuels, H21 and H18 refit validation: PASS\n")
