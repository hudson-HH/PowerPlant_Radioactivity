# Recover exact panel-specific sample sizes and numerically verify Figure 6.
required_packages(c("data.table","fixest","lme4"))
suppressPackageStartupMessages({library(data.table);library(fixest)})
d <- fread(file.path(MODEL_DATA_DIR,"filter_period_primary_2001_2024.csv.gz"))
d[, time_fe := interaction(collect_start,collect_end,drop=TRUE,lex.order=TRUE)]
env <- c("c_prop","sunspot_monthly","pm_model_centered","inv_hpbl","rhum","soilm","vel","temp_centered","temp2","temp3","sunspot_monthly:lat_diff_35","surface_u_ppm_article:soilm")
mixed_env <- c(env,"surface_u_ppm_article","lat_diff_35","year_model","I(year_model^2)")
opts <- list(optimizer="bobyqa",optCtrl=list(maxfun=200000),check.conv.singular="ignore",check.rankX="message+drop.cols")
if("autoscale" %in% names(formals(lme4::lmerControl))) opts$autoscale<-TRUE
control <- do.call(lme4::lmerControl,opts)
extract_slope <- function(fit,term,wind,mixed=FALSE){
  b<-if(mixed)lme4::fixef(fit) else coef(fit); v<-as.matrix(vcov(fit))
  ti<-paste0(term,":vel"); a<-b[term]+wind*b[ti]
  se<-sqrt(v[term,term]+wind^2*v[ti,ti]+2*wind*v[term,ti])
  data.table(estimate=unname(a),ci_low=unname(a-1.96*se),ci_high=unname(a+1.96*se),N_model=nobs(fit))
}
grid<-CJ(fuel=c("coal","gas","oil","biomass"),metric=c("count","heat"))
out<-list(); legacy_out<-list()
for(i in seq_len(nrow(grid))){
  fuel<-grid$fuel[i];metric<-grid$metric[i]
  ex<-paste0("fpmean_campd_",fuel,"_upwind_",if(metric=="count")"n_plants" else "heat_observed","_20km_90deg")
  x<-copy(d[is.finite(get(ex))]);x[,exposure:=get(ex)/if(metric=="heat")10000 else 1]
  wind<-median(x$vel)
  main<-feols(as.formula(paste0("gross_beta_mBq_m3 ~ exposure + exposure:vel + ",paste(env,collapse=" + ")," | monitor_year_fe")),data=x,vcov=~loc_num,notes=FALSE)
  strict<-feols(as.formula(paste0("gross_beta_mBq_m3 ~ exposure + exposure:vel + ",paste(env,collapse=" + ")," | monitor_year_fe + time_fe")),data=x,vcov=~loc_num,notes=FALSE)
  nopm<-feols(as.formula(paste0("gross_beta_mBq_m3 ~ exposure + exposure:vel + ",paste(setdiff(env,"pm_model_centered"),collapse=" + ")," | monitor_year_fe")),data=x,vcov=~loc_num,notes=FALSE)
  a<-rbind(extract_slope(main,"exposure",wind)[,model:="Monitor × year FE"],extract_slope(strict,"exposure",wind)[,model:="Monitor × year + common-time FE"])[,panel:="A"]
  cc<-rbind(extract_slope(main,"exposure",wind)[,model:="PM₂.₅ adjusted"],extract_slope(nopm,"exposure",wind)[,model:="No PM₂.₅ adjustment"])[,panel:="C"]
  x[,loc_num_factor:=factor(loc_num)]
  original<-lme4::lmer(as.formula(paste0("gross_beta_mBq_m3 ~ exposure + exposure:vel + ",paste(mixed_env,collapse=" + ")," + (1 | loc_num_factor)")),data=x,REML=TRUE,na.action=na.fail,control=control)
  x[,exposure_between:=mean(exposure),by=monitor_year_fe];x[,exposure_within:=exposure-exposure_between]
  mundlak<-lme4::lmer(as.formula(paste0("gross_beta_mBq_m3 ~ exposure_within + exposure_between + exposure_within:vel + exposure_between:vel + ",paste(mixed_env,collapse=" + ")," + (1 | loc_num_factor)")),data=x,REML=TRUE,na.action=na.fail,control=control)
  bb<-rbind(extract_slope(original,"exposure",wind,TRUE)[,model:="Original random-intercept"],extract_slope(mundlak,"exposure_within",wind,TRUE)[,model:="Within monitor-year"],extract_slope(mundlak,"exposure_between",wind,TRUE)[,model:="Between monitor-years"])[,panel:="B"]
  r<-rbind(a,bb,cc);r[,`:=`(fuel=fuel,metric=metric,N_pre=nrow(x),N_monitors=uniqueN(x$loc_num),wind_ms=wind,
     wind_reference="Shared fuel-by-metric baseline sample median")];out[[i]]<-r
  legacy_b<-rbind(extract_slope(original,"exposure",1,TRUE)[,model:="Original random-intercept"],extract_slope(mundlak,"exposure_within",1,TRUE)[,model:="Within monitor-year"],extract_slope(mundlak,"exposure_between",1,TRUE)[,model:="Between monitor-years"])[,panel:="B"]
  legacy_r<-rbind(a,legacy_b,cc);legacy_r[,`:=`(fuel=fuel,metric=metric,N_pre=nrow(x),
    N_monitors=uniqueN(x$loc_num),wind_ms=ifelse(panel=="B",1,wind))];legacy_out[[i]]<-legacy_r
}
results<-rbindlist(out)
frozen<-rbindlist(lapply(c("A","B","C"),function(p){z<-read_csv(file.path(SOURCE_DIR,"main_figures",paste0("Figure6",p,".csv")));z[,panel:=p];z[,.(panel,fuel,metric,model,frozen_estimate=estimate,frozen_low=ci_low,frozen_high=ci_high)]}))
cmp<-merge(results,frozen,by=c("panel","fuel","metric","model"),all.x=TRUE)
assert_true(nrow(cmp)==56L && all(is.finite(cmp$frozen_estimate)),"Figure6 reference row mismatch.")
cmp[,max_abs_diff:=pmax(abs(estimate-frozen_estimate),abs(ci_low-frozen_low),abs(ci_high-frozen_high))]
# Mixed-model optimizer tolerances are separate from exact OLS reproduction.
cmp[,tolerance:=ifelse(panel=="B",1e-5,1e-7)]
cmp[,status:=ifelse(max_abs_diff<=tolerance,"PASS","FAIL")]
fwrite(cmp,file.path(VALIDATION_DIR,"Figure6_refit_vs_frozen.csv"))
fwrite(results,file.path(OUTPUT_DIR,"tables","Figure6_panel_sample_sizes_and_estimates.csv"))
assert_true(all(cmp$status=="PASS"),"Figure6 numerical comparison failed; inspect validation report.")
legacy_results<-rbindlist(legacy_out)
legacy_frozen<-rbindlist(lapply(c("A","B","C"),function(p){z<-read_csv(file.path(ROOT,"provenance","revision_v3","legacy_1ms",paste0("Figure6",p,".csv")));z[,panel:=p];z[,.(panel,fuel,metric,model,frozen_estimate=estimate,frozen_low=ci_low,frozen_high=ci_high)]}))
legacy_cmp<-merge(legacy_results,legacy_frozen,by=c("panel","fuel","metric","model"),all.x=TRUE)
assert_true(nrow(legacy_cmp)==56L && all(is.finite(legacy_cmp$frozen_estimate)),"Legacy Figure6 reference mismatch.")
legacy_cmp[,max_abs_diff:=pmax(abs(estimate-frozen_estimate),abs(ci_low-frozen_low),abs(ci_high-frozen_high))]
legacy_cmp[,tolerance:=ifelse(panel=="B",1e-5,1e-7)]
legacy_cmp[,status:=ifelse(max_abs_diff<=tolerance,"PASS","FAIL")]
fwrite(legacy_cmp,file.path(VALIDATION_DIR,"Figure6_legacy_1ms_same_fit_comparison.csv"))
assert_true(all(legacy_cmp$status=="PASS"),"Legacy anchor validation failed; model fits changed beyond wind evaluation.")
cat("Figure6 refit and panel sample-size recovery: PASS\n")
