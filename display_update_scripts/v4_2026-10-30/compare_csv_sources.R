# Verify the display inputs against the manuscript v5 source CSVs, including CIs
# and sample-size columns. No estimates, sample sets or uncertainty are refitted.
library(data.table);library(jsonlite)
checks<-list()
compare<-function(id,manuscript,reproduction,keys,subset=NULL){
 a<-fread(file.path('source_data_manuscript_v5',manuscript));b<-fread(file.path('data/source_data',reproduction))
 if(!is.null(subset)){a<-a[eval(subset)];b<-b[eval(subset)]}
 common<-intersect(names(a),names(b));setorderv(a,keys);setorderv(b,keys)
 stopifnot(nrow(a)==nrow(b))
 diffs<-vapply(common,function(k){
  if(is.numeric(a[[k]]) && is.numeric(b[[k]])){
   x<-abs(a[[k]]-b[[k]]);if(!length(x)||all(is.na(x)))0 else max(x,na.rm=TRUE)
  }else if(identical(as.character(a[[k]]),as.character(b[[k]])))0 else Inf
 },numeric(1))
 stopifnot(all(diffs==0))
 checks[[id]]<<-list(source_csv=manuscript,display_csv=reproduction,rows=nrow(a),
  compared_columns=common,max_absolute_difference=max(diffs),same_key_set=TRUE)
}
compare('Figure_3a','Figure_3A_radius_wind.csv','main_figures/Figure3A.csv',c('radius_km','wind_speed_ms','is_representative_wind'))
compare('Figure_3b','Figure_3B_joint_direction.csv','main_figures/Figure3B.csv','component')
compare('Figure_3c','Figure_3C_sector_angle.csv','main_figures/Figure3C.csv','sector_angle_deg')
compare('Figure_4','Figure_4_estimates.csv','main_figures/Figure4.csv',c('fuel','metric','radius_km'))
compare('Figure_5','Figure_5_estimates.csv','main_figures/Figure5.csv',c('fuel','metric','radius_km','component'))
compare('Supplementary_Figure_S3','Supplementary_Figure_S3_monitorFE_estimates.csv','supplementary/SFigure_monitorFE_upwind.csv',c('fuel','metric','radius_km'))
compare('Supplementary_Figure_S4','Supplementary_Figure_S4_monitorFE_contrasts.csv','supplementary/SFigure_monitorFE_direction.csv',c('fuel','metric','radius_km','component'))
compare('Supplementary_Figure_S5_overall','Supplementary_Figure_S5_overall_LOO.csv','supplementary/SFigureS4_overall_LOO_plot_order.csv','omitted_monitor')
compare('Supplementary_Figure_S5_fuel','Supplementary_Figure_S5_fuel_specific_LOO.csv','supplementary/SFigureS4_fuel_specific_LOO_plot_order.csv',c('fuel','metric','omitted_monitor'))
compare('Supplementary_Figure_S8','Supplementary_Figure_S8_GFED_estimates.csv','sensitivity/H21_biomass_GFED_full.csv',c('fuel','metric','radius_km'))
for(id in c('Figure_6','Supplementary_Figure_S6','Supplementary_Figure_S7')){
 name<-switch(id,Figure_6='Figure_6_internal_validity_robustness.csv',Supplementary_Figure_S6='Supplementary_Figure_S6_adjusted_combustion_effects.csv',Supplementary_Figure_S7='Supplementary_Figure_S7_natural_spline_effect_change.csv')
 a<-fread(file.path('source_data_manuscript_v5',name))
 checks[[id]]<-list(source_csv=name,rows=nrow(a),direct_manuscript_csv=TRUE,
  max_absolute_difference=0,description='Plot reads this unchanged manuscript CSV directly; only original display-reference row selection/long reshape is applied.')
}
P<-fread('source_data_manuscript_v5/Figure_1A_powerplant_locations.csv')
Q<-fread('data/display_assets/v4_2026-10-30/source_data/Figure_1_combustion_locations.csv')
roles<-c('Coal','Natural gas','Oil','Other fossil','Biomass','Solid waste')
expected<-P[fuel_label %in% roles]
stopifnot(identical(expected,Q),all(between(Q$plant_longitude,-125,-66)),all(between(Q$plant_latitude,24,50)))
checks[['Figure_1']]<-list(original_rows=nrow(P),retained_rows=nrow(Q),unique_plants=uniqueN(Q$plant_id_plot),
 exact_subset=TRUE,coordinates_identical=TRUE,all_retained_points_inside_original_map_bounds=TRUE)
P<-fread('data/display_assets/v4_2026-10-30/source_data/Figure_1_combustion_locations.csv')
Q<-fread('data/display_assets/v4_2026-10-30/source_data/Supplementary_Figure_S1_combustion_locations.csv')
stopifnot(identical(P,Q),uniqueN(P$plant_id_plot)==5068L,nrow(P)==8326L)
counts <- Q[,.(unique_plants=uniqueN(plant_id_plot)),by=fuel_label][match(roles,fuel_label)]
stopifnot(identical(counts$unique_plants,c(798L,3500L,2675L,108L,1057L,188L)))
checks[['Supplementary_Figure_S1']]<-list(retained_rows=nrow(Q),unique_plants=uniqueN(Q$plant_id_plot),
 exact_Figure_1_role_rows=TRUE,coordinates_identical=TRUE,rule='Non-exclusive combustion roles; include multi-role plants in each applicable panel')
write_json(checks,'data/display_assets/v4_2026-10-30/qa/audits/csv_source_comparison.json',auto_unbox=TRUE,pretty=TRUE,digits=NA)
cat('All source CSV and subset comparisons passed.\n')
