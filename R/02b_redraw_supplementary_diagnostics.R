# Presentation-only redraws of frozen estimates; no model refits here.
# Preserve every supplied row and confidence interval.
plot_monitor_alternative <- function(file, title, joint=FALSE) {
  d <- read_csv(file.path(SOURCE_DIR,"supplementary",file))
  if(joint)d<-d[component=="upwind_minus_downwind"]
  assert_true(nrow(d)==56L,"Unexpected monitor-FE display rows")
  d[, `:=`(fuel_label=factor(fuel_labels[fuel],levels=unname(fuel_labels)),
           metric_label=factor(metric_labels[metric],levels=unname(metric_labels)))]
  p<-ggplot(d,aes(radius_km,estimate))+
    geom_hline(yintercept=0,linetype=2,linewidth=.3)+
    geom_errorbar(aes(ymin=ci95_low,ymax=ci95_high),width=.7,linewidth=.4)+
    geom_line(linewidth=.45)+geom_point(size=1.2)+
    facet_wrap(vars(metric_label,fuel_label),ncol=4,scales="free_y")+
    scale_x_continuous(breaks=seq(20,50,10))+
    labs(title=title,subtitle="Monitor fixed effects; linear and quadratic calendar time; sample-median wind",
         x="Cumulative radius (km)",y=if(joint)"Upwind − downwind (mBq/m³)" else "Gross-beta change (mBq/m³)")+
    theme_paper()
  harmonize_facet_axes(p)
}
ps3<-plot_monitor_alternative("SFigure_monitorFE_upwind.csv","Supplementary Figure S3. Upwind associations: monitor fixed effects")
ps4<-plot_monitor_alternative("SFigure_monitorFE_direction.csv","Supplementary Figure S4. Directional contrasts: monitor fixed effects",TRUE)
save_both(ps3,"Supplementary_Figure_S3_monitorFE",11,6.7)
save_both(ps4,"Supplementary_Figure_S4_monitorFE_direction",11,6.7)

lo<-read_csv(file.path(SOURCE_DIR,"supplementary","SFigureS4_overall_LOO_plot_order.csv"))
lf<-read_csv(file.path(SOURCE_DIR,"supplementary","SFigureS4_fuel_specific_LOO_plot_order.csv"))
assert_true(nrow(lo)==140L && nrow(lf)==8L*133L,"Unexpected leave-one-monitor-out rows")
assert_true(all(lo$status=="PASS") && all(lf$status=="PASS"),"Missing frozen LOO fits")
lo[, `:=`(panel="A  Overall inventory count",full_estimate=fa[radius_km==20 & is_representative_wind==TRUE,estimate_mBq_m3])]
panel_key<-data.table(fuel=rep(c("coal","gas","oil","biomass"),2),metric=rep(c("count","heat"),each=4),
  panel=c("B  Coal count","C  Natural gas count","D  Oil count","E  Biomass count",
          "F  Coal heat input","G  Natural gas heat input","H  Oil heat input","I  Biomass heat input"))
lf<-merge(lf,panel_key,by=c("fuel","metric"),sort=FALSE)
loo_display<-rbind(lo[,.(panel,ordered_index,estimate,ci95_low,ci95_high,full_estimate)],
                  lf[,.(panel,ordered_index,estimate,ci95_low,ci95_high,full_estimate)])
loo_display[,panel:=factor(panel,levels=c("A  Overall inventory count",panel_key$panel))]
ps5<-ggplot(loo_display,aes(estimate,ordered_index))+
  geom_vline(xintercept=0,linetype=2,linewidth=.3)+
  geom_vline(data=unique(loo_display[,.(panel,full_estimate)]),aes(xintercept=full_estimate),colour="grey45",linewidth=.3)+
  geom_segment(aes(x=ci95_low,xend=ci95_high,yend=ordered_index),colour="grey55",linewidth=.15)+
  geom_point(colour="#A51C30",size=.45)+
  facet_wrap(~panel,ncol=3,scales="free_x")+
  scale_y_continuous(breaks=c(1,35,70,105,140))+
  labs(title="Supplementary Figure S5. Leave-one-monitor-out influence",
       subtitle="20 km, 90°; 140 inventory omissions and 133 per fuel–metric; full-sample wind references",
       x="Gross-beta change (mBq/m³; count or heat-input units)",y="Omitted monitor, ordered by estimate")+
  theme_paper()+theme(panel.grid.major=element_blank(),plot.subtitle=element_text(size=8))
ps5<-harmonize_facet_axes(ps5)
save_both(ps5,"Supplementary_Figure_S5_leave_one_monitor_out",11,9)
