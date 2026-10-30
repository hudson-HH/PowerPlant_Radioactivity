# Figure 3 is redrawn from its numerical exports, not from historical pixels.
# Sourced by 02_reproduce_figures_tables.R after theme_paper and save_both.
# This step does not refit the quasi-daily mixed models.
fa <- read_csv(file.path(SOURCE_DIR,"Figure_3A_radius_wind.csv"))
fb <- read_csv(file.path(SOURCE_DIR,"Figure_3B_joint_direction.csv"))
fc <- read_csv(file.path(SOURCE_DIR,"Figure_3C_sector_angle.csv"))
assert_true(nrow(fa)==66L && nrow(fb)==3L && nrow(fc)==3L,"Unexpected Figure 3 source rows")
assert_true(all(is.finite(fa$estimate_mBq_m3)) && all(is.finite(fb$estimate_mBq_m3)) && all(is.finite(fc$estimate_mBq_m3)),"Missing Figure 3 estimates")
fa[, wind_label:=ifelse(is_representative_wind,"Mean wind (2.986 m/s)",paste0(wind_speed_ms," m/s"))]
fa[, wind_label:=factor(wind_label,levels=c(paste0(1:5," m/s"),"Mean wind (2.986 m/s)"))]
wind_cols<-c("#BCD2E8","#8CAECD","#638AB2","#3D668E","#173B61","#111111")
names(wind_cols)<-levels(fa$wind_label)
p3a<-ggplot(fa,aes(radius_km,estimate_mBq_m3,colour=wind_label,group=wind_label))+
  geom_hline(yintercept=0,linetype=2,linewidth=.3,colour="grey40")+
  geom_ribbon(data=fa[is_representative_wind==TRUE],aes(ymin=ci95_low_mBq_m3,ymax=ci95_high_mBq_m3),fill="grey70",alpha=.18,colour=NA)+
  geom_line(aes(linewidth=is_representative_wind)) + geom_point(size=1.1)+
  scale_linewidth_manual(values=c(`FALSE`=.35,`TRUE`=.8),guide="none")+
  scale_colour_manual(values=wind_cols)+scale_x_continuous(breaks=c(20,50,100,150))+
  labs(title="a  Radius and wind",x="Cumulative radius (km)",y="Conditional change (mBq/m³ per plant)",colour=NULL)+
  theme_paper()+theme(legend.position="bottom",legend.text=element_text(size=8),plot.title=element_text(size=11))+
  guides(colour=guide_legend(nrow=3,byrow=TRUE))
fb[,component_label:=factor(component,levels=c("downwind","upwind","difference"),labels=c("Downwind","Upwind","Upwind minus\ndownwind"))]
p3b<-ggplot(fb,aes(component_label,estimate_mBq_m3))+
  geom_hline(yintercept=0,linetype=2,linewidth=.3)+
  geom_errorbar(aes(ymin=ci95_low_mBq_m3,ymax=ci95_high_mBq_m3),width=.16,linewidth=.5)+
  geom_point(size=2,colour="#234A84")+
  labs(title="b  Joint direction",x=NULL,y=NULL)+theme_paper()+theme(plot.title=element_text(size=11),axis.text.x=element_text(size=8))
p3c<-ggplot(fc,aes(factor(sector_angle_deg),estimate_mBq_m3))+
  geom_hline(yintercept=0,linetype=2,linewidth=.3)+
  geom_errorbar(aes(ymin=ci95_low_mBq_m3,ymax=ci95_high_mBq_m3),width=.16,linewidth=.5)+
  geom_point(size=2,colour="#234A84")+
  labs(title="c  Sector width",x="Sector width (degrees)",y=NULL)+theme_paper()+theme(plot.title=element_text(size=11))
p3_baseline<-p3a+p3b+p3c+plot_layout(nrow=1,widths=c(1.35,1,1))+
  plot_annotation(title="Overall plant inventory and ambient gross-beta activity",
    subtitle="140 monitors; 558,915 quasi-daily rows; directional comparisons evaluated at mean wind 2.986 m/s")
# This legend belongs only to a; b and c are the original plot objects.
p3a_inside <- p3a + theme(legend.position='inside',legend.position.inside=c(.94,.98),
  legend.justification=c(1,1),legend.background=element_blank(),
  legend.margin=margin(1,1,1,1),legend.key.height=unit(3.2,'mm'),
  legend.key.width=unit(5,'mm'),legend.spacing.y=unit(0,'mm')) +
  guides(colour=guide_legend(ncol=1,byrow=TRUE))
p3 <- p3a_inside+p3b+p3c+plot_layout(nrow=1,widths=c(1.35,1,1))
BASELINES[['Figure_3']] <- p3_baseline
# Optional curator hook: alignment QA before export in the release environment.
# Ordinary reproduction does not require a local skill installation.
qa_hook<-Sys.getenv("RADNET_FIGURE3_ALIGNMENT_HOOK",unset="")
if(nzchar(qa_hook))source(qa_hook,local=FALSE)
save_both(p3,"Figure3_overall_inventory",11,5.8)
