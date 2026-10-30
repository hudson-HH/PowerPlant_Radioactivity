source('display_update_scripts/v4_2026-10-30/export_helpers.R')
library(magick)
before<-list.files('data/display_assets/frozen_v3/figures',pattern='png$',full.names=TRUE)
for(old in before){
 bits<-strsplit(basename(old),'_',fixed=TRUE)[[1]]
 id<-paste(head(bits,if(bits[1]=='Supplementary')3 else 2),collapse='_')
 new<-file.path('data/display_assets/v4_2026-10-30/figures',paste0(figure_stem(id),'.png'))
 a<-image_scale(image_read(old),'1000x');b<-image_scale(image_read(new),'1000x')
 ah<-image_info(a)$height;bh<-image_info(b)$height;h<-max(ah,bh)+38
 a<-image_extent(a,paste0('1020x',h),gravity='south',color='white')
 b<-image_extent(b,paste0('1020x',h),gravity='south',color='white')
 a<-image_annotate(a,'Before: frozen',size=22,gravity='northwest',location='+12+5',color='black')
 b<-image_annotate(b,'After: presentation edits',size=22,gravity='northwest',location='+12+5',color='black')
 joined<-image_append(image_join(list(a,b)),stack=FALSE)
 image_write(joined,path=file.path('data/display_assets/v4_2026-10-30/qa',paste0(id,'_before_after.png')),format='png')
 if(id=='Figure_2')image_write(image_scale(image_read(new),'1000x'),path=file.path('data/display_assets/v4_2026-10-30/qa/previews',paste0(id,'_after.png')),format='png')
}
cat('All 14 before/after comparisons written.\n')
