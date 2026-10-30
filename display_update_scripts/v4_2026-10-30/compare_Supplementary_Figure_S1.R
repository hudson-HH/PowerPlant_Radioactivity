# Update only S1's established frozen/new comparison using the R image workflow.
library(magick)
old <- 'data/display_assets/frozen_v3/figures/Supplementary_Figure_S1_individual_fuel_type_maps.png'
new <- 'data/display_assets/v4_2026-10-30/figures/Supplementary_Figure_S1_individual_fuel_type_maps.png'
a <- image_scale(image_read(old),'1000x')
b <- image_scale(image_read(new),'1000x')
h <- max(image_info(a)$height,image_info(b)$height)+38
a <- image_extent(a,paste0('1020x',h),gravity='center',color='white')
b <- image_extent(b,paste0('1020x',h),gravity='center',color='white')
a <- image_annotate(a,'Before: frozen',size=22,gravity='northwest',location='+12+5',color='black')
b <- image_annotate(b,'After: non-exclusive roles',size=22,gravity='northwest',location='+12+5',color='black')
image_write(image_append(image_join(list(a,b)),stack=FALSE),
            path='data/display_assets/v4_2026-10-30/qa/Supplementary_Figure_S1_before_after.png',format='png')
