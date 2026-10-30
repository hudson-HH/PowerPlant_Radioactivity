FIGURE_NAMES <- unlist(jsonlite::read_json('display_update_scripts/v4_2026-10-30/figure_names.json'))
figure_stem <- function(id) unname(FIGURE_NAMES[[id]])
options(stringsAsFactors=FALSE,warn=1,scipen=999)
options(device=function(...)grDevices::pdf(file='data/display_assets/v4_2026-10-30/qa/audits/R_default_device.pdf',...))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(grid);library(jsonlite)})
BASELINES <- new.env(parent=emptyenv())
NUMERICAL <- if(file.exists('data/display_assets/v4_2026-10-30/qa/audits/numerical_invariance.json'))read_json('data/display_assets/v4_2026-10-30/qa/audits/numerical_invariance.json',simplifyVector=FALSE)else list()
EXPORTS <- if(file.exists('data/display_assets/v4_2026-10-30/qa/audits/export_dimensions.json'))read_json('data/display_assets/v4_2026-10-30/qa/audits/export_dimensions.json',simplifyVector=FALSE)else list()
required_packages <- function(x)stopifnot(all(vapply(x,requireNamespace,logical(1),quietly=TRUE)))
assert_true <- function(x,message)if(!isTRUE(x))stop(message)
read_csv <- function(path)data.table::fread(path,showProgress=FALSE)
fuel_labels <- c(coal='Coal',gas='Natural gas',oil='Oil',biomass='Biomass')
metric_labels <- c(count='CAMPD plant count',heat='Heat input')

plot_signature <- function(p) {
 if(inherits(p,'patchwork'))return(lapply(seq_len(length(p)),function(i)plot_signature(p[[i]])))
 b <- ggplot_build(p)
 list(layers=lapply(b$data,function(d)as.data.frame(d)),
      panels=lapply(b$layout$panel_params,function(q)list(
       x_range=q$x.range,y_range=q$y.range,x_limits=q$x$limits,y_limits=q$y$limits,
       x_breaks=q$x$breaks,y_breaks=q$y$breaks,
       x_labels=q$x$get_labels(),y_labels=q$y$get_labels())),layout=as.data.frame(b$layout$layout))
}

draw_plot <- function(p) {
 if(inherits(p,'grob')){grid.newpage();grid.draw(p)}else print(p)
}

measure_panels <- function(p,id,w,h) {
 grDevices::cairo_pdf(file.path('data/display_assets/v4_2026-10-30/qa/audits',paste0(id,'_geometry.pdf')),width=w,height=h,family='Helvetica')
 draw_plot(p);grid.force()
 vps <- grid.ls(viewports=TRUE,grobs=FALSE,print=FALSE)
 hits <- which(vps$type=='vpListing' & grepl('^panel(?:-[0-9]+)*[.]',vps$name,perl=TRUE))
 boxes <- lapply(hits,function(j){
  path <- c(strsplit(vps$vpPath[j],'::',fixed=TRUE)[[1]],vps$name[j]);path <- path[path!='ROOT' & nzchar(path)]
  upViewport(0);downViewport(do.call(vpPath,as.list(path)))
  a<-deviceLoc(unit(0,'npc'),unit(0,'npc'),valueOnly=TRUE); b<-deviceLoc(unit(1,'npc'),unit(1,'npc'),valueOnly=TRUE)
  c(a$x,a$y,b$x,b$y)*72
 })
 upViewport(0);dev.off()
 stopifnot(length(boxes)>0)
 rect <- unique(do.call(rbind,boxes));rect <- rect[rect[,3]>rect[,1]+1 & rect[,4]>rect[,2]+1,,drop=FALSE]
 rect <- rect[order(-round(rect[,4],2),rect[,1]),,drop=FALSE]
 # Shared rows/columns inferred from the actually drawn rectangles, including fixed-aspect maps.
 ids <- paste0('p',seq_len(nrow(rect)))
 rows <- split(ids,round(rect[,4],1));cols <- split(ids,round(rect[,1],1))
 rows <- unname(Filter(function(x)length(x)>1,rows));cols <- unname(Filter(function(x)length(x)>1,cols))
 exemptions <- list()
 if(id=='Figure_3')exemptions <- list(list(panels=ids,checks='panel-width',reason='Preserve supplied 1.35:1:1 layout; panel a carries radius and wind curves.'))
 # S2 deliberately retains original geographic bounds and fixed map aspects.
 if(id=='Supplementary_Figure_S2'){
  rows<-list(ids)
  exemptions <- list(list(panels=ids,checks=c('panel-width','horizontal-gutter'),
   reason='Original regional maps use distinct geographic bounds and fixed aspects; preserve those original map scales.'),
   list(panels=ids[3],checks='row',reason='Preserve the original regional geographic aspect ratio; Eastern CONUS is centered in its fixed layout cell.'))
 }
 man <- list(schema_version=1,backend='r-forced-cairo-PDF-viewport',figure=list(width_pt=w*72,height_pt=h*72),
  panels=lapply(seq_len(nrow(rect)),function(i)list(id=ids[i],bbox_pt=rect[i,])),
  row_groups=rows,column_groups=cols,exemptions=exemptions)
 write_json(man,file.path('data/display_assets/v4_2026-10-30/qa/audits',paste0(id,'.alignment-layout.json')),auto_unbox=TRUE,pretty=TRUE,digits=NA)
 py<-Sys.getenv('QA_PYTHON',unset='python3')
 status<-system2(py,c('display_update_scripts/v4_2026-10-30/qa_helpers/audit_panel_alignment.py',
  file.path('data/display_assets/v4_2026-10-30/qa/audits',paste0(id,'.alignment-layout.json')),'--json-out',
  file.path('data/display_assets/v4_2026-10-30/qa/audits',paste0(id,'.alignment.json')),'--strict'))
 if(status!=0)stop(id,' rendered panel alignment did not pass')
 invisible(man)
}

export_figure <- function(p,id,w,h,baseline=NULL) {
 cat('Exporting',id,w,h,'inches\n')
 original <- attr(p,'baseline_plot');new <- attr(p,'new_plot')
 if(!is.null(original)){baseline<-original}else new<-p
 if(!is.null(baseline)){
  oldsig<-plot_signature(baseline);newsig<-plot_signature(new)
  equal <- isTRUE(all.equal(oldsig,newsig,tolerance=0,check.attributes=TRUE))
  if(!equal)stop(id,' plotted-value/axis invariance failed: ',paste(all.equal(oldsig,newsig),collapse='; '))
  saveRDS(list(before=oldsig,after=newsig),file.path('data/display_assets/v4_2026-10-30/qa/audits',paste0(id,'_plot_signatures.rds')))
  NUMERICAL[[id]] <<- list(plotted_values_identical=equal,axis_limits_breaks_labels_identical=equal,
   layer_colours_shapes_lines_identical=equal,comparison='All ggplot-built layer columns, facet layout, numeric ranges, breaks and tick labels; exact tolerance 0')
 }
 man <- measure_panels(p,id,w,h)
 stem <- file.path('data/display_assets/v4_2026-10-30/figures',figure_stem(id))
 grDevices::cairo_pdf(paste0(stem,'.pdf'),width=w,height=h,family='Helvetica',bg='white',onefile=FALSE)
 draw_plot(p);dev.off()
 svglite::svglite(paste0(stem,'.svg'),width=w,height=h,bg='white');draw_plot(p);dev.off()
 ragg::agg_png(paste0(stem,'.png'),width=w,height=h,units='in',res=600,background='white');draw_plot(p);dev.off()
 ragg::agg_tiff(paste0(stem,'.tiff'),width=w,height=h,units='in',res=600,compression='lzw',background='white');draw_plot(p);dev.off()
 # Separate screen-review preview; submission PNG/TIFF above are always 600 dpi.
 preview_dpi<-110
 ragg::agg_png(file.path('data/display_assets/v4_2026-10-30/qa/previews',paste0(id,'_after.png')),width=w,height=h,units='in',res=preview_dpi,background='white');draw_plot(p);dev.off()
 EXPORTS[[id]] <<- list(width_in=w,height_in=h,dpi=600,panels=length(man$panels))
 write_json(NUMERICAL,'data/display_assets/v4_2026-10-30/qa/audits/numerical_invariance.json',auto_unbox=TRUE,pretty=TRUE,digits=NA)
 write_json(EXPORTS,'data/display_assets/v4_2026-10-30/qa/audits/export_dimensions.json',auto_unbox=TRUE,pretty=TRUE,digits=NA)
 invisible(p)
}

save_both <- function(plot,stem,width,height){
 map <- c(Figure3_overall_inventory='Figure_3',Figure4_upwind_actual_operation_monitorYearFE='Figure_4',
 Figure5_upwind_minus_downwind_monitorYearFE='Figure_5',Figure_GFED_biomass_sensitivity='Supplementary_Figure_S8',
 Supplementary_Figure_S3_monitorFE='Supplementary_Figure_S3',Supplementary_Figure_S4_monitorFE_direction='Supplementary_Figure_S4',
 Supplementary_Figure_S5_leave_one_monitor_out='Supplementary_Figure_S5')
 id<-unname(map[stem]);stopifnot(!is.na(id))
 baseline<-if(exists(id,envir=BASELINES,inherits=FALSE))BASELINES[[id]]else NULL
 export_figure(plot,id,width,height,baseline)
}
