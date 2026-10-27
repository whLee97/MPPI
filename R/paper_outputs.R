# Paper figures retain the manuscript's panel layout, colours and comparison methods.
mppi_export_paper <- function(summary_dir, dxa_dir, output_dir) {
  required <- c("simulation_summary.csv", "mean_test_curves.csv", "regression_test_curves.csv", "normal_diagnostics.csv")
  mppi_assert(all(file.exists(file.path(summary_dir,required))), "Missing simulation summaries; run summarize_results.R first.")
  out_pdf <- file.path(output_dir,"figures"); out_tab <- file.path(output_dir,"tables")
  dir.create(out_pdf,recursive=TRUE,showWarnings=FALSE)
  dir.create(out_tab,recursive=TRUE,showWarnings=FALSE)
  summary <- read.csv(file.path(summary_dir,"simulation_summary.csv"),check.names=FALSE)
  curves <- read.csv(file.path(summary_dir,"mean_test_curves.csv"),check.names=FALSE)
  method_order <- c('MPPI','EW','Classic','PPI source 1','PPI source 2','PPI source 3',
                    'PPI++ source 1','PPI++ source 2','PPI++ source 3')
  method_labels <- c('MPPI','EW','Classic','PPI (Source 1)','PPI (Source 2)','PPI (Source 3)',
                     'PPI++ (Source 1)','PPI++ (Source 2)','PPI++ (Source 3)')
  table_method_order <- c('MPPI','Classic','PPI source 1','PPI source 2','PPI source 3',
                          'PPI++ source 1','PPI++ source 2','PPI++ source 3','EW')
  method_map <- setNames(method_labels, method_order)
  colors <- c('#D73027','#2A9D6F','#222222','#8ECAE6','#3696C7','#126782','#E9A3C9','#E76F51','#7B61A8')
  point_chars <- c(15, 4, 19, 18, 18, 18, 17, 17, 17)
  line_widths <- c(2.7, 1.8, 1.9, rep(1.55, 6))

  panel_plot <- function(dat, title) {
    xlim <- range(dat$theta0); ylim <- c(0, 1.05)
    plot(NA, xlim=xlim, ylim=ylim, xlab=expression(theta[0]), ylab='Empirical rejection probability',
         axes=FALSE, main=title, cex.main=1.18, cex.lab=1.05, font.main=2,
         xaxs='i', yaxs='i')
    rect(xlim[1], ylim[1], xlim[2], ylim[2], col='white', border=NA)
    xt <- pretty(xlim, 5); yt <- seq(0,1,0.2)
    abline(v=xt, h=yt, col='#E9ECEF', lwd=0.7)
    truth <- dat$theta0[which.min(dat$distance)]
    abline(h=0.05, v=truth, lty=2, col='#333333', lwd=1.1)
    for (j in seq_along(method_order)) {
      z <- dat[dat$method == method_order[j],]
      z <- z[order(z$theta0),]
      lines(z$theta0, z$rejection_frequency, col=colors[j], lwd=line_widths[j])
      idx <- unique(round(seq(1,nrow(z),length.out=10)))
      points(z$theta0[idx], z$rejection_frequency[idx], col=colors[j], pch=point_chars[j], cex=0.48)
    }
    axis(1, at=xt, las=1, cex.axis=0.94, col='#666666', col.axis='#222222')
    axis(2, at=yt, las=1, cex.axis=0.94, col='#666666', col.axis='#222222')
    box(col='#777777', lwd=0.8)
    text(xlim[2] - 0.035*diff(xlim), 0.077, expression(alpha == 0.05),
         adj=c(1,0), cex=1.12)
    label_x <- if (grepl('Domain shift: nonlinear', title)) xlim[1] + 0.045*diff(xlim) else truth + 0.045*diff(xlim)
    text(label_x, 0.91, expression(theta[0] == theta^'*'),
         adj=c(0,0.5), cex=1.15)
  }

  make_shift_figure <- function(scenario_ids, label, filename) {
    pdf(file.path(out_pdf, filename), width=11.2, height=4.5, pointsize=14,
        family='Helvetica', bg='white', useDingbats=FALSE)
    op <- par(no.readonly=TRUE)
    on.exit({par(op); dev.off()}, add=TRUE)
    layout(matrix(1:3, nrow=1), widths=c(1,1,0.42))
    par(mar=c(4.25,4.55,2.65,0.75), mgp=c(2.55,0.72,0), tcl=-0.22,
        bg='white', family='Helvetica')
    for (k in seq_along(c('linear','nonlinear'))) {
      dgp <- c('linear','nonlinear')[k]
      dat <- curves[curves$scenario_id == scenario_ids[k],]
      panel_plot(dat, sprintf('%s: %s DGP', label, dgp))
    }
    par(mar=c(0,0,0,0))
    plot.new()
    legend('center', bty='n', title='Method', legend=method_labels, col=colors,
           lwd=line_widths, pch=point_chars, pt.cex=0.78, cex=1.06,
           y.intersp=1.05, x.intersp=0.82, seg.len=2.2)
  }

  if (all(c(1L,7L) %in% curves$scenario_id)) make_shift_figure(c(1L,7L),'Homogeneous','mean_rejection_homogeneous.pdf')
  if (all(c(2L,8L) %in% curves$scenario_id)) make_shift_figure(c(2L,8L),'Covariate shift','mean_rejection_covariate.pdf')
  # Domain-shift mean results use the final equal-width main panels with linear raw zooms.

  reg_curves <- read.csv(file.path(summary_dir, 'regression_test_curves.csv'), check.names=FALSE)
  normal_diag <- read.csv(file.path(summary_dir, 'normal_diagnostics.csv'), check.names=FALSE)

  regression_panel <- function(dat, title) {
    xlim <- range(dat$radius); ylim <- c(0, 1.05)
    plot(NA, xlim=xlim, ylim=ylim, xlab=expression(r == paste("||", theta[0] - theta^"*", "||")[2]),
         ylab='Empirical rejection probability', axes=FALSE, main=title, cex.main=1.18,
         cex.lab=1.05, font.main=2, xaxs='i', yaxs='i')
    rect(xlim[1], ylim[1], xlim[2], ylim[2], col='white', border=NA)
    xt <- seq(0,0.1,0.02); yt <- seq(0,1,0.2)
    abline(v=xt, h=yt, col='#E9ECEF', lwd=0.7)
    abline(h=0.05, lty=2, col='#333333', lwd=1.1)
    for (j in seq_along(method_order)) {
      z <- dat[dat$method == method_order[j],]
      z <- z[order(z$radius),]
      lines(z$radius, z$rejection_frequency, col=colors[j], lwd=line_widths[j])
      points(z$radius, z$rejection_frequency, col=colors[j], pch=point_chars[j], cex=0.48)
    }
    axis(1, at=xt, cex.axis=0.94, col='#666666', col.axis='#222222')
    axis(2, at=yt, las=1, cex.axis=0.94, col='#666666', col.axis='#222222')
    box(col='#777777', lwd=0.8)
    text(xlim[2] - 0.035*diff(xlim), 0.077, expression(alpha == 0.05),
         adj=c(1,0), cex=1.12)
  }

  make_regression_figure <- function(scenario_ids, label, filename) {
    pdf(file.path(out_pdf, filename), width=11.2, height=4.5, pointsize=14,
        family='Helvetica', bg='white', useDingbats=FALSE)
    op <- par(no.readonly=TRUE)
    on.exit({par(op); dev.off()}, add=TRUE)
    layout(matrix(1:3, nrow=1), widths=c(1,1,0.42))
    par(mar=c(4.25,4.55,2.65,0.75), mgp=c(2.55,0.72,0), tcl=-0.22,
        bg='white', family='Helvetica')
    for (k in seq_along(c('linear','nonlinear'))) {
      dgp <- c('linear','nonlinear')[k]
      dat <- reg_curves[reg_curves$scenario_id == scenario_ids[k],]
      regression_panel(dat, sprintf('%s: %s DGP', label, dgp))
    }
    par(mar=c(0,0,0,0)); plot.new()
    legend('center', bty='n', title='Method', legend=method_labels, col=colors,
           lwd=line_widths, pch=point_chars, pt.cex=0.78, cex=1.06,
           y.intersp=1.05, x.intersp=0.82, seg.len=2.2)
  }

  if (all(c(4L,10L) %in% reg_curves$scenario_id)) make_regression_figure(c(4L,10L),'Homogeneous','regression_rejection_homogeneous.pdf')
  if (all(c(5L,11L) %in% reg_curves$scenario_id)) make_regression_figure(c(5L,11L),'Covariate shift','regression_rejection_covariate.pdf')
  if (all(c(3L,9L) %in% curves$scenario_id) && all(c(6L,12L) %in% reg_curves$scenario_id))
    mppi_export_domain_raw_zoom(summary_dir, out_pdf)

  normal_mppi <- normal_diag[normal_diag$target=='mean' & normal_diag$method=='MPPI' & normal_diag$coordinate==1L,]
  normal_ids <- matrix(c(1L,2L,3L,7L,8L,9L), nrow=2, byrow=TRUE)
  normal_titles <- matrix(c('Homogeneous, linear','Covariate shift, linear','Domain shift, linear',
                            'Homogeneous, nonlinear','Covariate shift, nonlinear','Domain shift, nonlinear'), nrow=2, byrow=TRUE)
  if (all(normal_ids %in% normal_mppi$scenario_id)) {
  pdf(file.path(out_pdf,'normal_diag_mean_MPPI_all_settings.pdf'), width=10.5, height=6.8, family='Helvetica', useDingbats=FALSE)
  op <- par(mfrow=c(2,3), mar=c(3.6,3.7,2.3,0.7), oma=c(0.2,0.2,0.2,0.2), mgp=c(2.25,0.65,0), tcl=-0.22)
  for (i in seq_len(nrow(normal_ids))) for (j in seq_len(ncol(normal_ids))) {
    z <- normal_mppi$z[normal_mppi$scenario_id==normal_ids[i,j]]
    if (!length(z)) {plot.new(); next}
    extent <- max(4, ceiling(max(abs(z),na.rm=TRUE)/0.4)*0.4 + 0.4)
    h <- hist(z, breaks=seq(-extent,extent,by=0.4), plot=FALSE)
    ymax <- max(h$density, dnorm(0))*1.08
    plot(h, freq=FALSE, xlim=c(-4,4), ylim=c(0,ymax), col='#D9EAF4', border='white',
         main=normal_titles[i,j], xlab=expression(Z[b]), ylab='Density', axes=FALSE)
    abline(v=seq(-4,4,2), h=pretty(c(0,ymax),4), col='#E5E5E5', lwd=0.7)
    plot(h, freq=FALSE, add=TRUE, col=adjustcolor('#78B4D4',0.55), border='white')
    curve(dnorm(x), from=-4, to=4, add=TRUE, col='#E63946', lwd=2)
    axis(1, at=seq(-4,4,2)); axis(2, las=1); box(col='#B0B0B0')
    legend('topright', legend=c('Standard normal'), col='#E63946', lwd=2, bty='n', cex=0.72)
  }
  par(op); dev.off()
  }

  # Export numeric tables, without rounding the underlying reference CSVs.
  write_table <- function(dat, name) mppi_atomic_csv(dat,file.path(out_tab,paste0(name,".csv")))
  for (target in c("mean","regression")) for (shift in c("homogeneous","covariate","domain")) {
    dat <- summary[summary$target==target & summary$shift==shift,]
    if (!nrow(dat)) next
    dat <- dat[order(match(dat$method,table_method_order),match(dat$dgp,c("linear","nonlinear"))),]
    write_table(dat[,c("method","dgp","coverage","coverage_mcse","mean_volume","B_expected","B_valid","B_failed","B_missing")],paste(target,shift,"acp_volume",sep="_"))
    weights <- dat[dat$method=="MPPI",c("dgp","mean_w0","mean_w1","mean_w2","mean_w3")]
    write_table(weights,paste(target,shift,"mppi_weights",sep="_"))
  }
  dxa_path <- file.path(dxa_dir,"dxa_results.csv")
  if (file.exists(dxa_path)) {
    dxa <- read.csv(dxa_path,check.names=FALSE)
    main <- (dxa$Group=="Male 40-59" & dxa$Threshold==30) | (dxa$Group=="Female 40-59" & dxa$Threshold==40)
    write_table(dxa[main,],"dxa_main")
    write_table(dxa[!main,],"dxa_threshold_sensitivity")
  }
  pvalue_path <- file.path(dxa_dir,"dxa_pvalue_curves.csv")
  if (file.exists(pvalue_path)) {
    pvalues <- read.csv(pvalue_path,check.names=FALSE)
    pdf(file.path(out_pdf,"dxa_pvalue_curves.pdf"),width=11.2,height=4.5,pointsize=14,family="Helvetica",useDingbats=FALSE)
    op <- par(no.readonly=TRUE)
    layout(matrix(1:3,nrow=1),widths=c(1,1,0.42))
    par(mar=c(4.25,4.55,2.65,0.75),mgp=c(2.55,0.72,0),tcl=-0.22)
    for (group in c("Male 40-59","Female 40-59")) {
      dat <- pvalues[pvalues$Group==group,]
      plot(NA,xlim=c(0,1),ylim=c(0,1.05),xlab=expression(theta[0]),ylab="Wald p-value",
           main=paste0(group,", threshold ",unique(dat$Threshold)),axes=FALSE,xaxs="i",yaxs="i",cex.main=1.18)
      abline(v=seq(0,1,.2),h=seq(0,1,.2),col="#E9ECEF",lwd=.7)
      abline(h=.05,lty=2,col="#333333")
      for(j in seq_along(method_order)) {
        z <- dat[dat$Method==method_order[j],];z <- z[order(z$theta0),]
        lines(z$theta0,z$p_value,col=colors[j],lwd=line_widths[j])
      }
      axis(1,at=seq(0,1,.2));axis(2,at=seq(0,1,.2),las=1);box(col="#777777")
    }
    par(mar=c(0,0,0,0));plot.new()
    legend("center",bty="n",title="Method",legend=method_labels,col=colors,lwd=line_widths,cex=1.06,seg.len=2.2)
    par(op);dev.off()
  }
  invisible(output_dir)
}
