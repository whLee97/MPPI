# Domain-shift plots with equal-width main panels and two raw-probability linear zooms.
mppi_export_domain_raw_zoom <- function(summary_dir, figure_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("cowplot", quietly = TRUE))
    stop("Domain zoom figures require ggplot2 and cowplot. See README.md.", call. = FALSE)
  g <- asNamespace("ggplot2"); c <- asNamespace("cowplot")
  mean_curves <- read.csv(file.path(summary_dir, "mean_test_curves.csv"), check.names = FALSE)
  reg_curves <- read.csv(file.path(summary_dir, "regression_test_curves.csv"), check.names = FALSE)
  order <- c("MPPI", "EW", "Classic", paste0("PPI source ", 1:3), paste0("PPI++ source ", 1:3))
  labels <- c("MPPI", "EW", "Classic", paste0("PPI (Source ", 1:3, ")"), paste0("PPI++ (Source ", 1:3, ")"))
  colors <- setNames(c("#E31A1C", "#A6DBA0", "#000000", "#A6CEE3", "#6BAED6", "#1F78B4",
                       "#FCC2D7", "#FBB4AE", "#D1B3FF"), order)
  shapes <- setNames(c(15, 7, 16, 18, 18, 18, 17, 17, 17), order)
  as_methods <- function(dat) { dat$method <- factor(dat$method, levels = order); dat }
  curve_aes <- function(xvar, points = FALSE) {
    if (xvar == "theta0") {
      if (points) g$aes(x = theta0, y = rejection_frequency, color = method, shape = method, group = method)
      else g$aes(x = theta0, y = rejection_frequency, color = method, group = method)
    } else {
      if (points) g$aes(x = radius, y = rejection_frequency, color = method, shape = method, group = method)
      else g$aes(x = radius, y = rejection_frequency, color = method, group = method)
    }
  }
  scale_methods <- function() list(g$scale_color_manual(values = colors, breaks = order, labels = labels, drop = FALSE),
                                   g$scale_shape_manual(values = shapes, breaks = order, labels = labels, drop = FALSE))
  theme_main <- function() g$theme_minimal(base_size = 12) + g$theme(
    panel.grid.minor = g$element_blank(), panel.grid.major = g$element_line(color = "#E5E5E5", linewidth = 0.32),
    plot.title = g$element_text(size = 14, hjust = 0.5, margin = g$margin(b = 5)),
    axis.title.x = g$element_text(size = 13, margin = g$margin(t = 5)),
    axis.title.y = g$element_text(size = 13, margin = g$margin(r = 5)),
    axis.text = g$element_text(size = 11), plot.margin = g$margin(4, 4, 4, 4), legend.position = "none")
  theme_zoom <- function() g$theme_minimal(base_size = 9.5) + g$theme(
    legend.position = "none", panel.grid.minor = g$element_blank(),
    panel.grid.major = g$element_line(color = "#E5E5E5", linewidth = 0.27),
    plot.title = g$element_text(size = 10.5, hjust = 0.5, margin = g$margin(b = 3)),
    axis.title.x = g$element_text(size = 9), axis.title.y = g$element_text(size = 8.5),
    axis.text = g$element_text(size = 8), plot.margin = g$margin(2, 2, 2, 2))

  make_plot <- function(dat, xvar, title, zoom = NULL, specs = NULL) {
    dat <- as_methods(dat); aes_line <- curve_aes(xvar); aes_point <- curve_aes(xvar, TRUE)
    others <- dat[dat$method != "MPPI", , drop = FALSE]
    mppi <- dat[dat$method == "MPPI", , drop = FALSE]
    small <- !is.null(zoom)
    if (small) {
      dat <- dat[dat[[xvar]] >= zoom$xlim[1] & dat[[xvar]] <= zoom$xlim[2], , drop = FALSE]
      others <- dat[dat$method != "MPPI", , drop = FALSE]
      mppi <- dat[dat$method == "MPPI", , drop = FALSE]
    }
    p <- g$ggplot() + g$geom_line(data = others, mapping = aes_line, linewidth = if (small) 0.18 else 0.22, alpha = 0.82) +
      g$geom_point(data = others, mapping = aes_point, size = 0.28, alpha = 0.80) +
      g$geom_line(data = mppi, mapping = aes_line, linewidth = if (small) 0.27 else 0.31) +
      g$geom_point(data = mppi, mapping = aes_point, size = if (small) 0.40 else 0.38) +
      g$geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = if (small) 0.38 else 0.45) + scale_methods()
    if (xvar == "theta0") {
      p <- p + g$geom_vline(xintercept = 1, linetype = "dashed", linewidth = if (small) 0.38 else 0.45)
      x_label <- expression(theta[0])
    } else x_label <- if (small) expression(r) else expression(r == paste("||", theta[0] - theta^"*", "||")[2])
    if (small) return(p + g$coord_cartesian(xlim = zoom$xlim, ylim = zoom$ylim) +
                      g$labs(title = title, x = x_label, y = "Rejection probability") + theme_zoom())
    p <- p + g$scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      g$labs(title = title, x = x_label, y = "Empirical rejection probability") + theme_main() +
      g$annotate("text", x = max(dat[[xvar]]), y = 0.085, label = "alpha == 0.05", parse = TRUE, hjust = 1, size = 4)
    if (xvar == "theta0") p <- p + g$annotate("text", x = 1.012, y = 0.92,
                                               label = "theta[0] == theta^'*'", parse = TRUE, hjust = 0, size = 4.1)
    if (!is.null(specs)) for (s in specs)
      p <- p + g$annotate("rect", xmin = s$xlim[1], xmax = s$xlim[2], ymin = s$ylim[1],
                           ymax = s$ylim[2], fill = NA, color = "#888888", linewidth = 0.20)
    p
  }

  make_zoom <- function(dat, xvar, xlim, title, min_pad) {
    z <- dat[dat[[xvar]] >= xlim[1] & dat[[xvar]] <= xlim[2], , drop = FALSE]
    stopifnot(nrow(z) > 0)
    rr <- range(c(z$rejection_frequency, if (title == "Near null") 0.05), na.rm = TRUE)
    pad <- max(min_pad, 0.12 * diff(rr))
    spec <- list(xlim = xlim, ylim = c(max(0, rr[1] - pad), min(1, rr[2] + pad)))
    list(plot = make_plot(dat, xvar, title, zoom = spec), spec = spec)
  }
  make_legend <- function(dat, xvar) {
    dat <- as_methods(dat)
    p <- g$ggplot(dat) + g$geom_line(mapping = curve_aes(xvar), linewidth = 0.36) +
      g$geom_point(mapping = curve_aes(xvar, TRUE), size = 0.95) +
      scale_methods() + g$labs(color = "Method", shape = "Method") +
      g$guides(color = g$guide_legend(nrow = 3, byrow = TRUE, override.aes = list(linewidth = 0.43, size = 1.4)),
               shape = g$guide_legend(nrow = 3, byrow = TRUE)) + g$theme_minimal() +
      g$theme(legend.position = "bottom", legend.title = g$element_text(size = 10.5),
              legend.text = g$element_text(size = 10),
              legend.key.width = grid::unit(1.15, "cm"), legend.key.height = grid::unit(0.40, "cm"))
    c$get_legend(p)
  }
  make_figure <- function(curves, ids, xvar, file) {
    linear <- curves[curves$scenario_id == ids[1], , drop = FALSE]
    nonlinear <- curves[curves$scenario_id == ids[2], , drop = FALSE]
    if (!nrow(linear) || !nrow(nonlinear)) return(invisible(NULL))
    ranges <- if (xvar == "theta0") list(c(0.985, 1.015), c(0.94, 0.96)) else
      list(c(0, 0.03), c(0.04, 0.06))
    pads <- if (xvar == "theta0") c(0.004, 0.008) else c(0.005, 0.008)
    zoom1 <- make_zoom(linear, xvar, ranges[[1]], "Near null", pads[1])
    zoom2 <- make_zoom(linear, xvar, ranges[[2]], "Moderate alternative", pads[2])
    linear_main <- make_plot(linear, xvar, "DGP-linear under Domain shift", specs = list(zoom1$spec, zoom2$spec))
    nonlinear_main <- make_plot(nonlinear, xvar, "DGP-nonlinear under Domain shift")
    zoom_col <- c$plot_grid(zoom1$plot, zoom2$plot, ncol = 1, rel_heights = c(1, 1))
    body <- c$plot_grid(linear_main, zoom_col, nonlinear_main, ncol = 3, rel_widths = c(1, 0.58, 1),
                        align = "hv", axis = "tb")
    final <- c$plot_grid(body, make_legend(linear, xvar), ncol = 1, rel_heights = c(1, 0.15))
    g$ggsave(file.path(figure_dir, file), final, width = 16, height = 6.3)
  }
  make_figure(mean_curves, c(3L, 9L), "theta0", "mean_rejection_domain_raw_zoom.pdf")
  make_figure(reg_curves, c(6L, 12L), "radius", "regression_rejection_domain_raw_zoom.pdf")
  invisible(figure_dir)
}
