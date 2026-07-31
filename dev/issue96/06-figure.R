## Figure for the issue #96 write-up: mean bias in the recovered slopes under
## the two re-draw policies, by estimator and test length.
##
##   cd dev/issue96 && Rscript 06-figure.R

library(ggplot2)

pc <- readRDS('policy-comparison.rds')
d <- pc$tab[pc$tab$policy %in% c('listwise', 'marginal'), ]

d$Policy <- factor(ifelse(d$policy == 'listwise',
                          'Re-draw on failure', 'Keep failures'),
                   levels = c('Re-draw on failure', 'Keep failures'))
d$Estimator <- factor(d$estimator, levels = c('FIML', 'DWLS'),
                      labels = c('FIML (mirt)', 'DWLS (lavaan)'))
d$Items <- factor(paste(d$nitems, 'items'), levels = c('10 items', '20 items'))

SURFACE <- '#fcfcfb'
INK     <- '#0b0b0b'
INK2    <- '#52514e'
SERIES  <- c('Re-draw on failure' = '#2a78d6', 'Keep failures' = '#eb6834')

lab <- d[d$Estimator == 'DWLS (lavaan)' & d$Items == '10 items' &
         d$sample_size == 50, ]

p <- ggplot(d, aes(sample_size, bias, colour = Policy, shape = Policy)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.4, fill = SURFACE, stroke = 0.9) +
    geom_text(data = lab, aes(label = Policy), hjust = -0.14,
              vjust = c(-0.9, 1.8), size = 3.1, show.legend = FALSE) +
    facet_grid(Estimator ~ Items) +
    scale_x_log10(breaks = c(50, 100, 250, 500, 1000),
                  labels = c('50', '100', '250', '500', '1000'),
                  expand = expansion(mult = c(0.08, 0.08))) +
    scale_y_continuous(limits = c(0, 0.22), breaks = seq(0, 0.2, 0.05),
                       expand = expansion(mult = c(0.02, 0.04))) +
    scale_colour_manual(values = SERIES) +
    scale_shape_manual(values = c(16, 17)) +
    labs(x = 'Sample size (log scale)',
         y = 'Mean bias in recovered slopes',
         title = 'Re-drawing on convergence failure shifts the reported bias',
         subtitle = paste('Vignette conditions start at N = 250. The gap is',
                          'confined to the estimator that did not cause the failures.'),
         caption = paste('SimDesign 2.26.1 - 2000 replications per cell -',
                         'philchalmers/SimDesign#96')) +
    theme_minimal(base_size = 11) +
    theme(
        plot.background   = element_rect(fill = SURFACE, colour = NA),
        panel.background  = element_rect(fill = SURFACE, colour = NA),
        panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(colour = '#e6e5e1', linewidth = 0.3),
        panel.spacing     = unit(1.1, 'lines'),
        strip.text        = element_text(colour = INK, face = 'bold', size = 9.5),
        axis.text         = element_text(colour = INK2, size = 8.5),
        axis.title        = element_text(colour = INK2, size = 9),
        plot.title        = element_text(colour = INK, face = 'bold', size = 12.5),
        plot.subtitle     = element_text(colour = INK2, size = 9, margin = margin(b = 8)),
        plot.caption      = element_text(colour = INK2, size = 7.5, hjust = 0),
        legend.position   = 'top',
        legend.title      = element_blank(),
        legend.text       = element_text(colour = INK2, size = 9),
        legend.key        = element_rect(fill = SURFACE, colour = NA)
    )

ggsave('redraw-effect.png', p, width = 7.8, height = 5.2, dpi = 160, bg = SURFACE)
cat('wrote redraw-effect.png\n')
