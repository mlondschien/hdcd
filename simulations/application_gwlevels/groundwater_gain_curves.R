# PLOT GROUNDWATER GAIN CURVES OVER EACH OTHER

# This file generates the plots for the gain curves
library(cowplot)
library(hdcdwithmissingvalues)
directory <- './simulations/application_gwlevels'

data <- read.csv(file.path(directory, 'groundwater_levels.csv'))

res_pair_50 <- readRDS(file.path(directory, 'gw_4_big_lambda_line_pair_50_0.025'))
res_average_50 <- readRDS(file.path(directory, 'gw_4_big_lambda_line_average_50_0.025'))
res_loh_50 <- readRDS(file.path(directory, 'gw_4_big_lambda_line_loh_50_0.025'))

res_pair_full <- readRDS(file.path(directory, 'gw_4_big_lambda_line_pair_full_0.025'))
res_average_full <- readRDS(file.path(directory, 'gw_4_big_lambda_line_average_full_0.025'))
res_loh_full <- readRDS(file.path(directory, 'gw_4_big_lambda_line_loh_full_0.025'))


# utility function
f <- function(i){
  if(is.null(i)){
    rep(NA, 753)
  } else {
    i
  }
}

data_plot_full <- data.frame(x = 1 : nrow(res_pair_full$data),
                   y_pair_0 = res_pair_full$tree$gain,
                   y_pair_left = f(res_pair_full$tree$children[[1]]$gain), y_pair_right = f(res_pair_full$tree$children[[2]]$gain),
                   y_pair_left_left = f(res_pair_full$tree$children[[1]]$children[[1]]$gain), y_pair_left_right = f(res_pair_full$tree$children[[1]]$children[[2]]$gain),
                   y_pair_right_left =  f(res_pair_full$tree$children[[2]]$children[[1]]$gain), y_pair_right_right = f(res_pair_full$tree$children[[2]]$children[[2]]$gain),

                   y_loh_0 = f(res_loh_full$tree$gain),
                   y_loh_left = f(res_loh_full$tree$children[[1]]$gain), y_loh_right =  f(res_loh_full$tree$children[[2]]$gain),
                   y_loh_left_left = f(res_loh_full$tree$children[[1]]$children[[1]]$gain), y_loh_left_right = f(res_loh_full$tree$children[[1]]$children[[2]]$gain),
                   y_loh_right_left =  f(res_loh_full$tree$children[[2]]$children[[1]]$gain), y_loh_right_right = f(res_loh_full$tree$children[[2]]$children[[2]]$gain),

                   y_average_0 = f(res_average_full$tree$gain),
                   y_average_left = f(res_average_full$tree$children[[1]]$gain), y_average_right =  f(res_average_full$tree$children[[2]]$gain),
                   y_average_left_left = f(res_average_full$tree$children[[1]]$children[[1]]$gain), y_average_left_right = f(res_average_full$tree$children[[1]]$children[[2]]$gain),
                   y_average_right_left =  f(res_average_full$tree$children[[2]]$children[[1]]$gain), y_average_right_right = f(res_average_full$tree$children[[2]]$children[[2]]$gain)
)

data_plot_50 <- data.frame(x = 1 : nrow(res_pair_50$data),
                             y_pair_0 = f(res_pair_50$tree$gain),
                             y_pair_left = f(res_pair_50$tree$children[[1]]$gain), y_pair_right =  f(res_pair_50$tree$children[[2]]$gain),
                             y_pair_left_left = f(res_pair_50$tree$children[[1]]$children[[1]]$gain), y_pair_left_right = f(res_pair_50$tree$children[[1]]$children[[2]]$gain),
                             y_pair_right_left =  f(res_pair_50$tree$children[[2]]$children[[1]]$gain), y_pair_right_right = f(res_pair_50$tree$children[[2]]$children[[2]]$gain),

                             y_loh_0 = f(res_loh_50$tree$gain),
                             y_loh_left = f(res_loh_50$tree$children[[1]]$gain), y_loh_right =  f(res_loh_50$tree$children[[2]]$gain),
                             y_loh_left_left = f(res_loh_50$tree$children[[1]]$children[[1]]$gain), y_loh_left_right = f(res_loh_50$tree$children[[1]]$children[[2]]$gain),
                             y_loh_right_left =  f(res_loh_50$tree$children[[2]]$children[[1]]$gain), y_loh_right_right = f(res_loh_50$tree$children[[2]]$children[[2]]$gain),

                             y_average_0 = f(res_average_50$tree$gain),
                             y_average_left = f(res_average_50$tree$children[[1]]$gain), y_average_right =  f(res_average_50$tree$children[[2]]$gain),
                             y_average_left_left = f(res_average_50$tree$children[[1]]$children[[1]]$gain), y_average_left_right = f(res_average_50$tree$children[[1]]$children[[2]]$gain),
                             y_average_right_left =  f(res_average_50$tree$children[[2]]$children[[1]]$gain), y_average_right_right = f(res_average_50$tree$children[[2]]$children[[2]]$gain)
)

library(ggplot2)
m0 <- c(0,0,0,0)

first_splits_full <- c(
  res_pair_full$tree$split_point,
  res_average_full$tree$split_point,
  res_loh_full$tree$split_point
)

first_splits_50 <- c(
  res_pair_50$tree$split_point,
  res_average_50$tree$split_point,
  res_loh_50$tree$split_point
)

second_splits_50 <- c(
  res_pair_50$tree$children[[1]]$split_point,
  res_pair_50$tree$children[[2]]$split_point,
  res_average_50$tree$children[[1]]$split_point,
  res_average_50$tree$children[[2]]$split_point,
  res_loh_50$tree$children[[1]]$split_point,
  res_loh_50$tree$children[[2]]$split_point
)

second_splits_full <- c(
  res_pair_full$tree$children[[1]]$split_point,
  res_pair_full$tree$children[[2]]$split_point,
  res_average_full$tree$children[[1]]$split_point,
  res_average_full$tree$children[[2]]$split_point,
  res_loh_full$tree$children[[1]]$split_point,
  res_loh_full$tree$children[[2]]$split_point
)

third_splits_full <- c(
  res_pair_full$tree$children[[1]]$children[[1]]$split_point,
  res_pair_full$tree$children[[1]]$children[[2]]$split_point,
  res_pair_full$tree$children[[2]]$children[[1]]$split_point,
  res_pair_full$tree$children[[2]]$children[[2]]$split_point,
  
  res_average_full$tree$children[[1]]$children[[1]]$split_point,
  res_average_full$tree$children[[1]]$children[[2]]$split_point,
  res_average_full$tree$children[[2]]$children[[1]]$split_point,
  res_average_full$tree$children[[2]]$children[[2]]$split_point,
  
  res_loh_full$tree$children[[1]]$children[[1]]$split_point,
  res_loh_full$tree$children[[1]]$children[[2]]$split_point,
  res_loh_full$tree$children[[2]]$children[[1]]$split_point,
  res_loh_full$tree$children[[2]]$children[[2]]$split_point
)

third_splits_50 <- c(
  res_pair_50$tree$children[[1]]$children[[1]]$split_point,
  res_pair_50$tree$children[[1]]$children[[2]]$split_point,
  res_pair_50$tree$children[[2]]$children[[1]]$split_point,
  res_pair_50$tree$children[[2]]$children[[2]]$split_point,
  
  res_average_50$tree$children[[1]]$children[[1]]$split_point,
  res_average_50$tree$children[[1]]$children[[2]]$split_point,
  res_average_50$tree$children[[2]]$children[[1]]$split_point,
  res_average_50$tree$children[[2]]$children[[2]]$split_point,
  
  res_loh_50$tree$children[[1]]$children[[1]]$split_point,
  res_loh_50$tree$children[[1]]$children[[2]]$split_point,
  res_loh_50$tree$children[[2]]$children[[1]]$split_point,
  res_loh_50$tree$children[[2]]$children[[2]]$split_point
)

colors = c('green', 'red', 'blue')
linetypes <-  c('longdash', 'twodash', 'dashed')

p1_full <- ggplot(data_plot_full) +
  geom_line(aes(x = x, y = y_pair_0), color = 'green') +
  geom_line(aes(x = x, y = y_loh_0), color = 'blue') +
  geom_line(aes(x = x, y = y_average_0), color = 'red') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(5,5,0,0), 'pt')) +
  labs(y = 'gain') +
  scale_x_continuous(expand=c(0,0), breaks = c(1, 108 + 120 * 0 : 5), labels = as.character(c(1951, 1960 + 10 * 0 : 5)), limits = c(1, 753)) +
  geom_vline(xintercept = first_splits_full, color = colors, lty = linetypes)

p2_full <- ggplot(data_plot_full) +
  geom_line(aes(x = x, y = y_pair_left), color = 'green') +
  geom_line(aes(x = x, y = y_pair_right), color = 'green') +
  geom_line(aes(x = x, y = y_average_left), color = 'red')+
  geom_line(aes(x = x, y = y_average_right), color = 'red') +
  geom_line(aes(x = x, y = y_loh_left), color = 'blue')+
  geom_line(aes(x = x, y = y_loh_right), color = 'blue') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(5,0,0,0), 'pt')) +
  labs(y = 'gain') +
  scale_x_continuous(expand=c(0,0), breaks = c(1, 108 + 120 * 0 : 5), labels = as.character(c(1951, 1960 + 10 * 0 : 5))) +
  geom_vline(xintercept = second_splits_full, color = rep(colors, each = 2), lty = rep(linetypes, each = 2))

p3_full <-  ggplot(data_plot_full) +
  geom_line(aes(x = x, y = y_pair_left_left), color = 'green') +
  geom_line(aes(x = x, y = y_pair_left_right), color = 'green') +
  geom_line(aes(x = x, y = y_pair_right_left), color = 'green') +
  #geom_line(aes(x = x, y = y_pair_right_right), color = 'green') +
  geom_line(aes(x = x, y = y_average_left_left), color = 'red') +
  geom_line(aes(x = x, y = y_average_left_right), color = 'red') +
  geom_line(aes(x = x, y = y_average_right_left), color = 'red') +
  geom_line(aes(x = x, y = y_average_right_right), color = 'red') +
  geom_line(aes(x = x, y = y_loh_left_left), color = 'blue') +
  geom_line(aes(x = x, y = y_loh_left_right), color = 'blue') +
  geom_line(aes(x = x, y = y_loh_right_left), color = 'blue') +
  geom_line(aes(x = x, y = y_loh_right_right), color = 'blue') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(5,0,0,0), 'pt')) +
  labs(y = 'gain')+
  scale_x_continuous(expand=c(0,0), breaks = c(1, 108 + 120 * 0 : 5), labels = as.character(c(1951, 1960 + 10 * 0 : 5))) +
  geom_vline(xintercept = third_splits_full, color = rep(colors, each = 4), lty = rep(linetypes, each = 4))
 
p4_full <- plot_missingness_structure(res_pair_full$data, order = T) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = 'white'),
        axis.text.y = element_text(size = 15, color = 'white'),
        axis.ticks.y =  element_line(color = 'white'),
        plot.margin = unit(m0, 'pt')) +
  scale_x_continuous(expand = c(0,0), breaks = c(1, 108 + 120 * 0 : 5), labels = as.character(c(1951, 1960 + 10 * 0 : 5)))


p1_50 <- ggplot(data_plot_50) +
  geom_line(aes(x = x, y = y_pair_0), color = 'green') +
  geom_line(aes(x = x, y = y_loh_0), color = 'blue') +
  geom_line(aes(x = x, y = y_average_0), color = 'red') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.title.y = element_text(vjust = -5),  plot.margin = unit(c(5,0,0,0), 'pt')) +
  labs(y = 'gain') +
  scale_x_continuous(expand=c(0,0), breaks = c(1, 108 + 120 * 0 : 5), labels = as.character(c(1951, 1960 + 10 * 0 : 5))) +
  geom_vline(xintercept = first_splits_50, color = rep(colors, each = 1), lty = linetypes)

p2_50 <- ggplot(data_plot_50) +
  geom_line(aes(x = x, y = y_pair_left), color = 'green') +
  geom_line(aes(x = x, y = y_pair_right), color = 'green') +
  geom_line(aes(x = x, y = y_average_left), color = 'red')+
  geom_line(aes(x = x, y = y_average_right), color = 'red') +
  geom_line(aes(x = x, y = y_loh_left), color = 'blue')+
  geom_line(aes(x = x, y = y_loh_right), color = 'blue') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.title.y = element_text(vjust = -5), plot.margin = unit(m0, 'pt')) +
  labs(y = 'gain') +
  scale_x_continuous(expand=c(0,0), breaks = c(1, 108 + 120 * 0 : 5), labels = as.character(c(1951, 1960 + 10 * 0 : 5))) +
  geom_vline(xintercept = second_splits_50, color = rep(colors, each = 2), lty = rep(linetypes, each = 2))

p3_50 <-  ggplot(data_plot_50) +
  geom_line(aes(x = x, y = y_pair_left_left), color = 'green') +
  geom_line(aes(x = x, y = y_pair_left_right), color = 'green') +
  geom_line(aes(x = x, y = y_pair_right_left), color = 'green') +
  geom_line(aes(x = x, y = y_pair_right_right), color = 'green') +
  geom_line(aes(x = x, y = y_average_left_left), color = 'red') +
  geom_line(aes(x = x, y = y_average_left_right), color = 'red') +
  geom_line(aes(x = x, y = y_average_right_left), color = 'red') +
  geom_line(aes(x = x, y = y_average_right_right), color = 'red') +
  geom_line(aes(x = x, y = y_loh_left_left), color = 'blue') +
  geom_line(aes(x = x, y = y_loh_left_right), color = 'blue') +
  geom_line(aes(x = x, y = y_loh_right_left), color = 'blue') +
  geom_line(aes(x = x, y = y_loh_right_right), color = 'blue') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(vjust = -5), plot.margin = unit(m0, 'pt')) +
  labs(y = 'gain') +
  scale_x_continuous(expand=c(0,0), breaks = c(1, 108 + 120 * 0 : 5), labels = as.character(c(1951, 1960 + 10 * 0 : 5))) +
  geom_vline(xintercept = third_splits_50, color = rep(colors, each = 4), lty = rep(linetypes, each = 4))

p4_50 <- plot_missingness_structure(res_pair_50$data, order = T) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = 'white'),
        axis.text.y = element_text(size = 15, color = 'white'),
        axis.ticks.y =  element_line(color = 'white'),
        plot.margin = unit(m0, 'pt')) +
  scale_x_continuous(expand = c(0,0), breaks = c(1, 108 + 120 * 0 : 5), labels = as.character(c(1951, 1960 + 10 * 0 : 5)))


plot_grid(p1_50, p1_full, p2_50, p2_full, p3_50, p3_full, p4_50, p4_full, nrow = 4, align = 'vh')

