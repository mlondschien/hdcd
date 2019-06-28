data <- read.csv2('./simulations/application_SnP500/SnP500_logreturns.csv')
SNP_loh <- readRDS('./simulations/application_SnP500/snp_2_section_loh_0.025_0.01')
SNP_pair <- readRDS('./simulations/application_SnP500/snp_2_section_pair_0.025_0.01')
SNP_average <- readRDS('./simulations/application_SnP500/snp_2_section_average_0.025_0.01')

SNP_loh <- readRDS('./simulations/application_SnP500/snp_3_comp_section_loh_0.025_0.01')
SNP_pair <- readRDS('./simulations/application_SnP500/snp_3_comp_section_pair_0.025_0.01')
SNP_average <- readRDS('./simulations/application_SnP500/snp_3_comp_section_average_0.025_0.01')


conversion <- data.frame(
  id = get_change_points_from_tree(SNP_loh),
  date = data[hdcdwithmissingvalues::get_change_points_from_tree(SNP_loh), 1]
)

library(hdcdwithmissingvalues)
library(cowplot)

yearly_changes <- c(1, which(substr(as.character(data[, 1]), start = 1, stop = 4) != c(substr(as.character(data[, 1]), start = 1, stop = 4)[-1], '2019')))

f <- function(i){
  if(is.null(i)){
    rep(NA, nrow(data))
  } else {
    i
  }
}

first_splits <- c(
  SNP_pair$split_point,
  SNP_average$split_point,
  SNP_loh$split_point
)

second_splits <- c(
  SNP_pair$children[[1]]$split_point,
  SNP_pair$children[[2]]$split_point,
  SNP_average$children[[1]]$split_point,
  SNP_average$children[[2]]$split_point,
  SNP_loh$children[[1]]$split_point,
  SNP_loh$children[[2]]$split_point
)

third_splits <- c(
  SNP_pair$children[[1]]$children[[1]]$split_point,
  SNP_pair$children[[1]]$children[[2]]$split_point,
  SNP_pair$children[[2]]$children[[1]]$split_point,
  SNP_pair$children[[2]]$children[[2]]$split_point,
  
  SNP_average$children[[1]]$children[[1]]$split_point,
  SNP_average$children[[1]]$children[[2]]$split_point,
  SNP_average$children[[2]]$children[[1]]$split_point,
  SNP_average$children[[2]]$children[[2]]$split_point,
  
  SNP_loh$children[[1]]$children[[1]]$split_point,
  SNP_loh$children[[1]]$children[[2]]$split_point,
  SNP_loh$children[[2]]$children[[1]]$split_point,
  SNP_loh$children[[2]]$children[[2]]$split_point
)


data_plot_SNP <- data.frame(x = 1 : nrow(data),
                             y_pair_0 = SNP_pair$gain,
                             y_pair_left = f(SNP_pair$children[[1]]$gain), y_pair_right = f(SNP_pair$children[[2]]$gain),
                             y_pair_left_left = f(SNP_pair$children[[1]]$children[[1]]$gain), y_pair_left_right = f(SNP_pair$children[[1]]$children[[2]]$gain),
                             y_pair_right_left =  f(SNP_pair$children[[2]]$children[[1]]$gain), y_pair_right_right = f(SNP_pair$children[[2]]$children[[2]]$gain),
                             
                             y_loh_0 = f(SNP_loh$gain),
                             y_loh_left = f(SNP_loh$children[[1]]$gain), y_loh_right =  f(SNP_loh$children[[2]]$gain),
                             y_loh_left_left = f(SNP_loh$children[[1]]$children[[1]]$gain), y_loh_left_right = f(SNP_loh$children[[1]]$children[[2]]$gain),
                             y_loh_right_left =  f(SNP_loh$children[[2]]$children[[1]]$gain), y_loh_right_right = f(SNP_loh$children[[2]]$children[[2]]$gain),
                             
                             y_average_0 = f(SNP_average$gain),
                             y_average_left = f(SNP_average$children[[1]]$gain), y_average_right =  f(SNP_average$children[[2]]$gain),
                             y_average_left_left = f(SNP_average$children[[1]]$children[[1]]$gain), y_average_left_right = f(SNP_average$children[[1]]$children[[2]]$gain),
                             y_average_right_left =  f(SNP_average$children[[2]]$children[[1]]$gain), y_average_right_right = f(SNP_average$children[[2]]$children[[2]]$gain)
)

colors = c('green', 'red', 'blue')
linetypes <-  c('longdash', 'twodash', 'dashed')


p1 <- ggplot() +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_pair_0), ], aes(x = x, y = y_pair_0), color = 'green') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_loh_0), ], aes(x = x, y = y_loh_0), color = 'blue') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_average_0), ], aes(x = x, y = y_average_0), color = 'red') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(5,5,0,0), 'pt')) +
  labs(y = 'gain') +
  scale_x_continuous(breaks = yearly_changes, labels = 1990 : 2019) +
  geom_vline(xintercept = first_splits, color = colors, lty = linetypes)

p2 <- ggplot(data_plot_SNP) +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_pair_left), ], aes(x = x, y = y_pair_left), color = 'green') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_pair_right), ],aes(x = x, y = y_pair_right), color = 'green') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_average_left), ],aes(x = x, y = y_average_left), color = 'red')+
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_average_right), ],aes(x = x, y = y_average_right), color = 'red') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_loh_left), ],aes(x = x, y = y_loh_left), color = 'blue')+
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_loh_right), ],aes(x = x, y = y_loh_right), color = 'blue') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(5,0,0,0), 'pt')) +
  labs(y = 'gain') +
  scale_x_continuous(breaks = yearly_changes, labels = 1990 : 2019) +
  geom_vline(xintercept = second_splits, color = rep(colors, each = 2), lty = rep(linetypes, each = 2))

p3 <-  ggplot(data_plot_SNP) +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_pair_left_left), ], aes(x = x, y = y_pair_left_left), color = 'green') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_pair_left_right), ],aes(x = x, y = y_pair_left_right), color = 'green') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_pair_right_left), ], aes(x = x, y = y_pair_right_left), color = 'green') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_pair_right_right), ],aes(x = x, y = y_pair_right_right), color = 'green') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_average_left_left), ], aes(x = x, y = y_average_left_left), color = 'red') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_average_left_right), ],aes(x = x, y = y_average_left_right), color = 'red') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_average_right_left), ], aes(x = x, y = y_average_right_left), color = 'red') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_average_right_right), ],aes(x = x, y = y_average_right_right), color = 'red') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_loh_left_left), ], aes(x = x, y = y_loh_left_left), color = 'blue') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_loh_left_right), ],aes(x = x, y = y_loh_left_right), color = 'blue') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_loh_right_left), ], aes(x = x, y = y_loh_right_left), color = 'blue') +
  geom_line(data = data_plot_SNP[!is.na(data_plot_SNP$y_loh_right_right), ],aes(x = x, y = y_loh_right_right), color = 'blue') +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(5,0,0,0), 'pt')) +
  labs(y = 'gain')+
  scale_x_continuous(breaks = yearly_changes, labels = 1990 : 2019) +
  geom_vline(xintercept = third_splits, color = rep(colors, each = 4), lty = rep(linetypes, each = 4))

p4 <- plot_missingness_structure(data, order = T) +
  theme(legend.position = "none",
        axis.title.x = element_text(color = 'black', size = 15),
        axis.title.y = element_text(color = 'white'),
        axis.text.y = element_text(size = 15, color = 'white'),
        axis.ticks.y =  element_line(color = 'white'),
        plot.margin = unit(c(5,0,0,0), 'pt')) +
        scale_x_continuous(breaks = yearly_changes, labels = 1990 : 2019)

library(cowplot)
plot_grid(p1, p2, p3, p4, nrow = 4, align = 'vh')




# 2009-11-09
# 2003-06-24

FUN <- function(node) {
  node <- data.tree::Navigate(node, "..") # We want to prune the parent tree where the split would've occured

  if (is.null(node)) {
    TRUE
  } else {
    node[['cv_train_improvement']] - 2 * node[['cv_train_improvement_sd']] > 0
  }
}

GetSignificantChangepoints <- function(x) {
  stopifnot(is(x, "bs_tree"))
  GetChangePointsFromLeafs(data.tree::Clone(x, pruneFun = FUN))
}


library(ggplot2)
library(scales)
ggplot(data = data.frame(x = as.Date(SnP500_logreturns[, 1]), y = res_SNP_line$gain)) +
  geom_line(aes(x = x, y = y)) +
  scale_x_date(labels = date_format("%m-%Y")) +
  geom_point(data = data.frame(x = as.Date(SnP500_logreturns[GetSignificantChangepoints(res_SNP_OBS_loh), 1]),
                           y = res_SNP_line$gain[GetSignificantChangepoints(res_SNP_OBS_loh)]),
              aes(x = x, y = y)) +
  geom_text(data = data.frame(x = as.Date(SnP500_logreturns[GetSignificantChangepoints(res_SNP_OBS_loh), 1]),
                              y = res_SNP_line$gain[GetSignificantChangepoints(res_SNP_OBS_loh)]),
            aes(x = x, y = y, label = x))



plot(res_SNP_line$gain, type = 'n', xlab = SnP500_logreturns[, 1])
lines(res_SNP_line$gain)
