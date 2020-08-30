library(hdcdwithmissingvalues)

set.seed(1)
x <- simulate_from_model(model <- create_model(500, 100, c(120, 240, 310), RandomNetwork))
x_del_blockwise <- delete_values(x, 0.3, 'blockwise')
x_del_mcar <- delete_values(x, 0.3, 'mcar')

lambda <- 0.1
curve_complete <- sapply(5 : 495, get_glasso_gain_function(hdcd_control(segment_loss_min_points = 5, n_obs = 500))(x, start = 0, end = 500, lambda = lambda))

curve_loh_mcar <- sapply(5 : 495, get_glasso_gain_function(hdcd_control(glasso_NA_method = 'loh', segment_loss_min_points = 5, n_obs = 500))(x_del_mcar, start = 0, end = 500, lambda = lambda))
curve_pair_mcar <- sapply(5 : 495, get_glasso_gain_function(hdcd_control(glasso_NA_method = 'pair', segment_loss_min_points = 5, n_obs = 500))(x_del_mcar, start = 0, end = 500, lambda = lambda))
curve_av_mcar <- sapply(5 : 495, get_glasso_gain_function(hdcd_control(glasso_NA_method = 'av', segment_loss_min_points = 5, n_obs = 500))(x_del_mcar, start = 0, end = 500, lambda = lambda))

curve_loh_blockwise <- sapply(5 : 495, get_glasso_gain_function(hdcd_control(glasso_NA_method = 'loh', segment_loss_min_points = 5, n_obs = 500))(x_del_blockwise, start = 0, end = 500, lambda = lambda))
curve_pair_blockwise <- sapply(5 : 495, get_glasso_gain_function(hdcd_control(glasso_NA_method = 'pair', segment_loss_min_points = 5, n_obs = 500))(x_del_blockwise, start = 0, end = 500, lambda = lambda))
curve_av_blockwise <- sapply(5 : 495, get_glasso_gain_function(hdcd_control(glasso_NA_method = 'av', segment_loss_min_points = 5, n_obs = 500))(x_del_blockwise, start = 0, end = 500, lambda = lambda))
                          

library(cowplot)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

red <- 'red' # cbbPalette[7]
green <- 'green' # cbbPalette[4]
blue <- 'blue' # cbbPalette[3]

p1 <- ggplot(data = data.frame(x = 5 : 495,
                               curve_complete = curve_complete,
                               av = curve_av_mcar,
                               loh = curve_loh_mcar,
                               pair = curve_pair_mcar)) +
  scale_y_continuous(breaks = 0 : 2, labels = 0 : 2, expand = c(0, 0.02)) +
  geom_line(aes(x = x, y = curve_complete), col = 'black') +
  #geom_line(aes(x = x, y = curve_true), col = 'black', lty = 'dashed') +
  geom_line(aes(x = x, y = av), col = red) +
  geom_line(aes(x = x, y = loh), col = blue) +
  geom_line(aes(x = x, y = pair), col = green) +
  theme(plot.margin = unit(c(0, 5,0,0), 'pt'),
        axis.title.x = element_blank(), axis.text.x = element_text(color = 'white', size = 15),
        axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15)) +
  labs(y = 'gain')

p2 <- ggplot(data = data.frame(x = 5 : 495,
                               curve_complete = curve_complete,
                               av = curve_av_blockwise,
                               loh = curve_loh_blockwise,
                               pair = curve_pair_blockwise)) +
  scale_y_continuous(breaks = 0 : 2, labels = 0 : 2, expand = c(0, 0.02)) +
  geom_line(aes(x = x, y = curve_complete), col = 'black') +
  #geom_line(aes(x = x, y = curve_true), col = 'black', lty = 'dashed') +
  geom_line(aes(x = x, y = av), col = red) +
  geom_line(aes(x = x, y = loh), col = blue) +
  geom_line(aes(x = x, y = pair), col = green) +
  theme(plot.margin = unit(c(0, 5,0,0), 'pt'),
        axis.title.x = element_blank(), axis.text.x = element_text(color = 'white', size = 15),
        axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15)) +
  labs(y = 'gain')

library(cowplot)
plot_grid(p1, p2, nrow = 2)
