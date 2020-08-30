#' Copy of the function get_glasso_gain_function from rfcd package, but with additional parameter
#' smoothen. If FALSE, we ignore the steps in the accompaniying paper section 3.2
gain <- function(control, smoothen = TRUE){
  
  function(x, start, end, lambda) {
    
    if(is.null(lambda)){
      stop('Please supply a value for lambda to the glasso gain function')
    }
    
    control$n_obs <- nrow(x)
    
    fit_global <- get_glasso_fit(x[(start + 1) : end, , drop = F], lambda = lambda, control = control)
    
    function(split_point, evaluate_all = FALSE){
      
      fit_left <- get_glasso_fit(x[(start + 1) : split_point, , drop = F], lambda = lambda, control = control)
      fit_right <- get_glasso_fit(x[(split_point + 1) : end, , drop = F], lambda = lambda, control = control)
      
      
    
   
      if(smoothen){
        (sum(loglikelihood(x[(start + 1) : split_point, fit_left$inds, drop = F],
                           fit_global$mu[fit_left$inds, drop = F],
                           #fit_global$w[fit_left$inds, fit_left$inds, drop = F],
                           fit_global$wi[fit_left$inds, fit_left$inds, drop = F])) +
           sum(loglikelihood(x[(split_point + 1) : end, fit_right$inds, drop = F],
                             fit_global$mu[fit_right$inds, drop = F],
                             #fit_global$w[fit_right$inds, fit_right$inds, drop = F],
                             fit_global$wi[fit_right$inds, fit_right$inds, drop = F])) -
           sum(loglikelihood(x[(start + 1) : split_point, fit_left$inds, drop = F],
                             fit_left$mu[fit_left$inds, drop = F],
                             #fit_left$w[fit_left$inds, fit_left$inds, drop = F],
                             fit_left$wi[fit_left$inds, fit_left$inds, drop = F])) -
           sum(loglikelihood(x[(split_point + 1) : end, fit_right$inds, drop = F],
                             fit_right$mu[fit_right$inds, drop = F],
                             #fit_right$w[fit_right$inds, fit_right$inds, drop = F],
                             fit_right$wi[fit_right$inds, fit_right$inds, drop = F]))) / nrow(x)
      } else {
        (sum(loglikelihood(x[(start + 1) : end, fit_global$inds, drop = F],
                           fit_global$mu[fit_global$inds, drop = F],
                          # fit_global$w[fit_global$inds, fit_global$inds, drop = F],
                           fit_global$wi[fit_global$inds, fit_global$inds, drop = F])) -
           sum(loglikelihood(x[(start + 1) : split_point, fit_left$inds, drop = F],
                             fit_left$mu[fit_left$inds, drop = F],
                          #   fit_left$w[fit_left$inds, fit_left$inds, drop = F],
                             fit_left$wi[fit_left$inds, fit_left$inds, drop = F])) -
           sum(loglikelihood(x[(split_point + 1) : end, fit_right$inds, drop = F],
                             fit_right$mu[fit_right$inds, drop = F],
                          #   fit_right$w[fit_right$inds, fit_right$inds, drop = F],
                             fit_right$wi[fit_right$inds, fit_right$inds, drop = F]))) / nrow(x)
      }
    }
  }
}

library('hdcdwithmissingvalues')
library(data.table)

set.seed(12)
x <- simulate_from_model(model <- create_model(200, 100, c(70), RandomNetwork))
x_del <- delete_values(x, 0.3, 'blockwise')
curve_complete <- sapply(5 : 195, gain(hdcd_control())(x, start = 0, end = nrow(x), lambda = 2))
curve_loh_blockwise_smoothened <- sapply(5 : 195, gain(hdcd_control(glasso_NA_method = 'loh', segment_loss_min_points = 5), smoothen = T)(x_del, start = 0, end = nrow(x), lambda = 2))
curve_loh_blockwise_non_smoothened <- sapply(5 : 195, gain(hdcd_control(glasso_NA_method = 'loh', segment_loss_min_points = 5), smoothen = F)(x_del, start = 0, end = nrow(x), lambda = 2))

library(cowplot)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
red <- 'orange'
blue <- 'blue'

p1 <- ggplot(data = data.frame(x = 1 : 200,
                               curve_complete = c(rep(NA, 4), curve_complete, rep(NA, 5)),
                               curve_1 = c(rep(NA, 4), curve_loh_blockwise_non_smoothened, rep(NA, 5)),
                               curve_2 = c(rep(NA, 4), curve_loh_blockwise_smoothened, rep(NA, 5)))) +
  scale_y_continuous(breaks = c(), labels = c(), expand = c(0, 0.02)) +
  scale_x_continuous(expand = c(0,0)) +
  geom_line(aes(x = x, y = curve_complete), col = 'black') +
  geom_line(aes(x = x, y = curve_1), col = red) +
  geom_line(aes(x = x, y = curve_2), col = blue) +
  theme(plot.margin = unit(c(0, 5,0,0), 'pt'),
        axis.title.x = element_blank(), axis.text.x = element_text(color = 'white', size = 15),
        axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15)) +
  labs(y = 'gain') +
  geom_vline(xintercept = 70, linetype = 2)
set.seed(2) # Ordering of variables random
p2 <- plot_missingness_structure(x_del, order = T) + 
  theme(plot.margin = unit(c(0, 5,0,0), 'pt'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  scale_x_continuous(expand = c(0,0), labels = c())

  
  
plot_grid(p1, p2, labels = c('', ''), nrow = 2)


