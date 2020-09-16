# install.packages("ggplot2")
# install.packages("TSP")
# install.packages("cowplot")

#' Copy of the function get_glasso_gain_function from rfcd package, but with additional parameter
#' smoothen. If FALSE, we ignore the steps in the accompaniying paper section 3.2
gain <- function(control, smoothen = TRUE){
  
  function(x, start, end, lambda) {
    
    if(is.null(lambda)){
      stop('Please supply a value for lambda to the glasso gain function')
    }
    
    control$n_obs <- nrow(x)
    
    fit_global <- hdcd::get_glasso_fit(x[(start + 1) : end, , drop = F], lambda = lambda, control = control)
    
    function(split_point, evaluate_all = FALSE){
      
      fit_left <- hdcd::get_glasso_fit(x[(start + 1) : split_point, , drop = F], lambda = lambda, control = control)
      fit_right <- hdcd::get_glasso_fit(x[(split_point + 1) : end, , drop = F], lambda = lambda, control = control)
      
      
    
   
      if(smoothen){
        (sum(hdcd::loglikelihood(x[(start + 1) : split_point, fit_left$inds, drop = F],
                           fit_global$mu[fit_left$inds, drop = F],
                           #fit_global$w[fit_left$inds, fit_left$inds, drop = F],
                           fit_global$wi[fit_left$inds, fit_left$inds, drop = F])) +
           sum(hdcd::loglikelihood(x[(split_point + 1) : end, fit_right$inds, drop = F],
                             fit_global$mu[fit_right$inds, drop = F],
                             #fit_global$w[fit_right$inds, fit_right$inds, drop = F],
                             fit_global$wi[fit_right$inds, fit_right$inds, drop = F])) -
           sum(hdcd::loglikelihood(x[(start + 1) : split_point, fit_left$inds, drop = F],
                             fit_left$mu[fit_left$inds, drop = F],
                             #fit_left$w[fit_left$inds, fit_left$inds, drop = F],
                             fit_left$wi[fit_left$inds, fit_left$inds, drop = F])) -
           sum(hdcd::loglikelihood(x[(split_point + 1) : end, fit_right$inds, drop = F],
                             fit_right$mu[fit_right$inds, drop = F],
                             #fit_right$w[fit_right$inds, fit_right$inds, drop = F],
                             fit_right$wi[fit_right$inds, fit_right$inds, drop = F]))) / nrow(x)
      } else {
        (sum(hdcd::loglikelihood(x[(start + 1) : end, fit_global$inds, drop = F],
                           fit_global$mu[fit_global$inds, drop = F],
                          # fit_global$w[fit_global$inds, fit_global$inds, drop = F],
                           fit_global$wi[fit_global$inds, fit_global$inds, drop = F])) -
           sum(hdcd::loglikelihood(x[(start + 1) : split_point, fit_left$inds, drop = F],
                             fit_left$mu[fit_left$inds, drop = F],
                          #   fit_left$w[fit_left$inds, fit_left$inds, drop = F],
                             fit_left$wi[fit_left$inds, fit_left$inds, drop = F])) -
           sum(hdcd::loglikelihood(x[(split_point + 1) : end, fit_right$inds, drop = F],
                             fit_right$mu[fit_right$inds, drop = F],
                          #   fit_right$w[fit_right$inds, fit_right$inds, drop = F],
                             fit_right$wi[fit_right$inds, fit_right$inds, drop = F]))) / nrow(x)
      }
    }
  }
}

set.seed(12)
x <- hdcd::simulate_from_model(model <- hdcd::create_model(200, 100, c(70), hdcd::RandomNetwork))
x_del <- hdcd::delete_values(x, 0.3, 'blockwise')
curve_complete <- sapply(5 : 195, gain(hdcd::hdcd_control())(x, start = 0, end = nrow(x), lambda = 2))
curve_loh_blockwise_smoothened <- sapply(5 : 195, gain(hdcd::hdcd_control(glasso_NA_method = 'loh'), smoothen = T)(x_del, start = 0, end = nrow(x), lambda = 2))
curve_loh_blockwise_non_smoothened <- sapply(5 : 195, gain(hdcd::hdcd_control(glasso_NA_method = 'loh'), smoothen = F)(x_del, start = 0, end = nrow(x), lambda = 2))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
red <- 'orange'
blue <- 'blue'

p1 <- ggplot2::ggplot(data = data.frame(x = 1 : 200,
                               curve_complete = c(rep(NA, 4), curve_complete, rep(NA, 5)),
                               curve_1 = c(rep(NA, 4), curve_loh_blockwise_non_smoothened, rep(NA, 5)),
                               curve_2 = c(rep(NA, 4), curve_loh_blockwise_smoothened, rep(NA, 5)))) +
  ggplot2::scale_y_continuous(breaks = c(), labels = c(), expand = c(0, 0.02)) +
  ggplot2::scale_x_continuous(expand = c(0,0)) +
  ggplot2::geom_line(ggplot2::aes(x = x, y = curve_complete), col = 'black') +
  ggplot2::geom_line(ggplot2::aes(x = x, y = curve_1), col = red) +
  ggplot2::geom_line(ggplot2::aes(x = x, y = curve_2), col = blue) +
  ggplot2::theme(plot.margin = unit(c(0, 5,0,0), 'pt'),
        axis.title.x = element_blank(), axis.text.x = element_text(color = 'white', size = 15),
        axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15)) +
  ggplot2::labs(y = 'gain') +
  ggplot2::geom_vline(xintercept = 70, linetype = 2)
set.seed(2) # The ordering of variables is (semi-)random
p2 <- hdcd::plot_missingness_structure(x_del, order = T) + 
  ggplot2::theme(plot.margin = unit(c(0, 5,0,0), 'pt'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  ggplot2::scale_x_continuous(expand = c(0,0), labels = c())


cowplot::plot_grid(p1, p2, labels = c('', ''), nrow = 2)


