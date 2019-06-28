# Main simulation file run on external cores

args <- commandArgs(FALSE)

# args <- c(0,0,0,'test', 500, 100, 'ChainNetwork', 'average', 'section', 'blockwise', '0.2', '2', '1', '4', '70', '120', '120', '190')
# args <- c(0,0,0,'test', 100, 100, 'RandomNetwork', 'pair', 'line', 'mcar', '0.3', '2', '1', '1', '100')
comment <- args[4]
n <- as.integer(args[5])
p <- as.integer(args[6])
Network <- args[7]
NA_method <-  args[8]
optimizer <- args[9]
deletion <-  args[10]
delete_fraction <- as.numeric(args[11])
seed <- as.integer(args[12])
nrep <- as.integer(args[13])
n_of_segments <- as.integer(args[14])
segment_lengths <- as.integer(args[14 + seq_len(n_of_segments)])

permute_segments <- TRUE

cat('Running Simulation with comment = ', comment, ', seed = ', seed,', Network Structure = ', Network, ', NA_method = ', NA_method, ', optimizer = ',optimizer, ', deletion method = ', deletion,
    ', deletion_fraction = ', delete_fraction, ', nrep = ', nrep, ', n = ', n, ', p = ', p, ' and ', n_of_segments, ' segments of length ', paste(segment_lengths, collapse = ' '),
    sep = '')


stopifnot(sum(segment_lengths) == n)

# load library
library('hdcdwithmissingvalues', lib.loc = '../R/x86_64-slackware-linux-gnu-library/')
library(data.table)

# evaluate function
Network_fun <- eval(parse(text = Network))

# if we delete no values, we scale down the segment lengths and total number of observations (for comparison)
if(deletion == 'none'){
  segment_lengths <- segment_lengths * (1 - delete_fraction)
  n <-  n *  (1 - delete_fraction)
}
# initialize dataset to save information later
data <- NULL

# repeat
for (i in 1 : nrep){
  set.seed(i + 10000 * seed)
  
  # if neccessary permute the order of the segment lengths
  if(permute_segments){
    segment_lengths <- segment_lengths[sample(1 : n_of_segments)]
  }
  alpha <- cumsum(segment_lengths)[-n_of_segments] # set of change points
  
  # draw training dataset
  x <- hdcdwithmissingvalues::simulate_from_model(hdcdwithmissingvalues::create_model(n, p, alpha, Network_fun))

  # delete values from training dataset
  x_del <- hdcdwithmissingvalues::delete_values(x, delete_fraction, deletion)
  
  tree <- hdcdwithmissingvalues::hdcd(x_del, method = 'glasso', lambda = 0.1, delta = 0.1, optimizer = optimizer,
             control = hdcdwithmissingvalues::hdcd_control(cv_inner = T, cv_inner_search_lambda = T, segment_loss_min_points = 5, section_search_min_points = 5, glasso_NA_method = NA_method))
  
  # save the true change points in a list
  true_cpts <- list(alpha)
  # extract the found change points from the tree
  found_cpts <- list(hdcdwithmissingvalues::get_change_points_from_tree(tree, variable = 'cv_improvement'))
  
  # get adj. Rand Index and other performance measures
  performance_measure <- hdcdwithmissingvalues::compare_change_points(alpha, found_cpts[[1]], n)
  
  data <- rbind(data,
                data.table::data.table(seed = 10000 * seed + i,
                                       Network = Network,
                                       NA_method = NA_method,
                                       optimizer = optimizer,
                                       deletion_method = deletion,
                                       delete_fraction = delete_fraction,
                                       n = n,
                                       p = p,
                                       n_of_segments = n_of_segments,
                                       segment_lengths = list(segment_lengths[order(segment_lengths)]),
                                       true_cpts = true_cpts,
                                       found_cpts = found_cpts,
                                       adj_Rand = performance_measure[2]
                ))
}


saveRDS(data, file = paste('simulation_output', paste(args[-(1:3)], collapse = '_'), sep = '_'))
