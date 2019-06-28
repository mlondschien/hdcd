args <- commandArgs(FALSE)
#args <- c(0, 0, 0, 'test', 'sect', 'pair', '0.05', '0.2', '0')
comment <- args[4]
optimizer <- args[5]
NA_method <- args[6]
delta <- as.numeric(args[7])
section_search_stepsize <-  as.numeric(args[8])

x <- read.csv2('SnP500_logreturns.csv')
#x <- read.csv2('./simulations/application_SnP500/SnP500_logreturns.csv')
# for numerical stability
dt <- as.matrix(x[, -1]) * 100

dt <- dt[, apply(!is.na(dt), 2, all)]

library('hdcdwithmissingvalues', lib.loc = '../R/x86_64-slackware-linux-gnu-library/')

# get an initial guess for lambda. This will be selected using inner cv anyway. Better choose a small lambda for
# initialization.
cov_mat <- cov(as.matrix(dt), use = 'pairwise')
lambda <- max(abs(cov_mat[upper.tri(cov_mat)]), na.rm  = T)

res <- hdcd(dt[], method = 'glasso', lambda = lambda, optimizer = optimizer, delta = delta,
            control = hdcd_control(glasso_NA_method = NA_method, cv_inner = T,  cv_inner_search_lambda = T, segment_loss_min_points = 10, section_search_min_points = 20,
                                   section_search_stepsize = section_search_stepsize))

saveRDS(res, file = paste(args[-c(1,2,3)], collapse = '_'))
