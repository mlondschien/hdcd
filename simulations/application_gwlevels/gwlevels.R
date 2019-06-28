args <- commandArgs(FALSE)

comment <- args[4]
optimizer <- args[5]
NA_method <- args[6]
part <- args[7]
delta <- as.numeric(args[8])

dt <- read.csv('groundwater_levels.csv')
dt <- read.csv('./simulations/application_gwlevels/groundwater_levels.csv')

library('hdcdwithmissingvalues', lib.loc = '../R/x86_64-slackware-linux-gnu-library/')
#library(hdcdwithmissingvalues)

dt_processed <- data.frame(dt[, -c(22,1)]) #delete time and totally missing values

LogSpace <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}

# we remove the seasonality (yearly and monthly trend) for each column seperately
for (j in 1 : ncol(dt_processed)){
  fit <- lm(y ~ trend_yearly + trend_total, data= data.frame(y = dt_processed[, j], trend_total = 1 : nrow(dt), trend_yearly = as.factor(rep(1 : 12, ceiling(nrow(dt_processed) / 12))[1 : nrow(dt_processed)]) ))
  dt_processed[!is.na(dt_processed[, j]), j] <- fit$residuals * 10
}

# possibly only consider variables with more than 50% of observations available
if(part == 'full'){
  dt_processed <- dt_processed
} else if (part == 'cut'){
  dt_processed <- dt_processed[-50 : - 1, ]
} else if (part == '50'){
  available_obs <- apply(dt_processed, 2, function(z) sum(!is.na(z)))
  dt_processed <- dt_processed[, available_obs > nrow(dt_processed) / 2]
}

# lambda for initialisation
cov_mat <- cov(as.matrix(dt_processed), use = 'pairwise')
lambda <- max(abs(cov_mat[upper.tri(cov_mat)]), na.rm  = T)
#lambda <- log_space(1, 10, 10)

res <- list(
  tree = hdcd(dt_processed, lambda = lambda[1], method = 'glasso', optimizer = optimizer,  delta = delta, control = hdcd_control(
    cv_inner_search_lambda = TRUE, cv_inner = T, cv_inner_lambda = lambda, segment_loss_min_points = 5, section_search_min_points = 10, glasso_NA_method = NA_method)),
  data = dt_processed,
  part = part,
  comment = comment,
  NA_method = NA_method,
  delta = delta,
  optimizer = optimizer
)

saveRDS(res, paste(args[-c(1,2,3)], collapse = '_'))
