# Read outputs of single simulations and combine into one data.table

name <- 'simulation_overfitting_1'
directory <- file.path('./simulations', name)
library(data.table)

data <- NULL
i <- 0
fail_files <- NULL
for(file in list.files(path = directory)[]){ 
  i <- i + 1
  tryCatch({
    cat('Reading file ', i, ' / ', length(list.files(path = directory)), ' : ', file, '...', sep = '')
    dt <- readRDS(file.path(directory, file))
    cat('Done \n')
    data <- rbind(data, data.table(dt[,], file = file))
    },
  error = function(e){
    cat('\n Error: Could not read file', file, '\n', sep = ' ')
  })
}

# save into RDS file
# saveRDS(data, file = file.path(directory, paste('data_', name, sep = '')))

# print overfitting behaviour, First number is proportion of simulations without any found 
# change points, second with one and third with more than one found change point
dcast(
  data[,
       paste(
         format(round(mean(0 == lapply(found_cpts, length)),2), nsmall = 2),
         format(round(mean(1 == lapply(found_cpts, length)),2), nsmall = 2),
         format(round(mean(1 < lapply(found_cpts, length)),2), nsmall = 2)
       ),
       by = c('Network', 'NA_method', 'deletion_method', 'delete_fraction', 'optimizer')],
  Network + NA_method + optimizer + deletion_method ~ delete_fraction
)

data[, n_cpts_found := unlist(lapply(found_cpts, length))]
data[, sum(n_cpts_found), by = .(NA_method, deletion_method)] # none found for MCAR missingness
data[deletion_method == 'blockwise', sum(n_cpts_found), by = .(NA_method, delete_fraction)]