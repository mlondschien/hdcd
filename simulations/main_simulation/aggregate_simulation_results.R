# Read outputs of single simulations and combine into one data.table

name <- 'main_simulation'
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

# cast adjusted Rand Indices into more managable form
dcast(
  data[,
       paste(format(round(mean(adj_Rand), digits = 3), nsmall = 3),
             ' (', format(round(sd(adj_Rand), digits = 3), nsmall = 3), ')',
             sep = ''),
       by = c('Network', 'NA_method', 'deletion_method', 'delete_fraction', 'optimizer')],
  Network + deletion_method +NA_method + optimizer  ~ delete_fraction
)

# prepare generation of output for latex table
library(xtable)
data_table <- data
data_table[Network == 'RandomNetwork', NW := 'RN']
data_table[Network == 'ChainNetwork', NW := 'CN']
data_table[deletion_method == 'blockwise', MISSINGNESS := 'block']
data_table[deletion_method == 'mcar', MISSINGNESS := 'MCAR']
data_table[optimizer == 'line', OPT := 'BS']
data_table[optimizer == 'section', OPT := 'OBS']
data_table[NA_method == 'loh', NA_method := 'LW']

# Print the main table
print(xtable(
  dcast(
    data_table[!is.na(MISSINGNESS),
               paste(format(round(mean(adj_Rand), digits = 3), nsmall = 3),
                     ' (', format(round(sd(adj_Rand), digits = 3), nsmall = 3), ')',
                     sep = ''),
               by = c('NW', 'OPT', 'MISSINGNESS', 'NA_method', 'delete_fraction')],
    NW + MISSINGNESS+  NA_method +  OPT ~ delete_fraction
  )
), include.rownames=FALSE)

# Print the table with the reduced amount of observations for comparison
print(xtable(
  dcast(
    data_table[deletion_method == 'none',
               paste(format(round(mean(adj_Rand), digits = 3), nsmall = 3),
                     ' (', format(round(sd(adj_Rand), digits = 3), nsmall = 3), ')',
                     sep = ''),
               by = c('NW', 'OPT', 'delete_fraction')],
    NW +  OPT ~ delete_fraction
  )
), include.rownames=FALSE)
