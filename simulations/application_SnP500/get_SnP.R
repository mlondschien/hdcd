library(quantmod)
library(data.table)
library(timeDate)

library(rvest)

# get list of SP500 tickers
url <- "https://en.wikipedia.org/wiki/List_of_S%26P_500_companies"
SP500 <- url %>%
  html() %>%
  html_nodes(xpath='//*[@id="mw-content-text"]/div/table[1]') %>%
  html_table()

names <- SP500[[1]]$Symbol
names <- gsub('\\.', '', names) # remove dots

# download returns until ...
start_date <- "1990-01-01"

#' get SnP500 data
#' 
#' Loads returns of SnP500 stocks with tickers in \code{name} from \code{start_date} until
#' \code{end_date}.
#' 
get_SnP <- function(names, start_date, end_date = NULL){

  if (is.null(end_date)) end_date <- substr(Sys.time(), start = 1, stop = 10)

  dates <- timeSequence(as.Date(start_date), as.Date(end_date))
  # dates <- dates[isWeekday(dates)]
  data <- data.table(date = as.character(dates))
  fail <- NULL

  for (ticker in names){
    cat("Downloading time series for symbol '", ticker, "' ...\n", sep = "")
    tryCatch({
      S <- as.data.table(getSymbols(ticker, auto.assign=F, from = as.Date(start_date)))[, 1 : 2]
      S[, index := as.character(index)]
      data <- merge(data, S, by.x = 'date', by.y = 'index', all.x = T)},
      error = function(e){
        cat('Failed to download time series for symbol ', ticker, '...\n', sep = '')
        fail <- c(fail, ticker)})
  }
  if(length(fail) >0) cat('Failed to download time series for symbols', fail, sep = ' ' )
  data
}

start_date <- "1990-01-01"
SnP500 <- get_SnP(names, start_date)

# save
write.csv2(SnP500, file = 'SnP500.csv', row.names = FALSE)

# remove weekends
weekend_id <- apply(is.na(SnP500[, -1]), 1, all)
SnP500 <- SnP500[!weekend_id, ]

# we are interested in the logreturns of the stocks
log_returns <- data.frame(date = SnP500[-1, date], log(SnP500[-1, -1]) - log(SnP500[-nrow(SnP500), -1]))

# Remove entries where the return was greater than 2 or smaller than 1/2, as these come from untypical results
# and distort the sassumption of the data following a multivariate gaussian
log_returns[, -1][!is.na(log_returns[, -1]) & abs(log_returns[, -1]) > log(2)] <- 0

# some quick check: how many entries lie 10 sd's away from the mean.
# Note that this can be a high amount since we only assume the time series to be piecewise gaussian
n <- sum(!is.na(log_returns[, -1]))
sigma <- apply(log_returns[, -1], 1, sd, na.rm = T)
mu <- apply(log_returns[, -1], 1, mean, na.rm = T)
sum(abs(log_returns[, -1] - mu) > 10 * sigma, na.rm = T) / n

write.csv2(log_returns, 'SnP500_logreturns.csv', row.names = F)