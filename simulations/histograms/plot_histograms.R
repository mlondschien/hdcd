library(ggplot2)

data <- data.table::data.table(readRDS('./simulations/histograms/data_histograms'))

### PLOT HISTOGRAMS
NW1 <- 'RandomNetwork'
missingness1 <- 'blockwise'
frac1 <- 0.1
frac2 <- 0.2
frac3 <- 0.3
frac4 <- 0.4
frac5 <- 0.5
method1 <- 'loh'
optimizer1 <- 'section'
true_cpts1 <- c(120, 240, 310)

data_hist1 <-  data.frame(changepoints = data[Network == NW1 & NA_method == method1 & optimizer == optimizer1 & deletion_method == missingness1 & delete_fraction == frac1, unlist(found_cpts)])
rand1 <- mean(data[Network == NW1 & NA_method == method1 & optimizer == optimizer1 & deletion_method == missingness1 & delete_fraction == frac1, adj_Rand])
data_hist2 <-  data.frame(changepoints = data[Network == NW1 & NA_method == method1 & optimizer == optimizer1 & deletion_method == missingness1 & delete_fraction == frac2, unlist(found_cpts)])
rand2 <- mean(data[Network == NW1 & NA_method == method1 & optimizer == optimizer1 & deletion_method == missingness1 & delete_fraction == frac2, adj_Rand])
data_hist3 <-  data.frame(changepoints = data[Network == NW1 & NA_method == method1 & optimizer == optimizer1 & deletion_method == missingness1 & delete_fraction == frac3, unlist(found_cpts)])
rand3 <- mean(data[Network == NW1 & NA_method == method1 & optimizer == optimizer1 & deletion_method == missingness1 & delete_fraction == frac3, adj_Rand])
data_hist4 <-  data.frame(changepoints = data[Network == NW1 & NA_method == method1 & optimizer == optimizer1 & deletion_method == missingness1 & delete_fraction == frac4, unlist(found_cpts)])
rand4 <- mean(data[Network == NW1 & NA_method == method1 & optimizer == optimizer1 & deletion_method == missingness1 & delete_fraction == frac4, adj_Rand])
data_hist5 <-  data.frame(changepoints = data[Network == NW1 & NA_method == method1 & optimizer == optimizer1 & deletion_method == missingness1 & delete_fraction == frac5, unlist(found_cpts)])
rand5 <- mean(data[Network == NW1 & NA_method == method1 & optimizer == optimizer1 & deletion_method == missingness1 & delete_fraction == frac5, adj_Rand])


trans <- 'log10'

library(cowplot)
y_max <- 400
p1 <- ggplot(data_hist1) +
  geom_histogram(aes(x = changepoints), binwidth = 1) +
  scale_y_continuous(expand = c(0,0), trans=trans, limits = c(1, 400)) +
  xlim(1, 500) +theme( axis.title.x = element_blank(),axis.text.x =element_text(color = 'white'), axis.title.y = element_blank()) +
  geom_text(data = data.frame(x = 430, y = 100, label = paste('adj. Rand = ', round(rand1, 3))), aes(x = x, y = y, label = label))

p2 <- ggplot(data_hist2) +
  geom_histogram(aes(x = changepoints), binwidth = 1) +
  scale_y_continuous(expand = c(0,0), trans=trans, limits = c(1, 400)) +
  xlim(1, 500) +theme(axis.title.x = element_blank(),axis.text.x = element_text(color = 'white'), axis.title.y = element_blank()) +
  geom_text(data = data.frame(x = 430, y = 100, label = paste('adj. Rand = ', round(rand2, 3))), aes(x = x, y = y, label = label))


p3 <- ggplot(data_hist3) +
  geom_histogram(aes(x = changepoints), binwidth = 1) +
  scale_y_continuous(expand = c(0,0), trans=trans, limits = c(1, 400)) +
  xlim(1, 500) +theme(axis.title.x = element_blank(),axis.text.x = element_text(color = 'white'), axis.title.y = element_blank()) +
  geom_text(data = data.frame(x = 430, y = 100, label = paste('adj. Rand = ', round(rand3, 3))), aes(x = x, y = y, label = label))


p4 <- ggplot(data_hist4) +
  geom_histogram(aes(x = changepoints), binwidth = 1) +
  scale_y_continuous(expand = c(0,0),trans=trans, limits = c(1, 400)) +
  xlim(1, 500) +theme(axis.title.x = element_blank(),axis.text.x = element_text(color = 'white'), axis.title.y = element_blank()) +
  geom_text(data = data.frame(x = 430, y = 100, label = paste('adj. Rand = ', round(rand4, 3))), aes(x = x, y = y, label = label))


p5 <- ggplot(data_hist5) +
  geom_histogram(aes(x = changepoints), binwidth = 1) +
  scale_y_continuous(expand = c(0,0), trans = trans, limits = c(1, 400)) +
  xlim(1, 500) +theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_text(data = data.frame(x = 430, y = 100, label = paste('adj. Rand = ', round(rand5, 3))), aes(x = x, y = y, label = label))

plot_grid(p1, p2, p3, p4, p5, labels = c(), nrow = 5, align = 'vh')


data[, n_cpts_found := unlist(lapply(found_cpts, length))]
data[deletion_method == 'blockwise' & delete_fraction == 0.3, sum(n_cpts_found)] # 1343
sum(data[deletion_method == 'blockwise' & delete_fraction == 0.3, mapply(function(i,j){sum(sapply(i, function(k) min(abs(k - j))) >= 6)}, found_cpts, true_cpts)]) # 264
sum(data[deletion_method == 'blockwise' & delete_fraction == 0.3, mapply(function(i,j){sum(sapply(i, function(k) min(abs(k - j))) >= 11)}, found_cpts, true_cpts)]) # 153

data[deletion_method == 'blockwise' & delete_fraction == 0.4, sum(n_cpts_found)] # 1118
sum(data[deletion_method == 'blockwise' & delete_fraction == 0.4, mapply(function(i,j){sum(sapply(i, function(k) min(abs(k - j))) >= 6)}, found_cpts, true_cpts)]) # 364
sum(data[deletion_method == 'blockwise' & delete_fraction == 0.4, mapply(function(i,j){sum(sapply(i, function(k) min(abs(k - j))) >= 11)}, found_cpts, true_cpts)]) # 236

data[deletion_method == 'blockwise' & delete_fraction == 0.5, sum(n_cpts_found)] # 724
sum(data[deletion_method == 'blockwise' & delete_fraction == 0.5, mapply(function(i,j){sum(sapply(i, function(k) min(abs(k - j))) >= 6)}, found_cpts, true_cpts)]) # 371
sum(data[deletion_method == 'blockwise' & delete_fraction == 0.5, mapply(function(i,j){sum(sapply(i, function(k) min(abs(k - j))) >= 11)}, found_cpts, true_cpts)]) # 277
