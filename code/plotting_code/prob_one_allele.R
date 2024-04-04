### This code is used to plot the probability of observing at least one allele
# based on sample size and AF using a binomial distribution, mainly to determine
# what default reference sample size to use


library(dplyr)
library(tidyr)
library(ggplot2)

dir_out =  paste0('C:/Users/sagee/Documents/GitHub/masters_project/Results/manuscript_plots/')

N = seq(500, 10000, by=500)
AF = c(1/100, 1/1000, 1/2000, 1/5000, 1/10000)

p1 = pbinom(1, 2*N, 1/100, lower.tail = FALSE)
p2 = pbinom(1, 2*N, 1/1000, lower.tail = FALSE)
p3 = pbinom(1, 2*N, 1/2000, lower.tail = FALSE)
p4 = pbinom(1, 2*N, 1/5000, lower.tail = FALSE)
p5 = pbinom(1, 2*N, 1/10000, lower.tail = FALSE)

df <- data.frame(AF = rep(AF, each=length(N)), N=rep(N, times=length(AF)), prob=c(p1, p2, p3, p4, p5))
df$AF = factor(df$AF, levels = AF)

cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 = c("#999999", "#BC9F4C", "#56B4E9", "#009E73", "#0072B2")
af_cols = c("#CC79A7", "#D55E00", "#BC9F4C", "#009E73", "#56B4E9")

p1 <- ggplot(df, aes(N, prob, color=AF)) +
  geom_line(linewidth = 1) +
  geom_point(aes(color=AF)) +
  geom_vline(xintercept = 900, linetype = 2) +
  scale_color_manual(values = af_cols) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_x_continuous(breaks = c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)) +
  labs(y='Probability', x='Sample Size', title = 'Probability of observing at least one allele by allele frequency (AF)') 
p1

ggsave(file = paste0(dir_out, 'prob_by_AF_and_sample_size_more.jpg'), plot = p1, height = 8, width = 10, units = 'in')



