library(tidyverse)
library(ggplot2)

file <- list.files('1.read-length-data/',pattern = '.txt')

# bacth read
map_df(file,function(x){
  tmp = read.table(paste('1.read-length-data/',x,sep = ''))
  colnames(tmp) <- c('readlength','numbers')
  tmp$percent <- tmp$numbers/sum(tmp$numbers)
  tmp$sample <- sapply(strsplit(x,split = '\\.'),'[',1)
  # filter length 25-35 nt
  tmp <- tmp %>% filter(readlength >= 25 & readlength <= 35)
  return(tmp)
}) -> len

# plot 6x15
ggplot(len,aes(x = factor(readlength),y = percent)) +
  geom_col() +
  theme_bw(base_size = 14) +
  xlab('') +
  facet_wrap(~sample,scales = 'free',ncol = 4)
