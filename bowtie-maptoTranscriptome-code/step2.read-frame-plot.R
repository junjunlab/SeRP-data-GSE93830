library(tidyverse)
library(data.table)
library(ggsci)
library(ggplot2)
library(scales)

file <- list.files('2.length-frame-data/',pattern = '.txt')

# bacth read
map_df(file,function(x){
  tmp = fread(paste('2.length-frame-data/',x,sep = ''))
  colnames(tmp) <- c('readlength','frame','rel2st','rel2sp','numbers')
  tmp$sample <- sapply(strsplit(x,split = '\\.'),'[',1)
  # filter length 25-35 nt
  tmp <- tmp %>% filter(readlength >= 25 & readlength <= 35)
  return(tmp)
}) -> frame_len

###########################################
# summarise all frame type
all_frame <- frame_len %>% group_by(sample,frame) %>% 
  summarise(count = sum(numbers)) %>% 
  mutate(percent = count/sum(count))

# show colors
show_col(pal_lancet("lanonc")(9))

# plot
ggplot(all_frame,aes(x = factor(frame),y = percent)) +
  geom_col(aes(fill = factor(frame)),width = 0.6) +
  theme_classic() +
  scale_fill_manual(values = c('#ED0000FF','#0099B4FF','#42B540FF'),
                    name = '') +
  theme(strip.background = element_rect(colour = NA,fill = 'grey')) +
  facet_wrap(~sample,scales = 'free',ncol = 4)

###########################################
# summarise all frame type
all_frame_bylen <- frame_len %>% group_by(sample,readlength,frame) %>% 
  summarise(count = sum(numbers))

# plot
ggplot(all_frame_bylen,aes(x = readlength,y = count)) +
  geom_col(aes(fill = factor(frame)),
           position = position_dodge2()) +
  theme_classic() +
  scale_fill_manual(values = c('#ED0000FF','#0099B4FF','#42B540FF'),
                    name = '') +
  theme(strip.background = element_rect(colour = NA,fill = 'grey')) +
  facet_wrap(~sample,scales = 'free',ncol = 4)

###########################################
tost <- frame_len %>% group_by(sample,readlength,rel2st,frame) %>% 
  summarise(count = sum(numbers)) %>%
  filter(rel2st >= -40 & rel2st <= 20)

# plot
ggplot(tost %>% filter(readlength == 29),
       aes(x = rel2st,y = count)) +
  geom_col(aes(fill = factor(frame)),
           position = position_dodge2()) +
  theme_classic() +
  scale_fill_manual(values = c('#ED0000FF','#0099B4FF','#42B540FF'),
                    name = '') +
  xlab('') +
  theme(strip.background = element_rect(colour = NA,fill = 'grey')) +
  facet_wrap(~sample,scales = 'free',ncol = 4)
  # facet_grid(sample~readlength,scales = 'free')

tosp <- frame_len %>% group_by(sample,readlength,rel2sp,frame) %>% 
  summarise(count = sum(numbers)) %>%
  filter(rel2sp >= -40 & rel2sp <= 20)

# plot
ggplot(tosp %>% filter(readlength == 31),
       aes(x = rel2sp,y = count)) +
  geom_col(aes(fill = factor(frame)),
           position = position_dodge2()) +
  theme_classic() +
  scale_fill_manual(values = c('#ED0000FF','#0099B4FF','#42B540FF'),
                    name = '') +
  xlab('') +
  theme(strip.background = element_rect(colour = NA,fill = 'grey')) +
  facet_wrap(~sample,scales = 'free',ncol = 4)

