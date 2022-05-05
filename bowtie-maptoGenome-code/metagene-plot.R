library(ggplot2)
library(tidyverse)
library(ggsci)
library(data.table)

###################################################################### start codon
name <- list.files(pattern = 'MetaDist2StartCodon.txt')

map_df(1:length(name), function(x){
  tmp = read.table(name[x])
  colnames(tmp) <- c('pos','density')
  tmp$sample <- sapply(strsplit(name[x],split = '\\.'),'[',1)
  tmp$type <- sapply(strsplit(name[x],split = '\\-'),'[',1)
  tmp$exp <- sapply(strsplit(name[x],split = '\\-'),'[',2)
  return(tmp)
}) %>% data.table() -> df_st      

# single group
ggplot(df_st,aes(x = pos,y = density)) +
  # geom_line(aes(color = exp)) +
  stat_summary(geom = 'line',fun = 'mean',aes(color = exp)) +
  geom_hline(yintercept = 1,lty = 'dashed',size = 1,color = 'grey50') +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  theme(legend.position = 'top',
        strip.background = element_rect(color = NA,fill = 'grey')) +
  facet_grid(exp~type) +
  ylab('Mean read density [AU]') +
  xlab('Distance from start codon (nt)')

########################## codon position transform

sq <- seq(-51,1500,3);rg <- c(-17:500)

map_df(unique(df_st$sample),function(x){
  tmp = df_st[sample %in% x] %>% arrange(pos)
  map_df(1:length(sq),function(z){
    tmp1 = tmp[pos >= sq[z] & pos <= sq[z] + 2
    ][,.(mean_norm_density = mean(density)),by = .(type,exp,sample)
    ][,`:=`(codon_pos = rg[z])]
    return(tmp1)
  }) -> res
}) -> codon_meta_st

# codon plot
ggplot(codon_meta_st,aes(x = codon_pos,y = mean_norm_density)) +
  geom_line(aes(color = exp)) +
  geom_hline(yintercept = 1,lty = 'dashed',size = 1,color = 'grey50') +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  theme(legend.position = 'top',
        strip.background = element_rect(color = NA,fill = 'grey')) +
  facet_wrap(~type) +
  ylab('Mean read density [AU]') +
  xlab('Codons / amino acids')

# calculate replicates sd
mean_sd <- codon_meta_st %>% group_by(type,exp,codon_pos) %>% 
  summarise(mean_density = mean(mean_norm_density),
            sd = sd(mean_norm_density))

# mean-sd codon plot
ggplot(mean_sd,aes(x = codon_pos,y = mean_density)) +
  geom_line(aes(color = exp)) +
  geom_ribbon(aes(ymin = mean_density - sd,
                  ymax = mean_density + sd,
                  fill = exp),
              alpha = 0.5) +
  geom_hline(yintercept = 1,lty = 'dashed',size = 1,color = 'grey50') +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  theme(legend.position = 'top',
        strip.background = element_rect(color = NA,fill = 'grey')) +
  facet_wrap(~type) +
  ylab('Mean read density [AU]') +
  xlab('Codons / amino acids')


###################################################################### stop codon
name <- list.files(pattern = 'MetaDist2StopCodon.txt')

map_df(1:length(name), function(x){
  tmp = read.table(name[x])
  colnames(tmp) <- c('pos','density')
  tmp$sample <- sapply(strsplit(name[x],split = '\\.'),'[',1)
  tmp$type <- sapply(strsplit(name[x],split = '\\-'),'[',1)
  tmp$exp <- tmp$exp <- sapply(strsplit(name[x],split = '\\-'),'[',2)
  return(tmp)
}) %>% data.table() -> df_sp      

# single group
ggplot(df_sp,aes(x = pos,y = density)) +
  # geom_line(aes(color = exp)) +
  stat_summary(geom = 'line',fun = 'mean',aes(color = exp)) +
  geom_hline(yintercept = 1,lty = 'dashed',size = 1,color = 'grey50') +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  theme(legend.position = 'top',
        strip.background = element_rect(color = NA,fill = 'grey')) +
  facet_grid(exp~type) +
  ylab('Mean read density [AU]') +
  xlab('Distance from start codon (nt)')

########################## codon position transform

sq <- seq(-1501,50,3);rg <- c(-500:17)

map_df(unique(df_sp$sample),function(x){
  tmp = df_sp[sample %in% x] %>% arrange(pos)
  map_df(1:length(sq),function(z){
    tmp1 = tmp[pos >= sq[z] & pos <= sq[z] + 2
    ][,.(mean_norm_density = mean(density)),by = .(type,exp,sample)
    ][,`:=`(codon_pos = rg[z])]
    return(tmp1)
  }) -> res
}) -> codon_meta_sp

# codon plot
ggplot(codon_meta_sp,aes(x = codon_pos,y = mean_norm_density)) +
  geom_line(aes(color = exp)) +
  geom_hline(yintercept = 1,lty = 'dashed',size = 1,color = 'grey50') +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  theme(legend.position = 'top',
        strip.background = element_rect(color = NA,fill = 'grey')) +
  facet_wrap(~type) +
  ylab('Mean read density [AU]') +
  xlab('Codons / amino acids')

# calculate replicates sd
mean_sd <- codon_meta_sp %>% group_by(type,exp,codon_pos) %>% 
  summarise(mean_density = mean(mean_norm_density),
            sd = sd(mean_norm_density))

# mean-sd codon plot
ggplot(mean_sd,aes(x = codon_pos,y = mean_density)) +
  geom_line(aes(color = exp)) +
  geom_hline(yintercept = 1,lty = 'dashed',size = 1,color = 'grey50') +
  geom_ribbon(aes(ymin = mean_density - sd,
                  ymax = mean_density + sd,
                  fill = exp),
              alpha = 0.5) +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  theme(legend.position = 'top',
        strip.background = element_rect(color = NA,fill = 'grey')) +
  facet_wrap(~type) +
  ylab('Mean read density [AU]') +
  xlab('Codons / amino acids')
