library(ggplot2)
library(tidyverse)
library(ggsci)
library(Rmisc)
library(data.table)

###################################################################### start codon
name <- list.files(path = '4.metagene-data/',pattern = '.txt')

map_df(1:length(name), function(x){
  tmp = read.table(paste('4.metagene-data/',name[x],sep = ''))
  colnames(tmp) <- c('pos','density')
  tmp$sample <- sapply(strsplit(name[x],split = '\\.'),'[',1)
  tmp$type <- sapply(strsplit(name[x],split = '\\-'),'[',1)
  tmp$exp <- sapply(strsplit(name[x],split = '\\-'),'[',2)
  return(tmp)
}) %>% data.table() -> df_st

# single group
ggplot(df_st,aes(x = pos,y = density)) +
  # geom_line(aes(color = exp)) +
  stat_summary(geom = 'line',fun = 'mean',aes(color = exp),size = 1) +
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
  # geom_line(aes(color = exp)) +
  stat_summary(geom = 'line',fun = 'mean',aes(color = exp),size = 1) +
  geom_hline(yintercept = 1,lty = 'dashed',size = 1,color = 'grey50') +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  theme(legend.position = 'top',
        strip.background = element_rect(color = NA,fill = 'grey')) +
  facet_wrap(~type) +
  ylab('Mean read density [AU]') +
  xlab('Codons / amino acids')

# calculate replicates 95% CI
mean_sd_ci <- codon_meta_st %>% 
  dplyr::group_by(type,exp,codon_pos) %>% 
  dplyr::summarise(mean_density = mean(mean_norm_density),
            sd = sd(mean_norm_density),
            upper = CI(mean_norm_density,ci = 0.95)[1],
            lower = CI(mean_norm_density,ci = 0.95)[3])

# mean ci 0.95 codon plot
ggplot(mean_sd_ci,aes(x = codon_pos,y = mean_density)) +
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = exp),
              alpha = 0.3) +
  geom_line(aes(color = exp),size = 1) +
  geom_hline(yintercept = 1,lty = 'dashed',size = 1,color = 'grey50') +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  theme(legend.position = 'top',
        strip.background = element_rect(color = NA,fill = 'grey')) +
  facet_wrap(~type) +
  ylab('Mean read density [AU]') +
  xlab('Codons / amino acids')

# mean sd
ggplot(mean_sd_ci,aes(x = codon_pos,y = mean_density)) +
  geom_ribbon(aes(ymin = mean_density - sd,
                  ymax = mean_density + sd,
                  fill = exp),
              alpha = 0.5) +
  geom_line(aes(color = exp),size = 1) +
  geom_hline(yintercept = 1,lty = 'dashed',size = 1,color = 'grey50') +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  theme(legend.position = 'top',
        strip.background = element_rect(color = NA,fill = 'grey')) +
  facet_wrap(~type) +
  ylab('Mean read density [AU]') +
  xlab('Codons / amino acids')


#############################################################
# enrichment ratio

InputFile <- list.files(path = '4.metagene-data/',pattern = 'trans')
IpFile <- list.files(path = '4.metagene-data/',pattern = 'inter')

map_df(1:length(InputFile),function(x){
  # 1.load ip file
  input = read.table(paste('4.metagene-data/',InputFile[x],sep = ''))
  colnames(input) <- c('pos','density')
  # 2.load input file
  ip = read.table(paste('4.metagene-data/',IpFile[x],sep = ''))
  colnames(ip) <- c('pos','density')
  # 3.merge
  mer_df <- merge(input,ip,by = 'pos')
  # 4.calculate ratio IP/Input
  mer_df$ratio <- ifelse(mer_df$density.x == 0,0,mer_df$density.y/mer_df$density.x)
  # add sample info
  id = sapply(strsplit(InputFile[x],split = '\\.'),'[',1)
  mer_df$type <- sapply(strsplit(id,split = '\\-'),'[',1)
  mer_df$group <- paste(sapply(strsplit(id,split = '\\-'),'[',1),
                        sapply(strsplit(id,split = '\\-'),'[',3),
                        sep = '-')
  mer_df <- mer_df %>% select(type,group,pos,ratio)
  return(mer_df)
}) %>% data.table() -> enrich_df

########################## codon position transform

sq <- seq(-51,1500,3);rg <- c(-17:500)

map_df(unique(enrich_df$group),function(x){
  tmp = enrich_df[group %in% x] %>% arrange(pos)
  map_df(1:length(sq),function(z){
    tmp1 = tmp[pos >= sq[z] & pos <= sq[z] + 2
    ][,.(mean_ratio = mean(ratio)),by = .(type,group)
    ][,`:=`(codon_pos = rg[z])]
    return(tmp1)
  }) -> res
}) -> codon_enrich_df


# calculate replicates 95% CI
mean_sd_ci_enrichdf <- codon_enrich_df %>% 
  dplyr::group_by(type,codon_pos) %>% 
  dplyr::summarise(rep_mean_ratio = mean(mean_ratio),
                   sd = sd(mean_ratio),
                   upper = CI(mean_ratio,ci = 0.95)[1],
                   lower = CI(mean_ratio,ci = 0.95)[3])


# mean ci 0.95 codon plot
ggplot(mean_sd_ci_enrichdf,aes(x = codon_pos,y = rep_mean_ratio)) +
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = type),
              alpha = 0.3) +
  geom_line(aes(color = type),size = 1) +
  geom_hline(yintercept = 1,lty = 'dashed',size = 1,color = 'grey50') +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  theme(legend.position = 'top',
        strip.background = element_rect(color = NA,fill = 'grey')) +
  facet_wrap(~type) +
  ylab('Mean enrichment [AU]') +
  xlab('Codons / amino acids')

# mean sd
ggplot(mean_sd_ci_enrichdf,aes(x = codon_pos,y = rep_mean_ratio)) +
  geom_ribbon(aes(ymin = rep_mean_ratio - sd,
                  ymax = rep_mean_ratio + sd,
                  fill = type),
              alpha = 0.5) +
  geom_line(aes(color = type),size = 1) +
  geom_hline(yintercept = 1,lty = 'dashed',size = 1,color = 'grey50') +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  theme(legend.position = 'top',
        strip.background = element_rect(color = NA,fill = 'grey')) +
  facet_wrap(~type) +
  ylab('Mean enrichment [AU]') +
  xlab('Codons / amino acids')
