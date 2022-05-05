library(ggplot2)
library(tidyverse)
library(ggsci)
library(data.table)
library(patchwork)

###################################################################### plus strand
pname <- list.files(pattern = 'Penrichment.txt')

map_df(1:length(pname), function(x){
  tmp = fread(pname[x])
  colnames(tmp) <- c('pos','density')
  tmp$sample <- sapply(strsplit(pname[x],split = '\\.'),'[',1)
  tmp$type <- sapply(strsplit(pname[x],split = '\\-'),'[',1)
  return(tmp)
}) %>% data.table() -> df_enrichment      

# load gene
gene <- read.table('target-gene.txt')
colnames(gene) <- c('gene_name','start','end')
gene
#   gene_name  start    end
# 1      PMT1 287059 289512
# 2     CDC37 790328 791848
# 3      CCT3 407558 409162
# 4      SSC1 519638 521602
# 5      PDI1  48653  50221

# filter gene
map_df(1:3,function(x){
  ginfo <- gene[x,]
  # filter gene region
  tmp = df_enrichment %>% filter(pos >= ginfo$start & pos <= ginfo$end)
  tmp$gene <- ginfo$gene_name
  # get replicates mean values
  mean_tmp <- tmp %>% group_by(gene,type,pos) %>% 
    summarise(mean_density = mean(density),
              sd = sd(density)) %>% data.table()
  # add codon position
  sq <- seq(ginfo$start,ginfo$end,3);rg <- c(1:length(sq))
  map_df(1:length(sq),function(x){
    tmp <- mean_tmp[pos >= sq[x] & pos <= sq[x] + 2
    ][,.(mean_density = mean(mean_density),mean_sd = mean(sd)),by = .(gene,type)
    ][,`:=`(codon_pos = rg[x])]
  }) -> res
  return(res)
}) -> my_geneinfo

# check
head(my_geneinfo,3)
#    gene type mean_density    mean_sd codon_pos
# 1: PMT1 Ssb1    0.7860251 0.15900585         1
# 2: PMT1 Ssb2    0.7021540 0.00389072         1
# 3: PMT1 Ssb1    0.7508004 0.13647891         2

# PMT1
PMT1 <- ggplot(my_geneinfo %>% filter(gene == 'PMT1'),
       aes(x = codon_pos,y = mean_density)) +
  geom_line(aes(color = type),size = 1) +
  geom_hline(yintercept = 1.5,lty = 'dashed',color = 'grey50',size = 1) +
  geom_ribbon(aes(ymin = mean_density - mean_sd,
                  ymax = mean_density + mean_sd,
                  fill = type),
              alpha = 0.5) +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  theme(legend.position = c(0.85,0.9),
        legend.background = element_blank(),
        strip.background = element_rect(color = NA,fill = 'grey')) +
  ylab('Enrichment [AU]') +
  xlab('Codons / amino acids') +
  ggtitle('PMT1')

CDC37 <- ggplot(my_geneinfo %>% filter(gene == 'CDC37'),
               aes(x = codon_pos,y = mean_density)) +
  geom_line(aes(color = type),size = 1) +
  geom_hline(yintercept = 1.5,lty = 'dashed',color = 'grey50',size = 1) +
  geom_ribbon(aes(ymin = mean_density - mean_sd,
                  ymax = mean_density + mean_sd,
                  fill = type),
              alpha = 0.5) +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  theme(legend.position = c(0.85,0.9),
        legend.background = element_blank(),
        strip.background = element_rect(color = NA,fill = 'grey')) +
  ylab('Enrichment [AU]') +
  xlab('Codons / amino acids') +
  ggtitle('CDC37')

CCT3 <- ggplot(my_geneinfo %>% filter(gene == 'CCT3'),
                aes(x = codon_pos,y = mean_density)) +
  geom_line(aes(color = type),size = 1) +
  geom_hline(yintercept = 1.5,lty = 'dashed',color = 'grey50',size = 1) +
  geom_ribbon(aes(ymin = mean_density - mean_sd,
                  ymax = mean_density + mean_sd,
                  fill = type),
              alpha = 0.5) +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  theme(legend.position = c(0.85,0.9),
        legend.background = element_blank(),
        strip.background = element_rect(color = NA,fill = 'grey')) +
  ylab('Enrichment [AU]') +
  xlab('Codons / amino acids') +
  ggtitle('CCT3')

# 3x12
PMT1 + CDC37 + CCT3


###################################################################### minus strand
mname <- list.files(pattern = 'Menrichment.txt')

map_df(1:length(pname), function(x){
  tmp = fread(mname[x])
  colnames(tmp) <- c('pos','density')
  tmp$sample <- sapply(strsplit(mname[x],split = '\\.'),'[',1)
  tmp$type <- sapply(strsplit(mname[x],split = '\\-'),'[',1)
  return(tmp)
}) %>% data.table() -> df_enrichment      

# load gene
gene <- read.table('target-gene.txt')
colnames(gene) <- c('gene_name','start','end')

# filter gene
map_df(4:5,function(x){
  ginfo <- gene[x,]
  tmp = df_enrichment %>% filter(pos >= ginfo$start & pos <= ginfo$end)
  tmp$gene <- ginfo$gene_name
  # get replicates mean values
  mean_tmp <- tmp %>% group_by(gene,type,pos) %>% 
    summarise(mean_density = mean(density),
              sd = sd(density)) %>% data.table()
  # add codon position
  sq <- seq(ginfo$start,ginfo$end,3);rg <- c(length(sq):1)
  # sq <- seq(ginfo$start,ginfo$end,3);rg <- c(1:length(sq))
  map_df(1:length(sq),function(x){
    tmp <- mean_tmp[pos >= sq[x] & pos <= sq[x] + 2
    ][,.(mean_density = mean(mean_density),mean_sd = mean(sd)),by = .(gene,type)
    ][,`:=`(codon_pos = rg[x])]
  }) -> res
  return(res)
}) -> my_geneinfo

# SSC1
SSC1 <- ggplot(my_geneinfo %>% filter(gene == 'SSC1'),
               aes(x = codon_pos,y = mean_density)) +
  geom_line(aes(color = type),size = 1) +
  geom_hline(yintercept = 1.5,lty = 'dashed',color = 'grey50',size = 1) +
  geom_ribbon(aes(ymin = mean_density - mean_sd,
                  ymax = mean_density + mean_sd,
                  fill = type),
              alpha = 0.5) +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  theme(legend.position = c(0.85,0.9),
        legend.background = element_blank(),
        strip.background = element_rect(color = NA,fill = 'grey')) +
  ylab('Enrichment [AU]') +
  xlab('Codons / amino acids') +
  ggtitle('SSC1')

# PDI1
PDI1 <- ggplot(my_geneinfo %>% filter(gene == 'PDI1'),
               aes(x = codon_pos,y = mean_density)) +
  geom_line(aes(color = type),size = 1) +
  geom_hline(yintercept = 1.5,lty = 'dashed',color = 'grey50',size = 1) +
  geom_ribbon(aes(ymin = mean_density - mean_sd,
                  ymax = mean_density + mean_sd,
                  fill = type),
              alpha = 0.5) +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  theme(legend.position = c(0.85,0.9),
        legend.background = element_blank(),
        strip.background = element_rect(color = NA,fill = 'grey')) +
  ylab('Enrichment [AU]') +
  xlab('Codons / amino acids') +
  ggtitle('PDI1')

SSC1 + PDI1
