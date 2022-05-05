library(ggplot2)
library(tidyverse)
library(ggsci)
library(Rmisc)
library(data.table)

name <- list.files(path = '5.enrichment-data/',pattern = '.txt')

map_df(1:length(name), function(x){
  tmp <- fread(paste('5.enrichment-data/',name[x],sep = ''),sep = '\t')
  colnames(tmp) <- c('id','pos','ip_denisty','input_density','ratio')
  tmp <- tmp %>%  select(id,pos,ratio)
  # add gene name and filter target gene
  tmp <- tmp[, c("gene_name") := tstrsplit(id, "|", fixed=TRUE)[1]] %>% 
    filter(gene_name %in% c('PMT1','CDC37','CCT3','SSC1','PDI1'))
  
  # loop for every gene to process data
  map_df(unique(tmp$gene_name),function(y){
    tmp1 <- tmp[gene_name == y]
    start <- sapply(strsplit(tmp1$id[1],split = '\\|'),'[',3) %>% as.numeric()
    end <- sapply(strsplit(tmp1$id[1],split = '\\|'),'[',4) %>% as.numeric()
    
    # 1.add nt pos
    tmp1$nt_pos <- tmp1$pos - start + 1
    
    # 2.transform to codon pos
    sq <- seq(1,(end - start + 1),3);rg <- c(1:length(sq))
    map_df(1:length(sq),function(z){
      tmp2 = tmp1[nt_pos >= sq[z] & nt_pos <= sq[z] + 2
      ][,.(mean_ratio = mean(ratio)),by = .(id,gene_name)
      ][,`:=`(codon_pos = rg[z])]
      return(tmp2)
    }) -> codon_res
    
    # 3.add to continues codon positions
    map_df(rg,function(s){
      if(s %in% codon_res$codon_pos){
        tmp3 <- as.vector(codon_res[codon_pos == s])
      }else{
        tmp3 <- data.table(id = codon_res$id[1],gene_name = codon_res$gene_name[1],
                           mean_ratio = 0,codon_pos = s)
      }
      return(tmp3)
    }) -> continues_codon_res
    return(continues_codon_res)
  }) -> final_res
  
  # add sample info
  final_res$sample <- sapply(strsplit(name[x],split = '\\.'),'[',1)
  final_res$type <- sapply(strsplit(name[x],split = '\\-'),'[',1)
  final_res$exp <- sapply(strsplit(name[x],split = '\\-'),'[',2)
  return(final_res)
}) %>% data.table() -> df_ratio

##################################################

# mean for replicates
merge_rep <- df_ratio %>% 
  dplyr::group_by(type,gene_name,codon_pos) %>% 
  dplyr::summarise(mean_rep_ratio = mean(mean_ratio),
            mean_sd = sd(mean_ratio))

# check
head(merge_rep,3)
# # A tibble: 3 x 5
# # Groups:   type, gene_name [1]
#   type  gene_name codon_pos mean_rep_ratio mean_sd
#  <chr> <chr>          <int>          <dbl>   <dbl>
# 1 ssb1  CCT3              1          0.906  0.272 
# 2 ssb1  CCT3              2          0.830  0.0120
# 3 ssb1  CCT3              3          0.553  0.149 

###################################################
merge_rep$gene_name <- factor(merge_rep$gene_name,
                              levels = c('PMT1','CDC37','CCT3','SSC1','PDI1'))

# plot
ggplot(merge_rep,aes(x = codon_pos,y = mean_rep_ratio)) +
  geom_line(aes(color = type)) +
  geom_hline(yintercept = 1.5,lty = 'dashed',color = 'red',size = 1) +
  geom_ribbon(aes(ymin = mean_rep_ratio - mean_sd,
                  ymax = mean_rep_ratio + mean_sd,
                  fill = type),
              alpha = 0.4) +
  theme_classic(base_size = 16) +
  scale_color_d3(name = '') +
  scale_fill_d3(name = '') +
  theme(legend.background = element_blank(),
        strip.background = element_rect(color = NA,fill = 'grey')) +
  ylab('Mean enrichment [AU] \n (co-IP/total)') +
  xlab('Ribosome position \n (Codons/amino acids)') +
  facet_wrap(~gene_name,scales = 'free',ncol = 3)
  