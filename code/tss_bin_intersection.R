library(tidyverse)
library(vroom)
library(GenomicRanges)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-------------------------------
hic_dat_in<-function(dat_file){
  chr_dat<-vroom(dat_file,delim = "\t",col_names = F,trim_ws = T,escape_double = F)
  return(chr_dat%>%
           # ensure HiC score variable is a numeric
           mutate(X3=as.numeric(X3))%>%
           # Filter self-interaction
           filter(!(is.na(X3)))%>%filter(X1!=X2) %>% 
           dplyr::rename(bin_a = X1,
                  bin_b=X2,
                  score=X3))
  
}
data_tbl_load_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}

#-------------------------------
CAGE_enh_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_enh_GM12878_Grange.Rda"
CAGE_tss_GRange_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/GRanges/CAGE_tss_GM12878_Grange.Rda"
HiC_data_folder <- "~/Documents/multires_bhicect/data/GM12878/5kb/"
hub_file<-"~/Documents/multires_bhicect/BootTHiC_snakemake_v2/data/results/final_cluster_table/GM12878_final_tbl.tsv"
spec_res_folder<-"~/Documents/multires_bhicect/data/GM12878/determinate_spec_res/"
chromo <- "chr22"

data_file<-paste0(HiC_data_folder,chromo,".txt")
spec_res_file <- paste0(spec_res_folder,chromo,"_spec_res.Rda")

hub_tbl<-read_tsv(hub_file)
chr_hub_tbl <- hub_tbl %>% 
  filter(chr == chromo) %>% 
  mutate(bins=chr_spec_res$cl_member[node]) %>% 
  mutate(bins=map(bins,as.numeric)) 
hub_long_tbl <- chr_hub_tbl %>% 
  unnest(cols = c(bins)) %>% 
  mutate(end = bins + res_num[res] - 1)

chr_hub_Grange<-   GenomicRanges::reduce(GRanges(seqnames=hub_long_tbl$chr,
                                                  ranges = IRanges(start=as.numeric(hub_long_tbl$bins),
                                                                   end=hub_long_tbl$end
                                                  )))


hic_dat_tbl<-hic_dat_in(data_file)

CAGE_tss_GRange <- data_tbl_load_fn(CAGE_tss_GRange_file)
CAGE_enh_GRange <- data_tbl_load_fn(CAGE_enh_GRange_file)

chr_spec_res<-data_tbl_load_fn(spec_res_file)

chr_bins <- unique(c(hic_dat_tbl$bin_a,hic_dat_tbl$bin_b))

chr_bin_Grange<-  GRanges(seqnames=chromo,
                          ranges = IRanges(start=as.numeric(chr_bins),
                                           end=chr_bins + 4999
                                                  ))

hub_bins <- chr_bins[unique(queryHits(findOverlaps(chr_bin_Grange,chr_hub_Grange)))]
tss_bins <- chr_bins[unique(queryHits(findOverlaps(chr_bin_Grange,CAGE_tss_GRange)))]
enh_bins <- chr_bins[unique(queryHits(findOverlaps(chr_bin_Grange,CAGE_enh_GRange)))]

sum(tss_bins %in% hub_bins)/length(tss_bins)
sum(enh_bins %in% hub_bins)/length(enh_bins)
sum(countOverlaps(CAGE_tss_GRange,chr_hub_Grange)>0)/length(CAGE_tss_GRange[seqnames(CAGE_tss_GRange) == chromo])
sum(countOverlaps(CAGE_enh_GRange,chr_hub_Grange)>0)/length(CAGE_enh_GRange[seqnames(CAGE_enh_GRange) == chromo])

prom_enh_hic_tbl<-bind_rows(
  hic_dat_tbl %>% 
    filter(bin_a %in% tss_bins & bin_b %in% enh_bins),
  hic_dat_tbl %>% 
    filter(bin_a %in% enh_bins & bin_b %in% tss_bins)  
) %>% 
  mutate(type="enh_prom")

prom_enh_hic_tbl %>% 
  mutate(d=abs(bin_a-bin_b)) %>% 
  summarise(m=max(d))

