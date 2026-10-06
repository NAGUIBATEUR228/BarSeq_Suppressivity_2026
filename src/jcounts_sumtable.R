library(tidyverse)
library("ggrepel")
library('patchwork')

"%+%" <- function(...){
  paste0(...)
}

remove_excess <-function(string, excess){
  Reduce(function(x, y) gsub(y, "", x), excess, init = string)
}

folder_with_fastq<-'2026_article_suppr'
#Set working directory
workingDIR <- "C:/Users/zokmi/space/2026_article_suppr/"
setwd(workingDIR)

#write down a folder with fastq
folder <- 'seq_data/'
#output folder
output_folder <- 'output/'
#summary
sum_folder <- 'data/'
#reference tables
ref_folder <- 'ref/'
excess <- c('_R', '_L00')

seqfiles<-unlist(map(output_folder, ~paste(.x,dir(.x)%>%str_subset("output_count.csv$"),sep='')))

tablenames<-seqfiles%>%
  str_match_all("\\s*(.*?)\\s*_output_count.csv")%>%
  map_chr(function(x){x[,2]})

#table with other synonyms of gene names. https://www.yeastgenome.org/
yeast_genes <- read_tsv(ref_folder%+%'yeast_gene_names.txt')
#function converts all gene names either to systematic or standard with to_std=T option. 
convert_names <- function(vals, to_std = F){
  if (length(vals) == 0){return(vals)}
  l <- tibble(syn = vals, id = (1:length(vals)))%>%
    #gene names, combined with '|' are converted separately and returned back joined
    separate_longer_delim(syn, '|')%>%
    left_join(yeast_genes, by = 'syn')
  if (to_std) {
    return(l%>%
             mutate(res = ifelse(is.na(std_name),
                                 syn,
                                 std_name))%>%
             group_by(id)%>%
             summarise(res = str_flatten(res%>%unique, collapse = '|'))%>%
             pull(res))
  }
  return(l%>%
           mutate(res = ifelse(is.na(sys_name),
                               syn,
                               sys_name))%>%
           group_by(id)%>%
           summarise(res = str_flatten(res%>%unique, collapse = '|'))%>%
           pull(res))
}

tighten <- function(x){
  x%>%
    mutate(std_name = convert_names(Confirmed_deletion, to_std = T))%>%
    transmute(sys_name = Confirmed_deletion, std_name, count)
}

table_joining <- function(x){
  exp_name<-ifelse(str_ends(x, '_output_count.csv'),
                   str_match_all(x, output_folder%+%"\\s*(.*?)\\s*_output_count.csv$")%>%map_chr(function(x){x[,2]}),
                   ifelse(str_ends(x, '_output_dm_count.csv'),
                          str_match_all(x, output_folder%+%"\\s*(.*?)\\s*_output_dm_count.csv$")%>%map_chr(function(x){x[,2]}),
                          ifelse(str_ends(x, '_blasted_dm.csv'),
                                 str_match_all(x, output_folder%+%"\\s*(.*?)\\s*_blasted_dm.csv$")%>%map_chr(function(x){x[,2]}),
                                 str_match_all(x, output_folder%+%"\\s*(.*?)\\s*_blasted.csv$")%>%map_chr(function(x){x[,2]}))))
  
  read_csv(x, show_col_types = F)%>%
    tighten%>%
    mutate(exp = exp_name)
}

tbs<-c(tablenames %+% '_output_count.csv',
       tablenames %+% '_output_dm_count.csv',
       tablenames %+% '_blasted_dm.csv',
       tablenames %+% '_blasted.csv')
#making jcounts table
counts <- purrr::reduce(map(tbs, table_joining), full_join)
wide_counts<-counts%>%
  mutate(exp = remove_excess(exp, excess)%>%tolower())%>%
  group_by(sys_name, std_name, exp)%>%
  summarise(count=sum(count))%>%
  ungroup()%>%
  pivot_wider(names_from = exp, values_from = count)%>%
  mutate(across(where(is.numeric),
                ~ifelse(is.na(.x), 0, .x)))%>%
  mutate(name = ifelse(is.na(std_name), sys_name, std_name))
colnames(wide_counts)

jcounts <- wide_counts%>%dplyr::select(name, where(is.numeric))

jcounts <- jcounts %>% pivot_longer(- name, 
                                    values_to = "count", 
                                    names_to = "Conditions") %>%
  mutate(Experiments = str_sub(Conditions, 1, -2)) %>%
  group_by(name, Experiments) %>%
  summarise(count = sum(count)) %>% 
  ungroup()%>%
  pivot_wider(names_from = "Experiments", values_from= "count")
write_csv(jcounts,sum_folder%+%'jcounts_'%+%folder_with_fastq%+%'.csv')

jcounts<-read_csv(sum_folder%+%'jcounts_'%+%folder_with_fastq%+%'.csv')
reference<-read_csv(ref_folder%+%'reference.txt')
sumtsh<-read_tsv(sum_folder%+%'sumtable_short.txt')
sb<-read_csv(sum_folder%+%'sumblast.csv')
sbd<-read_csv(sum_folder%+%'sumblastdm.csv')
sd<-read_csv(sum_folder%+%'sumdm.csv')

sumtable<-tibble(exp='tmp',
                 total_reads=0,
                 barcoded=0,
                 barcoded_dm=0,
                 matched=0,
                 matched_dm=0,
                 blasted=0,
                 blasted_dm=0,
                 final_mean=0,
                 final_max=0,
                 
                 final_unique=0,
                 final_unique_ref=0,
                 final_unique_10=0,
                 final_unique_ref_10=0,
                 
                 barcoded_mean_qual=0,
                 dm_mean_qual=0,
                 nmdm_mean_qual=0,
                 p_barcoded=0,
                 p_barcoded_dm=0,
                 matched_in_barcoded=0,
                 matched_in_barcoded_dm=0,
                 blasted_in_barcoded=0,
                 blasted_in_barcoded_dm=0,
                 ORFs_in_Puddu=0,
                 ORFs_in_YKOC=0,
                 final_ORF_collision_count=0,
                 reads_from_Puddu=0,
                 reads_from_YKOC=0
)

make_row <-function(file){
  
  ref<-get('reference',envir=.GlobalEnv)
  r<-get('sumtsh',envir=.GlobalEnv)
  sb<-get('sb',envir=.GlobalEnv)
  sbd<-get('sbd',envir=.GlobalEnv)
  sd<-get('sd',envir=.GlobalEnv)
  jcounts<-get('jcounts',envir=.GlobalEnv)
  
  a<-read_csv(paste(file,'output_count',sep='_')%+%'.csv')
  dm<-read_csv(paste(file,'output_dm_count',sep='_')%+%'.csv')
  b<-read_csv(paste(file,'blasted',sep='_')%+%'.csv')
  bdm<-read_csv(paste(file,'blasted_dm',sep='_')%+%'.csv')
  nm<-read_csv(paste(file,'not_matched',sep='_')%+%'.csv')
  nmdm<-read_csv(paste(file,'not_matched_dm',sep='_')%+%'.csv')
  q<-read_csv(paste(file,'barcoded_qual_count',sep='_')%+%'.csv')
  nb<-read_csv(paste(file,'not_barcoded_raw_qual_count',sep='_')%+%'.csv')
  
  unique_count<-dm%>%
    select(barcode,count)%>%
    mutate(n=str_count(barcode,"\\|")+1)%>%
    separate_longer_delim(barcode,'|')%>%
    mutate(count=count/n)%>%
    select(barcode, count)%>%
    add_row(q%>%select(barcode,count))%>%
    add_row(nmdm%>%select(barcode,count))%>%
    group_by(barcode)%>%
    summarise(count=sum(count))%>%
    ungroup%>%
    mutate(greater10=(count>10),
           inref=(barcode%in%ref$UPTAG_seqs))%>%
    select(greater10, inref)%>%
    table
  #reads from original (not mutated)
  rfo<-map_dfr(list(a,dm,b%>%mutate(barcode=original_barcode),bdm%>%mutate(barcode=original_barcode)),function(x){x%>%
      filter(barcode %in%
               (reference%>%
                  filter(str_detect(UPTAG_notes,'YKOC'))%>%
                  pull(UPTAG_seqs)))%>%
      dplyr::select(count)})%>%
    sum
  #reads from both
  rfb<-map_dfr(list(a,dm,b%>%mutate(barcode=original_barcode),bdm%>%mutate(barcode=original_barcode)),function(x){x%>%
      filter(barcode %in%
               (reference%>%
                  filter(str_detect(UPTAG_notes,'\\|'))%>%
                  pull(UPTAG_seqs)))%>%
      dplyr::select(count)})%>%
    sum
  
  exp_name<-str_match(file, output_folder%+%"\\s*(.*?)\\s*$")[,2]
  
  sumrow<-tibble(exp=exp_name, 
                 total_reads= (r%>%filter(str_detect(name, exp_name%+%'$')))$total_count,
                 barcoded=sum(q$count),
                 barcoded_dm=sum(dm$count)+sum(nmdm$count),
                 matched=sum(a$count),
                 matched_dm=sum(dm$count),
                 blasted=sum(b$count),
                 blasted_dm=sum(bdm$count),
                 final_mean=mean(jcounts%>%pull(exp_name%>%remove_excess(excess)%>%tolower()%>%str_sub(1,-2))),
                 final_max=max(jcounts%>%pull(exp_name%>%remove_excess(excess)%>%tolower()%>%str_sub(1,-2))),
                 
                 final_unique=sum(unique_count[,]),
                 final_unique_ref=sum(unique_count["TRUE",]),
                 final_unique_10=sum(unique_count[,"TRUE"]),
                 final_unique_ref_10=unique_count["TRUE","TRUE"],
                 
                 barcoded_mean_qual=mean(q$qual),
                 dm_mean_qual=mean(nb$qual),
                 nmdm_mean_qual=mean(nmdm$qual),
                 p_barcoded=0,
                 p_barcoded_dm=0,
                 matched_in_barcoded=0,
                 matched_in_barcoded_dm=0,
                 blasted_in_barcoded=0,
                 blasted_in_barcoded_dm=0,
                 ORFs_in_Puddu=map_dfr(
                   list(a,dm,b,bdm),
                   function(x){x%>%
                       filter(str_detect(notes,'\\|')|str_detect(notes,'puddu'))%>%
                       dplyr::select(Confirmed_deletion)})%>%
                   unique%>%nrow,
                 ORFs_in_YKOC=map_dfr(
                   list(a,dm,b,bdm),
                   function(x){x%>%filter(str_detect(notes,'YKOC'))%>%
                       dplyr::select(Confirmed_deletion)})%>%
                   unique%>%nrow,
                 final_ORF_collision_count=sum(jcounts%>%filter(str_detect(name,'\\|'))%>%dplyr::select(sym(exp_name%>%remove_excess(excess)%>%tolower()%>%str_sub(1,-2)))),
                 reads_from_Puddu=rfb,#temporary value, see below
                 reads_from_YKOC=rfo
  )
  return(sumrow)
}

for (i in tablenames){
  sumtable<-sumtable%>%
    add_row(make_row(i))
}

sumtable<-sumtable%>%
  filter(exp!='tmp')%>%
  mutate(
    p_barcoded=barcoded/total_reads,
    p_barcoded_dm=barcoded_dm/total_reads,
    matched_in_barcoded=matched/barcoded,
    matched_in_barcoded_dm=matched_dm/barcoded_dm,
    blasted_in_barcoded=blasted/barcoded,
    blasted_in_barcoded_dm=blasted_dm/barcoded_dm,
    reads_from_Puddu=matched+matched_dm-reads_from_YKOC+reads_from_Puddu
  )%>%
  rowwise%>%
  mutate(
   
    barcoded_total = barcoded + barcoded_dm,
    assigned       = matched + blasted + matched_dm + blasted_dm,
    yield = assigned /total_reads,
    
    not_matched    = barcoded    - matched    - blasted,
    not_matched_dm = barcoded_dm - matched_dm - blasted_dm,
    
    p_assigned_of_barcoded = assigned / barcoded_total,
    
    dm_share_of_barcoded   = barcoded_dm / barcoded_total,
    dm_share_of_assigned   = (matched_dm + blasted_dm) / assigned,
    
    p_ref_of_deep = final_unique_ref_10 / final_unique_10,
    p_deep_of_ref = final_unique_ref_10 / final_unique_ref,
    
    
    )

write_csv(sumtable,sum_folder%+%'sumtable.csv')

