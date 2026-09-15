library(tidyverse)
library("ggrepel")
library('patchwork')
library(DESeq2)

#Set working directory
workingDIR <- "C:/Users/zokmi/space/2026_article_suppr/"#/
setwd(workingDIR)

#write down a folder with fastq
folder <- 'seq_data/'
#output folder
output_folder <- 'output/'
#summary
sum_folder <- 'data/'
#reference tables
ref_folder <- 'ref/'

"%+%" <- function(...){
  paste0(...)
}

jcounts<-read_tsv(sum_folder%+%'jcounts_final.txt')

jc <- jc <- jcounts
n = jc%>%dplyr::select(!name)%>%colSums%>%mean
jc<-jc%>%mutate(across(!name, ~ .x / sum(.x) * n))

noMtDNA <- read_csv2(ref_folder%+%'noMtDNA.csv')

#table with other synonyms of gene names. https://www.yeastgenome.org/
yeast_genes <- read_tsv(ref_folder%+%'yeast_gene_names.txt')
annot<-read_tsv(ref_folder%+%"gene_annotation.txt")
go<-read_tsv(ref_folder%+%"go_annotation.txt")
all_go_descriptions<-read_tsv(ref_folder%+%"go_description.txt")

#####
# FUNCTIONS
#####

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

t_test <- function(mean.x, mean.y, sigma.x, sigma.y, n.x, n.y) {
  # Стандартная ошибка разности средних (для неравных дисперсий)
  se <- sqrt(sigma.x^2 / n.x + sigma.y^2 / n.y)
  
  # t-статистика
  t <- (mean.x - mean.y) / se
  
  # Степени свободы (Уэлча-Саттертуэйта)
  df <- (sigma.x^2 / n.x + sigma.y^2 / n.y)^2 / 
    ((sigma.x^2 / n.x)^2 / (n.x - 1) + (sigma.y^2 / n.y)^2 / (n.y - 1))
  
  # Двустороннее p-value
  p <- 2 * pt(abs(t), df, lower.tail = FALSE)
  
  return(p)
}

reads_z_test <- function(data, t, c){
  data%>%
    dplyr::select(name, starts_with(t), starts_with(c))%>%
    pivot_longer(cols = !name, names_to = 'exp')%>%
    mutate(rep = str_sub(exp,-1,-1)%>%as.numeric,
           exp = str_sub(exp,1,-2))%>%
    mutate(exp = exp%>%factor(levels=c(t,c)))%>%
    arrange(name,exp,rep)%>%
    # filter(-Inf != value)%>%
    group_by(name, exp)%>%
    mutate(m = mean(value, na.rm=T), s=sd(value,na.rm=T))%>%
    na.omit%>%
    ungroup -> to_test
  
  to_test%>%
    group_by(name,exp)%>%
    summarise(m=max(m), n=n(), s=max(s),.groups = 'drop')%>%
    ungroup%>%
    pivot_wider(
      names_from = exp,
      values_from = c(m, n, s),  # разделяем средние и размеры
      names_sep = "_"
    )  %>%
    mutate(delta=!!sym("m_"%+%t)-!!sym("m_"%+%c),
           f = (!!sym("m_"%+%t)-!!sym("m_"%+%c))/abs(!!sym("m_"%+%c)))%>%
    na.omit%>%
    mutate(f_WT=f-median(f),
           delta_WT=delta-median(delta),
           delta_WT_norm=(delta-median(delta))/abs(median(delta)))%>%
    mutate(
      pvalue_WT = t_test(
        mean.x = delta_WT,
        mean.y = 0,
        sigma.x = !!sym("s_"%+%t),
        sigma.y = !!sym("s_"%+%c),
        n.x = !!sym("n_"%+%t),
        n.y = !!sym("n_"%+%c)
      )
    )%>%
    mutate(
      pvalue = t_test(
        mean.x = delta,
        mean.y = 0,
        sigma.x = !!sym("s_"%+%t),
        sigma.y = !!sym("s_"%+%c),
        n.x = !!sym("n_"%+%t),
        n.y = !!sym("n_"%+%c)
      )
    )%>%
    na.omit%>%
    mutate(control_val = !!sym(paste0("m_", t)))%>%
    transmute(name,LogFC=!!sym("m_"%+%t) - !!sym("m_"%+%c),"mean_{t}" :=control_val, mean_control=!!sym("m_"%+%c),f,f_WT, delta, delta_WT, delta_WT_norm, pvalue,pvalue_WT)%>%
    mutate(padj=p.adjust(pvalue, method='BH'),
           padj_WT=p.adjust(pvalue_WT, method='BH'))
}

difexp<-function(jc,t,u){
  
  cd <- jc %>% select(!name)%>%select(starts_with(t)|starts_with(u))
  
  colnames <- data.frame(experiment = colnames(cd)[1:ncol(cd)])
  rownames(colnames) <- colnames[,1]
  colnames <- colnames %>% mutate(condition = str_sub(experiment, 1,-2)) %>%
    mutate(rep = str_sub(experiment, -1, -1))  
  
  ### DDS
  
  dds <- DESeqDataSetFromMatrix(countData = cd,
                                colData = colnames,
                                design = ~ condition)
  dds
  keep <- rowSums(counts(dds)) >= 100
  dds <- dds[keep,]
  
  dds <- DESeq(dds)
  
  dds$condition
  res <- results(dds, contrast=c("condition",t,u))
  return(res)
}
# res <- results(dds, contrast=c("condition","treated","untreated"))

#######

noMtDNA1 <- convert_names(na.omit(noMtDNA$noMtDNA_Puddu))
noMtDNA2 <- convert_names(na.omit(noMtDNA$Pet_Stenger))
all_noMtDNA <- unique(c(noMtDNA1, noMtDNA2))

no_mating<-go%>%filter(go_term=='GO:0000747')%>%pull(gene_id)

names(jcounts)[-1]%>%str_sub(1,-2)%>%unique
conds<-c('pd','pg','rd','rg','yd','yg')
# cond_pairs<-str_split_1("pg_pd,yg_yd,rg_rd,pd_rd,pg_rg,pd_yd,pg_yg,rd_yd,rg_yg",",")
cond_pairs<-str_split_1("pg_pd,pg_rg",",")

logd <- 'filter/'
dir.create(workingDIR %+% sum_folder %+% logd, showWarnings = FALSE)
logd<- sum_folder%+% logd


#####
# first approach: jcounts -> DESeq2
#####
apr_name <- 'f1'

### res <- results(dds, contrast=c("condition","treated","untreated"))

jc<-column_to_rownames(jcounts,'name')%>%
  mutate(name=jcounts$name)
for (pair in cond_pairs){
  t<-str_split_1(pair,"_")[1]
  u <-str_split_1(pair,"_")[2]
  filename <- paste(pair, apr_name, sep='_')
  res<-jc%>%
    difexp(t,u)
  res%>%
    as.data.frame%>%
    rownames_to_column('name')%>%
    as_tibble%>%
    transmute(name,LogFC=log2FoldChange,pvalue,padj)%>%
    assign(x=filename,value=.,envir = .GlobalEnv)
  get(filename)%>%
    write_csv(logd%+% filename %+% '.csv')
}
cond_list<-paste(cond_pairs, apr_name, sep='_')
assign(apr_name, purrr::reduce(map(cond_list,
                                   function(x){
                                     get(x)%>%
                                       mutate(exp=str_sub(x,1,5))}),
                               ~add_row(.x,.y)))
write_csv(get(apr_name), logd%+%apr_name%+%'.csv')


#####
# log reads
#####

jc <- jcounts %>%
  filter(if_any(where(is.numeric), ~.x!=0))
jc%>%
  dplyr::select(name, starts_with('yd'))%>%
  pivot_longer(cols = !name, names_to = 'exp')%>%
  group_by(name)%>%
  mutate(s = sum(value))%>%
  filter(s >= 100)%>%
  transmute(name)%>%
  ungroup()%>%
  distinct()%>%
  left_join(jc)->jc
n = jc %>% dplyr::select(!name) %>% colSums %>% median
jc <- jc %>% mutate(across(!name, ~ (.x+1) / sum(.x+1) * n))


jc%>%
  mutate(across(!name,~ifelse(.x==0,NA,log2(.x))))->data
write_csv(data,logd%+%'log_reads.csv')
sum(unique(convert_names(c(noMtDNA2,noMtDNA1), to_std=T))%in%data$name)
data<-read_csv(logd%+%'log_reads.csv')


#####
# second approach
#####
apr_name <- 'f2'

#'data' is already recorded

for (pair in cond_pairs){
  t<-str_split_1(pair,"_")[1]
  c <-str_split_1(pair,"_")[2]
  filename <- paste(pair, apr_name, sep='_')
  assign(x=filename,value=reads_z_test(data, t, c),envir = .GlobalEnv)
  get(filename)%>%
    write_csv(logd%+% filename %+% '.csv')
}
cond_list<-paste(cond_pairs, apr_name, sep='_')
assign(apr_name, purrr::reduce(map(cond_list,
                                   function(x){
                                     get(x)%>%
                                       mutate(exp=str_sub(x,1,5))}),
                               ~add_row(.x,.y)))


write_csv(get(apr_name), logd%+%apr_name%+%'.csv')


######
## normalise
######


#cond1 - cond2 vs cond3 - cond4
# cond<-c("pg", "rg", "rd", "yd")
cond<-c("pg", "rg", "pd", "rd")

cond_numbers <- setNames(1:4, cond)

res_cont<-data %>% #mutate(across(!name,~log(exp(.x)+1)))%>%
  pivot_longer(cols = -name, names_to = "variable", values_to = "value") %>%
  mutate(
    condition = str_sub(variable, 1, 2),   # первые два символа: pd, pg, rd, rg
    rep = str_sub(variable, 3)              # остальное — номер реплики
  ) %>%
  filter(condition %in% cond) %>%   # оставляем только нужные условия
  # filter(is.finite(value)) %>%                            # убираем -Inf, NA и т.п.
  left_join(jcounts%>%pivot_longer(!name, names_to = "variable", values_to = "reads"))%>%
  group_by(name)%>%
  filter(sum(reads)>50)%>%
  dplyr::select(name, condition, rep, value)%>%
  group_by(name, condition) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val = sd(value, na.rm = TRUE),
    n_val = n(),
    .groups = "drop"
  )%>%
  mutate(condition = cond_numbers[condition]) %>%
  pivot_wider(
    id_cols = name,
    names_from = condition,
    values_from = c(mean_val, sd_val, n_val),
    names_sep = "_"
  ) %>%
  mutate(
    # разности для супрессивного и нормального фона
    diff_t = mean_val_1 - mean_val_2,
    diff_u = mean_val_3 - mean_val_4,
    # интересующая разность разностей
    delta = diff_t - diff_u,
    
    # дисперсии разностей (квадраты стандартных ошибок)
    var_t = sd_val_1^2 / n_val_1 + sd_val_2^2 / n_val_2,
    var_u = sd_val_3^2 / n_val_3 + sd_val_4^2 / n_val_4,
    # стандартная ошибка для delta
    SE_delta = sqrt(var_t + var_u),
    
    # t-статистика
    t_stat = delta / SE_delta,
    
    # степени свободы по формуле Саттертуэйта для линейной комбинации четырёх средних
    # var_contrib_i = (коэффициент)^2 * (sd_i^2 / n_i)
    # коэффициенты: pg: +1, pd: -1, rg: -1, rd: +1
    var_contrib_1 = (sd_val_1^2 / n_val_1),
    var_contrib_2 = (sd_val_2^2 / n_val_2),
    var_contrib_3 = (sd_val_3^2 / n_val_3),
    var_contrib_4 = (sd_val_4^2 / n_val_4),
    
    df = (SE_delta^4) / (
      (var_contrib_1^2 / (n_val_1 - 1)) +
        (var_contrib_2^2 / (n_val_2 - 1)) +
        (var_contrib_3^2 / (n_val_3 - 1)) +
        (var_contrib_4^2 / (n_val_4 - 1))
    ),
    
    # двустороннее p-value
    pvalue = 2 * pt(abs(t_stat), df, lower.tail = FALSE)
  ) %>%
  dplyr::select(name, delta, SE_delta, t_stat, df, pvalue) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(pvalue)

res_cont%>%
  write_csv(sum_folder%+%'norm_pd.txt')
