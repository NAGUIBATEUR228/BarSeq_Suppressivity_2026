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


#########
## volcanos
#########
volcano_plot <- function(data,
                         x_col,
                         y_col,
                         xlabel='effect',
                         ylabel=expression(-Log[10] ~ (P[adj])),#expression(Log[2] ~ FC)
                         gene_col = "name",
                         center = NULL,
                         FCcutoff = 0.2,
                         pCutoff = 0.01,
                         xlim = NULL,
                         ylim = NULL,
                         pointSize = 4,
                         point_alpha=0.6,
                         labSize = 6,
                         drawConnectors = TRUE,
                         widthConnectors = 0.75,
                         max_overlaps = 20,
                         show_legend = FALSE,
                         title = NULL,
                         show_center_label = FALSE) {
  
  # Подготовка данных
  df <- data %>%
    transmute(
      x = .data[[x_col]],
      y = -log10(.data[[y_col]]),
      gene = .data[[gene_col]]
    ) %>%
    filter(!is.na(x) & !is.na(y) & is.finite(y))
  
  # Центр (по умолчанию медиана)
  if (is.null(center)) {
    center <- median(df$x, na.rm = TRUE)
    message("Центр (медиана x): ", round(center, 4))
  }
  
  # Пороги
  x_left  <- center - FCcutoff
  x_right <- center + FCcutoff
  y_cutoff <- -log10(pCutoff)
  
  # Категории точек (как в EnhancedVolcano)
  df <- df %>%
    mutate(
      category = case_when(
        y < y_cutoff & abs(x - center) < FCcutoff ~ "NS",
        y >= y_cutoff & abs(x - center) < FCcutoff ~ "p-sig",
        y < y_cutoff & abs(x - center) >= FCcutoff ~ "FC-sig",
        y >= y_cutoff & abs(x - center) >= FCcutoff ~ "both-sig",
        TRUE ~ "NS"
      ),
      category = factor(category, levels = c("NS", "p-sig", "FC-sig", "both-sig"))
    )
  
  # Цвета
  colors <- c(
    "NS"       = "grey30",
    "p-sig"    = "green3",
    "FC-sig"   = "royalblue",
    "both-sig" = "red2"
  )
  
  # Подписи только для «both-sig»
  label_df <- df %>% filter(category == "both-sig")
  
  # Базовый график
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_hline(yintercept = y_cutoff, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = c(x_left, x_right), linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_point(aes(color = category), size = pointSize, alpha = point_alpha) +
    scale_color_manual(values = colors) +
    labs(
      title = if (!is.null(title)) title else "Volcano plot",
      x = xlabel,
      y = ylabel
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # Добавляем подписи, если нужно
  if (nrow(label_df) > 0 && drawConnectors) {
    p <- p +
      geom_text_repel(
        data = label_df,
        aes(label = gene),
        size = labSize / 2.5,
        segment.color = "black",
        segment.size = widthConnectors,
        segment.alpha = 0.5,
        arrow = arrow(length = unit(0.01, "npc"), type = "closed"),
        box.padding = 0.25,
        point.padding = 0.3,
        max.overlaps = max_overlaps,
        min.segment.length = 0.2,
        color = "black",
        bg.color = "white",
        bg.r = 0.15
      )
  }
  
  # Подпись центра (только если запрошено)
  if (show_center_label) {
    p <- p +
      annotate("text",
               x = center,
               y = if (nrow(df) > 0) max(df$y, na.rm = TRUE) * 0.95 else 0,
               label = paste("Центр:", round(center, 3)),
               color = "darkred", size = 3, hjust = 0.5, fontface = "italic")
  }
  
  # Устанавливаем пределы осей, если заданы
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)
  }
  if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  p
}

genes <- str_split_1("SDS3,EFG1",",")

norm_pd<-read_csv(sum_folder%+%'norm_pd.txt')
f2<-read_csv(logd%+%'f2.csv')
f2_pg_pd<-filter(f2,exp=='pg_pd')
f2_pg_rg<-filter(f2,exp=='pg_rg')

# filter(!(convert_names(name)%in% no_mating))

for(table in c('norm_pd')){#,'norm_yd',#'f2_pg_pd','f2_pg_rg'
  m<-'delta'
  pvalue<-'pvalue'
  to_volcano<-get(table) %>%
    # mutate(!!sym(m):=100*!!sym(m))%>%
    # filter(exp==i)%>%
    transmute(
      name = name,
      logfc = !!sym(m),
      pval = !!sym(pvalue)
    ) %>%
    filter(!is.na(logfc) & !is.na(pval))
  
  to_volcano%>%
    filter(!(convert_names(name)%in%no_mating))%>%
    volcano_plot(
      data = .,
      ylabel=expression(-Log[10] ~ (P[value])),
      xlabel = expression(Log[2] ~ (italic(FoldChange))),#expression(frac(Delta*mu - median(Delta*mu), median(Delta*mu))~", %"),
      x_col = "logfc",
      y_col = "pval",
      center = mean(to_volcano$logfc, na.rm = TRUE),#0
      FCcutoff = 0.3,
      pCutoff = 0.05,
      pointSize = 2,
      point_alpha = 0.5,
      # xlim = c(-2.5, 2.5),
      title = toupper(table)  # добавляем заголовок
    )#+
    # geom_point(
    #   data = to_volcano %>% filter(convert_names(name) %in% no_mating),
    #   aes(x = logfc, y = -log10(pval)),
    #   color = "aquamarine",
    #   size = 3,
    #   stroke = 1,
    #   shape = 21  # круг с обводкой
    # ) +
    
    # geom_point(
    #   data = to_volcano %>% filter(name %in% genes),
    #   aes(x = logfc, y = -log10(pval)),
    #   color = "yellow",
    #   size = 3,
    #   stroke = 1,
    #   shape = 21  # круг с обводкой
    # ) +
    # ggrepel::geom_label_repel(
    #   data = to_volcano %>% filter(name %in% genes),
    #   aes(x = logfc, y = -log10(pval), label = name),
    #   size = 3,
    #   color = "black",
    #   fill = "yellow",
    #   point.padding = 0.3,
    #   force=100,
    #   direction = "both",
    #   min.segment.length = 0
    # )
  
  ggsave('figures/'%+%paste(table,'volcano','.png',sep='_'),width = 1280, height = 1280, units = 'px', scale=1.5)
}

res_cont<-read_tsv(sum_folder%+%'norm_pd.txt')
res_cont%>%
  filter(!(convert_names(name)%in%no_mating))%>%
  filter(abs(delta)>0.3, pvalue<0.05)%>%
  arrange(delta>0,pvalue)%>%
  mutate(across(where(is.numeric), ~round(.x,digits=3)))%>%
  transmute(name,delta,pvalue, gene_id=convert_names(name))%>%
  left_join(annot%>%select(gene_id,description),by='gene_id')%>%
  select(!gene_id)%>%write_tsv(sum_folder%+%'stats_hits.txt')
