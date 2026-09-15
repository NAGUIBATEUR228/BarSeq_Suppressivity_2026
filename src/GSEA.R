# ==============================================================================
# GSEA_final.R
# Анализ обогащения GO терминов для данных по супрессорам
# ==============================================================================

# 1. Загрузка библиотек ------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(org.Sc.sgd.db)

# 2. Настройка окружения и путей ---------------------------------------
workingDIR <- "C:/Users/zokmi/space/2026_article_suppr/"
setwd(workingDIR)

folder <- 'seq_data/'
output_folder <- 'output/'
sum_folder <- 'data/'
ref_folder <- 'ref/'

# Оператор конкатенации строк
"%+%" <- function(...){ paste0(...) }

# 3. Загрузка данных ------------------------------------------------
# Загружаем финальные нормализованные данные
jcounts <- read_tsv(sum_folder %+% 'jcounts_final.txt')
jc <- jcounts

# Нормализация по медиане общего количества ридов (как в оригинальном анализе)
n = jc %>% dplyr::select(!name) %>% colSums %>% median
jc <- jc %>% mutate(across(!name, ~ .x / sum(.x) * n))

# Загрузка аннотаций
yeast_genes <- read_tsv(ref_folder %+% 'yeast_gene_names.txt')
annot <- read_tsv(ref_folder %+% "gene_annotation.txt")
go <- read_tsv(ref_folder %+% "go_annotation.txt")
all_go_descriptions <- read_tsv(ref_folder %+% "go_description.txt")

# 4. Вспомогательные функции ----------------------------------------
# Функция конвертации имен генов
convert_names <- function(vals, to_std = F){
  if (length(vals) == 0){return(vals)}
  l <- tibble(syn = vals, id = (1:length(vals))) %>%
    separate_longer_delim(syn, '|') %>%
    left_join(yeast_genes, by = 'syn')
  if (to_std) {
    return(l %>%
             mutate(res = ifelse(is.na(std_name), syn, std_name)) %>%
             group_by(id) %>%
             summarise(res = str_flatten(res %>% unique, collapse = '|')) %>%
             pull(res))
  }
  return(l %>%
           mutate(res = ifelse(is.na(sys_name), syn, sys_name)) %>%
           group_by(id) %>%
           summarise(res = str_flatten(res %>% unique, collapse = '|')) %>%
           pull(res))
}

# 5. Подготовка данных для GSEA -------------------------------------
# Исключаем гены спаривания (как в последнем блоке)
no_mating <- go %>% filter(go_term == 'GO:0000747') %>% pull(gene_id)

# Создаем таблицу `to_go` для анализа

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

to_go <- data %>%
  # filter(if_all(starts_with('yd'), ~ .x > 10)) %>%
  pivot_longer(!name, names_to = 'exp') %>%
  mutate(exp = str_sub(exp, 1, -2)) %>%
  group_by(name, exp) %>%
  # mutate(value=log2(value))%>%
  mutate(value = mean(value, na.rm = T)) %>%
  ungroup() %>%
  distinct() %>%
  pivot_wider(names_from = 'exp')

# 6. Основной цикл GSEA ---------------------------------------------
# Инициализируем пустой датафрейм для результатов
GO_norm <- tibble(pvalue = numeric(), 
                  Description = character(), 
                  condition = character(), 
                  NES = numeric(), 
                  geneID = character(),
                  coverage = numeric(),
                  setSize = numeric(), qvalue = numeric())

# Запускаем GSEA для пар условий
# Выполняем для 'pg_rg' и 'pg_pd'
for (i in c('pg_rg', 'pg_pd')) {
  condition <- i
  t <- str_split_1(condition, '_')[1]  # первое условие (например, 'pg')
  u <- str_split_1(condition, '_')[2]  # второе условие (например, 'rg')
  
  # Подготовка данных для GSEA
  res <- to_go %>%
    dplyr::select(name, starts_with(t) | starts_with(u)) %>%
    mutate(log2FoldChange = !!sym(t) - !!sym(u)) %>%
    mutate(gene = name) %>%
    dplyr::select(gene, log2FoldChange)
  
  # Исключаем гены с "|" (коллизии)
  cut <- res %>% filter(!str_detect(gene, '\\|'))
  
  # Создаем ранжированный список генов
  gene_list <- cut$log2FoldChange
  names(gene_list) <- cut$gene %>% convert_names
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Запускаем GSEA
  # Важно: pvalueCutoff = 1, pAdjustMethod = "none", как в истории
  norm_gsea2 <- clusterProfiler::gseGO(
    geneList = gene_list,
    OrgDb = org.Sc.sgd.db,
    pvalueCutoff = 0.05,
    keyType = "ENSEMBL",
    ont = "BP",      # ← Оставляем только значимые
    pAdjustMethod = "none",       # ← Сразу коррекция
    minGSSize = 10,             # ← Минимальный размер термина
    maxGSSize = 200,            # ← Максимальный размер термина
  )
  
  norm_gsea2@result$gene_count <- map_int(norm_gsea2@result$core_enrichment, ~ {
    if(is.na(.x)) return(0)
    str_count(.x, "/") + 1
  })
  

  norm_gsea2@result$coverage <- (norm_gsea2@result$gene_count / 
                                   norm_gsea2@result$setSize) * 100
  
  norm_gsea<-clusterProfiler::simplify(
    norm_gsea2,  # Ваш результат gseGO
    cutoff = 0.7,  # Порог сходства (0.7 = умеренный)
    by = "pvalue",
    select_fun = min  
  )
  
  # Сохраняем результаты, если они есть
  if ("core_enrichment" %in% names(norm_gsea@result)) {
    gene_clust <- norm_gsea@result$core_enrichment %>%
      map(~ convert_names(str_split(., '\\/') %>% unlist(), to_std = TRUE)) %>%
      map(~ str_flatten(., collapse = '/')) %>% as.character()
    
    GO_norm <- norm_gsea@result %>%
      # mutate(
      #   # Количество генов из вашего списка в термине
      #   gene_count = map_int(core_enrichment, ~ str_count(.x, "/") + 1),
      #   
      #   # Общее количество генов в GO термине
      #   total_genes = setSize,
      #   
      #   # Насыщенность (доля)
      #   coverage = gene_count / total_genes
      # ) %>%
      dplyr::select(Description, pvalue, coverage, setSize, NES, qvalue) %>%
      mutate(condition = condition, geneID = gene_clust) %>%
      add_row(GO_norm, .)
  }
}


GO_all%>%
  mutate(padj=p.adjust(pvalue,'BH'))%>%
  dplyr::select(Description, condition,padj)%>%
  right_join(GO_norm)%>%write_tsv(sum_folder%+%"GO_simple_unified.txt")

GO_all%>%
  mutate(padj=p.adjust(pvalue,'BH'))%>%
  write_tsv(sum_folder%+%'GO_no_correction_unified.txt')

GO_norm<-read_tsv(sum_folder%+%'GOsimple.txt')%>%
  filter(condition=='pg_pd')

# Упорядочиваем термины по NES
GO_norm %>%
  arrange(desc(NES)) %>%
  pull(Description) %>%
  unique() -> lvls


p<-to_pic %>%
  ggplot(aes(y = factor(Description, levels = rev(lvls)),
             x = NES,
             fill = ifelse(NES > 0, "Обогащен", "Обеднен"))) +
  geom_col(width = 0.7) +
  labs(y = "GO-terms",
       x = "NES (Normalized Enrichment Score)",
       title = "GO-term enrichment") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_fill_manual(name = "Направление",
                    values = c("Обогащен" = "steelblue", "Обеднен" = "salmon"),
                    labels = c("Обогащен" = "Положительный NES", "Обеднен" = "Отрицательный NES")) +
  coord_fixed(ratio = 0.5)

# Сохраняем график
ggsave('figures/GSEA_simple_pg_pd.png', plot = p, scale = 2.5)
# ggsave('figures/GSEA_pg_pd.svg', plot = p, scale = 2.5)