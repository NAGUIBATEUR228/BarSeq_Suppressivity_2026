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
n = jc %>% dplyr::select(!name) %>% colSums %>% mean
jc <- jc %>% mutate(across(!name, ~ .x / sum(.x) * n))

# Загрузка аннотаций
yeast_genes <- read_tsv(ref_folder %+% 'yeast_gene_names.txt')
annot <- read_tsv(ref_folder %+% "gene_annotation.txt")
go <- read_tsv(ref_folder %+% "go_annotation.txt")
all_go_descriptions <- read_tsv(ref_folder %+% "go_description.txt")

# # Создаем полную аннотацию: ген → GO термин
# gene_to_go <- go %>%
#   left_join(all_go_descriptions, by = c("go_term" = "go_id")) %>%
#   filter(!is.na(go_description))  # Убираем NA
# # Теперь у вас есть:
# # gene_id → go_term (GO ID) → go_description (название термина) + ontology
# # Для каждого GO термина получаем список генов
# go_genes <- gene_to_go %>%
#   group_by(go_description) %>%
#   summarise(
#     genes = list(unique(gene_id)),
#     n_genes = n()
#   ) %>%
#   filter(n_genes >= 3)  # Минимальный размер
# # Функция проверки, является ли термин A подмножеством B
# is_subset <- function(genes_A, genes_B) {
#   all(genes_A %in% genes_B)
# }

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
to_go <- jc %>%
  # Оставляем гены с суммой ридов по репликам yd > 10 (как в последней версии)
  filter(if_all(starts_with('yd'), ~ .x > 10)) %>%
  pivot_longer(!name, names_to = 'exp') %>%
  mutate(exp = str_sub(exp, 1, -2)) %>%
  group_by(name, exp) %>%
  filter(value != 0) %>%
  mutate(value = mean(log2(value + 1), na.rm = T)) %>%
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


GO_norm%>%
  # filter(
  #   # 1. Минимальный размер термина (3-5 генов)
  #   setSize >= 5,
  #   
  #   # 2. Максимальный размер (слишком общие термины)
  #   setSize <= 200,
  #   coverage*setSize>=4
  # )%>%
  mutate(padj=p.adjust(pvalue, method='BH'))%>%
  filter(condition=='pg_pd')%>%
  view()

# 7. Сохранение результатов -----------------------------------------
# Сохраняем результаты для 'pg_pd' с p-value < 0.01 (как в финальной версии)
# GO_norm %>%
#   filter(pvalue < 0.01, condition == 'pg_pd') %>%
#   write_tsv(sum_folder %+% 'GOpgpd.txt')

# 8. Создание графика (фильтрация и визуализация) --------------------
# Загружаем сохраненные данные для pg_pd
# GO_norm <- read_tsv(sum_folder %+% 'GOpgpd.txt')
GO_norm%>%write_tsv(sum_folder%+%'GOsimple.txt')

# # Функция для фильтрации избыточных GO терминов (по подмножеству генов)
# genesets <- GO_norm$geneID %>% str_split('/')
# 
# # Оставляем только термины, которые НЕ являются надмножеством для других
# to_keep <- sapply(1:nrow(GO_norm), function(i) {
#   !any(sapply(1:nrow(GO_norm), function(j) {
#     i != j &&
#       all(genesets[[j]] %in% genesets[[i]]) &&  # j подмножество i
#       length(genesets[[j]]) < length(genesets[[i]]) # j строго меньше
#   }))
# })
# GO_norm <- GO_norm[to_keep, ]

GO_norm<-read_tsv(sum_folder%+%'GOsimple.txt')%>%
  filter(condition=='pg_pd')

# Упорядочиваем термины по NES
GO_norm %>%
  arrange(desc(NES)) %>%
  pull(Description) %>%
  unique() -> lvls

# Строим финальный график
p <- GO_norm %>%
  # mutate(pvalue = p.adjust(pvalue, method = 'BH')) %>% # В финальной версии этот шаг пропущен
  # filter(pvalue < 0.01) %>% # В финальной версии фильтр по p-value не применяется на графике
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