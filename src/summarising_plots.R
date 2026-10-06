# ==============================================================================
# Barplots for total_reads/barcoded, boxplots for derived metrics,
# repeat correlations per condition, ORF detection barplot.
# ==============================================================================

library(tidyverse)
library(patchwork)

"%+%" <- function(...){ paste0(...) }

workingDIR <- "C:/Users/zokmi/space/2026_article_suppr/"
setwd(workingDIR)

sum_folder <- 'data/'
ref_folder <- 'ref/'
fig_folder <- 'figures/'
dir.create(fig_folder, showWarnings = FALSE)

# -------- read
sumtable <- read_csv(sum_folder %+% 'sumtable.csv', show_col_types = FALSE)


# jcounts with replicates
jcounts_rep <- read_csv(sum_folder %+% 'jcounts_2026_article_suppr.csv',
                        show_col_types = FALSE)

condition_map <- c(
  "yd" = "1n-all (YPD)",
  "yg" = "1n-Grande (YPGly)",
  "rd" = "rho0->2n-all (YPD)",
  "rg" = "rho0->2n-Grande (YPGly)",
  "pd" = "rho- ->2n-all (YPD)",
  "pg" = "rho- ->2n-Grande (YPGly)"
)

# ==============================================================================
# 1. Barplot: total_reads with total_barcoded inside, per sample
# ==============================================================================

bar_reads <- sumtable %>%
  mutate(
    barcoded_all = barcoded + barcoded_dm,
    other        = total_reads - barcoded_all
  ) %>%
  dplyr::select(exp, total_reads, barcoded_all, other) %>%
  pivot_longer(c(barcoded_all, other),
               names_to = 'class', values_to = 'reads') %>%
  mutate(class = factor(class,
                        levels = c('other', 'barcoded_all'),
                        labels = c('not barcoded', 'barcoded')))


combos <- expand.grid(
  Var4 = c(" F", " R"),        
  Var2 = " rep ",
  Var3 = c(1, 2, 3),
  Var1 = unname(condition_map),
  stringsAsFactors = FALSE
)

my_order <- paste0(combos$Var1, combos$Var2, combos$Var3, combos$Var4)

bar_reads <- bar_reads %>%
  mutate(exp=str_remove(exp,'_L00_R')%>%tolower(),
         condition_code = str_sub(exp, 1, 2),
         condition_label = recode(condition_code, !!!condition_map),
         exp=condition_label%+%" rep "%+%str_sub(exp,-2,-2)%+%ifelse(
           str_sub(exp,-1,-1)=='1'," F", " R"
         ))%>%
  mutate(exp=factor(exp, levels = rev(my_order)))

p_reads <- ggplot(bar_reads, aes(x = exp, y = reads, fill = class)) +
  geom_col(color = 'black', linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = c('not barcoded' = 'grey70',
                               'barcoded'     = '#F06C98')) +
  labs(x = 'Sample', y = 'Reads', fill = NULL,
       title = 'Total reads and barcoded fraction per sample') +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    legend.position = 'top'
  )+
  scale_y_continuous(labels=function(x)paste0(x/1e6,'M'))+
  coord_flip()

ggsave(fig_folder %+% 'sumtable_reads_bar.svg', p_reads,
       width = 12, height = 6, dpi = 300)
# ==============================================================================
# 2. Boxplots of derived metrics
# ==============================================================================

cor_mat <- jcounts_rep %>%
  dplyr::select(-name) %>%
  cor(method = 'spearman')

repeat_cor <- cor_mat %>%
  as.data.frame() %>%
  rownames_to_column('sample_a') %>%
  pivot_longer(-sample_a, names_to = 'sample_b', values_to = 'rho') %>%
  mutate(condition_a = str_sub(sample_a, 1, 2),
         condition_b = str_sub(sample_b, 1, 2)) %>%
  filter(condition_a == condition_b,
         sample_a < sample_b) %>%          
  transmute(condition = condition_a,
            metric = 'repeat correlation (Spearman)',
            value = rho)
# --- derived metrics per sample ---
metrics_plot <- sumtable %>%
  mutate(condition = str_sub(exp, 1, 2)) %>%
  transmute(
    condition,
    `yield`                          = yield,
    `assigned / barcoded_total`      = p_assigned_of_barcoded
  ) %>%
  pivot_longer(-condition, names_to = 'metric', values_to = 'value')

# --- combine ---
all_metrics <- bind_rows(metrics_plot, repeat_cor) %>%
  mutate(metric = factor(metric,
                         levels = c('yield',
                                    'assigned / barcoded_total',
                                    'repeat correlation (Spearman)')))

p_metrics <- ggplot(all_metrics, aes(x = metric, y = value, fill = metric)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = 'Value', title = 'Quality metrics') +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold'),

    legend.position = 'none',
    axis.text.x = element_text(size = 18),
    axis.text.y=element_text(size = 18),
    axis.title.y = element_blank()
  )+
  scale_x_discrete(labels=c('yield',
                            'assigned /\nbarcoded_total',
                            'repeat correlation\n(Spearman)'))

ggsave(fig_folder %+% 'sumtable_metrics_box.svg', p_metrics,
       width = 8, height = 5, dpi = 300)



# ==============================================================================
# 3. ORF detection barplot (count != 0 in jcounts)
#    Horizontal boxes, condition labels mapped to readable names.
# ==============================================================================

condition_map <- c(
  "yd" = "1n-all (YPD)",
  "yg" = "1n-Grande (YPGly)",
  "rd" = "rho0->2n-all (YPD)",
  "rg" = "rho0->2n-Grande (YPGly)",
  "pd" = "rho- ->2n-all (YPD)",
  "pg" = "rho- ->2n-Grande (YPGly)"
)

# (a) per replicate
orf_per_rep <- jcounts_rep %>%
  mutate(across(!name, ~as.integer(.x > 0))) %>%
  dplyr::select(!name) %>%
  colSums() %>%
  enframe(name = 'sample', value = 'n_detected') %>%
  mutate(condition = str_sub(sample, 1, 2),
         condition = recode(condition, !!!condition_map))

p_orf_rep <- ggplot(orf_per_rep, aes(x = n_detected, y = condition, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  geom_jitter(height = 0.15, size = 1.5, alpha = 0.7) +
  labs(x = 'ORFs detected (count > 0)', y = NULL,
       title = 'Detected ORFs per replicate') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
        legend.position = 'none')

# (b) detected in >= 2 of 3 replicates
orf_consensus <- jcounts_rep %>%
  mutate(across(!name, ~as.integer(.x > 0))) %>%
  pivot_longer(-name, names_to = 'sample', values_to = 'detected') %>%
  mutate(condition = str_sub(sample, 1, 2)) %>%
  group_by(condition, name) %>%
  summarise(n_rep = sum(detected), .groups = 'drop') %>%
  filter(n_rep >= 2) %>%
  count(condition, name = 'n_detected') %>%
  mutate(condition = recode(condition, !!!condition_map))

p_orf_cons <- ggplot(orf_consensus, aes(x = n_detected, y = condition, fill = condition)) +
  geom_col(color = 'black', linewidth = 0.3, width = 0.7) +

  geom_text(aes(label = n_detected), hjust = -0.2, size = 6) +

  labs(x = 'ORFs detected (>= 2 of 3 reps)', y = NULL,
       title = 'Consensus detected ORFs') +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
        legend.position = 'none',
        axis.text.y = element_text(size=18,color = 'black'),
        axis.text.x = element_text(color='black'))+
  coord_cartesian(xlim=c(0,3900))

ggsave(fig_folder %+% 'jcounts_orf_consensus.svg', p_orf_cons,
       width = 10, height = 8, dpi = 300)

p_orf <- p_orf_rep / p_orf_cons

ggsave(fig_folder %+% 'jcounts_orf_detection.svg', p_orf,
       width = 10, height = 9, dpi = 300)

# ==============================================================================
## yeast_phenome_correlation
# ==============================================================================

ypg<-read_tsv(ref_folder%+%'ypg16489.txt')[,3:4]

jcounts<-jcounts_rep

#reference NPVs

npvise<-function(col){
  col<-ifelse(col!=0,log(col),NA)
  dens<-density(col,bw=0.25,kernel='gaussian',na.rm=T)
  mode<-dens$x[which.max(dens$y)]
  dif<-col-mode
  sq<-dif*dif
  NPV<-ifelse(is.na(col),NA,dif/sqrt(1/nrow(na.omit(jcounts))*sum(sq,na.rm=T)))
  return(NPV)
}


colnames(ypg)<-c('name','yg')
jcounts%>%dplyr::select(name, starts_with('y'))->g
n <- mean(colSums(g%>%dplyr::select(where(is.numeric))))
g <- g %>%
  mutate(across(!name, ~ .x / sum(.x) * n))

library(multimode)

g %>%
  pivot_longer( - name, names_to = 'exp')%>%
  group_by(name)%>%
  mutate(average = mean(value))%>%
  pivot_wider(names_from = 'exp')%>%
  ungroup()%>%
  mutate(mode=density(average,bw=0.25,kernel='gaussian',na.rm=T)$x[which.max(density(average,bw=0.25,kernel='gaussian',na.rm=T)$y)])%>%
  mutate(dif=average - mode)%>%
  mutate(sq = sqrt(sum(dif*dif, na.rm = T) / nrow(ypg)))%>%
  full_join(ypg)%>%
  mutate(yg = yg * sq + mode)%>%
  na.omit%>%
  pull('yg')%>%hist(breaks = 100)

(g%>%
    pivot_longer( - name, names_to = 'exp')%>%
    group_by(name)%>%
    mutate(average = mean(value, na.rm = T))%>%
    pivot_wider(names_from = 'exp')%>%
    ungroup()%>%full_join(ypg)%>%
    mutate(across(!name & !yg, npvise))%>%
    dplyr::select(average, yg)%>%
    na.omit%>%
    with(cor.test(pull(., average), pull(., yg), method = 'p')))%>%names()
g%>%
  pivot_longer( - name, names_to = 'exp')%>%
  group_by(name)%>%
  mutate(average = mean(value, na.rm = T))%>%
  pivot_wider(names_from = 'exp')%>%
  ungroup()%>%full_join(ypg)%>%
  mutate(across(!name & !yg, npvise))%>%
  dplyr::select(average, yg)%>%
  na.omit%>%ggplot()+
  aes(x = average, y = yg)+
  geom_point(color = 'red3', alpha = 0.2)+
  theme_bw()+
  labs(y = "reference YKOC YPGly fitness, normalised", x = 'log( average count ) on YPGly, normalised')
ggsave(fig_folder%+%'yph_cor.svg', height = 1080, width = 1920, units = 'px')



aa_test<-g%>%
  pivot_longer( - name, names_to = 'exp')%>%
  group_by(name)%>%
  mutate(average = mean(value, na.rm = T))%>%
  pivot_wider(names_from = 'exp')%>%
  ungroup()%>%full_join(ypg)%>%
  mutate(across(!name & !yg, npvise))%>%
  dplyr::select(name,average, yg)%>%
  na.omit%>%
  mutate(ref_group = ifelse(yg< -1,'low','high'))

testres<-t.test(filter(aa_test, ref_group=='low')$average, filter(aa_test, ref_group=='high')$average)

ggplot(aa_test, aes(x = ref_group, y = average, fill = ref_group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  scale_fill_manual(values = c('low' = '#E64B35', 'high' = '#4DBBD5')) +
  scale_x_discrete(labels = c('low' = 'fitness < -1', 'high' = 'fitness > -1')) +
  labs(x = NULL, y = 'log( average count ) on YPGly') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(size=18,color='black'))+
  annotate('text', x = 1.5, y = max(aa_test$average, na.rm = TRUE) * 1.05,
           label = paste0('Welch t-test, p = ', format.pval(testres$p.value, digits = 2)))


ggsave(fig_folder%+%'yph_test.svg', height = 1080, width = 1920, units = 'px')

