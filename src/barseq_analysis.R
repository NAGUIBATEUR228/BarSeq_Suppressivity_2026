#https://www.bioconductor.org/
library(tidyverse)
library(Biostrings)
library(ShortRead)
# for faster processing
library(data.table)


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

#write down the extension of fastq files
typefq <- '.gz'
#write down the substring to be removed to get shorter output filenames if needed
excess <- c('.fq', '202304172125_211101003_2P221122075US2S2691BX_A_428_Knorre_PE75_17_04_2023_', 'L00_')

"%+%" <- function(...){
  paste0(...)
}

remove_excess <-function(string, excess){
  Reduce(function(x, y) gsub(y, "", x), excess, init = string)
}

# reverse_complement_string
rev_compl_string <- function(x) {
  as.character(reverseComplement(DNAStringSet(x)))
}

#read file with reference barcodes of YKOC (S. genome deletion project)
ref_nm <- read_tsv(ref_folder %+% "sgdp.txt") %>%
  transmute(Confirmed_deletion = ORF,
            UPTAG_seqs = UPTAG_seq,
            UPTAG_notes = 'YKOC')

#read file with Reference reads file generated from Puddu et al. 2019
reference <- read_tsv(ref_folder%+%"puddu.txt")%>%
  select("Confirmed_deletion","UPTAG_seqs","UPTAG_seqs1",'UPTAG_notes')%>%
  filter(!is.na(UPTAG_seqs)|!is.na(UPTAG_seqs1))%>%
  unique()%>%
  pivot_longer(starts_with("UPTAG_seqs"),values_to = "UPTAG_seqs")%>%
  select(!'name')%>%
  filter(!is.na(UPTAG_seqs))%>%unique()%>% #filter(if_any(everything(),~str_detect(.x,"\\|")))
  full_join(.,ref_nm)%>%
  pivot_longer(cols=c(Confirmed_deletion,UPTAG_notes))%>%
  pivot_wider(values_fn = ~str_flatten(.x%>%unique,collapse = '|'))%>%
  mutate(UPTAG_notes=ifelse(is.na(UPTAG_notes),'puddu',UPTAG_notes))

write_csv(reference,ref_folder%+%"reference.txt")

u1 <- DNAString("GATGTCCACGAGGTCTCT")
u1_rc <- reverseComplement(u1)

u2 <- DNAString("CGTACGCTGCAGGTCGAC")
u2_rc <- reverseComplement(u2)

fastqdirectories <- list.dirs(path = workingDIR%+%folder, recursive = T)

fastqfiles <- unlist(map(fastqdirectories, 
                         ~ paste(.x, dir(.x), sep='')))%>%
  str_subset(typefq%+%'$')
sumtable_short<-tibble(name=character(), total_count = numeric())

for (j in fastqfiles){
  setwd(workingDIR)
  strm <- FastqStreamer(j)
  
  cat(format(Sys.time(), "%H:%M:%S"), j, '\n')
  
  # to collect all data in one data.table
  all_data <- NULL
  
  repeat {
    fq <- yield(strm)
    if (length(fq) == 0) {break}
    
    # estimate mean quality for each read
    quals <- rowMeans(as(quality(fq), "matrix"))
    seqs <- as.character(sread(fq))
    
    # make temp data.table
    dt <- data.table(
      seq = seqs,
      qual = quals
    )
    
    if (is.null(all_data)) {
      all_data <- dt
    } else {
      all_data <- rbindlist(list(all_data, dt))
    }
  }
  
  # group and aggregate after data processing
  out_table <- all_data[, .(
    count = .N,
    qual = mean(qual)
  ), by = seq]
  
  raw_reads_count <- out_table[order(-count)]
  
  
  
  raw_reads_count <- raw_reads_count%>% as_tibble()
  
  dir.create(workingDIR %+% output_folder, showWarnings = FALSE)
  setwd(workingDIR %+% output_folder)
  i <- '../' %+% str_remove(j, workingDIR)
  name_for_files <- i%>%
    remove_excess(c(excess,typefq))%>%
    str_replace(folder, output_folder)
  
  total_counts <- sum(raw_reads_count$count)
  
  barcoded_reads <- raw_reads_count%>% 
    filter(
      str_detect(seq, u1%+%"\\s*(.*?)\\s*"%+%u2)|
        str_detect(seq, u2_rc%+%"\\s*(.*?)\\s*"%+%u1_rc)) %>%
    mutate(barcode = str_match(seq, u1%+%"\\s*(.*?)\\s*"%+%u2)[,2]) %>%
    mutate(barcode = ifelse(
      is.na(barcode),
      str_match(seq,u2_rc%+%"\\s*(.*?)\\s*"%+%u1_rc)[,2],
      barcode))%>%
    filter(!is.na(barcode))%>%#if u1u2? with no barcode
    mutate(barcode = ifelse(str_detect(seq, "GATGTCCACGAGGTCTCT"),
                            barcode, 
                            rev_compl_string(barcode)))
  
  write_csv(raw_reads_count%>% 
              filter(!(
                str_detect(seq, u1%+%"\\s*(.*?)\\s*"%+%u2)|str_detect(seq, u2_rc%+%"\\s*(.*?)\\s*"%+%u1_rc))), 
            paste(name_for_files, 
                  "not_barcoded_raw_qual_count.csv", 
                  sep = "_"))
  
  final_result <- barcoded_reads %>% 
    mutate(total_qual=count*qual)%>%
    group_by(barcode) %>%
    summarise(n = n(), count = sum(count),qual=sum(total_qual)/sum(count)) %>%
    ungroup()%>%
    arrange(desc(count))#column n means number of unique reads with barcode
  
  #table with information about quality of barcodes
  write_csv(final_result, 
            paste(name_for_files, 
                  "barcoded_qual_count.csv", 
                  sep = "_"))
  
  
  output <- reference %>% 
    left_join(.,final_result, by = c("UPTAG_seqs" = "barcode"))%>% 
    filter(!is.na(count))%>%#some barcodes from final result are not matched
    group_by(Confirmed_deletion)%>%
    summarise(barcode=str_flatten(unlist(str_split(UPTAG_seqs,'\\|'))%>%unique,collapse = '|'),n=sum(n),count=sum(count),notes=str_flatten(unlist(str_split(UPTAG_notes,'\\|'))%>%unique,collapse = '|'))%>%
    ungroup()%>%
    arrange(desc(count))
  
  write_csv(final_result%>%
              filter(!(barcode %in% reference$UPTAG_seqs))%>%
              arrange(desc(count)),
            paste(name_for_files, 
                  "not_matched.csv", sep = "_"))
  
  hist<-output %>% 
    ggplot() +
    aes(x = count)+
    geom_histogram(
      binwidth = 10,
      fill='#F06C98',
      color='black')+
    theme_bw()+
    labs(x = "Barcode read counts", y="Number of barcode with such read count")+
    scale_y_log10()
  
  ggsave(paste(name_for_files, 
               "hist.png", sep = "_"),
         hist,
         width = 1280,
         height = 720,
         units = "px",scale=2)
  
  write_csv(output,
            paste(name_for_files, 
                  "output_count.csv", sep = "_"))
  
  sumtable_short<-sumtable_short%>%add_row(name = name_for_files, total_count = total_counts)
}
write_csv(sumtable_short,workingDIR%+%sum_folder%+%'sumtable_short.txt')

setwd(workingDIR%+%'..')
system('python3 src/barcode_dm.py')
system('python3 src/barcode_blast.py')