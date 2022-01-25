```{r}
library(stringr)
library(tidyverse)
library(ggpubr)
library(foreach)
library(openxlsx)
library(doParallel)
library(readxl)
library(VennDiagram)
libs = c('stringr','tidyverse','openxlsx','foreach')

cl <- makeCluster(16)
registerDoParallel(cl) 
```

---

## Downdoad annotation files

```{r}

Mouse_miRBase_annotation = read.table('mouse/mmu.gff3',stringsAsFactors = F)
Mouse_miRBase_annotation$ID = foreach(i = 1:nrow(Mouse_miRBase_annotation),.combine = 'c') %do% {
  str_remove(strsplit(Mouse_miRBase_annotation$V9[i],';')[[1]][3],'Name=')
}
Mouse_piRNAdb_annotation = read.table('mouse/pirnadb.v1_7_6.mm10.gtf',stringsAsFactors = F)
Mouse_GtRNAdb_annotation = read.table('mouse/mm39-tRNAs.bed',stringsAsFactors = F)
Mouse_rRNA_ucsc_annotation = read.table('mouse/mouse_rRNA_ucsc.bed',stringsAsFactors = F)
Mouse_trfdb_raw_table = rbind.data.frame(
  read.xlsx('mouse/trfdb_excel_result_1.xlsx'),
  read.xlsx('mouse/trfdb_excel_result_2.xlsx'),
  read.xlsx('mouse/trfdb_excel_result_3.xlsx')
  )
write.table(
  unique(str_c('>',Mouse_trfdb_raw_table$Type,'_',Mouse_trfdb_raw_table$tRF.ID,'\n',Mouse_trfdb_raw_table$tRF.Sequence)),'mouse/tRFdb.fasta',row.names = F,col.names = F
)

```


---
## Build bed-files for micro, piwi, trf RNA -> liftover -> getfasta

```{r}
write.table(Mouse_miRBase_annotation[,c(1,4,5,10,7)],'mouse/Mouse_miRBase_mm10.bed', col.names = F,
            row.names = F, sep = '\t')
write.table(Mouse_piRNAdb_annotation[1:1000000,c(1,4,5,10,7)],'mouse/Mouse_piRNAdb_mm10_1.bed', col.names = F,
            row.names = F, sep = '\t')
write.table(Mouse_piRNAdb_annotation[1000001:2000000,c(1,4,5,10,7)],'mouse/Mouse_piRNAdb_mm10_2.bed', col.names = F,
            row.names = F, sep = '\t')
write.table(Mouse_piRNAdb_annotation[2000001:nrow(Mouse_piRNAdb_annotation),c(1,4,5,10,7)],'mouse/Mouse_piRNAdb_mm10_3.bed', col.names = F,
            row.names = F, sep = '\t')
Mouse_trfdb_annotation = foreach(i = 1:nrow(Mouse_trfdb_raw_table),.combine = 'rbind.data.frame') %do% {
  p = unlist(c(strsplit(Mouse_trfdb_raw_table$`tRNA.Gene.Co-ordinates`[i],'-')[[1]],
  Mouse_trfdb_raw_table[i,c(1,3,8)]))
  x = p[2]
  if (as.numeric(p[2]) < as.numeric(p[3])) {
    p[2] = as.numeric(x) - 2 + as.numeric(strsplit(str_remove(p[6],' Start:'),' -')[[1]][1])
  p[3] = as.numeric(x) -1 + as.numeric(strsplit(str_remove(p[6],' '),' - End:')[[1]][2])
    o = c(p[1:5],'+')}
  else {
    p[2] = as.numeric(x) + 1 - as.numeric(strsplit(str_remove(p[6],' Start:'),' -')[[1]][1])
  p[3] = as.numeric(x) - as.numeric(strsplit(str_remove(p[6],' '),' - End:')[[1]][2])
    o = c(p[c(1,3,2,4,5)],'-')}
  names(o) = c('chr','start','end','ID','type','stand')
  t(as.data.frame(o))
}
Mouse_trfdb_annotation$chr = str_remove(Mouse_trfdb_annotation$chr,' ')

  
write.table(Mouse_trfdb_annotation
,'mouse/Mouse_tRFdb_mm9.bed', col.names = F,
            row.names = F, sep = '\t')


Mouse_miRBase_mm39_annotation = read.table('mouse/Mouse_miRBase_mm39.bed',stringsAsFactors = F)
Mouse_tRFdb_mm39_annotation = read.table('mouse/Mouse_tRFdb_mm39.bed',stringsAsFactors = F)
Mouse_piRNAdb_mm39_annotation = rbind.data.frame(
  read.table('mouse/Mouse_piRNAdb_mm38_1.bed',stringsAsFactors = F),
  read.table('mouse/Mouse_piRNAdb_mm38_2.bed',stringsAsFactors = F),
  read.table('mouse/Mouse_piRNAdb_mm38_3.bed',stringsAsFactors = F)
)
write.table(Mouse_piRNAdb_mm38_annotation,'mouse/Mouse_piRNAdb_mm38.bed', col.names = F,
            row.names = F, sep = '\t') 

Mouse_miRBase_getfasta = read.table('mouse/Mouse_miRBase_mm39.getfasta.tsv',stringsAsFactors = F)
colnames(Mouse_miRBase_getfasta) = c('Loci','Loci_seq')
Mouse_GtRNAdb_getfasta = read.table('mouse/mm39-tRNAs.getfasta.tsv',stringsAsFactors = F)
colnames(Mouse_GtRNAdb_getfasta) = c('Loci','Loci_seq')

Mouse_piRNAdb_getfasta_fa = read.table('mouse/Mouse_piRNAdb_mm39.getfasta',stringsAsFactors = F)
Mouse_piRNAdb_getfasta = 
  data.frame(Loci = Mouse_piRNAdb_getfasta_fa$V1[seq(1,length(Mouse_piRNAdb_getfasta_fa$V1), by = 2)],
             Loci_seq = Mouse_piRNAdb_getfasta_fa$V1[seq(2,length(Mouse_piRNAdb_getfasta_fa$V1), by = 2)],stringsAsFactors = F)

Mouse_piRNAdb_getfasta$Loci = str_remove_all(Mouse_piRNAdb_getfasta$Loci,'>')

Mouse_rRNA_ucsc_getfasta = read.table('mouse/mouse_rRNA_ucsc.getfasta.tsv',stringsAsFactors = F)
colnames(Mouse_rRNA_ucsc_getfasta) = c('Loci','Loci_seq')
Mouse_tRFdb_getfasta = read.table('mouse/Mouse_tRFdb_mm39.getfasta.tsv',stringsAsFactors = F)
colnames(Mouse_tRFdb_getfasta) = c('Loci','Loci_seq')

```

## Build loci databases

```{r}

mouse_databases = list()
mouse_databases$miRBase = data.frame(
  ID = Mouse_miRBase_mm39_annotation$V4,
  Loci = str_c(Mouse_miRBase_mm39_annotation$V1,':',Mouse_miRBase_mm39_annotation$V2,'-',Mouse_miRBase_mm39_annotation$V3),
  Loci_chr = Mouse_miRBase_mm39_annotation$V1,
  Loci_start = Mouse_miRBase_mm39_annotation$V2,
  Loci_end = Mouse_miRBase_mm39_annotation$V3,
  Loci_len = as.numeric(Mouse_miRBase_mm39_annotation$V3) - as.numeric(Mouse_miRBase_mm39_annotation$V2),
  Loci_strand = Mouse_miRBase_mm39_annotation$V5,
  stringsAsFactors = F
)
mouse_databases$piRNAdb = data.frame(
  ID = Mouse_piRNAdb_mm39_annotation$V4,
  Loci = str_c(Mouse_piRNAdb_mm39_annotation$V1,':',Mouse_piRNAdb_mm39_annotation$V2,'-',Mouse_piRNAdb_mm39_annotation$V3),
  Loci_chr = Mouse_piRNAdb_mm39_annotation$V1,
  Loci_start = Mouse_piRNAdb_mm39_annotation$V2,
  Loci_end = Mouse_piRNAdb_mm39_annotation$V3,
  Loci_len = as.numeric(Mouse_piRNAdb_mm39_annotation$V3) - as.numeric(Mouse_piRNAdb_mm39_annotation$V2),
  Loci_strand = Mouse_piRNAdb_mm39_annotation$V5,
  stringsAsFactors = F
)
mouse_databases$GtRNAdb = data.frame(
  ID = Mouse_GtRNAdb_annotation$V4,
  Loci = str_c(Mouse_GtRNAdb_annotation$V1,':',Mouse_GtRNAdb_annotation$V2,'-',Mouse_GtRNAdb_annotation$V3),
  Loci_chr = Mouse_GtRNAdb_annotation$V1,
  Loci_start = Mouse_GtRNAdb_annotation$V2,
  Loci_end = Mouse_GtRNAdb_annotation$V3,
  Loci_len = as.numeric(Mouse_GtRNAdb_annotation$V3) - as.numeric(Mouse_GtRNAdb_annotation$V2),
  Loci_strand = Mouse_GtRNAdb_annotation$V6,
  stringsAsFactors = F
)
mouse_databases$rRNA_ucsc = data.frame(
  ID = Mouse_rRNA_ucsc_annotation$V4,
  Loci = str_c(Mouse_rRNA_ucsc_annotation$V1,':',Mouse_rRNA_ucsc_annotation$V2,'-',Mouse_rRNA_ucsc_annotation$V3),
  Loci_chr = Mouse_rRNA_ucsc_annotation$V1,
  Loci_start = Mouse_rRNA_ucsc_annotation$V2,
  Loci_end = Mouse_rRNA_ucsc_annotation$V3,
  Loci_len = as.numeric(Mouse_rRNA_ucsc_annotation$V3) - as.numeric(Mouse_rRNA_ucsc_annotation$V2),
  Loci_strand = Mouse_rRNA_ucsc_annotation$V6,
  stringsAsFactors = F
)
mouse_databases$tRFdb = data.frame(
  ID = str_c(Mouse_tRFdb_mm39_annotation$V5,'_',Mouse_tRFdb_mm39_annotation$V4),
  Loci = str_c(Mouse_tRFdb_mm39_annotation$V1,':',Mouse_tRFdb_mm39_annotation$V2,'-',Mouse_tRFdb_mm39_annotation$V3),
  Loci_chr = Mouse_tRFdb_mm39_annotation$V1,
  Loci_start = as.numeric(as.vector(Mouse_tRFdb_mm39_annotation$V2)),
  Loci_end = as.numeric(as.vector(Mouse_tRFdb_mm39_annotation$V3)),
  Loci_len = as.numeric(as.vector(Mouse_tRFdb_mm39_annotation$V3)) - as.numeric(as.vector(Mouse_tRFdb_mm39_annotation$V2)),
  Loci_strand = Mouse_tRFdb_mm39_annotation$V6,
  stringsAsFactors = F
)

```

---
## Add getfasta sequences to loci -> filtering -> reverse compliment sequences

```{r}

mouse_databases$miRBase$Loci_seq = foreach(l = mouse_databases$miRBase$Loci,.combine = 'c') %dopar% {
    Mouse_miRBase_getfasta[which(l == Mouse_miRBase_getfasta$Loci)[1],'Loci_seq']}

mouse_databases$GtRNAdb$Loci_seq = foreach(l = mouse_databases$GtRNAdb$Loci,.combine = 'c') %dopar% {
    Mouse_GtRNAdb_getfasta[which(l == Mouse_GtRNAdb_getfasta$Loci)[1],'Loci_seq']}
mouse_databases$GtRNAdb = filter(mouse_databases$GtRNAdb,!is.na(Loci_seq))

mean(mouse_databases$piRNAdb$Loci == Mouse_piRNAdb_getfasta$Loci)
mouse_databases$piRNAdb$Loci_seq = Mouse_piRNAdb_getfasta$Loci_seq

mouse_databases$rRNA_ucsc$Loci_seq = foreach(l = mouse_databases$rRNA_ucsc$Loci,.combine = 'c') %do% {
    Mouse_rRNA_ucsc_getfasta[which(l == Mouse_rRNA_ucsc_getfasta$Loci)[1],'Loci_seq']}
mouse_databases$rRNA_ucsc = filter(mouse_databases$rRNA_ucsc,!is.na(Loci_seq))

mouse_databases$tRFdb$Loci_seq = foreach(l = mouse_databases$tRFdb$Loci,.combine = 'c') %do% {
    Mouse_tRFdb_getfasta[which(l == Mouse_tRFdb_getfasta$Loci)[1],'Loci_seq']}

mouse_databases$miRBase$Loci_seq = str_to_upper(mouse_databases$miRBase$Loci_seq)
mouse_databases$piRNAdb$Loci_seq = str_to_upper(mouse_databases$piRNAdb$Loci_seq)
mouse_databases$GtRNAdb$Loci_seq = str_to_upper(mouse_databases$GtRNAdb$Loci_seq)
mouse_databases$rRNA_ucsc$Loci_seq = str_to_upper(mouse_databases$rRNA_ucsc$Loci_seq)
mouse_databases$tRFdb$Loci_seq = str_to_upper(mouse_databases$tRFdb$Loci_seq)

rev_comp = function(f) {
  return(foreach(x = f,.combine = 'c',
                 .packages = c('stringr','foreach')) %dopar% {
    forward = strsplit(x,'')[[1]]
    foreach(i = length(forward):1,.combine = 'str_c') %do% {
    o = ''
      if (forward[i] == 'A'){o = 'T'}
      if (forward[i] == 'T'){o = 'A'}
      if (forward[i] == 'G'){o = 'C'}
      if (forward[i] == 'C'){o = 'G'}
    o
    }})
}

for (x in names(mouse_databases)[5]) {
  print(x)
  mouse_databases[[x]] = mutate(mouse_databases[[x]],Loci_rev_seq = rev_comp(Loci_seq))
}

```

## Download fasta and sam files

```{r}
Mouse_miRBase_fasta = rbind.data.frame(
  read.table('mouse/miRBase_hairpin.fasta.tsv',stringsAsFactors = F,sep = '\t'),
  read.table('mouse/miRBase_mature.fasta.tsv',stringsAsFactors = F,sep = '\t'))
colnames(Mouse_miRBase_fasta) = c('ID','Fasta_seq')
Mouse_miRBase_fasta$ID = foreach(i = 1:nrow(Mouse_miRBase_fasta),.combine = 'c') %do% {
  strsplit(Mouse_miRBase_fasta$ID[i],' ')[[1]][1]
}

Mouse_GtRNAdb_fasta = read.table('mouse/mm39-tRNAs.fa.tsv',stringsAsFactors = F,sep = '\t')
colnames(Mouse_GtRNAdb_fasta) = c('ID','Fasta_seq')
Mouse_GtRNAdb_fasta_annotation = foreach(i = 1:nrow(Mouse_GtRNAdb_fasta),.combine = 'c') %do% {
  strsplit(Mouse_GtRNAdb_fasta$ID[i],' ')[[1]][11]
}
Mouse_GtRNAdb_fasta$ID = foreach(i = 1:nrow(Mouse_GtRNAdb_fasta),.combine = 'c') %do% {
  str_remove(strsplit(Mouse_GtRNAdb_fasta$ID[i],' ')[[1]][1],'Mus_musculus_')
}

Mouse_piRNAdb_fasta = read.table('mouse/piRNAdb.mmu.v1_7_6.fa.tsv',stringsAsFactors = F)
colnames(Mouse_piRNAdb_fasta) = c('ID','Fasta_seq')

Mouse_rRNA_ucsc_fasta = read.table('mouse/mouse_rRNA_ucsc.fasta.tsv',stringsAsFactors = F)
Mouse_rRNA_ucsc_fasta_annotation = str_remove_all(Mouse_rRNA_ucsc_fasta$V2,'range=')
Mouse_rRNA_ucsc_fasta = Mouse_rRNA_ucsc_fasta[,c(1,7)]
colnames(Mouse_rRNA_ucsc_fasta) = c('ID','Fasta_seq')

Mouse_miRBase_sam_stats = rbind.data.frame(
  read.table('mouse/miRBase_mature_sam_stats.txt',stringsAsFactors = F),
  read.table('mouse/miRBase_hairpin_sam_stats.txt',stringsAsFactors = F))
colnames(Mouse_miRBase_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Mouse_miRBase_sam_stats$Fasta_start = Mouse_miRBase_sam_stats$Fasta_start - 1
Mouse_miRBase_sam_stats$Fasta_mapped = foreach(i = Mouse_miRBase_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

Mouse_GtRNAdb_sam_stats = read.table('mouse/GtRNAdb_sam_stats.txt',stringsAsFactors = F)
colnames(Mouse_GtRNAdb_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Mouse_GtRNAdb_sam_stats$Fasta_start = Mouse_GtRNAdb_sam_stats$Fasta_start - 1
Mouse_GtRNAdb_sam_stats$ID = foreach(i = Mouse_GtRNAdb_sam_stats$ID,.combine = 'c') %do% {
  strsplit(i,'_')[[1]][3]
}
Mouse_GtRNAdb_sam_stats$Fasta_mapped = foreach(i = Mouse_GtRNAdb_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]}

Mouse_piRNAdb_sam_stats = read.table('mouse/piRNAdb_sam_stats.txt',stringsAsFactors = F)
colnames(Mouse_piRNAdb_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                           'Fasta_len','Fasta_mapped')
Mouse_piRNAdb_sam_stats$Fasta_start = Mouse_piRNAdb_sam_stats$Fasta_start - 1
Mouse_piRNAdb_sam_stats$Fasta_mapped = foreach(i = Mouse_piRNAdb_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

Mouse_rRNA_sam_stats = read.table('mouse/rRNA_ucsc_sam_stats.txt',stringsAsFactors = F)
colnames(Mouse_rRNA_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Mouse_rRNA_sam_stats$Fasta_start = Mouse_rRNA_sam_stats$Fasta_start - 1
Mouse_rRNA_sam_stats$ID = str_remove_all(Mouse_rRNA_sam_stats$ID,'rn7_rmsk_')
Mouse_rRNA_sam_stats$Fasta_mapped = foreach(i = Mouse_rRNA_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]}

Mouse_tRFdb_sam_stats = read.table('mouse/tRFdb_sam_stats.txt',stringsAsFactors = F)
colnames(Mouse_rRNA_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Mouse_rRNA_sam_stats$Fasta_start = Mouse_rRNA_sam_stats$Fasta_start - 1
Mouse_rRNA_sam_stats$ID = str_remove_all(Mouse_rRNA_sam_stats$ID,'rn7_rmsk_')
Mouse_rRNA_sam_stats$Fasta_mapped = foreach(i = Mouse_rRNA_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]}

```

## Add fasta sequences to loci -> check equal seq

```{r}
mouse_databases$miRBase$Fasta_seq = NA
foreach(i = 1:nrow(mouse_databases$miRBase),.combine = 'c') %do% {
  mouse_databases$miRBase$Fasta_seq[i] = 
    Mouse_miRBase_fasta[which(Mouse_miRBase_fasta$ID == mouse_databases$miRBase$ID[i])[1],2]
}
  
mouse_databases$GtRNAdb$Fasta_seq = NA
foreach(i = 1:nrow(mouse_databases$GtRNAdb),.combine = 'c') %do% {
  mouse_databases$GtRNAdb$Fasta_seq[i] = 
    Mouse_GtRNAdb_fasta[which(Mouse_GtRNAdb_fasta$ID == mouse_databases$GtRNAdb$ID[i])[1],2]
}

mouse_databases$piRNAdb$Fasta_seq = NA
mouse_databases$piRNAdb$Fasta_seq = foreach(i = 1:nrow(mouse_databases$piRNAdb),.combine = 'c') %do% {
  if (i %% 100000 == 0) {print(i)}
  Mouse_piRNAdb_fasta[which(Mouse_piRNAdb_fasta$ID == mouse_databases$piRNAdb$ID[i])[1],2]
}

mouse_databases$rRNA_ucsc$Fasta_seq = foreach(i = 1:nrow(mouse_databases$rRNA_ucsc),.combine = 'c') %do% {
  Mouse_rRNA_ucsc_fasta[which(
    Mouse_rRNA_ucsc_fasta_annotation == 
      str_c(mouse_databases$rRNA_ucsc$Loci_chr[i],':',
            mouse_databases$rRNA_ucsc$Loci_start[i]+1,'-',
            mouse_databases$rRNA_ucsc$Loci_end[i]))[1],2]
}

mouse_databases$tRFdb$Fasta_seq = Mouse_trfdb_raw_table$tRF.Sequence

mouse_databases$miRBase$Fasta_len = foreach(i = strsplit(mouse_databases$miRBase$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
mouse_databases$GtRNAdb$Fasta_len = foreach(i = strsplit(mouse_databases$GtRNAdb$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
mouse_databases$piRNAdb$Fasta_len = foreach(i = strsplit(mouse_databases$piRNAdb$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
mouse_databases$rRNA_ucsc$Fasta_len = foreach(i = strsplit(mouse_databases$rRNA_ucsc$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
mouse_databases$tRFdb$Fasta_seq = str_remove_all(mouse_databases$tRFdb$Fasta_seq,' ')
mouse_databases$tRFdb$Fasta_len = foreach(i = strsplit(mouse_databases$tRFdb$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}

mouse_databases$miRBase$equal_seq = 
  mouse_databases$miRBase$Fasta_seq == mouse_databases$miRBase$Loci_seq
mouse_databases$GtRNAdb$equal_seq = 
  mouse_databases$GtRNAdb$Fasta_seq == mouse_databases$GtRNAdb$Loci_seq
mouse_databases$piRNAdb$equal_seq = 
  mouse_databases$piRNAdb$Fasta_seq == mouse_databases$piRNAdb$Loci_seq
mouse_databases$rRNA_ucsc$equal_seq = 
  mouse_databases$rRNA_ucsc$Fasta_seq == mouse_databases$rRNA_ucsc$Loci_seq
mouse_databases$tRFdb$equal_seq = 
  mouse_databases$tRFdb$Fasta_seq == mouse_databases$tRFdb$Loci_seq

mouse_databases$miRBase$equal_rev_seq = 
  mouse_databases$miRBase$Fasta_seq == mouse_databases$miRBase$Loci_rev_seq
mouse_databases$GtRNAdb$equal_rev_seq = 
  mouse_databases$GtRNAdb$Fasta_seq == mouse_databases$GtRNAdb$Loci_rev_seq
mouse_databases$piRNAdb$equal_rev_seq = 
  mouse_databases$piRNAdb$Fasta_seq == mouse_databases$piRNAdb$Loci_rev_seq
mouse_databases$rRNA_ucsc$equal_rev_seq = 
  mouse_databases$rRNA_ucsc$Fasta_seq == mouse_databases$rRNA_ucsc$Loci_rev_seq
mouse_databases$tRFdb$equal_rev_seq = 
  mouse_databases$tRFdb$Fasta_seq == mouse_databases$tRFdb$Loci_rev_seq

```

## Stats of united database

```{r}

write.xlsx(data.frame(
  Unique_ID_amount = c(
length(unique(mouse_databases$miRBase$ID)),
length(unique(mouse_databases$piRNAdb$ID)),
length(unique(mouse_databases$GtRNAdb$ID)),
length(unique(mouse_databases$rRNA_ucsc$ID)),
length(unique(mouse_databases$tRFdb$ID))
),
  Fasta_seqs_amount = c(
sum(!is.na(mouse_databases$miRBase$Fasta_seq)),
sum(!is.na(mouse_databases$piRNAdb$Fasta_seq)),
sum(!is.na(mouse_databases$GtRNAdb$Fasta_seq)),
sum(!is.na(mouse_databases$rRNA_ucsc$Fasta_seq)),
sum(!is.na(mouse_databases$tRFdb$Fasta_seq))
),
  Equal_Fasta_Loci_seqs = c(
sum(mouse_databases$miRBase$equal_seq,na.rm = T) + sum(mouse_databases$miRBase$equal_rev_seq,na.rm = T),
sum(mouse_databases$piRNAdb$equal_seq,na.rm = T) + sum(mouse_databases$piRNAdb$equal_rev_seq,na.rm = T),
sum(mouse_databases$GtRNAdb$equal_seq,na.rm = T) + sum(mouse_databases$GtRNAdb$equal_rev_seq,na.rm = T),
sum(mouse_databases$rRNA_ucsc$equal_seq,na.rm = T) + sum(mouse_databases$rRNA_ucsc$equal_rev_seq,na.rm = T),
sum(mouse_databases$tRFdb$equal_seq,na.rm = T) + sum(mouse_databases$tRFdb$equal_rev_seq,na.rm = T)
), row.names = names(mouse_databases)),
'mouse/smRNA_compamouseive_table_stats.xlsx',rowNames = T)

write.xlsx(mouse_databases,'mouse/smRNA_compamouseive_table_mouse.xlsx',rowNames = T)

```

## Pictures

```{r} 

venn.diagram(list(unique(mouse_databases$miRBase$ID), unique(Mouse_miRBase_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'mouse/Venn_mouse_miRBase_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red','blue'), cat.dist = c(0.05,0.065),
             main = 'Unique ID of miRBase', fill = c('red','blue'))

venn.diagram(list(unique(mouse_databases$piRNAdb$ID), unique(Mouse_piRNAdb_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'mouse/Venn_mouse_piRNAdb_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red','blue'), cat.dist = c(0.05,0.04),
             main = 'Unique ID of piRNAdb', fill = c('red','blue'))

venn.diagram(list(unique(mouse_databases$GtRNAdb$ID), unique(Mouse_GtRNAdb_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'mouse/Venn_mouse_GtRNAdb_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red','blue'), cat.dist = c(0.05,0.04),
             main = 'Unique ID of GtRNAdb', fill = c('red','blue'))

ggplot(mouse_databases$miRBase) + aes(Loci_len,Fasta_len) +
  theme_classic() + geom_point() + geom_abline() +
  ggtitle('Comparison of length of miRBase loci and fasta seqs')
ggsave('mouse/Comparison of mouse length of microRNA.tiff')

ggplot(mouse_databases$tRFdb) + aes(Loci_len,Fasta_len) +
  theme_classic() + geom_point() + geom_abline() +
  ggtitle('Comparison of length of tRFdb loci and fasta seqs')
ggsave('mouse/Comparison of mouse length of tRNA-derived.tiff')



```


## Building united database

```{r loci&fa} 
Mouse_miRBase_precursors_loci_fa = filter(mouse_databases$miRBase,equal_seq+equal_rev_seq == 1,
                                    str_detect(ID,'mmu-mir')) #is null
Mouse_miRBase_mature_loci_fa = filter(mouse_databases$miRBase,equal_seq+equal_rev_seq == 1,
                                    str_detect(ID,'mmu-miR')) #is null
Mouse_piRNAdb_loci_fa = filter(mouse_databases$piRNAdb,equal_seq+equal_rev_seq == 1)
Mouse_GtRNAdb_loci_fa = filter(mouse_databases$GtRNAdb,equal_seq+equal_rev_seq == 1)
Mouse_rRNA_ucsc_loci_fa = filter(mouse_databases$rRNA_ucsc,equal_seq+equal_rev_seq == 1)
Mouse_tRFdb_loci_fa = filter(mouse_databases$tRFdb,equal_seq+equal_rev_seq == 1)

mouse_smRNA_DATABASE = rbind.data.frame(
data.frame(type = 'tRNA-derived',
           ID = Mouse_tRFdb_loci_fa$ID,
           Annotation_consensus = Mouse_tRFdb_loci_fa$Loci,
           Annotation_chr = Mouse_tRFdb_loci_fa$Loci_chr,
           Annotation_start = Mouse_tRFdb_loci_fa$Loci_start,
           Annotation_end = Mouse_tRFdb_loci_fa$Loci_end,
           Sequence_consensus = Mouse_tRFdb_loci_fa$Fasta_seq,
           Length_consensus = Mouse_tRFdb_loci_fa$Loci_len,
           strand = Mouse_tRFdb_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'tRNA',
           ID = Mouse_GtRNAdb_loci_fa$ID,
           Annotation_consensus = Mouse_GtRNAdb_loci_fa$Loci,
           Annotation_chr = Mouse_GtRNAdb_loci_fa$Loci_chr,
           Annotation_start = Mouse_GtRNAdb_loci_fa$Loci_start,
           Annotation_end = Mouse_GtRNAdb_loci_fa$Loci_end,
           Sequence_consensus = Mouse_GtRNAdb_loci_fa$Fasta_seq,
           Length_consensus = Mouse_GtRNAdb_loci_fa$Loci_len,
           strand = Mouse_GtRNAdb_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'rRNA',
           ID = Mouse_rRNA_ucsc_loci_fa$ID,
           Annotation_consensus = Mouse_rRNA_ucsc_loci_fa$Loci,
           Annotation_chr = Mouse_rRNA_ucsc_loci_fa$Loci_chr,
           Annotation_start = Mouse_rRNA_ucsc_loci_fa$Loci_start,
           Annotation_end = Mouse_rRNA_ucsc_loci_fa$Loci_end,
           Sequence_consensus = Mouse_rRNA_ucsc_loci_fa$Fasta_seq,
           Length_consensus = Mouse_rRNA_ucsc_loci_fa$Loci_len,
           strand = Mouse_rRNA_ucsc_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
)
)

```

```{r loci only}

Mouse_miRBase_precursors_loci_only = rbind.data.frame(
  filter(mouse_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'mmu-mir')),
  filter(mouse_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'mmu-let')))
Mouse_miRBase_precursors_loci_only = filter(Mouse_miRBase_precursors_loci_only,
                                            !(str_detect(ID,'p')))
Mouse_miRBase_precursors_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(Mouse_miRBase_precursors_loci_only),.combine = 'c') %do% {
    if (Mouse_miRBase_precursors_loci_only$Loci_strand[i] == '-') {out = Mouse_miRBase_precursors_loci_only$Loci_seq[i]}
    if (Mouse_miRBase_precursors_loci_only$Loci_strand[i] == '+') {out = Mouse_miRBase_precursors_loci_only$Loci_rev_seq[i]}
    out
  }

Mouse_miRBase_mature_loci_only = unique(rbind.data.frame(
  filter(mouse_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'mmu-miR')),
  filter(mouse_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'p'))))
Mouse_miRBase_mature_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(Mouse_miRBase_mature_loci_only),.combine = 'c') %do% {
    if (Mouse_miRBase_mature_loci_only$Loci_strand[i] == '-') {out = Mouse_miRBase_mature_loci_only$Loci_seq[i]}
    if (Mouse_miRBase_mature_loci_only$Loci_strand[i] == '+') {out = Mouse_miRBase_mature_loci_only$Loci_rev_seq[i]}
    out
  }
Mouse_piRNAdb_loci_only = filter(mouse_databases$piRNAdb,is.na(Fasta_seq)) # is null

Mouse_GtRNAdb_loci_only = filter(mouse_databases$GtRNAdb,is.na(Fasta_seq))
Mouse_GtRNAdb_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(Mouse_GtRNAdb_loci_only),.combine = 'c') %do% {
    if (Mouse_GtRNAdb_loci_only$Loci_strand[i] == '-') {out = Mouse_GtRNAdb_loci_only$Loci_seq[i]}
    if (Mouse_GtRNAdb_loci_only$Loci_strand[i] == '+') {out = Mouse_GtRNAdb_loci_only$Loci_rev_seq[i]}
    out
  }
Mouse_tRFdb_loci_only = filter(mouse_databases$tRFdb,is.na(Fasta_seq)) # is null
Mouse_rRNA_ucsc_loci_only = filter(mouse_databases$rRNA_ucsc,is.na(Fasta_seq)) # is null

mouse_smRNA_DATABASE = rbind.data.frame(mouse_smRNA_DATABASE,
  data.frame(type = 'precursor_miRNA',
             ID = Mouse_miRBase_precursors_loci_only$ID,
             Annotation_consensus = Mouse_miRBase_precursors_loci_only$Loci,
             Annotation_chr = Mouse_miRBase_precursors_loci_only$Loci_chr,
             Annotation_start = Mouse_miRBase_precursors_loci_only$Loci_start,
             Annotation_end = Mouse_miRBase_precursors_loci_only$Loci_end,
             Sequence_consensus = Mouse_miRBase_precursors_loci_only$Fasta_seq,
             Length_consensus = Mouse_miRBase_precursors_loci_only$Loci_len,
             strand = Mouse_miRBase_precursors_loci_only$Loci_strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ),
  data.frame(type = 'mature_miRNA',
             ID = Mouse_miRBase_mature_loci_only$ID,
             Annotation_consensus = Mouse_miRBase_mature_loci_only$Loci,
             Annotation_chr = Mouse_miRBase_mature_loci_only$Loci_chr,
             Annotation_start = Mouse_miRBase_mature_loci_only$Loci_start,
             Annotation_end = Mouse_miRBase_mature_loci_only$Loci_end,
             Sequence_consensus = Mouse_miRBase_mature_loci_only$Fasta_seq,
             Length_consensus = Mouse_miRBase_mature_loci_only$Loci_len,
             strand = Mouse_miRBase_mature_loci_only$Loci_strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ),
  data.frame(type = 'tRNA',
             ID = Mouse_GtRNAdb_loci_only$ID,
             Annotation_consensus = Mouse_GtRNAdb_loci_only$Loci,
             Annotation_chr = Mouse_GtRNAdb_loci_only$Loci_chr,
             Annotation_start = Mouse_GtRNAdb_loci_only$Loci_start,
             Annotation_end = Mouse_GtRNAdb_loci_only$Loci_end,
             Sequence_consensus = Mouse_GtRNAdb_loci_only$Fasta_seq,
             Length_consensus = Mouse_GtRNAdb_loci_only$Loci_len,
             strand = Mouse_GtRNAdb_loci_only$Loci_strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ))

```

```{r fa only}

Mouse_miRBase_precursors_fa_only = filter(Mouse_miRBase_sam_stats, !(ID %in% mouse_databases$miRBase$ID), str_detect(ID,'mmu-mir'), Fasta_chr != '*')
Mouse_miRBase_mature_fa_only = filter(Mouse_miRBase_sam_stats, !(ID %in% mouse_databases$miRBase$ID),Fasta_chr != '*', str_detect(ID,'rno-miR')) # is null
Mouse_piRNAdb_fa_only = filter(Mouse_piRNAdb_sam_stats,!(ID %in% mouse_databases$piRNAdb$ID),
                         Fasta_chr != '*')
Mouse_GtRNAdb_fa_only = filter(Mouse_GtRNAdb_sam_stats,!(ID %in% mouse_databases$GtRNAdb$ID),
                         Fasta_chr != '*') # is null
Mouse_rRNA_ucsc_fa_only = filter(Mouse_rRNA_sam_stats,!(ID %in% mouse_databases$rRNA_ucsc$ID),
                         Fasta_chr != '*') # is null

mouse_smRNA_DATABASE = rbind.data.frame(mouse_smRNA_DATABASE,
                                  data.frame(type = 'precursor_miRNA',
                                             ID = Mouse_miRBase_precursors_fa_only$ID,
        Annotation_consensus = str_c(Mouse_miRBase_precursors_fa_only$Fasta_chr,':',
                                     Mouse_miRBase_precursors_fa_only$Fasta_start,'-',
                                     as.character(Mouse_miRBase_precursors_fa_only$Fasta_start+Mouse_miRBase_precursors_fa_only$Fasta_len)),
                                             Annotation_chr = Mouse_miRBase_precursors_fa_only$Fasta_chr,
                                             Annotation_start = Mouse_miRBase_precursors_fa_only$Fasta_start,
                                             Annotation_end = Mouse_miRBase_precursors_fa_only$Fasta_start+
          Mouse_miRBase_precursors_fa_only$Fasta_len,
                                             Sequence_consensus = Mouse_miRBase_precursors_fa_only$Fasta_seq,
                                             Length_consensus = Mouse_miRBase_precursors_fa_only$Fasta_len,
        strand = '+',
        strand_from_annotation = F,
                                             Source = 'fasta',
                                             stringsAsFactors = F
                                  ),
                                  data.frame(type = 'piRNA',
                                             ID = Mouse_piRNAdb_fa_only$ID,
                                             Annotation_consensus = str_c(Mouse_piRNAdb_fa_only$Fasta_chr,':',Mouse_piRNAdb_fa_only$Fasta_start,'-',
                                                                          as.character(Mouse_piRNAdb_fa_only$Fasta_start+Mouse_piRNAdb_fa_only$Fasta_len)),
                                             Annotation_chr = Mouse_piRNAdb_fa_only$Fasta_chr,
                                             Annotation_start = Mouse_piRNAdb_fa_only$Fasta_start,
                                             Annotation_end = Mouse_piRNAdb_fa_only$Fasta_start+Mouse_piRNAdb_fa_only$Fasta_len,
                                             Sequence_consensus = Mouse_piRNAdb_fa_only$Fasta_seq,
                                             Length_consensus = Mouse_piRNAdb_fa_only$Fasta_len,
                                             strand = '+',
                                             strand_from_annotation = F,
                                             Source = 'fasta',
                                             stringsAsFactors = F
                                  ))

```

```{r count H-dist & suff-pref intersect}

Mouse_miRBase_precursors_loci_not_fa = rbind.data.frame(
  filter(mouse_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'mmu-mir')),
  filter(mouse_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'mmu-let')))
Mouse_miRBase_precursors_loci_not_fa = filter(Mouse_miRBase_precursors_loci_not_fa,
                                              !(str_detect(ID,'p')))
Mouse_miRBase_mature_loci_not_fa = unique(rbind.data.frame(
  filter(mouse_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'mmu-miR')),
  filter(mouse_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'p'))))
Mouse_GtRNAdb_loci_not_fa = filter(mouse_databases$GtRNAdb, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq)) # is null
Mouse_rRNA_ucsc_loci_not_fa = filter(mouse_databases$rRNA_ucsc, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq)) # is null
Mouse_piRNAdb_loci_not_fa = filter(mouse_databases$piRNAdb, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq))
Mouse_tRFdb_loci_not_fa = filter(mouse_databases$tRFdb, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq))

hamming_dist_unequal = function(x,y) {
  library(foreach)
  x = strsplit(x,'')[[1]]
  y = strsplit(y,'')[[1]]
  baseseq = list(x,y)[[which.max(c(length(x),length(y)))]]
  moveseq = list(x,y)[[-which.max(c(length(x),length(y)))]]
  return(min(foreach(i = 1:(length(baseseq)-length(moveseq)+1),.combine = 'c') %dopar% {
    sum(baseseq[i:(i+length(moveseq)-1)] != moveseq)
  }))
}
intersequences = function(x,y,mismatch = 0) {
  x = strsplit(x,'')[[1]]
  y = strsplit(y,'')[[1]]
  baseseq = list(x,y)[[which.min(c(length(x),length(y)))]]
  moveseq = list(x,y)[[-which.min(c(length(x),length(y)))]]
  tab = foreach(i = 1:length(baseseq),.combine = 'rbind.data.frame') %dopar% {
    pref = 0
    suff = 0
    if (sum(baseseq[1:i] == moveseq[(length(moveseq)+1-i):length(moveseq)]) >= i-mismatch) 
    {pref = i}
    if (sum(moveseq[1:i] == baseseq[(length(baseseq)+1-i):length(baseseq)]) >= i-mismatch) 
    {suff = i}
    c(max(pref,suff),c(NA,0,1)[which.max(c(0.1,pref,suff))])
  }
  out = tab[which.max(tab[,1]),]
  colnames(out) = c('inter_len','inter_seq')
  if (is.na(out[,2])) {return(out)}
  else if (out[,2] == '1') {out[,2] = str_c(moveseq[1:out[1,1]], collapse = '')
  return(out)}
  else if (out[,2] == '0') {out[,2] = str_c(baseseq[1:out[1,1]], collapse = '')
  return(out)}
  
}

Mouse_loci_not_fa = foreach(d = list(Mouse_miRBase_precursors_loci_not_fa,
Mouse_miRBase_mature_loci_not_fa,
Mouse_tRFdb_loci_not_fa)) %do% {
  print(nrow(d))
d = cbind.data.frame(
  d,
  foreach(i = 1:nrow(d),.combine = 'rbind',.packages = libs) %do% {
    if (i %% 100 == 0) {print(i)}
  fseq = d$Fasta_seq[i]
  lfseq = d$Loci_seq[i]
  lrseq = d$Loci_rev_seq[i]
  forw_split = strsplit(fseq,lfseq, fixed = T)[[1]]
  rev_split = strsplit(fseq,lrseq, fixed = T)[[1]]
  if (forw_split[1] != fseq) {
    if (forw_split[1] == ''){forw_split[1] = NA}
      o = data.frame(p1 = 'loci_forward_in_fa',
      p2 = length(na.omit(forw_split[1])),p3 = length(na.omit(forw_split[2])),stringsAsFactors = F)
  }
  else if (rev_split[1] != fseq) {
    if (rev_split[1] == ''){rev_split[1] = NA}
      o = data.frame(p1 = 'loci_reverse_in_fa',
      p2 = length(na.omit(rev_split[1])),p3 = length(na.omit(rev_split[2])),stringsAsFactors = F)
  }
  else if (strsplit(lfseq,fseq, fixed = T)[[1]] != lfseq) {
    o = data.frame(p1 = 'fa_in_loci',p2 = 'forward',p3 = NA,stringsAsFactors = F)
  }
  else if (strsplit(lrseq,fseq, fixed = T)[[1]] != lrseq) {
    o = data.frame(p1 = 'fa_in_loci',p2 = 'reverse',p3 = NA,stringsAsFactors = F)
  }
  else {
    inter = rbind(intersequences(d$Loci_seq[i],d$Fasta_seq[i]),
    intersequences(d$Loci_rev_seq[i],d$Fasta_seq[i]))
    colnames(inter) = c('p2','p3')
    o = data.frame(
      p1 = min(hamming_dist_unequal(d$Loci_seq[i],d$Fasta_seq[i]),
          hamming_dist_unequal(d$Loci_rev_seq[i],d$Fasta_seq[i])),
      inter[which.max(inter[,1]),],stringsAsFactors = F
    )
  }
  o
  }
)
d
}
names(Mouse_loci_not_fa) = c('miRBase_precursors','miRBase_mature','tRFdb')

Mouse_loci_not_fa$piRNAdb = Mouse_piRNAdb_loci_not_fa
d = Mouse_loci_not_fa$piRNAdb
d = unique(d[,8:10])

d = cbind.data.frame(d,
foreach(i = 1:nrow(d),.combine = 'rbind',.packages = libs) %do% {
    if (i %% 1000 == 0) {print(i)}
  fseq = d$Fasta_seq[i]
  lfseq = d$Loci_seq[i]
  lrseq = d$Loci_rev_seq[i]
  forw_split = strsplit(fseq,lfseq, fixed = T)[[1]]
  rev_split = strsplit(fseq,lrseq, fixed = T)[[1]]
  if (forw_split[1] != fseq) {
      o = data.frame(p1 = 'loci_forward_in_fa',
      p2 = length(forw_split[1]),p3 = length(na.omit(forw_split[2])),stringsAsFactors = F)
  }
  else if (rev_split[1] != fseq) {
      o = data.frame(p1 = 'loci_reverse_in_fa',
      p2 = length(rev_split[1]),p3 = length(na.omit(rev_split[2])),stringsAsFactors = F)
  }
  else if (strsplit(lfseq,fseq, fixed = T)[[1]][1] != lfseq) {
    o = data.frame(p1 = 'fa_in_loci',p2 = 'forward',p3 = NA,stringsAsFactors = F)
  }
  else if (strsplit(lrseq,fseq, fixed = T)[[1]][1] != lrseq) {
    o = data.frame(p1 = 'fa_in_loci',p2 = 'reverse',p3 = NA,stringsAsFactors = F)
  }
  else {
    inter = rbind(intersequences(d$Loci_seq[i],d$Fasta_seq[i]),
    intersequences(d$Loci_rev_seq[i],d$Fasta_seq[i]))
    colnames(inter) = c('p2','p3')
    o = data.frame(
      p1 = min(hamming_dist_unequal(d$Loci_seq[i],d$Fasta_seq[i]),
          hamming_dist_unequal(d$Loci_rev_seq[i],d$Fasta_seq[i])),
      inter[which.max(inter[,1]),],stringsAsFactors = F
    )
  }
  o
}
)
d$ID = str_c(d$Loci_seq,'_',d$Fasta_seq)
length(d$ID) == length(unique(d$ID))
Mouse_loci_not_fa$piRNAdb$ID = str_c(Mouse_loci_not_fa$piRNAdb$Loci_seq,'_',Mouse_loci_not_fa$piRNAdb$Fasta_seq)
Mouse_loci_not_fa$piRNAdb = 
  cbind.data.frame(Mouse_loci_not_fa$piRNAdb,
                   foreach(i = 1:nrow(Mouse_loci_not_fa$piRNAdb),.combine = 'rbind') %do% {
                     if (i %% 10000 == 0) {print(i)}
                     k = Mouse_loci_not_fa$piRNAdb$ID[i]
                     d[which(d$ID == k),4:6]
                   }
                   )

```

```{r fix loci not fa}
Mouse_loci_not_fa$miRBase_precursors$type = 'precursor_miRNA'
Mouse_loci_not_fa$miRBase_mature$type = 'mature_miRNA'
Mouse_loci_not_fa$tRFdb$type = 'tRNA-derived'
Mouse_loci_not_fa$piRNAdb$type = 'piRNA'

Mouse_loci_not_fa$piRNAdb$ID = Mouse_piRNAdb_loci_not_fa$ID


for (i in 1:3) {
  d = Mouse_loci_not_fa[[i]]
  mouse_smRNA_DATABASE = rbind.data.frame(mouse_smRNA_DATABASE, 
    foreach(k = 1:nrow(d),.combine = 'rbind',.packages = libs) %dopar% {
    if (k %% 10000 == 0) {print(k)}
      Annotation_start = 0
      Annotation_end = 0
      Sequence_consensus = '0'
      Length_consensus = 0
      Source = 'trash'
    if (d$p1[k] %in% c('loci_forward_in_fa')) {
      Annotation_start = d$Loci_start[k]-as.numeric(d$p2[k])
      Annotation_end = d$Loci_end[k]+as.numeric(d$p3[k])
      Sequence_consensus = d$Fasta_seq[k]
      Length_consensus = d$Fasta_len[k]
      Source = 'wide_annotation&fasta'
    }
    if (d$p1[k] %in% c('loci_reverse_in_fa')) {
      Annotation_start = d$Loci_start[k]-as.numeric(d$p3[k])
      Annotation_end = d$Loci_end[k]+as.numeric(d$p2[k])
      Sequence_consensus = d$Fasta_seq[k]
      Length_consensus = d$Fasta_len[k]
      Source = 'wide_annotation&fasta'
    }
    if (d$p1[k] == 'fa_in_loci') {
      Annotation_start = d$Loci_start[k]
      Annotation_end = d$Loci_end[k]
      if (d$p2[k] == 'forward') {Sequence_consensus = d$Loci_seq[k]}
      else if (d$p2[k] == 'reverse') {Sequence_consensus = d$Loci_rev_seq[k]}
      
      Length_consensus = d$Fasta_len[k]
      Source = 'annotation&wide_fasta'
    }
    else if (d$Fasta_len[k] - as.numeric(d$p2[k]) == 1) {
      Length_consensus = d$Loci_len[k]+1
      Source = 'wide_annotation&wide_fasta'
      ssplit = strsplit(d$Fasta_seq[k],d$p3[k])[[1]]
      if (strsplit(d$Loci_seq[k],d$p3[k])[[1]][1] != d$Loci_seq[k]) {
        Sequence_consensus = d$Loci_seq[k]
        if (length(ssplit) == 1) {
          Sequence_consensus = str_c(ssplit, Sequence_consensus)
          Annotation_start = d$Loci_start[k]-1
          Annotation_end = d$Loci_end[k]
          
      }
      else if (length(ssplit) == 2 && ssplit[1] == '') {
          Sequence_consensus = str_c(Sequence_consensus, ssplit[2])
          Annotation_start = d$Loci_start[k]
          Annotation_end = d$Loci_end[k]+1
        }
      }
      else {
        Sequence_consensus = d$Loci_rev_seq[k]
        if (length(ssplit) == 1) {
          Sequence_consensus = str_c(ssplit, Sequence_consensus)
          Annotation_start = d$Loci_start[k]
          Annotation_end = d$Loci_end[k]+1
      }
      else if (length(ssplit) == 2 && ssplit[1] == '') {
          Sequence_consensus = str_c(Sequence_consensus, ssplit[2])
          Annotation_start = d$Loci_start[k]-1
          Annotation_end = d$Loci_end[k]
        }
      }
      
      
    }
      
    data.frame(type = d$type[k],ID = d$ID[k],
               Annotation_consensus = str_c(d$Loci_chr[k],':',Annotation_start,'-', Annotation_end),
               Annotation_chr = d$Loci_chr[k], Annotation_start = Annotation_start,
               Annotation_end = Annotation_end, Sequence_consensus = Sequence_consensus,
               Length_consensus = Length_consensus,
               strand = d$Loci_strand[k], strand_from_annotation = T,
               Source = Source,stringsAsFactors = F
               )
  })
}
mouse_smRNA_DATABASE = filter(mouse_smRNA_DATABASE,Source != 'trash')

d = Mouse_loci_not_fa[[4]]
d$diff_len = d$Fasta_len - d$Loci_len
d$fix_p2 = d$Fasta_len - d$Loci_len == as.numeric(d$p3)
d$p2[which(d$fix_p2 == T)] = '0'

wide_annot = filter(d,p1 %in% c('loci_forward_in_fa'))
Annotation_start = wide_annot$Loci_start-as.numeric(wide_annot$p2)
Annotation_end = wide_annot$Loci_end+as.numeric(wide_annot$p3)
Sequence_consensus = wide_annot$Fasta_seq
Length_consensus = wide_annot$Fasta_len
Source = 'wide_annotation&fasta'
mouse_smRNA_DATABASE = rbind.data.frame(mouse_smRNA_DATABASE, 
      data.frame(type = wide_annot$type,ID = wide_annot$ID,
               Annotation_consensus = str_c(wide_annot$Loci_chr,':',Annotation_start,'-', Annotation_end),
               Annotation_chr = wide_annot$Loci_chr, Annotation_start = Annotation_start,
               Annotation_end = Annotation_end, Sequence_consensus = Sequence_consensus,
               Length_consensus = Length_consensus, 
               strand = wide_annot$Loci_strand, strand_from_annotation = T,
               Source = Source,stringsAsFactors = F
               ))
wide_annot = filter(d,p1 %in% c('loci_reverse_in_fa'))
Annotation_start = wide_annot$Loci_start-as.numeric(wide_annot$p3)
Annotation_end = wide_annot$Loci_end+as.numeric(wide_annot$p2)
Sequence_consensus = wide_annot$Fasta_seq
Length_consensus = wide_annot$Fasta_len
Source = 'wide_annotation&fasta'
mouse_smRNA_DATABASE = rbind.data.frame(mouse_smRNA_DATABASE, 
      data.frame(type = wide_annot$type,ID = wide_annot$ID,
               Annotation_consensus = str_c(wide_annot$Loci_chr,':',Annotation_start,'-', Annotation_end),
               Annotation_chr = wide_annot$Loci_chr, Annotation_start = Annotation_start,
               Annotation_end = Annotation_end, Sequence_consensus = Sequence_consensus,
               Length_consensus = Length_consensus, 
               strand = wide_annot$Loci_strand, strand_from_annotation = T,
               Source = Source,stringsAsFactors = F
               ))


wide_fasta = filter(d,p1 == 'fa_in_loci')
Annotation_start = wide_fasta$Loci_start
Annotation_end = wide_fasta$Loci_end
Sequence_consensus = foreach(k = 1:nrow(wide_fasta),.combine = 'c') %do% {
  if (wide_fasta$p2[k] == 'forward') {o = wide_fasta$Loci_seq[k]}
  else if (wide_fasta$p2[k] == 'reverse') {o = wide_fasta$Loci_rev_seq[k]}
  o
}
Length_consensus = wide_fasta$Loci_len
Source = 'annotation&wide_fasta'  
mouse_smRNA_DATABASE = rbind.data.frame(mouse_smRNA_DATABASE, 
      data.frame(type = wide_fasta$type,ID = wide_fasta$ID,
               Annotation_consensus = str_c(wide_fasta$Loci_chr,':',Annotation_start,'-', Annotation_end),
               Annotation_chr = wide_fasta$Loci_chr, Annotation_start = Annotation_start,
               Annotation_end = Annotation_end, Sequence_consensus = Sequence_consensus,
               Length_consensus = Length_consensus, 
               strand = wide_fasta$Loci_strand, strand_from_annotation = T,
               Source = Source,stringsAsFactors = F
               ))

wide_annot_fasta = filter(d,Fasta_len - as.numeric(p2) == 1, diff_len == 0)
d = wide_annot_fasta
mouse_smRNA_DATABASE = rbind.data.frame(mouse_smRNA_DATABASE, 
    foreach(k = 1:nrow(d),.combine = 'rbind',.packages = libs) %dopar% {
      Annotation_start = 0
      Annotation_end = 0
      Sequence_consensus = '0'
      Length_consensus = 0
      Source = 'trash'
    if (d$Fasta_len[k] - as.numeric(d$p2[k]) == 1) {
      Length_consensus = d$Loci_len[k]+1
      Source = 'wide_annotation&wide_fasta'
      ssplit = strsplit(d$Fasta_seq[k],d$p3[k])[[1]]
      if (strsplit(d$Loci_seq[k],d$p3[k])[[1]][1] != d$Loci_seq[k]) {
        Sequence_consensus = d$Loci_seq[k]
        if (length(ssplit) == 1) {
          Sequence_consensus = str_c(ssplit, Sequence_consensus)
          Annotation_start = d$Loci_start[k]-1
          Annotation_end = d$Loci_end[k]
          
      }
      else if (length(ssplit) == 2 && ssplit[1] == '') {
          Sequence_consensus = str_c(Sequence_consensus, ssplit[2])
          Annotation_start = d$Loci_start[k]
          Annotation_end = d$Loci_end[k]+1
        }
      }
      else {
        Sequence_consensus = d$Loci_rev_seq[k]
        if (length(ssplit) == 1) {
          Sequence_consensus = str_c(ssplit, Sequence_consensus)
          Annotation_start = d$Loci_start[k]
          Annotation_end = d$Loci_end[k]+1
      }
      else if (length(ssplit) == 2 && ssplit[1] == '') {
          Sequence_consensus = str_c(Sequence_consensus, ssplit[2])
          Annotation_start = d$Loci_start[k]-1
          Annotation_end = d$Loci_end[k]
        }
      }
      
      
    }
      
    data.frame(type = d$type[k],ID = d$ID[k],
               Annotation_consensus = str_c(d$Loci_chr[k],':',Annotation_start,'-', Annotation_end),
               Annotation_chr = d$Loci_chr[k], Annotation_start = Annotation_start,
               Annotation_end = Annotation_end, Sequence_consensus = Sequence_consensus,
               Length_consensus = Length_consensus,
               strand = d$Loci_strand[k], strand_from_annotation = T,
               Source = Source,stringsAsFactors = F
               )
  })

```


## Interceptions

```{r write bed}

write.table(str_c(filter(mouse_smRNA_DATABASE,type == 'mature_miRNA')$Annotation_chr,filter(mouse_smRNA_DATABASE,type == 'mature_miRNA')$Annotation_start,filter(mouse_smRNA_DATABASE,type == 'mature_miRNA')$Annotation_end,filter(mouse_smRNA_DATABASE,type == 'mature_miRNA')$ID,sep = '\t'),'mouse/mouse_smRNA_DATABASE_miRNA.bed',col.names = F, row.names = F)
write.table(str_c(filter(mouse_smRNA_DATABASE,type == 'piRNA')$Annotation_chr,filter(mouse_smRNA_DATABASE,type == 'piRNA')$Annotation_start,filter(mouse_smRNA_DATABASE,type == 'piRNA')$Annotation_end,filter(mouse_smRNA_DATABASE,type == 'piRNA')$ID,sep = '\t'),'mouse/mouse_smRNA_DATABASE_piRNA.bed',col.names = F, row.names = F)
write.table(str_c(filter(mouse_smRNA_DATABASE,type == 'tRNA')$Annotation_chr,filter(mouse_smRNA_DATABASE,type == 'tRNA')$Annotation_start,filter(mouse_smRNA_DATABASE,type == 'tRNA')$Annotation_end,filter(mouse_smRNA_DATABASE,type == 'tRNA')$ID,sep = '\t'),'mouse/mouse_smRNA_DATABASE_tRNA.bed',col.names = F, row.names = F)
write.table(str_c(filter(mouse_smRNA_DATABASE,type == 'rRNA')$Annotation_chr,filter(mouse_smRNA_DATABASE,type == 'rRNA')$Annotation_start,filter(mouse_smRNA_DATABASE,type == 'rRNA')$Annotation_end,filter(mouse_smRNA_DATABASE,type == 'rRNA')$ID,sep = '\t'),'mouse/mouse_smRNA_DATABASE_rRNA.bed',col.names = F, row.names = F)

```

```{r filter by intersept}

mouse_inters = rbind.data.frame(
  unique(read.table('mouse/miRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('mouse/miRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('mouse/miRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
  unique(read.table('mouse/piRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('mouse/piRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('mouse/piRNA_inter_in.txt',stringsAsFactors = F,sep = '\t')[,8]),1:4]),
unique(read.table('mouse/rRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('mouse/rRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('mouse/rRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('mouse/tRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('mouse/tRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('mouse/tRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),

unique(read.table('mouse/piRNA_inter_out.txt',stringsAsFactors = F)[,1:4]),
unique(read.table('mouse/tRNA_inter_out.txt',stringsAsFactors = F)[,1:4]),
unique(read.table('mouse/rRNA_inter_out.txt',stringsAsFactors = F)[,1:4])
)
mouse_inters = unique(mouse_inters)
mouse_inters_vector = str_c(mouse_inters$V1,mouse_inters$V2,mouse_inters$V3,mouse_inters$V4)

mouse_smRNA_DATABASE$filtered = foreach(i = 1:nrow(mouse_smRNA_DATABASE),.combine = 'c',.packages = libs) %dopar% {
  str_c(mouse_smRNA_DATABASE$Annotation_chr[i],
        mouse_smRNA_DATABASE$Annotation_start[i],
        mouse_smRNA_DATABASE$Annotation_end[i],mouse_smRNA_DATABASE$ID[i]) %in% mouse_inters_vector
}

write.table(filter(mouse_smRNA_DATABASE, type %in% c('mature_miRNA', 'piRNA', 'tRNA', 'rRNA'),filtered == F),'mouse/mouse_smRNA_DATABASE_mipitr_filtered.tsv',sep = '\t', row.names = F)

write.table(str_c('>',filter(mouse_smRNA_DATABASE,type %in% c('mature_miRNA','piRNA','rRNA','tRNA'), filtered == F)$ID,'\n',filter(mouse_smRNA_DATABASE,type %in% c('mature_miRNA','piRNA','rRNA','tRNA'), filtered == F)$Sequence_consensus,collapse = '\n'),'mouse/mouse_smRNA_DATABASE_filtered.fasta',row.names = F,col.names = F)

```

```{r}
mouse_inters_in = rbind.data.frame(
  unique(read.table('mouse/miRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('mouse/miRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('mouse/miRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
  unique(read.table('mouse/piRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('mouse/piRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('mouse/piRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('mouse/rRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('mouse/rRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('mouse/rRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('mouse/tRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('mouse/tRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('mouse/tRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]))

mouse_inters_in = unique(mouse_inters_in)
mouse_inters_in_vector = str_c(mouse_inters_in$V1,mouse_inters_in$V2,mouse_inters_in$V3,mouse_inters_in$V4)


mouse_smRNA_DATABASE$filtered_in = F
mouse_smRNA_DATABASE$filtered_in[which(str_c(mouse_smRNA_DATABASE$Annotation_chr,
        mouse_smRNA_DATABASE$Annotation_start,
        mouse_smRNA_DATABASE$Annotation_end,mouse_smRNA_DATABASE$ID) %in% mouse_inters_in_vector)] = T


for (i in unique(mouse_smRNA_DATABASE$type)) {
  write.xlsx(filter(mouse_smRNA_DATABASE,filtered_in == F,type == i),str_c('mouse/mouse_smRNA_DATABASE_',i,'.xlsx'))
}
```

```{r Unique rRNA}

mouse_smRNA_DATABASE$ID[which(mouse_smRNA_DATABASE$type == 'rRNA')] = str_remove(filter(mouse_smRNA_DATABASE, type == 'rRNA')$ID,'_[[:digit:]]+')
for (i in 1:length(unique(filter(mouse_smRNA_DATABASE,type == 'rRNA')[,7]))){
  mouse_smRNA_DATABASE$ID[which(mouse_smRNA_DATABASE$Sequence_consensus == unique(filter(mouse_smRNA_DATABASE,type == 'rRNA')[,7])[i] )] = str_c(mouse_smRNA_DATABASE$ID[which(mouse_smRNA_DATABASE$Sequence_consensus == unique(filter(mouse_smRNA_DATABASE,type == 'rRNA')[,7])[i] )],'_seq_',i)
}

```

```{r write DATABASE}
Mouse_Source_df = mouse_smRNA_DATABASE %>% group_by(type) %>% 
  summarise(n = table(Source),
            Source = names(table(Source)))
write.xlsx(pivot_wider(Mouse_Source_df, names_from = type, values_from = n),'mouse/mouse_smRNA_DATABASE_stats.xlsx')
write.xlsx(mouse_smRNA_DATABASE,'mouse/mouse_smRNA_DATABASE.xlsx')
write.table(str_c('>',mouse_smRNA_DATABASE$ID[which(mouse_smRNA_DATABASE$type %in% c('mature_miRNA','piRNA','rRNA','tRNA'))],'\n',mouse_smRNA_DATABASE$Sequence_consensus[which(mouse_smRNA_DATABASE$type %in% c('mature_miRNA','piRNA','rRNA','tRNA'))],collapse = '\n'),'mouse/mouse_smRNA_DATABASE.fasta',row.names = F,col.names = F)
for (i in unique(mouse_smRNA_DATABASE$type)) {
  write.xlsx(filter(mouse_smRNA_DATABASE,type == i),str_c('mouse/mouse_smRNA_DATABASE_',i,'.xlsx'))
}

```
