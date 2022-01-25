
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

Celegans_miRBase_annotation = read.table('celegans/cel.gff3.txt',stringsAsFactors = F)
Celegans_miRBase_annotation$ID = foreach(i = 1:nrow(Celegans_miRBase_annotation),.combine = 'c') %do% {
  str_remove(strsplit(Celegans_miRBase_annotation$V9[i],';')[[1]][3],'Name=')
}
Celegans_piRNAdb_annotation = read.table('celegans/pirnadb.v1_7_6.ce11.gtf',stringsAsFactors = F)
Celegans_GtRNAdb_annotation = read.table('celegans/ce11-tRNAs.bed',stringsAsFactors = F)
Celegans_rRNA_ucsc_annotation = read.table('celegans/celegans_rRNA_ucsc.bed',stringsAsFactors = F)
Celegans_trfdb_raw_table = rbind.data.frame(
  read.xlsx('celegans/trfdb_excel_result_1.xlsx'),
  read.xlsx('celegans/trfdb_excel_result_2.xlsx')
  )
write.table(
  unique(str_c('>',Celegans_trfdb_raw_table$Type,'_',Celegans_trfdb_raw_table$tRF.ID,'\n',Celegans_trfdb_raw_table$tRF.Sequence)),'celegans/tRFdb.fasta',row.names = F,col.names = F
)

```


---
## Build bed-files for trf RNA -> liftover -> getfasta

```{r}

Celegans_trfdb_annotation = foreach(i = 1:nrow(Celegans_trfdb_raw_table),.combine = 'rbind.data.frame') %do% {
  p = unlist(c(strsplit(Celegans_trfdb_raw_table$`tRNA.Gene.Co-ordinates`[i],'-')[[1]],
  Celegans_trfdb_raw_table[i,c(1,3,8)]))
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
Celegans_trfdb_annotation$chr = str_remove(Celegans_trfdb_annotation$chr,' ')

write.table(Celegans_trfdb_annotation
,'celegans/Celegans_tRFdb_ce6.bed', col.names = F,
            row.names = F, sep = '\t')

Celegans_tRFdb_ce11_annotation = read.table('celegans/Celegans_tRFdb_ce11.bed',stringsAsFactors = F)
Celegans_piRNAdb_annotation$V4[which(Celegans_piRNAdb_annotation$V7 == '+')] = Celegans_piRNAdb_annotation$V4[which(Celegans_piRNAdb_annotation$V7 == '+')] - 1
Celegans_piRNAdb_annotation$V5[which(Celegans_piRNAdb_annotation$V7 == '-')] = Celegans_piRNAdb_annotation$V5[which(Celegans_piRNAdb_annotation$V7 == '-')] + 1
write.table(Celegans_piRNAdb_annotation[,c(1,4,5,7,10)],'celegans/Celegans_piRNAdb_ce11.bed', col.names = F,row.names = F, sep = '\t')
Celegans_miRBase_annotation$V4[which(Celegans_miRBase_annotation$V7 == '-')] = Celegans_miRBase_annotation$V4[which(Celegans_miRBase_annotation$V7 == '-')] - 1
Celegans_miRBase_annotation$V5[which(Celegans_miRBase_annotation$V7 == '+')] = Celegans_miRBase_annotation$V5[which(Celegans_miRBase_annotation$V7 == '+')] + 1
write.table(Celegans_miRBase_annotation[,c(1,4,5,7,10)],'celegans/Celegans_miRBase_ce11.bed', col.names = F,row.names = F, sep = '\t')

Celegans_miRBase_getfasta = read.table('celegans/Celegans_miRBase_ce11.getfasta.tsv',stringsAsFactors = F)
colnames(Celegans_miRBase_getfasta) = c('Loci','Loci_seq')
Celegans_GtRNAdb_getfasta = read.table('celegans/ce11-tRNAs.getfasta.tsv',stringsAsFactors = F)
colnames(Celegans_GtRNAdb_getfasta) = c('Loci','Loci_seq')
Celegans_piRNAdb_getfasta = read.table('celegans/Celegans_piRNAdb_ce11.getfasta.tsv',stringsAsFactors = F)
colnames(Celegans_piRNAdb_getfasta) = c('Loci','Loci_seq')
Celegans_rRNA_ucsc_getfasta = read.table('celegans/celegans_rRNA_ucsc.getfasta.tsv',stringsAsFactors = F)
colnames(Celegans_rRNA_ucsc_getfasta) = c('Loci','Loci_seq')
Celegans_tRFdb_getfasta = read.table('celegans/Celegans_tRFdb_ce11.getfasta.tsv',stringsAsFactors = F)
colnames(Celegans_tRFdb_getfasta) = c('Loci','Loci_seq')

```

## Build loci databases

```{r}

celegans_databases = list()
celegans_databases$miRBase = data.frame(
  ID = Celegans_miRBase_annotation$ID,
  Loci = str_c(Celegans_miRBase_annotation$V1,':',Celegans_miRBase_annotation$V4,'-',Celegans_miRBase_annotation$V5),
  Loci_chr = Celegans_miRBase_annotation$V1,
  Loci_start = Celegans_miRBase_annotation$V4,
  Loci_end = Celegans_miRBase_annotation$V5,
  Loci_len = as.numeric(Celegans_miRBase_annotation$V5) - as.numeric(Celegans_miRBase_annotation$V4),
  Loci_strand = Celegans_miRBase_annotation$V7,
  stringsAsFactors = F
)
celegans_databases$piRNAdb = data.frame(
  ID = Celegans_piRNAdb_annotation$V10,
  Loci = str_c(Celegans_piRNAdb_annotation$V1,':',Celegans_piRNAdb_annotation$V4,'-',Celegans_piRNAdb_annotation$V5),
  Loci_chr = Celegans_piRNAdb_annotation$V1,
  Loci_start = Celegans_piRNAdb_annotation$V4,
  Loci_end = Celegans_piRNAdb_annotation$V5,
  Loci_len = as.numeric(Celegans_piRNAdb_annotation$V5) - as.numeric(Celegans_piRNAdb_annotation$V4),
  Loci_strand = Celegans_piRNAdb_annotation$V7,
  stringsAsFactors = F
)
celegans_databases$GtRNAdb = data.frame(
  ID = Celegans_GtRNAdb_annotation$V4,
  Loci = str_c(Celegans_GtRNAdb_annotation$V1,':',Celegans_GtRNAdb_annotation$V2,'-',Celegans_GtRNAdb_annotation$V3),
  Loci_chr = Celegans_GtRNAdb_annotation$V1,
  Loci_start = Celegans_GtRNAdb_annotation$V2,
  Loci_end = Celegans_GtRNAdb_annotation$V3,
  Loci_len = as.numeric(Celegans_GtRNAdb_annotation$V3) - as.numeric(Celegans_GtRNAdb_annotation$V2),
  Loci_strand = Celegans_GtRNAdb_annotation$V6,
  stringsAsFactors = F
)
celegans_databases$rRNA_ucsc = data.frame(
  ID = Celegans_rRNA_ucsc_annotation$V4,
  Loci = str_c(Celegans_rRNA_ucsc_annotation$V1,':',Celegans_rRNA_ucsc_annotation$V2,'-',Celegans_rRNA_ucsc_annotation$V3),
  Loci_chr = Celegans_rRNA_ucsc_annotation$V1,
  Loci_start = Celegans_rRNA_ucsc_annotation$V2,
  Loci_end = Celegans_rRNA_ucsc_annotation$V3,
  Loci_len = as.numeric(Celegans_rRNA_ucsc_annotation$V3) - as.numeric(Celegans_rRNA_ucsc_annotation$V2),
  Loci_strand = Celegans_rRNA_ucsc_annotation$V6,
  stringsAsFactors = F
)
celegans_databases$tRFdb = data.frame(
  ID = str_c(Celegans_tRFdb_ce11_annotation$V5,'_',Celegans_tRFdb_ce11_annotation$V4),
  Loci = str_c(Celegans_tRFdb_ce11_annotation$V1,':',Celegans_tRFdb_ce11_annotation$V2,'-',Celegans_tRFdb_ce11_annotation$V3),
  Loci_chr = Celegans_tRFdb_ce11_annotation$V1,
  Loci_start = as.numeric(as.vector(Celegans_tRFdb_ce11_annotation$V2)),
  Loci_end = as.numeric(as.vector(Celegans_tRFdb_ce11_annotation$V3)),
  Loci_len = as.numeric(as.vector(Celegans_tRFdb_ce11_annotation$V3)) - as.numeric(as.vector(Celegans_tRFdb_ce11_annotation$V2)),
  Loci_strand = Celegans_tRFdb_ce11_annotation$V6,
  stringsAsFactors = F
)

```

---
## Add getfasta sequences to loci -> filtering -> reverse compliment sequences

```{r}

celegans_databases$miRBase$Loci_seq = foreach(l = celegans_databases$miRBase$Loci,.combine = 'c') %dopar% {
    Celegans_miRBase_getfasta[which(l == Celegans_miRBase_getfasta$Loci)[1],'Loci_seq']}

celegans_databases$GtRNAdb$Loci_seq = foreach(l = celegans_databases$GtRNAdb$Loci,.combine = 'c') %dopar% {
    Celegans_GtRNAdb_getfasta[which(l == Celegans_GtRNAdb_getfasta$Loci)[1],'Loci_seq']}
celegans_databases$GtRNAdb = filter(celegans_databases$GtRNAdb,!is.na(Loci_seq))

mean(celegans_databases$piRNAdb$Loci == Celegans_piRNAdb_getfasta$Loci)
celegans_databases$piRNAdb$Loci_seq = Celegans_piRNAdb_getfasta$Loci_seq

celegans_databases$rRNA_ucsc$Loci_seq = foreach(l = celegans_databases$rRNA_ucsc$Loci,.combine = 'c') %do% {
    Celegans_rRNA_ucsc_getfasta[which(l == Celegans_rRNA_ucsc_getfasta$Loci)[1],'Loci_seq']}
celegans_databases$rRNA_ucsc = filter(celegans_databases$rRNA_ucsc,!is.na(Loci_seq))

celegans_databases$tRFdb$Loci_seq = foreach(l = celegans_databases$tRFdb$Loci,.combine = 'c') %do% {
    Celegans_tRFdb_getfasta[which(l == Celegans_tRFdb_getfasta$Loci)[1],'Loci_seq']}

celegans_databases$miRBase$Loci_seq = str_to_upper(celegans_databases$miRBase$Loci_seq)
celegans_databases$piRNAdb$Loci_seq = str_to_upper(celegans_databases$piRNAdb$Loci_seq)
celegans_databases$GtRNAdb$Loci_seq = str_to_upper(celegans_databases$GtRNAdb$Loci_seq)
celegans_databases$rRNA_ucsc$Loci_seq = str_to_upper(celegans_databases$rRNA_ucsc$Loci_seq)
celegans_databases$tRFdb$Loci_seq = str_to_upper(celegans_databases$tRFdb$Loci_seq)

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

for (x in names(celegans_databases)) {
  print(x)
  celegans_databases[[x]] = na.omit(celegans_databases[[x]])
  celegans_databases[[x]] = mutate(celegans_databases[[x]],Loci_rev_seq = rev_comp(Loci_seq))
}

```

## Download fasta and sam files

```{r}
Celegans_miRBase_fasta = rbind.data.frame(
  read.table('celegans/miRBase_hairpin.fasta.tsv',stringsAsFactors = F,sep = '\t'),
  read.table('celegans/miRBase_mature.fasta.tsv',stringsAsFactors = F,sep = '\t'))
colnames(Celegans_miRBase_fasta) = c('ID','Fasta_seq')
Celegans_miRBase_fasta$ID = foreach(i = 1:nrow(Celegans_miRBase_fasta),.combine = 'c') %do% {
  strsplit(Celegans_miRBase_fasta$ID[i],' ')[[1]][1]
}

Celegans_GtRNAdb_fasta = read.table('celegans/ce11-tRNAs.fa.tsv',stringsAsFactors = F,sep = '\t')
colnames(Celegans_GtRNAdb_fasta) = c('ID','Fasta_seq')
Celegans_GtRNAdb_fasta$ID = foreach(i = Celegans_GtRNAdb_fasta$ID,.combine = 'c') %do% {
  strsplit(str_remove(i,'Caenorhabditis_elegans_'),' ')[[1]][1]
} 

Celegans_GtRNAdb_fasta_annotation = foreach(i = 1:nrow(Celegans_GtRNAdb_fasta),.combine = 'c') %do% {
  strsplit(Celegans_GtRNAdb_fasta$ID[i],' ')[[1]][11]
}

Celegans_piRNAdb_fasta = read.table('celegans/piRNAdb.cel.v1_7_6.fa.tsv',stringsAsFactors = F)
colnames(Celegans_piRNAdb_fasta) = c('ID','Fasta_seq')

Celegans_rRNA_ucsc_fasta = read.table('celegans/celegans_rRNA_ucsc.fasta.tsv',stringsAsFactors = F)
Celegans_rRNA_ucsc_fasta_annotation = str_remove_all(Celegans_rRNA_ucsc_fasta$V2,'range=')
Celegans_rRNA_ucsc_fasta = Celegans_rRNA_ucsc_fasta[,c(1,7)]
colnames(Celegans_rRNA_ucsc_fasta) = c('ID','Fasta_seq')

Celegans_miRBase_sam_stats = rbind.data.frame(
  read.table('celegans/miRBase_mature_sam_stats.txt',stringsAsFactors = F),
  read.table('celegans/miRBase_hairpin_sam_stats.txt',stringsAsFactors = F))
colnames(Celegans_miRBase_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Celegans_miRBase_sam_stats$Fasta_start = Celegans_miRBase_sam_stats$Fasta_start - 1
Celegans_miRBase_sam_stats$Fasta_mapped = foreach(i = Celegans_miRBase_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

Celegans_GtRNAdb_sam_stats = read.table('celegans/GtRNAdb_sam_stats.txt',stringsAsFactors = F)
colnames(Celegans_GtRNAdb_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Celegans_GtRNAdb_sam_stats$Fasta_start = Celegans_GtRNAdb_sam_stats$Fasta_start - 1
Celegans_GtRNAdb_sam_stats$ID = foreach(i = Celegans_GtRNAdb_sam_stats$ID,.combine = 'c') %do% {
  strsplit(i,'_')[[1]][3]
}
Celegans_GtRNAdb_sam_stats$Fasta_mapped = foreach(i = Celegans_GtRNAdb_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]}

Celegans_piRNAdb_sam_stats = read.table('celegans/piRNAdb_sam_stats.txt',stringsAsFactors = F)
colnames(Celegans_piRNAdb_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                           'Fasta_len','Fasta_mapped')
Celegans_piRNAdb_sam_stats$Fasta_start = Celegans_piRNAdb_sam_stats$Fasta_start - 1
Celegans_piRNAdb_sam_stats$Fasta_mapped = foreach(i = Celegans_piRNAdb_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

Celegans_rRNA_sam_stats = read.table('celegans/rRNA_ucsc_sam_stats.txt',stringsAsFactors = F)
colnames(Celegans_rRNA_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Celegans_rRNA_sam_stats$Fasta_start = Celegans_rRNA_sam_stats$Fasta_start - 1
Celegans_rRNA_sam_stats$ID = str_remove_all(Celegans_rRNA_sam_stats$ID,'rn7_rmsk_')
Celegans_rRNA_sam_stats$Fasta_mapped = foreach(i = Celegans_rRNA_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]}

Celegans_tRFdb_sam_stats = read.table('celegans/tRFdb_sam_stats.txt',stringsAsFactors = F)
colnames(Celegans_rRNA_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Celegans_rRNA_sam_stats$Fasta_start = Celegans_rRNA_sam_stats$Fasta_start - 1
Celegans_rRNA_sam_stats$ID = str_remove_all(Celegans_rRNA_sam_stats$ID,'rn7_rmsk_')
Celegans_rRNA_sam_stats$Fasta_mapped = foreach(i = Celegans_rRNA_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]}

```

## Add fasta sequences to loci -> check equal seq

```{r} 
celegans_databases$miRBase$Fasta_seq = NA
foreach(i = 1:nrow(celegans_databases$miRBase),.combine = 'c') %do% {
  celegans_databases$miRBase$Fasta_seq[i] = 
    Celegans_miRBase_fasta[which(Celegans_miRBase_fasta$ID == celegans_databases$miRBase$ID[i])[1],2]
}
  
celegans_databases$GtRNAdb$Fasta_seq = NA
foreach(i = 1:nrow(celegans_databases$GtRNAdb),.combine = 'c') %do% {
  celegans_databases$GtRNAdb$Fasta_seq[i] = 
    Celegans_GtRNAdb_fasta[which(Celegans_GtRNAdb_fasta$ID == celegans_databases$GtRNAdb$ID[i])[1],2]
}

celegans_databases$piRNAdb$Fasta_seq = NA
celegans_databases$piRNAdb$Fasta_seq = foreach(i = 1:nrow(celegans_databases$piRNAdb),.combine = 'c') %do% {
  if (i %% 100000 == 0) {print(i)}
  Celegans_piRNAdb_fasta[which(Celegans_piRNAdb_fasta$ID == celegans_databases$piRNAdb$ID[i])[1],2]
}

celegans_databases$rRNA_ucsc$Fasta_seq = foreach(i = 1:nrow(celegans_databases$rRNA_ucsc),.combine = 'c') %do% {
  Celegans_rRNA_ucsc_fasta[which(
    Celegans_rRNA_ucsc_fasta_annotation == 
      str_c(celegans_databases$rRNA_ucsc$Loci_chr[i],':',
            celegans_databases$rRNA_ucsc$Loci_start[i]+1,'-',
            celegans_databases$rRNA_ucsc$Loci_end[i]))[1],2]
}

celegans_databases$tRFdb$Fasta_seq = Celegans_trfdb_raw_table$tRF.Sequence

celegans_databases$miRBase$Fasta_len = foreach(i = strsplit(celegans_databases$miRBase$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
celegans_databases$GtRNAdb$Fasta_len = foreach(i = strsplit(celegans_databases$GtRNAdb$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
celegans_databases$piRNAdb$Fasta_len = foreach(i = strsplit(celegans_databases$piRNAdb$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
celegans_databases$rRNA_ucsc$Fasta_len = foreach(i = strsplit(celegans_databases$rRNA_ucsc$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
celegans_databases$tRFdb$Fasta_seq = str_remove_all(celegans_databases$tRFdb$Fasta_seq,' ')
celegans_databases$tRFdb$Fasta_len = foreach(i = strsplit(celegans_databases$tRFdb$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}

celegans_databases$miRBase$equal_seq = 
  celegans_databases$miRBase$Fasta_seq == celegans_databases$miRBase$Loci_seq
celegans_databases$GtRNAdb$equal_seq = 
  celegans_databases$GtRNAdb$Fasta_seq == celegans_databases$GtRNAdb$Loci_seq
celegans_databases$piRNAdb$equal_seq = 
  celegans_databases$piRNAdb$Fasta_seq == celegans_databases$piRNAdb$Loci_seq
celegans_databases$rRNA_ucsc$equal_seq = 
  celegans_databases$rRNA_ucsc$Fasta_seq == celegans_databases$rRNA_ucsc$Loci_seq
celegans_databases$tRFdb$equal_seq = 
  celegans_databases$tRFdb$Fasta_seq == celegans_databases$tRFdb$Loci_seq

celegans_databases$miRBase$equal_rev_seq = 
  celegans_databases$miRBase$Fasta_seq == celegans_databases$miRBase$Loci_rev_seq
celegans_databases$GtRNAdb$equal_rev_seq = 
  celegans_databases$GtRNAdb$Fasta_seq == celegans_databases$GtRNAdb$Loci_rev_seq
celegans_databases$piRNAdb$equal_rev_seq = 
  celegans_databases$piRNAdb$Fasta_seq == celegans_databases$piRNAdb$Loci_rev_seq
celegans_databases$rRNA_ucsc$equal_rev_seq = 
  celegans_databases$rRNA_ucsc$Fasta_seq == celegans_databases$rRNA_ucsc$Loci_rev_seq
celegans_databases$tRFdb$equal_rev_seq = 
  celegans_databases$tRFdb$Fasta_seq == celegans_databases$tRFdb$Loci_rev_seq

```

## Stats of united database

```{r}

write.xlsx(data.frame(
  Unique_ID_amount = c(
length(unique(celegans_databases$miRBase$ID)),
length(unique(celegans_databases$piRNAdb$ID)),
length(unique(celegans_databases$GtRNAdb$ID)),
length(unique(celegans_databases$rRNA_ucsc$ID)),
length(unique(celegans_databases$tRFdb$ID))
),
  Fasta_seqs_amount = c(
sum(!is.na(celegans_databases$miRBase$Fasta_seq)),
sum(!is.na(celegans_databases$piRNAdb$Fasta_seq)),
sum(!is.na(celegans_databases$GtRNAdb$Fasta_seq)),
sum(!is.na(celegans_databases$rRNA_ucsc$Fasta_seq)),
sum(!is.na(celegans_databases$tRFdb$Fasta_seq))
),
  Equal_Fasta_Loci_seqs = c(
sum(celegans_databases$miRBase$equal_seq,na.rm = T) + sum(celegans_databases$miRBase$equal_rev_seq,na.rm = T),
sum(celegans_databases$piRNAdb$equal_seq,na.rm = T) + sum(celegans_databases$piRNAdb$equal_rev_seq,na.rm = T),
sum(celegans_databases$GtRNAdb$equal_seq,na.rm = T) + sum(celegans_databases$GtRNAdb$equal_rev_seq,na.rm = T),
sum(celegans_databases$rRNA_ucsc$equal_seq,na.rm = T) + sum(celegans_databases$rRNA_ucsc$equal_rev_seq,na.rm = T),
sum(celegans_databases$tRFdb$equal_seq,na.rm = T) + sum(celegans_databases$tRFdb$equal_rev_seq,na.rm = T)
), row.names = names(celegans_databases)),
'celegans/smRNA_comparative_table_stats.xlsx',rowNames = T)

write.xlsx(celegans_databases,'celegans/smRNA_comparative_table_celegans.xlsx',rowNames = T)

```

## Pictures

```{r} 

venn.diagram(list(unique(celegans_databases$miRBase$ID), unique(Celegans_miRBase_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'celegans/Venn_celegans_miRBase_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red','blue'), cat.dist = c(0.05,0.065),
             main = 'Unique ID of miRBase', fill = c('red','blue'))

venn.diagram(list(unique(celegans_databases$piRNAdb$ID), unique(Celegans_piRNAdb_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'celegans/Venn_celegans_piRNAdb_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red','blue'), cat.dist = c(0.05,0.04),
             main = 'Unique ID of piRNAdb', fill = c('red','blue'))

venn.diagram(list(unique(celegans_databases$GtRNAdb$ID), unique(Celegans_GtRNAdb_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'celegans/Venn_celegans_GtRNAdb_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red','blue'), cat.dist = c(0.05,0.04),
             main = 'Unique ID of GtRNAdb', fill = c('red','blue'))

ggplot(celegans_databases$miRBase) + aes(Loci_len,Fasta_len) +
  theme_classic() + geom_point() + geom_abline() +
  ggtitle('Comparison of length of miRBase loci and fasta seqs')
ggsave('celegans/Comparison of celegans length of microRNA.tiff')

ggplot(celegans_databases$tRFdb) + aes(Loci_len,Fasta_len) +
  theme_classic() + geom_point() + geom_abline() +
  ggtitle('Comparison of length of tRFdb loci and fasta seqs')
ggsave('celegans/Comparison of celegans length of tRNA-derived.tiff')



```


## Building united database

```{r loci&fa} 
Celegans_miRBase_precursors_loci_fa = filter(celegans_databases$miRBase,equal_seq+equal_rev_seq == 1,
                                    str_detect(ID,'cel-mir'))
Celegans_miRBase_mature_loci_fa = filter(celegans_databases$miRBase,equal_seq+equal_rev_seq == 1,
                                    str_detect(ID,'cel-miR')) 
Celegans_piRNAdb_loci_fa = filter(celegans_databases$piRNAdb,equal_seq+equal_rev_seq == 1)
Celegans_GtRNAdb_loci_fa = filter(celegans_databases$GtRNAdb,equal_seq+equal_rev_seq == 1)
Celegans_rRNA_ucsc_loci_fa = filter(celegans_databases$rRNA_ucsc,equal_seq+equal_rev_seq == 1)
Celegans_tRFdb_loci_fa = filter(celegans_databases$tRFdb,equal_seq+equal_rev_seq == 1)

celegans_smRNA_DATABASE = rbind.data.frame(
  data.frame(type = 'mature_miRNA',
           ID = Celegans_miRBase_mature_loci_fa$ID,
           Annotation_consensus = Celegans_miRBase_mature_loci_fa$Loci,
           Annotation_chr = Celegans_miRBase_mature_loci_fa$Loci_chr,
           Annotation_start = Celegans_miRBase_mature_loci_fa$Loci_start,
           Annotation_end = Celegans_miRBase_mature_loci_fa$Loci_end,
           Sequence_consensus = Celegans_miRBase_mature_loci_fa$Fasta_seq,
           Length_consensus = Celegans_miRBase_mature_loci_fa$Loci_len,
           strand = Celegans_miRBase_mature_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
  data.frame(type = 'precursor_miRNA',
           ID = Celegans_miRBase_precursors_loci_fa$ID,
           Annotation_consensus = Celegans_miRBase_precursors_loci_fa$Loci,
           Annotation_chr = Celegans_miRBase_precursors_loci_fa$Loci_chr,
           Annotation_start = Celegans_miRBase_precursors_loci_fa$Loci_start,
           Annotation_end = Celegans_miRBase_precursors_loci_fa$Loci_end,
           Sequence_consensus = Celegans_miRBase_precursors_loci_fa$Fasta_seq,
           Length_consensus = Celegans_miRBase_precursors_loci_fa$Loci_len,
           strand = Celegans_miRBase_precursors_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'piRNA',
           ID = Celegans_piRNAdb_loci_fa$ID,
           Annotation_consensus = Celegans_piRNAdb_loci_fa$Loci,
           Annotation_chr = Celegans_piRNAdb_loci_fa$Loci_chr,
           Annotation_start = Celegans_piRNAdb_loci_fa$Loci_start,
           Annotation_end = Celegans_piRNAdb_loci_fa$Loci_end,
           Sequence_consensus = Celegans_piRNAdb_loci_fa$Fasta_seq,
           Length_consensus = Celegans_piRNAdb_loci_fa$Loci_len,
           strand = Celegans_piRNAdb_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'tRNA-derived',
           ID = Celegans_tRFdb_loci_fa$ID,
           Annotation_consensus = Celegans_tRFdb_loci_fa$Loci,
           Annotation_chr = Celegans_tRFdb_loci_fa$Loci_chr,
           Annotation_start = Celegans_tRFdb_loci_fa$Loci_start,
           Annotation_end = Celegans_tRFdb_loci_fa$Loci_end,
           Sequence_consensus = Celegans_tRFdb_loci_fa$Fasta_seq,
           Length_consensus = Celegans_tRFdb_loci_fa$Loci_len,
           strand = Celegans_tRFdb_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'tRNA',
           ID = Celegans_GtRNAdb_loci_fa$ID,
           Annotation_consensus = Celegans_GtRNAdb_loci_fa$Loci,
           Annotation_chr = Celegans_GtRNAdb_loci_fa$Loci_chr,
           Annotation_start = Celegans_GtRNAdb_loci_fa$Loci_start,
           Annotation_end = Celegans_GtRNAdb_loci_fa$Loci_end,
           Sequence_consensus = Celegans_GtRNAdb_loci_fa$Fasta_seq,
           Length_consensus = Celegans_GtRNAdb_loci_fa$Loci_len,
           strand = Celegans_GtRNAdb_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'rRNA',
           ID = Celegans_rRNA_ucsc_loci_fa$ID,
           Annotation_consensus = Celegans_rRNA_ucsc_loci_fa$Loci,
           Annotation_chr = Celegans_rRNA_ucsc_loci_fa$Loci_chr,
           Annotation_start = Celegans_rRNA_ucsc_loci_fa$Loci_start,
           Annotation_end = Celegans_rRNA_ucsc_loci_fa$Loci_end,
           Sequence_consensus = Celegans_rRNA_ucsc_loci_fa$Fasta_seq,
           Length_consensus = Celegans_rRNA_ucsc_loci_fa$Loci_len,
           strand = Celegans_rRNA_ucsc_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
)
)

```

```{r loci only}

Celegans_miRBase_precursors_loci_only = rbind.data.frame(
  filter(celegans_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'cel-mir')),
  filter(celegans_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'cel-l')))
Celegans_miRBase_precursors_loci_only = filter(Celegans_miRBase_precursors_loci_only,
                                            !(str_detect(ID,'p')))
Celegans_miRBase_precursors_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(Celegans_miRBase_precursors_loci_only),.combine = 'c') %do% {
    if (Celegans_miRBase_precursors_loci_only$Loci_strand[i] == '-') {out = Celegans_miRBase_precursors_loci_only$Loci_seq[i]}
    if (Celegans_miRBase_precursors_loci_only$Loci_strand[i] == '+') {out = Celegans_miRBase_precursors_loci_only$Loci_rev_seq[i]}
    out
  }

Celegans_miRBase_mature_loci_only = unique(rbind.data.frame(
  filter(celegans_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'cel-miR')),
  filter(celegans_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'p'))))
Celegans_miRBase_mature_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(Celegans_miRBase_mature_loci_only),.combine = 'c') %do% {
    if (Celegans_miRBase_mature_loci_only$Loci_strand[i] == '-') {out = Celegans_miRBase_mature_loci_only$Loci_seq[i]}
    if (Celegans_miRBase_mature_loci_only$Loci_strand[i] == '+') {out = Celegans_miRBase_mature_loci_only$Loci_rev_seq[i]}
    out
  }
Celegans_piRNAdb_loci_only = filter(celegans_databases$piRNAdb,is.na(Fasta_seq)) # is null

Celegans_GtRNAdb_loci_only = filter(celegans_databases$GtRNAdb,is.na(Fasta_seq))
Celegans_GtRNAdb_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(Celegans_GtRNAdb_loci_only),.combine = 'c') %do% {
    if (Celegans_GtRNAdb_loci_only$Loci_strand[i] == '-') {out = Celegans_GtRNAdb_loci_only$Loci_seq[i]}
    if (Celegans_GtRNAdb_loci_only$Loci_strand[i] == '+') {out = Celegans_GtRNAdb_loci_only$Loci_rev_seq[i]}
    out
  }
Celegans_tRFdb_loci_only = filter(celegans_databases$tRFdb,is.na(Fasta_seq)) # is null
Celegans_rRNA_ucsc_loci_only = filter(celegans_databases$rRNA_ucsc,is.na(Fasta_seq)) # is null

celegans_smRNA_DATABASE = rbind.data.frame(celegans_smRNA_DATABASE,
  data.frame(type = 'precursor_miRNA',
             ID = Celegans_miRBase_precursors_loci_only$ID,
             Annotation_consensus = Celegans_miRBase_precursors_loci_only$Loci,
             Annotation_chr = Celegans_miRBase_precursors_loci_only$Loci_chr,
             Annotation_start = Celegans_miRBase_precursors_loci_only$Loci_start,
             Annotation_end = Celegans_miRBase_precursors_loci_only$Loci_end,
             Sequence_consensus = Celegans_miRBase_precursors_loci_only$Fasta_seq,
             Length_consensus = Celegans_miRBase_precursors_loci_only$Loci_len,
             strand = Celegans_miRBase_precursors_loci_only$Loci_strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ),
  data.frame(type = 'mature_miRNA',
             ID = Celegans_miRBase_mature_loci_only$ID,
             Annotation_consensus = Celegans_miRBase_mature_loci_only$Loci,
             Annotation_chr = Celegans_miRBase_mature_loci_only$Loci_chr,
             Annotation_start = Celegans_miRBase_mature_loci_only$Loci_start,
             Annotation_end = Celegans_miRBase_mature_loci_only$Loci_end,
             Sequence_consensus = Celegans_miRBase_mature_loci_only$Fasta_seq,
             Length_consensus = Celegans_miRBase_mature_loci_only$Loci_len,
             strand = Celegans_miRBase_mature_loci_only$Loci_strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ),
  data.frame(type = 'tRNA',
             ID = Celegans_GtRNAdb_loci_only$ID,
             Annotation_consensus = Celegans_GtRNAdb_loci_only$Loci,
             Annotation_chr = Celegans_GtRNAdb_loci_only$Loci_chr,
             Annotation_start = Celegans_GtRNAdb_loci_only$Loci_start,
             Annotation_end = Celegans_GtRNAdb_loci_only$Loci_end,
             Sequence_consensus = Celegans_GtRNAdb_loci_only$Fasta_seq,
             Length_consensus = Celegans_GtRNAdb_loci_only$Loci_len,
             strand = Celegans_GtRNAdb_loci_only$Loci_strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ))

```

```{r fa only}

Celegans_miRBase_precursors_fa_only = filter(Celegans_miRBase_sam_stats, !(ID %in% celegans_databases$miRBase$ID), str_detect(ID,'cel-mir'), Fasta_chr != '*') # is null
Celegans_miRBase_mature_fa_only = filter(Celegans_miRBase_sam_stats, !(ID %in% celegans_databases$miRBase$ID),Fasta_chr != '*', str_detect(ID,'rno-miR')) # is null
Celegans_piRNAdb_fa_only = filter(Celegans_piRNAdb_sam_stats,!(ID %in% celegans_databases$piRNAdb$ID),
                         Fasta_chr != '*') # is null
Celegans_GtRNAdb_fa_only = filter(Celegans_GtRNAdb_sam_stats,!(ID %in% celegans_databases$GtRNAdb$ID),
                         Fasta_chr != '*') # is null
Celegans_rRNA_ucsc_fa_only = filter(Celegans_rRNA_sam_stats,!(ID %in% celegans_databases$rRNA_ucsc$ID),
                         Fasta_chr != '*')

celegans_smRNA_DATABASE = rbind.data.frame(celegans_smRNA_DATABASE,
                                  data.frame(type = 'rRNA',
                                             ID = Celegans_rRNA_ucsc_fa_only$ID,
        Annotation_consensus = str_c(Celegans_rRNA_ucsc_fa_only$Fasta_chr,':',
                                     Celegans_rRNA_ucsc_fa_only$Fasta_start,'-',
                                     as.character(Celegans_rRNA_ucsc_fa_only$Fasta_start+Celegans_rRNA_ucsc_fa_only$Fasta_len)),
                                             Annotation_chr = Celegans_rRNA_ucsc_fa_only$Fasta_chr,
                                             Annotation_start = Celegans_rRNA_ucsc_fa_only$Fasta_start,
                                             Annotation_end = Celegans_rRNA_ucsc_fa_only$Fasta_start+
          Celegans_rRNA_ucsc_fa_only$Fasta_len,
                                             Sequence_consensus = Celegans_rRNA_ucsc_fa_only$Fasta_seq,
                                             Length_consensus = Celegans_rRNA_ucsc_fa_only$Fasta_len,
        strand = '+',
        strand_from_annotation = F,
                                             Source = 'fasta',
                                             stringsAsFactors = F
                                  ))

```

```{r count H-dist & suff-pref intersect}

Celegans_miRBase_precursors_loci_not_fa = rbind.data.frame(
  filter(celegans_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'cel-mir')),
  filter(celegans_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'cel-let')),
  filter(celegans_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'cel-ban')))
Celegans_miRBase_precursors_loci_not_fa = filter(Celegans_miRBase_precursors_loci_not_fa,
                                              !(str_detect(ID,'p')))
Celegans_miRBase_mature_loci_not_fa = unique(rbind.data.frame(
  filter(celegans_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'cel-miR')),
  filter(celegans_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'p'))))
Celegans_GtRNAdb_loci_not_fa = filter(celegans_databases$GtRNAdb, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq)) # is null
Celegans_rRNA_ucsc_loci_not_fa = filter(celegans_databases$rRNA_ucsc, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq)) # is null
Celegans_piRNAdb_loci_not_fa = filter(celegans_databases$piRNAdb, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq))
Celegans_tRFdb_loci_not_fa = filter(celegans_databases$tRFdb, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq))

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

Celegans_loci_not_fa = foreach(d = list(Celegans_miRBase_precursors_loci_not_fa,
Celegans_miRBase_mature_loci_not_fa,
Celegans_tRFdb_loci_not_fa)) %do% {
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
names(Celegans_loci_not_fa) = c('miRBase_precursors','miRBase_mature','tRFdb')

Celegans_loci_not_fa$piRNAdb = Celegans_piRNAdb_loci_not_fa
d = Celegans_loci_not_fa$piRNAdb
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
Celegans_loci_not_fa$piRNAdb$ID = str_c(Celegans_loci_not_fa$piRNAdb$Loci_seq,'_',Celegans_loci_not_fa$piRNAdb$Fasta_seq)
Celegans_loci_not_fa$piRNAdb = 
  cbind.data.frame(Celegans_loci_not_fa$piRNAdb,
                   foreach(i = 1:nrow(Celegans_loci_not_fa$piRNAdb),.combine = 'rbind') %do% {
                     if (i %% 10000 == 0) {print(i)}
                     k = Celegans_loci_not_fa$piRNAdb$ID[i]
                     d[which(d$ID == k),4:6]
                   }
                   )

```

```{r fix loci not fa}
Celegans_loci_not_fa$miRBase_precursors$type = 'precursor_miRNA'
Celegans_loci_not_fa$miRBase_mature$type = 'mature_miRNA'
Celegans_loci_not_fa$tRFdb$type = 'tRNA-derived'
Celegans_loci_not_fa$piRNAdb$type = 'piRNA'

Celegans_loci_not_fa$piRNAdb$ID = Celegans_piRNAdb_loci_not_fa$ID


for (i in 1:3) {
  d = Celegans_loci_not_fa[[i]]
  celegans_smRNA_DATABASE = rbind.data.frame(celegans_smRNA_DATABASE, 
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
celegans_smRNA_DATABASE = filter(celegans_smRNA_DATABASE,Source != 'trash')

d = Celegans_loci_not_fa[[4]]
d$diff_len = d$Fasta_len - d$Loci_len
d$fix_p2 = d$Fasta_len - d$Loci_len == as.numeric(d$p3)
d$p2[which(d$fix_p2 == T)] = '0'

wide_annot = filter(d,p1 %in% c('loci_forward_in_fa'))
Annotation_start = wide_annot$Loci_start-as.numeric(wide_annot$p2)
Annotation_end = wide_annot$Loci_end+as.numeric(wide_annot$p3)
Sequence_consensus = wide_annot$Fasta_seq
Length_consensus = wide_annot$Fasta_len
Source = 'wide_annotation&fasta'
celegans_smRNA_DATABASE = rbind.data.frame(celegans_smRNA_DATABASE, 
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
celegans_smRNA_DATABASE = rbind.data.frame(celegans_smRNA_DATABASE, 
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
celegans_smRNA_DATABASE = rbind.data.frame(celegans_smRNA_DATABASE, 
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
celegans_smRNA_DATABASE = rbind.data.frame(celegans_smRNA_DATABASE, 
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

write.table(str_c(filter(celegans_smRNA_DATABASE,type == 'mature_miRNA')$Annotation_chr,filter(celegans_smRNA_DATABASE,type == 'mature_miRNA')$Annotation_start,filter(celegans_smRNA_DATABASE,type == 'mature_miRNA')$Annotation_end,filter(celegans_smRNA_DATABASE,type == 'mature_miRNA')$ID,sep = '\t'),'celegans/celegans_smRNA_DATABASE_miRNA.bed',col.names = F, row.names = F)
write.table(str_c(filter(celegans_smRNA_DATABASE,type == 'piRNA')$Annotation_chr,filter(celegans_smRNA_DATABASE,type == 'piRNA')$Annotation_start,filter(celegans_smRNA_DATABASE,type == 'piRNA')$Annotation_end,filter(celegans_smRNA_DATABASE,type == 'piRNA')$ID,sep = '\t'),'celegans/celegans_smRNA_DATABASE_piRNA.bed',col.names = F, row.names = F)
write.table(str_c(filter(celegans_smRNA_DATABASE,type == 'tRNA')$Annotation_chr,filter(celegans_smRNA_DATABASE,type == 'tRNA')$Annotation_start,filter(celegans_smRNA_DATABASE,type == 'tRNA')$Annotation_end,filter(celegans_smRNA_DATABASE,type == 'tRNA')$ID,sep = '\t'),'celegans/celegans_smRNA_DATABASE_tRNA.bed',col.names = F, row.names = F)
write.table(str_c(filter(celegans_smRNA_DATABASE,type == 'rRNA')$Annotation_chr,filter(celegans_smRNA_DATABASE,type == 'rRNA')$Annotation_start,filter(celegans_smRNA_DATABASE,type == 'rRNA')$Annotation_end,filter(celegans_smRNA_DATABASE,type == 'rRNA')$ID,sep = '\t'),'celegans/celegans_smRNA_DATABASE_rRNA.bed',col.names = F, row.names = F)

```

```{r filter by intersept}

celegans_inters = rbind.data.frame(
  unique(read.table('celegans/miRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('celegans/miRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('celegans/miRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
  unique(read.table('celegans/piRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('celegans/piRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('celegans/piRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('celegans/rRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('celegans/rRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('celegans/rRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('celegans/tRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('celegans/tRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('celegans/tRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('celegans/miRNA_inter_out.txt',stringsAsFactors = F)[,1:4]),
unique(read.table('celegans/piRNA_inter_out.txt',stringsAsFactors = F)[,1:4]),
unique(read.table('celegans/tRNA_inter_out.txt',stringsAsFactors = F)[,1:4]),
unique(read.table('celegans/rRNA_inter_out.txt',stringsAsFactors = F)[,1:4])
)
celegans_inters = unique(celegans_inters)
celegans_inters_vector = str_c(celegans_inters$V1,celegans_inters$V2,celegans_inters$V3,celegans_inters$V4)

celegans_smRNA_DATABASE$filtered = foreach(i = 1:nrow(celegans_smRNA_DATABASE),.combine = 'c',.packages = libs) %dopar% {
  str_c(celegans_smRNA_DATABASE$Annotation_chr[i],
        celegans_smRNA_DATABASE$Annotation_start[i],
        celegans_smRNA_DATABASE$Annotation_end[i],celegans_smRNA_DATABASE$ID[i]) %in% celegans_inters_vector
}

write.xlsx(filter(celegans_smRNA_DATABASE, type %in% c('mature_miRNA', 'piRNA', 'tRNA'),filtered == F),'Celegans/celegans_smRNA_DATABASE_mipitr_filtered.xlsx')

write.table(str_c('>',filter(celegans_smRNA_DATABASE,type %in% c('mature_miRNA','piRNA','rRNA','tRNA'), filtered == F,filtered_in == F)$ID,'\n',filter(celegans_smRNA_DATABASE,type %in% c('mature_miRNA','piRNA','rRNA','tRNA'), filtered == F,filtered_in == F)$Sequence_consensus,collapse = '\n'),'Celegans/celegans_smRNA_DATABASE_filtered.fasta',row.names = F,col.names = F)

```

```{r write DATABASE}

Celegans_Source_df = celegans_smRNA_DATABASE %>% group_by(type) %>% 
  summarise(n = table(Source),
            Source = names(table(Source)))
write.xlsx(pivot_wider(Celegans_Source_df, names_from = type, values_from = n),'celegans/celegans_smRNA_DATABASE_stats.xlsx')
write.xlsx(celegans_smRNA_DATABASE,'celegans/celegans_smRNA_DATABASE.xlsx')
write.table(str_c('>',celegans_smRNA_DATABASE$ID[which(celegans_smRNA_DATABASE$type %in% c('mature_miRNA','piRNA','rRNA','tRNA'))],'\n',celegans_smRNA_DATABASE$Sequence_consensus[which(celegans_smRNA_DATABASE$type %in% c('mature_miRNA','piRNA','rRNA','tRNA'))],collapse = '\n'),'celegans/smRNA_DATABASE.fasta',row.names = F,col.names = F)

celegans_inters_in = rbind.data.frame(
  unique(read.table('celegans/miRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('celegans/miRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('celegans/miRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
  unique(read.table('celegans/piRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('celegans/piRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('celegans/piRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('celegans/rRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('celegans/rRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('celegans/rRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('celegans/tRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('celegans/tRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('celegans/tRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]))
celegans_inters_in = unique(celegans_inters_in)
celegans_inters_in_vector = str_c(celegans_inters_in$V1,celegans_inters_in$V2,celegans_inters_in$V3,celegans_inters_in$V4)
celegans_smRNA_DATABASE$filtered_in = foreach(i = 1:nrow(celegans_smRNA_DATABASE),.combine = 'c',.packages = libs) %dopar% {
  o = str_c(celegans_smRNA_DATABASE$Annotation_chr[i],
        celegans_smRNA_DATABASE$Annotation_start[i],
        celegans_smRNA_DATABASE$Annotation_end[i],celegans_smRNA_DATABASE$ID[i]) %in% celegans_inters_in_vector
  if (celegans_smRNA_DATABASE$type[i] == 'rRNA' && !(str_detect(celegans_smRNA_DATABASE$ID[i],'rRNA'))) {o = T}
  o
}

celegans_smRNA_DATABASE_rRNA = filter(celegans_smRNA_DATABASE,type == 'rRNA', str_detect(ID,'rRNA'),Source == 'annotation&fasta')
write.xlsx(celegans_smRNA_DATABASE_rRNA,'Celegans/celegans_smRNA_DATABASE_rRNA.xlsx')
for (i in unique(celegans_smRNA_DATABASE$type)) {
  write.xlsx(filter(celegans_smRNA_DATABASE, filtered_in == F, type == i), str_c('Celegans/celegans_smRNA_DATABASE_',i,'.xlsx'))
}
```
