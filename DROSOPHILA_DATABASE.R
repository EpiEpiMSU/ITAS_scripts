---

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

Drosophila_miRBase_annotation = read.table('drosophila/dme.gff3',stringsAsFactors = F)
Drosophila_miRBase_annotation$ID = foreach(i = 1:nrow(Drosophila_miRBase_annotation),.combine = 'c') %do% {
  str_remove(strsplit(Drosophila_miRBase_annotation$V9[i],';')[[1]][3],'Name=')
}
Drosophila_piRNAdb_annotation = read.table('drosophila/pirnadb.v1_7_6.dm6.gtf',stringsAsFactors = F)
Drosophila_GtRNAdb_annotation = read.table('drosophila/dm6-tRNAs.bed',stringsAsFactors = F)
Drosophila_rRNA_ucsc_annotation = read.table('drosophila/drosophila_rRNA_ucsc.bed',stringsAsFactors = F)
Drosophila_trfdb_raw_table = rbind.data.frame(
  read.xlsx('drosophila/trfdb_excel_result_1.xlsx'),
  read.xlsx('drosophila/trfdb_excel_result_2.xlsx'),
  read.xlsx('drosophila/trfdb_excel_result_3.xlsx')
  )
write.table(
  unique(str_c('>',Drosophila_trfdb_raw_table$Type,'_',Drosophila_trfdb_raw_table$tRF.ID,'\n',Drosophila_trfdb_raw_table$tRF.Sequence)),'drosophila/tRFdb.fasta',row.names = F,col.names = F
)

```


---
## Build bed-files for trf RNA -> liftover -> getfasta

```{r}

Drosophila_trfdb_annotation = foreach(i = 1:nrow(Drosophila_trfdb_raw_table),.combine = 'rbind.data.frame') %do% {
  p = unlist(c(strsplit(Drosophila_trfdb_raw_table$`tRNA.Gene.Co-ordinates`[i],'-')[[1]],
  Drosophila_trfdb_raw_table[i,c(1,3,8)]))
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
Drosophila_trfdb_annotation$chr = str_remove(Drosophila_trfdb_annotation$chr,' ')

  
write.table(Drosophila_trfdb_annotation
,'drosophila/Drosophila_tRFdb_dm3.bed', col.names = F,
            row.names = F, sep = '\t')

Drosophila_tRFdb_dm6_annotation = read.table('drosophila/Drosophila_tRFdb_dm6.bed',stringsAsFactors = F)
Drosophila_piRNAdb_annotation$V4[which(Drosophila_piRNAdb_annotation$V7 == '+')] = Drosophila_piRNAdb_annotation$V4[which(Drosophila_piRNAdb_annotation$V7 == '+')] - 1
Drosophila_piRNAdb_annotation$V5[which(Drosophila_piRNAdb_annotation$V7 == '-')] = Drosophila_piRNAdb_annotation$V5[which(Drosophila_piRNAdb_annotation$V7 == '-')] + 1
write.table(Drosophila_piRNAdb_annotation[,c(1,4,5,7,10)],'drosophila/Drosophila_piRNAdb_dm6.bed', col.names = F,row.names = F, sep = '\t')
Drosophila_miRBase_annotation$V4[which(Drosophila_miRBase_annotation$V7 == '-')] = Drosophila_miRBase_annotation$V4[which(Drosophila_miRBase_annotation$V7 == '-')] - 1
Drosophila_miRBase_annotation$V5[which(Drosophila_miRBase_annotation$V7 == '+')] = Drosophila_miRBase_annotation$V5[which(Drosophila_miRBase_annotation$V7 == '+')] + 1
write.table(Drosophila_miRBase_annotation[,c(1,4,5,7,10)],'drosophila/Drosophila_miRBase_dm6.bed', col.names = F,row.names = F, sep = '\t')

Drosophila_miRBase_getfasta = read.table('drosophila/Drosophila_miRBase_dm6.getfasta.tsv',stringsAsFactors = F)
colnames(Drosophila_miRBase_getfasta) = c('Loci','Loci_seq')
Drosophila_GtRNAdb_getfasta = read.table('drosophila/dm6-tRNAs.getfasta.fasta.tsv',stringsAsFactors = F)
colnames(Drosophila_GtRNAdb_getfasta) = c('Loci','Loci_seq')
Drosophila_piRNAdb_getfasta = read.table('drosophila/Drosophila_piRNAdb_dm6.getfasta.tsv',stringsAsFactors = F)
colnames(Drosophila_piRNAdb_getfasta) = c('Loci','Loci_seq')
Drosophila_rRNA_ucsc_getfasta = read.table('drosophila/drosophila_rRNA_ucsc.getfasta.fasta.tsv',stringsAsFactors = F)
colnames(Drosophila_rRNA_ucsc_getfasta) = c('Loci','Loci_seq')
Drosophila_tRFdb_getfasta = read.table('drosophila/Drosophila_tRFdb_dm6.getfasta.fasta.tsv',stringsAsFactors = F)
colnames(Drosophila_tRFdb_getfasta) = c('Loci','Loci_seq')

```

## Build loci databases

```{r}

drosophila_databases = list()
drosophila_databases$miRBase = data.frame(
  ID = Drosophila_miRBase_annotation$ID,
  Loci = str_c(Drosophila_miRBase_annotation$V1,':',Drosophila_miRBase_annotation$V4,'-',Drosophila_miRBase_annotation$V5),
  Loci_chr = Drosophila_miRBase_annotation$V1,
  Loci_start = Drosophila_miRBase_annotation$V4,
  Loci_end = Drosophila_miRBase_annotation$V5,
  Loci_len = as.numeric(Drosophila_miRBase_annotation$V5) - as.numeric(Drosophila_miRBase_annotation$V4),
  Loci_strand = Drosophila_miRBase_annotation$V7,
  stringsAsFactors = F
)
drosophila_databases$piRNAdb = data.frame(
  ID = Drosophila_piRNAdb_annotation$V10,
  Loci = str_c(Drosophila_piRNAdb_annotation$V1,':',Drosophila_piRNAdb_annotation$V4,'-',Drosophila_piRNAdb_annotation$V5),
  Loci_chr = Drosophila_piRNAdb_annotation$V1,
  Loci_start = Drosophila_piRNAdb_annotation$V4,
  Loci_end = Drosophila_piRNAdb_annotation$V5,
  Loci_len = as.numeric(Drosophila_piRNAdb_annotation$V5) - as.numeric(Drosophila_piRNAdb_annotation$V4),
  Loci_strand = Drosophila_piRNAdb_annotation$V7,
  stringsAsFactors = F
)
drosophila_databases$GtRNAdb = data.frame(
  ID = Drosophila_GtRNAdb_annotation$V4,
  Loci = str_c(Drosophila_GtRNAdb_annotation$V1,':',Drosophila_GtRNAdb_annotation$V2,'-',Drosophila_GtRNAdb_annotation$V3),
  Loci_chr = Drosophila_GtRNAdb_annotation$V1,
  Loci_start = Drosophila_GtRNAdb_annotation$V2,
  Loci_end = Drosophila_GtRNAdb_annotation$V3,
  Loci_len = as.numeric(Drosophila_GtRNAdb_annotation$V3) - as.numeric(Drosophila_GtRNAdb_annotation$V2),
  Loci_strand = Drosophila_GtRNAdb_annotation$V6,
  stringsAsFactors = F
)
drosophila_databases$rRNA_ucsc = data.frame(
  ID = Drosophila_rRNA_ucsc_annotation$V4,
  Loci = str_c(Drosophila_rRNA_ucsc_annotation$V1,':',Drosophila_rRNA_ucsc_annotation$V2,'-',Drosophila_rRNA_ucsc_annotation$V3),
  Loci_chr = Drosophila_rRNA_ucsc_annotation$V1,
  Loci_start = Drosophila_rRNA_ucsc_annotation$V2,
  Loci_end = Drosophila_rRNA_ucsc_annotation$V3,
  Loci_len = as.numeric(Drosophila_rRNA_ucsc_annotation$V3) - as.numeric(Drosophila_rRNA_ucsc_annotation$V2),
  Loci_strand = Drosophila_rRNA_ucsc_annotation$V6,
  stringsAsFactors = F
)
drosophila_databases$tRFdb = data.frame(
  ID = str_c(Drosophila_tRFdb_dm6_annotation$V5,'_',Drosophila_tRFdb_dm6_annotation$V4),
  Loci = str_c(Drosophila_tRFdb_dm6_annotation$V1,':',Drosophila_tRFdb_dm6_annotation$V2,'-',Drosophila_tRFdb_dm6_annotation$V3),
  Loci_chr = Drosophila_tRFdb_dm6_annotation$V1,
  Loci_start = as.numeric(as.vector(Drosophila_tRFdb_dm6_annotation$V2)),
  Loci_end = as.numeric(as.vector(Drosophila_tRFdb_dm6_annotation$V3)),
  Loci_len = as.numeric(as.vector(Drosophila_tRFdb_dm6_annotation$V3)) - as.numeric(as.vector(Drosophila_tRFdb_dm6_annotation$V2)),
  Loci_strand = Drosophila_tRFdb_dm6_annotation$V6,
  stringsAsFactors = F
)

```

---
## Add getfasta sequences to loci -> filtering -> reverse compliment sequences

```{r}

drosophila_databases$miRBase$Loci_seq = foreach(l = drosophila_databases$miRBase$Loci,.combine = 'c') %dopar% {
    Drosophila_miRBase_getfasta[which(l == Drosophila_miRBase_getfasta$Loci)[1],'Loci_seq']}

drosophila_databases$GtRNAdb$Loci_seq = foreach(l = drosophila_databases$GtRNAdb$Loci,.combine = 'c') %dopar% {
    Drosophila_GtRNAdb_getfasta[which(l == Drosophila_GtRNAdb_getfasta$Loci)[1],'Loci_seq']}
drosophila_databases$GtRNAdb = filter(drosophila_databases$GtRNAdb,!is.na(Loci_seq))

mean(drosophila_databases$piRNAdb$Loci == Drosophila_piRNAdb_getfasta$Loci)
drosophila_databases$piRNAdb$Loci_seq = Drosophila_piRNAdb_getfasta$Loci_seq

drosophila_databases$rRNA_ucsc$Loci_seq = foreach(l = drosophila_databases$rRNA_ucsc$Loci,.combine = 'c') %do% {
    Drosophila_rRNA_ucsc_getfasta[which(l == Drosophila_rRNA_ucsc_getfasta$Loci)[1],'Loci_seq']}
drosophila_databases$rRNA_ucsc = filter(drosophila_databases$rRNA_ucsc,!is.na(Loci_seq))

drosophila_databases$tRFdb$Loci_seq = foreach(l = drosophila_databases$tRFdb$Loci,.combine = 'c') %do% {
    Drosophila_tRFdb_getfasta[which(l == Drosophila_tRFdb_getfasta$Loci)[1],'Loci_seq']}

drosophila_databases$miRBase$Loci_seq = str_to_upper(drosophila_databases$miRBase$Loci_seq)
drosophila_databases$piRNAdb$Loci_seq = str_to_upper(drosophila_databases$piRNAdb$Loci_seq)
drosophila_databases$GtRNAdb$Loci_seq = str_to_upper(drosophila_databases$GtRNAdb$Loci_seq)
drosophila_databases$rRNA_ucsc$Loci_seq = str_to_upper(drosophila_databases$rRNA_ucsc$Loci_seq)
drosophila_databases$tRFdb$Loci_seq = str_to_upper(drosophila_databases$tRFdb$Loci_seq)

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

for (x in names(drosophila_databases)) {
  print(x)
  drosophila_databases[[x]] = na.omit(drosophila_databases[[x]])
  drosophila_databases[[x]] = mutate(drosophila_databases[[x]],Loci_rev_seq = rev_comp(Loci_seq))
}

```

## Download fasta and sam files

```{r}
Drosophila_miRBase_fasta = rbind.data.frame(
  read.table('drosophila/miRBase_hairpin.fasta.tsv',stringsAsFactors = F,sep = '\t'),
  read.table('drosophila/miRBase_mature.fasta.tsv',stringsAsFactors = F,sep = '\t'))
colnames(Drosophila_miRBase_fasta) = c('ID','Fasta_seq')
Drosophila_miRBase_fasta$ID = foreach(i = 1:nrow(Drosophila_miRBase_fasta),.combine = 'c') %do% {
  strsplit(Drosophila_miRBase_fasta$ID[i],' ')[[1]][1]
}

Drosophila_GtRNAdb_fasta = read.table('drosophila/dm6-tRNAs.fa.tsv',stringsAsFactors = F,sep = '\t')
colnames(Drosophila_GtRNAdb_fasta) = c('ID','Fasta_seq')
Drosophila_GtRNAdb_fasta$ID = foreach(i = Drosophila_GtRNAdb_fasta$ID,.combine = 'c') %do% {
  strsplit(str_remove(i,'Drosophila_melanogaster_'),' ')[[1]][1]
} 

Drosophila_GtRNAdb_fasta_annotation = foreach(i = 1:nrow(Drosophila_GtRNAdb_fasta),.combine = 'c') %do% {
  strsplit(Drosophila_GtRNAdb_fasta$ID[i],' ')[[1]][11]
}

Drosophila_piRNAdb_fasta = read.table('drosophila/piRNAdb.dme.v1_7_6.fa.tsv',stringsAsFactors = F)
colnames(Drosophila_piRNAdb_fasta) = c('ID','Fasta_seq')

Drosophila_rRNA_ucsc_fasta = read.table('drosophila/drosophila_rRNA_ucsc.fasta.tsv',stringsAsFactors = F)
Drosophila_rRNA_ucsc_fasta_annotation = str_remove_all(Drosophila_rRNA_ucsc_fasta$V2,'range=')
Drosophila_rRNA_ucsc_fasta = Drosophila_rRNA_ucsc_fasta[,c(1,7)]
colnames(Drosophila_rRNA_ucsc_fasta) = c('ID','Fasta_seq')

Drosophila_miRBase_sam_stats = rbind.data.frame(
  read.table('drosophila/miRBase_mature_sam_stats.txt',stringsAsFactors = F),
  read.table('drosophila/miRBase_hairpin_sam_stats.txt',stringsAsFactors = F))
colnames(Drosophila_miRBase_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Drosophila_miRBase_sam_stats$Fasta_start = Drosophila_miRBase_sam_stats$Fasta_start - 1
Drosophila_miRBase_sam_stats$Fasta_mapped = foreach(i = Drosophila_miRBase_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

Drosophila_GtRNAdb_sam_stats = read.table('drosophila/GtRNAdb_sam_stats.txt',stringsAsFactors = F)
colnames(Drosophila_GtRNAdb_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Drosophila_GtRNAdb_sam_stats$Fasta_start = Drosophila_GtRNAdb_sam_stats$Fasta_start - 1
Drosophila_GtRNAdb_sam_stats$ID = foreach(i = Drosophila_GtRNAdb_sam_stats$ID,.combine = 'c') %do% {
  strsplit(i,'_')[[1]][3]
}
Drosophila_GtRNAdb_sam_stats$Fasta_mapped = foreach(i = Drosophila_GtRNAdb_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]}

Drosophila_piRNAdb_sam_stats = read.table('drosophila/piRNAdb_sam_stats.txt',stringsAsFactors = F)
colnames(Drosophila_piRNAdb_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                           'Fasta_len','Fasta_mapped')
Drosophila_piRNAdb_sam_stats$Fasta_start = Drosophila_piRNAdb_sam_stats$Fasta_start - 1
Drosophila_piRNAdb_sam_stats$Fasta_mapped = foreach(i = Drosophila_piRNAdb_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

Drosophila_rRNA_sam_stats = read.table('drosophila/rRNA_ucsc_sam_stats.txt',stringsAsFactors = F)
colnames(Drosophila_rRNA_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Drosophila_rRNA_sam_stats$Fasta_start = Drosophila_rRNA_sam_stats$Fasta_start - 1
Drosophila_rRNA_sam_stats$ID = str_remove_all(Drosophila_rRNA_sam_stats$ID,'rn7_rmsk_')
Drosophila_rRNA_sam_stats$Fasta_mapped = foreach(i = Drosophila_rRNA_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]}

Drosophila_tRFdb_sam_stats = read.table('drosophila/tRFdb_sam_stats.txt',stringsAsFactors = F)
colnames(Drosophila_rRNA_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Drosophila_rRNA_sam_stats$Fasta_start = Drosophila_rRNA_sam_stats$Fasta_start - 1
Drosophila_rRNA_sam_stats$ID = str_remove_all(Drosophila_rRNA_sam_stats$ID,'rn7_rmsk_')
Drosophila_rRNA_sam_stats$Fasta_mapped = foreach(i = Drosophila_rRNA_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]}

```

## Add fasta sequences to loci -> check equal seq

```{r} 
drosophila_databases$miRBase$Fasta_seq = NA
foreach(i = 1:nrow(drosophila_databases$miRBase),.combine = 'c') %do% {
  drosophila_databases$miRBase$Fasta_seq[i] = 
    Drosophila_miRBase_fasta[which(Drosophila_miRBase_fasta$ID == drosophila_databases$miRBase$ID[i])[1],2]
}
  
drosophila_databases$GtRNAdb$Fasta_seq = NA
foreach(i = 1:nrow(drosophila_databases$GtRNAdb),.combine = 'c') %do% {
  drosophila_databases$GtRNAdb$Fasta_seq[i] = 
    Drosophila_GtRNAdb_fasta[which(Drosophila_GtRNAdb_fasta$ID == drosophila_databases$GtRNAdb$ID[i])[1],2]
}

drosophila_databases$piRNAdb$Fasta_seq = NA
drosophila_databases$piRNAdb$Fasta_seq = foreach(i = 1:nrow(drosophila_databases$piRNAdb),.combine = 'c') %do% {
  if (i %% 100000 == 0) {print(i)}
  Drosophila_piRNAdb_fasta[which(Drosophila_piRNAdb_fasta$ID == drosophila_databases$piRNAdb$ID[i])[1],2]
}

drosophila_databases$rRNA_ucsc$Fasta_seq = foreach(i = 1:nrow(drosophila_databases$rRNA_ucsc),.combine = 'c') %do% {
  Drosophila_rRNA_ucsc_fasta[which(
    Drosophila_rRNA_ucsc_fasta_annotation == 
      str_c(drosophila_databases$rRNA_ucsc$Loci_chr[i],':',
            drosophila_databases$rRNA_ucsc$Loci_start[i]+1,'-',
            drosophila_databases$rRNA_ucsc$Loci_end[i]))[1],2]
}

drosophila_databases$tRFdb$Fasta_seq = Drosophila_trfdb_raw_table$tRF.Sequence

drosophila_databases$miRBase$Fasta_len = foreach(i = strsplit(drosophila_databases$miRBase$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
drosophila_databases$GtRNAdb$Fasta_len = foreach(i = strsplit(drosophila_databases$GtRNAdb$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
drosophila_databases$piRNAdb$Fasta_len = foreach(i = strsplit(drosophila_databases$piRNAdb$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
drosophila_databases$rRNA_ucsc$Fasta_len = foreach(i = strsplit(drosophila_databases$rRNA_ucsc$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
drosophila_databases$tRFdb$Fasta_seq = str_remove_all(drosophila_databases$tRFdb$Fasta_seq,' ')
drosophila_databases$tRFdb$Fasta_len = foreach(i = strsplit(drosophila_databases$tRFdb$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}

drosophila_databases$miRBase$equal_seq = 
  drosophila_databases$miRBase$Fasta_seq == drosophila_databases$miRBase$Loci_seq
drosophila_databases$GtRNAdb$equal_seq = 
  drosophila_databases$GtRNAdb$Fasta_seq == drosophila_databases$GtRNAdb$Loci_seq
drosophila_databases$piRNAdb$equal_seq = 
  drosophila_databases$piRNAdb$Fasta_seq == drosophila_databases$piRNAdb$Loci_seq
drosophila_databases$rRNA_ucsc$equal_seq = 
  drosophila_databases$rRNA_ucsc$Fasta_seq == drosophila_databases$rRNA_ucsc$Loci_seq
drosophila_databases$tRFdb$equal_seq = 
  drosophila_databases$tRFdb$Fasta_seq == drosophila_databases$tRFdb$Loci_seq

drosophila_databases$miRBase$equal_rev_seq = 
  drosophila_databases$miRBase$Fasta_seq == drosophila_databases$miRBase$Loci_rev_seq
drosophila_databases$GtRNAdb$equal_rev_seq = 
  drosophila_databases$GtRNAdb$Fasta_seq == drosophila_databases$GtRNAdb$Loci_rev_seq
drosophila_databases$piRNAdb$equal_rev_seq = 
  drosophila_databases$piRNAdb$Fasta_seq == drosophila_databases$piRNAdb$Loci_rev_seq
drosophila_databases$rRNA_ucsc$equal_rev_seq = 
  drosophila_databases$rRNA_ucsc$Fasta_seq == drosophila_databases$rRNA_ucsc$Loci_rev_seq
drosophila_databases$tRFdb$equal_rev_seq = 
  drosophila_databases$tRFdb$Fasta_seq == drosophila_databases$tRFdb$Loci_rev_seq

```

## Stats of united database

```{r}

write.xlsx(data.frame(
  Unique_ID_amount = c(
length(unique(drosophila_databases$miRBase$ID)),
length(unique(drosophila_databases$piRNAdb$ID)),
length(unique(drosophila_databases$GtRNAdb$ID)),
length(unique(drosophila_databases$rRNA_ucsc$ID)),
length(unique(drosophila_databases$tRFdb$ID))
),
  Fasta_seqs_amount = c(
sum(!is.na(drosophila_databases$miRBase$Fasta_seq)),
sum(!is.na(drosophila_databases$piRNAdb$Fasta_seq)),
sum(!is.na(drosophila_databases$GtRNAdb$Fasta_seq)),
sum(!is.na(drosophila_databases$rRNA_ucsc$Fasta_seq)),
sum(!is.na(drosophila_databases$tRFdb$Fasta_seq))
),
  Equal_Fasta_Loci_seqs = c(
sum(drosophila_databases$miRBase$equal_seq,na.rm = T) + sum(drosophila_databases$miRBase$equal_rev_seq,na.rm = T),
sum(drosophila_databases$piRNAdb$equal_seq,na.rm = T) + sum(drosophila_databases$piRNAdb$equal_rev_seq,na.rm = T),
sum(drosophila_databases$GtRNAdb$equal_seq,na.rm = T) + sum(drosophila_databases$GtRNAdb$equal_rev_seq,na.rm = T),
sum(drosophila_databases$rRNA_ucsc$equal_seq,na.rm = T) + sum(drosophila_databases$rRNA_ucsc$equal_rev_seq,na.rm = T),
sum(drosophila_databases$tRFdb$equal_seq,na.rm = T) + sum(drosophila_databases$tRFdb$equal_rev_seq,na.rm = T)
), row.names = names(drosophila_databases)),
'drosophila/smRNA_comparative_table_stats.xlsx',rowNames = T)

write.xlsx(drosophila_databases,'drosophila/smRNA_comparative_table_drosophila.xlsx',rowNames = T)

```

## Pictures

```{r} 

venn.diagram(list(unique(drosophila_databases$miRBase$ID), unique(Drosophila_miRBase_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'drosophila/Venn_drosophila_miRBase_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red','blue'), cat.dist = c(0.05,0.065),
             main = 'Unique ID of miRBase', fill = c('red','blue'))

venn.diagram(list(unique(drosophila_databases$piRNAdb$ID), unique(Drosophila_piRNAdb_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'drosophila/Venn_drosophila_piRNAdb_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red','blue'), cat.dist = c(0.05,0.04),
             main = 'Unique ID of piRNAdb', fill = c('red','blue'))

venn.diagram(list(unique(drosophila_databases$GtRNAdb$ID), unique(Drosophila_GtRNAdb_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'drosophila/Venn_drosophila_GtRNAdb_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red','blue'), cat.dist = c(0.05,0.04),
             main = 'Unique ID of GtRNAdb', fill = c('red','blue'))

ggplot(drosophila_databases$miRBase) + aes(Loci_len,Fasta_len) +
  theme_classic() + geom_point() + geom_abline() +
  ggtitle('Comparison of length of miRBase loci and fasta seqs')
ggsave('drosophila/Comparison of drosophila length of microRNA.tiff')

ggplot(drosophila_databases$tRFdb) + aes(Loci_len,Fasta_len) +
  theme_classic() + geom_point() + geom_abline() +
  ggtitle('Comparison of length of tRFdb loci and fasta seqs')
ggsave('drosophila/Comparison of drosophila length of tRNA-derived.tiff')



```


## Building united database

```{r loci&fa} 
Drosophila_miRBase_precursors_loci_fa = filter(drosophila_databases$miRBase,equal_seq+equal_rev_seq == 1,
                                    str_detect(ID,'dme-mir'))
Drosophila_miRBase_mature_loci_fa = filter(drosophila_databases$miRBase,equal_seq+equal_rev_seq == 1,
                                    str_detect(ID,'dme-miR')) 
Drosophila_piRNAdb_loci_fa = filter(drosophila_databases$piRNAdb,equal_seq+equal_rev_seq == 1)
Drosophila_GtRNAdb_loci_fa = filter(drosophila_databases$GtRNAdb,equal_seq+equal_rev_seq == 1)
Drosophila_rRNA_ucsc_loci_fa = filter(drosophila_databases$rRNA_ucsc,equal_seq+equal_rev_seq == 1)
Drosophila_tRFdb_loci_fa = filter(drosophila_databases$tRFdb,equal_seq+equal_rev_seq == 1)

drosophila_smRNA_DATABASE = rbind.data.frame(
  data.frame(type = 'mature_miRNA',
           ID = Drosophila_miRBase_mature_loci_fa$ID,
           Annotation_consensus = Drosophila_miRBase_mature_loci_fa$Loci,
           Annotation_chr = Drosophila_miRBase_mature_loci_fa$Loci_chr,
           Annotation_start = Drosophila_miRBase_mature_loci_fa$Loci_start,
           Annotation_end = Drosophila_miRBase_mature_loci_fa$Loci_end,
           Sequence_consensus = Drosophila_miRBase_mature_loci_fa$Fasta_seq,
           Length_consensus = Drosophila_miRBase_mature_loci_fa$Loci_len,
           strand = Drosophila_miRBase_mature_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
  data.frame(type = 'precursor_miRNA',
           ID = Drosophila_miRBase_precursors_loci_fa$ID,
           Annotation_consensus = Drosophila_miRBase_precursors_loci_fa$Loci,
           Annotation_chr = Drosophila_miRBase_precursors_loci_fa$Loci_chr,
           Annotation_start = Drosophila_miRBase_precursors_loci_fa$Loci_start,
           Annotation_end = Drosophila_miRBase_precursors_loci_fa$Loci_end,
           Sequence_consensus = Drosophila_miRBase_precursors_loci_fa$Fasta_seq,
           Length_consensus = Drosophila_miRBase_precursors_loci_fa$Loci_len,
           strand = Drosophila_miRBase_precursors_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'piRNA',
           ID = Drosophila_piRNAdb_loci_fa$ID,
           Annotation_consensus = Drosophila_piRNAdb_loci_fa$Loci,
           Annotation_chr = Drosophila_piRNAdb_loci_fa$Loci_chr,
           Annotation_start = Drosophila_piRNAdb_loci_fa$Loci_start,
           Annotation_end = Drosophila_piRNAdb_loci_fa$Loci_end,
           Sequence_consensus = Drosophila_piRNAdb_loci_fa$Fasta_seq,
           Length_consensus = Drosophila_piRNAdb_loci_fa$Loci_len,
           strand = Drosophila_piRNAdb_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'tRNA-derived',
           ID = Drosophila_tRFdb_loci_fa$ID,
           Annotation_consensus = Drosophila_tRFdb_loci_fa$Loci,
           Annotation_chr = Drosophila_tRFdb_loci_fa$Loci_chr,
           Annotation_start = Drosophila_tRFdb_loci_fa$Loci_start,
           Annotation_end = Drosophila_tRFdb_loci_fa$Loci_end,
           Sequence_consensus = Drosophila_tRFdb_loci_fa$Fasta_seq,
           Length_consensus = Drosophila_tRFdb_loci_fa$Loci_len,
           strand = Drosophila_tRFdb_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'tRNA',
           ID = Drosophila_GtRNAdb_loci_fa$ID,
           Annotation_consensus = Drosophila_GtRNAdb_loci_fa$Loci,
           Annotation_chr = Drosophila_GtRNAdb_loci_fa$Loci_chr,
           Annotation_start = Drosophila_GtRNAdb_loci_fa$Loci_start,
           Annotation_end = Drosophila_GtRNAdb_loci_fa$Loci_end,
           Sequence_consensus = Drosophila_GtRNAdb_loci_fa$Fasta_seq,
           Length_consensus = Drosophila_GtRNAdb_loci_fa$Loci_len,
           strand = Drosophila_GtRNAdb_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'rRNA',
           ID = Drosophila_rRNA_ucsc_loci_fa$ID,
           Annotation_consensus = Drosophila_rRNA_ucsc_loci_fa$Loci,
           Annotation_chr = Drosophila_rRNA_ucsc_loci_fa$Loci_chr,
           Annotation_start = Drosophila_rRNA_ucsc_loci_fa$Loci_start,
           Annotation_end = Drosophila_rRNA_ucsc_loci_fa$Loci_end,
           Sequence_consensus = Drosophila_rRNA_ucsc_loci_fa$Fasta_seq,
           Length_consensus = Drosophila_rRNA_ucsc_loci_fa$Loci_len,
           strand = Drosophila_rRNA_ucsc_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
)
)

```

```{r loci only}

Drosophila_miRBase_precursors_loci_only = rbind.data.frame(
  filter(drosophila_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'dme-mir')),
  filter(drosophila_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'dme-ban')),
  filter(drosophila_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'dme-let')))
Drosophila_miRBase_precursors_loci_only = filter(Drosophila_miRBase_precursors_loci_only,
                                            !(str_detect(ID,'p')))
Drosophila_miRBase_precursors_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(Drosophila_miRBase_precursors_loci_only),.combine = 'c') %do% {
    if (Drosophila_miRBase_precursors_loci_only$Loci_strand[i] == '-') {out = Drosophila_miRBase_precursors_loci_only$Loci_seq[i]}
    if (Drosophila_miRBase_precursors_loci_only$Loci_strand[i] == '+') {out = Drosophila_miRBase_precursors_loci_only$Loci_rev_seq[i]}
    out
  }

Drosophila_miRBase_mature_loci_only = unique(rbind.data.frame(
  filter(drosophila_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'dme-miR')),
  filter(drosophila_databases$miRBase, is.na(Fasta_seq),str_detect(ID,'p'))))
Drosophila_miRBase_mature_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(Drosophila_miRBase_mature_loci_only),.combine = 'c') %do% {
    if (Drosophila_miRBase_mature_loci_only$Loci_strand[i] == '-') {out = Drosophila_miRBase_mature_loci_only$Loci_seq[i]}
    if (Drosophila_miRBase_mature_loci_only$Loci_strand[i] == '+') {out = Drosophila_miRBase_mature_loci_only$Loci_rev_seq[i]}
    out
  }
Drosophila_piRNAdb_loci_only = filter(drosophila_databases$piRNAdb,is.na(Fasta_seq)) # is null

Drosophila_GtRNAdb_loci_only = filter(drosophila_databases$GtRNAdb,is.na(Fasta_seq))
Drosophila_GtRNAdb_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(Drosophila_GtRNAdb_loci_only),.combine = 'c') %do% {
    if (Drosophila_GtRNAdb_loci_only$Loci_strand[i] == '-') {out = Drosophila_GtRNAdb_loci_only$Loci_seq[i]}
    if (Drosophila_GtRNAdb_loci_only$Loci_strand[i] == '+') {out = Drosophila_GtRNAdb_loci_only$Loci_rev_seq[i]}
    out
  }
Drosophila_tRFdb_loci_only = filter(drosophila_databases$tRFdb,is.na(Fasta_seq)) # is null
Drosophila_rRNA_ucsc_loci_only = filter(drosophila_databases$rRNA_ucsc,is.na(Fasta_seq)) # is null

drosophila_smRNA_DATABASE = rbind.data.frame(drosophila_smRNA_DATABASE,
  data.frame(type = 'precursor_miRNA',
             ID = Drosophila_miRBase_precursors_loci_only$ID,
             Annotation_consensus = Drosophila_miRBase_precursors_loci_only$Loci,
             Annotation_chr = Drosophila_miRBase_precursors_loci_only$Loci_chr,
             Annotation_start = Drosophila_miRBase_precursors_loci_only$Loci_start,
             Annotation_end = Drosophila_miRBase_precursors_loci_only$Loci_end,
             Sequence_consensus = Drosophila_miRBase_precursors_loci_only$Fasta_seq,
             Length_consensus = Drosophila_miRBase_precursors_loci_only$Loci_len,
             strand = Drosophila_miRBase_precursors_loci_only$Loci_strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ),
  data.frame(type = 'mature_miRNA',
             ID = Drosophila_miRBase_mature_loci_only$ID,
             Annotation_consensus = Drosophila_miRBase_mature_loci_only$Loci,
             Annotation_chr = Drosophila_miRBase_mature_loci_only$Loci_chr,
             Annotation_start = Drosophila_miRBase_mature_loci_only$Loci_start,
             Annotation_end = Drosophila_miRBase_mature_loci_only$Loci_end,
             Sequence_consensus = Drosophila_miRBase_mature_loci_only$Fasta_seq,
             Length_consensus = Drosophila_miRBase_mature_loci_only$Loci_len,
             strand = Drosophila_miRBase_mature_loci_only$Loci_strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ),
  data.frame(type = 'tRNA',
             ID = Drosophila_GtRNAdb_loci_only$ID,
             Annotation_consensus = Drosophila_GtRNAdb_loci_only$Loci,
             Annotation_chr = Drosophila_GtRNAdb_loci_only$Loci_chr,
             Annotation_start = Drosophila_GtRNAdb_loci_only$Loci_start,
             Annotation_end = Drosophila_GtRNAdb_loci_only$Loci_end,
             Sequence_consensus = Drosophila_GtRNAdb_loci_only$Fasta_seq,
             Length_consensus = Drosophila_GtRNAdb_loci_only$Loci_len,
             strand = Drosophila_GtRNAdb_loci_only$Loci_strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ))

```

```{r fa only}

Drosophila_miRBase_precursors_fa_only = filter(Drosophila_miRBase_sam_stats, !(ID %in% drosophila_databases$miRBase$ID), str_detect(ID,'dme-mir'), Fasta_chr != '*') # is null
Drosophila_miRBase_mature_fa_only = filter(Drosophila_miRBase_sam_stats, !(ID %in% drosophila_databases$miRBase$ID),Fasta_chr != '*', str_detect(ID,'rno-miR')) # is null
Drosophila_piRNAdb_fa_only = filter(Drosophila_piRNAdb_sam_stats,!(ID %in% drosophila_databases$piRNAdb$ID),
                         Fasta_chr != '*') # is null
Drosophila_GtRNAdb_fa_only = filter(Drosophila_GtRNAdb_sam_stats,!(ID %in% drosophila_databases$GtRNAdb$ID),
                         Fasta_chr != '*') # is null
Drosophila_rRNA_ucsc_fa_only = filter(Drosophila_rRNA_sam_stats,!(ID %in% drosophila_databases$rRNA_ucsc$ID),
                         Fasta_chr != '*')

drosophila_smRNA_DATABASE = rbind.data.frame(drosophila_smRNA_DATABASE,
                                  data.frame(type = 'rRNA',
                                             ID = Drosophila_rRNA_ucsc_fa_only$ID,
        Annotation_consensus = str_c(Drosophila_rRNA_ucsc_fa_only$Fasta_chr,':',
                                     Drosophila_rRNA_ucsc_fa_only$Fasta_start,'-',
                                     as.character(Drosophila_rRNA_ucsc_fa_only$Fasta_start+Drosophila_rRNA_ucsc_fa_only$Fasta_len)),
                                             Annotation_chr = Drosophila_rRNA_ucsc_fa_only$Fasta_chr,
                                             Annotation_start = Drosophila_rRNA_ucsc_fa_only$Fasta_start,
                                             Annotation_end = Drosophila_rRNA_ucsc_fa_only$Fasta_start+
          Drosophila_rRNA_ucsc_fa_only$Fasta_len,
                                             Sequence_consensus = Drosophila_rRNA_ucsc_fa_only$Fasta_seq,
                                             Length_consensus = Drosophila_rRNA_ucsc_fa_only$Fasta_len,
        strand = '+',
        strand_from_annotation = F,
                                             Source = 'fasta',
                                             stringsAsFactors = F
                                  ))

```

```{r count H-dist & suff-pref intersect}

Drosophila_miRBase_precursors_loci_not_fa = rbind.data.frame(
  filter(drosophila_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'dme-mir')),
  filter(drosophila_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'dme-let')),
  filter(drosophila_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'dme-ban')))
Drosophila_miRBase_precursors_loci_not_fa = filter(Drosophila_miRBase_precursors_loci_not_fa,
                                              !(str_detect(ID,'p')))
Drosophila_miRBase_mature_loci_not_fa = unique(rbind.data.frame(
  filter(drosophila_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'dme-miR')),
  filter(drosophila_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'p'))))
Drosophila_GtRNAdb_loci_not_fa = filter(drosophila_databases$GtRNAdb, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq)) # is null
Drosophila_rRNA_ucsc_loci_not_fa = filter(drosophila_databases$rRNA_ucsc, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq)) # is null
Drosophila_piRNAdb_loci_not_fa = filter(drosophila_databases$piRNAdb, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq))
Drosophila_tRFdb_loci_not_fa = filter(drosophila_databases$tRFdb, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq))

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

Drosophila_loci_not_fa = foreach(d = list(Drosophila_miRBase_precursors_loci_not_fa,
Drosophila_miRBase_mature_loci_not_fa,
Drosophila_tRFdb_loci_not_fa)) %do% {
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
names(Drosophila_loci_not_fa) = c('miRBase_precursors','miRBase_mature','tRFdb')

Drosophila_loci_not_fa$piRNAdb = Drosophila_piRNAdb_loci_not_fa
d = Drosophila_loci_not_fa$piRNAdb
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
Drosophila_loci_not_fa$piRNAdb$ID = str_c(Drosophila_loci_not_fa$piRNAdb$Loci_seq,'_',Drosophila_loci_not_fa$piRNAdb$Fasta_seq)
Drosophila_loci_not_fa$piRNAdb = 
  cbind.data.frame(Drosophila_loci_not_fa$piRNAdb,
                   foreach(i = 1:nrow(Drosophila_loci_not_fa$piRNAdb),.combine = 'rbind') %do% {
                     if (i %% 10000 == 0) {print(i)}
                     k = Drosophila_loci_not_fa$piRNAdb$ID[i]
                     d[which(d$ID == k),4:6]
                   }
                   )

```

```{r fix loci not fa}
Drosophila_loci_not_fa$miRBase_precursors$type = 'precursor_miRNA'
Drosophila_loci_not_fa$miRBase_mature$type = 'mature_miRNA'
Drosophila_loci_not_fa$tRFdb$type = 'tRNA-derived'
Drosophila_loci_not_fa$piRNAdb$type = 'piRNA'

Drosophila_loci_not_fa$piRNAdb$ID = Drosophila_piRNAdb_loci_not_fa$ID


for (i in 1:3) {
  d = Drosophila_loci_not_fa[[i]]
  drosophila_smRNA_DATABASE = rbind.data.frame(drosophila_smRNA_DATABASE, 
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
drosophila_smRNA_DATABASE = filter(drosophila_smRNA_DATABASE,Source != 'trash')

d = Drosophila_loci_not_fa[[4]]
d$diff_len = d$Fasta_len - d$Loci_len
d$fix_p2 = d$Fasta_len - d$Loci_len == as.numeric(d$p3)
d$p2[which(d$fix_p2 == T)] = '0'

wide_annot = filter(d,p1 %in% c('loci_forward_in_fa'))
Annotation_start = wide_annot$Loci_start-as.numeric(wide_annot$p2)
Annotation_end = wide_annot$Loci_end+as.numeric(wide_annot$p3)
Sequence_consensus = wide_annot$Fasta_seq
Length_consensus = wide_annot$Fasta_len
Source = 'wide_annotation&fasta'
drosophila_smRNA_DATABASE = rbind.data.frame(drosophila_smRNA_DATABASE, 
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
drosophila_smRNA_DATABASE = rbind.data.frame(drosophila_smRNA_DATABASE, 
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
drosophila_smRNA_DATABASE = rbind.data.frame(drosophila_smRNA_DATABASE, 
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
drosophila_smRNA_DATABASE = rbind.data.frame(drosophila_smRNA_DATABASE, 
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

write.table(str_c(filter(drosophila_smRNA_DATABASE,type == 'mature_miRNA')$Annotation_chr,filter(drosophila_smRNA_DATABASE,type == 'mature_miRNA')$Annotation_start,filter(drosophila_smRNA_DATABASE,type == 'mature_miRNA')$Annotation_end,filter(drosophila_smRNA_DATABASE,type == 'mature_miRNA')$ID,sep = '\t'),'drosophila/drosophila_smRNA_DATABASE_miRNA.bed',col.names = F, row.names = F)
write.table(str_c(filter(drosophila_smRNA_DATABASE,type == 'piRNA')$Annotation_chr,filter(drosophila_smRNA_DATABASE,type == 'piRNA')$Annotation_start,filter(drosophila_smRNA_DATABASE,type == 'piRNA')$Annotation_end,filter(drosophila_smRNA_DATABASE,type == 'piRNA')$ID,sep = '\t'),'drosophila/drosophila_smRNA_DATABASE_piRNA.bed',col.names = F, row.names = F)
write.table(str_c(filter(drosophila_smRNA_DATABASE,type == 'tRNA')$Annotation_chr,filter(drosophila_smRNA_DATABASE,type == 'tRNA')$Annotation_start,filter(drosophila_smRNA_DATABASE,type == 'tRNA')$Annotation_end,filter(drosophila_smRNA_DATABASE,type == 'tRNA')$ID,sep = '\t'),'drosophila/drosophila_smRNA_DATABASE_tRNA.bed',col.names = F, row.names = F)
write.table(str_c(filter(drosophila_smRNA_DATABASE,type == 'rRNA')$Annotation_chr,filter(drosophila_smRNA_DATABASE,type == 'rRNA')$Annotation_start,filter(drosophila_smRNA_DATABASE,type == 'rRNA')$Annotation_end,filter(drosophila_smRNA_DATABASE,type == 'rRNA')$ID,sep = '\t'),'drosophila/drosophila_smRNA_DATABASE_rRNA.bed',col.names = F, row.names = F)

```

```{r filter by intersept}

drosophila_inters = rbind.data.frame(
  unique(read.table('drosophila/miRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('drosophila/miRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('drosophila/miRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
  unique(read.table('drosophila/piRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('drosophila/piRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('drosophila/piRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('drosophila/rRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('drosophila/rRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('drosophila/rRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('drosophila/tRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('drosophila/tRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('drosophila/tRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('drosophila/miRNA_inter_out.txt',stringsAsFactors = F)[,1:4]),
unique(read.table('drosophila/piRNA_inter_out.txt',stringsAsFactors = F)[,1:4]),
unique(read.table('drosophila/tRNA_inter_out.txt',stringsAsFactors = F)[,1:4]),
unique(read.table('drosophila/rRNA_inter_out.txt',stringsAsFactors = F)[,1:4])
)
drosophila_inters = unique(drosophila_inters)
drosophila_inters_vector = str_c(drosophila_inters$V1,drosophila_inters$V2,drosophila_inters$V3,drosophila_inters$V4)



drosophila_smRNA_DATABASE$filtered = foreach(i = 1:nrow(drosophila_smRNA_DATABASE),.combine = 'c',.packages = libs) %dopar% {
  str_c(drosophila_smRNA_DATABASE$Annotation_chr[i],
        drosophila_smRNA_DATABASE$Annotation_start[i],
        drosophila_smRNA_DATABASE$Annotation_end[i],drosophila_smRNA_DATABASE$ID[i]) %in% drosophila_inters_vector
}

write.xlsx(filter(drosophila_smRNA_DATABASE, type %in% c('mature_miRNA', 'piRNA', 'tRNA', 'rRNA'),filtered == F, filtered_in == F),'drosophila/drosophila_smRNA_DATABASE_mipitr_filtered.xlsx')

write.table(str_c('>',filter(drosophila_smRNA_DATABASE,type %in% c('mature_miRNA','piRNA','rRNA','tRNA'), filtered == F)$ID,'\n',filter(drosophila_smRNA_DATABASE,type %in% c('mature_miRNA','piRNA','rRNA','tRNA'), filtered == F)$Sequence_consensus,collapse = '\n'),'drosophila/drosophila_smRNA_DATABASE_filtered.fasta',row.names = F,col.names = F)

```

```{r Unique rRNA}

drosophila_smRNA_DATABASE$ID[which(drosophila_smRNA_DATABASE$type == 'rRNA')] = str_remove(filter(drosophila_smRNA_DATABASE, type == 'rRNA')$ID,'_[[:digit:]]+')
for (i in 1:length(unique(filter(drosophila_smRNA_DATABASE,type == 'rRNA')[,7]))){
  drosophila_smRNA_DATABASE$ID[which(drosophila_smRNA_DATABASE$Sequence_consensus == unique(filter(drosophila_smRNA_DATABASE,type == 'rRNA')[,7])[i] )] = str_c(drosophila_smRNA_DATABASE$ID[which(drosophila_smRNA_DATABASE$Sequence_consensus == unique(filter(drosophila_smRNA_DATABASE,type == 'rRNA')[,7])[i] )],'_seq_',i)
}

```

```{r write DATABASE}
Drosophila_Source_df = drosophila_smRNA_DATABASE %>% group_by(type) %>% 
  summarise(n = table(Source),
            Source = names(table(Source)))
write.xlsx(pivot_wider(Drosophila_Source_df, names_from = type, values_from = n),'drosophila/drosophila_smRNA_DATABASE_stats.xlsx')
write.xlsx(drosophila_smRNA_DATABASE,'drosophila/drosophila_smRNA_DATABASE.xlsx')
write.table(str_c('>',drosophila_smRNA_DATABASE$ID[which(drosophila_smRNA_DATABASE$type %in% c('mature_miRNA','piRNA','rRNA','tRNA'))],'\n',drosophila_smRNA_DATABASE$Sequence_consensus[which(drosophila_smRNA_DATABASE$type %in% c('mature_miRNA','piRNA','rRNA','tRNA'))],collapse = '\n'),'drosophila/smRNA_DATABASE.fasta',row.names = F,col.names = F)

drosophila_inters_in = rbind.data.frame(
  unique(read.table('drosophila/miRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('drosophila/miRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('drosophila/miRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
  unique(read.table('drosophila/piRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('drosophila/piRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('drosophila/piRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('drosophila/rRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('drosophila/rRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('drosophila/rRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('drosophila/tRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('drosophila/tRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('drosophila/tRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]))

drosophila_inters_in = unique(drosophila_inters_in)
drosophila_inters_in_vector = str_c(drosophila_inters_in$V1,drosophila_inters_in$V2,drosophila_inters_in$V3,drosophila_inters_in$V4)


drosophila_smRNA_DATABASE$filtered_in = foreach(i = 1:nrow(drosophila_smRNA_DATABASE),.combine = 'c',.packages = libs) %dopar% {
  str_c(drosophila_smRNA_DATABASE$Annotation_chr[i],
        drosophila_smRNA_DATABASE$Annotation_start[i],
        drosophila_smRNA_DATABASE$Annotation_end[i],drosophila_smRNA_DATABASE$ID[i]) %in% drosophila_inters_in_vector
}

for (i in unique(drosophila_smRNA_DATABASE$type)) {
  write.xlsx(filter(drosophila_smRNA_DATABASE,filtered_in == F,type == i),str_c('drosophila/drosophila_smRNA_DATABASE_',i,'.xlsx'))
}
```
