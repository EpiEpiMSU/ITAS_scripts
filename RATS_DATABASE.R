```{r}
library(stringr)
library(tidyverse)
library(ggpubr)
library(foreach)
library(openxlsx)
library(doParallel)
library(VennDiagram)
libs = c('stringr','tidyverse','openxlsx','foreach')

cl <- makeCluster(6)
registerDoParallel(cl) 
```

---

## Downdoad annotation files

```{r}

Rat_miRBase_annotation = read.table('rat/rno.gff3',stringsAsFactors = F)
Rat_miRBase_annotation$ID = foreach(i = 1:nrow(Rat_miRBase_annotation),.combine = 'c') %do% {
  str_remove(strsplit(Rat_miRBase_annotation$V9[i],';')[[1]][3],'Name=')
}
Rat_piRNAdb_annotation = read.table('rat/pirnadb.v1_7_6.rn6.gff3',stringsAsFactors = F)
Rat_GtRNAdb_annotation = read.table('rat/rn7-tRNAs.bed',stringsAsFactors = F)
Rat_rRNA_ucsc_annotation = read.table('rat/rRNA_ucsc_rat.gff',stringsAsFactors = F)

write.table(Rat_rRNA_ucsc_annotation[,c(1,4,5,13,7)],'rat/rat_rRNA.bed', col.names = F,
            row.names = F, sep = '\t')

```


---
## Build bed-files for micro, piwi RNA -> liftover -> getfasta

```{r}
write.table(Rat_miRBase_annotation[,c(1,4,5,10,7)],'rat/rat_miRBase_rn6.bed', col.names = F,
            row.names = F, sep = '\t')
write.table(Rat_piRNAdb_annotation[1:1000000,c(1,4,5,9,7)],'rat/rat_piRNAdb_rn6_1.bed', col.names = F,
            row.names = F, sep = '\t')
write.table(Rat_piRNAdb_annotation[1000001:2000000,c(1,4,5,9,7)],'rat/rat_piRNAdb_rn6_2.bed', col.names = F,
            row.names = F, sep = '\t')
write.table(Rat_piRNAdb_annotation[2000001:nrow(Rat_piRNAdb_annotation),c(1,4,5,9,7)],'rat/rat_piRNAdb_rn6_3.bed', col.names = F,
            row.names = F, sep = '\t')

Rat_miRBase_rn7_annotation = read.table('rat/rat_miRBase_rn7.bed',stringsAsFactors = F)
Rat_piRNAdb_rn7_annotation = read.table('rat/rat_piRNAdb_rn7.bed',stringsAsFactors = F)

Rat_miRBase_getfasta = read.table('rat/rat_miRBase_getfasta.tsv',stringsAsFactors = F)
colnames(Rat_miRBase_getfasta) = c('Loci','Loci_seq')
Rat_GtRNAdb_getfasta = read.table('rat/rat_GtRNAdb_getfasta.tsv',stringsAsFactors = F)
colnames(Rat_GtRNAdb_getfasta) = c('Loci','Loci_seq')

Rat_piRNAdb_getfasta_fa = read.table('rat/rat_piRNAdb_getfasta.fasta',stringsAsFactors = F)
Rat_piRNAdb_getfasta = 
  data.frame(Loci = Rat_piRNAdb_getfasta_fa$V1[seq(1,length(Rat_piRNAdb_getfasta_fa$V1), by = 2)],
             Loci_seq = Rat_piRNAdb_getfasta_fa$V1[seq(2,length(Rat_piRNAdb_getfasta_fa$V1), by = 2)],stringsAsFactors = F)

Rat_piRNAdb_getfasta$Loci = str_remove_all(Rat_piRNAdb_getfasta$Loci,'>')

Rat_rRNA_ucsc_getfasta = read.table('rat/rat_rRNA_ucsc_getfasta.tsv',stringsAsFactors = F)
colnames(Rat_rRNA_ucsc_getfasta) = c('Loci','Loci_seq')

```

---
## Build loci databases

```{r}

rats_databases = list()
rats_databases$miRBase = data.frame(
  ID = Rat_miRBase_rn7_annotation$V4,
  Loci = str_c(Rat_miRBase_rn7_annotation$V1,':',Rat_miRBase_rn7_annotation$V2,'-',Rat_miRBase_rn7_annotation$V3),
  Loci_chr = Rat_miRBase_rn7_annotation$V1,
  Loci_start = Rat_miRBase_rn7_annotation$V2,
  Loci_end = Rat_miRBase_rn7_annotation$V3,
  Loci_len = as.numeric(Rat_miRBase_rn7_annotation$V3) - as.numeric(Rat_miRBase_rn7_annotation$V2),
  Loci_strand = Rat_miRBase_rn7_annotation$V5,
  stringsAsFactors = F
)
rats_databases$piRNAdb = data.frame(
  ID = Rat_piRNAdb_rn7_annotation$V4,
  Loci = str_c(Rat_piRNAdb_rn7_annotation$V1,':',Rat_piRNAdb_rn7_annotation$V2,'-',Rat_piRNAdb_rn7_annotation$V3),
  Loci_chr = Rat_piRNAdb_rn7_annotation$V1,
  Loci_start = Rat_piRNAdb_rn7_annotation$V2,
  Loci_end = Rat_piRNAdb_rn7_annotation$V3,
  Loci_len = as.numeric(Rat_piRNAdb_rn7_annotation$V3) - as.numeric(Rat_piRNAdb_rn7_annotation$V2),
  Loci_strand = Rat_piRNAdb_rn7_annotation$V5,
  stringsAsFactors = F
)
rats_databases$GtRNAdb = data.frame(
  ID = Rat_GtRNAdb_annotation$V4,
  Loci = str_c(Rat_GtRNAdb_annotation$V1,':',Rat_GtRNAdb_annotation$V2,'-',Rat_GtRNAdb_annotation$V3),
  Loci_chr = Rat_GtRNAdb_annotation$V1,
  Loci_start = Rat_GtRNAdb_annotation$V2,
  Loci_end = Rat_GtRNAdb_annotation$V3,
  Loci_len = as.numeric(Rat_GtRNAdb_annotation$V3) - as.numeric(Rat_GtRNAdb_annotation$V2),
  Loci_strand = Rat_GtRNAdb_annotation$V6,
  stringsAsFactors = F
)
rats_databases$rRNA_ucsc = data.frame(
  ID = Rat_rRNA_ucsc_annotation$V10,
  Loci = str_c(Rat_rRNA_ucsc_annotation$V1,':',Rat_rRNA_ucsc_annotation$V4,'-',Rat_rRNA_ucsc_annotation$V5),
  Loci_chr = Rat_rRNA_ucsc_annotation$V1,
  Loci_start = Rat_rRNA_ucsc_annotation$V4,
  Loci_end = Rat_rRNA_ucsc_annotation$V5,
  Loci_len = as.numeric(Rat_rRNA_ucsc_annotation$V5) - as.numeric(Rat_rRNA_ucsc_annotation$V4),
  Loci_strand = Rat_rRNA_ucsc_annotation$V7,
  stringsAsFactors = F
)

```

---
## Add getfasta sequences to loci -> filtering -> reverse compliment sequences

```{r}

rats_databases$miRBase$Loci_seq = foreach(l = rats_databases$miRBase$Loci,.combine = 'c') %dopar% {
    Rat_miRBase_getfasta[which(l == Rat_miRBase_getfasta$Loci)[1],'Loci_seq']}
rats_databases$GtRNAdb$Loci_seq = foreach(l = rats_databases$GtRNAdb$Loci,.combine = 'c') %dopar% {
    Rat_GtRNAdb_getfasta[which(l == Rat_GtRNAdb_getfasta$Loci)[1],'Loci_seq']}
rats_databases$GtRNAdb = filter(rats_databases$GtRNAdb,!is.na(Loci_seq))

Rat_piRNAdb_getfasta = arrange(Rat_piRNAdb_getfasta,Loci)
rats_databases$piRNAdb = filter(rats_databases$piRNAdb,Loci %in% Rat_piRNAdb_getfasta$Loci)
rats_databases$piRNAdb = arrange(rats_databases$piRNAdb,Loci)
rats_databases$piRNAdb$Loci_seq = Rat_piRNAdb_getfasta$Loci_seq

rats_databases$rRNA_ucsc$Loci_seq = foreach(l = rats_databases$rRNA_ucsc$Loci,.combine = 'c') %do% {
    Rat_rRNA_ucsc_getfasta[which(l == Rat_rRNA_ucsc_getfasta$Loci)[1],'Loci_seq']}
rats_databases$rRNA_ucsc = filter(rats_databases$rRNA_ucsc,!is.na(Loci_seq))

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

for (x in 1:4) {
  print(x)
  rats_databases[[x]] = mutate(rats_databases[[x]],Loci_rev_seq = rev_comp(Loci_seq))
}

Rat_piRNAdb_Loci_seq = data.frame(Loci_seq = unique(rats_databases$piRNAdb$Loci_seq),
                                  stringsAsFactors = F)

Rat_piRNAdb_Loci_seq$Loci_rev_seq = foreach(i = Rat_piRNAdb_Loci_seq$Loci_seq,.combine = 'c') %dopar% {
  rev_comp(i)
}
Rat_piRNAdb_Loci_seq = arrange(Rat_piRNAdb_Loci_seq, Loci_seq)
rats_databases$piRNAdb = arrange(rats_databases$piRNAdb, Loci_seq)

Rat_piRNAdb_Loci_rev_seq = foreach(i = 1:nrow(Rat_piRNAdb_Loci_seq),.combine = 'c') %do% {
  if (i %% 1000 == 0) {print(i)}
  rep(Rat_piRNAdb_Loci_seq$Loci_rev_seq[i],
      nrow(filter(rats_databases$piRNAdb, Loci_seq == Rat_piRNAdb_Loci_seq$Loci_seq[i])))
}

rats_databases$piRNAdb$Loci_rev_seq = Rat_piRNAdb_Loci_rev_seq

```

---

## Download fasta and sam files
```{r}
Rat_miRBase_fasta = rbind.data.frame(
  read.table('rat/miRBase_hairpin.fasta.tsv',stringsAsFactors = F,sep = '\t'),
  read.table('rat/miRBase_mature.fasta.tsv',stringsAsFactors = F,sep = '\t'))
colnames(Rat_miRBase_fasta) = c('ID','Fasta_seq')
Rat_miRBase_fasta$ID = foreach(i = 1:nrow(Rat_miRBase_fasta),.combine = 'c') %do% {
  strsplit(Rat_miRBase_fasta$ID[i],' ')[[1]][1]
}

Rat_GtRNAdb_fasta = read.table('rat/rn7-tRNAs.fa.tsv',stringsAsFactors = F,sep = '\t')
Rat_GtRNAdb_fasta$Fasta_seq = str_c(Rat_GtRNAdb_fasta$V2,Rat_GtRNAdb_fasta$V3)
Rat_GtRNAdb_fasta = Rat_GtRNAdb_fasta[,-(2:3)]
colnames(Rat_GtRNAdb_fasta) = c('ID','Fasta_seq')
Rat_GtRNAdb_fasta$ID = str_remove(Rat_GtRNAdb_fasta$ID,'Rattus_norvegicus_')
Rat_GtRNAdb_fasta$ID = t(data.frame(strsplit(Rat_GtRNAdb_fasta$ID,' '),stringsAsFactors = F))[,1]

Rat_piRNAdb_fasta = read.table('rat/piRNAdb.rno.v1_7_6.fa.tsv',stringsAsFactors = F)
colnames(Rat_piRNAdb_fasta) = c('ID','Fasta_seq')

Rat_rRNA_ucsc_fasta = read.table('rat/rRNA_ucsc_rat.fasta.tsv',stringsAsFactors = F)
Rat_rRNA_ucsc_fasta = Rat_rRNA_ucsc_fasta[,c(1,2,7)]
colnames(Rat_rRNA_ucsc_fasta) = c('ID','Loci','Fasta_seq')
Rat_rRNA_ucsc_fasta$ID = str_remove_all(Rat_rRNA_ucsc_fasta$ID,'rn7_rmsk_')
Rat_rRNA_ucsc_fasta$Loci = str_remove_all(Rat_rRNA_ucsc_fasta$Loci,'range=')
Rat_rRNA_ucsc_fasta = filter(Rat_rRNA_ucsc_fasta,Loci %in% rats_databases$rRNA_ucsc$Loci)

Rat_miRBase_sam_stats = rbind.data.frame(
  read.table('rat/miRBase_mature_sam_stats.txt',stringsAsFactors = F),
  read.table('rat/miRBase_hairpin_sam_stats.txt',stringsAsFactors = F))
colnames(Rat_miRBase_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Rat_miRBase_sam_stats$Fasta_start = Rat_miRBase_sam_stats$Fasta_start - 1
Rat_miRBase_sam_stats$Fasta_mapped = foreach(i = Rat_miRBase_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

Rat_GtRNAdb_sam_stats = read.table('rat/tRNA_sam_stats.txt',stringsAsFactors = F)
colnames(Rat_GtRNAdb_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Rat_GtRNAdb_sam_stats$Fasta_start = Rat_GtRNAdb_sam_stats$Fasta_start - 1
Rat_GtRNAdb_sam_stats$ID = foreach(i = Rat_GtRNAdb_sam_stats$ID,.combine = 'c') %do% {
  strsplit(i,'_')[[1]][3]
}
Rat_GtRNAdb_sam_stats$Fasta_mapped = foreach(i = Rat_GtRNAdb_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]}

Rat_piRNAdb_sam_stats = read.table('rat/piRNAdb_sam_stats.txt',stringsAsFactors = F)
colnames(Rat_piRNAdb_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                           'Fasta_len','Fasta_mapped')
Rat_piRNAdb_sam_stats$Fasta_start = Rat_piRNAdb_sam_stats$Fasta_start - 1
Rat_piRNAdb_sam_stats$Fasta_mapped = foreach(i = Rat_piRNAdb_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

Rat_rRNA_sam_stats = read.table('rat/rRNA_sam_stats.txt',stringsAsFactors = F)
colnames(Rat_rRNA_sam_stats) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
Rat_rRNA_sam_stats$Fasta_start = Rat_rRNA_sam_stats$Fasta_start - 1
Rat_rRNA_sam_stats$ID = str_remove_all(Rat_rRNA_sam_stats$ID,'rn7_rmsk_')
Rat_rRNA_sam_stats$Fasta_mapped = foreach(i = Rat_rRNA_sam_stats$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]}

```

---
## Add fasta sequences to loci -> check equal seq

```{r}
rats_databases$miRBase$Fasta_seq = NA
foreach(i = 1:nrow(rats_databases$miRBase),.combine = 'c') %do% {
  rats_databases$miRBase$Fasta_seq[i] = 
    Rat_miRBase_fasta[which(Rat_miRBase_fasta$ID == rats_databases$miRBase$ID[i])[1],2]
}
  
rats_databases$GtRNAdb$Fasta_seq = NA
foreach(i = 1:nrow(rats_databases$GtRNAdb),.combine = 'c') %do% {
  rats_databases$GtRNAdb$Fasta_seq[i] = 
    Rat_GtRNAdb_fasta[which(Rat_GtRNAdb_fasta$ID == rats_databases$GtRNAdb$ID[i])[1],2]
}

rats_databases$piRNAdb$Fasta_seq = NA
rats_databases$piRNAdb$Fasta_seq = foreach(i = 1:nrow(rats_databases$piRNAdb),.combine = 'c') %do% {
  if (i %% 100000 == 0) {print(i)}
  Rat_piRNAdb_fasta[which(Rat_piRNAdb_fasta$ID == rats_databases$piRNAdb$ID[i])[1],2]
}

rats_databases$rRNA_ucsc$Fasta_seq = Rat_rRNA_ucsc_fasta$Fasta_seq

rats_databases$miRBase$Fasta_len = foreach(i = strsplit(rats_databases$miRBase$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
rats_databases$GtRNAdb$Fasta_len = foreach(i = strsplit(rats_databases$GtRNAdb$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
rats_databases$piRNAdb$Fasta_len = foreach(i = strsplit(rats_databases$piRNAdb$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}
rats_databases$rRNA_ucsc$Fasta_len = foreach(i = strsplit(rats_databases$rRNA_ucsc$Fasta_seq,''),.combine = 'c') %do% {
  length(i)
}

rats_databases$miRBase$equal_seq = 
  rats_databases$miRBase$Fasta_seq == rats_databases$miRBase$Loci_seq
rats_databases$GtRNAdb$equal_seq = 
  rats_databases$GtRNAdb$Fasta_seq == rats_databases$GtRNAdb$Loci_seq
rats_databases$piRNAdb$equal_seq = 
  rats_databases$piRNAdb$Fasta_seq == rats_databases$piRNAdb$Loci_seq
rats_databases$rRNA_ucsc$equal_seq = 
  rats_databases$rRNA_ucsc$Fasta_seq == rats_databases$rRNA_ucsc$Loci_seq

rats_databases$miRBase$equal_rev_seq = 
  rats_databases$miRBase$Fasta_seq == rats_databases$miRBase$Loci_rev_seq
rats_databases$GtRNAdb$equal_rev_seq = 
  rats_databases$GtRNAdb$Fasta_seq == rats_databases$GtRNAdb$Loci_rev_seq
rats_databases$piRNAdb$equal_rev_seq = 
  rats_databases$piRNAdb$Fasta_seq == rats_databases$piRNAdb$Loci_rev_seq
rats_databases$rRNA_ucsc$equal_rev_seq = 
  rats_databases$rRNA_ucsc$Fasta_seq == rats_databases$rRNA_ucsc$Loci_rev_seq

```

---

## Sam stats add

```{r}

Rat_sam_stats = list(Rat_miRBase_sam_stats,Rat_piRNAdb_sam_stats,Rat_GtRNAdb_sam_stats,
                     Rat_rRNA_sam_stats)
names(Rat_sam_stats) = names(rats_databases)

for (x in c(1,3)) {
  print(x)
  d = rats_databases[[x]]
  
  t = foreach(pos = 1:nrow(d),.combine = 'rbind') %do% {
    if (pos %% 1000 == 0){print(pos)}
      #if (pos > nrow(d)) {break}
  #    d$mapped2loci[pos] = 
    data.frame(
       mapped2loci =  mean(str_c(d$Loci_chr[pos],d$Loci_start[pos]) == str_c(
          Rat_sam_stats[[x]]$Fasta_chr[which(Rat_sam_stats[[x]]$ID == d$ID[pos])],
          Rat_sam_stats[[x]]$Fasta_start[which(Rat_sam_stats[[x]]$ID == d$ID[pos])]
          )) > 0,
      #d$Fasta_mapped[pos] = 
        Fasta_mapped = Rat_sam_stats[[x]]$Fasta_mapped[which(Rat_sam_stats[[x]]$ID == d$ID[pos])][1],stringsAsFactors = F
    )
  }
  
  d = cbind.data.frame(d,t)
  rats_databases[[x]] = d
}

```


---

## Stats of united database

```{r}
table(filter(smRNA_DATABASE,type == 'tRNA-derived')$Annotation_chr)
write.xlsx(data.frame(table(smRNA_DATABASE$Annotation_chr)),'table_chr_human.xlsx')

write.xlsx(data.frame(table(rbind.data.frame(rats_databases$miRBase,
                                       rats_databases$piRNAdb, rats_databases$GtRNAdb)$Loci_chr)),'table_chr_rats_loci.xlsx')

length(unique(rats_databases$miRBase$ID))
length(unique(rats_databases$piRNAdb$ID))
length(unique(rats_databases$GtRNAdb$ID))
length(unique(rats_databases$rRNA_ucsc$ID))

sum(!is.na(rats_databases$miRBase$Fasta_seq))
sum(!is.na(rats_databases$piRNAdb$Fasta_seq))
sum(!is.na(rats_databases$GtRNAdb$Fasta_seq))
sum(!is.na(rats_databases$rRNA_ucsc$Fasta_seq))

sum(rats_databases$piRNAdb$equal_seq,na.rm = T) + sum(rats_databases$piRNAdb$equal_rev_seq,na.rm = T)
sum(rats_databases$rRNA_ucsc$equal_seq,na.rm = T) + sum(rats_databases$rRNA_ucsc$equal_rev_seq,na.rm = T)

for (i in 1:4) {
  print(names(rats_databases[i]))
  print(sum(rats_databases[[i]]$mapped2loci,na.rm = T))
  print(sum(rats_databases[[i]]$Fasta_mapped == '1',na.rm = T))
}

write.xlsx(rats_databases,'rat/smRNA_databases_rats.xlsx')

```

---

## Pictures

```{r} 

venn.diagram(list(unique(rats_databases$miRBase$ID), unique(Rat_miRBase_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'rat/Venn_rat_miRBase_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('orange','blue'), cat.dist = c(-0.1,0.02),
             main = 'Unique ID of miRBase', fill = c('red','blue'))

venn.diagram(list(unique(rats_databases$piRNAdb$ID), unique(Rat_piRNAdb_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'rat/Venn_rats_piRNAdb_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red3','blue'), cat.dist = c(0.04,0.03),
             main = 'Unique ID of piRNAdb', fill = c('red','blue'))

venn.diagram(list(unique(rats_databases$GtRNAdb$ID), unique(Rat_GtRNAdb_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'rat/Venn_rat_GtRNAdb_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red3','blue'), cat.dist = c(0.04,0.03),
             main = 'Unique ID of GtRNAdb', fill = c('red','blue'))

ggplot(rats_databases$rRNA_ucsc) + aes(Loci_len,Fasta_len, color = Fasta_mapped) +
  theme_classic() + geom_point() + geom_abline() +
  ggtitle('Comparison of length of rRNA ucsc loci and fasta seqs')
ggsave('rat/Comparison of rat length of rRNA.tiff')

ggplot(rats_databases$rRNA_ucsc) + aes(Loci_len - Fasta_len) +
  theme_classic() + geom_histogram(binwidth = 1) +
  ggtitle('Difference of length of rRNA ucsc loci and fasta')
ggsave('rat/Difference of length of rat rRNA.tiff',height = 10)

```

---

## Building united database

```{r loci&fa} 
Rat_miRBase_precursors_loci_fa = filter(rats_databases$miRBase,equal_seq+equal_rev_seq == 1,
                                    str_detect(ID,'rno-mir')) #is null
Rat_miRBase_mature_loci_fa = filter(rats_databases$miRBase,equal_seq+equal_rev_seq == 1,
                                    str_detect(ID,'rno-miR')) #is null
Rat_piRNAdb_loci_fa = filter(rats_databases$piRNAdb,equal_seq+equal_rev_seq == 1)
Rat_GtRNAdb_loci_fa = filter(rats_databases$GtRNAdb,equal_seq+equal_rev_seq == 1)
Rat_rRNA_ucsc_loci_fa = filter(rats_databases$rRNA_ucsc,equal_seq+equal_rev_seq == 1)

rats_smRNA_DATABASE = rbind.data.frame(
data.frame(type = 'piRNA',
           ID = Rat_piRNAdb_loci_fa$ID,
           Annotation_consensus = Rat_piRNAdb_loci_fa$Loci,
           Annotation_chr = Rat_piRNAdb_loci_fa$Loci_chr,
           Annotation_start = Rat_piRNAdb_loci_fa$Loci_start,
           Annotation_end = Rat_piRNAdb_loci_fa$Loci_end,
           Sequence_consensus = Rat_piRNAdb_loci_fa$Fasta_seq,
           Length_consensus = Rat_piRNAdb_loci_fa$Loci_len,
           mapped2loci = Rat_piRNAdb_loci_fa$mapped2loci,
           strand = Rat_piRNAdb_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'tRNA',
           ID = Rat_GtRNAdb_loci_fa$ID,
           Annotation_consensus = Rat_GtRNAdb_loci_fa$Loci,
           Annotation_chr = Rat_GtRNAdb_loci_fa$Loci_chr,
           Annotation_start = Rat_GtRNAdb_loci_fa$Loci_start,
           Annotation_end = Rat_GtRNAdb_loci_fa$Loci_end,
           Sequence_consensus = Rat_GtRNAdb_loci_fa$Fasta_seq,
           Length_consensus = Rat_GtRNAdb_loci_fa$Loci_len,
           mapped2loci = Rat_GtRNAdb_loci_fa$mapped2loci,
           strand = Rat_GtRNAdb_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'rRNA',
           ID = Rat_rRNA_ucsc_loci_fa$ID,
           Annotation_consensus = Rat_rRNA_ucsc_loci_fa$Loci,
           Annotation_chr = Rat_rRNA_ucsc_loci_fa$Loci_chr,
           Annotation_start = Rat_rRNA_ucsc_loci_fa$Loci_start,
           Annotation_end = Rat_rRNA_ucsc_loci_fa$Loci_end,
           Sequence_consensus = Rat_rRNA_ucsc_loci_fa$Fasta_seq,
           Length_consensus = Rat_rRNA_ucsc_loci_fa$Loci_len,
           mapped2loci = Rat_rRNA_ucsc_loci_fa$mapped2loci,
           strand = Rat_rRNA_ucsc_loci_fa$Loci_strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
)
)

```

```{r loci only}

Rat_miRBase_precursors_loci_only = filter(rats_databases$miRBase, is.na(Fasta_seq),
                                    str_detect(ID,'rno-mir'))
Rat_miRBase_precursors_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(Rat_miRBase_precursors_loci_only),.combine = 'c') %do% {
    if (Rat_miRBase_precursors_loci_only$Loci_strand[i] == '-') {out = Rat_miRBase_precursors_loci_only$Loci_seq[i]}
    if (Rat_miRBase_precursors_loci_only$Loci_strand[i] == '+') {out = Rat_miRBase_precursors_loci_only$Loci_rev_seq[i]}
    out
  }

Rat_miRBase_mature_loci_only = filter(rats_databases$miRBase,is.na(Fasta_seq),
                                str_detect(ID,'rno-miR'))
Rat_miRBase_mature_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(Rat_miRBase_mature_loci_only),.combine = 'c') %do% {
    if (Rat_miRBase_mature_loci_only$Loci_strand[i] == '-') {out = Rat_miRBase_mature_loci_only$Loci_seq[i]}
    if (Rat_miRBase_mature_loci_only$Loci_strand[i] == '+') {out = Rat_miRBase_mature_loci_only$Loci_rev_seq[i]}
    out
  }
Rat_piRNAdb_loci_only = filter(rats_databases$piRNAdb,is.na(Fasta_seq)) # is null

Rat_GtRNAdb_loci_only = filter(rats_databases$GtRNAdb,is.na(Fasta_seq))
Rat_GtRNAdb_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(Rat_GtRNAdb_loci_only),.combine = 'c') %do% {
    if (Rat_GtRNAdb_loci_only$Loci_strand[i] == '-') {out = Rat_GtRNAdb_loci_only$Loci_seq[i]}
    if (Rat_GtRNAdb_loci_only$Loci_strand[i] == '+') {out = Rat_GtRNAdb_loci_only$Loci_rev_seq[i]}
    out
  }

Rat_rRNA_ucsc_loci_only = filter(rats_databases$rRNA_ucsc,is.na(Fasta_seq)) # is null

rats_smRNA_DATABASE = rbind.data.frame(rats_smRNA_DATABASE,
  data.frame(type = 'precursor_miRNA',
             ID = Rat_miRBase_precursors_loci_only$ID,
             Annotation_consensus = Rat_miRBase_precursors_loci_only$Loci,
             Annotation_chr = Rat_miRBase_precursors_loci_only$Loci_chr,
             Annotation_start = Rat_miRBase_precursors_loci_only$Loci_start,
             Annotation_end = Rat_miRBase_precursors_loci_only$Loci_end,
             Sequence_consensus = Rat_miRBase_precursors_loci_only$Fasta_seq,
             Length_consensus = Rat_miRBase_precursors_loci_only$Loci_len,
             mapped2loci = Rat_miRBase_precursors_loci_only$mapped2loci,
             strand = Rat_miRBase_precursors_loci_only$Loci_strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ),
  data.frame(type = 'mature_miRNA',
             ID = Rat_miRBase_mature_loci_only$ID,
             Annotation_consensus = Rat_miRBase_mature_loci_only$Loci,
             Annotation_chr = Rat_miRBase_mature_loci_only$Loci_chr,
             Annotation_start = Rat_miRBase_mature_loci_only$Loci_start,
             Annotation_end = Rat_miRBase_mature_loci_only$Loci_end,
             Sequence_consensus = Rat_miRBase_mature_loci_only$Fasta_seq,
             Length_consensus = Rat_miRBase_mature_loci_only$Loci_len,
             mapped2loci = Rat_miRBase_mature_loci_only$mapped2loci,
             strand = Rat_miRBase_mature_loci_only$Loci_strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ),
  data.frame(type = 'tRNA',
             ID = Rat_GtRNAdb_loci_only$ID,
             Annotation_consensus = Rat_GtRNAdb_loci_only$Loci,
             Annotation_chr = Rat_GtRNAdb_loci_only$Loci_chr,
             Annotation_start = Rat_GtRNAdb_loci_only$Loci_start,
             Annotation_end = Rat_GtRNAdb_loci_only$Loci_end,
             Sequence_consensus = Rat_GtRNAdb_loci_only$Fasta_seq,
             Length_consensus = Rat_GtRNAdb_loci_only$Loci_len,
             mapped2loci = Rat_GtRNAdb_loci_only$mapped2loci,
             strand = Rat_GtRNAdb_loci_only$Loci_strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ))

```

```{r fa only}

Rat_miRBase_precursors_fa_only = filter(Rat_miRBase_sam_stats, !(ID %in% rats_databases$miRBase$ID), str_detect(ID,'rno-mir'), Fasta_chr != '*')
Rat_miRBase_mature_fa_only = filter(Rat_miRBase_sam_stats, !(ID %in% rats_databases$miRBase$ID),Fasta_chr != '*', str_detect(ID,'rno-miR'))
Rat_piRNAdb_fa_only = filter(Rat_piRNAdb_sam_stats,!(ID %in% rats_databases$piRNAdb$ID),
                         Fasta_chr != '*')
Rat_GtRNAdb_fa_only = filter(Rat_GtRNAdb_sam_stats,!(ID %in% rats_databases$GtRNAdb$ID),
                         Fasta_chr != '*') # is null
Rat_rRNA_ucsc_fa_only = filter(Rat_rRNA_sam_stats,!(ID %in% rats_databases$rRNA_ucsc$ID),
                         Fasta_chr != '*') # is null

rats_smRNA_DATABASE = rbind.data.frame(rats_smRNA_DATABASE,
                                  data.frame(type = 'precursor_miRNA',
                                             ID = Rat_miRBase_precursors_fa_only$ID,
        Annotation_consensus = str_c(Rat_miRBase_precursors_fa_only$Fasta_chr,':',
                                     Rat_miRBase_precursors_fa_only$Fasta_start,'-',
                                     as.character(Rat_miRBase_precursors_fa_only$Fasta_start+Rat_miRBase_precursors_fa_only$Fasta_len)),
                                             Annotation_chr = Rat_miRBase_precursors_fa_only$Fasta_chr,
                                             Annotation_start = Rat_miRBase_precursors_fa_only$Fasta_start,
                                             Annotation_end = Rat_miRBase_precursors_fa_only$Fasta_start+
          Rat_miRBase_precursors_fa_only$Fasta_len,
                                             Sequence_consensus = Rat_miRBase_precursors_fa_only$Fasta_seq,
                                             Length_consensus = Rat_miRBase_precursors_fa_only$Fasta_len,
                                             mapped2loci = NA,
        strand = '+',
        strand_from_annotation = F,
                                             Source = 'fasta',
                                             stringsAsFactors = F
                                  ),
                                  data.frame(type = 'mature_miRNA',
                                             ID = Rat_miRBase_mature_fa_only$ID,
                                             Annotation_consensus = str_c(Rat_miRBase_mature_fa_only$Fasta_chr,':', Rat_miRBase_mature_fa_only$Fasta_start,'-',
as.character(Rat_miRBase_mature_fa_only$Fasta_start + Rat_miRBase_mature_fa_only$Fasta_len)),
                                             Annotation_chr = Rat_miRBase_mature_fa_only$Fasta_chr,
                                             Annotation_start = Rat_miRBase_mature_fa_only$Fasta_start,
                                             Annotation_end = Rat_miRBase_mature_fa_only$Fasta_start+
                                               Rat_miRBase_mature_fa_only$Fasta_len,
                                             Sequence_consensus = Rat_miRBase_mature_fa_only$Fasta_seq,
                                             Length_consensus = Rat_miRBase_mature_fa_only$Fasta_len,
                                             mapped2loci = NA,
                                             strand = '+',
                                             strand_from_annotation = F,
                                             Source = 'fasta',
                                             stringsAsFactors = F
                                  ),
                                  data.frame(type = 'piRNA',
                                             ID = Rat_piRNAdb_fa_only$ID,
                                             Annotation_consensus = str_c(Rat_piRNAdb_fa_only$Fasta_chr,':',Rat_piRNAdb_fa_only$Fasta_start,'-',
                                                                          as.character(Rat_piRNAdb_fa_only$Fasta_start+Rat_piRNAdb_fa_only$Fasta_len)),
                                             Annotation_chr = Rat_piRNAdb_fa_only$Fasta_chr,
                                             Annotation_start = Rat_piRNAdb_fa_only$Fasta_start,
                                             Annotation_end = Rat_piRNAdb_fa_only$Fasta_start+Rat_piRNAdb_fa_only$Fasta_len,
                                             Sequence_consensus = Rat_piRNAdb_fa_only$Fasta_seq,
                                             Length_consensus = Rat_piRNAdb_fa_only$Fasta_len,
                                             mapped2loci = NA,
                                             strand = '+',
                                             strand_from_annotation = F,
                                             Source = 'fasta',
                                             stringsAsFactors = F
                                  ))

```

```{r count H-dist & suff-pref intersect}

Rat_miRBase_precursors_loci_not_fa = filter(rats_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'rno-mir'))
Rat_miRBase_mature_loci_not_fa = filter(rats_databases$miRBase, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq), str_detect(ID,'rno-miR'))
Rat_GtRNAdb_loci_not_fa = filter(rats_databases$GtRNAdb, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq)) # is null
Rat_rRNA_ucsc_loci_not_fa = filter(rats_databases$rRNA_ucsc, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq))
Rat_piRNAdb_loci_not_fa = filter(rats_databases$piRNAdb, equal_seq+equal_rev_seq == 0, !is.na(Fasta_seq))

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

Rat_loci_not_fa = foreach(d = list(Rat_miRBase_precursors_loci_not_fa,
Rat_miRBase_mature_loci_not_fa,
Rat_rRNA_ucsc_loci_not_fa)) %do% {
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
names(Rat_loci_not_fa) = c('miRBase_precursors','miRBase_mature','rRNA_ucsc')

Rat_loci_not_fa$piRNAdb = Rat_piRNAdb_loci_not_fa
d = Rat_loci_not_fa$piRNAdb
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
Rat_loci_not_fa$piRNAdb$ID = str_c(Rat_loci_not_fa$piRNAdb$Loci_seq,'_',Rat_loci_not_fa$piRNAdb$Fasta_seq)
Rat_loci_not_fa$piRNAdb = 
  cbind.data.frame(Rat_loci_not_fa$piRNAdb,
                   foreach(i = 1:nrow(Rat_loci_not_fa$piRNAdb),.combine = 'rbind') %do% {
                     if (i %% 10000 == 0) {print(i)}
                     k = Rat_loci_not_fa$piRNAdb$ID[i]
                     d[which(d$ID == k),4:6]
                   }
                   )

```

```{r H dist & inter pictures}
for (i in 1:3) {
  d = Rat_loci_not_fa[[i]]
  d$source = d$p1
  d$source[which(!(d$source %in% c('loci_forward_in_fa','loci_reverse_in_fa')))] = 'others'
  
  ggplot(d) + aes(source) + theme_classic() + geom_histogram(stat="count") +
    ggtitle(str_c('Source of difference loci & fa of ',names(Rat_loci_not_fa)[i]), subtitle = str_c('n = ',nrow(d)))
  ggsave(str_c('rat/Source of difference loci & fa of ',names(Rat_loci_not_fa)[i],'.tiff'))
  
  ggplot(d) + aes(Fasta_len - Loci_len) + theme_classic() + geom_histogram(binwidth = 0.5) +
    ggtitle(str_c('Difference fasta & loci length of ',names(Rat_loci_not_fa)[i]), subtitle = str_c('n = ',nrow(d)))
  ggsave(str_c('rat/Difference fasta & loci length of ',names(Rat_loci_not_fa)[i],'.tiff'))
  df = filter(d,source == 'others')
  
  ggplot(df) + aes(as.numeric(p1),p2) + theme_classic() + geom_point() + 
    scale_x_continuous(name = 'Hamming distance') +
    scale_y_continuous(name = 'Suff-pref intersection length') +
    ggtitle(str_c('Distribution of Hamming distance and intersection length of ',names(Rat_loci_not_fa)[i]), subtitle = str_c('n = ',nrow(df)))
  ggsave(str_c('rat/Distribution of Hamming distance and intersection length of ',names(Rat_loci_not_fa)[i],'.tiff'))
}


```

```{r fix loci not fa}
Rat_loci_not_fa$miRBase_precursors$type = 'precursor_miRNA'
Rat_loci_not_fa$miRBase_mature$type = 'mature_miRNA'
Rat_loci_not_fa$rRNA_ucsc$type = 'rRNA'
Rat_loci_not_fa$piRNAdb$type = 'piRNA'

Rat_loci_not_fa$piRNAdb$ID = Rat_piRNAdb_loci_not_fa$ID


for (i in 1:3) {
  d = Rat_loci_not_fa[[i]]
  rats_smRNA_DATABASE = rbind.data.frame(rats_smRNA_DATABASE, 
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
               Length_consensus = Length_consensus, mapped2loci = d$mapped2loci[k],
               strand = d$Loci_strand[k], strand_from_annotation = T,
               Source = Source,stringsAsFactors = F
               )
  })
}
rats_smRNA_DATABASE = filter(rats_smRNA_DATABASE,Source != 'trash')

d = Rat_loci_not_fa[[4]]
d$diff_len = d$Fasta_len - d$Loci_len
d$fix_p2 = d$Fasta_len - d$Loci_len == as.numeric(d$p3)
d$p2[which(d$fix_p1 == T)] = '0'

wide_annot = filter(d,p1 %in% c('loci_forward_in_fa'))
Annotation_start = wide_annot$Loci_start-as.numeric(wide_annot$p2)
Annotation_end = wide_annot$Loci_end+as.numeric(wide_annot$p3)
Sequence_consensus = wide_annot$Fasta_seq
Length_consensus = wide_annot$Fasta_len
Source = 'wide_annotation&fasta'
rats_smRNA_DATABASE = rbind.data.frame(rats_smRNA_DATABASE, 
      data.frame(type = wide_annot$type,ID = wide_annot$ID,
               Annotation_consensus = str_c(wide_annot$Loci_chr,':',Annotation_start,'-', Annotation_end),
               Annotation_chr = wide_annot$Loci_chr, Annotation_start = Annotation_start,
               Annotation_end = Annotation_end, Sequence_consensus = Sequence_consensus,
               Length_consensus = Length_consensus, mapped2loci = wide_annot$mapped2loci,
               strand = wide_annot$Loci_strand, strand_from_annotation = T,
               Source = Source,stringsAsFactors = F
               ))
wide_annot = filter(d,p1 %in% c('loci_reverse_in_fa'))
Annotation_start = wide_annot$Loci_start-as.numeric(wide_annot$p3)
Annotation_end = wide_annot$Loci_end+as.numeric(wide_annot$p2)
Sequence_consensus = wide_annot$Fasta_seq
Length_consensus = wide_annot$Fasta_len
Source = 'wide_annotation&fasta'
rats_smRNA_DATABASE = rbind.data.frame(rats_smRNA_DATABASE, 
      data.frame(type = wide_annot$type,ID = wide_annot$ID,
               Annotation_consensus = str_c(wide_annot$Loci_chr,':',Annotation_start,'-', Annotation_end),
               Annotation_chr = wide_annot$Loci_chr, Annotation_start = Annotation_start,
               Annotation_end = Annotation_end, Sequence_consensus = Sequence_consensus,
               Length_consensus = Length_consensus, mapped2loci = wide_annot$mapped2loci,
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
rats_smRNA_DATABASE = rbind.data.frame(rats_smRNA_DATABASE, 
      data.frame(type = wide_fasta$type,ID = wide_fasta$ID,
               Annotation_consensus = str_c(wide_fasta$Loci_chr,':',Annotation_start,'-', Annotation_end),
               Annotation_chr = wide_fasta$Loci_chr, Annotation_start = Annotation_start,
               Annotation_end = Annotation_end, Sequence_consensus = Sequence_consensus,
               Length_consensus = Length_consensus, mapped2loci = wide_fasta$mapped2loci,
               strand = wide_fasta$Loci_strand, strand_from_annotation = T,
               Source = Source,stringsAsFactors = F
               ))

wide_annot_fasta = filter(d,Fasta_len - as.numeric(p2) == 1, diff_len == 0)
d = wide_annot_fasta
rats_smRNA_DATABASE = rbind.data.frame(rats_smRNA_DATABASE, 
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
               Length_consensus = Length_consensus, mapped2loci = d$mapped2loci[k],
               strand = d$Loci_strand[k], strand_from_annotation = T,
               Source = Source,stringsAsFactors = F
               )
  })

```



----

## Intersection in mature miRNA, piRNA, tRNA, rRNA base

```{r write bed -> inters}

Rat_mipitrrib_data = unique(filter(rats_smRNA_DATABASE,type %in% c('tRNA','piRNA','mature_miRNA','rRNA')))

write.table(cbind.data.frame(Rat_mipitrrib_data[,4:6],1:nrow(Rat_mipitrrib_data)),sep = '\t', row.names = F, col.names = F,'rat/Rat_mipitrrib.bed')

Rat_mipitrrib_inter = read.table('rat/Rat_mi_pi_tr_rib_inter.bed',stringsAsFactors = F)
inter_record = table(Rat_mipitrrib_inter$V4)
inter_record = inter_record[which(inter_record == 1)]

Rat_mipitrrib_data_no_inters = Rat_mipitrrib_data[as.numeric(names(inter_record)),]
Rat_mipitrrib_data_no_inters_check = table(foreach(i = 1:nrow(Rat_mipitrrib_data_no_inters),.combine = 'c') %do% {
  if (i %% 100000 == 0) {print(i)}
  str_c(Rat_mipitrrib_data_no_inters$Annotation_chr[i],'_',
        (Rat_mipitrrib_data_no_inters$Annotation_start[i]+1):Rat_mipitrrib_data_no_inters$Annotation_end[i])
})
unique(Rat_mipitrrib_data_no_inters_check)

write.xlsx(pivot_wider(Rat_mipitrrib_data_no_inters %>% group_by(type) %>% 
  summarise(n = table(Source),
            Source = names(table(Source))), names_from = type, values_from = n),'rat/rats_miR_pi_rib_tRNA_DATABASE_filter_stats.xlsx')

#fix rRNA####
for (i in 1:length(unique(filter(Rat_mipitrrib_data_no_inters,type == 'rRNA')[,7]))){
  Rat_mipitrrib_data_no_inters$ID[which(Rat_mipitrrib_data_no_inters$Sequence_consensus == unique(filter(Rat_mipitrrib_data_no_inters,type == 'rRNA')[,7])[i] )] = str_c(Rat_mipitrrib_data_no_inters$ID[which(Rat_mipitrrib_data_no_inters$Sequence_consensus == unique(filter(Rat_mipitrrib_data_no_inters,type == 'rRNA')[,7])[i] )],'_seq_',i)
}
  
write.xlsx(Rat_mipitrrib_data_no_inters,'rat/Rat_smRNA_DATABASE_filter_miR_pi_rib_tRNA.xlsx')

```


```{r intersection with %}
Rat_positions_smRNA_DATABASE = foreach(i = 1:nrow(rats_smRNA_DATABASE)) %do% {
  if (i %% 100000 == 0) {print(i)}
  str_c(rats_smRNA_DATABASE$Annotation_chr[i],'_',
        rats_smRNA_DATABASE$Annotation_start[i]:rats_smRNA_DATABASE$Annotation_end[i])
}
names(Rat_positions_smRNA_DATABASE) = foreach(i = 1:length(Rat_positions_smRNA_DATABASE),
                                          .combine = 'c')%do%{
  if (i %% 100000 == 0) {print(i)}
  Rat_positions_smRNA_DATABASE[[i]][1]
                                          }

Rat_positions_mipitrrib_DATABASE = Rat_positions_smRNA_DATABASE[
  str_c(Rat_mipitrrib_data$Annotation_chr,'_',Rat_mipitrrib_data$Annotation_start)]
d_chr_rat = foreach(chr = unique(Rat_mipitrrib_data$Annotation_chr),.combine = 'rbind') %do% {
  print(chr)
  d = filter(Rat_mipitrrib_data,Annotation_chr == chr)
  p = Rat_positions_mipitrrib_DATABASE[
    str_detect(names(Rat_positions_mipitrrib_DATABASE),str_c(chr,'_'))]
  positions = table(unlist(p))
  positions = positions[which(positions >1)]
  n = names(positions)
  d$interseptions = foreach(i = p,.combine = 'c') %dopar% {
    sum(i %in% n) 
  }/(d$Length_consensus+1)
  d
}

```

---

## Count repeats in smRNA_DATABASE

```{r}

table_rats_smRNA_DATABASE_ID = table(rats_smRNA_DATABASE$ID)

rats_smRNA_DATABASE$repeats = foreach(i = rats_smRNA_DATABASE$ID[1:1000],.combine = 'c') %do% {
  table_rats_smRNA_DATABASE_ID[i]
}
for (i in names(table(rats_smRNA_DATABASE$type))) {
  ggplot(filter(rats_smRNA_DATABASE, type == i, repeats > 1)) + aes(y = repeats) + theme_classic() +
    geom_boxplot() +
    ggtitle('Number of loci for one ID', subtitle = str_c('type = ',i))
  ggsave(str_c('rat/Number of loci for one ID for rat ',i,'.tiff'))
}

ggplot(data.frame(num = as.vector(table(rats_smRNA_DATABASE$repeats)),
                  similarity = as.numeric(names(table(rats_smRNA_DATABASE$repeats))),
                  n = as.vector(table(rats_smRNA_DATABASE$repeats))/
                    as.numeric(names(table(rats_smRNA_DATABASE$repeats))))[2:50,]) +
  theme_classic() + aes(similarity, n) + geom_point() +
  ggtitle('Distribution of number of repets')
ggsave('rat/Distribution of number of repets rats.tiff') 
```



## Interceptions

```{r write bed}
rats_smRNA_DATABASE = read.xlsx('Rat/rats_smRNA_DATABASE.xlsx')

write.table(str_c(filter(rats_smRNA_DATABASE,type == 'mature_miRNA')$Annotation_chr,filter(rats_smRNA_DATABASE,type == 'mature_miRNA')$Annotation_start,filter(rats_smRNA_DATABASE,type == 'mature_miRNA')$Annotation_end,filter(rats_smRNA_DATABASE,type == 'mature_miRNA')$ID,sep = '\t'),'rat/rats_smRNA_DATABASE_miRNA.bed',col.names = F, row.names = F)
write.table(str_c(filter(rats_smRNA_DATABASE,type == 'piRNA')$Annotation_chr,filter(rats_smRNA_DATABASE,type == 'piRNA')$Annotation_start,filter(rats_smRNA_DATABASE,type == 'piRNA')$Annotation_end,filter(rats_smRNA_DATABASE,type == 'piRNA')$ID,sep = '\t'),'rat/rats_smRNA_DATABASE_piRNA.bed',col.names = F, row.names = F)
write.table(str_c(filter(rats_smRNA_DATABASE,type == 'tRNA')$Annotation_chr,filter(rats_smRNA_DATABASE,type == 'tRNA')$Annotation_start,filter(rats_smRNA_DATABASE,type == 'tRNA')$Annotation_end,filter(rats_smRNA_DATABASE,type == 'tRNA')$ID,sep = '\t'),'rat/rats_smRNA_DATABASE_tRNA.bed',col.names = F, row.names = F)
write.table(str_c(filter(rats_smRNA_DATABASE,type == 'rRNA')$Annotation_chr,filter(rats_smRNA_DATABASE,type == 'rRNA')$Annotation_start,filter(rats_smRNA_DATABASE,type == 'rRNA')$Annotation_end,filter(rats_smRNA_DATABASE,type == 'rRNA')$ID,sep = '\t'),'rat/rats_smRNA_DATABASE_rRNA.bed',col.names = F, row.names = F)

```

```{r}
rats_inters_in = rbind.data.frame(
  unique(read.table('rat/miRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('rat/miRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('rat/miRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
  unique(read.table('rat/piRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('rat/piRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('rat/piRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('rat/rRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('rat/rRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('rat/rRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]),
unique(read.table('rat/tRNA_inter_in.txt',stringsAsFactors = F)[which(read.table('rat/tRNA_inter_in.txt',stringsAsFactors = F)[,4] != 
  read.table('rat/tRNA_inter_in.txt',stringsAsFactors = F)[,8]),1:4]))
rats_inters_in = unique(rats_inters_in)
rats_inters_in_vector = str_c(rats_inters_in$V1,rats_inters_in$V2,rats_inters_in$V3,rats_inters_in$V4)

rats_smRNA_DATABASE$filtered_in = F
rats_smRNA_DATABASE$filtered_in[which(str_c(rats_smRNA_DATABASE$Annotation_chr,
        rats_smRNA_DATABASE$Annotation_start,
        rats_smRNA_DATABASE$Annotation_end,rats_smRNA_DATABASE$ID) %in% rats_inters_in_vector)] = T

for (i in unique(rats_smRNA_DATABASE$type)) {
  write.xlsx(filter(rats_smRNA_DATABASE, filtered_in == F, type == i), str_c('rat/rats_smRNA_DATABASE_',i,'.xlsx'))
}


```

```{r write DATABASE}
Rat_Source_df = rats_smRNA_DATABASE %>% group_by(type) %>% 
  summarise(n = table(Source),
            Source = names(table(Source)))
write.xlsx(pivot_wider(Rat_Source_df, names_from = type, values_from = n),'rat/rats_smRNA_DATABASE_stats.xlsx')
write.xlsx(rats_smRNA_DATABASE,'rat/rats_smRNA_DATABASE.xlsx')
```

## write fasta smRNA_database

```{r}
Rat_mipitrrib_data_no_inters_unique = unique(Rat_mipitrrib_data_no_inters[,c(2,7)])
for (i in names(table(Rat_mipitrrib_data_no_inters_unique$ID)[
  which(table(Rat_mipitrrib_data_no_inters_unique$ID) >1)]) ) {
  k = nrow(filter(Rat_mipitrrib_data_no_inters_unique,ID == i))
  Rat_mipitrrib_data_no_inters_unique[
    Rat_mipitrrib_data_no_inters_unique$ID == i,]$ID = str_c(
      Rat_mipitrrib_data_no_inters_unique[
    Rat_mipitrrib_data_no_inters_unique$ID == i,]$ID,'_seq_',1:k
    )
  }

Rat_mipitrrib_data_no_inters_unique_CCA_tRNA = 
  filter(Rat_mipitrrib_data_no_inters_unique,startsWith(ID,'tR'))
Rat_mipitrrib_data_no_inters_unique_CCA_tRNA$ID = str_c(
  Rat_mipitrrib_data_no_inters_unique_CCA_tRNA$ID,'-CCA')
Rat_mipitrrib_data_no_inters_unique_CCA_tRNA$Sequence_consensus = str_c(
  Rat_mipitrrib_data_no_inters_unique_CCA_tRNA$Sequence_consensus,'CCA')

Rat_mipitrrib_data_no_inters_unique = rbind.data.frame(Rat_mipitrrib_data_no_inters_unique,
Rat_mipitrrib_data_no_inters_unique_CCA_tRNA)

Rat_mipitrrib_data_no_inters_fasta = foreach(i = 1:nrow(Rat_mipitrrib_data_no_inters_unique),.combine = 'rbind') %do% {
  d = rbind.data.frame(str_c('>',Rat_mipitrrib_data_no_inters_unique$ID[i]),
                       Rat_mipitrrib_data_no_inters_unique$Sequence_consensus[i])
  colnames(d) = 'f'
  d
                                         }

write.table(Rat_mipitrrib_data_no_inters_fasta,
            'rat/Rat_smRNA_DATABASE_filter_miR_pi_rib_tRNA_tRNA-CCA.fasta',
            col.names = F,row.names = F)

for (i in unique(rats_smRNA_DATABASE$type)) {
  write.xlsx(filter(rats_smRNA_DATABASE,type == i),str_c('rat/rats_smRNA_DATABASE_',i,'.xlsx'))
}
```

