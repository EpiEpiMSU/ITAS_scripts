library(stringr)
library(tidyverse)
library(ggpubr)
library(foreach)
library(openxlsx)
library(doParallel)
library(VennDiagram)
library(UpSetR)
library(reshape)

cl <- makeCluster(6)
registerDoParallel(cl) 

#MINTbase -> bed####
MINTbase = read.table('MINTbase.tsv',stringsAsFactors = F)
MINTbase_out = data.frame(
  chr = str_c('chr',MINTbase$V7),
  start = MINTbase$V9,
  end = MINTbase$V10,
  name = str_c(MINTbase$V1,MINTbase$V2,MINTbase$V5,MINTbase$V6,sep = '_'),
  stringsAsFactors = F
)
for (i in 1:nrow(MINTbase_out)) {
  if (is.na(as.numeric(MINTbase_out$start[i]))) {
    MINTbase_out$start[i] = as.integer(strsplit(MINTbase_out$start[i],'-')[[1]][1])
  }
  if (is.na(as.numeric(MINTbase_out$end[i]))) {
    MINTbase_out$end[i] = as.integer(strsplit(MINTbase_out$end[i],'\\+')[[1]][1])
  }
  else {
    MINTbase_out$start[i] = as.integer(MINTbase_out$start[i])
    MINTbase_out$end[i] = as.integer(MINTbase_out$end[i])
  }
  
}
write.table(MINTbase_out,'MINTbase.bed',row.names = F,sep = '\t',col.names = F)

#Table for databases loci####
MINTbase = read.table('MINTBase_hg38.bed')
databases = list()
databases$MINTbase = data.frame(
  ID = MINTbase$V4,
  Loci = str_c(MINTbase$V1,':',MINTbase$V2,'-',MINTbase$V3),
  Loci_chr = MINTbase$V1,
  Loci_start = MINTbase$V2,
  Loci_end = MINTbase$V3,
  Loci_len = as.numeric(MINTbase$V3) - as.numeric(MINTbase$V2),
  stringsAsFactors = F
)
miRBase = read.table('miRBase.bed',stringsAsFactors = F)
databases$miRBase = data.frame(
  ID = foreach(i = strsplit(miRBase$V10,';'),.combine = 'c') %do% {
    str_remove(i[3],'Name=')
  },
  Loci = str_c(miRBase$V1,':',miRBase$V2,'-',miRBase$V3),
  Loci_chr = miRBase$V1,
  Loci_start = miRBase$V2,
  Loci_end = miRBase$V3,
  Loci_len = as.numeric(miRBase$V3) - as.numeric(miRBase$V2),
  stringsAsFactors = F
)
databases$piRNAdb = data.frame(
  ID = piRNAdb_bedcov$V11,
  Loci = str_c(piRNAdb_bedcov$V1,':',piRNAdb_bedcov$V2,'-',piRNAdb_bedcov$V3),
  Loci_chr = piRNAdb_bedcov$V1,
  Loci_start = piRNAdb_bedcov$V2,
  Loci_end = piRNAdb_bedcov$V3,
  Loci_len = as.numeric(piRNAdb_bedcov$V3) - as.numeric(piRNAdb_bedcov$V2),
  stringsAsFactors = F
)
databases$GtRNAdb = data.frame(
  ID = GtRNAdb_bedcov$V4,
  Loci = str_c(GtRNAdb_bedcov$V1,':',GtRNAdb_bedcov$V2,'-',GtRNAdb_bedcov$V3),
  Loci_chr = GtRNAdb_bedcov$V1,
  Loci_start = GtRNAdb_bedcov$V2,
  Loci_end = GtRNAdb_bedcov$V3,
  Loci_len = as.numeric(GtRNAdb_bedcov$V3) - as.numeric(GtRNAdb_bedcov$V2),
  stringsAsFactors = F
)
#getfasta add####
MINTbase_getfasta = read.table('MINTBase_hg38_getfasta.tsv',stringsAsFactors = F)
colnames(MINTbase_getfasta) = c('Loci','Loci_seq')
GtRNAdb_getfasta = read.table('hg38-tRNAs_getfasta.tsv',stringsAsFactors = F)
colnames(GtRNAdb_getfasta) = c('Loci','Loci_seq')
miRBase_getfasta = read.table('miRBase_getfasta.tsv',stringsAsFactors = F)
colnames(miRBase_getfasta) = c('Loci','Loci_seq')
pirnadb_getfasta_fa = read.table('pirnadb_getfasta.fasta',stringsAsFactors = F)
pirnadb_getfasta = data.frame(Loci = vector(length = nrow(pirnadb_getfasta_fa)/2),
                              Loci_seq = vector(length = nrow(pirnadb_getfasta_fa)/2))
foreach(i = 1:(nrow(pirnadb_getfasta_fa)/2),.combine = 'rbind') %do% {
  if (i %% 50000 == 0) {print(i)}
  pirnadb_getfasta$Loci[i] = str_remove(pirnadb_getfasta_fa[(i*2)-1,1],'>')
  pirnadb_getfasta$Loci_seq[i] = pirnadb_getfasta_fa[(i*2),1]
}
getfasta = list(MINTbase_getfasta,miRBase_getfasta,pirnadb_getfasta,GtRNAdb_getfasta)
databases[["miRBase"]]$Loci_seq = foreach(l = databases[["miRBase"]]$Loci,.combine = 'c') %do% {
    miRBase_getfasta[which(l == miRBase_getfasta$Loci)[1],'Loci_seq']}
databases[["MINTbase"]]$Loci_seq = foreach(l = databases[["MINTbase"]]$Loci,.combine = 'c') %do% {
  MINTbase_getfasta[which(l == MINTbase_getfasta$Loci)[1],'Loci_seq']}
databases[["GtRNAdb"]]$Loci_seq = foreach(l = databases[["GtRNAdb"]]$Loci,.combine = 'c') %do% {
  GtRNAdb_getfasta[which(l == GtRNAdb_getfasta$Loci)[1],'Loci_seq']}
databases[["piRNAdb"]]$Loci_seq = ''
foreach(l = 1:nrow(databases[["piRNAdb"]]),.combine = 'c') %do% {
  if (l %% 50000 == 0){print(l)}
  databases[["piRNAdb"]]$Loci_seq[l] = 
    pirnadb_getfasta[which(databases[["piRNAdb"]]$Loci[l] == 
                           pirnadb_getfasta$Loci)[1],'Loci_seq']}

#Fasta seq add####
MINTbase_data = read.table('MINTbase.tsv',stringsAsFactors = F)

MINTbase_fasta = data.frame(ID_seq = MINTbase_data$V2,
                            Fasta_seq = MINTbase_data$V11,
                            stringsAsFactors = F)
MINTbase_ID_seq = t(data.frame(strsplit(databases$MINTbase$ID,'_'),stringsAsFactors = F))[,2]
which(!(MINTbase_fasta$ID_seq %in% MINTbase_ID_seq))
databases$MINTbase$Fasta_seq = NA
foreach(i = 1:nrow(databases$MINTbase),.combine = 'c') %do% {
  databases$MINTbase$Fasta_seq[i] = 
    MINTbase_fasta[which(MINTbase_fasta$ID_seq == MINTbase_ID_seq[i])[1],2]
}

GtRNAdb_fasta = read.table('hg38-tRNAs.fa.tsv',stringsAsFactors = F,sep = '\t')
colnames(GtRNAdb_fasta) = c('ID','Fasta_seq')
GtRNAdb_fasta$ID = str_remove(GtRNAdb_fasta$ID,'Homo_sapiens_')
GtRNAdb_fasta$ID = t(data.frame(strsplit(GtRNAdb_fasta$ID,' '),stringsAsFactors = F))[,1]
which(!(GtRNAdb_fasta$ID %in% databases$GtRNAdb$ID))
databases$GtRNAdb$Fasta_seq = NA
foreach(i = 1:nrow(databases$GtRNAdb),.combine = 'c') %do% {
  databases$GtRNAdb$Fasta_seq[i] = 
    GtRNAdb_fasta[which(GtRNAdb_fasta$ID == databases$GtRNAdb$ID[i])[1],2]
}

miRBase_fasta = read.table('miRBase.fasta.tsv',stringsAsFactors = F,sep = '\t')
colnames(miRBase_fasta) = c('ID','Fasta_seq')
miRBase_fasta$ID = foreach(i = 1:nrow(miRBase_fasta),.combine = 'c') %do% {
  strsplit(miRBase_fasta$ID[i],' ')[[1]][1]
}
which(!(miRBase_fasta$ID %in% databases$miRBase$ID))
databases$miRBase$Fasta_seq = NA
foreach(i = 1:nrow(databases$miRBase),.combine = 'c') %do% {
  databases$miRBase$Fasta_seq[i] = 
    miRBase_fasta[which(miRBase_fasta$ID == databases$miRBase$ID[i])[1],2]
}

pirnadb_fasta = read.table('piRNAdb.hsa.v1_7_6.fa.tsv',stringsAsFactors = F)
colnames(pirnadb_fasta) = c('ID','Fasta_seq')
which(!(pirnadb_fasta$ID %in% databases$piRNAdb$ID))
databases$piRNAdb$Fasta_seq = NA
foreach(i = 1:nrow(databases$piRNAdb),.combine = 'c') %do% {
  if (i %% 50000 == 0) {print(i)}
  databases$piRNAdb$Fasta_seq[i] = 
    pirnadb_fasta[which(pirnadb_fasta$ID == databases$piRNAdb$ID[i])[1],2]
}
databases$MINTbase$ID = MINTbase_ID_seq

for (f in 1:4) {
  databases[[f]]$equal_seq = databases[[f]]$Loci_seq == databases[[f]]$Fasta_seq
}


#Sam-file stats####
miRBase_sam = rbind.data.frame(read.table('miRBase_mature_sam_stats.txt',
                                          stringsAsFactors = F,sep = '\t'),
                               read.table('miRBase_hairpin_sam_stats.txt',
                                          stringsAsFactors = F,sep = '\t'))
colnames(miRBase_sam) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
miRBase_sam$Fasta_start = miRBase_sam$Fasta_start - 1
miRBase_sam$Fasta_mapped = foreach(i = miRBase_sam$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

GtRNAdb_sam = read.table('hg38-tRNAs_sam_stats.txt',stringsAsFactors = F,sep = '\t')
colnames(GtRNAdb_sam) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
GtRNAdb_sam$Fasta_start = GtRNAdb_sam$Fasta_start - 1
GtRNAdb_sam$ID = foreach(i = GtRNAdb_sam$ID,.combine = 'c') %do% {
  strsplit(i,'_')[[1]][3]
}
GtRNAdb_sam$Fasta_mapped = foreach(i = GtRNAdb_sam$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

MINTbase_sam = read.table('MINTbase_sam_stats.tsv',stringsAsFactors = F)
colnames(MINTbase_sam) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
MINTbase_sam$Fasta_start = MINTbase_sam$Fasta_start - 1
MINTbase_sam$ID = foreach(i = MINTbase_sam$ID,.combine = 'c') %do% {
  strsplit(i,'_')[[1]][2]
}
MINTbase_sam$Fasta_mapped = foreach(i = MINTbase_sam$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

pirnadb_sam = read.table('pirnadb_sam_stats.txt',stringsAsFactors = F,sep = '\t')
colnames(pirnadb_sam) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                           'Fasta_len','Fasta_mapped')
pirnadb_sam$Fasta_start = pirnadb_sam$Fasta_start - 1
pirnadb_sam$Fasta_mapped = foreach(i = pirnadb_sam$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

sam_stats = list(MINTbase_sam,miRBase_sam,pirnadb_sam,GtRNAdb_sam)
names(sam_stats) = names(databases)

for (x in 1:4) {
  print(x)
  d = databases[[x]]
  d$mapped2loci = vector(length = nrow(d))
  d$equal_IDs = vector(length = nrow(d))
  d$Fasta_mapped = vector(length = nrow(d))
  for (j in 1:nrow(d)) {
      pos = j
      if (pos > nrow(d)) {break}
      d$mapped2loci[pos] = 
        mean(str_c(d$Loci_chr[pos],d$Loci_start[pos]) == str_c(
          sam_stats[[x]]$Fasta_chr[which(sam_stats[[x]]$ID == d$ID[pos])],
          sam_stats[[x]]$Fasta_start[which(sam_stats[[x]]$ID == d$ID[pos])]
          )) > 0
      d$Fasta_mapped[pos] = sam_stats[[x]]$Fasta_mapped[which(sam_stats[[x]]$ID == d$ID[pos])][1]
      d$equal_IDs[pos] = length(which(sam_stats[[x]]$ID == d$ID[pos]))
      
    
  }
  databases[[x]] = d
}

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



#reverse compliment seq####
for (x in 1:4) {
  print(x)
  databases[[x]] = mutate(databases[[x]],Loci_rev_seq = rev_comp(Loci_seq))
}
for (x in 1:4) {
  databases[[x]] = mutate(databases[[x]],equal_rev_seq = Loci_rev_seq == Fasta_seq)
}
databases$MINTbase$equal_rev_seq = foreach(i = 1:nrow(MINTbase),.combine = 'c') %do% {
  str_c(strsplit(databases$MINTbase$Fasta_seq[i],databases$MINTbase$Loci_rev_seq[i])[[1]],
        collapse = '') != 
    databases$MINTbase$Fasta_seq[i]
}

#Fix MINTbase$equal_seq####
databases$MINTbase$equal_seq = foreach(i = 1:nrow(MINTbase),.combine = 'c') %do% {
  str_c(strsplit(databases$MINTbase$Fasta_seq[i],databases$MINTbase$Loci_seq[i])[[1]],
        collapse = '') != 
    databases$MINTbase$Fasta_seq[i]
}
databases$MINTbase$Fasta_len = foreach(i = 1:nrow(MINTbase),.combine = 'c') %do% {
  length(strsplit(databases$MINTbase$Fasta_seq[i],'')[[1]])
}
databases_cutted = foreach(i = 1:4) %do% {
  filter(databases[[i]], equal_seq == T)
}
names(databases_cutted) = names(databases)

#Value####
for (x in 1:4) {
  databases[[x]] = mutate(databases[[x]], value = 
           (equal_seq == T || equal_rev_seq == T) + (mapped2loci == T) + (Fasta_mapped == '1'))
}
write.xlsx(databases[1:4],'smRNA_databases_full.xlsx')

#rRNA add####

rRNA_ucsc_loci = read.table('rRNA_annot_ucsc.gtf',stringsAsFactors = F)
write.table(cbind(rRNA_ucsc_loci[,c(1,4,5)],row.names(rRNA_ucsc_loci)), 'rRNA_ucsc_loci.bed',sep = '\t',row.names = F)
rRNA_ucsc_getfasta = read.table('rRNA_ucsc_loci_getfasta.tsv',stringsAsFactors = F)
colnames(rRNA_ucsc_getfasta) = c('Loci','Loci_seq')
rRNA_ucsc_getfasta$Loci_rev_seq = rev_comp(rRNA_ucsc_getfasta$Loci_seq)

rRNA_Silva_fasta = read.table('rRNA_Silva_v138_1.tsv',stringsAsFactors = F)
rRNA_Silva_sam = read.table('rRNA_Silva_v138_1_sam_stats.tsv',stringsAsFactors = F)
colnames(rRNA_Silva_sam) = c('ID','Fasta_chr','Fasta_start','Fasta_seq',
                          'Fasta_len','Fasta_mapped')
rRNA_Silva_sam$Fasta_start = rRNA_Silva_sam$Fasta_start - 1
rRNA_Silva_sam$Fasta_mapped = foreach(i = rRNA_Silva_sam$Fasta_mapped,.combine = 'c') %do% {
  strsplit(i,'NH:i:')[[1]][2]
}

databases$rRNA_ucsc_Silva = data.frame(
  gene_ID = rRNA_ucsc_loci$V13,
  Loci = str_c(rRNA_ucsc_loci$V1,':',rRNA_ucsc_loci$V4,'-',rRNA_ucsc_loci$V5),
  Loci_chr = rRNA_ucsc_loci$V1,
  Loci_start = rRNA_ucsc_loci$V4,
  Loci_end = rRNA_ucsc_loci$V5,
  Loci_len = as.numeric(rRNA_ucsc_loci$V5) - as.numeric(rRNA_ucsc_loci$V4),
  stringsAsFactors = F
)
databases$rRNA_ucsc_Silva = cbind.data.frame(databases$rRNA_ucsc_Silva,
      foreach(i = 1:nrow(databases$rRNA_ucsc_Silva),.combine = 'rbind') %do% {
        filter(rRNA_ucsc_getfasta,Loci == databases$rRNA_ucsc_Silva$Loci[i])[1,2:3]
      })

databases$rRNA_ucsc_Silva = cbind.data.frame(databases$rRNA_ucsc_Silva,
      foreach(i = 1:nrow(databases$rRNA_ucsc_Silva),.combine = 'rbind') %do% {                             
        filter(rRNA_Silva_sam,
   strsplit(Fasta_seq,databases$rRNA_ucsc_Silva$Loci_seq[i])[[1]][1] != 
                 Fasta_seq[i])[1,]
      })       


rRNA_Silva_sam = mutate(rRNA_Silva_sam, Fasta_end = Fasta_start + Fasta_len)
write.table(filter(rRNA_Silva_sam,Fasta_chr != '*')[,c(2,3,7)], 'rRNA_Silva_sam.bed',sep = '\t',row.names = F)

rRNA_inter = read.table('rRNA_inter.bed',stringsAsFactors = F)
colnames(rRNA_inter) = c('Fasta_chr','Fasta_start','Fasta_end',
                         'Loci_chr','Loci_start','Loci_end','Loci_num')

rRNA_inter_ucsc_Silva = foreach(i = unique(rRNA_inter$Loci_num),.combine = 'rbind') %do% {
  Fasta_ID = filter(rRNA_Silva_sam, Fasta_chr %in% filter(rRNA_inter,Loci_num == i)$Fasta_chr,
         Fasta_start %in% filter(rRNA_inter,Loci_num == i)$Fasta_start)$ID
  foreach(f = Fasta_ID,.combine = 'rbind') %do% {
    data.frame(Fasta_ID = Fasta_ID,
                    databases$rRNA_ucsc_Silva[i,],
         stringsAsFactors = F)
  }
}
rRNA_inter_ucsc_Silva = unique(rRNA_inter_ucsc_Silva)
rRNA_inter_ucsc_Silva$Fasta_len = foreach(i = 1:nrow(rRNA_inter_ucsc_Silva),.combine = 'c') %do% {
  filter(rRNA_Silva_sam,ID == rRNA_inter_ucsc_Silva$Fasta_ID[i])$Fasta_len[1]
}
rRNA_inter_ucsc_Silva$Fasta_mapped = foreach(i = 1:nrow(rRNA_inter_ucsc_Silva),.combine = 'c') %do% {
  filter(rRNA_Silva_sam,ID == rRNA_inter_ucsc_Silva$Fasta_ID[i])$Fasta_mapped[1]
}

write.xlsx(databases$rRNA_ucsc_Silva, 'rRNA_ucsc_loci.xlsx')
write.xlsx(rRNA_inter_ucsc_Silva,'rRNA_inter_ucsc_Silva.xlsx')

#posthoc stats####
databases_stats = data.frame(
  Type = names(databases),
  num_loci = foreach(i = 1:4,.combine = 'c') %do% {nrow(databases[[i]])},
  uniqueIDs = foreach(i = 1:4,.combine = 'c') %do% {length(unique(databases[[i]]$ID))},
  Loci_Fasta = foreach(i = 1:4,.combine = 'c') %do% {sum(!is.na(databases[[i]]$equal_seq))},
  Loci_Fasta_equal = foreach(i = 1:4,.combine = 'c') %do% 
    {sum(c(databases[[i]]$equal_seq, databases[[i]]$equal_rev_seq),na.rm = T)},
  Fasta_mapped2Loci = foreach(i = 1:4,.combine = 'c') %do% 
    {sum(databases[[i]]$mapped2loci,na.rm = T)},
  Loci_Fasta_mapped_once = foreach(i = 1:4,.combine = 'c') %do% 
    {sum(databases[[i]]$Fasta_mapped == '1',na.rm = T)},
  value_1 = foreach(i = 1:4,.combine = 'c') %do% 
    {sum(databases[[i]]$value == 1,na.rm = T)},
  value_2 = foreach(i = 1:4,.combine = 'c') %do% 
    {sum(databases[[i]]$value == 2,na.rm = T)},
  value_3 = foreach(i = 1:4,.combine = 'c') %do% 
    {sum(databases[[i]]$value == 3,na.rm = T)},
  
  stringsAsFactors = F
)

write.xlsx(databases_stats,'smRNA_databases_stats.xlsx')

#Pictures####

venn.diagram(list(unique(databases$MINTbase$ID), unique(MINTbase_fasta$ID_seq)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'Venn_MINTbase_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('orange','blue'), cat.dist = c(-0.1,0.02),
             main = 'Unique ID of MINTbase', fill = c('red','blue'))

venn.diagram(list(unique(databases$miRBase$ID), unique(miRBase_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'Venn_miRBase_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red3','blue'), cat.dist = c(0.04,0.03),
             main = 'Unique ID of miRBase', fill = c('red','blue'))

venn.diagram(list(unique(databases$piRNAdb$ID), unique(pirnadb_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'Venn_piRNAdb_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red3','blue'), cat.dist = c(0.04,0.03),
             main = 'Unique ID of piRNAdb', fill = c('red','blue'))

venn.diagram(list(unique(databases$GtRNAdb$ID), unique(GtRNAdb_fasta$ID)), 
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'Venn_GtRNAdb_uniqueID.tiff'
             ,category.names = c('Loci_ID', 'Fasta_ID'),
             cat.col = c('red3','blue'), cat.dist = c(0.04,0.03),
             main = 'Unique ID of GtRNAdb', fill = c('red','blue'))

ggplot(rRNA_inter_ucsc_Silva) + aes(Loci_len,Fasta_len, color = Fasta_mapped) +
  theme_classic() + geom_point() + geom_abline() +
  ggtitle('Comparison of length of rRNA ucsc loci and Silva fasta seqs')
ggsave('Comparison of length of rRNA.tiff')

ggplot(databases$MINTbase) + aes(Loci_len,Fasta_len, color = as.factor(value)) +
  theme_classic() + geom_point() + geom_abline() +
  ggtitle('Comparison of length of tRNA-derived MINTbase loci and fasta')
ggsave('Comparison of length of MINTbase.tiff')

ggplot(databases$MINTbase) + aes(Loci_len - Fasta_len, fill = as.factor(value)) +
  theme_classic() + geom_histogram(binwidth = 1) + geom_abline() +
  ggtitle('Difference of length of tRNA-derived MINTbase loci and fasta')
ggsave('Difference of length of MINTbase.tiff',height = 10)

#####
databases$miRBase$strand = miRBase_gff3$V7
databases$piRNAdb$strand = piRNAdb_bedcov$V6
databases$GtRNAdb$strand = GtRNAdb_bedcov$V6

#Building united databases####

#Construct loci&fa data####
miRBase_precursors_loci_fa = filter(databases$miRBase,equal_seq+equal_rev_seq == 1,
                                    str_detect(ID,'hsa-mir'))
miRBase_mature_loci_fa = filter(databases$miRBase,equal_seq+equal_rev_seq == 1,
                                    str_detect(ID,'hsa-miR'))
piRNAdb_loci_fa = filter(databases$piRNAdb,equal_seq+equal_rev_seq == 1)
GtRNAdb_loci_fa = filter(databases$GtRNAdb,equal_seq+equal_rev_seq == 1)

smRNA_DATABASE = rbind.data.frame(
  data.frame(type = 'precursor_miRNA',
             ID = miRBase_precursors_loci_fa$ID,
             Annotation_consensus = miRBase_precursors_loci_fa$Loci,
             Annotation_chr = miRBase_precursors_loci_fa$Loci_chr,
             Annotation_start = miRBase_precursors_loci_fa$Loci_start,
             Annotation_end = miRBase_precursors_loci_fa$Loci_end,
             Sequence_consensus = miRBase_precursors_loci_fa$Fasta_seq,
             Length_consensus = miRBase_precursors_loci_fa$Loci_len,
             mapped2loci = miRBase_precursors_loci_fa$mapped2loci,
             strand = miRBase_precursors_loci_fa$strand,
             strand_from_annotation = T,
             Source = 'annotation&fasta',
             stringsAsFactors = F
             ),
  data.frame(type = 'mature_miRNA',
             ID = miRBase_mature_loci_fa$ID,
             Annotation_consensus = miRBase_mature_loci_fa$Loci,
             Annotation_chr = miRBase_mature_loci_fa$Loci_chr,
             Annotation_start = miRBase_mature_loci_fa$Loci_start,
             Annotation_end = miRBase_mature_loci_fa$Loci_end,
             Sequence_consensus = miRBase_mature_loci_fa$Fasta_seq,
             Length_consensus = miRBase_mature_loci_fa$Loci_len,
             mapped2loci = miRBase_mature_loci_fa$mapped2loci,
             strand = miRBase_mature_loci_fa$strand,
             strand_from_annotation = T,
             Source = 'annotation&fasta',
             stringsAsFactors = F
  ),
data.frame(type = 'piRNA',
           ID = piRNAdb_loci_fa$ID,
           Annotation_consensus = piRNAdb_loci_fa$Loci,
           Annotation_chr = piRNAdb_loci_fa$Loci_chr,
           Annotation_start = piRNAdb_loci_fa$Loci_start,
           Annotation_end = piRNAdb_loci_fa$Loci_end,
           Sequence_consensus = piRNAdb_loci_fa$Fasta_seq,
           Length_consensus = piRNAdb_loci_fa$Loci_len,
           mapped2loci = piRNAdb_loci_fa$mapped2loci,
           strand = piRNAdb_loci_fa$strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
),
data.frame(type = 'tRNA',
           ID = GtRNAdb_loci_fa$ID,
           Annotation_consensus = GtRNAdb_loci_fa$Loci,
           Annotation_chr = GtRNAdb_loci_fa$Loci_chr,
           Annotation_start = GtRNAdb_loci_fa$Loci_start,
           Annotation_end = GtRNAdb_loci_fa$Loci_end,
           Sequence_consensus = GtRNAdb_loci_fa$Fasta_seq,
           Length_consensus = GtRNAdb_loci_fa$Loci_len,
           mapped2loci = GtRNAdb_loci_fa$mapped2loci,
           strand = GtRNAdb_loci_fa$strand,
           strand_from_annotation = T,
           Source = 'annotation&fasta',
           stringsAsFactors = F
))

#Construct loci only data####
miRBase_gff3 = read.table('E:/Cygwin/home/Виталий/RSF_Project/databases/check_databases/miRBase.gff3',
                          stringsAsFactors = F)
miRBase_gff3$strand = miRBase_gff3$V7 == '-'
miRBase_gff3$ID = foreach(i = 1:nrow(miRBase_gff3),.combine = 'c') %do% {
  str_remove_all(strsplit(miRBase_gff3$V9[i],';')[[1]][3],'Name=')
}

miRBase_precursors_loci_only = filter(databases$miRBase, is.na(Fasta_seq),
                                    str_detect(ID,'hsa-mir'))
miRBase_precursors_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(miRBase_precursors_loci_only),.combine = 'c') %do% {
    seq = mean(filter(miRBase_gff3,ID == miRBase_precursors_loci_only$ID[i],
                V1 == miRBase_precursors_loci_only$Loci_chr[i],
                V4-1 == miRBase_precursors_loci_only$Loci_start[i],
                V5 == miRBase_precursors_loci_only$Loci_end[i],
           )$strand)
    if (seq == 1) {out = miRBase_precursors_loci_only$Loci_seq[i]}
    if (seq == 0) {out = miRBase_precursors_loci_only$Loci_rev_seq[i]}
    out
  }
miRBase_mature_loci_only = filter(databases$miRBase,is.na(Fasta_seq),
                                str_detect(ID,'hsa-miR'))
miRBase_mature_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(miRBase_mature_loci_only),.combine = 'c') %do% {
    seq = mean(filter(miRBase_gff3,ID == miRBase_mature_loci_only$ID[i],
                      V1 == miRBase_mature_loci_only$Loci_chr[i],
                      V4-1 == miRBase_mature_loci_only$Loci_start[i],
                      V5 == miRBase_mature_loci_only$Loci_end[i],
    )$strand)
    if (seq == 1) {out = miRBase_mature_loci_only$Loci_seq[i]}
    if (seq == 0) {out = miRBase_mature_loci_only$Loci_rev_seq[i]}
    out
  }
piRNAdb_loci_only = filter(databases$piRNAdb,is.na(Fasta_seq)) # is null

GtRNAdb_loci_only = filter(databases$GtRNAdb,is.na(Fasta_seq))
GtRNAdb_loci_only$Fasta_seq = 
  foreach(i = 1:nrow(GtRNAdb_loci_only),.combine = 'c') %do% {
    seq = mean(filter(GtRNAdb_bedcov,V4 == GtRNAdb_loci_only$ID[i],
                      V1 == GtRNAdb_loci_only$Loci_chr[i],
                      V2 == GtRNAdb_loci_only$Loci_start[i],
                      V3 == GtRNAdb_loci_only$Loci_end[i],
    )$V6 == '-')
    if (seq == 0) {out = GtRNAdb_loci_only$Loci_seq[i]}
    if (seq == 1) {out = GtRNAdb_loci_only$Loci_rev_seq[i]}
    out
  }

smRNA_DATABASE = rbind.data.frame(smRNA_DATABASE,
  data.frame(type = 'precursor_miRNA',
             ID = miRBase_precursors_loci_only$ID,
             Annotation_consensus = miRBase_precursors_loci_only$Loci,
             Annotation_chr = miRBase_precursors_loci_only$Loci_chr,
             Annotation_start = miRBase_precursors_loci_only$Loci_start,
             Annotation_end = miRBase_precursors_loci_only$Loci_end,
             Sequence_consensus = miRBase_precursors_loci_only$Fasta_seq,
             Length_consensus = miRBase_precursors_loci_only$Loci_len,
             mapped2loci = miRBase_precursors_loci_only$mapped2loci,
             strand = miRBase_precursors_loci_only$strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ),
  data.frame(type = 'mature_miRNA',
             ID = miRBase_mature_loci_only$ID,
             Annotation_consensus = miRBase_mature_loci_only$Loci,
             Annotation_chr = miRBase_mature_loci_only$Loci_chr,
             Annotation_start = miRBase_mature_loci_only$Loci_start,
             Annotation_end = miRBase_mature_loci_only$Loci_end,
             Sequence_consensus = miRBase_mature_loci_only$Fasta_seq,
             Length_consensus = miRBase_mature_loci_only$Loci_len,
             mapped2loci = miRBase_mature_loci_only$mapped2loci,
             strand = miRBase_mature_loci_only$strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ),
  data.frame(type = 'tRNA',
             ID = GtRNAdb_loci_only$ID,
             Annotation_consensus = GtRNAdb_loci_only$Loci,
             Annotation_chr = GtRNAdb_loci_only$Loci_chr,
             Annotation_start = GtRNAdb_loci_only$Loci_start,
             Annotation_end = GtRNAdb_loci_only$Loci_end,
             Sequence_consensus = GtRNAdb_loci_only$Fasta_seq,
             Length_consensus = GtRNAdb_loci_only$Loci_len,
             mapped2loci = GtRNAdb_loci_only$mapped2loci,
             strand = GtRNAdb_loci_only$strand,
             strand_from_annotation = T,
             Source = 'annotation',
             stringsAsFactors = F
  ))




#Construct fa only data####

miRBase_precursors_fa_only = filter(miRBase_sam, !(ID %in% databases$miRBase$ID),
                                      str_detect(ID,'hsa-mir'), Fasta_chr != '*')
miRBase_mature_fa_only = filter(miRBase_sam, !(ID %in% databases$miRBase$ID),
                                Fasta_chr != '*', str_detect(ID,'hsa-miR'))
piRNAdb_fa_only = filter(pirnadb_sam,!(ID %in% databases$piRNAdb$ID),
                         Fasta_chr != '*')
GtRNAdb_fa_only = filter(GtRNAdb_sam,!(ID %in% databases$GtRNAdb$ID),
                         Fasta_chr != '*')

smRNA_DATABASE = rbind.data.frame(smRNA_DATABASE,
                                  data.frame(type = 'precursor_miRNA',
                                             ID = miRBase_precursors_fa_only$ID,
        Annotation_consensus = str_c(miRBase_precursors_fa_only$Fasta_chr,':',
                                     miRBase_precursors_fa_only$Fasta_start,'-',
                                     as.character(miRBase_precursors_fa_only$Fasta_start+
                                                    miRBase_precursors_fa_only$Fasta_len)),
                                             Annotation_chr = miRBase_precursors_fa_only$Fasta_chr,
                                             Annotation_start = miRBase_precursors_fa_only$Fasta_start,
                                             Annotation_end = miRBase_precursors_fa_only$Fasta_start+
          miRBase_precursors_fa_only$Fasta_len,
                                             Sequence_consensus = miRBase_precursors_fa_only$Fasta_seq,
                                             Length_consensus = miRBase_precursors_fa_only$Fasta_len,
                                             mapped2loci = NA,
        strand = '+',
        strand_from_annotation = T,
                                             Source = 'fasta',
                                             stringsAsFactors = F
                                  ),
                                  data.frame(type = 'mature_miRNA',
                                             ID = miRBase_mature_fa_only$ID,
                                             Annotation_consensus = str_c(miRBase_mature_fa_only$Fasta_chr,':',
                                                   miRBase_mature_fa_only$Fasta_start,'-',
                                                   as.character(miRBase_mature_fa_only$Fasta_start+
                                                                  miRBase_mature_fa_only$Fasta_len)),
                                             Annotation_chr = miRBase_mature_fa_only$Fasta_chr,
                                             Annotation_start = miRBase_mature_fa_only$Fasta_start,
                                             Annotation_end = miRBase_mature_fa_only$Fasta_start+
                                               miRBase_mature_fa_only$Fasta_len,
                                             Sequence_consensus = miRBase_mature_fa_only$Fasta_seq,
                                             Length_consensus = miRBase_mature_fa_only$Fasta_len,
                                             mapped2loci = NA,
                                             strand = '+',
                                             strand_from_annotation = T,
                                             Source = 'fasta',
                                             stringsAsFactors = F
                                  ),
                                  data.frame(type = 'piRNA',
                                             ID = piRNAdb_fa_only$ID,
                                             Annotation_consensus = str_c(piRNAdb_fa_only$Fasta_chr,':',
                                                                          piRNAdb_fa_only$Fasta_start,'-',
                                                                          as.character(piRNAdb_fa_only$Fasta_start+
                                                                                         piRNAdb_fa_only$Fasta_len)),
                                             Annotation_chr = piRNAdb_fa_only$Fasta_chr,
                                             Annotation_start = piRNAdb_fa_only$Fasta_start,
                                             Annotation_end = piRNAdb_fa_only$Fasta_start+
                                               piRNAdb_fa_only$Fasta_len,
                                             Sequence_consensus = piRNAdb_fa_only$Fasta_seq,
                                             Length_consensus = piRNAdb_fa_only$Fasta_len,
                                             mapped2loci = NA,
                                             strand = '+',
                                             strand_from_annotation = T,
                                             Source = 'fasta',
                                             stringsAsFactors = F
                                  ),
                                  data.frame(type = 'tRNA',
                                             ID = GtRNAdb_fa_only$ID,
                                             Annotation_consensus = str_c(GtRNAdb_fa_only$Fasta_chr,':',
                                                                          GtRNAdb_fa_only$Fasta_start,'-',
                                                                          as.character(GtRNAdb_fa_only$Fasta_start+
                                                                                         GtRNAdb_fa_only$Fasta_len)),
                                             Annotation_chr = GtRNAdb_fa_only$Fasta_chr,
                                             Annotation_start = GtRNAdb_fa_only$Fasta_start,
                                             Annotation_end = GtRNAdb_fa_only$Fasta_start+
                                               GtRNAdb_fa_only$Fasta_len,
                                             Sequence_consensus = GtRNAdb_fa_only$Fasta_seq,
                                             Length_consensus = GtRNAdb_fa_only$Fasta_len,
                                             mapped2loci = NA,
                                             strand = '+',
                                             strand_from_annotation = T,
                                             Source = 'fasta',
                                             stringsAsFactors = F
                                  ))







#Construct MINTbase data#####

MINTbase = read.table('MINTbase_strand.bed',stringsAsFactors = F)
databases$MINTbase$strand = MINTbase$V5

MINTbase_loci_fa = filter(databases$MINTbase,equal_seq+equal_rev_seq == 1)
MINTbase_loci_fa$diff_len = MINTbase_loci_fa$Fasta_len - MINTbase_loci_fa$Loci_len

MINTbase_loci_only = filter(databases$MINTbase,is.na(Fasta_seq)) #null
MINTbase_fa_only = filter(MINTbase_sam,!(ID %in% databases$MINTbase$ID),
                         Fasta_chr != '*')

MINTbase_loci_fa = cbind.data.frame(MINTbase_loci_fa,
       foreach(i = 1:nrow(MINTbase_loci_fa),.combine = 'rbind') %dopar% {
         if (MINTbase_loci_fa$equal_seq[i] == T) {
           widestr = strsplit(MINTbase_loci_fa$Fasta_seq[i],MINTbase_loci_fa$Loci_seq[i])[[1]]
         }
         if (MINTbase_loci_fa$equal_rev_seq[i] == T) {
           widestr = strsplit(MINTbase_loci_fa$Fasta_seq[i],MINTbase_loci_fa$Loci_rev_seq[i])[[1]]
         }
         if (is.na(widestr[2])) {widestr[2] = ''}
         data.frame(left_wide = length(strsplit(widestr[1],'')[[1]]),
                    right_wide = length(strsplit(widestr[2],'')[[1]]))
       }) 

MINTbase_loci_fa$start_wide = MINTbase_loci_fa$Loci_start - MINTbase_loci_fa$left_wide
MINTbase_loci_fa$end_wide = MINTbase_loci_fa$Loci_end + MINTbase_loci_fa$right_wide


smRNA_DATABASE = rbind.data.frame(smRNA_DATABASE,
                                  data.frame(type = 'tRNA-derived',
                                             ID = MINTbase_fa_only$ID,
                                             Annotation_consensus = str_c(MINTbase_fa_only$Fasta_chr,':',
                                                                          MINTbase_fa_only$Fasta_start,'-',
                                                                          as.character(MINTbase_fa_only$Fasta_start+
                                                                                         MINTbase_fa_only$Fasta_len)),
                                             Annotation_chr = MINTbase_fa_only$Fasta_chr,
                                             Annotation_start = MINTbase_fa_only$Fasta_start,
                                             Annotation_end = MINTbase_fa_only$Fasta_start+
                                               MINTbase_fa_only$Fasta_len,
                                             Sequence_consensus = MINTbase_fa_only$Fasta_seq,
                                             Length_consensus = MINTbase_fa_only$Fasta_len,
                                             mapped2loci = NA,
                                             strand = '+',
                                             strand_from_annotation = F,
                                             Source = 'fasta',
                                             stringsAsFactors = F
                                  ),
                                  data.frame(type = 'tRNA-derived',
                                             ID = MINTbase_loci_fa$ID,
                                             Annotation_consensus = str_c(MINTbase_loci_fa$Loci_chr,':',
                                                                          MINTbase_loci_fa$start_wide,'-',
                                                                          MINTbase_loci_fa$end_wide),
                                             Annotation_chr = MINTbase_loci_fa$Loci_chr,
                                             Annotation_start = MINTbase_loci_fa$start_wide,
                                             Annotation_end = MINTbase_loci_fa$end_wide,
                                             Sequence_consensus = MINTbase_loci_fa$Fasta_seq,
                                             Length_consensus = MINTbase_loci_fa$Fasta_len,
                                             mapped2loci = MINTbase_loci_fa$mapped2loci,
                                             strand = MINTbase_loci_fa$strand,
                                             strand_from_annotation = T,
                                             Source = 'wide_annotation&fasta',
                                             stringsAsFactors = F
                                  ))


#Construct loci fa conflict data (piRNAdb & MINTbase)#####
piRNAdb_loci_not_fa = filter(databases$piRNAdb,equal_seq+equal_rev_seq == 0)

smRNA_DATABASE = rbind.data.frame(smRNA_DATABASE,
                                  data.frame(type = 'piRNA',
                                             ID = piRNAdb_loci_not_fa$ID,
                                             Annotation_consensus = piRNAdb_loci_not_fa$Loci,
                                             Annotation_chr = piRNAdb_loci_not_fa$Loci_chr,
                                             Annotation_start = piRNAdb_loci_not_fa$Loci_start,
                                             Annotation_end = piRNAdb_loci_not_fa$Loci_end,
                                             Sequence_consensus = piRNAdb_loci_not_fa$Fasta_seq,
                                             Length_consensus = piRNAdb_loci_not_fa$Loci_len,
                                             mapped2loci = piRNAdb_loci_not_fa$mapped2loci,
                                             strand = piRNAdb_loci_not_fa$strand,
                                             strand_from_annotation = T,
                                             Source = 'annotation&fasta',
                                             stringsAsFactors = F
                                  ))


MINTbase_loci_not_fa = filter(databases$MINTbase,equal_seq+equal_rev_seq == 0)
MINTbase_loci_not_fa$diff_len = MINTbase_loci_not_fa$Fasta_len - MINTbase_loci_not_fa$Loci_len
intersequences = function(x,y,mismatch = 0) {
  x = strsplit(x,'')[[1]]
  y = strsplit(y,'')[[1]]
  baseseq = list(x,y)[[which.min(c(length(x),length(y)))]]
  moveseq = list(x,y)[[which.max(c(length(x),length(y)))]]
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
  if (is.na(out[,2])) {return(out)}
  else if (out[,2] == '1') {out[,2] = str_c(moveseq[1:out[1,1]], collapse = '')
  return(out)}
  else if (out[,2] == '0') {out[,2] = str_c(baseseq[1:out[1,1]], collapse = '')
  return(out)}
  
}

MINTbase_loci_not_fa = cbind.data.frame(MINTbase_loci_not_fa,
      foreach(i = 1:nrow(MINTbase_loci_not_fa),.combine = 'rbind',
              .packages = c('foreach','stringr')) %dopar% {
        inter = cbind(
          intersequences(MINTbase_loci_not_fa$Fasta_seq[i],MINTbase_loci_not_fa$Loci_seq[i]),
          intersequences(MINTbase_loci_not_fa$Fasta_seq[i],MINTbase_loci_not_fa$Loci_rev_seq[i]))
        data.frame(inter_length = max(inter[1,1],inter[1,3]),
                   inter_seq = inter[1,which.max(c(inter[1,1],inter[1,3]))*2],
                   stringsAsFactors = F)
      })

MINTbase_loci_not_fa = cbind.data.frame(MINTbase_loci_not_fa,
                                        foreach(i = 1:nrow(MINTbase_loci_not_fa),.combine = 'rbind',
                                                .packages = c('foreach','stringr')) %dopar% {
                                                  inter = cbind(
                                                    intersequences(MINTbase_loci_not_fa$Fasta_seq[i],
                                                                   MINTbase_loci_not_fa$Loci_seq[i], mismatch = 2),
                                                    intersequences(MINTbase_loci_not_fa$Fasta_seq[i],
                                                                   MINTbase_loci_not_fa$Loci_rev_seq[i], mismatch = 2))
                                                  data.frame(inter_mis2_length = max(inter[1,1],inter[1,3]),
                                                             inter_mis2_seq = inter[1,which.max(c(inter[1,1],inter[1,3]))*2],
                                                             stringsAsFactors = F)
                                                })

MINTbase_loci_not_fa = cbind.data.frame(MINTbase_loci_not_fa,
                                        foreach(i = 1:nrow(MINTbase_loci_not_fa),.combine = 'rbind',
                                                .packages = c('foreach','stringr')) %dopar% {
                                                  inter = cbind(
                                                    intersequences(MINTbase_loci_not_fa$Fasta_seq[i],
                                                                   MINTbase_loci_not_fa$Loci_seq[i], mismatch = 2999),
                                                    intersequences(MINTbase_loci_not_fa$Fasta_seq[i],
                                                                   MINTbase_loci_not_fa$Loci_rev_seq[i], mismatch = 2999))
                                                  data.frame(inter_mis_length = max(inter[1,1],inter[1,3]),
                                                             stringsAsFactors = F)
                                                })

ggplot(MINTbase_loci_not_fa) + aes(inter_length) + theme_classic() +
  geom_histogram(binwidth = 0.5) + scale_x_continuous(breaks = 0:30) +
  scale_y_continuous(breaks = (0:5)*100) +
  ggtitle('Intersection of Fasta and Annotation sequences for MINTbase')
ggsave('Intersection of Fasta and Annotation sequences for MINTbase.tiff')
ggplot(MINTbase_loci_not_fa) + aes(inter_mis2_length) + theme_classic() +
  geom_histogram(binwidth = 0.5) + scale_x_continuous(breaks = 0:50) +
  scale_y_continuous(breaks = (0:5)*100) +
  ggtitle('Intersection of Fasta and Annotation sequences for MINTbase',
          subtitle = '2 mismatches allowed')
ggsave('Intersection of Fasta and Annotation sequences for MINTbase 2 mismatches.tiff')
ggplot(MINTbase_loci_not_fa) + aes(inter_mis2_length, inter_length) + theme_classic() +
  geom_point() + scale_x_continuous(breaks = 0:50) +
  scale_y_continuous(breaks = (0:50)) + geom_abline() +
  ggtitle('Intersection of Fasta and Annotation sequences for MINTbase',
          subtitle = '0 vs 2 mismatches allowed')
ggsave('Intersection of Fasta and Annotation sequences for MINTbase 0vs2 mismatches.tiff')
MINTbase_loci_not_fa$min_len = foreach(i = 1:nrow(MINTbase_loci_not_fa),.combine = 'c') %do% {
  min(MINTbase_loci_not_fa$Fasta_len[i],MINTbase_loci_not_fa$Loci_len[i])
}

MINTbase_loci_more_fa = filter(MINTbase_loci_not_fa, diff_len < 0)

#fa in loci####
MINTbase_loci_not_fa$fa_in_loci = foreach(i = 1:nrow(MINTbase_loci_not_fa),.combine = 'c') %dopar% {
  strsplit(MINTbase_loci_not_fa$Loci_seq[i],MINTbase_loci_not_fa$Fasta_seq[i])[[1]] != 
    MINTbase_loci_not_fa$Loci_seq[i] ||
  strsplit(MINTbase_loci_not_fa$Loci_rev_seq[i],MINTbase_loci_not_fa$Fasta_seq[i])[[1]] != 
    MINTbase_loci_not_fa$Loci_rev_seq[i]
}
MINTbase_fa_in_loci = filter(MINTbase_loci_not_fa, inter_length == min_len)

MINTbase_fa_in_loci$wide_fasta_seq = foreach(i = 1:nrow(MINTbase_fa_in_loci),.combine = 'rbind') %dopar% {
        widestr = strsplit(MINTbase_fa_in_loci$Loci_seq[i],MINTbase_fa_in_loci$Fasta_seq[i])[[1]]
        if (widestr == MINTbase_fa_in_loci$Loci_seq[i]) {
           o = MINTbase_fa_in_loci$Loci_rev_seq[i]
        }
        else{o = MINTbase_fa_in_loci$Loci_seq[i]}
        o
} 
#suff-pref data####
ggarrange(
ggplot(MINTbase_loci_not_fa) + aes(Mod(diff_len), inter_length, 
                                   color = min_len - inter_length <= 1) +
  theme_classic() + geom_point() + 
  scale_y_continuous(breaks = 0:30, name = 'Intersection length') +
  scale_x_continuous(breaks = 0:30, name = 'Difference of length') + 
  ggtitle('Comparison of fasta-getfasta intersection length with sequenses length'),
ggplot(MINTbase_loci_not_fa) + aes(min_len, inter_length,
                                   color = min_len - inter_length <= 1) +
  theme_classic() + geom_point() + geom_abline() + 
  scale_y_continuous(breaks = 0:30, name = 'Intersection length') +
  scale_x_continuous(breaks = 0:100, name = 'Min length'),
ncol = 1)
ggsave('Comparison_intersection_length_MINTbase.tiff',height = 8, width = 10)

MINTbase_loci_fa_suffpref = filter(MINTbase_loci_not_fa, min_len - inter_length == 1)
MINTbase_loci_fa_suffpref$left_wide_fasta = 
  unlist(strsplit(MINTbase_loci_fa_suffpref$Loci_rev_seq, 
                  MINTbase_loci_fa_suffpref$inter_seq))
MINTbase_loci_fa_suffpref$right_wide_fasta = 
  t(as.data.frame(strsplit(MINTbase_loci_fa_suffpref$Fasta_seq, 
                  MINTbase_loci_fa_suffpref$inter_seq))[2,])
  
MINTbase_loci_fa_suffpref$left_wide = MINTbase_loci_fa_suffpref$Loci_start -
  unlist(rapply(strsplit(MINTbase_loci_fa_suffpref$left_wide_fasta,''), length, how="list"))
MINTbase_loci_fa_suffpref$right_wide = MINTbase_loci_fa_suffpref$Loci_end +
  unlist(rapply(strsplit(MINTbase_loci_fa_suffpref$right_wide_fasta,''), length, how="list"))

MINTbase_loci_fa_suffpref$wide_fasta = str_c(MINTbase_loci_fa_suffpref$left_wide_fasta,
                                             MINTbase_loci_fa_suffpref$Fasta_seq,
                                as.character(MINTbase_loci_fa_suffpref$right_wide_fasta))


smRNA_DATABASE = rbind.data.frame(smRNA_DATABASE,
                                  data.frame(type = 'tRNA-derived',
                                             ID = MINTbase_fa_in_loci$ID,
                                             Annotation_consensus = MINTbase_fa_in_loci$Loci,
                                             Annotation_chr = MINTbase_fa_in_loci$Loci_chr,
                                             Annotation_start = MINTbase_fa_in_loci$Loci_start,
                                             Annotation_end = MINTbase_fa_in_loci$Loci_end,
                                             Sequence_consensus = MINTbase_fa_in_loci$wide_fasta_seq,
                                             Length_consensus = MINTbase_fa_in_loci$Loci_len,
                                             mapped2loci = MINTbase_fa_in_loci$mapped2loci,
                                             strand = MINTbase_fa_in_loci$strand,
                                             strand_from_annotation = T,
                                             Source = 'annotation&wide_fasta',
                                             stringsAsFactors = F
                                  ),
                                  data.frame(type = 'tRNA-derived',
                                             ID = MINTbase_loci_fa_suffpref$ID,
                                             Annotation_consensus = str_c(MINTbase_loci_fa_suffpref$Loci_chr,':',
                                                                          MINTbase_loci_fa_suffpref$left_wide,'-',
                                                                          MINTbase_loci_fa_suffpref$right_wide),
                                             Annotation_chr = MINTbase_loci_fa_suffpref$Loci_chr,
                                             Annotation_start = MINTbase_loci_fa_suffpref$left_wide,
                                             Annotation_end = MINTbase_loci_fa_suffpref$right_wide,
                                             Sequence_consensus = MINTbase_loci_fa_suffpref$wide_fasta,
                                             Length_consensus = MINTbase_loci_fa_suffpref$right_wide - 
                                               MINTbase_loci_fa_suffpref$left_wide,
                                             mapped2loci = MINTbase_loci_fa_suffpref$mapped2loci,
                                             strand = MINTbase_loci_fa_suffpref$strand,
                                             strand_from_annotation = T,
                                             Source = 'wide_annotation&wide_fasta',
                                             stringsAsFactors = F
                                  ))

#Count unions#####

Counts_smRNA_databases = foreach(i = names(table(smRNA_DATABASE$type))) %do% {
  table(filter(smRNA_DATABASE,type == i)$Source)
}
names(Counts_smRNA_databases) = names(table(smRNA_DATABASE$type))

#rRNA add####

databases$rRNA_ucsc_Silva$strand = rRNA_ucsc_loci$V7
databases$rRNA_ucsc_Silva = na.omit(databases$rRNA_ucsc_Silva)

smRNA_DATABASE = rbind.data.frame(smRNA_DATABASE,
                                  data.frame(type = 'rRNA',
                                             ID = databases$rRNA_ucsc_Silva$gene_ID,
                                             Annotation_consensus = databases$rRNA_ucsc_Silva$Loci,
                                             Annotation_chr = databases$rRNA_ucsc_Silva$Loci_chr,
                                             Annotation_start = databases$rRNA_ucsc_Silva$Loci_start,
                                             Annotation_end = databases$rRNA_ucsc_Silva$Loci_end,
                                             Sequence_consensus = 
      foreach(i = 1:nrow(databases$rRNA_ucsc_Silva),.combine = 'c') %do% {
        if (databases$rRNA_ucsc_Silva$strand[i] == '+'){
          o = databases$rRNA_ucsc_Silva$Loci_seq[i]
        }
        else {o = databases$rRNA_ucsc_Silva$Loci_rev_seq[i]}
        o
      },
                                             Length_consensus = databases$rRNA_ucsc_Silva$Loci_len,
                                             mapped2loci = NA,
      strand = databases$rRNA_ucsc_Silva$strand,
      strand_from_annotation = T,
                                             Source = 'annotation',
                                             stringsAsFactors = F
                                  ))

smRNA_DATABASE = filter(smRNA_DATABASE,!(is.na(Sequence_consensus)))
smRNA_DATABASE$repeats = repeats
write.xlsx(smRNA_DATABASE,'smRNA_DATABASE.xlsx')

for (i in 1:length(unique(filter(mipitrrib_data_no_inters,type == 'rRNA')[,7]))){
  mipitrrib_data_no_inters$ID[which(mipitrrib_data_no_inters$Sequence_consensus == unique(filter(mipitrrib_data_no_inters,type == 'rRNA')[,7])[i] )] = str_c(mipitrrib_data_no_inters$ID[which(mipitrrib_data_no_inters$Sequence_consensus == unique(filter(mipitrrib_data_no_inters,type == 'rRNA')[,7])[i] )],'_seq_',i)
}

#1718 MINTbase loci#####

MINTbase_loci_not_fa_unsorted = 
  filter(MINTbase_loci_not_fa, min_len - inter_length > 1)

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

MINTbase_loci_not_fa_unsorted$similarity = 
  foreach(i = 1:nrow(MINTbase_loci_not_fa_unsorted),.combine = 'c',
          .packages = c('foreach','stringr')) %dopar% {
    min(hamming_dist_unequal(MINTbase_loci_not_fa_unsorted$Loci_seq[i],
                         MINTbase_loci_not_fa_unsorted$Fasta_seq[i]),
    hamming_dist_unequal(MINTbase_loci_not_fa_unsorted$Loci_rev_seq[i],
                         MINTbase_loci_not_fa_unsorted$Fasta_seq[i]))
  }
ggplot(MINTbase_loci_not_fa_unsorted) + aes(similarity) + theme_classic() +
  geom_histogram(binwidth = 0.5) + scale_x_continuous(breaks = 0:30) +
  scale_y_continuous(breaks = (0:5)*100) +
  ggtitle('Similarity of sequences for MINTbase',subtitle = 'n = 1718')
ggsave('Similarity of sequences for MINTbase.tiff')

#Count repeats in smRNA_DATABASE#####

table_smRNA_DATABASE_ID = table(smRNA_DATABASE$ID)
smRNA_DATABASE$repeats = foreach(i = smRNA_DATABASE$ID,.combine = 'c') %dopar% {
  table_smRNA_DATABASE_ID[i]
}
for (i in names(table(smRNA_DATABASE$type))) {
  ggplot(filter(smRNA_DATABASE, type == i, repeats > 1)) + aes(y = repeats) + theme_classic() +
    geom_boxplot() +
    ggtitle('Number of loci for one ID', subtitle = str_c('type = ',i))
  ggsave(str_c('Number of loci for one ID for ',i,'.tiff'))
}

ggplot(data.frame(num = as.vector(table(smRNA_DATABASE$repeats)),
                  similarity = as.numeric(names(table(smRNA_DATABASE$repeats))),
                  n = as.vector(table(smRNA_DATABASE$repeats))/
                    as.numeric(names(table(smRNA_DATABASE$repeats))))[2:50,]) +
  theme_classic() + aes(similarity, n) + geom_point() +
  ggtitle('Distribution of number of repets')
ggsave('Distribution of number of repets.tiff')  

smRNA_DATABASE = smRNA_DATABASE[-(which(is.na(smRNA_DATABASE$Sequence_consensus))),]

#inter in mat mirna,pirna,trna, rRNA####

names(positions_smRNA_DATABASE) = foreach(i = 1:length(positions_smRNA_DATABASE),
                                          .combine = 'c')%do%{
  if (i %% 100000 == 0) {print(i)}
  positions_smRNA_DATABASE[[i]][1]
}

mipitrrib_data = filter(smRNA_DATABASE,type %in% c('tRNA','piRNA','mature_miRNA','rRNA'))
positions_mipitrrib_DATABASE = positions_smRNA_DATABASE[
  str_c(mipitrrib_data$Annotation_chr,'_',mipitrrib_data$Annotation_start)]
d_chr = foreach(chr = unique(mipitrrib_data$Annotation_chr),.combine = 'rbind') %do% {
  print(chr)
  d = filter(mipitrrib_data,Annotation_chr == chr)
  p = positions_mipitrrib_DATABASE[
    str_detect(names(positions_mipitrrib_DATABASE),str_c(chr,'_'))]
  positions = table(unlist(p))
  positions = positions[which(positions >1)]
  n = names(positions)
  
  d$interseptions = foreach(i = p,.combine = 'c') %dopar% {
    sum(i %in% n) 
  }/(d$Length_consensus+1)
  d
}

#interseptions ggplot#####
for (i in unique(d_chr$type)) {
  ggplot(filter(d_chr,type == i)) + aes(interseptions, fill = type) + theme_classic() +
    geom_histogram(binwidth = 0.05) + 
    scale_x_continuous(breaks = (0:10)/10, name = 'Part of sequence in interseption',
                       limits = c(-0.05,1.05)) +
    ggtitle('Number of records with interseptions in mature microRNA,piRNA,tRNA,rRNA database for', 
            subtitle = i)
  ggsave(str_c('interseptions_mi_pi_rib_tr/Number of records with interseptions in mi_pi_rib_tRNA database for ',i,'.tiff'))
}
ggplot(d_chr) + aes(interseptions, fill = type) + theme_classic() +
  geom_histogram(binwidth = 0.05) + 
  scale_x_continuous(breaks = (0:10)/10, name = 'Part of sequence in interseption',
                     limits = c(-0.05,1.05)) +
  ggtitle('Number of records with interseptions in mature microRNA,piRNA,tRNA,rRNA database for', 
          subtitle = 'mature miRNA,piRNA,tRNA')
ggsave(str_c('interseptions_mi_pi_rib_tr/Number of records with interseptions in mi_pi_rib_tRNA  database for mi_pi_tRNA.tiff'))

mipitrrib_data = d_chr
mipitrrib_data_no_inters = filter(mipitrrib_data,interseptions == 0)
table(mipitrrib_data_no_inters$type)
write.xlsx(mipitrrib_data_no_inters,'smRNA_DATABASE_filter_miR_pi_rib_tRNA.xlsx')

mipitrrib_data_no_inters_check = table(foreach(i = 1:nrow(mipitrrib_data_no_inters),.combine = 'c') %do% {
  if (i %% 100000 == 0) {print(i)}
  str_c(mipitrrib_data_no_inters$Annotation_chr[i],'_',
        mipitrrib_data_no_inters$Annotation_start[i]:mipitrrib_data_no_inters$Annotation_end[i])
})

#####
smRNA_DATABASE$x = 1:nrow(smRNA_DATABASE)
inter_df_by_types = cbind.data.frame(inter_df_by_types,
foreach(n = c(2),.combine = 'cbind') %do% {
  print(names(positions_smRNA_DATABASE_by_types[n]))
  f = positions_smRNA_DATABASE_by_types[[n]]
  
  foreach(i = 1:length(positions_smRNA_DATABASE),.combine = 'c') %do% {
    if (i %% 10000 == 0) {print(i)}
    mean(positions_smRNA_DATABASE[[i]] %in% 
           filter(f,chr == strsplit(positions_smRNA_DATABASE[[i]][1],'_')[[1]][1])$coord)
  }
})
colnames(inter_df_by_types) = c("tRNA","tRNA-derived",
                                "mature_miRNA","precursor_miRNA","rRNA")
inter_df_by_types$type = smRNA_DATABASE$type

inter_df_by_types$type_inter = foreach(i = 1:nrow(inter_df_by_types),.combine = 'c') %do% {
  if (i %% 100000 == 0) {print(i)}
  str_c(names(inter_df_by_types)[c(1:5,9)]
        [which(inter_df_by_types[i,c(1:5,9)] > 0)],collapse = '&')
}

f = unique(unlist(positions_smRNA_DATABASE[which(smRNA_DATABASE$ID %in%
                                                   pi_inter_transcripts)]))
f = str_split_fixed(f,'_',2) 
f = cbind.data.frame(f,unique(unlist(positions_smRNA_DATABASE[
  which(smRNA_DATABASE$ID %in% pi_inter_transcripts)])))
colnames(f) = c('chr','pos','coord')
f$coord = as.character(f$coord)
f$chr = as.character(f$chr)
f$pos = as.numeric(f$pos)

inter_df_by_types$piRNA = sapply(positions_smRNA_DATABASE, function(x) {
  mean(x %in% filter(f,chr == strsplit(x[1],'_')[[1]][1])$coord)
})
inter_df_by_types$piRNA[which(inter_df_by_types$type == 'piRNA')] = 1
write.xlsx(inter_df_by_types,'interseptions_between_bases.xlsx')

#write bed for bases from smRNA_DATABASE####
for (i in unique(smRNA_DATABASE$type)) {
  write.table(filter(smRNA_DATABASE,type == i)[,c(4:6,2)],str_c(i,'_DATABASE.bed'),
              sep = '\t',row.names = F,col.names = F)
}

inter_list_piRNA = foreach(i = unique(smRNA_DATABASE$type)[-3]) %do% {
  d = read.table(str_c('pi_',i,'.bed'),stringsAsFactors = F)
  unique(data.frame(Loci = str_c(d[,1],':',d[,2],'-',d[,3]),
             ID = d[,4],stringsAsFactors = F))
}
names(inter_list_piRNA) = unique(smRNA_DATABASE$type)[-3]
Num_of_pi_inter_transcripts = foreach(i = 1:length(inter_list_piRNA),.combine = 'c') %do% {
  length(unique(inter_list_piRNA[[i]]$ID))
}
pi_inter_transcripts = unique(foreach(i = 1:length(inter_list_piRNA),.combine = 'c') %do% {
  inter_list_piRNA[[i]]$ID
})

#barplot inter transcripts####
type_inter_bar = foreach(i = unique(smRNA_DATABASE$type),.combine = 'cbind.data.frame') %do% {
  print(i)
  foreach(k = unique(inter_df_by_types$type_inter),.combine = 'c') %do% {
    mean(filter(inter_df_by_types,type == i)$type_inter == k)
  }
}
colnames(type_inter_bar) = unique(smRNA_DATABASE$type)
rownames(type_inter_bar) = unique(inter_df_by_types$type_inter)
type_inter_bar = melt(t(type_inter_bar))
type_inter_bar$X2 = as.character(type_inter_bar$X2)
type_inter_bar$X1 = as.character(type_inter_bar$X1)

type_inter_bar$X2[str_detect(type_inter_bar$X2,'&.+&')] = '>2databases'
ggplot(type_inter_bar) + aes(X1,value,fill = X2) + theme_classic() +
  geom_bar(stat = "identity") + ggtitle('Intersections between databases') +
  scale_fill_manual(values = c('black',colorRampPalette(c('lightblue','red',
                                                          'green', "yellow"))(12)),
                    name = NULL) +
  scale_x_discrete(name = 'small RNA type')
ggsave('Intersections_barplot.tiff',width = 10)

#Venn inter transcripts####
inter_df_by_types$x = smRNA_DATABASE$ID
venn.diagram(
  list(mature_miRNA = filter(inter_df_by_types,type %in% c('mature_miRNA','precursor_miRNA'),
      mature_miRNA > 0)$x,
       precursor_miRNA = filter(inter_df_by_types,type %in% c('mature_miRNA','precursor_miRNA'),
                                mature_miRNA > 0)$x),
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'Venn_microRNA_types_transcripts_id.tiff'
             ,category.names = c('mature_miRNA','precursor_miRNA'),
              cat.col = c('red','blue'),
             main = 'Transcripts loci of types of human microRNA', 
             fill = c('red','blue'))
venn.diagram(
  list(tRNA = filter(inter_df_by_types,type %in% c('tRNA','tRNA-derived'),
                     tRNA > 0)$x,
       tRNAderived = filter(inter_df_by_types,type %in% c('tRNA','tRNA-derived'),
                            `tRNA-derived` > 0)$x),
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'Venn_tRNA_types_transcripts.tiff'
             ,category.names = c('tRNA','tRNA-derived'),
              cat.col = c('red','blue'),
             main = 'Transcripts loci of types of human tRNA', 
             fill = c('red','blue'))
venn.diagram(
  list(rRNA = filter(inter_df_by_types,type %in% c('rRNA','tRNA-derived','precursor_miRNA'),
                     rRNA > 0)$x,
       tRNAderived = filter(inter_df_by_types,type %in% c('rRNA','tRNA-derived','precursor_miRNA'),
           `tRNA-derived` > 0)$x,
       precursor_miRNA = filter(inter_df_by_types,type %in% c('mature_miRNA','precursor_miRNA'),
                                mature_miRNA > 0)$x),
             width = 4500, height = 4000,imagetype = 'tiff',
             filename = 'Venn_tsRNA_rRNA_premiRNA_transcripts.tiff'
             ,category.names = c('rRNA','tRNA-derived','precursor_miRNA'),
              cat.col = c('red','blue','green'), cat.pos = 0,
             main = 'Transcripts loci of human tRNA-derived, rRNA, miRNA', 
             fill = c('red','blue','green'))

inter_list_piRNA_id = foreach(i = inter_list_piRNA) %do% {
  unique(i$ID)
}
names(inter_list_piRNA_id) = names(inter_list_piRNA)

venn.diagram(inter_list_piRNA_id,
  width = 4500, height = 4000,imagetype = 'tiff',
  filename = 'Venn_piRNA_intersected_transcripts.tiff'
  ,category.names = names(inter_list_piRNA_id),
  cat.col = c('red','blue','green','violet','grey'), cat.pos = c(0,0,180,180,0),
  main = 'Transcripts of piRNA intersected other types', 
  fill = c('red','blue','green','violet','grey'))

#add hsa-let to DATABASE####

hsalet_chr_start_end = 
  foreach(i = 1:nrow(databases$miRBase[str_detect(databases$miRBase$ID,'hsa-let'),])) %do% {
            str_c(databases$miRBase[str_detect(databases$miRBase$ID,'hsa-let'),]$Loci_chr[i],'_',
                  databases$miRBase[str_detect(databases$miRBase$ID,'hsa-let'),]$Loci_start[i]:
                  databases$miRBase[str_detect(databases$miRBase$ID,'hsa-let'),]$Loci_end[i])          
  }
names(hsalet_chr_start_end) = databases$miRBase$ID[str_detect(databases$miRBase$ID,'hsa-let')]
positions_smRNA_DATABASE_unique = unlist(positions_smRNA_DATABASE_by_types)
hsalet_chr_start_end = hsalet_chr_start_end[
foreach(i = hsalet_chr_start_end,.combine = 'c') %do% {
  mean(i %in% positions_smRNA_DATABASE_unique) == 0
}]
hsalet = filter(databases$miRBase,ID %in% names(hsalet_chr_start_end))
hsalet$source = is.na(hsalet$Fasta_seq)
hsalet$mature = hsalet$Loci_len < 50
smRNA_DATABASE = rbind.data.frame(
  smRNA_DATABASE,
  foreach(i = 1:nrow(hsalet),.combine = 'rbind') %do% {
    if (hsalet$mature[i] == T) {type = 'mature_miRNA'}
    else {type = 'precursor_miRNA'}
    if (hsalet$source[i] == T) {
      Source = 'annotation'
      if (hsalet$strand[i] == '-'){Sequence_consensus = hsalet$Loci_seq[i]}
      else {Sequence_consensus = hsalet$Loci_rev_seq[i]}
      }
    else {
      Source = 'annotation&fasta'
      Sequence_consensus = hsalet$Fasta_seq[i]}
    data.frame(
      type = type,
      ID = hsalet$ID[i],Annotation_consensus = hsalet$Loci[i],
      Annotation_chr = hsalet$Loci_chr[i],Annotation_start = hsalet$Loci_start[i],
      Annotation_end = hsalet$Loci_end[i],Sequence_consensus = Sequence_consensus,
      Length_consensus = hsalet$Loci_len[i],mapped2loci = hsalet$mapped2loci[i],
      strand = hsalet$strand[i],strand_from_annotation = T,
      Source = Source,repeats = 1,interseptions = 0,
      x = nrow(smRNA_DATABASE)+i,stringsAsFactors = F
    )
  }
  )


#write fasta smRNA_database#####
mipitrrib_data_no_inters = read.xlsx('smRNA_DATABASE_filter_miR_pi_rib_tRNA.xlsx')

mipitrrib_data_no_inters_unique = unique(mipitrrib_data_no_inters[,c(2,7)])

mipitrrib_data_no_inters_unique_CCA_tRNA = 
  filter(mipitrrib_data_no_inters_unique,startsWith(ID,'tR'))
mipitrrib_data_no_inters_unique_CCA_tRNA$ID = str_c(
  mipitrrib_data_no_inters_unique_CCA_tRNA$ID,'-CCA')
mipitrrib_data_no_inters_unique_CCA_tRNA$Sequence_consensus = str_c(
  mipitrrib_data_no_inters_unique_CCA_tRNA$Sequence_consensus,'CCA')

mipitrrib_data_no_inters_unique = rbind.data.frame(mipitrrib_data_no_inters_unique,
                                                   mipitrrib_data_no_inters_unique_CCA_tRNA)
mipitrrib_data_no_inters_fasta = foreach(i = 1:nrow(mipitrrib_data_no_inters_unique),
                                         .combine = 'rbind') %do% {
  d = rbind.data.frame(str_c('>',mipitrrib_data_no_inters_unique$ID[i]),
                       mipitrrib_data_no_inters_unique$Sequence_consensus[i])
  colnames(d) = 'f'
  d
                                         }

write.table(mipitrrib_data_no_inters_fasta,
            'smRNA_DATABASE_filter_miR_pi_rib_tRNA_CCA_added.fasta',
            col.names = F,row.names = F)

which(table(mipitrrib_data_no_inters_unique$ID) >1)

