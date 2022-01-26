## Libraries

```{r}
library(tximportData)
library(tximport)
library(DESeq2)
library(stringr)
library(tidyverse)
library(limma)
library(pheatmap)
library(edgeR)
library(foreach)
library(doParallel)
library(openxlsx)
library(VennDiagram)
libs = c('tximportData',
'tximport',
'DESeq2',
'stringr',
'tidyverse',
'limma',
'pheatmap',
'edgeR',
'foreach',
'doParallel',
'openxlsx')

cl <- makeCluster(16)
registerDoParallel(cl) 

```

## Sports

```{r prepare counts data}

Ingerslev_Sports_default = data.frame()
for (i in 0:23) {
  d = read.table(str_c('sports_output and counts/all_counts/Ingerslev_all_counts/un_sra_data',i,'_summary.txt'),header = T, stringsAsFactors = F)
  Ingerslev_Sports_default[d$Sub_Class,as.character(i)] = d$Reads
}
Ingerslev_Sports_default = Ingerslev_Sports_default[-(1:3),]
Ingerslev_Sports_default_smRNA = Ingerslev_Sports_default[!(endsWith(
  row.names(Ingerslev_Sports_default), 'end'
)),]
Ingerslev_Sports_default_fragments = Ingerslev_Sports_default[endsWith(
  row.names(Ingerslev_Sports_default), 'end'
),]

write.xlsx(list(Ingerslev_Sports_default_smRNA,Ingerslev_Sports_default_fragments),'sports_output and counts/Ingerslev_Sports.xlsx', rowNames = T)

Donkin_Sports_default = data.frame()
for (i in 0:22) {
  d = read.table(str_c('sports_output and counts/all_counts/Donkin_all_counts/',i,'_summary.txt'),header = T, stringsAsFactors = F)
  Donkin_Sports_default[d$Sub_Class,as.character(i)] = d$Reads
}
Donkin_Sports_default = Donkin_Sports_default[-(1:3),]
Donkin_Sports_default_smRNA = Donkin_Sports_default[!(endsWith(
  row.names(Donkin_Sports_default), 'end'
)),]
Donkin_Sports_default_fragments = Donkin_Sports_default[endsWith(
  row.names(Donkin_Sports_default), 'end'
),]

write.xlsx(list(Donkin_Sports_default_smRNA,Donkin_Sports_default_fragments),'sports_output and counts/Donkin_Sports.xlsx', rowNames = T)

Hua_Sports_default = data.frame()
for (i in 902:988) {
  d = read.table(str_c('sports_output and counts/all_counts_our_database/Hua_all_counts/un_SRR8543',i,'_summary.txt'),header = T, stringsAsFactors = F)
  Hua_Sports_default[d$Sub_Class,as.character(i)] = d$Reads
}
Hua_Sports_default = Hua_Sports_default[-(1:3),]
Hua_Sports_default_smRNA = Hua_Sports_default[!(endsWith(
  row.names(Hua_Sports_default), 'end'
)),]
Hua_Sports_default_fragments = Hua_Sports_default[endsWith(
  row.names(Hua_Sports_default), 'end'
),]

write.xlsx(list(Hua_Sports_default_smRNA,Hua_Sports_default_fragments),'sports_output and counts/Hua_Sports.xlsx', rowNames = T)

Ingerslev_Sports_default_smRNA_filter =
  na.omit(Ingerslev_Sports_default_smRNA[rowMeans(Ingerslev_Sports_default_smRNA) >= 10,])
Donkin_Sports_default_smRNA_filter = 
  na.omit(Donkin_Sports_default_smRNA[rowMeans(Donkin_Sports_default_smRNA) >= 10,])
Hua_Sports_default_smRNA_filter = 
  na.omit(Hua_Sports_default_smRNA[rowMeans(Hua_Sports_default_smRNA) >= 10,])
Ingerslev_Sports_default_fragments_filter =
  na.omit(Ingerslev_Sports_default_fragments[rowMeans(Ingerslev_Sports_default_fragments) >= 10,])
Donkin_Sports_default_fragments_filter = 
  na.omit(Donkin_Sports_default_fragments[rowMeans(Donkin_Sports_default_fragments) >= 10,])
Hua_Sports_default_fragments_filter = 
  na.omit(Hua_Sports_default_fragments[rowMeans(Hua_Sports_default_fragments) >= 10,])
Sports_default_smRNA_filter_list = list(
  Ingerslev = Ingerslev_Sports_default_smRNA_filter,
  Donkin = Donkin_Sports_default_smRNA_filter,
  Hua = Hua_Sports_default_smRNA_filter
)
Sports_default_fragments_filter_list = list(
  Ingerslev = Ingerslev_Sports_default_fragments_filter,
  Donkin = Donkin_Sports_default_fragments_filter,
  Hua = Hua_Sports_default_fragments_filter
)

```

```{r prepare covariates}

Covariates_list = list(
  Ingerslev = data.frame(cov = c(rep(c('Untrained','Trained','Detrained'),2),
                rep(c('Untrained','Trained'),2),
                rep(c('Untrained','Trained','Detrained'),1),
                rep(c('Untrained','Trained'),1),
                rep(c('Untrained','Trained','Detrained'),3)),
                row.names = colnames(Ingerslev_Sports_default_filter)),
  Donkin = data.frame(cov = c(rep('Lean',13),rep('Obese',10)),
                      row.names = colnames(Donkin_Sports_default_filter)),
  Hua = data.frame(cov = c(rep('H-GQE',23),rep('L-GQE',64)),
                   row.names = colnames(Hua_Sports_default_filter))
)
write.xlsx(Covariates_list,'Covariates_literature.xlsx',rowNames = T)

```

```{r DESeq2 and heatmap}
for (i in 1:3) {
  
dds = DESeqDataSetFromMatrix(round(Sports_default_smRNA_filter_list[[i]]),
        colData = Covariates_list[[i]],design = ~cov)
dds = DESeq(dds)
res_df = results(dds)
res_df = as.data.frame(res_df)
res_df = arrange(res_df,pvalue)
write.xlsx(res_df, str_c('sports_output and counts/',names(Sports_default_smRNA_filter_list)[i],'_Sports_filtered_results_smRNA.xlsx'),row.names = T)

#row.names(filter(res_df,pvalue < 0.01))[which(row.names(filter(res_df,pvalue < 0.01)) %!in% str_remove_all(str_replace_all(row.names(Literature_Rsubread$Ingerslev),'miR','mir'),'-[[:digit:]]p'))]

ggsave(str_c('sports_output and counts/',names(Sports_default_smRNA_filter_list)[i],'_heatmap_DESeq2_counts_Sports_pval005_smRNA.tiff'),
       pheatmap(assay(varianceStabilizingTransformation(dds))[
  row.names(filter(res_df,pvalue< 0.05)),],
  annotation_col = Covariates_list[[i]],
  main = 'VST counts for pvalue < 0.05 for Sports with our database'))

ggplot(data.frame(Name = rep(row.names(filter(res_df,pvalue< 0.05)),each = nrow(Covariates_list[[i]])),
                    Counts = as.double(t(assay(varianceStabilizingTransformation(dds))[
  row.names(filter(res_df,pvalue< 0.05)),])),
  Groups = rep(Covariates_list[[i]]$cov,nrow(filter(res_df,pvalue< 0.05))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() +
    theme_classic() + ggtitle('VST counts for Sports with our database for pvalue < 0.05') + theme(axis.text.x = element_text(angle = 90))
ggsave(str_c('sports_output and counts/',names(Sports_default_smRNA_filter_list)[i],'_boxplot_DESeq2_filtered_counts_Sports_DESeq2_smRNA.tiff'),width = 12,height = 10)
}

dds = DESeqDataSetFromMatrix(round(Sports_default_fragments_filter_list[[i]]),
        colData = Covariates_list[[i]],design = ~cov)
dds = DESeq(dds)
res_df = results(dds)
res_df = as.data.frame(res_df)
res_df = arrange(res_df,pvalue)
write.xlsx(res_df, str_c('sports_output and counts/',names(Sports_default_fragments_filter_list)[i],'_Sports_our_database_filtered_results_fragments.xlsx'),row.names = T)

ggsave(str_c('sports_output and counts/',names(Sports_default_fragments_filter_list)[i],'_heatmap_DESeq2_counts_Sports_pval001_fragments.tiff'),
       pheatmap(assay(varianceStabilizingTransformation(dds))[
  row.names(filter(res_df,pvalue< 0.01)),],
  annotation_col = Covariates_list[[i]],
  main = 'VST counts for pvalue < 0.01 for Sports with our database'))

ggplot(data.frame(Name = rep(row.names(filter(res_df,pvalue< 0.01)),each = nrow(Covariates_list[[i]])),
                    Counts = as.double(t(assay(varianceStabilizingTransformation(dds))[
  row.names(filter(res_df,pvalue< 0.01)),])),
  Groups = rep(Covariates_list[[i]]$cov,nrow(filter(res_df,pvalue< 0.01))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() +
    theme_classic() + ggtitle('VST counts for Sports with our database for pvalue < 0.01') + theme(axis.text.x = element_text(angle = 90))
ggsave(str_c('sports_output and counts/',names(Sports_default_fragments_filter_list)[i],'_boxplot_DESeq2_filtered_counts_Sports_DESeq2_fragments.tiff'),width = 12,height = 10)
}

```

## MINT

```{r prepare counts data}

MINT_gene_list = read.table('Human/abundanceMINT.tsv',header = T,stringsAsFactors = F)[,1]

Samples_literature = list(
  Ingerslev = str_c('filt_tRNA_un_sra_data',colnames(Sports_default_filter_list$Ingerslev)),
  Donkin = str_c('filt_tRNA_',colnames(Sports_default_filter_list$Donkin)),
  Hua = str_c('filt_tRNA_un_SRR8543',colnames(Sports_default_filter_list$Hua))
)

Literature_MINT_kallisto = foreach(i = names(Samples_literature)) %do% {
  files = file.path('Human',str_c(i,'_kallisto'),str_c(Samples_literature[[i]],'_MINT_11'),"abundance.h5")
  names(files) = Samples_literature[[i]]
  txi <- tximport(files, type = "kallisto", 
                             tx2gene = cbind.data.frame(MINT_gene_list,
                                                        MINT_gene_list), 
                             ignoreAfterBar = TRUE)
  txi
}
names(Literature_MINT_kallisto) = names(Samples_literature)

write.xlsx(foreach(i = Human_txi_kallisto_kmers) %do% {i$counts},
           'Human/Literature_MINT_kallisto_counts.xlsx',row.names = T)

Literature_MINT_kallisto_filtered = foreach(i = Literature_MINT_kallisto) %do% {
  i$abundance = i$abundance[rowMeans(i$counts) >= 10,]
  i$length = i$length[rowMeans(i$counts) >= 10,]
  i$counts = i$counts[rowMeans(i$counts) >= 10,]
  i
}
names(Literature_MINT_kallisto_filtered) = names(Literature_MINT_kallisto)
write.xlsx(foreach(i = Literature_MINT_kallisto_filtered) %do% {i$counts},
           'Human/Literature_MINT_kallisto_counts_filtered.xlsx',row.names = T)

```

```{r DESeq2 and heatmap}

for (i in 1:3) {
  
dds = DESeqDataSetFromTximport(txi = Literature_MINT_kallisto_filtered[[i]],
        colData = Covariates_list[[i]],design = ~cov)
dds = DESeq(dds)
res_df = results(dds)
res_df = as.data.frame(res_df)
res_df = arrange(res_df,padj)
write.xlsx(res_df, str_c('Human/kallisto_DESeq2_MINT/',names(Literature_MINT_kallisto_filtered)[i],'_MINT_filtered_results.xlsx'),row.names = T)

ggsave(str_c('Human/kallisto_DESeq2_MINT/',names(Literature_MINT_kallisto_filtered)[i],'_heatmap_DESeq2_counts_MINT_padj005.tiff'),
       pheatmap(assay(varianceStabilizingTransformation(dds))[
  row.names(filter(res_df,padj< 0.05)),],
  annotation_col = Covariates_list[[i]],
  main = 'VST counts for padj < 0.05 for Sports'),width = 7,height = 10)

ggplot(data.frame(Name = rep(row.names(filter(res_df,padj < 0.05)),each = nrow(Covariates_list[[i]])),
                    Counts = as.double(t(assay(varianceStabilizingTransformation(dds))[
  row.names(filter(res_df,padj < 0.05)),])),
  Groups = rep(Covariates_list[[i]]$cov,nrow(filter(res_df,padj < 0.05))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() +
    theme_classic() + ggtitle('VST counts for Sports for padj < 0.05') + theme(axis.text.x = element_text(angle = 90))
ggsave(str_c('Human/kallisto_DESeq2_MINT/',names(Literature_MINT_kallisto_filtered)[i],'_boxplot_DESeq2_filtered_counts_MINT_DESeq2.tiff'),width = 12,height = 10)
}

```

## Rsubread

```{r prepare counts data}

Literature_Rsubread = list(
  Ingerslev = read.table('Ingerslev_Rsubread_allRNA.csv',stringsAsFactors = F,header = T),
  Donkin = read.table('Donkin_Rsubread_allRNA.csv',stringsAsFactors = F,header = T),
  Hua = read.table('Hua_Rsubread_allRNA.csv',stringsAsFactors = F,header = T)
)
for (i in 1:3) {
  row.names(Literature_Rsubread[[i]]) = Literature_Rsubread[[i]]$Gene
  Literature_Rsubread[[i]] = Literature_Rsubread[[i]][,-1]
}
Literature_Rsubread$Ingerslev = Literature_Rsubread$Ingerslev[,c(24,11,16:23,1:10,12:15)]
write.xlsx(Literature_Rsubread,'Human/Rsubread_literature/Rsubread_counts.xlsx')

```

```{r DESeq2 and heatmap}

for (i in 1:3) {
  
dds = DESeqDataSetFromMatrix(Literature_Rsubread[[i]],
        colData = Covariates_list[[i]],design = ~cov)
dds = DESeq(dds)
res_df = results(dds)
res_df = as.data.frame(res_df)
res_df = arrange(res_df,pvalue)
write.xlsx(res_df, str_c('Human/Rsubread_literature/',names(Sports_default_filter_list)[i],'_Rsubread_filtered_results.xlsx'),row.names = T)

ggsave(str_c('Human/Rsubread_literature/',names(Sports_default_filter_list)[i],'_heatmap_DESeq2_counts_Rsubread_pval001.tiff'),
       pheatmap(assay(varianceStabilizingTransformation(dds))[
  row.names(filter(res_df,pvalue< 0.01)),],
  annotation_col = Covariates_list[[i]],
  main = 'VST counts for pvalue < 0.01 for Sports'))

ggplot(data.frame(Name = rep(row.names(filter(res_df,pvalue< 0.01)),each = nrow(Covariates_list[[i]])),
                    Counts = as.double(t(assay(varianceStabilizingTransformation(dds))[
  row.names(filter(res_df,pvalue< 0.01)),])),
  Groups = rep(Covariates_list[[i]]$cov,nrow(filter(res_df,pvalue< 0.01))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() +
    theme_classic() + ggtitle('VST counts for Sports for pvalue < 0.01') + theme(axis.text.x = element_text(angle = 90))
ggsave(str_c('Human/Rsubread_literature/',names(Sports_default_filter_list)[i],'_boxplot_DESeq2_filtered_counts_Rsubread_DESeq2.tiff'),width = 10,height = 7)
}

```


## Pictures

```{r Donkin smRNA}

Donkin_res_smRNA = list(Sports_def = read.xlsx('sports_output and counts/Donkin_Sports_filtered_results_smRNA.xlsx',rowNames = T),
                  Sports_our = read.xlsx('sports_output and counts/Donkin_Sports_our_database_filtered_results_smRNA.xlsx',rowNames = T),
                  Rsubread = read.xlsx('Human/Rsubread_literature/Donkin_Rsubread_filtered_results.xlsx',rowNames = T))

ggsave('D_Sdef_h.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(na.omit(round(read.xlsx('sports_output and counts/Donkin_Sports_default.xlsx',rowNames = T)))))[
  row.names(filter(Donkin_res_smRNA$Sports_def,pvalue< 0.05)),],
  annotation_col = Covariates_list$Donkin,
  main = 'Heatmap with vst counts for Donkin et.al\npvalue < 0.05, Sports default analysis'))

ggplot(data.frame(Name = rep(row.names(filter(Donkin_res_smRNA$Sports_def,pvalue< 0.05)),each = nrow(Covariates_list$Donkin)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(na.omit(round(read.xlsx('sports_output and counts/Donkin_Sports_default.xlsx',rowNames = T)))))[
  row.names(filter(Donkin_res_smRNA$Sports_def,pvalue< 0.05)),])),
  Groups = rep(Covariates_list$Donkin$cov,nrow(filter(Donkin_res_smRNA$Sports_def,pvalue< 0.05))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = 
      str_remove(sort(row.names(filter(Donkin_res_smRNA$Sports_def,pvalue< 0.05))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 12, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('D_Sdef_b.tiff')

ggsave('D_Sour_h.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(na.omit(round(read.xlsx('sports_output and counts/Donkin_Sports_our_database.xlsx',rowNames = T)))))[
  row.names(filter(Donkin_res_smRNA$Sports_our,pvalue< 0.05)),],
  annotation_col = Covariates_list$Donkin,
  main = 'Heatmap with vst counts for Donkin et.al\npvalue < 0.05, Sports analysis with our database'))

ggplot(data.frame(Name = rep(row.names(filter(Donkin_res_smRNA$Sports_our,pvalue< 0.05)),each = nrow(Covariates_list$Donkin)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(na.omit(round(read.xlsx('sports_output and counts/Donkin_Sports_our_database.xlsx',rowNames = T)))))[
  row.names(filter(Donkin_res_smRNA$Sports_our,pvalue< 0.05)),])),
  Groups = rep(Covariates_list$Donkin$cov,nrow(filter(Donkin_res_smRNA$Sports_our,pvalue< 0.05))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = str_remove(sort(row.names(filter(Donkin_res_smRNA$Sports_our,pvalue< 0.05))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 12, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('D_Sour_b.tiff')

ggsave('D_Rs_h.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(Literature_Rsubread$Donkin))[
  row.names(filter(Donkin_res_smRNA$Rsubread,pvalue< 0.01)),],
  annotation_col = Covariates_list$Donkin,
  main = 'Heatmap with vst counts for Donkin et.al\npvalue < 0.01, Rsubread analysis with our database'))

ggplot(data.frame(Name = rep(row.names(filter(Donkin_res_smRNA$Rsubread,pvalue< 0.01)),each = nrow(Covariates_list$Donkin)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(Literature_Rsubread$Donkin))[
  row.names(filter(Donkin_res_smRNA$Rsubread,pvalue< 0.01)),])),
  Groups = rep(Covariates_list$Donkin$cov,nrow(filter(Donkin_res_smRNA$Rsubread,pvalue< 0.01))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = str_remove(sort(row.names(filter(Donkin_res_smRNA$Rsubread,pvalue< 0.01))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 12, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('D_Rs_b.tiff')



```

```{r Donkin fragments}

Donkin_res_fragments = list(Sports_def = read.xlsx('sports_output and counts/Donkin_Sports_filtered_results_fragments.xlsx',rowNames = T),
                  Sports_our = read.xlsx('Human/kallisto_DESeq2_MINT/Donkin_MINT_filtered_results.xlsx',rowNames = T))

ggsave('D_Sdef_h_ts.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(na.omit(round(Donkin_Sports_default_fragments))))[
  row.names(filter(Donkin_res_fragments$Sports_def,padj< 0.1)),],
  annotation_col = Covariates_list$Donkin,
  main = 'Heatmap with vst counts for Donkin et.al\npadj < 0.1, Sports default analysis'))

ggplot(data.frame(Name = rep(row.names(filter(Donkin_res_fragments$Sports_def,padj< 0.1)),each = nrow(Covariates_list$Donkin)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(na.omit(round(Donkin_Sports_default_fragments))))[
  row.names(filter(Donkin_res_fragments$Sports_def,padj< 0.1)),])),
  Groups = rep(Covariates_list$Donkin$cov,nrow(filter(Donkin_res_fragments$Sports_def,padj< 0.1))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = 
      str_remove(sort(row.names(filter(Donkin_res_fragments$Sports_def,padj< 0.1))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 12, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('D_Sdef_b_ts.tiff')

ggsave('D_Sour_h_ts.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(na.omit(round(Literature_MINT_kallisto_filtered$Donkin$counts))))[
  row.names(filter(Donkin_res_fragments$Sports_our,padj < 0.01)),],
  annotation_col = Covariates_list$Donkin,
  main = 'Heatmap with vst counts for Donkin et.al\npadj < 0.01, Analysis with MINTbase'))

ggplot(data.frame(Name = rep(row.names(filter(Donkin_res_fragments$Sports_our,padj < 0.01)),each = nrow(Covariates_list$Donkin)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(na.omit(round(Literature_MINT_kallisto_filtered$Donkin$counts))))[
  row.names(filter(Donkin_res_fragments$Sports_our,padj < 0.01)),])),
  Groups = rep(Covariates_list$Donkin$cov,nrow(filter(Donkin_res_fragments$Sports_our,padj < 0.01))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = str_remove(sort(row.names(filter(Donkin_res_fragments$Sports_our,padj < 0.01))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 12, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('D_Sour_b_ts.tiff')

```

```{r Ingerslev smRNA}

Ingerslev_res_smRNA = list(Sports_def = read.xlsx('sports_output and counts/Ingerslev_Sports_filtered_results_smRNA.xlsx',rowNames = T),
                  Sports_our = read.xlsx('sports_output and counts/Ingerslev_Sports_our_database_filtered_results_smRNA.xlsx',rowNames = T),
                  Rsubread = read.xlsx('Human/Rsubread_literature/Ingerslev_Rsubread_filtered_results.xlsx',rowNames = T))

ggsave('I_Sdef_h.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(na.omit(round(read.xlsx('sports_output and counts/Ingerslev_Sports_default.xlsx',rowNames = T)))))[
  row.names(filter(Ingerslev_res_smRNA$Sports_def,pvalue< 0.01)),],
  annotation_col = Covariates_list$Ingerslev,
  main = 'Heatmap with vst counts for Ingerslev et.al\npvalue < 0.01, Sports default analysis'))

ggplot(data.frame(Name = rep(row.names(filter(Ingerslev_res_smRNA$Sports_def,pvalue< 0.01)),each = nrow(Covariates_list$Ingerslev)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(na.omit(round(read.xlsx('sports_output and counts/Ingerslev_Sports_default.xlsx',rowNames = T)))))[
  row.names(filter(Ingerslev_res_smRNA$Sports_def,pvalue< 0.01)),])),
  Groups = rep(Covariates_list$Ingerslev$cov,nrow(filter(Ingerslev_res_smRNA$Sports_def,pvalue< 0.01))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = 
      str_remove(sort(row.names(filter(Ingerslev_res_smRNA$Sports_def,pvalue< 0.01))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 12, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('I_Sdef_b.tiff')

ggsave('I_Sour_h.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(na.omit(round(read.xlsx('sports_output and counts/Ingerslev_Sports_our_database.xlsx',rowNames = T)))))[
  row.names(filter(Ingerslev_res_smRNA$Sports_our,pvalue< 0.01)),],
  annotation_col = Covariates_list$Ingerslev,
  main = 'Heatmap with vst counts for Ingerslev et.al\npvalue < 0.01, Sports analysis with our database'))

ggplot(data.frame(Name = rep(row.names(filter(Ingerslev_res_smRNA$Sports_our,pvalue< 0.01)),each = nrow(Covariates_list$Ingerslev)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(na.omit(round(read.xlsx('sports_output and counts/Ingerslev_Sports_our_database.xlsx',rowNames = T)))))[
  row.names(filter(Ingerslev_res_smRNA$Sports_our,pvalue< 0.01)),])),
  Groups = rep(Covariates_list$Ingerslev$cov,nrow(filter(Ingerslev_res_smRNA$Sports_our,pvalue< 0.01))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = str_remove(sort(row.names(filter(Ingerslev_res_smRNA$Sports_our,pvalue< 0.01))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 12, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('I_Sour_b.tiff')

ggsave('I_Rs_h.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(Literature_Rsubread$Ingerslev))[
  row.names(filter(Ingerslev_res_smRNA$Rsubread,pvalue< 0.005)),],
  annotation_col = Covariates_list$Ingerslev,
  main = 'Heatmap with vst counts for Ingerslev et.al\npvalue < 0.005, Rsubread analysis with our database'))

ggplot(data.frame(Name = rep(row.names(filter(Ingerslev_res_smRNA$Rsubread,pvalue< 0.005)),each = nrow(Covariates_list$Ingerslev)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(Literature_Rsubread$Ingerslev))[
  row.names(filter(Ingerslev_res_smRNA$Rsubread,pvalue< 0.005)),])),
  Groups = rep(Covariates_list$Ingerslev$cov,nrow(filter(Ingerslev_res_smRNA$Rsubread,pvalue< 0.005))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = str_remove(sort(row.names(filter(Ingerslev_res_smRNA$Rsubread,pvalue< 0.005))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 12, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('I_Rs_b.tiff')



```

```{r Ingerslev fragments}

Ingerslev_res_fragments = list(Sports_def = read.xlsx('sports_output and counts/Ingerslev_Sports_filtered_results_fragments.xlsx',rowNames = T),
                  Sports_our = read.xlsx('Human/kallisto_DESeq2_MINT/Ingerslev_MINT_filtered_results.xlsx',rowNames = T))

ggsave('I_Sdef_h_ts.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(na.omit(round(Ingerslev_Sports_default_fragments))))[
  row.names(filter(Ingerslev_res_fragments$Sports_def,padj< 0.1)),],
  annotation_col = Covariates_list$Ingerslev,
  main = 'Heatmap with vst counts for Ingerslev et.al\npadj < 0.1, Sports default analysis'))

ggplot(data.frame(Name = rep(row.names(filter(Ingerslev_res_fragments$Sports_def,padj< 0.1)),each = nrow(Covariates_list$Ingerslev)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(na.omit(round(Ingerslev_Sports_default_fragments))))[
  row.names(filter(Ingerslev_res_fragments$Sports_def,padj< 0.1)),])),
  Groups = rep(Covariates_list$Ingerslev$cov,nrow(filter(Ingerslev_res_fragments$Sports_def,padj< 0.1))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = 
      str_remove(sort(row.names(filter(Ingerslev_res_fragments$Sports_def,padj< 0.1))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 12, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('I_Sdef_b_ts.tiff')

ggsave('I_Sour_h_ts.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(na.omit(round(Literature_MINT_kallisto_filtered$Ingerslev$counts))))[
  row.names(filter(Ingerslev_res_fragments$Sports_our,padj < 0.01)),],
  annotation_col = Covariates_list$Ingerslev,
  main = 'Heatmap with vst counts for Ingerslev et.al\npadj < 0.01, Analysis with MINTbase'))

ggplot(data.frame(Name = rep(row.names(filter(Ingerslev_res_fragments$Sports_our,padj < 0.01)),each = nrow(Covariates_list$Ingerslev)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(na.omit(round(Literature_MINT_kallisto_filtered$Ingerslev$counts))))[
  row.names(filter(Ingerslev_res_fragments$Sports_our,padj < 0.01)),])),
  Groups = rep(Covariates_list$Ingerslev$cov,nrow(filter(Ingerslev_res_fragments$Sports_our,padj < 0.01))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = str_remove(sort(row.names(filter(Ingerslev_res_fragments$Sports_our,padj < 0.01))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 8, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('I_Sour_b_ts.tiff')

```

```{r Hua smRNA}

Hua_res_smRNA = list(Sports_def = read.xlsx('sports_output and counts/Hua_Sports_filtered_results_smRNA.xlsx',rowNames = T),
                  Sports_our = read.xlsx('sports_output and counts/Hua_Sports_our_database_filtered_results_smRNA.xlsx',rowNames = T),
                  Rsubread = read.xlsx('Human/Rsubread_literature/Hua_Rsubread_filtered_results.xlsx',rowNames = T))

ggsave('H_Sdef_h.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(na.omit(round(read.xlsx('sports_output and counts/Hua_Sports_default.xlsx',rowNames = T)))))[
  row.names(filter(Hua_res_smRNA$Sports_def,pvalue< 0.01)),],
  annotation_col = Covariates_list$Hua,
  main = 'Heatmap with vst counts for Hua et.al\npvalue < 0.01, Sports default analysis'))

ggplot(data.frame(Name = rep(row.names(filter(Hua_res_smRNA$Sports_def,pvalue< 0.01)),each = nrow(Covariates_list$Hua)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(na.omit(round(read.xlsx('sports_output and counts/Hua_Sports_default.xlsx',rowNames = T)))))[
  row.names(filter(Hua_res_smRNA$Sports_def,pvalue< 0.01)),])),
  Groups = rep(Covariates_list$Hua$cov,nrow(filter(Hua_res_smRNA$Sports_def,pvalue< 0.01))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = 
      str_remove(sort(row.names(filter(Hua_res_smRNA$Sports_def,pvalue< 0.01))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 12, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('H_Sdef_b.tiff')

ggsave('H_Sour_h.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(na.omit(round(read.xlsx('sports_output and counts/Hua_Sports_our_database.xlsx',rowNames = T)))))[
  row.names(filter(Hua_res_smRNA$Sports_our,pvalue< 0.05)),],
  annotation_col = Covariates_list$Hua,
  main = 'Heatmap with vst counts for Hua et.al\npvalue < 0.05, Sports analysis with our database'))

ggplot(data.frame(Name = rep(row.names(filter(Hua_res_smRNA$Sports_our,pvalue< 0.05)),each = nrow(Covariates_list$Hua)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(na.omit(round(read.xlsx('sports_output and counts/Hua_Sports_our_database.xlsx',rowNames = T)))))[
  row.names(filter(Hua_res_smRNA$Sports_our,pvalue< 0.05)),])),
  Groups = rep(Covariates_list$Hua$cov,nrow(filter(Hua_res_smRNA$Sports_our,pvalue< 0.05))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = str_remove(sort(row.names(filter(Hua_res_smRNA$Sports_our,pvalue< 0.05))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 12, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('H_Sour_b.tiff')

ggsave('H_Rs_h.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(Literature_Rsubread$Hua))[
  row.names(filter(Hua_res_smRNA$Rsubread,padj< 0.001)),],
  annotation_col = Covariates_list$Hua,
  main = 'Heatmap with vst counts for Hua et.al\npadj < 0.001, Rsubread analysis with our database'))

ggplot(data.frame(Name = rep(row.names(filter(Hua_res_smRNA$Rsubread,padj< 0.001)),each = nrow(Covariates_list$Hua)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(Literature_Rsubread$Hua))[
  row.names(filter(Hua_res_smRNA$Rsubread,padj< 0.001)),])),
  Groups = rep(Covariates_list$Hua$cov,nrow(filter(Hua_res_smRNA$Rsubread,padj< 0.001))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = str_remove(sort(row.names(filter(Hua_res_smRNA$Rsubread,padj< 0.001))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 12, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('H_Rs_b.tiff')



```

```{r Hua fragments}

Hua_res_fragments = list(
                  Sports_our = read.xlsx('Human/kallisto_DESeq2_MINT/Hua_MINT_filtered_results.xlsx',rowNames = T))


ggsave('H_Sour_h_ts.tiff',
pheatmap(
  varianceStabilizingTransformation(as.matrix(na.omit(round(Literature_MINT_kallisto_filtered$Hua$counts))))[
  row.names(filter(Hua_res_fragments$Sports_our,padj < 0.01)),],
  annotation_col = Covariates_list$Hua,
  main = 'Heatmap with vst counts for Hua et.al\npadj < 0.01, Analysis with MINTbase'))

ggplot(data.frame(Name = rep(row.names(filter(Hua_res_fragments$Sports_our,padj < 0.01)),each = nrow(Covariates_list$Hua)),
                    Counts = as.double(t(varianceStabilizingTransformation(as.matrix(na.omit(round(Literature_MINT_kallisto_filtered$Hua$counts))))[
  row.names(filter(Hua_res_fragments$Sports_our,padj < 0.01)),])),
  Groups = rep(Covariates_list$Hua$cov,nrow(filter(Hua_res_fragments$Sports_our,padj < 0.01))),stringsAsFactors = F)) +
    aes(Name,Counts,fill = Groups) + geom_boxplot() + 
    scale_y_continuous(name = NULL) +
    scale_x_discrete(name = NULL, labels = str_remove(sort(row.names(filter(Hua_res_fragments$Sports_our,padj < 0.01))), fixed('mature-'))) +
    theme_classic() + theme(axis.text.x = element_text(size = 8, face = 'bold',angle = 90),
                            axis.text.y = element_text(size = 14, face = 'bold'))
ggsave('H_Sour_b_ts.tiff')

```

```{r}

print(xtable(filter(Hua_res_smRNA$Rsubread, padj < 0.001), type = "latex"), file = "Hua_res_smRNA_Rsubread.tex")
print(xtable(filter(Hua_res_smRNA$Sports_def, pvalue < 0.01), type = "latex"), file = "Hua_res_smRNA_Sports.tex")
print(xtable(filter(Hua_res_fragments$Sports_our, padj < 0.01), type = "latex"), file = "Hua_res_fragments_MINT.tex")
print(xtable(filter(Ingerslev_res_smRNA$Sports_def, pvalue < 0.01), type = "latex"), file = "Ingerslev_res_smRNA_Sports.tex")
print(xtable(filter(Ingerslev_res_smRNA$Rsubread, pvalue < 0.005), type = "latex"), file = "Ingerslev_res_smRNA_Rsubread.tex")
print(xtable(filter(Ingerslev_res_fragments$Sports_def, padj < 0.1), type = "latex"), file = "Ingerslev_res_fragments_Sports.tex")
print(xtable(filter(Ingerslev_res_fragments$Sports_our, padj < 0.01), type = "latex"), file = "Ingerslev_res_fragments_MINT.tex")
print(xtable(filter(Donkin_res_smRNA$Rsubread, pvalue < 0.01), type = "latex"), file = "Donkin_res_smRNA_Rsubread.tex")
print(xtable(filter(Donkin_res_smRNA$Sports_def, pvalue < 0.05), type = "latex"), file = "Donkin_res_smRNA_Sports.tex")
print(xtable(filter(Donkin_res_fragments$Sports_def, padj < 0.1), type = "latex"), file = "Donkin_res_fragments_Sports.tex")
print(xtable(filter(Donkin_res_fragments$Sports_our, padj < 0.01), type = "latex"), file = "Donkin_res_fragments_MINT.tex")

```


