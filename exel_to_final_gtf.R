library(readxl)
library(magrittr)
library(dplyr)
library(stringr)

"%+%" <- function(...){
  paste0(...)
}

#list of species (human, mouse, celegans, drosophila, rats)
species_list <- read.table("/home/ivan/Desktop/Grant2021/aligment_approach/annot_sum/other_species/species_list.txt", header = FALSE, sep = "\n")
species_list <- split(species_list, seq(nrow(species_list)))

#list of snRNA's types (all, mature_miRNA, piRNA, tRNA, rRNA, tRNA_derived) 
types_list <- read.table("/home/ivan/Desktop/Grant2021/aligment_approach/annot_sum/other_species/types_list.txt", header = FALSE, sep = "\n")
types_list <- split(types_list, seq(nrow(types_list)))

#cycle for converting exel-table to final gtf
for (specie in species_list){
  specie_directory = "/home/ivan/Desktop/Grant2021/aligment_approach/annot_sum/other_species/result/" %+% specie
  dir.create(specie_directory)
  for (type in types_list) {
    specie_type_file = "/home/ivan/Desktop/Grant2021/aligment_approach/annot_sum/other_species/" %+% specie %+% "/" %+% specie %+% "_smRNA_DATABASE_" %+% type %+% ".xlsx"
    output_directory = "/home/ivan/Desktop/Grant2021/aligment_approach/annot_sum/other_species/result/" %+% specie %+% "/" %+% specie %+% "_" %+% type %+% ".gtf"
    t3 <- read_excel(specie_type_file)
    
    colnames(t3)[2] <- "ID"
    colnames(t3)[1] <- "source_db"
    
    t3 <- t3 %>% group_by(ID) %>% mutate(num = row_number())
    t3$feauture <- "exon"
    t3$transcript <- "transcript_id \""
    t3$transcript_copy <- "\"; transcript_copy_id \""
    t3$line <- "_"
    t3$dots_end <- "\"; "
    t3$dots <- "."
    t3$seq_name <- "sequence \""
    t3$quotes <- "\""
    t3$group <- str_c(t3$transcript, t3$ID, t3$transcript_copy, t3$ID, t3$line, t3$num, t3$dots_end, t3$seq_name, t3$Sequence_consensus, t3$quotes)
    
    t3$source_db <- gsub("mature_miRNA", "miRBase_v22_1_mature_miRNA", t3$source_db)
    t3$source_db <- gsub("piRNA", "pirnadb_v1_7_6", t3$source_db)
    t3$source_db <- gsub("tRNA", "GtRNAdb_v2_0", t3$source_db)
    t3$source_db <- gsub("rRNA", "UCSC", t3$source_db)
    
    all_table_gtf <- t3[, c("Annotation_chr", "source_db", "feauture", "Annotation_start", "Annotation_end", "dots", "strand", "dots", "group")]
    write.table(all_table_gtf, file = output_directory, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}
    