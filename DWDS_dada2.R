###########################################
### V4 processing for Kimbell et al, 2020
### Lou LaMartina, finalized April 27, 2020
###########################################


# PURPOSE
# Sequence V4 regions of 16S rRNA genes of microbial
# communities from a decommissioned pipe from a 
# drinking water distributions system (DWDS).

# PCR amplicons from DW pipe were sequenced on 
# Illumina MiSeq. Samples from an anerobic digester
# were sequenced on the same plate. Initial investigations
# of processed data showed evidence that the DW pipe
# samples were contaminated with anaerobic digester sequences.
# Also, there was a fair amount of contamination from the 
# no-template/negative control (extraction & PCR blanks).

# Therefore, following filtering/merging/error-correcting/
# chimera-checking with DADA2, additional steps were added:

# 1. remove exact sequence matches from the mock community
# 2. remove sequences believed to be from the negative control:
#    remove sequences whose mean relative abundance was 3x more 
#    in the control than the rest of the dataset
# 3. remove any ASVs classified as mitochondria, eukaryota,
#    or chloroplast.
# 4. remove ASVs that were likely from the anaerobic digester
#    (ANDIG) dataset: if mean relative abundance of an ASV was 
#    10x greater in DWDS than ANDIG, then it was considered DWDS-derived.


library(dada2)
library(phyloseq)


# set working directory
setwd("~/Desktop/Lab/Projects/Kimbell")




###################################
### prepare data for processing ###
###################################


# set file paths
path <- "./cutadapt"
pathF <- "./cutadapt/fastqF"
pathR <- "./cutadapt/fastqR"


# set paths for filtered reads
filtered_pathF <- "./cutadapt/fastqF/Filtered"
filtered_pathR <- "./cutadapt/fastqR/Filtered"


# sort forward and reverse reads into file paths
fastqFs <- sort(list.files(pathF, pattern = "_R1.fastq.bz2", full.names = TRUE))
fastqRs <- sort(list.files(pathR, pattern = "_R2.fastq.bz2", full.names = TRUE))


# extract file names
Sample_names <- sapply(strsplit(basename(fastqFs), "_"), '[', 1)


# remove "MKE-" from sample names
Sample_names[grepl("MKE", Sample_names) == TRUE] <- sapply(strsplit(Sample_names[grepl("MKE", Sample_names) == TRUE], "-"), '[', 2)




################################
### Inspect sequence quality ###
################################

# visualize quality of reads
plotQualityProfile(fastqFs[1:4])
plotQualityProfile(fastqRs[1:4])


# check: if there are not the same number of F and R files, stop.
if(length(fastqFs) != length(fastqRs)) 
  stop("Forward and reverse files do not match.")




#######################
### Quality control ###
#######################

# give filtered files new names and paths
filtered_Fs <- file.path(filtered_pathF, paste0(Sample_names, "_F_filt.fastq.bz2"))
filtered_Rs <- file.path(filtered_pathR, paste0(Sample_names, "_R_filt.fastq.bz2"))


# filter based on quality and read length (chose 200bp because that's the median length, according to FASTQC)
filtered_out <- filterAndTrim(fastqFs, filtered_Fs, fastqRs, filtered_Rs, 
                              truncLen = 200, maxEE = 2, maxN = 0, rm.phix = TRUE, 
                              truncQ = 10, compress = TRUE, verbose = TRUE, multithread = TRUE)


# inspect how many reads were filtered out of each sample
filtered_out
(1 - (filtered_out[,2] / filtered_out[,1])) * 100
mean((1 - (filtered_out[,2] / filtered_out[,1])) * 100)
# [1] 31.65811 % reads removed from each sample on average.
# that's a little high


# set sample names to the ID only
names(filtered_Fs) <- Sample_names
names(filtered_Rs) <- Sample_names




############################
### Learning error rates ###
############################

# learn error rates of the F and R reads
error_F <- learnErrors(filtered_Fs, multithread = TRUE)
error_R <- learnErrors(filtered_Rs, multithread = TRUE)


# visualize error frequencies
plotErrors(error_F, nominalQ = TRUE, err_in = TRUE, err_out = TRUE)
plotErrors(error_R, nominalQ = TRUE, err_in = TRUE, err_out = TRUE)




################################
### Merging paired-end reads ###
################################

# create list of merged reads
Mergers <- vector("list", length(Sample_names))
names(Mergers) <- Sample_names


# sample inference and merging paired-end reads
for(sample in Sample_names) 
{
  cat ("Processing:", sample, "\n")
  derep_F <- derepFastq(filtered_Fs[[sample]])
  dada_F <- dada(derep_F, err = error_F, multithread = TRUE)
  derep_R <- derepFastq(filtered_Rs[[sample]])
  dada_R <- dada(derep_R, err = error_R, multithread = TRUE)
  Merger <- mergePairs(dada_F, derep_F, dada_R, derep_R)
  Mergers[[sample]] <- Merger
}


# removing dereps to save memory
rm(derep_F, derep_R)


# contruct a sequence table
Sequence_table <- makeSequenceTable(Mergers)


# dimensions: num rows x num columns
dim(Sequence_table)
# [1]   26 3938





##################################
### Quality control: processed ###
##################################


########
### trim

# inspect sequence distribution
seq_distribution <- data.frame(table(nchar(getSequences(Sequence_table)))) 
colnames(seq_distribution) <- c("Length", "Freq")


# remove reads of non target length, 5% above and below the median 
median(nchar(getSequences(Sequence_table)))
# [1] 253 perfect!

min_len <- floor(median(nchar(getSequences(Sequence_table))) * 0.95)
max_len <- ceiling(median(nchar(getSequences(Sequence_table))) * 1.05)


# modify sequence table with new guidelines
Seq_table_trimmed <- Sequence_table[ ,nchar(colnames(Sequence_table)) 
                                    %in% seq(min_len, max_len)]

dim(Sequence_table)
# [1]   26 3938


dim(Seq_table_trimmed)
# [1]   26 3806



###################
### Remove chimeras

# removing chimeras with denovo screening
Sequence_table_nochim <- removeBimeraDenovo(Seq_table_trimmed, method = "consensus",
                                            multithread = TRUE, verbose = TRUE)

# how many unique sequences were moved?
ncol(Sequence_table) - ncol(Sequence_table_nochim)
# [1] 960


# what percentage of reads were identified as chimeras?
(1 - sum(Sequence_table_nochim) / sum(Sequence_table)) * 100
# [1] 3.950415 awesome




#######################
### Assign taxonomy ###
#######################

Sequence_taxa_table <- assignTaxonomy(Sequence_table_nochim, 
                                      "~/Desktop/Lab/Projects/Misc/Decontam/silva_nr_v132_train_set.fa.gz", 
                                      multithread = TRUE)



# done with DADA2

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# now clean up data




##################################
### remove mock community ASVs ###
##################################

# save FASTA
uniquesToFasta(Sequence_table_nochim, ids = paste0("nochim", 1:ncol(Sequence_table_nochim)),
                                                   fout = "./RData/Seqtabnochim.fasta")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
### used BLAST in command line ###
# makeblastdb -dbtype nucl -in Seqtabnochim.fasta -input_type fasta

# blastn -db Seqtabnochim.fasta -query 16S_mock_all_trimmed.fasta \
#  -task blastn -perc_identity 100 -outfmt 6 -out mock_align.txt
# # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # #


# load results (organized & modified .txt to .csv in excel)
mock_align <- read.csv("./RData/nochim_mock_align.csv")


# add taxonomy
Sequence_taxa_table.df <- data.frame(FASTA = rownames(Sequence_taxa_table),
                                     Sequence_taxa_table)

mock_align <- merge(mock_align, Sequence_table_nochim.fastas, by = "sseqid")
mock_align <- merge(mock_align, Sequence_taxa_table.df, by = "FASTA")


# what were they?
as.character(unique(mock_align$qseqid))
# [1] "Acinetobacter"     "Pseudomonas"       "Bacteroides"       "Rhodobacter"      
# [5] "Escherichia"       "Deinococcus"       "Actinomyces"       "Neisseria"        
# [9] "Propionibacterium" "Streptococcus"     "Bacillus"          "Staphylococcus"   
# [13] "Clostridium"       "Enterococcus"      "Lactobacillus"     "Listeria" 


# What the mock community file says:
# Listeria
# Neisseria
# Propionibacterium
# Pseudomonas
# Rhodobacter
# Staphylococcus
# Streptococcus
# Acinetobacter
# Actinomyces
# Bacillus
# Bacteroides
# Clostridium (same as Clostridium_sensu_stricto_1)
# Deinococcus
# Enterococcus
# Escherichia
# Helicobacter
# Lactobacillus


# save their FASTAS
mock_sseqids <- as.character(unique(mock_align$FASTA))


# remove them
Sequence_table_nomock <- Sequence_table_nochim[, ! colnames(Sequence_table_nochim) %in% mock_sseqids]

ncol(Sequence_table_nochim)
# [1] 2978

ncol(Sequence_table_nomock)
# [1] 2957


# remove the mock sample
rownames(Sequence_table_nomock)[25]
Sequence_table_nomock <- Sequence_table_nomock[-25,]
rownames(Sequence_table_nomock)


# subset the taxonomy table
Taxa_table_nomock <- Sequence_taxa_table[rownames(Sequence_taxa_table) %in% colnames(Sequence_table_nomock),]
identical(rownames(Taxa_table_nomock), colnames(Sequence_table_nomock))
# [1] TRUE


# save it as a data frame
Taxonomy_all <- data.frame(Sequence_taxa_table, FASTA = rownames(Sequence_taxa_table))


# keep the mock ASVs
Mock_ASVs <- subset(subset(Taxonomy_all, FASTA %in% mock_sseqids))
Mock_ASVs$Source <- "Mock"


# subset the good taxa
Taxonomy_all <- data.frame(Taxa_table_nomock, FASTA = rownames(Taxa_table_nomock))




#######################################
### remove no template control ASVs ###
#######################################

# subset ASVs from the mock sample
NTC_ASVs <- Sequence_table_nomock["NTC",]
NTC_ASVs <- NTC_ASVs[NTC_ASVs > 0]
length(NTC_ASVs)
# [1] 276


# calculate the mean abundance of NTC ASVs in the dataset (rounded up),
# turn into a data frame
NTC_ASVs <- data.frame(NTC_abund = NTC_ASVs, 
                   data_mean = ceiling(colMeans(Sequence_table_nomock[-25, colnames(Sequence_table_nomock) %in% names(NTC_ASVs)])))


# ASV mean ≥ 3*NTC = keep
# ASV mean < 3*NTC = contaminant
NTC_ASVs$Contaminant <- NTC_ASVs$data_mean < (3 * NTC_ASVs$NTC_abund)


# what taxa are they?
NTC_ASVs$FASTA <- rownames(NTC_ASVs)
NTC_ASVs <- merge(NTC_ASVs, Taxonomy_all, by = "FASTA")


# save
#write.csv(NTC_ASVs, "./RData/NTC_contaminants.csv", row.names = F, na = "")


# remove them
NTC_ASVs <- subset(NTC_ASVs, Contaminant == TRUE)
Sequence_table_nomockNTC <- Sequence_table_nomock[, ! colnames(Sequence_table_nomock) %in% NTC_ASVs$FASTA]


ncol(Sequence_table_nomock)
# [1] 2957


ncol(Sequence_table_nomockNTC)
# [1] 2693


# remove the NTC sample
rownames(Sequence_table_nomockNTC)[25]
Sequence_table_nomockNTC <- Sequence_table_nomockNTC[-25,]
rownames(Sequence_table_nomockNTC)




###############################
### remove nonspecific ASVs ###
###############################

# create phyloseq object, remove them that way
DWDS_object <- phyloseq(otu_table(Sequence_table_nomockNTC, taxa_are_rows = FALSE),
                             tax_table(Sequence_taxa_table))

euk_ASVs <- taxa_names(subset_taxa(DWDS_object, Kingdom == "Eukaryota"))
chloro_ASVs <- taxa_names(subset_taxa(DWDS_object, Order == "Chloroplast"))
mito_ASVs <- taxa_names(subset_taxa(DWDS_object, Family == "Mitochondria"))


# combine them
NTC_ASVs$Source <- "NTC"

euk_ASVs <- subset(Taxonomy_all, FASTA %in% euk_ASVs)
euk_ASVs$Source <- "Eukaryota"

chloro_ASVs <- subset(Taxonomy_all, FASTA %in% chloro_ASVs)
chloro_ASVs$Source <- "Chloroplasts"

mito_ASVs <- subset(Taxonomy_all, FASTA %in% mito_ASVs)
mito_ASVs$Source <- "Mitochondria"

Contaminants <- rbind(NTC_ASVs[c(5:10,1,11)], Mock_ASVs, euk_ASVs, chloro_ASVs, mito_ASVs)
Data_means <- data.frame(Data_mean = ceiling(colMeans(Sequence_table_nochim[,colnames(Sequence_table_nochim) %in% Contaminants$FASTA])))
Data_means$FASTA <- rownames(Data_means)
Contaminants <- merge(Data_means, Contaminants, by = "FASTA")


# save
#write.csv(Contaminants, "./RData/DWDS_contaminants.csv", row.names = FALSE, na = "")


# remove bad ASVs from object
DWDS_object
# 2693 taxa and 24 samples

DWDS_object <- subset_taxa(DWDS_object, ! taxa_names(DWDS_object) %in% Contaminants$FASTA)

DWDS_object
# 2670 taxa and 24 samples




#########################
### Organize and save ### 
#########################

# update sample info
Sample_info <- read.csv("./RData/DWDS_sample_info.csv")
rownames(Sample_info) <- Sample_info$Sample_name
sample_data(DWDS_object) <- Sample_info


# add FASTA variable to taxa table
temp <- data.frame(DWDS_object@tax_table@.Data)
temp$FASTA <- rownames(temp)
tax_table(DWDS_object) <- as.matrix(temp)


# change taxa names from FASTA sequences to something else
taxa_names(DWDS_object) <- paste0("DWDS", 1:ntaxa(DWDS_object))


# save
#saveRDS(DWDS_object, "./RData/DWDS_phyloseq_object.RData")




############################
### remove ANDIG ASVs ###
############################

# 1. find shared ASVs between ANDIG and DWDS datasets
# 2. find mean relative abundances of ASVs in each dataset
# 3. if mean relative abundance of an ASV is 10x greater in DWDS
#     than ANDIG, then it is DWDS's


# load ANDIG data (see McNamara_dada2.R from 3/26/20)
ANDIG_object <- readRDS("~/Desktop/Lab/Projects/McNamara/RData/McNamara_phyloseq_object.RData")


# which ASVs are shared? (using FASTA sequences, since they have unique taxa names)
Shared_ASVs <- intersect(tax_table(ANDIG_object)[,7], tax_table(DWDS_object)[,7])
length(Shared_ASVs)
# [1] 716


# subset abundance tables to those
Shared_DWDS_counts <- data.frame(subset_taxa(DWDS_object, FASTA %in% Shared_ASVs)@otu_table@.Data)
Shared_ANDIG_counts <- data.frame(subset_taxa(ANDIG_object, FASTA %in% Shared_ASVs)@otu_table@.Data)


# find mean relative abundances of those ASVs, and combine
Shared_DWDS_means <- data.frame(DWDS_ASVs = colnames(Shared_DWDS_counts),
                                   DWDS_means = colMeans(Shared_DWDS_counts / rowSums(Shared_DWDS_counts)),
                                   FASTA = tax_table(subset_taxa(DWDS_object, FASTA %in% Shared_ASVs))[,7])

Shared_ANDIG_means <- data.frame(ANDIG_ASVs = colnames(Shared_ANDIG_counts),
                                   ANDIG_means = colMeans(Shared_ANDIG_counts / rowSums(Shared_ANDIG_counts)),
                                   FASTA = tax_table(subset_taxa(ANDIG_object, FASTA %in% Shared_ASVs))[,7])


# did i do that right?
identical(as.character(Shared_DWDS_means$FASTA), 
                       as.character(tax_table(subset_taxa(DWDS_object, FASTA %in% Shared_ASVs))[,7]))
identical(as.character(Shared_DWDS_means$DWDS_ASVs), 
          as.character(taxa_names(subset_taxa(DWDS_object, FASTA %in% Shared_ASVs))))

identical(as.character(Shared_ANDIG_means$FASTA), 
          as.character(tax_table(subset_taxa(ANDIG_object, FASTA %in% Shared_ASVs))[,7]))
identical(as.character(Shared_ANDIG_means$ANDIG_ASVs), 
          as.character(taxa_names(subset_taxa(ANDIG_object, FASTA %in% Shared_ASVs))))


# merge
Shared_means <- merge(Shared_DWDS_means, Shared_ANDIG_means, by = "FASTA")


# which are Lee's (10x greater in DWDS data than ANDIG data)?
Shared_means$DWDS_final <- Shared_means$DWDS_means > (Shared_means$ANDIG_means * 10)


# how many must be removed?
length(subset(Shared_means, DWDS_final == FALSE)$DWDS_ASVs)
# [1] 646


# add taxonomy
Shared_means <- merge(Shared_means, Taxonomy_all, by = "FASTA")


# save
#write.csv(Shared_means, "./RData/ANDIG_contamination.csv", na = "", row.names = F)


# final dataset
# 2661 taxa and 24 samples
DWDS_object <- subset_taxa(DWDS_object, ! taxa_names(DWDS_object) %in% 
                      subset(Shared_means, DWDS_final == FALSE)$DWDS_ASVs)
DWDS_object
# 2024 taxa and 24 samples


# change DWDS to ASV
taxa_names(DWDS_object) <- paste0("ASV", 1:ntaxa(DWDS_object))


# save
DWDS_counts <- data.frame(t(DWDS_object@otu_table@.Data))
DWDS_tax <- data.frame(DWDS_object@tax_table@.Data)
DWDS_relabun <- DWDS_counts / colSums(DWDS_counts)
DWDS_relabun$Mean_relative_abundance <- rowMeans(DWDS_relabun)
identical(rownames(DWDS_tax), rownames(DWDS_relabun))
DWDS_relabun <- cbind(DWDS_relabun[25], DWDS_tax)



# save
#saveRDS(DWDS_object, "./RData/DWDS_phyloseq_object.RData")
#write.csv(DWDS_relabun, "./RData/Final_DWDS_ASVs.csv", na = "")





# # # # # # # # # # #
# save R environment
save.image("./RData/DWDS_dada2_environment2.RData")





