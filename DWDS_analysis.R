##########################################
### Analysis for Kimbell et al, 2020
### Lou LaMartina, finalized April 28, 2020
##########################################


# set working directory
setwd("~/Desktop/Lab/Projects/Kimbell")


library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(OTUtable)
library(reshape2)
library(ape)
library(ggrepel)
library(indicspecies)
library(ggbeeswarm)


# load phyloseq object (see DWDS_dada2.R from 3/26/20)
DWDS_object <- readRDS("./RData/Kimbell_phyloseq_object.RData")


# convert to relative abundance
DWDS_object_relabun <- transform_sample_counts(DWDS_object, function(x) x/sum(x))


# extract relative abundance info
DWDS_relabun <- data.frame(DWDS_object_relabun@otu_table@.Data)


# remove empty ASVs
DWDS_relabun <- DWDS_relabun[, colSums(DWDS_relabun) > 0]


# load sample info
Sample_info <- read.csv("./RData/Kimbell_sample_info.csv")
rownames(Sample_info) <- Sample_info$Sample_name


# extract taxonomy
Taxonomy_all <- data.frame(DWDS_object@tax_table@.Data)
Taxonomy_all$ASV <- rownames(Taxonomy_all)


# add sample IDs to sample info
# order samples (first by type, then by location, then by distance)
info_order <- c("P30", "P14", "P36", "P27", "P11", "P17", "P3", "P9", 
                "P15", "P6", "P12", "P18", "P23", "P7", "P13", "P26", 
                "P4", "P10", "P19", "P21", "P20", "P44", "P45", "P46")
Sample_info <- Sample_info[match(info_order, Sample_info$Sample_name),]
Sample_info$Sample_num <-  1:nrow(Sample_info)




#######################
### alpha diversity ###
#######################

# measure alpha diversity
shannon <- data.frame(Sample_name = rownames(DWDS_relabun), shannon = diversity(DWDS_relabun, "shannon"))


# add sample info
alpha.div <- merge(shannon, Sample_info, by = "Sample_name")


# are t-tests appropriate to use? (is the data normally distributed?)
for(i in unique(alpha.div$Sample_type)){
  print(i)
  print(shapiro.test(subset(alpha.div[c(7,2)], Sample_type == i)$shannon))
}
# p > 0.05 so yep


# yes, so i can use ANOVA and tests
# are there alpha diversity differences btwn sample types?
summary(aov(shannon ~ Sample_type, alpha.div))
# p = 0.237, no


# pipe location?
summary(aov(shannon ~ Pipe_location, alpha.div))
# p = 0.527


# distance into pipe?
summary(aov(shannon ~ Distance_into_pipe, alpha.div))
# p = 0.955


# between tubercles and scrapes? (looks most different)
t.test(subset(alpha.div, Sample_type == "Scrape")$shannon, 
       subset(alpha.div, Sample_type == "Tubercle")$shannon)
# p = 0.1079


### figure 2 ###
alpha <-
  ggplot(alpha.div, aes(x = Sample_type, y = shannon, color = Sample_type)) +
  geom_point(size = 1.5) +
  geom_text_repel(data = alpha.div, size = 1.8, segment.size = 0.25,
                  aes(label = Sample_num)) +
  theme_classic() +
  scale_color_manual(values = rev(brewer.pal(4, "Set1"))) +
  scale_x_discrete(labels = c("Scrape", "Sub-surface\nswab",
                              "Surface\nswab", "Tubercle")) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_text(size = 6, color = "black", face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25)) +
  guides(color = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1),
         shape = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1)) +
  labs(y = "Shannon alpha diversity score", x = "Sample type")
alpha

#ggsave("./Plots/alpha.pdf", plot = alpha, device = "pdf", width = 2.5, height = 2.5, units = "in")


# save sample information for figure
#write.csv(Sample_info, "./Plots/plotpoints.csv", row.names = FALSE)




##################
### ordination ###
##################

# stat
DWDS.pcoa <- pcoa(vegdist(DWDS_relabun, method = "bray"))


# extract values
DWDS.pcoa.df <- data.frame(DWDS.pcoa$vectors[,1:2])


# add info
DWDS.pcoa.df$Sample_name <- rownames(DWDS.pcoa.df)
DWDS.pcoa.df <- merge(DWDS.pcoa.df, Sample_info, by = "Sample_name")


### figure 3A ###
pcoa <-
  ggplot(DWDS.pcoa.df, aes(x = Axis.1, y = Axis.2, color = Sample_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_text_repel(data = DWDS.pcoa.df, size = 1.8, segment.size = 0.25,
                  aes(label = Sample_num)) +
  geom_point(size = 2) +
  theme_classic() +
  scale_color_manual(values = rev(brewer.pal(4, "Set1")),
                     labels = c("Scrape", "Sub-surface swab",
                                "Surface swab", "Tubercle")) +
  theme(axis.text.x = element_text(size = 5, color = "black"),
        axis.text.y = element_text(size = 5, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_text(size = 6, color = "black", face = "bold"),
        legend.text = element_text(size = 5, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25)) +
  guides(color = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1),
         shape = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1)) +
  labs(y = "Axis 2\n(13.2%)", x = "Axis 1\n(23.5%)", color = "Sample type")
pcoa

#ggsave("./Plots/pcoa.pdf", plot = pcoa, device = "pdf", width = 4.86, height = 2.6, units = "in")




#############################################
### combine taxa to lowest classification ###
#############################################

# get ASVs that were not classified to kingdom level
NAkingdom <- subset(Taxonomy_all, is.na(Taxonomy_all$Kingdom) == TRUE)


# get ASVs that were not classified to phylum level
NAphylum <- subset(Taxonomy_all, is.na(Taxonomy_all$Phylum) == TRUE & 
                     ! Taxonomy_all$ASV %in% NAkingdom$ASV)


# get ASVs that were not classified to class level
NAclass <- subset(Taxonomy_all, is.na(Taxonomy_all$Class) == TRUE & 
                    ! Taxonomy_all$ASV %in% NAkingdom$ASV &
                    ! Taxonomy_all$ASV %in% NAphylum$ASV)


# get ASVs that were not classified to order level
NAorder <- subset(Taxonomy_all, is.na(Taxonomy_all$Order) == TRUE & 
                    ! Taxonomy_all$ASV %in% NAkingdom$ASV &
                    ! Taxonomy_all$ASV %in% NAphylum$ASV &
                    ! Taxonomy_all$ASV %in% NAclass$ASV)


# get ASVs that were not classified to family level
NAfamily <- subset(Taxonomy_all, is.na(Taxonomy_all$Family) == TRUE & 
                     ! Taxonomy_all$ASV %in% NAkingdom$ASV &
                     ! Taxonomy_all$ASV %in% NAphylum$ASV &
                     ! Taxonomy_all$ASV %in% NAclass$ASV &
                     ! Taxonomy_all$ASV %in% NAorder$ASV)


# get ASVs that were not classified to genus level
NAgenus <- subset(Taxonomy_all, is.na(Taxonomy_all$Genus) == TRUE & 
                    ! Taxonomy_all$ASV %in% NAkingdom$ASV &
                    ! Taxonomy_all$ASV %in% NAphylum$ASV &
                    ! Taxonomy_all$ASV %in% NAclass$ASV &
                    ! Taxonomy_all$ASV %in% NAorder$ASV &
                    ! Taxonomy_all$ASV %in% NAfamily$ASV)


# get ASV that were classified to genus level
yesgenus <- subset(Taxonomy_all, is.na(Taxonomy_all$Genus) == FALSE)


# give them a new variable name
NAkingdom_names <- paste0(rep("Unknown", nrow(NAkingdom)))
names(NAkingdom_names) <- NAkingdom$ASV

NAphylum_names <- paste0(rep("Kingdom:", nrow(NAkingdom)), NAphylum$Kingdom)
names(NAphylum_names) <- NAphylum$ASV

NAclass_names <- paste0(rep("Phylum:", nrow(NAclass)), NAclass$Phylum)
names(NAclass_names) <- NAclass$ASV

NAorder_names <- paste0(rep("Class:", nrow(NAorder)), NAorder$Class)
names(NAorder_names) <- NAorder$ASV

NAfamily_names <- paste0(rep("Order:", nrow(NAfamily)), NAfamily$Order)
names(NAfamily_names) <- NAfamily$ASV

NAgenus_names <- paste0(rep("Family:", nrow(NAgenus)), NAgenus$Family)
names(NAgenus_names) <- NAgenus$ASV

genus_names <- paste0(rep("Genus:", nrow(yesgenus)), yesgenus$Genus)
names(genus_names) <- yesgenus$ASV


# combine
lowest_class <- c(NAkingdom_names, NAphylum_names, NAclass_names, NAorder_names, NAfamily_names, NAgenus_names, genus_names)
lowest_class.df <- data.frame(Tax = lowest_class, ASV = names(lowest_class))


# want this in phyloseq object, so can glom
Taxonomy_all2 <- merge(Taxonomy_all, lowest_class.df, by = "ASV")
rownames(Taxonomy_all2) <- Taxonomy_all2$ASV
Taxonomy_all2$Genus <- Taxonomy_all2$Tax
tax_table(DWDS_object_relabun) <- as.matrix(Taxonomy_all2[2:7])


# glom
DWDS_object_relabun
# 2024 taxa

DWDS_tax_glom <- speedyseq::tax_glom(DWDS_object_relabun, "Genus")
DWDS_tax_glom
# 701 taxa




##########################
### stacked bar charts ###
##########################

# extract 11 most abundant
DWDS_top_object <- prune_taxa(names(sort(taxa_sums(DWDS_tax_glom), decreasing = TRUE))[1:11], 
                                 DWDS_tax_glom)


# extract that info
DWDS_top_abundance <- data.frame(DWDS_top_object@otu_table@.Data)
DWDS_top_taxa <- data.frame(DWDS_top_object@tax_table@.Data)


# make column names of abundance table the genus names
colnames(DWDS_top_abundance) <- DWDS_top_taxa$Genus


# want an "other"
DWDS_top_abundance$Other <- 1 - rowSums(DWDS_top_abundance)


# make sample name column
DWDS_top_abundance$Sample_name <- rownames(DWDS_top_abundance)


# melt
Top_genera.m <- melt(DWDS_top_abundance, variable.name = "Genus", value.name = "Relabun")


# rank by axis 1 score
axis_order <- DWDS.pcoa.df[order(DWDS.pcoa.df$Axis.1, decreasing = FALSE), "Sample_name"]


# add axis 1 score
Top_genera.m <- merge(Top_genera.m, DWDS.pcoa.df[1:2], by = "Sample_name")


# add space after colon
toptaxa <- sub(":", ": ", colnames(DWDS_top_abundance[-ncol(DWDS_top_abundance)]))


### figure 3B ###
bars <-
  ggplot(Top_genera.m, aes(x = Sample_name, y = Relabun, fill = Genus, color = Genus)) +
  geom_bar(stat = "identity", width = 0.75, size = 0.1) +
  scale_fill_manual(values = c(brewer.pal(11, "Paired"), "grey90"), labels = toptaxa) +
  scale_color_manual(values = c(brewer.pal(11, "Paired"), "grey90"), guide = FALSE) +
  theme_classic() +
  scale_x_discrete(limits = axis_order, labels = Sample_info[match(axis_order, rownames(Sample_info)), "Sample_num"]) + 
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y = element_text(size = 5, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_text(size = 6, color = "black", face = "bold"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        legend.text = element_text(size = 5, color = "black", face = "italic"),
        panel.border = element_blank(), 
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25)) +
  guides(fill = guide_legend(keyheight = 0.5, keywidth = 0.4, units = "in")) +
  labs(y = "Relative abundances of taxa", x = "Sample", 
       fill = "Lowest\ntaxonomic\nclassification")
bars

#ggsave("./Plots/stackedbars.pdf", plot = bars, device = "pdf", width = 5, height = 2.5, units = "in")


# save data frame
#write.csv(DWDS_top_abundance, "./RData/DWDS_top20_relabun.csv")




##################
### dendrogram ###
##################

# extract info
DWDS_tax_abundance <- data.frame(DWDS_tax_glom@otu_table@.Data)
DWDS_tax <- data.frame(DWDS_tax_glom@tax_table@.Data)
DWDS_tax$ASV <- rownames(DWDS_tax)


# remove ASVs with less than 2% as their maximum relative abundance
maxtax <- apply(DWDS_tax_abundance, 2, max)
mintax <- names(which(maxtax < 0.02))
DWDS_relabun_filt <- DWDS_tax_abundance[, -which(colnames(DWDS_tax_abundance) %in% mintax)]
dim(DWDS_relabun_filt)
# [1] 24 53


identical(rownames(subset(DWDS_tax, rownames(DWDS_tax) %in% colnames(DWDS_relabun_filt))), colnames(DWDS_relabun_filt))
colnames(DWDS_relabun_filt) <- subset(DWDS_tax, rownames(DWDS_tax) %in% colnames(DWDS_relabun_filt))$Genus


# convert to z score  = (x-μ)/σ
DWDS_relabun.z <- t(zscore(t(DWDS_relabun_filt)))


# euclidian distances on z score matrices
ASV_dist <- vegdist(t(DWDS_relabun.z), method = "euclidian")
smp_dist <- vegdist(DWDS_relabun.z, method = "euclidian")


# cluster ASVs and samples based on euclidian distances
ASV_clus <- hclust(ASV_dist, method = "average")
smp_clus <- hclust(smp_dist, method = "average")


# define colors: blue is low, red is high
heat_colors <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(299))


# default heatmap
heatmap(as.matrix(t(DWDS_relabun.z)),
        Rowv = as.dendrogram(ASV_clus),
        Colv = as.dendrogram(smp_clus),
        margins = c(2, 12),
        col = heat_colors)
# save for dendrograms


# add taxa info
ASVs.df <- data.frame(ASV_order = 1:length(colnames(DWDS_relabun.z)), Genus = colnames(DWDS_relabun.z))


# extract order of ASVs
ASV_clus_order.df <- data.frame(Genus = ASV_clus$labels[c(ASV_clus$order)])
ASV_clus_order.df$ASV_order <- 1:nrow(ASV_clus_order.df)


# extract order of samples, add sample info
smp_clus_order.df <- data.frame(Sample_name = smp_clus$labels[c(smp_clus$order)])
smp_clus_order.df$smp_order <- 1:nrow(smp_clus_order.df)
smp_clus_order.df <- merge(smp_clus_order.df, Sample_info, by = "Sample_name")
smp_clus_order.df <- smp_clus_order.df[order(smp_clus_order.df$smp_order),]


# make same heatmap in ggplot2, so i can actually modify it
DWDS_relabun.z.df <- data.frame(DWDS_relabun.z)
DWDS_relabun.z.df$Sample_name <- rownames(DWDS_relabun.z.df)
DWDS_relabun.z.df.m <- melt(DWDS_relabun.z.df, id.vars = "Sample_name", variable.name = "Genus", value.name = "Z")


# need to keep order of ASVs in dendrogram
DWDS_relabun.z.df.m$Genus <- sub("\\.", ": ", DWDS_relabun.z.df.m$Genus)
ASV_clus_order.df$Genus <- sub(":", ": ", ASV_clus_order.df$Genus)
DWDS_relabun.z.df.m <- merge(DWDS_relabun.z.df.m, ASV_clus_order.df, by = "Genus", all.x = TRUE)
DWDS_relabun.z.df.m <- DWDS_relabun.z.df.m[order(DWDS_relabun.z.df.m$ASV_order), ]


# make ASVs an ordered factor
DWDS_relabun.z.df.m$Genus <- factor(DWDS_relabun.z.df.m$Genus, levels = unique(DWDS_relabun.z.df.m$Genus))


# define colors
heat.colors <- c(brewer.pal(11, "RdBu")[9], "white", brewer.pal(11, "RdBu")[5:4])
heat.values <- scales::rescale(c(min(DWDS_relabun.z.df.m$Z), 0, 2, max(DWDS_relabun.z.df.m$Z)))


# ggplot
dendro <- 
  ggplot(DWDS_relabun.z.df.m, aes(x = Sample_name, y = Genus, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colors = heat.colors, values = heat.values) +
  theme_classic() +
  scale_x_discrete(limits = as.character(smp_clus_order.df$Sample_name),
                   labels = as.character(smp_clus_order.df$Sample_num)) +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(size = 5, color = "black"),
        axis.text.y = element_text(size = 4, color = "black", face = "italic"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_text(size = 6, color = "black", face = "bold"),
        legend.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_blank(),
        axis.ticks = element_line(size = 0.25),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) +
  guides(fill = guide_colorbar(barheight = 2, barwidth = 0.3, units = "in")) +
  labs(x = "Sample", y = "Lowest taxonomic classification", fill = "Z score")
dendro

#ggsave("./Plots/heatmap.pdf", plot = dendro, device = "pdf", width = 5, height = 5.5, units = "in")




################################
### determine core community ###
################################

# stat
Sample_info <- Sample_info[order(rownames(Sample_info)),]
identical(rownames(DWDS_relabun), rownames(Sample_info))
ARG.cor <- cbind(DWDS_relabun, Sample_info[10:15])
ARG.cor <- cor(ARG.cor, method = "spearman")


# how do ASV abundances correlate with ARGs?
ARG.cor.m <- melt(ARG.cor, value.name = "rho")


# keep only comparisons of ASVs and ARGs
ARG.cor.m <- subset(ARG.cor.m, Var1 %in% colnames(DWDS_relabun) &
                      Var2 %in% colnames(Sample_info[10:15]))


# which genera are present in how many samples?
DWDS_relabun_binary <- as.matrix((DWDS_relabun > 0) + 0)
Freq <- data.frame(freq = colSums(DWDS_relabun_binary), Var1 = colnames(DWDS_relabun))


# combine
ARG.cor.m <- merge(ARG.cor.m, Freq, by = "Var1")


# get rid of freq = 1, those are not real
ARG.cor.m <- subset(ARG.cor.m, freq > 1)


# any correlation to frequency and rho score?
# using poisson because frequency/counts against continuous rho
summary(glm(freq ~ rho, family = "poisson", data = ARG.cor.m))
# p = < 2e-16 ***



glm <- 
  ggplot(ARG.cor.m, aes(x = freq, y = rho)) +
  geom_smooth(formula = rho ~ freq, color = "black", fill = "grey90", size = 0.3) +
  theme_classic() +
  scale_x_continuous(breaks = NULL) +
  theme(axis.text.x = element_text(size = 5, color = "black"),
        axis.text.y = element_text(size = 5, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25)) +
  labs(x = "Frequency of ASV in dataset\nout of 24 samples", 
       y = "Spearman rho score\ncorrelating ASV to ARG abundance")
glm

#ggsave("./Plots/glm.pdf", plot = glm, device = "pdf", width = 5, height = 2.5, units = "in")



### figure 5 ###
means <- aggregate(rho ~ freq, mean, data = ARG.cor.m)

rho <-
  ggplot(ARG.cor.m, aes(x = as.factor(freq), y = rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_quasirandom(size = 0.3, width = 0.1, alpha = 0.5, color = "grey30") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 5, color = "black"),
        axis.text.y = element_text(size = 5, color = "black"),
        axis.title.y = element_text(size = 6, color = "black", face = "bold"),
        axis.title.x = element_text(size = 6, color = "black", face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25)) +
  geom_crossbar(data = means, aes(ymin = rho, ymax = rho), size = 0.2, width = 0.4, color = "black") +
  labs(x = "Frequency of ASV in dataset\nout of 24 samples", 
       y = "Spearman rho score\ncorrelating ASV to ARG abundance")
rho

#ggsave("./Plots/rho.pdf", plot = rho, device = "pdf", width = 5, height = 3.5, units = "in")




#################################
### correlations among genera ###
#################################

# extract 20 most abundant
DWDS_top_object <- prune_taxa(names(sort(taxa_sums(DWDS_tax_glom), decreasing = TRUE))[1:20], 
                                 DWDS_tax_glom)


# extract that info
DWDS_top_abundance <- data.frame(DWDS_top_object@otu_table@.Data)
DWDS_top_taxa <- data.frame(DWDS_top_object@tax_table@.Data)


# make column names of abundance table the genus names
colnames(DWDS_top_abundance) <- DWDS_top_taxa$Genus


# loop and get rho and p values for each top 20 genus
rho.dm = matrix(NA, 20, 20)
p.dm = matrix(NA, 20, 20)

rownames(rho.dm) <- colnames(DWDS_top_abundance)
colnames(rho.dm) <- colnames(DWDS_top_abundance)
rownames(p.dm) <- colnames(DWDS_top_abundance)
colnames(p.dm) <- colnames(DWDS_top_abundance)

for(var1 in colnames(DWDS_top_abundance)){
  for(var2 in colnames(DWDS_top_abundance)){
    rho.dm[var1,var2] <- cor.test(DWDS_top_abundance[[var1]], DWDS_top_abundance[[var2]],
                                  method = "spearman")$estimate
    p.dm[var1,var2] <- cor.test(DWDS_top_abundance[[var1]], DWDS_top_abundance[[var2]],
                                method = "spearman")$p.value
  }
}


# melt
rho.m <- melt(data.frame(Var1 = rownames(rho.dm), rho.dm), variable.name = "Var2", value.name = "rho")
p.m <- melt(data.frame(Var1 = rownames(p.dm), p.dm), variable.name = "Var2", value.name = "p")


# combine
rho.results <- cbind(rho.m, p = p.m[,3])


# now get rid of duplicates
rho.results <- subset(rho.results, rho < 0.99)
rho.results <- rho.results[ ! duplicated(rho.results$rho),]


# great! now sub the . for : and -
rho.results$Var2 <- sub("\\.", ":", rho.results$Var2)
rho.results[grep("Burkholderia", rho.results$Var2), "Var2"] <- "Genus:Burkholderia-Caballeronia-Paraburkholderia"


# save
#write.csv(rho.results, "./RData/Genus_corr_results.csv", row.names = FALSE)


### figure S7 ###
rho2 <-
  ggplot(rho.results, aes(y = Var2, x = Var1, fill = rho)) +
  geom_tile() +
  scale_fill_gradientn(values = scales::rescale(c(min(rho.results$rho), -0.3, 0, 0.3, max(rho.results$rho))),
                       colors = c(brewer.pal(5, "Spectral")[1:2], "white", brewer.pal(5, "Spectral")[4:5]),
                       breaks = c(-0.76, -0.35, 0.0, 0.45, 0.93)) +
  geom_text(data = subset(rho.results, p < 0.05), 
            aes(label = round(rho, 2)), size = 1.5, color = "white") +
  theme_classic() +
  scale_y_discrete(position = "right", limits = as.character(colnames(DWDS_top_abundance))[-20]) +
  scale_x_discrete(limits = as.character(colnames(DWDS_top_abundance))[-1]) +
  theme(axis.text.x = element_text(size = 5, color = "black", angle = 90, hjust = 1, face = "italic"),
        axis.text.y = element_text(size = 5, color = "black", face = "italic"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 5, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_blank(),
        axis.ticks = element_line(size = 0.25),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) +
  guides(fill = guide_colorbar(barheight = 3, barwidth = 0.5, units = "in")) +
  labs(y = "Genus correlations")
rho2

#ggsave("./Plots/rho2.pdf", plot = rho2, device = "pdf", width = 6, height = 5, units = "in")




######################################
### correlations of genera to ARGs ###
######################################

# stat
identical(rownames(DWDS_top_abundance), rownames(Sample_info))
ARG_cor_data <- cbind(DWDS_top_abundance, Sample_info[10:15])



# loop and get rho and p values for each top 20 genus
ARG_rho.dm = matrix(NA, 26, 26)
ARG_p.dm = matrix(NA, 26, 26)

rownames(ARG_rho.dm) <- colnames(ARG_cor_data)
colnames(ARG_rho.dm) <- colnames(ARG_cor_data)
rownames(ARG_p.dm) <- colnames(ARG_cor_data)
colnames(ARG_p.dm) <- colnames(ARG_cor_data)

for(var1 in colnames(ARG_cor_data)){
  for(var2 in colnames(ARG_cor_data)){
    ARG_rho.dm[var1,var2] <- cor.test(ARG_cor_data[[var1]], ARG_cor_data[[var2]],
                                  method = "spearman")$estimate
    ARG_p.dm[var1,var2] <- cor.test(ARG_cor_data[[var1]], ARG_cor_data[[var2]],
                                method = "spearman")$p.value
  }
}


# melt
ARG_rho.m <- melt(data.frame(Var1 = rownames(ARG_rho.dm), ARG_rho.dm), variable.name = "Var2", value.name = "rho")
ARG_p.m <- melt(data.frame(Var1 = rownames(ARG_p.dm), ARG_p.dm), variable.name = "Var2", value.name = "p")


# combine
ARG_rho.results <- cbind(ARG_rho.m, p = ARG_p.m[,3])


# only want comparisons of genus ~ ARG
ARG_rho.results <- ARG_rho.results[grepl("copies", ARG_rho.results$Var1),]
ARG_rho.results <- ARG_rho.results[! grepl("copies", ARG_rho.results$Var2),]
ARG_rho.results$Var1 <- droplevels(ARG_rho.results$Var1)
ARG_rho.results$Var2 <- droplevels(ARG_rho.results$Var2)


# now get rid of duplicates
ARG_rho.results <- subset(ARG_rho.results, rho < 1)
ARG_rho.results <- ARG_rho.results[ ! duplicated(ARG_rho.results$rho),]


# great! now sub the . for : and -
ARG_rho.results$Var2 <- sub("\\.", ":", ARG_rho.results$Var2)
ARG_rho.results[grep("Burkholderia", ARG_rho.results$Var2), "Var2"] <- "Genus:Burkholderia-Caballeronia-Paraburkholderia"


### figure S8 ###
rho3 <-
  ggplot(ARG_rho.results, aes(y = Var1, x = Var2, fill = rho)) +
  geom_tile() +
  scale_fill_gradientn(values = scales::rescale(c(min(ARG_rho.results$rho), -0.3, 0, 0.3, max(ARG_rho.results$rho))),
                       colors = c(brewer.pal(5, "Spectral")[1:2], "white", brewer.pal(5, "Spectral")[4:5]),
                       breaks = c(-0.54, -0.25, 0.0, 0.25, 0.51)) +
  geom_text(data = subset(ARG_rho.results, p < 0.05), 
            aes(label = round(rho, 2)), size = 1.5, color = "white") +
  theme_classic() +
  scale_y_discrete(position = "right", labels = rev(c("sul1", "intI1", "czcD", "copA", "blaTEM", "blaSHV"))) +
  theme(axis.text.x = element_text(size = 5, color = "black", angle = 90, hjust = 1, face = "italic"),
        axis.text.y = element_text(size = 5, color = "black", face = "italic"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 5, color = "black"),
        legend.title = element_text(size = 6, color = "black", face = "bold"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_blank(),
        axis.ticks = element_line(size = 0.25),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) +
  guides(fill = guide_colorbar(barheight = 3, barwidth = 0.5, units = "in")) +
  labs(y = "Genus correlations")
rho3

#ggsave("./Plots/rho3.pdf", plot = rho3, device = "pdf", width = 5, height = 3, units = "in")




####################
### indicspecies ###
####################

# add taxa names
indic_data <- data.frame(ASV = colnames(DWDS_relabun), t(DWDS_relabun))
indic_data <- merge(lowest_class.df, indic_data, by = "ASV")
rownames(indic_data) <- paste0(indic_data$ASV, "__", indic_data$Tax)
indic_data <- t(indic_data[-c(1:2)])


# indicpsecies analysis
Sample_info <- Sample_info[order(rownames(Sample_info)),]
identical(rownames(indic_data), rownames(Sample_info))
DWDS_indic <- multipatt(indic_data, as.vector(Sample_info$Sample_type), control = how(nperm = 999))
summary(DWDS_indic, indvalcomp = TRUE)
# copy and pasted into excel




#####################################
### explain community composition ###
#####################################

# subset only relevant sample info, including if something is a tubercle or not
variables <- Sample_info[, c("Sample_type", "Pipe_location", "Distance_into_pipe", "X16S_copies")]
identical(rownames(indic_data), rownames(variables))
variables$Tubercle <- FALSE
variables$Tubercle[variables$Sample_type == "Tubercle"] <- TRUE


# environmental fit
DWDS.cca <- cca(indic_data ~ ., variables, na.action = na.exclude)
DWDS.fit <- envfit(DWDS.cca, variables, permu = 999, na.rm = TRUE)
DWDS.fit

#                CCA1    CCA2     r2 Pr(>r)  
# X16S_copies 0.88300 0.46937 0.2468  0.042 *

#                        r2 Pr(>r)    
# Sample_type        0.3243  0.009 ** 
# Pipe_location      0.1058  0.096 .  
# Distance_into_pipe 0.3596  0.001 ***
# Tubercle           0.2811  0.002 ** 




#################
### questions ###
#################

# how many ASVs belonged mycobacterium?
nrow(subset(Taxonomy_all, Genus == "Mycobacterium")) / nrow(Taxonomy_all) * 100
# [1] 22 / 2024 = 1.086957 %


# how many reads belong to mycobacterium?
sum(subset_taxa(DWDS_object, Genus == "Mycobacterium")@otu_table@.Data) /
  sum(DWDS_object@otu_table@.Data) * 100
# [1] 32.00618 %


# how many genera were there across the dataset?
DWDS_tax_glom # combined to lowest possible classification
# 701 taxa
length(unique(Taxonomy_all$Genus)) # genus only
# [1] 469


# How many were assigned to the “core” community?
# (using frequency) -
data.frame(table(Freq$freq))


# compare envfit results PERMANOVA
identical(rownames(indic_data), rownames(Sample_info))
type.adonis <- adonis(indic_data ~ Sample_info$Sample_type, method = "bray", perm = 999)
type.adonis # 0.029 *
DWDS.fit # 0.001 ***

loc.adonis <- adonis(indic_data ~ Sample_info$Pipe_location, method = "bray", perm = 999)
loc.adonis # 0.712
DWDS.fit # 0.010 **

dist.adonis <- adonis(indic_data ~ Sample_info$Distance_into_pipe, method = "bray", perm = 999)
dist.adonis # 0.471
DWDS.fit # 0.017 *

tub.adonis <- adonis(indic_data ~ variables$Tubercle, method = "bray", perm = 999)
tub.adonis # 0.004 **
DWDS.fit # 0.001 ***


# compare envfit results PERMANOVA
type.ano <- anosim(indic_data, Sample_info$Sample_type, distance = "bray", perm = 9999)
type.ano # 0.0275

loc.ano <- anosim(indic_data, Sample_info$Pipe_location, distance = "bray", perm = 9999)
loc.ano # 0.4746 

dist.ano <- anosim(indic_data, Sample_info$Distance_into_pipe, distance = "bray", perm = 9999)
dist.ano # 0.4356

tub.ano <- anosim(indic_data, variables$Tubercle, distance = "bray", perm = 9999)
tub.ano # 0.0174

