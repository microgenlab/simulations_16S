library(phyloseq)
library(ggplot2)
library(stringr)
library(dplyr)
library(Metrics)
library(viridis)
library(dplyr)
library(Metrics)
library(RColorBrewer)
library(tidyr)
library(ggsignif)

setwd("./")

df <- read.csv("raw/BVBRC_genome.csv", header = T)
head(df)

dff <- select(df, Species, Genome.ID, GenBank.Accessions)

h <- read.table("raw/headers_human.txt", 
                sep = "\t", header = F)
head(h)

h$V1 <- gsub("sim_16S_barrnap_accn_", "", h$V1)
h$V1 <- gsub(".fasta.fastq", "", h$V1)

colnames(h) <- c("GenBank.Accessions")

meta <-merge(dff, h, by = "GenBank.Accessions")
meta <- meta[-which(meta$Species == "null"),]
dim(meta)

retrive <- as.data.frame(meta$GenBank.Accessions)

dir.create("output_tables")
write.table(retrive, "output_tables/retrive.txt",
            sep = "\t", row.names = F, quote = F)



####################################
##########porefile##################
####################################

#Matching step


setwd("raw/results_silva_polishON_human/")
list.files()

c1 <- as.data.frame(read.table("COUNTS.tsv", sep = "\t", header = T))
t1 <- as.data.frame(read.table("TAXCLA.tsv", sep = "\t", header = T))
name <- as.data.frame(colnames(c1))
df <- as.data.frame(str_split_fixed(name$`colnames(c1)`, "_", 3))
df1 <- as.data.frame(cbind(name, df$V2))
colnames(df1) <- c("file", "id")
head(df1)

meta1 <- df1
meta1$round <- c("match")
row.names(meta1) <- meta1$file
head(meta1)
dim(meta1)

c1 <- as.matrix(read.table("COUNTS.tsv", sep = "\t", header = T))
t1 <- as.matrix(read.table("TAXCLA.tsv", sep = "\t", header = T))

OTU1 = otu_table(c1, taxa_are_rows = TRUE)
TAX1 = tax_table(t1)
samples = sample_data(meta1)
OTU1
TAX1

ps1 = phyloseq(OTU1, TAX1, samples)
ps1

top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:10000]
ps.top20 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

graph1 <- psmelt(ps.top20)
head(graph1)
tab2 <- graph1
tab2$step <- c("match")

coul <- brewer.pal(8, "Spectral")
coul <- colorRampPalette(coul)(100)

target <- meta$Species
pos <- filter(graph1, Species %in% target)
min(pos$Abundance)
neg <- graph1[!graph1$Species %in% pos$Species, ]
mean(neg$Abundance)


########################################################################
#Polishing step

count <- read.table("NanoPlots/summary.tsv", sep = "\t", header = T)
c1 <- as.data.frame(read.table("COUNTS_polished.tsv", sep = "\t", header = T))
t1 <- as.data.frame(read.table("TAXCLA_polished.tsv", sep = "\t", header = T))
name <- as.data.frame(colnames(c1))
df <- as.data.frame(str_split_fixed(name$`colnames(c1)`, "_", 3))
df1 <- as.data.frame(cbind(name, df$V2))
colnames(df1) <- c("file", "id")
head(df1)

meta1 <- df1
meta1$round <- c("polish")
row.names(meta1) <- meta1$file
head(meta1)
dim(meta1)

c1 <- as.matrix(read.table("COUNTS_polished.tsv", sep = "\t", header = T))
t1 <- as.matrix(read.table("TAXCLA_polished.tsv", sep = "\t", header = T))

OTU1 = otu_table(c1, taxa_are_rows = TRUE)
TAX1 = tax_table(t1)
samples = sample_data(meta1)
OTU1
TAX1

ps1 = phyloseq(OTU1, TAX1, samples)
ps1

top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:10000]
ps.top20 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

graph2 <- psmelt(ps.top20)

target <- meta$Species
pos <- filter(graph2, Species %in% target)
dim(pos)
min(pos$Abundance)
neg <- graph2[!graph2$Species %in% pos$Species, ]
mean(neg$Abundance)

pos <- pos %>% group_by(id, Species) %>%
  summarise(sum = sum(Abundance))

oth <- neg %>% group_by(id) %>%
  summarise(sum = sum(Abundance))
oth$Species <- c("Others")

oth <- select(oth, id, Species, sum) #porefile have 6.8% of misclassifications after the polishing step at Species-level
oth$sum*100

all <- as.data.frame(rbind(pos, oth))

p <- ggplot(all, aes(x = id, y = sum, fill = Species))
g2 <- p +  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(values = coul) +
  theme(legend.position="right", axis.line = element_line(colour = "black"), 
        text = element_text(size = 10), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 10, angle = 90),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(fill=guide_legend(ncol = 5)) + ggtitle("porefile") +
  xlab("") + ylab("Relative Abundance") + labs(fill = "Species") 
g2

metadata <- meta
dim(metadata)
metadata$expected <- as.numeric(("0.99"))

joined <- left_join(metadata, all, by = "Species")
joined$observed <- joined$sum*100
head(joined)
joined$original <- c("retrived")
joined$id <- joined$id %>% replace_na('not_retrived')
joined$status <- joined$original == joined$id
joined$observed[is.na(joined$observed)] <- 0
head(joined)
table(joined$status) #90 out of 101 HCC components are detected with porefile

pf <-joined[which(joined$status == "FALSE"),]


porefile <- select(joined, Species, expected, observed)
porefile$observed <- round(porefile$observed, 1)
porefile$method <- c("porefile")
porefile[is.na(porefile)] <- 0
porefile$rmse <- rmse(porefile$expected, porefile$observed)
mean(porefile$observed)

#################################
##############EMU################
#################################

emu <- read.table("../emu_all_retrived_human/all_retrived_human_sim_default_rel-abundance.tsv", sep = "\t", header = T)
emu$abundance <- as.numeric(emu$abundance)
emu$abundance <- round(emu$abundance, 4)
emu$abundance[is.na(emu$abundance)] <- 0

emu$id <- c("human")

target <- meta$Species
dim(emu)
pos <- filter(emu, species %in% target)
dim(pos)
round(min(pos$abundance), 5)

neg <- emu[!emu$species %in% pos$species, ]
dim(neg)
mean(neg$abundance)
neg$abundance <- as.numeric(neg$abundance)
neg$abundance <- round(neg$abundance, 4)

pos <- pos %>% group_by(id, species) %>%
  summarise(sum = sum(abundance))

oth <- neg %>% group_by(id) %>%
  summarise(sum = sum(as.numeric(abundance)))
oth$species <- c("Others") #11.2% of misclasifications with EMU
oth$sum*100

oth <- select(oth, id, species, sum)

all <- as.data.frame(rbind(pos, oth))

p <- ggplot(all, aes(x = id, y = sum, fill = species))
g2 <- p +  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(values = coul) +
  theme(legend.position="right", axis.line = element_line(colour = "black"), 
        text = element_text(size = 10), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 10, angle = 90),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(fill=guide_legend(ncol = 5)) + ggtitle("porefile") +
  xlab("") + ylab("Relative Abundance") + labs(fill = "Species") 
g2

metadata <- meta
metadata$expected <- as.numeric(("0.99"))
head(all)
colnames(all)[2] <- c("Species")

joined <- left_join(metadata, all, by = "Species")
joined$observed <- joined$sum*100
head(joined)
joined$id <- joined$id %>% replace_na('not_retrived')
joined$original <- c("human")
joined$status <- joined$original == joined$id
joined$observed[is.na(joined$observed)] <- 0
head(joined)
table(joined$status) #92 out of 101 components were detected with EMU

e <- joined[which(joined$status == "FALSE"),]

intersect(pf$Species, e$Species)
setdiff(pf$Species, e$Species)

emu <- select(joined, Species, expected, observed)
emu$observed <- round(emu$observed, 1)
emu$method <- c("EMU")
emu[is.na(emu)] <- 0
emu$rmse <- rmse(emu$expected, emu$observed)
mean(emu$observed)


##########################################
#########wf-metagenomics##################
##########################################

setwd("../output_wf_human/")
df <- read.table("wf-metagenomics-counts.tsv", sep = "\t", header = T)
df$barcode01 <- as.numeric(df$barcode01)
colnames(df) <- c("Species", "count")
head(df)
s <- sum(df$count)
df$abundance <- df$count/s
df$observed <- df$abundance*100
df$observed <- round(df$observed, 2)
head(df)
df$id <- c("human")

target <- meta$Species
pos <- filter(df, Species %in% target)
dim(pos)
max(pos$abundance)
neg <- df[!df$Species %in% pos$Species, ]
mean(neg$abundance)
neg$abundance <- as.numeric(neg$abundance)

pos <- pos %>% group_by(id, Species) %>%
  summarise(sum = sum(abundance))

oth <- neg %>% group_by(id) %>%
  summarise(sum = sum(abundance))
oth$Species <- c("Others")
oth <- select(oth, id, Species, sum)
oth$sum*100 #25.7% of misclassifications with wf-metagenomics

all <- as.data.frame(rbind(pos, oth))
head(all)

metadata <- meta
metadata$expected <- as.numeric(("0.99"))
head(all)
colnames(all)[2] <- c("Species")

joined <- left_join(metadata, all, by = "Species")
joined$observed <- joined$sum*100
head(joined)
dim(joined)
joined$id <- joined$id %>% replace_na('not_retrived')
joined$original <- c("human")
joined$status <- joined$original == joined$id
joined$observed[is.na(joined$observed)] <- 0
head(joined)
table(joined$status) #94 out of 101 components were detected with wf-metagenomics

wm <- joined[which(joined$status == "FALSE"),]

wf <- select(joined, Species, expected, observed)
wf$observed <- round(wf$observed, 2)
wf$method <- c("wf-metagenomics")
wf[is.na(wf)] <- 0
wf$rmse <- rmse(wf$expected, wf$observed)
mean(wf$observed)


########################################
#############combine output#############
########################################

tab <- as.data.frame(rbind(porefile, emu, wf))
tab[is.na(tab)] <- 0
head(tab)
tab$method_f = factor(tab$method, levels=c("Expected", "porefile", "EMU", "wf-metagenomics"))

tab$host <- c("human")
tab1 <- tab
human <- tab
head(tab1)

p <- ggplot(tab, aes(x=method_f, y=observed)) +
  geom_boxplot(alpha = 0) +
  theme(legend.position="none", 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 10),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  ylab("Abundance (%)") + 
  geom_signif(comparisons = list(c("porefile", "EMU")), map_signif_level = T, y_position = 2) +
  geom_signif(comparisons = list(c("porefile", "wf-metagenomics")), map_signif_level = T, y_position = 2.5) +
  geom_signif(comparisons = list(c("EMU", "wf-metagenomics")), map_signif_level = T, y_position = 2.2) +
  scale_fill_manual(values = c("#d62828", "#f77f00","#fcbf49")) + 
  geom_hline(yintercept = 0.99, linetype = "dotted") +
  geom_jitter(size = 5, pch=21, aes(fill = method_f, alpha = 0.3)) 
p

library(ggpubr)

ls <- ggscatter(tab, x = "expected", y = "observed", color = "Species", size = 5, alpha = 0.5,
                add = "reg.line",  add.params = list(color = "black", fill = "lightgray"),
                conf.int = TRUE) + stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),
                                            label.x = 1) +
  geom_point(pch=21, size = 5) +
  theme(text = element_text(size = 12),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(color=guide_legend(nrow=25,byrow=TRUE)) +
  scale_color_viridis(option = "F", direction = -1, discrete = T) +
  ylab("Observed") + xlab("Expected") + 
  scale_x_continuous(breaks=seq(0,1.5,by=0.5)) +
  geom_hline(yintercept = 0.99, linetype = "dotted") +
  facet_grid(~method_f, scales = "free_x") +
  ggtitle("")
ls

head(tab)
tab$dis <- tab$expected - tab$observed
tab[is.na(tab)] <- 0

ptab <- select(tab, Species, expected)
ptab$method <- c("Expected")
colnames(ptab)[2] <- c("abundance")
head(ptab) 

ctab <- select(tab, Species, observed, method)
colnames(ctab)[2] <- c("abundance")
head(ctab)
c <- as.data.frame(rbind(ptab, ctab))
head(c)
c <- c %>% filter(abundance > 0)

c$method_f <- factor(c$method, levels=c("Expected", "porefile", "EMU", "wf-metagenomics"))
c$method <- NULL

wide <- reshape(c, idvar = "Species", timevar = "method_f", direction = "wide")
wide[is.na(wide)] <- 0
set.seed(335)
w <- wide[sample(nrow(wide), 50), ]

library(data.table)
long <- melt(setDT(w), id.vars = c("Species"), variable.name = "method")
colnames(long) <- c("Species", "method", "value")
long$method <- gsub("abundance.", "", long$method)
long$method_f <- factor(long$method, levels=c("Expected", "porefile", "EMU", "wf-metagenomics"))
long <- long %>% filter(value > 0)


coul <- brewer.pal(8, "Spectral")
coul <- colorRampPalette(coul)(10)

q <- ggplot(long, aes(x=method_f, y=Species, fill=value))+
  geom_tile(colour="white", size=0.5) +
  scale_fill_gradient2(low = "#f6f4d2", mid = "#f6bd60", high = "#e63946") +
  theme(legend.position="right", axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(face = "italic"),
        text = element_text(size = 12), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 10, angle = 70, hjust = 1, vjust = 1),
        strip.background = element_rect(color="black", fill="white", size=.1, linetype="solid")) +
  labs(fill='Abundance (%)') 
q


#################################################################################################

#######################
#####Environment#######
#######################


library(dplyr)
library(Metrics)

setwd("../")

df <- read.csv("BVBRC_genome_env.csv", header = T)
head(df)

dff <- select(df, Species, Genome.ID, GenBank.Accessions)


#######################################
###########porefile####################
#######################################

#Matching agains SILVA
setwd("results_silva_polishON_env/")

c1 <- as.data.frame(read.table("COUNTS.tsv", sep = "\t", header = T))
t1 <- as.data.frame(read.table("TAXCLA.tsv", sep = "\t", header = T))
name <- as.data.frame(colnames(c1))
df <- as.data.frame(str_split_fixed(name$`colnames(c1)`, "_", 3))
df1 <- as.data.frame(cbind(name, df$V2))
colnames(df1) <- c("file", "id")
head(df1)

meta1 <- df1
meta1$round <- c("match")
row.names(meta1) <- meta1$file
head(meta1)
dim(meta1)

c1 <- as.matrix(read.table("COUNTS.tsv", sep = "\t", header = T))
t1 <- as.matrix(read.table("TAXCLA.tsv", sep = "\t", header = T))

OTU1 = otu_table(c1, taxa_are_rows = TRUE)
TAX1 = tax_table(t1)
samples = sample_data(meta1)
OTU1
TAX1

ps1 = phyloseq(OTU1, TAX1, samples)
ps1

top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:110]
ps.top20 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)


graph1 <- psmelt(ps.top20)
head(graph1)
tab2 <- graph1
tab2$step <- c("match")

coul <- brewer.pal(8, "Spectral")
coul <- colorRampPalette(coul)(100)

target <- dff$Species
pos <- filter(graph1, Species %in% target)
min(pos$Abundance)
neg <- graph1[!graph1$Species %in% pos$Species, ]
mean(neg$Abundance)


p <- ggplot(graph1, aes(x = id, y = Abundance, fill = Species))
g1 <- p +  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(values = coul) +
  theme(legend.position="right", axis.line = element_line(colour = "black"), 
        text = element_text(size = 10), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 10, angle = 90),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(fill=guide_legend(ncol = 5)) + ggtitle("porefile") +
  xlab("") + ylab("Relative Abundance") + labs(fill = "Species") 
g1


########################################################################
count <- read.table("NanoPlots/summary.tsv", sep = "\t", header = T)
c1 <- as.data.frame(read.table("COUNTS_polished.tsv", sep = "\t", header = T))
t1 <- as.data.frame(read.table("TAXCLA_polished.tsv", sep = "\t", header = T))
name <- as.data.frame(colnames(c1))
df <- as.data.frame(str_split_fixed(name$`colnames(c1)`, "_", 3))
df1 <- as.data.frame(cbind(name, df$V2))
colnames(df1) <- c("file", "id")
head(df1)

meta1 <- df1
meta1$round <- c("polish")
row.names(meta1) <- meta1$file
head(meta1)
dim(meta1)

c1 <- as.matrix(read.table("COUNTS_polished.tsv", sep = "\t", header = T))
t1 <- as.matrix(read.table("TAXCLA_polished.tsv", sep = "\t", header = T))

OTU1 = otu_table(c1, taxa_are_rows = TRUE)
TAX1 = tax_table(t1)
samples = sample_data(meta1)
OTU1
TAX1

ps1 = phyloseq(OTU1, TAX1, samples)
ps1

top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:10000]
ps.top20 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
graph1 <- psmelt(ps.top20)


target <- dff$Species
pos <- filter(graph1, Species %in% target)
dim(pos)
min(pos$Abundance)
neg <- graph1[!graph1$Species %in% pos$Species, ]
mean(neg$Abundance)

pos <- pos %>% group_by(id, Species) %>%
  summarise(sum = sum(Abundance))

oth <- neg %>% group_by(id) %>%
  summarise(sum = sum(Abundance))
oth$Species <- c("Others")

oth <- select(oth, id, Species, sum)

all <- as.data.frame(rbind(pos, oth))

p <- ggplot(all, aes(x = id, y = sum, fill = Species))
g2 <- p +  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(values = coul) +
  theme(legend.position="right", axis.line = element_line(colour = "black"), 
        text = element_text(size = 10), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 10, angle = 90),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(fill=guide_legend(ncol = 5)) + ggtitle("porefile") +
  xlab("") + ylab("Relative Abundance") + labs(fill = "Species") 
g2

metadata <- dff

metadata$expected <- as.numeric(("8.3"))

joined <- left_join(metadata, all, by = "Species")
joined$observed <- joined$sum*100
head(joined)

porefile <- select(joined, Species, expected, observed)
porefile$observed <- round(porefile$observed, 1)
porefile$method <- c("porefile")
porefile[is.na(porefile)] <- 0
porefile$rmse <- rmse(porefile$expected, porefile$observed)
mean(porefile$observed)
joined$original <- c("retrived")
joined$id <- joined$id %>% replace_na('not_retrived')
joined$status <- joined$original == joined$id
joined$observed[is.na(joined$observed)] <- 0
head(joined)
table(joined$status)

pf <-joined[which(joined$status == "FALSE"),]

################################
#############EMU################
################################

emu <- read.table("../emu_all_retrived_env/all_retrived_env_sim_default_rel-abundance.tsv", sep = "\t", header = T)
dim(emu)
emu$abundance <- as.numeric(emu$abundance)
emu$abundance <- round(emu$abundance, 4)
emu$abundance[is.na(emu$abundance)] <- 0
emu$id <- c("environment")

target <- dff$Species
pos <- filter(emu, species %in% target)
dim(pos)
sum(pos$abundance)
head(emu)
neg <- emu[!emu$species %in% pos$species, ]
mean(neg$abundance)
neg$abundance <- as.numeric(neg$abundance)

pos <- pos %>% group_by(id, species) %>%
  summarise(sum = sum(abundance))

oth <- neg %>% group_by(id) %>%
  summarise(sum = sum(abundance))
oth$species <- c("Others")

oth <- select(oth, id, species, sum)

all <- as.data.frame(rbind(pos, oth))

p <- ggplot(all, aes(x = id, y = sum, fill = species))
g2 <- p +  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(values = coul) +
  theme(legend.position="right", axis.line = element_line(colour = "black"), 
        text = element_text(size = 10), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 10, angle = 90),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(fill=guide_legend(ncol = 5)) + ggtitle("porefile") +
  xlab("") + ylab("Relative Abundance") + labs(fill = "Species") 
g2

metadata <- dff

metadata$expected <- as.numeric(("8.3"))
head(all)
colnames(all)[2] <- c("Species")

joined <- left_join(metadata, all, by = "Species")
joined$observed <- joined$sum*100
head(joined)
joined$id <- joined$id %>% replace_na('not_retrived')
joined$original <- c("environment")
joined$status <- joined$original == joined$id
joined$observed[is.na(joined$observed)] <- 0
head(joined)
table(joined$status)

e <- joined[which(joined$status == "FALSE"),]

emu <- select(joined, Species, expected, observed)
emu$observed <- round(emu$observed, 1)
emu$method <- c("EMU")
emu[is.na(emu)] <- 0
emu$rmse <- rmse(emu$expected, emu$observed)
mean(emu$observed)

###########################
#######wf-metagenomics#####
###########################

setwd("../output_wf_env/")
df <- read.table("wf-metagenomics-counts.tsv", sep = "\t", header = T)
df$barcode01 <- as.numeric(df$barcode01)
colnames(df) <- c("Species", "count")
head(df)
s <- sum(df$count)
df$abundance <- df$count/s
df$observed <- df$abundance*100
df$observed <- round(df$observed, 2)
head(df)
df$id <- c("environment")

target <- dff$Species
pos <- filter(df, Species %in% target)
dim(pos)
max(pos$abundance)
neg <- df[!df$Species %in% pos$Species, ]
mean(neg$abundance)
neg$abundance <- as.numeric(neg$abundance)

pos <- pos %>% group_by(id, Species) %>%
  summarise(sum = sum(abundance))

oth <- neg %>% group_by(id) %>%
  summarise(sum = sum(abundance))
oth$Species <- c("Others")
oth <- select(oth, id, Species, sum)

all <- as.data.frame(rbind(pos, oth))

joined <- left_join(metadata, all, by = "Species")
joined$observed <- joined$sum*100
head(joined)
dim(joined)
joined$id <- joined$id %>% replace_na('not_retrived')
joined$original <- c("environment")
joined$status <- joined$original == joined$id
joined$observed[is.na(joined$observed)] <- 0
head(joined)
table(joined$status)

wm <- joined[which(joined$status == "FALSE"),]

wf <- select(joined, Species, expected, observed)
wf$observed <- round(wf$observed, 2)
wf$method <- c("wf-metagenomics")
wf[is.na(wf)] <- 0
wf$rmse <- rmse(wf$expected, wf$observed)
mean(wf$observed)


tab <- as.data.frame(rbind(porefile, emu, wf))
tab[is.na(tab)] <- 0
head(tab)

tab$method_f = factor(tab$method, levels=c("Expected", "porefile", "EMU", "wf-metagenomics"))
tab$host <- c("environment")
tab2 <- tab
envir <- tab 

p <- ggplot(tab, aes(x=method_f, y=observed)) +
  geom_boxplot(alpha = 0) +
  theme(legend.position="none", 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 12),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  ylab("Abundance (%)") + 
  geom_signif(comparisons = list(c("porefile", "EMU")), map_signif_level = T, y_position = 25) +
  geom_signif(comparisons = list(c("porefile", "wf-metagenomics")), map_signif_level = T, y_position = 26) +
  geom_signif(comparisons = list(c("EMU", "wf-metagenomics")), map_signif_level = T, y_position = 27) +
  scale_fill_manual(values = c("#d62828", "#f77f00","#fcbf49")) + 
  geom_hline(yintercept = 8.3, linetype = "dotted") +
  geom_jitter(size = 5, pch=21, aes(fill = method_f, alpha = 0.3)) 
p


ls <- ggscatter(tab, x = "expected", y = "observed", color = "Species", size = 5, alpha = 0.5,
                add = "reg.line",  add.params = list(color = "black", fill = "lightgray"),
                conf.int = TRUE) + stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),
                                            label.x = 1) +
  geom_point(pch=21, size = 5)+
  theme(text = element_text(size = 12),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(color=guide_legend(nrow=25,byrow=TRUE)) +
  scale_color_viridis(option = "F", direction = -1, discrete = T) +
  ylab("Observed") + xlab("Expected") + 
  scale_x_continuous(breaks=seq(0,1.5,by=0.5)) +
  geom_hline(yintercept = 8.3, linetype = "dotted") +
  facet_grid(~method_f, scales = "free_x") +
  ggtitle("")
ls

head(tab)
tab$dis <- tab$expected - tab$observed
tab[is.na(tab)] <- 0

ptab <- select(tab, Species, expected)
ptab$method <- c("Expected")
colnames(ptab)[2] <- c("abundance")
head(ptab) 

ctab <- select(tab, Species, observed, method)
colnames(ctab)[2] <- c("abundance")
head(ctab)
c <- as.data.frame(rbind(ptab, ctab))
head(c)
c <- c %>% filter(abundance > 0)

c$method_f <- factor(c$method, levels=c("Expected", "porefile", "EMU", "wf-metagenomics"))
c$method <- NULL
head(c)

coul <- brewer.pal(8, "Spectral")
coul <- colorRampPalette(coul)(10)

q <- ggplot(c, aes(x=method_f, y=Species, fill=abundance))+
  geom_tile(colour="white", size=0.5) +
  scale_fill_gradient2(low = "#f6f4d2", mid = "#f6bd60", high = "#e63946") +
  theme(legend.position="right", axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(face = "italic"),
        text = element_text(size = 12), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 10, angle = 70, hjust = 1, vjust = 1),
        strip.background = element_rect(color="black", fill="white", size=.1, linetype="solid")) +
  labs(fill='Abundance (%)') 
q


##########################################################
#########combine human and environment results############
##########################################################

tab <- as.data.frame(rbind(human, envir))
p <- ggplot(tab, aes(x=method_f, y=observed)) +
  geom_boxplot(alpha = 0) +
  theme(legend.position="bottom",
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_text(size = 28),
        axis.title.y = element_text(size = 34),
        strip.text = element_text(size = 28),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  ylab("Abundance (%)") + 
  labs(fill="Source") +
  geom_signif(comparisons = list(c("porefile", "EMU")), map_signif_level = T, y_position = 25, size = 1, textsize = 10) +
  geom_signif(comparisons = list(c("porefile", "wf-metagenomics")), map_signif_level = T, y_position = 29, size = 1, textsize = 10) +
  geom_signif(comparisons = list(c("EMU", "wf-metagenomics")), map_signif_level = T, y_position = 27, size = 1, textsize = 10) +
  scale_fill_manual(values = c("#5c8001", "#f77f00")) + 
  geom_hline(yintercept = 8.3, linetype = "dotted", size = 3) +
  geom_hline(yintercept = 0.99, linetype = "dotted", size = 3) +
  geom_jitter(size = 10, pch=21, aes(fill = host)) 
p



ls <- ggscatter(tab, x = "expected", y = "observed", color = "host", size = 10, alpha = 1,
                palette = c("#5c8001", "#f77f00"),
                add = "reg.line",  add.params = list(color = "black", fill = "lightgray"),
                conf.int = TRUE) + stat_cor(size = 14, aes(label = paste(..rr.label.., sep = "~`,`~")),
                                            label.x = 1) +
  geom_point(pch=21, size = 10)+
  theme(text = element_text(size = 32),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 28),
        strip.text = element_text(size = 28),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  labs(color="Source") +
  ylab("Observed") + xlab("Expected") + 
  scale_x_continuous(breaks=seq(0,30,by=10)) +
  geom_hline(yintercept = c(0.99, 8.3), linetype = "dotted") +
  facet_grid(~method_f, scales = "free_x") 
ls


head(tab)
tab$dis <- tab$expected - tab$observed
tab[is.na(tab)] <- 0

ptab <- select(tab, Species, expected, host)
ptab$method <- c("Expected")
colnames(ptab)[2] <- c("abundance")
head(ptab) 

ctab <- select(tab, Species, observed, method, host)
colnames(ctab)[2] <- c("abundance")
head(ctab)
c <- as.data.frame(rbind(ptab, ctab))
head(c)
c <- c %>% filter(abundance > 0)

c$method_f <- factor(c$method, levels=c("Expected", "porefile", "EMU", "wf-metagenomics"))
c$method <- NULL

wide <- reshape(c, idvar = c("Species", "host"), timevar = "method_f", direction = "wide")
wide[is.na(wide)] <- 0

abutab <- wide
write.table(abutab, "../../output_tables/compare_detection_abundance.tsv", sep = "\t", row.names = F, quote = F)

abutab %>% filter(abundance.porefile == 0)
abutab %>% filter(abundance.EMU == 0)
abutab %>% filter(`abundance.wf-metagenomics` == 0)

set.seed(227)
w <- wide[sample(nrow(wide), 20), ]

library(data.table)
long <- melt(setDT(w), id.vars = c("Species", "host"), variable.name = "method")
colnames(long) <- c("Species", "host", "method", "value")
long$method <- gsub("abundance.", "", long$method)
long$method_f <- factor(long$method, levels=c("Expected", "porefile", "EMU", "wf-metagenomics"))
long <- long %>% filter(value > 0)
long <- as.data.frame(long)

coul <- brewer.pal(8, "Spectral")
coul <- colorRampPalette(coul)(10)

head(long)
long$value <- as.numeric(long$value)

q <- ggplot(long, aes(x=method_f, y=Species, fill=value))+
  geom_tile(colour="white", size=0.5) +
  scale_fill_gradient2(low = "#f6f4d2", mid = "#f6bd60", high = "#e63946") +
  theme(legend.position="bottom", 
        legend.text=element_text(size=14),
        legend.title = element_text(size = 24),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(size = 32, face = "italic"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 28, angle = 50, hjust = 1, vjust = 1),
        strip.background = element_rect(color="black", fill="white", size=.1, linetype="solid")) +
  labs(fill='Abundance (%)') 
q

g <- ggarrange(ls, p, nrow = 2, common.legend = T, legend = "bottom", labels = c("B", "C"),
               font.label = list(size = 40, color = "black"))

fig2 <- ggarrange(q, g, ncol = 2, heights = c(0.5, 0.5), widths = c(0.5, 0.5), labels = c("A", ""),
                  font.label = list(size = 40, color = "black"))

dir.create("../../fig2")
png("../../fig2/Fig2.png", res = 600, height = 50, width = 65, units = "cm")
fig2
dev.off()

########
##RMSE##
########

p <- as.data.frame(rmse(wide$abundance.Expected, wide$abundance.porefile))
row.names(p) <- c("porefile")
colnames(p) <- c("RMSE")
e <- as.data.frame(rmse(wide$abundance.Expected, wide$abundance.EMU))
row.names(e) <- c("EMU")
colnames(e) <- c("RMSE")
w <- as.data.frame(rmse(wide$abundance.Expected, wide$`abundance.wf-metagenomics`))
row.names(w) <- c("wf-metagenomics")
colnames(w) <- c("RMSE")
rmse <- as.data.frame(rbind(p, e, w))
rmse$method <- row.names(rmse)
row.names(rmse) <- NULL

write.table(rmse, "../../output_tables/methods_rmse_HCC.tsv", sep = "\t", row.names = F, quote = F)

