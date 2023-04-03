library(phyloseq)
library(ggplot2)
library(stringr)
library(dplyr)
library(Metrics)
library(RColorBrewer)

########################################
######Output pre-processing#############
########################################

setwd("./raw/results_silva_zymo_noyacrd_simdefault_polishtop5/")

list.files()

############
##porefile##
############

#generate metadata from file names
c1 <- as.data.frame(read.table("COUNTS_polished.tsv", sep = "\t", header = T))
t1 <- as.data.frame(read.table("TAXCLA_polished.tsv", sep = "\t", header = T))
name <- as.data.frame(colnames(c1))
df <- as.data.frame(str_split_fixed(name$`colnames(c1)`, "_", 3))
df1 <- as.data.frame(cbind(name, df$V1))
colnames(df1) <- c("file", "id")
head(df1)

meta1 <- df1
meta1$round <- c('polish')
row.names(meta1) <- meta1$file
head(meta1)
dim(meta1)

#load count table and taxonomy to a phyloseq object
c1 <- as.matrix(read.table("COUNTS_polished.tsv", sep = "\t", header = T))
t1 <- as.matrix(read.table("TAXCLA_polished.tsv", sep = "\t", header = T))

OTU1 = otu_table(c1, taxa_are_rows = TRUE)
TAX1 = tax_table(t1)
samples = sample_data(meta1)
OTU1
TAX1

ps1 = phyloseq(OTU1, TAX1, samples) #phyloseq object
ps1

rel <- psmelt(ps1) #get the read count per taxa

#add column with an alternative ID to the sample
rel <- rel %>%
  mutate(id2 = case_when(
    endsWith(id, "mock1") ~ "A",
    endsWith(id, "mock2") ~ "B",
    endsWith(id, "mock3") ~ "C",
    endsWith(id, "mock4") ~ "D",
    endsWith(id, "mock5") ~ "E",
    endsWith(id, "mock6") ~ "F",
    endsWith(id, "mock7") ~ "G",
    endsWith(id, "mock8") ~ "H",
    endsWith(id, "mock9") ~ "I",
    endsWith(id, "mock10") ~ "J"))
head(rel)
unique(rel$Species)

rel <- select(rel, id, Abundance, Species) #reduce the table 

dir.create("../../inputs")
write.table(rel, "../../inputs/prefile_relative_abundance.tsv", row.names = F, quote = F, sep = "\t")

###################
##wf-metagenomics##
###################

setwd("../output_wf/")

df <- read.table("wf-metagenomics-counts.tsv", sep = "\t", header = T)
s <- colSums(df[,-1])
dff <- df[,-1]/s
row.names(dff) <- df$X
dff$Species <- row.names(dff)
row.names(dff) <- NULL
head(dff)

library(data.table)
dff = melt(df, id.vars = c("X"),
           measure.vars = c("barcode01", "barcode02", "barcode03",
                            "barcode04", "barcode05", "barcode06",
                            "barcode07", "barcode08", "barcode09", "barcode10"))
head(dff)
dff$variable <- gsub("barcode", "mock", dff$variable)
colnames(dff) <- c("Species", "id", "observed")
dff$id <- gsub("barcode", "mock", dff$id)
dff$Species <- gsub("Bacillus spizizenii", "Bacillus subtilis", dff$Species)
dff <- select(dff, id, observed, Species)
colnames(dff)[2] <- c("Abundance")

write.table(dff, "../../inputs/wf-metagneomics_abundance_tab.tsv", sep = "\t", row.names = F, quote = F)


#######
##EMU##
#######

library(fs) 
library(magrittr)
library(vroom)
library(data.table) 
library(formattable) 
library(tidyverse)  
library(tidyselect)
library(writexl)

setwd("../emu_output/")
list.files()
# read file path
all_paths <-
  list.files(path = "./",
             pattern = "*.tsv",
             full.names = TRUE)
# read file content
all_content <-
  all_paths %>%
  lapply(read.table,
         header = TRUE,
         sep = "\t",
         encoding = "UTF-8")


# read file name
all_filenames <- all_paths %>%
  basename() %>%
  as.list()

# combine file content list and file name list
all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)

# unlist all lists and change column name
all_result <- rbindlist(all_lists, fill = T)
dim(all_result)
# change column name
names(all_result)[15] <- "File.Path"
all_result$File.Path <- gsub("_rel-abundance.tsv", "", all_result$File.Path)
head(all_result)

emu <- select(all_result, File.Path, abundance, species)
emu$abundance <- gsub("abundance", "", emu$abundance)
colnames(emu) <- c("id", "Abundance", "Species")
emu$Abundance <- as.numeric(emu$Abundance)
emu$id <- gsub("filt_", "", emu$id)
emu <- as.data.frame(emu)
emu <- emu[complete.cases(emu),]
emu <- select(emu, id, Abundance, Species)

write.table(emu, "../../inputs/emu_abundance_tab.tsv", sep = "\t", row.names = F, quote = F)

######################################################################################################
#######################Process porefile, EMU and wf-output############################################
######################################################################################################

#porefile
setwd("../../inputs")
porefile <- read.table("prefile_relative_abundance.tsv", sep = "\t", header = T) 
porefile$method <- c("porefile")
head(porefile)
pr <- porefile %>% group_by(id) %>%
  summarise(count = sum(Abundance))
porefile <- merge(porefile, pr, by = "id")
porefile$abundance <- porefile$Abundance/porefile$count #convert to relative abundance 
porefile <- select(porefile, id, abundance, Species, method)
colnames(porefile)[2] <- c("Abundance") 

#EMU
emu <- read.table("emu_abundance_tab.tsv", sep = "\t", header = T)
emu$method <- c("EMU")
head(emu)

#wf-metagenomics
wf <- read.table("wf-metagneomics_abundance_tab.tsv", sep = "\t", header = T)
w <- wf %>% group_by(id) %>%
  summarise(count = sum(Abundance))
wf <- merge(wf, w, by = "id")
wf$abundance <- wf$Abundance/wf$count #convert to relative abundance 
wf$method <- c("wf-metagenomics")
wf <- select(wf, id, abundance, Species, method)
colnames(wf)[2] <- c("Abundance") 

#combine outputs
head(porefile)
head(wf)
head(emu)
tab <- rbind(porefile, emu, wf)
tab$Species <- gsub("Bacillus spizizenii", "Bacillus subtilis", tab$Species) #in the SILVA databse Bacillus subtilis subsp. spizizenii for this reason it is replaced by Bacillus subtilis 
tabb <- tab
tabb <- tabb %>% group_by(id, Species, method) %>%
  summarise(Abundance = sum(Abundance))


#Theoretical composition of the LCC
z <- read.table("../zmock.tsv", sep = "\t", header = T)
z$expected <- as.numeric(gsub(",", ".", z$expected))
merged <- merge(z, tab, by = "Species")

merged <- merged %>% group_by(id, Species, expected, method) %>%
  summarise(sum = sum(Abundance))
merged$observed <- merged$sum*100
head(merged)

m <- merged

target <- c("Listeria monocytogenes", "Pseudomonas aeruginosa", "Bacillus subtilis", "Escherichia coli",
            "Salmonella enterica","Limosilactobacillus fermentum", "Enterococcus faecalis", "Staphylococcus aureus")

int <- filter(merged, Species %in% target)
neg <- tab[!tab$Species %in% merged$Species, ]
neg <- neg %>% group_by(id, method) %>%
  summarise(sum = sum(Abundance))
neg$Species <- c("Others")
neg$expected <- as.numeric(c("0"))
neg$observed <- neg$sum*100
neg <- select(neg, id, Species, expected, method, sum, observed)

int <- as.data.frame(rbind(int, neg))
unique(int$id)
int <- int%>%
  mutate(id2 = case_when(
    endsWith(id, "mock1") ~ "A",
    endsWith(id, "mock2") ~ "B",
    endsWith(id, "mock3") ~ "C",
    endsWith(id, "mock4") ~ "D",
    endsWith(id, "mock5") ~ "E",
    endsWith(id, "mock6") ~ "F",
    endsWith(id, "mock7") ~ "G",
    endsWith(id, "mock8") ~ "H",
    endsWith(id, "mock9") ~ "I",
    endsWith(id, "mock01") ~ "A",
    endsWith(id, "mock02") ~ "B",
    endsWith(id, "mock03") ~ "C",
    endsWith(id, "mock04") ~ "D",
    endsWith(id, "mock05") ~ "E",
    endsWith(id, "mock06") ~ "F",
    endsWith(id, "mock07") ~ "G",
    endsWith(id, "mock08") ~ "H",
    endsWith(id, "mock09") ~ "I",
    endsWith(id, "mock10") ~ "J"))
head(int)
dir.create("../output_tables")
write.table(int, "../output_tables/figure_AB.table.tsv", sep = "\t", row.names = F, quote = F)

########################
#############

coul <- brewer.pal(8, "Spectral")
coul <- colorRampPalette(coul)(12)

p <- ggplot(int, aes(x = id2, y = sum, fill = Species))
g <- p +  geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = coul) +
  facet_grid(method~id2, scales = "free_x") +
  theme(legend.position="bottom",
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40),
        legend.text = element_text(size = 36),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  guides(fill=guide_legend(ncol = 2)) + ggtitle("") +
  xlab("") + ylab("Relative Abundance") + labs(fill = "Species")
g


library(ggpubr)
int <- int[-which(int$Species == "Others"),] #remove other classification than the LCC components

ls <- ggscatter(int, x = "expected", y = "observed", color = "Species", size = 10, alpha = 1, 
                add = "reg.line", cor.coef.size = 25, add.params = list(color = "black", fill = "lightgray"),
                conf.int = TRUE) + stat_cor(size = 14, aes(label = paste(..rr.label.., sep = "~`,`~")),
                                            label.x = 1) +
  geom_point(pch=21, size = 10)+
  theme(text = element_text(size = 46),
        legend.position = "bottom",
        legend.text = element_text(size = 36),
        legend.title = element_blank(),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  scale_color_manual(values = coul) +
  guides(color=guide_legend(ncol = 2)) +
  ylab("Observed") + xlab("Expected") +
  facet_wrap(~method, scales = "free_x", ncol = 5) +
  ggtitle("")

ls

fig1AB <- ggarrange(g, ls, ncol = 2, labels = c("A", "B"), font.label = list(size = 40, color = "black"), 
                    heights = c(0.5,0.5), widths = c(0.5,0.5))

###############################################################################################
###################################################################################

#porefile read count
p_rc <- read.table("../raw/results_silva_zymo_noyacrd_simdefault_polishtop5/NanoPlots/summary.tsv", sep = "\t", header = T)
porefile <- read.table("prefile_relative_abundance.tsv", sep = "\t", header = T) 
porefile$Species <- gsub("Bacillus spizizenii", "Bacillus subtilis", porefile$Species)
porefile$method <- c("porefile")
head(porefile)

p <- porefile %>% group_by(id) %>%
  summarise(sum = sum(Abundance))

p$FN <- 40000-p$sum

target <- c("Listeria monocytogenes", "Pseudomonas aeruginosa", "Bacillus subtilis", "Escherichia coli",
            "Salmonella enterica","Limosilactobacillus fermentum", "Enterococcus faecalis", "Staphylococcus aureus")

int <- filter(porefile, Species %in% target)
TP <- int %>% group_by(id) %>%
  summarise(sum = sum(Abundance))
colnames(TP) <- c("id", "TP")
head(TP)
p <- merge(p, TP, by = "id")
head(p)

neg <- porefile[!porefile$Species %in% int$Species,] #no false positives with porefile
p$FP <- as.numeric(c("0"))
p$method <- c("porefile")


#emu read count
e_rc <- all_result
head(e_rc)
e_rc <- select(e_rc, File.Path, estimated.counts)
e_rc <- as.data.frame(e_rc)
e_rc <- e_rc[complete.cases(e_rc),]
colnames(e_rc) <- c("id", "count")
e_rc <- e_rc %>% group_by(id) %>%
  summarise(count = sum(count))
e_rc <- as.data.frame(e_rc)

emu <- read.table("emu_abundance_tab.tsv", sep = "\t", header = T) 
emu$Species <- gsub("Bacillus spizizenii", "Bacillus subtilis", emu$Species)
emu$method <- c("EMU")

em <- merge(emu, e_rc, by = "id")
head(em)
em$reads <- em$Abundance*em$count
em$reads <- round(em$reads, 0)
em$count <- NULL

e <- em %>% group_by(id) %>%
  summarise(sum = sum(reads))
e$FN <- 40000-e$sum

int <- filter(em, Species %in% target)
TP <- int %>% group_by(id) %>%
  summarise(sum = sum(reads))
colnames(TP) <- c("id", "TP")
e <- merge(e, TP, by = "id")

neg <- em[!em$Species %in% int$Species,]
FP <- neg %>% group_by(id) %>%
  summarise(sum = sum(reads))
colnames(FP) <- c("id", "FP")
e <- merge(e, FP, id = "id")
e$method <- c("EMU")


#wf-metagenomics read count
wf <- read.table("wf-metagneomics_abundance_tab.tsv", sep = "\t", header = T)
wf$method <- c("wf-metagenomics")
head(wf)

wf$Species <- gsub("Bacillus spizizenii", "Bacillus subtilis", wf$Species)

w <- wf %>% group_by(id) %>%
  summarise(sum = sum(Abundance))

w$FN <- 40000-w$sum

int <- filter(wf, Species %in% target)
TP <- int %>% group_by(id) %>%
  summarise(sum = sum(Abundance))
colnames(TP) <- c("id", "TP")
w <- merge(w, TP, by = "id")

neg <- wf[!wf$Species %in% int$Species,]
FP <- neg %>% group_by(id) %>%
  summarise(sum = sum(Abundance))
colnames(FP) <- c("id", "FP")
w <- merge(w, FP, by = "id")
w$method <- c("wf-metagenomics")

tab <- as.data.frame(rbind(p, e, w))

tab$precision <- tab$TP/(tab$TP+tab$FP)*100
tab$recall <- tab$TP/(tab$TP+tab$FN)*100

emup <- tab[which(tab$method == "EMU"),]
mean(emup$precision)
sd(emup$precision)
mean(emup$recall)
sd(emup$recall)

porp <- tab[which(tab$method == "porefile"),]
mean(porp$precision)
sd(porp$precision)
mean(porp$recall)
sd(porp$recall)

wfp <- tab[which(tab$method == "wf-metagenomics"),]
mean(wfp$precision)
sd(wfp$precision)
mean(wfp$recall)
sd(wfp$recall)


write.table(tab, "../output_tables/precision_recall_table.tsv", sep = "\t", row.names = F, quote = F)

coul <- brewer.pal(3, "Spectral")

prec <- ggplot(tab, aes(x=method, y=precision, fill = method)) +
  geom_boxplot(alpha = 0) +
  theme(legend.position="none", 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 45, angle = 30, hjust = 0.8, vjust = 0.8, face = "bold"),
        strip.text = element_text(size = 32, face = "bold"),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  ylab("PRECISION (%)") + scale_y_continuous(limits = c(0, 120), breaks = c(0, 20, 40, 60, 80, 100)) +
  geom_signif(comparisons = list(c("porefile", "EMU")), map_signif_level = T, y_position = 110, size = 1, textsize = 15) +
  geom_signif(comparisons = list(c("porefile", "wf-metagenomics")), map_signif_level = T, y_position = 105, size = 1, textsize = 15) +
  scale_fill_manual(values = coul) + guides(color=guide_legend(title="")) +
  geom_jitter(size = 10, pch=21, aes(fill = method, alpha = 0.5)) 
prec

rec <- ggplot(tab, aes(x=method, y=recall)) +
  geom_boxplot(alpha = 0) +
  theme(legend.position="none", 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 40),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 45, angle = 30, hjust = 0.8, vjust = 0.8, face = "bold"),
        strip.text = element_text(size = 32, face = "bold"),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  ylab("RECALL (%)") + scale_y_continuous(limits = c(0, 120), breaks = c(0, 20, 40, 60, 80, 100)) +
  geom_signif(comparisons = list(c("porefile", "EMU")), map_signif_level = T, y_position = 90, size = 1, textsize = 15) +
  geom_signif(comparisons = list(c("porefile", "wf-metagenomics")), map_signif_level = T, y_position = 105, size = 1, textsize = 15) +
  scale_fill_manual(values = coul) + guides(color=guide_legend(title="")) +
  geom_jitter(size = 10, pch=21, aes(fill = method, alpha = 0.5)) 
rec

figc <- ggarrange(prec, rec, ncol = 2)

###############################################################################################################

#porefile mock components
m1 <- m[which(m$method == "porefile"),]

bs <- m1[which(m1$Species == "Bacillus subtilis"),]
r_bs <- rmse(bs$expected, bs$observed)
ef <- m1[which(m1$Species == "Enterococcus faecalis"),]
r_ef <-rmse(ef$expected, ef$observed)
ec <- m1[which(m1$Species == "Escherichia coli"),]
r_ec <- rmse(ec$expected, ec$observed)
lf <- m1[which(m1$Species == "Limosilactobacillus fermentum"),]
r_lf <- rmse(lf$expected, lf$observed)
pa <- m1[which(m1$Species == "Pseudomonas aeruginosa"),]
r_pa <- rmse(pa$expected, pa$observed)
se <- m1[which(m1$Species == "Salmonella enterica"),]
r_se <- rmse(se$expected, se$observed)
sa <- m1[which(m1$Species == "Staphylococcus aureus"),]
r_sa <- rmse(sa$expected, sa$observed)
lm <- m1[which(m1$Species == "Listeria monocytogenes"),]
r_lm <- rmse(lm$expected, lm$observed)

r_sp <- as.data.frame(rbind(r_bs, r_ef, r_ec, r_lf, r_pa, r_se, r_sa, r_lm))
r_sp$method <- c("porefile")
colnames(r_sp) <- c("RMSE")
r_sp$sample <- row.names(r_sp)
row.names(r_sp) <- NULL

#emu mock components
m2 <- m[which(m$method == "EMU"),]

bs <- m2[which(m2$Species == "Bacillus subtilis"),]
r_bs <- rmse(bs$expected, bs$observed)
ef <- m2[which(m2$Species == "Enterococcus faecalis"),]
r_ef <-rmse(ef$expected, ef$observed)
ec <- m2[which(m2$Species == "Escherichia coli"),]
r_ec <- rmse(ec$expected, ec$observed)
lf <- m2[which(m2$Species == "Limosilactobacillus fermentum"),]
r_lf <- rmse(lf$expected, lf$observed)
pa <- m2[which(m2$Species == "Pseudomonas aeruginosa"),]
r_pa <- rmse(pa$expected, pa$observed)
se <- m2[which(m2$Species == "Salmonella enterica"),]
r_se <- rmse(se$expected, se$observed)
sa <- m2[which(m2$Species == "Staphylococcus aureus"),]
r_sa <- rmse(sa$expected, sa$observed)
lm <- m2[which(m2$Species == "Listeria monocytogenes"),]
r_lm <- rmse(lm$expected, lm$observed)

r_sp2 <- as.data.frame(rbind(r_bs, r_ef, r_ec, r_lf, r_pa, r_se, r_sa, r_lm))
r_sp2$method <- c("EMU")
colnames(r_sp2) <- c("RMSE")
r_sp2$sample <- row.names(r_sp2)
row.names(r_sp2) <- NULL

#wf-metagenomics components
m3 <- m[which(m$method == "wf-metagenomics"),]

bs <- m3[which(m3$Species == "Bacillus subtilis"),]
r_bs <- rmse(bs$expected, bs$observed)
ef <- m3[which(m3$Species == "Enterococcus faecalis"),]
r_ef <-rmse(ef$expected, ef$observed)
ec <- m3[which(m3$Species == "Escherichia coli"),]
r_ec <- rmse(ec$expected, ec$observed)
lf <- m3[which(m3$Species == "Limosilactobacillus fermentum"),]
r_lf <- rmse(lf$expected, lf$observed)
pa <- m3[which(m3$Species == "Pseudomonas aeruginosa"),]
r_pa <- rmse(pa$expected, pa$observed)
se <- m3[which(m3$Species == "Salmonella enterica"),]
r_se <- rmse(se$expected, se$observed)
sa <- m3[which(m3$Species == "Staphylococcus aureus"),]
r_sa <- rmse(sa$expected, sa$observed)
lm <- m3[which(m3$Species == "Listeria monocytogenes"),]
r_lm <- rmse(lm$expected, lm$observed)

r_sp3 <- as.data.frame(rbind(r_bs, r_ef, r_ec, r_lf, r_pa, r_se, r_sa, r_lm))
r_sp3$method <- c("wf-metagenomics")
colnames(r_sp3) <- c("RMSE")
r_sp3$sample <- row.names(r_sp3)
row.names(r_sp3) <- NULL

list <- as.data.frame(rbind(r_sp, r_sp2, r_sp3))
head(list)
colnames(list)[2] <- c("method")
list <- list %>%
  mutate(Species = case_when(
    endsWith(sample, "r_bs") ~ "Bacillus subtilis",
    endsWith(sample, "r_ef") ~ "Enterococcus faecalis",
    endsWith(sample, "r_ec") ~ "Escherichia coli",
    endsWith(sample, "r_lf") ~ "Limosilactobacillus fermentum",
    endsWith(sample, "r_pa") ~ "Pseudomonas aeruginosa",
    endsWith(sample, "r_se") ~ "Salmonella enterica",
    endsWith(sample, "r_sa") ~ "Staphylococcus aureus",
    endsWith(sample, "r_lm") ~ "Listeria monocytogenes"))

write.table(list, "../output_tables/rmsa_species_table.tsv", sep = "\t", row.names = F, quote = F)


coul <- brewer.pal(4, "Spectral")


q <- ggplot(list, aes(x=method, y=Species, fill=RMSE)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_gradient2(high = "#d62828", mid = "#fcbf49", low = "#eae2b7", midpoint = 8) +
  labs(fill="RMSE") +
  theme(legend.position="right", axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(face = "italic"),
        text = element_text(size = 40), 
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 22),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 32, face = "bold", angle = 30, hjust = 1, vjust = 1),
        strip.background = element_rect(color="black", fill="white", size=.1, linetype="solid"))
q


l <- list[order(list$RMSE, decreasing = T),]
ll <- l[which(l$method == "porefile"),]
mean(ll$RMSE)
sd(ll$RMSE)

lll <- l[which(l$method == "EMU"),]
mean(lll$RMSE)
sd(lll$RMSE)

w <- l[which(l$method == "wf-metagenomics"),]
mean(w$RMSE)
sd(w$RMSE)


figCD <- ggarrange(figc, q, ncol = 2, labels = c("C", "D"), font.label = list(size = 40, color = "black"), 
                   heights = c(0.6, 0.4), widths = c(0.6, 0.4))

fig1 <- ggarrange(fig1AB, figCD, nrow = 2, heights = c(0.6, 0.4), widths = c(0.6, 0.4))

dir.create("../fig1")
png("../fig1/Fig1.png", res = 600, height = 60, width = 78, units = "cm")
fig1
dev.off()

####################################
####Supplementary Figure 2##########
####################################


head(tabb)
tabb <- select(tabb, Species, Abundance, method)
z <- read.table("../zmock.tsv", sep = "\t", header = T)
z$expected <- gsub(",", ".", z$expected)
z$expected <- as.numeric(z$expected)
z$Genus <- NULL
z$Family <- NULL
z$Abundance <- z$expected/100
z$method <- c("Expected")
z$expected <- NULL
z$id <- c("mix")
head(z)
z <- select(z, id, Species, Abundance, method)
list <- as.data.frame(rbind(z, tabb))

target <- c("Listeria monocytogenes", "Pseudomonas aeruginosa", "Bacillus subtilis", "Escherichia coli",
            "Salmonella enterica","Limosilactobacillus fermentum", "Enterococcus faecalis", "Staphylococcus aureus")

list2 <- filter(list, Species %in% target)

list2 <-list2 %>%
  mutate(paired = case_when(
    endsWith(Species, "Bacillus subtilis") ~ "1",
    endsWith(Species, "Enterococcus faecalis") ~ "2",
    endsWith(Species, "Escherichia coli") ~ "3",
    endsWith(Species, "Limosilactobacillus fermentum") ~ "4",
    endsWith(Species, "Pseudomonas aeruginosa") ~ "5",
    endsWith(Species, "Salmonella enterica") ~ "6",
    endsWith(Species, "Staphylococcus aureus") ~ "7",
    endsWith(Species, "Listeria monocytogenes") ~ "8"))
head(list2)


list2 <-list2 %>%  group_by(Species, method, paired) %>%
  summarise(mean = mean(Abundance))
list2$observed <- list2$mean*100


list2$method2 <- factor(list2$method, levels=c("Expected", "porefile", "EMU", "wf-metagenomics"))

library(viridis)
figs2a <- 
  list2 %>%
  ggplot(aes(method2, observed, color=Species)) +
  geom_dotplot(alpha = 0) +
  geom_point(size = 10)+ 
  geom_line(aes(group=paired)) +
  scale_color_viridis(option = "F", discrete = T) +
  theme(legend.position = "bottom", 
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(size = 18),
        text = element_text(size = 32), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        legend.key=element_blank(),
        axis.text.x = element_text(size = 18, face = "bold")) +
  guides(color=guide_legend(ncol=2)) +
  ylab("Abundance (%)") 
figs2a

dir.create("../sup1")
png("../sup1/FigS2.png", res = 600, height = 40, width = 40, units = "cm")
figs2a
dev.off()

