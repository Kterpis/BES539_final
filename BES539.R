library(tidyverse)
library(ggplot2)

Orthogroups <- read_tsv("~/Desktop/Statistics_PerSpecies.tsv")

first_10_rows <- Orthogroups %>% slice(1:10)

changed<- as.data.frame(t(first_10_rows))

names(changed)<- changed %>% slice(1) %>% unlist()
changed <- changed %>% slice(-1)
rownames(changed)[rownames(changed) == "1" ] = "Aplano"
rownames(changed)[rownames(changed) == "2" ] = "Cafe"
rownames(changed)[rownames(changed) == "3" ] = "Cylindro"
rownames(changed)[rownames(changed) == "4" ] = "Dinobryon"
species <- c("Aplano", "Cafe", "Cylindro", "Dinobryon")
changed$species=species


factorToNumeric <- function(f) as.numeric(levels(f))[as.integer(f)]
cols <- c(1, 2:ncol(changed))
changed[cols] <- lapply(changed[cols], factorToNumeric)
changed$species=species

ggplot(changed, aes(x=species, y=`Number of genes`, fill=species)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  scale_fill_manual(values = mycols) +
  theme_bw()+
  labs(title = "Number of genes")
       

num_of_genes <- changed %>% arrange(desc(species)) %>% mutate(lab.ypos = cumsum(`Number of genes in species-specific orthogroups` - 0.5*`Number of genes in species-specific orthogroups`))
num_of_genes

mycols <- c("#d55e00", "#009E73", "#0072b2", "#f0e442")

ggplot(num_of_genes, aes(x=species, y=`Number of genes in species-specific orthogroups`, fill = species)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  scale_fill_manual(values = mycols) +
  theme_bw()+
  labs(title = "Number of genes in species specific orthogroups")


num_in_ortho <- changed %>% arrange(desc(species)) %>% mutate(lab.ypos = cumsum(`Number of genes in orthogroups` - 0.5*`Number of genes in orthogroups`))
num_in_ortho

mycols <- c("#d55e00", "#009E73", "#0072b2", "#f0e442")

ggplot(num_in_ortho, aes(x=species, y=`Number of genes in orthogroups`, fill = species)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  scale_fill_manual(values = mycols) +
  theme_bw() +
  labs(title = "Number of genes in orthogroups")


num_of_spe_spec <- changed %>% arrange(desc(species)) %>% mutate(lab.ypos = cumsum(`Number of species-specific orthogroups` - 0.5*`Number of species-specific orthogroups`))
num_of_spe_spec

num_of_spe_spec<- c("#d55e00", "#009E73", "#0072b2", "#f0e442")


mycols <- c("#d55e00", "#009E73", "#0072b2", "#f0e442")

ggplot(num_unassign, aes(x=species, y=`Number of unassigned genes`, fill = species)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  scale_fill_manual(values = mycols) + theme_bw() +
  labs(title = "Number of Unassigned genes")

#percent_genes_ortho <- changed %>% arrange(desc(species)) %>% mutate(lab.ypos = cumsum(`Percentage of genes in orthogroups` - 0.5*`Percentage of genes in orthogroups`))
#percent_genes_ortho

mycols <- c("#d55e00", "#009E73", "#0072b2", "#f0e442")

ggplot(changed, aes(x="", y=`Number of species-specific orthogroups`, fill = species)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  geom_text(aes(label=`Number of species-specific orthogroups`) , position = position_stack(vjust = 0.5)) +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = mycols) +
  theme_void() +
  labs(title = "Number of species-specific orthogroups")

ggplot(percent_genes_ortho, aes(x="", y=`Percentage of genes in orthogroups`, fill = species)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = mycols) +
  theme_void() +
  labs(title = "Percent of genes in orthogroups")

percent_unassign_genes <- changed %>% arrange(desc(species)) %>% mutate(lab.ypos = cumsum(`Percentage of unassigned genes` - 0.5*`Percentage of unassigned genes`))
percent_unassign_genes

mycols <- c("#f0e442", "#009E73", "#0072b2", "#d55e00")


##This is a bad plot
#ggplot(percent_unassign_genes, aes(x="", y=`Percentage of unassigned genes`, fill = species)) +
  #geom_bar(width = 1, stat = "identity", color = "black") +
  #geom_text(aes(label=`Percentage of unassigned genes`) , position = position_stack(vjust = 0.5)) +
  #coord_polar("y", start = 0)+
  #scale_fill_manual(values = mycols) +
  #theme_void() +
  #labs(title = "Percent of unassigned genes")
    