
library("phyloseq")
library("microbiomeutilities") # BiocManager::install("microbiome") AND #devtools::install_github("microsud/microbiomeutilities")


###############################################################
########## preprocessing for the mock experiments #############
###############################################################

# Import the raw phyloseq mock_raw object 

Runs_raw <- readRDS("R_objects/Runs_raw.rds")

library("tidyverse")
taxM <- data.frame(phyloseq::tax_table(Runs_raw));taxM
tax.cleanM <- data.frame(row.names = row.names(taxM),
                         Kingdom = str_replace(taxM[,1], "k__",""),
                         Phylum = str_replace(taxM[,2], "p__",""),
                         Class = str_replace(taxM[,3], "c__",""),
                         Order = str_replace(taxM[,4], "o__",""),
                         Family = str_replace(taxM[,5], "f__",""),
                         Genus = str_replace(taxM[,6], "g__",""),
                         Species = str_replace(taxM[,7], "s__",""),
                         stringsAsFactors = FALSE)
for (i in 1:7){ tax.cleanM[,i] <- as.character(tax.cleanM[,i])}

tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "g__", " ")) %>% 
  mutate(Species=str_replace(Species, "g__", " "))%>% 
  mutate(Genus=str_replace(Genus, "f__", " "))%>% 
  mutate(Family=str_replace(Family, "o__", " "))%>% 
  mutate(Species=str_replace_all(Species, "_", " "));tax.cleanM


  phyloseq::tax_table(Runs_raw) <- as.matrix(tax.cleanM) ; phyloseq::tax_table(Runs_raw)


################################################################################
######### Summary stats ########

library(microbiome) # BiocManager::install("microbiome")
library(microbiomeutilities) #remotes::install_github("microsud/microbiomeutilities")

summarize_phyloseq(Runs_raw)

Dep_raw<-plot_read_distribution(Runs_raw, groups = "YEAR", 
                                plot.type = "histogram")+
  theme_biome_utils()+
  scale_x_continuous(trans='log10',limits=c(1, 10000000))+
  scale_fill_manual(values=c("#111111"))+
  geom_vline(xintercept = 2500, colour = "black", linetype="dashed")+
  theme(legend.position="none")+
  labs(x = "", y = "Count");Dep_raw

Runs_rawP <- prune_samples(sample_sums(Runs_raw) >= 3100, Runs_raw)

summarize_phyloseq(Runs_rawP)

Dep_P<-plot_read_distribution(Runs_rawP, groups = "YEAR", 
                              plot.type = "histogram")+
  theme_biome_utils()+
  scale_x_continuous(trans='log10', limits=c(1, 1000000))+
  scale_fill_manual(values=c("#111111"))+ 
  geom_vline(xintercept = 2500, colour = "black", linetype="dashed")+
  theme(legend.position="none")+
  labs(x = "Reads per samples", y = "Count");Dep_P


library("cowplot")
depth<-plot_grid(Dep_raw+theme(legend.position="none"),
                 Dep_P+theme(legend.position="none"), 
                 align="vh",
                 labels = c("A", "B"),
                 hjust = -1,
                 vjust= 2,
                 nrow = 2);depth

depth_final<-plot_grid(depth, ncol = 1, rel_heights = c(0.8, .05))
depth_final

library("ggpubr")
ggsave(file="Figures/FigS3_depth_runs.pdf", 
       width=8, height=5, units="in", dpi=300)

#save the basic Phyloseq object


#Runs_rawP = prune_taxa(taxa_sums(Runs_rawP) > 9, Runs_rawP); Runs_rawP

saveRDS(Runs_rawP, file = "R_objects/Runs_rawP.rds")

################################################################################

#For the mock experiments, the data were first rarefied
set.seed(1235)
Runs_R = rarefy_even_depth(Runs_rawP,
                           sample.size=10000,
                           rngseed=1024,
                           replace=FALSE)


summarize_phyloseq(Runs_R)

RunsRP = prune_taxa(taxa_sums(Runs_R) > 9, Runs_R); RunsRP

summarize_phyloseq(RunsRP)

# Because singletons still remains after rarefaction, taxa not supported with 
# at least 10 ASVs, were removed from the final object. 

saveRDS(RunsRP, "R_objects/RunsR_R.rds")



