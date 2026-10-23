
library("readxl")      
library("tibble")
library("vegan")
library("ggstar")
library("patchwork")
library("ggpubr")
library("cowplot")
library("phyloseq") #BiocManager::install("phyloseq")

Runs_PRF <- readRDS("R_objects/RunsR_R.rds")


library(MicrobiotaProcess)
classtaxa <- get_taxadf(obj=Runs_PRF, taxlevel=2)

pclass2 <- ggbartax(obj=classtaxa, count=TRUE, plotgroup=TRUE, facetNames="YEAR") +
  xlab(NULL) +
  ylab("count reads") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
pclass2


##################################################################
################# compute and plot stacked bar charts  ###########
##################################################################

### Extract the data table from the phyloseq object ##############
Runs_PRF_dat <- psmelt(Runs_PRF)

#write.csv(Runs_PRF_dat, "Runs_PRF_dat.csv")

str(Runs_PRF_dat)

Runs_PRF_dat$DOY<-as.numeric( as.numeric(as.character(Runs_PRF_dat$DOY)))

head(Runs_PRF_dat)
str(Runs_PRF_dat)

library("forcats")

########## filter the taxonomy to keep the most abundant only ####
str(Runs_PRF_dat)

predefined_genus = c("Alternaria",
                     "Botrytis",
                     "Blumeria",
                     "Cladosporium",
                     "Epicoccum",
                     "Fusarium",
                     "Ganoderma",
                     "Hyaloperonospora",
                     "Mycosphaerella",
                     "Neoascochyta",
                     "Neosetophoma",
                     "Peniophora",
                     "Peronospora",
                     "Pseudopithomyces",
                     "Pyrenophora",
                     "Rhodotorula",
                     "Septoriella",
                     "Stemphylium",
                     "Vishniacozyma"
                      );predefined_genus

#################################### 

library("tidyverse")
library("ggtext")
library("dplyr")

Runs_PRF_dat2 <-  Runs_PRF_dat %>%
  mutate(Genus = case_when(
    Genus %in% 
      predefined_genus ~ Genus, 
    TRUE ~ "Other"
  ));Runs_PRF_dat2

tail(Runs_PRF_dat2)
head(Runs_PRF_dat2)

################# Filter the relative abundance by year and DOY ####### 

otu_rel_abund <- Runs_PRF_dat2 %>%
  group_by(YEAR, DOY) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  select(-Abundance) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
               names_to="level",
               values_to="taxon");head(otu_rel_abund)

##################  #################

run_abund<-otu_rel_abund %>%
  filter(level=="Genus") %>%
  group_by(YEAR, DOY, taxon) %>%
  dplyr::summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(YEAR, DOY, taxon) %>%
  dplyr::summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = factor(taxon,
                        levels=c("Alternaria",
                                "Botrytis",
                                "Blumeria",
                                "Cladosporium",
                                "Epicoccum",
                                "Fusarium",
                                "Ganoderma",
                                "Hyaloperonospora",
                                "Mycosphaerella",
                                "Neoascochyta",
                                "Neosetophoma",
                                "Peniophora",
                                "Peronospora",
                                "Pseudopithomyces",
                                "Pyrenophora",
                                "Rhodotorula",
                                "Septoriella",
                                "Stemphylium",
                                "Vishniacozyma",
                                "Other")))

# plot the stacked bar chart 

run_stacked <- ggplot(data = run_abund, 
                      aes(x = DOY, 
                          y = mean_rel_abund, 
                          fill = taxon)) +
  geom_col(colour = "black", width = 2, linewidth = 0.1) +
  facet_wrap(vars(YEAR), nrow = 2) +
  labs(x = "Day of the year",
       y = "Relative abundance (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits = c(150, 231), n.breaks = 15, expand = c(0, 0)) +
  theme(legend.text = element_text(face = "italic"),
        legend.position = "top") +
  guides(
    fill = guide_legend(
      title.position = "top",
      keywidth = 0.6, 
      keyheight = 0.7, 
      ncol = 6
    )
  ) +
  theme(axis.text.x = element_text(angle = 60, 
                                   vjust = 0.5, 
                                   hjust = 0.4))
library(wesanderson)
library(dplyr)
library(ggnewscale)
library(RColorBrewer)

str(Runs_PRF_dat)

# 1. Un seul BBCH par site/date
bbch_site <- Runs_PRF_dat %>%
  group_by(YEAR, DOY, Site) %>%
  summarise(BBCH = unique(BBCH), .groups = "drop")

# 2. Le mode BBCH parmi les sites pour chaque DOY
bbch_strip <- bbch_site %>%
  group_by(YEAR, DOY) %>%
  summarise(BBCH = as.character(names(sort(table(BBCH), decreasing=TRUE)[1])), .groups="drop")


print(bbch_strip, n = 70)

# Palette et levels pour l'affichage
bbch_levels <- sort(unique(as.numeric(as.character(bbch_strip$BBCH)))) 
bbch_strip$BBCH <- factor(bbch_strip$BBCH, levels=bbch_levels, ordered=TRUE)
n_stages <- length(bbch_levels)

palette_bbch <- c(
  wesanderson::wes_palette("Zissou1", 4, "discrete"),
  wesanderson::wes_palette("Darjeeling1", 5, "discrete"),
  wesanderson::wes_palette("GrandBudapest1", 4, "discrete")
)



run_stacked_C <- run_stacked +
  scale_fill_manual(
    values = c("#9E0142", 
               "#D53E4F", 
               "#FDAE61", 
               "#FEE08B", 
               "#FFFFBF", 
               "#E6F598", 
               "#ABDDA4", 
               "#66C2A5", 
               "#3288BD", 
               "#5E4FA2", 
               "#2D004B", 
               "#01665E", 
               "#E0E0E0",              
               "#35978F",                
               "#F46D43",   
               "#B9F9F5",
               "#DE77AE", 
               "#80CDC1", 
               "#A7EAE5", 
               "#a9a9a9"

               ),
    name = "Genus",
    guide = guide_legend(order = 1)
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = bbch_strip,
    aes(x = DOY, y = 106, fill = BBCH, group=YEAR),
    width = 1.8, height = 4,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = palette_bbch,
    name = "Growth stage (BBCH)",
    guide = guide_legend(order = 2, title.position = "top")
  )+
  guides(
    fill = guide_legend(
      title.position = "top",
      keywidth = 0.6, 
      keyheight = 0.7, 
      ncol = 3)
    );run_stacked_C




################################################################################
############### with funguild ##################################################

Runs_PRF_G_melt <- psmelt(Runs_PRF)
head(Runs_PRF_G_melt)

FunGuildData <- read.csv("Data/runs/FunGuildData.csv", sep=";")%>% 
  mutate(Species=taxon)

head(FunGuildData)

# Nettoyer la colonne Species dans Runs_PRF_G_melt
Runs_PRF_G_melt <- Runs_PRF_G_melt %>%
  mutate(Species = str_replace(Species, " cluster$", ""))

Runs_PRF_G_melt_2 <- merge(Runs_PRF_G_melt, FunGuildData, by = "Species") %>% 
  select(-c("taxon"))
head(Runs_PRF_G_melt_2)

unique(data.frame(Species = subset(Runs_PRF_G_melt_2, Phenotype == "Saprotroph")$Species))
# 
# unique(data.frame(Species = subset(Runs_PRF_G_melt_2, is.na(simpleFunct))$Species))

#write.csv(Runs_PRF_G_melt_2, "Data/runs/melt_guild.csv")

################# Filter the relative abundance by year and DOY ####### 

otu_rel_abund_guild_2021 <- Runs_PRF_G_melt_2 %>%
  subset(YEAR==2021) %>% 
  group_by(BBCH_grouped) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  select(-Abundance) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", 
                 "Species", "growthForm", "Function", 
                 "Guild_1", "Phenotype", "trophic_mode"),
               names_to="level",
               values_to="taxon");head(otu_rel_abund_guild_2021)


##################  #################

run_abund_guild_2021 <- otu_rel_abund_guild_2021 %>%
  filter(level=="Phenotype") %>%
  group_by(BBCH_grouped, taxon) %>%
  dplyr::summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(BBCH_grouped, taxon) %>%
  dplyr::summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") 

# plot the stacked bar chart 

PIE_Guild_2021<-run_abund_guild_2021 %>% 
  group_by(BBCH_grouped) %>%
  mutate(BBCH_label = paste("BBCH", BBCH_grouped)) %>% 
  mutate(ypos = sum(mean_rel_abund) / 2) %>% 
  ggplot(aes(x = 2, y = mean_rel_abund, fill = taxon)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 200) +
  xlim(0.5, 2.5) +
  theme_void() +
  theme(strip.text = element_blank())+
  facet_wrap(vars(BBCH_grouped), nrow = 1)+
  labs(x="", y="") +
  scale_fill_manual(values=c( '#D9E6EBFF',
                              "#E9967A",
                              '#E37E00FF',
                              '#53869F',
                              '#A0527DFF',
                              '#9ACD32',
                              '#2D6D66FF',
                              '#CAC27EFF',
                              "#0A9F9D",
                              '#68228B',
                              '#FFD200FF',
                              '#C10534FF',
                              "#3F5151",
                              "#1A477FFF",
                              "#90353BFF",
                              '#228B22',
                              "#39312F", 
                              "#0F0D0E",
                              "#E6A0C4",
                              "#446455",
                              "#D3D4D8",
                              "#9D9C85",
                              "#5BBCD6",
                              "#BDC881",
                              "#AF4E24",
                              "#9986A5",
                              "#6C8645",
                              "grey"
                             ),
                   name = "Phenotype")+
  geom_text(aes(x = 0.5, y = ypos, label = str_wrap(BBCH_label, width = 5)), 
            color = "black", size = 2.5)+
  theme(legend.position="right")+
  theme(legend.box.margin = margin(6, 10, 6, 10))+
  guides(fill= guide_legend(keywidth = 0.6, 
                            keyheight = 0.7, 
                            ncol=2));PIE_Guild_2021


################################################################################
otu_rel_abund_guild_2022 <- Runs_PRF_G_melt_2 %>%
  subset(YEAR==2022) %>% 
  group_by(BBCH_grouped) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  select(-Abundance) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", 
                 "Species", "growthForm", "Function", 
                 "Guild_1", "Phenotype"),
               names_to="level",
               values_to="taxon");head(otu_rel_abund_guild_2022)

##################  #################

run_abund_guild_2022 <- otu_rel_abund_guild_2022 %>%
  filter(level=="Phenotype") %>%
  group_by(BBCH_grouped, taxon) %>%
  dplyr::summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(BBCH_grouped, taxon) %>%
  dplyr::summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") 

# plot the stacked bar chart 

PIE_Guild_2022<-run_abund_guild_2022 %>% 
  group_by(BBCH_grouped) %>%
  mutate(BBCH_label = paste("BBCH", BBCH_grouped)) %>% 
  mutate(ypos = sum(mean_rel_abund) / 2) %>% 
  ggplot(aes(x = 2, y = mean_rel_abund, fill = taxon)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 200) +
  xlim(0.5, 2.5) +
  theme_void() +
  theme(strip.text = element_blank())+
  facet_wrap(vars(BBCH_grouped), nrow = 1)+
  labs(x="", y="") +
  scale_fill_manual(values=c('#D9E6EBFF',
                              "#E9967A",
                              '#E37E00FF',
                              '#53869F',
                              '#A0527DFF',
                              '#9ACD32',
                              '#2D6D66FF',
                              '#CAC27EFF',
                              "#0A9F9D",
                              '#68228B',
                              '#FFD200FF',
                              '#C10534FF',
                              "#3F5151",
                              "#1A477FFF",
                              "#90353BFF",
                              '#228B22',
                              "#39312F",
                              "#0F0D0E",
                              "#E6A0C4",
                              "#446455",
                              "#D3D4D8",
                              "#9D9C85",
                              "#5BBCD6",
                              "#BDC881",
                              "#AF4E24",
                              "#9986A5",
                              "#6C8645",
                              "grey"),
                     name = "Phenotype")+
  geom_text(aes(x = 0.5, y = ypos, label = str_wrap(BBCH_label, width = 5)), 
            color = "black", size = 2.5)+
  theme(legend.position="right")+
  theme(legend.box.margin = margin(6, 10, 6, 10))+
  guides(fill= guide_legend(keywidth = 0.6, 
                            keyheight = 0.7, 
                            ncol=2));PIE_Guild_2022

#### arrange the 3 plots in a single figure. 

ggarrange(ggarrange(run_stacked_C, labels = "A"), 
          nrow = 2, 
          heights = c(2,1),
          ggarrange
          (PIE_Guild_2021, PIE_Guild_2022, nrow = 2, 
            labels = c("B", "C"), 
            common.legend = TRUE, 
            legend="right"))


##### save the figure 

ggsave(file="Figures/Fig3_Run_stacked.pdf", width=8.5, 
       height=10.5, units="in", dpi=900)
