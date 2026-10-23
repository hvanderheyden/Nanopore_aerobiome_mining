
library(phyloseq)
library(tidyverse)
library(ggtext)
library(dplyr)
library(vegan)
library(pheatmap)
library(ggplot2)

my_cols <- c("#899DA4", '#90353BFF', "#3B9AB2", "#74A089", 
             "#BDC881", "#C10534FF", "#F2AD00", "#000000", 
             "#C7B19C", "#446455", "#91BAB6", "#D69C4E", "grey")


mock_rare <- readRDS("R_objects/mock_R_rare.rds")

phyloseq::tax_table(mock_rare)

# Predefined species 

predefined_species <- c(
  "Botrytis porri",
  "Alternaria sect. alternata",
  "Botrytis cinerea cluster",
  "Peronospora variabilis",
  "Epicoccum nigrum",
  "Stemphylium vesicarium cluster",
  "Peronospora destructor",
  "Sclerotinia sclerotiorum cluster",
  "Alternaria sect. brassicicola",
  "Botrytis squamosa",
  "Fusarium oxysporum",
  "Hyaloperonospora brassicae"
)

mock_dat <- psmelt(mock_rare)

mock_dat <- mock_dat %>%
  mutate(Species = case_when(
    is.na(Species) | Species == "" ~ "Other",
    Species %in% predefined_species ~ Species,
    TRUE ~ "Other"
  ))


mock_abund <- mock_dat %>%
  group_by(Name, Description, Sample, Species) %>%
  dplyr::summarize(Abundance = sum(Abundance)) %>%
  ungroup() %>%
  group_by(Name, Description, Sample) %>%
  mutate(rel_abund = 100 * Abundance / sum(Abundance)) %>%
  ungroup() %>%
  group_by(Name, Description, Species) %>%
  dplyr::summarize(mean_rel_abund = mean(rel_abund)) %>%
  ungroup() %>%
  mutate(Description = factor(Description, 
                              levels = c("Theoretical", 
                                         "Minimap2",
                                         "Emu",
                                         "Bracken",
                                         "Qiime2",
                                         "Kraken")),
         Species = factor(Species, 
                          levels = c("Botrytis cinerea cluster",
                                     "Botrytis porri",
                                     "Botrytis squamosa",
                                     "Sclerotinia sclerotiorum cluster",
                                     "Alternaria sect. alternata",
                                     "Alternaria sect. brassicicola",
                                     "Peronospora variabilis",
                                     "Peronospora destructor",
                                     "Hyaloperonospora brassicae",
                                     "Epicoccum nigrum",
                                     "Stemphylium vesicarium cluster",
                                     "Fusarium oxysporum",
                                     "Other")))


write.csv(mock_abund, "figures/summary_mock.csv", row.names = FALSE)

# Plot barplot stacked - créer un plot pour chaque mock
plot_barplot <- function(data, mock_name) {
  ggplot(data = data %>% filter(Name == mock_name), 
         aes(x = Description, 
             y = mean_rel_abund, 
             fill = Species)) +
    geom_col(colour = "black", width = 0.8, linewidth = 0.1) +
    theme(legend.title = element_blank()) +
    labs(x = NULL, y = "Relative Abundance (%)") +
    theme(legend.text = element_text(face = "italic", size = 7.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = my_cols) +
    guides(fill = guide_legend(ncol = 2, keyheight = 0.4, keywidth = 0.4)) +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust = 0.4),
          plot.margin = margin(t = 10, r = 20, b = 0, l = 25, unit = "pt"))
}

mock_a <- plot_barplot(mock_abund, "mock")
mock_b <- plot_barplot(mock_abund, "mock2")
mock_c <- plot_barplot(mock_abund, "mock3")



# Préparer les données pour comparaison classifiers vs theoretical
classifiers <- c("Minimap2",
                 "Emu",
                 "Bracken",
                 "Qiime2",
                 "Kraken")

# PARAMÈTRE: Choisir la méthode de distance
# ============================================================================
# Note: La distance est calculée comme Hellinger transformation + Bray-Curtis
# (Hellinger est appliquée pour normaliser, puis Bray-Curtis pour la distance)

comparison_data <- mock_abund %>%
  pivot_wider(
    id_cols = c(Name, Species),
    names_from = Description,
    values_from = mean_rel_abund
  )

# Calculer les métriques par classifier et par groupe
metrics_by_group <- data.frame()

for (mock_name in c("mock", "mock2", "mock3")) {
  
  mock_data <- mock_abund %>%
    filter(Name == mock_name) %>%
    pivot_wider(
      id_cols = Species,
      names_from = Description,
      values_from = mean_rel_abund,
      values_fill = 0
    )
  
  for (clf in classifiers) {
    
    # Récupérer les données
    theo <- mock_data$Theoretical
    obs <- mock_data[[clf]]
    
    # Calculer les métriques
    rmse <- sqrt(mean((theo - obs)^2, na.rm = TRUE))
    # Pourcentage de reads mal classifiés = % Other
    # Récupérer la valeur de "Other" pour ce classifier depuis la ligne correspondante
    other_row <- mock_data[mock_data$Species == "Other", ]
    pct_misclassified <- as.numeric(other_row[[clf]])
    # Proportion of reads classified as Other
    proportion_other <- pct_misclassified
    
    # Distance: Hellinger transformation + Bray-Curtis
    # Normaliser les vecteurs à proportions (0-1)
    theo_norm <- theo / sum(theo)
    obs_norm <- obs / sum(obs)
    # Appliquer transformation Hellinger (sqrt(x/sum(x)))
    combined <- rbind(theo_norm, obs_norm)
    hell_transformed <- decostand(combined, method = "hellinger")
    # Puis calculer Bray-Curtis sur données transformées
    dist_value <- vegdist(hell_transformed, method = "bray")[1]
    one_minus_dist <- 1 - dist_value 
    
    metrics_by_group <- rbind(metrics_by_group, data.frame(
      Group = mock_name,
      Classifier = clf,
      RMSE = rmse,
      Proportion_Other = proportion_other,
      Distance_Similarity = as.numeric(one_minus_dist)
    ))
  }
}

# Calculer moyenne et SE par classifier (tous groupes confondus)
metrics_summary <- metrics_by_group %>%
  group_by(Classifier) %>%
  dplyr::summarize(
    Mean_RMSE = mean(RMSE),
    RMSE_se = sd(RMSE) / sqrt(n()),
    Mean_Proportion_Other = mean(Proportion_Other)/100,
    Proportion_Other_se = (sd(Proportion_Other) / sqrt(n()))/100,
    Mean_Distance_Similarity = mean(Distance_Similarity),
    Distance_Similarity_se = sd(Distance_Similarity) / sqrt(n()),
    .groups = "drop"
  )

print(metrics_summary)

# Créer les 4 plots avec style horizontal (barres)
library(gridExtra)

# Ordonner les classifiers pour affichage
metrics_summary$Classifier <- factor(metrics_summary$Classifier, 
                                     levels = rev(c("Minimap2",
                                                    "Emu",
                                                    "Bracken",
                                                    "Qiime2",
                                                    "Kraken"
                                                    )))


# Plot 1: Proportion Other (horizontal) avec errorbars

p1 <- ggplot(metrics_summary, aes(x = Mean_Proportion_Other, 
                                  y = Classifier, 
                                  fill = Classifier)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(xmin = Mean_Proportion_Other, 
                    xmax = Mean_Proportion_Other + Proportion_Other_se),
                width=0.3, color = "black", linewidth = 0.6) +
  geom_text(aes(label = paste0(round(Mean_Proportion_Other, 3))), 
            hjust = -0.95, vjust = 0.5, fontface = "bold", size = 2.5) +
  scale_fill_manual(values = c("Minimap2" ="#3A9AB2",
                               "Emu"= "#A5C2A3",
                               "Bracken" = "#E3B710",
                               "Qiime2" = "#E79805",
                               "Kraken" = "#5B1A18"),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 1.08), 
                     breaks = seq(0, 1, by = 0.2), 
                     expand = c(0, 0)) +
  theme_minimal()+
  labs(x = "False Positive Proportion")+    
  theme(axis.title.y = element_blank(),
        plot.margin = margin(t = 10, r = 30, b = 10, l = 12, unit = "pt")) 

# Plot 3: RMSE (horizontal) avec errorbars
p4 <- ggplot(metrics_summary, aes(x = Mean_RMSE, 
                                  y = Classifier,
                                  fill = Classifier)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(xmin = Mean_RMSE, 
                    xmax = Mean_RMSE + RMSE_se),
                width=0.3, color = "black", linewidth = 0.6) +
  geom_text(aes(label = round(Mean_RMSE, 3), 
                x = Mean_RMSE + RMSE_se + 0.005), 
            hjust = -0.6, vjust = 0.5, fontface = "bold", size = 2.5) +
  scale_fill_manual(values = c("Minimap2" ="#3A9AB2",
                               "Emu"= "#A5C2A3",
                               "Bracken" = "#E3B710",
                               "Qiime2" = "#E79805",
                               "Kraken" = "#5B1A18"),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 25), 
                     breaks = seq(0, 25, by = 5), 
                     expand = c(0, 0)) +
  theme_minimal()+
  labs(x = "RMSE")+
  theme(axis.title.y = element_blank(),
        plot.margin = margin(t = 10, r = 30, b = 10, l = 12, unit = "pt")) 

# Plot 4: Bray-Curtis Similarity (horizontal) avec errorbars
p3 <- ggplot(metrics_summary, aes(x = Mean_Distance_Similarity, 
                                  y = Classifier, 
                                  fill = Classifier)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(xmin = Mean_Distance_Similarity, 
                    xmax = Mean_Distance_Similarity + Distance_Similarity_se),
                width=0.3, color = "black", linewidth = 0.6) +
  geom_text(aes(label = round(Mean_Distance_Similarity, 3)), 
            hjust =-0.6, vjust = 0.5, fontface = "bold", size = 2.5) +
  scale_fill_manual(values = c("Minimap2" ="#3A9AB2",
                               "Emu"= "#A5C2A3",
                               "Bracken" = "#E3B710",
                               "Qiime2" = "#E79805",
                               "Kraken" = "#5B1A18"),
                    guide = "none") +
  scale_x_continuous(limits = c(0, 1.05), 
                     breaks = seq(0, 1, by = 0.2), 
                     expand = c(0, 0)) +
  theme_minimal()+
  labs(x = "Similarity (1 - Bray)")+
  theme(axis.title.y = element_blank(),
        plot.margin = margin(t = 10, r = 30, b = 10, l = 12, unit = "pt")) 

# Afficher tous les plots ensemble avec ggarrange
library(ggpubr)

# Combiner les 3 mocks avec légende commune
mocks <- ggarrange(mock_a, mock_b, mock_c,
                   nrow = 3,
                   labels = c("A", "B", "C"),
                   hjust = -0.6,
                   vjust = 2,
                   common.legend = TRUE,
                   legend = "bottom",
                   font.label = list(size = 12, face = "bold"))

# Combiner les 3 metrics
metrics <- ggarrange(p1, p4, p3,
                     nrow = 3,
                     labels = c("D", "E", "F"),
                     hjust = 0.4,
                     vjust = 2.2,
                     font.label = list(size = 12, face = "bold"))

# Combiner les deux groupes côte à côte
ggarrange(mocks, metrics,
          ncol = 2,
          widths = c(1, 0.7),
          heights = c(1, 0.7))

ggsave(file="Figures/mock_final.pdf", width=7.1, height=9.5, units="in", dpi=1200)

