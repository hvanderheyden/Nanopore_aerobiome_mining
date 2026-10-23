

library(tidyverse)
library(phyloseq)
library(Hmisc)
library(ComplexHeatmap)
library(grid)
library(circlize)
library(stringr)
library(ggplot2)
library(forcats)

top_n <- 25    # <--- Set here the number of most abundant OTUs you want to keep

# --- Charger les données ---
Runs_PRF <- readRDS("R_objects/RunsR_R.rds")

Runs_PRF <- subset_taxa(Runs_PRF, Species!="unclassified")

Runs_PRF_melt <- psmelt(Runs_PRF) %>%
  mutate(Species_clean = str_replace(Species, " cluster$", ""))

FunGuildData <- read.csv("Data/runs/FunGuildData.csv", sep=";") %>%
  mutate(Species = taxon) %>%
  select(-c("taxon", "Genus_2", "Species_2", 
            "confidenceRanking",  "growthForm", 
            "citationSource", "host", "substrate", 
            "Function"))

Runs_PRF_G_melt <- Runs_PRF_melt %>%
  left_join(FunGuildData, by = c("Species_clean" = "Species")) %>%
  select(-Species_clean)

head(Runs_PRF_G_melt)

Runs_PRF_G_melt_patho <- Runs_PRF_G_melt %>%
filter(Guild_1 == "Plant Pathogen")

str(Runs_PRF_G_melt_patho)

Runs_PRF_G_melt_patho_corr <- Runs_PRF_G_melt %>%
  select(-c("OTU",  "Mean_BLB", "Mean_SLB", "DM", "Sample", "Site", "LON", "LAT", "date", 
            "BBCH_grouped", "BBCH_bin", "Month", "Week", "DOY", 
            "SampleType", "P_destructor", "B_squamosa", 
            "Kingdom", "Phylum", "Class", "OM", "Thrips",
            "Order", "Family", "Genus",
            "X.", "Guild_1", "simpleFunct", "trophic_mode"))

str(Runs_PRF_G_melt_patho_corr)

# --- Sélection des top N OTUs pathogènes les plus abondantes ---
top_otus <- Runs_PRF_G_melt_patho_corr %>%
  group_by(Species) %>%
  summarise(total_abd = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(total_abd)) %>%
  slice_head(n = top_n) %>%
  pull(Species)

Runs_PRF_G_melt_patho_corr_top <- Runs_PRF_G_melt_patho_corr %>%
  filter(Species %in% top_otus)

# --- Agrégation par variables environnementales et par espèce ---
env_vars <- c("BBCH", "YEAR",
              "Precipitation", "Solar_radiation", "Vapor_pressure", "Wind_speed",
              "Wind_direction", "Mean_temp", "Degree_days", "Mean_temp_3",
              "Mean_temp_5", "Mean_temp_7", "Cum_rain_3", "Cum_rain_5",
              "Cum_rain_7", "Relative_humidity")

# On agrège par espèce et variables environnementales (moyenne d'abondance)
data_corr <- Runs_PRF_G_melt_patho_corr_top %>%
  group_by(Species, across(all_of(env_vars))) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

# On pivote pour avoir une matrice espèce x variables environnementales
mat_for_corr <- data_corr %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0)

# --- Séparer les colonnes OTUs et environnementales ---
otu_cols <- top_otus
env_cols <- env_vars

# Vérifie que les colonnes existent bien
otu_cols <- otu_cols[otu_cols %in% colnames(mat_for_corr)]
env_cols <- env_cols[env_cols %in% colnames(mat_for_corr)]

# Matrices pour corrélation
mat_otu <- as.matrix(mat_for_corr[, otu_cols])
mat_env <- as.matrix(mat_for_corr[, env_cols])

# --- Calcul des corrélations (Spearman) ---
cor_mat <- matrix(NA, nrow = length(otu_cols), ncol = length(env_cols),
                  dimnames = list(otu_cols, env_cols))
p_mat <- matrix(NA, nrow = length(otu_cols), ncol = length(env_cols),
                dimnames = list(otu_cols, env_cols))

for (i in seq_along(otu_cols)) {
  for (j in seq_along(env_cols)) {
    test <- cor.test(mat_otu[, i], mat_env[, j], method = "spearman")
    cor_mat[i, j] <- test$estimate
    p_mat[i, j] <- test$p.value
  }
}

# --- Ajustement des p-values et étoiles ---
p_adj <- matrix(p.adjust(as.vector(p_mat), method = "none"),
                nrow = nrow(p_mat), ncol = ncol(p_mat),
                dimnames = dimnames(p_mat))
stars_mat <- matrix("", nrow=nrow(p_adj), ncol=ncol(p_adj), dimnames=dimnames(p_adj))
stars_mat[p_adj < 0.05 & p_adj >= 0.01]  <- "*"
stars_mat[p_adj < 0.01 & p_adj >= 0.001] <- "**"
stars_mat[p_adj < 0.001]                 <- "***"

# --- Heatmap ---
col_fun <- colorRamp2(c(-1, -0.25, 0, 0.25, 1), c("darkblue", "lightblue",  "#D9E9EBFF", "orange", "red"))

set.seed(1065)
ht <- Heatmap(
  cor_mat,
  name = "Spearman\nrho",
  col = col_fun,
  row_names_gp = grid::gpar(fontface = "italic", fontsize = 8),
  column_names_gp = grid::gpar(fontface = "plain", fontsize = 8),
  row_labels = rownames(cor_mat),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "grey", fill = NA, lwd = 1))
    grid.text(stars_mat[i, j], x, y, gp = gpar(fontsize = 10, fontface = "bold", col="black"))
  },
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  row_split = 7,
  column_split = 4,
  row_title = NULL,
  column_title = NULL
); ht

ht_grob <- grid::grid.grabExpr(
  ComplexHeatmap::draw(ht),
  width = 6,
  height = 4
)

ggplot2::ggsave(
  filename = "Figures/corrplot.pdf",
  plot = ht_grob,
  width = 8.6,
  height = 5.6,
  units = "in"
)

################################################################################

library(ggplot2)
library(forcats)
library(tidyr)
library(dplyr)
library(stringr)

# Filtrer les pathogènes de plante et sélectionner les top N OTUs
Runs_PRF_G_melt_patho_top <- Runs_PRF_G_melt %>%
  group_by(Species) %>%
  summarise(total_abd = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(total_abd)) %>%
  slice_head(n = top_n) %>%
  pull(Species)


# Création de otu_long avec Site inclus
otu_long <- Runs_PRF_G_melt %>%
  filter(Species %in% Runs_PRF_G_melt_patho_top) %>%
  mutate(
    Week = as.integer(Week),
    Species_clean = str_remove(Species, " cluster$")
  ) %>%
  group_by(YEAR, Site, Week, Species) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(YEAR, Site, Week) %>%
  mutate(Abundance_pct = 100 * Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    plot_label = ifelse(Abundance_pct >= 0.01, round(Abundance_pct, 0.1), ""),
    Abundance_pct_plot = ifelse(Abundance_pct == 0, NA, Abundance_pct)
  )

# Ordre des espèces par abondance moyenne
mean_order <- otu_long %>%
  group_by(Species) %>%
  summarise(mean_abd = mean(Abundance_pct, na.rm = TRUE)) %>%
  arrange(desc(mean_abd)) %>%
  pull(Species)

otu_long <- otu_long %>%
  mutate(DisplayLabel = factor(Species, levels = mean_order))


# Séparer l'année et le site pour le faceting
otu_long <- otu_long %>%
  mutate(Year = as.factor(str_extract(Site, "2021|2022")),
         Site_simple = str_remove(Site, "_2021|_2022"))

heat_map <- ggplot(otu_long, aes(x = factor(Week),
                     y = fct_rev(DisplayLabel),
                     fill = Abundance_pct_plot)) +
  geom_tile(color = "grey80", size = 0.7, width = 1, height = 1) +
  geom_text(aes(label = plot_label),
            size = 2, color = "black", family = "sans") +
  scale_fill_gradientn(
    colors = c("blue", "lightblue", "#c6ecff", "#D9E9EBFF", "#fee391", "orange", "red"),
    trans = "log10",
    values = scales::rescale(c(0, 10, 20, 30, 50, 75, 100)),
    limits = c(0.1, max(otu_long$Abundance_pct, na.rm = TRUE)),
    name = "% Read\nAbundance",
    na.value = "grey") +
  facet_grid(Year ~ Site_simple, switch = "y") +
  theme_minimal() +
  theme(legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 1),
    axis.text.y = element_text(face = "italic", size = 7),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    strip.placement = "outside") +
  labs(x = "Week", y = NULL)


save_plot(
  filename = "Figures/heat_map.pdf",  # Output file name
  plot = heat_map,                  # The ggplot object
  base_height = 6.5,           # Height in inches
  base_width = 11.3             # Width in inches
)
