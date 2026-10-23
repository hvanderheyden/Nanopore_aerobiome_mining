

library("phyloseq")
library("tidyverse")
library("tibble")

#######################################################################
############## build the phyloseq object for the runs #################
#######################################################################


## load the biom_table, taxonomy and metadata ####
Apogee_biom<- read.csv("Data/runs/Runs_otu.Csv", 
                       header=TRUE, 
                       sep=";")
head(Apogee_biom)
str(Apogee_biom)

# Apogee_biom<-Apogee_biom %>% select(-363)
# head(Apogee_biom)

Apogee_taxo<- read.csv("Data/runs/Runs_taxo.csv", 
                       header=TRUE, 
                       sep=";")
head(Apogee_taxo)

Apogee_meta <- read.csv("Data/runs/Meta_full.tsv", 
                        header=TRUE, 
                        sep="\t")
head(Apogee_meta)
str(Apogee_meta)

Apogee_meta$Year<-as.integer(Apogee_meta$Year)
str(Apogee_meta)

################################################################################
##################### adding weather data to the metadata table ################

meteo <- read.csv ("Data/Meteo/Meteo.csv")
str(meteo)

Apogee_meta <- left_join(
  Apogee_meta,
  meteo,
  by = c("Site", "DOY")
)

head(Apogee_meta)
str(Apogee_meta)

####################### Cleaning the metadata table ############################

Apogee_meta <- Apogee_meta %>% 
  select(-c(ID_S, 
            run_name, 
            flowcell,
            Field,
            barcode, 
            Site_A,
            Name,
            X,
            Year,
            T_Max,
            T_Min,
            Day_length,
            SWE,
            Date))

head(Apogee_meta)

################################################################################

Apogee_meta <- Apogee_meta %>% 
  tibble::column_to_rownames("Sample_ID")

library("phyloseq")
samplesS = sample_data(Apogee_meta)

# define the row names from the otu column ####

Apogee_biom <- Apogee_biom %>%
  tibble::column_to_rownames("otu") 

Apogee_taxo <- Apogee_taxo %>%
  tibble::column_to_rownames("otu") 

# Transform into matrixes otu and tax tables (sample table can be left as data frame) ####

Apogee_biom <- as.matrix(Apogee_biom)
Apogee_taxo <- as.matrix(Apogee_taxo)

class(Apogee_biom)
class(Apogee_taxo)

# convert to phyloseq objects ####

Apogee_OTU = otu_table(Apogee_biom, taxa_are_rows = TRUE)
Apogee_TAX = phyloseq::tax_table(Apogee_taxo)


Apogee_run_raw <- phyloseq(Apogee_OTU, 
                           Apogee_TAX, 
                           samplesS)
Apogee_run_raw

sample_variables(Apogee_run_raw)

#create and merge 
library("ape")
random_tree = rtree(ntaxa(Apogee_run_raw), 
                    rooted=TRUE, 
                    tip.label=taxa_names(Apogee_run_raw))

Apogee_run_raw <- phyloseq(Apogee_OTU, 
                           Apogee_TAX, 
                           samplesS, 
                           random_tree)
Apogee_run_raw

saveRDS(Apogee_run_raw, file= "R_objects/Runs_raw.rds")



