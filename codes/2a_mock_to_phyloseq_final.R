

library("phyloseq")
library("tidyverse")
library("tibble")
#######################################################################
############## build the phyloseq object for the mocks ################
#######################################################################

#######################################
########## Bracken ####################
#######################################

## load the biom_table and taxonomy ####
biom_bracken<- read.csv("Data/new_mock/bracken_otu_raw.biom", 
                        header=TRUE, sep=",")
head(biom_bracken)

taxo_bracken<- read.csv("Data/new_mock/bracken_taxo_raw.tsv", 
                        header=TRUE, sep="\t")
head(taxo_bracken)

# define the row names from the otu column ####

biom_bracken <- biom_bracken %>%
  tibble::column_to_rownames("otu") 

taxo_bracken <- taxo_bracken %>%
  tibble::column_to_rownames("otu") 

# Transform into matrices otu and tax tables ####

biom_bracken <- as.matrix(biom_bracken)
taxo_bracken <- as.matrix(taxo_bracken)

class(taxo_bracken)
class(biom_bracken)

# convert to phyloseq objects ####
OTU_bracken = otu_table(biom_bracken, taxa_are_rows = TRUE)
TAX_bracken = phyloseq::tax_table(taxo_bracken)

phylo_bracken <- phyloseq(OTU_bracken, TAX_bracken)
phylo_bracken

#######################################
########## Kracken ####################
#######################################


## load the biom_table and taxonomy ####
biom_kracken<- read.csv("Data/new_mock/kraken_otu_raw.biom", 
                        header=TRUE, sep=",")
head(biom_kracken)

taxo_kracken<- read.csv("Data/new_mock/kraken_taxo_raw.tsv", 
                        header=TRUE, sep="\t")
head(taxo_kracken)

# define the row names from the otu column ####

biom_kracken <- biom_kracken %>%
  tibble::column_to_rownames("otu") 

taxo_kracken <- taxo_kracken %>%
  tibble::column_to_rownames("otu") 

# Transform into matrices otu and tax tables ####

biom_kracken <- as.matrix(biom_kracken)
taxo_kracken <- as.matrix(taxo_kracken)

class(taxo_kracken)
class(biom_kracken)

# convert to phyloseq objects ####
OTU_kracken = otu_table(biom_kracken, taxa_are_rows = TRUE)
TAX_kracken = phyloseq::tax_table(taxo_kracken)


phylo_kracken <- phyloseq(OTU_kracken, TAX_kracken)
phylo_kracken


###############################
## Qiime2 ####
###############################

## load the biom_table and taxonomy ####
biom_qiime2<- read.csv("Data/new_mock/qiime2_otu.tsv", 
                       header=TRUE, sep="\t")
head(biom_qiime2)

taxo_qiime2<- read.csv("Data/new_mock/qiime2_taxo.tsv", 
                       header=TRUE, sep="\t")
head(taxo_qiime2)

# define the row names from the otu column ####

biom_qiime2 <- biom_qiime2 %>%
  tibble::column_to_rownames("otu") 

taxo_qiime2 <- taxo_qiime2 %>%
  tibble::column_to_rownames("otu") 

# Transform into matrixes otu and tax tables ####

biom_qiime2 <- as.matrix(biom_qiime2)
taxo_qiime2 <- as.matrix(taxo_qiime2)

class(biom_qiime2)
class(taxo_qiime2)

# convert to phyloseq objects ####
OTU_qiime2 = otu_table(biom_qiime2, taxa_are_rows = TRUE)
TAX_qiime2 = phyloseq::tax_table(taxo_qiime2)


phylo_qiime2 <- phyloseq(OTU_qiime2, TAX_qiime2)
phylo_qiime2
#######################################
########## APOGEE #####################
#######################################


## load the biom_table and taxonomy ####
biom_minimap2<- read.csv("Data/new_mock/APOGEE_otu_mock.csv", 
                         header=TRUE, sep=";")
head(biom_minimap2)

taxo_minimap2<- read.csv("Data/new_mock/APOGEE_taxo_mock.csv", 
                         header=TRUE, sep=";")
head(taxo_minimap2)

# define the row names from the otu column ####

biom_minimap2 <- biom_minimap2 %>%
  tibble::column_to_rownames("otu") 

taxo_minimap2 <- taxo_minimap2 %>%
  tibble::column_to_rownames("otu") 

# Transform into matrices otu and tax tables ####

biom_minimap2 <- as.matrix(biom_minimap2)
taxo_minimap2 <- as.matrix(taxo_minimap2)

class(taxo_minimap2)
class(biom_minimap2)

# convert to phyloseq objects ####
OTU_minimap2 = otu_table(biom_minimap2, taxa_are_rows = TRUE)
TAX_minimap2 = phyloseq::tax_table(taxo_minimap2)


phylo_minimap2 <- phyloseq(OTU_minimap2, TAX_minimap2)
phylo_minimap2

#######################################
########## Emu ########################
#######################################


## load the biom_table and taxonomy ####
biom_emu<- read.csv("Data/new_mock/Emu_otu.csv", 
                    header=TRUE, sep=";")
head(biom_emu)

taxo_emu<- read.csv("Data/new_mock/Emu_taxo.csv", 
                    header=TRUE, sep=";")
head(taxo_emu)


# define the row names from the otu column ####

biom_emu <- biom_emu %>%
  tibble::column_to_rownames("otu") 

taxo_emu <- taxo_emu %>%
  tibble::column_to_rownames("otu") 

# Transform into matrices otu and tax tables ####

biom_emu <- as.matrix(biom_emu)
taxo_emu <- as.matrix(taxo_emu)

class(taxo_emu)
class(biom_emu)

# convert to phyloseq objects ####
OTU_emu = phyloseq::otu_table(biom_emu, taxa_are_rows = TRUE)
TAX_emu = phyloseq::tax_table(taxo_emu)

taxa_names(TAX_emu)

phylo_emu <- phyloseq(OTU_emu, TAX_emu)
phylo_emu

###################################
###########   Theoretical #########
###################################

## load the biom_table and taxonomy ####
biom_theo<- read.csv("Data/new_mock/theor_otu_table.csv", 
                     header=TRUE, sep=";")
head(biom_theo)

taxo_theo<- read.csv("Data/new_mock/theor_taxo.csv", 
                     header=TRUE, sep=";")
head(taxo_theo)

# define the row names from the otu column ####

biom_theo <- biom_theo %>%
  tibble::column_to_rownames("otu") 

taxo_theo <- taxo_theo %>%
  tibble::column_to_rownames("otu") 

# Transform into matrices otu and tax tables ####

biom_theo <- as.matrix(biom_theo)
taxo_theo <- as.matrix(taxo_theo)

class(taxo_theo)
class(biom_theo)

# convert to phyloseq objects ####
OTU_theo = otu_table(biom_theo, taxa_are_rows = TRUE)
TAX_theo = phyloseq::tax_table(taxo_theo)


phylo_theo <- phyloseq(OTU_theo, TAX_theo)
phylo_theo

###################################
# Merge the 5 PS objects  #
###################################

# extract the tax tables from PS objects 1 to 6
tax1  = phyloseq::tax_table(phylo_bracken)
tax2 = phyloseq::tax_table(phylo_kracken)
tax3 = phyloseq::tax_table(phylo_qiime2)
tax4  = phyloseq::tax_table(phylo_minimap2)
tax5  = phyloseq::tax_table(phylo_theo)
tax6 = phyloseq::tax_table(phylo_emu)

# merge the tax tables 
tax <- merge_phyloseq(tax1, tax2, tax3, tax4, tax5, tax6)
str(tax)

head(tax)

tax <- as.matrix(tax)

tax = phyloseq::tax_table(tax)

# extract the otu tables from PS objects 1 to 6
otu1 = otu_table(phylo_bracken)
otu2 = otu_table(phylo_kracken)
otu3 = otu_table(phylo_qiime2)
otu4  = otu_table(phylo_minimap2)
otu5  = otu_table(phylo_theo)
otu6 = otu_table(phylo_emu)

# merge the otu tables 
otu <- merge_phyloseq(otu1, otu2, otu3, otu4, otu5, otu6)

sample_names(phylo_bracken)
sample_names(phylo_kracken)
sample_names(phylo_qiime2)
rank_names(phylo_qiime2)
sample_names(phylo_minimap2)
sample_names(phylo_emu)

# import the metadata file 
meta <- read.csv("Data/new_mock/meta.csv", 
                 header=TRUE, sep=";")
head(meta)

meta <- meta %>% 
  tibble::column_to_rownames("Sample_ID")

samples = sample_data(meta)

# merge the otu, tax and meta objects into a final PS object: 
mock_raw <- merge_phyloseq(otu, tax, samples)

mock_raw

# inspect your  phyloseq object: 
sample_names(mock_raw)
rank_names(mock_raw)
sample_variables(mock_raw)

get_taxa(mock_raw)

saveRDS(mock_raw, file= "R_objects/mock_raw.rds")
