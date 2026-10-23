
library("phyloseq")
library("microbiomeutilities")

###############################################################
########## preprocessing for the mock experiments #############
###############################################################

# Import the raw phyloseq mock_raw object 

mock_raw <- readRDS("R_objects/mock_raw.rds")

library("tidyverse")

###############################################################
# Because APOGEE perform species clustering we apply the species 
# clustering for the other pipelines here 
###############################################################

# split APOGEE from the other data sets 
mock_raw_APOGEE <- subset_samples(mock_raw, Description == "Minimap2")
mock_raw_others <- subset_samples(mock_raw, Description != "Minimap2")

# Prune taxa to keep only present in each group 
mock_raw_APOGEE <- prune_taxa(taxa_sums(mock_raw_APOGEE) > 0, mock_raw_APOGEE)
mock_raw_others <- prune_taxa(taxa_sums(mock_raw_others) > 0, mock_raw_others)

###############################################################
# Clean taxa for aPOGEE 
###############################################################

taxM_APOGEE <- data.frame(phyloseq::tax_table(mock_raw_APOGEE))
tax.cleanM_APOGEE <- data.frame(row.names = row.names(taxM_APOGEE),
                        Kingdom = str_replace(taxM_APOGEE[,1], "k__",""),
                        Phylum = str_replace(taxM_APOGEE[,2], "p__",""),
                        Class = str_replace(taxM_APOGEE[,3], "c__",""),
                        Order = str_replace(taxM_APOGEE[,4], "o__",""),
                        Family = str_replace(taxM_APOGEE[,5], "f__",""),
                        Genus = str_replace(taxM_APOGEE[,6], "g__",""),
                        Species = str_replace(taxM_APOGEE[,7], "s__",""),
                        stringsAsFactors = FALSE)
for (i in 1:7){ tax.cleanM_APOGEE[,i] <- as.character(tax.cleanM_APOGEE[,i])}

tax.cleanM_APOGEE <- tax.cleanM_APOGEE %>% 
  mutate(Species=str_replace_all(Species, "_", " "))

phyloseq::tax_table(mock_raw_APOGEE) <- as.matrix(tax.cleanM_APOGEE)

###############################################################
# Clean taxa and 
###############################################################

taxM <- data.frame(phyloseq::tax_table(mock_raw_others))
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
mutate(Species=str_replace_all(Species, "_", " "))

####################################
###### Phytophthora species#########
####################################

# Actualize taxonomy for certain accessions 

tax.cleanM <- tax.cleanM %>% 
  #####################################
########## OOMYCETES ################
#####################################

# cluster 1
  mutate(Species=str_replace(Species, "Achlya americana", "Achlya americana cluster"))%>%
  mutate(Species=str_replace(Species, "Achlya aquatica", "Achlya americana cluster"))%>%
  mutate(Species=str_replace(Species, "Achlya dubia", "Achlya americana cluster"))%>%
  mutate(Species=str_replace(Species, "Achlya glomerata", "Achlya americana cluster"))%>%
  
  # cluster 2
  mutate(Species=str_replace(Species, "Achlya bisexualis", "Achlya bisexualis cluster"))%>%
  mutate(Species=str_replace(Species, "Achlya heterosexualis", "Achlya bisexualis cluster"))%>%
  
  # cluster 3
  mutate(Species=str_replace(Species, "Achlya caroliniana", "Achlya caroliniana cluster"))%>%
  mutate(Species=str_replace(Species, "Achlya conspicua", "Achlya caroliniana cluster"))%>%
  mutate(Species=str_replace(Species, "Achlya flagellata", "Achlya caroliniana cluster"))%>%
  
  # cluster 4
  mutate(Species=str_replace(Species, "Achlya colorata", "Achlya colorata cluster"))%>%
  mutate(Species=str_replace(Species, "Achlya racemosa", "Achlya colorata cluster"))%>%
  mutate(Species=str_replace(Species, "Achlya sparrowii", "Achlya colorata cluster"))%>%
  
  # cluster 16
  mutate(Species=str_replace(Species, "Bremia centaureae", "Bremia centaureae cluster"))%>%
  mutate(Species=str_replace(Species, "Bremia tulasnei", "Bremia centaureae cluster"))%>%
  
  # cluster 19
  mutate(Species=str_replace(Species, "Bremia lactucae", "Bremia lactucae cluster"))%>%
  mutate(Species=str_replace(Species, "Bremia elliptica", "Bremia lactucae cluster"))%>%
  
  # cluster 27
  mutate(Species=str_replace(Species, "Brevilegnia gracilis", "Brevilegnia gracilis cluster"))%>%
  mutate(Species=str_replace(Species, "Globisporangium ultimum", "Brevilegnia gracilis cluster"))%>%
  
  # cluster 35
  mutate(Species=str_replace(Species, "Globisporangium sylvaticum", "Globisporangium sylvaticum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium terrestris", "Globisporangium sylvaticum cluster"))%>%
  
  # cluster 49
  mutate(Species=str_replace(Species, "Peronospora alchemillae", "Peronospora alchemillae cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora potentillae-anserinae", "Peronospora alchemillae cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora potentillae-reptantis", "Peronospora alchemillae cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora sanguisorbae", "Peronospora alchemillae cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora sparsa", "Peronospora alchemillae cluster"))%>%
  
  # cluster 50
  mutate(Species=str_replace(Species, "Peronospora alpicola", "Peronospora alpicola cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora pulveracea", "Peronospora alpicola cluster"))%>%
  
  # cluster 52
  mutate(Species=str_replace(Species, "Peronospora aparines", "Peronospora aparines cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora calotheca", "Peronospora aparines cluster"))%>%
  
  # cluster 54
  mutate(Species=str_replace(Species, "Peronospora arborescens", "Peronospora arborescens cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora argemones", "Peronospora arborescens cluster"))%>%
  
  # cluster 56
  mutate(Species=str_replace(Species, "Peronospora atriplicis-hastatae", "Peronospora atriplicis-hastatae cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora litoralis", "Peronospora atriplicis-hastatae cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora minor", "Peronospora atriplicis-hastatae cluster"))%>%
  
  # cluster 57
  mutate(Species=str_replace(Species, "Peronospora belbahrii", "Peronospora belbahrii cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora glechomae", "Peronospora belbahrii cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora salviae-officinalis", "Peronospora belbahrii cluster"))%>%
  
  # cluster 59
  mutate(Species=str_replace(Species, "Peronospora bulbocapni", "Peronospora bulbocapni cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora chrysosplenii", "Peronospora bulbocapni cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora eranthidis", "Peronospora bulbocapni cluster"))%>%
  
  # cluster 66
  mutate(Species=str_replace(Species, "Peronospora cyparissiae", "Peronospora cyparissiae cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora esulae", "Peronospora cyparissiae cluster"))%>%
  
  # cluster 68
  mutate(Species=str_replace(Species, "Peronospora effusa", "Peronospora effusa cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora farinosa", "Peronospora effusa cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora lepigoni", "Peronospora effusa cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora obovata", "Peronospora effusa cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora polycarpi", "Peronospora effusa cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora rumicis", "Peronospora effusa cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora schachtii", "Peronospora effusa cluster"))%>%
  
  # cluster 83
  mutate(Species=str_replace(Species, "Peronospora pseudostellariae", "Peronospora pseudostellariae cluster"))%>%
  mutate(Species=str_replace(Species, "Peronospora stellariae-uliginosae", "Peronospora pseudostellariae cluster"))%>%
  
  # cluster 99
  mutate(Species=str_replace(Species, "Phytophthora abietivora", "Phytophthora abietivora cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora europaea", "Phytophthora abietivora cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora flexuosa", "Phytophthora abietivora cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora tyrrhenica", "Phytophthora abietivora cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora uliginosa", "Phytophthora abietivora cluster"))%>%
  
  # cluster 101
  mutate(Species=str_replace(Species, "Phytophthora acerina", "Phytophthora acerina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora caryae", "Phytophthora acerina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora citricola", "Phytophthora acerina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora inflata", "Phytophthora acerina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora pachypleura", "Phytophthora acerina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora pini", "Phytophthora acerina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora plurivora", "Phytophthora acerina cluster"))%>%
  
  # cluster 103
  mutate(Species=str_replace(Species, "Phytophthora agathidicida", "Phytophthora agathidicida cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora castaneae", "Phytophthora agathidicida cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora cocois", "Phytophthora agathidicida cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora katsurae", "Phytophthora agathidicida cluster"))%>%
  
  # cluster 104
  mutate(Species=str_replace(Species, "Phytophthora aleatoria", "Phytophthora cactorum cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora alpina", "Phytophthora cactorum cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora cactorum", "Phytophthora cactorum cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora hedraiandra", "Phytophthora cactorum cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora idaei", "Phytophthora cactorum cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora pseudotsugae", "Phytophthora cactorum cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora x serendipita", "Phytophthora cactorum cluster"))%>%
  
  # cluster 105
  mutate(Species=str_replace(Species, "Phytophthora alni", "Phytophthora alni cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora attenuata", "Phytophthora alni cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora cambivora", "Phytophthora alni cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora formosa", "Phytophthora alni cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora intricata", "Phytophthora alni cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora uniformis", "Phytophthora alni cluster"))%>%
  
  # cluster 106
  mutate(Species=str_replace(Species, "Phytophthora alticola", "Phytophthora alticola cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora arenaria", "Phytophthora alticola cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora boodjera", "Phytophthora alticola cluster"))%>%
  
  # cluster 107
  mutate(Species=str_replace(Species, "Phytophthora amaranthi", "Phytophthora amaranthi cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora capsici", "Phytophthora amaranthi cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora gloveri", "Phytophthora amaranthi cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora mexicana", "Phytophthora amaranthi cluster"))%>%
  
  # cluster 109
  mutate(Species=str_replace(Species, "Phytophthora andina", "Phytophthora andina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora infestans", "Phytophthora andina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora ipomoeae", "Phytophthora andina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora mirabilis", "Phytophthora andina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora phaseoli", "Phytophthora andina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora urerae", "Phytophthora andina cluster"))%>%
  
  # cluster 112
  mutate(Species=str_replace(Species, "Phytophthora arecae", "Phytophthora arecae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora heterospora", "Phytophthora arecae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora palmivora", "Phytophthora arecae cluster"))%>%
  
  # cluster 113
  mutate(Species=str_replace(Species, "Phytophthora asiatica", "Phytophthora asiatica cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora niederhauseri", "Phytophthora asiatica cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora pisi", "Phytophthora asiatica cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora sojae", "Phytophthora asiatica cluster"))%>%
  
  # cluster 116
  mutate(Species=str_replace(Species, "Phytophthora aysenensis", "Phytophthora aysenensis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora mengei", "Phytophthora aysenensis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora siskiyouensis", "Phytophthora aysenensis cluster"))%>%
  
  # cluster 117
  mutate(Species=str_replace(Species, "Phytophthora balyanboodja", "Phytophthora balyanboodja cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora chesapeakensis", "Phytophthora balyanboodja cluster"))%>%
  
  # cluster 122
  mutate(Species=str_replace(Species, "Phytophthora borealis", "Phytophthora borealis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora chlamydospora", "Phytophthora borealis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora mississippiae", "Phytophthora borealis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora ornamentata", "Phytophthora borealis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora x stagnum", "Phytophthora borealis cluster"))%>%
  
  # cluster 123
  mutate(Species=str_replace(Species, "Phytophthora botryosa", "Phytophthora botryosa cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora citrophthora", "Phytophthora botryosa cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora himalsilva", "Phytophthora botryosa cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora occultans", "Phytophthora botryosa cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora terminalis", "Phytophthora botryosa cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora x vanyenensis", "Phytophthora botryosa cluster"))%>%
  
  # cluster 124
  mutate(Species=str_replace(Species, "Phytophthora brassicae", "Phytophthora brassicae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora dauci", "Phytophthora brassicae cluster"))%>%
  
  # cluster 125
  mutate(Species=str_replace(Species, "Phytophthora cacuminis", "Phytophthora cacuminis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora captiosa", "Phytophthora cacuminis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora fallax", "Phytophthora cacuminis cluster"))%>%
  
  # cluster 126
  mutate(Species=str_replace(Species, "Phytophthora cajani", "Phytophthora cajani cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora vignae", "Phytophthora cajani cluster"))%>%
  
  # cluster 127
  mutate(Species=str_replace(Species, "Phytophthora capensis", "Phytophthora capensis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora emzansi", "Phytophthora capensis cluster"))%>%
  
  # cluster 128
  mutate(Species=str_replace(Species, "Phytophthora castanetorum", "Phytophthora castanetorum cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora quercina", "Phytophthora castanetorum cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora versiformis", "Phytophthora castanetorum cluster"))%>%
  
  # cluster 130
  mutate(Species=str_replace(Species, "Phytophthora celebensis", "Phytophthora celebensis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora javanensis", "Phytophthora celebensis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora multiglobulosa", "Phytophthora celebensis cluster"))%>%
  
  # cluster 131
  mutate(Species=str_replace(Species, "Phytophthora chilensis", "Phytophthora chilensis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora kernoviae", "Phytophthora chilensis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora pseudochilensis", "Phytophthora chilensis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora pseudokernoviae", "Phytophthora chilensis cluster"))%>%
  
  # cluster 136
  mutate(Species=str_replace(Species, "Phytophthora colocasiae", "Phytophthora colocasiae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora insulinativitatica", "Phytophthora colocasiae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora meadii", "Phytophthora colocasiae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora multibullata", "Phytophthora colocasiae cluster"))%>%
  
  # cluster 137
  mutate(Species=str_replace(Species, "Phytophthora condilina", "Phytophthora condilina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora humicola", "Phytophthora condilina cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora inundata", "Phytophthora condilina cluster"))%>%
  
  # cluster 139
  mutate(Species=str_replace(Species, "Phytophthora cooljarloo", "Phytophthora cooljarloo cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora kwongonina", "Phytophthora cooljarloo cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora oreophila", "Phytophthora cooljarloo cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora pseudorosacearum", "Phytophthora cooljarloo cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora rosacearum", "Phytophthora cooljarloo cluster"))%>%
  
  # cluster 140
  mutate(Species=str_replace(Species, "Phytophthora crassamura", "Phytophthora crassamura cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora megasperma", "Phytophthora crassamura cluster"))%>%
  
  # cluster 141
  mutate(Species=str_replace(Species, "Phytophthora cryptogea", "Phytophthora cryptogea cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora erythroseptica", "Phytophthora cryptogea cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora himalayensis", "Phytophthora cryptogea cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora kelmanii", "Phytophthora cryptogea cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora pseudocryptogea", "Phytophthora cryptogea cluster"))%>%
  
  # cluster 148
  mutate(Species=str_replace(Species, "Phytophthora fragariae", "Phytophthora fragariae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora rubi", "Phytophthora fragariae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora x heterohybrida", "Phytophthora fragariae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora x incrassata", "Phytophthora fragariae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora x multiformis", "Phytophthora fragariae cluster"))%>%
  
  # cluster 149
  mutate(Species=str_replace(Species, "Phytophthora fragariaefolia", "Phytophthora fragariaefolia cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora nagaii", "Phytophthora fragariaefolia cluster"))%>%
  
  # cluster 151
  mutate(Species=str_replace(Species, "Phytophthora gallica", "Phytophthora gallica cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora subarctica", "Phytophthora gallica cluster"))%>%
  
  # cluster 153
  mutate(Species=str_replace(Species, "Phytophthora gibbosa", "Phytophthora gibbosa cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora gregata", "Phytophthora gibbosa cluster"))%>%
  
  # cluster 159
  mutate(Species=str_replace(Species, "Phytophthora ilicis", "Phytophthora ilicis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora nemorosa", "Phytophthora ilicis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora pseudosyringae", "Phytophthora ilicis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora psychrophila", "Phytophthora ilicis cluster"))%>%
  
  # cluster 164
  mutate(Species=str_replace(Species, "Phytophthora lactucae", "Phytophthora lactucae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora pseudolactucae", "Phytophthora lactucae cluster"))%>%
  
  # cluster 165
  mutate(Species=str_replace(Species, "Phytophthora lacustris", "Phytophthora lacustris cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora riparia", "Phytophthora lacustris cluster"))%>%
  
  # cluster 172
  mutate(Species=str_replace(Species, "Phytophthora macrochlamydospora", "Phytophthora macrochlamydospora cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora richardiae", "Phytophthora macrochlamydospora cluster"))%>%
  
  # cluster 178
  mutate(Species=str_replace(Species, "Phytophthora melonis", "Phytophthora melonis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora sinensis", "Phytophthora melonis cluster"))%>%
  
  # cluster 183
  mutate(Species=str_replace(Species, "Phytophthora nicotianae", "Phytophthora nicotianae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora tabaci", "Phytophthora nicotianae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora x pelgrandis", "Phytophthora nicotianae cluster"))%>%
  
  # cluster 186
  mutate(Species=str_replace(Species, "Phytophthora panamensis", "Phytophthora panamensis cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora quercetorum", "Phytophthora panamensis cluster"))%>%
  
  # cluster 193
  mutate(Species=str_replace(Species, "Phytophthora polonica", "Phytophthora polonica cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora pseudopolonica", "Phytophthora polonica cluster"))%>%
  
  # cluster 194
  mutate(Species=str_replace(Species, "Phytophthora porri", "Phytophthora porri cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora primulae", "Phytophthora porri cluster"))%>%
  
  # cluster 200
  mutate(Species=str_replace(Species, "Phytophthora rhizophorae", "Phytophthora rhizophorae cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora virginiana", "Phytophthora rhizophorae cluster"))%>%
  
  # cluster 220
  mutate(Species=str_replace(Species, "Phytopythium boreale", "Phytopythium boreale cluster"))%>%
  mutate(Species=str_replace(Species, "Phytopythium megacarpum", "Phytopythium boreale cluster"))%>%
  
  # cluster 236
  mutate(Species=str_replace(Species, "Pythium acrogynum", "Pythium acrogynum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium hypogynum", "Pythium acrogynum cluster"))%>%
  
  # cluster 238
  mutate(Species=str_replace(Species, "Pythium afertile", "Pythium afertile cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium kashmirense", "Pythium afertile cluster"))%>%
  
  # cluster 239
  mutate(Species=str_replace(Species, "Pythium aff. diclinum", "Pythium aff. diclinum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium aff. dictyosporum", "Pythium aff. diclinum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium coloratum", "Pythium aff. diclinum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium diclinum", "Pythium aff. diclinum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium dissotocum", "Pythium aff. diclinum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium lutarium", "Pythium aff. diclinum cluster"))%>%
  
  # cluster 242
  mutate(Species=str_replace(Species, "Pythium amasculinum", "Pythium amasculinum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium hydnosporum", "Pythium amasculinum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium lycopersicum", "Pythium amasculinum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium oligandrum", "Pythium amasculinum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium ornamentatum", "Pythium amasculinum cluster"))%>%
  
  # cluster 248
  mutate(Species=str_replace(Species, "Pythium aquatile", "Pythium aquatile cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium sukuiense", "Pythium aquatile cluster"))%>%
  
  # cluster 249
  mutate(Species=str_replace(Species, "Pythium aristosporum", "Pythium aristosporum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium arrhenomanes", "Pythium aristosporum cluster"))%>%
  
  # cluster 251
  mutate(Species=str_replace(Species, "Pythium buismaniae", "Pythium buismaniae cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium megalacanthum", "Pythium buismaniae cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium polymastum", "Pythium buismaniae cluster"))%>%
  
  # cluster 254
  mutate(Species=str_replace(Species, "Pythium capillosum", "Pythium capillosum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium flevoense", "Pythium capillosum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium pectinolyticum", "Pythium capillosum cluster"))%>%
  
  # cluster 256
  mutate(Species=str_replace(Species, "Pythium catenulatum", "Pythium catenulatum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium folliculosum", "Pythium catenulatum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium rhizo-oryzae", "Pythium catenulatum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium torulosum", "Pythium catenulatum cluster"))%>%
  
  # cluster 257
  mutate(Species=str_replace(Species, "Pythium chondricola", "Pythium chondricola cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium porphyrae", "Pythium chondricola cluster"))%>%
  
  # cluster 258
  mutate(Species=str_replace(Species, "Pythium conidiophorum", "Pythium conidiophorum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium salpingophorum", "Pythium conidiophorum cluster"))%>%
  
  # cluster 259
  mutate(Species=str_replace(Species, "Pythium contiguanum", "Pythium contiguanum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium dissimile", "Pythium contiguanum cluster"))%>%
  
  # cluster 264
  mutate(Species=str_replace(Species, "Pythium debaryanum", "Pythium debaryanum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium viniferum", "Pythium debaryanum cluster"))%>%
  
  # cluster 266
  mutate(Species=str_replace(Species, "Pythium dimorphum", "Pythium dimorphum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium undulatum", "Pythium dimorphum cluster"))%>%
  
  # cluster 267
  mutate(Species=str_replace(Species, "Pythium echinulatum", "Pythium echinulatum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium erinaceum", "Pythium echinulatum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium ornacarpum", "Pythium echinulatum cluster"))%>%
  
  # cluster 269
  mutate(Species=str_replace(Species, "Pythium glomeratum", "Pythium glomeratum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium heterothallicum", "Pythium glomeratum cluster"))%>%
  
  # cluster 270
  mutate(Species=str_replace(Species, "Pythium graminicola", "Pythium graminicola cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium periilum", "Pythium graminicola cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium tardicrescens", "Pythium graminicola cluster"))%>%
  
  # cluster 278
  mutate(Species=str_replace(Species, "Pythium kunmingense", "Pythium kunmingense cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium spinosum", "Pythium kunmingense cluster"))%>%
  
  # cluster 288
  mutate(Species=str_replace(Species, "Pythium minus", "Pythium minus cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium pleroticum", "Pythium minus cluster"))%>%
  
  # cluster 291
  mutate(Species=str_replace(Species, "Pythium myriotylum", "Pythium myriotylum cluster"))%>%
  mutate(Species=str_replace(Species, "Pythium zingiberis", "Pythium myriotylum cluster"))%>%
  
  # cluster 326
  mutate(Species=str_replace(Species, "Saprolegnia delica", "Saprolegnia delica cluster"))%>%
  mutate(Species=str_replace(Species, "Saprolegnia parasitica", "Saprolegnia delica cluster"))%>%
  
  # cluster 327
  mutate(Species=str_replace(Species, "Saprolegnia diclina", "Saprolegnia diclina cluster"))%>%
  mutate(Species=str_replace(Species, "Saprolegnia ferax", "Saprolegnia diclina cluster"))%>%
  mutate(Species=str_replace(Species, "Saprolegnia lapponica", "Saprolegnia diclina cluster"))%>%
  mutate(Species=str_replace(Species, "Saprolegnia mixta", "Saprolegnia diclina cluster"))%>%
  mutate(Species=str_replace(Species, "Saprolegnia unispora", "Saprolegnia diclina cluster"))%>%
  
  # cluster 337
  mutate(Species=str_replace(Species, "Thraustotheca clavata", "Thraustotheca clavata cluster"))%>%
  mutate(Species=str_replace(Species, "Thraustotheca terrestris", "Thraustotheca clavata cluster"))

#####################################
########## Fungi ################
#####################################

tax.cleanM <- tax.cleanM %>% 
  
  #####################################
###### Sclerotinia species ##########
#####################################

# S. sclerotiorum cluster
mutate(Species=str_replace(Species, "Sclerotinia sclerotiorum", "Sclerotinia sclerotiorum cluster"))%>%
  mutate(Species=str_replace(Species, "Sclerotinia minor", "Sclerotinia sclerotiorum cluster"))%>%
  mutate(Species=str_replace(Species, "Sclerotinia trifoliorum", "Sclerotinia sclerotiorum cluster"))%>%
  
  #####################################
###### Botrytis species #############
#####################################

# B. cinerea cluster
mutate(Species=str_replace(Species, "Botrytis cinerea", "Botrytis cinerea cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis fabae", "Botrytis cinerea cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis pelargonii", "Botrytis cinerea cluster")) %>%
  mutate(Species=str_replace(Species, "Botrytis eucalypti", "Botrytis cinerea cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis pseudocinerea", "Botrytis cinerea cluster")) %>%
  mutate(Species=str_replace(Species, "Botrytis polyphyllae", "Botrytis cinerea cluster")) %>%
  
  # B. calthae cluster 
  mutate(Species=str_replace(Species, "Botrytis calthae", "Botrytis calthae cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis californica", "Botrytis calthae cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis sinoviticola", "Botrytis calthae cluster")) %>%
  
  # B. byssoidea cluster 
  mutate(Species=str_replace(Species, "Botrytis byssoidea", "Botrytis byssoidea cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis globosa", "Botrytis byssoidea cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis sphaerosperma", "Botrytis byssoidea cluster")) %>%
  mutate(Species=str_replace(Species, "Botrytis polyblastis", "Botrytis byssoidea cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis tulipae", "Botrytis byssoidea cluster")) %>%
  mutate(Species=str_replace(Species, "Botrytis gladiolorum", "Botrytis byssoidea cluster")) %>%
  mutate(Species=str_replace(Species, "Botrytis narcissicola", "Botrytis byssoidea cluster")) %>%
  
  # B. croci cluster 
  mutate(Species=str_replace(Species, "Botrytis croci", "Botrytis croci cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis hyacinthi", "Botrytis croci cluster"))%>%
  
  # B. caroliniana cluster 
  mutate(Species=str_replace(Species, "Botrytis caroliniana", "Botrytis caroliniana cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis fabiopsis", "Botrytis caroliniana cluster")) %>%
  
  # B. fragariae cluster 
  mutate(Species=str_replace(Species, "Botrytis fragariae", "Botrytis fragariae cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis prunorum", "Botrytis fragariae cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis convoluta", "Botrytis fragariae cluster")) %>%
  mutate(Species=str_replace(Species, "Botrytis paeoniae", "Botrytis fragariae cluster")) %>%
  
  # B. aclada cluster 
  mutate(Species=str_replace(Species, "Botrytis aclada", "Botrytis aclada cluster"))%>%
  mutate(Species=str_replace(Species, "Botrytis euroamericana", "Botrytis aclada cluster"))%>%
  
  # B. squamosa cluster 
  # mutate(Species=str_replace(Species, "Botrytis squamosa", "Botrytis squamosa cluster"))%>%
  #  mutate(Species=str_replace(Species, "Botrytis elliptica", "Botrytis squamosa cluster"))%>%
  #  mutate(Species=str_replace(Species, "Botrytis deweyae", "Botrytis squamosa cluster"))%>%
  #  mutate(Species=str_replace(Species, "Botrytis ficarium", "Botrytis squamosa cluster"))%>%
  #  mutate(Species=str_replace(Species, "Botrytis ranunculi", "Botrytis squamosa cluster"))%>%
  #  mutate(Species=str_replace(Species, "Botrytis sinoallii", "Botrytis squamosa cluster"))%>%  
  
  ####################################
###### Stemphylium species #########
####################################

# Stemphylium vesicarium cluster
mutate(Species=str_replace(Species, "Stemphylium vesicarium", "Stemphylium vesicarium cluster"))%>%
  mutate(Species=str_replace(Species, "Stemphylium botryosum", "Stemphylium vesicarium cluster"))%>%
  mutate(Species=str_replace(Species, "Stemphylium eturmiunum", "Stemphylium vesicarium cluster"))%>%
  
  # Stemphylium trifolii cluster
  mutate(Species=str_replace(Species, "Stemphylium trifolii", "Stemphylium trifolii cluster"))%>%
  mutate(Species=str_replace(Species, "Stemphylium triglochinicola", "Stemphylium trifolii cluster"))%>%
  
  #####################################
######### Fusarium species ##########
#####################################

# Fusarium acaciae-mearnsii cluster
mutate(Species=str_replace(Species, "Fusarium acaciae-mearnsii", "Fusarium acaciae-mearnsii cluster"))%>%
  mutate(Species=str_replace(Species, "Fusarium falsibabinda", "Fusarium acaciae-mearnsii cluster"))%>%
  
  # Fusarium acutatum cluster
  mutate(Species=str_replace(Species, "Fusarium acutatum", "Fusarium acutatum cluster"))%>%
  mutate(Species=str_replace(Species, "Fusarium coeruleum", "Fusarium acutatum cluster"))%>%
  mutate(Species=str_replace(Species, "Fusarium commune", "Fusarium acutatum cluster"))%>%
  
  # Fusarium bambusarum cluster
  mutate(Species=str_replace(Species, "Fusarium bambusarum", "Fusarium bambusarum cluster"))%>%
  mutate(Species=str_replace(Species, "Fusarium polyphialidicum", "Fusarium bambusarum cluster"))%>%
  
  # Fusarium cicatricum cluster
  mutate(Species=str_replace(Species, "Fusarium cicatricum", "Fusarium cicatricum cluster"))%>%
  mutate(Species=str_replace(Species, "Fusarium zealandicum", "Fusarium cicatricum cluster"))%>%
  
  # Fusarium nygamai cluster
  mutate(Species=str_replace(Species, "Fusarium nygamai", "Fusarium nygamai cluster"))%>%
  mutate(Species=str_replace(Species, "Fusarium xylarioides", "Fusarium nygamai cluster"))%>%
  
  # Fusarium verticillioides
  mutate(Species=str_replace(Species, "Fusarium verticillioides", "Fusarium verticillioides cluster"))%>%
  mutate(Species=str_replace(Species, "Fusarium pseudoanthophilum", "Fusarium verticillioides cluster"))%>%
  
  #####################################
######### Alternaria species ######## 
#####################################

# Section Alternantherae 
mutate(Species=str_replace(Species, "Alternaria alternantherae", "Alternaria sect. alternantherae"))%>%
  mutate(Species=str_replace(Species, "Alternaria celosiicola", "Alternaria sect. alternantherae"))%>%
  mutate(Species=str_replace(Species, "Alternaria gomphrenae", "Alternaria sect. alternantherae"))%>%
  mutate(Species=str_replace(Species, "Alternaria perpunctulata", "Alternaria sect. alternantherae"))%>%
  
  # Section Alternata 
  mutate(Species=str_replace(Species, "Alternaria alternata", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria angustiovoidea", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria arborescens", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria burnsii", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria cerealis", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria citriarbusti", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria citrimaculari", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria colombiana", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria daucifollii", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria destruens", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria dumosa", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria gaisen", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria herbiphorbicola", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria limoniasperae", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria longipes", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria perangusta", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria postmessia", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria tangelonis", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria tenuissima", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria toxicogenica", "Alternaria sect. alternata"))%>%
  mutate(Species=str_replace(Species, "Alternaria turkisafria", "Alternaria sect. alternata"))%>%
  
  # Section Brassicicola 
  mutate(Species=str_replace(Species, "Alternaria brassicicola", "Alternaria sect. brassicicola"))%>%
  mutate(Species=str_replace(Species, "Alternaria conoidea", "Alternaria sect. brassicicola"))%>%
  mutate(Species=str_replace(Species, "Alternaria mimicula", "Alternaria sect. brassicicola"))%>%
  mutate(Species=str_replace(Species, "Alternaria septorioides", "Alternaria sect. brassicicola"))%>%
  mutate(Species=str_replace(Species, "Alternaria solidaccana", "Alternaria sect. brassicicola"))%>%
  
  # Section Chalastospora 
  mutate(Species=str_replace(Species, "Alternaria obclavata", "Alternaria sect. chalastospora"))%>%
  mutate(Species=str_replace(Species, "Alternaria abundans", "Alternaria sect. chalastospora"))%>%
  mutate(Species=str_replace(Species, "Alternaria armoraciae", "Alternaria sect. chalastospora"))%>%
  mutate(Species=str_replace(Species, "Alternaria cetera", "Alternaria sect. chalastospora"))%>%
  mutate(Species=str_replace(Species, "Alternaria malorum", "Alternaria sect. chalastospora"))%>%
  
  # Section Cheiranthus 
  mutate(Species=str_replace(Species, "Alternaria cheiranthi", "Alternaria sect. cheiranthus"))%>%
  mutate(Species=str_replace(Species, "Alternaria indefessa", "Alternaria sect. cheiranthus"))%>%
  
  # Section Crivellia 
  mutate(Species=str_replace(Species, "Alternaria crivellia", "Alternaria sect. crivellia"))%>%
  mutate(Species=str_replace(Species, "Alternaria papavericola", "Alternaria sect. crivellia"))%>%
  mutate(Species=str_replace(Species, "Alternaria penicillata", "Alternaria sect. crivellia"))%>%
  
  # Section Dianthicola 
  mutate(Species=str_replace(Species, "Alternaria dianthicola", "Alternaria sect. dianthicola"))%>%
  mutate(Species=str_replace(Species, "Alternaria elegans", "Alternaria sect. dianthicola"))%>%
  mutate(Species=str_replace(Species, "Alternaria simsimi", "Alternaria sect. dianthicola"))%>%
  
  # Section Embellisia 
  mutate(Species=str_replace(Species, "Alternaria embellisia", "Alternaria sect. embellisia"))%>%
  mutate(Species=str_replace(Species, "Alternaria chlamydosporigena", "Alternaria sect. embellisia"))%>%
  mutate(Species=str_replace(Species, "Alternaria tellustris", "Alternaria sect. embellisia"))%>%
  
  # Section Embellisioides
  mutate(Species=str_replace(Species, "Alternaria embellisioides", "Alternaria sect. embellisioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria botryospora", "Alternaria sect. embellisioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria hyacinthi", "Alternaria sect. embellisioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria lolii", "Alternaria sect. embellisioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria planifunda", "Alternaria sect. embellisioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria proteae", "Alternaria sect. embellisioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria tumida", "Alternaria sect. embellisioides"))%>%
  
  # Section Eureka
  mutate(Species=str_replace(Species, "Alternaria eureka", "Alternaria sect. eureka"))%>%
  mutate(Species=str_replace(Species, "Alternaria anigozanthi", "Alternaria sect. eureka"))%>%
  mutate(Species=str_replace(Species, "Alternaria cumini", "Alternaria sect. eureka"))%>%
  mutate(Species=str_replace(Species, "Alternaria geniostomatis", "Alternaria sect. eureka"))%>%
  mutate(Species=str_replace(Species, "Alternaria leptinellae", "Alternaria sect. eureka"))%>%
  mutate(Species=str_replace(Species, "Alternaria triglochinicola", "Alternaria sect. eureka"))%>%
  
  # Section Gypsophilae
  mutate(Species=str_replace(Species, "Alternaria gypsophilae", "Alternaria sect. gypsophilae"))%>%
  mutate(Species=str_replace(Species, "Alternaria axiaeriisporifera", "Alternaria sect. gypsophilae"))%>%
  mutate(Species=str_replace(Species, "Alternaria ellipsoidea", "Alternaria sect. gypsophilae"))%>%
  mutate(Species=str_replace(Species, "Alternaria juxtiseptata", "Alternaria sect. gypsophilae"))%>%
  mutate(Species=str_replace(Species, "Alternaria nobilis", "Alternaria sect. gypsophilae"))%>%
  mutate(Species=str_replace(Species, "Alternaria saponariae", "Alternaria sect. gypsophilae"))%>%
  mutate(Species=str_replace(Species, "Alternaria vaccariae", "Alternaria sect. gypsophilae"))%>%
  mutate(Species=str_replace(Species, "Alternaria vaccariicola", "Alternaria sect. gypsophilae"))%>%
  
  # Section Infectoriae
  mutate(Species=str_replace(Species, "Alternaria infectoria", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria alternarina", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria arbusti", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria caespitosa", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria californica", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria conjuncta", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria daucicaulis", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria ethzedia", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria frumenti", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria graminicola", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria hordeiaustralica", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria hordeicola", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria humuli", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria incomplexa", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria intercepta", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria merytae", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria metachromatica", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria novae-zelandiae", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria oregonensis", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria slovaca", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria triticimaculans", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria triticina", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria ventricosa", "Alternaria sect. infectoriae"))%>%
  mutate(Species=str_replace(Species, "Alternaria viburni", "Alternaria sect. infectoriae"))%>%
  
  # Section Japonicae
  mutate(Species=str_replace(Species, "Alternaria japonica", "Alternaria sect. japonicae"))%>%
  mutate(Species=str_replace(Species, "Alternaria nepalensis", "Alternaria sect. japonicae"))%>%
  
  # Section Nimbya
  mutate(Species=str_replace(Species, "Alternaria nimbya", "Alternaria sect. nimbya"))%>%
  mutate(Species=str_replace(Species, "Alternaria caricis", "Alternaria sect. nimbya"))%>%
  mutate(Species=str_replace(Species, "Alternaria scirpinfestans", "Alternaria sect. nimbya"))%>%
  mutate(Species=str_replace(Species, "Alternaria scirpivora", "Alternaria sect. nimbya"))%>%
  
  
  # Section Panax
  mutate(Species=str_replace(Species, "Alternaria panax", "Alternaria sect. panax"))%>%
  mutate(Species=str_replace(Species, "Alternaria avenicola", "Alternaria sect. panax"))%>%
  mutate(Species=str_replace(Species, "Alternaria calycipyricola", "Alternaria sect. panax"))%>%
  mutate(Species=str_replace(Species, "Alternaria eryngii", "Alternaria sect. panax"))%>%
  mutate(Species=str_replace(Species, "Alternaria photistica", "Alternaria sect. panax"))%>%
  
  # Section Phragmosporae
  mutate(Species=str_replace(Species, "Alternaria phragmospora", "Alternaria sect. phragmosporae"))%>%
  mutate(Species=str_replace(Species, "Alternaria chlamydospora", "Alternaria sect. phragmosporae"))%>%
  mutate(Species=str_replace(Species, "Alternaria didymospora", "Alternaria sect. phragmosporae"))%>%
  mutate(Species=str_replace(Species, "Alternaria limaciformis", "Alternaria sect. phragmosporae"))%>%
  mutate(Species=str_replace(Species, "Alternaria molesta", "Alternaria sect. phragmosporae"))%>%
  mutate(Species=str_replace(Species, "Alternaria mouchaccae", "Alternaria sect. phragmosporae"))%>%
  
  # Section Porri
  mutate(Species=str_replace(Species, "Alternaria porri", "Alternaria sect. porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria acalyphicola", "Alternaria sect. porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria agerati", "Alternaria sect. porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria agripestis", "Alternaria sect. porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria anagallidis", "Alternaria sect. porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria aragakii", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria argyroxiphii", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria bataticola", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria blumeae", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria calendulae", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria capsici", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria carthami", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria cassiae", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria cichorii", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria cirsinoxia", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria crassa", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria cretica", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria cucumerina", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria cyphomandrae", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria danida", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria dauci", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria dichondrae", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria euphorbiicola", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria grandis", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria hawaiiensis", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria limicola", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria linicola", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria macrospora", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria multirostrata", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria nitrimali", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria passiflorae", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria poonensis", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria protenta", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria pseudorostrata", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria ricini", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria rostellata", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria scorzonerae", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria sesami", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria solani", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria solani-nigri", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria steviae", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria subcylindrica", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria tagetica", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria tomatophila", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria tropica", "Alternaria sect. Porri"))%>%
  mutate(Species=str_replace(Species, "Alternaria zinniae", "Alternaria sect. Porri"))%>%
  
  # Section Pseudoulocladium
  mutate(Species=str_replace(Species, "Alternaria phragmospora", "Alternaria sect. pseudoulocladium"))%>%
  mutate(Species=str_replace(Species, "Alternaria aspera", "Alternaria sect. pseudoulocladium"))%>%
  mutate(Species=str_replace(Species, "Alternaria chartarum", "Alternaria sect. pseudoulocladium"))%>%
  mutate(Species=str_replace(Species, "Alternaria septospora", "Alternaria sect. pseudoulocladium"))%>%
  
  # Section Radicina
  mutate(Species=str_replace(Species, "Alternaria radicina ", "Alternaria sect. radicina"))%>%
  mutate(Species=str_replace(Species, "Alternaria carotiincultae ", "Alternaria sect. radicina"))%>%
  mutate(Species=str_replace(Species, "Alternaria petroselini", "Alternaria sect. radicina"))%>%
  mutate(Species=str_replace(Species, "Alternaria selini", "Alternaria sect. radicina"))%>%
  mutate(Species=str_replace(Species, "Alternaria smyrnii", "Alternaria sect. radicina"))%>%
  
  # Section Sonchi
  mutate(Species=str_replace(Species, "Alternaria sonchi", "Alternaria sect. sonchi"))%>%
  mutate(Species=str_replace(Species, "Alternaria cinerariae", "Alternaria sect. sonchi"))%>%
  
  # Section Teretispora
  mutate(Species=str_replace(Species, "Alternaria leucanthemi", "Alternaria sect. teretispora"))%>%
  
  # Section Ulocladioides
  mutate(Species=str_replace(Species, "Alternaria ulocladioides ", "Alternaria sect. ulocladioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria atra", "Alternaria sect. ulocladioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria brassicae-pekinensis", "Alternaria sect. ulocladioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria cantlous", "Alternaria sect. ulocladioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria consortialis", "Alternaria sect. ulocladioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria cucurbitae", "Alternaria sect. ulocladioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria heterospora", "Alternaria sect. ulocladioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria multiformis", "Alternaria sect. ulocladioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria obovoidea", "Alternaria sect. ulocladioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria subcucurbitae", "Alternaria sect. ulocladioides"))%>%
  mutate(Species=str_replace(Species, "Alternaria terricola", "Alternaria sect. ulocladioides"))%>%
  
  # Section Ulocladium
  mutate(Species=str_replace(Species, "Alternaria ulocladium ", "Alternaria sect. ulocladium"))%>%
  mutate(Species=str_replace(Species, "Alternaria alternariae", "Alternaria sect. ulocladium"))%>%
  mutate(Species=str_replace(Species, "Alternaria botrytis", "Alternaria sect. ulocladium"))%>%
  mutate(Species=str_replace(Species, "Alternaria capsici-annui", "Alternaria sect. ulocladium"))%>%
  mutate(Species=str_replace(Species, "Alternaria oudemansii", "Alternaria sect. ulocladium"))%>%
  
  # Section Undifilum
  mutate(Species=str_replace(Species, "Alternaria bornmuelleri", "Alternaria sect. undifilum"))%>%
  mutate(Species=str_replace(Species, "Alternaria cinerea", "Alternaria sect. undifilum"))%>%
  mutate(Species=str_replace(Species, "Alternaria fulva", "Alternaria sect. undifilum"))%>%
  mutate(Species=str_replace(Species, "Alternaria oxytropis", "Alternaria sect. undifilum"))%>%
  
  #####################################
######### Neosetophoma species ######
#####################################

# Neosetophoma buxi cluster
mutate(Species=str_replace(Species, "Neosetophoma buxi", "Neosetophoma buxi cluster"))%>%
  mutate(Species=str_replace(Species, "Neosetophoma phragmitis", "Neosetophoma buxi cluster"))%>%
  
  # Neosetophoma clematidis cluster
  mutate(Species=str_replace(Species, "Neosetophoma clematidis", "Neosetophoma clematidis cluster"))%>%
  mutate(Species=str_replace(Species, "Neosetophoma lunariae", "Neosetophoma clematidis cluster"))%>%
  
  #####################################
######### Peniophora species ######
#####################################
# Peniophora cinerea cluster
mutate(Species=str_replace(Species, "Peniophora cinerea", "Peniophora cinerea cluster"))%>%
  mutate(Species=str_replace(Species, "Peniophora piceae", "Peniophora cinerea cluster"))%>%
  
  #####################################
######### Pyrenophora species ######
#####################################
# Pyrenophora avenicola cluster
mutate(Species=str_replace(Species, "Pyrenophora avenicola", "Pyrenophora avenicola cluster"))%>%
  mutate(Species=str_replace(Species, "Pyrenophora lolii", "Pyrenophora avenicola cluster"))%>%
  
  # Pyrenophora trichostoma cluster
  mutate(Species=str_replace(Species, "Pyrenophora trichostoma", "Pyrenophora trichostoma cluster"))%>%
  mutate(Species=str_replace(Species, "Pyrenophora tritici-repentis", "Pyrenophora trichostoma cluster"))%>%
  
  #####################################
######### Rhodotorula species ######
#####################################
# Rhodotorula araucariae cluster
mutate(Species=str_replace(Species, "Rhodotorula araucariae", "Rhodotorula araucariae cluster"))%>%
  mutate(Species=str_replace(Species, "Rhodotorula babjevae", "Rhodotorula araucariae cluster"))%>%
  mutate(Species=str_replace(Species, "Rhodotorula graminis", "Rhodotorula araucariae cluster"))%>%
  
  # Rhodotorula dairenensis cluster
  mutate(Species=str_replace(Species, "Rhodotorula dairenensis", "Rhodotorula dairenensis cluster"))%>%
  mutate(Species=str_replace(Species, "Rhodotorula mucilaginosa", "Rhodotorula dairenensis cluster"))%>%
  
  # Rhodotorula diobovata cluster
  mutate(Species=str_replace(Species, "Rhodotorula diobovata", "Rhodotorula diobovata cluster"))%>%
  mutate(Species=str_replace(Species, "Rhodotorula kratochvilovae", "Rhodotorula diobovata cluster"))%>%
  
  #####################################
######### Ganoderma species ######
#####################################

# Ganoderma bambusicola
mutate(Species=str_replace(Species, "Ganoderma bambusicola", "Ganoderma bambusicola cluster"))%>%
  mutate(Species=str_replace(Species, "Ganoderma neojaponicum", "Ganoderma bambusicola cluster"))%>%
  
  # Ganoderma carnosum
  mutate(Species=str_replace(Species, "Ganoderma carnosum", "Ganoderma carnosum cluster"))%>%
  mutate(Species=str_replace(Species, "Ganoderma lucidum", "Ganoderma carnosum cluster"))%>%
  
  # Ganoderma cupreolaccatum
  mutate(Species=str_replace(Species, "Ganoderma cupreolaccatum", "Ganoderma cupreolaccatum cluster"))%>%
  mutate(Species=str_replace(Species, "Ganoderma mutabile", "Ganoderma cupreolaccatum cluster"))%>%
  mutate(Species=str_replace(Species, "Ganoderma pfeifferi", "Ganoderma cupreolaccatum cluster"))%>%
  
  # Ganoderma curtisii
  mutate(Species=str_replace(Species, "Ganoderma curtisii", "Ganoderma curtisii cluster"))%>%
  mutate(Species=str_replace(Species, "Ganoderma meredithae", "Ganoderma curtisii cluster"))%>%
  
  # Ganoderma destructans
  mutate(Species=str_replace(Species, "Ganoderma destructans", "Ganoderma destructans cluster"))%>%
  mutate(Species=str_replace(Species, "Ganoderma steyaertanum", "Ganoderma destructans cluster"))%>%
  
  # Ganoderma sichuanense
  mutate(Species=str_replace(Species, "Ganoderma sichuanense", "Ganoderma sichuanense cluster"))%>%
  mutate(Species=str_replace(Species, "Ganoderma tropicum", "Ganoderma sichuanense cluster"))%>%
  
  #####################################
######### Erysiphe species ######
#####################################

# Erysiphe abeliana
mutate(Species=str_replace(Species, "Erysiphe abeliana", "Erysiphe abeliana cluster"))%>%
  mutate(Species=str_replace(Species, "Erysiphe chifengensis", "Erysiphe abeliana cluster"))%>%
  
  # Erysiphe flexibilis
  mutate(Species=str_replace(Species, "Erysiphe flexibilis", "Erysiphe flexibilis cluster"))%>%
  mutate(Species=str_replace(Species, "Erysiphe symphoricarpi", "Erysiphe flexibilis cluster"))%>%
  
  # Erysiphe havrylenkoana
  mutate(Species=str_replace(Species, "Erysiphe havrylenkoana", "Erysiphe havrylenkoana cluster"))%>%
  mutate(Species=str_replace(Species, "Erysiphe nothofagi", "Erysiphe havrylenkoana cluster"))%>%
  
  # Erysiphe kissiana
  mutate(Species=str_replace(Species, "Erysiphe kissiana", "Erysiphe kissiana cluster"))%>%
  mutate(Species=str_replace(Species, "Erysiphe monoperidiata", "Erysiphe kissiana cluster"))%>%
  
  # Erysiphe miyabei
  mutate(Species=str_replace(Species, "Erysiphe miyabei", "Erysiphe miyabei cluster"))%>%
  mutate(Species=str_replace(Species, "Erysiphe verniciferae", "Erysiphe miyabei cluster"))%>%
  
  #####################################
######### Elsinoe species ######
#####################################

# Elsinoe caleae cluster
mutate(Species=str_replace(Species, "Elsinoe caleae", "Elsinoe caleae cluster"))%>%
  mutate(Species=str_replace(Species, "Elsinoe diospyri", "Elsinoe caleae cluster"))%>%
  
  #####################################
######### Penicillium species ######
#####################################

# Penicillium alogum cluster
mutate(Species=str_replace(Species, " alogum", "Penicillium alogum cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium Penicilliumsubarcticum", "Penicillium alogum cluster"))%>%
  
  # Penicillium amphipolaria cluster
  mutate(Species=str_replace(Species, "Penicillium amphipolaria", "Penicillium amphipolaria cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium singorense", "Penicillium amphipolaria cluster"))%>%
  
  # Penicillium austrosinicum cluster
  mutate(Species=str_replace(Species, "Penicillium austrosinicum", "Penicillium austrosinicum cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium maximae", "Penicillium austrosinicum cluster"))%>%
  
  # Penicillium bissettii cluster
  mutate(Species=str_replace(Species, "Penicillium bissettii", "Penicillium bissettii cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium mariae-crucis", "Penicillium bissettii cluster"))%>%
  
  # Penicillium brefeldianum cluster
  mutate(Species=str_replace(Species, "Penicillium brefeldianum", "Penicillium brefeldianum cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium limosum", "Penicillium brefeldianum cluster"))%>%
  
  # Penicillium viticola cluster"
  mutate(Species=str_replace(Species, "Penicillium viticola", "Penicillium viticola cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium jacksonii", "Penicillium viticola cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium sclerotiorum", "Penicillium viticola cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium cainii", "Penicillium viticola cluster"))%>%
  
  # Penicillium christenseniae cluster
  mutate(Species=str_replace(Species, "Penicillium christenseniae", "Penicillium christenseniae cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium manginii", "Penicillium christenseniae cluster"))%>%
  
  # Penicillium chrysogenum
  mutate(Species=str_replace(Species, "Penicillium chrysogenum", "Penicillium chrysogenum cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium fluviserpens", "Penicillium chrysogenum cluster"))%>%
  
  # Penicillium dokdoense cluster
  mutate(Species=str_replace(Species, "Penicillium dokdoense", "Penicillium dokdoense cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium terrigenum", "Penicillium dokdoense cluster"))%>%
  
  # Penicillium echinulatum cluster
  mutate(Species=str_replace(Species, "Penicillium echinulatum", "Penicillium echinulatum cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium infrapurpureum", "Penicillium echinulatum cluster"))%>%
  
  # Penicillium labradoris cluster
  mutate(Species=str_replace(Species, "Penicillium labradoris", "Penicillium labradoris cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium parvum", "Penicillium labradoris cluster"))%>%
  
  # Penicillium lagena cluster
  mutate(Species=str_replace(Species, "Penicillium lagena", "Penicillium lagena cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium oregonense", "Penicillium lagena cluster"))%>%
  
  # Penicillium lemhiflumine cluster
  mutate(Species=str_replace(Species, "Penicillium lemhiflumine", "Penicillium lemhiflumine cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium salmoniflumine", "Penicillium lemhiflumine cluster"))%>%
  
  # Penicillium levitum cluster
  mutate(Species=str_replace(Species, "Penicillium levitum", "Penicillium levitum cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium meloforme", "Penicillium levitum cluster"))%>%
  
  # Penicillium malacosphaerulum cluster
  mutate(Species=str_replace(Species, "Penicillium malacosphaerulum", "Penicillium malacosphaerulum cluster"))%>%
  mutate(Species=str_replace(Species, "Penicillium reticulisporum", "Penicillium malacosphaerulum cluster"))%>%
  
  #####################################
######### Puccinia species ######
#####################################

# Puccinia allii cluster
mutate(Species=str_replace(Species, "Puccinia allii", "Puccinia allii cluster"))%>%
  mutate(Species=str_replace(Species, "Puccinia mixta", "Puccinia allii cluster"))%>%
  
  # Puccinia xanthii cluster
  mutate(Species=str_replace(Species, "Puccinia xanthii", "Puccinia xanthii cluster"))%>%
  mutate(Species=str_replace(Species, "Puccinia melampodii", "Puccinia xanthii cluster"))%>%
  
  #####################################
######### Diaporthe species ######
#####################################

# Diaporthe cinnamomi cluster
mutate(Species=str_replace(Species, "Diaporthe cinnamomi", "Diaporthe cinnamomi cluster"))%>%
  mutate(Species=str_replace(Species, "Diaporthe discoidispora", "Diaporthe cinnamomi cluster"))%>%
  
  # Diaporthe cynaroidis cluster
  mutate(Species=str_replace(Species, "Diaporthe cynaroidis", "Diaporthe cynaroidis cluster"))%>%
  mutate(Species=str_replace(Species, "Diaporthe italiana", "Diaporthe cynaroidis cluster"))%>%
  
  #####################################
######### Aspergillus species ######
#####################################

# Aspergillus aeneus cluster
mutate(Species=str_replace(Species, "Aspergillus aeneus", "Aspergillus aeneus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus karnatakaensis", "Aspergillus aeneus cluster"))%>%
  
  # Aspergillus aureolatus cluster
  mutate(Species=str_replace(Species, "Aspergillus aureolatus", "Aspergillus aureolatus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus viridicatenatus", "Aspergillus aureolatus cluster"))%>%
  
  # Aspergillus terreus cluster
  mutate(Species=str_replace(Species, "Aspergillus terreus", "Aspergillus terreus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus aureoterreus", "Aspergillus terreus cluster"))%>%
  
  # Aspergillus biplanus cluster
  mutate(Species=str_replace(Species, "Aspergillus biplanus", "Aspergillus biplanus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus diversus", "Aspergillus biplanus cluster"))%>%
  
  # Aspergillus brevipes cluster
  mutate(Species=str_replace(Species, "Aspergillus brevipes", "Aspergillus brevipes cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus felis", "Aspergillus brevipes cluster"))%>%
  
  # Aspergillus caatingaensis cluster
  mutate(Species=str_replace(Species, "Aspergillus caatingaensis", "Aspergillus caatingaensis cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus shendaweii", "Aspergillus caatingaensis cluster"))%>%
  
  # Aspergillus caesiellus cluster
  mutate(Species=str_replace(Species, "Aspergillus caesiellus", "Aspergillus caesiellus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus pachycaulis", "Aspergillus caesiellus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus restrictus", "Aspergillus caesiellus cluster"))%>%
  
  # Aspergillus curvatus cluster
  mutate(Species=str_replace(Species, "Aspergillus curvatus", "Aspergillus curvatus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus flocculosus", "Aspergillus curvatus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus insulicola", "Aspergillus curvatus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus ochraceopetaliformis", "Aspergillus curvatus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus pseudoelegans", "Aspergillus curvatus cluster"))%>%
  
  # Aspergillus elegans cluster
  mutate(Species=str_replace(Species, "Aspergillus elegans", "Aspergillus elegans cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus steynii", "Aspergillus elegans cluster"))%>%
  
  # Aspergillus elsenburgensis cluster
  mutate(Species=str_replace(Species, "Aspergillus elsenburgensis", "Aspergillus elsenburgensis cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus galapagensis", "Aspergillus elsenburgensis cluster"))%>%
  
  # Aspergillus fumigatiaffinis cluster
  mutate(Species=str_replace(Species, "Aspergillus fumigatiaffinis", "Aspergillus fumigatiaffinis cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus laciniosus", "Aspergillus fumigatiaffinis cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus novofumigatus", "Aspergillus fumigatiaffinis cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus udagawae", "Aspergillus fumigatiaffinis cluster"))%>%
  
  # Aspergillus germanicus cluster
  mutate(Species=str_replace(Species, "Aspergillus germanicus", "Aspergillus germanicus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus ustus", "Aspergillus germanicus cluster"))%>%
  
  # Aspergillus hiratsukae cluster
  mutate(Species=str_replace(Species, "Aspergillus hiratsukae", "Aspergillus hiratsukae cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus marvanovae", "Aspergillus hiratsukae cluster"))%>%
  
  # Aspergillus iranicus cluster
  mutate(Species=str_replace(Species, "Aspergillus iranicus", "Aspergillus iranicus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus neoindicus", "Aspergillus iranicus cluster"))%>%
  
  # Aspergillus lanosus cluster
  mutate(Species=str_replace(Species, "Aspergillus lanosus", "Aspergillus lanosus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus vandermerwei", "Aspergillus lanosus cluster"))%>%
  
  # Aspergillus multicolor cluster
  mutate(Species=str_replace(Species, "Aspergillus multicolor", "Aspergillus multicolor cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus mulundensis", "Aspergillus multicolor cluster"))%>%
  
  # Aspergillus nidulans cluster
  mutate(Species=str_replace(Species, "Aspergillus nidulans", "Aspergillus nidulans cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus rugulosus", "Aspergillus nidulans cluster"))%>%
  
  # Aspergillus nishimurae cluster
  mutate(Species=str_replace(Species, "Aspergillus nishimurae", "Aspergillus nishimurae cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus viridinutans", "Aspergillus nishimurae cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus waksmanii", "Aspergillus nishimurae cluster"))%>%
  
  # Aspergillus phialiformis cluster
  mutate(Species=str_replace(Species, "Aspergillus phialiformis", "Aspergillus phialiformis cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus phialosimplex", "Aspergillus phialiformis cluster"))%>%
  
  # Aspergillus piperis cluster
  mutate(Species=str_replace(Species, "Aspergillus piperis", "Aspergillus piperis cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus welwitschiae", "Aspergillus piperis cluster"))%>%
  
  # Aspergillus tasmanicus cluster
  mutate(Species=str_replace(Species, "Aspergillus tasmanicus", "Aspergillus tasmanicus cluster"))%>%
  mutate(Species=str_replace(Species, "Aspergillus tatenoi", "Aspergillus tasmanicus cluster"))%>%
  
  #####################################
######### Curvularia species ######
#####################################

# Curvularia affinis cluster
mutate(Species=str_replace(Species, "Curvularia affinis", "Curvularia affinis cluster"))%>%
  mutate(Species=str_replace(Species, "Curvularia verruciformis", "Curvularia affinis cluster"))%>%
  
  # Curvularia akaii cluster
  mutate(Species=str_replace(Species, "Curvularia akaii", "Curvularia akaii cluster"))%>%
  mutate(Species=str_replace(Species, "Curvularia bothriochloae", "Curvularia akaii cluster"))%>%
  
  # Curvularia boeremae cluster
  mutate(Species=str_replace(Species, "Curvularia boeremae", "Curvularia boeremae cluster"))%>%
  mutate(Species=str_replace(Species, "Curvularia sesuvii", "Curvularia boeremae cluster"))%>%
  
  # Curvularia fallax cluster
  mutate(Species=str_replace(Species, "Curvularia fallax", "Curvularia fallax cluster"))%>%
  mutate(Species=str_replace(Species, "Curvularia lycopersici", "Curvularia fallax cluster"))%>%
  
  # Curvularia heteropogonis cluster
  mutate(Species=str_replace(Species, "Curvularia heteropogonis", "Curvularia heteropogonis cluster"))%>%
  mutate(Species=str_replace(Species, "Curvularia tsudae", "Curvularia heteropogonis cluster"))%>%
  
  # Curvularia pseudorobusta cluster
  mutate(Species=str_replace(Species, "Curvularia pseudorobusta", "Curvularia pseudorobusta cluster"))%>%
  mutate(Species=str_replace(Species, "Curvularia robusta", "Curvularia pseudorobusta cluster"))%>%
  
  #####################################
######### Ramularia species ######
#####################################

# Ramularia citricola cluster
mutate(Species=str_replace(Species, "Ramularia citricola", "Ramularia citricola cluster"))%>%
  mutate(Species=str_replace(Species, "Ramularia hydrangeae-macrophyllae", "Ramularia citricola cluster"))%>%
  
  # Ramularia diervillae cluster
  mutate(Species=str_replace(Species, "Ramularia diervillae", "Ramularia diervillae cluster"))%>%
  mutate(Species=str_replace(Species, "Ramularia geraniicola", "Ramularia diervillae cluster"))%>%
  
  #####################################
######### Pseudocercospora species ######
#####################################

# Pseudocercospora araliae cluster
mutate(Species=str_replace(Species, "Pseudocercospora araliae", "Pseudocercospora araliae cluster"))%>%
  mutate(Species=str_replace(Species, "Pseudocercospora arecacearum", "Pseudocercospora araliae cluster"))%>%
  
  # Pseudocercospora leandrae-fragilis cluster
  mutate(Species=str_replace(Species, "Pseudocercospora leandrae-fragilis", "Pseudocercospora leandrae-fragilis cluster"))%>%
  mutate(Species=str_replace(Species, "Pseudocercospora sphaerulinae", "Pseudocercospora leandrae-fragilis cluster"))%>%
  
  #####################################
######### Cercospora species ######
#####################################

# Cercospora capsici cluster
mutate(Species=str_replace(Species, "Cercospora capsici", "Cercospora capsici cluster"))%>%
  mutate(Species=str_replace(Species, "Cercospora coniogrammes", "Cercospora capsici cluster"))%>%
  
  #####################################
######### Podosphaera species ######
#####################################

# Podosphaera clandestina cluster
mutate(Species=str_replace(Species, "Podosphaera clandestina", "Podosphaera clandestina cluster"))%>%
  mutate(Species=str_replace(Species, "Podosphaera prunicola", "Podosphaera clandestina cluster"))%>%
  
  # Podosphaera euphorbiae cluster
  mutate(Species=str_replace(Species, "Podosphaera euphorbiae", "Podosphaera euphorbiae cluster"))%>%
  mutate(Species=str_replace(Species, "Podosphaera lini", "Podosphaera euphorbiae cluster"))%>%
  
  #####################################
######### Golovinomyces species ######
#####################################

# Golovinomyces glandulariae cluster
mutate(Species=str_replace(Species, "Golovinomyces glandulariae", "Golovinomyces glandulariae cluster"))%>%
  mutate(Species=str_replace(Species, "Golovinomyces magnicellulatus", "Golovinomyces glandulariae cluster"))%>%
  
  #####################################
######### Colletotrichum species (need to double check C.acutatum) ######
#####################################

# Colletotrichum acerbum cluster
mutate(Species=str_replace(Species, "Colletotrichum acerbum", "Colletotrichum acerbum cluster"))%>%
  mutate(Species=str_replace(Species, "Colletotrichum rhombiforme", "Colletotrichum acerbum cluster"))%>%
  
  # Colletotrichum cariniferi cluster
  mutate(Species=str_replace(Species, "Colletotrichum cariniferi", "Colletotrichum cariniferi cluster"))%>%
  mutate(Species=str_replace(Species, "Colletotrichum coelogynes", "Colletotrichum cariniferi cluster"))%>%
  
  # Colletotrichum dracaenophilum cluster
  mutate(Species=str_replace(Species, "Colletotrichum dracaenophilum", "Colletotrichum dracaenophilum cluster"))%>%
  mutate(Species=str_replace(Species, "Colletotrichum sinuatum", "Colletotrichum dracaenophilum cluster"))%>%
  
  #####################################
######### Bipolaris species ######
#####################################

# Bipolaris brachiariae cluster
mutate(Species=str_replace(Species, "Bipolaris brachiariae", "Bipolaris brachiariae cluster"))%>%
  mutate(Species=str_replace(Species, "Bipolaris saccharicola", "Bipolaris brachiariae cluster"))%>%
  
  # Bipolaris microstegii cluster
  mutate(Species=str_replace(Species, "Bipolaris microstegii", "Bipolaris microstegii cluster"))%>%
  mutate(Species=str_replace(Species, "Bipolaris urochloae", "Bipolaris microstegii cluster"))%>%
  
  # Bipolaris stenospila cluster
  mutate(Species=str_replace(Species, "Bipolaris stenospila", "Bipolaris stenospila cluster"))%>%
  mutate(Species=str_replace(Species, "Bipolaris yamadai", "Bipolaris stenospila cluster"))%>%
  
  #####################################
######### Microdochium species ######
#####################################

# Microdochium poae cluster
mutate(Species=str_replace(Species, "Microdochium poae", "Microdochium poae cluster"))%>%
  mutate(Species=str_replace(Species, "Microdochium colombiense", "Microdochium poae cluster"))%>%
  
  #####################################
######### Venturia species ######
#####################################

# Venturia martianoffiana cluster
mutate(Species=str_replace(Species, "Venturia martianoffiana", "Venturia martianoffiana cluster"))%>%
  mutate(Species=str_replace(Species, "Venturia phaeosepta", "Venturia martianoffiana cluster"))%>%
  
  # Venturia populina cluster
  mutate(Species=str_replace(Species, "Venturia populina", "Venturia populina cluster"))%>%
  mutate(Species=str_replace(Species, "Venturia saliciperda", "Venturia populina cluster"))%>%
  
  #####################################
######### Trichoderma species ######
#####################################

# Trichoderma alni cluster
mutate(Species=str_replace(Species, "Trichoderma alni", "Trichoderma alni cluster"))%>%
  mutate(Species=str_replace(Species, "Trichoderma rufobrunneum", "Trichoderma alni cluster"))%>%
  
  # Trichoderma hypoxylon cluster
  mutate(Species=str_replace(Species, "Trichoderma hypoxylon", "Trichoderma hypoxylon cluster"))%>%
  mutate(Species=str_replace(Species, "Trichoderma taxi", "Trichoderma hypoxylon cluster"))%>%
  
  # Trichoderma protopulvinatum cluster
  mutate(Species=str_replace(Species, "Trichoderma protopulvinatum", "Trichoderma protopulvinatum cluster"))%>%
  mutate(Species=str_replace(Species, "Trichoderma pulvinatum", "Trichoderma protopulvinatum cluster"))%>%
  
  #####################################
######### Uromyces species ######
#####################################

# Uromyces laburni cluster
mutate(Species=str_replace(Species, "Uromyces laburni", "Uromyces laburni cluster"))%>%
  mutate(Species=str_replace(Species, "Uromyces punctatus", "Uromyces laburni cluster"))%>%
  
  #####################################
######### Ustilago species ######
#####################################

# Ustilago airae-caespitosae cluster
mutate(Species=str_replace(Species, "Ustilago airae-caespitosae", "Ustilago airae-caespitosae cluster"))%>%
  mutate(Species=str_replace(Species, "Ustilago nunavutica", "Ustilago airae-caespitosae cluster"))%>%
  
  # Ustilago avenae cluster
  mutate(Species=str_replace(Species, "Ustilago avenae", "Ustilago avenae cluster"))%>%
  mutate(Species=str_replace(Species, "Ustilago hordei", "Ustilago avenae cluster"))%>%
  
  #####################################
######### [Candida] species ######
#####################################

# [Candida] glaebosa cluster
mutate(Species=str_replace(Species, "[Candida] glaebosa", "[Candida] glaebosa cluster"))%>%
  mutate(Species=str_replace(Species, "[Candida] pseudoglaebosa", "[Candida] glaebosa cluster"))%>%
  
  # [Candida] membranifaciens cluster
  mutate(Species=str_replace(Species, "[Candida] membranifaciens", "[Candida] membranifaciens cluster"))%>%
  mutate(Species=str_replace(Species, "Candida friedrichii", "[Candida] membranifaciens cluster"))%>%
  
  #####################################
######### Mucor species ######
#####################################

# Mucor atramentarius cluster
mutate(Species=str_replace(Species, "Mucor atramentarius", "Mucor atramentarius cluster"))%>%
  mutate(Species=str_replace(Species, "Mucor bainieri", "Mucor atramentarius cluster"))%>%
  
  # Mucor irregularis cluster
  mutate(Species=str_replace(Species, "Mucor irregularis", "Mucor irregularis cluster"))%>%
  mutate(Species=str_replace(Species, "Mucor takensis", "Mucor irregularis cluster"))%>%
  
  # Mucor koreanus cluster
  mutate(Species=str_replace(Species, "Mucor koreanus", "Mucor koreanus cluster"))%>%
  mutate(Species=str_replace(Species, "Mucor piriformis", "Mucor koreanus cluster"))%>%
  
  #####################################
######### Tilletia species ######
#####################################

# Tilletia holci cluster
mutate(Species=str_replace(Species, "Tilletia holci", "Tilletia holci cluster"))%>%
  mutate(Species=str_replace(Species, "Tilletia olida", "Tilletia holci cluster"))%>%
  
  #####################################
######### Pestalotiopsis species ####
#####################################

# Pestalotiopsis arceuthobii cluster
mutate(Species=str_replace(Species, "Pestalotiopsis arceuthobii", "Pestalotiopsis arceuthobii cluster"))%>%
  mutate(Species=str_replace(Species, "Pestalotiopsis diversiseta", "Pestalotiopsis arceuthobii cluster"))%>%
  
  # Pestalotiopsis iberica cluster
  mutate(Species=str_replace(Species, "Pestalotiopsis iberica", "Pestalotiopsis iberica cluster"))%>%
  mutate(Species=str_replace(Species, "Pestalotiopsis rhododendri", "Pestalotiopsis iberica cluster"))%>%
  
  #####################################
######### Ophiostoma species ####
#####################################

# Ophiostoma ainoae cluster
mutate(Species=str_replace(Species, "Ophiostoma ainoae", "Ophiostoma ainoae cluster"))%>%
  mutate(Species=str_replace(Species, "Ophiostoma brevipilosi", "Ophiostoma ainoae cluster"))%>%
  
  # Ophiostoma arduennense cluster
  mutate(Species=str_replace(Species, "Ophiostoma arduennense", "Ophiostoma arduennense cluster"))%>%
  mutate(Species=str_replace(Species, "Ophiostoma taphrorychi", "Ophiostoma arduennense cluster"))%>%
  
  # Ophiostoma sparsiannulatum cluster
  mutate(Species=str_replace(Species, "Ophiostoma sparsiannulatum", "Ophiostoma sparsiannulatum cluster"))%>%
  mutate(Species=str_replace(Species, "Ophiostoma subannulatum", "Ophiostoma sparsiannulatum cluster"))%>%
  
  #####################################
######### Mortierella species ####
#####################################

# Mortierella biramosa cluster
mutate(Species=str_replace(Species, "Mortierella biramosa", "Mortierella biramosa cluster"))%>%
  mutate(Species=str_replace(Species, "Mortierella nantahalensis", "Mortierella biramosa cluster"))%>%
  
  #####################################
######### Metarhizium species ####
#####################################

# Metarhizium chaiyaphumense cluster
mutate(Species=str_replace(Species, "Metarhizium chaiyaphumense", "Metarhizium chaiyaphumense cluster"))%>%
  mutate(Species=str_replace(Species, "Metarhizium cylindrosporae", "Metarhizium chaiyaphumense cluster"))%>%
  
  # Metarhizium flavoviride cluster
  mutate(Species=str_replace(Species, "Metarhizium flavoviride", "Metarhizium flavoviride cluster"))%>%
  mutate(Species=str_replace(Species, "Metarhizium gaoligongense", "Metarhizium flavoviride cluster"))%>%
  
  #####################################
######### Cortinarius species ####
#####################################

# cluster 2
mutate(Species=str_replace(Species, "Cortinarius abiegnus", "Cortinarius abiegnus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius millaresensis", "Cortinarius abiegnus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius uxorum", "Cortinarius abiegnus cluster"))%>%
  
  # cluster 3
  mutate(Species=str_replace(Species, "Cortinarius acutispissipes", "Cortinarius acutispissipes cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius alboviolaceus", "Cortinarius acutispissipes cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius obliquus", "Cortinarius acutispissipes cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius paralbocyaneus", "Cortinarius acutispissipes cluster"))%>%
  
  # cluster 9
  mutate(Species=str_replace(Species, "Cortinarius aff. orellanus", "Cortinarius aff. orellanus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius renidens", "Cortinarius aff. orellanus cluster"))%>%
  
  # cluster 10
  mutate(Species=str_replace(Species, "Cortinarius aff. torvus", "Cortinarius aff. torvus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius multicolor", "Cortinarius aff. torvus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subionophyllus", "Cortinarius aff. torvus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius ultimiionophyllus", "Cortinarius aff. torvus cluster"))%>%
  
  # cluster 12
  mutate(Species=str_replace(Species, "Cortinarius aff. infractus", "Cortinarius aff. infractus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius infractiflavus", "Cortinarius aff. infractus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius infractus", "Cortinarius aff. infractus cluster"))%>%
  
  # cluster 13
  mutate(Species=str_replace(Species, "Cortinarius ahsii", "Cortinarius ahsii cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius vernalisierraensis", "Cortinarius ahsii cluster"))%>%
  
  # cluster 21
  mutate(Species=str_replace(Species, "Cortinarius albobrunnescens", "Cortinarius albobrunnescens cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius murinascens", "Cortinarius albobrunnescens cluster"))%>%
  
  # cluster 33
  mutate(Species=str_replace(Species, "Cortinarius allutus", "Cortinarius allutus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius fulminoides", "Cortinarius allutus cluster"))%>%
  
  # cluster 35
  mutate(Species=str_replace(Species, "Cortinarius alpinus", "Cortinarius alpinus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius favrei", "Cortinarius alpinus cluster"))%>%
  
  # cluster 36
  mutate(Species=str_replace(Species, "Cortinarius alticaudus", "Cortinarius alticaudus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius glaucopus", "Cortinarius alticaudus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius glaucopus", "Cortinarius alticaudus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius pansa", "Cortinarius alticaudus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subfoetens", "Cortinarius alticaudus cluster"))%>%
  
  # cluster 39
  mutate(Species=str_replace(Species, "Cortinarius ammophiloides", "Cortinarius ammophiloides cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius caliginosus", "Cortinarius ammophiloides cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius gallurae", "Cortinarius ammophiloides cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius punctatiformis", "Cortinarius ammophiloides cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius sagacitas", "Cortinarius ammophiloides cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subobtusobrunneus", "Cortinarius ammophiloides cluster"))%>%
  
  # cluster 42
  mutate(Species=str_replace(Species, "Cortinarius andreae", "Cortinarius andreae cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius pinisquamulosus", "Cortinarius andreae cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius roseocastaneus", "Cortinarius andreae cluster"))%>%
  
  # cluster 51
  mutate(Species=str_replace(Species, "Cortinarius anomalopacificus", "Cortinarius anomalopacificus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius tabularis", "Cortinarius anomalopacificus cluster"))%>%
  
  # cluster 75
  mutate(Species=str_replace(Species, "Cortinarius auroripes", "Cortinarius auroripes cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius campester", "Cortinarius auroripes cluster"))%>%
  
  # cluster 79
  mutate(Species=str_replace(Species, "Cortinarius azureus", "Cortinarius azureus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius nymphatus", "Cortinarius azureus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius pseudobulliardioides", "Cortinarius azureus cluster"))%>%
  
  # cluster 80
  mutate(Species=str_replace(Species, "Cortinarius badioflavidus", "Cortinarius badioflavidus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius hinnuleocervinus", "Cortinarius badioflavidus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius hinnuloides", "Cortinarius badioflavidus cluster"))%>%
  
  # cluster 89
  mutate(Species=str_replace(Species, "Cortinarius balteatus", "Cortinarius balteatus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius pseudonaevosus", "Cortinarius balteatus cluster"))%>%
  
  # cluster 99
  mutate(Species=str_replace(Species, "Cortinarius biriensis", "Cortinarius biriensis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius bonachei", "Cortinarius biriensis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius ferrugineovelatus", "Cortinarius biriensis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius griseovioleipes", "Cortinarius biriensis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius intempestivus", "Cortinarius biriensis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius miwok", "Cortinarius biriensis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius roseobasilis", "Cortinarius biriensis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subpuellaris", "Cortinarius biriensis cluster"))%>%
  
  # cluster 100
  mutate(Species=str_replace(Species, "Cortinarius bivelosimilis", "Cortinarius bivelosimilis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius bivelus", "Cortinarius bivelosimilis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius brunneifolius", "Cortinarius bivelosimilis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius cf.", "Cortinarius bivelosimilis cluster"))%>%
  mutate(Species=str_replace(Species, "Ectomycorrhizal root", "Cortinarius bivelosimilis cluster"))%>%
  
  # cluster 101
  mutate(Species=str_replace(Species, "Cortinarius borealis", "Cortinarius borealis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius selinolens", "Cortinarius borealis cluster"))%>%
  
  # cluster 105
  mutate(Species=str_replace(Species, "Cortinarius boreotrichus", "Cortinarius boreotrichus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius expallens", "Cortinarius boreotrichus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subpaleaceus", "Cortinarius boreotrichus cluster"))%>%
  
  # cluster 112
  mutate(Species=str_replace(Species, "Cortinarius bridgei", "Cortinarius bridgei cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius flavobasilis", "Cortinarius bridgei cluster"))%>%
  
  # cluster 114
  mutate(Species=str_replace(Species, "Cortinarius brunneoaffinis", "Cortinarius brunneoaffinis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius centrirufus", "Cortinarius brunneoaffinis cluster"))%>%
  
  # cluster 126
  mutate(Species=str_replace(Species, "Cortinarius cadi-aguirrei", "Cortinarius cadi-aguirrei cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius fulvescens", "Cortinarius cadi-aguirrei cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius fulvescentoideus", "Cortinarius cadi-aguirrei cluster"))%>%
  
  # cluster 133
  mutate(Species=str_replace(Species, "Cortinarius cagei", "Cortinarius cagei cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius evernius", "Cortinarius cagei cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius paranomalus", "Cortinarius cagei cluster"))%>%
  
  # cluster 146
  mutate(Species=str_replace(Species, "Cortinarius casimiri", "Cortinarius casimiri cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius decipientoides", "Cortinarius casimiri cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius hirtus", "Cortinarius casimiri cluster"))%>%
  
  # cluster 148
  mutate(Species=str_replace(Species, "Cortinarius castaneopallidus", "Cortinarius castaneopallidus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius quercoconicus", "Cortinarius castaneopallidus cluster"))%>%
  
  # cluster 151
  mutate(Species=str_replace(Species, "Cortinarius cepistipes", "Cortinarius cepistipes cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius pescolanensis", "Cortinarius cepistipes cluster"))%>%
  
  # cluster 153
  mutate(Species=str_replace(Species, "Cortinarius cf.", "Cortinarius cf. cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius cinnamomeus", "Cortinarius cf. cluster"))%>%
  
  # cluster 155
  mutate(Species=str_replace(Species, "Cortinarius chevassutii", "Cortinarius chevassutii cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius crassisporus", "Cortinarius chevassutii cluster"))%>%
  
  # cluster 157
  mutate(Species=str_replace(Species, "Cortinarius chrysma", "Cortinarius chrysma cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius icterinoides", "Cortinarius chrysma cluster"))%>%
  
  # cluster 171
  mutate(Species=str_replace(Species, "Cortinarius clarisordidus", "Cortinarius clarisordidus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius clarobrunneus", "Cortinarius clarisordidus cluster"))%>%
  
  # cluster 190
  mutate(Species=str_replace(Species, "Cortinarius cremeolina", "Cortinarius cremeolina cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius rapaceus", "Cortinarius cremeolina cluster"))%>%
  
  # cluster 201
  mutate(Species=str_replace(Species, "Cortinarius danili", "Cortinarius danili cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius lindstroemii", "Cortinarius danili cluster"))%>%
  
  # cluster 205
  mutate(Species=str_replace(Species, "Cortinarius decipiens", "Cortinarius decipiens cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subturibulosus", "Cortinarius decipiens cluster"))%>%
  
  # cluster 211
  mutate(Species=str_replace(Species, "Cortinarius diasemospermus", "Cortinarius diasemospermus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius geraniolens", "Cortinarius diasemospermus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius paezii", "Cortinarius diasemospermus cluster"))%>%
  
  # cluster 215
  mutate(Species=str_replace(Species, "Cortinarius diosmus", "Cortinarius diosmus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius urbicus", "Cortinarius diosmus cluster"))%>%
  
  # cluster 217
  mutate(Species=str_replace(Species, "Cortinarius dolabratoides", "Cortinarius dolabratoides cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius turgidipes", "Cortinarius dolabratoides cluster"))%>%
  
  # cluster 232
  mutate(Species=str_replace(Species, "Cortinarius elegantio-montanus", "Cortinarius elegantio-montanus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius elegantio-occidentalis", "Cortinarius elegantio-montanus cluster"))%>%
  
  # cluster 240
  mutate(Species=str_replace(Species, "Cortinarius exsularis", "Cortinarius exsularis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius vinaceogrisescens", "Cortinarius exsularis cluster"))%>%
  
  # cluster 242
  mutate(Species=str_replace(Species, "Cortinarius ferrei", "Cortinarius ferrei cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius roseivelatus", "Cortinarius ferrei cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subcastaneus", "Cortinarius ferrei cluster"))%>%
  
  # cluster 257
  mutate(Species=str_replace(Species, "Cortinarius fraudulosus", "Cortinarius fraudulosus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius patrickensis", "Cortinarius fraudulosus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subfraudulosus", "Cortinarius fraudulosus cluster"))%>%
  
  # cluster 258
  mutate(Species=str_replace(Species, "Cortinarius fructuodorus", "Cortinarius fructuodorus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius leucophaeatus", "Cortinarius fructuodorus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius tigrinipes", "Cortinarius fructuodorus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius venustus", "Cortinarius fructuodorus cluster"))%>%
  
  # cluster 261
  mutate(Species=str_replace(Species, "Cortinarius fulvopaludosus", "Cortinarius fulvopaludosus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius umbrinolens", "Cortinarius fulvopaludosus cluster"))%>%
  
  # cluster 263
  mutate(Species=str_replace(Species, "Cortinarius fuscescens", "Cortinarius fuscescens cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius lobatus", "Cortinarius fuscescens cluster"))%>%
  
  # cluster 265
  mutate(Species=str_replace(Species, "Cortinarius fuscobovinaster", "Cortinarius fuscobovinaster cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius fuscobovinus", "Cortinarius fuscobovinaster cluster"))%>%
  
  # cluster 266
  mutate(Species=str_replace(Species, "Cortinarius fuscodiscus", "Cortinarius fuscodiscus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius pseudofusisporus", "Cortinarius fuscodiscus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius rusticellus", "Cortinarius fuscodiscus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius semivestitus", "Cortinarius fuscodiscus cluster"))%>%
  
  # cluster 282
  mutate(Species=str_replace(Species, "Cortinarius griseascens", "Cortinarius griseascens cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius herculinus", "Cortinarius griseascens cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius luridus", "Cortinarius griseascens cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius roseonudipes", "Cortinarius griseascens cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius semiodoratus", "Cortinarius griseascens cluster"))%>%
  
  # cluster 283
  mutate(Species=str_replace(Species, "Cortinarius griseocoeruleus", "Cortinarius griseocoeruleus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius leonicolor", "Cortinarius griseocoeruleus cluster"))%>%
  
  # cluster 285
  mutate(Species=str_replace(Species, "Cortinarius gualalaensis", "Cortinarius gualalaensis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius tuolumnensis", "Cortinarius gualalaensis cluster"))%>%
  
  # cluster 290
  mutate(Species=str_replace(Species, "Cortinarius heatherae", "Cortinarius heatherae cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subbrunneus", "Cortinarius heatherae cluster"))%>%
  
  # cluster 291
  mutate(Species=str_replace(Species, "Cortinarius helodes", "Cortinarius helodes cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius lacorum", "Cortinarius helodes cluster"))%>%
  
  # cluster 295
  mutate(Species=str_replace(Species, "Cortinarius heterocycloideus", "Cortinarius heterocycloideus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subheterocyclus", "Cortinarius heterocycloideus cluster"))%>%
  
  # cluster 296
  mutate(Species=str_replace(Species, "Cortinarius heterodepressus", "Cortinarius heterodepressus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius pelargoniostriatulus", "Cortinarius heterodepressus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius roseomyceliosus", "Cortinarius heterodepressus cluster"))%>%
  
  # cluster 302
  mutate(Species=str_replace(Species, "Cortinarius huddartensis", "Cortinarius huddartensis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius lebretonii", "Cortinarius huddartensis cluster"))%>%
  
  # cluster 304
  mutate(Species=str_replace(Species, "Cortinarius humboldtensis", "Cortinarius humboldtensis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius phoeniceus", "Cortinarius humboldtensis cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius tinctorum", "Cortinarius humboldtensis cluster"))%>%
  
  # cluster 310
  mutate(Species=str_replace(Species, "Cortinarius illuminoides", "Cortinarius illuminoides cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius illuminus", "Cortinarius illuminoides cluster"))%>%
  
  # cluster 313
  mutate(Species=str_replace(Species, "Cortinarius indotatus", "Cortinarius indotatus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius luteostriatulus", "Cortinarius indotatus cluster"))%>%
  
  # cluster 319
  mutate(Species=str_replace(Species, "Cortinarius jonimitchelliae", "Cortinarius jonimitchelliae cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius pastoralis", "Cortinarius jonimitchelliae cluster"))%>%
  
  # cluster 321
  mutate(Species=str_replace(Species, "Cortinarius kaimanawa", "Cortinarius kaimanawa cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius rhipiduranus", "Cortinarius kaimanawa cluster"))%>%
  
  # cluster 325
  mutate(Species=str_replace(Species, "Cortinarius kranabetteri", "Cortinarius kranabetteri cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius perrotensis", "Cortinarius kranabetteri cluster"))%>%
  
  # cluster 346
  mutate(Species=str_replace(Species, "Cortinarius luteo-ornatus", "Cortinarius luteo-ornatus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius rossicioenochelis", "Cortinarius luteo-ornatus cluster"))%>%
  
  # cluster 353
  mutate(Species=str_replace(Species, "Cortinarius majorinus", "Cortinarius majorinus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius tenebricus", "Cortinarius majorinus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius vernus", "Cortinarius majorinus cluster"))%>%
  
  # cluster 357
  mutate(Species=str_replace(Species, "Cortinarius megacystidiosus", "Cortinarius megacystidiosus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius pseudoflabellus", "Cortinarius megacystidiosus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius sagarum", "Cortinarius megacystidiosus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius umbrinolutescens", "Cortinarius megacystidiosus cluster"))%>%
  
  # cluster 358
  mutate(Species=str_replace(Species, "Cortinarius melanotus", "Cortinarius melanotus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius viridans", "Cortinarius melanotus cluster"))%>%
  
  # cluster 370
  mutate(Species=str_replace(Species, "Cortinarius minusculus", "Cortinarius minusculus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius paululus", "Cortinarius minusculus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius rubrocinctus", "Cortinarius minusculus cluster"))%>%
  
  # cluster 371
  mutate(Species=str_replace(Species, "Cortinarius misermontii", "Cortinarius misermontii cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius scaurocaninus", "Cortinarius misermontii cluster"))%>%
  
  # cluster 388
  mutate(Species=str_replace(Species, "Cortinarius neofallax", "Cortinarius neofallax cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius pseudofallax", "Cortinarius neofallax cluster"))%>%
  
  # cluster 389
  mutate(Species=str_replace(Species, "Cortinarius neosanguineus", "Cortinarius neosanguineus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius sierraensis", "Cortinarius neosanguineus cluster"))%>%
  
  # cluster 409
  mutate(Species=str_replace(Species, "Cortinarius olididisjungendus", "Cortinarius olididisjungendus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius orasericeus", "Cortinarius olididisjungendus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius piceidisjungendus", "Cortinarius olididisjungendus cluster"))%>%
  
  # cluster 428
  mutate(Species=str_replace(Species, "Cortinarius paludosaniosus", "Cortinarius paludosaniosus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subsaniosus", "Cortinarius paludosaniosus cluster"))%>%
  
  # cluster 430
  mutate(Species=str_replace(Species, "Cortinarius panellus", "Cortinarius panellus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius raphanoides", "Cortinarius panellus cluster"))%>%
  
  # cluster 453
  mutate(Species=str_replace(Species, "Cortinarius perplexus", "Cortinarius perplexus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius pseudogracilior", "Cortinarius perplexus cluster"))%>%
  
  # cluster 455
  mutate(Species=str_replace(Species, "Cortinarius persplendidus", "Cortinarius persplendidus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius sinapivelus", "Cortinarius persplendidus cluster"))%>%
  
  # cluster 467
  mutate(Species=str_replace(Species, "Cortinarius porphyroideus", "Cortinarius porphyroideus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius violaceovolvatus", "Cortinarius porphyroideus cluster"))%>%
  
  # cluster 480
  mutate(Species=str_replace(Species, "Cortinarius pseudophlegma", "Cortinarius pseudophlegma cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subcarneinatus", "Cortinarius pseudophlegma cluster"))%>%
  
  # cluster 481
  mutate(Species=str_replace(Species, "Cortinarius pseudotorvus", "Cortinarius pseudotorvus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius torvus", "Cortinarius pseudotorvus cluster"))%>%
  
  # cluster 484
  mutate(Species=str_replace(Species, "Cortinarius puellaris", "Cortinarius puellaris cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius vinaceobrunneus", "Cortinarius puellaris cluster"))%>%
  
  # cluster 512
  mutate(Species=str_replace(Species, "Cortinarius rubrosanguineus", "Cortinarius rubrosanguineus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius sanguineus", "Cortinarius rubrosanguineus cluster"))%>%
  
  # cluster 526
  mutate(Species=str_replace(Species, "Cortinarius scaurotraganoides", "Cortinarius scaurotraganoides cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subrimosus", "Cortinarius scaurotraganoides cluster"))%>%
  
  # cluster 534
  mutate(Species=str_replace(Species, "Cortinarius serratissimus", "Cortinarius serratissimus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius subserratissimus", "Cortinarius serratissimus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius venustissimus", "Cortinarius serratissimus cluster"))%>%
  
  # cluster 573
  mutate(Species=str_replace(Species, "Cortinarius tacitus", "Cortinarius tacitus cluster"))%>%
  mutate(Species=str_replace(Species, "Cortinarius torvoides", "Cortinarius tacitus cluster"))%>%
  
  #####################################
######### Amanita species ####
#####################################

# cluster 102
mutate(Species=str_replace(Species, "Amanita cf.", "Amanita cf. cluster"))%>%
  mutate(Species=str_replace(Species, "Amanita ovoidea", "Amanita cf. cluster"))%>%
  
  # cluster 106
  mutate(Species=str_replace(Species, "Amanita chiui", "Amanita chiui cluster"))%>%
  mutate(Species=str_replace(Species, "Amanita griseofusca", "Amanita chiui cluster"))%>%
  
  # cluster 125
  mutate(Species=str_replace(Species, "Amanita djarilmari", "Amanita djarilmari cluster"))%>%
  mutate(Species=str_replace(Species, "Amanita gardneri", "Amanita djarilmari cluster"))%>%
  mutate(Species=str_replace(Species, "Amanita millsii", "Amanita djarilmari cluster"))%>%
  
  # cluster 153
  mutate(Species=str_replace(Species, "Amanita griseopantherina", "Amanita griseopantherina cluster"))%>%
  mutate(Species=str_replace(Species, "Amanita pakistanica", "Amanita griseopantherina cluster"))%>%
  
  # cluster 181
  mutate(Species=str_replace(Species, "Amanita mairei", "Amanita mairei cluster"))%>%
  mutate(Species=str_replace(Species, "Amanita mansehraensis", "Amanita mairei cluster"))%>%
  
  # cluster 193
  mutate(Species=str_replace(Species, "Amanita murrilliana", "Amanita murrilliana cluster"))%>%
  mutate(Species=str_replace(Species, "Amanita tuza", "Amanita murrilliana cluster"))%>%
  
  # cluster 228
  mutate(Species=str_replace(Species, "Amanita pruittii", "Amanita pruittii cluster"))%>%
  mutate(Species=str_replace(Species, "Amanita singeri", "Amanita pruittii cluster"))%>%
  
  #####################################
######### Inocybe species ####
#####################################

# cluster 439
mutate(Species=str_replace(Species, "Inocybe kohistanensis", "Inocybe kohistanensis cluster"))%>%
  mutate(Species=str_replace(Species, "Inocybe piceae", "Inocybe kohistanensis cluster"))%>%
  
  # cluster 475
  mutate(Species=str_replace(Species, "Inocybe multifolia", "Inocybe multifolia cluster"))%>%
  mutate(Species=str_replace(Species, "Inocybe subdestricta", "Inocybe multifolia cluster"))%>%
  
  #####################################
######### Acrodontium species ####
#####################################

# cluster 104
mutate(Species=str_replace(Species, "Acrodontium crateriforme", "Acrodontium crateriforme cluster"))%>%
  mutate(Species=str_replace(Species, "Acrodontium neolitseae", "Acrodontium crateriforme cluster"))%>%
  
  #####################################
######### Agaricus species ####
#####################################

# cluster 185
mutate(Species=str_replace(Species, "Agaricus crocodilinus", "Agaricus crocodilinus cluster"))%>%
  mutate(Species=str_replace(Species, "Agaricus julius", "Agaricus crocodilinus cluster"))%>%
  
  # cluster 187
  mutate(Species=str_replace(Species, "Agaricus cupressicola", "Agaricus cupressicola cluster"))%>%
  mutate(Species=str_replace(Species, "Agaricus cupressophilus", "Agaricus cupressicola cluster"))%>%
  
  # cluster 198
  mutate(Species=str_replace(Species, "Agaricus endoxanthus", "Agaricus endoxanthus cluster"))%>%
  mutate(Species=str_replace(Species, "Agaricus xanthosarcus", "Agaricus endoxanthus cluster"))%>%
  
  # cluster 232
  mutate(Species=str_replace(Species, "Agaricus lanipedisimilis", "Agaricus lanipedisimilis cluster"))%>%
  mutate(Species=str_replace(Species, "Agaricus planipileus", "Agaricus lanipedisimilis cluster"))%>%
  
  #####################################
######### Apiospora species ####
#####################################

# cluster 178
mutate(Species=str_replace(Species, "Apiospora gaoyouensis", "Apiospora gaoyouensis cluster"))%>%
  mutate(Species=str_replace(Species, "Arthrinium marii", "Apiospora gaoyouensis cluster"))%>%
  
  # cluster 180
  mutate(Species=str_replace(Species, "Apiospora hydei", "Apiospora hydei cluster"))%>%
  mutate(Species=str_replace(Species, "Arthrinium aureum", "Apiospora hydei cluster"))%>%
  
  #####################################
######### Aplosporella species ####
#####################################

# cluster 208
mutate(Species=str_replace(Species, "Aplosporella artocarpi", "Aplosporella artocarpi cluster"))%>%
  mutate(Species=str_replace(Species, "Aplosporella thailandica", "Aplosporella artocarpi cluster"))%>%
  
  #####################################
######### Arthroderma species ####
#####################################

# cluster 303
mutate(Species=str_replace(Species, "Arthroderma gertleri", "Arthroderma gertleri cluster"))%>%
  mutate(Species=str_replace(Species, "Arthroderma gloriae", "Arthroderma gertleri cluster"))%>%
  
  #####################################
######### Ascosphaera species ####
#####################################

# cluster 351
mutate(Species=str_replace(Species, "Ascosphaera atra", "Ascosphaera atra cluster"))%>%
  mutate(Species=str_replace(Species, "Ascosphaera duoformis", "Ascosphaera atra cluster"))%>%
  
  #####################################
######### Antherospora species ####
#####################################

# cluster 502
mutate(Species=str_replace(Species, "Antherospora hortensis", "Antherospora hortensis cluster"))%>%
  mutate(Species=str_replace(Species, "Antherospora muscari-botryoidis", "Antherospora hortensis cluster"))%>%
  mutate(Species=str_replace(Species, "Antherospora vaillantii", "Antherospora hortensis cluster"))%>%
  
  # cluster 503
  mutate(Species=str_replace(Species, "Antherospora scillae", "Antherospora scillae cluster"))%>%
  mutate(Species=str_replace(Species, "Antherospora vindobonensis", "Antherospora scillae cluster"))%>%
  
  #####################################
######### Apiotrichum species ####
#####################################

# cluster 552
mutate(Species=str_replace(Species, "Apiotrichum dehoogii", "Apiotrichum dehoogii cluster"))%>%
  mutate(Species=str_replace(Species, "Apiotrichum sporotrichoides", "Apiotrichum dehoogii cluster"))%>%
  
  #####################################
######### Armillaria species ####
#####################################

# cluster 564
mutate(Species=str_replace(Species, "Armillaria fumosa", "Armillaria fumosa cluster"))%>%
  mutate(Species=str_replace(Species, "Armillaria pallidula", "Armillaria fumosa cluster"))%>%
  
  #####################################
######### Caloplaca species ####
#####################################

# cluster 25
mutate(Species=str_replace(Species, "Caloplaca confusa", "Caloplaca confusa cluster"))%>%
  mutate(Species=str_replace(Species, "Caloplaca microthallina", "Caloplaca confusa cluster"))%>%
  
  #####################################
######### Calvitimela species ####
#####################################

# cluster 91
mutate(Species=str_replace(Species, "Calvitimela armeniaca", "Calvitimela armeniaca cluster"))%>%
  mutate(Species=str_replace(Species, "Calvitimela melaleuca", "Calvitimela armeniaca cluster"))%>%
  
  #####################################
######### Capnodium species ####
#####################################

# cluster 149
mutate(Species=str_replace(Species, "Capnodium alfenasii", "Capnodium alfenasii cluster"))%>%
  mutate(Species=str_replace(Species, "Capnodium blackwelliae", "Capnodium alfenasii cluster"))%>%
  
  #####################################
######### Celoporthe species ####
#####################################

# cluster 200
mutate(Species=str_replace(Species, "Celoporthe borbonica", "Celoporthe borbonica cluster"))%>%
  mutate(Species=str_replace(Species, "Celoporthe cerciana", "Celoporthe borbonica cluster"))%>%
  mutate(Species=str_replace(Species, "Celoporthe eucalypti", "Celoporthe borbonica cluster"))%>%
  
  #####################################
######### Cephalotrichum species ####
#####################################

# cluster 215
mutate(Species=str_replace(Species, "Cephalotrichum guizhouense", "Cephalotrichum guizhouense cluster"))%>%
  mutate(Species=str_replace(Species, "Cephalotrichum lignatile", "Cephalotrichum guizhouense cluster"))%>%
  
  #####################################
######### Ceratocystiopsis species ##
#####################################

# cluster 221
mutate(Species=str_replace(Species, "Ceratocystiopsis collifera", "Ceratocystiopsis collifera cluster"))%>%
  mutate(Species=str_replace(Species, "Ceratocystiopsis ranaculosa", "Ceratocystiopsis collifera cluster"))%>%
  
  # cluster 235
  mutate(Species=str_replace(Species, "Ceratocystis caryae", "Ceratocystis caryae cluster"))%>%
  mutate(Species=str_replace(Species, "Ceratocystis smalleyi", "Ceratocystis caryae cluster"))%>%
  
  # cluster 243
  mutate(Species=str_replace(Species, "Ceratocystis ecuadoriana", "Ceratocystis ecuadoriana cluster"))%>%
  mutate(Species=str_replace(Species, "Ceratocystis neglecta", "Ceratocystis ecuadoriana cluster"))%>%
  
  # cluster 247
  mutate(Species=str_replace(Species, "Ceratocystis fimbriatomima", "Ceratocystis fimbriatomima cluster"))%>%
  mutate(Species=str_replace(Species, "Ceratocystis manginecans", "Ceratocystis fimbriatomima cluster"))%>%
  
  # cluster 248
  mutate(Species=str_replace(Species, "Ceratocystis fujiensis", "Ceratocystis fujiensis cluster"))%>%
  mutate(Species=str_replace(Species, "Ceratocystis laricicola", "Ceratocystis fujiensis cluster"))%>%
  mutate(Species=str_replace(Species, "Ceratocystis polonica", "Ceratocystis fujiensis cluster"))%>%
  
  # cluster 249
  mutate(Species=str_replace(Species, "Ceratocystis huliohia", "Ceratocystis huliohia cluster"))%>%
  mutate(Species=str_replace(Species, "Ceratocystis uchidae", "Ceratocystis huliohia cluster"))%>%
  
  #####################################
######### Chaetomium species ####
#####################################

# cluster 342
mutate(Species=str_replace(Species, "Chaetomium cervicicola", "Chaetomium cervicicola cluster"))%>%
  mutate(Species=str_replace(Species, "Chaetomium grande", "Chaetomium cervicicola cluster"))%>%
  
  # cluster 373
  mutate(Species=str_replace(Species, "Chaetothyrina guttulata", "Chaetothyrina guttulata cluster"))%>%
  mutate(Species=str_replace(Species, "Chaetothyrina spondiadis", "Chaetothyrina guttulata cluster"))%>%
  
  #####################################
######### Cladonia species ####
#####################################

# cluster 522
mutate(Species=str_replace(Species, "Cladonia digitata", "Cladonia digitata cluster"))%>%
  mutate(Species=str_replace(Species, "Cladonia sulphurina", "Cladonia digitata cluster"))%>%
  
  # cluster 538
  mutate(Species=str_replace(Species, "Cladonia isabellina", "Cladonia isabellina cluster"))%>%
  mutate(Species=str_replace(Species, "Cladonia verticillata", "Cladonia isabellina cluster"))%>%
  
  #####################################
######### ChaetomClarireediaium species ####
#####################################

# cluster 640
mutate(Species=str_replace(Species, "Clarireedia jacksonii", "Clarireedia jacksonii cluster"))%>%
  mutate(Species=str_replace(Species, "Clarireedia paspali", "Clarireedia jacksonii cluster"))%>%
  
  # cluster 662
  mutate(Species=str_replace(Species, "Claviceps setariicola", "Claviceps setariicola cluster"))%>%
  mutate(Species=str_replace(Species, "Claviceps tenuispora", "Claviceps setariicola cluster"))%>%
  
  #####################################
######### Chroogomphus species ####
#####################################

# cluster 820
mutate(Species=str_replace(Species, "Chroogomphus fulmineus", "Chroogomphus fulmineus cluster"))%>%
  mutate(Species=str_replace(Species, "Chroogomphus subfulmineus", "Chroogomphus fulmineus cluster"))%>%
  
  
  # cluster 83
  mutate(Species=str_replace(Species, "Collophorina badensis", "Collophorina badensis cluster"))%>%
  mutate(Species=str_replace(Species, "Collophorina germanica", "Collophorina badensis cluster"))%>%
  
  # cluster 101
  mutate(Species=str_replace(Species, "Coniella pseudokoreana", "Coniella pseudokoreana cluster"))%>%
  mutate(Species=str_replace(Species, "Coniella quercicola", "Coniella pseudokoreana cluster"))%>%
  
  # cluster 103
  mutate(Species=str_replace(Species, "Coniella vitis", "Coniella vitis cluster"))%>%
  mutate(Species=str_replace(Species, "Coniella diplodiella", "Coniella vitis cluster"))%>%
  
  # cluster 108
  mutate(Species=str_replace(Species, "Coniocessia maxima", "Coniocessia maxima cluster"))%>%
  mutate(Species=str_replace(Species, "Coniocessia minima", "Coniocessia maxima cluster"))%>%
  
  # cluster 154
  mutate(Species=str_replace(Species, "Cordyceps blackwelliae", "Cordyceps blackwelliae cluster"))%>%
  mutate(Species=str_replace(Species, "Cordyceps cf.", "Cordyceps blackwelliae cluster"))%>%
  
  # cluster 184
  mutate(Species=str_replace(Species, "Corynascus novoguineensis", "Corynascus novoguineensis cluster"))%>%
  mutate(Species=str_replace(Species, "Corynascus verrucosus", "Corynascus novoguineensis cluster"))%>%
  
  # cluster 204
  mutate(Species=str_replace(Species, "Coryneum sinense", "Coryneum sinense cluster"))%>%
  mutate(Species=str_replace(Species, "Coryneum suttonii", "Coryneum sinense cluster"))%>%
  
  # cluster 234
  mutate(Species=str_replace(Species, "Cryphonectria quercicola", "Cryphonectria quercicola cluster"))%>%
  mutate(Species=str_replace(Species, "Cryphonectria quercus", "Cryphonectria quercicola cluster"))%>%
  
  # cluster 258
  mutate(Species=str_replace(Species, "Cryptosporella alni-rubrae", "Cryptosporella alni-rubrae cluster"))%>%
  mutate(Species=str_replace(Species, "Cryptosporella jaklitschii", "Cryptosporella alni-rubrae cluster"))%>%
  
  # cluster 329
  mutate(Species=str_replace(Species, "Cylindrocladiella pseudoparva", "Cylindrocladiella pseudoparva cluster"))%>%
  mutate(Species=str_replace(Species, "Cylindrocladiella stellenboschensis", "Cylindrocladiella pseudoparva cluster"))%>%
  
  # cluster 642
  mutate(Species=str_replace(Species, "Coriolopsis dendriformis", "Coriolopsis dendriformis cluster"))%>%
  mutate(Species=str_replace(Species, "Coriolopsis retropicta", "Coriolopsis dendriformis cluster"))%>%
  
  # cluster 734
  mutate(Species=str_replace(Species, "Cryptococcus aciditolerans", "Cryptococcus aciditolerans cluster"))%>%
  mutate(Species=str_replace(Species, "Cryptococcus metallitolerans", "Cryptococcus aciditolerans cluster"))%>%
  
  # cluster 788
  mutate(Species=str_replace(Species, "Cyanosporus caesius", "Cyanosporus caesius cluster"))%>%
  mutate(Species=str_replace(Species, "Cyanosporus submicroporus", "Cyanosporus caesius cluster"))%>%
  
  # cluster 790
  mutate(Species=str_replace(Species, "Cyanosporus fusiformis", "Cyanosporus fusiformis cluster"))%>%
  mutate(Species=str_replace(Species, "Cyanosporus glauca", "Cyanosporus fusiformis cluster"))%>%
  
  
  
  # cluster 55
  mutate(Species=str_replace(Species, "Debaryomyces robertsiae", "Debaryomyces robertsiae cluster"))%>%
  mutate(Species=str_replace(Species, "Debaryomyces udenii", "Debaryomyces robertsiae cluster"))%>%
  
  # cluster 79
  mutate(Species=str_replace(Species, "Dendrostoma castaneae", "Dendrostoma castaneae cluster"))%>%
  mutate(Species=str_replace(Species, "Dendrostoma luteum", "Dendrostoma castaneae cluster"))%>%
  
  # cluster 93
  mutate(Species=str_replace(Species, "Deniquelata barringtoniae", "Deniquelata barringtoniae cluster"))%>%
  mutate(Species=str_replace(Species, "Deniquelata hypolithi", "Deniquelata barringtoniae cluster"))%>%
  
  # cluster 500
  mutate(Species=str_replace(Species, "Dentipellis parmastoi", "Dentipellis parmastoi cluster"))%>%
  mutate(Species=str_replace(Species, "Dentipellis tropicalis", "Dentipellis parmastoi cluster"))%>%
  
  # cluster 522
  mutate(Species=str_replace(Species, "Derxomyces boekhoutii", "Derxomyces boekhoutii cluster"))%>%
  mutate(Species=str_replace(Species, "Derxomyces pseudoboekhoutii", "Derxomyces boekhoutii cluster"))%>%
  
  # cluster 649
  mutate(Species=str_replace(Species, "Entoloma aurorae-borealis", "Entoloma aurorae-borealis cluster"))%>%
  mutate(Species=str_replace(Species, "Entoloma serpens", "Entoloma aurorae-borealis cluster"))%>%
  
  # cluster 709
  mutate(Species=str_replace(Species, "Entoloma hirtipes", "Entoloma hirtipes cluster"))%>%
  mutate(Species=str_replace(Species, "Entoloma tibiicystidiatum", "Entoloma hirtipes cluster"))%>%
  
  # cluster 821
  mutate(Species=str_replace(Species, "Entyloma arnicale", "Entyloma arnicale cluster"))%>%
  mutate(Species=str_replace(Species, "Entyloma cosmi", "Entyloma arnicale cluster"))%>%
  mutate(Species=str_replace(Species, "Entyloma dahliae", "Entyloma arnicale cluster"))%>%
  
  # cluster 827
  mutate(Species=str_replace(Species, "Entyloma jolantae", "Entyloma jolantae cluster"))%>%
  mutate(Species=str_replace(Species, "Entyloma savchenkoi", "Entyloma jolantae cluster"))%>%
  
  # cluster 74
  mutate(Species=str_replace(Species, "Exserohilum curvatum", "Exserohilum curvatum cluster"))%>%
  mutate(Species=str_replace(Species, "Exserohilum holmii", "Exserohilum curvatum cluster"))%>%
  
  # cluster 175
  mutate(Species=str_replace(Species, "Fuscopannaria protensa", "Fuscopannaria protensa cluster"))%>%
  mutate(Species=str_replace(Species, "Fuscopannaria sorediata", "Fuscopannaria protensa cluster"))%>%
  
  # cluster 202
  mutate(Species=str_replace(Species, "Gaeumannomyces fusiformis", "Gaeumannomyces fusiformis cluster"))%>%
  mutate(Species=str_replace(Species, "Gaeumannomyces oryzinus", "Gaeumannomyces fusiformis cluster"))%>%
  
  # cluster 270
  mutate(Species=str_replace(Species, "Geopyxis majalis", "Geopyxis majalis cluster"))%>%
  mutate(Species=str_replace(Species, "Geopyxis rehmii", "Geopyxis majalis cluster"))%>%
  
  # cluster 277
  mutate(Species=str_replace(Species, "Geosmithia luteobrunnea", "Geosmithia luteobrunnea cluster"))%>%
  mutate(Species=str_replace(Species, "Geosmithia radiata", "Geosmithia luteobrunnea cluster"))%>%
  
  # cluster 305
  mutate(Species=str_replace(Species, "Gliocephalotrichum bacillisporum", "Gliocephalotrichum bacillisporum cluster"))%>%
  mutate(Species=str_replace(Species, "Gliocephalotrichum bulbilium", "Gliocephalotrichum bacillisporum cluster"))%>%
  
  # cluster 381
  mutate(Species=str_replace(Species, "Graphilbum furuicola", "Graphilbum furuicola cluster"))%>%
  mutate(Species=str_replace(Species, "Graphilbum sexdentatum", "Graphilbum furuicola cluster"))%>%
  
  # cluster 453
  mutate(Species=str_replace(Species, "Gyromitra esculenta", "Gyromitra esculenta cluster"))%>%
  mutate(Species=str_replace(Species, "Gyromitra tasmanica", "Gyromitra esculenta cluster"))%>%
  
  # cluster 609
  mutate(Species=str_replace(Species, "Fomitiporia maxonii", "Fomitiporia maxonii cluster"))%>%
  mutate(Species=str_replace(Species, "Fomitiporia sonorae", "Fomitiporia maxonii cluster"))%>%
  
  # cluster 667
  mutate(Species=str_replace(Species, "Fuscoporia ramulicola", "Fuscoporia ramulicola cluster"))%>%
  mutate(Species=str_replace(Species, "Fuscoporia subferrea", "Fuscoporia ramulicola cluster"))%>%
  
  # cluster 712
  mutate(Species=str_replace(Species, "Ganoderma curtisii", "Ganoderma curtisii cluster"))%>%
  mutate(Species=str_replace(Species, "Ganoderma meredithae", "Ganoderma curtisii cluster"))%>%
  
  # cluster 798
  mutate(Species=str_replace(Species, "Geastrum smardae", "Geastrum smardae cluster"))%>%
  mutate(Species=str_replace(Species, "Geastrum triplex", "Geastrum smardae cluster"))%>%
  
  # cluster 818
  mutate(Species=str_replace(Species, "Gerronema nemorale", "Gerronema nemorale cluster"))%>%
  mutate(Species=str_replace(Species, "Gerronema subclavatum", "Gerronema nemorale cluster"))%>%
  
  # cluster 1032
  mutate(Species=str_replace(Species, "Gongronella lacrispora", "Gongronella lacrispora cluster"))%>%
  mutate(Species=str_replace(Species, "Gongronella orasabula", "Gongronella lacrispora cluster"))%>%
  
  # cluster 10
  mutate(Species=str_replace(Species, "Hamigera avellanea", "Hamigera avellanea cluster"))%>%
  mutate(Species=str_replace(Species, "Hamigera paravellanea", "Hamigera avellanea cluster"))%>%
  
  # cluster 42
  mutate(Species=str_replace(Species, "Harknessia pilularis", "Harknessia pilularis cluster"))%>%
  mutate(Species=str_replace(Species, "Harknessia rhabdosphaera", "Harknessia pilularis cluster"))%>%
  
  # cluster 97
  mutate(Species=str_replace(Species, "Helvella dovrensis", "Helvella dovrensis cluster"))%>%
  mutate(Species=str_replace(Species, "Helvella philonotis", "Helvella dovrensis cluster"))%>%
  
  # cluster 167
  mutate(Species=str_replace(Species, "Heterosphaeria linariae", "Heterosphaeria linariae cluster"))%>%
  mutate(Species=str_replace(Species, "Heterosphaeria patella", "Heterosphaeria linariae cluster"))%>%
  
  # cluster 195
  mutate(Species=str_replace(Species, "Holocryphia sp.", "Holocryphia sp. cluster"))%>%
  mutate(Species=str_replace(Species, "Holocryphia sp.", "Holocryphia sp. cluster"))%>%
  
  # cluster 302
  mutate(Species=str_replace(Species, "Hebeloma nanum", "Hebeloma nanum cluster"))%>%
  mutate(Species=str_replace(Species, "Hebeloma pungens", "Hebeloma nanum cluster"))%>%
  
  # cluster 339
  mutate(Species=str_replace(Species, "Hemistropharia albocrenulata", "Hemistropharia albocrenulata cluster"))%>%
  mutate(Species=str_replace(Species, "Hemistropharia subalbocrenulata", "Hemistropharia albocrenulata cluster"))%>%
  
  # cluster 378
  mutate(Species=str_replace(Species, "Hohenbuehelia bonii", "Hohenbuehelia bonii cluster"))%>%
  mutate(Species=str_replace(Species, "Hohenbuehelia leiospora", "Hohenbuehelia bonii cluster"))%>%
  
  # cluster 455
  mutate(Species=str_replace(Species, "Hydnum berkeleyanum", "Hydnum berkeleyanum cluster"))%>%
  mutate(Species=str_replace(Species, "Hydnum ventricosum", "Hydnum berkeleyanum cluster"))%>%
  
  # cluster 456
  mutate(Species=str_replace(Species, "Hydnum boreorepandum", "Hydnum boreorepandum cluster"))%>%
  mutate(Species=str_replace(Species, "Hydnum repandum", "Hydnum boreorepandum cluster"))%>%
  
  # cluster 351
  mutate(Species=str_replace(Species, "Hygrophorus morrisii", "Hygrophorus morrisii cluster"))%>%
  mutate(Species=str_replace(Species, "Hygrophorus odoratus", "Hygrophorus morrisii cluster"))%>%
  
  # cluster 104
  mutate(Species=str_replace(Species, "Kionochaeta microspora", "Kionochaeta microspora cluster"))%>%
  mutate(Species=str_replace(Species, "Kionochaeta spissa", "Kionochaeta microspora cluster"))%>%
  
  # cluster 233
  mutate(Species=str_replace(Species, "Kockovaella calophylli", "Kockovaella calophylli cluster"))%>%
  mutate(Species=str_replace(Species, "Kockovaella vietnamensis", "Kockovaella calophylli cluster"))%>%
  
  # cluster 332
  mutate(Species=str_replace(Species, "Lactarius aspideoides", "Lactarius aspideoides cluster"))%>%
  mutate(Species=str_replace(Species, "Lactarius salicis-reticulatae", "Lactarius aspideoides cluster"))%>%
  
  # cluster 378
  mutate(Species=str_replace(Species, "Lactarius giennensis", "Lactarius giennensis cluster"))%>%
  mutate(Species=str_replace(Species, "Lactarius subgiennensis", "Lactarius giennensis cluster"))%>%
  
  # cluster 436
  mutate(Species=str_replace(Species, "Lactarius populicola", "Lactarius populicola cluster"))%>%
  mutate(Species=str_replace(Species, "Lactarius soehneri", "Lactarius populicola cluster"))%>%
  
  # cluster 143
  mutate(Species=str_replace(Species, "Lecanosticta acicola", "Lecanosticta acicola cluster"))%>%
  mutate(Species=str_replace(Species, "Lecanosticta brevispora", "Lecanosticta acicola cluster"))%>%
  mutate(Species=str_replace(Species, "Lecanosticta guatemalensis", "Lecanosticta acicola cluster"))%>%
  mutate(Species=str_replace(Species, "Lecanosticta jani", "Lecanosticta acicola cluster"))%>%
  mutate(Species=str_replace(Species, "Lecanosticta longispora", "Lecanosticta acicola cluster"))%>%
  
  # cluster 183
  mutate(Species=str_replace(Species, "Lectera longa", "Lectera longa cluster"))%>%
  mutate(Species=str_replace(Species, "Lectera nordwiniana", "Lectera longa cluster"))%>%
  
  # cluster 201
  mutate(Species=str_replace(Species, "Lentithecium clioninum", "Lentithecium clioninum cluster"))%>%
  mutate(Species=str_replace(Species, "Lentithecium pseudoclioninum", "Lentithecium clioninum cluster"))%>%
  
  # cluster 223
  mutate(Species=str_replace(Species, "Lepraria finkii", "Lepraria finkii cluster"))%>%
  mutate(Species=str_replace(Species, "Lepraria lobificans", "Lepraria finkii cluster"))%>%
  
  # cluster 379
  mutate(Species=str_replace(Species, "Laetiporus ailaoshanensis", "Laetiporus ailaoshanensis cluster"))%>%
  mutate(Species=str_replace(Species, "Laetiporus zonatus", "Laetiporus ailaoshanensis cluster"))%>%
  
  # cluster 434
  mutate(Species=str_replace(Species, "Lentinula boryana", "Lentinula boryana cluster"))%>%
  mutate(Species=str_replace(Species, "Lentinula detonsa", "Lentinula boryana cluster"))%>%
  
  # cluster 481
  mutate(Species=str_replace(Species, "Lepiota nigrescentipes", "Lepiota nigrescentipes cluster"))%>%
  mutate(Species=str_replace(Species, "Lepiota parvannulata", "Lepiota nigrescentipes cluster"))%>%
  
  # cluster 44
  mutate(Species=str_replace(Species, "Lindgomyces angustiascus", "Lindgomyces angustiascus cluster"))%>%
  mutate(Species=str_replace(Species, "Lindgomyces lemonweirensis", "Lindgomyces angustiascus cluster"))%>%
  
  # cluster 77
  mutate(Species=str_replace(Species, "Lobaria hengduanensis", "Lobaria hengduanensis cluster"))%>%
  mutate(Species=str_replace(Species, "Lobaria irrugulosa", "Lobaria hengduanensis cluster"))%>%
  
  # cluster 101
  mutate(Species=str_replace(Species, "Lobariella flynniana", "Lobariella flynniana cluster"))%>%
  mutate(Species=str_replace(Species, "Lobariella sandwicensis", "Lobariella flynniana cluster"))%>%
  
  # cluster 139
  mutate(Species=str_replace(Species, "Lophiostoma caryophyllacearum", "Lophiostoma caryophyllacearum cluster"))%>%
  mutate(Species=str_replace(Species, "Lophiostoma montanae", "Lophiostoma caryophyllacearum cluster"))%>%
  
  # cluster 197
  mutate(Species=str_replace(Species, "Macroconia bulbipes", "Macroconia bulbipes cluster"))%>%
  mutate(Species=str_replace(Species, "Macroconia papilionacearum", "Macroconia bulbipes cluster"))%>%
  
  # cluster 314
  mutate(Species=str_replace(Species, "Leucopaxillus gentianeus", "Leucopaxillus gentianeus cluster"))%>%
  mutate(Species=str_replace(Species, "Leucopaxillus laterarius", "Leucopaxillus gentianeus cluster"))%>%
  
  # cluster 342
  mutate(Species=str_replace(Species, "Lopharia cinerascens", "Lopharia cinerascens cluster"))%>%
  mutate(Species=str_replace(Species, "Lopharia mirabilis", "Lopharia cinerascens cluster"))%>%
  
  # cluster 362
  mutate(Species=str_replace(Species, "Lyomyces fimbriatus", "Lyomyces fimbriatus cluster"))%>%
  mutate(Species=str_replace(Species, "Lyomyces leptocystidiatus", "Lyomyces fimbriatus cluster"))%>%
  
  # cluster 385
  mutate(Species=str_replace(Species, "Macalpinomyces langdonii", "Macalpinomyces langdonii cluster"))%>%
  mutate(Species=str_replace(Species, "Macalpinomyces terrae-australis", "Macalpinomyces langdonii cluster"))%>%
  
  # cluster 432
  mutate(Species=str_replace(Species, "Luteoperenniporia mopanshanensis", "Luteoperenniporia mopanshanensis cluster"))%>%
  mutate(Species=str_replace(Species, "Luteoperenniporia yinggelingensis", "Luteoperenniporia mopanshanensis cluster"))%>%
  
  # cluster 3
  mutate(Species=str_replace(Species, "Maublancomyces korfii", "Maublancomyces korfii cluster"))%>%
  mutate(Species=str_replace(Species, "Maublancomyces ticinianus", "Maublancomyces korfii cluster"))%>%
  
  # cluster 16
  mutate(Species=str_replace(Species, "Mariannaea macrochlamydospora", "Mariannaea macrochlamydospora cluster"))%>%
  mutate(Species=str_replace(Species, "Mariannaea submersa", "Mariannaea macrochlamydospora cluster"))%>%
  
  # cluster 42
  mutate(Species=str_replace(Species, "Massaria mediterranea", "Massaria mediterranea cluster"))%>%
  mutate(Species=str_replace(Species, "Massaria platanoidea", "Massaria mediterranea cluster"))%>%
  
  # cluster 95
  mutate(Species=str_replace(Species, "Melanelixia ahtii", "Melanelixia ahtii cluster"))%>%
  mutate(Species=str_replace(Species, "Melanelixia glabroides", "Melanelixia ahtii cluster"))%>%
  
  # cluster 237
  mutate(Species=str_replace(Species, "Microascus alveolaris", "Microascus alveolaris cluster"))%>%
  mutate(Species=str_replace(Species, "Microascus campaniformis", "Microascus alveolaris cluster"))%>%
  
  # cluster 245
  mutate(Species=str_replace(Species, "Microascus cleistocarpus", "Microascus cleistocarpus cluster"))%>%
  mutate(Species=str_replace(Species, "Microascus hyalinus", "Microascus cleistocarpus cluster"))%>%
  
  # cluster 375
  mutate(Species=str_replace(Species, "Megacollybia rodmani", "Megacollybia rodmani cluster"))%>%
  mutate(Species=str_replace(Species, "Megacollybia subfurfuracea", "Megacollybia rodmani cluster"))%>%
  
  # cluster 381
  mutate(Species=str_replace(Species, "Megalocystidium pellitum", "Megalocystidium pellitum cluster"))%>%
  mutate(Species=str_replace(Species, "Megalocystidium perticatum", "Megalocystidium pellitum cluster"))%>%
  
  # cluster 420
  mutate(Species=str_replace(Species, "Melanoleuca exscissa", "Melanoleuca exscissa cluster"))%>%
  mutate(Species=str_replace(Species, "Melanoleuca substrictipes", "Melanoleuca exscissa cluster"))%>%
  
  # cluster 438
  mutate(Species=str_replace(Species, "Melanoleuca porphyropoda", "Melanoleuca porphyropoda cluster"))%>%
  mutate(Species=str_replace(Species, "Melanoleuca stepposa", "Melanoleuca porphyropoda cluster"))%>%
  mutate(Species=str_replace(Species, "Melanoleuca stridula", "Melanoleuca porphyropoda cluster"))%>%
  
  # cluster 465
  mutate(Species=str_replace(Species, "Microbotryum coronariae", "Microbotryum coronariae cluster"))%>%
  mutate(Species=str_replace(Species, "Microbotryum silenes-inflatae", "Microbotryum coronariae cluster"))%>%
  
  # cluster 470
  mutate(Species=str_replace(Species, "Microbotryum lychnidis-dioicae", "Microbotryum lychnidis-dioicae cluster"))%>%
  mutate(Species=str_replace(Species, "Microbotryum minuartiae", "Microbotryum lychnidis-dioicae cluster"))%>%
  
  # cluster 46
  mutate(Species=str_replace(Species, "Neofusicladium eucalypti", "Neofusicladium eucalypti cluster"))%>%
  mutate(Species=str_replace(Species, "Neofusicladium eucalypticola", "Neofusicladium eucalypti cluster"))%>%
  
  # cluster 121
  mutate(Species=str_replace(Species, "Neosartorya quadricincta", "Neosartorya quadricincta cluster"))%>%
  mutate(Species=str_replace(Species, "Neosartorya stramenia", "Neosartorya quadricincta cluster"))%>%
  
  # cluster 186
  mutate(Species=str_replace(Species, "Niebla eburnea", "Niebla eburnea cluster"))%>%
  mutate(Species=str_replace(Species, "Niebla homalea", "Niebla eburnea cluster"))%>%
  
  # cluster 191
  mutate(Species=str_replace(Species, "Niebla podetiaforma", "Niebla podetiaforma cluster"))%>%
  mutate(Species=str_replace(Species, "Niebla sorocarpia", "Niebla podetiaforma cluster"))%>%
  
  # cluster 288
  mutate(Species=str_replace(Species, "Ogataea ganodermae", "Ogataea ganodermae cluster"))%>%
  mutate(Species=str_replace(Species, "Ogataea neopini", "Ogataea ganodermae cluster"))%>%
  
  # cluster 390
  mutate(Species=str_replace(Species, "Nigrohirschioporus confertus", "Nigrohirschioporus confertus cluster"))%>%
  mutate(Species=str_replace(Species, "Nigrohirschioporus sector", "Nigrohirschioporus confertus cluster"))%>%
  
  # cluster 194
  mutate(Species=str_replace(Species, "Pachyphlodes conglomerata", "Pachyphlodes conglomerata cluster"))%>%
  mutate(Species=str_replace(Species, "Pachyphlodes melanoxantha", "Pachyphlodes conglomerata cluster"))%>%
  
  # cluster 294
  mutate(Species=str_replace(Species, "Paralulworthia gigaspora", "Paralulworthia gigaspora cluster"))%>%
  mutate(Species=str_replace(Species, "Paralulworthia mediterranea", "Paralulworthia gigaspora cluster"))%>%
  
  # cluster 328
  mutate(Species=str_replace(Species, "Paraphoma pye", "Paraphoma pye cluster"))%>%
  mutate(Species=str_replace(Species, "Paraphoma vinacea", "Paraphoma pye cluster"))%>%
  
  # cluster 30
  mutate(Species=str_replace(Species, "Phialophora chinensis", "Phialophora chinensis cluster"))%>%
  mutate(Species=str_replace(Species, "Phialophora verrucosa", "Phialophora chinensis cluster"))%>%
  
  # cluster 92
  mutate(Species=str_replace(Species, "Phyllactinia antarctica", "Phyllactinia antarctica cluster"))%>%
  mutate(Species=str_replace(Species, "Phyllactinia parietariae", "Phyllactinia antarctica cluster"))%>%
  
  # cluster 43
  mutate(Species=str_replace(Species, "Plagiostoma barriae", "Plagiostoma barriae cluster"))%>%
  mutate(Species=str_replace(Species, "Plagiostoma dilatatum", "Plagiostoma barriae cluster"))%>%
  
  # cluster 72
  mutate(Species=str_replace(Species, "Plenodomus agnitus", "Plenodomus agnitus cluster"))%>%
  mutate(Species=str_replace(Species, "Plenodomus fallaciosus", "Plenodomus agnitus cluster"))%>%
  
  # cluster 101
  mutate(Species=str_replace(Species, "Pleopunctum ellipsoideum", "Pleopunctum ellipsoideum cluster"))%>%
  mutate(Species=str_replace(Species, "Pleopunctum pseudoellipsoideum", "Pleopunctum ellipsoideum cluster"))%>%
  
  # cluster 112
  mutate(Species=str_replace(Species, "Pleurophoma italica", "Pleurophoma italica cluster"))%>%
  mutate(Species=str_replace(Species, "Pleurophoma ossicola", "Pleurophoma italica cluster"))%>%
  
  # cluster 317
  mutate(Species=str_replace(Species, "Pluteus alniphilus", "Pluteus alniphilus cluster"))%>%
  mutate(Species=str_replace(Species, "Pluteus spegazzinianus", "Pluteus alniphilus cluster"))%>%
  
  # cluster 389
  mutate(Species=str_replace(Species, "Podoscypha bolleana", "Podoscypha bolleana cluster"))%>%
  mutate(Species=str_replace(Species, "Podoscypha simulans", "Podoscypha bolleana cluster"))%>%
  
  # cluster 433
  mutate(Species=str_replace(Species, "Porodaedalea himalayensis", "Porodaedalea himalayensis cluster"))%>%
  mutate(Species=str_replace(Species, "Porodaedalea microsperma", "Porodaedalea himalayensis cluster"))%>%
  
  # cluster 50
  mutate(Species=str_replace(Species, "Pseudocosmospora eutypellae", "Pseudocosmospora eutypellae cluster"))%>%
  mutate(Species=str_replace(Species, "Pseudocosmospora rogersonii", "Pseudocosmospora eutypellae cluster"))%>%
  
  # cluster 65
  mutate(Species=str_replace(Species, "Pseudocyphellaria epiflavoides", "Pseudocyphellaria epiflavoides cluster"))%>%
  mutate(Species=str_replace(Species, "Pseudocyphellaria hawaiiensis", "Pseudocyphellaria epiflavoides cluster"))%>%
  
  # cluster 235
  mutate(Species=str_replace(Species, "Pyrenochaetopsis kuksensis", "Pyrenochaetopsis kuksensis cluster"))%>%
  mutate(Species=str_replace(Species, "Pyrenochaetopsis leptospora", "Pyrenochaetopsis kuksensis cluster"))%>%
  
  # cluster 83
  mutate(Species=str_replace(Species, "Ramalina crispatula", "Ramalina crispatula cluster"))%>%
  mutate(Species=str_replace(Species, "Ramalina webbii", "Ramalina crispatula cluster"))%>%
  
  # cluster 167
  mutate(Species=str_replace(Species, "Readeriella callista", "Readeriella callista cluster"))%>%
  mutate(Species=str_replace(Species, "Readeriella pseudocallista", "Readeriella callista cluster"))%>%
  
  # cluster 186
  mutate(Species=str_replace(Species, "Requienella fraxini", "Requienella fraxini cluster"))%>%
  mutate(Species=str_replace(Species, "Requienella seminuda", "Requienella fraxini cluster"))%>%
  
  # cluster 441
  mutate(Species=str_replace(Species, "Rhodotorula babjevae", "Rhodotorula babjevae cluster"))%>%
  mutate(Species=str_replace(Species, "Rhodotorula graminis", "Rhodotorula babjevae cluster"))%>%
  
  # cluster 3
  mutate(Species=str_replace(Species, "Rhynchosporium commune", "Rhynchosporium commune cluster"))%>%
  mutate(Species=str_replace(Species, "Rhynchosporium secalis", "Rhynchosporium commune cluster"))%>%
  
  # cluster 12
  mutate(Species=str_replace(Species, "Rhytidhysteron rufulum", "Rhytidhysteron rufulum cluster"))%>%
  mutate(Species=str_replace(Species, "Rhytidhysteron thailandicum", "Rhytidhysteron rufulum cluster"))%>%
  
  # cluster 137
  mutate(Species=str_replace(Species, "Saccharomyces bayanus", "Saccharomyces bayanus cluster"))%>%
  mutate(Species=str_replace(Species, "Saccharomyces kudriavzevii", "Saccharomyces bayanus cluster"))%>%
  
  # cluster 193
  mutate(Species=str_replace(Species, "Sarcostroma diversiseptatum", "Sarcostroma diversiseptatum cluster"))%>%
  mutate(Species=str_replace(Species, "Sarcostroma paragrevilleae", "Sarcostroma diversiseptatum cluster"))%>%
  
  # cluster 232
  mutate(Species=str_replace(Species, "Savoryella aquatica", "Savoryella aquatica cluster"))%>%
  mutate(Species=str_replace(Species, "Savoryella verrucosa", "Savoryella aquatica cluster"))%>%
  
  # cluster 278
  mutate(Species=str_replace(Species, "Rubroboletus haematinus", "Rubroboletus haematinus cluster"))%>%
  mutate(Species=str_replace(Species, "Rubroboletus pulcherrimus", "Rubroboletus haematinus cluster"))%>%
  
  # cluster 298
  mutate(Species=str_replace(Species, "Russula atrovirens", "Russula atrovirens cluster"))%>%
  mutate(Species=str_replace(Species, "Russula atroviridis", "Russula atrovirens cluster"))%>%
  
  # cluster 313
  mutate(Species=str_replace(Species, "Russula catillus", "Russula catillus cluster"))%>%
  mutate(Species=str_replace(Species, "Russula indocatillus", "Russula catillus cluster"))%>%
  
  # cluster 1
  mutate(Species=str_replace(Species, "Scheffersomyces lignicola", "Scheffersomyces lignicola cluster"))%>%
  mutate(Species=str_replace(Species, "Scheffersomyces queiroziae", "Scheffersomyces lignicola cluster"))%>%
  
  # cluster 232
  mutate(Species=str_replace(Species, "Sphaerulina pelargonii", "Sphaerulina pelargonii cluster"))%>%
  mutate(Species=str_replace(Species, "Sphaerulina westendorpii", "Sphaerulina pelargonii cluster"))%>%
  
  # cluster 27
  mutate(Species=str_replace(Species, "Stagonosporopsis actaeae", "Stagonosporopsis actaeae cluster"))%>%
  mutate(Species=str_replace(Species, "Stagonosporopsis oculi-hominis", "Stagonosporopsis actaeae cluster"))%>%
  mutate(Species=str_replace(Species, "Stagonosporopsis trachelii", "Stagonosporopsis actaeae cluster"))%>%
  
  # cluster 41
  mutate(Species=str_replace(Species, "Staphylotrichum acaciicola", "Staphylotrichum acaciicola cluster"))%>%
  mutate(Species=str_replace(Species, "Staphylotrichum koreanum", "Staphylotrichum acaciicola cluster"))%>%
  
  # cluster 120
  mutate(Species=str_replace(Species, "Sticta aff. peltigerella", "Sticta aff. peltigerella cluster"))%>%
  mutate(Species=str_replace(Species, "Sticta fuliginoides", "Sticta aff. peltigerella cluster"))%>%
  
  # cluster 307
  mutate(Species=str_replace(Species, "Synnemapestaloides juniperi", "Synnemapestaloides juniperi cluster"))%>%
  mutate(Species=str_replace(Species, "Synnemapestaloides rhododendri", "Synnemapestaloides juniperi cluster"))%>%
  
  # cluster 340
  mutate(Species=str_replace(Species, "Stephanospora flava", "Stephanospora flava cluster"))%>%
  mutate(Species=str_replace(Species, "Stephanospora occidentiaustralis", "Stephanospora flava cluster"))%>%
  
  # cluster 368
  mutate(Species=str_replace(Species, "Strobilomyces alpinus", "Strobilomyces alpinus cluster"))%>%
  mutate(Species=str_replace(Species, "Strobilomyces longistipitatus", "Strobilomyces alpinus cluster"))%>%
  
  # cluster 92
  mutate(Species=str_replace(Species, "Tubakia americana", "Tubakia americana cluster"))%>%
  mutate(Species=str_replace(Species, "Tubakia dryinoides", "Tubakia americana cluster"))%>%
  
  # cluster 228
  mutate(Species=str_replace(Species, "Tremella cinnabarina", "Tremella cinnabarina cluster"))%>%
  mutate(Species=str_replace(Species, "Tremella flava", "Tremella cinnabarina cluster"))%>%
  
  # cluster 312
  mutate(Species=str_replace(Species, "Tricholoma olivaceotinctum", "Tricholoma olivaceotinctum cluster"))%>%
  mutate(Species=str_replace(Species, "Tricholoma squarrulosum", "Tricholoma olivaceotinctum cluster"))%>%
  
  # cluster 355
  mutate(Species=str_replace(Species, "Trichosporon aquatile", "Trichosporon aquatile cluster"))%>%
  mutate(Species=str_replace(Species, "Trichosporon ovoides", "Trichosporon aquatile cluster"))%>%
  
  # cluster 399
  mutate(Species=str_replace(Species, "Tubaria serrulata", "Tubaria serrulata cluster"))%>%
  mutate(Species=str_replace(Species, "Tubaria similis", "Tubaria serrulata cluster"))%>%
  
  # cluster 481
  mutate(Species=str_replace(Species, "Typhrasa polycystis", "Typhrasa polycystis cluster"))%>%
  mutate(Species=str_replace(Species, "Typhrasa rugocephala", "Typhrasa polycystis cluster"))%>%
  
  # cluster 3
  mutate(Species=str_replace(Species, "Ulocladium chartarum", "Ulocladium chartarum cluster"))%>%
  mutate(Species=str_replace(Species, "Ulocladium cucurbitae", "Ulocladium chartarum cluster"))%>%
  
  # cluster 21
  mutate(Species=str_replace(Species, "Umbilicaria haplocarpa", "Umbilicaria haplocarpa cluster"))%>%
  mutate(Species=str_replace(Species, "Umbilicaria leprosa", "Umbilicaria haplocarpa cluster"))%>%
  
  # cluster 67
  mutate(Species=str_replace(Species, "Usnea barbata", "Usnea barbata cluster"))%>%
  mutate(Species=str_replace(Species, "Usnea filipendula", "Usnea barbata cluster"))%>%
  mutate(Species=str_replace(Species, "Usnea intermedia", "Usnea barbata cluster"))%>%
  
  # cluster 77
  mutate(Species=str_replace(Species, "Usnea cornuta", "Usnea cornuta cluster"))%>%
  mutate(Species=str_replace(Species, "Usnea subpectinata", "Usnea cornuta cluster"))%>%
  
  # cluster 135
  mutate(Species=str_replace(Species, "Uwebraunia australiensis", "Uwebraunia australiensis cluster"))%>%
  mutate(Species=str_replace(Species, "Uwebraunia musae", "Uwebraunia australiensis cluster"))%>%
  
  # cluster 324
  mutate(Species=str_replace(Species, "Wickerhamiella tropicalis", "Wickerhamiella tropicalis cluster"))%>%
  mutate(Species=str_replace(Species, "Wickerhamiella verensis", "Wickerhamiella tropicalis cluster"))%>%
  
  # cluster 458
  mutate(Species=str_replace(Species, "Wallemia canadensis", "Wallemia canadensis cluster"))%>%
  mutate(Species=str_replace(Species, "Wallemia tropicalis", "Wallemia canadensis cluster"))%>%
  
  # cluster 78
  mutate(Species=str_replace(Species, "Xanthoria digitata", "Xanthoria digitata cluster"))%>%
  mutate(Species=str_replace(Species, "Xanthoria papillifera", "Xanthoria digitata cluster"))%>%
  
  # cluster 340
  mutate(Species=str_replace(Species, "Xeromphalina brevipes", "Xeromphalina brevipes cluster"))%>%
  mutate(Species=str_replace(Species, "Xeromphalina utricularis", "Xeromphalina brevipes cluster"))%>%
  
  # cluster 352
  mutate(Species=str_replace(Species, "Xerula longipes", "Xerula longipes cluster"))%>%
  mutate(Species=str_replace(Species, "Xerula sinopudens", "Xerula longipes cluster"))


phyloseq::tax_table(mock_raw_others) <- as.matrix(tax.cleanM) ; phyloseq::tax_table(mock_raw_others)

###############################################################
# Merger mock_raw_others et mock_raw_APOGEE 
# (sans arbres pour éviter conflits de tips)
###############################################################

# Merger
mock_raw <- merge_phyloseq(mock_raw_others, mock_raw_APOGEE)

# Reconstruire un arbre aléatoire
library("ape")
random_tree <- rtree(ntaxa(mock_raw), rooted=TRUE, tip.label=taxa_names(mock_raw))
mock_raw <- phyloseq::phyloseq(otu_table(mock_raw), 
                                tax_table(mock_raw), 
                                sample_data(mock_raw),
                                random_tree)


##### save the un-rarefied object 

# Taxa not supported with at least 10 occurrences, were removed from the final "mock" object. 
mockR = prune_taxa(taxa_sums(mock_raw) > 09, mock_raw); mockR

saveRDS(mockR, "R_objects/mock_R.rds")


##### save the un-rarefied object 

mock_raw_rare <- rarefy_even_depth(mock_raw,
                                rngseed=1515,
                                sample.size=10000,
                                replace=FALSE)


# Taxa not supported with at least 10 occurrences, were removed from the final "mock" object. 
mock_raw_rare_F = prune_taxa(taxa_sums(mock_raw_rare) > 09, mock_raw_rare); mock_raw_rare_F


saveRDS(mock_raw_rare_F, "R_objects/mock_R_rare.rds")