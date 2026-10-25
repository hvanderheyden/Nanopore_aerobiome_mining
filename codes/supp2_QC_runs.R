
library(forcats)
library(tidyverse)

run_summary<-read.csv("Data/QC/runs_summary.txt", sep="\t")
head(run_summary)

run_summary$filt <- c("Raw")

run_final_summary<-read.csv("Data/QC/runs_summary_final.txt", sep="\t")
head(run_final_summary)

run_final_summary$filt <- c("Filtered")

Runs<-rbind(run_summary, run_final_summary)

library("plyr")
library("stringfish")

Runs$run_id<-sf_gsub(Runs$run_id, "be1e6bb5a7b953e819d788c46f30b62f081f2c86", "A", nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "04ebfd995e04463e2236f45b7eb19a4d01725e35", "B",  nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "e9a5112705c1b4eba6c171de29880da37664b528", "C", nthreads = 8);gc()

Runs$run_id<-sf_gsub(Runs$run_id, "5e268682ec04d421801f51aeb2f95f263bf93347", "D", nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "00c90becb47f13eb0d56c578377e24fb4f850b7c", "F",  nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "6f3e8dd7148c0f3f6b15f68cc3bbc83e36caadcd", "G", nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "a21874c067d5979e08f49a63e5075993bfe166d7", "G", nthreads = 8);gc()

Runs$run_id<-sf_gsub(Runs$run_id, "5aa25e7487617709b5f7665b63eb92c9fb012d8e", "H", nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "ab2090c144afc6c86cc56f886323e23648e7ce6a", "I",  nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "3cb61b505bda66acd48cc439907dfc34891990d9", "J", nthreads = 8);gc()

Runs$run_id<-sf_gsub(Runs$run_id, "25229685dfba4d12977f2098cc1aba89cabbf9a3", "K", nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "b028b1050684717631972d0f22c5e5e7651ab958", "L", nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "100751dda22e1d7db084cb563d096baccb5c97e2", "M",  nthreads = 8);gc()

Runs$run_id<-sf_gsub(Runs$run_id, "e3a0ea18f8a773e7a4f6a521f7f2318ec8a2ecac", "N", nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "dccd8e29b6563fc3c3bad7039ac31b59f9d1b763", "O",  nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "2deddfc3491835138fd09ab63e2a97d8bf7d9794", "P", nthreads = 8);gc()

Runs$run_id<-sf_gsub(Runs$run_id, "08b541dcf0c98e672d5042fdd3f74ded9097ba25", "Q", nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "af2c258224dd7222ae544c06f4f8e35abc7a56ac", "R",  nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "d8ff8d72561494a9d48c8be6d66d2d58688e11a7", "S", nthreads = 8);gc()

Runs$run_id<-sf_gsub(Runs$run_id, "56fc8c0acdeb71b50c0a7bbd72c93c1d57d38da4", "T", nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "e197af606aa8a7f242ac0f9a293093dbd79c9788", "U",  nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "5141fd5506dec6195e5baa03b975161ab939c748", "V", nthreads = 8);gc()

Runs$run_id<-sf_gsub(Runs$run_id, "8727616b3c6c7863250d3cc5cd90d2c23c46d072", "W", nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "b2a0ace5558a6a9427faf86e4bded3702597d149", "X",  nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "a49a6e89327ed1ee0c09762ce1b51d6349b02f49", "Y", nthreads = 8);gc()
Runs$run_id<-sf_gsub(Runs$run_id, "7215d1c8de62dc81e7dfe5cc2b031c4768ff2bff", "Y", nthreads = 8);gc()

Runs$run_id<-sf_gsub(Runs$run_id, "c7648e4555e47bfd0471d15c4daecd98e1ed47c4", "Negatives",nthreads = 8);gc()

Runs$filt <- factor(Runs$filt , levels=c("Raw", "Filtered"))

head(Runs)
tail(Runs)

library("plyr")
summary_r <- ddply(Runs, c("filt", "run_id"), summarise,
               N    = length(mean_qscore_template),
               median = median(mean_qscore_template),
               mean = mean(mean_qscore_template),
               sd   = sd(mean_qscore_template),
               se   = sd / sqrt(N));summary_r

library(ggplot2)
library(hrbrthemes)

# curves for q-scores  ##

D<-Runs %>%
ggplot(aes(x=mean_qscore_template, group=filt, fill=filt)) +
  geom_density(adjust=1.5, alpha=.4)+
  #facet_grid(.~filt)+
  xlim(6,30)+
  theme(
    legend.title = element_blank(), 
    legend.position = c(0.9, 0.8),
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8))+
  xlab("Q-score")+
  ylab("Density");D

# curves for read lenght  ##

summary2_r <- ddply(Runs, c("filt", "run_id"), summarise,
                 N    = length(sequence_length_template),
                 median = median(sequence_length_template),
                 mean = mean(sequence_length_template),
                 sd   = sd(sequence_length_template),
                 se   = sd / sqrt(N));summary2_r

E<-Runs %>%
  ggplot(aes(x=filt, y=sequence_length_template, group=filt, fill=filt)) +
  geom_violin(adjust=1.5)+
  #facet_grid(. ~ filt, scales = "free")+
  scale_y_continuous(trans='log10', limits = c(100, 7000))+
  theme(
    legend.title = element_blank(), 
    legend.position = "none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8))+
  xlab("")+
  ylab("Read lenght (bp)");E

# curves for GC content  ##

summary3_r <- ddply(Runs, c("filt", "run_id"), summarise,
                  N    = length(gc_percent),
                  median = median(gc_percent),
                  mean = mean(gc_percent),
                  sd   = sd(gc_percent),
                  se   = sd / sqrt(N));summary3_r

F<-Runs %>%
  ggplot(aes(x=filt, y=gc_percent, group=filt, fill=filt)) +
  geom_violin(adjust=1.5)+
  #facet_grid(. ~ filt, scales = "free")+
  ylim(30, 60)+
  theme(
    legend.title = element_blank(), 
    legend.position = "none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8))+
  ylab("GC content (%)")+
  xlab("");F


library(ggpubr)
ggarrange(D,                                                 
          ggarrange(E, F, ncol = 2, labels = c("B", "C")),
          nrow = 2, 
          labels = "A")


ggsave(file="Figures/FigS2_summary_Runs.pdf", width=5.5, height=6, units="in", dpi=900)


