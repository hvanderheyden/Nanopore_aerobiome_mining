
library(forcats)
library(tidyverse)

########################################
############## QC Mocs #################
########################################


mock_summary<-read.csv("Data/QC/mock_summary.txt", sep="\t")
head(mock_summary)

mock_summary$filt <- c("Raw")

mock_final_summary<-read.csv("Data/QC/mock_summary_final.txt", sep="\t")
head(mock_final_summary)

mock_final_summary$filt <- c("Filtered")

mock<-rbind(mock_summary, mock_final_summary)

library(stringfish)

mock$run_id<-sf_gsub(mock$run_id, "2e592184e983f643460022d81e7659e55cc8e59f", "Mock", nthreads = 24)
mock$run_id<-sf_gsub(mock$run_id, "c7648e4555e47bfd0471d15c4daecd98e1ed47c4", "Mock2",  nthreads = 24)
mock$run_id<-sf_gsub(mock$run_id,"dcd81b2b714b74a5070320859195b98658ca9204", "Mock3", nthreads = 24)


head(mock)
tail(mock)

mock$filt <- factor(mock$filt , levels=c("Raw", "Filtered"))


library("plyr")
summary <- ddply(mock, c("filt","run_id"), summarise,
               N    = length(mean_qscore_template),
               median = median(mean_qscore_template),
               mean = mean(mean_qscore_template),
               sd   = sd(mean_qscore_template),
               se   = sd / sqrt(N));summary



library(ggplot2)
library(hrbrthemes)

# curves for q-scores  ##

mA<-mock %>%
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
  ylab("Density");mA

# curves for read lenght  ##

library("plyr")
summary2 <- ddply(mock, c("filt", "run_id"), summarise,
                 N    = length(sequence_length_template),
                 median = median(sequence_length_template),
                 mean = mean(sequence_length_template),
                 sd   = sd(sequence_length_template),
                 se   = sd / sqrt(N));summary2

mB<-mock %>%
  ggplot(aes(x=filt, y=sequence_length_template, group=filt, fill=filt)) +
  geom_violin(adjust=1.5)+
  #facet_grid(. ~ filt, scales = "free")+
  scale_y_continuous(trans='log10', limits = c(10, 7000))+
  theme(
    legend.title = element_blank(), 
    legend.position = "none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8))+
  xlab("")+
  ylab("Read lenght (bp)");mB

# curves for GC content  ##

mC<-mock %>%
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
  xlab("");mC


library(ggpubr)
ggarrange(mA,                                                 
          ggarrange(mB, mC, ncol = 2, labels = c("B", "C")),
          nrow = 2, 
          labels = "A")


ggsave(file="Figures/FigS1_summary_mock.pdf", width=5.5, height=6, units="in", dpi=900)


