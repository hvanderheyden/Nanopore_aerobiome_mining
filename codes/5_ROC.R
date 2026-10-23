

################################################################################
########################### ROC Analysis #######################################
################################################################################

####################### P. destructor ##########################################

library("phyloseq")
library("tidyverse")
library("pROC")

Runs_PRF <- readRDS("R_objects/Runs_rawP.rds")

p_dest <- subset_taxa(Runs_PRF, Species=="Peronospora destructor")

p_dest_dat <- psmelt(p_dest)

p_dest_dat$P_destructor <- as.numeric(p_dest_dat$P_destructor)

p_dest_dat<-p_dest_dat%>%
  mutate(NoT=ifelse(p_dest_dat$P_destructor > 0, 1, 0))

shapiro.test(p_dest_dat$P_destructor) # => p < 2.2e-16
shapiro.test(p_dest_dat$Abundance) # => p < 2.2e-16

corr_PD <- cor.test(log(p_dest_dat$P_destructor+1), 
                    log(p_dest_dat$Abundance+1), 
                    method = "pearson");corr_PD

ggplot(
  p_dest_dat %>% 
    filter(!is.na(Abundance), !is.na(P_destructor), Abundance > 0, P_destructor > 0),
  aes(x = log(Abundance), y = log(P_destructor))
) +
  geom_point() +
  stat_smooth(method = "lm") +
  ylim(0, 7)

library(cutpointr)

set.seed(1324)
cp1 <- cutpointr(p_dest_dat, 
                Abundance, 
                NoT,
                method = maximize_metric,
                metric = youden, 
                na.rm = TRUE)

summary(cp1)

library(plotROC)
set.seed(1324)
P1<-ggplot(p_dest_dat, aes(d = NoT, m = Abundance )) + 
  geom_roc(n.cuts = 9, labelsize = 2, labelround = 2) + 
  style_roc(xlab = "1-Specificity (FPR)",
            ylab = "Sensitivity (TPR)")+ 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="darkgrey", linetype="dashed")

P1<-P1+
  annotate('text', x = 0.61, y = .25,
           label="paste(italic(r)==0.800)", 
           parse=TRUE)+
  annotate('text', x = .85, y = .25,
           label="(italic(P)<0.0001)", 
           parse=TRUE)+
  annotate("text", x = 0.67, y = .17, 
           label = paste("AUC =", 
                         round(calc_auc(P1)$AUC, 3)))+
  annotate('text', x = 0.65, y = .10,
           label="ACC == 0.860", 
           parse=TRUE);P1



####################### B. squamosa ############################################

B_sq <- subset_taxa(Runs_PRF, Species=="Botrytis squamosa")

B_sq_dat<-psmelt(B_sq)

#write.csv(B_sq_dat, "Data/squamosa")

B_sq_dat$B_squamosa <- as.numeric(B_sq_dat$B_squamosa)

B_sq_dat<-B_sq_dat%>%
  mutate(NoT=ifelse(B_sq_dat$B_squamosa > 0, 1, 0))%>%
  mutate(Abund=ifelse(B_sq_dat$Abundance > 0, Abundance, 0))

#str(B_sq_dat)

shapiro.test(B_sq_dat$B_squamosa) # => p < 2.2e-16
shapiro.test(B_sq_dat$Abund) # => p < 2.2e-16

corr_bs <- cor.test(log(B_sq_dat$B_squamosa+1), 
                    log(B_sq_dat$Abund+1), 
                    method = "pearson");corr_bs

ggplot(
  B_sq_dat %>% 
    filter(!is.na(Abund), !is.na(B_squamosa), Abund > 0, B_squamosa > 0),
  aes(x = log(Abund), y = log(B_squamosa))
) +
  geom_point() +
  stat_smooth(method = "lm") +
  ylim(0, 8)


library(Epi)
library(cutpointr)

set.seed(1324)
cp <- cutpointr(B_sq_dat, 
                Abund, 
                NoT,
                method = maximize_metric,
                metric = youden, 
                na.rm = TRUE)
summary(cp)

library(plotROC)
set.seed(1324)
P2 <- ggplot(B_sq_dat, aes(d=NoT, m=Abund)) + 
  geom_roc(n.cuts = 9, labelsize = 2, labelround = 2)+
  style_roc(
    xlab = "1-Specificity (FPR)",
    ylab = "Sensitivity (TPR)"
  )+ 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="darkgrey", linetype="dashed")


P2<-P2+
  annotate('text', x = 0.61, y = .25,
           label="paste(italic(r)==0.652)", 
           parse=TRUE)+
  annotate('text', x = .85, y = .25,
           label="(italic(P)<0.0001)", 
           parse=TRUE)+
  annotate("text", x = 0.65, y = .17, 
           label = paste("AUC =", 
                         round(calc_auc(P2)$AUC, 3)))+
  annotate('text', x = 0.64, y = .10,
           label="ACC == 0.790", 
           parse=TRUE);P2


######### 

library(ggpubr)
ggarrange(P1, P2, ncol = 2, labels = c("A", "B"))

ggsave(file="Figures/Fig4_ROC.pdf", 
       width=8, height=4, units="in", dpi=900)










