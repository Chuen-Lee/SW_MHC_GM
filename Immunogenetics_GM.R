# Immunogenetics in the Seychelles warblers - 16s + shotgun metagenomics
# Host immunogenetics and Microbiome 
# How do we do this? What do we want to compare?
# TLR and MHC, controlling for genome wide heterozygosity, age, sex, season,sample year, catch time, and number of days the samples were stored in the fridge.

# libraries ----
library("tidyverse")
library("microbiome")
library("speedyseq")
library("arm")
library("vegan")
library("lmerTest")
library("DHARMa")
library("car")
library("ggeffects")
library("ggthemes")
library("MuMIn")
library("ANCOMBC")
library("gllvm")
library("ggpubr")
library("ggrepel")
#library("ALDEx2")
library("ggVennDiagram")
#library("decontam")

setwd("~/Documents/PhD/R_analysis/HostGenetics/")

## functions ----
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

alpha_estimate <- function(physeq) {
  ttt <- transform_sample_counts(physeq, function(x) round(x, 0))
  richnessEstRare<-estimate_richness(ttt, split=TRUE, measures= c("Chao1", "Shannon", "Observed"))
  RareMeta <- as.data.frame(sample_data(ttt))
  RareMeta$Chao1 <- richnessEstRare$Chao1
  RareMeta$Shannon <- richnessEstRare$Shannon
  RareMeta$Observed <- richnessEstRare$Observed
  AlphaMeta <- as_tibble(RareMeta) %>% mutate(TerminalYear = as.factor(TerminalYear)) %>% mutate(Timeinfridge = as.numeric(Timeinfridge)) %>% 
    mutate(CatchTime = tidyr::replace_na(CatchTime, 43292)) %>% mutate(CatchTime = as.numeric(CatchTime))
  return(AlphaMeta)
}
# TLR3 and MHC from Charli
TL3_MHC_20240530 <- read_csv("InputTables/TL3_MHC_20240530.csv") %>% mutate_at(vars(contains("Ase")), funs(as.factor)) %>% mutate(TLR3new =case_when(TLR3 %in% "A:A" ~ "HomA", TLR3 %in% "C:C" ~ "HomC", TLR3 %in% "A:C" ~ "Het"))
TL3_MHC_20240530full <- TL3_MHC_20240530 %>% dplyr::select(-Comment) %>% na.omit
TL3_MHC_20240530 %>% distinct(TLR3new)

# heterozygosity
SWhet <- read.csv("Genhet_heterozygosity/SWheterozygosity_17_05_22.txt",sep="\t")
MHChet <- merge(TL3_MHC_20240530full,SWhet, by="BirdID")
ggplot(MHChet, aes(x=MHC2_Diversity,y=Hs_obs,group=MHC2_Diversity)) + geom_boxplot()
ggplot(MHChet, aes(x=MHC1_Diversity,y=Hs_obs,group=MHC1_Diversity)) + geom_boxplot()

MHC1hetlm <- lm(MHC1_Diversity ~ Hs_obs, MHChet)
summary(MHC1hetlm)
ggplot(MHChet, aes(MHC1_Diversity)) + geom_histogram()
ggplot(MHChet, aes(Hs_obs)) + geom_histogram()
ggplot(MHChet, aes(MHC1_Diversity,Hs_obs)) + geom_point(position = "jitter") + geom_smooth() + xlim(0,8)
MHC2hetlm <- lm(MHC2_Diversity ~ Hs_obs, MHChet)
summary(MHC2hetlm)


# 16s ----
physeq3TLR <- read_rds("InputTables/physeq3TLR.rds")

#mhc 1 diversity metrics
physeq3TLR %>% summarise(meanMHC1 = mean(MHC1_Diversity), sd=(sd(MHC1_Diversity)), se= sd(MHC1_Diversity/sqrt(n())))
physeq3TLR %>% summarise(minMHC1 = min(MHC1_Diversity))
physeq3TLR %>% summarise(maxMHC1 = max(MHC1_Diversity))

#mhc 2 diversity metrics
physeq3TLR %>% summarise(meanMHC2 = mean(MHC2_Diversity), sd=(sd(MHC2_Diversity)), se= sd(MHC2_Diversity/sqrt(n())))
physeq3TLR %>% summarise(minMHC2 = min(MHC2_Diversity))
physeq3TLR %>% summarise(maxMHC2 = max(MHC2_Diversity))

ggplot(physeq3TLR, aes(MHC1_Diversity)) + geom_histogram()
ggplot(physeq3TLR, aes(MHC2_Diversity)) + geom_histogram()

ggplot(physeq3TLR, aes(Ageclass,MHC1_Diversity, group=Ageclass)) + geom_point() + geom_boxplot()
ggplot(physeq3TLR, aes(Ageclass,MHC2_Diversity, group=Ageclass)) + geom_point() + geom_boxplot()

physeq3TLRp <- phyloseq(tax_table(physeq4A),otu_table(physeq4A),sample_data(physeq3TLR))
physeq4ataxa <- data.frame(tax_table(physeq4A)) %>% rownames_to_column("taxon")


## 16s alpha diversity ----
physeq3TLRprarefy <- rarefy_even_depth(physeq3TLRp, rngseed=88,sample.size = 8000, replace=F )

physeq3TLRpalpha <- alpha_estimate(physeq3TLRprarefy)%>% mutate_at(vars(contains('Ase')), as.factor) %>% mutate(SampleYear = as.factor(SampleYear)) %>% mutate(cutage = pmin(SamplingAge,10))

ggplot(physeq3TLRpalpha, aes(x=Observed)) + geom_histogram()
physeqobs <- glmer.nb(Observed ~ Hs_obs + Asedab3 + Asedab4 + Asedab5 + Aseua1+Aseua3+Aseua4+Aseua5+Aseua6+Aseua7+ Aseua8+ Aseua9 + Aseua11 + cutage + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled+ SampleYear + (1|BirdID), data=physeq3TLRpalpha,nAGQ=0 )
summary(physeqobs)
car::Anova(physeqobs,type="III")
vif(physeqobs)
simulateResiduals(physeqobs,plot=T)
r.squaredGLMM(physeqobs)

physeqobs2 <- glmer.nb(Observed ~ Hs_obs + MHC1_Diversity +MHC2_Diversity + + cutage+season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled  + SampleYear + (1|BirdID), data=physeq3TLRpalpha,nAGQ=0 )
summary(physeqobs2)
car::Anova(physeqobs2,type="III")

physeqobs3 <- glmer.nb(Observed ~ Hs_obs + I(Hs_obs^2)+ MHC1_Diversity + I(MHC1_Diversity^2) +MHC2_Diversity  + I(MHC2_Diversity^2) +cutage + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + SampleYear + (1|BirdID), data=physeq3TLRpalpha ,nAGQ=0)
summary(physeqobs3)
car::Anova(physeqobs3,type="III")

ggplot(physeq3TLRpalpha, aes(x=Shannon)) + geom_histogram()
physeqsha <- lmer(Shannon ~ Hs_obs + Asedab3 + Asedab4 + Asedab5 + Aseua1+Aseua3+Aseua4+Aseua5+Aseua6+Aseua7+ Aseua8+ Aseua9 + Aseua11 + cutage + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled  + SampleYear + (1|BirdID), data=physeq3TLRpalpha)
summary(physeqsha)
car::Anova(physeqsha,type="III")
vif(physeqsha)
simulateResiduals(physeqsha,plot=T)
r.squaredGLMM(physeqsha)

physeqsha2 <- lmer(Shannon ~ Hs_obs+ MHC1_Diversity +MHC2_Diversity + cutage + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + SampleYear + (1|BirdID), data=physeq3TLRpalpha)
summary(physeqsha2)
car::Anova(physeqsha2,type="III")

physeqsha3 <- lmer(Shannon ~ Hs_obs + I(Hs_obs^2)+ MHC1_Diversity + I(MHC1_Diversity^2) +MHC2_Diversity  + I(MHC2_Diversity^2) +SamplingAge + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + SampleYear + (1|BirdID), data=physeq3TLRpalpha)
summary(physeqsha3)


## 16s beta ----
physeq3TLRprare <- transform(physeq3TLRp,"compositional") %>% core_members(., detection = 0.001, prevalence = 0.05)
physeq3TLRpCLR <- prune_taxa(physeq3TLRprare,physeq3TLRp) %>% transform(.,"clr") %>% mutate_sample_data(MHC1Dbin = as.factor(case_when(MHC1_Diversity < 4 ~ "<4", MHC1_Diversity <7 ~ "4-6",MHC1_Diversity >= 7 ~ "7")),MHC2Dbin = as.factor(case_when(MHC2_Diversity < 2 ~ "1", MHC2_Diversity <= 4 ~ "2-4",MHC2_Diversity >= 5 ~ "5")), Hsbin = as.factor(case_when(Hs_obs < 1 ~ "<1", Hs_obs < 1.1 ~ "1-1.1", Hs_obs >= 1.1 ~ ">1.1")))

physeq3TLRpCLRord <- phyloseq::ordinate(physeq3TLRpCLR, method="RDA",distance = "euclidean")

Aseua5_16s_PCA <- plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="Aseua5") + 
  stat_ellipse() + 
  scale_colour_manual(name = expression(italic("Ase-ua5")),values = c("deepskyblue2","orangered3")) + 
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))
Aseua7_16s_PCA <-plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="Aseua7") + 
  stat_ellipse()+ 
  scale_colour_manual(name = expression(italic("Ase-ua7")),values = c("deepskyblue2","orangered3")) + 
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))
Aseua9_16s_PCA <-plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="Aseua9") + 
  stat_ellipse()+ 
  scale_colour_manual(name = expression(italic("Ase-ua9")),values = c("deepskyblue2","orangered3")) + 
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

#mhc1bin
mhc1binpcadf <- plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="MHC1Dbin",justDF = T) %>% group_by(MHC1Dbin) %>% dplyr::summarise(meanPC1 = mean(PC1), meanPC2 = mean(PC2))
mhc1binpcaplot <- plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="MHC1Dbin") + 
  geom_point(data=mhc1binpcadf,aes(x=meanPC1,y=meanPC2,fill=MHC1Dbin), shape=23, stroke=1, size=5, colour="black") +
  scale_colour_manual("MHC-I Diversity",values = c("olivedrab2","deepskyblue2","orangered3"))+
  scale_fill_manual("MHC-I Diversity",values = c("olivedrab2","deepskyblue2","orangered3"))+ 
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

#mhc2bin
mhc2binpcadf <- plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="MHC2Dbin",justDF = T)%>% group_by(MHC2Dbin) %>% dplyr::summarise(meanPC1 = mean(PC1), meanPC2 = mean(PC2))
mhc2binpcaplot <- plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="MHC2Dbin") +
  geom_point(data=mhc2binpcadf,aes(x=meanPC1,y=meanPC2,fill=MHC2Dbin), shape=23, stroke=1, size=5, colour="black") + 
  scale_colour_manual("MHC-II Diversity",values = c("olivedrab2","deepskyblue2","orangered3")) +
  scale_fill_manual("MHC-II Diversity",values = c("olivedrab2","deepskyblue2","orangered3"))+
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

#hsbin 
hsobsbinpcadf <- plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="Hsbin",justDF = T)%>% group_by(Hsbin) %>% dplyr::summarise(meanPC1 = mean(PC1), meanPC2 = mean(PC2))
plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="Hsbin") + 
  stat_ellipse(level = 0.3,type = "euclid")+ 
  geom_point(data=hsobsbinpcadf,aes(x=meanPC1,y=meanPC2,fill=Hsbin), shape=23, stroke=1, size=5, colour="black")+
  scale_colour_manual(values = c("olivedrab2","deepskyblue2","orangered3")) +
  scale_fill_manual(values = c("olivedrab2","deepskyblue2","orangered3")) +
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

hsobs_16s_pca <- plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="Hs_obs")+ 
  scale_color_gradient2(
    low = "olivedrab2", 
    mid = "skyblue2", 
    high = "orangered3", 
    midpoint = 1.1)+ 
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

mhc1_16s_pca <- plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="MHC1_Diversity")+ 
  scale_color_gradient2(
    low = "olivedrab2", 
    mid = "skyblue2", 
    high = "orangered1", 
    midpoint = 4.5)+ 
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

mhc2_16s_pca <- plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="MHC2_Diversity")+ 
  scale_color_gradient2(
    low = "olivedrab2", 
    mid = "skyblue2", 
    high = "orangered1", 
    midpoint = 3)+ 
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

tiff("Output/MHC_16s_PCA.tiff", res =400, units = "in", width = 16, height = 9, bg = "transparent")
ggarrange(ggarrange( mhc1binpcaplot,mhc2binpcaplot, nrow = 1, labels = c("A","B")), ggarrange(Aseua5_16s_PCA,Aseua7_16s_PCA,Aseua9_16s_PCA, nrow=1, labels = c("C","D","E")), nrow = 2, ncol = 1)
dev.off()

mhc1pcs <- plot_ordination(physeq3TLRpCLR,physeq3TLRpCLRord,axes = 1:2 ,color="MHC1_Diversity",justDF = T)

mhc1pcslm <- lmer(PC1 ~ MHC1_Diversity + I(MHC1_Diversity^2)  + (1),mhc1pcs)
summary(mhc1pcslm)

physeq3TLRpCLR_mat <- vegan_otu(physeq3TLRpCLR)
physeq3TLRpCLR_st <- as(sample_data(physeq3TLRpCLR),"data.frame") %>% group_by(BirdID) %>% dplyr::select(SamplingAge,TerminalYear, season, SampleYear,SexEstimate, Timeinfridge_scaled,CatchTime_scaled,BirdID,contains("Ase"),MHC1_Diversity,MHC2_Diversity,Hs_obs,TQcorrected_scaled,survive) %>% mutate_at(vars(contains('Ase')), as.factor) %>% mutate(SampleYear = as.factor(SampleYear), cutage = pmin(SamplingAge,10))

perm <- how(nperm = 9999, blocks = physeq3TLRpCLR_st$BirdID)
set.seed(888)
physeq3TLRpCLRperm<- adonis2(physeq3TLRpCLR_mat ~ Hs_obs + Asedab3+Asedab4+Asedab5 + Aseua1+Aseua3+Aseua4+Aseua5+Aseua6+Aseua7+ Aseua8+Aseua9+ Aseua11  +cutage+ season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled    , data=physeq3TLRpCLR_st, permutations = perm, method = "euclidean", by= "margin")
physeq3TLRpCLRperm
AICc_permanova2(physeq3TLRpCLRperm)

physeq3TLRpCLRperm_S<- adonis2(physeq3TLRpCLR_mat ~ Hs_obs + Aseua5+Aseua7+Aseua9  + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled , data=physeq3TLRpCLR_st, permutations = perm, method = "euclidean", by= "margin")
physeq3TLRpCLRperm_S

physeq3TLRpCLRperm2<- adonis2(physeq3TLRpCLR_mat ~ survive + Hs_obs + MHC1_Diversity + I(MHC1_Diversity^2) + MHC2_Diversity + I(MHC2_Diversity^2) + cutage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled   , data=physeq3TLRpCLR_st, permutations = perm, method = "euclidean", by= "margin")
physeq3TLRpCLRperm2

physeq3TLRpCLRperm2big<- adonis2(physeq3TLRpCLR_mat ~ Hs_obs + I(Hs_obs^2) + MHC1_Diversity + I(MHC1_Diversity^2) + MHC2_Diversity + I(MHC2_Diversity^2) + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled   , data=physeq3TLRpCLR_st, permutations = perm, method = "euclidean", by= "margin")
physeq3TLRpCLRperm2big

physeq3TLRpCLRperm3<- adonis2(physeq3TLRpCLR_mat ~ Hs_obs + MHC1_Diversity + I(MHC1_Diversity^2) + MHC2_Diversity+ SamplingAge + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled +TQcorrected_scaled   , data=physeq3TLRpCLR_st, permutations = perm, method = "euclidean", by= "margin")
physeq3TLRpCLRperm3

# mpa ----
mpaTLR <- readRDS("InputTables/mpaTLR.rds")

## mpa alpha ----
mpaTLRrarefy <- rarefy_even_depth(mpaTLR, rngseed=88,sample.size = 5500, replace=F)

mpaTLRalpha <- alpha_estimate(mpaTLRrarefy)%>% mutate(SamplingAge2 = round(SamplingAge)) %>% mutate_at(vars(contains('Ase')), as.factor) %>% mutate(cutage = pmin(SamplingAge,10))
n_distinct(mpaTLRalpha$BirdID) # 57 birds 
ggplot(mpaTLRalpha, aes(x=Observed)) + geom_histogram()

mpaobs <- glmer.nb(Observed ~ Aseua11 + cutage+ season + SexEstimate + Timeinfridge_scaled + CatchTime_scaled + SampleYear  + (1|BirdID), data=mpaTLRalpha,nAGQ=0)
summary(mpaobs)
vif(mpaobs)
anova(mpaobs)
car::Anova(mpaobs, type = 3)
simulateResiduals(mpaobs,plot=T)

mpaobs2 <- glmer.nb(Observed ~ Hs_obs + MHC1_Diversity + MHC2_Diversity + cutage + season + SexEstimate + Timeinfridge_scaled + CatchTime_scaled + SampleYear + (1|BirdID), data=mpaTLRalpha,nAGQ=0 )
summary(mpaobs2)
car::Anova(mpaobs2, type = 3)
simulateResiduals(mpaobs2,plot=T)

ggplot(mpaTLRalpha, aes(x=Shannon)) + geom_histogram()
mpashan <- lmer(Shannon ~ Aseua11+ cutage+ season + SexEstimate + Timeinfridge_scaled + CatchTime_scaled + SampleYear   + (1|BirdID), data=mpaTLRalpha )
summary(mpashan)
vif(mpashan)
Anova(mpashan, type = "III")
simulateResiduals(mpashan,plot=T)

mpashan2 <- lmer(Shannon ~ Hs_obs + MHC1_Diversity + MHC2_Diversity +cutage +season + SexEstimate + Timeinfridge_scaled + CatchTime_scaled + SampleYear + (1|BirdID), data=mpaTLRalpha )
summary(mpashan2)
simulateResiduals(mpashan2,plot=T)
car::Anova(mpashan2, type="III")

ggplot(mpaTLRalpha,aes(log(Chao1))) + geom_histogram()
mpachao <- lmer(log(Chao1) ~ Asedab3+Asedab4+Asedab5 + Aseua1+Aseua3+Aseua4+Aseua5+Aseua6+Aseua7+ Aseua8+ Aseua11 +SamplingAge2 +TerminalYear + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + SampleYear + (1|BirdID), data=mpaTLRalpha )
summary(mpachao)
simulateResiduals(mpachao,plot=T)

mpachao2 <- lmer(log(Chao1) ~ MHC1_Diversity + MHC2_Diversity +SamplingAge2 +TerminalYear + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + SampleYear + (1|BirdID), data=mpaTLRalpha )
summary(mpachao2)

## mpa beta ----
mpaTLRrare <- transform(mpaTLR,"compositional") %>% core_members(., detection = 0.001, prevalence = 0.05)
mpaTLRCLR <- prune_taxa(mpaTLRrare,mpaTLR) %>% transform(.,"clr") %>% mutate_sample_data(MHC2Dbin = as.factor(case_when(MHC2_Diversity < 2 ~ "1", MHC2_Diversity <= 4 ~ "2-4",MHC2_Diversity >= 5 ~ "5")), Hsbin = as.factor(case_when(Hs_obs < 1 ~ "<1", Hs_obs < 1.1 ~ "1-1.1", Hs_obs >= 1.1 ~ ">1.1")))

mpaTLRCLRord <- ordinate(mpaTLRCLR, method="RDA",distance = "euclidean")

Aseua7_mpa_pca <- plot_ordination(mpaTLRCLR,mpaTLRCLRord,axes = 1:2 ,color="Aseua7") + stat_ellipse()+ 
  scale_colour_manual(name = expression(italic("Ase-ua7")),values = c("deepskyblue2","orangered3")) + 
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

hsobs_mpa_pca <-plot_ordination(mpaTLRCLR,mpaTLRCLRord,axes = 1:2 ,color="Hs_obs")+ 
  scale_color_gradient2(name="Heterozygosity",
                        low = "olivedrab2", 
                        mid = "deepskyblue2", 
                        high = "orangered3", 
                        midpoint = 1.1)+ 
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

mhc2_mpa_pca <-plot_ordination(mpaTLRCLR,mpaTLRCLRord,axes = 1:2 ,color="MHC2_Diversity")+ 
  scale_color_gradient2(
    low = "olivedrab2", 
    mid = "skyblue2", 
    high = "orangered1", 
    midpoint = 3)+ 
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

#mhc2bin
mhc2binmpadf <- plot_ordination(mpaTLRCLR,mpaTLRCLRord,axes = 1:2 ,color="MHC2Dbin",justDF = T)%>% group_by(MHC2Dbin) %>% dplyr::summarise(meanPC1 = mean(PC1), meanPC2 = mean(PC2))
mhc2binmpaplot <- plot_ordination(mpaTLRCLR,mpaTLRCLRord,axes = 1:2 ,color="MHC2Dbin") +
  geom_point(data=mhc2binmpadf,aes(x=meanPC1,y=meanPC2,fill=MHC2Dbin), shape=23, stroke=1, size=5, colour="black") + 
  scale_colour_manual("MHC-II Diversity",values = c("olivedrab2","deepskyblue2","orangered3")) +
  scale_fill_manual("MHC-II Diversity",values = c("olivedrab2","deepskyblue2","orangered3"))+
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

#hsbin 
hsobsbinmpaadf <- plot_ordination(mpaTLRCLR,mpaTLRCLRord,axes = 1:2 ,color="Hsbin",justDF = T)%>% group_by(Hsbin) %>% dplyr::summarise(meanPC1 = mean(PC1), meanPC2 = mean(PC2))
hsobsbinmpaaplot <- plot_ordination(mpaTLRCLR,mpaTLRCLRord,axes = 1:2 ,color="Hsbin") + 
  geom_point(data=hsobsbinmpaadf,aes(x=meanPC1,y=meanPC2,fill=Hsbin), shape=23, stroke=1, size=5, colour="black")+
  scale_colour_manual(name = "Heterozygosity",values = c("olivedrab2","deepskyblue2","orangered3")) +
  scale_fill_manual(name = "Heterozygosity",values = c("olivedrab2","deepskyblue2","orangered3")) +
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

tiff("Output/mpa_PCA.tiff", res =400, units = "in", width = 10, height = 6, bg = "transparent")
ggarrange(hsobsbinmpaaplot,mhc2binmpaplot,Aseua7_mpa_pca, labels = c("A","B","C"))
dev.off()

mpaTLRCLRmat <-vegan_otu(mpaTLRCLR)
mpaTLRCLR_st <- as(sample_data(mpaTLRCLR),"data.frame") %>% 
  mutate(survive = as.factor(survive)) %>% group_by(BirdID) %>%
  mutate(deltaAge = capage - mean(capage), meanAge = mean(capage),TerminalYearbird = case_when(any(TerminalYear %in% "1") ~ "1", TRUE ~ "0")) %>% dplyr::select(SamplingAge,TerminalYear, season, SampleYear,SexEstimate, Timeinfridge_scaled,CatchTime_scaled,TQcorrected_scaled,BirdID,contains("Ase"),MHC1_Diversity,MHC2_Diversity,Hs_obs) %>% mutate(cutage = pmin(SamplingAge,10))

perm <- how(nperm = 9999, blocks = mpaTLRCLR_st$BirdID)
set.seed(888)

mpaTLRCLRperm2S<- adonis2(mpaTLRCLRmat ~ Aseua5 + Aseua7+Aseua9 + cutage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled, data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
mpaTLRCLRperm2S

mpaTLRCLRperm4<- adonis2(mpaTLRCLRmat ~Hs_obs + MHC1_Diversity + MHC2_Diversity+ cutage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled   , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
mpaTLRCLRperm4



## mpa DAA ----
mpaTLRcoremembers <- transform(mpaTLR, "compositional") %>% core_members(., detection = 0.00001, prevalence = 0.1)
mpaTLRcore <- prune_taxa(mpaTLRcoremembers, mpaTLR)

mpaTLRcoreclrmelt <- psmelt(mpaTLRcoreclr)
mpaTLRcoreclrmelt92 <- mpaTLRcoreclrmelt %>% filter(OTU %in% "OTU92")
ggplot(mpaTLRcoreclrmelt92, aes(Hs_obs,Abundance)) + geom_point() + geom_smooth(method = "lm")
mpaTLRcoreclrmelt92lm <- lmer(Abundance ~ Hs_obs+ Aseua7 + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + (1|BirdID),data=mpaTLRcoreclrmelt92)
summary(mpaTLRcoreclrmelt92lm)

mpaTLRcoremelt <- psmelt(transform(mpaTLRcore,"compositional"))
mpaTLRcoremelt409 <- mpaTLRcoremelt %>% filter(OTU %in% "OTU409")
ggplot(mpaTLRcoremelt409, aes(Hs_obs,Abundance)) + geom_point() + geom_smooth(method = "lm")
mpaTLRcoremelt409lm <- lmer(Abundance ~ Hs_obs+ Aseua7 + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + (1|BirdID),data=mpaTLRcoremelt409)
summary(mpaTLRcoremelt409lm)
mpaTLRcoremelt409lmdata <- ggpredict(mpaTLRcoremelt409lm, terms="Hs_obs")
plot(mpaTLRcoremelt409lmdata, show_data = T)

mpaTLRcoreclrmelt409 <- mpaTLRcoreclrmelt %>% filter(OTU %in% "OTU409")
ggplot(mpaTLRcoreclrmelt409, aes(Hs_obs,Abundance)) + geom_point() + geom_smooth(method = "lm")
mpaTLRcoreclrmelt409lm <- lmer(Abundance ~ Hs_obs+ Aseua7 + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + (1|BirdID),data=mpaTLRcoreclrmelt409)
summary(mpaTLRcoreclrmelt409lm)
mpaTLRcoreclrmelt409lmdata <- ggpredict(mpaTLRcoreclrmelt409lm, terms="Hs_obs")
plot(mpaTLRcoreclrmelt409lmdata, show_data = T)


### aldex2 mpa ----
library(ALDEx2)
mpataxo7spec <- mpataxo7 %>% rownames_to_column("taxonsig")
mpaTLRcoremat <- otu_table(mpaTLRcore)
mpaTLRcorest <- data.frame(sample_data(mpaTLRcore))
mm <- model.matrix(~ Hs_obs + MHC1_Diversity + MHC2_Diversity + Aseua7 + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled, mpaTLRcorest) 
x.aldex <- ALDEx2::aldex(mpaTLRcoremat, mm, mc.samples=200, test="glm", effect=TRUE, include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL)

x.aldexhsobs <- x.aldex %>% dplyr::select(contains("Hs_obs")) %>% rownames_to_column("OTU") %>% rename(pval = Hs_obs.pval,) %>% mutate(direct = case_when(Hs_obs.Est < 0 & pval < 0.05 ~ "p05", Hs_obs.Est > 0 & pval < 0.05 ~ "p05", TRUE ~ "p > 0.05"),cog = case_when(pval < 0.05 ~ OTU)) %>% mutate(taxonsig = case_when(pval < 0.05 ~ OTU)) %>% merge(.,mpataxo7spec, by = "taxonsig", all.x=T)%>% mutate(Species2=str_extract(Species, "(?<=__).*"))

hsobsmpaaldexplot <- ggplot(x.aldexhsobs, aes(x= Hs_obs.Est,y=  -log(pval),color=direct, label=Species2)) +
  geom_point(size =3) +
  geom_text_repel(size=3, colour="black") +
  ggtitle("A. Heterozygosity")+
  xlab("Log fold change Heterozygosity") +
  ylab("-log(p-value)")+
  scale_color_manual("Significance",values = c("p05"="olivedrab4","p > 0.05"="seashell3"))  + 
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1)) 

x.aldexmhc2 <- x.aldex %>% dplyr::select(contains("MHC2_Diversity")) %>% rownames_to_column("OTU") %>% rename(pval = MHC2_Diversity.pval,) %>% mutate(direct = case_when(MHC2_Diversity.Est < 0 & pval < 0.05 ~ "p05", MHC2_Diversity.Est > 0 & pval < 0.05 ~ "p05", TRUE ~ "p > 0.05"),cog = case_when(pval < 0.05 ~ OTU)) %>% mutate(taxonsig = case_when(pval < 0.05 ~ OTU)) %>% merge(.,mpataxo7spec, by = "taxonsig", all.x=T)%>% mutate(Species2=str_extract(Species, "(?<=__).*"))

mhc2mpaaldexplot <- ggplot(x.aldexmhc2, aes(x= MHC2_Diversity.Est,y=  -log(pval),color=direct, label=Species2)) +
  geom_point(size =3) +
  geom_text_repel(size=3, colour="black") +
  ggtitle("C. MHC-II Diversity")+
  xlab("Log fold change MHC-II Diversity") +
  ylab("-log(p-value)")+
  scale_color_manual("Significance",values = c("p05"="olivedrab4","p > 0.05"="seashell3"),labels = c("p > 0.05", "p < 0.05")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

x.aldexaseua7 <- x.aldex %>% dplyr::select(contains("Aseua7")) %>% rownames_to_column("OTU") %>% rename(pval = Aseua71.pval,) %>% mutate(direct = case_when(Aseua71.Est < 0 & pval < 0.05 ~ "p05", Aseua71.Est > 0 & pval < 0.05 ~ "p05", TRUE ~ "p > 0.05"),cog = case_when(pval < 0.05 ~ OTU)) %>% mutate(taxonsig = case_when(pval < 0.05 ~ OTU)) %>% merge(.,mpataxo7spec, by = "taxonsig", all.x=T) %>% mutate(Species2=str_extract(Species, "(?<=__).*"))

Aseua7mpaaldexplot <- ggplot(x.aldexaseua7, aes(x= Aseua71.Est,y=  -log(pval),color=direct, label=Species2)) +
  geom_point(size =3) +
  geom_text_repel(size=3, colour="black") +
  ggtitle("D. Aseua7")+
  xlab("Log fold change Aseua7") +
  ylab("-log(p-value)")+
  scale_color_manual("Significance",values = c("p05"="olivedrab4","p > 0.05"="seashell3"),labels = c("p > 0.05", "p < 0.05")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))



plotlegend = get_legend(Aseua7mpaaldexplot)

x.aldexmhc1 <- x.aldex %>% dplyr::select(contains("MHC1_Diversity")) %>% rownames_to_column("OTU") %>% rename(pval = MHC1_Diversity.pval,) %>% mutate(direct = case_when(MHC1_Diversity.Est < 0 & pval < 0.05 ~ "p05", MHC1_Diversity.Est > 0 & pval < 0.05 ~ "p05", TRUE ~ "p > 0.05"),cog = case_when(pval < 0.05 ~ OTU)) %>% mutate(taxonsig = case_when(pval < 0.05 ~ OTU)) %>% merge(.,mpataxo7spec, by = "taxonsig", all.x=T) %>% mutate(Species2=str_extract(Species, "(?<=__).*"))

mhc1mpaaldexplot <- ggplot(x.aldexmhc1, aes(x= MHC1_Diversity.Est,y=  -log(pval),color=direct, label=Species2)) +
  geom_point(size =3) +
  geom_text_repel(size=3, colour="black") +
  ggtitle("B. MHC-I Diversity")+
  xlab("Log fold change MHC-I Diversity") +
  ylab("-log(p-value)")+
  scale_color_manual("Significance",values = c("p05"="olivedrab4","p > 0.05"="seashell3"),labels = c("p > 0.05", "p < 0.05")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

tiff("Output/mpa_ALDEx2.tiff", res =400, units = "in", width = 10, height = 10, bg = "transparent")
ggarrange(hsobsmpaaldexplot,mhc1mpaaldexplot,mhc2mpaaldexplot,Aseua7mpaaldexplot, common.legend = T, legend.grob = plotlegend, legend = "bottom")
dev.off()


tiff("Output/mpa_ANCOM_GLLVM_ALDEx2.tiff", res =400, units = "in", width = 15, height = 9, bg = "transparent")
ggarrange(hsobsmpaancomplot,hsobsmpagllvmplot,hsobsmpaaldexplot,aseua7mpaancomplot,Aseua7mpagllvmplot,Aseua7mpaaldexplot,common.legend = T,legend = "right")
dev.off()



# NOG ----
NOGTLR <- readRDS("InputTables/NOGTLR.rds")

## NOG alpha ----
NOGrarefy <- rarefy_even_depth(NOGTLR, rngseed=88,sample.size = 100000, replace=F)

NOGTLRalpha <- alpha_estimate(NOGrarefy) %>% mutate_at(vars(contains('Ase')), as.factor) %>% mutate(cutage = pmin(SamplingAge,10))

ggplot(NOGTLRalpha, aes(x=MHC1_Diversity)) + geom_histogram()
ggplot(NOGTLRalpha, aes(x=MHC2_Diversity)) + geom_histogram()
ggplot(NOGTLRalpha, aes(x=MHC1_Diversity, y= Lifespan)) + geom_point() + geom_smooth()

#AsePC <- prcomp(~ as.numeric(Asedab3)+as.numeric(Asedab4)+as.numeric(Asedab5) + as.numeric(Aseua1)+as.numeric(Aseua3)+as.numeric(Aseua4)+as.numeric(Aseua5)+as.numeric(Aseua6)+as.numeric(Aseua7)+ as.numeric(Aseua8)+ as.numeric(Aseua11)  + as.numeric(Aseua10), data=NOGTLRalpha)
#plot(AsePC)
#summary(AsePC)
#AsePCdf <- as.data.frame(AsePC[["x"]]) %>% cbind(NOGTLRalpha)


ggplot(NOGTLRalpha, aes(x=exp(arm::rescale(Observed)))) + geom_histogram()
NOGobs <- lmer(exp(arm::rescale(Observed)) ~  Aseua11  + cutage+ season + SexEstimate + Timeinfridge_scaled + CatchTime_scaled + SampleYear  + (1|BirdID), data=NOGTLRalpha )
summary(NOGobs)
car::Anova(NOGobs,type="III")
vif(NOGobs)
simulateResiduals(NOGobs,plot=T)
r.squaredGLMM(NOGobs)


NOGobs2 <- lmer(exp(arm::rescale(Observed)) ~ Hs_obs + MHC1_Diversity + MHC2_Diversity + cutage+ season + SexEstimate + Timeinfridge_scaled + CatchTime_scaled + SampleYear  + (1|BirdID), data=NOGTLRalpha )
summary(NOGobs2)
r.squaredGLMM(NOGobs2)
car::Anova(NOGobs2, type="III",test.statistic = "Chisq")

# NOGobs3 <- lmer(exp(arm::rescale(Observed)) ~ PC1 + PC2 + PC3 + PC4 +SamplingAge2  + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + SampleYear + (1|BirdID), data=AsePCdf )
# summary(NOGobs3)
# r.squaredGLMM(NOGobs3)

NOGobs2data <- ggpredict(NOGobs2, terms = c("MHC1_Diversity"))
plot(NOGobs2data)

ggplot(NOGTLRalpha, aes(x=Shannon)) + geom_histogram()
NOGshan <- lmer(exp(Shannon) ~ Aseua11 + cutage+ season + SexEstimate + Timeinfridge_scaled + CatchTime_scaled + SampleYear  + (1|BirdID), data=NOGTLRalpha )
summary(NOGshan)
vif(NOGshan)
simulateResiduals(NOGshan,plot=T)
car::Anova(NOGshan,type="III")

# NOGshandata <- ggpredict(NOGshan,terms = c("SamplingAge2","Asedab4")) 
# plot(NOGshandata)

NOGshan2 <- lmer(exp(Shannon) ~ Hs_obs + MHC1_Diversity + MHC2_Diversity + cutage+ season + SexEstimate + Timeinfridge_scaled + CatchTime_scaled + SampleYear + + (1|BirdID), data=NOGTLRalpha )
summary(NOGshan2)
simulateResiduals(NOGshan2,plot = T)
car::Anova(NOGshan2,type="III")

NOGTLRalpha %>% group_by(MHC2_Diversity) %>% summarise(minage = min(SamplingAge2),maxage = max(SamplingAge2),count=n())

## NOG beta ----

NOGTLRDB <- aggregate_taxa(NOGTLR,"Class")
plot_composition(NOGTLRDB,group_by = "Aseua11")

NOGTLRrare <- transform(NOGTLR,"compositional") %>% core_members(., detection = 0.0001, prevalence = 0.05)
NOGTLRCLR <- prune_taxa(NOGTLRrare,NOGTLR) %>% transform(.,"clr") %>% mutate_sample_data(MHC1Dbin = as.factor(case_when(MHC1_Diversity < 4 ~ "<4", MHC1_Diversity <7 ~ "4-6",MHC1_Diversity >= 7 ~ "7")))

NOGTLRCLRord <- ordinate(NOGTLRCLR, method="RDA",distance = "euclidean")

plot_ordination(NOGTLRCLR,NOGTLRCLRord, color="Aseua11") + stat_ellipse()
plot_ordination(NOGTLRCLR,NOGTLRCLRord, color="MHC1Dbin") + stat_ellipse(level = 0.95,type = "euclid")

plot_ordination(NOGTLRCLR,NOGTLRCLRord, axes = 1:2, color="MHC1_Diversity") +
  scale_color_gradient2(
    low = "olivedrab2", 
    mid = "deepskyblue2", 
    high = "orangered3", 
    midpoint = 4.5) + 
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))

#mhc1bin
mhc1binNOGdf <- plot_ordination(NOGTLRCLR,NOGTLRCLRord,axes = 1:2 ,color="MHC1Dbin",justDF = T) %>% group_by(MHC1Dbin) %>% dplyr::summarise(meanPC1 = mean(PC1), meanPC2 = mean(PC2))
mhc1binNOGplot <- plot_ordination(NOGTLRCLR,NOGTLRCLRord,axes = 1:2 ,color="MHC1Dbin") + 
  geom_point(data=mhc1binNOGdf,aes(x=meanPC1,y=meanPC2,fill=MHC1Dbin), shape=23, stroke=1, size=5, colour="black") +
  scale_colour_manual("MHC-I Diversity",values = c("olivedrab2","deepskyblue2","orangered3"))+
  scale_fill_manual("MHC-I Diversity",values = c("olivedrab2","deepskyblue2","orangered3"))+ 
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1))


NOGTLRCLRmat <-vegan_otu(NOGTLRCLR)
NOGTLRCLR_st <- as(sample_data(NOGTLRCLR),"data.frame") %>% 
  mutate(survive = as.factor(survive)) %>% group_by(BirdID) %>%
  mutate(deltaAge = capage - mean(capage), meanAge = mean(capage),TerminalYearbird = case_when(any(TerminalYear %in% "1") ~ "1", TRUE ~ "0")) %>% dplyr::select(SamplingAge,TerminalYear, season, SampleYear,SexEstimate, Timeinfridge_scaled,CatchTime_scaled,TQcorrected_scaled,BirdID,contains("Ase"), MHC1_Diversity,MHC2_Diversity,Hs_obs,capage) %>% mutate(SamplingAge2 = round(SamplingAge), cutage = pmin(SamplingAge,10))

ggplot(NOGTLRCLR_st, aes(x=SamplingAge,y=MHC1_Diversity)) + geom_point()

perm <- how(nperm = 9999, blocks = NOGTLRCLR_st$BirdID)
set.seed(888)

NOGTLRCLRpermS<- adonis2(NOGTLRCLRmat ~ Aseua5+Aseua7+Aseua9 + cutage+ season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
NOGTLRCLRpermS

NOGTLRCLRperm3<- adonis2(NOGTLRCLRmat ~ Hs_obs+MHC1_Diversity +  MHC2_Diversity + cutage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
NOGTLRCLRperm3

## NOG DAA ----
Chuen_new_COG_tax <- read.csv("~/Downloads/cog_all_results.tsv",sep="\t") %>% distinct(COG,.keep_all = T)
# nog daa but with core cog
NOGTLRcoremembers <- transform(NOGTLR,"compositional") %>% core_members(.,detection = 0.001, prevalence = 0.5)
NOGTLRCore <- prune_taxa(NOGTLRcoremembers,NOGTLR) 
# aldex2 nog
eggNOG_Cat_description <- read_csv("InputTables/eggNOG_Cat_description.csv")
NOGTLRcoremat <- round(data.frame(otu_table(NOGTLRCore)))
NOGTLRcorest <- data.frame(sample_data(NOGTLRCore))
mmnog <- model.matrix(~  Hs_obs + MHC1_Diversity + MHC2_Diversity + Timeinfridge_scaled, NOGTLRcorest) 
x.aldexnog <- ALDEx2::aldex(NOGTLRcoremat, mmnog, mc.samples=200, test="glm", effect=TRUE, include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL)

x.aldexnogMHC1 <- x.aldexnog %>% dplyr::select(contains("MHC1_Diversity")) %>% rownames_to_column("OTU") %>% rename(pval = MHC1_Diversity.pval,) %>% mutate(direct = case_when(MHC1_Diversity.Est < 0 & pval < 0.05 ~ "Negative_LFC", MHC1_Diversity.Est > 0 & pval < 0.05 ~ "Positive_LFC", TRUE ~ "p > 0.05"),cog = case_when(pval < 0.05 ~ OTU)) %>% mutate(taxonsig = case_when(pval < 0.1 ~ OTU)) %>% merge(.,Chuen_new_COG_tax,by.x="cog",by.y="COG",all.x=T ) 

x.aldexnogMHC1 %>% filter(pval < 0.05)

x.aldexnogMHC1sig <- x.aldexnogMHC1 %>% filter(pval < 0.05) %>% dplyr::select(cog, direct, Annotation, Cat)

write_csv(x.aldexnogMHC1sig, "Output/aldexNOG_MHC1_sig.csv")

mhc1nogaldexplot  <- ggplot(x.aldexnogMHC1, aes(x= MHC1_Diversity.Est,y=  -log(pval),color=direct, label=cog)) +
  geom_point(size =3) +
  geom_text_repel(size=3, colour="black") +
  ggtitle("")+
  xlab("Log fold change MHC-I Diversity") +
  ylab("-log(p-value)")+
  scale_color_manual("Direction",values = c("Positive_LFC"="deepskyblue2","Negative_LFC"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1),plot.title = element_text(hjust = 0.5))

# cat nog
x.aldexnogMHC1_cat <- x.aldexnogMHC1 %>% dplyr::select(cog,Cat,direct) %>% mutate(Cat = case_when(is.na(Cat)~"`",Cat == "" ~ "`" ,TRUE ~ Cat)) %>% separate_longer_position(Cat,1) %>% filter(grepl("^[A-Za-z]", Cat)) %>% group_by(Cat,direct) %>% summarise(countaldex2 = n()) %>% mutate(realcount = case_when(direct %in% "Negative_LFC" ~ -countaldex2, direct %in% "Positive_LFC" ~countaldex2 )) %>% merge(.,eggNOG_Cat_description,by="Cat",all.x=T) %>% arrange(desc(realcount)) %>% mutate(groupid =case_when(Cat %in% "V" ~ 1,Cat %in% "M" ~ 2,Cat %in% "T" ~ 3,Cat %in% "R" ~ 4,Cat %in% "K" ~ 5,Cat %in% "L" ~ 6,Cat %in% "O" ~ 7,Cat %in% "J" ~ 8,Cat %in% "H" ~ 9,Cat %in% "P" ~ 10,Cat %in% "I" ~ 11,Cat %in% "E" ~ 12,Cat %in% "G" ~ 13 ))

mhc1nogaldexplotB <- ggplot(x.aldexnogMHC1_cat, aes(x=realcount,y=reorder(COG_Description, -groupid), fill = direct)) + 
  geom_bar(stat = "identity") + 
  xlab("Number of functional genes") + 
  ylab("Functional categories") +
  ggtitle("")+
  scale_fill_manual("Direction",values = c("Positive_LFC"="deepskyblue2","Negative_LFC"="orangered3"))+
  theme_tufte(base_size = 15, base_family = "Arial")+ 
  theme(axis.line = element_line(colour = "black", linetype=1),plot.title = element_text(hjust = 1))

tiff("Output/ALDEx2_MHC1.tiff", res =400, units = "in", width = 10, height = 10, bg = "transparent")
ggarrange(ggarrange(mhc1nogaldexplot,NULL, widths = c(4,0.5), legend = "top"),mhc1nogaldexplotB, nrow = 2, labels = c("A.","B."), legend =F)
dev.off()

