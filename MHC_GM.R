# Host Genetics and Microbiome 
# How do we do this? What do we want to compare?
# TLR, MHC, Inbreeding, Telomere length, Methylation, 
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

# Functional eggNOG phyloseq
otuNOG <- read.delim("~/Documents/PhD/R_analysis/Phyloseq/InputTables/MFF_Feb2024/emapper/EM.NOGL0.txt", row.names=1, sep = "\t")
#load tax table
taxNOG <- read.csv("~/Documents/PhD/R_analysis/Phyloseq/InputTables/NOG.Chuen.annotations.csv",row.names=1)
#load sample table (needs editing)
st <- read.csv("~/Documents/PhD/R_analysis/SampleMetadata/InputTables/st.csv")
# sample data 
st$TubeNumber <- sub("^","SW",st$TubeNumber)
# add TubeNumber as rownames 
stA <- st %>% filter(!duplicated(TubeNumber)) %>% mutate(deathdate = LatestEligibleDate) %>% mutate(SampleYear = format(as.Date(SampleDate, format="%Y-%m-%d"), "%Y")) 
st3<-as.data.frame(stA)
st4 <- st3[,-1]
rownames(st4) <- st3[,1]
stdata <- sample_data(st4)
st4 #from taxo

#convert otu and tax tables to matrix
NOGotumat <- as.matrix(otuNOG)
NOGtaxmat <- as.matrix(taxNOG)

NOGOTU <- otu_table(NOGotumat, taxa_are_rows = TRUE)
NOGTAX = tax_table(NOGtaxmat)
stdata <- sample_data(st4)
NOGp = phyloseq(NOGOTU, NOGTAX, stdata)

NOGpo <- NOGp %>% mutate_sample_data(Type = case_when(SampleType %in% c("F","f") ~ "F", TRUE ~ "C"))
NOGpf <- NOGpo %>% filter_sample_data(Type=="F") %>% filter_sample_data(SamplingAge > 0.5) %>% mutate_sample_data(SamplingAge2 = round(SamplingAge)) %>% mutate_sample_data(AgeClass = case_when(SamplingAge2 <= 5 ~ "A", SamplingAge2 <= 10 ~ "B", SamplingAge2 > 10 ~ "C"), capage=pmin(SamplingAge,12)) %>% mutate_sample_data(Age_scaled = arm::rescale(SamplingAge,"center"),Timeinfridge_scaled = arm::rescale(Timeinfridge, "center"),CatchTime_scaled = arm::rescale(CatchTime, "center"),TQcorrected_scaled = arm::rescale(TQcorrected, "center")) 

NOGpfsampledata <- data.frame(sample_data(NOGpf)) %>% rownames_to_column("Rownames")

NOGpfTLR <- merge(NOGpfsampledata,TL3_MHC_20240530full, by="BirdID") %>% merge(.,SWhet, by="BirdID",all.x=T) %>% column_to_rownames("Rownames")

NOGTLR <- phyloseq(tax_table(NOGpf),otu_table(NOGpf),sample_data(NOGpfTLR)) 

# NOG alpha ----
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

# NOGshan2data <- ggpredict(NOGshan2, terms=c("SamplingAge2 [all]","MHC2_Diversity [all]"), back_transform = T)
# plot(NOGshan2data)
# NOGshan2data2 <- NOGshan2data %>% dplyr::rename(capage=x,Shannon=predicted,MHC2_Diversity=group) %>% filter((MHC2_Diversity %in% "1" & capage >= 1 & capage <= 5) | (MHC2_Diversity %in% "2" & capage >= 1 & capage <= 9)| (MHC2_Diversity %in% "3" & capage >= 1 & capage <= 15)| (MHC2_Diversity %in% "4" & capage >= 1 & capage <= 16)| (MHC2_Diversity %in% "5" & capage >= 1& capage <= 9)) 
# NOGshan2data3 <- as.data.frame(NOGshan2data2) %>% mutate(Shannon=log(Shannon),conf.low=log(abs(conf.low)),conf.high=log(conf.high))

#ggplot(NOGshan2data3,aes(x=capage,y=Shannon,colour=as.character(MHC2_Diversity)) ) + geom_line() + geom_ribbon(aes(ymin=conf.low,ymax=conf.high, fill =as.character(MHC2_Diversity) ),alpha = 0.1, colour=FALSE) + geom_point(data=NOGTLRalpha,aes(x=capage,y=Shannon,colour=as.character(MHC2_Diversity)), inherit.aes = F) + geom_line(data=NOGTLRalpha,aes(x=capage,y=Shannon,colour=as.character(MHC2_Diversity),group=BirdID), inherit.aes = F, colour="grey",linetype=2) + xlab("Age (years)") + scale_colour_colorblind(name = "MHC 2 Diversity") +scale_fill_colorblind(name = "MHC 2 Diversity") + theme_tufte(base_size = 15, base_family = "Arial") + theme(axis.line = element_line(colour = "black", linetype=1)) + facet_wrap(vars(as.character(MHC2_Diversity))) + geom_hline(yintercept = 6)

#ggplot(NOGshan2data2,aes(x=capage,y=Shannon,colour=as.character(MHC2_Diversity)) ) + geom_line() + geom_ribbon(aes(ymin=conf.low,ymax=conf.high, fill =as.character(MHC2_Diversity) ),alpha = 0.1, colour=FALSE) + geom_point(data=NOGTLRalpha,aes(x=SamplingAge2,y=exp(Shannon),colour=as.character(MHC2_Diversity)), inherit.aes = F) + geom_line(data=NOGTLRalpha,aes(x=SamplingAge2,y=exp(Shannon),colour=as.character(MHC2_Diversity),group=BirdID), inherit.aes = F, colour="grey",linetype=2) + xlab("Age (years)")+ ylab("Exponential transformed Shannon") + scale_colour_colorblind(name = "MHC 2 Diversity") +scale_fill_colorblind(name = "MHC 2 Diversity") + theme_tufte(base_size = 15, base_family = "Arial") + theme(axis.line = element_line(colour = "black", linetype=1))

# NOGshan3 <- lmer(exp(Shannon) ~PC1 + PC2 + PC3 + PC4 +SamplingAge2  + season + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + SampleYear + (1|BirdID), data=AsePCdf )
# summary(NOGshan3)
# r.squaredGLMM(NOGshan3)

# NOG beta ----

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
#NOGTLRCLRperm0<- adonis2(NOGTLRCLRmat ~ SamplingAge + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled  , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
#NOGTLRCLRperm0

#NOGTLRCLRperm<- adonis2(NOGTLRCLRmat ~ Aseua5+Aseua7+Aseua9 + SamplingAge + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled  , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
#NOGTLRCLRperm

#NOGTLRCLRperm2<- adonis2(NOGTLRCLRmat ~Hs_obs+Aseua1+ Aseua5+Aseua7+Aseua9 + SamplingAge + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled  , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
#NOGTLRCLRperm2

NOGTLRCLRpermS<- adonis2(NOGTLRCLRmat ~ Aseua5+Aseua7+Aseua9 + cutage+ season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
NOGTLRCLRpermS

# NOGTLRCLRperm2 <- adonis2(NOGTLRCLRmat ~ Aseua7+Aseua9 + SamplingAge2 + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled  , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# NOGTLRCLRperm2

# NOGTLRCLRperm2<- adonis2(NOGTLRCLRmat ~ capage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + Asedab3+Asedab4+Asedab5 + Aseua1+Aseua3+Aseua4+Aseua5+Aseua6+Aseua7+ Aseua8+ Aseua11  , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# NOGTLRCLRperm2

NOGTLRCLRperm3<- adonis2(NOGTLRCLRmat ~ Hs_obs+MHC1_Diversity +  MHC2_Diversity + cutage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
NOGTLRCLRperm3

NOGTLRCLRperm4<- adonis2(NOGTLRCLRmat ~Hs_obs+ MHC1_Diversity + I(MHC1_Diversity^2) +  MHC2_Diversity + SamplingAge + season + SampleYear  + SexEstimate+ Timeinfridge_scaled +CatchTime_scaled , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
NOGTLRCLRperm4

# NOGTLRCLRperm3A<- adonis2(NOGTLRCLRmat ~ MHC1_Diversity +  MHC2_Diversity + Timeinfridge_scaled , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# NOGTLRCLRperm3A
# 
# NOGTLRCLRperm4A<- adonis2(NOGTLRCLRmat ~ MHC1_Diversity + I(MHC1_Diversity^2) +  MHC2_Diversity + Timeinfridge_scaled  , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# NOGTLRCLRperm4A

NOGTLRCLRperm5<- adonis2(NOGTLRCLRmat ~ Hs_obs+MHC1_Diversity + I(MHC1_Diversity^2) +  MHC2_Diversity+ I(MHC2_Diversity^2) + season + SampleYear  + SexEstimate+ Timeinfridge_scaled   , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
NOGTLRCLRperm5

# NOGTLRCLRperm4<- adonis2(NOGTLRCLRmat ~ capage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + MHC1_Diversity*SexEstimate + MHC2_Diversity*SexEstimate  , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# NOGTLRCLRperm4
# 
# NOGTLRCLRperm5<- adonis2(NOGTLRCLRmat ~ Aseua4 + Aseua8+Aseua9+ Aseua11 + SamplingAge2 + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled  , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# NOGTLRCLRperm5
# 
# NOGMHCdiversityCLRperm<- adonis2(NOGTLRCLRmat ~ MHC1_Diversity + MHC2_Diversity + capage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled  , data=NOGTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# NOGMHCdiversityCLRperm
# 
# NOGTLRCLR_stmatrix <- as(sample_data(NOGTLRCLR),"data.frame")  %>% dplyr::select(contains("Ase"))
# d = dist(NOGTLRCLR_stmatrix, method = "binary")
# NOGMHCmatrix <- as.matrix(d)
# 
# NOGmicrobiomematrix <- as.matrix(dist(t(data.frame(otu_table(NOGTLRCLR))), method = "euclidean"))
# 
# mantel(NOGMHCmatrix,NOGmicrobiomematrix)


### NOG DAA ----
Chuen_new_COG_tax <- read.csv("~/Downloads/cog_all_results.tsv",sep="\t") %>% distinct(COG,.keep_all = T)
# nog daa but with core cog
NOGTLRcoremembers <- transform(NOGTLR,"compositional") %>% core_members(.,detection = 0.001, prevalence = 0.5)
NOGTLRCore <- prune_taxa(NOGTLRcoremembers,NOGTLR) 

#NOGcoreancomid <- ancombc2(NOGTLRCore,fix_formula = "MHC1_Diversity + Timeinfridge_scaled",rand_formula = "(1|BirdID)" , p_adj_method = "BH")
#NOGcoreancomidres <- NOGcoreancomid$res

#NOGcoreancomidresmhc1 <- NOGcoreancomidres %>% dplyr::select(taxon, contains("MHC1")) %>% mutate(direct = case_when(lfc_MHC1_Diversity < 0 & p_MHC1_Diversity < 0.05 ~ "Negative_LFC", lfc_MHC1_Diversity > 0 & p_MHC1_Diversity < 0.05 ~ "Positive_LFC", TRUE ~ "p > 0.05"), cog = case_when(p_MHC1_Diversity < 0.05 ~ taxon)) %>% merge(.,Chuen_new_COG_tax,by.x="cog",by.y="COG",all.x=T ) 

# NOGcoreancomidresmhc1plot <- ggplot(NOGcoreancomidresmhc1, aes(x=lfc_MHC1_Diversity,y= reorder(taxon, +lfc_MHC1_Diversity),color=direct)) +
#   geom_point() +
#   geom_errorbar(aes(xmin = lfc_MHC1_Diversity - 1.96*se_MHC1_Diversity, xmax = lfc_MHC1_Diversity + 1.96*se_MHC1_Diversity)) +
#   geom_vline(xintercept = 0, linetype = 3) +
#   ggtitle("A. ANCOMBC2 MHC-I Diversity")+
#   ylab("") +
#   xlab("Log fold change with MHC-I Diversity") +
#   scale_color_manual(values = c("Positive_LFC"="olivedrab4","Negative_LFC"="orangered3","p > 0.05"="seashell3"), name="Significance") +
#   theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

#mhc1nogancomplot <- ggplot(NOGcoreancomidresmhc1, aes(x= lfc_MHC1_Diversity,y=  -log(p_MHC1_Diversity),color=direct, label=cog)) +
  # geom_point(size=3) +
  # geom_text_repel(size=3, colour="black") +
  # ggtitle("A. ANCOMBC2 MHC-I Diversity")+
  # xlab("Log fold change MHC-I Diversity") +
  # ylab("-log(p-value)") +
  # scale_color_manual("Direction",values = c("Positive_LFC"="olivedrab","Negative_LFC"="orangered3","p > 0.05"="seashell3")) +
  # theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

#NOGcoreancomidresmhc1_cat <- NOGcoreancomidresmhc1 %>% dplyr::select(cog,Cat) %>% mutate(Cat = case_when(is.na(Cat)~"`",Cat == "" ~ "`" ,TRUE ~ Cat)) %>% separate_longer_position(Cat,1) %>% filter(grepl("^[A-Za-z]", Cat)) %>% group_by(Cat) %>% summarise(countancom = n())

# NOGcoreancomidreshs <- NOGcoreancomidres %>% dplyr::select(taxon, contains("Hs_obs")) %>% mutate(direct = case_when(lfc_Hs_obs < 0 & p_Hs_obs < 0.05 ~ "Negative_LFC", lfc_Hs_obs > 0 & p_Hs_obs < 0.05 ~ "Positive_LFC", TRUE ~ "p > 0.05"), cog = case_when(p_Hs_obs < 0.05 ~ taxon))
# 
# hsobsnogancomplot <- ggplot(NOGcoreancomidreshs, aes(x= lfc_Hs_obs,y=  -log(p_Hs_obs),color=direct, label=cog)) +
#   geom_point() +
#   geom_text_repel(size=3, colour="black") +
#   ggtitle("A. ANCOMBC2 Heterozygosity")+
#   xlab("Log fold change with Heterozygosity") +
#   ylab("-log(p-value)") +
#   scale_color_manual("Direction",values = c("Positive_LFC"="olivedrab4","Negative_LFC"="orangered3","p > 0.05"="seashell3")) +
#   theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))
# 
# NOGcoreancomid2 <- ancombc2(NOGTLRCore,fix_formula = "Aseua7+Aseua9 + SamplingAge2 + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled",rand_formula = "(1|BirdID)" , p_adj_method = "BH")
# NOGcoreancomidres2 <- NOGcoreancomid2$res

# NOG DAA gllvm
# NOGTLRCoreclr <- transform(NOGTLRCore,"clr")
# NOGTLRCoreclrmat <- vegan_otu(NOGTLRCoreclr)
# NOGTLRCoreclrst <- as(sample_data(NOGTLRCoreclr),"data.frame") %>% 
#   mutate(TerminalYear = as.factor(TerminalYear)) %>% 
#   mutate(capage=pmin(SamplingAge,12))%>%     
#   dplyr::select(MHC1_Diversity, MHC2_Diversity,TerminalYear,BirdID,season,SampleYear,SamplingAge2, SexEstimate,CatchTime_scaled,Timeinfridge_scaled,TQcorrected_scaled,Hs_obs) 
# NOGTLRCoreclrgllvm <- gllvm::gllvm(NOGTLRCoreclrmat, NOGTLRCoreclrst, formula = ~  MHC1_Diversity + Timeinfridge_scaled , family = "gaussian",num.lv = 1,sd.errors = TRUE, row.eff = ~(1|BirdID), seed = 8,starting.val="zero")
# gllvm::coefplot.gllvm(NOGTLRCoreclrgllvm,which.Xcoef = 1:2)
# summary(NOGTLRCoreclrgllvm)
# NOGTLRCoreclrgllvmsum <- summary(NOGTLRCoreclrgllvm)
# NOGTLRCoreclrgllvmcoef <- as.data.frame(NOGTLRCoreclrgllvmsum[["Coef.tableX"]]) %>% mutate(terms= row.names(.)) %>% separate_wider_regex(cols = terms, c(terms = ".*", ":", OTU = ".*"))
# gllvm::AICc.gllvm(NOGTLRCoreclrgllvm)
# 
# NOGXgllvmcoefmhc1 <- NOGTLRCoreclrgllvmcoef %>% filter(terms %in% "MHC1_Diversity") %>% rename(pval = 'Pr(>|z|)', se = 'Std. Error') %>% mutate(direct = case_when(Estimate < 0 & pval < 0.05 ~ "Negative_LFC", Estimate > 0 & pval < 0.05 ~ "Positive_LFC", TRUE ~ "p > 0.05"),cog = case_when(pval < 0.05 ~ OTU)) %>% merge(.,Chuen_new_COG_tax,by.x="cog",by.y="COG",all.x=T ) 
# 
# mhc1noggllvmplot <- ggplot(NOGXgllvmcoefmhc1, aes(x= Estimate,y=  -log(pval),color=direct, label=cog)) +
#   geom_point(size =3) +
#   geom_text_repel(size=3, colour="black") +
#   ggtitle("B. GLLVM MHC-I Diversity")+
#   xlab("Log fold change MHC-I Diversity") +
#   ylab("-log(p-value)")+
#   scale_color_manual("Direction",values = c("Positive_LFC"="olivedrab4","Negative_LFC"="orangered3","p > 0.05"="seashell3")) +
#   theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))
# 
# ggarrange(mhc1nogancomplot,mhc1noggllvmplot, labels = c("A","B"),common.legend = T,legend = "right")

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

NOGXgllvmcoefmhc1_cat <- NOGXgllvmcoefmhc1 %>% dplyr::select(cog,Cat) %>% mutate(Cat = case_when(is.na(Cat)~"`",Cat == "" ~ "`" ,TRUE ~ Cat)) %>% separate_longer_position(Cat,1) %>% filter(grepl("^[A-Za-z]", Cat)) %>% group_by(Cat) %>% summarise(countgllvm = n())

catmhc1 <- merge(NOGcoreancomidresmhc1_cat,NOGXgllvmcoefmhc1_cat,by="Cat",all=T) %>% merge(.,x.aldexnogMHC1_cat,all=T) %>% pivot_longer(cols = starts_with("count"), names_to = "method") %>% group_by(Cat) %>% mutate(sumvalue = sum(value,na.rm = T))

catmhc1plot <- ggplot(catmhc1, aes(x=value,y=reorder(Cat, -sumvalue), fill=method)) + 
  geom_bar(stat = "identity") + 
  xlab("Number of eggNOG categories") + 
  ylab("eggNOG categories")+
  ggtitle("E. eggNOG categories with MHC-I Diversity")+
  scale_fill_manual("Method",values = c("countaldex2"="skyblue","countancom"="olivedrab4","countgllvm"="orangered3"), labels = c("ALDEx2","ANCOMBC2","GLLVM"), name = "Method") +
  theme_tufte(base_size = 15, base_family = "Arial")+ 
  theme(axis.line = element_line(colour = "black", linetype=1))

NOGXgllvmcoefhs <- NOGTLRCoreclrgllvmcoef %>% filter(terms %in% "Hs_obs") %>% rename(pval = 'Pr(>|z|)', se = 'Std. Error') %>% mutate(direct = case_when(Estimate < 0 & pval < 0.05 ~ "Negative_LFC", Estimate > 0 & pval < 0.05 ~ "Positive_LFC", TRUE ~ "p > 0.05"),cog = case_when(pval < 0.05 ~ OTU))

hsnoggllvmplot <- ggplot(NOGXgllvmcoefhs, aes(x= Estimate,y=  -log(pval),color=direct, label=cog)) +
  geom_point() +
  geom_text_repel(size=3, colour="black") +
  ggtitle("B. GLLVM Heterozygosity")+
  xlab("Log fold change Heterozygosity") +
  ylab("-log(p-value)")+
  scale_color_manual("Direction",values = c("Positive_LFC"="olivedrab4","Negative_LFC"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

# venn diagram of DAA
NOGXgllvmcoefmhc1_data <- NOGXgllvmcoefmhc1 %>% dplyr::select(direct,cog) %>% mutate(method="GLLVM") %>% filter(!is.na(cog))

NOGcoreancomidresmhc1_data <- NOGcoreancomidresmhc1 %>% dplyr::select(direct,cog) %>% mutate(method="ANCOMBC2") %>% filter(!is.na(cog)) %>% data.frame()

nogaldex2_data <- x.aldexnogMHC1 %>% dplyr::select(direct,cog) %>% mutate(method="ALDEx2") %>% filter(!is.na(cog)) %>% data.frame()

noggllvmancom <- list(ANCOMBC2=NOGcoreancomidresmhc1_data$cog,GLLVM=NOGXgllvmcoefmhc1_data$cog, ALDEx2=nogaldex2_data$cog)



nogmhc1venn <- ggVennDiagram(noggllvmancom,set_color = "black") + scale_fill_gradient2(low="white",mid="white",high = "white")  +
  ggtitle("D. ANCOMBC2 vs GLLVM vs ALDEx2 MHC-I Diversity")+
  theme_tufte(base_size = 15, base_family = "Arial")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.1),legend.position="none")

tiff("Output/NOG_ANCOM_GLLVM_ALDEx2_MHC1.tiff", res =400, units = "in", width = 15, height = 9, bg = "transparent")
ggarrange(ggarrange(mhc1nogancomplot,mhc1noggllvmplot,mhc1nogaldexplot, common.legend = T, legend = "right",nrow = 1), ggarrange(nogmhc1venn,NULL,catmhc1plot, nrow = 1, widths = c(1,0.3,1)), ncol = 1)
dev.off()

noggllvmancom2 <- merge(NOGXgllvmcoefmhc1_data,NOGcoreancomidresmhc1_data,by=c("direct","cog"),all=T) %>% mutate(tests = case_when(!is.na(method.x) & !is.na(method.y) ~ "both", !is.na(method.x) & is.na(method.y) ~ method.x, is.na(method.x) & !is.na(method.y) ~ method.y, TRUE ~ "none"))

noggllvmancom3 <- noggllvmancom2 %>% filter(tests %in% "both") %>% merge(.,Chuen_new_COG_tax,by.x="cog",by.y="COG",all.x=T )
noggllvmancom3$cog

# hs obs cat

NOGcoreancomidreshs_cat <- NOGcoreancomidreshs %>% merge(.,Chuen_new_COG_tax,by.x="cog",by.y="COG",all.x=T ) %>% dplyr::select(cog,Cat) %>% mutate(Cat = case_when(is.na(Cat)~"`",Cat == "" ~ "`" ,TRUE ~ Cat)) %>% separate_longer_position(Cat,1) %>% filter(grepl("^[A-Za-z]", Cat)) %>% group_by(Cat) %>% summarise(countancom = n())

NOGXgllvmcoefhs_cat <- NOGXgllvmcoefhs %>% merge(.,Chuen_new_COG_tax,by.x="cog",by.y="COG",all.x=T ) %>% dplyr::select(cog,Cat) %>% mutate(Cat = case_when(is.na(Cat)~"`",Cat == "" ~ "`" ,TRUE ~ Cat)) %>% separate_longer_position(Cat,1) %>% filter(grepl("^[A-Za-z]", Cat)) %>% group_by(Cat) %>% summarise(countgllvm = n())

cathsobs <- merge(NOGcoreancomidreshs_cat,NOGXgllvmcoefhs_cat,by="Cat",all=T) %>% pivot_longer(cols = starts_with("count"), names_to = "method") %>% group_by(Cat) %>% mutate(sumvalue = sum(value,na.rm = T))

ggplot(cathsobs, aes(x=value,y=reorder(Cat, -sumvalue), fill=method)) + 
  geom_bar(stat = "identity") + 
  xlab("Number of eggNOG categories") + 
  ylab("eggNOG categories")+
  ggtitle("A. Differentially abundant eggNOG categories with genome-wide heterozygosity")+
  scale_fill_manual("Method",values = c("countancom"="olivedrab4","countgllvm"="orangered")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ 
  theme(axis.line = element_line(colour = "black", linetype=1))

# hs obs venn
NOGXgllvmcoefhsobs_data <- NOGXgllvmcoefhs %>% dplyr::select(direct,cog) %>% mutate(method="GLLVM") %>% filter(!is.na(cog))

NOGcoreancomidreshsobs_data <- NOGcoreancomidreshs %>% dplyr::select(direct,cog) %>% mutate(method="ANCOMBC2") %>% filter(!is.na(cog)) %>% data.frame()

noggllvmancom_hs <- list(ANCOMBC2=NOGcoreancomidreshsobs_data$cog,GLLVM=NOGXgllvmcoefhsobs_data$cog)

noghsvenn <- ggVennDiagram(noggllvmancom_hs,set_color = "black") + scale_fill_gradient2(low="white",mid="white",high = "white") +
  coord_flip() +
  ggtitle("C. ANCOMBC2 vs GLLVM MHC-I Diversity")+
  theme_tufte(base_size = 15, base_family = "Arial")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.1),legend.position="none")

tiff("Output/NOG_ANCOM_GLLVM_hsobs.tiff", res =400, units = "in", width = 15, height = 9, bg = "transparent")
ggarrange(hsobsnogancomplot,hsnoggllvmplot, noghsvenn, heights = c(1.5,1))
dev.off()

noggllvmancom2 <- merge(NOGXgllvmcoefmhc1_data,NOGcoreancomidresmhc1_data,by=c("direct","cog"),all=T) %>% mutate(tests = case_when(!is.na(method.x) & !is.na(method.y) ~ "both", !is.na(method.x) & is.na(method.y) ~ method.x, is.na(method.x) & !is.na(method.y) ~ method.y, TRUE ~ "none"))

noggllvmancom3 <- noggllvmancom2 %>% filter(tests %in% "both") %>% merge(.,Chuen_new_COG_tax,by.x="cog",by.y="COG",all.x=T )




# mpa ----
mpamat <- read.csv("~/Documents/PhD/R_analysis/Phyloseq/InputTables/MFF_Feb2024/merged_reads_abundance.txt", sep = "\t") %>% pivot_wider(names_from = TubeNo, values_from = Reads) %>% dplyr::mutate(OTU = paste0("OTU", 1:nrow(.))) %>% dplyr::rename(SW370=SW370_, SW29=SW29_,SW131=SW131_,SW55=SW55_)
# taxa table
mpataxa <- mpamat %>% dplyr::select(OTU,OUT) %>% separate_wider_delim(OUT, delim = "|", names = c("Kingdom","Phylum","Order","Class","Family","Genus","Species","SGB"), too_few = "align_start") %>% column_to_rownames(., var="OTU")
mpataxo7 <- mpataxa %>% filter(!is.na(SGB))
# otu table
mpaL7 <- mpamat  %>% dplyr::select(-OUT) %>% column_to_rownames(., var="OTU") %>% replace(is.na(.), 0)
mpaL7 <- mpaL7[rownames(mpaL7) %in% row.names(mpataxo7), ]

#matrix and phyloseq
mpaOTUL7 <- otu_table(as.matrix(mpaL7), taxa_are_rows = TRUE)

mpaTAXL7 = tax_table(as.matrix(mpataxo7))

mpapL7 <- phyloseq(mpaOTUL7, mpaTAXL7, stdata) 
#mpapL7 %>% filter(Type %in% "C") %>% distinct(Sample)
mpapL7sd <- data.frame(sample_data(mpapL7)) %>% rownames_to_column("TubeNumber") %>% mutate(negpos = case_when(TubeNumber %in% c("SW2021","SW2022","SW2017","SW2019","SW2020","SW1421") ~ TRUE, TRUE ~ FALSE)) %>% column_to_rownames("TubeNumber")

mpaL72 <- phyloseq(tax_table(mpapL7),otu_table(mpapL7),sample_data(mpapL7sd))
contamdf.prev <- isContaminant(mpaL72, method="prevalence", neg="negpos",threshold=0.4)
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))
ps.pa <- transform_sample_counts(mpaL72, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$negpos == TRUE, ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$negpos == FALSE, ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

mpapL7.bac <- subset_taxa(mpapL7, Kingdom %in% "k__Bacteria") %>% filter_sample_data(Type%in%"F") %>% filter_sample_data(SamplingAge > 0.5) %>% filter_sample_data(Type %in% "F") %>% mutate_sample_data(Age_scaled = arm::rescale(SamplingAge,"center"),Timeinfridge_scaled = arm::rescale(Timeinfridge, "center"),CatchTime_scaled = arm::rescale(CatchTime, "center"),TQcorrected_scaled = arm::rescale(TQcorrected, "center"), capage=pmin(SamplingAge,12)) #%>% filter_sample_data(!(TubeNo %in% c("SW1168","SW1183")))
mpapL7samplesums <- data.frame(sample_sums(mpapL7.bac))

##different levels
# mpa1p <- tax_glom(mpapL7.bac, "Phylum")
# mpa2p <- tax_glom(mpapL7.bac, "Class")
# mpa3p <- tax_glom(mpapL7.bac, "Order")
# mpa4p <- tax_glom(mpapL7.bac, "Family")
# mpa5p <- tax_glom(mpapL7.bac, "Genus")
# mpa6p <- tax_glom(mpapL7.bac, "Species")


mpapL7.bacpfsampledata <- data.frame(sample_data(mpapL7.bac)) %>% rownames_to_column("Rownames")

mpal7TLR <- merge(mpapL7.bacpfsampledata,TL3_MHC_20240530full, by="BirdID") %>% merge(.,SWhet, by="BirdID",all.x=T) %>% column_to_rownames("Rownames")

mpal7TLR %>% filter(!is.na(MHC1_Diversity),!is.na(Aseua1)) %>% nrow(.)
mpal7TLR %>% filter(!is.na(MHC1_Diversity),!is.na(Aseua1)) %>% {n_distinct(.$BirdID)}

mpal7TLR %>% group_by(Aseua7) %>% summarise(n=n())

mpaTLR <- phyloseq(tax_table(mpapL7.bac),otu_table(mpapL7.bac),sample_data(mpal7TLR)) %>% mutate_sample_data(SamplingAge2 = round(SamplingAge))

# difference between mpa and eggNOG - why is it different?
NOGsamplenames <- NOGpfTLR %>% rownames_to_column("nognames") %>% dplyr::select(nognames) %>% arrange(nognames)

mpasamplenames <- mpal7TLR %>% rownames_to_column("mpanames") %>% dplyr::select(mpanames) %>% arrange(mpanames)

dplyr::symdiff(NOGsamplenames$nognames,mpasamplenames$mpanames)



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

# mpa beta 
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
#mpaTLRCLRperm0<- adonis2(mpaTLRCLRmat ~ SamplingAge + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled  +TQcorrected_scaled   , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
#mpaTLRCLRperm0

#mpaTLRCLRperm<- adonis2(mpaTLRCLRmat ~  Aseua1+Aseua5+Aseua7+Aseua9 + SamplingAge + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled  +TQcorrected_scaled   , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
#mpaTLRCLRperm
#AICc_permanova2(mpaTLRCLRperm)

#mpaTLRCLRperm2<- adonis2(mpaTLRCLRmat ~ Aseua5 + Aseua7+Aseua9 + SamplingAge + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled +TQcorrected_scaled    , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
#mpaTLRCLRperm2

mpaTLRCLRperm2S<- adonis2(mpaTLRCLRmat ~ Aseua5 + Aseua7+Aseua9 + cutage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled, data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
mpaTLRCLRperm2S

# mpaTLRCLRperm2A<- adonis2(mpaTLRCLRmat ~ Aseua7+Aseua9  + season + SampleYear  + SexEstimate + CatchTime_scaled    , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# mpaTLRCLRperm2A


# mpaTLRCLRperm2<- adonis2(mpaTLRCLRmat ~ capage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + Asedab3+Asedab4+Asedab5 + Aseua1+Aseua3+Aseua4+Aseua5+Aseua6+Aseua7+ Aseua8+ Aseua11   , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# mpaTLRCLRperm2
# 
# mpaTLRCLRperm3<- adonis2(mpaTLRCLRmat ~ capage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + Aseua11   , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# mpaTLRCLRperm3

mpaTLRCLRperm4<- adonis2(mpaTLRCLRmat ~Hs_obs + MHC1_Diversity + MHC2_Diversity+ cutage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled   , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
mpaTLRCLRperm4

mpaTLRCLRperm5<- adonis2(mpaTLRCLRmat ~ Hs_obs +MHC1_Diversity + I(MHC1_Diversity^2)+ MHC2_Diversity+  SamplingAge + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled   , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
mpaTLRCLRperm5

mpaTLRCLRperm4S<- adonis2(mpaTLRCLRmat ~Hs_obs + MHC1_Diversity + MHC2_Diversity + season + SampleYear  + SexEstimate +Timeinfridge_scaled+ CatchTime_scaled, data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
mpaTLRCLRperm4S

# mpaTLRCLRperm4A<- adonis2(mpaTLRCLRmat ~ MHC1_Diversity + MHC2_Diversity + + season + SampleYear  + SexEstimate + CatchTime_scaled    , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# mpaTLRCLRperm4A
# 
# mpaTLRCLRperm5A<- adonis2(mpaTLRCLRmat ~ MHC1_Diversity + I(MHC1_Diversity^2)+ MHC2_Diversity + season + SampleYear  + SexEstimate + CatchTime_scaled    , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# mpaTLRCLRperm5A


mpaTLRCLRperm6<- adonis2(mpaTLRCLRmat ~ Hs_obs +MHC1_Diversity + I(MHC1_Diversity^2)+ MHC2_Diversity + I(MHC2_Diversity^2)+  SamplingAge + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled   , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
mpaTLRCLRperm6

# mpaTLRCLRperm5<- adonis2(mpaTLRCLRmat ~ capage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + MHC1_Diversity*SexEstimate + MHC2_Diversity*SexEstimate   , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# mpaTLRCLRperm5
# 
# mpaTLRCLRperm6<- adonis2(mpaTLRCLRmat ~ capage + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled + MHC1_Diversity + I(MHC1_Diversity^2) + MHC2_Diversity + I(MHC2_Diversity^2)    , data=mpaTLRCLR_st, permutations = perm, method = "euclidean", by= "margin")
# mpaTLRCLRperm6

# mpaTLRCLR_stmatrix <- as(sample_data(mpaTLRCLR),"data.frame")  %>% dplyr::select(Aseua11,Aseua8) %>% arrange(row.names(.)) %>% mutate_all(as.numeric) %>% vegdist(., method="jaccard")
# 
# mpaMHCmatrix <- as.matrix(mpaTLRCLR_stmatrix)
# 
# mpamicrobiomematrix <- data.frame(t(data.frame(otu_table(mpaTLRCLR)))) %>% arrange(row.names(.)) %>% vegdist(., method="euclidean")
# mpamicrobiomematrix2 <- as.matrix(mpamicrobiomematrix)
# 
# mpaTLRCLR_stnew <- data.frame(sample_data(mpaTLRCLR))%>% dplyr::select(SampleYear,SexEstimate)%>% mutate_all(as.character) %>% mutate_all(as.numeric) %>% arrange(row.names(.)) %>% vegdist(.)
# mpaTLRCLR_stnewmat <- as.matrix(mpaTLRCLR_stnew)
# 
# mantel(mpamicrobiomematrix2,mpaMHCmatrix, method="pearson")
# 
# cor(c(mpamicrobiomematrix2),c(mpaMHCmatrix))
# cxy <- cancor(mpamicrobiomematrix2,mpaTLRCLR_stnewmat)
# summary(cxy)
# cor.test(c(mpamicrobiomematrix2),c(mpaTLRCLR_stnewmat),method="pearson")
# cor.test(c(mpamicrobiomematrix2),c(mpaMHCmatrix),method="pearson")

# mpa DAA
mpaTLRcoremembers <- transform(mpaTLR, "compositional") %>% core_members(., detection = 0.00001, prevalence = 0.1)
mpaTLRcore <- prune_taxa(mpaTLRcoremembers, mpaTLR)
# mpa DAA ancombc2 
mpaTLRancomdiv <- ancombc2(mpaTLRcore,fix_formula = "Hs_obs + Aseua7 + season + SampleYear  + SexEstimate + CatchTime_scaled",rand_formula = "(1|BirdID)" , p_adj_method = "BH", pseudo_sens = FALSE)
mpaTLRancomdivres <- mpaTLRancomdiv$res

mpataxo7spec <- mpataxo7 %>% rownames_to_column("taxonsig")
mpaTLRancomdivreshsobs <- mpaTLRancomdivres %>% dplyr::select(taxon, contains("Hs_obs")) %>% mutate(direct = case_when(lfc_Hs_obs < 0 & p_Hs_obs < 0.05 ~ "Negative_LFC", lfc_Hs_obs > 0 & p_Hs_obs < 0.05 ~ "Positive_LFC", TRUE ~ "p > 0.05"), cog = case_when(p_Hs_obs < 0.05 ~ taxon)) %>% mutate(taxonsig = case_when(p_Hs_obs < 0.1 ~ taxon)) %>% merge(.,mpataxo7spec, by = "taxonsig", all.x=T)

hsobsmpaancomplot <- ggplot(mpaTLRancomdivreshsobs, aes(x= lfc_Hs_obs,y=  -log(p_Hs_obs),color=direct, label=Species)) +
  geom_point(size=3) +
  geom_text_repel(size=3, colour="black") +
  ggtitle("A. ANCOMBC2 Heterozygosity")+
  xlab("Log fold change Heterozygosity") +
  ylab("-log(p-value)") +
  scale_color_manual("Direction",values = c("Positive_LFC"="olivedrab","Negative_LFC"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

mpaTLRancomdivresaseua7 <- mpaTLRancomdivres %>% dplyr::select(taxon, contains("Aseua71")) %>% mutate(direct = case_when(lfc_Aseua71 < 0 & p_Aseua71 < 0.05 ~ "Negative_LFC", lfc_Aseua71 > 0 & p_Aseua71 < 0.05 ~ "Positive_LFC", TRUE ~ "p > 0.05"), cog = case_when(p_Aseua71 < 0.05 ~ taxon)) %>% mutate(taxonsig = case_when(p_Aseua71 < 0.1 ~ taxon)) %>% merge(.,mpataxo7spec, by = "taxonsig", all.x=T)

aseua7mpaancomplot <- ggplot(mpaTLRancomdivresaseua7, aes(x= lfc_Aseua71,y=  -log(p_Aseua71),color=direct,label =Species)) +
  geom_point(size=3) +
  geom_text_repel(size=3, colour="black") +
  ggtitle("C. ANCOMBC2 Aseua 7")+
  xlab("Log fold change Aseua 7") +
  ylab("-log(p-value)") +
  scale_color_manual("Direction",values = c("Positive_LFC"="olivedrab","Negative_LFC"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

#mpa DAA gllvm
mpaTLRcoreclr <- transform(mpaTLRcore,"clr")

mpaTLRcoreclrmat <- vegan_otu(mpaTLRcoreclr)
mpaTLRcoreclrst <- as(sample_data(mpaTLRcoreclr),"data.frame") %>% 
  mutate(TerminalYear = as.factor(TerminalYear)) %>% 
  mutate(capage=pmin(SamplingAge,12))%>%     
  dplyr::select(MHC1_Diversity, MHC2_Diversity,TerminalYear,BirdID,season,SampleYear,SamplingAge2, SexEstimate,CatchTime_scaled,Timeinfridge_scaled,TQcorrected_scaled,Aseua7, Aseua11,Hs_obs) 
mpaTLRcoreclrgllvm <- gllvm::gllvm(mpaTLRcoreclrmat, mpaTLRcoreclrst, formula = ~ Hs_obs + Aseua7 + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled , family = "gaussian",num.lv = 1,sd.errors = TRUE, row.eff = ~(1|BirdID), seed = 8,starting.val="zero")
gllvm::coefplot.gllvm(mpaTLRcoreclrgllvm,which.Xcoef = 1:2)
summary(mpaTLRcoreclrgllvm)

mpaTLRCoreclrgllvmsum <- summary(mpaTLRcoreclrgllvm)
mpaTLRCoreclrgllvmcoef <- as.data.frame(mpaTLRCoreclrgllvmsum[["Coef.tableX"]]) %>% mutate(terms= row.names(.)) %>% separate_wider_regex(cols = terms, c(terms = ".*", ":", OTU = ".*"))

mpaTLRCoreclrgllvmhsobs <- mpaTLRCoreclrgllvmcoef %>% filter(terms %in% "Hs_obs") %>% rename(pval = 'Pr(>|z|)', se = 'Std. Error') %>% mutate(direct = case_when(Estimate < 0 & pval < 0.05 ~ "Negative_LFC", Estimate > 0 & pval < 0.05 ~ "Positive_LFC", TRUE ~ "p > 0.05"),cog = case_when(pval < 0.05 ~ OTU)) %>% mutate(taxonsig = case_when(pval < 0.1 ~ OTU)) %>% merge(.,mpataxo7spec, by = "taxonsig", all.x=T)

hsobsmpagllvmplot <- ggplot(mpaTLRCoreclrgllvmhsobs, aes(x= Estimate,y=  -log(pval),color=direct, label=Species)) +
  geom_point(size =3) +
  geom_text_repel(size=3, colour="black") +
  ggtitle("B. GLLVM Heterozygosity")+
  xlab("Log fold change Heterozygosity") +
  ylab("-log(p-value)")+
  scale_color_manual("Direction",values = c("Positive_LFC"="olivedrab4","Negative_LFC"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

mpaTLRCoreclrgllvmaseua7 <- mpaTLRCoreclrgllvmcoef %>% filter(terms %in% "Aseua71") %>% rename(pval = 'Pr(>|z|)', se = 'Std. Error') %>% mutate(direct = case_when(Estimate < 0 & pval < 0.05 ~ "Negative_LFC", Estimate > 0 & pval < 0.05 ~ "Positive_LFC", TRUE ~ "p > 0.05"),cog = case_when(pval < 0.05 ~ OTU)) %>% mutate(taxonsig = case_when(pval < 0.1 ~ OTU)) %>% merge(.,mpataxo7spec, by = "taxonsig", all.x=T)

Aseua7mpagllvmplot <- ggplot(mpaTLRCoreclrgllvmaseua7, aes(x= Estimate,y=  -log(pval),color=direct, label=Species)) +
  geom_point(size =3) +
  geom_text_repel(size=3, colour="black") +
  ggtitle("D. GLLVM Aseua 7")+
  xlab("Log fold change Aseua 7") +
  ylab("-log(p-value)")+
  scale_color_manual("Direction",values = c("Positive_LFC"="olivedrab4","Negative_LFC"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))



tiff("Output/mpa_ANCOM_GLLVM.tiff", res =400, units = "in", width = 15, height = 9, bg = "transparent")
ggarrange(hsobsmpaancomplot,hsobsmpagllvmplot,aseua7mpaancomplot,Aseua7mpagllvmplot,common.legend = T,legend = "right")
dev.off()

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


# aldex2 mpa
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

# physeq mhc ----
physeq4 <- readRDS("InputTables/physeq4.rds")
physeq4A <- subset_samples(physeq4, Ageclass %in% c("A"))
physeq3sampledata <- data.frame(sample_data(physeq4A)) %>% rownames_to_column("Rownames")

physeq3TLR <- merge(physeq3sampledata,TL3_MHC_20240530full, by="BirdID") %>% dplyr::rename(season=Season,Timeinfridge=TimeAt4Degrees,CatchTime=MinutesSinceSunrise, SamplingAge = AgeYears) %>% mutate(Timeinfridge_scaled = arm::rescale(Timeinfridge),CatchTime_scaled = arm::rescale(CatchTime), TQcorrected_scaled = arm::rescale(TQcorrected) ) %>% mutate(TQcorrected_scaled = replace_na(TQcorrected_scaled, mean(TQcorrected_scaled,na.rm=T))) %>% mutate_at(vars(contains('Ase')), as.factor) %>% mutate(SampleYear = as.factor(SampleYear), season = as.factor(season), SexEstimate = as.factor(SexEstimate)) %>% merge(.,SWhet, by="BirdID",all.x=T) %>% mutate(survive = case_when(FieldPeriodID == FieldPeriodID_LastSeen ~ "died", TRUE ~ "survived")) %>% column_to_rownames("Rownames")

n_distinct(physeq3TLR$BirdID)
n_distinct(physeq3TLR$TubeNumber)
dplyr::distinct(physeq3TLR,Ageclass)
physeq3TLR %>% group_by(Aseua7) %>% summarise(n=n())

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


# 16s alpha diversity ----
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

#16s ancom
physeq3TLRcoremembers <- transform(physeq3TLRp,"compositional") %>% core_members(.,detection = 0.001, prevalence = 0.2)
physeq3TLRCore <- prune_taxa(physeq3TLRcoremembers,physeq3TLRp) 

physeq3coreancomid <- ancombc2(physeq3TLRCore,fix_formula = "Hs_obs + MHC1_Diversity  + MHC2_Diversity  + SamplingAge + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled",rand_formula = "(1|BirdID)" , p_adj_method = "BH")
physeq3coreancomidres <- physeq3coreancomid$res %>% merge(.,physeq4ataxa, by="taxon")

physeq3coreancomidresmhc1 <- physeq3coreancomidres %>% dplyr::select(taxon, contains("MHC1"),Species,Genus,Family,Phylum,Order) %>% mutate(direct = case_when(lfc_MHC1_Diversity < 0 & p_MHC1_Diversity < 0.05 ~ "Negative_LFC", lfc_MHC1_Diversity > 0 & p_MHC1_Diversity < 0.05 ~ "Positive_LFC", TRUE ~ "p > 0.05"), cog = case_when(p_MHC1_Diversity < 0.05 ~ taxon)) 

physeq3coreancomidresmhc2 <- physeq3coreancomidres %>% dplyr::select(taxon, contains("MHC2"),Species,Genus,Family,Phylum,Order) %>% mutate(direct = case_when(lfc_MHC2_Diversity < 0 & p_MHC2_Diversity < 0.05 ~ "Negative_LFC", lfc_MHC2_Diversity > 0 & p_MHC2_Diversity < 0.05 ~ "Positive_LFC", TRUE ~ "p > 0.05"), cog = case_when(p_MHC2_Diversity < 0.05 ~ taxon)) 

physeq3coreancomidreshsobs <- physeq3coreancomidres %>% dplyr::select(taxon, contains("Hs_obs"),Species,Genus,Family,Phylum,Order) %>% mutate(direct = case_when(lfc_Hs_obs < 0 & p_Hs_obs < 0.05 ~ "Negative_LFC", lfc_Hs_obs > 0 & p_Hs_obs < 0.05 ~ "Positive_LFC", TRUE ~ "p > 0.05"), cog = case_when(p_Hs_obs < 0.05 ~ taxon)) 

physeq3coreancomidresdiv <- physeq3coreancomidres %>% dplyr::select(taxon, contains("MHC"), contains("Hs_obs"),Species,Genus,Family,Phylum,Order) %>% filter(p_MHC1_Diversity < 0.05 | p_MHC2_Diversity < 0.05 | p_Hs_obs < 0.05 )


physeq3coreancomid2 <- ancombc2(physeq3TLRCore,fix_formula = "Aseua1 + Aseua5 + Aseua7 + Aseua9 + SamplingAge + season + SampleYear  + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled",rand_formula = "(1|BirdID)" , p_adj_method = "BH")
physeq3coreancomid2res <- physeq3coreancomid2$res %>% merge(.,physeq4ataxa, by="taxon")

physeq3coreancomid2resase <- physeq3coreancomid2res %>% dplyr::select(contains("Ase"), taxon,Species,Genus,Family,Phylum,Order) %>% filter(p_Aseua11 < 0.05 | p_Aseua51 < 0.05 | p_Aseua71 < 0.05 | p_Aseua91 < 0.05)

# 16s gllvm 
physeq3TLRCoreclr <- transform(physeq3TLRCore,"clr")
physeq3TLRCoreclrmat <- vegan_otu(physeq3TLRCoreclr)
physeq3TLRCoreclrst <- as(sample_data(physeq3TLRCoreclr),"data.frame") %>% 
  mutate(TerminalYear = as.factor(TerminalYear)) %>% 
  mutate(capage=pmin(SamplingAge,12))%>%     
  dplyr::select(MHC1_Diversity, MHC2_Diversity,TerminalYear,BirdID,season,SampleYear,SamplingAge, SexEstimate,CatchTime_scaled,Timeinfridge_scaled,TQcorrected_scaled,Hs_obs, contains("Ase")) 
physeq3TLRCoreclrgllvm <- gllvm(physeq3TLRCoreclrmat, physeq3TLRCoreclrst, formula = ~  Hs_obs+MHC1_Diversity  + MHC2_Diversity  + SamplingAge + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled , family = "gaussian",num.lv = 4,sd.errors = TRUE, row.eff = ~(1|BirdID), seed = 8,starting.val="random")
coefplot.gllvm(physeq3TLRCoreclrgllvm,which.Xcoef = 1:3)
coefplot.gllvm(physeq3TLRCoreclrgllvm,which.Xcoef = 4:5)
summary(physeq3TLRCoreclrgllvm)
physeq3TLRCoreclrgllvmsum <- summary(physeq3TLRCoreclrgllvm)
physeq3TLRCoreclrgllvmcoef <- as.data.frame(physeq3TLRCoreclrgllvmsum[["Coef.tableX"]]) %>% mutate(terms= row.names(.)) %>% separate_wider_regex(cols = terms, c(terms = ".*", ":", OTU = ".*")) %>% merge(.,physeq4ataxa,by.x="OTU", by.y="taxon",all.x=T)
gllvm::AICc(physeq3TLRCoreclrgllvm)

physeq3TLRCoreclrgllvmcoefdiv <- physeq3TLRCoreclrgllvmcoef %>% rename(pcoef = 'Pr(>|z|)') %>% filter(terms %in% c("MHC1_Diversity","MHC2_Diversity","Hs_obs") & pcoef < 0.05) 

#16s hs obs
physeq3coreanhsdiv <- physeq3coreancomidreshsobs %>% filter(p_Hs_obs < 0.05) %>% dplyr::select(lfc_Hs_obs,se_Hs_obs,Order,Genus) %>% mutate(method = "ancombc2")
physeq3TLRCoreclrgllvmhsdiv <- physeq3TLRCoreclrgllvmcoefdiv %>% filter(terms %in% "Hs_obs") %>% rename(lfc_Hs_obs=Estimate,se_Hs_obs = 'Std. Error')%>% dplyr::select(lfc_Hs_obs,se_Hs_obs,Order,Genus) %>% mutate(method = "gllvm")
physeq3corehs <- rbind(physeq3coreanhsdiv,physeq3TLRCoreclrgllvmhsdiv)

hs16s <- ggplot(physeq3corehs, aes(x=lfc_Hs_obs, y= Order, colour=method)) + geom_point() + geom_errorbarh(aes(xmin=lfc_Hs_obs-se_Hs_obs, xmax=lfc_Hs_obs+se_Hs_obs),height=0.2)+
  scale_color_manual("DAA method",values = c("ancombc2"="olivedrab4","gllvm"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

#16s mhc1 diversity
physeq3coreanmhc1div <- physeq3coreancomidresmhc1 %>% filter(p_MHC1_Diversity < 0.05) %>% dplyr::select(lfc_MHC1_Diversity,se_MHC1_Diversity,Order,Genus) %>% mutate(method = "ancombc2")
physeq3TLRCoreclrgllvmmhc1div <- physeq3TLRCoreclrgllvmcoefdiv %>% filter(terms %in% "MHC1_Diversity") %>% rename(lfc_MHC1_Diversity=Estimate,se_MHC1_Diversity = 'Std. Error')%>% dplyr::select(lfc_MHC1_Diversity,se_MHC1_Diversity,Order,Genus) %>% mutate(method = "gllvm")
physeq3coremhc1 <- rbind(physeq3coreanmhc1div,physeq3TLRCoreclrgllvmmhc1div)

mhc16s <- ggplot(physeq3coremhc1, aes(x=lfc_MHC1_Diversity, y= Order, colour=method)) + geom_point() + geom_errorbarh(aes(xmin=lfc_MHC1_Diversity-se_MHC1_Diversity, xmax=lfc_MHC1_Diversity+se_MHC1_Diversity),height=0.2)+
  scale_color_manual("DAA method",values = c("ancombc2"="olivedrab4","gllvm"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

#16s mhc2 diversity
physeq3coreanmhc2div <- physeq3coreancomidresmhc2 %>% filter(p_MHC2_Diversity < 0.05) %>% dplyr::select(lfc_MHC2_Diversity,se_MHC2_Diversity,Order,Genus) %>% mutate(method = "ancombc2")
physeq3TLRCoreclrgllvmmhc2div <- physeq3TLRCoreclrgllvmcoefdiv %>% filter(terms %in% "MHC2_Diversity") %>% rename(lfc_MHC2_Diversity=Estimate,se_MHC2_Diversity = 'Std. Error')%>% dplyr::select(lfc_MHC2_Diversity,se_MHC2_Diversity,Order,Genus) %>% mutate(method = "gllvm")
physeq3coremhc2 <- rbind(physeq3coreanmhc2div,physeq3TLRCoreclrgllvmmhc2div)

mhc216s<-ggplot(physeq3coremhc2, aes(x=lfc_MHC2_Diversity, y= Order, colour=method)) + geom_point() + geom_errorbarh(aes(xmin=lfc_MHC2_Diversity-se_MHC2_Diversity, xmax=lfc_MHC2_Diversity+se_MHC2_Diversity),height=0.2)+
  scale_color_manual("DAA method",values = c("ancombc2"="olivedrab4","gllvm"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

# 16s daa gllvm 2
physeq3TLRCoreclrgllvm2 <- gllvm(physeq3TLRCoreclrmat, physeq3TLRCoreclrst, formula = ~  Aseua1 + Aseua5 + Aseua7 + Aseua9  + SamplingAge + season + SampleYear + SexEstimate+ Timeinfridge_scaled + CatchTime_scaled + TQcorrected_scaled , family = "gaussian",num.lv = 2,sd.errors = TRUE, row.eff = ~(1|BirdID), seed = 8,starting.val="random")
coefplot.gllvm(physeq3TLRCoreclrgllvm2,which.Xcoef = 1:3)

physeq3TLRCoreclrgllvm2sum <- summary(physeq3TLRCoreclrgllvm2)
physeq3TLRCoreclrgllvm2sum2 <-  as.data.frame(physeq3TLRCoreclrgllvm2sum[["Coef.tableX"]]) %>% mutate(terms= row.names(.)) %>% separate_wider_regex(cols = terms, c(terms = ".*", ":", OTU = ".*")) %>% merge(.,physeq4ataxa,by.x="OTU", by.y="taxon",all.x=T)  %>% rename(pcoef = 'Pr(>|z|)') 

# 16s aseua1
physeq3anase1 <- physeq3coreancomid2res %>% dplyr::select(contains("Aseua1"),Order,Genus) %>% filter(p_Aseua11 < 0.05) %>% dplyr::select(lfc_Aseua11,se_Aseua11,Order,Genus) %>% mutate(method = "ancombc2")
physeq3TLRCoreclrgllvase1 <- physeq3TLRCoreclrgllvm2sum2 %>% filter(terms %in% "Aseua11" & pcoef < 0.05) %>% rename(lfc_Aseua11=Estimate,se_Aseua11 = 'Std. Error')%>% dplyr::select(lfc_Aseua11,se_Aseua11,Order,Genus) %>% mutate(method = "gllvm")
physeq3coreaseua1 <- rbind(physeq3anase1,physeq3TLRCoreclrgllvase1)

aseua116s <- ggplot(physeq3coreaseua1, aes(x=lfc_Aseua11, y= Order, colour=method)) + geom_point() + geom_errorbarh(aes(xmin=lfc_Aseua11-se_Aseua11, xmax=lfc_Aseua11+se_Aseua11),height=0.2)+
  scale_color_manual("DAA method",values = c("ancombc2"="olivedrab4","gllvm"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

# 16s aseua5
physeq3anase5 <- physeq3coreancomid2res %>% dplyr::select(contains("Aseua5"),Order,Genus) %>% filter(p_Aseua51 < 0.05) %>% dplyr::select(lfc_Aseua51,se_Aseua51,Order,Genus) %>% mutate(method = "ancombc2")
physeq3TLRCoreclrgllvase5 <- physeq3TLRCoreclrgllvm2sum2 %>% filter(terms %in% "Aseua51" & pcoef < 0.05) %>% rename(lfc_Aseua51=Estimate,se_Aseua51 = 'Std. Error')%>% dplyr::select(lfc_Aseua51,se_Aseua51,Order,Genus) %>% mutate(method = "gllvm")
physeq3coreaseua5 <- rbind(physeq3anase5,physeq3TLRCoreclrgllvase5)

aseua516s <-ggplot(physeq3coreaseua5, aes(x=lfc_Aseua51, y= Order, colour=method)) + geom_point() + geom_errorbarh(aes(xmin=lfc_Aseua51-se_Aseua51, xmax=lfc_Aseua51+se_Aseua51),height=0.2)+
  scale_color_manual("DAA method",values = c("ancombc2"="olivedrab4","gllvm"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

# 16s aseua7
physeq3anase7 <- physeq3coreancomid2res %>% dplyr::select(contains("Aseua7"),Order,Genus) %>% filter(p_Aseua71 < 0.05) %>% dplyr::select(lfc_Aseua71,se_Aseua71,Order,Genus) %>% mutate(method = "ancombc2")
physeq3TLRCoreclrgllvase7 <- physeq3TLRCoreclrgllvm2sum2 %>% filter(terms %in% "Aseua71" & pcoef < 0.05) %>% rename(lfc_Aseua71=Estimate,se_Aseua71 = 'Std. Error')%>% dplyr::select(lfc_Aseua71,se_Aseua71,Order,Genus) %>% mutate(method = "gllvm")
physeq3coreaseua7 <- rbind(physeq3anase7,physeq3TLRCoreclrgllvase7)

aseua716s <-ggplot(physeq3coreaseua7, aes(x=lfc_Aseua71, y= Order, colour=method)) + geom_point() + geom_errorbarh(aes(xmin=lfc_Aseua71-se_Aseua71, xmax=lfc_Aseua71+se_Aseua71),height=0.2)+
  scale_color_manual("DAA method",values = c("ancombc2"="olivedrab4","gllvm"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

# 16s aseua9
physeq3anase9 <- physeq3coreancomid2res %>% dplyr::select(contains("Aseua9"),Order,Genus) %>% filter(p_Aseua91 < 0.05) %>% dplyr::select(lfc_Aseua91,se_Aseua91,Order,Genus) %>% mutate(method = "ancombc2")
physeq3TLRCoreclrgllvase9 <- physeq3TLRCoreclrgllvm2sum2 %>% filter(terms %in% "Aseua91" & pcoef < 0.05) %>% rename(lfc_Aseua91=Estimate,se_Aseua91 = 'Std. Error')%>% dplyr::select(lfc_Aseua91,se_Aseua91,Order,Genus) %>% mutate(method = "gllvm")
physeq3coreaseua9 <- rbind(physeq3anase9,physeq3TLRCoreclrgllvase9)

aseua916s <-ggplot(physeq3coreaseua9, aes(x=lfc_Aseua91, y= Order, colour=method)) + geom_point() + geom_errorbarh(aes(xmin=lfc_Aseua91-se_Aseua91, xmax=lfc_Aseua91+se_Aseua91),height=0.2)+
  scale_color_manual("DAA method",values = c("ancombc2"="olivedrab4","gllvm"="orangered3","p > 0.05"="seashell3")) +
  theme_tufte(base_size = 15, base_family = "Arial")+ theme(axis.line = element_line(colour = "black", linetype=1))

ggarrange(hs16s,mhc16s,mhc216s,aseua116s,aseua516s,aseua716s,aseua916s, common.legend = T)


# GO terms
EMGO <- read_tsv("~/Documents/PhD/R_analysis/Phyloseq/InputTables/MFF_Feb2024/emapper/EM.GOL0.txt")
