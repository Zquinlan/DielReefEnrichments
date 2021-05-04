--- 
layout: posts 
classes: wide 
sidebar:
  nav: "content" 
---
```
## Script Written by Zach Quinlan 06/19/19
# Re-organization of DORCIERR_FCM_fDOM.R because it needs to be cleaner 07/15/2019
# Only working on daytime exudation and remineralization
# Rewritten with changes to RR3 starting pipeline 10/11/2019
# 16s rRNA amplicon seque3nce data added in a upaded 10/18/2019
# Added in ClassyFire annotations from Inchi_keys 10/07/2019
# This has been rewritten for the new pipeline 12/18/2019
# This was split to just look at metabolite trends on 05/12/2020


# LOADING -- packages -------------------------------------------------------
#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)
library(readxl)
library(multcomp)
library(CHNOSZ)
library(furrr)
library(future)
library(biclustermd)
library(webchem)
library(classyfireR)
library(randomForest)
library(ggpubr)
library(rsq)

#PCoA, PERMANOVA
library(vegan)
library(ape)

#Visualizations
library(wesanderson)
library(RColorBrewer)
library(gplots)
library(gtable)
library(ggnewscale)

#Defining functions and removing issues with overlapping function calls
map <- purrr::map
select <- dplyr::select
tidy <- broom::tidy
rename <- dplyr::rename
mutate <- dplyr::mutate

zscore <- function(x) {
  (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}


# LOADING -- Dataframes ---------------------------------------------------
## FCM
dorc_fcm_fdom <- read_xlsx("~/Documents/GitHub/DORCIERR/data/raw/DOC_fDOM_FCM/DORCIERR_fDOM_FCM.xlsx")%>%
  rename(sample_name =`Sample Name of DORCIERR_FCM_Final`)%>%
  rename('DayNight' = 'Experiment')%>%
  rename(Organism = 'Organsim')%>%
  mutate(Organism = case_when(Organism == "Water" ~ "Water control",
                              TRUE ~ as.character(Organism)))

# 16s rRNA sequences
microbe_abundance_raw <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/microbes/MCR2017.16S.Nelson.Pipeline.October2019/abundance_table_100.shared.tsv")
microbe_taxonomy <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/microbes/MCR2017.16S.Nelson.Pipeline.October2019/annotations_100.taxonomy.tsv")

# PRE-STATS CLEANING -- Microbes and pre-filtering OTUs for abundance-----------------------------------------------
microbe_combined <- microbe_abundance_raw%>%
  dplyr::select(-1)%>%
  mutate(Group = case_when(Group == "Dorcierr_D_DT_1_TFD" ~ "D_DT_1_TFD",
                           Group == "DORCIERR_D_WA_2_TFN" ~ "D_WA_2_TFN",
                           Group == "D_PV_2_TFN_SA504_SC704" ~ "D_PV_2_TFN",
                           Group == "D_PV_2_TFN_SA503_SC704" ~ "D_PV_3_TFN",
                           Group == "D_WA_4_TFN_SA503_SC703" ~ "D_PL_3_TFN",
                           Group == "D_WA_4_TFN_SA504_SC703" ~ "D_WA_4_TFN",
                           TRUE ~ as.character(Group)))%>%
  filter(Group %like% "%D_%")%>%
  rename(sample_code = Group)%>%
  gather(OTU, reads, 3:ncol(.))%>%
  # left_join(., microbe_taxonomy, by = "OTU")%>%
  # dplyr::select(-Size)%>%
  # separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";", remove = FALSE)%>%
  separate(sample_code, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_", remove = FALSE)%>%
  dplyr::mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                       Experiment == "M" ~ "mordor",
                                       Experiment == "R" ~ "RR3",
                                       TRUE ~ as.character(Experiment)))%>%
  dplyr::mutate(Organism = case_when(Organism == "CC" ~ "CCA",
                                     Organism == "DT" ~ "Dictyota",
                                     Organism == "PL" ~ "Porites lobata",
                                     Organism == "PV" ~ "Pocillopora verrucosa",
                                     Organism == "TR" ~ "Turf",
                                     Organism == "WA" ~ "Water control",
                                     Organism == "IN" ~ "Influent",
                                     Organism == "OF" ~ "Offshore",
                                     TRUE ~ as.character(Organism)))%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))%>%
  filter(Organism != "Offshore",
         Organism != "Influent")%>%
  mutate(reads = case_when(reads == 0 ~ 1/numOtus,
                           TRUE ~ as.numeric(reads)))%>%
  # left_join(fcm_wdf%>%
  #             rename(sample_code = sample_name)%>%
  #             select(c(`Cells µL-1`, sample_code)), by = "sample_code")%>%
  group_by(sample_code)%>%
  mutate(ra = reads/sum(reads))%>%
  # cell_abun = ra*`Cells µL-1`,
  # log10 = log10(cell_abun + 0.001))%>%
  ungroup()%>%
  group_by(OTU)%>%
  mutate(abundant = case_when(max(ra) > 0.01 | sum(ra > 0.001) >= 3 ~ "abundant",
                              TRUE ~ "rare"))

microbe_no_rare <- microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU, DayNight, Organism)%>%
  filter(abundant == "abundant")%>%
  mutate(asin = asin(sqrt(ra)))%>%
  select(-c(ra, abundant))


# PRE-STATS CLEANING -- microbe RA data  --------------------------------------------
ra_bigger_TF <- microbe_combined%>%
  filter(abundant == "abundant",
         Timepoint == 'TF')%>%
  select(-c(abundant, reads, numOtus, sample_code, Replicate))%>%
  group_by(Organism, DayNight, OTU)%>%
  summarize_if(is.numeric, mean)%>%
  spread(Organism, ra)%>%
  gather(Organism, ra, 3:7)%>%
  mutate(difference = ra - `Water control`)%>%
  filter(difference > 0)%>%
  select(DayNight, OTU, Organism)

average_ra <- microbe_no_rare%>%
  select(-c(Experiment, Replicate))%>%
  group_by(OTU, Organism, DayNight, Timepoint)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()

# PRE-STATS CLEAINING -- FCM Stats prep---------------------------------------------------------------------
## This should calculate mean cells per hour per µL for each organism
## Just looks at log10(TF) - mean(log10(T0))
## Will result in log10(cells*µL-1)*hr^-1
fcm_wdf <- dorc_fcm_fdom%>%
  dplyr::select(c(1:8,34:ncol(.)))

fcm_th_t0_prep <-fcm_wdf%>%
  filter(Timepoint %like any% c('T0', 'T4'),
         !Organism == 'Influent',
         !Organism == 'Offshore')%>%
  mutate(log_10_cellsµL = log10(`Cells µL-1`))%>%
  dplyr::select(c(5:8, 15))%>%
  spread(Timepoint, log_10_cellsµL)

fcm_t0 <- fcm_th_t0_prep%>%
  dplyr::select(-c(T4, Replicate))%>%
  group_by(Organism, DayNight)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)

fcm_rate_th_t0 <- fcm_th_t0_prep%>%
  dplyr::select(-T0)%>%
  left_join(., fcm_t0, by = c("Organism", "DayNight"))%>%
  add_column(change_per_hour = (.$T4 - .$T0)/24)

fcm_t7 <- fcm_wdf%>%
  mutate(log_10_cellsµL = log10(`Cells µL-1`))%>%
  dplyr::select(c(5:8, 15))%>%
  filter(Timepoint == c('TF'),
         !Organism == 'Influent')%>%
  spread(Timepoint, log_10_cellsµL)

## This is the actual dataframe to use for running FCM stats
fcm_stats_df <- left_join(fcm_t7, fcm_rate_th_t0, by = c("Organism", "Replicate", "DayNight"))%>%
  dplyr::select(-c(5,6))%>%
  gather(test, log_value, 4:5)


# SET SEED ----------------------------------------------------------------
set.seed(2005)

# STATS ANOVA -- Microbe TWO-Way ------------------------------------------
aov_microbe <- microbe_no_rare%>%
  group_by(OTU)%>%
  nest()%>%
  mutate(anova = map(data, ~ aov(asin ~ Organism*DayNight, .x)%>%
                       tidy()%>%
                       filter(!term == "Residuals")%>%
                       dplyr::select(term, p.value)))%>%
  dplyr::select(-data)%>%
  unnest(anova)

anova_microbe_pvalues <- aov_microbe%>%
  ungroup()%>%  
  filter(p.value < 0.05)%>%
  add_column(FDR = p.adjust(.$p.value, method = "BH"))%>%
  filter(FDR < 0.05)

organism_significant_microbes <- as.vector(anova_microbe_pvalues%>%
                                             filter(term == "Organism"))$OTU

DayNight_significant_microbes <- as.vector(anova_microbe_pvalues%>%
                                             filter(term != "Organism"))$OTU
# PRE-POST-HOC CLEANING -- Microbe Dunnetts and DayNight anova -------------------------------
mic_organism_post_hoc <- microbe_no_rare%>%
  filter(OTU %in% organism_significant_microbes)

daynight_microbe_post_hoc <- microbe_no_rare%>%
  filter(OTU %in% DayNight_significant_microbes)

# STATS POST-HOC -- FCM Tukeys ----------------------------------------------------------
# Tukey growth rates for the first half of the expierment
tukey_model_fcm <-fcm_stats_df%>%
  group_by(test, DayNight)%>%
  nest()%>%
  mutate(data = map(data, ~aov(log_value ~ Organism, .x)%>%
                      TukeyHSD(p.adjust.methods = "BH")%>%
                      tidy()))%>%
  unnest(data)%>%
  filter(adj.p.value < 0.05)

# STATS POST-HOC -- MICROBES Dunnetts -----------------------------
organism_order_micro <- as.factor(mic_organism_post_hoc$Organism)%>%
  relevel("Water control")%>%
  levels()%>%
  as.vector()

dunnett_microbe_pvals <- mic_organism_post_hoc%>%
  group_by(DayNight, OTU)%>%
  mutate(sum = sum(asin))%>%
  filter(sum != 0)%>%
  dplyr::select(-sum)%>%
  mutate(Organism = factor(Organism))%>%
  mutate(Organism = fct_relevel(Organism, organism_order_micro))%>%
  nest()%>%
  mutate(dunnett = map(data, ~ aov(asin ~ Organism, .x)%>%
                         glht(linfct = mcp(Organism = "Dunnett"))),
         dunnett_summary = map(dunnett, ~summary(.x)%>%
                                 tidy()))%>%
  dplyr::select(-c(data,dunnett))%>%
  unnest(dunnett_summary)%>%
  dplyr::select(-c(4:7))%>%
  mutate(lhs = gsub(" - Water control", "", lhs))%>%
  rename("Organism" = "lhs")%>%
  ungroup()%>%   
  add_column(FDR = p.adjust(.$p.value, method = "BH"))%>%
  filter(FDR < 0.05)

# STATS POST-HOC -- DayNight MICROBES t-test organism --------------------
daynight_microbe_pvals <- mic_organism_post_hoc%>%
  group_by(Organism, OTU)%>%
  mutate(sum = sum(asin))%>%
  filter(sum != 0)%>%
  dplyr::select(-sum)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(asin ~ DayNight, .x)%>%
                      tidy()))%>%
  unnest(data)%>%
  dplyr::select(-c(4:7))%>%
  filter(term != "Residuals")%>%
  ungroup()%>%   
  add_column(FDR = p.adjust(.$p.value, method = "BH"))%>%
  filter(FDR < 0.05)


# META-STATS -- microbes --------------------------------------------------
dunnett_micro_analysis <- dunnett_microbe_pvals%>%
  dplyr::select(-p.value)%>%
  spread(Organism, FDR)%>%
  add_column(number_exudate_organisms = rowSums(.[3:ncol(.)] >= 0, na.rm = TRUE))%>%
  mutate(microbe_organism = case_when(is.na(CCA) == FALSE & 
                                        is.na(Dictyota) &
                                        is.na(Turf) &
                                        is.na(`Pocillopora verrucosa`)  &
                                        is.na(`Porites lobata`) ~ "CCA",
                                      is.na(CCA)  & 
                                        is.na(Dictyota) == FALSE &
                                        is.na(Turf) &
                                        is.na(`Pocillopora verrucosa`)  &
                                        is.na(`Porites lobata`) ~ "Dictyota",
                                      is.na(CCA)  & 
                                        is.na(Dictyota) &
                                        is.na(Turf) == FALSE &
                                        is.na(`Pocillopora verrucosa`)  &
                                        is.na(`Porites lobata`) ~ "Turf",
                                      is.na(CCA) & 
                                        is.na(Dictyota) &
                                        is.na(Turf) &
                                        is.na(`Pocillopora verrucosa`) == FALSE &
                                        is.na(`Porites lobata`) ~ "Pocillopora verrucosa",
                                      is.na(CCA) & 
                                        is.na(Dictyota) &
                                        is.na(Turf) &
                                        is.na(`Pocillopora verrucosa`)  &
                                        is.na(`Porites lobata`) == FALSE ~ "Porites lobata",
                                      is.na(`Pocillopora verrucosa`) == FALSE &
                                        is.na(CCA) == FALSE & 
                                        is.na(Dictyota) &
                                        is.na(Turf)  |
                                        is.na(`Porites lobata`) == FALSE &
                                        is.na(CCA) == FALSE &
                                        is.na(Dictyota) &
                                        is.na(Turf) ~ "Corraline",
                                      is.na(`Pocillopora verrucosa`) &
                                        is.na(CCA) == FALSE & 
                                        is.na(Dictyota) == FALSE &
                                        is.na(`Porites lobata`) |
                                        is.na(CCA) == FALSE &
                                        is.na(Turf) == FALSE &
                                        is.na(`Porites lobata`) &
                                        is.na(`Pocillopora verrucosa`) ~ "Algae",
                                      is.na(`Dictyota`) &
                                        is.na(CCA)  &
                                        is.na(Turf) &
                                        is.na(`Porites lobata`) == FALSE &
                                        is.na(`Pocillopora verrucosa`) == FALSE ~ "Coral",
                                      is.na(`Dictyota`) == FALSE &
                                        is.na(CCA)  &
                                        is.na(Turf) == FALSE &
                                        is.na(`Porites lobata`) &
                                        is.na(`Pocillopora verrucosa`) ~ "Fleshy Algae",
                                      is.na(CCA)  & 
                                        is.na(Dictyota) &
                                        is.na(Turf) == FALSE &
                                        is.na(`Pocillopora verrucosa`)  &
                                        is.na(`Porites lobata`) ~ "Turf",
                                      number_exudate_organisms > 3 ~ "Primary Producers",
                                      TRUE ~ "Cosmo"))%>%
  left_join(microbe_taxonomy, by = "OTU")

micro_sigs_vector <- as.vector(dunnett_micro_analysis$OTU)


# META-STATS -- Hierarchical cluster matrix----------------------------------------
hc_microbe <- mic_organism_post_hoc%>%
  # filter(OTU %in% rf_microbe_sigs)%>%
  ungroup()%>%
  group_by(OTU)%>%
  mutate(zscore = zscore(log10))%>%
  ungroup()%>%
  dplyr::select(c(Organism, DayNight, Replicate, OTU, zscore))%>%
  unite(sample, c("Organism", "DayNight", "Replicate"), sep = "_")%>%
  left_join(microbe_taxonomy, by = "OTU")%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  unite(OFGO, c("Order", "Family", "Genus", "OTU"), sep = ";")%>%
  select(-c(Size, OTU_id, Kingdom, Phylum, Class))%>%
  spread(OFGO, zscore)

write_csv(hc_microbe%>%
            select(everything(), sample), "~/Documents/GitHub/DORCIERR/data/plots/microbe_hc_df.csv")

# META-STATS -- RA significant dunnetts -----------------------------------
ra_dunnetts <- dunnett_microbe_pvals%>%
  left_join(average_ra, by = c('OTU', 'Organism', "DayNight"))%>%
  left_join(microbe_taxonomy, by = 'OTU')%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "sp"), sep = ";")%>%
  unite(Tax_plot, c("Order", "Family", "Genus", "OTU"), sep = " ", remove = FALSE)

top5_microbes <- ra_dunnetts%>%
  group_by(DayNight, Organism)%>%
  nest()%>%
  mutate(data = map(data, ~ top_n(.x, 5, ra)))%>%
  filter(DayNight == 'Day')%>%
  unnest(data)

top5_microbes%>%
  ggplot(aes(Organism, ra, fill = Tax_plot)) +
  geom_bar(stat = 'identity')

pdf("./plots/org_microbe_ra.pdf", width = 15, height = 10)
average_ra%>%
  left_join(microbe_taxonomy, by = 'OTU')%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "sp"), sep = ";")%>%
  unite(Tax_plot, c("Order", "Family", "Genus", "OTU"), sep = " ", remove = FALSE)%>%
  filter(DayNight == 'Day')%>%
  ggplot(aes(x="", y=ra, fill = Class)) +
  geom_bar(width = 1, size = 0.5, color = "white", stat = "identity") +
  labs(x =NULL, y = NULL, title = "") +
  theme_classic() +
  facet_wrap(~Organism) +
  # scale_fill_manual(values = c('#78B7C5', '#EBCC2A', "#00A08A")) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  coord_polar("y", start = 0)
dev.off()

```