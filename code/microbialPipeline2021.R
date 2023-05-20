## Script Written by Zach Quinlan 06/19/19
# Re-organization of DORCIERR_FCM_fDOM.R because it needs to be cleaner 07/15/2019
# Only working on daytime exudation and remineralization
# Rewritten with changes to RR3 starting pipeline 10/11/2019
# 16s rRNA amplicon seque3nce data added in a upaded 10/18/2019
# Added in ClassyFire annotations from Inchi_keys 10/07/2019
# This has been rewritten for the new pipeline 12/18/2019
# This was split to just look at metabolite trends on 05/12/2020
# Edited and made for final publication merge 08/23/2021


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
library(ggpubr)
library(rsq)

library(emojifont)
load.emojifont('OpenSansEmoji.ttf')

#PCoA, PERMANOVA, permutational anova
library(vegan)
library(ape)
# library(permuco)
library(lmPerm)


#Visualizations
library(wesanderson)
library(RColorBrewer)
library(gplots)
library(gtable)
library(ggnewscale)
library(ggrepel)

#Defining functions and removing issues with overlapping function calls
map <- purrr::map
select <- dplyr::select
tidy <- broom::tidy
rename <- dplyr::rename
mutate <- dplyr::mutate


zscore <- function(x) {
  (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

gen_theme <-  function(x){
  theme(plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 25),
        axis.text.y = element_text(size = 25),
        plot.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
        panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"))
}

nameSamples <- function(x) { 
  separate(x, sample, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_", remove = FALSE)%>%
    mutate(Organism = case_when(Organism == "CC" ~ "CCA",
                                Organism == "DT" ~ "Dictyota",
                                Organism == "PL" ~ "Porites lobata",
                                Organism == "PV" ~ "Pocillopora verrucosa",
                                Organism == "TR" ~ "Turf",
                                Organism == "WA" ~ "Water control",
                                Organism == "IN" ~ "Influent",
                                Organism == "OF" ~ "Offshore",
                                TRUE ~ as.character(Organism)))%>%
    filter(!Organism %in% c('Influent', 'Offshore'))%>%
    separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
    mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                                TRUE ~ "Night"))
}


# LOADING -- Dataframes ---------------------------------------------------
# 16s rRNA sequences
# microbe_abundance_raw <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/microbes/MCR2017.16S.Nelson.Pipeline.October2019/abundance_table_100.shared.tsv")
# microbe_taxonomy <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/microbes/MCR2017.16S.Nelson.Pipeline.October2019/annotations_100.taxonomy.tsv")

unifracRaw <- read_tsv('~/Documents/GitHub/DORCIERR/data/raw/microbes/otu_repr_100.tre1.tsv')

microbe2022 <- read_tsv('~/Documents/GitHub/DORCIERR/data/raw/microbes/052022_abundance_table_100.tsv')%>%
  select(c(OTUConTaxonomy, OTUNumber, everything()))%>%
  gather(sample, ra, 3:ncol(.))%>%
  filter(sample %like any% c('D_%', 'DORC%'))

# PRE-STATS CLEANING -- Microbes and pre-filtering OTUs for abundance-----------------------------------------------
microbe_combined <- microbe2022%>%
  # select(-1)%>%
  mutate(sample = case_when(sample == "Dorcierr_D_DT_1_TFD" ~ "D_DT_1_TFD",
                            sample == "DORCIERR_D_WA_2_TFN" ~ "D_WA_2_TFN",
                           # Group == "D_PV_2_TFN_SA504_SC704" ~ "D_PV_2_TFN",
                           sample == "D_PV_2_TFN_SA503_SC704" ~ "D_PL_3_TFN",
                           sample == "D_WA_4_TFN_SA503_SC703" ~ "D_PV_3_TFN",
                           sample == "D_WA_4_TFN_SA504_SC703" ~ "D_WA_4_TFN",
                           TRUE ~ as.character(sample)))%>%
  filter(sample %like% "%D_%",
         sample != 'D_PV_2_TFN_SA504_SC704',
         sample != 'D_PV_4_TFD')%>%
  rename(sample_code = sample)%>%
  separate(sample_code, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_", remove = FALSE)%>%
  filter(Experiment == 'D')%>%
  mutate(Organism = case_when(Organism == "CC" ~ "CCA",
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
  # mutate(reads = case_when(reads == 0 ~ 1/numOtus,
  #                          TRUE ~ as.numeric(reads)))%>%
  # group_by(sample_code)%>%
  # mutate(ra = reads/sum(reads))%>%
  # ungroup()%>%
  group_by(OTUNumber)%>%
  mutate(ra = as.numeric(ra),
         abundant = case_when(max(ra) > 0.05 | sum(ra > 0.001) >= 3 ~ "abundant",
                              TRUE ~ "rare"))%>%
  unite(OTU, c('OTUConTaxonomy', 'OTUNumber'), sep = '')

# PRE-STATS -- qqplots ----------------------------------------------------
pdf('plots/qqplots.pdf', width = 15, height = 10)
microbe_combined%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU, DayNight, Organism)%>%
  filter(abundant == "abundant")%>% 
  ggplot(aes(ra)) + 
  geom_histogram(bins = 100) +
  scale_y_log10() + 
  labs(title = 'Distribution of relative abundances', y = 'Count', x = 'Relative Abundance') +
  gen_theme()

car::qqPlot((microbe_combined%>%
               filter(Timepoint == "TF")%>%
               group_by(OTU, DayNight, Organism)%>%
               filter(abundant == "abundant"))$ra, 
            ylab = "Relative Abundance", xlab = "Normal quantiles",
            main = 'QQ-plot: Relative Abundance') 

car::qqPlot((microbe_combined%>%
               filter(Timepoint == "TF")%>%
               group_by(OTU, DayNight, Organism)%>%
               filter(abundant == "abundant")%>%
               mutate(asin = asin(sqrt(ra))))$asin, 
            ylab = "Relative Abundance", xlab = "Normal quantiles",
            main = 'QQ-plot: Angular Transformation')

# car::qqPlot((microbe_combined%>%
#                filter(Timepoint == "TF")%>%
#                group_by(OTUNumber, DayNight, Organism)%>%
#                filter(abundant == "abundant")%>%
#                mutate(log = log10(ra)))$log, 
#             ylab = "Relative Abundance", xlab = "Normal quantiles",
#             main = 'QQ-plot: log10')
# 
# car::qqPlot((microbe_combined%>%
#                filter(Timepoint == "TF")%>%
#                group_by(OTUNumber, DayNight, Organism)%>%
#                filter(abundant == "abundant")%>%
#                mutate(log = Logit(ra)))$log, 
#             ylab = "Relative Abundance", xlab = "Normal quantiles",
#             main = 'QQ-plot: Logit')


car::qqPlot((microbe_combined%>%
               filter(Timepoint == "TF")%>%
               group_by(OTU, DayNight, Organism)%>%
               filter(abundant == "abundant")%>%
               mutate(clr = compositions::clr(ra)))$clr, 
            ylab = "Relative Abundance", xlab = "Normal quantiles",
            main = 'QQ-plot: clr')
dev.off()


# PRE-STATS -- removing non-rare asvs and transforming --------------------
microbe_no_rare <- microbe_combined%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU, DayNight, Organism)%>%
  filter(abundant == "abundant")%>%
  select(-abundant)


greater_than70 <- microbe_combined%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU, DayNight, Organism)%>%
  filter(abundant == "abundant")%>%
  summarize_if(is.numeric, mean)%>%
  arrange(-ra)%>%
  ungroup()%>%
  group_by(DayNight, Organism)%>%
  nest()%>%
  mutate(data = map(data, ~ mutate(.x, cum = cumsum(ra))))%>%
  unnest(data)%>%
  ungroup()%>%
  filter(ra >= 0.02)%>%
  # filter(cum <= 0.7)%>%
  select(OTU)%>%
  unique()


# PRE-STATS CLEANING -- microbe RA data  --------------------------------------------
ra_bigger_TF <- microbe_combined%>%
  filter(abundant == "abundant",
         Timepoint == 'TF')%>%
  select(-c(sample_code, Replicate))%>%
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

# SET SEED and color values----------------------------------------------------------------
set.seed(2005)

org_colors_no_water <- c("#A30029","#669900", "#FF850A", "#9900FF", "#33CC33")

# STATS permutation ANOVA -- Microbe TWO-Way ------------------------------------------
aov_microbe <- microbe_combined%>%
  filter(abundant == "abundant")%>%
  select(-abundant)%>%
  inner_join(ra_bigger_TF%>%
               ungroup()%>%
               select(OTU)%>%
               unique(), by = 'OTU')%>%
  spread(Timepoint, ra)%>%
  group_by(Organism, DayNight, OTU)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         delta = TF-T0,
         delta = case_when(delta < 0 ~ 0,
                           TRUE ~ delta))%>%
  ungroup()%>%
  group_by(OTU)%>%
  nest()%>%
  mutate(anova = map(data, ~ aovp(delta ~ Organism*DayNight, data = .x)%>%
                       tidy()%>%
                       filter(!term == "Residuals")%>%
                       select(term, `Pr.Prob.`)))%>%
  select(-data)%>%
  unnest(anova)


anova_microbe_pvalues <- aov_microbe%>%
  ungroup()%>%  
  # filter(p.value < 0.05)%>%
  mutate(FDR = p.adjust(.$`Pr.Prob.`, method = "BH"),
         term = gsub('1', '', term))%>%
  filter(FDR < 0.05)

organism_significant_microbes <- (anova_microbe_pvalues%>%
                                    filter(term == "Organism" | term == 'Organism:DayNight'))$OTU%>%
  unique()%>%
  as.vector()

DayNight_significant_microbes <- (anova_microbe_pvalues%>%
                                    filter(term == "Organism:DayNight"))$OTU%>%
  as.vector()

different_microbes <- anova_microbe_pvalues$OTU%>%
  unique()%>%
  as.vector()


# VIZUALIZATIONS -- asvPlots ------------------------------------------
asvPlots <- microbe_combined%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU)%>%
  filter(max(ra) > 0.02)%>%
  ungroup()%>%
  separate(OTU, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";", remove = FALSE)%>%
  unite(lab, c('Family', 'Genus'), sep = ' ')%>%
  group_by(OTU, lab)%>%
  nest()%>%
  mutate(data = map(data, ~ggplot(.x, aes(Organism, ra, fill = Organism)) + 
                                 geom_boxplot() +
                                 geom_point() +
                                 scale_fill_manual(values = c(org_colors_no_water, 'blue')) +
                                 scale_color_manual(values = c('grey', 'grey', 'grey', 'grey', 'grey', 'grey')) +
                                 facet_wrap(~DayNight, nrow = 1) +
                                 labs(y = 'Relative Abundance', x = 'Organism', title = lab) +
                                 # ylim(0.0, 0.8) +
                                 gen_theme()))

pdf('plots/asvDiel.pdf', width = 15, height = 10)
asvPlots$data
dev.off()


# VIZUALIZATIONS -- Asv summary plot -------------------------------------
aovMicrobeSummary <- anova_microbe_pvalues%>%
  select(-Pr.Prob.)%>%
  spread(term, FDR)

asvSummaryPlot <- microbe_combined%>%
  filter(Timepoint == 'TF', 
         OTU %in% different_microbes)%>%
  inner_join(ra_bigger_TF, by = c('OTU', 'DayNight', 'Organism'))%>%
  group_by(OTU, DayNight, Organism)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  filter(Organism != 'CCA')%>%
  group_by(OTU)%>%
  mutate(maxRa = max(ra),
         size = case_when(maxRa <= 0.05 ~ 'ra <= 0.05',
                          maxRa > 0.05 & maxRa <= 0.2 ~ '0.5 < ra <= 0.2',
                          maxRa > 0.2 & maxRa <= 0.3 ~ '0.2 < ra <= 0.3',
                          maxRa > 0.3  ~ 'ra > 0.3'),
         size = as.factor(size),
         # waRa = ra - `Water control`,
         # waRa = case_when(ra < 0 ~ 0,
         #                TRUE ~ waRa),
         coral = case_when(Organism == 'Pocillopora verrucosa' | Organism == 'Porites lobata' ~ ra,
                           TRUE ~ NA_real_),
         algae = case_when(Organism == 'Dictyota' | Organism == 'Turf' ~ ra,
                           TRUE ~ NA_real_),
         night = case_when(DayNight == 'Night' ~ ra,
                           TRUE ~ NA_real_),
         day = case_when(DayNight == 'Day' ~ ra,
                         TRUE ~ NA_real_),
         #max values
         maxCoral = max(coral, na.rm = TRUE),
         maxCoral = as.numeric(gsub(-Inf, 0, maxCoral)),
         maxAlgae = max(algae, na.rm = TRUE),
         maxAlgae = as.numeric(gsub(-Inf, 0, maxAlgae)),
         maxNight = max(night, na.rm = TRUE),
         maxNight = as.numeric(gsub(-Inf, 0, maxNight)),
         maxDay = max(day, na.rm = TRUE),
         maxDay = as.numeric(gsub(-Inf, 0, maxDay)),
         maxVal = max(ra),
         maxOrg = case_when(maxVal == ra ~ Organism,
                            TRUE ~ NA_character_))%>%
  select(OTU, maxRa, size, maxCoral:maxOrg)%>%
  unique()%>%
  filter(!is.na(maxOrg))%>%
  separate(OTU, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";", remove = FALSE)%>%
  mutate(xVal = (maxAlgae - maxCoral)/maxVal,
         yVal = (maxDay - maxNight)/maxVal,
         label = case_when(maxRa > 0.3 ~ Genus,
                           maxRa > 0.01 & abs(xVal) > 0.5 ~ Genus,
                           TRUE ~ NA_character_),
         label = case_when(label == 'uncultured' ~ paste('Uncultured', Family),
                           TRUE ~ label),
         hjust = case_when(xVal < 0 ~ -0,
                           TRUE ~ 0),
         vjust = case_when(yVal < 0 ~ -0.01,
                           TRUE ~ 0.01),
         shape = case_when(yVal < 0 ~ emoji('crescent_moon'),
                           TRUE ~ '\u2600'))%>%
  left_join(aovMicrobeSummary, by = 'OTU')%>%
  mutate(maxOrg = case_when(is.na(Organism) & is.na(`Organism:DayNight`) ~ NA_character_,
                            TRUE ~ maxOrg),
         shape = case_when(is.na(DayNight) & is.na(`Organism:DayNight`) ~ '\u2022',
                           TRUE ~ shape))
         # label = case_when(label %like% '%Rhodobacter'))

asvSummarySigs <- asvSummaryPlot%>%
  filter(shape != '\u2022')

asvSummaryNonDayNightSigs <- asvSummaryPlot%>%
  filter(shape == '\u2022')

  
asvNonSigSummary <- microbe_combined%>%
  filter(Timepoint == 'TF', 
         !OTU %in% different_microbes)%>%
  group_by(OTU, DayNight, Organism)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  filter(Organism != 'CCA')%>%
  group_by(OTU)%>%
  mutate(maxRa = max(ra),
         size = case_when(maxRa <= 0.05 ~ 'ra <= 0.05',
                          maxRa > 0.05 & maxRa <= 0.2 ~ '0.5 < ra <= 0.2',
                          maxRa > 0.2 & maxRa <= 0.3 ~ '0.2 < ra <= 0.3',
                          maxRa > 0.3  ~ 'ra > 0.3'),
         size = as.factor(size),
         coral = case_when(Organism == 'Pocillopora verrucosa' | Organism == 'Porites lobata' ~ ra,
                           TRUE ~ NA_real_),
         algae = case_when(Organism == 'Dictyota' | Organism == 'Turf' ~ ra,
                           TRUE ~ NA_real_),
         night = case_when(DayNight == 'Night' ~ ra,
                           TRUE ~ NA_real_),
         day = case_when(DayNight == 'Day' ~ ra,
                         TRUE ~ NA_real_),
         #max values
         maxCoral = max(coral, na.rm = TRUE),
         maxCoral = as.numeric(gsub(-Inf, 0, maxCoral)),
         maxAlgae = max(algae, na.rm = TRUE),
         maxAlgae = as.numeric(gsub(-Inf, 0, maxAlgae)),
         maxNight = max(night, na.rm = TRUE),
         maxNight = as.numeric(gsub(-Inf, 0, maxNight)),
         maxDay = max(day, na.rm = TRUE),
         maxDay = as.numeric(gsub(-Inf, 0, maxDay)),
         maxVal = max(ra),
         maxOrg = 'Non significant')%>%
  select(OTU, maxRa, size, maxCoral:maxOrg)%>%
  unique()%>%
  filter(!is.na(maxOrg),
         maxRa > 0)%>%
  separate(OTU, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";", remove = FALSE)%>%
  mutate(xVal = (maxAlgae - maxCoral)/maxVal,
         yVal = (maxDay - maxNight)/maxVal,
         shape = case_when(yVal < 0 ~ emoji('crescent_moon'),
                           TRUE ~ '\u2600'))



pdf('plots/asvSummary.pdf', width = 15, height = 12)
asvSummarySigs%>%
  ggplot(aes(xVal, yVal, color = maxOrg, size = maxRa)) +
  geom_text(aes(label = shape), family= 'OpenSansEmoji') +
  geom_text_repel(aes(label = label), nudge_y = 0.04, nudge_x = 0.05, cex = 5, max.overlaps = 20) +
  geom_point(data = asvSummaryNonDayNightSigs, aes(xVal, yVal, color = maxOrg, size = maxRa)) +
  geom_point(data = asvNonSigSummary, aes(xVal, yVal, color = 'Water control', size = maxRa, alpha = 0.7)) +
  scale_color_manual(values = c("#669900", "#FF850A", "#9900FF", "#33CC33", 'Grey')) +
  scale_size_continuous(range = c(3,15), breaks = c(0.005, 0.01, 0.1, 0.3, 0.5)) +
  guides(size = guide_legend(title = 'Relative Abundance', override.aes = aes(label = '\u2022')),
         color  = guide_legend(title = 'Organism', override.aes = aes(label = '\u2022'))) +
  labs(color = 'Organism', 
       x = bquote(atop(Organism ~specificity, ~((Algal ~max- ~Coral ~max)/(~Overall ~max)))),
       y= bquote(atop(Diel ~specificity, ~((Day ~max- ~Night ~max)/(~Overall ~max))))) +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic"),
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25))
dev.off()


# STATS -- Diversity and Richness -----------------------------------------
asvDiversity <- microbe_combined%>%
  filter(Timepoint == 'TF', 
         OTU %in% different_microbes)%>%
  inner_join(ra_bigger_TF, by = c('OTU', 'DayNight', 'Organism'))%>%
  spread(OTU, ra)%>%
  gather(OTU, ra, 8:ncol(.))%>%
  mutate(ra = replace_na(ra, 0))%>%
  spread(OTU, ra)%>%
  select(Organism, DayNight, 8:ncol(.))%>%
  group_by(Organism, DayNight)%>%
  nest()%>%
  mutate(shannon = map(data, ~ diversity(.x, index = 'shannon')),
         specN = map(data, ~ rowSums(.x > 0, na.rm = TRUE)))%>%
  unnest(c(shannon, specN))%>%
  mutate(eveness = shannon/log(specN))

asvDiversity%>% 
  select(-data)%>% 
  mutate(orgName = case_when(Organism %like% 'Poc%' ~ 'Coral',
                             Organism %like% 'Por%' ~ 'Coral',
                             Organism %like% 'Dic%' ~ 'FAlgae',
                             Organism %like% 'Turf%' ~ 'FAlgae'))%>%
  group_by(orgName, DayNight)%>% 
  mutate(sd = sd(eveness))%>%
  summarize_if(is.numeric, mean)

waterMeans <- microbe_combined%>%
  filter(Timepoint == 'TF', 
         OTU %in% different_microbes)%>%
  inner_join(ra_bigger_TF%>%
               ungroup()%>%
               select(OTU)%>%
               unique(), by = 'OTU')%>%
  filter(Organism == 'Water control')%>%
  select(OTU, DayNight, ra)%>%
  group_by(OTU, DayNight)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  rename(waterRa = ra)

asvMean <- microbe_combined%>%
  filter(Timepoint == 'TF', 
         OTU %in% different_microbes)%>%
  inner_join(ra_bigger_TF, by = c('OTU', 'DayNight', 'Organism'))%>%
  filter(ra > 0)%>%
  left_join(waterMeans, by = c('DayNight', 'OTU'))%>%
  mutate(waterSub = ra - waterRa)

pdf('plots/MicrobeDiversity.pdf', width = 15, height = 10)
asvDiversity%>%
  ggplot(aes(DayNight, shannon, color = Organism)) +
  geom_boxplot() +
  facet_wrap(~Organism, nrow = 1) +
  gen_theme() +
  labs(y = 'Shannon Diveristy', x = 'Diel Cycle') +
  scale_color_manual(values = org_colors_no_water)

asvDiversity%>%
  ggplot(aes(DayNight, eveness, color = Organism)) +
  geom_boxplot() +
  facet_wrap(~Organism, nrow = 1) +
  gen_theme() +
  labs(y = "Pielou's Eveneess", x = 'Diel Cycle') +
  scale_color_manual(values = org_colors_no_water)
dev.off()

asvMean%>%
  separate(OTU, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";", remove = FALSE)%>%
  unite(microbe, c('Genus', 'OTU_id'), sep = ' ')%>%
  ggplot(aes(DayNight, ra, color = Organism)) +
  geom_violin() +
  geom_point(size = 4) + 
  facet_wrap(~Organism) +
  gen_theme() +
  labs(y = "Relative Abundance") +
  scale_color_manual(values = org_colors_no_water)
  
evenessStats <- asvDiversity%>%
  group_by(Organism)%>%
  nest()%>%
  mutate(data = map(data, ~ lm(eveness~DayNight, data = .x)%>%
                      tidy()%>%
                      filter(!term %like% '%Intercept%')))%>%
  unnest(data)
 

# VIZUALIZATIONS -- Hierarchical cluster matrix and asv table----------------------------------------
hc_microbe <- microbe_combined%>%
  # select(-c(sample_code, reads, numOtus))%>%
  # ungroup()%>%
  filter(Timepoint == "TF")%>%
  select(c(Organism, DayNight, Replicate, OTU, ra))%>%
  group_by(OTU)%>%
  filter(max(ra) > 0.01)%>%
  unite(sample, c("DayNight", "Organism", "Replicate"), sep = "_")%>%
  separate(OTU, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  unite(OFGO, c("Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  select(-c(Kingdom, Phylum, Class))%>%
  group_by(OFGO)%>%
  mutate(zscore = zscore(ra))

# write_csv(hc_microbe%>%
#             select(sample, everything()), "analysis/microbe_hc_df.csv")

pdf('plots/bacterioplanktonHeatmap.pdf', width = 15, height = 10)
hc_microbe%>%
  # gather(OFGO, zscore, 2:ncol(.))%>%
  ggplot(aes(sample, OFGO, fill = zscore
             )) +
  geom_tile() +
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  # scale_y_discrete(limits = rev(levels(depletion_levels))) +
  # scale_x_discrete(limits = levels(organism_levels)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

microbe_table <- hc_microbe%>%
  separate(sample, c("DayNight", "Organism", "Replicate"), sep = "_")%>%
  group_by(DayNight, Organism, OFGO)%>%
  summarize_if(is.numeric, mean)%>%
  unite(sample, c('DayNight', 'Organism'), sep = '_')%>%
  separate(OFGO, c("Order", "Family", "Genus", "OTU"), sep = ";")%>%
  select(-zscore)%>%
  spread(sample, ra)%>%
  inner_join(aov_microbe%>%
               ungroup()%>%
               separate(OTU, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"), sep = ";")%>%
               select(-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))%>%
               mutate(FDR = p.adjust(`Pr.Prob.`, method = "BH"),
                      term = gsub('1', '', term))%>%
               select(-`Pr.Prob.`)%>%
               group_by(OTU)%>%
               filter(min(FDR) < 0.05)%>%
               ungroup()%>%
               spread(term, FDR), by = 'OTU')%>%
  select(Order, Family, Genus, OTU, Organism, DayNight, `Organism:DayNight`, everything())
  
write_csv(microbe_table, 'analysis/microbesBiggerThan2percent.csv')


# VIZUALIZATIONS -- PCOA --------------------------------------------------
unifracHead <- unifracRaw%>%
  rename(sample = Groups,
         uW = WScore)%>%
  filter(sample %like% "D_%",
         sample %like% "%-D%")%>%
  rename(sample_code = sample)%>%
  separate(sample_code, c('sample1', 'sample2'), sep = '-')%>%
  filter(!sample1 %like% '%D_PV_2_TFN_SA504_SC704%',
         !sample1 %like% '%D_PV_4_TFD%',
         !sample2 %like% '%D_PV_2_TFN_SA504_SC704%',
         !sample2 %like% '%D_PV_4_TFD%',
         !sample2 %like% '%Dorcierr_unlabeled%')%>%
  mutate(sample1 = case_when(sample1 == "Dorcierr_D_DT_1_TFD" ~ "D_DT_1_TFD",
                             sample1 == "DORCIERR_D_WA_2_TFN" ~ "D_WA_2_TFN",
                             # Group == "D_PV_2_TFN_SA504_SC704" ~ "D_PV_2_TFN",
                             sample1 == "D_PV_2_TFN_SA503_SC704" ~ "D_PL_3_TFN",
                             sample1 == "D_WA_4_TFN_SA503_SC703" ~ "D_PV_3_TFN",
                             sample1 == "D_WA_4_TFN_SA504_SC703" ~ "D_WA_4_TFN",
                             TRUE ~ as.character(sample1)),
         sample2 = case_when(sample2 == "Dorcierr_D_DT_1_TFD" ~ "D_DT_1_TFD",
                             sample2 == "DORCIERR_D_WA_2_TFN" ~ "D_WA_2_TFN",
                             # Group == "D_PV_2_TFN_SA504_SC704" ~ "D_PV_2_TFN",
                             sample2 == "D_PV_2_TFN_SA503_SC704" ~ "D_PL_3_TFN",
                             sample2 == "D_WA_4_TFN_SA503_SC703" ~ "D_PV_3_TFN",
                             sample2 == "D_WA_4_TFN_SA504_SC703" ~ "D_WA_4_TFN",
                             TRUE ~ as.character(sample2)))%>%
  select(sample1, sample2, uW)

unifracTail <- unifracHead%>%
  rename(sample1 = 2,
         sample2 = 1)%>%
  bind_rows(unifracHead)%>%
  spread(sample2, uW)%>%
  column_to_rownames(var = 'sample1')%>%
  replace(is.na(.), 0)%>%
  as.matrix()

microbePcoa <- pcoa(unifracTail)

pcoaPlots <- microbePcoa$vectors%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample")%>%
  nameSamples()%>%
  group_by(Organism, Timepoint, DayNight)%>%
  mutate(xerr = sd(Axis.1), yerr = sd(Axis.2))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(shape = case_when(DayNight == 'Day' ~ '\u2600',
                           DayNight == 'Night' ~ emoji('crescent_moon')))

pdf('plots/pcoaMicrobes.pdf', width = 15, height = 10)
pcoaPlots%>%
  filter(DayNight == 'Day')%>%
  ggplot(aes(x = Axis.1, y = Axis.2, color = Organism)) +
  geom_text(aes(label = shape), cex = 13, family= 'OpenSansEmoji') +
  geom_errorbar(aes(ymin = Axis.2 - yerr, ymax = Axis.2 + yerr), size = 1) +
  geom_errorbarh(aes(xmin = Axis.1 - xerr, xmax = Axis.1 + xerr), size = 1) +
  geom_line(aes(group = Organism)) +
  scale_color_manual(values = c(org_colors_no_water, 'blue')) +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic"),
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25)) +
  xlab(str_c("Axis 1", " (", round(microbePcoa$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ylab(str_c("Axis 2", " (", round(microbePcoa$values$Relative_eig[2], digits = 4)*100, "%)", sep = ""))

pcoaPlots%>%
  filter(DayNight == 'Night')%>%
  ggplot(aes(x = Axis.1, y = Axis.2, color = Organism)) +
  geom_text(aes(label = shape), cex = 13, family= 'OpenSansEmoji') +
  geom_errorbar(aes(ymin = Axis.2 - yerr, ymax = Axis.2 + yerr), size = 1) +
  geom_errorbarh(aes(xmin = Axis.1 - xerr, xmax = Axis.1 + xerr), size = 1) +
  geom_line(aes(group = Organism)) +
  scale_color_manual(values = c(org_colors_no_water, 'blue')) +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic"),
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25)) +
  xlab(str_c("Axis 1", " (", round(microbePcoa$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ylab(str_c("Axis 2", " (", round(microbePcoa$values$Relative_eig[2], digits = 4)*100, "%)", sep = ""))
dev.off()

# STATS --Unifrac distances ------------------------------------------------------------
ylab <- bquote(atop(Weighted ~unifrac ~distances, ~(T[0]:T[F])))

unifracValues <- unifracTail%>%
  as.data.frame()%>%
  rownames_to_column(var = 'start')%>%
  filter(start %like% '%T0%')%>%
  gather(sample, unifracDissimilarity, 2:ncol(.))%>%
  nameSamples()%>%
  separate(start, c("startExperiment", "startOrganism", "startReplicate", "startTimepoint"), sep = "_", remove = FALSE)%>%
  mutate(startOrganism = case_when(startOrganism == "CC" ~ "CCA",
                                   startOrganism == "DT" ~ "Dictyota",
                                   startOrganism == "PL" ~ "Porites lobata",
                                   startOrganism == "PV" ~ "Pocillopora verrucosa",
                                   startOrganism == "TR" ~ "Turf",
                                   startOrganism == "WA" ~ "Water control",
                                   startOrganism == "IN" ~ "Influent",
                                   startOrganism == "OF" ~ "Offshore",
                              TRUE ~ as.character(startOrganism)))%>%
  filter(!startOrganism %in% c('Influent', 'Offshore'))%>%
  separate(startTimepoint, c("startTimepoint", "startDayNight"), sep = 2)%>%
  mutate(startDayNight = case_when(startDayNight == "D" ~ "Day",
                              TRUE ~ "Night"))%>%
  filter(Timepoint == 'TF',
         DayNight == startDayNight,
         startOrganism == Organism)%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(mannWhitneyCCA = map(data, ~ (wilcox.test(unifracDissimilarity ~ Organism, 
                                                  subset = Organism %in% c('Water control', 'CCA'), 
                                                  data = .x))$p.value),
         mannWhitneyDictyota = map(data, ~ (wilcox.test(unifracDissimilarity ~ Organism, 
                                                   subset = Organism %in% c('Water control', 'Dictyota'), 
                                                   data = .x))$p.value),
         mannWhitneyPoc = map(data, ~ (wilcox.test(unifracDissimilarity ~ Organism, 
                                                   subset = Organism %in% c('Water control', 'Pocillopora verrucosa'), 
                                                   data = .x))$p.value),
         mannWhitneyPor = map(data, ~ (wilcox.test(unifracDissimilarity ~ Organism, 
                                                   subset = Organism %in% c('Water control', 'Porites lobata'), 
                                                   data = .x))$p.value),
         mannWhitneyTurf = map(data, ~ (wilcox.test(unifracDissimilarity ~ Organism, 
                                                   subset = Organism %in% c('Water control', 'Turf'), 
                                                   data = .x))$p.value),
         plot = map(data, ~ ggplot(.x, aes(Organism, unifracDissimilarity))+
                      geom_boxplot()+
                      gen_theme()+
                      labs(y = ylab)))

pdf('plots/microbeUnifrac.pdf', width = 15, height = 10)
unifracValues$plot
dev.off()

# STATS -- PERMANOVAs --------------------------------------------------------
permanova_organism <- unifracTail%>%
  as.data.frame()%>%
  rownames_to_column(var = 'sample')%>%
  nameSamples()%>%
  column_to_rownames(var = 'sample')%>%
  unite(group, c('Organism', 'Timepoint', 'DayNight'), sep = '_', remove = FALSE)

permanova_matrix <- permanova_organism[7:ncol(permanova_organism)]%>%
  as.matrix()

pairwisePermanovaPVals <- permanova_matrix%>%
  pairwiseAdonis::pairwise.adonis(., permanova_organism$group, p.adjust.m = 'BH')%>%
  filter(!pairs %like% '%T0%')%>%
  mutate(fdr = p.adjust(p.value, method = 'BH'))

permanova_DayNight <- permanova_organism%>%
  adonis(.[7:ncol(.)] ~ DayNight, ., perm = 1000)

permanova_twoway <- permanova_organism%>%
  adonis(.[7:ncol(.)] ~ Organism*DayNight, ., perm = 1000)


