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
library(permuco)


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


# LOADING -- Dataframes ---------------------------------------------------
# 16s rRNA sequences
microbe_abundance_raw <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/microbes/MCR2017.16S.Nelson.Pipeline.October2019/abundance_table_100.shared.tsv")
microbe_taxonomy <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/microbes/MCR2017.16S.Nelson.Pipeline.October2019/annotations_100.taxonomy.tsv")

# PRE-STATS CLEANING -- Microbes and pre-filtering OTUs for abundance-----------------------------------------------
microbe_combined <- microbe_abundance_raw%>%
  select(-1)%>%
  mutate(Group = case_when(Group == "Dorcierr_D_DT_1_TFD" ~ "D_DT_1_TFD",
                           Group == "DORCIERR_D_WA_2_TFN" ~ "D_WA_2_TFN",
                           # Group == "D_PV_2_TFN_SA504_SC704" ~ "D_PV_2_TFN",
                           Group == "D_PV_2_TFN_SA503_SC704" ~ "D_PL_3_TFN",
                           Group == "D_WA_4_TFN_SA503_SC703" ~ "D_PV_3_TFN",
                           Group == "D_WA_4_TFN_SA504_SC703" ~ "D_WA_4_TFN",
                           TRUE ~ as.character(Group)))%>%
  filter(Group %like% "%D_%",
         Group != 'D_PV_2_TFN_SA504_SC704')%>%
  rename(sample_code = Group)%>%
  gather(OTU, reads, 3:ncol(.))%>%
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
  mutate(reads = case_when(reads == 0 ~ 1/numOtus,
                           TRUE ~ as.numeric(reads)))%>%
  group_by(sample_code)%>%
  mutate(ra = reads/sum(reads))%>%
  ungroup()%>%
  group_by(OTU)%>%
  mutate(abundant = case_when(max(ra) > 0.05 | sum(ra > 0.001) >= 3 ~ "abundant",
                              TRUE ~ "rare"))

# PRE-STATS -- qqplots ----------------------------------------------------
pdf('plots/qqplots.pdf', width = 15, height = 10)
microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU, DayNight, Organism)%>%
  filter(abundant == "abundant")%>% 
  ggplot(aes(ra)) + 
  geom_histogram(bins = 100) +
  labs(title = 'Distribution of relative abundances', y = 'Count', x = 'Relative Abundance') +
  gen_theme()

car::qqPlot((microbe_combined%>%
               select(-c(sample_code, reads, numOtus))%>%
               ungroup()%>%
               filter(Timepoint == "TF")%>%
               group_by(OTU, DayNight, Organism)%>%
               filter(abundant == "abundant"))$ra, 
            ylab = "Relative Abundance", xlab = "Normal quantiles",
            main = 'QQ-plot: Relative Abundance') 

car::qqPlot((microbe_combined%>%
               select(-c(sample_code, reads, numOtus))%>%
               ungroup()%>%
               filter(Timepoint == "TF")%>%
               group_by(OTU, DayNight, Organism)%>%
               filter(abundant == "abundant")%>%
               mutate(asin = asin(sqrt(ra))))$asin, 
            ylab = "Relative Abundance", xlab = "Normal quantiles",
            main = 'QQ-plot: Angular Transformation')

car::qqPlot((microbe_combined%>%
               select(-c(sample_code, reads, numOtus))%>%
               ungroup()%>%
               filter(Timepoint == "TF")%>%
               group_by(OTU, DayNight, Organism)%>%
               filter(abundant == "abundant")%>%
               mutate(log = log10(ra)))$log, 
            ylab = "Relative Abundance", xlab = "Normal quantiles",
            main = 'QQ-plot: log10')

car::qqPlot((microbe_combined%>%
               select(-c(sample_code, reads, numOtus))%>%
               ungroup()%>%
               filter(Timepoint == "TF")%>%
               group_by(OTU, DayNight, Organism)%>%
               filter(abundant == "abundant")%>%
               mutate(log = Logit(ra)))$log, 
            ylab = "Relative Abundance", xlab = "Normal quantiles",
            main = 'QQ-plot: Logit')


car::qqPlot((microbe_combined%>%
               select(-c(sample_code, reads, numOtus))%>%
               ungroup()%>%
               filter(Timepoint == "TF")%>%
               group_by(OTU, DayNight, Organism)%>%
               filter(abundant == "abundant")%>%
               mutate(clr = compositions::clr(ra)))$clr, 
            ylab = "Relative Abundance", xlab = "Normal quantiles",
            main = 'QQ-plot: clr')
dev.off()


# PRE-STATS -- removing non-rare asvs and transforming --------------------
microbe_no_rare <- microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU, DayNight, Organism)%>%
  filter(abundant == "abundant")%>%
  mutate(asin = asin(sqrt(ra)))%>%
  select(-c(ra, abundant))


greater_than70 <- microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
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

# SET SEED and color values----------------------------------------------------------------
set.seed(2005)

org_colors_no_water <- c("#A30029","#669900", "#FF850A", "#9900FF", "#33CC33")

# STATS permutation ANOVA -- Microbe TWO-Way ------------------------------------------
aov_microbe <- microbe_no_rare%>%
  inner_join(ra_bigger_TF%>%
               ungroup()%>%
               select(OTU)%>%
               unique(), by = 'OTU')%>%
  group_by(OTU)%>%
  nest()%>%
  mutate(anova = map(data, ~ aovperm(asin ~ Organism*DayNight, data = .x, np = 10000)%>%
                       summary()%>%
                       rownames_to_column(var = 'term')%>%
                       filter(!term == "Residuals")%>%
                       select(term, `permutation P(>F)`)))%>%
  select(-data)%>%
  unnest(anova)


anova_microbe_pvalues <- aov_microbe%>%
  ungroup()%>%  
  # filter(p.value < 0.05)%>%
  mutate(FDR = p.adjust(.$`permutation P(>F)`, method = "BH"))%>%
  filter(FDR < 0.05)

organism_significant_microbes <- (anova_microbe_pvalues%>%
                                    filter(term == "Organism"))$OTU%>%
  as.vector()

DayNight_significant_microbes <- (anova_microbe_pvalues%>%
                                    filter(term == "Organism:DayNight"))$OTU%>%
  as.vector()

different_microbes <- anova_microbe_pvalues$OTU%>%
  unique()%>%
  as.vector()


# STATS -- t-tests -----------------------------
otu_ttest <- microbe_no_rare%>%
  filter(OTU %in% different_microbes)%>%
  select(-Experiment)

otu_pairwise <- otu_ttest%>%
  mutate(Organism = factor(Organism),
         Organism = fct_relevel(Organism, "Water control"))%>%
  # filter(OTU %like any% c('Otu0005', 'Otu0003', 'Otu0006'))%>%
  group_by(OTU, Organism)%>%
  nest()%>%
  mutate(data = map(data, ~ pairwise.wilcox.test(.x$asin, .x$DayNight,
                                                 p.adjust.method = "BH")%>%
                      tidy()))%>%
  unnest(data)%>%
  ungroup()%>%
  filter(p.value < 0.05)


# STATS -- PERMANOVAs --------------------------------------------------------
permanova_df <- microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF")%>%
  inner_join(greater_than70, by = 'OTU')%>%
  select(c(Organism, DayNight, Replicate, OTU, ra))%>%
  # unite(sample, c("DayNight", "Organism", "Replicate"), sep = "_")%>%
  left_join(microbe_taxonomy, by = "OTU")%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  unite(OFGO, c("Order", "Family", "Genus", "OTU"), sep = ";")%>%
  # group_by(OFGO)%>%
  # mutate(zscore = zscore(ra))%>%
  # ungroup()%>%
  mutate(asin = asin(sqrt(ra)))%>%
  select(-c(Size, OTU_id, Kingdom, Phylum, Class, ra))%>%
  spread(OFGO, asin)

permanova_organism <- permanova_df%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(data = map(data, ~ pairwiseAdonis::pairwise.adonis(.x[3:ncol(.)], .x$Organism, sim.function = 'vegdist', 
                                                            sim.method =  'bray', p.adjust.m = 'BH')))%>%
  unnest(data)%>%
  ungroup()%>%
  mutate(fdr = p.adjust(p.value, method = 'BH'))

permanova_DayNight <- permanova_df%>%
  adonis(.[4:ncol(.)] ~ DayNight, ., perm = 1000, method = "bray", na.rm = TRUE)

# VIZUALIZATIONS -- Rhodobacter ------------------------------------------
pdf('plots/rhodobacterDiel.pdf', width = 15, height = 10)
microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU, DayNight, Organism)%>%
  filter(OTU == 'Otu0005')%>%
  ggplot(aes(DayNight, ra, fill = Organism)) + 
  geom_boxplot() +
  geom_point() +
  scale_fill_manual(values = c(org_colors_no_water, 'blue')) +
  scale_color_manual(values = c('grey', 'grey', 'grey', 'grey', 'grey', 'grey')) +
  facet_wrap(~Organism, nrow = 1) +
  labs(y = 'Relative Abundance', x = 'Diel cycle', title = 'Unclassified rhodobacteraceae') +
  ylim(0.0,0.7) +
  gen_theme()

microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU, DayNight, Organism)%>%
  filter(OTU == 'Otu0003')%>%
  ggplot(aes(DayNight, ra, fill = Organism)) + 
  geom_boxplot() +
  geom_point() +
  scale_fill_manual(values = c(org_colors_no_water, 'blue')) +
  scale_color_manual(values = c('grey', 'grey', 'grey', 'grey', 'grey', 'grey')) +
  facet_wrap(~Organism, nrow = 1) +
  labs(y = 'Relative Abundance', x = 'Diel cycle', title = 'Nautella sp.') +
  ylim(0.0,0.7) +
  gen_theme()

microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU, DayNight, Organism)%>%
  filter(OTU == 'Otu0006')%>%
  ggplot(aes(DayNight, ra, fill = Organism)) + 
  geom_boxplot() +
  geom_point() +
  scale_fill_manual(values = c(org_colors_no_water, 'blue')) +
  scale_color_manual(values = c('grey', 'grey', 'grey', 'grey', 'grey', 'grey')) +
  facet_wrap(~Organism, nrow = 1) +
  labs(y = 'Relative Abundance', x = 'Diel cycle', title = 'Thalassobious sp.') +
  ylim(0.0,0.7) +
  gen_theme()

microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU, DayNight, Organism)%>%
  filter(OTU == 'Otu0002')%>%
  ggplot(aes(DayNight, ra, fill = Organism)) + 
  geom_boxplot() +
  geom_point() +
  scale_fill_manual(values = c(org_colors_no_water, 'blue')) +
  scale_color_manual(values = c('grey', 'grey', 'grey', 'grey', 'grey', 'grey')) +
  facet_wrap(~Organism, nrow = 1) +
  labs(y = 'Relative Abundance', x = 'Diel cycle', title = 'HIMB11 sp.') +
  ylim(0.0,0.7) +
  gen_theme()

microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU, DayNight, Organism)%>%
  filter(OTU == 'Otu0004')%>%
  ggplot(aes(DayNight, ra, fill = Organism)) + 
  geom_boxplot() +
  geom_point() +
  scale_fill_manual(values = c(org_colors_no_water, 'blue')) +
  scale_color_manual(values = c('grey', 'grey', 'grey', 'grey', 'grey', 'grey')) +
  facet_wrap(~Organism, nrow = 1) +
  labs(y = 'Relative Abundance', x = 'Diel cycle', title = 'Flavobacter cryomorphaceae') +
  ylim(0.0,0.7) +
  gen_theme()

microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU, DayNight, Organism)%>%
  filter(OTU == 'Otu0010')%>%
  ggplot(aes(DayNight, ra, fill = Organism)) + 
  geom_boxplot() +
  geom_point() +
  scale_fill_manual(values = c(org_colors_no_water, 'blue')) +
  scale_color_manual(values = c('grey', 'grey', 'grey', 'grey', 'grey', 'grey')) +
  facet_wrap(~Organism, nrow = 1) +
  labs(y = 'Relative Abundance', x = 'Diel cycle', title = 'Thalassobious sp.') +
  ylim(0.0,0.7) +
  gen_theme()

dev.off()

# VIZUALIZATIONS -- Hierarchical cluster matrix ----------------------------------------
hc_microbe <- microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF")%>%
  inner_join(greater_than70, by = 'OTU')%>%
  select(c(Organism, DayNight, Replicate, OTU, ra))%>%
  unite(sample, c("DayNight", "Organism", "Replicate"), sep = "_")%>%
  left_join(microbe_taxonomy, by = "OTU")%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  unite(OFGO, c("Order", "Family", "Genus", "OTU"), sep = ";")%>%
  select(-c(Size, OTU_id, Kingdom, Phylum, Class))%>%
  group_by(OFGO)%>%
  mutate(zscore = zscore(ra))

# write_csv(hc_microbe%>%
#             select(sample, everything()), "analysis/microbe_hc_df.csv")

pdf('bacterioplanktonHeatmap.pdf', width = 15, height = 10)
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
  left_join(aov_microbe%>%
              ungroup()%>%  
              mutate(FDR = p.adjust(`permutation P(>F)`, method = "BH"))%>%
              select(-`permutation P(>F)`)%>%
              spread(term, FDR), by = 'OTU')%>%
  select(Order, Family, Genus, OTU, DayNight, Organism, `Organism:DayNight`, everything())
  
write_csv(microbe_table, 'analysis/microbesBiggerThan2percent.csv')

# VIZUALIZATIONS -- PCOA --------------------------------------------------
pcoaMicrobe <- microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF" | Timepoint == 'T0' & Organism == 'Water control',
         DayNight == 'Day')%>%
  group_by(OTU, DayNight, Organism, Timepoint)%>%
  filter(abundant == "abundant")%>%
  mutate(asin = asin(sqrt(ra)))%>%
  select(-c(ra, abundant))%>%
  # filter(OTU %in% different_microbes)%>%
  select(c(Organism, DayNight, Replicate, OTU, Timepoint, asin))%>%
  unite(sample, c("Organism", "DayNight", "Replicate", "Timepoint"), sep = "_")%>%
  left_join(microbe_taxonomy, by = "OTU")%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  unite(OFGO, c("Order", "Family", "Genus", "OTU"), sep = ";")%>%
  select(-c(Size, OTU_id, Kingdom, Phylum, Class))%>%
  # mutate(asin = asin - min(asin))%>% 
  spread(OFGO, asin)%>%
  column_to_rownames(var = "sample")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()

pdf('plots/pcoaMicrobesDay.pdf', width = 15, height = 12)
pcoaMicrobe$vectors%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample")%>%
  separate(sample, c('Organism', 'DayNight', 'Replicate', 'Timepoint'), sep = '_')%>%
  mutate(shape = case_when(DayNight == 'Day' ~ '\u2600',
                           TRUE ~ emoji('crescent_moon')))%>%
  ggplot(., aes(x = Axis.1, y = Axis.2, color = Organism)) +
  geom_text(aes(label = shape), cex = 13, family= 'OpenSansEmoji') +
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
  xlab(str_c("Axis 1", " (", round(pcoaMicrobe$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ylab(str_c("Axis 2", " (", round(pcoaMicrobe$values$Relative_eig[2], digits = 4)*100, "%)", sep = ""))
dev.off()


# STATS --Bray-Curtis distances ------------------------------------------------------------
ylab <- bquote(atop(Bray-curtis ~dissimilarity, ~(T[0] ~Water ~control:T[F] ~Organism)))

brayCurtis <- microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  filter(Timepoint == "TF" | Timepoint == 'T0' & Organism == 'Water control')%>%
  group_by(OTU, DayNight, Organism, Timepoint)%>%
  filter(abundant == "abundant")%>%
  mutate(asin = asin(sqrt(ra)))%>%
  select(-c(ra, abundant))%>%
  # filter(OTU %in% different_microbes)%>%
  select(c(Organism, DayNight, Replicate, OTU, Timepoint, asin))%>%
  unite(sample, c("Organism", "Replicate", "Timepoint"), sep = "_")%>%
  left_join(microbe_taxonomy, by = "OTU")%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  unite(OFGO, c("Order", "Family", "Genus", "OTU"), sep = ";")%>%
  select(-c(Size, OTU_id, Kingdom, Phylum, Class))%>%
  # mutate(asin = asin - min(asin))%>% 
  spread(OFGO, asin)%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(data = map(data, ~ column_to_rownames(.x, var = "sample")%>% 
                      vegdist(na.rm = TRUE)%>%
                      as.matrix()%>%
                      as.data.frame()%>%
                      rownames_to_column(var = 'waterStart')%>%
                      filter(waterStart %like% 'Water control%',
                             waterStart %like% '%T0')%>%
                      gather(sample, brayDissimilarity, 2:ncol(.))%>%
                      separate(sample, c("Organism", "Replicate", "Timepoint"), sep = "_")%>%
                      filter(Timepoint == 'TF')),
         mannWhitneyCCA = map(data, ~ (wilcox.test(brayDissimilarity ~ Organism, 
                                                  subset = Organism %in% c('Water control', 'CCA'), 
                                                  data = .x))$p.value),
         mannWhitneyDictyota = map(data, ~ (wilcox.test(brayDissimilarity ~ Organism, 
                                                   subset = Organism %in% c('Water control', 'Dictyota'), 
                                                   data = .x))$p.value),
         mannWhitneyPoc = map(data, ~ (wilcox.test(brayDissimilarity ~ Organism, 
                                                   subset = Organism %in% c('Water control', 'Pocillopora verrucosa'), 
                                                   data = .x))$p.value),
         mannWhitneyPor = map(data, ~ (wilcox.test(brayDissimilarity ~ Organism, 
                                                   subset = Organism %in% c('Water control', 'Porites lobata'), 
                                                   data = .x))$p.value),
         mannWhitneyTurf = map(data, ~ (wilcox.test(brayDissimilarity ~ Organism, 
                                                   subset = Organism %in% c('Water control', 'Turf'), 
                                                   data = .x))$p.value),
         plot = map(data, ~ ggplot(.x, aes(Organism, brayDissimilarity))+
                      geom_boxplot()+
                      gen_theme()+
                      labs(y = ylab)))

pdf('microbeBrayCurtis.pdf', width = 15, height = 10)
brayCurtis$plot
dev.off()






