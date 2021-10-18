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
library(VennDiagram)


#PCoA, PERMANOVA
library(vegan)
library(ape)

#Visualizations
library(wesanderson)
library(RColorBrewer)
library(gplots)
library(gtable)
library(ggnewscale)
library(ggforce)
library(emojifont)
load.emojifont('OpenSansEmoji.ttf')

#modeling
library(yardstick)
library(glmnet)
library(recipes)
library(foreach)

#Defining functions and removing issues with overlapping function calls
map <- purrr::map
select <- dplyr::select
tidy <- broom::tidy
rename <- dplyr::rename
mutate <- dplyr::mutate

zscore <- function(x) {
  (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

angular_transform <- function(x) {
  asin(sqrt(x))
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
                           Group == "D_PV_2_TFN_SA504_SC704" ~ "D_PV_2_TFN",
                           Group == "D_PV_2_TFN_SA503_SC704" ~ "D_PV_3_TFN",
                           Group == "D_WA_4_TFN_SA503_SC703" ~ "D_PL_3_TFN",
                           Group == "D_WA_4_TFN_SA504_SC703" ~ "D_WA_4_TFN",
                           TRUE ~ as.character(Group)))%>%
  filter(Group %like% "%D_%")%>%
  rename(sample_code = Group)%>%
  gather(OTU, reads, 3:ncol(.))%>%
  # left_join(., microbe_taxonomy, by = "OTU")%>%
  # select(-Size)%>%
  # separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";", remove = FALSE)%>%
  separate(sample_code, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_", remove = FALSE)%>%
  mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                       Experiment == "M" ~ "mordor",
                                       Experiment == "R" ~ "RR3",
                                       TRUE ~ as.character(Experiment)))%>%
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

# SET SEED ----------------------------------------------------------------
set.seed(2005)

# STATS ANOVA -- Microbe TWO-Way ------------------------------------------
aov_microbe <- microbe_no_rare%>%
  group_by(OTU)%>%
  nest()%>%
  mutate(anova = map(data, ~ aov(asin ~ Organism*DayNight, .x)%>%
                       tidy()%>%
                       filter(!term == "Residuals")%>%
                       select(term, p.value)))%>%
  select(-data)%>%
  unnest(anova)

anova_microbe_pvalues <- aov_microbe%>%
  ungroup()%>%  
  filter(p.value < 0.05)%>%
  add_column(FDR = p.adjust(.$p.value, method = "BH"))%>%
  filter(FDR < 0.05)

organism_significant_microbes <- (anova_microbe_pvalues%>%
  filter(term == "Organism"))$OTU%>%
  as.vector()

DayNight_significant_microbes <- (anova_microbe_pvalues%>%
  filter(term != "Organism"))$OTU%>%
  as.vector()

different_microbes <- anova_microbe_pvalues$OTU%>%
  unique()%>%
  as.vector()
  
# PRE-POST-HOC CLEANING -- Microbe Dunnetts and DayNight anova -------------------------------
mic_organism_post_hoc <- microbe_no_rare%>%
  filter(OTU %in% organism_significant_microbes)

daynight_microbe_post_hoc <- microbe_no_rare%>%
  filter(OTU %in% DayNight_significant_microbes)

# STATS POST-HOC -- MICROBES Dunnetts -----------------------------
organism_order_micro <- as.factor(mic_organism_post_hoc$Organism)%>%
  relevel("Water control")%>%
  levels()%>%
  as.vector()

dunnett_microbe_pvals <- mic_organism_post_hoc%>%
  group_by(DayNight, OTU)%>%
  mutate(sum = sum(asin))%>%
  filter(sum != 0)%>%
  select(-sum)%>%
  mutate(Organism = factor(Organism))%>%
  mutate(Organism = fct_relevel(Organism, organism_order_micro))%>%
  nest()%>%
  mutate(dunnett = map(data, ~ aov(asin ~ Organism, .x)%>%
                         glht(linfct = mcp(Organism = "Dunnett"))),
         dunnett_summary = map(dunnett, ~summary(.x)%>%
                                 tidy()))%>%
  select(-c(data,dunnett))%>%
  unnest(dunnett_summary)%>%
  select(-c(4:7))%>%
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
  select(-sum)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(asin ~ DayNight, .x)%>%
                      tidy()))%>%
  unnest(data)%>%
  select(-c(4:7))%>%
  filter(term != "Residuals")%>%
  ungroup()%>%   
  add_column(FDR = p.adjust(.$p.value, method = "BH"))%>%
  filter(FDR < 0.05)


# META-STATS -- microbes --------------------------------------------------
dunnett_micro_analysis <- dunnett_microbe_pvals%>%
  select(-p.value)%>%
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

# VIZUALIZATIONS -- RGB hex codes for orgs --------------------------------
org_colors_no_water <- c("#A30029","#669900", "#FF850A", "#9900FF", "#33CC33")

# META-STATS -- Hierarchical cluster matrix----------------------------------------
hc_microbe <- mic_organism_post_hoc%>%
  filter(OTU %in% different_microbes)%>%
  ungroup()%>%
  group_by(OTU)%>%
  mutate(zscore = zscore(asin))%>%
  ungroup()%>%
  select(c(Organism, DayNight, Replicate, OTU, zscore))%>%
  unite(sample, c("Organism", "DayNight", "Replicate"), sep = "_")%>%
  left_join(microbe_taxonomy, by = "OTU")%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  unite(OFGO, c("Order", "Family", "Genus", "OTU"), sep = ";")%>%
  select(-c(Size, OTU_id, Kingdom, Phylum, Class))%>%
  spread(OFGO, zscore)

write_csv(hc_microbe%>%
            select(sample, everything()), "analysis/microbe_hc_df.csv")

# VIZUALIZATIONS -- PCOA --------------------------------------------------
pcoaMicrobe <- mic_organism_post_hoc%>%
  # inner_join(dunnett_microbe_pvals%>%
  #              select(DayNight, OTU, Organism), by = c('DayNight', 'OTU', 'Organism'))%>%
  filter(OTU %in% different_microbes)%>%
  filter(Organism != 'Water control')%>%
  select(c(Organism, DayNight, Replicate, OTU, asin))%>%
  unite(sample, c("Organism", "DayNight", "Replicate"), sep = "_")%>%
  left_join(microbe_taxonomy, by = "OTU")%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_id"), sep = ";")%>%
  unite(OFGO, c("Order", "Family", "Genus", "OTU"), sep = ";")%>%
  select(-c(Size, OTU_id, Kingdom, Phylum, Class))%>%
  spread(OFGO, asin)%>%
  column_to_rownames(var = "sample")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()

pdf('plots/pcoaMicrobes.pdf', width = 15, height = 12)
pcoaMicrobe$vectors%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample")%>%
  separate(sample, c('Organism', 'DayNight', 'Replicate'), sep = '_')%>%
  mutate(shape = case_when(DayNight == 'Day' ~ emoji('sunny'),
                           TRUE ~ emoji('crescent_moon')))%>%
  ggplot(., aes(x = Axis.1, y = Axis.2, color = Organism)) +
  geom_text(aes(label = shape), cex = 13, family= 'OpenSansEmoji') +
  scale_color_manual(values = org_colors_no_water) +
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

