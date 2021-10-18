## Script Written by Zach Quinlan 06/19/19
# Re-organization of DORCIERR_FCM_fDOM.R because it needs to be cleaner 07/15/2019
# Only working on daytime exudation and remineralization
# Rewritten with changes to RR3 starting pipeline 10/11/2019
# 16s rRNA amplicon seque3nce data added in a upaded 10/18/2019
# Added in ClassyFire annotations from Inchi_keys 10/07/2019
# This has been rewritten for the new pipeline 12/18/2019
# This was split to just look at metabolite trends on 05/12/2020
# This was again split to focus on the figures for the paper and multiple linear regressions on 12/07/2020


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


# CORES -- setting processors available -----------------------------------
doParallel::registerDoParallel(cores = 12)

# LOADING -- dataframes  ------------------------------------------------------
## FCM and fDOM data
dorc_fcm_fdom <- read_xlsx("~/Documents/GitHub/DORCIERR/data/raw/DOC_fDOM_FCM/DORCIERR_fDOM_FCM.xlsx")%>%
  rename(sample_name =`Sample Name of DORCIERR_FCM_Final`)%>%
  rename('DayNight' = 'Experiment')%>%
  rename(Organism = 'Organsim')%>%
  mutate(Organism = case_when(Organism == "Water" ~ "Water control",
                              TRUE ~ as.character(Organism)))

#DOC data
moorea_doc <- read_xlsx("~/Documents/GitHub/DORCIERR/data/raw/DOC_fDOM_FCM/MO17_ExpSummary_DOC.2018.04.04.xlsx")%>%
  dplyr::select(1:2)%>%
  rename(sample_name = 1)%>%
  rename(DOC = 2)

# True hits and analog hits are exported CSVs from GNPS
# True hits are more strictly matched to the library
true_hits <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/metabolomics/Library-hits.tsv")%>%
  rename("feature_number" = '#Scan#')

analog_hits <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/metabolomics/Analog-hits.tsv")%>%
  rename("feature_number" = '#Scan#')

# Node info includes networking information about each feature
node_info <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/metabolomics/Node_info.tsv")%>%
  rename('feature_number' = 'cluster index',
         'network' = 'componentindex')

csi_finger_id <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/metabolomics/summary_csi_fingerid.tsv")%>%
  rename(feature_number = experimentName)

# Canopus tries to classify each feature
canopus_anotations <- read_csv("~/Documents/GitHub/DORCIERR/data/raw/metabolomics/Canopus_classes.csv")

chemont_anotations <- read_csv("~/Documents/GitHub/DORCIERR/data/raw/metabolomics/categories.canopus.strings.nelsonMarch2019.CSV")%>%
  rename('canopus_annotation' = 'name')

canopusSirius4 <- read_csv('~/Documents/SDSU_Scripps/EcoNet/investigations/canopus_summary.csv')

# Sirius and Zodiac both try to assign molecular formulas to all the features
sirius_zodiac_anotations <- read_csv("~/Documents/GitHub/DORCIERR/data/raw/metabolomics/SIRIUS_Zodiac_converted.csv")%>%
  rename(feature_number = 1)%>%
  dplyr::select(-c(14:ncol(.)))


# Feature table has all features found within the experiments and blanks
# The columns need to be changed to the actual experiment sample codes
# Feature_table_raw is the raw export from MZMine
feature_table_raw <- read_csv("~/Documents/GitHub/DORCIERR/data/raw/metabolomics/Morrea_Feayures-Table_all_Gap-Filled5.csv")%>%
  rename('feature_number' = 'row ID')

ms_sample_codes <- read_csv("~/Documents/GitHub/DORCIERR/data/raw/metabolomics/Mo'orea 2017 Mass spec sample codes - Sheet1.csv")%>%
  rename('run_code' = 'Sample ID',
         'sample_code' = 'Sample Name')

# 16s rRNA sequences
microbe_abundance_raw <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/microbes/MCR2017.16S.Nelson.Pipeline.October2019/abundance_table_100.shared.tsv")
microbe_taxonomy <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/microbes/MCR2017.16S.Nelson.Pipeline.October2019/annotations_100.taxonomy.tsv")

# NAP
nap_df <- read_tsv("~/Documents/GitHub/DORCIERR/data/raw/metabolomics/moorea2017_NAP.tsv")%>%
  rename("feature_number" = "cluster.index")

inchikey_df <- read_csv("~/Documents/SDSU_Scripps/Moorea_2017/csi_inchikey.csv")%>%
  mutate(feature_number = as.character(feature_number))

# EcoNet
<<<<<<< HEAD
ecoNet <- read_csv('~/Documents/Github/DORCIERR/data/raw/metabolomics/ecoNetConsensus.csv')%>%
  select(-1)
=======
ecoNet <- read_csv('~/Documents/Github/EcoNet/src/EcoNetMoorea_dereplicated.csv')
>>>>>>> 866b30c22f5339816f0bda22d3e41936737b507f

# Linda annotations
# molnet_class <- read_csv("raw/metabolomics/Dorcierr_POC_TopDepletolites_Classified.csv")%>%
#   select(FinalClass, feature_number, ClassString)%>%
#   mutate(Organism = "Pocillopora verrucosa")%>%
#   bind_rows(read_csv("raw/metabolomics/Dorcierr_DIC_TopDepletolitesClass.csv")%>%
#               rename(FinalClass = FinalClassify)%>%
#               select(FinalClass, feature_number, ClassString)%>%
#               mutate(Organism = "Dictyota"),
#             read_csv("raw/metabolomics/Dorcierr_CCA_TopDepletolitesClassified.csv")%>%
#               select(FinalClass, feature_number, ClassString)%>%
#               mutate(Organism = "CCA"))%>%
#   group_by(feature_number)%>%
#   filter(row_number(feature_number) == 1)%>%
#   ungroup()%>%
#   mutate(feature_number = as.character(feature_number))%>%
#   distinct(FinalClass, feature_number, ClassString, Organism)

# MolNetEnhancer
molnet_class <- read_csv("~/Documents/GitHub/DORCIERR/data/analysis/Moorea2017_MolNetEnhancer.csv")%>%
  rename(feature_number = "cluster index")%>%
  mutate(feature_number = as.character(feature_number))%>%
  select(feature_number, CF_kingdom, CF_class, CF_subclass, CF_superclass)%>%
  unite(molnet_string, c(CF_kingdom, CF_superclass, CF_class, CF_subclass), sep = ";")

#PPl extraction efficiency
extraction_efficiency <- read_csv("~/Documents/GitHub/DORCIERR/data/raw/metabolomics/Extraction_efficiency.csv")%>%
  select(-c(X6:X8))

# compare <- deplete_classified%>%
#   inner_join(molnet_class, by = "feature_number")%>%
#   left_join(true_hits%>%
#               select(feature_number, Compound_Name)%>%
#               mutate(feature_number = as.character(feature_number)), by = "feature_number")%>%
#   left_join(node_info%>%
#               select(feature_number, network)%>%
#               mutate(feature_number = as.character(feature_number)), by = "feature_number")
#   
# 
# write_csv(compare, "~/Downloads/dorcierr_annotation_comparison.csv")

# CLEANING -- SIRIUS_Zodiac elemental composition of molecular formulas -------------------------------------------
networking_elements <- sirius_zodiac_anotations%>%
  filter(!ZodiacMF == "not_explainable")%>%
  group_by(feature_number)%>% 
  do(., rownames_to_column(as.data.frame(makeup(.$ZodiacMF, multiplier = 1), var = "element")))%>%
  spread(rowname, 3)

networking_elements[is.na(networking_elements)] <- 0

#Filter NOSC and Zodiac
networking_energy <- networking_elements%>%
  dplyr::select(c(1, 'C', 'H', 'N', 'O', 'P', 'S'))%>%
  add_column(NOSC = (-((4*.$C + .$H - 3*.$N - 2*.$O + 5*.$P - 2*.$S)/.$C)+4))%>%
  add_column(dG = 60.3-28.5*.$NOSC)%>%
  filter(NOSC < 4 & NOSC > -4)%>%
  filter(P <= 2)

# CLEANING -- Canopus---------------------------

canopus_annotation_names <- canopus_anotations%>%
  gather(canopus_annotation, canopus_probability, 2:ncol(.))

canopus_chemonnt_tidy <- left_join(canopus_annotation_names, chemont_anotations, by = "canopus_annotation")


# CLEANING -- FCM ---------------------------------------------------------
fcm_wdf <- dorc_fcm_fdom%>%
  dplyr::select(c(1:8,34:ncol(.)))

fcm_T0_5 <- dorc_fcm_fdom%>%
  select(-c(1:4, 9:36, 38, 39))%>%
  # filter(Timepoint == 'T0' | Timepoint == 'T4')%>%
  spread(Timepoint, `Cells µL-1`)%>%
  mutate(hour = case_when(Organism == 'Turf' ~ 23,
                          Organism == 'Dictyota'  ~ 23,
                          Organism == 'CCA'  ~ 23,
                          Organism == 'Porites lobata' ~ 48,
                          Organism == 'Pocillopora verrucosa' ~ 37,
                          Organism == 'Water control' ~ 37),
         final_cells = case_when(Organism == 'Turf' ~ T4,
                                 Organism == 'Dictyota'  ~ T4,
                                 Organism == 'CCA'  ~ T4,
                                 Organism == 'Porites lobata' ~ TF,
                                 Organism == 'Pocillopora verrucosa' ~ T6,
                                 Organism == 'Water control' ~ T6))%>%
  group_by(Organism, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         cells_ul = (log(final_cells) - log(T0))/(hour))%>%
  select(-c(T0:TF, hour))

# SET -- CANOPUS filters --------------------------------------------------
# Canopus annotations which are above level 3 AND 80% probability
canopus_filtered_tidy <- canopus_chemonnt_tidy%>%
  rename('feature_number' = 'name')%>%
  group_by(feature_number)%>%
  filter(canopus_probability >= 0.80,
         level > 3,
         canopus_probability == max(canopus_probability),
         level == max(level),
         nchar(CLASS_STRING) == max(nchar(CLASS_STRING)),
         nchar(canopus_annotation) == max(nchar(canopus_annotation)))%>%
  ungroup()

write_csv(canopus_filtered_tidy, "~/Documents/GitHub/DORCIERR/data/analysis/canopus_filtered_tidy.csv")

# Combines canopus, sirus, and zodiac
super_computer_annotations <- full_join(full_join(canopus_filtered_tidy, 
                                                  sirius_zodiac_anotations, by = "feature_number"),
                                        networking_energy, by = "feature_number")%>%
  add_column(`characterization scores` = .$`quality`)%>%
  mutate(`characterization scores` = case_when(`characterization scores` != "Good" ~ "Bad quality",
                                               ZodiacScore < .98 ~ "Low probability",
                                               TRUE ~ as.character(`characterization scores`)))



# CLEANING -- METADATA and filter out bad samples --------------------------
## join library hits, analog hits and super computer predictions
metadata <- full_join(node_info, 
                      full_join(super_computer_annotations,
                                full_join(feature_table_raw[1:4],
                                          full_join(true_hits, analog_hits, by = "feature_number",
                                                    suffix = c("Library_", "Analog_")),
                                          by = "feature_number"),
                                by = "feature_number"),
                      by = "feature_number")%>%
  left_join(molnet_class%>%
              mutate(feature_number = as.numeric(feature_number)), by = "feature_number")%>%
  # left_join(nap_df, by = "feature_number", suffix = c("", "_nap"))%>%
  left_join(csi_finger_id, by = "feature_number", suffix = c("", "_csi"))%>%
  add_column(binary_ID = .$LibraryID, .before = 1)%>%
  mutate(binary_ID = case_when(binary_ID != "N/A" ~ "1",
                               Compound_NameAnalog_ != "NA" ~ "2",
                               !is.na(molnet_string) ~ "3",
                               TRUE ~ as.character(binary_ID)))%>%
  add_column(combined_ID = .$LibraryID, .before = 1)%>%
  mutate(combined_ID = case_when(binary_ID == "1" ~ LibraryID,
                                 binary_ID == "3" ~ molnet_string,
                                 binary_ID == "2" ~ Compound_NameAnalog_,
                                 binary_ID == "N/A" ~ canopus_annotation,
                                 TRUE ~ as.character(binary_ID)))%>%
  mutate(inchi_binary = case_when(!is.na(INCHILibrary_) ~ "1",
                                  !is.na(INCHIAnalog_) & is.na(INCHILibrary_) | 
                                    !is.na(INCHIAnalog_) & INCHILibrary_ == "N/A" ~ "2",
                                  TRUE ~ "3"))%>%
  mutate(inchi_combined = case_when(inchi_binary == "1" ~ INCHILibrary_,
                                    inchi_binary == "2" ~ INCHIAnalog_,
                                    TRUE ~ "N/A"))


metadata$feature_number <- as.character(metadata$feature_number)

networking <- metadata%>%
  select(c(feature_number, network,combined_ID, binary_ID,  SmilesLibrary_, SmilesAnalog_,
           canopus_annotation:CLASS_STRING, ZodiacMF, `characterization scores`,
           C:dG, inchi_binary, inchi_combined))%>%
  left_join(molnet_class, by = 'feature_number')%>%
  separate(molnet_string, c("CF_kingdom", "CF_superclass", "CF_class", "CF_subclass"), sep = ";")%>%
  mutate(c_temp = case_when(C > 0 ~ "C",
                            TRUE ~ "_"),
         o_temp = case_when(O > 0 ~ "O",
                            TRUE ~ "_"),
         h_temp = case_when(H > 0 ~ "H",
                            TRUE~ "_"),
         n_temp = case_when(N > 0 ~ "N",
                            TRUE ~ "_"),
         p_temp = case_when(P > 0 ~ "P",
                            TRUE ~ "_"),
         s_temp = case_when(S > 0 ~ "S",
                            TRUE ~ "_"))%>%
  unite(simplified_makeup, c("c_temp", "h_temp", "o_temp", "n_temp", "p_temp", "s_temp"), sep = "")%>%
  mutate(simplified_makeup = gsub("_","", simplified_makeup),
         simplified_makeup = case_when(`characterization scores` != "Good" ~ "uncharacterized",
                                       simplified_makeup == "" ~ "uncharacterized",
                                       TRUE ~ as.character(simplified_makeup)))%>%
  left_join(inchikey_df, by = "feature_number", suffix = c("", "_inchi"))


## making feature table so we can remove blanks
# have to change the MS codes for sample codes
# dplyr::selecting for RR3 
feature_table_temp <- feature_table_raw%>%
  dplyr::select(-X326)%>%
  dplyr::select(-c(2:4))%>%
  gather(run_code, ion_charge, 2:ncol(.))%>%
  spread(feature_number, ion_charge)

feature_table_temp$run_code <- feature_table_temp$run_code%>%
  gsub(".mzXML Peak area", "", .)%>%
  gsub("_MSMS", "", .)

feature_table_dirty <- left_join(ms_sample_codes, feature_table_temp, by = "run_code")%>%
  dplyr::select(-run_code)%>%
  filter(sample_code %like any% c("%Blank%","R_%", "D_%", "M_%", "SPIFFy_%"))%>%
  filter(!sample_code %like any% c("%C18%", "%XAD%"))%>%
  gather(feature_number, ion_charge, 2:ncol(.))%>%
  spread(sample_code, ion_charge)

# DEFINING -- Samples and blanks ------------------------------------------
## defining what columns the samples are in
ions_samples <- 10:259

## defining different blanks
ions_blanks <- c(2:8, 260)


# FLAGGING -- BACKGROUND FEATURES flagging and removing ------------------------------------------------
# Background features are defined as features where max(blanks) >= 0.5*mean(samples)
ions <- feature_table_dirty

ions$max_blanks <- apply(ions[ions_blanks], 1, max)

ions$mean_samples <- apply(ions[ions_samples], 1, mean)


# This section finds the features where mean area under the curve across all samples is larger than 2 * max of the blanks
# The non background features are saved into a vector so that they can be filtered from the master database
# mean(samples)*0.5 > max_blanks
no_background <-ions%>%
  mutate(mean_samples = case_when(mean_samples*0.5 > max_blanks ~ "real",
                                  TRUE ~ "background"))%>%
  rename("background" = "mean_samples")

# FLAGGING -- TRANSIENT FEATURES flagging and removing ---------------------------------------------
# Transient features are defined as features who's area under the curve is not more than 2E5 in at least 3 samples
# This was determined by comparing gap filled and non-gap filled data and selecting for the lowest peak area in non-gap filled data
# This gives the assumption that anything lower than 2E5 is noise. See supplemental figure
feature_table_no_background_trans_finder <- feature_table_dirty%>%
  gather(sample, xic, 2:ncol(.))%>%
  separate(sample, "experiment", sep = "_", remove = FALSE, extra = "drop")%>%
  filter(!experiment %like% "%Blank%")%>%
  filter(!sample %like% "%Blank%")%>%
  mutate(experiment = case_when(experiment == "D" ~ "Dorcierr_transient",
                                experiment == "M" ~ "Mordor_transient",
                                experiment == "R" ~ "RR3_transient",
                                TRUE ~ "SPIFFy_transient"))%>%
  group_by(experiment)%>%
  nest()%>%
  mutate(data = map(data, ~ spread(.x, sample, xic)%>%
                      add_column(trans_feature_finder = rowSums(.[3:ncol(.)] > 2E5), .before = 2)%>%
                      mutate(transient = case_when(trans_feature_finder >= 3 ~ "real",
                                                   TRUE ~ "transient"))%>%
                      dplyr::select(c(feature_number, transient))))%>%
  unnest(data)%>%
  spread(experiment, transient)


feature_table_no_back_trans_filter <- full_join(feature_table_no_background_trans_finder, no_background, by = "feature_number")%>%
  dplyr::select(feature_number, background, everything())

# FILTERING -- out background and transient features ----------------------
dorcierr_real_features <- as.vector(feature_table_no_back_trans_filter%>%
                                      filter(background == "real")%>%
                                      filter(Dorcierr_transient == "real"))$feature_number

feature_table_no_back_trans <- feature_table_dirty%>%
  gather(sample, val, 2:ncol(.))%>%
  spread(feature_number, val)%>%
  filter(sample %like any% c("D_%", "%Blank%", "%OF%"))%>%
  dplyr::select(c(1, dorcierr_real_features))%>%
  gather(feature_number, val, 2:ncol(.))%>%
  spread(sample, val)


# FILTERING -- LOG2 bottleneck --------------------------------------------
## Everything has to double from T0 to TF (1 > log2(TF/T0) < 1)
log2_features_clean <- function(x) {
  new <- x%>%
    gather(sample_name, xic, 2:ncol(.))%>%
    ungroup()%>%
    separate(sample_name, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
    filter(!Experiment %like% "%Blank%",
           !Organism %like% "%Blank")%>%
    mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                  Experiment == "M" ~ "mordor",
                                  Experiment == "R" ~ "RR3",
                                  TRUE ~ as.character(Experiment)),
           Organism = case_when(Organism == "CC" ~ "CCA",
                                Organism == "DT" ~ "Dictyota",
                                Organism == "PL" ~ "Porites lobata",
                                Organism == "PV" ~ "Pocillopora verrucosa",
                                Organism == "TR" ~ "Turf",
                                Organism == "WA" ~ "Water control",
                                TRUE ~ as.character(Organism)))%>%
    separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
    mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                                TRUE ~ "Night"))%>%
    spread(Timepoint, xic)%>%
    group_by(Organism, DayNight, feature_number)%>%
    summarize_if(is.numeric, mean, na.rm = TRUE)%>%
    ungroup()%>%
    mutate(log2_change = log2(TF/T0),
           log2_change = case_when(TF > 0 & T0 == 0 ~ 6.6,
                                   TF == 0 & T0 > 0 ~ -6.6,
                                   TF == 0 & T0 == 0 ~ 0,
                                   TRUE ~ as.numeric(log2_change)))%>%
    filter(log2_change >= 1 | log2_change <= -1)
}

log2_features <- feature_table_no_back_trans%>%
  log2_features_clean()%>%
  select(-c(Organism, T0, TF))%>%
  group_by(DayNight, feature_number)%>%
  summarize_if(is.numeric, mean)


major_deplete_features <- feature_table_no_back_trans%>%
  log2_features_clean()%>%
  group_by(DayNight, feature_number, Organism)%>%
  filter(min(log2_change) <= -3.3)%>%
  select(-c(T0, TF, log2_change))%>%
  ungroup()%>%
  group_by(DayNight, feature_number)%>%
  summarize_if(is.numeric, mean)

log2_change_vals <- feature_table_no_back_trans%>%
  gather(sample_name, xic, 2:ncol(.))%>%
  ungroup()%>%
  separate(sample_name, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
  filter(!Experiment %like% "%Blank%",
         !Organism %like% "%Blank")%>%
  mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                Experiment == "M" ~ "mordor",
                                Experiment == "R" ~ "RR3",
                                TRUE ~ as.character(Experiment)),
         Organism = case_when(Organism == "CC" ~ "CCA",
                              Organism == "DT" ~ "Dictyota",
                              Organism == "PL" ~ "Porites lobata",
                              Organism == "PV" ~ "Pocillopora verrucosa",
                              Organism == "TR" ~ "Turf",
                              Organism == "WA" ~ "Water control",
                              TRUE ~ as.character(Organism)))%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))%>%
  spread(Timepoint, xic)%>%
  group_by(Organism, DayNight, feature_number)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         log2_change = log2(TF/T0))%>%
  ungroup()%>%
  mutate(log2_change = case_when(TF > 0 & T0 == 0 ~ 6.6,
                                 TF == 0 & T0 > 0 ~ -6.6,
                                 TF == 0 & T0 == 0 ~ 0,
                                 TRUE ~ as.numeric(log2_change)))



# FILTERING -- Log2 Exometabolites ----------------------------------------
exometabolite_features <- feature_table_no_back_trans%>%
  gather(sample_name, xic, 2:ncol(.))%>%
  ungroup()%>%
  separate(sample_name, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
  filter(!Experiment %like% "%Blank%",
         !Organism %like% "%Blank")%>%
  mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                Experiment == "M" ~ "mordor",
                                Experiment == "R" ~ "RR3",
                                TRUE ~ as.character(Experiment)),
         Organism = case_when(Organism == "CC" ~ "CCA",
                              Organism == "DT" ~ "Dictyota",
                              Organism == "PL" ~ "Porites lobata",
                              Organism == "PV" ~ "Pocillopora verrucosa",
                              Organism == "TR" ~ "Turf",
                              Organism == "WA" ~ "Water control",
                              TRUE ~ as.character(Organism)))%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))%>%
  group_by(DayNight, Timepoint, feature_number, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  spread(Organism, xic)%>%
  gather(Organism, xic, `CCA`:`Turf`)%>%
  mutate(log_org = log2(xic/`Water control`))%>%
  ungroup()%>%
  filter(log_org > 3.3)

day_exometabolites <- exometabolite_features%>%
  filter(DayNight == 'Day')%>%
  select(feature_number)%>%
  unique()

org_exometabolites <- exometabolite_features%>%
  select(feature_number, Organism, DayNight)%>%
  unique()

overlapping_exometabolites <- org_exometabolites%>% 
  mutate(num_organisms = 1)%>% 
  ungroup()%>% 
  select(-Organism)%>%
  group_by(feature_number)%>%
  summarise_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(num_organisms)%>%
  mutate(num_features = 1)%>%
  summarize_if(is.numeric, sum)

unique_benthic_metabolites <- org_exometabolites%>%
  mutate(num_organisms = 1)%>% 
  ungroup()%>% 
  group_by(feature_number)%>%
  mutate(num_organisms = sum(num_organisms))%>%
  filter(num_organisms == 1)%>%
  select(-num_organisms)


benthic_produced_exometabolites <- exometabolite_features%>%
  filter(DayNight == 'Day',
         Timepoint == 'T0')%>%
  select(feature_number, Organism)%>%
  mutate(bin = 'yes',
         dorc_prd = "produced")%>%
  group_by(feature_number)%>%
  add_tally(name = 'num_organism')%>%
  unite(Organism, c('Organism', 'dorc_prd'), sep = '_')%>%
  spread(Organism, bin)%>%
  mutate(dorcierr = 'yes')

benthic_produced_exometabolites[is.na(benthic_produced_exometabolites)] <- 'no'

write_csv(benthic_produced_exometabolites, './analysis/dorcierr_day_benthic_exometabolite_features.csv')


# FILTERING -- raw XIC minimum --------------------------------------------
min_filter_pre <- feature_table_no_back_trans%>%
  gather(sample_name, xic, 2:ncol(.))%>%
  ungroup()%>%
  separate(sample_name, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
  filter(!Experiment %like% "%Blank%",
         !Organism %like% "%Blank")%>%
  mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                Experiment == "M" ~ "mordor",
                                Experiment == "R" ~ "RR3",
                                TRUE ~ as.character(Experiment)),
         Organism = case_when(Organism == "CC" ~ "CCA",
                              Organism == "DT" ~ "Dictyota",
                              Organism == "PL" ~ "Porites lobata",
                              Organism == "PV" ~ "Pocillopora verrucosa",
                              Organism == "TR" ~ "Turf",
                              Organism == "WA" ~ "Water control",
                              TRUE ~ as.character(Organism)))%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))

#Building filters and dataframes for plotting
min_filter <- min_filter_pre%>%
  group_by(feature_number, Organism, DayNight, Timepoint)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  group_by(feature_number, Organism, DayNight)%>%
  filter(max(xic) >= 1*10^6)%>%
  left_join(networking%>%
              select(network, feature_number), by = 'feature_number')

# FILTERING -- Effect of all three filters --------------------------------
min_filter_feature_tag <- min_filter%>%
  ungroup()%>%
  filter(DayNight == 'Day',
         network %in% c('131','21','7','55','107','198','466', '165', '619', '627', '756', '924', '1314','80','141', '249','65','336','355','346'))%>%
  select(feature_number, Organism)%>%
  mutate(xic_min_filter = 1)

no_min_filter <- min_filter_pre%>%
  left_join(networking%>%
              select(network, feature_number), by = 'feature_number')

log2_filter <- min_filter_pre%>%
  left_join(networking%>%
              select(network, feature_number), by = 'feature_number')%>%
  inner_join(log2_features%>%
               select(DayNight, feature_number), by = c('DayNight', 'feature_number'))%>%
  inner_join(min_filter%>%
               select(feature_number, DayNight, Organism), by = c('DayNight', 'Organism', 'feature_number'))

all_filter <- min_filter_pre%>%
  left_join(networking%>%
              select(network, feature_number), by = 'feature_number')%>%
  inner_join(log2_features%>%
               select(DayNight, feature_number), by = c('DayNight', 'feature_number'))%>%
  inner_join(min_filter%>%
               select(feature_number, DayNight, Organism), by = c('DayNight', 'Organism', 'feature_number'))%>%
  inner_join(org_exometabolites, by = c('DayNight', 'Organism', 'feature_number'))

#COunting number of features in each network
num_features_no_filter <- no_min_filter%>%
  select(feature_number, network)%>%
  unique()%>%
  group_by(network)%>%
  mutate(num_features = 1)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()

num_features_min_filter <- min_filter%>%
  ungroup()%>%
  select(feature_number, network)%>%
  unique()%>%
  group_by(network)%>%
  mutate(num_features = 1)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()

num_features_log2_filter <- log2_filter%>%
  ungroup()%>%
  select(feature_number, network)%>%
  unique()%>%
  group_by(network)%>%
  mutate(num_features = 1)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()

num_features_all_filter <- all_filter%>%
  ungroup()%>%
  select(feature_number, network)%>%
  unique()%>%
  group_by(network)%>%
  mutate(num_features = 1)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()

num_features_incremental <- no_min_filter%>%
  ungroup()%>%
  filter(DayNight == 'Day')%>%
  nest(data = everything())%>%
  mutate(log2 = map(data, ~inner_join(.x, log2_features%>%
                                        select(DayNight, feature_number), by = c('DayNight', 'feature_number'))%>%
                      select(feature_number)%>%
                      unique()%>%
                      nrow()),
         exometabolite = map(data, ~inner_join(.x, log2_features%>%
                                                 select(DayNight, feature_number), by = c('DayNight', 'feature_number'))%>%
                               inner_join(min_filter%>%
                                            select(feature_number, DayNight, Organism), by = c('DayNight', 'Organism', 'feature_number'))%>%
                               select(feature_number, Organism)%>%
                               unique()%>%
                               group_by(Organism)%>%
                               mutate(count = 1)%>%
                               summarize_if(is.numeric, sum)),
         min = map(data, ~ inner_join(.x, log2_features%>%
                                        select(DayNight, feature_number), by = c('DayNight', 'feature_number'))%>%
                     inner_join(min_filter%>%
                                  select(feature_number, DayNight, Organism), by = c('DayNight', 'Organism', 'feature_number'))%>%
                     inner_join(org_exometabolites, by = c('DayNight', 'Organism', 'feature_number'))%>%
                     select(feature_number, Organism)%>%
                     unique()%>%
                     group_by(Organism)%>%
                     mutate(count = 1)%>%
                     summarize_if(is.numeric, sum))
         )%>%
  select(-data)



#Making the plots
plots_no_filter <- no_min_filter%>%
  group_by(feature_number, Timepoint, Organism)%>%
  filter(DayNight == 'Day',
         network %in% c('131','21','7','55','107','198','466', '165', '619', '627', '756', '924', '1314','80','141', '249','65','336','355','346'))%>%
  left_join(min_filter_feature_tag, by = c('feature_number', 'Organism'))%>%
  mutate(mean = mean(xic),
         sd = sd(xic),
         fill_color = case_when(is.na(xic_min_filter) ~ "#F98400",
                                TRUE ~ "#006658"),
         fill_name = case_when(is.na(xic_min_filter) ~ 'b_Noise',
                               TRUE ~ 'a_Real'))%>%
  select(-c(xic, Replicate, xic_min_filter))%>%
  unique()%>%
  ungroup()%>%
  left_join(num_features_no_filter, by = 'network')%>%
  arrange(network)%>%
  mutate(color_b = 'black',
         label = 'No filters')
# group_by(network, num_features, label)%>%
# nest()%>%
# mutate(plots = map(data, ~ ggplot(.x, aes(Timepoint, mean, fill = fill_bl, color = feature_number)) +
#                      geom_bar(stat = 'identity', size = 0.0001) +
#                      facet_wrap(~Organism) +
#                      scale_fill_manual(values = c("#F2AD00", "#006658")) +
#                      ylab("Mean feature XIC") +
#                      scale_color_manual(values = .x$color_b) +
#                      ggtitle(sprintf("Net:%4.0f (%3.0f features)", network, num_features)) +
#                      theme(legend.position = 'none',
#                            plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
#                            axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
#                            axis.text.y = element_text(size = 20),
#                            plot.title = element_text(size = 15, face = "bold"),
#                            panel.background = element_rect(fill = "transparent"), # bg of the panel
#                            plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#                            panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
#                            panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"))))%>%
# select(-data)%>%
# mutate(label = 'No filters')


plots_log2_filter <- log2_filter%>%
  group_by(feature_number, Timepoint, Organism)%>%
  filter(DayNight == 'Day',
         network %in% c('131','21','7','55','107','198','466', '165', '619', '627', '756', '924', '1314','80','141', '249','65','336','355','346'))%>%
  mutate(mean = mean(xic),
         sd = sd(xic))%>%
  select(-c(xic, Replicate))%>%
  unique()%>%
  ungroup()%>%
  left_join(num_features_log2_filter, by = 'network')%>%
  arrange(network)%>%
  mutate(color_b = 'black',
         fill_color = "#006658",
         fill_name = "black",
         label = 'Timepoint filters')
# group_by(network, num_features)%>%
# nest()%>%
# mutate(plots = map(data, ~ ggplot(.x, aes(Timepoint, mean, color = feature_number, fill = fill_bl)) +
#                      geom_bar(stat = 'identity', size= 0.0001) +
#                      facet_wrap(~Organism) +
#                      ylab("Mean feature XIC") +
#                      scale_color_manual(values = .x$color_b) +
#                      scale_fill_manual(values = .x$fill_bl) +
#                      ggtitle(sprintf("Net:%4.0f (%3.0f features)", network, num_features)) +
#                      theme(legend.position = 'none',
#                            plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
#                            axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
#                            axis.text.y = element_text(size = 20),
#                            plot.title = element_text(size = 15, face = "bold"),
#                            panel.background = element_rect(fill = "transparent"), # bg of the panel
#                            plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#                            panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
#                            panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"))))%>%
# select(-data)%>%
# mutate(label = 'Timepoint filter')

plots_all_filter <- all_filter%>%
  group_by(feature_number, Timepoint, Organism)%>%
  filter(DayNight == 'Day',
         network %in% c('131','21','7','55','107','198','466', '165', '619', '627', '756', '924', '1314','80','141','249','65','336','355','346'))%>%
  mutate(mean = mean(xic),
         sd = sd(xic))%>%
  select(-c(xic, Replicate))%>%
  unique()%>%
  ungroup()%>%
  left_join(num_features_all_filter, by = 'network')%>%
  arrange(network)%>%
  mutate(color_b = 'black',
         fill_color = "#006658",
         fill_name = "black",
         label = 'All filters')%>%
  bind_rows(plots_no_filter, plots_log2_filter)%>%
  group_by(network, Organism, Timepoint)%>%
  mutate(y_max = sum(mean))%>%
  ungroup()%>%
  group_by(network)%>%
  mutate(y_max = max(y_max))%>%
  ungroup()%>%
  group_by(network, num_features, label)%>%
  nest()%>%
  mutate(plots = map(data, ~ ggplot(.x, aes(Timepoint, mean, color = feature_number, fill = fill_name)) +
                       geom_bar(stat = 'identity', size= 0.0001) +
                       facet_wrap(~Organism) +
                       ylab("Mean feature XIC") +
                       scale_color_manual(values = .x$color_b) +
                       scale_fill_manual(values = c("#006658", "#F2AD00")) +
                       ggtitle(sprintf("Net:%4.0f (%3.0f features)", network, num_features)) +
                       ylim(0, min(.x$y_max)) +
                       theme(legend.position = 'none',
                             plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
                             axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
                             axis.text.y = element_text(size = 20),
                             plot.title = element_text(size = 15, face = "bold"),
                             panel.background = element_rect(fill = "transparent"), # bg of the panel
                             plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                             panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
                             panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"))))%>%
  select(-data)


# pdf("plots/networks_xic_notfiltered.pdf", width = 15, height = 10)
# plots_no_filter$plots
# dev.off()
# 
# pdf("plots/networks_xic_filtered.pdf", width = 15, height = 10)
# plots_min_filter$plots
# dev.off()


# DORCIERR feature_table --------------------------------------------------
feature_table_relnorm <- feature_table_no_back_trans%>%
  gather(sample_name, xic, 2:ncol(.))%>%
  ungroup()%>%
  separate(sample_name, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_", remove = FALSE)%>%
  filter(!Experiment %like% "%Blank%",
         !Organism %like% "%Blank")%>%
  mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                Experiment == "M" ~ "mordor",
                                Experiment == "R" ~ "RR3",
                                TRUE ~ as.character(Experiment)),
         Organism = case_when(Organism == "CC" ~ "CCA",
                              Organism == "DT" ~ "Dictyota",
                              Organism == "PL" ~ "Porites lobata",
                              Organism == "PV" ~ "Pocillopora verrucosa",
                              Organism == "TR" ~ "Turf",
                              Organism == "WA" ~ "Water control",
                              TRUE ~ as.character(Organism)))%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))%>%
  select(-c("Experiment", "Organism", "Replicate", "Timepoint", "DayNight"))


feature_table_combined <- right_join(metadata, feature_table_relnorm, by = "feature_number")

dorcierr_features_wdf  <- feature_table_combined%>%
  dplyr::select(c(feature_number, everything()))



# PRE-CLEANING -- Making Dorcierr working data frame for stats -----------------------------
# write_csv(dorcierr_features_wdf, "~/Documents/GitHub/DORCIERR/data/analysis/Dorcierr_feature_table_master_post_filtered.csv")

blanks_wdf <- dorcierr_features_wdf%>%
  filter(sample_name %like% "%Blank%",
         !sample_name == 'D_Blank_DI')%>%
  mutate(sample_name = case_when(sample_name == "Blank_Lot_6350565_01"~ "Blank_Blank_635056501",
                                 sample_name == "Blank_SD_01_A" ~ "Blank_Blank_SD01A",
                                 sample_name == "Blank_SD_01_B" ~ "Blank_Blank_SD01B",
                                 sample_name == "Blank_SD_LoRDI" ~ "Blank_Blank_SDLoRDI",
                                 sample_name == "Blank_SD_PPL" ~ "Blank_Blank_SDPPL",
                                 sample_name == "Blank? look up name on PPL" ~ "Blank_Blank_unknown",
                                 sample_name == "D_Blank" ~ "Blank_Blank_D",
                                 TRUE ~ "Blank_Blank_Spiffy"))%>%
  separate(sample_name, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")

dorc_wdf <- dorcierr_features_wdf%>%
  separate(sample_name, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
  filter(!Experiment %like% "%Blank%",
         !Organism %like% "%Blank")%>%
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
                                     TRUE ~ as.character(Organism)))%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))

# PRE-STATS CLEANING -- DOM-stats - RELATIVIZATION AND NORMALIZATION -- xic_log10 -----------------------------------------
dom_stats_wdf<- dorc_wdf%>%
  filter(!Organism == "Influent",
         !Organism == "Offshore",
         Timepoint == "T0" | Timepoint == "TF")%>%
  select(feature_number, Organism:ncol(.))%>%
  inner_join(log2_features%>%
               select(feature_number, DayNight), by = c('DayNight', 'feature_number'))%>%
  inner_join(min_filter%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('Organism', 'DayNight', 'feature_number'))

feature_stats_wdf <- dom_stats_wdf%>%
  inner_join(org_exometabolites, by = c('feature_number', 'Organism', 'DayNight'))%>%
  group_by(Organism, Replicate, Timepoint, DayNight)%>%
  mutate(log10 = log10(xic + 1),
         ra = xic/sum(xic),
         asin = asin(ra))

# SET SEED ----------------------------------------------------------------
set.seed(2005)


# STATS  -- T-TEST Network level ------------------------------------------
net_test <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = "feature_number")%>%
  filter(network != '-1')%>%
  group_by(network, DayNight)%>%
  nest()%>%
  mutate(greater = map(data, ~ t.test(log10 ~ Timepoint, .x, alternative = "greater")),
         lesser = map(data, ~ t.test(log10 ~ Timepoint, .x, alternative = "less")))%>%
  select(-data)%>%
  ungroup()%>%
  mutate(net_act = network)

singlenode_test <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = "feature_number")%>%
  filter(network == '-1')%>%
  group_by(feature_number, DayNight)%>%
  nest()%>%
  mutate(greater = map(data, ~ t.test(log10 ~ Timepoint, .x, alternative = "greater")),
         lesser = map(data, ~ t.test(log10 ~ Timepoint, .x, alternative = "less")))%>%
  select(-data)%>%
  ungroup()%>%
  mutate(net_act = -as.numeric(feature_number))


all_activity <- net_test%>%
  select(-network)%>%
  bind_rows(singlenode_test%>%
              select(-feature_number))%>%
  mutate(greater = map(greater, ~ .x["p.value"][[1]]))%>%
  mutate(lesser = map(lesser, ~ .x["p.value"][[1]]))%>%
  ungroup()%>%
  mutate(greater = as.numeric(greater),
         lesser = as.numeric(lesser),
         FDR_greater = p.adjust(greater, method = "BH"),
         FDR_lesser = p.adjust(lesser, method = "BH"))%>%
  mutate(activity = case_when(FDR_greater < 0.05 ~ "depletolite",
                              FDR_lesser < 0.05 ~ "accumolite",
                              TRUE ~ "recalcitrant"))%>%
  select(net_act, DayNight, activity)


# STATS - multiple regressions --------------------------------------------
mul_reg <- log2_change_vals%>%
  inner_join(feature_stats_wdf%>%
               ungroup()%>%
               select(feature_number, Organism, DayNight)%>%
               unique(), 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network, N, P, C, O, H, NOSC),
            by = 'feature_number')%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  filter(DayNight == 'Day',
         T0 != 0)%>%
  mutate(n_presence = case_when(N > 0 ~ 1,
                                TRUE ~ 1))%>%
  mutate(nc = N/C,
         oc = O/C,
         hc = H/C,
         log_nc = log10(nc),
         log_oc = log10(oc),
         log_hc = log10(hc),
         logxic = log10(T0),
         log_snc = log10(nc*T0),
         log_soc = log10(oc*T0),
         log_shc = log10(hc*T0))%>%
  mutate(Replicate = as.character(Replicate))%>%
  left_join(feature_stats_wdf%>%
              ungroup()%>%
              filter(Timepoint == 'T0')%>%
              select(feature_number, Organism, Replicate, DayNight, ra),
            by = c('feature_number', 'Organism', 'Replicate', 'DayNight'))%>%
  group_by(feature_number, Organism, DayNight)%>%
  mutate(ra_mean = mean(ra, na.rm = TRUE))%>%
  ungroup()

pdf("./plots/correlation_verify.pdf", width = 15, height = 10)
corr_verify <- mul_reg%>%
  filter(nc > 0 & oc > 0,
         NOSC < 0)%>%
  group_by(feature_number, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  select(log2_change, `row m/z`, NOSC, log_snc, log_soc, log_shc, log_nc, log_oc, log_hc, logxic)%>%
  cor()%>%
  corrplot::corrplot()
dev.off()


## N model
n_mulreg <- mul_reg%>%
  filter(N > 0,
         O > 0,
         NOSC < 0)%>%
  group_by(feature_number, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()


lm_test_mulreg <- lm(log2_change ~ `row m/z`+ NOSC + log_nc + log_oc + log_hc + logxic, data = n_mulreg)

#Step AIC
sink('./analysis/aic_model_selection.txt')
n_model_select <- stepAIC(lm_test_mulreg, direction = 'forward')
sink()

#Dredge
options(na.action = na.fail)
dredge_n_select <- MuMIn::dredge(lm_test_mulreg)
sink('./analysis/dredge_model_selection.txt')
head(dredge_n_select)
sink()

#organism random effect and residuals
lm_orgrand <- lme4::lmer(log2_change ~ log_nc + log_oc + log_hc + logxic + NOSC + (1|Organism), data = n_mulreg)

resids <- n_mulreg%>%
  mutate(residuals = (lm(log2_change ~ log_nc + log_oc + log_hc + logxic + NOSC, data = n_mulreg))$residuals)

resids_model <- lm(residuals ~ log_nc + log_oc + log_hc + logxic + NOSC, data = resids)

sink('./analysis/org_random_residsmodels.txt')
summary(lm_orgrand)
summary(resids_model)
sink()

#grouped by organism
lm_n_mulreg <- n_mulreg%>%
  group_by(Organism)%>%
  nest()%>%
  mutate(data = map(data, ~lm(log2_change ~ log_nc + log_oc + log_hc + logxic + NOSC, data = .x)),
         adj.r2 = map(data, ~ summary(.x)[['adj.r.squared']]),
         model_p = map(data, ~ summary(.x)[["fstatistic"]]%>%
                         as.data.frame()%>%
                         mutate(p = pf(.[1], .[2], .[3], lower.tail = FALSE))),
         data = map(data, ~ summary(.x)[["coefficients"]]%>%
                      as.data.frame()%>%
                      rename(coefficient = 1,
                             pval = 4)%>%
                      select(1,4)%>%
                      rownames_to_column(var = 'variable')))%>%
  unnest(c(data, adj.r2))%>%
  mutate(coefficient = case_when(pval >= 0.05 ~ NA_real_,
                                 TRUE ~ as.numeric(coefficient)))

lm_n_coefficient_table <- lm_n_mulreg%>%  
  select(-c(pval, adj.r2, model_p))%>%
  spread(variable, coefficient)

lm_n_r2_table <-lm_n_mulreg%>%
  select(-c(pval, coefficient, variable, model_p))%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)

lm_n_model_p <- lm_n_mulreg%>%
  select(-c(pval, coefficient, variable, adj.r2))%>%
  unnest(model_p)%>%
  select(-2)%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)

n_mean_coefficients <- lm_n_coefficient_table%>%
  ungroup()%>%
  select(-Organism)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)
<<<<<<< HEAD

lm_n_combined_table <- lm_n_model_p%>%
  left_join(lm_n_r2_table, by = 'Organism')%>%
  left_join(lm_n_coefficient_table, by = 'Organism')

sink('./analysis/linear_model_dredged_variabled.txt')
lm_n_combined_table
sink()

=======

lm_n_combined_table <- lm_n_model_p%>%
  left_join(lm_n_r2_table, by = 'Organism')%>%
  left_join(lm_n_coefficient_table, by = 'Organism')

sink('./analysis/linear_model_dredged_variabled.txt')
lm_n_combined_table
sink()

>>>>>>> 866b30c22f5339816f0bda22d3e41936737b507f
## Absent model
absent_mulreg <- mul_reg%>%
  filter(N == 0,
         NOSC < 0)%>%
  group_by(feature_number, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()

lm_test_absent <- lm(log2_change ~ NOSC + logxic, data = absent_mulreg)

dredge_absent_select <- MuMIn::dredge(lm_test_absent)

head(dredge_absent_select)

#Coefficients
lm_absent_coefficients <- absent_mulreg%>%
  group_by(Organism)%>%
  nest()%>%
  mutate(n = map(data, ~ mutate(.x, n = 1,
                                n = sum(n))%>%
                   select(n)),
         data = map(data, ~ lm(log2_change ~ NOSC + logxic, data = .x)),
         adj.r2 = map(data, ~ summary(.x)[['adj.r.squared']]),
         model_p = map(data, ~ summary(.x)[["fstatistic"]]%>%
                         as.data.frame()%>%
                         mutate(p = pf(.[1], .[2], .[3], lower.tail = FALSE))),
         data = map(data, ~ summary(.x)[["coefficients"]]%>%
                      as.data.frame()%>%
                      rename(coefficient = 1,
                             pval = 4)%>%
                      select(1,4)%>%
                      rownames_to_column(var = 'variable')))%>%
  unnest(c(data, adj.r2))%>%
  mutate(coefficient = case_when(pval >= 0.05 ~ NA_real_,
                                 TRUE ~ as.numeric(coefficient)))

lm_absent_coefficient_table <- lm_absent_coefficients%>%  
  select(-c(n, pval, adj.r2, model_p))%>%
  spread(variable, coefficient)

lm_absent_r2_table <-lm_absent_coefficients%>%
  select(-c(pval, coefficient, variable, model_p))%>%
  unnest(n)%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)

lm_absent_model_p <- lm_absent_coefficients%>%
  select(-c(pval, coefficient, variable, adj.r2))%>%
  unnest(model_p)%>%
  select(-2)%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)

lm_absent_combined_table <- lm_absent_model_p%>%
  left_join(lm_absent_r2_table, by = 'Organism')%>%
  left_join(lm_absent_coefficient_table, by = 'Organism')

mean_absent_coefficients <- lm_absent_coefficient_table%>%
  ungroup()%>%
  select(-Organism)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)
<<<<<<< HEAD

#absent org rand
lm_orgrand <- lme4::lmer(log2_change ~ logxic + NOSC + (1|Organism), data = absent_mulreg)%>%
  tidy()

=======

#absent org rand
lm_orgrand <- lme4::lmer(log2_change ~ logxic + NOSC + (1|Organism), data = absent_mulreg)%>%
  tidy()

>>>>>>> 866b30c22f5339816f0bda22d3e41936737b507f

# sink("./analysis/multiple_reg_coefficients_whole_metabolome_aic_models.txt")
# summary(lm_n_mulreg)
# summary(o_mulreg)
# sink()

# VIZUALIZATIONS -- RGB hex codes for orgs --------------------------------
org_colors_no_water <- c("#A30029","#669900", "#FF850A", "#9900FF", "#33CC33")

# Supplemental -- Checking for gaussian distribution in NOSC, Mass, NC and PC --------
# LilliFors and other tests for normallity will not be effective due to large sample size. 
# QQ-plots used to establish how 'normal' the data appear
glm_df <- log2_change_vals%>%
  left_join(networking%>%
              select(feature_number, network, NOSC), by = "feature_number")%>%
  filter(network != "-1",
         NOSC <= 0)%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'recalcitrant',
         DayNight == 'Day')%>%
  mutate(activity2 = activity)%>% 
  select(-c(T0,TF))%>%
  gather(response_var, value, NOSC:`row m/z`)

log_reg <- glm_df%>%
  inner_join(min_filter%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  inner_join(org_exometabolites, 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  inner_join(log2_features%>%
               select(feature_number, DayNight), 
             by = c('DayNight', 'feature_number'))

nc_logreg <- log2_change_vals%>%
  inner_join(feature_stats_wdf%>%
               ungroup()%>%
               select(feature_number, Organism, DayNight)%>%
               unique(), 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network, C:dG)%>%
              filter(network != '-1'),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  group_by(Organism, DayNight, Replicate)%>%
  mutate(nc = N/C,
         pc = P/C,
         sample_nc = T0*nc,
         sample_pc = T0*pc,
         activity2 = activity)%>%
  ungroup()%>%
  mutate(log_snc = log(sample_nc),
         log_spc = log(sample_pc))

log2_check <- (log_reg%>%
                 filter(response_var == 'NOSC'))$log2_change

#N:C
nc_check <- (mul_reg%>%
               filter(nc > 0,
                      !is.na(nc)))$nc

lnc_check <- (mul_reg%>%
                filter(nc > 0,
                       !is.na(nc)))$log_nc
#O:C
oc_check <- (mul_reg%>%
               filter(oc > 0,
                      !is.na(oc)))$oc

loc_check <- (mul_reg%>%
                filter(oc > 0,
                       !is.na(oc)))$log_oc
#H:C
hc_check <- (mul_reg%>%
               filter(hc > 0,
                      !is.na(hc)))$hc

lhc_check <- (mul_reg%>%
                filter(hc > 0,
                       !is.na(hc)))$log_hc

#NOSC and m/z
nosc_check <- (log_reg%>%
                 filter(response_var == 'NOSC')%>%
                 mutate(NOSC = value))$NOSC

mass_check <- (log_reg%>%
                 filter(response_var == 'row m/z')%>%
                 mutate(`row m/z` = value))$`row m/z`

# OC and HC are normal NC is not
pdf("./plots/gaussian_distribution_non_transformed.pdf", width = 15, height = 10)
car::qqPlot(nc_check, 
            ylab = "N:C quantiles", xlab = "Normal quantiles",
            main = 'QQ-plot:N:C')
car::qqPlot(lnc_check, 
            ylab = bquote(Log[10]~(N:C ~quantiles)), xlab = "Normal quantiles",
            main = 'QQ-plot: Log10 transformed N:C')
car::qqPlot(oc_check, 
            ylab = "O:C quantiles", xlab = "Normal quantiles",
            main = 'QQ-plot: O:C')
car::qqPlot(loc_check,
            ylab = bquote(Log[10]~(O:C ~quantiles)), xlab = "Normal quantiles",
            main = 'QQ-plot: Log10 transformed O:C')
car::qqPlot(hc_check, 
            ylab = "H:C quantiles", xlab = "Normal quantiles",
            main = 'QQ-plot: H:C')
car::qqPlot(lhc_check,
            ylab = bquote(Log[10]~(H:C ~quantiles)), xlab = "Normal quantiles",
            main = 'QQ-plot: Log10 transformed H:C')

car::qqPlot(nosc_check, 
            ylab = 'Nominal Oxidation State of Carbon (NOSC) quantiles', xlab = "Normal quantiles",
            main = 'QQ-plot: NOSC')
car::qqPlot(mass_check, 
            ylab = "Mass/charge quantiles", xlab = "Normal quantiles",
            main = 'QQ-plot: Mass/Charge')
dev.off()

pdf("./plots/gaussian_distribution_log2change.pdf", width = 15, height = 10)
car::qqPlot(log2_check, 
            ylab = bquote(Log[2] ~change), xlab = "Normal quantiles",
            main = bquote(QQ-plot: ~Log[2] ~change))
car::qqPlot(log10(log2_check + 20),
            ylab = bquote(Log[2] ~change), xlab = "Normal quantiles",
            main = bquote(QQ-plot: ~Log[10](Log[2] ~change)))
dev.off()

# VIZUALIZATIONS -- Org Pie charts----------------------------------------------------------
org_log2_ra <- feature_stats_wdf%>%
  filter(Timepoint == "T0",
         DayNight == "Day")%>%
  left_join(log2_change_vals%>%
              select(-c(Replicate, T0, TF))%>%
              group_by(feature_number, Organism, DayNight)%>%
              summarize_if(is.numeric, mean),
            by = c("feature_number", "Organism", "DayNight"))%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))

org_pie <- org_log2_ra%>%
  group_by(Organism, DayNight, Replicate, Timepoint)%>%
  select(-c(log10:asin))%>%
  mutate(ra = xic/sum(xic),
         count = 1)%>%
  ungroup()%>%
  group_by(Organism, Replicate, activity)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(Organism, activity)%>%
  mutate(err = sd(ra))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  group_by(Organism)%>%
  mutate(end = 2 * pi * cumsum(ra)/sum(ra),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1),
         explosive = case_when(activity == 'depletolite' ~ 0.2,
                               TRUE ~ 0),
         texty = case_when(activity == 'depletolite' & Organism == 'Turf' ~ 0.8*cos(middle),
                           activity == 'depletolite' & Organism == 'Porites lobata' ~ 0.95*cos(middle),
                           activity == 'depletolite' & Organism == 'CCA' ~ 0.95*cos(middle),
                           activity == 'recalcitrant' & Organism == 'Pocillopora verrucosa' ~ 1.6*cos(middle),
                           activity == 'recalcitrant' & Organism == 'Dictyota' ~ 1.05*cos(middle),
                           activity == 'recalcitrant' & Organism == 'Turf' ~ 1.15*cos(middle),
                           activity == 'depletolite' & Organism == 'Dictyota' ~ 0.95*cos(middle),
                           activity == 'depletolite' ~ 1.2*cos(middle),
                           activity == 'accumolite' ~ 1.2*cos(middle),
                           TRUE ~ 1.2*cos(middle)),
         textx = case_when(activity == 'accumolite' & Organism == 'Pocillopora verrucosa' ~ -20*sin(middle),
                           activity == 'recalcitrant' & Organism == 'Pocillopora verrucosa' ~ -1*sin(middle),  # Without Error
                           activity == 'accumolite' & Organism == 'CCA' ~ -5*sin(middle),
                           activity == 'depletolite' & Organism == 'CCA' ~ 1.6*sin(middle),
                           # activity == 'accumolite' & Organism == 'Dictyota' ~ -200*sin(middle),   #with Error
                           activity == 'accumolite' & Organism == 'Dictyota' ~ -50*sin(middle),    #Without error
                           activity == 'depletolite' & Organism == 'Dictyota' ~ 2*sin(middle),
                           activity == 'accumolite' & Organism == 'Porites lobata' ~ sin(middle),
                           activity == 'depletolite' & Organism == 'Porites lobata' ~ 1.6*sin(middle),
                           activity == 'recalcitrant' & Organism == 'Porites lobata' ~ -1*sin(middle), # Without error
                           # activity == 'accumolite' & Organism == 'Turf' ~ -60*sin(middle),  #with error
                           activity == 'accumolite' & Organism == 'Turf' ~ -45*sin(middle),   #without error
                           activity == 'depletolite' & Organism == 'Turf' ~ 1.8*sin(middle),
                           activity == 'depletolite' ~ 1.4*sin(middle),
                           activity == 'recalcitrant' ~ -3*sin(middle),
                           TRUE ~ 1.05*sin(middle)),
         n = as.character(sum(count)))%>%
  ungroup()

# org_pie_vis <-org_pie%>%
#   group_by(Organism, n)%>%
#   nest()%>%
#   mutate(data = map(data, ~ggplot(.x) +
#                       geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
#                                        start = start, end = end, 
#                                        fill = activity, explode = explosive, linetype = NA)) +
#                       geom_text(aes(x = textx, y = texty, 
#                                     label =  paste0(round(ra, digits = 4)*100, "%", " ± ", round(err, digits = 5)*100, "%"),
#                                     color = activity,
#                                     hjust = hjust, vjust = vjust), size = 20, show.legend = FALSE) +
#                       coord_fixed() +
#                       scale_x_continuous(limits = c(-3, 4),  # Adjust so labels are not cut off
#                                          name = "", breaks = NULL, labels = NULL) +
#                       scale_y_continuous(limits = c(-1.2, 1.5),      # Adjust so labels are not cut off
#                                          name = "", breaks = NULL, labels = NULL) +
#                       theme_classic() +
#                       ggtitle(paste0(Organism, " (n = ", n, ")")) + 
#                       # facet_wrap(~ Organism) +
#                       scale_fill_manual(values = c('#78B7C5', '#EBCC2A', "#006658")) +
#                       scale_color_manual(values = c('#78B7C5', '#EBCC2A', "#006658")) + 
#                       labs(fill = 'Network activity', color = 'Network activity') +
#                       theme(axis.line = element_blank(),
#                             axis.text = element_text(size = 20),
#                             axis.ticks = element_blank(),
#                             legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#                             legend.box.background = element_blank(),
#                             panel.background = element_rect(fill = "transparent"), # bg of the panel
#                             plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#                             strip.text = element_text(size=20),
#                             legend.text = element_text(size = 20),
#                             legend.title = element_text(size = 20),
#                             legend.position = 'none',
#                             title = element_text(size = 40, hjust = 0.5),
#                             plot.title = element_text(hjust = 0.5, vjust = 3),
#                             strip.background = element_blank())))

org_pie_noerr <-org_pie%>%
  group_by(Organism, n)%>%
  nest()%>%
  mutate(data = map(data, ~ggplot(.x) +
                      geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                                       start = start, end = end, 
                                       fill = activity, explode = explosive, linetype = NA)) +
                      geom_text(aes(x = textx, y = texty, 
                                    label =  paste0(round(ra, digits = 4)*100, "%"),
                                    color = activity,
                                    hjust = hjust, vjust = vjust), size = 20, show.legend = FALSE) +
                      coord_fixed() +
                      scale_x_continuous(limits = c(-3, 4),  # Adjust so labels are not cut off
                                         name = "", breaks = NULL, labels = NULL) +
                      scale_y_continuous(limits = c(-1.2, 1.5),      # Adjust so labels are not cut off
                                         name = "", breaks = NULL, labels = NULL) +
                      theme_classic() +
                      ggtitle(paste0(Organism, " (n = ", n, ")")) + 
                      # facet_wrap(~ Organism) +
                      scale_fill_manual(values = c('#78B7C5', '#EBCC2A', "#006658")) +
                      scale_color_manual(values = c('#78B7C5', '#EBCC2A', "#006658")) + 
                      labs(fill = 'Network activity', color = 'Network activity') +
                      theme(axis.line = element_blank(),
                            axis.text = element_text(size = 20),
                            axis.ticks = element_blank(),
                            legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                            legend.box.background = element_blank(),
                            panel.background = element_rect(fill = "transparent"), # bg of the panel
                            plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                            strip.text = element_text(size=20),
                            legend.text = element_text(size = 20),
                            legend.title = element_text(size = 20),
                            legend.position = 'none',
                            title = element_text(size = 40, hjust = 0.5),
                            plot.title = element_text(hjust = 0.5, vjust = 3),
                            strip.background = element_blank())))

pdf("./plots/org_pie.pdf", width = 17, height = 12)
# org_pie_vis$data
org_pie_noerr$data
dev.off()



# VIZUALIZATIONS -- Organism comparisons log2_change ----------------------
org_t0tf <- feature_stats_wdf%>%
  filter(DayNight == "Day")%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  select(-c(log10, ra, asin))%>%
  spread(Timepoint,xic)%>%
  group_by(Organism, Replicate, activity)%>%
  mutate(count = 1)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  gather(Timepoint, xic, T0:TF)%>%
  group_by(Organism, Timepoint, activity)%>%
  mutate(err = sd(xic, na.rm = TRUE))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  filter(activity == 'depletolite')%>% 
  mutate(count = case_when(Timepoint == 'TF' ~ ' ',
                           TRUE ~ as.character(count)))

pdf("./plots/org_depletoliteebar_netactivity.pdf", width = 15, height = 12)
org_t0tf%>%
  ggplot(aes(Organism, xic)) +
  geom_bar(aes(fill = Timepoint), stat = 'identity', position = position_dodge2(width= 1)) +
  geom_errorbar(aes(min = xic - err, max = xic + err), position = position_dodge2(width= 1)) +
  scale_fill_manual(values = c('#EBCC2A', '#EBCC2A' 
                               # '#78B7C5', '#EBCC2A', "#006658"
  )) +
  # geom_text(aes(label = paste0(count), vjust = -1.5), size = 8) +
  labs(y = 'Depletolites intensity (xic)', fill = 'Timepoint: ') +
  ylim(0, 1.3e+10) +
  theme(
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    axis.text.x = element_text(size = 25, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_line('grey'), # get rid of major grid
    panel.grid.minor = element_line('grey'), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    legend.position = 'top',
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )
dev.off()

boxplot_df <- feature_stats_wdf%>%
  filter(Timepoint == "TF",
         DayNight == "Day")%>%
  left_join(log2_change_vals%>%
              select(-c(T0, TF))%>%
              select(feature_number, Organism, DayNight, Replicate, log2_change),
            by = c("feature_number", "Organism", "DayNight", "Replicate"))%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))


org_log2_anova <- boxplot_df%>%
  filter(activity == 'depletolite')%>%
  group_by(activity)%>%
  nest()%>%
  mutate(anova = map(data, ~ aov(log2_change ~ Organism, data = .x)%>%
                       tidy()),
         tukey = map(data, ~ aov(log2_change ~ Organism, data = .x)%>%
                       TukeyHSD(p.adjust.methods = 'BH')%>%
                       tidy()))%>%
  select(-data)%>%
  unnest(tukey)

sink('./analysis/tukey_log2.txt')
org_log2_anova$anova

org_log2_anova%>%
  select(-anova)
sink()


pdf("./plots/org_log2.pdf", width = 15, height = 10)
boxplot_df%>%
  filter(activity == 'depletolite')%>%
  group_by(Organism, Replicate)%>%
  summarize_if(is.numeric, sum)%>%
  mutate(err = sd(log2_change))%>%
  ungroup()%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  ggplot(aes(Organism, log2_change)) +
  # geom_boxplot(aes(color ='#EBCC2A'), lwd = 2) +
  geom_bar(stat = 'identity', fill = '#EBCC2A') +
  geom_errorbar(aes(min = log2_change - err, max = log2_change + err)) +
  labs(y = bquote(atop(Total ~microbial ~depletetion ~of, ~depletolites ~(log[2] ~change)))) +
  scale_fill_manual(values = '#EBCC2A') +
  theme(
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    axis.text.x = element_text(size = 25, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_line('grey'), # get rid of major grid
    panel.grid.minor = element_line('grey'), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    legend.position = 'top'
  )

boxplot_df%>%
  filter(activity == 'depletolite')%>%
  ggplot(aes(Organism, log2_change)) +
  geom_boxplot(aes(color ='#EBCC2A'), lwd = 2) +
  # geom_bar(stat = 'identity', fill = '#EBCC2A') +
  # geom_errorbar(aes(min = log2_change - err, max = log2_change + err)) +
  labs(y = bquote(atop(Average ~depletolite, ~log[2] ~change))) +
  scale_color_manual(values = '#EBCC2A') +
  theme(
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    axis.text.x = element_text(size = 25, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_line('grey'), # get rid of major grid
    panel.grid.minor = element_line('grey'), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    legend.position = 'top'
  )
dev.off()


<<<<<<< HEAD
# VIZUALIZATIONS -- activity by classifcation -----------------------------
lability_classes <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(DayNight == 'Day')%>%
  left_join(ecoNet%>%
            rename(feature_number = scan)%>%
            mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  select(-c(xic, log10, asin))%>%
  spread(Timepoint, ra)%>%
  mutate(val = case_when(activity == 'recalcitrant' ~ mean(c(T0, TF), na.rm = TRUE),
                         activity == 'depletolite' ~ T0,
                         activity == 'accumolite' ~ TF))%>%
  group_by(activity)%>%
  nest()%>%
  mutate(data = map(data, ~filter(.x, !is.na(val))%>%
                      # group_by(ecoNetConsensus, superclass_consensus, class_consensus, subclass_consensus, Replicate)%>%
                      group_by(superclass_consensus, Replicate)%>%
                      summarize_if(is.numeric, sum, na.rm = TRUE)%>%
                      ungroup()
                      # group_by(ecoNetConsensus, superclass_consensus, class_consensus, subclass_consensus)%>%
                      # mutate(total_production = mean(val, na.rm = TRUE),
                      #        std = sd(val, na.rm = TRUE))%>%
                      # select(-c(Replicate:net_act, numberOfNodes:val, ecoNetConsensusScore))%>%
                      # unique()
                    ))%>%
    unnest(data)

pdf('plots/superclassActivityCompare.pdf', width = 15, height = 10)
lability_classes%>%
  # filter(superclass_consensus != 'NA')%>%
  # unite(superclass_class, c(superclass_consensus, class_consensus), sep = '_')%>%
  # filter(ecoNetConsensus != 'NA')%>%
  # mutate(superclass_consensus = case_when(is.na(superclass_consensus) ~ 'No Consensus',
  #                                         TRUE ~ superclass_consensus))%>%
  ggplot(aes(superclass_consensus, val, color = activity)) +
  geom_point(stat = 'identity', size = 4) +
  coord_flip() +
  scale_color_manual(labels = c('Accumolite', 'Depletolite', 'Recalcitrant'), values = c('#78B7C5', '#EBCC2A', "#006658")) +
  # geom_bar(stat = 'identity') +
  # geom_errorbar(aes(ymax = total_production + std, ymin = total_production - std)) +
  # facet_wrap(~activity, labeller = labeller(activity = label_wrap_gen()), scales = 'free_y', nrow = 3) +
  labs(y = bquote(Sum ~Intensity ~(log[10] ~XIC)), x = 'Putative Superclass Annotation', color = 'Activity') +
  # scale_y_log10() +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        axis.text.y = element_text(size = 14))
dev.off()

lability_elements <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(DayNight == 'Day')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  left_join(networking%>%
              select(feature_number, C, H, N, O, P, S, NOSC), by = 'feature_number')%>%
  filter(Replicate == 1,
         Timepoint == 'T0',
         DayNight == 'Day')%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  mutate(nc = N/C,
         pc = P/C,
         oc = O/C)

lability_elements%>%
  ggplot(aes(superclass_consensus, pc, color = activity)) +
  geom_boxplot() +
  # geom_point(stat = 'identity', size = 4) +
  coord_flip() +
  scale_color_manual(labels = c('Accumolite', 'Depletolite', 'Recalcitrant'), values = c('#78B7C5', '#EBCC2A', "#006658")) +
  labs(x = 'Putative Superclass Annotation', color = 'Activity') +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        axis.text.y = element_text(size = 14))

# VIZUALIZATIONS -- Venn Diagram ------------------------------------------
venn <- boxplot_df%>%
  ungroup()%>%
  filter(activity == 'depletolite')%>%
  select(feature_number, Organism)%>%
  unique()

venn_list <- venn$feature_number%>%
  as.vector()

venn_split <- split(venn_list, venn$Organism)

venn.diagram(venn_split, 
             './plots/venn_depletolites.png',
             category.names = c('CCA', 'Dic', 'Poc', 'Por', 'T'),
             imagetype = 'png',  
             height = 1500, 
             width = 1800,
             resolution = 300,
             fill = org_colors_no_water)

venn_features <- calculate.overlap(venn_split)[c('a1', 'a2', 'a3', 'a4', 'a5', 'a13', 'a31', 'a10')]

names(venn_features) <- c('CCA', 'Dictyota', 'Pocillopora verrucosa', 'Porites lobata', 'Turf', 'coral', 'Primary producers', 'fleshy algae')

venn_df <- venn_features%>%
  unlist()%>%
  as.data.frame()%>%
  rename('feature_number' = 1)%>%
  rownames_to_column(var = 'Intersection')%>%
  mutate(Intersection = gsub('([0-9])', '', Intersection))%>%
  left_join(metadata%>% 
              select(feature_number, network), 
            by= 'feature_number')%>%
  left_join(ecoNet%>%
              select(-network)%>%
              mutate(scan = as.character(scan))%>%
              rename('feature_number' = 'scan'), by = 'feature_number')

# write_csv(venn_depletion, 'analysis/vennDepleteion.csv')

# VIZUALIZATIONS -- Venn class depletion ----------------------------------
venn_depletion <- venn_df%>%
  rename('Organism' = 'Intersection')%>%
  left_join(log2_change_vals%>%
              filter(DayNight == 'Day')%>%
              ungroup(), by = c('feature_number', 'Organism'))%>%
  separate(ecoNetConsensus, c('superclass', 'class', 'subclass'), sep = ";", remove = FALSE)
  

percent_reduction <- venn_depletion%>%
  group_by(ecoNetConsensus, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  group_by(ecoNetConsensus)%>%
  mutate(label = sum(TF, na.rm = TRUE)/sum(T0, na.rm = TRUE))%>%
  select(ecoNetConsensus, label)%>%
  unique()%>%
  mutate(label = round((1-label)*100, digits = 2))

pdf('plots/depletoliteClassificationsSuperclass.pdf', width = 15, height = 13)
venn_depletion%>%
  gather(variable, value, T0:ncol(.))%>%
  filter(!Organism %in% c('coral', 'Primary producers', 'fleshy algae'),
         !variable %in% c('percent_reduction', 'log2_change'),
         !is.na(value))%>%
  left_join(percent_reduction, by = 'ecoNetConsensus')%>%
  group_by(ecoNetConsensus)%>%
  mutate(yLabel = sum(value))%>%
  ungroup()%>%
  mutate(ecoNetConsensusLevel = factor(ecoNetConsensusLevel),
         ecoNetConsensusLevel = fct_relevel(ecoNetConsensusLevel, c('superclass', 'class','subclass')))%>%
  filter(ecoNetConsensus != 'NA',
         superclass != 'NA')%>%
  ggplot(aes(variable, value, fill = Organism)) +
  geom_bar(stat = 'summary', fun.y = 'sum', position = 'stack') +
  scale_fill_manual(values = org_colors_no_water) +
  facet_wrap(~superclass, labeller = labeller(superclass = label_wrap_gen()), scales = 'free_y') +
  geom_text(aes(y = yLabel*1.1, x = 1.5, label = paste0(label, '%')), size = 6, check_overlap = TRUE, inherit.aes = FALSE) +
  labs(x = 'Timepoint', y = 'Sum Intensity (XIC)') +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13))
dev.off()  

pdf('plots/depletoliteClassificationsOrganismFacet.pdf', width = 20, height = 15)
venn_depletion%>%
  gather(variable, value, T0:ncol(.))%>%
  filter(!Organism %in% c('coral', 'Primary producers', 'fleshy algae'),
         !variable %in% c('percent_reduction', 'log2_change'),
         !is.na(value))%>%
  left_join(percent_reduction, by = 'ecoNetConsensus')%>%
  filter(ecoNetConsensus != 'NA',
         variable != 'T0',
         Replicate != c(3,4))%>%
  group_by(ecoNetConsensus, ecoNetConsensusLevel, Replicate, Organism)%>%
  summarize_if(is.numeric, sum, na.rm  = TRUE)%>%
  ungroup()%>%
  group_by(ecoNetConsensus, Organism)%>%
  mutate(meanVal = mean(value),
         std = sd(value))%>%
  select(ecoNetConsensus, Organism, ecoNetConsensusLevel, meanVal, std)%>%
  unique()%>%
  ggplot(aes(ecoNetConsensus, meanVal, fill = Organism)) +
  # geom_boxplot() +
  geom_bar(stat = 'summary', fun.y = 'mean', position = 'stack') +
  geom_errorbar(aes(ymin = meanVal - std, ymax = meanVal + std)) +
  scale_fill_manual(values = org_colors_no_water) +
  facet_wrap(~Organism, labeller = labeller(Organism = label_wrap_gen()), nrow = 5) +
  labs(x = 'Putative Consensus Annotation', y = 'Sum Intensity (XIC)') +
  scale_y_log10() +
  # coord_flip() +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 7)) 
dev.off()  

boxplot_df%>%
  ungroup()%>%
  filter(activity == 'depletolite')%>%
  select(feature_number, Organism)%>%
  unique()%>%
  left_join(metadata%>% 
              select(feature_number, network), 
            by= 'feature_number')%>%
  left_join(ecoNet%>%
              select(-network)%>%
              mutate(scan = as.character(scan))%>%
              rename('feature_number' = 'scan'), by = 'feature_number')%>%
  left_join(log2_change_vals%>%
              filter(DayNight == 'Day')%>%
              ungroup(), by = c('feature_number', 'Organism'))%>%
  select(-c(TF, log2_change, network, numberOfNodes, log2_change, ecoNetConsensusScore))%>%
  group_by(ecoNetConsensus, ecoNetConsensusLevel, Replicate, Organism)%>%
  summarize_if(is.numeric, sum, na.rm  = TRUE)%>%
  ungroup()%>%
  group_by(ecoNetConsensus, Organism)%>%
  mutate(meanVal = mean(T0),
         std = sd(T0))%>%
  select(ecoNetConsensus, Organism, ecoNetConsensusLevel, meanVal, std)%>%
  unique()%>%
  ggplot(aes(ecoNetConsensus, meanVal, fill = Organism)) +
  # geom_boxplot() +
  geom_bar(stat = 'summary', fun.y = 'mean', position = 'stack') +
  geom_errorbar(aes(ymin = meanVal - std, ymax = meanVal + std)) +
  scale_fill_manual(values = org_colors_no_water) +
  facet_wrap(~Organism, labeller = labeller(Organism = label_wrap_gen()), nrow = 5) +
  labs(x = 'Putative Consensus Annotation', y = 'Sum Intensity (XIC)') +
  scale_y_log10() +
  # coord_flip() +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 7)) 


depletion_hc <- boxplot_df%>%
  ungroup()%>%
  filter(activity == 'depletolite')%>%
  select(feature_number, Organism)%>%
  unique()%>%
  left_join(metadata%>% 
              select(feature_number, network), 
            by= 'feature_number')%>%
  left_join(ecoNet%>%
              select(-network)%>%
              mutate(scan = as.character(scan))%>%
              rename('feature_number' = 'scan'), by = 'feature_number')%>%
  left_join(log2_change_vals%>%
              filter(DayNight == 'Day')%>%
              ungroup(), by = c('feature_number', 'Organism'))%>%
  select(-c(TF, log2_change, network, numberOfNodes, log2_change, ecoNetConsensusScore))%>%
  group_by(ecoNetConsensus, Replicate, Organism)%>%
  summarize_if(is.numeric, sum, na.rm  = TRUE)%>%
  ungroup()%>%
  group_by(ecoNetConsensus, Organism)%>%
  summarize_if(is.numeric, mean, na.rm  = TRUE)%>%
  ungroup()%>%
  spread(Organism, T0)%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  gather(Organism, T0, 2:ncol(.))%>%
  group_by(ecoNetConsensus)%>%
  mutate(T0 = zscore(T0))
  # spread(Organism, T0)
  
depletion_levels <- depletion_hc$ecoNetConsensus%>% 
  as.factor()
organism_levels <- depletion_hc$Organism%>%
  as.factor()%>%
  fct_relevel(c('CCA', 'Turf', 'Dictyota'))

pdf('plots/depletoliteHeatmap.pdf', width = 20, height = 15)
depletion_hc%>%
  ggplot(aes(Organism, ecoNetConsensus, fill = T0)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_y_discrete(limits = rev(levels(depletion_levels))) +
  scale_x_discrete(limits = levels(organism_levels))
dev.off()

write_csv(depletion_hc%>% 
            spread(Organism, T0), './analysis/depletion_hc.csv')

# VIZUALIZATIONS -- Nutrition/lability vs microbial community change -------
lability_val <- log2_change_vals%>%
  inner_join(feature_stats_wdf%>%
               ungroup()%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))

pdf('./plots/weighted_log2_boxplot.pdf', width = 12, height = 10)
lability_val%>%
  filter(activity == 'depletolite')%>%
=======
# VIZUALIZATIONS -- Venn Diagram ------------------------------------------
venn <- boxplot_df%>%
  ungroup()%>%
  filter(activity == 'depletolite')%>%
  select(feature_number, Organism)%>%
  unique()

venn_list <- venn$feature_number%>%
  as.vector()

venn_split <- split(venn_list, venn$Organism)

venn.diagram(venn_split, 
             './plots/venn_depletolites.png',
             category.names = c('CCA', 'Dic', 'Poc', 'Por', 'T'),
             imagetype = 'png',  
             height = 1500, 
             width = 1800,
             resolution = 300,
             fill = org_colors_no_water)

venn_features <- calculate.overlap(venn_split)[c('a1', 'a2', 'a3', 'a4', 'a5', 'a13', 'a31', 'a10')]

names(venn_features) <- c('CCA', 'Dictyota', 'Poc', 'Por', 'Turf', 'coral', 'Primary producers', 'fleshy algae')

venn_df <- venn_features%>%
  unlist()%>%
  as.data.frame()%>%
  rename('feature_number' = 1)%>%
  rownames_to_column(var = 'Intersection')%>%
  mutate(Intersection = gsub('([0-9])', '', Intersection))%>%
  left_join(metadata%>% 
              select(combined_ID, binary_ID, feature_number, network), 
            by= 'feature_number')%>%
  left_join(ecoNet%>%
              select(-scan)%>%
              filter(network != '-1')%>%
              unique(), by = 'network')
  


# VIZUALIZATIONS -- Nutrition/lability vs microbial community change -------
lability_val <- log2_change_vals%>%
  inner_join(feature_stats_wdf%>%
               ungroup()%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))

pdf('./plots/weighted_log2_boxplot.pdf', width = 12, height = 10)
lability_val%>%
  filter(activity == 'depletolite')%>%
  group_by(Organism, Replicate)%>%
  mutate(log10 = log10(T0 + 1),
         weighted_log2 = spatstat::weighted.median(log2_change, log10))%>%
  # mutate(weighted_lability = mean(modeled_lab))%>%
  group_by(Organism)%>%
  mutate(x_err = sd(weighted_log2))%>%
  select(Organism, weighted_log2, x_err)%>%
  unique()%>%
  ggplot() +
  geom_boxplot(aes(Organism, weighted_log2, color = '#EBCC2A'), lwd = 2) +
  labs(y = bquote(atop(Metabolite ~pool ~xic, weighted ~fold ~change))) +
  scale_color_manual(values = '#EBCC2A') +
  theme(
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    axis.text.x = element_text(size = 25, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_line('grey'), # get rid of major grid
    panel.grid.minor = element_line('grey'), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    legend.position = 'top'
  )
dev.off()

# VIZUALIZATIONS -- Pooled Weighted log2 change by xic---------------------------
sum_xic_x_log2 <- lability_val%>%
  group_by(Organism, Replicate)%>%
  mutate(log10 = log10(T0 + 1),
         weighted_log2 = spatstat::weighted.median(log2_change, log10))%>%
  # mutate(weighted_lability = mean(modeled_lab))%>%
  group_by(Organism)%>%
  mutate(x_err = sd(weighted_log2))%>%
  select(Organism, weighted_log2, x_err)%>%
  unique()%>%
  summarize_if(is.numeric, mean)%>%
  left_join(fcm_T0_5%>%
              filter(DayNight == 'Day')%>%
              group_by(Organism)%>%
              mutate(y_err = sd(cells_ul))%>%
              summarize_if(is.numeric, mean),
            by = c('Organism'))

sum_xic_x_log2_stats <- lability_val%>%
  group_by(Organism, Replicate)%>%
  mutate(log10 = log10(T0 + 1),
         weighted_log2 = spatstat::weighted.median(log2_change, log10))%>%
  # mutate(weighted_lability = mean(modeled_lab))%>%
  group_by(Organism)%>%
  mutate(x_err = sd(weighted_log2))%>%
  select(Organism, Replicate, weighted_log2, x_err)%>%
  unique()%>%
  left_join(fcm_T0_5%>%
              filter(DayNight == 'Day'),
            by = c('Organism', 'Replicate'))


weight_lability_lm <- sum_xic_x_log2_stats%>%
  lm(cells_ul ~ weighted_log2, data = .)

xic_p <- (weight_lability_lm%>% 
            tidy()%>% 
            filter(term == 'weighted_log2'))$p.value

xic_f <- (weight_lability_lm%>% 
            tidy()%>% 
            filter(term == 'weighted_log2'))$statistic

xic_slope <- weight_lability_lm$coefficients["weighted_log2"]
xic_intercept <- weight_lability_lm$coefficients["(Intercept)"]


xic_r2 <- summary(weight_lability_lm)$adj.r.squared

pdf("./plots/xicweighted_log2_fcm.pdf", width = 12, height = 10)
sum_xic_x_log2%>%
  ggplot(aes(weighted_log2, cells_ul)) +
  geom_point(aes(color = Organism), stat = 'identity', size = 5) +
  geom_errorbar(aes(ymin = cells_ul - y_err, ymax = cells_ul + y_err)) +
  geom_errorbarh(aes(xmin = weighted_log2 - x_err, xmax = weighted_log2 + x_err))+
  scale_color_manual(values = org_colors_no_water) +
  geom_smooth(method = 'lm') +
  labs(y = bquote(Specific ~growth ~rate ~('Cells'~µL^-1 ~hr^-1)), x = "Metabolite pool xic weighted fold change") +
  geom_text(aes(x = -1, y = 0.09,
                label = paste("p-value: ", xic_p%>%
                                formatC(format = "e", digits = 2), sep = "")), size = 9) +
  geom_text(aes(x = -1, y = 0.086,
                label = paste("F statistic: ", xic_f%>%
                                round(digits = 4), sep = "")), size = 9) +
  geom_text(aes(x = -1, y = 0.082,
                label = paste("r²: ", xic_r2%>%
                                round(digits = 4), sep = "")), size = 9) +
  geom_text(aes(x = -1, y = 0.078,
                label = paste("Cells µL^-1", " = ", xic_slope%>%
                                round(digits = 2), "*lability + ", xic_intercept%>%
                                round(digits = 2), sep = "")), size = 9) +
  theme(
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_line('grey'), # get rid of major grid
    panel.grid.minor = element_line('grey'), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    legend.position = 'top'
  )
dev.off()

# sum_xic_x_log2%>%
#   # filter(!is.na(activity))%>%
#   filter(DayNight == 'Day')%>%
#   ggplot(aes(change_x_xic, cells_ul)) +
#   geom_point(aes(color = Organism), stat = 'summary', fun.y = 'mean', size = 5) +
#   scale_color_manual(values = org_colors_no_water) + 
#   geom_smooth(method = 'lm') +
#   theme(
#     plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
#     axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
#     axis.text.y = element_text(size = 20),
#     axis.title = element_text(size = 20),
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     panel.grid.major = element_line('grey'), # get rid of major grid
#     panel.grid.minor = element_line('grey'), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     legend.box.background = element_rect(fill = "transparent"),
#     legend.position = 'top'
#   )
# 
# sum_xic_x_log2%>%
#   # filter(!is.na(activity))%>%
#   filter(DayNight == 'Day',
#          NOSC < 0)%>%
#   ggplot(aes(log2_change, cells_ul)) +
#   geom_point(aes(color = Organism), stat = 'summary', fun.y = 'mean', size = 5) +
#   scale_color_manual(values = org_colors_no_water) + 
#   geom_smooth(method = 'lm') +
#   theme(
#     plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
#     axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
#     axis.text.y = element_text(size = 20),
#     axis.title = element_text(size = 20),
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     panel.grid.major = element_line('grey'), # get rid of major grid
#     panel.grid.minor = element_line('grey'), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     legend.box.background = element_rect(fill = "transparent"),
#     legend.position = 'top'
#   )
# 
# sum_xic_x_log2%>%
#   # filter(!is.na(activity))%>%
#   filter(DayNight == 'Day')%>%
#   ggplot(aes(sum_xic, cells_ul)) +
#   geom_point(aes(color = Organism), stat = 'summary', fun.y = 'mean', size = 5) +
#   scale_color_manual(values = org_colors_no_water) + 
#   geom_smooth(method = 'lm') +
#   theme(
#     plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
#     axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
#     axis.text.y = element_text(size = 20),
#     axis.title = element_text(size = 20),
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     panel.grid.major = element_line('grey'), # get rid of major grid
#     panel.grid.minor = element_line('grey'), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     legend.box.background = element_rect(fill = "transparent"),
#     legend.position = 'top'
#   )




# VIZUALIZATIONS -- checking ratio variance -------------------------------
n_mulreg%>%
  select(Organism, log_nc, log_oc, log_hc)%>%
  gather(element, ratio, 2:ncol(.))%>%
  ggplot(aes(Organism, ratio, color = Organism)) +
  geom_boxplot() +
  facet_wrap(~element, scale = 'free_y', ncol = 1) +
  scale_color_manual(values = org_colors_no_water) +
  theme(
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_line('grey'), # get rid of major grid
    panel.grid.minor = element_line('grey'), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    legend.position = 'top'
  )

#Run ANOVAS in conjunction with these and this will be a supplemental figure

# VIZUALIZATIONS -- Lability value from multiple regressions --------------
weighted_lability <- feature_stats_wdf%>%
  filter(Timepoint == 'T0',
         DayNight == 'Day')%>%
  left_join(mul_reg%>%
              ungroup()%>%
              mutate(modeled_lability  = case_when(N > 0 & O > 0 ~ (log_nc*n_n_coe + log_oc*n_o_coe + log_hc*n_h_coe + NOSC*n_nosc_coe + logxic*n_xic_coe + n_intercept),
                                                   N == 0 & O > 0 ~ (logxic*o_xic_coe + NOSC*o_nosc_coe + o_intercept),
                                                   TRUE ~ NA_real_),
                     model_num = case_when(N > 0 & O > 0 ~ 'Nitrogen',
                                           N == 0 & O > 0 ~ 'Oxygen',
                                           TRUE ~ 'none'))%>%
              select(feature_number, Organism, modeled_lability)%>%
              unique(),
            by = c('feature_number', 'Organism'))%>%
  filter(!is.na(modeled_lability))%>%
>>>>>>> 866b30c22f5339816f0bda22d3e41936737b507f
  group_by(Organism, Replicate)%>%
  mutate(log10 = log10(T0 + 1),
         weighted_log2 = spatstat::weighted.median(log2_change, log10))%>%
  # mutate(weighted_lability = mean(modeled_lab))%>%
  group_by(Organism)%>%
  mutate(x_err = sd(weighted_log2))%>%
  select(Organism, weighted_log2, x_err)%>%
  unique()%>%
  ggplot() +
  geom_boxplot(aes(Organism, weighted_log2, color = '#EBCC2A'), lwd = 2) +
  labs(y = bquote(atop(Metabolite ~pool ~xic, weighted ~fold ~change))) +
  scale_color_manual(values = '#EBCC2A') +
  theme(
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    axis.text.x = element_text(size = 25, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_line('grey'), # get rid of major grid
    panel.grid.minor = element_line('grey'), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    legend.position = 'top'
  )
dev.off()

<<<<<<< HEAD
# VIZUALIZATIONS -- Pooled Weighted log2 change by xic---------------------------
sum_xic_x_log2 <- lability_val%>%
  group_by(Organism, Replicate)%>%
  mutate(log10 = log10(T0 + 1),
         weighted_log2 = spatstat::weighted.median(log2_change, log10))%>%
  # mutate(weighted_lability = mean(modeled_lab))%>%
  group_by(Organism)%>%
  mutate(x_err = sd(weighted_log2))%>%
  select(Organism, weighted_log2, x_err)%>%
  unique()%>%
  summarize_if(is.numeric, mean)%>%
  left_join(fcm_T0_5%>%
              filter(DayNight == 'Day')%>%
              group_by(Organism)%>%
              mutate(y_err = sd(cells_ul))%>%
              summarize_if(is.numeric, mean),
            by = c('Organism'))

sum_xic_x_log2_stats <- lability_val%>%
  group_by(Organism, Replicate)%>%
  mutate(log10 = log10(T0 + 1),
         weighted_log2 = spatstat::weighted.median(log2_change, log10))%>%
  # mutate(weighted_lability = mean(modeled_lab))%>%
  group_by(Organism)%>%
  mutate(x_err = sd(weighted_log2))%>%
  select(Organism, Replicate, weighted_log2, x_err)%>%
  unique()%>%
  left_join(fcm_T0_5%>%
              filter(DayNight == 'Day'),
            by = c('Organism', 'Replicate'))

=======
>>>>>>> 866b30c22f5339816f0bda22d3e41936737b507f

weight_lability_lm <- sum_xic_x_log2_stats%>%
  lm(cells_ul ~ weighted_log2, data = .)

xic_p <- (weight_lability_lm%>% 
            tidy()%>% 
            filter(term == 'weighted_log2'))$p.value

<<<<<<< HEAD
xic_f <- (weight_lability_lm%>% 
            tidy()%>% 
            filter(term == 'weighted_log2'))$statistic
=======
wlab_f <- (weight_lability_lm%>% 
             tidy()%>% 
             filter(term == 'weighted_lability'))$statistic
>>>>>>> 866b30c22f5339816f0bda22d3e41936737b507f

xic_slope <- weight_lability_lm$coefficients["weighted_log2"]
xic_intercept <- weight_lability_lm$coefficients["(Intercept)"]


xic_r2 <- summary(weight_lability_lm)$adj.r.squared

<<<<<<< HEAD
pdf("./plots/xicweighted_log2_fcm.pdf", width = 12, height = 10)
sum_xic_x_log2%>%
  ggplot(aes(weighted_log2, cells_ul)) +
=======
#Plotting
pdf("./plots/weighted_lability_nultiple_regressions_aic.pdf", width = 15, height = 10)
weighted_lability%>%
  ggplot(aes(weighted_lability, cells_ul)) +
>>>>>>> 866b30c22f5339816f0bda22d3e41936737b507f
  geom_point(aes(color = Organism), stat = 'identity', size = 5) +
  geom_errorbar(aes(ymin = cells_ul - y_err, ymax = cells_ul + y_err)) +
  geom_errorbarh(aes(xmin = weighted_log2 - x_err, xmax = weighted_log2 + x_err))+
  scale_color_manual(values = org_colors_no_water) +
  geom_smooth(method = 'lm') +
<<<<<<< HEAD
  labs(y = bquote(Specific ~growth ~rate ~('Cells'~µL^-1 ~hr^-1)), x = "Metabolite pool xic weighted fold change") +
  geom_text(aes(x = -1, y = 0.09,
                label = paste("p-value: ", xic_p%>%
                                formatC(format = "e", digits = 2), sep = "")), size = 9) +
  geom_text(aes(x = -1, y = 0.086,
                label = paste("F statistic: ", xic_f%>%
                                round(digits = 4), sep = "")), size = 9) +
  geom_text(aes(x = -1, y = 0.082,
                label = paste("r²: ", xic_r2%>%
                                round(digits = 4), sep = "")), size = 9) +
  geom_text(aes(x = -1, y = 0.078,
                label = paste("Cells µL^-1", " = ", xic_slope%>%
                                round(digits = 2), "*lability + ", xic_intercept%>%
                                round(digits = 2), sep = "")), size = 9) +
  theme(
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_line('grey'), # get rid of major grid
    panel.grid.minor = element_line('grey'), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    legend.position = 'top'
  )
dev.off()

# sum_xic_x_log2%>%
#   # filter(!is.na(activity))%>%
#   filter(DayNight == 'Day')%>%
#   ggplot(aes(change_x_xic, cells_ul)) +
#   geom_point(aes(color = Organism), stat = 'summary', fun.y = 'mean', size = 5) +
#   scale_color_manual(values = org_colors_no_water) + 
#   geom_smooth(method = 'lm') +
#   theme(
#     plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
#     axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
#     axis.text.y = element_text(size = 20),
#     axis.title = element_text(size = 20),
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     panel.grid.major = element_line('grey'), # get rid of major grid
#     panel.grid.minor = element_line('grey'), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     legend.box.background = element_rect(fill = "transparent"),
#     legend.position = 'top'
#   )
# 
# sum_xic_x_log2%>%
#   # filter(!is.na(activity))%>%
#   filter(DayNight == 'Day',
#          NOSC < 0)%>%
#   ggplot(aes(log2_change, cells_ul)) +
#   geom_point(aes(color = Organism), stat = 'summary', fun.y = 'mean', size = 5) +
#   scale_color_manual(values = org_colors_no_water) + 
#   geom_smooth(method = 'lm') +
#   theme(
#     plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
#     axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
#     axis.text.y = element_text(size = 20),
#     axis.title = element_text(size = 20),
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     panel.grid.major = element_line('grey'), # get rid of major grid
#     panel.grid.minor = element_line('grey'), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     legend.box.background = element_rect(fill = "transparent"),
#     legend.position = 'top'
#   )
# 
# sum_xic_x_log2%>%
#   # filter(!is.na(activity))%>%
#   filter(DayNight == 'Day')%>%
#   ggplot(aes(sum_xic, cells_ul)) +
#   geom_point(aes(color = Organism), stat = 'summary', fun.y = 'mean', size = 5) +
#   scale_color_manual(values = org_colors_no_water) + 
#   geom_smooth(method = 'lm') +
#   theme(
#     plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
#     axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
#     axis.text.y = element_text(size = 20),
#     axis.title = element_text(size = 20),
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     panel.grid.major = element_line('grey'), # get rid of major grid
#     panel.grid.minor = element_line('grey'), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     legend.box.background = element_rect(fill = "transparent"),
#     legend.position = 'top'
#   )




# VIZUALIZATIONS -- checking ratio variance -------------------------------
n_mulreg%>%
  select(Organism, log_nc, log_oc, log_hc)%>%
  gather(element, ratio, 2:ncol(.))%>%
  ggplot(aes(Organism, ratio, color = Organism)) +
  geom_boxplot() +
  facet_wrap(~element, scale = 'free_y', ncol = 1) +
  scale_color_manual(values = org_colors_no_water) +
=======
  labs(y = bquote(Specific ~growth ~rate ~('Cells'~µL^-1 ~hr^-1)), x = "Metabolite pool modeled lability") +
  # geom_text(aes(x = -1.15, y = 820,
  #               label = paste("p-value: ", wlab_p%>%
  #                               formatC(format = "e", digits = 2), sep = "")), size = 9) +
  # geom_text(aes(x = -1.15, y = 780,
  #               label = paste("F statistic: ", wlab_f%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = -1.15, y = 740,
  #               label = paste("r²: ", wlab_r2%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = -1.15, y = 700,
  #               label = paste("Cells µL^-1", " = ", wlab_slope%>%
#                               round(digits = 2), "*lability - ", -wlab_intercept%>%
#                               round(digits = 2), sep = "")), size = 9) +
theme(
  plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
  axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
  axis.text.y = element_text(size = 20),
  axis.title = element_text(size = 20),
  panel.background = element_rect(fill = "transparent"), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  panel.grid.major = element_line('grey'), # get rid of major grid
  panel.grid.minor = element_line('grey'), # get rid of minor grid
  legend.background = element_rect(fill = "transparent"), # get rid of legend bg
  legend.box.background = element_rect(fill = "transparent"),
  legend.position = 'top'
)
dev.off()



lability_val_check <- (mul_reg%>%
                         ungroup()%>%
                         mutate(multiple_reg_lability = case_when(N > 0 & O > 0 & NOSC < 0 ~ (NOSC*n_nosc_coe + log_nc*n_n_coe + oc*n_o_coe + hc*n_h_coe + n_intercept),
                                                                  N == 0 & NOSC < 0 | O == 0 & NOSC < 0 ~ (`row m/z`*no_mass_coe + hc*no_h_coe + no_intercept),
                                                                  is.na(NOSC) | NOSC >= 0 ~ (`row m/z`*m_mass + m_intercept),
                                                                  TRUE ~ NA_real_))%>%
                         filter(DayNight == 'Day'))$multiple_reg_lability


car::qqPlot(lability_val_check, 
            ylab = "Lability value", xlab = "Normal quantiles",
            main = 'QQ-plot: Lability value')




# VIZUALIZATIONS -- CANOPUS machine learning model ------------------------
cutoff_canopus <- canopus_anotations%>%
  gather(classification, value, 2:ncol(.))%>%
  group_by(classification)%>%
  filter(max(value) >= 0.4,
         sd(value) > 0.05)%>%
  spread(classification, value)

canopus_mulreg <- mul_reg%>%
  filter(NOSC < 0)%>%
  select(feature_number, Organism, Replicate, log2_change, `row m/z`, NOSC, logxic)%>%
  inner_join(cutoff_canopus%>%
               rename(feature_number = name)%>%
               mutate(feature_number = as.character(feature_number))%>%
               mutate_if(is.numeric, angular_transform), by = 'feature_number')%>%
  left_join(sirius_zodiac_anotations%>%
              select(feature_number, ZodiacScore)%>%
              mutate(feature_number = as.character(feature_number)), by = 'feature_number')%>%
  filter(ZodiacScore >= .98)%>%
  select(-ZodiacScore)%>%
  unite(row, c('feature_number', 'Organism', 'Replicate'), sep = "_")%>%
  column_to_rownames('row')

recipe <- canopus_mulreg%>%
  recipe(log2_change ~ .)%>%
  step_center(all_predictors())%>%
  step_scale(all_predictors())


##canopus training set testing

canopus_split <- rsample::initial_split(canopus_mulreg, 9/10)

canopus_training <- rsample::training(canopus_split)
canopus_test <- rsample::testing(canopus_split)

canopus_x <- canopus_training[2:ncol(canopus_training)]%>%
  as.matrix()

canopus_y <- canopus_training['log2_change']%>%
  as.matrix()

canopus_test_x <-  canopus_test[2:ncol(canopus_test)]%>%
  as.matrix()

canopus_test_y <- canopus_test['log2_change']%>%
  as.matrix()

#running 
glmnet_canopus <- cv.glmnet(canopus_x, canopus_y, parallel = TRUE)

lambda_se <- glmnet_canopus$lambda.1se

#running model with parameters after lasso
canopus_model <- glmnet(canopus_x, canopus_y, intercept = TRUE, lambda = lambda_se, standarize = TRUE)

canopus_coeff <- coef(canopus_model)%>%
  as.matrix()%>%
  as.data.frame()%>%
  rename(coefficients = 1)%>%
  rownames_to_column("variable")%>%
  filter(coefficients != 0)

intercept <- canopus_coeff[[1,2]]

predicted_training <- predict(canopus_model, newx = canopus_x)
predicted <- predict(canopus_model, newx = canopus_test_x)

rootmse <- rmse_vec(as.vector(canopus_test_y), as.vector(predicted))
rootmse_training <- rmse_vec(as.vector(canopus_y), as.vector(predicted_training))

canopus_model$dev.ratio

predicted%>%
  as.data.frame()%>%
  rename(predicted = 1)%>%
  add_column(actuals = canopus_test_y)%>%
  rownames_to_column('sample')%>%
  ggplot() +
  geom_point(aes(sample, actuals), color = 'blue') +
  geom_point(aes(sample, predicted), color = 'red')

sink('./analysis/model_rmse_r2.txt')
rootmse
rootmse_training
canopus_model$dev.ratio
sink()

canopus_weighted_lability <- feature_stats_wdf%>%
  filter(Timepoint == 'T0',
         DayNight == 'Day')%>%
  inner_join(canopus_mulreg%>%
               rownames_to_column("row")%>%
               separate(row, c('feature_number', 'Organism', 'Replicate'), sep = "_")%>%
               select(-log2_change)%>%
               gather(variable, value, 4:ncol(.))%>%
               inner_join(rf_canopus_coeff, by = 'variable')%>%
               mutate(modeled_lab = value*coefficients)%>%
               select(feature_number, Organism, Replicate, modeled_lab)%>%
               group_by(feature_number, Organism, Replicate)%>%
               summarize_if(is.numeric, sum)%>%
               ungroup()%>%
               mutate(modeled_lab = modeled_lab + rf_intercept),
             by = c('feature_number', 'Replicate', 'Organism'))%>%
  group_by(Organism, Replicate)%>%
  # mutate(weighted_lability = sum(modeled_lab))%>%
  mutate(weighted_lability = spatstat::weighted.median(modeled_lab, log10))%>%
  # mutate(weighted_lability = weighted.mean(modeled_lab, log10))%>%
  group_by(Organism)%>%
  mutate(x_err = sd(weighted_lability))%>%
  select(Organism, weighted_lability, x_err)%>%
  unique()%>%
  summarize_if(is.numeric, mean)%>%
  left_join(fcm_T0_5%>%
              filter(DayNight == 'Day')%>%
              group_by(Organism)%>%
              mutate(y_err = sd(cells_ul))%>%
              summarize_if(is.numeric, mean),
            by = c('Organism'))

# # dredge_n_select <- MuMIn::dredge(lm_test_canopus)
# # test <- AIC(lm_test_canopus)
# 

canopus_weighted_lability_stats <- feature_stats_wdf%>%
  filter(Timepoint == 'T0',
         DayNight == 'Day')%>%
  inner_join(canopus_mulreg%>%
               rownames_to_column("row")%>%
               separate(row, c('feature_number', 'Organism', 'Replicate'), sep = "_")%>%
               select(-log2_change)%>%
               gather(variable, value, 4:ncol(.))%>%
               inner_join(rf_canopus_coeff, by = 'variable')%>%
               mutate(modeled_lab = value*coefficients)%>%
               select(feature_number, Organism, modeled_lab)%>%
               group_by(feature_number, Organism)%>%
               summarize_if(is.numeric, sum)%>%
               ungroup()%>%
               mutate(modeled_lab = modeled_lab + rf_intercept),
             by = c('feature_number', 'Organism'))%>%
  group_by(Organism, Replicate)%>%
  mutate(weighted_lability = spatstat::weighted.median(modeled_lab, log10))%>%
  # mutate(weighted_lability = mean(modeled_lab))%>%
  group_by(Organism)%>%
  mutate(x_err = sd(weighted_lability))%>%
  select(Organism, weighted_lability, x_err)%>%
  unique()%>%
  left_join(fcm_T0_5%>%
              filter(DayNight == 'Day')%>%
              group_by(Organism)%>%
              mutate(y_err = sd(cells_ul))%>%
              summarize_if(is.numeric, mean),
            by = c('Organism'))

#Linear model
weight_lability_lm <- canopus_weighted_lability_stats%>%
  lm(cells_ul ~ weighted_lability, data = .)

wlab_p <- (weight_lability_lm%>% 
             tidy()%>% 
             filter(term == 'weighted_lability'))$p.value

wlab_f <- (weight_lability_lm%>% 
             tidy()%>% 
             filter(term == 'weighted_lability'))$statistic

wlab_slope <- weight_lability_lm$coefficients["weighted_lability"]
wlab_intercept <- weight_lability_lm$coefficients["(Intercept)"]


wlab_r2 <- summary(weight_lability_lm)$adj.r.squared

pdf('~/Documents/GitHub/DORCIERR/data/plots/canopus_weighted_lability.pdf', width = 12, height = 10)
canopus_weighted_lability%>%
  ggplot(aes(weighted_lability, cells_ul)) +
  geom_point(aes(color = Organism), stat = 'identity', size = 5) +
  geom_errorbar(aes(ymin = cells_ul - y_err, ymax = cells_ul + y_err)) +
  geom_errorbarh(aes(xmin = weighted_lability - x_err, xmax = weighted_lability + x_err))+
  scale_color_manual(values = org_colors_no_water) +
  geom_smooth(method = 'lm') +
  labs(y = bquote(Specific ~growth ~rate ~('Cells'~µL^-1 ~hr^-1)), x = "Metabolite pool modeled lability") +
  # geom_text(aes(x = -1.3, y = 0.09,
  #               label = paste("p-value: ", wlab_p%>%
  #                               formatC(format = "e", digits = 2), sep = "")), size = 9) +
  # geom_text(aes(x = -1.3, y = 0.086,
  #               label = paste("F statistic: ", wlab_f%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = -1.3, y = 0.082,
  #               label = paste("r²: ", wlab_r2%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = -1.3, y = 0.078,
  #               label = paste("Cells µL^-1", " = ", wlab_slope%>%
#                               round(digits = 2), "*lability + ", wlab_intercept%>%
#                               round(digits = 2), sep = "")), size = 9) +
theme(
  plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
  axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
  axis.text.y = element_text(size = 20),
  axis.title = element_text(size = 20),
  panel.background = element_rect(fill = "transparent"), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  panel.grid.major = element_line('grey'), # get rid of major grid
  panel.grid.minor = element_line('grey'), # get rid of minor grid
  legend.background = element_rect(fill = "transparent"), # get rid of legend bg
  legend.box.background = element_rect(fill = "transparent"),
  legend.position = 'top'
)
dev.off()

#checking gaussian
pdf("./plots/canopus_gaussian_check.pdf")
gaussian_check <- prep(recipe)%>%
  juice()%>%
  gather(variable, value)%>%
  inner_join(canopus_coeff%>% select(variable),
             by = 'variable')%>%
  group_by(variable)%>%
  nest()%>%
  mutate(non_transformed = map(data, ~ .x$value%>%
                                 na.omit()%>%
                                 car::qqPlot(ylab = variable, xlab = "Normal quantiles",
                                             main = paste(variable, " Centered-Scaled"))))
dev.off()


# VIZUALIZATIONS -- CANOPUS Random Forest ---------------------------------
names(canopus_training) <- make.names(names(canopus_training))

rf <- foreach(ntree = 100, .combine=randomForest::combine, .multicombine=TRUE, .packages='randomForest') %dopar% {
  randomForest(log2_change ~ ., canopus_training, importance = TRUE, proximity = TRUE, ntree = ntree)
}

mda_canopus <- rf$importance%>%
  as.data.frame()%>%
  rownames_to_column("feature")

# selecting variables to use in model
mean_se <- (mda_canopus%>%
              rename(inc_mse = `%IncMSE`)%>%
              filter(inc_mse >= mean(inc_mse) + sd(inc_mse)/sqrt(length(.))))$feature%>% 
  as.vector()

rf_xline <- length(mean_se)

ggplot(mda_canopus, aes(x= reorder(feature, -`%IncMSE`), y = `%IncMSE`)) +
  geom_point(stat = "identity") +
  geom_vline(xintercept = rf_xline + 0.5, color = 'red')



#model
training_rf_variables <- canopus_training%>%
  select('log2_change', all_of(mean_se))

test_rf_variable <- canopus_test

rf_can_model <- lm(log2_change ~ ., data = training_rf_variables)

rf_canopus_coeff <- coef(rf_can_model)%>%
  as.matrix()%>%
  as.data.frame()%>%
  rename(coefficients = 1)%>%
  rownames_to_column("variable")%>%
  filter(coefficients != 0)

rf_intercept = rf_canopus_coeff[[1,2]]

# VIZUALIZATIONS -- FCM ---------------------------------------------------
fcm_graphing <- fcm_wdf%>%
  filter(!Organism == 'Influent',
         !Organism == 'Offshore',
         DayNight == "Day")%>%
  mutate(Hours = case_when(Timepoint == "T0" ~ 0,
                           Timepoint == "T1" ~ 4,
                           Timepoint == "T2" ~ 13,
                           Timepoint == "T3" ~ 17,
                           Timepoint == "T4" ~ 23,
                           Timepoint == "T5" ~ 28,
                           Timepoint == "T6" ~ 37,
                           Timepoint == "TF" ~ 48))%>%
  group_by(Organism, Timepoint)%>%
  mutate(st_err = sd(`Cells µL-1`))%>%
  summarize_if(is.numeric, mean)

pdf("~/Documents/GitHub/DORCIERR/data/plots/FCM_day.pdf", width = 10, height = 10)
fcm_graphing%>%
  ggplot(aes(x= Hours, y = `Cells µL-1`, color = Organism))+
  geom_point(stat = "identity", size = 5) +
  geom_errorbar(aes(ymin = `Cells µL-1` - st_err, ymax = `Cells µL-1` + st_err)) +
  geom_line(aes(group = Organism)) +
  scale_color_manual(values = c(org_colors_no_water, "#3B9AB2")) +
  labs(y = bquote(Cells ~µL^-1)) +
  # facet_wrap(~ DayNight) +
  # scale_color_manual(values = c("darkorchid3", "#50A45C", "#AF814B", "#5BBCD6")) +
  scale_y_continuous(limits = c(0,900), breaks= seq(0, 900, 100)) +
  scale_x_continuous(limits = c(0,50), breaks = seq(0, 50, 5)) +
>>>>>>> 866b30c22f5339816f0bda22d3e41936737b507f
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.x = element_blank(), # get rid of major grid
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.position = 'top',
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )
dev.off()

<<<<<<< HEAD
#Run ANOVAS in conjunction with these and this will be a supplemental figure

# VIZUALIZATIONS -- Lability value from multiple regressions --------------
weighted_lability <- feature_stats_wdf%>%
  filter(Timepoint == 'T0',
         DayNight == 'Day')%>%
  left_join(mul_reg%>%
              ungroup()%>%
              mutate(modeled_lability  = case_when(N > 0 & O > 0 ~ (log_nc*n_n_coe + log_oc*n_o_coe + log_hc*n_h_coe + NOSC*n_nosc_coe + logxic*n_xic_coe + n_intercept),
                                                   N == 0 & O > 0 ~ (logxic*o_xic_coe + NOSC*o_nosc_coe + o_intercept),
                                                   TRUE ~ NA_real_),
                     model_num = case_when(N > 0 & O > 0 ~ 'Nitrogen',
                                           N == 0 & O > 0 ~ 'Oxygen',
                                           TRUE ~ 'none'))%>%
              select(feature_number, Organism, modeled_lability)%>%
              unique(),
            by = c('feature_number', 'Organism'))%>%
  filter(!is.na(modeled_lability))%>%
  group_by(Organism, Replicate)%>%
  mutate(weighted_lability = spatstat::weighted.median(modeled_lability, log10))%>%
  group_by(Organism)%>%
  mutate(x_err = sd(weighted_lability))%>%
  select(Organism, weighted_lability, x_err)%>%
  unique()%>%
  summarize_if(is.numeric, mean)%>%
  left_join(fcm_T0_5%>%
              filter(DayNight == 'Day')%>%
              group_by(Organism)%>%
              mutate(y_err = sd(cells_ul))%>%
              summarize_if(is.numeric, mean), 
            by = c('Organism'))


# Linear model weighted lability
weight_lability_lm <- weighted_lability%>%
  lm(cells_ul ~ weighted_lability, data = .)

wlab_p <- (weight_lability_lm%>% 
             tidy()%>% 
             filter(term == 'weighted_lability'))$p.value

wlab_f <- (weight_lability_lm%>% 
             tidy()%>% 
             filter(term == 'weighted_lability'))$statistic

wlab_slope <- weight_lability_lm$coefficients["weighted_lability"]
wlab_intercept <- weight_lability_lm$coefficients["(Intercept)"]


wlab_r2 <- summary(weight_lability_lm)$adj.r.squared

#Plotting
pdf("./plots/weighted_lability_nultiple_regressions_aic.pdf", width = 15, height = 10)
weighted_lability%>%
  ggplot(aes(weighted_lability, cells_ul)) +
  geom_point(aes(color = Organism), stat = 'identity', size = 5) +
  geom_errorbar(aes(ymin = cells_ul - y_err, ymax = cells_ul + y_err)) +
  geom_errorbarh(aes(xmin = weighted_lability - x_err, xmax = weighted_lability + x_err))+
  scale_color_manual(values = org_colors_no_water) + 
  geom_smooth(method = 'lm') +
  labs(y = bquote(Specific ~growth ~rate ~('Cells'~µL^-1 ~hr^-1)), x = "Metabolite pool modeled lability") +
  # geom_text(aes(x = -1.15, y = 820,
  #               label = paste("p-value: ", wlab_p%>%
  #                               formatC(format = "e", digits = 2), sep = "")), size = 9) +
  # geom_text(aes(x = -1.15, y = 780,
  #               label = paste("F statistic: ", wlab_f%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = -1.15, y = 740,
  #               label = paste("r²: ", wlab_r2%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = -1.15, y = 700,
  #               label = paste("Cells µL^-1", " = ", wlab_slope%>%
#                               round(digits = 2), "*lability - ", -wlab_intercept%>%
#                               round(digits = 2), sep = "")), size = 9) +
theme(
  plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
  axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
  axis.text.y = element_text(size = 20),
  axis.title = element_text(size = 20),
  panel.background = element_rect(fill = "transparent"), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  panel.grid.major = element_line('grey'), # get rid of major grid
  panel.grid.minor = element_line('grey'), # get rid of minor grid
  legend.background = element_rect(fill = "transparent"), # get rid of legend bg
  legend.box.background = element_rect(fill = "transparent"),
  legend.position = 'top'
)
dev.off()



lability_val_check <- (mul_reg%>%
                         ungroup()%>%
                         mutate(multiple_reg_lability = case_when(N > 0 & O > 0 & NOSC < 0 ~ (NOSC*n_nosc_coe + log_nc*n_n_coe + oc*n_o_coe + hc*n_h_coe + n_intercept),
                                                                  N == 0 & NOSC < 0 | O == 0 & NOSC < 0 ~ (`row m/z`*no_mass_coe + hc*no_h_coe + no_intercept),
                                                                  is.na(NOSC) | NOSC >= 0 ~ (`row m/z`*m_mass + m_intercept),
                                                                  TRUE ~ NA_real_))%>%
                         filter(DayNight == 'Day'))$multiple_reg_lability


car::qqPlot(lability_val_check, 
            ylab = "Lability value", xlab = "Normal quantiles",
            main = 'QQ-plot: Lability value')




# VIZUALIZATIONS -- CANOPUS machine learning model ------------------------
cutoff_canopus <- canopus_anotations%>%
  gather(classification, value, 2:ncol(.))%>%
  group_by(classification)%>%
  filter(max(value) >= 0.4,
         sd(value) > 0.05)%>%
  spread(classification, value)

canopus_mulreg <- mul_reg%>%
  filter(NOSC < 0)%>%
  select(feature_number, Organism, Replicate, log2_change, `row m/z`, NOSC, logxic)%>%
  inner_join(cutoff_canopus%>%
               rename(feature_number = name)%>%
               mutate(feature_number = as.character(feature_number))%>%
               mutate_if(is.numeric, angular_transform), by = 'feature_number')%>%
  left_join(sirius_zodiac_anotations%>%
              select(feature_number, ZodiacScore)%>%
              mutate(feature_number = as.character(feature_number)), by = 'feature_number')%>%
  filter(ZodiacScore >= .98)%>%
  select(-ZodiacScore)%>%
  unite(row, c('feature_number', 'Organism', 'Replicate'), sep = "_")%>%
  column_to_rownames('row')

recipe <- canopus_mulreg%>%
  recipe(log2_change ~ .)%>%
  step_center(all_predictors())%>%
  step_scale(all_predictors())


##canopus training set testing

canopus_split <- rsample::initial_split(canopus_mulreg, 9/10)

canopus_training <- rsample::training(canopus_split)
canopus_test <- rsample::testing(canopus_split)

canopus_x <- canopus_training[2:ncol(canopus_training)]%>%
  as.matrix()

canopus_y <- canopus_training['log2_change']%>%
  as.matrix()

canopus_test_x <-  canopus_test[2:ncol(canopus_test)]%>%
  as.matrix()

canopus_test_y <- canopus_test['log2_change']%>%
  as.matrix()

#running 
glmnet_canopus <- cv.glmnet(canopus_x, canopus_y, parallel = TRUE)

lambda_se <- glmnet_canopus$lambda.1se

#running model with parameters after lasso
canopus_model <- glmnet(canopus_x, canopus_y, intercept = TRUE, lambda = lambda_se, standarize = TRUE)

canopus_coeff <- coef(canopus_model)%>%
  as.matrix()%>%
  as.data.frame()%>%
  rename(coefficients = 1)%>%
  rownames_to_column("variable")%>%
  filter(coefficients != 0)

intercept <- canopus_coeff[[1,2]]

predicted_training <- predict(canopus_model, newx = canopus_x)
predicted <- predict(canopus_model, newx = canopus_test_x)

rootmse <- rmse_vec(as.vector(canopus_test_y), as.vector(predicted))
rootmse_training <- rmse_vec(as.vector(canopus_y), as.vector(predicted_training))

canopus_model$dev.ratio

predicted%>%
  as.data.frame()%>%
  rename(predicted = 1)%>%
  add_column(actuals = canopus_test_y)%>%
  rownames_to_column('sample')%>%
  ggplot() +
  geom_point(aes(sample, actuals), color = 'blue') +
  geom_point(aes(sample, predicted), color = 'red')

sink('./analysis/model_rmse_r2.txt')
rootmse
rootmse_training
canopus_model$dev.ratio
sink()

canopus_weighted_lability <- feature_stats_wdf%>%
  filter(Timepoint == 'T0',
         DayNight == 'Day')%>%
  inner_join(canopus_mulreg%>%
               rownames_to_column("row")%>%
               separate(row, c('feature_number', 'Organism', 'Replicate'), sep = "_")%>%
               select(-log2_change)%>%
               gather(variable, value, 4:ncol(.))%>%
               inner_join(rf_canopus_coeff, by = 'variable')%>%
               mutate(modeled_lab = value*coefficients)%>%
               select(feature_number, Organism, Replicate, modeled_lab)%>%
               group_by(feature_number, Organism, Replicate)%>%
               summarize_if(is.numeric, sum)%>%
               ungroup()%>%
               mutate(modeled_lab = modeled_lab + rf_intercept),
             by = c('feature_number', 'Replicate', 'Organism'))%>%
  group_by(Organism, Replicate)%>%
  # mutate(weighted_lability = sum(modeled_lab))%>%
  mutate(weighted_lability = spatstat::weighted.median(modeled_lab, log10))%>%
  # mutate(weighted_lability = weighted.mean(modeled_lab, log10))%>%
  group_by(Organism)%>%
  mutate(x_err = sd(weighted_lability))%>%
  select(Organism, weighted_lability, x_err)%>%
  unique()%>%
  summarize_if(is.numeric, mean)%>%
  left_join(fcm_T0_5%>%
              filter(DayNight == 'Day')%>%
              group_by(Organism)%>%
              mutate(y_err = sd(cells_ul))%>%
              summarize_if(is.numeric, mean),
            by = c('Organism'))

# # dredge_n_select <- MuMIn::dredge(lm_test_canopus)
# # test <- AIC(lm_test_canopus)
# 

canopus_weighted_lability_stats <- feature_stats_wdf%>%
  filter(Timepoint == 'T0',
         DayNight == 'Day')%>%
  inner_join(canopus_mulreg%>%
               rownames_to_column("row")%>%
               separate(row, c('feature_number', 'Organism', 'Replicate'), sep = "_")%>%
               select(-log2_change)%>%
               gather(variable, value, 4:ncol(.))%>%
               inner_join(rf_canopus_coeff, by = 'variable')%>%
               mutate(modeled_lab = value*coefficients)%>%
               select(feature_number, Organism, modeled_lab)%>%
               group_by(feature_number, Organism)%>%
               summarize_if(is.numeric, sum)%>%
               ungroup()%>%
               mutate(modeled_lab = modeled_lab + rf_intercept),
             by = c('feature_number', 'Organism'))%>%
  group_by(Organism, Replicate)%>%
  mutate(weighted_lability = spatstat::weighted.median(modeled_lab, log10))%>%
  # mutate(weighted_lability = mean(modeled_lab))%>%
  group_by(Organism)%>%
  mutate(x_err = sd(weighted_lability))%>%
  select(Organism, weighted_lability, x_err)%>%
  unique()%>%
  left_join(fcm_T0_5%>%
              filter(DayNight == 'Day')%>%
              group_by(Organism)%>%
              mutate(y_err = sd(cells_ul))%>%
              summarize_if(is.numeric, mean),
            by = c('Organism'))

#Linear model
weight_lability_lm <- canopus_weighted_lability_stats%>%
  lm(cells_ul ~ weighted_lability, data = .)

wlab_p <- (weight_lability_lm%>% 
             tidy()%>% 
             filter(term == 'weighted_lability'))$p.value

wlab_f <- (weight_lability_lm%>% 
             tidy()%>% 
             filter(term == 'weighted_lability'))$statistic

wlab_slope <- weight_lability_lm$coefficients["weighted_lability"]
wlab_intercept <- weight_lability_lm$coefficients["(Intercept)"]


wlab_r2 <- summary(weight_lability_lm)$adj.r.squared

pdf('~/Documents/GitHub/DORCIERR/data/plots/canopus_weighted_lability.pdf', width = 12, height = 10)
canopus_weighted_lability%>%
  ggplot(aes(weighted_lability, cells_ul)) +
  geom_point(aes(color = Organism), stat = 'identity', size = 5) +
  geom_errorbar(aes(ymin = cells_ul - y_err, ymax = cells_ul + y_err)) +
  geom_errorbarh(aes(xmin = weighted_lability - x_err, xmax = weighted_lability + x_err))+
  scale_color_manual(values = org_colors_no_water) +
  geom_smooth(method = 'lm') +
  labs(y = bquote(Specific ~growth ~rate ~('Cells'~µL^-1 ~hr^-1)), x = "Metabolite pool modeled lability") +
  # geom_text(aes(x = -1.3, y = 0.09,
  #               label = paste("p-value: ", wlab_p%>%
  #                               formatC(format = "e", digits = 2), sep = "")), size = 9) +
  # geom_text(aes(x = -1.3, y = 0.086,
  #               label = paste("F statistic: ", wlab_f%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = -1.3, y = 0.082,
  #               label = paste("r²: ", wlab_r2%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = -1.3, y = 0.078,
  #               label = paste("Cells µL^-1", " = ", wlab_slope%>%
#                               round(digits = 2), "*lability + ", wlab_intercept%>%
#                               round(digits = 2), sep = "")), size = 9) +
theme(
  plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
  axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
  axis.text.y = element_text(size = 20),
  axis.title = element_text(size = 20),
  panel.background = element_rect(fill = "transparent"), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  panel.grid.major = element_line('grey'), # get rid of major grid
  panel.grid.minor = element_line('grey'), # get rid of minor grid
  legend.background = element_rect(fill = "transparent"), # get rid of legend bg
  legend.box.background = element_rect(fill = "transparent"),
  legend.position = 'top'
)
dev.off()

#checking gaussian
pdf("./plots/canopus_gaussian_check.pdf")
gaussian_check <- prep(recipe)%>%
  juice()%>%
  gather(variable, value)%>%
  inner_join(canopus_coeff%>% select(variable),
             by = 'variable')%>%
  group_by(variable)%>%
  nest()%>%
  mutate(non_transformed = map(data, ~ .x$value%>%
                                 na.omit()%>%
                                 car::qqPlot(ylab = variable, xlab = "Normal quantiles",
                                             main = paste(variable, " Centered-Scaled"))))
dev.off()


# VIZUALIZATIONS -- CANOPUS Random Forest ---------------------------------
names(canopus_training) <- make.names(names(canopus_training))

rf <- foreach(ntree = 100, .combine=randomForest::combine, .multicombine=TRUE, .packages='randomForest') %dopar% {
  randomForest(log2_change ~ ., canopus_training, importance = TRUE, proximity = TRUE, ntree = ntree)
}

mda_canopus <- rf$importance%>%
  as.data.frame()%>%
  rownames_to_column("feature")

# selecting variables to use in model
mean_se <- (mda_canopus%>%
              rename(inc_mse = `%IncMSE`)%>%
              filter(inc_mse >= mean(inc_mse) + sd(inc_mse)/sqrt(length(.))))$feature%>% 
  as.vector()

rf_xline <- length(mean_se)

ggplot(mda_canopus, aes(x= reorder(feature, -`%IncMSE`), y = `%IncMSE`)) +
  geom_point(stat = "identity") +
  geom_vline(xintercept = rf_xline + 0.5, color = 'red')



#model
training_rf_variables <- canopus_training%>%
  select('log2_change', all_of(mean_se))

test_rf_variable <- canopus_test

rf_can_model <- lm(log2_change ~ ., data = training_rf_variables)

rf_canopus_coeff <- coef(rf_can_model)%>%
  as.matrix()%>%
  as.data.frame()%>%
  rename(coefficients = 1)%>%
  rownames_to_column("variable")%>%
  filter(coefficients != 0)

rf_intercept = rf_canopus_coeff[[1,2]]

# VIZUALIZATIONS -- FCM ---------------------------------------------------
fcm_graphing <- fcm_wdf%>%
  filter(!Organism == 'Influent',
         !Organism == 'Offshore',
         DayNight == "Day")%>%
  mutate(Hours = case_when(Timepoint == "T0" ~ 0,
                           Timepoint == "T1" ~ 4,
                           Timepoint == "T2" ~ 13,
                           Timepoint == "T3" ~ 17,
                           Timepoint == "T4" ~ 23,
                           Timepoint == "T5" ~ 28,
                           Timepoint == "T6" ~ 37,
                           Timepoint == "TF" ~ 48))%>%
  group_by(Organism, Timepoint)%>%
  mutate(st_err = sd(`Cells µL-1`))%>%
  summarize_if(is.numeric, mean)

pdf("~/Documents/GitHub/DORCIERR/data/plots/FCM_day.pdf", width = 10, height = 10)
fcm_graphing%>%
  ggplot(aes(x= Hours, y = `Cells µL-1`, color = Organism))+
  geom_point(stat = "identity", size = 5) +
  geom_errorbar(aes(ymin = `Cells µL-1` - st_err, ymax = `Cells µL-1` + st_err)) +
  geom_line(aes(group = Organism)) +
  scale_color_manual(values = c(org_colors_no_water, "#3B9AB2")) +
  labs(y = bquote(Cells ~µL^-1)) +
  # facet_wrap(~ DayNight) +
  # scale_color_manual(values = c("darkorchid3", "#50A45C", "#AF814B", "#5BBCD6")) +
  scale_y_continuous(limits = c(0,900), breaks= seq(0, 900, 100)) +
  scale_x_continuous(limits = c(0,50), breaks = seq(0, 50, 5)) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.x = element_blank(), # get rid of major grid
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.position = 'top',
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )
dev.off()

fcm_err <- fcm_T0_5%>%
  group_by(Organism, DayNight)%>%
  mutate(final_err = sd(final_cells),
         cells_err = sd(cells_ul))%>%
  summarize_if(is.numeric, mean)%>%
  filter(DayNight == 'Day',
         Organism != 'Offshore',
         Organism != 'Influent')%>%
  mutate(Replicate = 1)%>%
  # select(-c(final_cells, cells_ul))%>%
  mutate(Replicate = as.character(Replicate))

=======
fcm_err <- fcm_T0_5%>%
  group_by(Organism, DayNight)%>%
  mutate(final_err = sd(final_cells),
         cells_err = sd(cells_ul))%>%
  summarize_if(is.numeric, mean)%>%
  filter(DayNight == 'Day',
         Organism != 'Offshore',
         Organism != 'Influent')%>%
  mutate(Replicate = 1)%>%
  # select(-c(final_cells, cells_ul))%>%
  mutate(Replicate = as.character(Replicate))

>>>>>>> 866b30c22f5339816f0bda22d3e41936737b507f
organism_order_sgr <- as.factor(fcm_T0_5$Organism)%>%
  relevel(c('Turf', 'Dictyota', 'Pocillopora Verrucosa', 'CCA'))%>%
  levels()%>%
  as.vector()

sgr_colors <- c("#3B9AB2", "#9900FF", "#A30029", "#FF850A", "#669900", "#33CC33")

pdf("~/Documents/GitHub/DORCIERR/data/plots/SpecificGrowthRate.pdf", width = 7, height = 7)
fcm_T0_5%>%
  filter(DayNight == 'Day',
         Organism != 'Offshore',
         Organism != 'Influent')%>%
  # left_join(fcm_err, 
  #           by = c('Organism', 'Replicate'))%>%
  mutate(Organism = factor(Organism, levels = c('Water control', 'Porites lobata', 'CCA', 'Pocillopora verrucosa', 'Dictyota', 'Turf')))%>%
  ggplot(aes(cells_ul, final_cells, color = Organism)) +
  # geom_boxplot() +
  geom_point(aes(alpha = 0.2), stat = 'identity', size = 7) +
  geom_errorbarh(data = fcm_err, aes(xmin = cells_ul - cells_err, xmax = cells_ul + cells_err)) +
  geom_errorbar(data = fcm_err, aes(ymin = final_cells - final_err, ymax = final_cells + final_err)) +
  # coord_flip() +
  scale_color_manual(values = c(sgr_colors, sgr_colors)) +
  labs(x = bquote(Specific ~Growth ~Rate ~(Cells ~µL^-1 ~hr^-1)), y = bquote(Maximmum ~microbial ~load ~(cells ~µL^-1))) +
  scale_y_continuous(limits = c(0,900), breaks= seq(0, 900, 100)) +
  xlim(0.02, 0.09) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.x = element_blank(), # get rid of major grid
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.position = 'top',
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )
dev.off()

# Depletolite Networks ----------------------------------------------------
molnet_class_network <- molnet_class%>% 
  left_join(networking%>% 
              select(feature_number, network), by = 'feature_number')%>% 
  select(-feature_number)%>% 
  unique()%>%
  filter(network != -1)

deplete_nets <- mul_reg%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  ungroup()%>%
  filter(activity == 'depletolite')%>%
  select(-c(Replicate:TF, Experiment, log_nc:log_shc, ra, ra_mean))%>%
  group_by(feature_number, net_act, Organism, activity)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  left_join(canopus_mulreg%>%
              rownames_to_column("row")%>%
              separate(row, c('feature_number', 'Organism'), sep = "_")%>%
              select(-log2_change)%>%
              gather(variable, value, 3:ncol(.))%>%
              inner_join(canopus_coeff, by = 'variable')%>%
              mutate(modeled_lab = value*coefficients)%>%
              select(feature_number, Organism, modeled_lab)%>%
              group_by(feature_number, Organism)%>%
              summarize_if(is.numeric, sum)%>%
              ungroup()%>%
              mutate(modeled_lab = modeled_lab + intercept),
            by = c('Organism', 'feature_number'))%>%
  group_by(Organism, network, activity)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
<<<<<<< HEAD
  ungroup()%>%
  left_join(ecoNet%>%
              select(-scan)%>%
              filter(network != '-1')%>%
              unique(), by = 'network')%>%
  # left_join(molnet_class_network, by = 'network')%>%
  select(Organism, network, ecoNetConsensus, ecoNetConsensusLevel, ecoNetConsensusScore, everything())

write_csv(deplete_nets, './analysis/deplete_networks.csv')  

identify_depletes <- deplete_nets%>%
  select(Organism, network)%>%
  left_join(metadata%>%
              select(feature_number, network, combined_ID, binary_ID),
            by = 'network')
# filter(binary_ID != 3,
#        network != -1,
#        binary_ID == 1)


write_csv(identify_depletes, './analysis/identify_nets_libid.csv')

top_nets <- (deplete_nets%>%
               filter(log2_change < -5.64)%>% # log2(50)
               # group_by(Organism)%>%
               # nest()%>%
               # mutate(data = map(data, ~ top_n(.x, 5, -log2_change)%>%
               #                     select(network, log2_change)))%>%
               # unnest(data)%>%
               # ungroup()%>%
               select(network)%>%
               unique())$network%>%
  as.vector()

dendogram_df <- log2_change_vals%>%
  filter(DayNight == 'Day',
         Organism != 'Water control')%>%
  select(feature_number, Organism, Replicate, log2_change)%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  filter(network %in% top_nets)%>%
  group_by(Organism, Replicate, network)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(network)%>%
  mutate(zscore = zscore(-log2_change))%>%
  ungroup()%>%
  select(c(Organism, Replicate, network, zscore))%>%
  unite(sample, c('Organism', 'Replicate'), sep = ' ')%>%
  spread(network, zscore)

write_csv(dendogram_df, './analysis/hc_depletolites_df.csv')  

# VIZUALIZATIONS -- Cytoscape ---------------------------------------------
filtered_features <- feature_stats_wdf%>%
  ungroup()%>%
  select(feature_number, Organism, DayNight)%>%
  unique()%>%
  group_by(feature_number, DayNight)%>%
  summarize_all(paste0, collapse = ", ")%>%
  mutate(post_filtering = 1)

cyto_deplete <- dorc_wdf%>%
  filter(!Organism == "Influent",
         !Organism == "Offshore",
         !Timepoint == c("T1", "T2", "T3", "T4"))%>%
  select(feature_number, Organism:ncol(.))%>%
  left_join(log2_change_vals%>%
              select(feature_number, Organism, DayNight, Replicate, log2_change)%>%
              mutate(Replicate = as.character(Replicate)),
            by = c('feature_number', 'Organism', 'DayNight', 'Replicate'))%>%
  filter(Timepoint == 'T0')%>%
  group_by(feature_number, Organism, DayNight)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  group_by(feature_number, DayNight)%>%
  mutate(log2_sizing = min(log2_change))%>%
  ungroup()%>%
  gather(response_var, val, xic:log2_sizing)%>%
  unite(org_var, c(Organism, response_var), sep = '_')%>%
  spread(org_var, val)%>%
  left_join(filtered_features,
            by = c('feature_number', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  rename(net_activity = activity)%>%
  filter(DayNight == 'Day')%>%
  left_join(ecoNet%>%
              select(-scan)%>%
              filter(network != '-1')%>%
              unique(), by = 'network')
  

=======
  ungroup()%>%
  left_join(molnet_class_network, by = 'network')%>%
  select(Organism, network, molnet_string, modeled_lab, everything())

write_csv(deplete_nets, './analysis/deplete_networks.csv')  

identify_depletes <- deplete_nets%>%
  select(Organism, network)%>%
  left_join(metadata%>%
              select(feature_number, network, combined_ID, binary_ID),
            by = 'network')
# filter(binary_ID != 3,
#        network != -1,
#        binary_ID == 1)


write_csv(identify_depletes, './analysis/identify_nets_libid.csv')

top_nets <- (deplete_nets%>%
               filter(log2_change < -5.64)%>% # log2(50)
               # group_by(Organism)%>%
               # nest()%>%
               # mutate(data = map(data, ~ top_n(.x, 5, -log2_change)%>%
               #                     select(network, log2_change)))%>%
               # unnest(data)%>%
               # ungroup()%>%
               select(network)%>%
               unique())$network%>%
  as.vector()

dendogram_df <- log2_change_vals%>%
  filter(DayNight == 'Day',
         Organism != 'Water control')%>%
  select(feature_number, Organism, Replicate, log2_change)%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  filter(network %in% top_nets)%>%
  group_by(Organism, Replicate, network)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(network)%>%
  mutate(zscore = zscore(-log2_change))%>%
  ungroup()%>%
  select(c(Organism, Replicate, network, zscore))%>%
  unite(sample, c('Organism', 'Replicate'), sep = ' ')%>%
  spread(network, zscore)

write_csv(dendogram_df, './analysis/hc_depletolites_df.csv')  

# VIZUALIZATIONS -- Cytoscape ---------------------------------------------
filtered_features <- feature_stats_wdf%>%
  ungroup()%>%
  select(feature_number, Organism, DayNight)%>%
  unique()%>%
  group_by(feature_number, DayNight)%>%
  summarize_all(paste0, collapse = ", ")%>%
  mutate(post_filtering = 1)

cyto_deplete <- dorc_wdf%>%
  filter(!Organism == "Influent",
         !Organism == "Offshore",
         !Timepoint == c("T1", "T2", "T3", "T4"))%>%
  select(feature_number, Organism:ncol(.))%>%
  left_join(log2_change_vals%>%
              select(feature_number, Organism, DayNight, Replicate, log2_change)%>%
              mutate(Replicate = as.character(Replicate)),
            by = c('feature_number', 'Organism', 'DayNight', 'Replicate'))%>%
  filter(Timepoint == 'T0')%>%
  group_by(feature_number, Organism, DayNight)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  group_by(feature_number, DayNight)%>%
  mutate(log2_sizing = min(log2_change))%>%
  ungroup()%>%
  gather(response_var, val, xic:log2_sizing)%>%
  unite(org_var, c(Organism, response_var), sep = '_')%>%
  spread(org_var, val)%>%
  left_join(filtered_features,
            by = c('feature_number', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  rename(net_activity = activity)%>%
  filter(DayNight == 'Day')%>%
  left_join(ecoNet%>%
              select(-scan)%>%
              filter(network != '-1')%>%
              unique(), by = 'network')
  

>>>>>>> 866b30c22f5339816f0bda22d3e41936737b507f
write_csv(cyto_deplete, './analysis/cyto_depletes.csv')


# Linda df ----------------------------------------------------------------
linda_df <- feature_stats_wdf%>%
  ungroup()%>%
  filter(DayNight == "Day")%>%
  select(feature_number, Organism, Timepoint, Replicate, DayNight, xic)%>%
  group_by(feature_number, Organism, DayNight, Replicate, Timepoint)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  spread(Timepoint, xic)%>%
  left_join(log2_change_vals%>%
              select(-c(T0, TF))%>%
              group_by(feature_number, Organism, DayNight)%>%
              summarize_if(is.numeric, mean),
            by = c("feature_number", "Organism","DayNight"))%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  select(feature_number, network, activity, Organism, Replicate, DayNight, everything())%>%
  select(-net_act)

linda_sum <- linda_df%>%
  group_by(Organism, activity, Replicate, DayNight)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  select(-log2_change)%>%
  gather(Timepoint, xic, T0:TF)%>%
  group_by(Organism, activity, Timepoint, DayNight)%>%
  mutate(err = sd(xic))%>%
  # summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  unite(org_activity, c('Organism', 'activity'), sep = ' / ', remove = FALSE)%>%
  mutate(shape = case_when(Timepoint == 'T0' ~ 'T0',
                           Timepoint == 'TF' & activity == 'depletolite' ~ 'TF depletion',
                           Timepoint == 'TF' & activity == 'accumolite' ~ 'TF accumulation',
                           Timepoint == 'TF' & activity == 'recalcitrant' & Organism == 'Turf' ~ 'accumulation',
                           TRUE ~ 'TF depletion'))%>%
  group_by(Organism, Replicate, Timepoint)%>%
  mutate(ra = xic/sum(xic, na.rm = TRUE))%>%
  ungroup()


# png("./plots/change_xic_activity.png", width = 1500, height = 1000)
# linda_sum%>%
#   filter(activity != 'recalcitrant')%>%
#   mutate(log10 = log10(xic))%>%
#   ggplot(aes(ra, org_activity, color = Organism, shape = shape)) +
#   geom_point(stat = 'identity', size = 14, alpha = 0.45) +
#   # geom_errorbarh(aes(xmin = xic - err, xmax = xic + err)) +
#   geom_line(aes(group = org_activity)) +
#   scale_color_manual(values = org_colors_no_water) +
#   scale_shape_manual(values = c("\u25A0", "\u25C4", "\u25BA")) +
#   labs(x = 'Relative Abundance', y = 'Organism / Activity') +
#   # scale_x_log10() +
#   facet_wrap(~activity, nrow = 1, ) +
#   theme(axis.line = element_blank(),
#         axis.text = element_text(size = 20),
#         axis.ticks = element_blank(),
#         legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#         legend.box.background = element_blank(),
#         panel.background = element_rect(fill = "transparent"), # bg of the panel
#         plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#         strip.text = element_text(size=20),
#         legend.text = element_text(size = 20),
#         legend.title = element_text(size = 20),
#         legend.position = 'none',
#         title = element_text(size = 40, hjust = 0.5),
#         plot.title = element_text(hjust = 0.5, vjust = 3),
#         strip.background = element_blank())

linda_mean_diff <- linda_sum%>%
  select(-c(err:ra))%>%
  group_by(org_activity, Organism, activity)%>%
  mutate(x_val = max(xic, na.rm = TRUE))%>%
  ungroup()%>%
  group_by(org_activity, Organism, activity, Timepoint)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  spread(Timepoint, xic)%>%
  mutate(difference = TF-T0)%>%
  filter(activity != 'recalcitrant')

png("./plots/change_xic_activity.png", width = 2000, height = 1500)
linda_sum%>%
  filter(activity != 'recalcitrant')%>%
  mutate(log10 = log10(xic))%>%
  ggplot(aes(xic, Organism, color = Organism, shape = shape)) +
  geom_point(stat = 'identity', size = 14, alpha = 0.45) +
  # geom_errorbarh(aes(xmin = xic - err, xmax = xic + err)) +
  geom_line(aes(group = Organism)) +
  geom_text(data = linda_mean_diff, aes(x = x_val, y = Organism, 
                                        label = formatC(difference, format = 'e', digits = 2), 
                                        shape = NULL), vjust = -1, size = 7) +
  facet_wrap(~activity, nrow = 1, scales = 'free_x') +
  scale_color_manual(values = org_colors_no_water) +
  scale_shape_manual(values = c("\u25A0", "\u25BA", "\u25C4")) +
  labs(x = 'Sum Intensity (xic)', y = 'Organism', color = 'Organism: ', shape = 'Sample:') +
  # xlim(0,0.3) +
  # scale_x_log10() +
  # facet_wrap(~activity, nrow = 3) +
  theme(axis.line = element_blank(),
        axis.text = element_text(size = 20),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        strip.text = element_text(size= 25),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.position = 'top',
        title = element_text(size = 40, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, vjust = 3),
        strip.background = element_blank()) +
  guides(color = guide_legend(nrow=2, byrow=TRUE),
         shape = guide_legend(nrow=2, byrow=TRUE))
dev.off()

# linda_split <- linda_df%>%
#   group_by(Organism)%>%
#   split(.$Organism)

# list2env(linda_split, globalenv())
# 
# write_csv(`Pocillopora verrucosa`, '~/Documents/SDSU_Scripps/DORCIERR/Datasets/linda_features_pocillopora.csv')
# write_csv(`Porites lobata`, '~/Documents/SDSU_Scripps/DORCIERR/Datasets/linda_features_porites.csv')
# write_csv(`CCA`, '~/Documents/SDSU_Scripps/DORCIERR/Datasets/linda_features_cca.csv')
# write_csv(`Turf`, '~/Documents/SDSU_Scripps/DORCIERR/Datasets/linda_features_turf.csv')
# write_csv(`Dictyota`, '~/Documents/SDSU_Scripps/DORCIERR/Datasets/linda_features_dictyota.csv')

linda_pie_features <- org_log2_ra%>%
  filter(DayNight == 'Day')%>%
  group_by(feature_number, Organism, network, activity)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  left_join(metadata%>%
              select(combined_ID, binary_ID, feature_number, ZodiacMF, ZodiacScore, C, H, N, O, P, S, NOSC, dG),
            by = 'feature_number')

write_csv(linda_pie_features, '~/Documents/GitHub/DORCIERR/data/analysis/lindaPieChartMetadata.csv')  

linda_depletolites <- venn%>%
  left_join(metadata%>%
              select(combined_ID, binary_ID, SmilesLibrary_, SmilesAnalog_, feature_number, network, ZodiacMF, ZodiacScore, C, H, N, O, P, S, NOSC, dG),
            by = 'feature_number')%>%
  left_join(venn_df%>%
              select(Intersection, feature_number),
            by = 'feature_number')

write_csv(linda_depletolites, '~/Documents/GitHub/DORCIERR/data/analysis/lindaDepletolites.csv')
