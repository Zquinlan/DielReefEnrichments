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
# CORES -- setting processors available -----------------------------------
##Only used if future_mapping
# num_cores <- availableCores() -1
# don't murder your compututer and save your self a core
# this is the parellel planning step (changes global env so this is plan for all parellel 
# work unless specificed otherwise)
# plan(multiprocess, workers = num_cores) #defaults to sequential process, multiprocess is one option for parellel

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
molnet_class <- read_csv("analysis/Moorea2017_MolNetEnhancer.csv")%>%
  rename(feature_number = "cluster index")%>%
  mutate(feature_number = as.character(feature_number))%>%
  select(feature_number, CF_kingdom, CF_class, CF_subclass, CF_superclass)%>%
  unite(molnet_string, c(CF_kingdom, CF_superclass, CF_class, CF_subclass), sep = ";")

#PPl extraction efficiency
extraction_efficiency <- read_csv("./raw/metabolomics/Extraction_efficiency.csv")%>%
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

# SET -- CANOPUS filters --------------------------------------------------
# Canopus annotations which are above level 3 AND 80% probability

canopus_filtered_tidy <- canopus_chemonnt_tidy%>%
  rename('feature_number' = 'name')%>%
  group_by(feature_number)%>%
  nest()%>%
  mutate(data = map(data, ~filter(.x, canopus_probability >= 0.80)%>%
                      filter(level > 3)%>%
                      # filter(., level == max(level))%>%
                      filter(canopus_probability == max(canopus_probability))%>%
                      filter(level == max(level))%>%
                      filter(nchar(CLASS_STRING) == max(nchar(CLASS_STRING)))%>%
                      filter(nchar(canopus_annotation) == max(nchar(canopus_annotation)))))%>%
  unnest(data)

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
                                TRUE ~ "Night"),
           xic = case_when(xic == 0 ~ 1000,
                           TRUE ~ as.numeric(xic)))%>%
    group_by(Organism, Timepoint, DayNight, feature_number)%>%
    mutate(Replicate = as.numeric(Replicate),
           xic = case_when(sum(xic) == 4000 ~ xic + Replicate, 
                           sum(xic) == 2000 ~ xic + Replicate,
                           TRUE ~ as.numeric(xic)))%>%
    ungroup()%>%
    spread(Timepoint, xic)%>%
    group_by(Organism, DayNight, feature_number)%>%
    summarize_if(is.numeric, max, na.rm = TRUE)%>%
    mutate(log2_change = log2(TF/T0))%>%
    filter(log2_change >= 1 | log2_change <= -1)%>%
    ungroup()
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
                              TRUE ~ "Night"),
         xic = case_when(xic == 0 ~ 1000,
                         TRUE ~ as.numeric(xic)))%>%
  group_by(Organism, Timepoint, DayNight, feature_number)%>%
  mutate(Replicate = as.numeric(Replicate),
         xic = case_when(sum(xic) == 4000 ~ xic + Replicate, 
                         sum(xic) == 2000 ~ xic + Replicate,
                         TRUE ~ as.numeric(xic)))%>%
  ungroup()%>%
  spread(Timepoint, xic)%>%
  group_by(Organism, DayNight, feature_number)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         log2_change = log2(TF/T0),
         complete_removal = case_when(mean(TF) > 0 & mean(TF, na.rm = TRUE) == 1000 ~ "removed",
                                      mean(TF, na.rm = TRUE) > mean(T0, na.rm = TRUE) ~"accumolite",
                                      TRUE ~ "semi-removed"))%>%
  ungroup()


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
                              TRUE ~ "Night"),
         xic = case_when(xic == 0 ~ 1000,
                         TRUE ~ as.numeric(xic)))%>%
  group_by(Organism, Timepoint, DayNight, feature_number)%>%
  mutate(Replicate = as.numeric(Replicate),
         xic = case_when(sum(xic) == 4000 ~ xic + Replicate, 
                         sum(xic) == 2000 ~ xic + Replicate,
                         TRUE ~ as.numeric(xic)))%>%
  ungroup()%>%
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
  filter(DayNight == 'Day')%>%
  select(feature_number, Organism)%>%
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

# RELATIVIZATION AND NORMALIZATION -- xic_log10 -----------------
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
  mutate(xic = case_when(xic == 0 ~ 1000,
                         TRUE ~ as.numeric(xic)))%>%
  group_by(Organism, Timepoint, DayNight, feature_number)%>%
  mutate(Replicate = as.numeric(Replicate),
         xic = case_when(sum(xic) == 4000 ~ xic + Replicate,
                         sum(xic) == 2000 ~ xic + Replicate,
                         TRUE ~ as.numeric(xic)))%>%
  ungroup()%>%
  select(-c("Experiment", "Organism", "Replicate", "Timepoint", "DayNight"))%>%
  group_by(sample_name)%>%
  mutate(log10 = log10(xic),
         ra = xic/sum(xic),
         asin = asin(ra),
         feature_number = as.character(feature_number))%>%
  gather(transformation, values, xic:asin)%>%
  arrange(transformation)%>%
  unite(sample_transformed, c("sample_name", "transformation"), sep = "_")%>%
  spread(sample_transformed, values)

# DORCIERR feature_table --------------------------------------------------
feature_table_combined <- right_join(metadata, feature_table_relnorm, by = "feature_number")

dorcierr_table_wdf_temp <- feature_table_combined%>%
  dplyr::select(c(feature_number, everything()))

# CLEANING-- adding carbon normalized values to wdf ------------------
  carbon_normalized_xic_NOSC <- dorcierr_table_wdf_temp%>%
  filter(`characterization scores` == "Good")%>%
  dplyr::select(c(feature_number, C, ends_with("_xic")))%>%
  gather(sample_name, xic, 3:ncol(.))%>%
  group_by(sample_name)%>%
  mutate(ra = xic/sum(xic, na.rm = TRUE),
         percent_total_C = xic*C)%>%
  mutate(sum_c = sum(percent_total_C,  na.rm = TRUE),
         carbon_norm_temp = percent_total_C/sum_c)%>%
  ungroup()%>%
  right_join(metadata%>%
               select(c(feature_number, NOSC)),
             .,  by = "feature_number")%>%
  mutate(carbon_normalized_NOSC = carbon_norm_temp*NOSC,
         sample_name = gsub("_xic", "", sample_name))%>%
  select(c(feature_number, sample_name, carbon_normalized_NOSC))%>%
  separate(sample_name, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
  mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                       TRUE ~ as.character(Experiment)))%>%
  mutate(Organism = case_when(Organism == "CC" ~ "CCA",
                                     Organism == "DT" ~ "Dictyota",
                                     Organism == "PL" ~ "Porites lobata",
                                     Organism == "PV" ~ "Pocillopora verrucosa",
                                     Organism == "TR" ~ "Turf",
                                     Organism == "WA" ~ "Water control",
                                     TRUE ~ as.character(Organism)))%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))%>%
  group_by(feature_number, Organism, Timepoint, DayNight)%>%
  summarize_if(is.numeric, mean)%>%
  filter(!carbon_normalized_NOSC < -0.001)

dorcierr_features_wdf <-dorcierr_table_wdf_temp

write_csv(dorcierr_features_wdf, "~/Documents/GitHub/DORCIERR/data/analysis/Dorcierr_feature_table_master_post_filtered.csv")

# PRE-CLEANING -- Making Dorcierr working data frame for stats -----------------------------
dorc_transposed <- dorcierr_features_wdf%>%
  dplyr::select(c(feature_number, ends_with("_log10")))%>%
  gather(sample_ID, log, 2:ncol(.))%>%
  spread(feature_number, log)

dorc_transposed$sample_ID <- dorc_transposed$sample_ID%>%
  gsub("_log10", "", .)%>%
  gsub("-", "_", .)

blanks_wdf <- dorc_transposed%>%
  filter(sample_ID %like% "%Blank%",
         !sample_ID == 'D_Blank_DI')%>%
  mutate(sample_ID = case_when(sample_ID == "Blank_Lot_6350565_01"~ "Blank_Blank_635056501",
                               sample_ID == "Blank_SD_01_A" ~ "Blank_Blank_SD01A",
                               sample_ID == "Blank_SD_01_B" ~ "Blank_Blank_SD01B",
                               sample_ID == "Blank_SD_LoRDI" ~ "Blank_Blank_SDLoRDI",
                               sample_ID == "Blank_SD_PPL" ~ "Blank_Blank_SDPPL",
                               sample_ID == "Blank? look up name on PPL" ~ "Blank_Blank_unknown",
                               sample_ID == "D_Blank" ~ "Blank_Blank_D",
                               TRUE ~ "Blank_Blank_Spiffy"))%>%
  separate(sample_ID, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")

dorc_wdf <- dorc_transposed%>%
  separate(sample_ID, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
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

# PRE-STATS CLEANING -- DOM-stats -----------------------------------------
dom_stats_wdf<- dorc_wdf%>%
  # full_join(., fdom_doc_log10, by = c("Organism", "Timepoint", "Replicate", "DayNight"))%>%
  filter(!Organism == "Influent",
         !Organism == "Offshore",
         !Timepoint == c("T1", "T2", "T3", "T4"))%>%
  gather(feature_number, log, 6:ncol(.))%>%
  right_join(log2_features%>%
               select(-log2_change, -Replicate), by = c("feature_number", "DayNight"))

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

# PRE-STATS CLEANING -- fDOM and DOC --------------------------------------------------------------
fdom_log10 <-fdom_wdf%>%
  dplyr::select(c(1:5, 15:21))%>%
  gather(fluor, val, 6:ncol(.))%>%
  mutate(val = log10(val))%>%
  spread(fluor, val)


doc_log10 <- moorea_doc%>%
  filter(sample_name != "D_OF_1_T0N",
         sample_name != "D_IN_2_T0N",
         sample_name != "D_PL_3_TFN",
         sample_name != "D_TR_1_T0N",
         sample_name != "D_WA_2_T0D",
         sample_name != "D_WA_1_T0D",
         sample_name != "D_CC_1_T0D",
         sample_name != "D_CC_2_T0D")%>%
  separate(sample_name, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_", remove = FALSE)%>%
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
                              TRUE ~ "Night"),
         DOC_log = log10(DOC))

fdom_doc_log10 <- left_join(fdom_log10, doc_log10, by = "sample_name")%>%
  dplyr::select(-sample_name)
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
  mutate(abundant = case_when(max(ra) > 0.01 | sum(ra > 0.001) >=3 ~ "abundant",
                              TRUE ~ "rare"))

microbe_no_rare <- microbe_combined%>%
  select(-c(sample_code, reads, numOtus))%>%
  ungroup()%>%
  spread(Timepoint, ra)%>%
  group_by(OTU, DayNight, Organism)%>%
  mutate(log2_change = log2(TF/mean(T0, na.rm = TRUE)),
         log2_change = na_if(log2_change, "Inf"), 
         log2_change = na_if(log2_change, "-Inf"),
         log2_change = na_if(log2_change, "NaN"))%>%
  filter(abundant == "abundant")%>%
  select(-abundant)

# microbe_log2 <- microbe_combined%>%
#   select(-c(cell_abun, log10, numOtus, reads, sum, `Cells µL-1`))%>%
#   spread(Timepoint, ra)%>%
#   group_by(OTU, DayNight, Organism)%>%
#   mutate(log2_change = log2(TF/mean(T0, na.rm = TRUE)))%>%
#   filter(abundant == "abundant")%>%
#   select(-abundant)
  

# PRE-STATS CLEANING -- microbe RA data  --------------------------------------------
ra_bigger_TF <- microbe_combined%>%
  filter(abundant == "abundant")%>%
  select(-abundant)%>%
  select(-c(reads, numOtus, sample_code))%>%
  spread(Organism, ra)%>%
  gather(Organism, ra, 6:10)%>%
  mutate(difference = ra - `Water control`)%>%
  filter(difference > 0)%>%
  dplyr::select(c(DayNight, OTU, Organism))

average_ra <- microbe_no_rare%>%
  select(-c(Experiment, Replicate))%>%
  group_by(OTU, Organism, DayNight)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()

osm_ra_bigger_TF <- microbe_combined%>%
  filter(abundant == "abundant",
         Organism != "Turf",
         Organism != "Porites lobata")%>%
  select(-abundant)%>%
  select(-c(reads, numOtus, sample_code))%>%
  spread(Organism, ra)%>%
  gather(Organism, ra, 6:8)%>%
  mutate(difference = ra - `Water control`)%>%
  filter(difference > 0)%>%
  dplyr::select(c(DayNight, OTU, Organism))


# SET SEED ----------------------------------------------------------------
set.seed(2005)


# STATS  -- T-TEST Network level ------------------------------------------
net_test <- dom_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = "feature_number")%>%
  filter(network != "-1")%>%
  group_by(network, DayNight)%>%
  nest()%>%
  mutate(greater = map(data, ~ t.test(log ~ Timepoint, .x, alternative = "greater")),
         lesser = map(data, ~ t.test(log ~ Timepoint, .x, alternative = "less")))%>%
  select(-data)%>%
  mutate(greater = map(greater, ~ .x["p.value"][[1]]))%>%
  mutate(lesser = map(lesser, ~ .x["p.value"][[1]]))%>%
  ungroup()%>%
  mutate(greater = as.numeric(greater),
         lesser = as.numeric(lesser),
         FDR_greater = p.adjust(greater, method = "BH"),
         FDR_lesser = p.adjust(lesser, method = "BH"))%>%
  mutate(activity = case_when(FDR_greater < 0.05 ~ "depletolite",
                              FDR_lesser < 0.05 ~ "accumolite",
                              TRUE ~ "recalcitrant"))

net_activity <- net_test%>%
  group_by(network, activity)%>%
  summarize_if(is.numeric, mean)%>%
  select(network, activity)

# STATS - GLM Activity groupings -----------------------------------------------
#Difference within activity groupings
net_glm <- log2_change_vals%>%
  left_join(networking%>%
              select(c(feature_number, network, NOSC)), by = "feature_number")%>%
  filter(network != "-1")%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  left_join(net_activity, by = 'network')%>%
  inner_join(day_exometabolites, by = 'feature_number')%>%
  mutate(activity = case_when(is.na(activity) ~ 'recalcitrant',
                              TRUE ~ as.character(activity)))%>% #This line was added because entire networks did not pass the log2 bottleneck
  select(-c(T0,TF, complete_removal))%>%
  gather(response_var, value, NOSC:`row m/z`)%>%
  group_by(DayNight, activity, response_var)%>%
  nest()%>%
  mutate(model = map(data, ~ glm(log2_change ~ value, family = gaussian, .x)),
         p_vals = map(model, ~tidy(.x)%>%
                        filter(term == 'value')),
         r2 = map(model, ~with(summary(.x), 1-deviance/null.deviance)))%>%
  select(-c(data, model))%>%
  unnest(p_vals)%>%
  ungroup()%>%
  mutate(r2 = as.numeric(r2),
         r = sqrt(r2),
         FDR = p.adjust(p.value, method = 'BH'))

write_csv(net_glm, "./analysis/net_glm_07012020.csv")


org_glm <- log2_change_vals%>%
  left_join(networking%>%
              select(c(feature_number, network, NOSC)), by = "feature_number")%>%
  filter(network != "-1")%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  left_join(net_activity, by = 'network')%>%
  right_join(day_exometabolites, by = 'feature_number')%>%
  mutate(activity = case_when(is.na(activity) ~ 'recalcitrant',
                              TRUE ~ as.character(activity)))%>% #This line was added because entire networks did not pass the log2 bottleneck
  filter(activity == 'depletolite')%>%
  select(-c(T0,TF, complete_removal, activity))%>%
  gather(response_var, value, NOSC:`row m/z`)%>%
  group_by(DayNight, Organism, response_var)%>%
  nest()%>%
  mutate(model = map(data, ~ glm(log2_change ~ value, family = gaussian, .x)),
         p_vals = map(model, ~tidy(.x)%>%
                        filter(term == 'value')),
         r2 = map(model, ~with(summary(.x), 1-deviance/null.deviance)))%>%
  select(-c(data, model))%>%
  unnest(p_vals)%>%
  ungroup()%>%
  mutate(r2 = as.numeric(r2),
         r = sqrt(r2),
         FDR = p.adjust(p.value, method = 'BH'))


pdf("plots/glm_NOSC_Mass.pdf", width = 15, height = 10)
log2_change_vals%>%
  left_join(networking%>%
              select(feature_number, network, NOSC), by = "feature_number")%>%
  filter(network != "-1")%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  left_join(net_activity, by = 'network')%>%
  right_join(day_exometabolites, by = 'feature_number')%>%
  filter(activity != 'recalcitrant')%>%
  mutate(activity2 = activity)%>% 
  select(-c(T0,TF, complete_removal))%>%
  gather(response_var, value, NOSC:`row m/z`)%>%
  ggplot(aes(value, log2_change, color = activity)) +
  facet_wrap(~ response_var, scales = 'free_x') + 
  geom_point(stat = 'summary', fun.y = mean) +
  scale_color_manual(values = c('#78B7C5', '#EBCC2A')) +
  new_scale('color') + 
  geom_smooth(method = lm, aes(color = activity2)) +
  scale_color_manual(values = c('#3B9AB2', '#E1AF00')) +
  theme(
    # legend.position = "none",
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )
  dev.off()
  
  log2_change_vals%>%
    left_join(networking%>%
                select(feature_number, network, NOSC), by = "feature_number")%>%
    filter(network != "-1")%>%
    left_join(metadata%>%
                select(feature_number, `row m/z`), by = "feature_number")%>%
    left_join(net_activity, by = 'network')%>%
    inner_join(unique_benthic_metabolites, by = c('Organism', 'feature_number'))%>%
    filter(activity == 'depletolite')%>%
    mutate(Organism2 = Organism)%>% 
    select(-c(T0,TF, complete_removal))%>%
    gather(response_var, value, NOSC:`row m/z`)%>%
    ggplot(aes(value, log2_change, color = Organism)) +
    facet_wrap(~ response_var, scales = 'free_x') + 
    geom_point(stat = 'summary', fun.y = mean) +
    # scale_color_manual(values = c('#78B7C5', '#EBCC2A')) +
    new_scale('color') + 
    geom_smooth(method = lm, aes(color = Organism2)) +
    # scale_color_manual(values = c('#3B9AB2', '#E1AF00')) +
    theme(
      # legend.position = "none",
      # plot.margin = margin(2,.8,2,.8, "cm"),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent")
    )
  
nosc_quadrants <- log2_change_vals%>%
  left_join(networking%>%
              select(feature_number, network, NOSC), by = "feature_number")%>%
  filter(network != "-1")%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  left_join(net_activity, by = 'network')%>%
  filter(activity != 'recalcitrant')%>%
  group_by(feature_number, Organism, DayNight, network, activity)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(quadrant = case_when(log2_change > -3.3 & log2_change < 0 & NOSC > 0 ~ '6',
                              log2_change > -3.3 & log2_change < 0 & NOSC <= 0 ~ '5',
                              log2_change <= -3.3 & NOSC > 0 ~ '8',
                              log2_change <= -3.3 & NOSC <= 0 ~ '7',
                              log2_change > 0 & log2_change < 1 & NOSC <= 0 ~ '3',
                              log2_change > 0 & log2_change < 1 & NOSC > 0 ~ '4',
                              log2_change >= 1 & NOSC <= 0 ~ '1',
                              log2_change >= 1 & NOSC > 0 ~ '2'))%>%
  filter(!is.na(quadrant))
  
  

# STATS - ANOVA Activity groupings ----------------------------------------
# Difference between activity groupings
net_anova <- net_activity%>%
  left_join(networking%>%
              select(feature_number, network, NOSC), by = "network")%>%
  filter(network != "-1")%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  gather(response_var, value, NOSC:`row m/z`)%>%
  group_by(response_var)%>%
  nest()%>%
  mutate(anova = map(data, ~ aov(value ~ activity, .x)%>%
                      tidy()%>%
                      filter(!term == "Residuals")%>%
                      select(p.value)%>%
                      rename(anova = p.value)),
         tukey = map(data, ~ aov(value ~ activity, .x)%>%
                       TukeyHSD(p.adjust.methods = 'BH')%>%
                       tidy()))%>%
  select(-data)%>%
  unnest(c(anova, tukey))%>%
  ungroup()%>%
  filter(anova < 0.05,
         adj.p.value < 0.05)


net_activity%>%
  left_join(networking%>%
              select(feature_number, network, NOSC), by = "network")%>%
  filter(network != "-1")%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  gather(response_var, value, NOSC:`row m/z`)%>%
  ggplot(aes(activity, value)) +
  # geom_jitter() +
  geom_boxplot() +
  facet_wrap(~response_var, scale = 'free_y')



# STATS ANOVA -- Microbe TWO-Way ------------------------------------------
aov_microbe <- microbe_no_rare%>%
  group_by(OTU)%>%
  nest()%>%
  mutate(anova = map(data, ~ aov(log2_change ~ Organism*DayNight, .x)%>%
                       tidy()%>%
                       filter(!term == "Residuals")%>%
                       dplyr::select(term, p.value)))%>%
  dplyr::select(-data)%>%
  unnest(anova)

anova_microbe_pvalues <- aov_microbe%>%
  ungroup()%>%  
  filter(p.value < 0.05)
# add_column(FDR = p.adjust(.$p.value, method = "BH"))%>%
# filter(FDR < 0.05)

organism_significant_microbes <- as.vector(anova_microbe_pvalues%>%
                                             filter(term == "Organism"))$OTU

DayNight_significant_microbes <- as.vector(anova_microbe_pvalues%>%
                                             filter(term != "Organism"))$OTU

# STATS -- One-Sample T-test ----------------------------------------------
t_test <- dom_stats_wdf%>%
  group_by(feature_number, Organism, DayNight)%>%
  spread(Timepoint, log)%>%
  select(-c(Experiment, Replicate))%>%
  nest()%>%
  mutate(greater = map(data, ~ t.test(.$TF, mu = mean(.$T0, na.rm=TRUE), alternative = "greater")),
         lesser = map(data, ~ t.test(.$TF, mu = mean(.$T0, na.rm=TRUE), alternative = "less")))

t_pvals <- t_test%>%
  select(-data)%>%
  mutate(greater = map(greater, ~ .x["p.value"][[1]]))%>%
  mutate(lesser = map(lesser, ~ .x["p.value"][[1]]))%>%
  ungroup()%>%
  mutate(greater = as.numeric(greater),
         lesser = as.numeric(lesser),
         FDR_greater = p.adjust(greater, method = "BH"),
         FDR_lesser = p.adjust(lesser, method = "BH"))%>%
  # select(-Organism)%>%
  # group_by(feature_number, DayNight)%>%
  # nest()%>%
  # mutate(data = map(data, ~ mutate(.x, feature_number = as.character(feature_number))%>%
  #                     group_by(feature_number)%>%
  #                     summarize_if(is.numeric, min)%>%
  mutate(activity = case_when(FDR_greater < 0.05 ~ "accumolite",
                              FDR_lesser < 0.05 ~ "depletolite",
                              FDR_lesser >= 0.05 & FDR_greater >= 0.05 | 
                                is.na(FDR_lesser) & is.na(FDR_greater) ~ "recalcitrant"))

# STATS -- ANOVA metabolties Organism T0 MAJOR DEPLETEs-------------------------------------------------------
anova_dom_t0_df <- t_pvals%>%
  select(c(DayNight, feature_number, activity))%>%
  inner_join(major_deplete_features, by = c("DayNight", "feature_number"))%>%
  left_join(dom_stats_wdf, by = c("feature_number", "DayNight"))%>%
  filter(Timepoint == "T0")

anova_dom_t0 <- anova_dom_t0_df%>%
  group_by(feature_number, DayNight, activity)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(log ~ Organism, data = .x)%>%
                      tidy()%>%
                      filter(!term == "Residuals")%>%
                      dplyr::select(term, p.value)))%>%
  unnest(data)%>%
  ungroup()%>%
  mutate(FDR = p.adjust(p.value, method = "BH"),
         anova = case_when(FDR < 0.05 ~ "producer_specific",
                           !FDR < 0.05 | is.na(FDR) ~ "background"))

sigs_osm <- anova_dom_t0%>%
  filter(DayNight == "Day",
         activity == "depletolite",
         anova == "producer_specific")%>%
  select(feature_number, DayNight)
  

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
  mutate(sum = sum(log2_change))%>%
  filter(sum != 0)%>%
  dplyr::select(-sum)%>%
  mutate(Organism = factor(Organism))%>%
  mutate(Organism = fct_relevel(Organism, organism_order_micro))%>%
  nest()%>%
  mutate(dunnett = map(data, ~ aov(log2_change ~ Organism, .x)%>%
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
  mutate(sum = sum(log2_change))%>%
  filter(sum != 0)%>%
  dplyr::select(-sum)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(log2_change ~ DayNight, .x)%>%
                      tidy()))%>%
  unnest(data)%>%
  dplyr::select(-c(4:7))%>%
  filter(term != "Residuals")%>%
  ungroup()%>%   
  add_column(FDR = p.adjust(.$p.value, method = "BH"))%>%
  filter(FDR < 0.05)

# STATS POST-HOC -- DOM Dunnetts ------------------------------------------
dunnett_factors_dom <- as.factor(anova_dom_t0_df$Organism)%>%
  relevel("Water control")%>%
  levels()%>%
  as.vector()

dunnetts_dom <- anova_dom_t0_df%>%
  right_join(sigs_osm, by = c("feature_number", "DayNight"))%>%
  group_by(feature_number, DayNight)%>%
  mutate(sum = sum(log))%>%
  filter(sum != 0)%>%
  dplyr::select(-sum)%>%
  mutate(Organism = factor(Organism)%>%
           fct_relevel(dunnett_factors_dom))%>%
  nest()%>%
  mutate(dunnett = map(data, ~ aov(log ~ Organism, .x)%>%
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
  
# META-STATS --Compounds prevalance ---------------------------------------
compound_prevalance <- t_pvals%>%
  mutate(exudate_type = 1)%>%
  mutate(fdr_mixed = case_when(FDR_greater < 0.05 ~ FDR_greater,
                               FDR_lesser < 0.05 ~ FDR_lesser))%>%
  select(c(1:3, activity, fdr_mixed))%>%
  unite(org_act, c(Organism, activity), sep = "_")%>%
  spread(org_act, fdr_mixed)%>%
  group_by(DayNight, feature_number)%>%
  nest()%>%
  mutate(data = map(data, ~ mutate(.x, reactivity = paste(colnames(.)[!is.na(.) > 0], sep = ", ", collapse = ", "))))%>%
  unnest(data)%>%
  left_join(networking, by = "feature_number")%>%
  group_by(DayNight, reactivity)%>%
  nest()

# t_test_features <- compound_prevalance%>%
#   left_join(networking, by = "feature_number")
# 
# grouped_t_test_features <- t_test_features%>%
#   group_by(exudate_type, DayNight, activity)%>%
#   nest()

t_test_feature_inchi <- t_test_features%>%
  filter(activity != "recalcitrant")%>%
  ungroup()%>%
  select(feature_number, inchi_key)%>%
  unique()%>%
  filter(!is.na(inchi_key))

write_tsv(t_test_feature_inchi, "~/Documents/GitHub/DORCIERR/data/analysis/inchi_key_norecal.tsv")



# META-STATS -- Major Depleteolites ---------------------------------------
major_deplete <- feature_table_no_back_trans%>%
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
                              TRUE ~ "Night"),
         xic = case_when(xic == 0 ~ 1000,
                         TRUE ~ as.numeric(xic)))%>%
  spread(Timepoint, xic)%>%
  group_by(Organism, DayNight, feature_number)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  mutate(log2_change = log2(TF/T0))%>%
  ungroup()%>%
  select(-c(T0, TF))%>%
  spread(Organism, log2_change)%>%
  right_join(dunnetts_dom, by = c("DayNight", "feature_number"))%>%
  left_join(networking, by = "feature_number")%>%
  group_by(DayNight, feature_number, Organism)
  
major_deplete$max_log <- apply(major_deplete[3:8], 1, min)

major_depletolites <- major_deplete%>%
  filter(max_log < -3.3)

# META-STATS -- Network changes -------------------------------------------
network_means <- dom_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = "feature_number")%>%
  filter(network != "-1")%>%
  group_by(network, Organism, DayNight, Timepoint)%>%
  summarize_if(is.numeric, mean)%>%
  spread(Timepoint, log)%>%
  mutate(diff = TF-T0)%>%
  right_join(net_test, by = c('Organism', 'DayNight', 'network'))%>%
  filter(activity != 'recalcitrant')
  


# EXPORT -- Major depletolites for Classyfire annotations --------------------
classyfire_features <- major_depletolites$feature_number%>%
  unique()

classyfire_export <- networking%>%
  filter(feature_number %in% classyfire_features)

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

hc_compounds <- dom_stats_wdf%>%
  filter(feature_number %in% producer_specific_change)%>%
  mutate(log = case_when(is.na(log) ~ 1000,
                         TRUE ~ as.numeric(log)))%>%
  spread(Timepoint, log)%>%
  group_by(Organism, DayNight, feature_number)%>%
  mutate(mean_t0 = mean(T0, na.rm = TRUE),
         feature_difference = TF-mean_t0)%>%
  ungroup()%>%
  mutate(feature_difference = feature_difference + min(feature_difference))%>%
  group_by(feature_number)%>%
  mutate(zscore = zscore(feature_difference))%>%
  ungroup()%>%
  left_join(networking%>%
              select(feature_number, combined_ID), 
            by = "feature_number")%>%
  unite(sample, c("Organism", "DayNight", "Replicate"), sep = "_")%>%
  unite(feature, c(feature_number, combined_ID), sep = "_")%>%
  dplyr::select(-c(TF, T0, mean_t0, Experiment, feature_difference))%>%
  spread(feature, zscore)

hc_df <- hc_compounds%>%
  left_join(hc_microbe, by = "sample")

hc_compounds[is.na(hc_compounds)] <- 0

write_csv(hc_df%>%
            select(everything(), sample), "~/Documents/GitHub/DORCIERR/data/plots/combined_hc_df.csv")

write_csv(hc_microbe%>%
            select(everything(), sample), "~/Documents/GitHub/DORCIERR/data/plots/microbe_hc_df.csv")

write_csv(hc_compounds%>%
            select(everything(), sample), "~/Documents/GitHub/DORCIERR/data/plots/compounds_hc_df.csv")



# META-STATS -- Network stats summary sheet -------------------------------
net_summary <- net_activity%>%
  left_join(metadata%>%
              mutate(num_features = 1)%>%
              select(num_features, network)%>%
              group_by(network)%>%
              summarize_if(is.numeric, sum)%>%
              ungroup(), 
            by = 'network')%>%
  left_join(networking%>%
              group_by(feature_number)%>%
              mutate(min_nosc = min(NOSC, na.rm = TRUE),
                     max_nosc = max(NOSC, na.rm = TRUE),
                     nc = N/C,
                     pc = P/C)%>%
              select(feature_number, network, min_nosc, max_nosc, C:P, nc, pc)%>%
              unique(),
            by = 'network')%>%
  left_join(log2_change_vals%>%
              group_by(feature_number)%>%
              mutate(min_log2 = min(log2_change),
                     max_log2 = max(log2_change))%>%
              select(feature_number, min_log2, max_log2)%>%
              unique(),
            by = 'feature_number')%>%
  left_join(molnet_class%>%
              left_join(networking, by = 'feature_number')%>%
              select(c(network, molnet_string))%>%
              filter(network != '-1')%>%
              unique(),
            by = 'network')%>%
  left_join(nosc_quadrants%>%
              filter(network != '324', DayNight == "Day")%>%
              select(-c('Replicate', 'T0', 'TF', 'NOSC', 'row m/z'))%>%
              spread(Organism, log2_change)%>%
              select(-activity)%>%
              group_by(feature_number, network, DayNight)%>%
              summarize_all(paste0, collapse = ", ")%>%
              ungroup()%>%
              mutate_all(funs(gsub("NA, ", "", .)))%>%
              mutate_all(funs(gsub(", NA", "", .)))%>%
              na_if("NA")%>%
              mutate(network = as.numeric(network)),
            by = c('feature_number', 'network'))%>%
  left_join(networking%>%
              select(feature_number, combined_ID, binary_ID, NOSC), 
            by = 'feature_number')%>%
  select(feature_number, network, activity, num_features, DayNight, NOSC, quadrant, combined_ID, binary_ID, everything())%>%
  inner_join(day_exometabolites, by = 'feature_number')%>%
  left_join(benthic_produced_exometabolites, by = 'feature_number')%>%
  left_join(unique_benthic_metabolites%>%
              rename(unique_organism = Organism),
            by = 'feature_number')

write_csv(net_summary, './analysis/depletolite_net_summary.csv')



# META-STATS -- unexpected depletolites and accumolites -------------------
unexpected_features_deplet <- net_activity%>%
  left_join(log2_change_vals%>%
              group_by(DayNight, Organism, feature_number)%>%
              summarize_if(is.numeric, mean)%>%
              left_join(networking, by = 'feature_number'),
            by = 'network')%>%
  filter(log2_change > 1,
         activity == 'depletolite',
         DayNight == 'Day')
  
pdf("plots/depletolite_networks_unexpected_accum.pdf", width = 6, height = 5)
unexpected_features_deplet%>%
  ggplot(aes(CF_subclass, fill = Organism)) +
  geom_histogram(stat = 'count') + 
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic"))
dev.off()

unexpected_features_accum <- net_activity%>%
  left_join(log2_change_vals%>%
              group_by(DayNight, Organism, feature_number)%>%
              summarize_if(is.numeric, mean)%>%
              left_join(networking, by = 'feature_number'),
            by = 'network')%>%
  filter(log2_change < -1,
         activity == 'accumolite',
         DayNight == 'Day')

pdf("plots/accumolite_networks_unexpected_deplet.pdf", width = 6, height = 5)
unexpected_features_accum%>%
  ggplot(aes(CF_subclass, fill = Organism)) +
  geom_histogram(stat = 'count') + 
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic"))
dev.off()

print_df <- unexpected_features_accum%>%
  bind_rows(unexpected_features_deplet)%>%
  select(1:11)%>%
  filter(binary_ID == 1)

write_csv(print_df, '~/Downloads/print_df_dorc.csv')

networks_plotting <- net_activity%>%
  left_join(log2_change_vals%>%
              group_by(DayNight, Organism, feature_number)%>%
              summarize_if(is.numeric, mean)%>%
              left_join(networking, by = 'feature_number'),
            by = 'network')%>%
  mutate(coloring = case_when(log2_change < -3.3 ~ 'depletolite',
                              log2_change > -3.3 & log2_change < 3.3 ~ 'recalcitrant',
                              TRUE ~ 'accumolite'),
         expected = case_when(log2_change > 1 & activity == 'depletolite' ~ 'unexpected',
                              log2_change < -1 & activity == 'accumolite' ~ 'unexpected',
                              TRUE ~ 'expected'))

networks_plotting%>%
  filter(activity == 'accumolite',
         DayNight == 'Day')%>%
  ggplot(aes(NOSC, log2_change, color = coloring)) +
  geom_text(aes(label = network))

pdf("plots/depletolite_networks_NOSC_unexpectedlogs.pdf", width = 15, height = 10)
networks_plotting%>%
  filter(activity == 'depletolite',
         DayNight == 'Day')%>%
  ggplot(aes(NOSC, log2_change, color = coloring)) +
  geom_text(aes(label = network))

networks_plotting%>%
  group_by(network)%>%
  filter(activity == 'accumulite',
         DayNight == 'Day',
         mean(log2_change) >= 1)%>%
  ungroup()%>%
  arrange(network)%>%
  mutate(network = as.character(network))%>%
  ggplot(aes(network, fill = expected)) +
  geom_histogram(stat = 'count', position = 'dodge') +
  ylim(c(0,50)) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic"))
dev.off()


# GRAPHING -- N:C P:C ratios ----------------------------------------------
unique_nc <- net_summary%>%
  filter(num_organisms == 1)%>%
  left_join(log2_change_vals%>%
              filter(DayNight == 'Day')%>%
              select(feature_number, Organism, T0, TF, log2_change, Replicate)%>%
              rename(unique_organism = Organism),
            by = c('feature_number', 'unique_organism'))%>%
  group_by(unique_organism, Replicate)%>%
  mutate(sample_t0_tic = sum(T0, na.rm = TRUE),
         sample_c = sum(C, na.rm = TRUE),
         sample_n = sum(N, na.rm = TRUE),
         sample_p = sum(P, na.rm = TRUE),
         weighted_nc = N/(C*T0/(sample_c*sample_t0_tic)),
         weighted_pc = P/(C*T0/(sample_c*sample_t0_tic)),
         weighted_p = (P*T0/(sample_p*sample_t0_tic)),
         weighted_n = (N*T0/(sample_n*sample_t0_tic)))%>%
  ungroup()
  
unique_nc%>%
  group_by(unique_organism, Replicate)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  mutate(Replicate = as.factor(Replicate))%>%
  ggplot(aes(n, p, color = unique_organism, shape = Replicate)) +
  geom_point(stat = 'identity') +
  scale_color_manual(values = wes_palette('Darjeeling1', 5, type = 'continuous'))

# GRAPHING -- [OSM] Pvalues, Log2 Volcano plot -----------------------------------------
volcano <- feature_table_no_back_trans%>%
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
                              TRUE ~ "Night"),
         xic = case_when(xic == 0 ~ 1000,
                         TRUE ~ as.numeric(xic)))%>%
  spread(Timepoint, xic)%>%
  group_by(Organism, DayNight, feature_number)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  mutate(log2_change = log2(TF/T0))%>%
  ungroup()%>%
  left_join(t_pvals, by = c("Organism", "DayNight", "feature_number"))%>%
  mutate(fdr_combined = case_when(FDR_greater > FDR_lesser ~ FDR_lesser,
                                  FDR_lesser > FDR_greater ~ FDR_greater,
                                  TRUE ~ 1),
         fdr_color = case_when(fdr_combined <= 0.05 ~ "Significant",
                               TRUE ~ "Not Signifant"))

volcano_themes <- function(x) {
  ggplot(x, aes(log2_change, log10(T0), col = fdr_color)) +
  geom_point(stat = "identity", size = 2.5, shape = 1, alpha = 0.8) +
  scale_x_continuous(breaks= seq(-20, 18, 2)) +
  scale_y_continuous(breaks = seq(0, 10, 1)) +
  # scale_shape_manual(values = c(1,19)) +
  scale_color_manual(values = c("gray38", "darkred")) +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic"))
}

pdf("./plots/osm_volcano.pdf", width = 6, height = 5)
volcano%>%
  volcano_themes()

volcano%>%
  volcano_themes() +
  geom_vline(xintercept = -3.3, col = "red", linetype = "dashed") +
  geom_vline(xintercept = 3.3, col = "red", linetype = "dashed")

# volcano%>%
#   volcano_themes() +
#   geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
#   geom_vline(xintercept = 3.3, col = "red", linetype = "dashed") +
#   geom_vline(xintercept = -3.3, col = "red", linetype = "dashed")
dev.off()  

# GRAPHING -- [OSM] Major depletolites ------------------------------------------
osm_dunnett_features <- (dunnetts_dom%>%
                           filter(DayNight == "Day")%>%
                           select(-p.value)%>%
                           spread(Organism, FDR)%>%
                           group_by(DayNight, feature_number)%>%
                           add_column(num_cols = rowSums(.[3:ncol(.)] >= 0, na.rm = TRUE))%>%
                           filter(num_cols != 1)%>%
                           ungroup())$feature_number%>%
  as.vector()
  

osm_dunnetts <- dunnetts_dom%>%
  filter(feature_number %in% osm_dunnett_features)%>%
  left_join(feature_table_no_back_trans%>%
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
                                          TRUE ~ "Night"),
                     xic = case_when(xic == 0 ~ 1000,
                                     TRUE ~ as.numeric(xic)))%>%
              group_by(Organism, Timepoint, DayNight, feature_number)%>%
              mutate(Replicate = as.numeric(Replicate),
                     xic = case_when(sum(xic) == 4000 ~ xic + Replicate,
                                     sum(xic) == 2000 ~ xic + Replicate,
                                     TRUE ~ as.numeric(xic)))%>%
              ungroup()%>%
              spread(Timepoint, xic)%>%
              group_by(Organism, DayNight, feature_number)%>%
              summarize_if(is.numeric, mean, na.rm = TRUE)%>%
              mutate(log2_change = log2(TF/T0))%>%
              ungroup(), by = c("DayNight", "feature_number", "Organism"))%>%
  mutate(log2_change = as.numeric(as.character(log2_change)))%>%
  filter(log2_change <= -3.3)%>%
  left_join(molnet_class, by = c("feature_number"))%>%
  separate(molnet_string, c("CF_kingdom", "CF_superclass", "CF_class", "CF_subclass"), sep = ";")%>%
  gather(Timepoint, xic, T0:TF)

# ## Colors for heterocyclic
# colors_hetero <- c("#FF0000",
#                    "firebrick4",
#                    "#E69F00",
#                    "#32806E", 
#                    "#91A737",
#                    "olivedrab4",
#                    "olivedrab2",
#                    "darkorchid3", 
#                    "#F0E442", 
#                    "#0072B2",
#                    "steelblue3",
#                    "#5BBCD6")

## Making the PDF
pdf("~/Documents/GitHub/DORCIERR/data/plots/osm_compounds.pdf", height = 7, width = 13)
osm_dunnetts%>%
  filter(Timepoint == "T0")%>%
  # filter(`CF_class` != "NA",
  #        `CF_subclass` != "NA",
  #        # `CF_class` %like% "Lipids%" & !FinalClass %like% "%carnitine%"
  #        # `CF_class` %like% "Organohetero%"
  #        # FinalClass %like% "%carnitine%"
  #        `CF_class` %like% "Organic acids%"
  #        )%>%
  # rename(`Chemical Class` = `CF_subclass`)%>%
  # mutate(`Chemical Class` = case_when(`Chemical Class` %like% "%Carb%" ~ `Chemical Class`,
  #                                     TRUE ~ as.character(`Chemical Class`)))%>%
  # unite(compound, c("level 3", "FinalClass"), sep = " ")%>%
  ggplot(aes(Organism, xic, fill = CF_class)) +
  geom_bar(stat = "identity", position = "stack") +
  # scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "darkorchid3")) + # Colorblind pallette
  facet_wrap(~CF_superclass) +
  # coord_flip() +
  theme(
    # legend.position = "none",
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(angle = 60, size = 20, hjust = 1),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )

osm_dunnetts%>%
  filter(`CF_class` != "NA",
         `CF_superclass` %like% "Lipids%"
         # `CF_class` %like% "Organohetero%"
         # FinalClass %like% "%carnitine%"

  )%>%
  # unite(compound, c("level 3", "FinalClass"), sep = " ")%>%
  ggplot(aes(Timepoint, xic, fill = CF_subclass)) +
  geom_bar(stat = "identity", position = "stack") +
  # scale_fill_manual(values = "#0072B2") +
  facet_wrap(~Organism) +
  ggtitle("Lipids and Lipid-like compounds") +
  # coord_flip() +
  theme(
    # legend.position = "none",
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(angle = 60, size = 20,  hjust = 1),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )

osm_dunnetts%>%
  filter(`CF_class` != "NA",
         # `CF_class` %like% "Lipids%" & !FinalClass %like% "%carnitine%"
         `CF_superclass` %like% "Organohetero%"
  )%>%
  ggplot(aes(Timepoint, xic, fill = CF_subclass)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Organism) +
  ggtitle("Organoheterocyclic compounds") +
  theme(
    # legend.position = "none",
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(angle = 60, size = 20, hjust = 1),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )

osm_dunnetts%>%
  filter(`CF_class` != "NA",
         `CF_superclass` %like% "Organic acids%" 
         # & !FinalClass %like% "%carnitine%"

  )%>%
  ggplot(aes(Timepoint, xic, fill = CF_subclass)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Organism) +
  ggtitle("Organic acids and derivatives") + 
  theme(
    # legend.position = "none",
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(angle = 60, size = 20,  hjust = 1),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )

dev.off()

osm_dunnetts_hc <- dunnetts_dom%>%
  filter(feature_number %in% osm_dunnett_features)%>%
  left_join(feature_table_no_back_trans%>%
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
                                          TRUE ~ "Night"),
                     xic = case_when(xic == 0 ~ 1000,
                                     TRUE ~ as.numeric(xic)))%>%
              group_by(Organism, Timepoint, DayNight, feature_number)%>%
              mutate(Replicate = as.numeric(Replicate),
                     xic = case_when(sum(xic) == 4000 ~ xic + Replicate,
                                     sum(xic) == 2000 ~ xic + Replicate,
                                     TRUE ~ as.numeric(xic)))%>%
              ungroup()%>%
              spread(Timepoint, xic)%>%
              group_by(Organism, DayNight, feature_number)%>%
              mutate(log2_change = log2(TF/mean(T0, na.rm = TRUE)))%>%
              ungroup()%>%
              select(-c(T0, TF)), by = c("DayNight", "feature_number", "Organism"))%>%
  mutate(log2_change = as.numeric(as.character(log2_change)))%>%
  left_join(molnet_class%>%
              select(-Organism), by = c("feature_number"))%>%
  separate(molnet_string, c(CF_kingdom, CF_superclass, CF_class, CF_subclass), sep = ";")%>%
  unite(identifier, c("CF_superclass", "CF_class", "feature_number"), sep = " ")%>%
  unite(sample, c("Organism", "Replicate"), sep = "_")%>%
  group_by(sample)%>%
  mutate(zscore = zscore(log2_change + 17.6))%>%
  select(sample, identifier, zscore)%>%
  spread(identifier, zscore)

osm_dunnetts_hc[is.na(osm_dunnetts_hc)] <- 0

write_csv(osm_dunnetts_hc, "./analysis/osm_dom_heat.csv")


# GRAPHING -- [OSM] Important OTUs --------------------------------------------------------
osm_otus <- dunnett_microbe_pvals%>%
  filter(DayNight == "Day")%>%
  left_join(microbe_taxonomy, by = "OTU")%>%
  inner_join(osm_ra_bigger_TF, by = c("DayNight", "Organism", "OTU"))%>%
  left_join(average_ra, by = c("OTU", "Organism", "DayNight"))%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "sp"), sep = ";")%>%
  unite(Tax_plot, c("Order", "Family", "Genus", "OTU"), sep = " ", remove = FALSE)

# colors_otus <- c("#FF0000", #Altermonas Red
#                  "firebrick4",
#            "#32806E", #Flavobacters greenish
#            "#91A737", 
#            "olivedrab4", 
#            "olivedrab2",
#            "olivedrab2",
#            "olivedrab2",
#            "olivedrab2",
#            "turquoise3", # Rhodobacters Blue
#            "#5BBCD6", 
#            "steelblue3", 
#            "steelblue3",
#            "royalblue2", 
#            "slateblue4", 
#            "slateblue4", 
#            "slateblue4",
#            "slateblue4") 

osm_ots_heat <- osm_otus%>%
  filter(log2_change >= 1,
         TF >= 0.001)%>%
  select(Tax_plot, OTU)%>%
  distinct(Tax_plot, OTU)%>%
  left_join(microbe_no_rare, by = c("OTU"))%>%
  filter(DayNight == "Day",
         Organism != "Water control")%>%
  # group_by(Organism, Tax_plot, Replicate)%>%
  # summarize_if(is.numeric, sum)%>%
  # ungroup()%>%
  group_by(Tax_plot)%>%
  mutate(log2_change = zscore(log2_change + 1.82))%>%
  ungroup()%>%
  select(Tax_plot, Organism, Replicate, log2_change)%>%
  unite(sample, c("Organism", "Replicate"), sep = "_")%>%
  spread(Tax_plot, log2_change)

write_csv(osm_ots_heat, "./analysis/osm_heat.csv")

# pdf("~/Documents/GitHub/DORCIERR/data/plots/osm_otus.pdf", width = 7, height = 5)
# osm_otus%>%
#   filter(log2_change >= 1,
#          TF >= 0.001)%>%
#   ggplot(aes(x = Organism, y = ra, fill = Tax_plot)) +
#   geom_bar(stat = "summary", fun.y = "mean", position = "stack") +
#   facet_wrap(~DayNight) +
#   # scale_color_manual(values = "black") +
#   scale_fill_manual(values = colors_otus) +
#   theme(
#     # legend.position = "none",
#     # plot.margin = margin(2,.8,2,.8, "cm"),
#     axis.text.x = element_text(angle = 60, hjust = 1),
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     panel.grid.major = element_blank(), # get rid of major grid
#     panel.grid.minor = element_blank(), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     legend.box.background = element_rect(fill = "transparent"),
#     legend.position = "none"# get rid of legend panel bg
#   ) +
#   ylab("Relative Abundance") +
#   ggtitle("OTUs")
# 
# osm_otus%>%
#   filter(log2_change >= 1,
#          TF >= 0.001)%>%
#   ggplot(aes(x = Organism, y = log2_change, fill = Tax_plot)) +
#   geom_bar(stat = "summary", fun.y = "mean", position = "stack") +
#   facet_wrap(~DayNight) +
#   # scale_color_manual(values = "black") +
#   # scale_fill_manual(values = colors_otus) +
#   theme(
#     # legend.position = "none",
#     # plot.margin = margin(2,.8,2,.8, "cm"),
#     axis.text.x = element_text(angle = 60, hjust = 1),
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     panel.grid.major = element_blank(), # get rid of major grid
#     panel.grid.minor = element_blank(), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     legend.box.background = element_rect(fill = "transparent")
#     # legend.position = "none"# get rid of legend panel bg
#   ) +
#   ggtitle("OTUs")
# dev.off()

# GRAPHING -- [OSM] Correlation -------------------------------------------------------
osm_large_otu <- osm_otus%>%
  filter(log2_change >= 1,
         TF >= 0.001)%>%
  group_by(OTU)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  select(OTU)%>%
  left_join(microbe_no_rare, by = c("OTU"))%>%
  select(-c(T0, TF))%>%
  spread(OTU, log2_change)


#Depletolites
osm_features_corr <- osm_dunnetts$feature_number%>%
  as.vector()%>%
  unique()

osm_major_depletolites <- major_depletolites%>%
  # filter(feature_number %in% osm_features_corr)%>%
  group_by(feature_number)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  select(feature_number)%>%
  left_join(feature_table_no_back_trans%>%
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
                                          TRUE ~ "Night"),
                     xic = case_when(xic == 0 ~ 1000,
                                     TRUE ~ as.numeric(xic)))%>%
              group_by(Organism, Timepoint, DayNight, feature_number)%>%
              mutate(Replicate = as.numeric(Replicate),
                     xic = case_when(sum(xic) == 4000 ~ xic + Replicate, 
                                     sum(xic) == 2000 ~ xic + Replicate,
                                     TRUE ~ as.numeric(xic)),
                     Replicate = as.character(Replicate))%>%
              ungroup()%>%
              spread(Timepoint, xic)%>%
              group_by(Organism, DayNight, feature_number)%>%
              mutate(log2 = log2(TF/mean(T0, na.rm = TRUE)))%>%
              ungroup(), by = c("feature_number"))%>%
  select(-c("T0", "TF"))%>%
  spread(feature_number, log2)

osm_correlation <- osm_major_depletolites%>%
  left_join(osm_large_otu, by = c("Organism", "DayNight", "Replicate"))%>%
  gather(feature_number, log10, 5:685)%>%
  gather(OTU, log2_change, 6:31)%>%
  select(-c(contains("Experiment")))
  
osm_metabolite_metabolite_corr <- osm_major_depletolites%>%
  select(-Experiment)%>%
  left_join(osm_major_depletolites%>%
              select(-Experiment), by = c("Organism", "DayNight", "Replicate"))%>%
  gather(feature_number, log10, contains(".x"))%>%
  gather(feature_number.y, log10.y, contains(".y"))

osm_corr_test <- osm_correlation%>%
    group_by(OTU, feature_number)%>%
    nest()%>%
    mutate(data = map(data, ~ cor.test(.x$log10, .x$log2_change, method = "pearson")%>%
                               broom::tidy()))

# osm_metab_corr_test <- osm_metabolite_metabolite_corr%>%
#   group_by(feature_number, feature_number.y)%>%
#   nest()%>%
#   mutate(data = map(data, ~ cor.test(.x$log10, .x$log10.y, method = "pearson")%>%
#                       broom::tidy()))

write_csv(osm_metab_corr_test%>%
            unnest(data), "analysis/metab_corr_results_unensted.csv")

## Edge table
osm_corr_pvals <- osm_corr_test%>%
  unnest(data)%>%
  # left_join(networking, by = c("feature_number"))%>%
  mutate(fdr = p.adjust(p.value, method = "BH"))%>%
  filter(fdr < 0.0001)%>%
  select(1,2,fdr)

osm_metab_corr_pvals <- osm_metab_corr_test%>%
  unnest(data)%>%
  filter(feature_number != feature_number.y)%>%
  # left_join(networking, by = c("feature_number"))%>%
  mutate(fdr = p.adjust(p.value, method = "BH"))%>%
  filter(fdr < 0.001)%>%
  select(1,2,fdr)

osm_corr_combined <- osm_metab_corr_pvals%>%
  ungroup()%>%
  filter(fdr < 0.000000001)%>%
  mutate(feature_number = gsub(".x", "", feature_number),
         feature_number.y = gsub(".y", "", feature_number.y))%>%
  rename(node_2 = feature_number.y)%>%
  bind_rows(osm_corr_pvals%>%
              ungroup()%>%
              rename(node_2 = OTU))%>%
  filter(feature_number != node_2)
  
write_csv(osm_corr_combined, "./analysis/osm_cytoscape_correlations_edge.csv")

## Node table
corr_nodes <- osm_otus%>%
  select(-c(p.value:Class, Order, sp, T0, TF))%>%
  distinct(OTU, Tax_plot, Family, Genus, Organism, log2_change)%>%
  spread(Organism, log2_change)%>%
  mutate(sample = "ASV",
         OTU_code = OTU,
         label = Family,
         label = case_when(label == "Clade_I" ~ "SAR 11",
         Genus %like% "%NS5_%" ~ "NS5 Marine Group",
         Genus %like% "%Dong%" ~ "Donghicola",
         Genus %like% "Meso%" ~ "Mesoflavibacter",
         Genus %like% "Wino%" ~ "Winogradskyella",
         TRUE ~ as.character(label))
)%>%
  unite(color, c("label", "OTU_code"), sep = ";", remove = FALSE)%>%
  rename(shared_name = OTU)%>%
  bind_rows(molnet_class%>%
              mutate(sample = "metabolite")%>%
              select(-Organism)%>%
              separate(molnet_string, c(CF_kingdom, CF_superclass, CF_class, CF_subclass), sep = ";")%>%
              left_join(osm_dunnetts%>%
                          select(feature_number, Organism), by = "feature_number")%>%
              rename(shared_name = feature_number,
                     label = CF_subclass,
                     color = `CF_superclass`))

write_csv(corr_nodes, "analysis/osm_cyto_node.csv")


# GRAPHING -- data sleuthing Porites --------------------------------------
porites_depletes <- osm_dunnetts%>%
  filter(Organism == "Porites lobata")

porites_asv <- osm_otus%>%
  # filter(log2_change >= 1,
  #        TF >= 0.001)%>%
  select(Tax_plot, OTU)%>%
  distinct(Tax_plot, OTU)%>%
  left_join(microbe_no_rare, by = c("OTU"))%>%
  filter(DayNight == "Day",
         Organism == "Porites lobata")

# PLOT FOR CRAIG ----------------------------------------------------------
osm_correlation%>%
  filter(feature_number == "2734",
         OTU == "Otu0062")%>%
  ggplot(aes(log10, log2_change, col = Organism)) +
  geom_point()

pdf("./plots/cyromorphaceae.pdf", height = 7, width = 8)
osm_otus%>%
  filter(OTU == "Otu0062")%>%
  select(-log2_change)%>%
  gather(Timepoint, val, T0:TF)%>%
  ggplot(aes(Timepoint, val, fill == "gray")) +
  geom_bar(stat = "summary", fun.y = "mean") +
  # scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + # Colorblind pallette
  facet_wrap(~Organism) +
  # coord_flip() +
  theme(
    # legend.position = "none",
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(angle = 60, size = 20,  hjust = 1),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )
dev.off()

# GRAPHING -- [OSM] FCM data ----------------------------------------------------
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

pdf("~/Documents/GitHub/DORCIERR/data/plots/FCM_day.pdf", width = 7, height = 5)
fcm_graphing%>%
  ggplot(aes(x= Hours, y = `Cells µL-1`, color = Organism))+
  geom_point(stat = "identity") +
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
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.x = element_blank(), # get rid of major grid
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )
dev.off()


# GRAPHING -- [OSM] PCoAs Labile ------------------------------------------
osm_dom_pco <- dom_stats_wdf%>%
  spread(Timepoint, log)%>%
  filter(DayNight == "Day")%>%
  group_by(feature_number, Organism, DayNight)%>%
  mutate(mean_t0 = mean(T0, na.rm = TRUE))%>%
  mutate(change = TF - mean_t0)%>%
  ungroup()%>%
  mutate(zscore = ((change - mean(change, na.rm = TRUE))/sd(change, na.rm = TRUE)))%>%
  mutate(zscore = zscore + 78)%>%
  select(-c(TF, T0, mean_t0, change))%>%
  unite(sample, c(1:4), sep = "_")%>%
  column_to_rownames(var = "sample")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()

dom_stats_wdf%>%
  spread(Timepoint, log)%>%
  filter(DayNight == "Day")%>%
  group_by(feature_number, Organism, DayNight)%>%
  mutate(mean_t0 = mean(T0, na.rm = TRUE))%>%
  mutate(change = TF - mean_t0)%>%
  ungroup()%>%
  mutate(zscore = ((change - mean(change, na.rm = TRUE))/sd(change, na.rm = TRUE)))%>%
  mutate(zscore = zscore + 78)%>%
  select(-c(TF, T0, mean_t0, change))%>%
  spread(feature_number, zscore)%>%
  adonis(.[5:ncol(.)] ~ Organism, data = ., perm = 1000, method = 'bray', p.adjust.methods = "BH")

## Plot Eigenvalues
osm_dom_pco$values[1:10,]%>%
  as.data.frame()%>%
  rownames_to_column("Axis")%>%
  mutate(axis = as.numeric(Axis))%>%
  ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, color = "red", vjust = -0.5)

org_colors <- c('#ED220D', '#7DFB4C', '#F19E38', '#8B19F5', '#33330A', '#0F01C4')
## PCoA plot
pdf("./plots/osm_dom_pcoa.pdf", width = 6, height = 5)
osm_dom_pco$vectors%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample")%>%
  separate(sample, c("experiemnt", "Organism", "replicate", "DayNight"), sep = "_")%>%
  ggplot(., aes(x = Axis.1, y = Axis.2, color = Organism, shape = DayNight)) +
  geom_point(stat = "identity", aes(size = 0.2)) +
  scale_shape_manual(values = c(19)) +
  scale_color_manual(values = org_colors) +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic")) +
  xlab(str_c("Axis 1", " (", round(osm_dom_pco$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ylab(str_c("Axis 2", " (", round(osm_dom_pco$values$Relative_eig[2], digits = 4)*100, "%)", sep = "")) +
  ggtitle("Labile features")
dev.off()


# GRAPHING -- [OSM] PCoAs microbes ------------------------------------------
osm_pcoa_microbe <- microbe_no_rare%>%
  filter(DayNight == "Day")%>%
  unite(sample, c(Organism, DayNight, Replicate), sep = "_")%>%
  select(-c(Experiment, T0, TF))%>%
  group_by(sample, OTU)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  mutate(log2_change = log2_change + 19)%>%
  spread(OTU, log2_change)%>%
  column_to_rownames("sample")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()

## Plot Eigenvalues
osm_pcoa_microbe$values[1:10,]%>%
  as.data.frame()%>%
  rownames_to_column("Axis")%>%
  mutate(axis = as.numeric(Axis))%>%
  ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, color = "red", vjust = -0.5)

## PCoA plot
pdf("./plots/osm_microbes_pcoa.pdf", width = 6, height = 5)
osm_pcoa_microbe$vectors%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample")%>%
  separate(sample, c("Organism", "DayNight", "Replicate"), sep = "_")%>%
  ggplot(., aes(x = Axis.1, y = Axis.2, color = Organism, shape = DayNight)) +
  geom_point(stat = "identity", aes(size = 0.2)) +
  scale_shape_manual(values = c(19)) +
  # scale_color_manual(values = c("darkorchid3", "#50A45C", "#AF814B", "#5BBCD6")) +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic")) +
  xlab(str_c("Axis 1", " (", round(osm_pcoa_microbe$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ylab(str_c("Axis 2", " (", round(osm_pcoa_microbe$values$Relative_eig[2], digits = 4)*100, "%)", sep = "")) +
  ggtitle("microbial communities")
dev.off()



# GRAPHING —- PCoAs DOM --------------------------------------
#looking at exudate features
dom_pco <- dom_stats_wdf%>%
  spread(Timepoint, log)%>%
  group_by(feature_number, Organism, DayNight)%>%
  mutate(mean_t0 = mean(T0, na.rm = TRUE))%>%
  mutate(change = TF - mean_t0)%>%
  ungroup()%>%
  mutate(zscore = ((change - mean(change, na.rm = TRUE))/sd(change, na.rm = TRUE)))%>%
  mutate(zscore = zscore + 78)%>%
  select(-c(TF, T0, mean_t0, change))

dom_graphing <- dom_pco%>%
  spread(5,6)

dorc_all_pcoa <- dom_pco%>%
  unite(sample, c(1:4), sep = "_")%>%
  spread(2,3)%>%
  column_to_rownames(var = "sample")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()


#This plots Eigenvalues
dorc_all_pcoa$values[1:10,]%>%
  as.data.frame()%>%
  rownames_to_column("Axis")%>%
  mutate(axis = as.numeric(Axis))%>%
  ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, color = "red", vjust = -0.5)

# S3 method for pcoa
pco_scores_all <- dorc_all_pcoa$vectors%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample")%>%
  mutate(feature = "all")%>%
  separate(sample, c("experiemnt", "Organism", "replicate", "DayNight"), sep = "_")

# Dorc_all_labile_exudates plot
pdf("./plots/dom_pcoa.pdf", width = 6, height = 5)
dorc_all_pcoa$vectors%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample")%>%
  separate(sample, c("experiemnt", "Organism", "replicate", "DayNight"), sep = "_")%>%
  ggplot(., aes(x = Axis.1, y = Axis.2, color = Organism, shape = DayNight)) +
  geom_point(stat = "identity", aes(size = 0.2)) +
  scale_shape_manual(values = c(1,19)) +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic")) +
  xlab(str_c("Axis 1", " (", round(dorc_all_pcoa$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ylab(str_c("Axis 2", " (", round(dorc_all_pcoa$values$Relative_eig[2], digits = 4)*100, "%)", sep = "")) +
  ggtitle("All DOM features")
dev.off()

# GRAPHING -- PCoA's Microbes ---------------------------------------------
pcoa_microbe <- microbe_no_rare%>%
  # filter(OTU %in% micro_sigs_vector)%>%
  filter(Timepoint == "TF")%>%
  unite(sample, c(Organism, DayNight, Replicate), sep = "_")%>%
  select(-c(sample_code, Experiment, Timepoint,
            numOtus, sum, reads, ra))%>%
  group_by(sample, OTU)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  mutate(zscore = (log10 -mean(log10))/sd(log10))%>%
  mutate(zscore = zscore + 0.75)%>%
  select(-c(log10, `Cells µL-1`, cell_abun))%>%
  spread(OTU, zscore)%>%
  column_to_rownames("sample")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()

pcoa_microbe$values[1:10,]%>%
  as.data.frame()%>%
  rownames_to_column("Axis")%>%
  mutate(axis = as.numeric(Axis))%>%
  ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, color = "red", vjust = -0.5)

pdf("~/Documents/GitHub/DORCIERR/data/plots/microbe_pcoa.pdf", width = 6, height = 5)
pcoa_microbe$vectors%>%
  as.data.frame()%>%
  rownames_to_column("sample")%>%
  separate(sample, c("Organism", "DayNight", "Replicate"), sep = "_")%>%
  ggplot(., aes(x = Axis.1, y = Axis.2, color = Organism, shape = DayNight)) +
  geom_point(stat = "identity", aes(size = 0.2)) +
  scale_shape_manual(values = c(1,19)) +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic")) +
  xlab(str_c("Axis 1", " (", round(pcoa_microbe$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ylab(str_c("Axis 2", " (", round(pcoa_microbe$values$Relative_eig[2], digits = 4)*100, "%)", sep = ""))
dev.off()



# GRAPHING -- DOC and XIC -------------------------------------------------
pdf("plots/doc_xic_compare.pdf", width = 6, height = 5)
doc_log10%>%
  filter(DayNight == "Day",
         Organism != 'Influent',
         Organism != "Offshore")%>%
  ggplot(aes(Organism, DOC)) + 
  geom_boxplot() +
  facet_wrap(~Timepoint) + 
  theme(
    # legend.position = "none",
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(angle = 60, size = 20, hjust = 1),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )


dom_stats_wdf%>%
  filter(DayNight == "Day")%>%
  ggplot(aes(Organism, log)) + 
  geom_bar(stat = "summary", fun.y = "sum") +
  facet_wrap(~Timepoint) + 
  theme(
    # legend.position = "none",
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(angle = 60, size = 20, hjust = 1),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )
dev.off()


test <- dom_stats_wdf%>%
  group_by(Timepoint, Organism, DayNight)%>%
  summarize_if(is.numeric, mean)

write_csv(test, "test.csv")



# GRAPHING -- Extraction Efficiency ---------------------------------------
extraction_efficiency%>%
  filter(DayNight == "Day")%>%
  ggplot(aes(Organism, efficiency, color = Organism)) +
  geom_boxplot() +
  facet_wrap(~Timepoint) +
  theme(
    # legend.position = "none",
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(angle = 60, size = 20, hjust = 1),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )



# GRAPHING -- NOSC Vs Organism --------------------------------------------
org_nosc <- log2_change_vals%>%
  filter(DayNight == "Day")%>%
  group_by(feature_number, Organism, complete_removal)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  left_join(networking%>%
              select(feature_number, network, dG, NOSC, N, O, P, C, CLASS_STRING), by = "feature_number")%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`, `row retention time`), by = "feature_number")%>%
  left_join(net_activity, by = 'network')%>%
  inner_join(unique_benthic_metabolites, by = c('feature_number', 'Organism'))

org_nosc%>%
  mutate(Organism2 = Organism)%>%
  ggplot(aes(NOSC, log2_change, color = Organism)) +
  facet_wrap(~activity) +
  geom_point() +
  geom_smooth(method = lm, aes(color = Organism2))


# GRAPHING -- NOSC Vs. dG plot --------------------------------------------
nosc_plot <- log2_change_vals%>%
  filter(DayNight == "Day")%>%
  group_by(feature_number, Organism, complete_removal)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  filter(feature_number %in% osm_dunnett_features)%>%
  left_join(osm_dunnetts%>%
              select(feature_number, CF_superclass), by = "feature_number")%>%
  left_join(networking%>%
              select(feature_number, network, dG, NOSC, N, O), by = "feature_number")%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`, `row retention time`), by = "feature_number")%>%
  left_join(net_activity, by = 'network')%>%
  left_join(carbon_normalized_xic_NOSC%>%
              filter(Timepoint == "T0",
                     DayNight == "Day"), by = c("Organism", "feature_number"))

#Original NOSC plots
pdf('plots/lability_plots_052220.pdf', width = 6, height = 5)
nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'))%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(complete_removal, `row m/z`)) +
  geom_boxplot() +
  facet_wrap(~CF_superclass) +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red')

nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'))%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(complete_removal, log10(T0))) +
  geom_boxplot() +
  facet_wrap(~CF_superclass) +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red')

nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'))%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(complete_removal, log2_change)) +
  geom_boxplot() +
  facet_wrap(~CF_superclass) +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red')

nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'))%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(complete_removal, NOSC)) +
  geom_boxplot() +
  facet_wrap(~CF_superclass) +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red')

nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'))%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(complete_removal, carbon_normalized_NOSC)) +
  geom_boxplot() +
  facet_wrap(~CF_superclass) +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red')

nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'),
         complete_removal != "accumolite",
         activity == 'depletolite',)%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(`row m/z`, log2_change, color = carbon_normalized_NOSC)) +
  ggtitle('accumolites') +
  geom_point(stat = "identity") +
  geom_smooth(method = lm, color = 'grey') +
  xlim(c(125, 750)) +
  facet_wrap(~CF_superclass) +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red')

nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'),
         activity == 'depletolite',
         complete_removal == "semi-removed")%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(`row m/z`, log2_change, color = carbon_normalized_NOSC)) +
  geom_point(stat = "identity") +
  ggtitle('Not completely removed') +
  geom_smooth(method = lm, color = 'grey') +
  xlim(c(125, 750)) +
  facet_wrap(~CF_superclass) +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red')

nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'),
         activity == 'depletolite',
         complete_removal != "removed")%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(`row m/z`, log2_change, color = carbon_normalized_NOSC)) +
  ggtitle('Completely removed') +
  geom_point(stat = "identity") +
  geom_smooth(method = lm, color = 'grey') +
  xlim(c(125, 750)) +
  facet_wrap(~CF_superclass) +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red')

nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'),
         activity == 'depletolite',
         complete_removal == "semi-removed")%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(carbon_normalized_NOSC, log2_change, color = `row m/z`)) +
  geom_point(stat = "identity") +
  ggtitle('Not completely removed') +
  geom_smooth(method = lm, color = 'grey') +
  # xlim(c(125, 750)) +
  facet_wrap(~CF_superclass) +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red')

nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'),
         activity == 'depletolite',
         complete_removal == "semi-removed")%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(`row m/z`, NOSC, color = log2_change)) +
  geom_point(stat = "identity") +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red') +
  ggtitle('Not completely removed') +
  geom_smooth(method = lm, color = 'grey') +
  # xlim(c(125, 750)) +
  facet_wrap(~CF_superclass)

nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'),
         activity == 'depletolite',
         complete_removal == "removed")%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(`row m/z`, NOSC, color = log2_change)) +
  geom_point(stat = "identity") +
  ggtitle('Completely Removed') +
  geom_smooth(method = lm, color = 'grey') +
  # xlim(c(125, 750)) +x
  facet_wrap(~CF_superclass) +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red')

nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'),
         activity == 'depletolite',
         complete_removal == "semi-removed")%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(carbon_normalized_NOSC, log2_change, color = `row m/z`)) +
  geom_point(stat = "identity") +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red') +
  ggtitle('Not completely removed') +
  geom_smooth(method = lm, color = 'grey') +
  # xlim(c(125, 750)) +x
  facet_wrap(~CF_superclass)
  

nosc_plot%>%
  filter(CF_superclass %like any% c('Alkaloids%', 'Benzen%', 'Lipids%', 'Organic acid%', '%hetero%', '%propan%'),
         activity == 'depletolite',
         complete_removal == "semi-removed")%>%
  # ggplot(aes(dG, NOSC, color = log2_change)) +
  ggplot(aes(NOSC, log2_change, color = `row m/z`)) +
  geom_point(stat = "identity") +
  ggtitle('Not completely removed') +
  geom_smooth(method = lm, color = 'grey') +
  # xlim(c(125, 750)) +x
  facet_wrap(~CF_superclass) +
  scale_color_gradient2(low='#5011D1', mid = 'grey', high='red')


dev.off()


# SUMMARY -- Cytoscape ----------------------------------------------------
cyto <- networking%>%
  left_join(net_activity, by = 'network')%>%
  right_join(feature_table_relnorm%>%
               select(feature_number, contains('_xic'))%>%
               gather(sample_code, xic, 2:ncol(.))%>%
               separate(sample_code, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
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
               group_by(feature_number, Organism, Timepoint)%>%
               summarize_if(is.numeric, mean)%>%
               spread(Timepoint, xic), 
            by = 'feature_number')%>%
  mutate(xic_treat = case_when(activity == "depletolite" ~ T0,
                               TRUE ~ TF))%>%
  select(-c(canopus_annotation:CLASS_STRING, T0, TF))%>%
  spread(Organism, xic_treat)
  
write_csv(cyto, 'analysis/cytoscape_networks.csv')

# SUMMARY -- DOM-----------------------------------------------------------------
summary_no_background <- as.vector(feature_table_no_back_trans_filter%>%
                                     filter(background == "real"))$feature_number%>%
  length()

summary_no_back_trans <- as.vector(feature_table_no_back_trans_filter%>%
                                     filter(background == "real")%>%
                                     filter(Dorcierr_transient == "real"))$feature_number%>%
  length()                                    

summary_log2 <-log2_features%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(log2_change = map(data, ~ length(.x$feature_number)))%>%
  select(-data)

summary_ttest <- t_pvals%>%
  mutate(count = 1)%>%
  select(c(Organism, DayNight, activity, count))%>%
  group_by(Organism, DayNight, activity)%>%
  summarize_if(is.numeric, sum)%>%
  spread(DayNight, count)

summary_ttest_meta <- compound_prevalance%>%
  mutate(data = map(data, ~ mutate(.x, count = 1)%>%
                      select(count)%>%
                      summarize_if(is.numeric, sum, na.rm = TRUE)))%>%
  unnest(data)

summary_average_xic <- feature_table_no_back_trans%>%
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
                              TRUE ~ "Night"),
         xic = case_when(xic == 0 ~ 1000,
                         TRUE ~ as.numeric(xic)))%>%
  # filter(Timepoint == "T0")%>%
  # select(-Timepoint)%>%
  spread(Timepoint, xic)%>%
  group_by(Organism, DayNight, feature_number)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  mutate(log2_change = log2(TF/T0))%>%
  ungroup()%>%
  select(-c(T0, TF))%>%
  spread(Organism, log2_change)%>%
  right_join(t_pvals, by = c("DayNight", "feature_number"))

summary_dunnett <- dunnetts_dom%>%
  right_join(dom_stats_wdf%>%
              select(c(feature_number, Timepoint, DayNight, Organism, log))%>%
              group_by(feature_number,Timepoint, DayNight, Organism)%>%
              summarize_if(is.numeric, mean, na.rm = TRUE), by = c("feature_number", "DayNight", "Organism"))%>%
  left_join(networking, by = "feature_number")%>%
  select(-p.value)

summary_count_dunnett <- dunnetts_dom%>%
  select(activity, Organism, DayNight)%>%
  group_by(activity, Organism, DayNight)%>%
  mutate(count = 1)%>%
  summarize_if(is.numeric, sum)

# summary_ttest <- t_pvals%>%
#   mutate(depletolites = map(data, ~ filter(.x, activity == "depletolites")$feature_number%>%
#                               length()),
#          accumlites = map(data, ~ filter(.x, activity == "accumolites")$feature_number%>%
#                             length()),
#          recalcitrant = map(data, ~ filter(.x, activity == "recalcitrant")$feature_number%>%
#                               length()))%>%
#   select(-data)
# 
# summary_anova <- anova_dom_t0%>%
#   select(-c(feature_number, term, p.value, FDR))%>%
#   mutate(anova_groups_summary = 1)%>%
#   group_by(DayNight, activity, anova)%>%
#   summarize_if(is.numeric, sum)
  
write_csv(t_test_features%>%
            filter(`characterization scores` == "Good"), "~/Documents/GitHub/DORCIERR/data/analysis/t_test_features.csv")

write_csv(summary_dunnett%>%
            unite(sample, c(Organism, DayNight, Timepoint), sep = "_")%>%
            spread(sample, log), "~/Documents/GitHub/DORCIERR/data/analysis/cytoscape_dunnetts.csv")

# SUMMARY -- Organism chemical activity -----------------------------------
org_summary <- summary_average_xic%>%
  left_join(feature_table_no_back_trans%>%
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
                                          TRUE ~ "Night"),
                     xic = case_when(xic == 0 ~ 1000,
                                     TRUE ~ as.numeric(xic)))%>%
              filter(Timepoint == "T0")%>%
              select(-Timepoint)%>%
              # spread(Timepoint, xic)%>%
              group_by(Organism, DayNight, feature_number)%>%
              summarize_if(is.numeric, mean, na.rm = TRUE)%>%
              # mutate(log2_change = log2(TF/T0))%>%
              ungroup()%>%
              # select(-c(T0, TF))%>%
              spread(Organism, xic), 
            by = c("DayNight", "feature_number"), suffix = c("_log2", "_xic"))%>%
  left_join(networking, by = "feature_number")

poc_deplete <- org_summary%>%
  filter(Organism == "Pocillopora verrucosa",
         activity == "depletolite")

cca_deplete <- org_summary%>%
  filter(Organism == "CCA",
         activity == "depletolite")

dic_deplete <- org_summary%>%
  filter(Organism == "Dictyota",
         activity == "depletolite")

write_csv(poc_deplete, "./analysis/pocillopora_depletolites.csv")
write_csv(dic_deplete, "./analysis/dictyota_depletolites.csv")
write_csv(cca_deplete, "./analysis/cca_depletolites.csv")

org_summary%>%
  select(-c(contains("_log2"), contains("_xic")))%>%
  left_join(dom_stats_wdf%>%
              filter(Timepoint == "T0")%>%
              group_by(Organism, DayNight, feature_number)%>%
              summarize_if(is.numeric, mean)%>%
              ungroup(), by = c("DayNight", "feature_number", "Organism"))%>%
  filter(`characterization scores` == "Good",
         activity != "recalcitrant")%>%
  ggplot(aes(Organism, y = .$P/.$C, col = activity, size = log)) +
  scale_size_continuous(breaks = c(3,4,5,6,7,8,9,10), range = c(3,15)) +
  facet_wrap(~DayNight) +
  # geom_boxplot() +
  geom_point(stat = "identity") +
  theme(
    # legend.position = "none",
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
    # legend.position = "none"# get rid of legend panel bg
  )


# SUMMARY -- Microbe filtering  --------------------------------------------
microbe_summary <- microbe_combined%>%
  select(-c(Experiment, Organism, Replicate, Timepoint, DayNight, reads, asin, numOtus, sum))%>%
  spread(sample_code, ra)



# LINDAS DATASHEET  -------------------------------------------------------

# linda <- read_csv('~/Downloads/Dorcierr_depletolites_Reduced_July2020.csv')%>%
#   select(-Feat_MinLog2)
# 
# 
# linda_quadrant <- linda%>%
#   select(feature_number, quadrant, CCA:`Water control`)%>%
#   unique()%>%
#   group_by(feature_number)%>%
#   summarize_all(paste0, collapse = ", ")%>%
#   mutate_all(funs(gsub("NA, ", "", .)))%>%
#   mutate_all(funs(gsub(", NA", "", .)))%>%
#   na_if("NA")%>%
#   mutate(Feat_MinLog2 = apply(.[3:8], 1, max, na.rm = TRUE))%>%
#   left_join(linda%>%
#               select(-c(quadrant, CCA:`Water control`))%>%
#               unique()%>%
#               mutate(feature_number = as.character(feature_number)), 
#             by = 'feature_number')%>%
#   left_join(benthic_produced_exometabolites,
#             by = 'feature_number')
# 
# 
# write_csv(linda_quadrant, '~/Downloads/Dorcierr_depletolites_Reduced_July2020_zaq_changes.csv')

```