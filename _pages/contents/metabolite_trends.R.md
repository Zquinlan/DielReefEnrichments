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
library(ggforce)

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
    summarize_if(is.numeric, mean, na.rm = TRUE)%>%
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
         complete_removal = case_when(mean(T0, na.rm = TRUE) > 1002 & mean(TF, na.rm = TRUE) == 1002.5 ~ "removed",
                                      mean(TF, na.rm = TRUE) > mean(T0, na.rm = TRUE) & mean(TF, na.rm = TRUE) > 1003 ~"accumolite",
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
                              TRUE ~ "Night"),
         xic = case_when(xic == 0 ~ 1000,
                         TRUE ~ as.numeric(xic)))%>%
  group_by(Organism, Timepoint, DayNight, feature_number)%>%
  mutate(Replicate = as.numeric(Replicate),
         xic = case_when(sum(xic) == 4000 ~ xic + Replicate, 
                         sum(xic) == 2000 ~ xic + Replicate,
                         TRUE ~ as.numeric(xic)))%>%
  ungroup()

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




#Making the plots
plots_no_filter <- no_min_filter%>%
  group_by(feature_number, Timepoint, Organism)%>%
  filter(DayNight == 'Day',
         network %in% c('131','21','7','55','107','198','466', '165', '619', '627', '756', '924', '1314','80','141', '249','65','336','355','346'))%>%
  left_join(min_filter_feature_tag, by = c('feature_number', 'Organism'))%>%
  mutate(mean = mean(xic),
         sd = sd(xic),
         fill_color = case_when(is.na(xic_min_filter) ~ "#F98400",
                                TRUE ~ "#00A08A"),
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
#                      scale_fill_manual(values = c("#F2AD00", "#00A08A")) +
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
         fill_color = "#00A08A",
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
         fill_color = "#00A08A",
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
                       scale_fill_manual(values = c("#00A08A", "#F2AD00")) +
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
  mutate(xic = case_when(xic == 0 ~ 1000,
                         TRUE ~ as.numeric(xic)))%>%
  group_by(Organism, Timepoint, DayNight, feature_number)%>%
  mutate(Replicate = as.numeric(Replicate),
         xic = case_when(sum(xic) == 4000 ~ xic + Replicate,
                         sum(xic) == 2000 ~ xic + Replicate,
                         TRUE ~ as.numeric(xic)))%>%
  ungroup()%>%
  select(-c("Experiment", "Organism", "Replicate", "Timepoint", "DayNight"))


feature_table_combined <- right_join(metadata, feature_table_relnorm, by = "feature_number")

dorcierr_features_wdf  <- feature_table_combined%>%
  dplyr::select(c(feature_number, everything()))

# CLEANING-- adding carbon normalized values to wdf ------------------
# carbon_normalized_xic_NOSC <- dorcierr_table_wdf_temp%>%
#   filter(`characterization scores` == "Good")%>%
#   dplyr::select(c(feature_number, C, ends_with("_xic")))%>%
#   gather(sample_name, xic, 3:ncol(.))%>%
#   group_by(sample_name)%>%
#   mutate(ra = xic/sum(xic, na.rm = TRUE),
#          percent_total_C = xic*C)%>%
#   mutate(sum_c = sum(percent_total_C,  na.rm = TRUE),
#          carbon_norm_temp = percent_total_C/sum_c)%>%
#   ungroup()%>%
#   right_join(metadata%>%
#                select(c(feature_number, NOSC)),
#              .,  by = "feature_number")%>%
#   mutate(carbon_normalized_NOSC = carbon_norm_temp*NOSC,
#          sample_name = gsub("_xic", "", sample_name))%>%
#   select(c(feature_number, sample_name, carbon_normalized_NOSC))%>%
#   separate(sample_name, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
#   mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
#                                 TRUE ~ as.character(Experiment)))%>%
#   mutate(Organism = case_when(Organism == "CC" ~ "CCA",
#                               Organism == "DT" ~ "Dictyota",
#                               Organism == "PL" ~ "Porites lobata",
#                               Organism == "PV" ~ "Pocillopora verrucosa",
#                               Organism == "TR" ~ "Turf",
#                               Organism == "WA" ~ "Water control",
#                               TRUE ~ as.character(Organism)))%>%
#   separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
#   mutate(DayNight = case_when(DayNight == "D" ~ "Day",
#                               TRUE ~ "Night"))%>%
#   group_by(feature_number, Organism, Timepoint, DayNight)%>%
#   summarize_if(is.numeric, mean)%>%
#   filter(!carbon_normalized_NOSC < -0.001)



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
               select(feature_number, DayNight), by = c('DayNight', 'feature_number'))

feature_stats_wdf <- dom_stats_wdf%>%
  inner_join(min_filter%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('Organism', 'DayNight', 'feature_number'))%>%
  inner_join(org_exometabolites, by = c('feature_number', 'Organism', 'DayNight'))%>%
  group_by(Organism, Replicate, Timepoint, DayNight)%>%
  mutate(log10 = log10(xic),
         ra = xic/sum(xic),
         asin = asin(ra))

# SET SEED ----------------------------------------------------------------
set.seed(2005)


# STATS  -- T-TEST Network level ------------------------------------------
net_test <- dom_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = "feature_number")%>%
  filter(network != "-1")%>%
  group_by(network, DayNight)%>%
  nest()%>%
  mutate(greater = map(data, ~ t.test(xic ~ Timepoint, .x, alternative = "greater")),
         lesser = map(data, ~ t.test(xic ~ Timepoint, .x, alternative = "less")))%>%
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
  group_by(network, DayNight, activity)%>%
  summarize_if(is.numeric, mean)%>%
  select(network, DayNight, activity)

# STATS -- One-Sample T-test ----------------------------------------------
t_test <- feature_stats_wdf%>%
  select(-c(xic, ra, asin))%>%
  group_by(feature_number, Organism, DayNight)%>%
  spread(Timepoint, log10)%>%
  nest()%>%
  mutate(greater = map(data, ~ t.test(.$TF, mu = mean(.$T0, na.rm=TRUE), alternative = "greater")),
         lesser = map(data, ~ t.test(.$TF, mu = mean(.$T0, na.rm=TRUE), alternative = "less")))

t_pvals <- t_test%>%
  select(-data)%>%
  ungroup()%>%
  mutate(greater = map(greater, ~ .x["p.value"][[1]]),
         lesser = map(lesser, ~ .x["p.value"][[1]]),
         greater = as.numeric(greater),
         lesser = as.numeric(lesser),
         FDR_greater = p.adjust(greater, method = "BH"),
         FDR_lesser = p.adjust(lesser, method = "BH"),
         activity = case_when(FDR_greater < 0.05 ~ "accumolite",
                              FDR_lesser < 0.05 ~ "depletolite",
                              FDR_lesser >= 0.05 & FDR_greater >= 0.05 | 
                                is.na(FDR_lesser) & is.na(FDR_greater) ~ "recalcitrant"))

# STATS - OLD GLM Activity groupings -----------------------------------------------
#Difference within activity groupings
net_lm_df <- log2_change_vals%>%
  inner_join(min_filter%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  inner_join(org_exometabolites, by = c('feature_number', 'Organism', 'DayNight'))%>%
  inner_join(log2_features%>%
               select(feature_number, DayNight)%>%
               unique(), by = c('DayNight', 'feature_number'))%>%
  inner_join(t_pvals%>%
               filter(activity != 'recalcitrant')%>%
               select('feature_number', 'Organism', 'DayNight'),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(c(feature_number, network, NOSC)), by = "feature_number")%>%
  filter(network != "-1",
         NOSC <= 0)%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  left_join(net_activity, by = c('network', 'DayNight'))%>%
  mutate(activity = case_when(is.na(activity) ~ 'recalcitrant',
                              TRUE ~ as.character(activity)))%>% #This line was added because entire networks did not pass the log2 bottleneck
  select(-c(T0,TF))%>%
  gather(response_var, value, NOSC:`row m/z`)

net_glm <- net_lm_df%>%
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


org_glm <- net_lm_df%>%
  group_by(DayNight, activity, Organism, response_var)%>%
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



nosc_quadrants <- log2_change_vals%>%
  left_join(networking%>%
              select(feature_number, network, NOSC), by = "feature_number")%>%
  filter(network != "-1")%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  left_join(net_activity, by = c('network', 'DayNight'))%>%
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




# STATS - multiple regressions --------------------------------------------
mul_reg <- log2_change_vals%>%
  inner_join(feature_stats_wdf%>%
               ungroup()%>%
               select(feature_number, Organism, DayNight)%>%
               unique(), 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network, N, P, C, NOSC),
            by = 'feature_number')%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  filter(NOSC < 0,
         DayNight == 'Day')%>%
  mutate(n_presence = case_when(N > 0 ~ 1,
                                TRUE ~ 1))%>%
  mutate(nc = N/C,
         pc = P/C,
         sample_nc = T0*nc,  
         sample_pc = T0*pc, 
         log_snc = log(sample_nc),
         log_spc = log(sample_pc))

pdf("./plots/correlation_verify.pdf", width = 15, height = 10)
corr_verify <- mul_reg%>%
  group_by(feature_number, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  select(log2_change, `row m/z`, NOSC, sample_pc)%>%
  cor()%>%
  corrplot::corrplot()
dev.off()


#can just impute mean for NA values
# n_p_mulreg <- mul_reg%>%
#   filter(N > 0,
#          P > 0)%>%
#   group_by(feature_number, Organism)%>%
#   summarize_if(is.numeric, mean)%>%
#   ungroup()%>%
#   lm(log2_change ~ `row m/z`+NOSC+log_snc+log_spc, data = .)

# p_mulreg <- mul_reg%>%
#   filter(N == 0,
#          P > 0)%>%
#   group_by(feature_number, Organism)%>%
#   summarize_if(is.numeric, mean)%>%
#   ungroup()%>%
#   lm(log2_change ~ `row m/z`+NOSC+log_spc, data = .)

n_mulreg <- mul_reg%>%
  filter(N > 0)%>%
  group_by(feature_number, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  lm(log2_change ~ `row m/z`+NOSC+log_snc, data = .)

mass_mulreg <- mul_reg%>%
  filter(N == 0)%>%
  group_by(feature_number, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  lm(log2_change ~ `row m/z`+NOSC, data = .)

sink("./analysis/multiple_reg_coefficients_whole_metabolome.txt")
# summary(n_p_mulreg)
# summary(p_mulreg)
summary(n_mulreg)
summary(mass_mulreg)
sink()


# STATS -- Log-Lin regression NC ------------------------------------------
nc_loglin <- log2_change_vals%>%
  inner_join(t_pvals%>%
               # filter(activity != 'recalcitrant')%>%
               # rename(feat_activity = activity)%>%
               select(feature_number, Organism, DayNight)%>%
               unique(), 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network, C:dG)%>%
              filter(network != '-1'),
            by = 'feature_number')%>%
  left_join(net_activity, by = c('network', 'DayNight'))%>%
  group_by(Organism, DayNight, Replicate)%>%
  mutate(nc = N/C,
         activity2 = activity,
         sample_nc = T0*nc/sum(T0))%>%
  ungroup()%>%
  filter(activity != 'recalcitrant',
         nc > 0,
         !is.na(nc))%>%
  group_by(DayNight, activity)%>%
  nest()%>%
  mutate(model = map(data, ~ glm(log2_change ~ log(sample_nc), family = gaussian, .x)),
         p_vals = map(model, ~tidy(.x)%>%
                        filter(term == 'log(sample_nc)')),
         r2 = map(model, ~with(summary(.x), 1-deviance/null.deviance)))%>%
  select(-c(data, model))%>%
  unnest(p_vals)%>%
  ungroup()%>%
  mutate(r2 = as.numeric(r2),
         r = sqrt(r2),
         FDR = p.adjust(p.value, method = 'BH'))

pc_loglin <- log2_change_vals%>%
  inner_join(t_pvals%>%
               # filter(activity != 'recalcitrant')%>%
               # rename(feat_activity = activity)%>%
               select(feature_number, Organism, DayNight)%>%
               unique(), 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network, C:dG)%>%
              filter(network != '-1'),
            by = 'feature_number')%>%
  left_join(net_activity, by = c('network', 'DayNight'))%>%
  group_by(Organism, DayNight, Replicate)%>%
  mutate(pc = P/C,
         activity2 = activity,
         sample_pc = T0*pc/sum(T0))%>%
  ungroup()%>%
  filter(activity != 'recalcitrant',
         pc > 0,
         !is.na(pc))%>%
  group_by(DayNight, activity)%>%
  nest()%>%
  mutate(model = map(data, ~ glm(log2_change ~ log(sample_pc), family = gaussian, .x)),
         p_vals = map(model, ~tidy(.x)%>%
                        filter(term == 'log(sample_pc)')),
         r2 = map(model, ~with(summary(.x), 1-deviance/null.deviance)))%>%
  select(-c(data, model))%>%
  unnest(p_vals)%>%
  ungroup()%>%
  mutate(r2 = as.numeric(r2),
         r = sqrt(r2),
         FDR = p.adjust(p.value, method = 'BH'))

write_csv(nc_loglin, './analysis/nc_loglin_pvals.csv')
write_csv(pc_loglin, './analysis/pc_loglin_pvals.csv')

# STATS - ANOVA Activity groupings ----------------------------------------
# Difference between activity groupings
net_anova <- net_activity%>%
  left_join(networking%>%
              select(feature_number, network, NOSC), by = "network")%>%
  filter(network != "-1")%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  inner_join(min_filter%>%
               filter(DayNight == 'Day')%>%
               ungroup()%>%
               select(feature_number)%>%
               unique(),
             by = 'feature_number')%>%
  inner_join(org_exometabolites%>%
               filter(DayNight == 'Day')%>%
               ungroup()%>%
               select(feature_number)%>%
               unique(), by = c('feature_number'))%>%
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
  inner_join(min_filter%>%
               filter(DayNight == 'Day')%>%
               ungroup()%>%
               select(feature_number)%>%
               unique(),
             by = 'feature_number')%>%
  inner_join(org_exometabolites%>%
               filter(DayNight == 'Day')%>%
               ungroup()%>%
               select(feature_number)%>%
               unique(), by = c('feature_number'))%>%
  gather(response_var, value, NOSC:`row m/z`)%>%
  ggplot(aes(activity, value)) +
  # geom_jitter() +
  geom_boxplot() +
  facet_wrap(~response_var, scale = 'free_y')




# STATS -- ORGANISM ANOVA -------------------------------------------------
org_aov <- log2_change_vals%>%
  select(-c(T0,TF, complete_removal))%>%
  inner_join(t_pvals%>%
               filter(activity == 'depletolite')%>%
               select(feature_number, Organism, DayNight),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(anova = map(data, ~ aov(log2_change ~ Organism, .x)%>%
                           tidy()%>%
                           filter(!term == "Residuals")%>%
                           select(p.value)%>%
                           rename(anova = p.value)),
         tukey = map(data, ~ aov(log2_change ~ Organism, .x)%>%
                       TukeyHSD(p.adjust.methods = 'BH')%>%
                       tidy()))%>%
  select(-data)%>%
  unnest(c(anova, tukey))%>%
  ungroup()%>%
  filter(anova < 0.05,
         adj.p.value < 0.05)

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
            by = c('feature_number', 'DayNight'))

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



# VIZUALIZATIONS -- RGB hex codes for orgs --------------------------------
org_colors_no_water <- c("#CC0033","#669900", "#CC6600", "#9900FF", "#33CC33")

# VIZUALIZATIONS -- GLM ---------------------------------------------------
glm_df <- log2_change_vals%>%
  left_join(networking%>%
              select(feature_number, network, NOSC), by = "feature_number")%>%
  filter(network != "-1",
         NOSC <= 0)%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  left_join(net_activity, by = c('network', 'DayNight'))%>%
  filter(activity != 'recalcitrant',
         DayNight == 'Day')%>%
  mutate(activity2 = activity)%>% 
  select(-c(T0,TF))%>%
  gather(response_var, value, NOSC:`row m/z`)


glm_no_filt <-  glm_df%>%
  select(-network)%>%
  mutate(network = 1600,
         label = 'No filters',
         num_features = 3)

glm_log2_filt <- glm_df%>%
  inner_join(min_filter%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  inner_join(log2_features%>%
               select(feature_number, DayNight), by = c('DayNight', 'feature_number'))%>%
  select(-network)%>%
  mutate(network = 1600,
         label = 'Timepoint filters',
         num_features = 2)

glm_all_filt <- glm_df%>%
  inner_join(min_filter%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  inner_join(org_exometabolites, 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  inner_join(log2_features%>%
               select(feature_number, DayNight), 
             by = c('DayNight', 'feature_number'))%>%
  inner_join(t_pvals%>%
               filter(activity != 'recalcitrant')%>%
               select('feature_number', 'Organism', 'DayNight'),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  select(-network)%>%
  mutate(network = 1600,
         label = 'All filters',
         num_features = 1)%>%
  bind_rows(glm_no_filt, glm_log2_filt)%>%
  filter(complete_removal != 'removed')%>%
  group_by(network, num_features, label)%>%
  nest()%>%
  mutate(plots = map(data, ~ ggplot(.x, aes(value, log2_change, color = activity)) +
                       facet_wrap(~ response_var, scales = 'free_x') + 
                       geom_point(stat = 'identity') +
                       scale_color_manual(values = c('#78B7C5', '#EBCC2A')) +
                       new_scale('color') + 
                       geom_smooth(method = lm, aes(color = activity2)) +
                       scale_color_manual(values = c('#3B9AB2', '#E1AF00')) +
                       theme(
                         legend.position = "none",
                         plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
                         # plot.margin = margin(2,.8,2,.8, "cm"),
                         axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
                         axis.text.y = element_text(size = 20),
                         panel.background = element_rect(fill = "transparent"), # bg of the panel
                         plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                         panel.grid.major = element_blank(), # get rid of major grid
                         panel.grid.minor = element_blank(), # get rid of minor grid
                         legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                         legend.box.background = element_rect(fill = "transparent")
                       )))%>%
  select(-data)

filter_steps_na <- plots_all_filter%>% 
  select(-plots)%>% 
  spread(label, num_features)%>% 
  gather(label, num_features, 2:ncol(.))


pdf("./plots/nosc_trend_091520.pdf", width = 15, height = 10)
glm_all_filt$plots[1]
dev.off()


# VIZUALIZATIONS -- Logarithmic regression --------------------------------
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

pdf('./plots/glm_loglin_112320.pdf', width = 15, height = 10)
log_reg%>% 
  filter(response_var == 'NOSC',
         activity == 'depletolite')%>%
  mutate(NOSC = value)%>%
  ggplot(aes(NOSC, log2_change, color = activity)) +
  geom_point(stat = 'identity') +
  scale_color_manual(values = c('#EBCC2A')) + # '#78B7C5',  accumolites
  new_scale('color') +
  geom_smooth(method = 'lm', formula = y ~ x, aes(color = activity2)) +
  scale_color_manual(values = c('#E1AF00')) + # '#3B9AB2', accumolites
  new_scale('color') +
  scale_color_manual(values = c('#E1AF00')) +
  labs(x = 'NOSC', y = bquote(Log[2] ~change)) +
  geom_text(aes(x = -2.3, y = 12, color = 'depletolite',
                label = paste("FDR p-value: ", (net_glm%>%
                                                  filter(DayNight == 'Day',
                                                         activity == 'depletolite',
                                                         response_var == 'NOSC'))$FDR%>%
                                formatC(format = "e", digits = 2), sep = "")), size = 7) + 
  geom_text(aes(x = -2.3, y = 10.5, color = 'depletolite',
                label = paste("r²: ", (net_glm%>%
                                         filter(DayNight == 'Day',
                                                activity == 'depletolite',
                                                response_var == 'NOSC'))$r2%>%
                                round(digits = 4), sep = "")), size = 7) +
  # geom_text(aes(x = -2.285, y = 15, color = 'accumolite',
  #               label = paste("FDR p-value: ", (net_glm%>%
  #                                                 filter(DayNight == 'Day',
  #                                                        activity == 'accumolite',
  #                                                        response_var == 'NOSC'))$FDR%>%
  #                               formatC(format = "e", digits = 2), "*", sep = "")), size = 7) +
  # geom_text(aes(x = -2.3, y = 13.5, color = 'accumolite',
  #               label = paste("r²: ", (net_glm%>%
  #                                        filter(DayNight == 'Day',
  #                                               activity == 'accumolite',
  #                                               response_var == 'NOSC'))$r2%>%
  #                               round(digits = 4), sep = "")), size = 7) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 25, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )

log_reg%>% 
  filter(response_var == 'row m/z',
         activity == 'depletolite')%>%
  rename('row m/z' = value)%>%
  ggplot(aes(`row m/z`, log2_change, color = activity)) +
  geom_point(stat = 'identity') +
  scale_color_manual(values = c('#EBCC2A')) + # '#78B7C5', for accumolites
  new_scale('color') +
  geom_smooth(method = 'lm', formula = y ~ x, aes(color = activity2)) +
  scale_color_manual(values = c('#E1AF00')) + # '#3B9AB2', for accumolites
  new_scale('color') +
  scale_color_manual(values = c('#E1AF00')) +
  labs(x = 'm/z', y = bquote(Log[2] ~change)) +
  geom_text(aes(x = 251, y = 14, color = 'depletolite',
                label = paste("FDR p-value: ", (net_glm%>%
                                                  filter(DayNight == 'Day',
                                                         activity == 'depletolite',
                                                         response_var == 'row m/z'))$FDR%>%
                                formatC(format = "e", digits = 2), "***", sep = "")), size = 7) +
  geom_text(aes(x = 245, y = 12.5, color = 'depletolite',
                label = paste("r²: ", (net_glm%>%
                                                  filter(DayNight == 'Day',
                                                         activity == 'depletolite',
                                                         response_var == 'row m/z'))$r2%>%
                                round(digits = 4), sep = "")), size = 7) +
  # geom_text(aes(x = 245, y = 17, color = 'accumolite',
  #               label = paste("FDR p-value: ", (net_glm%>%
  #                                                 filter(DayNight == 'Day',
  #                                                        activity == 'accumolite',
  #                                                        response_var == 'row m/z'))$FDR%>%
  #                               formatC(format = "e", digits = 2), sep = "")), size = 7) +
  # geom_text(aes(x = 245, y = 15.5, color = 'accumolite',
  #               label = paste("r²: ", (net_glm%>%
  #                                        filter(DayNight == 'Day',
  #                                               activity == 'accumolite',
  #                                               response_var == 'row m/z'))$r2%>%
  #                               round(digits = 4), sep = "")), size = 7) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 25, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )
dev.off()

# GRAPHING -- BAD N:C P:C ratios ----------------------------------------------
unique_nc <- log2_change_vals%>%
  select(feature_number, Organism, T0, TF, log2_change, Replicate, DayNight)%>%
  inner_join(t_pvals%>%
               filter(activity != 'recalcitrant')%>%
               select(feature_number, Organism, DayNight)%>%
               unique(), 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network, C:dG)%>%
              filter(network != '-1'),
            by = 'feature_number')%>%
  left_join(net_activity,
            by = c('network', 'DayNight'))%>%
  group_by(Organism, Replicate)%>%
  mutate(sample_log2 = sum(log2_change, na.rm = TRUE),
         sample_c = sum(C, na.rm = TRUE),
         sample_n = sum(N, na.rm = TRUE),
         sample_p = sum(P, na.rm = TRUE),
         weighted_nc = N/C*(log2_change/(sample_c*sample_log2)),
         weighted_pc = P/C*(log2_change/(sample_c*sample_log2)),
         weighted_p = (P*T0/(sample_p*sample_log2)),
         weighted_n = (N*T0/(sample_n*sample_log2)))%>%
  ungroup()

pdf("./plots/nc_091520.pdf", width = 15, height = 10)
unique_nc%>%
  filter(activity != 'recalcitrant',
         !is.na(activity))%>%
  group_by(Organism, Replicate, activity)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  mutate(Replicate = as.factor(Replicate))%>%
  ggplot(aes(N/C, log2_change, color = Organism)) +
  geom_point(stat = 'identity') +
  scale_color_manual(values = org_colors_no_water) +
  facet_wrap(~ activity, scales = 'free_x') +
  theme(
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_line('grey'), # get rid of major grid
    panel.grid.minor = element_line('grey'), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    legend.position = 'top'
  )
dev.off()


  
nc_no_filt <- net_summary%>%
  select(-c(22:ncol(.)))%>%
  # gather(Organism, log2_change, 22:27)%>%
  # filter(!is.na(log2_change))%>%
  left_join(log2_change_vals%>%
              filter(DayNight == 'Day')%>%
              select(feature_number, Organism, T0, TF, log2_change, Replicate),
            by = c('feature_number'))%>%
  group_by(Organism, Replicate)%>%
  mutate(nc = N/C,
         sample_t0_tic = sum(T0, na.rm = TRUE),
         sample_c = sum(C, na.rm = TRUE),
         sample_n = sum(N, na.rm = TRUE),
         sample_p = sum(P, na.rm = TRUE),
         weighted_nc = N/(C*T0/(sample_c*sample_t0_tic)),
         weighted_pc = P/(C*T0/(sample_c*sample_t0_tic)),
         weighted_p = (P*T0/(sample_p*sample_t0_tic)),
         weighted_n = (N*T0/(sample_n*sample_t0_tic)))%>%
  ungroup()

nc_log_filt <- nc_no_filt%>%
  inner_join(log2_features%>%
               filter(DayNight == 'Day')%>%
               ungroup()%>%
               select(feature_number)%>%
               unique(),
             by = 'feature_number')%>%
  inner_join(min_filter%>%
               ungroup()%>%
               filter(DayNight == 'Day')%>%
               select(feature_number, Organism)%>%
               unique(),
             by = c('feature_number', 'Organism'))%>%
  mutate(label = 'Timepoint filters')

nc_all_filt <- nc_no_filt%>%
  inner_join(log2_features%>%
             filter(DayNight == 'Day')%>%
             ungroup()%>%
             select(feature_number)%>%
             unique(),
           by = 'feature_number')%>%
  inner_join(org_exometabolites%>%
               filter(DayNight == 'Day')%>%
               unique(), 
             by = c('feature_number', 'Organism'))%>%
  inner_join(min_filter%>%
               ungroup()%>%
               filter(DayNight == 'Day')%>%
               select(feature_number, Organism)%>%
               unique(),
             by = c('feature_number', 'Organism'))%>%
  mutate(label = 'All filters')

nc_filt_plots <- nc_no_filt%>%
  mutate(label = 'No filters')%>%
  bind_rows(nc_log_filt, nc_all_filt)%>%
  filter(activity != "recalcitrant")%>%
  mutate(activity2 = activity,
         network = 1700,
         num_features = 1)%>%
  group_by(network, num_features, label)%>%
  nest()%>%
  mutate(plots = map(data, ~ ggplot(.x, aes(nc, log2_change, color = activity)) +
                       geom_point(stat = 'identity') +
                       scale_color_manual(values = c('#78B7C5', '#EBCC2A')) +
                       new_scale('color') + 
                       geom_smooth(method = lm, aes(color = activity2)) +
                       scale_color_manual(values = c('#3B9AB2', '#E1AF00')) +
                       xlab("N:C") +
                       theme(
                         legend.position = "none",
                         plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
                         # plot.margin = margin(2,.8,2,.8, "cm"),
                         axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
                         axis.text.y = element_text(size = 20),
                         panel.background = element_rect(fill = "transparent"), # bg of the panel
                         plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                         panel.grid.major = element_blank(), # get rid of major grid
                         panel.grid.minor = element_blank(), # get rid of minor grid
                         legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                         legend.box.background = element_rect(fill = "transparent")
                       )))%>%
  select(-data)
    

# VIZUALIZATIONS -- NC log reg --------------------------------------------
nc_logreg <- log2_change_vals%>%
  inner_join(t_pvals%>%
               # filter(activity != 'recalcitrant')%>%
               # rename(feat_activity = activity)%>%
               select(feature_number, Organism, DayNight)%>%
               unique(), 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network, C:dG)%>%
              filter(network != '-1'),
            by = 'feature_number')%>%
  left_join(net_activity,
            by = c('network', 'DayNight'))%>%
  
  group_by(Organism, DayNight, Replicate)%>%
  mutate(nc = N/C,
         pc = P/C,
         sample_nc = T0*nc/sum(T0),
         sample_pc = T0*pc/sum(T0),
         activity2 = activity)%>%
  ungroup()%>%
  mutate(log_snc = log(sample_nc + 1),
         log_spc = log(sample_pc + 1))

pdf('./plots/nc_log2_logreg_112320.pdf', width = 15, height = 10)
nc_logreg%>%
  filter(activity == 'depletolite',
         nc > 0,
         !is.na(nc))%>%
  ggplot(aes(log_snc, log2_change, color = activity)) +
  geom_point(stat = 'identity') +
  scale_color_manual(values = c('#EBCC2A')) + # '#78B7C5', 
  new_scale('color') + 
  geom_smooth(method = 'lm', formula = y ~ x, aes(color = activity)) +
  scale_color_manual(values = c('#E1AF00')) + # '#3B9AB2', 
  labs(x = "Log(Weighted N:C)", y = bquote(Log[2] ~change)) +
  geom_text(aes(x = -6, y = 12, color = 'depletolite',
                label = paste("FDR p-value: ", (nc_loglin%>%
                                                  filter(DayNight == 'Day',
                                                         activity == 'depletolite'))$FDR%>%
                                formatC(format = "e", digits = 2), sep = "")), size = 7) +
  geom_text(aes(x = -6, y = 10.5, color = 'depletolite',
                label = paste("r²: ", (nc_loglin%>%
                                         filter(DayNight == 'Day',
                                                activity == 'depletolite'))$r2%>%
                                formatC(format = "e", digits = 2), sep = "")), size = 7) +
  # geom_text(aes(x = -6, y = 15, color = 'accumolite',
  #               label = paste("FDR p-value: ", (nc_loglin%>%
  #                                                 filter(DayNight == 'Day',
  #                                                        activity == 'accumolite'))$FDR%>%
  #                               formatC(format = "e", digits = 2), "***", sep = "")), size = 7) +
  # geom_text(aes(x = -6, y = 13.5, color = 'accumolite',
  #               label = paste("r²: ", (nc_loglin%>%
  #                                        filter(DayNight == 'Day',
  #                                               activity == 'accumolite'))$r2%>%
  #                               round(digits = 4), sep = "")), size = 7) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 25, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )

nc_logreg%>%
  filter(activity == 'depletolite',
         # pc > 0,
         # !is.na(pc)
         )%>%
  ggplot(aes(log_spc, log2_change, color = activity)) +
  geom_point(stat = 'identity') +
  scale_color_manual(values = c( '#EBCC2A')) + # '#78B7C5',
  new_scale('color') + 
  geom_smooth(method = 'lm', formula = y ~ x, aes(color = activity)) +
  scale_color_manual(values = c('#E1AF00')) + # '#3B9AB2', 
  labs(x = "Log(Weighted P:C)", y = bquote(Log[2] ~change)) +
  # geom_text(aes(x = -12, y = 12, color = 'depletolite',
  #               label = paste("FDR p-value: ", (pc_loglin%>%
  #                                                 filter(DayNight == 'Day',
  #                                                        activity == 'depletolite'))$FDR%>%
  #                               formatC(format = "e", digits = 2), "***", sep = "")), size = 7) +
  # geom_text(aes(x = -12, y = 10.5, color = 'depletolite',
  #               label = paste("r²: ", (pc_loglin%>%
  #                                        filter(DayNight == 'Day',
  #                                               activity == 'depletolite'))$r2%>%
  #                               round(digits = 4), sep = "")), size = 7) +
  # geom_text(aes(x = -12, y = 15, color = 'accumolite',
  #               label = paste("FDR p-value: ", (pc_loglin%>%
  #                                                 filter(DayNight == 'Day',
  #                                                        activity == 'accumolite'))$FDR%>%
  #                               formatC(format = "e", digits = 2), "***", sep = "")), size = 7) +
  # geom_text(aes(x = -12, y = 13.5, color = 'accumolite',
  #               label = paste("r²: ", (pc_loglin%>%
  #                                        filter(DayNight == 'Day',
  #                                               activity == 'accumolite'))$r2%>%
  #                               round(digits = 4), sep = "")), size = 7) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 25, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")
  )
dev.off()  

# VIZUALIZATIONS -- Checking for gaussian distribution in NOSC, Mass, NC and PC --------
# LilliFors and other tests for normallity will not be effective due to large sample size. 
# QQ-plots used to establish how 'normal' the data appear

log2_check <- (log_reg%>%
                 filter(response_var == 'NOSC',
                        activity != 'recalcitrant'))$log2_change

nc_check <- (nc_logreg%>%
               filter(activity != 'recalcitrant',
                      nc > 0,
                      !is.na(pc)))$sample_nc

lnc_check <- (nc_logreg%>%
               filter(activity != 'recalcitrant',
                      nc > 0,
                      !is.na(pc)))$log_snc

pc_check <- (nc_logreg%>%
               filter(activity != 'recalcitrant',
                      pc > 0,
                      !is.na(pc)))$sample_pc

lpc_check <- (nc_logreg%>%
                filter(activity == 'depletolite'))$log_spc


nosc_check <- (log_reg%>%
                 filter(response_var == 'NOSC',
                        activity != 'recalcitrant')%>%
                 mutate(NOSC = value))$NOSC

mass_check <- (log_reg%>%
                 filter(response_var == 'row m/z',
                        activity != 'recalcitrant')%>%
                 mutate(`row m/z` = value))$`row m/z`

pdf("./plots/gaussian_distribution_non_transformed.pdf", width = 15, height = 10)
car::qqPlot(nc_check, 
            ylab = "Weighted N:C quantiles", xlab = "Normal quantiles",
            main = 'QQ-plot: Weighted N:C')
car::qqPlot(lnc_check, 
            ylab = bquote(Log[10]~(Weighted ~N:C ~quantiles)), xlab = "Normal quantiles",
            main = 'QQ-plot: Log transformed Weighted N:C')
car::qqPlot(pc_check, 
            ylab = 'Weighted P:C quantiles', xlab = "Normal quantiles",
            main = 'QQ-plot: Weighted P:C')

car::qqPlot(lpc_check, 
            ylab = bquote(Log[10]~(Weighted ~P:C ~quantiles)), xlab = "Normal quantiles",
            main = 'QQ-plot: Log transformed Weighted P:C')

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
              select(-c(Replicate, T0, TF, complete_removal))%>%
              group_by(feature_number, Organism, DayNight)%>%
              summarize_if(is.numeric, mean),
            by = c("feature_number", "Organism", "DayNight"))%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  left_join(net_activity,
            by = c('network', 'DayNight'))

org_pie <- org_log2_ra%>%
  filter(network != -1)%>%
  group_by(Organism, DayNight, Replicate, Timepoint)%>%
  select(-c(log10:asin))%>%
  mutate(ra = xic/sum(xic),
         count = 1)%>%
  ungroup()%>%
  select(-Replicate)%>%
  group_by(Organism, feature_number, activity)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  group_by(Organism, activity)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(Organism)%>%
  mutate(end = 2 * pi * cumsum(ra)/sum(ra),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1),
         explosive = case_when(activity == 'depletolite' ~ 0.2,
                               TRUE ~ 0),
         texty = case_when(activity == 'depletolite' ~ 1.2*cos(middle),
                           activity == 'accumolite' ~ 1.2*cos(middle),
                           TRUE ~ 1.05*cos(middle)),
         textx = case_when(activity == 'accumolite' & Organism == 'Pocillopora verrucosa' ~ -10*sin(middle),
                           activity == 'accumolite' & Organism == 'CCA' ~ -3.5*sin(middle),
                           activity == 'accumolite' & Organism == 'Dictyota' ~ -5*sin(middle),
                           activity == 'accumolite' & Organism == 'Porites lobata' ~ -20*sin(middle),
                           activity == 'accumolite' & Organism == 'Turf' ~ -2.2*sin(middle),
                           activity == 'depletolite' & Organism == 'Turf' ~ 1.1*sin(middle),
                           activity == 'depletolite' ~ 1.4*sin(middle),
                           TRUE ~ 1.05*sin(middle)))%>%
  ungroup()

org_pie_vis <-org_pie%>%
  group_by(Organism)%>%
  nest()%>%
  mutate(data = map(data, ~ggplot(.x) +
                      geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                                       start = start, end = end, 
                                       fill = activity, explode = explosive, linetype = NA)) +
                      geom_text(aes(x = textx, y = texty, 
                                    label =  paste(round(ra, digits = 4)*100, "%", sep = ""),
                                    color = activity,
                                    hjust = hjust, vjust = vjust), size = 20, show.legend = FALSE) +
                      coord_fixed() +
                      scale_x_continuous(limits = c(-1.5, 1.5),  # Adjust so labels are not cut off
                                         name = "", breaks = NULL, labels = NULL) +
                      scale_y_continuous(limits = c(-1.2, 1.3),      # Adjust so labels are not cut off
                                         name = "", breaks = NULL, labels = NULL) +
                      theme_classic() +
                      ggtitle(Organism) + 
                      # facet_wrap(~ Organism) +
                      scale_fill_manual(values = c('#78B7C5', '#EBCC2A', "#00A08A")) +
                      scale_color_manual(values = c('#78B7C5', '#EBCC2A', "#00A08A")) + 
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
                            strip.background = element_blank())))

pdf("./plots/org_pie_111720.pdf", width = 17, height = 10)
org_pie_vis$data
  
dev.off()

# VIZUALIZATIONS -- Organism comparisons log2_change ----------------------
pdf("./plots/org_depletoliteebar_netactivity_120720.pdf", width = 15, height = 11)
org_log2_ra%>%
  inner_join(feature_stats_wdf%>%
               ungroup()%>%
               select(feature_number, Organism, DayNight),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  group_by(feature_number, Organism, DayNight, activity)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  group_by(Organism, DayNight, activity)%>%
  mutate(count = 1)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  filter(activity == 'depletolite')%>%
  ggplot(aes(Organism, xic)) +
  geom_bar(aes(fill = activity), stat = 'identity', position = 'stack') +
  scale_fill_manual(values = c('#EBCC2A'
    # '#78B7C5', '#EBCC2A', "#00A08A"
    )) +
  geom_text(aes(label = paste(count), vjust = -.2), size = 8) +
  labs(y = 'Sum depletolite intensity', fill = 'Network Activity: ') +
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


pdf("./plots/org_log2_120720.pdf", width = 15, height = 10)
org_log2_ra%>%
  filter(activity == 'depletolite')%>%
  group_by(feature_number, Organism, DayNight, activity)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  inner_join(feature_stats_wdf%>%
               ungroup()%>%
               select(feature_number, Organism, DayNight),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  ggplot(aes(Organism, log2_change)) +
  geom_boxplot(aes(color = '#EBCC2A')) +
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


# VIZUALIZATION -- Nutrition/lability vs microbial community change -------
fcm_T0_5 <- dorc_fcm_fdom%>%
  select(-c(1:4, 9:36, 38, 39))%>%
  filter(Timepoint == 'T0' | Timepoint == 'T5')%>%
  spread(Timepoint, `Cells µL-1`)%>%
  group_by(Organism, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
                   cells_ul = T5-T0)%>%
  select(-c(T0, T5))

lability_val <- log2_change_vals%>%
  inner_join(t_pvals%>%
               # filter(activity != 'depletolite')%>%
               select(feature_number, Organism, DayNight)%>%
               unique(), 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network, C:dG)%>%
              filter(network != '-1'),
            by = 'feature_number')%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`),
            by = 'feature_number')%>%
  left_join(net_activity,
            by = c('network', 'DayNight'))%>%
  left_join(net_glm%>%
              filter(DayNight == 'Day',
                     activity == 'depletolite',
                     response_var == 'row m/z')%>%
              select(DayNight, r2)%>%
              rename(mass_r2 = r2),
            by = 'DayNight')%>%
  left_join(nc_loglin%>%
              filter(activity == 'depletolite')%>%
              select(DayNight, r2)%>%
              rename(nc_r2 = r2),
            by = 'DayNight')%>%
  left_join(pc_loglin%>%
              filter(activity == 'depletolite')%>%
              select(DayNight, r2)%>%
              rename(pc_r2 = r2),
            by = 'DayNight')%>%
  left_join(fcm_T0_5%>%
              mutate(Replicate = as.numeric(Replicate)), 
            by = c('Organism', 'DayNight', 'Replicate'))%>%
  # filter(activity == 'depletolite')%>%
  select(-c(Experiment, TF))%>%
  group_by(Organism, DayNight, Replicate)%>%
  mutate(nc = N/C,
         pc = P/C,
         sample_nc = T0*nc,
         log_snc = log(sample_nc),
         sample_pc = T0*pc,
         log_spc = log(sample_pc),
         activity2 = activity,
         nc_ra = sample_nc/sum(T0),
         pc_ra = sample_pc/sum(T0),
         mass_x_r2 = `row m/z`/mean(`row m/z`)*mass_r2,
         nc_x_r2 = nc_ra*nc_r2,
         pc_x_r2 = pc_ra*pc_r2,
         lability = mass_x_r2  + pc_x_r2) #nc was taken out because not significant p value

sum_xic_x_log2 <- lability_val%>%
  filter(!is.na(activity))%>%
  # filter(activity == 'depletolite')%>%
  group_by(Organism, DayNight, Replicate)%>%
  summarize_if(is.numeric, median, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(Organism, Replicate)%>%
  mutate(sum_xic = sum(T0),
         change_x_xic = T0*log2_change)%>%
  ungroup()

# nutritional_stats <- lability_val%>%
#   select(-feature_number)%>%
#   filter(DayNight == 'Day')%>%
#   group_by(Organism, Replicate)%>%
#   summarize_if(is.numeric, mean, na.rm = TRUE)%>%
#   ungroup()%>%
#   mutate(exp = "d")%>%
#   group_by(exp)%>%
#   nest()%>%
#   mutate(model = map(data, ~ glm(lability ~ cells_ul, family = gaussian, .x)),
#          p_vals = map(model, ~tidy(.x)%>%
#                         filter(term == 'cells_ul')),
#          r2 = map(model, ~with(summary(.x), 1-deviance/null.deviance)))%>%
#   select(-c(data, model))%>%
#   unnest(p_vals)%>%
#   ungroup()%>%
#   mutate(r2 = as.numeric(r2),
#          r = sqrt(r2),
#          FDR = p.adjust(p.value, method = 'BH'))

# pdf("./plots/lability_nutrition_fcm.pdf", width = 15, height = 10)
# lability_val%>%
#   filter(DayNight == 'Day')%>%
#   ggplot(aes(cells_ul, lability)) +
#   geom_point(aes(color = Organism), stat = 'summary', fun.y = 'mean', size = 5) +
#   scale_color_manual(values = org_colors_no_water) + 
#   geom_smooth(method = 'lm') +
#   labs(x = bquote('Cells'~µL^-1), y = "Metabolite pool nutritional value") +
#   geom_text(aes(x = 425, y = 0.239,
#                 label = paste("p-value: ", nutritional_stats$p.value%>%
#                                 round(digits = 4), "*", sep = "")), size = 7) +
#   geom_text(aes(x = 425, y = 0.2375, 
#                 label = paste("r²: ", nutritional_stats$r2%>%
#                                 round(digits = 4), sep = "")), size = 7) +
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
# dev.off()


# VIZUALIZATIONS -- Lability value from multiple regressions --------------
# PN has been removed. Only N present or absent
# pn_mass_coe <- n_p_mulreg$coefficients["`row m/z`"]
# pn_nosc_coe <- n_p_mulreg$coefficients["NOSC"]
# pn_n_coe <- n_p_mulreg$coefficients["log_snc"]
# # pn_p_coe <- n_p_mulreg$coefficients["log_spc"]  not significant
# pn_intercept <- n_p_mulreg$coefficients["(Intercept)"]

# p_mass_coe <- p_mulreg$coefficients["`row m/z`"]
# p_nosc_coe <- p_mulreg$coefficients["NOSC"]
# p_p_coe <- p_mulreg$coefficients["log_spc"]
# p_intercept <- n_mulreg$coefficients["(Intercept)"] not significant 

# N present
n_mass_coe <- n_mulreg$coefficients["`row m/z`"]
n_nosc_coe <- n_mulreg$coefficients["NOSC"]
n_n_coe <- n_mulreg$coefficients["log_snc"]
n_intercept <- n_mulreg$coefficients["(Intercept)"]
# n_mano_coe <- n_mulreg$coefficients["`row m/z`:NOSC"]

#N absent
m_mass_coe <- mass_mulreg$coefficients["`row m/z`"]
m_nosc_coe <- mass_mulreg$coefficients["NOSC"]
m_intercept <- mass_mulreg$coefficients["(Intercept)"]

# Old four way model
# mul_reg_fcm <- lability_val%>%
#   ungroup()%>%
#   mutate(multiple_reg_lability = case_when(P > 0 & N > 0 ~ (`row m/z`*pn_mass_coe + NOSC*pn_nosc_coe + log_snc*pn_n_coe + pn_intercept),
#                                            P > 0 & N == 0 ~ (`row m/z`*p_mass_coe + NOSC*p_nosc_coe + log_spc*p_p_coe),
#                                            P == 0 & N > 0 ~ (`row m/z`*n_mass_coe + NOSC*n_nosc_coe + log_snc*n_n_coe + n_intercept),
#                                            P == 0 & N == 0 ~ (NOSC*m_nosc_coe + m_intercept),
#                                            TRUE ~ NA_real_))

mul_reg_fcm <- mul_reg%>%
  ungroup()%>%
  mutate(multiple_reg_lability = case_when(N > 0 ~ (`row m/z`*n_mass_coe + NOSC*n_nosc_coe + log_snc*n_n_coe + n_intercept),
                                           N == 0 ~ (`row m/z`*m_mass_coe + NOSC*m_nosc_coe + m_intercept),
                                           TRUE ~ NA_real_))%>%
  filter(DayNight == 'Day',
         NOSC < 0)%>%
  group_by(Organism, Replicate)%>%
  summarise_if(is.numeric, median, na.rm = TRUE)%>%
  ungroup()%>%
  left_join(fcm_T0_5%>%
              filter(DayNight == 'Day')%>%
              mutate(Replicate = as.numeric(Replicate)), 
            by = c('Organism', 'Replicate'))
  

#Linear model
lability_lm <- mul_reg_fcm%>%
  lm(cells_ul ~ multiple_reg_lability, data = .)

lab_p <- (lability_lm%>% 
            tidy()%>% 
            filter(term == 'multiple_reg_lability'))$p.value

lab_r2 <- summary(lability_lm)$adj.r.squared

#Plotting
pdf("./plots/lability_nultiple_regressions.pdf", width = 15, height = 10)
mul_reg_fcm%>%
  ggplot(aes(multiple_reg_lability, cells_ul)) +
  geom_point(aes(color = Organism), stat = 'identity', size = 5) +
  scale_color_manual(values = org_colors_no_water) + 
  geom_smooth(method = 'lm') +
  labs(y = bquote('Cells'~µL^-1), x = "Metabolite pool median modeled lability") +
  geom_text(aes(x = -1.37, y = 700,
                label = paste("p-value: ", lab_p%>%
                                formatC(format = "e", digits = 2), "***", sep = "")), size = 9) +
  geom_text(aes(x = -1.37, y = 670,
                label = paste("r²: ", lab_r2%>%
                                round(digits = 4), sep = "")), size = 9) +
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
  mutate(multiple_reg_lability = case_when(N > 0 ~ (`row m/z`*n_mass_coe + NOSC*n_nosc_coe + log_snc*n_n_coe + n_intercept),
                                           N == 0 ~ (`row m/z`*m_mass_coe + NOSC*m_nosc_coe + m_intercept),
                                           TRUE ~ NA_real_))%>%
  filter(DayNight == 'Day',
         NOSC < 0))$multiple_reg_lability


car::qqPlot(lability_val_check, 
            ylab = "Lability value", xlab = "Normal quantiles",
            main = 'QQ-plot: Lability value')

# Supplemental xic/log2 change dont explain fcm ---------------------------
pdf("./plots/supplemental_xiclog2_fcm.pdf", width = 15, height = 10)
sum_xic_x_log2%>%
  # filter(!is.na(activity))%>%
  filter(DayNight == 'Day')%>%
  ggplot(aes(change_x_xic, cells_ul)) +
  geom_point(aes(color = Organism), stat = 'summary', fun.y = 'mean', size = 5) +
  scale_color_manual(values = org_colors_no_water) + 
  geom_smooth(method = 'lm') +
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

sum_xic_x_log2%>%
  # filter(!is.na(activity))%>%
  filter(DayNight == 'Day',
         NOSC < 0)%>%
  ggplot(aes(log2_change, cells_ul)) +
  geom_point(aes(color = Organism), stat = 'summary', fun.y = 'mean', size = 5) +
  scale_color_manual(values = org_colors_no_water) + 
  geom_smooth(method = 'lm') +
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

sum_xic_x_log2%>%
  # filter(!is.na(activity))%>%
  filter(DayNight == 'Day')%>%
  ggplot(aes(sum_xic, cells_ul)) +
  geom_point(aes(color = Organism), stat = 'summary', fun.y = 'mean', size = 5) +
  scale_color_manual(values = org_colors_no_water) + 
  geom_smooth(method = 'lm') +
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

# VIZUALIZATIONS -- N/C to activity ~ Organism ----------------------------
unique_nc%>%
  # filter(activity == 'depletolite')%>%
  ggplot(aes(activity, N/C)) +
  # geom_point(stat = 'identity') +
  # geom_bar(stat = 'summary', fun.y = 'mean') +
  geom_boxplot() +
  # scale_color_manual(values = org_colors_no_water) +
  facet_wrap(~ Organism, scales = 'free_y') +
  theme(
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 15, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 20),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_line('grey'), # get rid of major grid
    panel.grid.minor = element_line('grey'), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    legend.position = 'top'
  )

nc_logreg%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), 
            by = 'feature_number')%>%
  filter(activity == 'depletolite')%>%
  filter(!is.na(N),
         P > 0)%>%
  ggplot(aes(Organism)) +
  geom_boxplot(aes(y = `row m/z`))
  # geom_boxplot(aes(y = log_snc))
  

# NC compare 21 and 131 ---------------------------------------------------
nc21_131 <- networking%>%
  filter(network == '131' | network == '21')%>%
  inner_join(t_pvals%>%
               filter(activity == 'depletolite',
                      DayNight == 'Day',
                      Organism == 'Pocillopora verrucosa' | Organism == 'Porites lobata')%>%
               select(feature_number, Organism),
             by = 'feature_number')%>%
  left_join(org_log2_ra%>%
              filter(Timepoint == 'T0',
                     DayNight == 'Day')%>%
              ungroup()%>%
              select(feature_number, Organism, xic, ra)%>%
              group_by(feature_number, Organism)%>%
              summarize_if(is.numeric, sum)%>%
              ungroup(),
            by = c('feature_number', 'Organism'))%>%
  group_by(Organism)%>%
  mutate(relative_nc = N/C*(xic/sum(xic)))

nc21_131%>%
  ggplot(aes(Organism, relative_nc)) +
  geom_boxplot()

# VIZUALIZATIONS -- Noise-removal effect ----------------------------------
noise_removal_pdf <- plots_all_filter%>%
  full_join(filter_steps_na, by = c('network', 'label', 'num_features'))%>%
  ungroup()%>%
  bind_rows(glm_all_filt, nc_filt_plots)%>%
  mutate(num_features = case_when(is.na(num_features) ~ 0,
                                  TRUE ~ as.numeric(num_features)),
         label = as.factor(label),
         label = fct_relevel(label, "No filters", "Timepoint filters"))%>%
  arrange(network, label)%>%
  ggarrange(plotlist = c(.$plots[2:nrow(.)], .$plots[1]), ncol = 3, nrow = 1,
            labels = .$label)

pdf("./plots/noise_removal_090120.pdf", width = 15, height = 10)
noise_removal_pdf
dev.off()

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
  left_join(net_activity, by = c('network', 'DayNight'))%>%
  inner_join(org_exometabolites, by = c('feature_number', 'Organism'))%>%
  inner_join(min_filter%>%
               ungroup()%>%
               filter(DayNight == 'Day')%>%
               select(feature_number, Organism)%>%
               unique(),
             by = c('feature_number', 'Organism'))

org_nosc%>%
  mutate(Organism2 = Organism)%>%
  ggplot(aes(NOSC, log2_change, color = Organism)) +
  facet_wrap(~activity) +
  geom_point() +
  geom_smooth(method = lm, aes(color = Organism2))


org_nosc%>%
  filter(activity != 'recalcitrant',
         !is.na(activity))%>%
  ggplot(aes(Organism, NOSC, color = activity)) +
  geom_boxplot()

org_nosc%>%
  filter(activity != 'recalcitrant',
         !is.na(activity))%>%
  ggplot(aes(activity, NOSC)) +
  geom_boxplot()

# VIZUALIZATIONS -- Cytoscape ---------------------------------------------
filtered_features <- feature_stats_wdf%>%
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
  left_join(net_activity,
            by = c('network', 'DayNight'))%>%
  rename(net_activity = activity)%>%
  left_join(t_pvals%>%
              select(feature_number, Organism, DayNight, activity)%>%
              rename(feature_activity = activity)%>%
              filter(feature_activity != 'recalcitrant')%>%
              group_by(feature_number, DayNight)%>%
              summarize_all(paste0, collapse = ", ")%>%
              ungroup()%>%
              rename(feature_organism = Organism),
            by = c('feature_number', 'DayNight'))%>%
  filter(DayNight == 'Day')%>%
  mutate(post_filtering = case_when(is.na(post_filtering) ~ 0,
                                    TRUE ~ as.numeric(post_filtering)),
         feature_activity = case_when(is.na(feature_activity) ~ 'recalcitrant',
                                            TRUE ~ as.character(feature_activity)),
         feature_organism = case_when(is.na(feature_organism) ~ '0',
                                      TRUE ~ as.character(feature_organism)))
  
write_csv(cyto_deplete, './analysis/cyto_depletes.csv')



# SUMMARY -- noise filtering steps ----------------------------------------
summary_no_background <- as.vector(feature_table_no_back_trans_filter%>%
                                     filter(background == "real"))$feature_number%>%
  length()

summary_no_back_trans <- as.vector(feature_table_no_back_trans_filter%>%
                                     filter(background == "real")%>%
                                     filter(Dorcierr_transient == "real"))$feature_number%>%
  length()                                    


filtering_feature_summary <- log2_change_vals%>%
  select(feature_number, Organism, DayNight)%>%
  unique()%>%
  mutate(DayNight2 = DayNight)%>%
  group_by(DayNight2)%>%
  nest()%>%
  mutate(log2_timepoint_bottleneck = map(data, ~ inner_join(.x, log2_features%>%
                                                              ungroup()%>%
                                                              select(feature_number, DayNight)%>%
                                                              unique(),
                                                            by = c('feature_number', 'DayNight'))%>%
                                           select(feature_number)%>%
                                           unique()%>%
                                           nrow()),
         xic_filter = map(data, ~ inner_join(.x, log2_features%>%
                                               ungroup()%>%
                                               select(feature_number, DayNight)%>%
                                               unique(),
                                             by = c('feature_number', 'DayNight'))%>%
                            inner_join(min_filter%>%
                                         ungroup()%>%
                                         select(feature_number, Organism, DayNight)%>%
                                         unique(),
                                       by = c('feature_number', 'Organism', 'DayNight'))%>%
                            mutate(xic_filter = 1)%>%
                            group_by(Organism)%>%
                            summarize_if(is.numeric, sum)),
         All_filters = map(data, ~ inner_join(.x, log2_features%>%
                                                ungroup()%>%
                                                select(feature_number, DayNight)%>%
                                                unique(),
                                              by = c('feature_number', 'DayNight'))%>%
                             inner_join(min_filter%>%
                                          ungroup()%>%
                                          select(feature_number, Organism, DayNight)%>%
                                          unique(),
                                        by = c('feature_number', 'Organism', 'DayNight'))%>%
                             inner_join(org_exometabolites, 
                                        by = c('feature_number', 'Organism', 'DayNight'))%>%
                             mutate(All_filters = 1)%>%
                             group_by(Organism)%>%
                             summarize_if(is.numeric, sum)%>%
                             add_row(Organism = "Water control", All_filters = 0)),
         one_s_ttest = map(data, ~ inner_join(.x, log2_features%>%
                                                ungroup()%>%
                                                select(feature_number, DayNight)%>%
                                                unique(),
                                              by = c('feature_number', 'DayNight'))%>%
                             inner_join(min_filter%>%
                                          ungroup()%>%
                                          select(feature_number, Organism, DayNight)%>%
                                          unique(),
                                        by = c('feature_number', 'Organism', 'DayNight'))%>%
                             inner_join(org_exometabolites, 
                                        by = c('feature_number', 'Organism', 'DayNight'))%>%
                             inner_join(t_pvals%>%
                                          filter(activity != 'recalcitrant')%>%
                                          select(feature_number, Organism, DayNight),
                                        by = c('feature_number', 'Organism', 'DayNight'))%>%
                             mutate(one_s_ttest = 1)%>%
                             group_by(Organism)%>%
                             summarize_if(is.numeric, sum)%>%
                             add_row(Organism = "Water control", one_s_ttest = 0)),
         filts = map2(xic_filter, All_filters, ~ left_join(.x, .y, by = 'Organism')),
         stats = map2(filts, one_s_ttest, ~ left_join(.x, .y, by = 'Organism')))%>%
  select(-c(data, xic_filter, All_filters, one_s_ttest, filts))%>%
  unnest(c(log2_timepoint_bottleneck, stats))%>%
  rename(DayNight = DayNight2)%>%
  mutate(Background_filter = 13843,
         Transient_filter = 10605)%>%
  select(DayNight, Organism, Background_filter, Transient_filter, everything())

write_csv(filtering_feature_summary, "./analysis/filtering_steps_summary.csv")


```