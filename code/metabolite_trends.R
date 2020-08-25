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
  group_by(feature_number, Organism, DayNight)%>%
  filter(max(xic) >= 5*10^6)%>%
  left_join(networking%>%
              select(network, feature_number), by = 'feature_number')

min_filter_feature_tag <- min_filter%>%
  ungroup()%>%
  filter(DayNight == 'Day',
         network %in% c('131','21','7','55','107','198','466', '165', '619', '627', '756', '924', '1314','80','141'))%>%
  select(feature_number, Organism)%>%
  mutate(xic_min_filter = 1)

no_min_filter <- min_filter_pre%>%
  left_join(networking%>%
              select(network, feature_number), by = 'feature_number')

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

#Making the plots
plots_no_filter <- no_min_filter%>%
  group_by(feature_number, Timepoint, Organism)%>%
  filter(DayNight == 'Day',
         network %in% c('131','21','7','55','107','198','466', '165', '619', '627', '756', '924', '1314','80','141'))%>%
  left_join(min_filter_feature_tag, by = c('feature_number', 'Organism'))%>%
  mutate(mean = mean(xic),
         sd = sd(xic),
         xic_min_filter = case_when(is.na(xic_min_filter) ~ 'Noise',
                                    TRUE ~ 'Real'))%>%
  select(-c(xic, Replicate))%>%
  unique()%>%
  ungroup()%>%
  left_join(num_features_no_filter, by = 'network')%>%
  arrange(network)%>%
  mutate(color_b = 'black')%>%
  group_by(network, num_features)%>%
  nest()%>%
  mutate(plots = map(data, ~ ggplot(.x, aes(Timepoint, mean, fill = xic_min_filter, color = feature_number)) +
                       geom_bar(stat = 'identity', size = 0.0001) +
                       facet_wrap(~Organism) +
                       ggtitle(sprintf('Network %3.0f          total unique features = %3.0f', network, num_features)) +
                       scale_fill_manual(values = c("#F2AD00", "#00A08A")) +
                       ylab("Mean feature XIC") +
                       scale_color_manual(values = .x$color_b) +
                       theme(legend.position = 'none',
                             axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
                             axis.text.y = element_text(size = 20),
                             plot.title = element_text(size = 20),
                             panel.background = element_rect(fill = "transparent"), # bg of the panel
                             plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                             panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
                             panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"))))%>%
  select(-data)


plots_min_filter <- min_filter%>%
  group_by(feature_number, Timepoint, Organism)%>%
  filter(DayNight == 'Day',
         network %in% c('131','21','7','55','107','198','466', '165', '619', '627', '756', '924', '1314','80','141'))%>%
  mutate(mean = mean(xic),
         sd = sd(xic))%>%
  select(-c(xic, Replicate))%>%
  unique()%>%
  ungroup()%>%
  left_join(num_features_min_filter, by = 'network')%>%
  arrange(network)%>%
  mutate(color_b = 'black',
         fill_bl = "#00A08A")%>%
  group_by(network, num_features)%>%
  nest()%>%
  mutate(plots = map(data, ~ ggplot(.x, aes(Timepoint, mean, color = feature_number, fill = fill_bl)) +
                       geom_bar(stat = 'identity', size= 0.0001) +
                       facet_wrap(~Organism) +
                       ggtitle(sprintf('Network %3.0f          total unique features = %3.0f', network, num_features)) +
                       ylab("Mean feature XIC") +
                       scale_color_manual(values = .x$color_b) +
                       scale_fill_manual(values = .x$fill_bl) +
                       theme(legend.position = 'none',
                             axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
                             axis.text.y = element_text(size = 20),
                             plot.title = element_text(size = 20),
                             panel.background = element_rect(fill = "transparent"), # bg of the panel
                             plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                             panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
                             panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"))))%>%
  select(-data)

pdf("plots/networks_xic_notfiltered.pdf", width = 15, height = 10)
plots_no_filter$plots
dev.off()

pdf("plots/networks_xic_filtered.pdf", width = 15, height = 10)
plots_min_filter$plots
dev.off()

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
  inner_join(min_filter%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
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
  inner_join(min_filter%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
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
  inner_join(min_filter%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
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
  inner_join(min_filter%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
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
  filter(num_organism == 1)%>%
  left_join(log2_change_vals%>%
              filter(DayNight == 'Day')%>%
              select(feature_number, Organism, T0, TF, log2_change, Replicate)%>%
              rename(unique_organism = Organism),
            by = c('feature_number', 'unique_organism'))%>%
  rename('Organism' = 'unique_organism')%>%
  inner_join(org_exometabolites, by = c('feature_number', 'Organism'))%>%
  inner_join(min_filter%>%
               ungroup()%>%
               filter(DayNight == 'Day')%>%
               select(feature_number, Organism)%>%
               unique(),
             by = c('feature_number', 'Organism'))%>%
  group_by(Organism, Replicate)%>%
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
  group_by(Organism, Replicate)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  mutate(Replicate = as.factor(Replicate))%>%
  ggplot(aes(N, P, color = Organism, shape = Replicate)) +
  geom_point(stat = 'identity') +
  scale_color_manual(values = wes_palette('Darjeeling1', 5, type = 'continuous'))

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


