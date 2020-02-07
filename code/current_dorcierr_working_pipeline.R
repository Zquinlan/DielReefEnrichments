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

#PCoA, PERMANOVA
library(vegan)
library(ape)

#Visualizations
library(wesanderson)
library(RColorBrewer)
library(gplots)
library(gtable)

#Defining functions and removing issues with overlapping function calls
map <- purrr::map
select <- dplyr::select
tidy <- broom::tidy
rename <- dplyr::rename

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

inchikey_df <- read_csv("~/Documents/SDSU/Moorea_2017/csi_inchikey.csv")%>%
  mutate(feature_number = as.character(feature_number))

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
  left_join(nap_df, by = "feature_number", suffix = c("", "_nap"))%>%
  left_join(csi_finger_id, by = "feature_number", suffix = c("", "_csi"))%>%
  add_column(binary_ID = .$LibraryID, .before = 1)%>%
  mutate(binary_ID = case_when(binary_ID != "N/A" ~ "1",
                               Compound_NameAnalog_ != "NA" ~ "2",
                               !is.na(ConsensusSC) ~ "3",
                               TRUE ~ as.character(binary_ID)))%>%
  add_column(combined_ID = .$LibraryID, .before = 1)%>%
  mutate(combined_ID = case_when(binary_ID == "1" ~ LibraryID,
                                 binary_ID == "3" ~ ConsensusSC,
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
  dplyr::select(c(feature_number, network,combined_ID, binary_ID, 
                  canopus_annotation:CLASS_STRING, ZodiacMF, `characterization scores`,
                  C:dG, inchi_binary, inchi_combined))%>%
  separate(CLASS_STRING, c("level 1", "level 2", "level 3",
                           "level 4", "level 5", "level 6", "level 7", "level 8"), sep = ";")%>%
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
log2_features <- feature_table_no_back_trans%>%
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
  filter(log2_change >= 1 | log2_change <= -1)%>%
  ungroup()%>%
  select(-c(Organism, T0, TF))%>%
  group_by(DayNight, feature_number)%>%
  summarize_if(is.numeric, mean)


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
               select(-log2_change), by = c("feature_number", "DayNight"))

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
  mutate_if(is.numeric, log10)%>%
  filter(sample_name != "D_OF_1_T0N",
         sample_name != "D_IN_2_T0N",
         sample_name != "D_PL_3_TFN",
         sample_name != "D_TR_1_T0N",
         sample_name != "D_WA_2_T0D",
         sample_name != "D_WA_1_T0D",
         sample_name != "D_CC_1_T0D",
         sample_name != "D_CC_2_T0D")

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
  left_join(fcm_wdf%>%
              rename(sample_code = sample_name)%>%
              select(c(`Cells µL-1`, sample_code)), by = "sample_code")%>%
  group_by(sample_code)%>%
  mutate(sum = sum(reads),
         ra = reads/sum,
         cell_abun = ra*`Cells µL-1`,
         log10 = log10(cell_abun + 0.001))%>%
  ungroup()%>%
  group_by(OTU)%>%
  mutate(abundant = case_when(max(ra) > 0.01 | sum(ra > 0.001) >=3 ~ "abundant",
                              TRUE ~ "rare"))


microbe_no_rare <- microbe_combined%>%
  filter(abundant == "abundant")%>%
  select(-abundant)


# PRE-STATS CLEANING -- microbe RA data  --------------------------------------------
ra_bigger_TF <- microbe_no_rare%>%
  select(-c(reads, sum, log10, numOtus, sample_code, `Cells µL-1`, cell_abun))%>%
  spread(Organism, ra)%>%
  gather(Organism, ra, 6:10)%>%
  mutate(difference = ra - `Water control`)%>%
  filter(difference > 0)%>%
  dplyr::select(c(DayNight, OTU, Organism))

average_ra <- microbe_no_rare%>%
  filter(Timepoint == "TF")%>%
  select(-c(sample_code, Experiment, Replicate, numOtus, reads))%>%
  group_by(OTU, Organism, DayNight)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()


# SET SEED ----------------------------------------------------------------
set.seed(2005)


# STATS ANOVA -- Microbe TWO-Way ------------------------------------------
aov_microbe <- microbe_no_rare%>%
  filter(Timepoint == "TF")%>%
  group_by(OTU)%>%
  nest()%>%
  mutate(anova = map(data, ~ aov(log10 ~ Organism*DayNight, .x)%>%
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
  mutate(activity = case_when(FDR_greater < 0.05 ~ "accumulite",
                              FDR_lesser < 0.05 ~ "depletolite",
                              FDR_lesser >= 0.05 & FDR_greater >= 0.05 | 
                                is.na(FDR_lesser) & is.na(FDR_greater) ~ "recalcitrant"))

# # STATS -- ANOVA metabolties Organism T0 -------------------------------------------------------
anova_dom_t0_df <- t_pvals%>%
  select(c(DayNight, feature_number, activity))%>%
  left_join(dom_stats_wdf, by = c("feature_number", "DayNight"))%>%
  filter(Timepoint == "T0")

anova_dom_t0 <- anova_dom_t0_df%>%
  group_by(feature_number, DayNight, activity)
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


# PRE-POST-HOC CLEANING -- Microbe Dunnetts and DayNight anova -------------------------------
mic_organism_post_hoc <- microbe_no_rare%>%
  filter(OTU %in% organism_significant_microbes)%>%
  filter(Timepoint == "TF")

daynight_microbe_post_hoc <- microbe_no_rare%>%
  filter(OTU %in% DayNight_significant_microbes)%>%
  filter(Timepoint == "TF")

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
  mutate(sum = sum(log10))%>%
  filter(sum != 0)%>%
  dplyr::select(-sum)%>%
  mutate(Organism = factor(Organism))%>%
  mutate(Organism = fct_relevel(Organism, organism_order_micro))%>%
  nest()%>%
  mutate(dunnett = map(data, ~ aov(log10 ~ Organism, .x)%>%
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
  mutate(sum = sum(log10))%>%
  filter(sum != 0)%>%
  dplyr::select(-sum)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(log10 ~ DayNight, .x)%>%
                      tidy()))%>%
  unnest(data)%>%
  dplyr::select(-c(4:7))%>%
  filter(term != "Residuals")%>%
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
  right_join(t_pvals, by = c("DayNight", "feature_number"))%>%
  left_join(networking, by = "feature_number")%>%
  filter(activity == "depletolite")%>%
  group_by(DayNight, feature_number, Organism)
  
major_deplete$max_log <- apply(major_deplete[3:8], 1, min)

major_depletolites <- major_deplete%>%
  filter(max_log < -3.3)%>%
  filter(Organism == "Pocillopora verrucosa" | Organism == "Dictyota" | Organism == "CCA")

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

# META-STATS Random Forest -- microbes ------------------------------------------
rf_microbe_prep <- microbe_no_rare%>%
  filter(Timepoint == "TF")%>%
  select(c(1:6, OTU, log10))%>%
  filter(OTU %in% micro_sigs_vector)%>%
  # group_by(Organism, Replicate, Timepoint, DayNight, OTU)%>%
  # summarize_if(is.numeric, sum)%>%
  # ungroup()%>%
  spread(OTU, log10)%>%
  select(Organism, DayNight, 7:ncol(.))%>%
  mutate(Organism = as.factor(Organism))

names(rf_microbe_prep) <- make.names(names(rf_microbe_prep))

rf_microbe <- rf_microbe_prep%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(data = map(data, ~randomForest(Organism ~ ., .x,
                                        importance = TRUE, proximity = TRUE,
                                        ntree = 5000, na.action=na.exclude)),
         mda = map(data, ~ .x$importance%>%
                     as.data.frame()%>%
                     rownames_to_column("otu")))

rf_microbe_mda <- rf_microbe%>%
  select(DayNight, mda)%>%
  unnest(mda)


pdf("~/Documents/GitHub/DORCIERR/data/plots/microbe_mda.pdf", width = 7, height = 5)
rf_microbe_mda%>%
  filter(DayNight == "Day")%>%
  ggplot(., aes(x= reorder(otu, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  ggtitle("Day") +
  geom_hline(yintercept = (top_n(rf_microbe_mda%>%
                                   filter(DayNight=="Day"), 30, MeanDecreaseAccuracy)%>%
                             arrange(-MeanDecreaseAccuracy))$MeanDecreaseAccuracy[30],
             col = "red") +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), 
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic")) +
  scale_x_discrete(breaks = seq(0, 568, 50))

rf_microbe_mda%>%
  filter(DayNight == "Night")%>%
  ggplot(., aes(x= reorder(otu, -MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_point(stat = "identity") +
  ggtitle("Night") +
  geom_hline(yintercept = (top_n(rf_microbe_mda%>%
                                   filter(DayNight=="Night"), 35, MeanDecreaseAccuracy)%>%
                             arrange(-MeanDecreaseAccuracy))$MeanDecreaseAccuracy[30],
             col = "red") +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), 
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    legend.text = element_text(face = "italic")) +
  scale_x_discrete(breaks = seq(0, 568, 50))
dev.off()

rf_microbe_sigs <- (rf_microbe_mda%>%
                      group_by(DayNight)%>%
                      nest()%>%
                      mutate(dic = map(data, ~ top_n(.x, 20, Dictyota)%>%
                                         select(otu)),
                             cca = map(data, ~ top_n(.x, 20, CCA)%>%
                                         select(otu)),
                             trf = map(data, ~ top_n(.x, 20, Turf)%>%
                                         select(otu)),
                             poc = map(data, ~ top_n(.x, 20, `Pocillopora verrucosa`)%>%
                                         select(otu)),
                             por = map(data, ~ top_n(.x, 20, `Porites lobata`)%>%
                                         select(otu)),
                             wat = map(data, ~ top_n(.x, 20, `Water control`)%>%
                                         select(otu)))%>%
                      select(-data)%>%
                      gather(species, importance, dic:wat)%>%
                      unnest(importance))$otu%>%
  as.vector()%>%
  unique()
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
  select(-c(T0, TF))%>%
  filter(Organism != "Turf",
         Organism != "Porites lobata",
         Organism != "Water control")%>%
  left_join(t_pvals, by = c("Organism", "DayNight", "feature_number"))%>%
  mutate(fdr_combined = case_when(FDR_greater > FDR_lesser ~ FDR_lesser,
                                  FDR_lesser > FDR_greater ~ FDR_greater,
                                  TRUE ~ 1))
volcano_themes <- function(x) {
  ggplot(x, aes(log2_change, -log10(fdr_combined))) +
  geom_point(stat = "identity", size = 2.5, shape = 1, alpha = 0.8) +
  scale_x_continuous(breaks= seq(-20, 18, 2)) +
  scale_y_continuous(breaks = seq(0, 10, 1)) +
  # scale_shape_manual(values = c(1,19)) +
  # scale_color_manual(values = c("#FF0000", "#50A45C", "#C49647", "#5BBCD6")) +
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

volcano%>%
  volcano_themes() +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  geom_vline(xintercept = 3.3, col = "red", linetype = "dashed") +
  geom_vline(xintercept = -3.3, col = "red", linetype = "dashed")
dev.off()  

# GRAPHING -- [OSM] Major depletolites ------------------------------------------
major_depletolites%>%
  # left_join(dom_stats_wdf%>%
  #             filter(Timepoint == "T0")
              # group_by(Organism, DayNight, feature_number)%>%
              # summarize_if(is.numeric, mean)%>%
              # ungroup()
            # , by = c("DayNight", "feature_number", "Organism"))%>%
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
              spread(Timepoint, xic)%>%
              group_by(Organism, DayNight, feature_number)%>%
              summarize_if(is.numeric, mean, na.rm = TRUE)%>%
              mutate(log2_change = log2(TF/T0))%>%
              ungroup()%>%
              select(-c(T0, TF)), by = c("DayNight", "feature_number", "Organism"))%>%
  mutate(log2_change = as.numeric(as.character(log2_change)))%>%
  filter(`characterization scores` == "Good")%>%
  ggplot(aes(simplified_makeup, y = -log2_change, fill = Organism)) +
  geom_bar(stat = "summary", fun.y = "mean") +
  # scale_color_manual(values = wes_palette("Zissou1", 30, type = "continuous")) +
  facet_wrap(~DayNight) +
  # geom_boxplot() +
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




# GRAPHING -- [OSM] Important OTUs --------------------------------------------------------
osm_otus <- dunnett_microbe_pvals%>%
  left_join(microbe_taxonomy, by = "OTU")%>%
  filter(Organism != "Turf",
         Organism != "Porites lobata")%>%
  left_join(average_ra, by = c("OTU", "Organism", "DayNight"))%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "sp"), sep = ";")%>%
  unite(Tax_plot, c("Family", "Genus", "OTU"), sep = " ", remove = FALSE)


osm_otus%>%
  filter(ra >= 0.015)%>%
  ggplot(aes(x = Organism, y = ra, col = Tax_plot, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~DayNight) +
  # scale_color_manual(values = "black") +
  scale_fill_manual(values = c("#FF0000", "#32806E", "#91A737", "#F49C00", "#C49647", "#5BBCD6")) +
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
  ) +
  ggtitle("OTUs")

# GRAPHING -- [OSM] FCM data ----------------------------------------------------
fcm_graphing <- fcm_wdf%>%
  filter(!Organism == 'Influent',
         !Organism == 'Offshore',
         Organism != "Turf",
         Organism != "Porites lobata")

pdf("~/Documents/GitHub/DORCIERR/data/plots/osm_fcm_DayNight.pdf", width = 7, height = 5)
fcm_graphing%>%
  ggplot(aes(x= Timepoint, y = `Cells µL-1`, color = Organism))+
  geom_point(stat = "summary", fun.y = "mean") +
  geom_line(aes(group = Organism), stat = "summary", fun.y = "mean") +
  facet_wrap(~ DayNight) +
  scale_color_manual(values = c("#FF0000", "#50A45C", "#F69100", "#5BBCD6")) +
  scale_y_continuous(limits = c(0,900), breaks= seq(0, 900, 100)) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.x = element_blank(), # get rid of major grid
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )
dev.off()


# GRAPHING -- [OSM] PCoAs Labile ------------------------------------------
osm_dom_pco <- dom_stats_wdf%>%
  spread(Timepoint, log)%>%
  filter(Organism != "Turf",
         Organism != "Porites lobata")%>%
  group_by(feature_number, Organism, DayNight)%>%
  mutate(mean_t0 = mean(T0, na.rm = TRUE))%>%
  mutate(change = TF - mean_t0)%>%
  ungroup()%>%
  mutate(zscore = ((change - mean(change, na.rm = TRUE))/sd(change, na.rm = TRUE)))%>%
  mutate(zscore = zscore + 78)%>%
  select(-c(TF, T0, mean_t0, change))%>%
  unite(sample, c(1:4), sep = "_")%>%
  spread(2,3)%>%
  column_to_rownames(var = "sample")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()

## Plot Eigenvalues
osm_dom_pco$values[1:10,]%>%
  as.data.frame()%>%
  rownames_to_column("Axis")%>%
  mutate(axis = as.numeric(Axis))%>%
  ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, color = "red", vjust = -0.5)

## PCoA plot
pdf("./plots/osm_dom_pcoa.pdf", width = 6, height = 5)
osm_dom_pco$vectors%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample")%>%
  separate(sample, c("experiemnt", "Organism", "replicate", "DayNight"), sep = "_")%>%
  ggplot(., aes(x = Axis.1, y = Axis.2, color = Organism, shape = DayNight)) +
  geom_point(stat = "identity", aes(size = 0.2)) +
  scale_shape_manual(values = c(1,19)) +
  scale_color_manual(values = c("#FF0000", "#50A45C", "#F69100", "#5BBCD6")) +
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
  filter(Timepoint == "TF",
         Organism != "Turf",
         Organism != "Porites lobata")%>%
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
  scale_shape_manual(values = c(1,19)) +
  scale_color_manual(values = c("#FF0000", "#50A45C", "#F69100", "#5BBCD6")) +
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
## split between day and night
# dorc_split <- dom_pco%>%
#   group_by(DayNight)%>%
#   nest()%>%
#   mutate(data = map(data, ~ unite(.x, sample, c(1:3), sep = "_")%>%
#                              spread(2,3)%>%
#                              column_to_rownames(var = "sample")%>%
#                              vegdist(na.rm = TRUE)%>%
#                              pcoa()),
#          vectors = map(data, ~ .x$vectors%>%
#                          as.data.frame()%>%
#                          rownames_to_column(var = "sample")))
# 
# dorc_split_vectors <- dorc_split%>%
#   dplyr::select(-data)%>%
#   unnest(vectors)%>%
#   mutate(feature = DayNight)%>%
#   separate(sample, c("experiment", "Organism", "replicate", "Timepoint"), sep = "_")
# 
# 
# pco_all <- bind_rows(pco_scores_all, dorc_split_vectors)%>%
#   mutate(dorc_shape = case_when(DayNight == "Day" ~ 1,
#                                 DayNight == "Night" ~ 19,
#                                 TRUE ~ 3))
# 
# pdf("pcoa_all.pdf", width = 6, height = 5)
# pco_all%>%
#   split(.$feature)%>%
#   map(~ ggplot(., aes(x = Axis.1, y = Axis.2, color = Organism, shape = DayNight)) +
#         geom_point(stat = "identity", aes(size = 0.2)) +
#         scale_shape_manual(values = .$dorc_shape) +
#         ggtitle(.$feature) +
#         theme(
#           panel.background = element_rect(fill = "transparent"), # bg of the panel
#           plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#           panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
#           panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
#           legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#           legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
#           legend.text = element_text(face = "italic")) +
#         xlab("Axis 1") +
#         ylab("Axis 2"))
# dev.off()


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

summary_dunnett <- dom_dunnetts%>%
  right_join(dom_stats_wdf%>%
              select(c(feature_number, Timepoint, DayNight, Organism, log))%>%
              group_by(feature_number,Timepoint, DayNight, Organism)%>%
              summarize_if(is.numeric, mean, na.rm = TRUE), by = c("feature_number", "DayNight", "Organism"))%>%
  left_join(networking, by = "feature_number")%>%
  select(-c(rhs:p.value))

summary_count_dunnett <- dom_dunnetts%>%
  select(activity, Organism, DayNight)%>%
  group_by(activity, Organism, DayNight)%>%
  mutate(count = 1)%>%
  summarize_if(is.numeric, sum)

# summary_ttest <- t_pvals%>%
#   mutate(depletolites = map(data, ~ filter(.x, activity == "depletolites")$feature_number%>%
#                               length()),
#          accumlites = map(data, ~ filter(.x, activity == "accumulites")$feature_number%>%
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

write_csv(poc_deplete, "./analysis/pocillopora_depletolites.csv")

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


# SUMMARY -- Microbe filtering --------------------------------------------
microbe_summary <- microbe_combined%>%
  select(-c(Experiment, Organism, Replicate, Timepoint, DayNight, reads, asin, numOtus, sum))%>%
  spread(sample_code, ra)


