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
                                TRUE ~ "Night"))%>%
    spread(Timepoint, xic)%>%
    group_by(Organism, DayNight, feature_number)%>%
    summarize_if(is.numeric, mean, na.rm = TRUE)%>%
    ungroup()%>%
    mutate(log2_change = log2(TF/T0),
           log2_change = case_when(TF > 0 & T0 == 0 ~ 6.6,
                                   TF == 0 & T0 > 0 ~ -6.6,
                                   TF == 0 & T0 == 0 ~ 0,
                                   log2_change < -6.6 ~ -6.6,
                                   log2_change > 6.6 ~ 6.6,
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
                                 log2_change < -6.6 ~ -6.6,
                                 log2_change > 6.6 ~ 6.6,
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
  filter(network != "-1")%>%
  group_by(network, DayNight)%>%
  nest()%>%
  mutate(greater = map(data, ~ t.test(log10 ~ Timepoint, .x, alternative = "greater")),
         lesser = map(data, ~ t.test(log10 ~ Timepoint, .x, alternative = "less")))%>%
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
         pc = P/C,
         oc = O/C,
         hc = H/C,
         log_nc = log10(nc),
         log_pc = log10(pc),
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
  select(log2_change, `row m/z`, NOSC, log_snc, log_soc, log_shc)%>%
  cor()%>%
  corrplot::corrplot()
dev.off()

n_mulreg <- mul_reg%>%
  filter(N > 0 & O > 0,
         NOSC < 0)%>%
  group_by(feature_number, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  lm(log2_change ~ `row m/z`+ NOSC + log_snc + log_shc + log_soc, data = .)

o_mulreg <- mul_reg%>%
  filter(N == 0 & O > 0,
         NOSC < 0)%>%
  group_by(feature_number, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  lm(log2_change ~ `row m/z` + NOSC + log_soc + log_shc, data = .)

# mass_mulreg <- mul_reg%>%
#   filter(N == 0, O == 0)%>%
#   group_by(feature_number, Organism)%>%
#   summarize_if(is.numeric, mean)%>%
#   ungroup()
#   lm(log2_change ~ `row m/z`, data = .)

sink("./analysis/multiple_reg_coefficients_whole_metabolome_N_O_models.txt")
summary(n_mulreg)
summary(o_mulreg)
# summary(mass_mulreg)
sink()

# VIZUALIZATIONS -- RGB hex codes for orgs --------------------------------
org_colors_no_water <- c("#CC0033","#669900", "#CC6600", "#9900FF", "#33CC33")

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
  left_join(net_activity, by = c('network', 'DayNight'))%>%
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
  left_join(net_activity,
            by = c('network', 'DayNight'))%>%
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
  inner_join(feature_stats_wdf%>%
               ungroup()%>%
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
         activity2 = activity)

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


# VIZUALIZATIONS -- Lability value from multiple regressions --------------
# N present  
n_mass_coe <- n_mulreg$coefficients["`row m/z`"]
# n_nosc_coe <- n_mulreg$coefficients["NOSC"]
n_n_coe <- n_mulreg$coefficients["log_nc"]
n_o_coe <- n_mulreg$coefficients["log_soc"]
n_h_coe <- n_mulreg$coefficients["log_shc"]
# n_xic_coe <- n_mulreg$coefficients["logxic"]
n_intercept <- n_mulreg$coefficients["(Intercept)"]


#N absent
o_mass_coe <- o_mulreg$coefficients["`row m/z`"]
o_nosc_coe <- o_mulreg$coefficients["NOSC"]
o_o_coe <- o_mulreg$coefficients["log_soc"]
o_h_coe <- o_mulreg$coefficients["log_shc"]
# o_xic_coe <- o_mulreg$coefficients["logxic"]
o_intercept <- o_mulreg$coefficients["(Intercept)"]

#nosc absent
# m_mass <- mass_mulreg$coefficients["`row m/z`"]
# m_intercept <- mass_mulreg$coefficients["(Intercept)"]


mul_reg_fcm <- mul_reg%>%
  ungroup()%>%
  filter(NOSC < 0)%>%
  mutate(multiple_reg_lability = case_when(N > 0 & O > 0 ~ (`row m/z`*n_mass_coe + log_snc*n_n_coe + log_soc*n_o_coe + log_shc*n_h_coe + n_intercept),
                                           N == 0 & O > 0 ~ (`row m/z`*o_mass_coe + NOSC*o_nosc_coe + log_soc*o_o_coe +log_shc*o_h_coe + o_intercept),
                                           TRUE ~ NA_real_),
         model_num = case_when(N > 0 & O > 0 ~ 'Nitrogen',
                               N == 0 & O > 0 ~ 'Oxygen',
                               TRUE ~ 'none'))%>%
  filter(DayNight == 'Day')%>%
  group_by(Organism, Replicate)%>%
  summarise_if(is.numeric, median, na.rm = TRUE)%>%
  ungroup()%>%
  left_join(fcm_T0_5%>%
              filter(DayNight == 'Day'), 
            by = c('Organism', 'Replicate'))


weighted_lability <- feature_stats_wdf%>%
  filter(Timepoint == 'T0',
         DayNight == 'Day')%>%
  left_join(mul_reg%>%
              ungroup()%>%
              mutate(modeled_lability = case_when(N > 0 & O > 0 ~ (`row m/z`*n_mass_coe + log_snc*n_n_coe + log_soc*n_o_coe + log_shc*n_h_coe + n_intercept),
                                                  N == 0 & O > 0 ~ (`row m/z`*o_mass_coe + NOSC*o_nosc_coe + log_soc*o_o_coe +log_shc*o_h_coe + o_intercept),
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

#Linear model
lability_lm <- mul_reg_fcm%>%
  lm(cells_ul ~ multiple_reg_lability, data = .)

lab_p <- (lability_lm%>% 
            tidy()%>% 
            filter(term == 'multiple_reg_lability'))$p.value

lab_r2 <- summary(lability_lm)$adj.r.squared

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
pdf("./plots/weighted_lability_nultiple_regressions.pdf", width = 15, height = 10)
weighted_lability%>%
  ggplot(aes(weighted_lability, cells_ul)) +
  geom_point(aes(color = Organism), stat = 'identity', size = 5) +
  geom_errorbar(aes(ymin = cells_ul - y_err, ymax = cells_ul + y_err)) +
  geom_errorbarh(aes(xmin = weighted_lability - x_err, xmax = weighted_lability + x_err))+
  scale_color_manual(values = org_colors_no_water) + 
  geom_smooth(method = 'lm') +
  labs(y = bquote('Cells'~µL^-1), x = "Metabolite pool modeled lability") +
  # geom_text(aes(x = -2.2, y = 820,
  #               label = paste("p-value: ", wlab_p%>%
  #                               formatC(format = "e", digits = 2), sep = "")), size = 9) +
  # geom_text(aes(x = -2.2, y = 780,
  #               label = paste("F statistic: ", wlab_f%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = -2.2, y = 740,
  #               label = paste("r²: ", wlab_r2%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = -2.2, y = 700,
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


# Depletolite Networks ----------------------------------------------------
molnet_class_network <- molnet_class%>% 
  left_join(networking%>% 
              select(feature_number, network), by = 'feature_number')%>% 
  select(-feature_number)%>% 
  unique()%>%
  filter(network != -1)

deplete_nets <- mul_reg%>%
  left_join(net_activity, by = c('DayNight', 'network'))%>%
  ungroup()%>%
  mutate(modeled_lability = case_when(N > 0 ~ (`row m/z`*n_mass_coe + NOSC*n_nosc_coe + log_snc*n_n_coe + n_intercept),
                                           N == 0 ~ (`row m/z`*m_mass_coe + NOSC*m_nosc_coe),
                                           TRUE ~ NA_real_))%>%
  filter(activity == 'depletolite')%>%
  select(-c(Replicate:TF, complete_removal, Experiment, `row m/z`:log_spc))%>%
  group_by(feature_number, network, Organism, activity)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  group_by(Organism, network, activity)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  left_join(molnet_class_network, by = 'network')%>%
  select(Organism, network, molnet_string, modeled_lability, everything())
  
write_csv(deplete_nets, './analysis/deplete_networks.csv')  


