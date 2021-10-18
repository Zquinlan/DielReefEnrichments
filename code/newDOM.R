# LOADING -- packages -------------------------------------------------------
#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)
library(readxl)
library(multcomp)
library(CHNOSZ)
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

# Canopus 2021
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

# EcoNet
ecoNet <- read_csv('~/Documents/Github/DORCIERR/data/raw/metabolomics/ecoNetConsensus.csv')%>%
  select(-1)

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
  select(c(1, 'C', 'H', 'N', 'O', 'P', 'S'))%>%
  mutate(NOSC = -1*((4*C + H - 3*N - 2*O + 5*P - 2*S)/C)+4,
         dG = 60.3-28.5*NOSC)%>%
  filter(NOSC < 4 & NOSC > -4)%>%
  filter(P <= 2)

super_computer_annotations <- canopusSirius4%>%
  rename(feature_number = scan)%>%
  full_join(sirius_zodiac_anotations, by = "feature_number")%>%
  full_join(networking_energy, by = "feature_number")%>%
  mutate(`characterization scores` = `quality`)%>%
  mutate(`characterization scores` = case_when(`characterization scores` != "Good" ~ "Bad quality",
                                               ZodiacScore < .98 ~ "Low probability",
                                               TRUE ~ as.character(`characterization scores`)))

# CLEANING -- METADATA and filter out bad samples --------------------------
## join library hits, analog hits and super computer predictions
metadata <- node_info%>%
  full_join(super_computer_annotations, by = "feature_number")%>%
  full_join(feature_table_raw[1:4], by = 'feature_number')%>%
  full_join(true_hits, by = 'feature_number')%>%
  full_join(analog_hits, by = "feature_number", suffix = c("Library_", "Analog_"))%>%
  add_column(binary_ID = .$LibraryID, .before = 1)%>%
  add_column(combined_ID = .$LibraryID, .before = 1)%>%
  mutate(binary_ID = case_when(binary_ID != "N/A" ~ "1",
                               Compound_NameAnalog_ != "NA" ~ "2",
                               !is.na(superclass) ~ "3",
                               TRUE ~ as.character(binary_ID)),
         combined_ID = case_when(binary_ID == "1" ~ LibraryID,
                                 binary_ID == "3" ~ superclass,
                                 binary_ID == "2" ~ Compound_NameAnalog_,
                                 binary_ID == "N/A" ~ superclass,
                                 TRUE ~ as.character(binary_ID)),
         inchi_binary = case_when(!is.na(INCHILibrary_) ~ "1",
                                  !is.na(INCHIAnalog_) & is.na(INCHILibrary_) | 
                                    !is.na(INCHIAnalog_) & INCHILibrary_ == "N/A" ~ "2",
                                  TRUE ~ "3"),
         inchi_combined = case_when(inchi_binary == "1" ~ INCHILibrary_,
                                    inchi_binary == "2" ~ INCHIAnalog_,
                                    TRUE ~ "N/A"))


metadata$feature_number <- as.character(metadata$feature_number)

networking <- metadata%>%
  select(c(feature_number, network,combined_ID, binary_ID,  SmilesLibrary_, SmilesAnalog_,
           subclass:superclass, ZodiacMF, `characterization scores`,
           C:dG, inchi_binary, inchi_combined))%>%
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
                                       TRUE ~ as.character(simplified_makeup)))


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
  # filter(DayNight == 'Day')%>%
  nest(data = everything())%>%
  mutate(log2 = map(data, ~inner_join(.x, log2_features%>%
                                        select(DayNight, feature_number), by = c('DayNight', 'feature_number'))%>%
                      select(feature_number, DayNight)%>%
                      unique()%>%
                      group_by(DayNight)%>%
                      mutate(count = 1)%>%
                      summarize_if(is.numeric, sum)),
         exometabolite = map(data, ~inner_join(.x, log2_features%>%
                                                 select(DayNight, feature_number), by = c('DayNight', 'feature_number'))%>%
                               inner_join(min_filter%>%
                                            select(feature_number, DayNight, Organism), by = c('DayNight', 'Organism', 'feature_number'))%>%
                               select(feature_number, DayNight)%>%
                               unique()%>%
                               group_by( DayNight)%>%
                               mutate(count = 1)%>%
                               summarize_if(is.numeric, sum)),
         min = map(data, ~ inner_join(.x, log2_features%>%
                                        select(DayNight, feature_number), by = c('DayNight', 'feature_number'))%>%
                     inner_join(min_filter%>%
                                  select(feature_number, DayNight, Organism), by = c('DayNight', 'Organism', 'feature_number'))%>%
                     inner_join(org_exometabolites, by = c('DayNight', 'Organism', 'feature_number'))%>%
                     select(feature_number, DayNight)%>%
                     unique()%>%
                     group_by(DayNight)%>%
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
             by = c('Organism', 'DayNight', 'feature_number'))%>%
    group_by(Organism, Replicate, Timepoint, DayNight)%>%
    mutate(log10 = log10(xic + 1),
           ra = xic/sum(xic),
           asin = asin(ra))

feature_stats_wdf <- dom_stats_wdf%>%
  inner_join(org_exometabolites, by = c('feature_number', 'Organism', 'DayNight'))%>%
  group_by(Organism, Replicate, Timepoint, DayNight)%>%
  mutate(log10 = log10(xic + 1),
         ra = xic/sum(xic),
         asin = asin(ra))

# SET SEED ----------------------------------------------------------------
set.seed(2005)


# STATS  -- T-TEST Network level ------------------------------------------
net_test <- dom_stats_wdf%>%
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

singlenode_test <- dom_stats_wdf%>%
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



# STATS -- Tukey production -----------------------------------------------
production_ttest <- feature_stats_wdf%>%
  filter(Timepoint == 'T0')%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  # separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  ungroup()%>%
  filter(activity == 'depletolite')

# NEWLY PRODUCED DOM -- Subtrcting Water --------------------------------------------------------
all_features <- feature_stats_wdf%>%
  ungroup()%>%
  select(feature_number, DayNight)%>%
  unique()

water_features <- dorc_wdf%>%
  filter(!Organism == "Influent",
         !Organism == "Offshore",
         Timepoint == "T0" | Timepoint == "TF")%>%
  select(feature_number, Organism:ncol(.))%>%
  inner_join(log2_features%>%
               select(feature_number, DayNight), by = c('DayNight', 'feature_number'))%>%
  filter(Organism == 'Water control')%>%
  right_join(all_features, by = c('feature_number', 'DayNight'))%>%
  select(-Organism)%>%
  mutate(xic = log10(xic + 1))%>%
  rename(water = xic)%>%
  mutate(water = case_when(is.infinite(water) ~ 0,
                           TRUE ~ water))

new_production <- feature_stats_wdf%>%
  select(-c(xic, ra:ncol(.)))%>%
  left_join(water_features, by = c('feature_number', 'Timepoint', 'DayNight', 'Replicate'))%>%
  mutate(netProduction = log10 - water,
         netProduction = case_when(netProduction < 0 ~ 0,
                                   TRUE ~ netProduction))


# VIZUALIZATIONS -- RGB hex codes for orgs --------------------------------
org_colors_no_water <- c("#A30029","#669900", "#FF850A", "#9900FF", "#33CC33")


# VIZUALIZATIONS --  Net Production pie charts -----------------------------------------------
new_prod_pie_df <- new_production%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(Timepoint == 'T0')%>%
  group_by(activity, Organism, Replicate, DayNight)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(activity, Organism, DayNight)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  group_by(Organism, DayNight)%>%
  mutate(ra = netProduction/sum(netProduction))

new_prod_pies <- new_prod_pie_df%>%
  # ungroup()%>%
  # group_by(Organism, Replicate, activity)%>%
  # summarize_if(is.numeric, sum)%>%
  # ungroup()%>%
  # group_by(Organism, Replicate)%>%
  # mutate(ra = netProduction/sum(netProduction))%>%
  # ungroup()%>%
  # group_by(Organism, activity)%>%
  # mutate(err = sd(ra))%>%
  # summarize_if(is.numeric, mean)%>%
  # ungroup()%>%
  # group_by(Organism)%>%
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
                           activity == 'recalcitrant' & Organism == 'Pocillopora verrucosa' ~ 0.5*cos(middle),
                           activity == 'recalcitrant' & Organism == 'Dictyota' ~ 1.05*cos(middle),
                           activity == 'recalcitrant' & Organism == 'Turf' ~ 0.5*cos(middle),
                           activity == 'depletolite' & Organism == 'Dictyota' ~ 0.95*cos(middle),
                           activity == 'depletolite' ~ 1.2*cos(middle),
                           activity == 'accumolite' ~ 1.2*cos(middle),
                           TRUE ~ 1.2*cos(middle)),
         textx = case_when(activity == 'accumolite' & Organism == 'Pocillopora verrucosa' ~ -20*sin(middle),
                           activity == 'recalcitrant' & Organism == 'Pocillopora verrucosa' ~ 1*sin(middle),  # Without Error
                           activity == 'accumolite' & Organism == 'CCA' ~ -5*sin(middle),
                           activity == 'depletolite' & Organism == 'CCA' ~ 1.6*sin(middle),
                           activity == 'recalcitrant' & Organism == 'CCA' ~ 1*sin(middle),
                           # activity == 'accumolite' & Organism == 'Dictyota' ~ -200*sin(middle),   #with Error
                           activity == 'accumolite' & Organism == 'Dictyota' ~ -5*sin(middle),    #Without error
                           activity == 'depletolite' & Organism == 'Dictyota' ~ 2*sin(middle),
                           activity == 'accumolite' & Organism == 'Porites lobata' ~ sin(middle),
                           activity == 'depletolite' & Organism == 'Porites lobata' ~ 1.6*sin(middle),
                           activity == 'recalcitrant' & Organism == 'Porites lobata' ~ 1.6*sin(middle), # Without error
                           # activity == 'accumolite' & Organism == 'Turf' ~ -60*sin(middle),  #with error
                           activity == 'accumolite' & Organism == 'Turf' ~ -10*sin(middle),   #without error
                           activity == 'depletolite' & Organism == 'Turf' ~ 1.8*sin(middle),
                           activity == 'recalcitrant' & Organism == 'Turf' ~ 1*sin(middle),
                           activity == 'depletolite' ~ 1.4*sin(middle),
                           activity == 'recalcitrant' ~ -3*sin(middle),
                           TRUE ~ 1.05*sin(middle)))%>%
  ungroup()

new_prod_noerr <- new_prod_pies%>%
  filter(DayNight == 'Day')%>%
  group_by(Organism)%>%
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
                      ggtitle(Organism) + 
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

pdf("plots/org_pie_newlyProduced.pdf", width = 17, height = 12)
new_prod_noerr$data
dev.off()

newLabilityStack <- new_production%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(Timepoint == 'T0')%>%
  group_by(activity, Organism, Replicate, DayNight)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  select(-c(log10, network, net_act))%>%
  gather(source, value, water:netProduction)%>%
  mutate(activity = case_when(source == 'water' ~ 'background',
                              TRUE ~ activity),
         activity = as.factor(activity),
         activity = fct_relevel(activity, 'accumolite', 'depletolite', 'recalcitrant', 'background'))
  # group_by(Replicate, DayNight)%>%
  # mutate(ra = netProduction/sum(netProduction))%>%
  # ungroup()%>%
  # group_by(activity, DayNight)%>%
  # mutate(std = sd(ra))%>%
  # summarize_if(is.numeric, mean)

average_newLability <- new_production%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(Timepoint == 'T0')%>%
  group_by(activity, Replicate, DayNight)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(Replicate, DayNight)%>%
  mutate(ra = netProduction/sum(netProduction))%>%
  ungroup()%>%
  group_by(activity, DayNight)%>%
  mutate(std = sd(ra))%>%
  summarize_if(is.numeric, mean)

pdf('plots/stackedLabilityBackground.pdf', width = 15, height = 10)
newLabilityStack%>%
  ggplot(aes(Organism, value, fill = activity)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_manual(values = c('#78B7C5','#EBCC2A', '#006658', 'grey')) +
  labs(y = bquote(Sum ~Intensity(log[10] ~XIC)), fill = 'Metabolite pool') +
  gen_theme()
dev.off()  

# Vizualizations -- Lability classifications of DOM -----------------------
lability_classes <- new_production%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  # filter(DayNight == 'Day')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  spread(Timepoint, netProduction)%>%
  mutate(val = case_when(activity == 'recalcitrant' ~ T0,
                         activity == 'depletolite' ~ T0,
                         activity == 'accumolite' ~ T0))%>%
  group_by(activity, DayNight)%>%
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

percent_lability_summary <- lability_classes%>%
  group_by(DayNight, superclass_consensus, Replicate)%>%
  mutate(relativeT0Prod = T0/sum(T0, na.rm = TRUE))%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, activity)%>%
  filter(sum(T0) != 0)%>%
  ungroup()

depletolite_percent <- percent_lability_summary%>%
  group_by(DayNight, superclass_consensus, activity)%>%
  mutate(std = sd(relativeT0Prod, na.rm = TRUE))%>%
  summarize_if(is.numeric, mean)%>%
  select(DayNight, superclass_consensus, activity, relativeT0Prod, std)%>%
  filter(activity == 'depletolite')

nitrogen_percent <- percent_lability_summary%>%
  filter(superclass_consensus %like% '%nitrogen%')%>%
  group_by(DayNight, activity, superclass_consensus)%>%
  mutate(std = sd(relativeT0Prod, na.rm = TRUE))%>%
  summarize_if(is.numeric, mean)%>%
  select(DayNight, activity, superclass_consensus, relativeT0Prod, std)

pdf('plots/superclassPercentProductionDayNight.pdf', width = 15, height = 10)
lability_classes%>%
  group_by(DayNight, superclass_consensus, Replicate)%>%
  mutate(relativeT0Prod = T0/sum(T0, na.rm = TRUE),
         relativeT0Prod = case_when(is.na(relativeT0Prod) ~ 0,
                                    TRUE ~ relativeT0Prod))%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, activity)%>%
  filter(sum(T0) != 0)%>%
  ungroup()%>%
  # filter(superclass_consensus != 'NA')%>%
  # unite(superclass_class, c(superclass_consensus, class_consensus), sep = '_')%>%
  # filter(ecoNetConsensus != 'NA')%>%
  # mutate(superclass_consensus = case_when(is.na(superclass_consensus) ~ 'No Consensus',
  #                                         TRUE ~ superclass_consensus))%>%
  ggplot(aes(superclass_consensus, relativeT0Prod*100, color = activity)) +
  geom_point(stat = 'identity', size = 4) +
  coord_flip() +
  scale_color_manual(labels = c('Accumolite', 'Depletolite', 'Recalcitrant'), values = c('#78B7C5', '#EBCC2A', "#006658")) +
  # geom_bar(stat = 'identity') +
  # geom_errorbar(aes(ymax = total_production + std, ymin = total_production - std)) +
  # facet_wrap(~activity, labeller = labeller(activity = label_wrap_gen()), scales = 'free_y', nrow = 3) +
  labs(y = 'Percent production (%)', x = 'Putative Superclass Annotation', color = 'Activity') +
  facet_wrap(~DayNight) +
  # scale_y_log10() +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        axis.text.y = element_text(size = 14))

dev.off()

org_new_prod <- new_production%>%
  select(-c(log10, water))%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(Timepoint == 'T0',
         activity == 'depletolite')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = 'feature_number')%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), sep = ';')%>%
  group_by(superclass_consensus, Organism, Replicate, DayNight)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  # filter(DayNight == 'Day')%>%
  group_by(superclass_consensus)%>%
  mutate(val = zscore(netProduction))%>%
  unite(sample, c('DayNight', 'Organism', 'Replicate'), sep = '_')

pdf('plots/organismDepletoliteSuperclassHeatmapDay.pdf', width = 15, height = 10)
org_new_prod%>%
  ggplot(aes(sample, superclass_consensus, fill = val)) +
  geom_tile() +
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major.y = element_line(colour = 'transparent'), # get rid of major grid
        panel.grid.major.x = element_line(colour = 'transparent'),
        axis.text = element_text(size = 25),
        axis.text.x = element_text(angle = 60, hjust = 1))


org_new_prod%>%
  separate(sample, c('DayNight', 'Organism', 'Replicate'), sep = '_')%>%
  group_by(DayNight, superclass_consensus)%>%
  filter(sum(netProduction) != 0)%>%
  mutate(std = sd(netProduction))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  ggplot(aes(superclass_consensus, netProduction)) +
  geom_errorbar(aes(ymax = netProduction + std, ymin = netProduction-1)) +
  geom_bar(stat = 'identity') +
  coord_flip() + 
  labs(y = bquote(Sum ~Intensity ~(log[10] ~XIC)), x = 'Putative Superclass Annotation') +
  facet_wrap(~DayNight, nrow = 1) +
  gen_theme()

org_new_prod%>%
  separate(sample, c('DayNight', 'Organism', 'Replicate'), sep = '_')%>%
  group_by(DayNight, superclass_consensus, Organism)%>%
  filter(sum(netProduction) != 0)%>%
  ungroup()%>%
  ggplot(aes(superclass_consensus, netProduction, color = Organism)) +
  # geom_errorbar(aes(ymax = netProduction + std, ymin = netProduction-1)) +
  geom_point(stat = 'identity', size = 4) +
  coord_flip() + 
  labs(y = bquote(Sum ~Intensity ~(log[10] ~XIC)), x = 'Putative Superclass Annotation') +
  facet_wrap(~DayNight, nrow = 1) +
  scale_color_manual(values = org_colors_no_water) +
  gen_theme() +
  theme(axis.text.y = element_text(size = 4))

org_new_prod%>%
  filter(superclass_consensus == 'Organic nitrogen compounds')%>%
  separate(sample, c('DayNight', 'Organism', 'Replicate'), sep = '_')%>%
  ggplot(aes(DayNight, netProduction)) +
  geom_boxplot() +
  labs(y = bquote(Sum ~nitrogen ~Intensity ~(log[10] ~XIC)), x = 'Diel production cycle') +
  gen_theme()

dev.off()


# VIZUALIZATIONS --  Lability Classes elemental ratios --------------------
elemental_classifications <- new_production%>%
  ungroup()%>%
  select(feature_number, Organism, DayNight)%>%
  unique()%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  # filter(activity == 'depletolite')%>%
  left_join(metadata%>%
              select(feature_number, C, H, N, O, P, NOSC), by = 'feature_number')%>%
  mutate(hc = H/C,
         oc = O/C,
         pc = P/C,
         nc = N/C)
  

pdf('plots/elementalRatiosClassifications.pdf', width = 15, height = 10)
elemental_classifications%>%
  filter(nc < 1,
         activity == 'depletolite')%>%
  select(-c(C:P))%>%
  gather(element, ratio, NOSC:nc)%>%
  filter(element != 'hc')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = 'feature_number')%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), sep = ';')%>%
  ggplot(aes(superclass_consensus, ratio, color = DayNight)) +
  # geom_boxplot() +
  geom_point(position = 'jitter') +
  coord_flip() +
  facet_wrap(~element, nrow = 1, scales = 'free_x') +
  scale_color_manual(values = c('orange', 'black')) +
  gen_theme()

elemental_classifications%>%
  filter(nc < 1)%>%
  ggplot(aes(DayNight, nc)) +
  geom_boxplot() +
  labs(x = 'Diel production cycle', y = 'N:C ratio') +
  gen_theme()

# elemental_classifications%>%
#   filter(nc < 1)%>%
#   select(-c(C:P))%>%
#   gather(element, ratio, NOSC:nc)%>%
#   filter(element != 'hc')%>%
#   left_join(ecoNet%>%
#               rename(feature_number = scan)%>%
#               mutate(feature_number = as.character(feature_number)), by = 'feature_number')%>%
#   separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), sep = ';')%>%
#   ggplot(aes(superclass_consensus, ratio, color = activity, shape = DayNight)) +
#   # geom_boxplot() +
#   geom_point(position = position_dodge2(width = 1)) +
#   coord_flip() +
#   facet_wrap(~element, nrow = 1, scales = 'free_x') +
#   scale_color_manual(values = c('#78B7C5','#EBCC2A', "#006658")) +
#   scale_shape_manual(values = c('\U263C', '\U263D')) +
#   gen_theme()
dev.off()

elementProduction_classifications <- new_production%>%
  ungroup()%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(Timepoint == 'T0')%>%
  left_join(metadata%>%
              select(feature_number, C, H, N, O, P, NOSC), by = 'feature_number')%>%
  mutate(oc = O/C,
         pc = P/C,
         nc = N/C,
         ocProd = netProduction*oc,
         pcProd = netProduction*pc,
         ncProd = netProduction*nc)%>%
  filter(nc < 1)%>%
  select(-c(C:nc))%>%
  gather(element, ratio, ocProd:ncProd)%>%
  filter(element != 'hc')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = 'feature_number')%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), sep = ';')
  
## Looks at approximate element production
elementProduction_classifications%>%
  filter(activity == 'depletolite')%>%
  ggplot(aes(superclass_consensus, ratio, color = DayNight)) +
  geom_boxplot() +
  # geom_point(position = 'jitter') +
  coord_flip() +
  facet_wrap(~element, nrow = 1, scales = 'free_x') +
  scale_color_manual(values = c('orange', 'black')) +
  gen_theme()

## elemental production compare acitivty
elementProduction_classifications%>%
  filter(DayNight == 'Day')%>%
  # filter(activity == 'depletolite')%>%
  ggplot(aes(superclass_consensus, ratio, color = activity)) +
  # geom_boxplot() +
  geom_point(position = 'jitter') +
  coord_flip() +
  facet_wrap(~element, nrow = 1, scales = 'free_x') +
  # scale_color_manual(values = c('orange', 'black')) +
  gen_theme()


# STATS -- Nitrogen and elemental ratio ttest -----------------------------
nitrogen_ttest <- org_new_prod%>%
  filter(superclass_consensus %like% '%nitrogen%')%>%
  separate(sample, c('DayNight', 'Organism', 'Replicate'), sep = '_')%>%
  ungroup()%>%
  t.test(netProduction ~ DayNight, data = .)

nitrogen_ratio_ttest <- elemental_classifications%>%
  filter(nc < 1)%>%
  t.test(nc ~ DayNight, data = .)
  

# VIZUALIZATIONS -- Metabolite pool XIC change ----------------------------
metabolitePool_xic_change <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  group_by(Organism, Timepoint, activity, Replicate, DayNight)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  select(-c(log10:asin))%>%
  group_by(Organism, activity, Timepoint, DayNight)%>%
  mutate(err = sd(xic))%>%
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

metabolitePool_meanChange <- metabolitePool_xic_change%>%
  select(-c(err:ra))%>%
  group_by(org_activity, Organism, activity, DayNight)%>%
  mutate(x_val = max(xic, na.rm = TRUE))%>%
  ungroup()%>%
  group_by(org_activity, Organism, activity, Timepoint, DayNight)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  spread(Timepoint, xic)%>%
  mutate(difference = TF-T0,
         x_val = case_when(Organism == 'Pocillopora verrucosa' ~ x_val - 0.05e10,
                           TRUE ~ x_val))%>%
  filter(activity == 'depletolite',
         DayNight == 'Night')


png("plots/change_xic_activityNight.png", width = 2000, height = 1200)
# cairo_pdf("plots/change_xic_activityDay.pdf", width = 20, height = 12)
metabolitePool_xic_change%>%
  filter(activity == 'depletolite',
         DayNight == 'Night')%>%
  mutate(log10 = log10(xic))%>%
  ggplot(aes(xic, Organism, color = Organism, shape = shape)) +
  geom_point(stat = 'identity', size = 18, alpha = 0.45) +
  # geom_errorbarh(aes(xmin = xic - err, xmax = xic + err)) +
  geom_line(aes(group = Organism)) +
  geom_text(data = metabolitePool_meanChange, aes(x = x_val, y = Organism, 
                                        label = formatC(difference, format = 'e', digits = 2), 
                                        shape = NULL), vjust = -1, size = 12) +
  # facet_wrap(~activity, nrow = 2) +
  scale_color_manual(values = org_colors_no_water) +
  scale_shape_manual(values = c("\u25A0", "\u25C4")) +
  labs(x = 'Sum Intensity (xic)', y = 'Organism', color = 'Organism: ', shape = 'Sample:') +
  # xlim(0,0.3) +
  # scale_x_log10() +
  # facet_wrap(~activity, nrow = 3) +
  theme(axis.line = element_blank(),
        axis.text = element_text(size = 30),
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
         shape = guide_legend(nrow=2, byrow=TRUE)) +
  xlim(0, 2.1e10)
dev.off()


# STATS -- depletolite pool xic mean change -------------------------------
xic_change_ttest <- metabolitePool_xic_change%>%
  filter(activity == 'depletolite')%>%
  select(-c(network:ra))%>%
  spread(Timepoint, xic)%>%
  group_by(Organism, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         difference = TF-T0,
         log10Difference = log10(abs(difference)))%>%
  ungroup()%>%
  group_by(Organism)%>%
  nest()%>%
  mutate(DayDepletion = map(data, ~ t.test(difference ~ DayNight, .x, alternative = "greater")["p.value"][[1]]),
         NightDepletion = map(data, ~ t.test(difference ~ DayNight, .x, alternative = "less")["p.value"][[1]]))%>%
  select(-data)%>%
  unnest(c(DayDepletion, NightDepletion))%>%
  gather(test, p, 2:3)%>%
  mutate(fdr = p.adjust(p, method = 'BH'))


# VIZUALIZATIONS -- FCM ---------------------------------------------------
fcm_graphing <- fcm_wdf%>%
  filter(!Organism == 'Influent',
         !Organism == 'Offshore')%>%
  mutate(Hours = case_when(Timepoint == "T0" ~ 0,
                           Timepoint == "T1" ~ 4,
                           Timepoint == "T2" ~ 8,
                           DayNight == 'Day' & Timepoint == "T3" ~ 16,
                           DayNight == 'Night' & Timepoint == "T3" ~ 17,
                           Timepoint == "T4" ~ 24,
                           DayNight == 'Day' & Timepoint == "T5" ~ 32,
                           DayNight == 'Night' & Timepoint == "T5" ~ 28,
                           DayNight == 'Day' & Timepoint == "T6" ~ 40,
                           DayNight == 'Night' & Timepoint == "T6" ~ 41,
                           Timepoint == "TF" ~ 48))%>%
  group_by(Organism, DayNight, Timepoint)%>%
  mutate(st_err = sd(`Cells µL-1`))%>%
  summarize_if(is.numeric, mean)

pdf("plots/FCM_Day.pdf", width = 15, height = 10)
fcm_graphing%>%
  filter(DayNight == 'Day')%>%
  ggplot(aes(x= Hours, y = `Cells µL-1`, color = Organism))+
  geom_point(stat = "identity", size = 5) +
  geom_errorbar(aes(ymin = `Cells µL-1` - st_err, ymax = `Cells µL-1` + st_err)) +
  geom_line(aes(group = Organism)) +
  scale_color_manual(values = c(org_colors_no_water, "#3B9AB2")) +
  labs(y = bquote(Cells ~µL^-1)) +
  facet_wrap(~ DayNight) +
  # scale_color_manual(values = c("darkorchid3", "#50A45C", "#AF814B", "#5BBCD6")) +
  scale_y_continuous(limits = c(0,900), breaks= seq(0, 900, 100)) +
  scale_x_continuous(limits = c(0,50), breaks = seq(0, 50, 5)) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 30),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.x = element_blank(), # get rid of major grid
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(), # get rid of minor grid
    # legend.position = 'top',
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

fcm_graphing%>%
  filter(DayNight == 'Night')%>%
  ggplot(aes(x= Hours, y = `Cells µL-1`, color = Organism))+
  geom_point(stat = "identity", size = 5) +
  geom_errorbar(aes(ymin = `Cells µL-1` - st_err, ymax = `Cells µL-1` + st_err)) +
  geom_line(aes(group = Organism)) +
  scale_color_manual(values = c(org_colors_no_water, "#3B9AB2")) +
  labs(y = bquote(Cells ~µL^-1)) +
  facet_wrap(~ DayNight) +
  # scale_color_manual(values = c("darkorchid3", "#50A45C", "#AF814B", "#5BBCD6")) +
  scale_y_continuous(limits = c(0,900), breaks= seq(0, 900, 100)) +
  scale_x_continuous(limits = c(0,50), breaks = seq(0, 50, 5)) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 30),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.x = element_blank(), # get rid of major grid
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(), # get rid of minor grid
    # legend.position = 'top',
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

# VIZUALIZATIONS -- Sum depletion growth rate first 24 hours --------------
metabolitePool_depletoliteChange <- metabolitePool_xic_change%>%
  filter(activity == 'depletolite')%>%
  select(-c(err, ra))%>%
  group_by(Organism, activity, Timepoint, DayNight, Replicate)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  spread(Timepoint, xic)%>%
  group_by(Organism, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         difference = TF-T0)

fcm_23 <- dorc_fcm_fdom%>%
  select(-c(1:4, 9:36, 38, 39))%>%
  spread(Timepoint, `Cells µL-1`)%>%
  mutate(hour = 24,
         final_cells = T4)%>%
  group_by(Organism, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         cells_ul = (log(final_cells) - log(T0))/(hour))%>%
  select(-c(T0:TF, hour))%>%
  left_join(metabolitePool_depletoliteChange,  by = c('Organism', 'DayNight', 'Replicate'))%>%
  mutate(log10 = log10(abs(difference)))

# Day lm
fcm_depletolite_lm <- fcm_23%>%
  filter(DayNight == 'Day')%>%
  lm(cells_ul ~ log10, data = .)

fcd_p <- (fcm_depletolite_lm%>% 
            tidy()%>% 
            filter(term == 'log10'))$p.value

fcd_f <- (fcm_depletolite_lm%>% 
            tidy()%>% 
            filter(term == 'log10'))$statistic

fcd_slope <- fcm_depletolite_lm$coefficients["log10"]
fcd_intercept <- fcm_depletolite_lm$coefficients["(Intercept)"]


fcd_r2 <- summary(fcm_depletolite_lm)$adj.r.squared

# Night lm
fcm_depletolite_lm_night <- fcm_23%>%
  filter(DayNight == 'Night')%>%
  lm(cells_ul ~ log10, data = .)

fcn_p <- (fcm_depletolite_lm_night%>% 
            tidy()%>% 
            filter(term == 'log10'))$p.value

fcn_f <- (fcm_depletolite_lm_night%>% 
            tidy()%>% 
            filter(term == 'log10'))$statistic

fcn_slope <- fcm_depletolite_lm_night$coefficients["log10"]
fcn_intercept <- fcm_depletolite_lm_night$coefficients["(Intercept)"]


fcn_r2 <- summary(fcm_depletolite_lm_night)$adj.r.squared


# Plotting
# Day
pdf('plots/fcm_sumDepletion.pdf', width = 15, height = 10)
fcm_23%>%
  filter(DayNight == 'Day',
         Organism != 'Influent',
         Organism != 'Water control',
         Organism != 'Offshore')%>%
  ggplot(aes(log10, cells_ul)) +
  geom_point(aes(color = Organism), size = 7) +
  geom_smooth(method = 'lm', se = FALSE, linetype = 3) +
  ylim(0.036,0.083) +
  # geom_text(aes(x = 9.6, y = 0.08,
  #               label = paste("p-value: ", fcd_p%>%
  #                               formatC(format = "e", digits = 2), sep = "")), size = 9) +
  # geom_text(aes(x = 9.6, y = 0.078,
  #               label = paste("F statistic: ", fcd_f%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = 9.6, y = 0.076,
  #               label = paste("r²: ", fcd_r2%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = 9.6, y = 0.074,
  #               label = paste("Cells µL^-1", " = ", fcd_slope%>%
  #                               round(digits = 2), "*lability + ", fcd_intercept%>%
  #                               round(digits = 2), sep = "")), size = 9) +
  scale_color_manual(values = org_colors_no_water) +
  gen_theme() +
  labs(y = 'Specific Growth Rate', x = bquote(Depletolite ~Total ~Depletion ~(log[10] ~XIC)))

#Night
fcm_23%>%
  filter(DayNight == 'Night',
         Organism != 'Influent',
         Organism != 'Water control',
         Organism != 'Offshore')%>%
  ggplot(aes(log10, cells_ul)) +
  geom_point(aes(color = Organism), size = 7) +
  geom_smooth(method = 'lm') +
  ylim(0.036,0.083) +
  # geom_text(aes(x = 9.6, y = 0.08,
  #               label = paste("p-value: ", fcn_p%>%
  #                               formatC(format = "e", digits = 2), sep = "")), size = 9) +
  # geom_text(aes(x = 9.6, y = 0.078,
  #               label = paste("F statistic: ", fcn_f%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = 9.6, y = 0.076,
  #               label = paste("r²: ", fcn_r2%>%
  #                               round(digits = 4), sep = "")), size = 9) +
  # geom_text(aes(x = 9.6, y = 0.074,
  #               label = paste("Cells µL^-1", " = ", fcn_slope%>%
  #                               round(digits = 2), "*SumIntensity + ", fcd_intercept%>%
  #                               round(digits = 2), sep = "")), size = 9) +
  scale_color_manual(values = org_colors_no_water) +
  gen_theme() +
  labs(y = 'Specific Growth Rate', x = bquote(Depletolite ~Total ~Depletion ~(log[10] ~XIC)))

dev.off()


# STATS -- Specific Growth Rate reduction ---------------------------------
microbialLoad_stats <- fcm_23%>%
  filter(!Organism %like any% c('Influent', 'Water control', 'Offshore'))%>%
  group_by(Organism)%>%
  nest()%>%
  mutate(data = map(data, ~t.test(final_cells ~ DayNight, .x)["p.value"][[1]]))%>%
  unnest(data)%>%
  mutate(test = 'load')

sgr_stats <- fcm_23%>%
  filter(!Organism %like any% c('Influent', 'Water control', 'Offshore'))%>%
  group_by(Organism)%>%
  nest()%>%
  mutate(data = map(data, ~t.test(cells_ul ~ DayNight, .x)["p.value"][[1]]))%>%
  unnest(data)%>%
  mutate(test = 'sgr')

microbialGrowth_fdr <- microbialLoad_stats%>%
  bind_rows(sgr_stats)%>%
  mutate(fdr = p.adjust(data, method = "BH"))


# VIZUALIZATIONS -- Venn Diagram ------------------------------------------
venn <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  ungroup()%>%
  filter(activity == 'depletolite')%>%
  select(feature_number, Organism, DayNight)%>%
  unique()

total_deplete <- venn%>%
  select(feature_number)%>%
  unique()%>%
  as.vector()

# Day venn
day_deplete <- venn%>%
  filter(DayNight == 'Day')%>%
  select(-DayNight)

day_list <- day_deplete$feature_number%>%
  as.vector()

day_split <- split(day_list, day_deplete$Organism)

venn.diagram(day_split, 
             'plots/venn_depletolitesDay.png',
             category.names = c('CCA', 'Dic', 'Poc', 'Por', 'T'),
             imagetype = 'png',  
             height = 1500, 
             width = 1800,
             resolution = 300,
             fill = org_colors_no_water)


# Night venn
night_deplete <- venn%>%
  filter(DayNight == 'Night')%>%
  select(-DayNight)

night_list <- night_deplete$feature_number%>%
  as.vector()

night_split <- split(night_list, night_deplete$Organism)

venn.diagram(night_split, 
            'plots/venn_depletolitesNight.png',
            category.names = c('CCA', 'Dic', 'Poc', 'Por', 'T'),
            imagetype = 'png',  
            height = 1500, 
            width = 1800,
            resolution = 300,
            fill = org_colors_no_water)

# Total venn
venn_list <- venn$feature_number%>%
  as.vector()

venn_split <- split(venn_list, venn$Organism)

venn.diagram(venn_split, 
             './plots/venn_depletolitesDayNight.png',
             category.names = c('CCA', 'Dic', 'Poc', 'Por', 'T'),
             imagetype = 'png',  
             height = 1500, 
             width = 1800,
             resolution = 300,
             fill = org_colors_no_water)

# Org Venn
org_venn <- venn%>%
  group_by(Organism)%>%
  nest()%>%
  mutate(list = map(data, ~ .x$feature_number%>%
                      as.vector()),
         dayNightList = map(data, ~ .x$DayNight),
         splits = map2(list, dayNightList, ~split(.x, .y)),
         diagram = map(splits, ~venn.diagram(.x, 
                                             paste0('plots/venn_depletolites', Organism, '.png'),
                                             category.names = c('', ''),
                                             imagetype = 'png',  
                                             height = 1500, 
                                             width = 1800,
                                             resolution = 300,
                                             cex = 3,
                                             fill = c('#FF9300', 'grey'))))


# VIZUALIZATIONS -- Depletolite heatmap -----------------------------------
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


depletion_hc <- feature_stats_wdf%>%
  filter(Timepoint == "T0")%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  ungroup()%>%
  filter(activity == 'depletolite')%>%
  left_join(ecoNet%>%
              select(-network)%>%
              mutate(scan = as.character(scan))%>%
              rename('feature_number' = 'scan'), by = 'feature_number')%>%
  group_by(ecoNetConsensus, DayNight, Replicate, Organism)%>%
  select(-c(xic, ra:activity, ecoNetConsensusScore:matchSource))%>%
  summarize_if(is.numeric, sum, na.rm  = TRUE)%>%
  ungroup()%>%
  unite(sample, c(DayNight, Organism, Replicate), sep = '  ')%>%
  spread(sample, log10)%>%
  # mutate_all(~replace(., is.na(.), 0))%>%
  gather(Organism, T0, 3:ncol(.))%>%
  group_by(ecoNetConsensus)%>%
  mutate(T0 = zscore(T0))

depletion_levels <- depletion_hc$ecoNetConsensus%>% 
  as.factor()

organism_levels <- depletion_hc$Organism%>%
  as.factor()
  # fct_relevel(c('CCA', 'Turf', 'Dictyota'))

pdf('plots/depletoliteHeatmap.pdf', width = 20, height = 15)
depletion_hc%>%
  # filter(DayNight == 'Day')
  ggplot(aes(Organism, ecoNetConsensus, fill = T0)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_y_discrete(limits = rev(levels(depletion_levels))) +
  scale_x_discrete(limits = levels(organism_levels))
dev.off()


# VIZUALIZATIONS -- Nitrogen Compounds ----------------------------------------
nitrogen_compounds <- feature_stats_wdf%>%
  filter(Timepoint == "T0")%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  ungroup()%>%
  # filter(activity == 'depletolite')%>%
  # group_by(feature_number, Organism, DayNight, network, activity)%>%
  # summarize_if(is.numeric, mean)%>%
  # ungroup()%>%
  left_join(ecoNet%>%
              select(-network)%>%
              mutate(scan = as.character(scan))%>%
              rename('feature_number' = 'scan'), by = 'feature_number')%>%
  filter(ecoNetConsensus %like% '%nitrogen%')

write_csv(nitrogen_compounds%>%
            group_by(feature_number, Organism, DayNight, network, activity, ecoNetConsensus)%>%
            summarize_if(is.numeric, mean)%>%
            ungroup()%>%
            left_join(true_hits%>%
                        mutate(feature_number = as.character(feature_number))%>%
                        select(feature_number, Compound_Name),
                      by = 'feature_number'), 'analysis/nitrogen_compounds.csv')

pdf('plots/aminesDifferDayNight.pdf', width = 15, height = 13)
nitrogen_compounds%>%
  select(-c(xic, ra, asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  spread(ecoNetConsensus, log10)%>%
  gather(ecoNetConsensus, log10, 9:12)%>%
  group_by(Replicate, DayNight, ecoNetConsensus, activity)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, ecoNetConsensus, activity)%>%
  mutate(std = sd(log10))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(std = case_when(log10 == 0 ~ NA_real_,
                         TRUE ~ std))%>%
  ggplot(aes(ecoNetConsensus, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10 - 2), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  facet_wrap(~DayNight, nrow = 2) +
  scale_fill_manual(values = c('#EBCC2A', '#006658')) + 
  # scale_fill_manual(values = c('orange', 'grey')) +
  gen_theme() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  labs(x = 'Putative Consensus Annotation', y = bquote(Total ~Production ~Intensity ~(log[10] ~XIC)))
dev.off()
  
pdf('plots/amineProductionOrganism.pdf', width = 15, height = 13)
nitrogen_compounds%>%
  filter(ecoNetConsensus %like% '%Amine%')%>%
  select(-c(xic, ra, asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  spread(Organism, log10)%>%
  gather(Organism, log10, CCA:Turf)%>%
  group_by(Organism, Replicate, DayNight, ecoNetConsensus, activity)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(Organism, DayNight, ecoNetConsensus, activity)%>%
  mutate(std = sd(log10, na.rm = TRUE))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(std = case_when(log10 == 0 ~ NA_real_,
                         TRUE ~ std))%>%
  ggplot(aes(Organism, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10 - 2), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  facet_wrap(~DayNight, nrow = 2) +
  scale_fill_manual(values = c('#EBCC2A', '#006658')) + 
  # coord_flip() + 
  # scale_fill_manual(values = c('orange', 'grey')) +
  gen_theme() +
  labs(x = 'Organism', y = bquote(atop(Total ~Amine ~Production, ~Intensity ~(log[10] ~XIC))))
dev.off()  

# Nitrogen overlap
nitrogenOverlap <- nitrogen_compounds%>%
  filter(ecoNetConsensus %like% '%Amine%')%>%
  select(feature_number, network, DayNight, activity)%>%
  unique()%>%
  filter(DayNight == 'Day' & activity == 'recalcitrant' | DayNight == 'Night' & activity == 'depletolite')

# Amine venn
nitrogen_list <- nitrogenOverlap$feature_number%>%
  as.vector()

nitrogen_split <- split(nitrogen_list, nitrogenOverlap$activity)
  
venn.diagram(nitrogen_split, 
             'plots/venn_amines.png',
             category.names = c('depletolite', 'recalcitrant'),
             imagetype = 'png',  
             height = 1500, 
             width = 1800,
             resolution = 300,
             fill = c('#EBCC2A', '#006658'))

pdf('plots/net2023production.pdf', width = 15, height = 10)
nitrogen_compounds%>%
  filter(ecoNetConsensus %like% '%Amine%',
         network == '2023'| feature_number %like any% c('5247', '7076'))%>%
  select(-c(xic, ra, asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  spread(Organism, log10)%>%
  gather(Organism, log10, CCA:Turf)%>%
  group_by(Organism, Replicate, DayNight, ecoNetConsensus, activity)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(Organism, DayNight, ecoNetConsensus, activity)%>%
  mutate(std = sd(log10, na.rm = TRUE))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(std = case_when(log10 == 0 ~ NA_real_,
                         TRUE ~ std))%>%
  ggplot(aes(Organism, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10 - 2), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  scale_fill_manual(values = c('#EBCC2A', '#006658')) + 
  facet_wrap(~DayNight, nrow = 1) +
  gen_theme() +
  labs(x = 'Organism', y = bquote(atop(Total ~Network ~2023 ~Production, ~Intensity ~(log[10] ~XIC))))

dev.off()

net2023Depletion <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  ungroup()%>%
  left_join(ecoNet%>%
              select(-network)%>%
              mutate(scan = as.character(scan))%>%
              rename('feature_number' = 'scan'), by = 'feature_number')%>%
  filter(network == '2023'| feature_number %like any% c('5247', '7076'))%>%
  group_by(Replicate, DayNight, activity, Timepoint)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, Timepoint, activity)%>%
  mutate(std = sd(log10, na.rm = TRUE))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()


net2023Depletion%>%
  ggplot(aes(Timepoint, xic, fill = activity)) +
  geom_errorbar(aes(ymax = xic + std, ymin = xic -3)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('#EBCC2A', '#006658')) +
  scale_y_log10()+
  facet_wrap(~activity) +
  gen_theme()


# VIZUALIZATIONS -- Lipid-like compounds ----------------------------------
lipid_compounds <- feature_stats_wdf%>%
  filter(Timepoint == "T0")%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  ungroup()%>%
  # filter(activity == 'depletolite')%>%
  # group_by(feature_number, Organism, DayNight, network, activity)%>%
  # summarize_if(is.numeric, mean)%>%
  # ungroup()%>%
  left_join(ecoNet%>%
              select(-network)%>%
              mutate(scan = as.character(scan))%>%
              rename('feature_number' = 'scan'), by = 'feature_number')%>%
  filter(ecoNetConsensus %like% '%ipid%')

pdf('plots/lipidLikeCompounds.pdf', width = 15, height = 10)
lipid_compounds%>%
  select(-c(xic, ra, asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  spread(ecoNetConsensus, log10)%>%
  gather(ecoNetConsensus, log10, 9:ncol(.))%>%
  group_by(Replicate, DayNight, ecoNetConsensus, activity)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, ecoNetConsensus, activity)%>%
  mutate(std = sd(log10))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(std = case_when(log10 == 0 ~ NA_real_,
                         TRUE ~ std))%>%
  ggplot(aes(ecoNetConsensus, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10 - 2), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  facet_wrap(~DayNight, nrow = 2) +
  scale_fill_manual(values = c('#78B7C5','#EBCC2A', "#006658")) +
  gen_theme() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  labs(x = 'Putative Consensus Annotation', y = bquote(Total ~Production ~Intensity ~(log[10] ~XIC)))
dev.off()

pdf('plots/fattyAcidEsterProductionOrganism.pdf', width = 15, height = 13)
lipid_compounds%>%
  filter(ecoNetConsensus %like% '%Fatty acid esters%')%>%
  select(-c(xic, ra, asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  spread(Organism, log10)%>%
  gather(Organism, log10, CCA:Turf)%>%
  group_by(Organism, Replicate, DayNight, ecoNetConsensus, activity)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(Organism, DayNight, ecoNetConsensus, activity)%>%
  mutate(std = sd(log10, na.rm = TRUE))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(std = case_when(log10 == 0 ~ NA_real_,
                         TRUE ~ std))%>%
  ggplot(aes(Organism, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10 - 2), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  facet_wrap(~DayNight, nrow = 2) +
  scale_fill_manual(values = c('#EBCC2A', '#006658')) + 
  # coord_flip() + 
  # scale_fill_manual(values = c('orange', 'grey')) +
  gen_theme() +
  labs(x = 'Organism', y = bquote(atop(Total ~fatty ~acid ~ester, ~production ~intensity ~(log[10] ~XIC))))
dev.off()

pdf('plots/carnitineProductionOrganism.pdf', width = 15, height = 13)
lipid_compounds%>%
  filter(network %like any% c('21', '1583', '1841'))%>%
  group_by(Organism, DayNight, activity, Replicate)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(Organism, DayNight, activity)%>%
  mutate(std = sd(log10))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  ggplot(aes(Organism, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10-3), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  facet_wrap(~DayNight, nrow = 2) +
  scale_fill_manual(values = c('#EBCC2A', '#006658')) +
  gen_theme() +
  labs(x = 'Organism', y = bquote(atop(Total ~carnitine ~production, ~intensity ~(log[10] ~XIC))))

lipid_compounds%>%
  filter(network %like any% c('21', '1583', '1841'))%>%
  group_by(Organism, DayNight, activity)%>%
  mutate(std = sd(log10))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  ggplot(aes(Organism, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10-3), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  facet_wrap(~DayNight, nrow = 2) +
  scale_fill_manual(values = c('#EBCC2A', '#006658')) +
  gen_theme() +
  labs(x = 'Organism', y = bquote(atop(Average ~carnitine ~production, ~intensity ~(log[10] ~XIC))))
dev.off()

fattyEsterTable <- lipid_compounds%>%
  filter(ecoNetConsensus %like% '%Fatty acid esters%')%>%
  select(network, net_act, activity)%>%
  unique()%>%
  left_join(networking%>%
              filter(binary_ID == 1)%>%
              mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                                         TRUE ~ network))%>%
              select(net_act, combined_ID)%>%
              unique(), by = 'net_act')%>%
  select(-net_act)

write_csv(fattyEsterTable, 'analysis/fattyEsterTable.csv')

# VIZUALIZATIONS -- Organic Acids -----------------------------------------
acid_compounds <- feature_stats_wdf%>%
  filter(Timepoint == "T0")%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  ungroup()%>%
  # filter(activity == 'depletolite')%>%
  # group_by(feature_number, Organism, DayNight, network, activity)%>%
  # summarize_if(is.numeric, mean)%>%
  # ungroup()%>%
  left_join(ecoNet%>%
              select(-network)%>%
              mutate(scan = as.character(scan))%>%
              rename('feature_number' = 'scan'), by = 'feature_number')%>%
  filter(ecoNetConsensus %like% 'Organic acid%')

pdf('plots/organicAcidCompounds.pdf', width = 15, height = 10)
acid_compounds%>%
  select(-c(xic, ra, asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  spread(ecoNetConsensus, log10)%>%
  gather(ecoNetConsensus, log10, 9:ncol(.))%>%
  group_by(Replicate, DayNight, ecoNetConsensus, activity)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, ecoNetConsensus, activity)%>%
  mutate(std = sd(log10))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(std = case_when(log10 == 0 ~ NA_real_,
                         TRUE ~ std))%>%
  ggplot(aes(ecoNetConsensus, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10 - 2), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  facet_wrap(~DayNight, nrow = 2) +
  scale_fill_manual(values = c('#78B7C5','#EBCC2A', "#006658")) +
  # scale_fill_manual(values = c('orange', 'grey')) +
  gen_theme() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  labs(x = 'Putative Consensus Annotation', y = bquote(Total ~Production ~Intensity ~(log[10] ~XIC)))
dev.off()


# VIZUALIZATIONS -- Organoheterocyclic -------------------------------------
organoheterocyclic_compounds <- feature_stats_wdf%>%
  filter(Timepoint == "T0")%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  ungroup()%>%
  # filter(activity == 'depletolite')%>%
  # group_by(feature_number, Organism, DayNight, network, activity)%>%
  # summarize_if(is.numeric, mean)%>%
  # ungroup()%>%
  left_join(ecoNet%>%
              select(-network)%>%
              mutate(scan = as.character(scan))%>%
              rename('feature_number' = 'scan'), by = 'feature_number')%>%
  filter(ecoNetConsensus %like% 'Organoheterocyclic%')

pdf('plots/organoheterocyclicCompounds.pdf', width = 15, height = 10)
organoheterocyclic_compounds%>%
  select(-c(xic, ra, asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  spread(ecoNetConsensus, log10)%>%
  gather(ecoNetConsensus, log10, 9:ncol(.))%>%
  group_by(Replicate, DayNight, ecoNetConsensus, activity)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, ecoNetConsensus, activity)%>%
  mutate(std = sd(log10))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(std = case_when(log10 == 0 ~ NA_real_,
                         TRUE ~ std))%>%
  ggplot(aes(ecoNetConsensus, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10 - 2), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  facet_wrap(~DayNight, nrow = 2) +
  scale_fill_manual(values = c('#EBCC2A', "#006658")) +
  # scale_fill_manual(values = c('orange', 'grey')) +
  gen_theme() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  labs(x = 'Putative Consensus Annotation', y = bquote(Total ~Production ~Intensity ~(log[10] ~XIC)))
dev.off()

pdf('plots/nAlkylindoleProductionOrganism.pdf', width = 15, height = 13)
organoheterocyclic_compounds%>%
  filter(ecoNetConsensus %like% '%N-alkyl%')%>%
  select(-c(xic, ra, asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  spread(Organism, log10)%>%
  gather(Organism, log10, CCA:Turf)%>%
  group_by(Organism, Replicate, DayNight, ecoNetConsensus, activity)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(Organism, DayNight, ecoNetConsensus, activity)%>%
  mutate(std = sd(log10, na.rm = TRUE))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(std = case_when(log10 == 0 ~ NA_real_,
                         TRUE ~ std))%>%
  ggplot(aes(Organism, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10 - 2), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  facet_wrap(~DayNight, nrow = 2) +
  scale_fill_manual(values = c('#EBCC2A', '#006658')) + 
  # coord_flip() + 
  # scale_fill_manual(values = c('orange', 'grey')) +
  gen_theme() +
  labs(x = 'Organism', y = bquote(atop(Total ~N-alkylindole ~Production, ~Intensity ~(log[10] ~XIC))))
dev.off()  

pdf('plots/triazineProductionOrganism.pdf', width = 15, height = 13)
organoheterocyclic_compounds%>%
  filter(ecoNetConsensus %like% '%triazines')%>%
  select(-c(xic, ra, asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  spread(Organism, log10)%>%
  gather(Organism, log10, CCA:Turf)%>%
  group_by(Organism, Replicate, DayNight, ecoNetConsensus, activity)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(Organism, DayNight, ecoNetConsensus, activity)%>%
  mutate(std = sd(log10, na.rm = TRUE))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(std = case_when(log10 == 0 ~ NA_real_,
                         TRUE ~ std))%>%
  ggplot(aes(Organism, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10 - 2), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  facet_wrap(~DayNight, nrow = 2) +
  scale_fill_manual(values = c('#EBCC2A', '#006658')) + 
  # coord_flip() + 
  # scale_fill_manual(values = c('orange', 'grey')) +
  gen_theme() +
  labs(x = 'Organism', y = bquote(atop(Total ~Triazine ~Production, ~Intensity ~(log[10] ~XIC))))
dev.off()  

# VIZUALIZATIONS -- Benzenoids compounds ----------------------------------
benzenoid_compounds <- feature_stats_wdf%>%
  filter(Timepoint == "T0")%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  ungroup()%>%
  # filter(activity == 'depletolite')%>%
  # group_by(feature_number, Organism, DayNight, network, activity)%>%
  # summarize_if(is.numeric, mean)%>%
  # ungroup()%>%
  left_join(ecoNet%>%
              select(-network)%>%
              mutate(scan = as.character(scan))%>%
              rename('feature_number' = 'scan'), by = 'feature_number')%>%
  filter(ecoNetConsensus %like% 'Benzenoids%')

pdf('plots/benzenoidsCompounds.pdf', width = 15, height = 10)
benzenoid_compounds%>%
  select(-c(xic, ra, asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  spread(ecoNetConsensus, log10)%>%
  gather(ecoNetConsensus, log10, 9:ncol(.))%>%
  group_by(Replicate, DayNight, ecoNetConsensus, activity)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, ecoNetConsensus, activity)%>%
  mutate(std = sd(log10))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(std = case_when(log10 == 0 ~ NA_real_,
                         TRUE ~ std))%>%
  ggplot(aes(ecoNetConsensus, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10 - 2), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  facet_wrap(~DayNight, nrow = 2) +
  scale_fill_manual(values = c('#78B7C5','#EBCC2A', "#006658")) +
  gen_theme() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  labs(x = 'Putative Consensus Annotation', y = bquote(Total ~Production ~Intensity ~(log[10] ~XIC)))
dev.off()

# VIZUALIZATIONS -- Organic Oxygen compounds ----------------------------------
oxygen_compounds <- feature_stats_wdf%>%
  filter(Timepoint == "T0")%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  ungroup()%>%
  # filter(activity == 'depletolite')%>%
  # group_by(feature_number, Organism, DayNight, network, activity)%>%
  # summarize_if(is.numeric, mean)%>%
  # ungroup()%>%
  left_join(ecoNet%>%
              select(-network)%>%
              mutate(scan = as.character(scan))%>%
              rename('feature_number' = 'scan'), by = 'feature_number')%>%
  filter(ecoNetConsensus %like% '%oxygen%')

pdf('plots/oxygenCompounds.pdf', width = 15, height = 10)
oxygen_compounds%>%
  select(-c(xic, ra, asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  spread(ecoNetConsensus, log10)%>%
  gather(ecoNetConsensus, log10, 9:ncol(.))%>%
  group_by(Replicate, DayNight, ecoNetConsensus, activity)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, ecoNetConsensus, activity)%>%
  mutate(std = sd(log10))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(std = case_when(log10 == 0 ~ NA_real_,
                         TRUE ~ std))%>%
  ggplot(aes(ecoNetConsensus, log10, fill = activity)) +
  geom_errorbar(aes(ymax = log10 + std, ymin = log10 - 2), position = 'dodge2') +
  geom_bar(stat = 'identity', position = 'dodge2') +
  facet_wrap(~DayNight, nrow = 2) +
  scale_fill_manual(values = c('#78B7C5','#EBCC2A', "#006658")) +
  gen_theme() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  labs(x = 'Putative Consensus Annotation', y = bquote(Total ~Production ~Intensity ~(log[10] ~XIC)))
dev.off()

# VIZUALIZATIONS -- Depletion heatmap -------------------------------------------------------
production_hc <- feature_stats_wdf%>%
  filter(Timepoint == 'T0')%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  # separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  ungroup()%>%
  filter(activity == 'depletolite')%>%
  # left_join(log2_change_vals%>%
  #             ungroup(), by = c('feature_number', 'Replicate', 'Organism', 'DayNight'))%>%
  select(-c(network, numberOfNodes, ecoNetConsensusScore, xic, ra, asin, net_act))%>%
  group_by(ecoNetConsensus, DayNight, Replicate, Organism)%>%
  summarize_if(is.numeric, sum, na.rm  = TRUE)%>%
  ungroup()%>%
  group_by(ecoNetConsensus, DayNight, Organism)%>%
  summarize_if(is.numeric, mean, na.rm  = TRUE)%>%
  ungroup()%>%
  unite(sample, c('DayNight', 'Organism'), sep = '  ')%>%
  spread(sample, log10)%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  gather(sample, log10, 2:ncol(.))

depletion_levels <- production_hc$ecoNetConsensus%>% 
  as.factor()

organism_levels <- production_hc$sample%>%
  as.factor()
  # fct_relevel(c('Day  CCA', 'Day  Turf', 'Day  Dictyota', 'Day  Pocillopora verrucosa', 'Day  Porites lobata', 'Night  CCA', 'Night  Turf', 'Night  Dictyota'))

pdf('plots/depletoliteHeatmapDayNight.pdf', width = 20, height = 15)
production_hc%>%
  group_by(ecoNetConsensus)%>%
  mutate(log10 = zscore(log10))%>%
  ggplot(aes(sample, ecoNetConsensus, fill = log10)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_y_discrete(limits = rev(levels(depletion_levels))) +
  scale_x_discrete(limits = levels(organism_levels)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()


# VIZUALIZATIONS -- PCOA of Consensus annotations -------------------------
# Have to load in EmojiFont here otherwise it messes with unicode points
library(emojifont)
load.emojifont('OpenSansEmoji.ttf')

pcoa <- feature_stats_wdf%>%
  filter(Timepoint == 'T0')%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  # separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  ungroup()%>%
  filter(activity == 'depletolite')%>%
  # left_join(log2_change_vals%>%
  #             ungroup(), by = c('feature_number', 'Replicate', 'Organism', 'DayNight'))%>%
  select(-c(network, numberOfNodes, ecoNetConsensusScore, xic, ra, asin, net_act))%>%
  group_by(ecoNetConsensus, DayNight, Replicate, Organism)%>%
  summarize_if(is.numeric, sum, na.rm  = TRUE)%>%
  ungroup()%>%
  group_by(ecoNetConsensus, DayNight, Organism, Replicate)%>%
  summarize_if(is.numeric, mean, na.rm  = TRUE)%>%
  ungroup()%>%
  unite(sample, c('DayNight', 'Organism', 'Replicate'), sep = '  ')%>%
  spread(sample, log10)%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  gather(sample, log10, 2:ncol(.))%>%
  spread(ecoNetConsensus, log10)%>%
  column_to_rownames(var = "sample")%>%
  vegdist(na.rm = TRUE)%>%
  pcoa()

## Plot Eigenvalues
pcoa$values[1:10,]%>%
  as.data.frame()%>%
  rownames_to_column("Axis")%>%
  mutate(axis = as.numeric(Axis))%>%
  ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, color = "red", vjust = -0.5)

pdf("plots/annotationPcoa.pdf", width = 15, height = 13)
pcoa$vectors%>%
  as.data.frame()%>%
  rownames_to_column(var = "sample")%>%
  separate(sample, c('DayNight', 'Organism', 'Replicate'), sep = '  ')%>%
  mutate(shape = case_when(DayNight == 'Day' ~ '\u2600',
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
  xlab(str_c("Axis 1", " (", round(pcoa$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
  ylab(str_c("Axis 2", " (", round(pcoa$values$Relative_eig[2], digits = 4)*100, "%)", sep = ""))
dev.off()

# Pairwise PERMANOVA
permanova_df <- feature_stats_wdf%>%
  filter(Timepoint == 'T0')%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  ungroup()%>%
  filter(activity == 'depletolite')%>%
  select(-c(network, numberOfNodes, ecoNetConsensusScore, xic, ra, asin, net_act))%>%
  group_by(ecoNetConsensus, DayNight, Replicate, Organism)%>%
  summarize_if(is.numeric, sum, na.rm  = TRUE)%>%
  ungroup()%>%
  spread(ecoNetConsensus, log10)%>%
  mutate_all(~replace(., is.na(.), 0))

permanova_organism <- permanova_df%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(data = map(data, ~adonis(.x[3:ncol(.x)] ~ Organism, .x, perm = 1000, method = "bray", na.rm = TRUE)))

permanova_DayNight <- permanova_df%>%
  adonis(.[4:ncol(.)] ~ DayNight, ., perm = 1000, method = "bray", na.rm = TRUE)

# VIZUALIZATIONS -- Cytoscape file ----------------------------------------
night_activity <- all_activity%>%
  filter(DayNight == 'Night')

day_activity <- all_activity%>%
  filter(DayNight == 'Day')

cyto_export <- feature_stats_wdf%>%
  filter(Timepoint == "T0")%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  ungroup()%>%
  left_join(ecoNet%>%
              select(-network)%>%
              mutate(scan = as.character(scan))%>%
              rename('feature_number' = 'scan'), by = 'feature_number')%>%
  select(feature_number, Organism, DayNight, net_act, ecoNetConsensus, log10)%>%
  group_by(feature_number, Organism, DayNight, ecoNetConsensus, net_act)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  unite(sample, c('DayNight', 'Organism'), sep = '_')%>%
  spread(sample, log10)%>%
  mutate(Daysum = rowSums(across(Day_CCA:Day_Turf), na.rm = TRUE),
         Nightsum = rowSums(across(Night_CCA:Night_Turf), na.rm = TRUE))%>%
  left_join(night_activity, by = 'net_act')%>%
  left_join(day_activity, by = 'net_act', suffix = c('_night', '_day'))



write_csv(cyto_export, 'analysis/cytoDiel.csv')
