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
library(pscl)
library(lme4)
library(processx)

#PCoA, PERMANOVA
library(vegan)
library(ape)

#Visualizations
library(wesanderson)
library(RColorBrewer)
library(ggnewscale)
library(ggpattern)
library(scales)

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

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
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
canopusSirius4 <- read_csv('~/Documents/SDSU_Scripps/ConCISE/investigations/canopus_summary.csv')

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

# canopusRaw <- read_csv('~/Documents/GitHub/conCISE/verification/dataset1_unmodified/canopus_summary.csv')
# # 
# net199 <- networking%>%
#   select(feature_number, network, ZodiacMF:NOSC)%>%
#   filter(network == '199')%>%
#   left_join(feature_table_raw%>%
#               select(feature_number, `row m/z`)%>%
#               mutate(feature_number = as.character(feature_number)), by = 'feature_number')%>%
#   left_join(ecoNet%>%
#               mutate(scan = as.character(scan)), by = c('feature_number' = 'scan', 'network' = 'network'))%>%
#   left_join(canopusRaw%>%
#               select(-c(molecularFormula:adduct))%>%
#               mutate(scan = as.character(scan)), by = c('feature_number' = 'scan'))
# 

# net199Features <- net199$feature_number
# write_lines(net199Features, '~/Downloads/Net199FeatureList.txt')

# write_csv(net199, '~/Downloads/net199DeeperLook.csv')

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
  mutate(hour = case_when(Organism == 'Turf' & DayNight == 'Day' ~ 23,
                          Organism == 'Dictyota'& DayNight == 'Day'  ~ 23,
                          Organism == 'CCA' & DayNight == 'Day' ~ 23,
                          Organism == 'Porites lobata'& DayNight == 'Day' ~ 48,
                          Organism == 'Pocillopora verrucosa'& DayNight == 'Day' ~ 37,
                          Organism == 'Water control'& DayNight == 'Day' ~ 37,
                          TRUE ~ 24),
         final_cells = case_when(Organism == 'Turf'& DayNight == 'Day' ~ T4,
                                 Organism == 'Dictyota'  & DayNight == 'Day'~ T4,
                                 Organism == 'CCA' & DayNight == 'Day' ~ T4,
                                 Organism == 'Porites lobata'& DayNight == 'Day' ~ TF,
                                 Organism == 'Pocillopora verrucosa' & DayNight == 'Day'~ T6,
                                 Organism == 'Water control' & DayNight == 'Day'~ T6,
                                 TRUE ~ T5))%>%
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

# fwrite(list(dorcierr_real_features), file = '~/Documents/SDSU_Scripps/DORCIERR/Datasets/SIRIUS4/shortFeatures.txt')

# FILTERING -- LOG2 change vals --------------------------------------------
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
  mutate(log2_change = case_when(TF > 0 & T0 == 0 ~ max(log2_change),
                                 TF == 0 & T0 > 0 ~ min(log2_change),
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

min_filtered <- min_filter_pre%>%
  left_join(networking%>%
              select(network, feature_number), by = 'feature_number')%>%
  inner_join(min_filter%>%
               select(feature_number, DayNight, Organism), by = c('DayNight', 'Organism', 'feature_number'))

org_filter <- min_filter_pre%>%
  left_join(networking%>%
              select(network, feature_number), by = 'feature_number')%>%
  inner_join(org_exometabolites, by = c('DayNight', 'Organism', 'feature_number'))

all_filter <- min_filter_pre%>%
  left_join(networking%>%
              select(network, feature_number), by = 'feature_number')%>%
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
  select(feature_number, DayNight)%>%
  unique()%>%
  group_by(DayNight)%>%
  mutate(num_features = 1)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()

num_features_org_filter <- org_filter%>%
  ungroup()%>%
  select(feature_number, DayNight)%>%
  unique()%>%
  group_by(DayNight)%>%
  mutate(num_features = 1)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()

num_features_all_filter <- all_filter%>%
  ungroup()%>%
  select(feature_number, DayNight)%>%
  unique()%>%
  group_by(DayNight)%>%
  mutate(num_features = 1)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()

num_features_incremental <- no_min_filter%>%
  ungroup()%>%
  nest(data = everything())%>%
  mutate(exometabolite = map(data, ~inner_join(.x, org_exometabolites, by = c('DayNight', 'Organism', 'feature_number'))%>%
                               # mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                               #                            TRUE ~ network))%>%
                               # filter(network != -1)%>%
                               select(feature_number, DayNight)%>%
                               unique()%>%
                               group_by(DayNight)%>%
                               mutate(count = 1)%>%
                               summarize_if(is.numeric, sum)),
         min = map(data, ~ inner_join(.x, min_filter%>%
                                  select(feature_number, DayNight, Organism), by = c('DayNight', 'Organism', 'feature_number'))%>%
                     # mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                     #                            TRUE ~ network))%>%
                     select(feature_number, DayNight)%>%
                     unique()%>%
                     group_by(DayNight)%>%
                     mutate(count = 1)%>%
                     summarize_if(is.numeric, sum)),
         all = map(data, ~ inner_join(.x, min_filter%>%
                                        select(feature_number, DayNight, Organism), by = c('DayNight', 'Organism', 'feature_number'))%>%
                     
                     inner_join(org_exometabolites, by = c('DayNight', 'Organism', 'feature_number'))%>%
                     mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                                                TRUE ~ network))%>%
                     # filter(net_act < 0)%>%
                     select(net_act)%>%
                     unique()%>%
                     # group_by(DayNight)%>%
                     mutate(count = 1)%>%
                     summarize_if(is.numeric, sum))
  )%>%
  select(-data)


numberClassifiedExometabolites <- no_min_filter%>%
  inner_join(min_filter%>%
               select(feature_number, DayNight, Organism), 
             by = c('DayNight', 'Organism', 'feature_number'))%>%
  inner_join(org_exometabolites%>%
               select(feature_number, DayNight, Organism),
             by = c('DayNight', 'Organism', 'feature_number'))%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                                                        TRUE ~ network))%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  filter(!is.na(ecoNetConsensus))%>%
  select(net_act)%>%
  unique()%>%
  nrow()
  


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
  # inner_join(log2_features%>%
  #              select(feature_number, DayNight), by = c('DayNight', 'feature_number'))%>%
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


# STATS  -- LMER -- subnetwork level ------------------------------------------
decreaseNetworks <- dom_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = "feature_number")%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  select(-c(xic, ra, asin))%>%
  spread(Timepoint, log10)%>%
  group_by(net_act, DayNight)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  filter(TF < T0)

increaseNetworks <- dom_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = "feature_number")%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  select(-c(xic, ra, asin))%>%
  spread(Timepoint, log10)%>%
  group_by(net_act, DayNight)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  filter(TF > T0)

decreaseNetworkVector <- decreaseNetworks$net_act
increaseNetworkVector <- increaseNetworks$net_act


lmer_test <- dom_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = "feature_number")%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  group_by(net_act, DayNight)%>%
  nest()%>%
  mutate(mixedModel = map(data, ~ mutate(.x, Organism = as.factor(Organism),
                                                 feature_number = as.factor(feature_number))%>%
                                     lmer(log10 ~ Timepoint + (1|Organism) + (1|feature_number) , data =., control =  lmerControl(check.nlev.gtr.1 = "ignore",
                                                                                                                                  check.conv.singular = 'ignore',
                                                                                                                                  check.nobs.vs.nRE = 'ignore'))%>%
                                     car::Anova()%>%
                                     .[['Pr(>Chisq)']]))%>%
  select(-data)%>%
  ungroup()

all_activity <- lmer_test%>%
  mutate(activity = p.adjust(mixedModel, method = 'BH'))%>%
  mutate(activity = case_when(mixedModel < 0.05 & net_act %in% decreaseNetworkVector ~ "labile",
                                mixedModel < 0.05 & net_act %in% increaseNetworkVector ~ 'accumolite',
                                           TRUE ~ "recalcitrant"))%>%
  select(net_act, DayNight, activity)



# Analaysis -- Counting number of abundant labile subetnworks -------------
all_activity%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number),
                     net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                                       TRUE ~ network))%>%
              select(net_act, ecoNetConsensus)%>%
              unique(),
            by = c('net_act'))%>%
  # filter(!is.na(ecoNetConsensus))%>%
  filter(activity == 'labile',
         DayNight == 'Day')%>%
  # net_act < 0)%>%
  select(-DayNight)%>%
  unique()%>%
  nrow()


feature_stats_wdf%>%
  ungroup()%>%
  select(feature_number, Organism, DayNight)%>%
  left_join(networking%>%
              select(feature_number, network), by = "feature_number")%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('DayNight', 'net_act'))%>%
  filter(activity == 'labile')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number))%>%
              select(ecoNetConsensus, feature_number, network),
                     by = c('feature_number', 'network'))%>%
  select(network, Organism, net_act, ecoNetConsensus, DayNight)%>%
  group_by(DayNight, Organism)%>%
  unique()%>%
  mutate(n = 1)%>%
  summarize_if(is.numeric, sum)
  # filter(!is.na(ecoNetConsensus),
  #        net_act < 0)%>%
  nrow()



# PRE-STATS -- affinity ----------------------------------------------------------------
fcm_affinity <- fcm_wdf%>%
  select(Organism, DayNight, Replicate, Timepoint, `Cells µL-1`)%>%
  spread(Timepoint, `Cells µL-1`)%>%
  group_by(Organism, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         biomass48 = TF,
         biomass0 = T0)%>%
  select(-c(T1:T6))%>%
  rename(cellsT0 = T0,
         cellsTF = TF)
  

affinity <- dom_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = "feature_number")%>%
  # filter(network != '-1')%>%
  select(-c(log10, ra, asin))%>%
  spread(Timepoint, xic)%>%
  group_by(feature_number, Organism, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE))%>%
  left_join(fcm_affinity, by = c('Organism', 'DayNight', 'Replicate'))%>%
  ungroup()%>%
  mutate(TF = case_when(TF == 0 ~ min(TF[.$TF != 0], na.rm = TRUE),
                        TRUE ~ as.numeric(TF)),
         T0 = case_when(T0 == 0 ~ min(T0[.$T0 != 0], na.rm = TRUE),
                        TRUE ~ as.numeric(T0)),
         numerator = log(TF/T0),
         r = log(cellsTF/cellsT0)/48,
         denominator = (cellsT0*exp(((r*48)/(r*3600)))),
         affinity = numerator/(-cellsTF*48),
         betterAffinity = -numerator/denominator)

affinity%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')%>%
  # filter(network == '437')%>%
  ggplot(aes(log10(T0), betterAffinity, color = activity)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = 'lm')
  # scale_color_manual(values = c("#FF850A", "lightBlue"))


# STATS -- Microbial load -------------------------------------------------
loadAnova <- fcm_T0_5%>%
  filter(Organism != 'Offshore',
         Organism != 'Influent')%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(anova = map(data, ~aov(final_cells~Organism, data = .x)%>%
                      tidy()),
         tukey = map(data, ~ aov(final_cells~ Organism, data = .x)%>%
                       TukeyHSD()%>%
                       tidy()))
  

# STATS -- glm stoichiometry ----------------------------------------------
mul_reg <- affinity%>%
  inner_join(dom_stats_wdf%>%
               ungroup()%>%
               select(feature_number, Organism, DayNight)%>%
               unique(), 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network, N, P, C, O, H, NOSC),
            by = c('feature_number', 'network'))%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  mutate(n_presence = case_when(N > 0 ~ 1,
                                TRUE ~ 1))%>%
  mutate(nc = N/C,
         oc = O/C,
         hc = H/C,
         log_nc = log10(nc +1),
         log_oc = log10(oc + 1),
         log_hc = log10(hc + 1),
         logxic = log10(T0),
         log_snc = log10(nc*T0),
         log_soc = log10(oc*T0),
         log_shc = log10(hc*T0),
         Replicate = as.character(Replicate))%>%
  filter(NOSC < 0,
         nc < 1)%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')%>%
  mutate(activity = as.factor(activity),
         network = as.factor(network))

craig <- mul_reg%>%
  left_join(ecoNet%>% 
              rename(feature_number = 1)%>%
              mutate(feature_number = as.character(feature_number),
                     network = as.factor(network)), 
            by = c('feature_number', 'network'))

lm_test_mulreg <- lm(betterAffinity ~ `row m/z`+ T0 + NOSC + log_oc + log_nc+ log_hc, data = mul_reg)

lm_grouped <- mul_reg%>%
  left_join(ecoNet%>%
              rename(feature_number = 1)%>%
              mutate(feature_number = as.character(feature_number),
                     network = as.factor(network)), by = c('feature_number', 'network'))%>%
  separate(ecoNetConsensus, c('superclass', 'class', 'subclass'), sep = ';')%>%
  group_by(superclass)%>%
  nest()%>%
  mutate(n = map(data, ~ nrow(.x)),
         data = map(data, ~ lm(betterAffinity ~ `row m/z`+ T0 + NOSC + log_nc + log_oc + log_hc, data = .x)),
         pVal = map(data, ~ lmp(.x)),
         rSquared = map(data, ~ summary(.x)['adj.r.squared']))%>%
  select(-data)%>%
  unnest(c(rSquared, pVal))%>%
  as.data.frame()

write_csv(lm_grouped, 'analysis/multipleRegressionSuperclassPvaluesRsq.csv')

## Network averaged stoichiometry
networkMeanMreg <- mul_reg%>%
  group_by(network, net_act, activity)%>%
  summarise_if(is.numeric, mean)%>%
  ungroup()


lm_activity <- glm(activity ~ `row m/z`+ T0 + NOSC + log_oc + log_nc + log_hc, family = 'binomial'(link ='logit'), data = networkMeanMreg)
nmod <- glm(activity~1, family = 'binomial'(link='logit'), data = networkMeanMreg) ##"null" mod
anova(nmod, lm_activity, test = 'Chisq')
pR2(lm_activity)
paste('Whole logistic Regression n =', nrow(networkMeanMreg))


lm_activity_grouped <- mul_reg%>%
  left_join(ecoNet%>%
              rename(feature_number = 1)%>%
              mutate(feature_number = as.character(feature_number),
                     network = as.factor(network)), by = c('feature_number', 'network'))%>%
  separate(ecoNetConsensus, c('superclass', 'class', 'subclass'), sep = ';')%>%
  group_by(net_act, activity, superclass)%>%
  summarise_if(is.numeric, mean)%>%
  ungroup()%>%
  group_by(superclass)%>%
  nest()%>%
  mutate(n = map(data, ~ nrow(.)))%>%
  filter(n > 2)%>%
  mutate(data = map(data, ~ glm(activity ~ `row m/z`+ T0 + NOSC + log_nc + log_oc + log_hc, family = 'binomial'(link ='logit'), data = .x)),
         # chi = map(data, ~ anova(.x, test="Chisq")),
         rSquared = map(data, ~pR2(.x)[['McFadden']]))%>%
  # select(-rSquared)%>%
  # select(-data)
  unnest(rSquared)



# STATS -- logit model network to log2 ------------------------------------
networkedFeatures <- mul_reg%>%
  select(network, feature_number)%>%
  unique()%>%
  filter(network != '-1')%>%
  group_by(network)%>%
  mutate(numberOfNodes = 1,
         numberOfNodes = sum(numberOfNodes))%>%
  ungroup()%>%
  select(network, numberOfNodes)

nodeLimits <- c(5,10,15,20,25,30)

networkAovDf <- mul_reg%>%
  left_join(networkedFeatures, by = 'network')%>%
  mutate(network = as.factor(network))%>%
  select(network, numberOfNodes, affinity, betterAffinity)

networkAov <- networkAovDf%>%
  mutate(minimumNode = 1)%>%
  bind_rows(networkAovDf%>%
              mutate(minimumNode = 2),
            networkAovDf%>%
              mutate(minimumNode = 4),
            networkAovDf%>%
              mutate(minimumNode = 9),
            networkAovDf%>%
              mutate(minimumNode = 14))%>%
  mutate(minimumNodes = minimumNode)%>%
  unique()%>%
  group_by(minimumNode)%>%
  nest()%>%
  mutate(data = map(data, ~ filter(.x, numberOfNodes > minimumNodes)),
         nNetworks = map(data, ~ unique(.x$network)%>%
                           length()),
         nRows = map(data, ~ nrow(.x)),
         anova = map(data, ~ aov(affinity ~ network, data = .x)%>%
                       tidy()),
         pVal = map(anova, ~filter(.x, !term == "Residuals")%>%
                      select(term, p.value)),
         rSquared = map(anova, ~ .x$sumsq[1]/ (.x$sumsq[1] + .x$sumsq[2])))%>%
  select(-c(data, anova))%>%
  unnest(c(nNetworks:rSquared))%>%
  as.data.frame()

write_csv(networkAov, 'analysis/anovaStructure.csv')

# Supplemental -- Checking for correlation and gaussian distribution in NOSC, Mass, NC and PC --------
# LilliFors and other tests for normallity will not be effective due to large sample size. 
# QQ-plots used to establish how 'normal' the data appear
glm_df <- affinity%>%
  left_join(networking%>%
              select(feature_number, NOSC), by = "feature_number")%>%
  filter(network != "-1",
         NOSC <= 0)%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')%>%
  mutate(activity2 = activity)%>% 
  select(-c(T0,TF))%>%
  gather(response_var, value, NOSC:`row m/z`)

log_reg <- glm_df%>%
  inner_join(min_filter%>%
               select(feature_number, Organism, DayNight)%>%
               unique(),
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  inner_join(org_exometabolites, 
             by = c('feature_number', 'Organism', 'DayNight'))

nc_logreg <- affinity%>%
  inner_join(feature_stats_wdf%>%
               ungroup()%>%
               select(feature_number, Organism, DayNight)%>%
               unique(), 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, C:dG),
              by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')%>%
  group_by(Organism, DayNight, Replicate)%>%
  mutate(nc = N/C,
         pc = P/C,
         sample_nc = T0*nc,
         sample_pc = T0*pc,
         activity2 = activity)%>%
  ungroup()%>%
  mutate(log_snc = log(sample_nc),
         log_spc = log(sample_pc))

affinityCheck <- (log_reg%>%
                 filter(response_var == 'NOSC'))$betterAffinity

#N:C
nc_check <- networkMeanMreg$nc

lnc_check <- networkMeanMreg$log_nc
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
networkMeanMreg%>%
  ggplot(aes(nc)) +
  geom_histogram(bins = 100) +
  gen_theme()

networkMeanMreg%>%
  ggplot(aes(log_nc)) +
  geom_histogram(bins = 100) +
  gen_theme()

networkMeanMreg%>%
  ggplot(aes(log_oc)) +
  geom_histogram(bins = 100) +
  gen_theme()

car::qqPlot(affinityCheck,
            ylab = 'Affinity', xlab = 'Normal quantiles')
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

pdf('plots/correlationPlot.pdf', width = 12, height = 10)
corr_verify <- mul_reg%>%
  filter(NOSC < 0,
         !is.na(betterAffinity))%>%
  group_by(feature_number, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  select(betterAffinity, `row m/z`, NOSC, log_nc, log_oc, log_hc)%>%
  cor()%>%
  corrplot::corrplot()
dev.off()


# VIZUALIZATIONS -- labile pie charts -------------------------------------
labilePieCharts <- feature_stats_wdf%>%
  ungroup()%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')

labileCounts <- labilePieCharts%>%
  select(net_act, activity)%>%
  unique()%>%
  group_by(net_act)%>%
  mutate(actBin = case_when(activity %like% 'reca%' ~ 0,
                            # activity %like% 'acc%' ~ 0,
                            activity %like% 'labile' ~ 1))%>%
  summarize_if(is.numeric, sum)%>%
  mutate(activity = case_when(actBin == 0 ~ 'recalcitrant',
                              actBin == 1 ~ 'labile'),
         count = 1)%>%
  left_join(ecoNet%>%
              mutate(net_act = case_when(network == -1 ~ -as.numeric(scan),
                                         TRUE ~ network))%>%
              select(-c(scan, network))%>%
              unique(), by = 'net_act')%>%
  group_by(net_act)%>%
  nest()%>%
  mutate(data = map(data, ~mutate(.x, bad = nrow(.),
                                  bad = case_when(bad == 2 & is.na(matchSource) ~ 1000,
                                                  TRUE ~ as.numeric(bad)))%>%
                      filter(bad != 1000)%>%
                      select(-bad)))%>%
  unnest(data)%>%
  ungroup()%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  select(superclass_consensus, activity, count)%>%
  group_by(superclass_consensus, activity)%>%
  summarize_if(is.numeric, sum)%>%
  rename(labels = activity,
         parents = superclass_consensus)%>%
  mutate(parents = case_when(is.na(parents) ~ 'No consensus<br>reached',
                             parents %like% 'Lipids%' ~ 'Lipids and lipid-like<br>molecules',
                             parents %like% '%acids%' ~ 'Organic acids',
                             parents %like% '%oxygen%' ~ 'Oxygen compounds',
                             parents %like% '%nitrogen%' ~ 'Nitrogen compounds',
                             TRUE ~ as.character(parents)))

classifications <- labileCounts%>%
  group_by(parents)%>%
  select(parents, count)%>%
  summarize_if(is.numeric, sum)%>%
  rename(labels = parents)%>%
  mutate(parents = "")

labileSunburst <- labileCounts%>%
  group_by(parents)%>%
  mutate(total = round(count/sum(count), digits = 2)*100,
         idEnds = '%)')%>%
  filter(labels != 'recalcitrant')%>%
  unite(labels, c(count, total), sep = ' labile subnetworks<br>(', remove = FALSE)%>%
  unite(labels, c(labels, idEnds), sep = '')%>%
  mutate(labels = case_when(parents %like% '%Nitrog%' ~ '44 labile subnetworks (46%)',
                            parents %like% '%hetero%' ~ '68 labile subnetworks (46%)',
                            parents %like% '%Benz%' ~ '76 labile subnetworks (37%)',
                            parents %like% '%Oxyg%' ~ '80 labile subnetworks (17%)',
                            TRUE ~ labels))%>%
  ungroup()%>%
  bind_rows(classifications)%>%
  unite(ids, c(parents, labels), sep = ' - ', remove = FALSE)%>%
  mutate(ids = case_when(parents == "" ~ labels,
                         TRUE ~ as.character(ids)))%>%
  filter(!labels %like% 'No con%',
         !parents %like% 'No con%')
  
  
  
p <- plotly::plot_ly(labileSunburst, ids = ~ids, labels = ~labels, 
                parents = ~parents, values = ~count, type = 'sunburst',
                branchvalues = 'total', size = 0.8)%>%
  plotly::layout(p, colorway = c("#aaaaaa", "#888888", "#666666", "#444444", "#333333"))

p

plotly::save_image(p, "plots/sunburstLabilityNoNoCon.png", scale = 6)


# VISUALIZATIONS -- Sunburst xic ------------------------------------------------------------
xicSunburst <- labilePieCharts%>%
  filter(Timepoint == 'T0',
         DayNight == 'Day')%>%
  left_join(ecoNet%>%
              mutate(net_act = case_when(network == -1 ~ -as.numeric(scan),
                                         TRUE ~ network))%>%
              select(-c(scan, network))%>%
              unique()%>%
              group_by(net_act)%>%
              nest()%>%
              mutate(data = map(data, ~mutate(.x, bad = nrow(.),
                                              bad = case_when(bad == 2 & is.na(matchSource) ~ 1000,
                                                              TRUE ~ as.numeric(bad)))%>%
                                  filter(bad != 1000)%>%
                                  select(-bad)))%>%
              unnest(data)%>%
              ungroup(), by = 'net_act')%>%
  select(activity, ecoNetConsensus, Replicate, xic)%>%
  group_by(activity, ecoNetConsensus, Replicate)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  group_by(superclass_consensus, activity)%>%
  # mutate(xic = log10(xic))%>%
  summarize_if(is.numeric, mean)%>%
  rename(labels = activity,
         parents = superclass_consensus)%>%
  mutate(parents = case_when(is.na(parents) ~ 'No consensus<br>reached',
                             parents %like% 'Lipids%' ~ 'Lipids and lipid-like<br>molecules',
                             parents %like% '%acids%' ~ 'Organic acids',
                             parents %like% '%oxygen%' ~ 'Oxygen compounds',
                             parents %like% '%nitrogen%' ~ 'Nitrogen compounds',
                             TRUE ~ as.character(parents)))

sunburstParents <- xicSunburst%>%
  # mutate(xic = xic/10^8)%>%
  group_by(parents)%>%
  select(parents, xic)%>%
  summarize_if(is.numeric, sum)%>%
  rename(labels = parents)%>%
  mutate(parents = "")

xicSunburstLabile <- xicSunburst%>%
  # mutate(xic = xic/10^8)%>%
  group_by(parents)%>%
  mutate(total = round(xic/sum(xic), digits = 2)*100,
         idEnds = '%)')%>%
  filter(labels != 'recalcitrant')%>%
  mutate(labels = as.character(total))%>%
  ungroup()%>%
  bind_rows(sunburstParents)%>%
  unite(ids, c(parents, labels), sep = ' - ', remove = FALSE)%>%
  mutate(ids = case_when(parents == "" ~ labels,
                         TRUE ~ as.character(ids)))%>%
  filter(!labels %like% 'No con%',
         !parents %like% 'No con%')

g <- plotly::plot_ly(xicSunburstLabile, ids = ~ids, labels = ~labels, 
                     parents = ~parents, values = ~count, type = 'sunburst',
                     branchvalues = 'total', size = 0.8)%>%
  layout(g, colorway = c("#aaaaaa", "#888888", "#666666", "#444444", "#333333"))

g

save_image(g, "plots/sunburstLabilityXic.png", scale = 6)


# VIZUALIZATIONS -- Supplemental Figure Histogram -------------------------
pdf('plots/xicHistogram.pdf', width = 15, height = 10)
feature_stats_wdf%>%
  mutate(xic = xic + 1)%>%
  ggplot(aes(xic)) +
  geom_histogram(bins = 100)+
  scale_x_log10() +
  geom_vline(xintercept = 10^6, color = 'black', alpha = 0.4) +
  gen_theme()
dev.off()


# VIZUALIZATIONS -- RGB hex codes for orgs --------------------------------
org_colors_no_water <- c("#A30029","#669900", "#FF850A", "#9900FF", "#33CC33")



# Vizualizations -- Lability classifications of DOM -----------------------
lability_classes <- feature_stats_wdf%>%
  select(-c(log10:asin))%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  spread(Timepoint, xic)%>%
  mutate(val = T0)


lability_table <- lability_classes%>%
  group_by(DayNight, Replicate, Organism)%>%
  mutate(ra = T0/sum(T0, na.rm = TRUE))%>%
  ungroup()%>%
  group_by(activity, DayNight, superclass_consensus, Replicate)%>%
  filter(!is.na(val))%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, Replicate)%>%
  mutate(relativeT0Prod = T0/sum(T0, na.rm = TRUE),
         relativeT0Prod = case_when(is.na(relativeT0Prod) ~ 0,
                                    TRUE ~ relativeT0Prod),
         prorotionRa = ra/sum(ra, na.rm = TRUE),
         TIC = sum(T0, na.rm = TRUE))%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, activity)%>%
  filter(sum(T0) != 0,
         activity == 'labile')%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus)%>%
  # mutate(std = sd(relativeT0Prod))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  gather(reporter, val, relativeT0Prod:TIC)%>%
  select(DayNight, superclass_consensus, reporter, val)%>%
  unite(spreader, c('reporter', 'DayNight'), sep = '_')%>%
  spread(spreader, val)


lability_features <- lability_classes%>%
  ungroup()%>%
  filter(!is.na(val),
         Replicate == 1,
         network == -1)%>%
  select(activity, Replicate, net_act)%>%
  # select(activity, DayNight, superclass_consensus, Replicate, network)%>%
  mutate(count = 1)%>%
  unique()%>%
  group_by(activity, Replicate)%>%
  # group_by(activity, DayNight, superclass_consensus, Replicate)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(Replicate)%>%
  # group_by(DayNight, superclass_consensus, Replicate)%>%
  mutate(percentFeatures = count/sum(count, na.rm = TRUE),
         sumFeature = sum(count, na.rm = TRUE))%>%
  ungroup()
  filter(activity == 'labile')%>%
  # group_by(superclass_consensus)%>%
  # group_by(DayNight, superclass_consensus)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()
  gather(reporter, val, count:sumFeature)%>%
  select(superclass_consensus, reporter, val)%>%
  spread(reporter, val)
  # select(DayNight, superclass_consensus, reporter, val)%>%
  # unite(spreader, c('reporter', 'DayNight'), sep = '_')%>%
  # spread(spreader, val)

lability_all <- lability_table%>%
  filter(!superclass_consensus %like% 'Lignan%')%>%
  select(superclass_consensus, relativeT0Prod_Day, relativeT0Prod_Night)%>%
  gather(proportion, val, 2:3)%>%
  separate(proportion, c('response', 'DayNight'), sep = '_')%>%
  left_join(lability_features%>%
              select(superclass_consensus, percentFeatures_Day, percentFeatures_Night)%>%
              gather(features, fVal, 2:3)%>%
              separate(features, c('response', 'DayNight'), sep = '_'), by = c('superclass_consensus', 'DayNight'))
  # filter(!superclass_consensus %like% '%nitrogen%' | DayNight != 'Night')

percentLm <- lability_all%>%
  lm(fVal ~ val, data = .)%>%
  summary()

percentLmNoNiNi <- lability_all%>%
  filter(!superclass_consensus %like% '%itrogen%' | DayNight != 'Night')%>%
  lm(fVal ~ val, data = .)%>%
  summary()


pdf('plots/percentTICvsFeatures.pdf', width = 15, height = 10)
lability_all%>%
  ggplot(aes(val, fVal)) +
  geom_point(aes(color = DayNight), stat = 'identity', size = 3) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = c('#F09837', 'grey')) +
  labs(x = 'Labile percent of TIC', y = 'Labile percent of features') +
  geom_text(aes(x = 0.2, y = 0.9), size = 18, label = 'Pvalue: 0.0009') +
  geom_text(aes(x = 0.2, y = 0.8), size = 18,  label = 'R^2: 0.528') +
  # ylim(0.15,0.43) +
  gen_theme()

lability_all%>%
  filter(!superclass_consensus %like% '%nitrogen%' | DayNight != 'Night')%>%
  ggplot(aes(val, fVal)) +
  geom_point(aes(color = DayNight), stat = 'identity', size = 3) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = c('#F09837', 'grey')) +
  labs(x = 'Labile percent of TIC', y = 'Labile percent of features') +
  geom_text(aes(x = 0.2, y = .9), size = 18, label = 'Pvalue: 0.0006') +
  geom_text(aes(x = 0.2, y = 0.8), size = 18,  label = 'R^2: 0.580') +
  # ylim(0.15,0.43) +
  gen_theme()
dev.off()


superclassNetCount <- lability_classes%>%
  ungroup()%>%
  mutate(count = 1)%>%
  select(superclass_consensus, net_act, DayNight, count)%>%
  unique()%>%
  group_by(superclass_consensus, DayNight)%>%
  # filter(!superclass_consensus %like any% c('Lignans%', 'Phenylprop%', 'Alkaloids%', 'Nucleosides%'))%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()

superclassXicTotal <- lability_classes%>%
  ungroup()%>%
  group_by(network, net_act, DayNight, superclass_consensus, Organism, activity, feature_number)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(superclass_consensus, DayNight)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight)%>%
  mutate(relT0 = T0/sum(T0))%>%
  filter(!superclass_consensus %like any% c('Lignans%', 'Phenylprop%', 'Alkaloids%', 'Nucleosides%'))
  

# VIZUALIZATIONS -- Organismal production of lability classes -------------
org_lability <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  filter(Timepoint == 'T0',
         activity == 'labile')%>%
  select(-c(log10:asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  group_by(DayNight, superclass_consensus, Organism, Replicate)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  spread(DayNight,xic)%>%
  gather(DayNight,xic, Day:Night)%>%
  group_by(Organism, DayNight, superclass_consensus)%>%
  mutate(n = 1,
         n = sum(n),
         nxic = xic,
         xic = log10(xic + 1),
         std = sd(xic)/sqrt(n))%>%
  summarize_if(is.numeric, mean)

org_proportional_lability <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  filter(Timepoint == 'T0')%>%
  select(-c(log10:asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  group_by(DayNight, activity, Organism)%>%
  mutate(proportion = sum(xic))%>%
  ungroup()%>%
  group_by(DayNight, Organism)%>%
  mutate(proportion = proportion/sum(xic),
         n = 1,
         n = sum(n),
         std = sd(proportion)/sqrt(n))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()


pdf('plots/OrganismProportionLabile.pdf', width = 15, height = 10)
org_proportional_lability%>%
  mutate(Organism = as.factor(Organism),
         Organism = fct_relevel(Organism, c('CCA', 'Turf', 'Dictyota', 'Pocillopora verrucosa', 'Porites lobata')))%>%
  ggplot(aes(Organism, proportion, fill = DayNight)) +
  # geom_boxplot() +
  geom_bar(stat = 'identity', position = 'dodge2') +
  geom_errorbar(aes(ymax = proportion + std, ymin = proportion - std), position = 'dodge2') +
  scale_fill_manual(values = c('#B0ADAD', '#4D4D4D')) +
  # scale_color_manual(values = c('#F09837', 'grey')) +
  labs(y = 'Proportion of TIC', fill = 'Diel Cycle') +
  scale_y_continuous(limits = c(0.45,0.65), oob = rescale_none) +
  gen_theme()
dev.off()  

org_tic <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  filter(!superclass_consensus %like any% c('Lignan%', 'Alkaloid%', 'Phenyl%'))%>%
  filter(Timepoint == 'T0',
         activity == 'labile')%>%
  select(-c(log10:asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(data = map(data, ~filter(.x, !is.na(xic))%>%
                      group_by(Organism, Replicate)%>%
                      summarize_if(is.numeric, sum, na.rm = TRUE)%>%
                      ungroup()
  ))%>%
  unnest(data)%>%
  spread(DayNight,xic)%>%
  gather(DayNight,xic, Day:Night)%>%
  left_join(org_proportional_lability%>%
              select(-c(xic,n,std)), by = c('DayNight', 'Organism'))%>%
  group_by(Organism, DayNight)%>%
  mutate(propXic = log10(xic*proportion +1),
         log10 = log10(xic +1),
         std = sd(xic))%>%
  summarize_if(is.numeric, mean)%>%
  mutate(xic = format(xic, scientific = TRUE))

reOrg_colors_no_water <- c("#A30029","#33CC33", "#669900", "#FF850A", "#9900FF")

minValues <- c(5, 8, 7, 6, 0, 6)
maxValues <- c(8.2, 10.2, 8, 8, 8.1, 8)

org_lability_limites <- org_lability%>%
  ungroup()%>%
  select(superclass_consensus)%>%
  filter(superclass_consensus %like any% c('Benz%', 'Lipi%', '%acids%', '%nitr%', '%ox%', '%hetero%'))%>%
  unique()%>%
  bind_cols(minValues, maxValues)%>%
  rename('min' = 2,
         'max' = 3)

mapSuperclasses <- function(x) {
  ggplot(x, aes(Organism, xic, fill = DayNight)) +
    # geom_bar_pattern(stat = 'identity', position = position_dodge(preserve = 'single')) +
    geom_bar(stat = 'identity', position = 'dodge2') + 
    geom_errorbar(aes(ymax = xic + std, ymin = xic - std), position = 'dodge2') +
    # facet_wrap(~superclass_consensus, ncol = 1, scales = 'free_y') +
    scale_fill_manual(values = c('#B0ADAD', '#4D4D4D')) +
    labs(y = 'Ion intensity (XIC)', color = 'Organism', lineType = 'Diel Cycle') +
    gen_theme() +
    theme(axis.text.y =element_text(size = 35))
}


pdf('plots/organismSuperclassLability.pdf', width = 15, height = 10)
org_lability%>%
  mutate(Organism = as.factor(Organism),
         Organism = fct_relevel(Organism, c('CCA', 'Turf', 'Dictyota', 'Pocillopora verrucosa', 'Porites lobata')))%>%
  # filter(!superclass_consensus %like any% c('Lignan%', 'Alkaloid%', 'Phenyl%'))%>%
  filter(superclass_consensus %like any% c('Benz%'))%>%
  left_join(org_lability_limites, by = 'superclass_consensus')%>%
  mapSuperclasses  +
  scale_y_continuous(limits = c(5,8.2), oob = rescale_none)

org_lability%>%
  mutate(Organism = as.factor(Organism),
         Organism = fct_relevel(Organism, c('CCA', 'Turf', 'Dictyota', 'Pocillopora verrucosa', 'Porites lobata')))%>%
  # filter(!superclass_consensus %like any% c('Lignan%', 'Alkaloid%', 'Phenyl%'))%>%
  filter(superclass_consensus %like any% c( 'Lipi%'))%>%
  left_join(org_lability_limites, by = 'superclass_consensus')%>%
  mapSuperclasses() +
  scale_y_continuous(limits = c(8,10.2), oob = rescale_none)

org_lability%>%
  mutate(Organism = as.factor(Organism),
         Organism = fct_relevel(Organism, c('CCA', 'Turf', 'Dictyota', 'Pocillopora verrucosa', 'Porites lobata')))%>%
  # filter(!superclass_consensus %like any% c('Lignan%', 'Alkaloid%', 'Phenyl%'))%>%
  filter(superclass_consensus %like any% c('%acids%'))%>%
  left_join(org_lability_limites, by = 'superclass_consensus')%>%
  mapSuperclasses() +
  scale_y_continuous(limits = c(7,10.2), oob = rescale_none)

org_lability%>%
  mutate(Organism = as.factor(Organism),
         Organism = fct_relevel(Organism, c('CCA', 'Turf', 'Dictyota', 'Pocillopora verrucosa', 'Porites lobata')))%>%
  # filter(!superclass_consensus %like any% c('Lignan%', 'Alkaloid%', 'Phenyl%'))%>%
  filter(superclass_consensus %like any% c('%nitr%'))%>%
  left_join(org_lability_limites, by = 'superclass_consensus')%>%
  mapSuperclasses() +
  scale_y_continuous(limits = c(6,9), oob = rescale_none)

org_lability%>%
  mutate(Organism = as.factor(Organism),
         Organism = fct_relevel(Organism, c('CCA', 'Turf', 'Dictyota', 'Pocillopora verrucosa', 'Porites lobata')))%>%
  # filter(!superclass_consensus %like any% c('Lignan%', 'Alkaloid%', 'Phenyl%'))%>%
  filter(superclass_consensus %like any% c('%ox%'))%>%
  left_join(org_lability_limites, by = 'superclass_consensus')%>%
  mapSuperclasses() +
  scale_y_continuous(limits = c(0,8.1), oob = rescale_none)

org_lability%>%
  mutate(Organism = as.factor(Organism),
         Organism = fct_relevel(Organism, c('CCA', 'Turf', 'Dictyota', 'Pocillopora verrucosa', 'Porites lobata')))%>%
  # filter(!superclass_consensus %like any% c('Lignan%', 'Alkaloid%', 'Phenyl%'))%>%
  filter(superclass_consensus %like any% c('%hetero%'))%>%
  left_join(org_lability_limites, by = 'superclass_consensus')%>%
  mapSuperclasses() +
  scale_y_continuous(limits = c(6,9.2), oob = rescale_none)
  
dev.off()

pdf('plots/OrganismLabileXIC.pdf', width = 15, height = 10)
org_tic%>%
  mutate(pattern = "add")%>%
  group_by(Organism, DayNight)%>%
  summarize_if(is.numeric, mean)%>%
  mutate(Organism = as.factor(Organism),
         Organism = fct_relevel(Organism, c('CCA', 'Turf', 'Dictyota', 'Pocillopora verrucosa', 'Porites lobata')))%>%
  ggplot(aes(Organism, xic, fill = DayNight, color = DayNight)) + 
  geom_bar(stat = 'identity', position = 'dodge2') + 
  geom_errorbar(aes(ymax = xic + std, ymin = xic - std), position = 'dodge2') +
  # facet_wrap(~superclass_consensus, scales = 'free_y') +
  # geom_bar(aes(Organism, propXic, fill = 'labile'), stat = 'identity', position = 'dodge2') +
  geom_col_pattern(aes(y=propXic), pattern_fill = 'black', position = 'dodge2') +
  # scale_pattern_manual(values = c('none', 'none', 'crosshatch')) +
  # geom_bar_pattern(aes(Organism, propXic), stat = 'identity', position = 'dodge2', color = "black", 
  #                  pattern_fill = "black",
  #                  pattern_angle = 45,
  #                  pattern_density = 0.1,
  #                  pattern_spacing = 0.025,
  #                  pattern_key_scale_factor = 0.6) + 
  # scale_pattern_manual(values = c(add = "stripe")) +
  scale_fill_manual(values = c('#B0ADAD', '#4D4D4D')) +
  scale_color_manual(values = c('#B0ADAD', '#4D4D4D')) +
  labs(y = 'Ion intensity (XIC)', fill = 'Diel Cycle') +
  scale_y_continuous(limits = c(9.3, 10.6), oob = rescale_none) +
  # scale_y_log10()+
  gen_theme()

dev.off()


OrganismalLabileFeatues <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity == 'labile')%>%
  ungroup()%>%
  select(feature_number)%>%
  unique()
  # group_by(DayNight)%>%
  # mutate(n = 1)%>%
  # summarize_if(is.numeric, sum)%>%
  # ungroup()
  

# STATS -- Org proportional lability between superclasses -----------------
org_proportional_lability_superclass <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  filter(Timepoint == 'T0')%>%
  select(-c(log10:asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  group_by(DayNight, Replicate,  activity, Organism, superclass_consensus)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  group_by(DayNight, Replicate, Organism, superclass_consensus)%>%
  mutate(proportion = xic/sum(xic))%>%
  filter(activity == 'labile')%>%
  select(-xic)%>%
  filter(!superclass_consensus %like any% c('Alk%', 'Lig%'))%>%
  spread(superclass_consensus, proportion)

proportionalPermanova <- org_proportional_lability_superclass%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(data = map(data, ~ adonis2(.x[4:ncol(.x)] ~ Organism, data = .x, perm = 100000, method = "bray", na.rm = TRUE)))


# STATS -- Organism TIC aov/tukey -----------------------------------------
orgTicStats <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  # filter(activity == 'labile')
  filter(activity != 'accumolite')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  filter(!superclass_consensus %like any% c('Lignan%', 'Alkaloid%', 'Phenyl%'))%>%
  filter(Timepoint == 'T0',
         activity == 'labile')%>%
  select(-c(log10:asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(data = map(data, ~filter(.x, !is.na(xic))%>%
                      group_by(Organism, Replicate)%>%
                      summarize_if(is.numeric, sum, na.rm = TRUE)%>%
                      ungroup()))%>%
  unnest(data)%>%
  spread(DayNight,xic)%>%
  gather(DayNight,xic, Day:Night)%>%
  left_join(org_proportional_lability%>%
              select(-c(xic, n,std)), by = c('DayNight', 'Organism'))%>%
  group_by(Organism, DayNight)%>%
  mutate(propXic = log10(xic*proportion +1),
         xic = log10(xic +1),
         std = sd(xic))%>%
  ungroup()%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'))%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(anova = map(data, ~aov(xic~Organism, .x)%>%
                       tidy()),
         tukey = map(data, ~aov(xic~Organism, .x)%>%
                       TukeyHSD()),
         proportionAnova = map(data, ~aov(proportion~Organism, .x)%>%
                                 tidy()),
         proportionTukey = map(data, ~aov(proportion~Organism, .x)%>%
                                 TukeyHSD()))


# STATS -- LMER -- community/lability change coral vs mixed vs algae --------------
labileChangeStats <- feature_stats_wdf%>%
  ungroup()%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  filter(activity == 'labile',
         Timepoint == 'T0')%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')


## Superclass and TIC
groupedStatsLabile <- labileChangeStats%>%
    filter(Organism != 'CCA')%>%
    mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'))%>%
    group_by(group, Organism, Replicate, DayNight, superclass_consensus)%>%
    summarize_if(is.numeric, sum, na.rm = TRUE)%>%
    ungroup()%>%
    mutate(xic = log10(xic + 1),
           test = 'superclass')

TICStats <- labileChangeStats%>%
  filter(Organism != 'CCA')%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'))%>%
  group_by(group, Organism, Replicate, DayNight)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  mutate(xic = log10(xic),
         test = 'TIC')

figSixstats <- groupedStatsLabile%>%
  filter(!superclass_consensus %like any% c('Alk%', 'Lig%', 'Pheny%'),
         !is.na(superclass_consensus))%>%
  bind_rows(TICStats)%>%
  group_by(test, superclass_consensus)%>%
  nest()%>%
  # mutate(data = map(data, ~lm(xic~group*DayNight, data = .x)%>%
  #                     tidy()))%>%
  # unnest(data)
  mutate(data = map(data, ~ lmer(xic~group*DayNight + (1|Organism), data =.x)%>%
                      car::Anova()%>%
                      rownames_to_column(var = 'var')%>%
                      select(var, `Pr(>Chisq)`)))%>%
  unnest(data)


selectSubclass <- labileChangeStats%>%
  filter(Organism != 'CCA')%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'))%>%
  # group_by(group, Organism, Replicate, DayNight, superclass_consensus,class_consensus, subclass_consensus, network, net_act)%>%
  # summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  # ungroup()%>%
  mutate(xic = log10(xic + 1),
         test = 'subclass')%>%
  filter(!superclass_consensus %like any% c('Alk%', 'Lig%', 'Pheny%'),
         !is.na(superclass_consensus))%>%
  distinct(superclass_consensus,class_consensus, subclass_consensus, group)%>%
  count(superclass_consensus, class_consensus, subclass_consensus)%>%
  filter(n > 1)%>%
  select(-n)
  # pull(uperclass_consensus, subclass_consensus)

subclassLmer <- labileChangeStats%>%
  filter(Organism != 'CCA')%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'))%>%
  group_by(group, Organism, Replicate, DayNight, superclass_consensus, class_consensus, subclass_consensus, feature_number, network, net_act)%>%
  # summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  # ungroup()%>%
  mutate(xic = log10(xic + 1),
         test = 'subclass')%>%
  filter(!superclass_consensus %like any% c('Alk%', 'Lig%', 'Pheny%'),
         !is.na(superclass_consensus),
         subclass_consensus %like any% c('Fatty acid es%', 'Fatty ami%', 'Terpene lact%', 'Tricarboxylic ac%', 'Carbohydrate%','1,3,5%') | 
           class_consensus %like% 'Fatty Acyl%' & is.na(subclass_consensus)|
           class_consensus %like% 'Carboxylic aci%' & is.na(subclass_consensus))%>%
  # inner_join(selectSubclass, by = c('superclass_consensus', 'class_consensus', 'subclass_consensus'))%>%
  # group_by(test, DayNight, superclass_consensus, subclass_consensus)%>%
  # mutate(numGroup = length(unique(group)),
  #        numOrg = length(unique(Organism)))%>%
  # filter(numGroup > 1,
  #        numOrg > 3)%>%
  mutate(Organism = as.factor(Organism),
         group = as.factor(group),
         DayNight = as.factor(DayNight))%>%
  ungroup()%>%
  group_by(test, superclass_consensus, class_consensus, subclass_consensus)%>%
  nest()
  mutate(data = map(data, ~ lmer(xic~group*DayNight + (1|Organism) + (1|feature_number), ., control = lmerControl(check.nlev.gtr.1 = "ignore",
                                                                                          check.conv.singular = 'ignore',
                                                                                          check.nobs.vs.nRE = 'ignore'))%>%
                      car::Anova()%>%
                      rownames_to_column(var = 'var')%>%
                      select(var, `Pr(>Chisq)`)))%>%
  unnest(data)

write_csv(subclassLmer, '~/Documents/GitHub/DORCIERR/data/analysis/subclassLMERResults.csv')

# Microbiak community dynamics
growthLoadStats <- fcm_23%>%
  filter(Organism != 'CCA')%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'),
         val = 'microbe')%>%
  group_by(val)%>%
  nest()%>%
  mutate(load = map(data, ~ lmer(final_cells~group*DayNight + (1|Organism), data =., control =  lmerControl(check.nlev.gtr.1 = "ignore",
                                                                                                    check.conv.singular = 'ignore',
                                                                                                    check.nobs.vs.nRE = 'ignore'))%>%
                      car::Anova()%>%
                      rownames_to_column(var = 'var')%>%
                      select(var, `Pr(>Chisq)`)%>%
                      mutate(test = 'load')),
         growth = map(data, ~ lmer(cells_ul~group*DayNight + (1|Organism), data =., control =  lmerControl(check.nlev.gtr.1 = "ignore",
                                                                                                              check.conv.singular = 'ignore',
                                                                                                              check.nobs.vs.nRE = 'ignore'))%>%
                        car::Anova()%>%
                        rownames_to_column(var = 'var')%>%
                        select(var, `Pr(>Chisq)`)%>%
                        mutate(test = 'growth')))

growthLoadPVals <- growthLoadStats%>%
  select(-c(data, growth))%>%
  unnest(load)%>%
  bind_rows(growthLoadStats%>%
              select(-c(data, load))%>%
              unnest(growth))%>%
  ungroup()%>%
  select(-val)
  # mutate(FDR = p.adjust(`Pr(>Chisq)`, method = 'BH'))


## Depletion ∆
depletionStats <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity == 'labile')%>%
  group_by(Organism, Timepoint, activity, Replicate, DayNight)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  select(-c(log10:asin, network, net_act, ra))%>%
  spread(Timepoint, xic)%>%
  filter(Organism != 'CCA')%>%
  group_by(Organism, DayNight)%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'),
         T0 = mean(T0, na.rm = TRUE),
         difference = abs(TF - T0),
         log10Difference = log10(difference),
         proportionalDiff = (TF-T0)/T0)%>%
  mutate(val = 'val')%>%
  group_by(val)%>%
  nest()%>%
  mutate(difference = map(data, ~ lmer(log10Difference~group*DayNight + (1|Organism), data =., control =  lmerControl(check.nlev.gtr.1 = "ignore",
                                                                               check.conv.singular = 'ignore',
                                                                               check.nobs.vs.nRE = 'ignore'))%>%
                      car::Anova()%>%
                      rownames_to_column(var = 'var')%>%
                      select(var, `Pr(>Chisq)`)%>%
                      mutate(test = 'depletion')),
         proportionalDifference = map(data, ~ lmer(proportionalDiff~group*DayNight + (1|Organism), data =., control =  lmerControl(check.nlev.gtr.1 = "ignore",
                                                                                                                      check.conv.singular = 'ignore',
                                                                                                                      check.nobs.vs.nRE = 'ignore'))%>%
                            car::Anova()%>%
                            rownames_to_column(var = 'var')%>%
                            select(var, `Pr(>Chisq)`)%>%
                            mutate(test = 'proportionalDepletion')))
                      
  

depletionStatsCombined <- depletionStats%>%
  select(-c(data, proportionalDifference))%>%
  unnest(difference)%>%
  bind_rows(depletionStats%>%
              select(-c(data, difference))%>%
              unnest(proportionalDifference))%>%
  ungroup()%>%
  select(-val)
  # mutate(FDR = p.adjust(`Pr(>Chisq)`, method = 'BH'))



# -- ANALYSIS - lmer proportionally higher --------------------------------
## Superclass and TIC
groupedStatsLabileLarger <- labileChangeStats%>%
  filter(Organism != 'CCA')%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'))%>%
  group_by(group, Organism, Replicate, DayNight, superclass_consensus)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  mutate(xic = log10(xic + 1),
         test = 'superclass')

#Superclassifications
superClassDayNight <- groupedStatsLabileLarger%>%
  select(superclass_consensus, DayNight, xic)%>% 
  group_by(superclass_consensus, DayNight)%>% 
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  spread(DayNight, xic)%>%
  mutate(larger = Day - Night,
         larger = ifelse(larger > 0, 'Day', 'Night'),
         var = 'DayNight',
         test = 'superclass')%>%
  select(-c(Day, Night))%>%
  filter(superclass_consensus %like any% c('Lipid%', '%nitrogen%', '%heterocy%'))

superClassInteraction <- groupedStatsLabileLarger%>%
  select(superclass_consensus, group, DayNight, xic)%>% 
  group_by(superclass_consensus, group, DayNight)%>% 
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  spread(group, xic)%>%
  mutate(larger = Coral - `Fleshy algae`,
         larger = ifelse(larger > 0, 'coral', 'Fleshy algae'),
         var = 'group:DayNight',
         test = 'superclass')%>%
  select(-c(Coral, `Fleshy algae`))%>%
  filter(superclass_consensus %like any% c('Lipid%', '%heterocy%'))%>%
  spread(DayNight, larger)

TICInteraction <- groupedStatsLabileLarger%>%
  filter(Organism != 'CCA')%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'))%>%
  group_by(group, Organism, Replicate, DayNight)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  mutate(xic = log10(xic),
         test = 'TIC')%>%
  select(test, group, DayNight, xic)%>%
  group_by(test, group, DayNight)%>% 
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  spread(group, xic)%>%
  mutate(larger = Coral - `Fleshy algae`,
         larger = ifelse(larger > 0, 'coral', 'Fleshy algae'),
         var = 'group:DayNight')%>%
  select(-c(Coral, `Fleshy algae`))%>%
  spread(DayNight, larger)

# Microbial community dynamics
growthLoadPostStats <- fcm_23%>%
  filter(Organism != 'CCA')%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'),
         val = 'microbe')%>%
  group_by(val)%>%
  nest()%>%
  mutate(loadInteraction = map(data, ~select(.x, group, DayNight, final_cells)%>%
                                 group_by(group, DayNight)%>% 
                                 summarize_if(is.numeric, mean, na.rm = TRUE)%>%
                                 spread(group, final_cells)%>%
                                 mutate(larger = Coral - `Fleshy algae`,
                                        larger = ifelse(larger > 0, 'coral', 'Fleshy algae'),
                                        var = 'group:DayNight')%>%
                                 select(-c(Coral, `Fleshy algae`))%>%
                                 spread(DayNight, larger)%>%
                                 mutate(test = 'load')),
         loadDayNight = map(data, ~select(.x, DayNight, final_cells)%>%
                              group_by(DayNight)%>% 
                              summarize_if(is.numeric, mean, na.rm = TRUE)%>%
                              spread(DayNight, final_cells)%>%
                              mutate(larger = Day - Night,
                                     larger = ifelse(larger > 0, 'Day', 'Night'),
                                     var = 'DayNight',
                                     test = 'load')%>%
                              select(-c(Day, Night))),
         growthDayNight = map(data, ~select(.x, DayNight, cells_ul)%>%
                                group_by(DayNight)%>% 
                                summarize_if(is.numeric, mean, na.rm = TRUE)%>%
                                spread(DayNight, cells_ul)%>%
                                mutate(larger = Day - Night,
                                       larger = ifelse(larger > 0, 'Day', 'Night'),
                                       var = 'DayNight',
                                       test = 'growth')%>%
                                select(-c(Day, Night))),
         growthInteraction = map(data, ~select(.x, group, DayNight, cells_ul)%>%
                                   group_by(group, DayNight)%>% 
                                   summarize_if(is.numeric, mean, na.rm = TRUE)%>%
                                   spread(group, cells_ul)%>%
                                   mutate(larger = Coral - `Fleshy algae`,
                                          larger = ifelse(larger > 0, 'coral', 'Fleshy algae'),
                                          var = 'group:DayNight')%>%
                                   select(-c(Coral, `Fleshy algae`))%>%
                                   spread(DayNight, larger)%>%
                                   mutate(test = 'growth')))
  

## Depletion ∆
depletionPostStats <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity == 'labile')%>%
  group_by(Organism, Timepoint, activity, Replicate, DayNight)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  select(-c(log10:asin, network, net_act, ra))%>%
  spread(Timepoint, xic)%>%
  filter(Organism != 'CCA')%>%
  group_by(Organism, DayNight)%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'),
         T0 = mean(T0, na.rm = TRUE),
         difference = abs(TF - T0),
         log10Difference = log10(difference),
         proportionalDiff = (TF-T0)/T0)%>%
  mutate(val = 'val')%>%
  group_by(val)%>%
  nest()%>%
  mutate(depletionInteraction = map(data, ~select(.x, group, DayNight, log10Difference)%>%
                           group_by(group, DayNight)%>% 
                           summarize_if(is.numeric, mean, na.rm = TRUE)%>%
                           spread(group, log10Difference)%>%
                           mutate(larger = Coral - `Fleshy algae`,
                                  larger = ifelse(larger > 0, 'coral', 'Fleshy algae'),
                                  var = 'group:DayNight')%>%
                           select(-c(Coral, `Fleshy algae`))%>%
                           spread(DayNight, larger)%>%
                           mutate(test = 'depletion')),
         proportionDepletionDayNight = map(data, ~select(.x, DayNight, proportionalDiff)%>%
                                             group_by(DayNight)%>% 
                                             summarize_if(is.numeric, mean, na.rm = TRUE)%>%
                                             spread(DayNight, proportionalDiff)%>%
                                             mutate(larger = Day - Night,
                                                    larger = ifelse(larger > 0, 'Day', 'Night'),
                                                    var = 'DayNight',
                                                    test = 'proportionalDepletion')%>%
                                             select(-c(Day, Night))))


lmerComparisons<- superClassDayNight%>%
  bind_rows(superClassInteraction,
            TICInteraction,
            growthLoadPostStats$loadInteraction,
            growthLoadPostStats$loadDayNight,
            growthLoadPostStats$growthDayNight,
            growthLoadPostStats$growthInteraction,
            depletionPostStats$depletionInteraction,
            depletionPostStats$proportionDepletionDayNight)%>%
  rename(DayNight = 'larger')


allLmerStats <- figSixstats %>%
  bind_rows(growthLoadPVals, depletionStatsCombined) %>%
  ungroup()%>%
  left_join(lmerComparisons, by = c('superclass_consensus', 'var', 'test'))
# mutate(fdr = p.adjust(`Pr(>Chisq)`, method = 'BH'))

write_csv(allLmerStats, '~/Documents/GitHub/DORCIERR/data/analysis/lmerStatistics.csv')


# STATS -- Org lability t-tests-------------------------------------------
org_lability_stats <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  filter(Timepoint == 'T0',
         activity == 'labile')%>%
  select(-c(log10:asin, network, net_act, ecoNetConsensusScore, numberOfNodes))%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(data = map(data, ~filter(.x, !is.na(xic))%>%
                      group_by(superclass_consensus, Organism, Replicate)%>%
                      summarize_if(is.numeric, sum, na.rm = TRUE)%>%
                      ungroup()
  ))%>%
  unnest(data)%>%
  spread(DayNight,xic)%>%
  gather(DayNight,xic, Day:Night)

org_lability_test <- org_lability_stats%>%
  filter(!superclass_consensus %like any% c('Lignan%', 'Alkaloid%', 'Phenyl%'))%>%
  mutate(xic = case_when(is.na(xic) ~ 0,
                         TRUE ~ xic),
         log10 = log10(xic + 1))%>%
  group_by(superclass_consensus, Organism)%>%
  nest()%>%
  mutate(data = map(data, ~ lm(log10 ~ DayNight, data = .x)%>%
                      tidy()%>%
                      filter(term == 'DayNightNight')%>%
                      select(p.value)))%>%
  unnest(data)%>%
  mutate(FDR = p.adjust(p.value, method = 'BH'))

org_lability_aov <- org_lability_stats%>%
  filter(!superclass_consensus %like any% c('Lignan%', 'Alkaloid%', 'Phenyl%'))%>%
  mutate(xic = case_when(is.na(xic) ~ 0,
                         TRUE ~ xic),
         log10 = log10(xic + 1))%>%
  group_by(superclass_consensus, DayNight)%>%
  nest()%>%
  mutate(aov = map(data, ~ aov(log10 ~ Organism, data = .x)%>%
                      tidy()%>%
                      filter(term == 'Organism')%>%
                      select(p.value)),
         tukey = map(data, ~aov(log10 ~ Organism, data = .x)%>%
                       TukeyHSD()%>%
                       tidy()%>%
                       select(comparison, adj.p.value)))%>%
  unnest(c(aov, tukey))%>%
  unnest(data)%>%
  mutate(aovFDR = p.adjust(p.value, method = 'BH'))
  filter(aovFDR < 0.05)

# VIZUALIZATIONS -- Metabolite pool XIC change ----------------------------
metabolitePool_xic_change <- feature_stats_wdf%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')%>%
  group_by(Organism, Timepoint, activity, Replicate, DayNight)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  select(-c(log10:asin))%>%
  group_by(Organism, activity, Timepoint, DayNight)%>%
  mutate(err = sd(xic))%>%
  ungroup()%>%
  unite(org_activity, c('Organism', 'activity'), sep = ' / ', remove = FALSE)%>%
  mutate(shape = case_when(Timepoint == 'T0' ~ 'T0',
                           Timepoint == 'TF' & activity == 'labile' ~ 'TF depletion',
                           # Timepoint == 'TF' & activity == 'accumolite' ~ 'TF accumulation',
                           Timepoint == 'TF' & activity == 'recalcitrant' & Organism == 'Turf' ~ 'accumulation',
                           TRUE ~ 'TF depletion'))%>%
  group_by(Organism, Replicate, Timepoint)%>%
  mutate(ra = xic/sum(xic, na.rm = TRUE))%>%
  ungroup()

metabolitePool_meanChangeNight <- metabolitePool_xic_change%>%
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
  filter(activity == 'labile',
         DayNight == 'Night')%>%
  unite(group, c('Organism', 'DayNight'), sep = '_', remove = FALSE)

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
  filter(activity == 'labile')%>%
  unite(group, c('Organism', 'DayNight'), sep = '_', remove = FALSE)

write_csv(metabolitePool_meanChange, 'analysis/labilePoolChange.csv')

png("plots/change_xic_activityDayNight.png", width = 2000, height = 1200)
# cairo_pdf("plots/change_xic_activityDay.pdf", width = 20, height = 12)
metabolitePool_xic_change%>%
  filter(activity == 'labile')%>%
  unite(group, c('Organism', 'DayNight'), sep = '_', remove = FALSE)%>%
  mutate(log10 = log10(xic + 1),
         group = factor(group, levels = c('Porites lobata_Night', 'Porites lobata_Day', 'Pocillopora verrucosa_Night', 'Pocillopora verrucosa_Day', 'Dictyota_Night', 'Dictyota_Day', 'Turf_Night', 'Turf_Day', 'CCA_Night', 'CCA_Day')))%>%
  ggplot(aes(log10, group, color = DayNight, shape = shape)) +
  geom_point(stat = 'identity', size = 18, alpha = 0.45) +
  geom_line(aes(group = group)) +
  scale_color_manual(values = c('#B0ADAD', '#4D4D4D')) +
  scale_shape_manual(values = c("\u25A0", "\u25C4")) +
  # geom_text(data = metabolitePool_meanChange, aes(x = x_val, y = group, 
  #                                                 label = formatC(difference, format = 'e', digits = 2), 
  #                                                 shape = NULL), vjust = -0.2, size = 20) +
  labs(x = 'Sum Intensity (xic)', y = 'Organism', color = 'Diel Cycle: ', shape = 'Sample:') +
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
         shape = guide_legend(nrow=2, byrow=TRUE))
  # xlim(0, 3.5e10) 
  # scale_x_log10()
dev.off()

# STATS -- labile pool xic mean change -------------------------------
xic_change_ttest <- metabolitePool_xic_change%>%
  filter(activity == 'labile')%>%
  select(-c(network:ra))%>%
  spread(Timepoint, xic)%>%
  group_by(Organism, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         difference = TF-T0,
         log10Difference = log10(abs(difference)))%>%
  ungroup()%>%
  group_by(Organism)%>%
  nest()%>%
  mutate(p = map(data, ~ lm(log10Difference ~ DayNight, .x)%>%
                              lmp()))%>%
  # select(-data)%>%
  # unnest(c(DayDepletion, NightDepletion))%>%
  # gather(test, p, 2:3)%>%
  mutate(fdr = p.adjust(p, method = 'BH'))

xicChangelmDF <- metabolitePool_xic_change%>%
  filter(activity == 'labile')%>%
  select(-c(network:ra))%>%
  spread(Timepoint, xic)%>%
  group_by(Organism, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         difference = TF-T0,
         log10Difference = log10(abs(difference)))%>%
  ungroup()

xicChangelm <- xicChangelmDF%>%
  lm(log10Difference~Organism*DayNight, .)%>%
  lmp()

xicChangeTukey <- xicChangelmDF%>%
  aov(log10Difference~Organism, .)%>%
  TukeyHSD()%>%
  tidy()


# xic_T0_ttest <- metabolitePool_xic_change%>%
#   filter(activity == 'labile')%>%
#   select(-c(network:ra))%>%
#   filter(Timepoint == 'T0')%>%
#   mutate(log10 = log10(xic))%>%
#   group_by(Organism)%>%
#   nest()%>%
#   mutate(p = map(data, ~ lm(log10 ~ DayNight, .x)%>%
#                    lmp()))%>%
#   # select(-data)%>%
#   # unnest(c(DayDepletion, NightDepletion))%>%
#   # gather(test, p, 2:3)%>%
#   mutate(fdr = p.adjust(p, method = 'BH'))
# 
# 
xic_tukey <- metabolitePool_xic_change%>%
  filter(activity == 'labile')%>%
  select(-c(network:ra))%>%
  spread(Timepoint, xic)%>%
  group_by(Organism, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         difference = TF-T0,
         log10Difference = log10(abs(difference)))%>%
  ungroup()%>%
  aov(log10Difference ~ Organism, .)%>%
  TukeyHSD()%>%
  tidy()
  # group_by(DayNight)%>%
  # nest()%>%
  # mutate(data = map(data, ~ aov(log10Difference ~ Organism, .x)%>%
  #                  TukeyHSD()%>%
  #                  tidy()))%>%
  # unnest(data)

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
  filter(Organism != 'Offshore',
         Organism != 'Influent')%>%
  mutate(Replicate = 1)%>%
  # select(-c(final_cells, cells_ul))%>%
  mutate(Replicate = as.character(Replicate))

# organism_order_sgr <- as.factor(fcm_T0_5$Organism)%>%
#   relevel(., c('Turf', 'Dictyota', 'Pocillopora Verrucosa', 'CCA'))%>%
#   levels()%>%
#   as.vector()

sgr_colors <- c(org_colors_no_water, "#3B9AB2")

pdf("~/Documents/GitHub/DORCIERR/data/plots/SpecificGrowthRate.pdf", width = 7, height = 7)
fcm_T0_5%>%
  filter(Organism != 'Offshore',
         Organism != 'Influent')%>%
  # left_join(fcm_err, 
  #           by = c('Organism', 'Replicate'))%>%
  mutate(Organism = factor(Organism, levels = c('Water control', 'Porites lobata', 'CCA', 'Pocillopora verrucosa', 'Dictyota', 'Turf')))%>%
  ggplot(aes(cells_ul, final_cells, color = Organism)) +
  # geom_boxplot() +
  geom_errorbarh(data = fcm_err, aes(xmin = cells_ul - cells_err, xmax = cells_ul + cells_err)) +
  geom_errorbar(data = fcm_err, aes(ymin = final_cells - final_err, ymax = final_cells + final_err)) +
  geom_point(aes(shape = DayNight), stat = 'identity', size = 7, alpha = 0.6) +
  # coord_flip() +
  scale_color_manual(values = c(sgr_colors)) +
  scale_shape_manual(values = c(16, 1)) +
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

# VIZUALIZATIONS -- Linear model: Sum depletion growth rate first 24 hours --------------
metabolitePool_labileChange <- metabolitePool_xic_change%>%
  filter(activity == 'labile')%>%
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
  left_join(metabolitePool_labileChange,  by = c('Organism', 'DayNight', 'Replicate'))%>%
  mutate(log10 = log10(abs(difference)))%>%
  filter(Organism != 'Influent',
         Organism != 'Water control',
         Organism != 'Offshore')

# Day lm
fcm_labile_lm <- fcm_23%>%
  filter(DayNight == 'Day')%>%
  lm(cells_ul ~ log10, data = .)

fcd_p <- (fcm_labile_lm%>% 
            tidy()%>% 
            filter(term == 'log10'))$p.value
# 
# fcd_f <- (fcm_labile_lm%>% 
#             tidy()%>% 
#             filter(term == 'log10'))$statistic
# 
# fcd_slope <- fcm_labile_lm$coefficients["log10"]
# fcd_intercept <- fcm_labile_lm$coefficients["(Intercept)"]
# 
# 
# fcd_r2 <- summary(fcm_labile_lm)$adj.r.squared

#Day Corals
fcm_labileCorals_lm <- fcm_23%>%
  filter(DayNight == 'Day',
         Organism %like any% c('Porites%', 'Poci%'))%>%
  lm(cells_ul ~ log10, data = .)

fcdCoral_p <- (fcm_labileCorals_lm%>% 
            tidy()%>% 
            filter(term == 'log10'))$p.value

fcdCoral_f <- (fcm_labileCorals_lm%>% 
            tidy()%>% 
            filter(term == 'log10'))$statistic

fcdCoral_slope <- fcm_labileCorals_lm$coefficients["log10"]
fcdCoral_intercept <- fcm_labileCorals_lm$coefficients["(Intercept)"]


fcdCoral_r2 <- summary(fcm_labileCorals_lm)$adj.r.squared

#Day algae
fcm_labileAlgae_lm <- fcm_23%>%
  filter(DayNight == 'Day',
         Organism %like any% c('Dict%', 'Turf'))%>%
  lm(cells_ul ~ log10, data = .)

fcdAlgae_p <- (fcm_labileAlgae_lm%>% 
                 tidy()%>% 
                 filter(term == 'log10'))$p.value

# fcdAlgae_f <- (fcm_labileAlgae_lm%>% 
#                  tidy()%>% 
#                  filter(term == 'log10'))$statistic
# 
# fcdAlgae_slope <- fcm_labileAlgae_lm$coefficients["log10"]
# fcdAlgae_intercept <- fcm_labileAlgae_lm$coefficients["(Intercept)"]
# 
# 
# fcdAlgae_r2 <- summary(fcm_labileAlgae_lm)$adj.r.squared


# Night lm
fcm_labile_lm_night <- fcm_23%>%
  filter(DayNight == 'Night')%>%
  lm(cells_ul ~ log10, data = .)

fcn_p <- (fcm_labile_lm_night%>% 
            tidy()%>% 
            filter(term == 'log10'))$p.value

fcn_f <- (fcm_labile_lm_night%>% 
            tidy()%>% 
            filter(term == 'log10'))$statistic

fcn_slope <- fcm_labile_lm_night$coefficients["log10"]
fcn_intercept <- fcm_labile_lm_night$coefficients["(Intercept)"]


fcn_r2 <- summary(fcm_labile_lm_night)$adj.r.squared

#Day Corals
fcm_NlabileCorals_lm <- fcm_23%>%
  filter(DayNight == 'Night',
         Organism %like any% c('Porites%', 'Poci%'))%>%
  lm(cells_ul ~ log10, data = .)

fcdNCoral_p <- (fcm_NlabileCorals_lm%>% 
                 tidy()%>% 
                 filter(term == 'log10'))$p.value

fcdNCoral_f <- (fcm_NlabileCorals_lm%>% 
                 tidy()%>% 
                 filter(term == 'log10'))$statistic

fcdNCoral_slope <- fcm_NlabileCorals_lm$coefficients["log10"]
fcdNCoral_intercept <- fcm_NlabileCorals_lm$coefficients["(Intercept)"]


fcdNCoral_r2 <- summary(fcm_NlabileCorals_lm)$adj.r.squared

#night algae
fcm_NlabileAlgae_lm <- fcm_23%>%
  filter(DayNight == 'Night',
         Organism %like any% c('Dict%', 'Turf'))%>%
  lm(cells_ul ~ log10, data = .)

fcdNAlgae_p <- (fcm_NlabileAlgae_lm%>% 
                 tidy()%>% 
                 filter(term == 'log10'))$p.value
# 
# fcdNAlgae_f <- (fcm_NlabileAlgae_lm%>% 
#                  tidy()%>% 
#                  filter(term == 'log10'))$statistic
# 
# fcdNAlgae_slope <- fcm_NlabileAlgae_lm$coefficients["log10"]
# fcdNAlgae_intercept <- fcm_NlabileAlgae_lm$coefficients["(Intercept)"]
# 
# 
# fcdNAlgae_r2 <- summary(fcm_NlabileAlgae_lm)$adj.r.squared


# Plotting
# Day
fcm23DayCorals <- fcm_23%>%
  filter(DayNight == 'Day',
         Organism %like any% c('Porites%', 'Pocil%'),
         Organism != 'Influent',
         Organism != 'Water control',
         Organism != 'Offshore')

fcm23DayAlgae <- fcm_23%>%
  filter(DayNight == 'Day',
         Organism %like any% c('Dict%', 'Turf%'),
         Organism != 'Influent',
         Organism != 'Water control',
         Organism != 'Offshore')

fcm23NightCorals <- fcm_23%>%
  filter(DayNight == 'Night',
         Organism %like any% c('Porites%', 'Pocil%'),
         Organism != 'Influent',
         Organism != 'Water control',
         Organism != 'Offshore')

fcm23NightAlgae <- fcm_23%>%
  filter(DayNight == 'Night',
         Organism %like any% c('Dict%', 'Turf%'),
         Organism != 'Influent',
         Organism != 'Water control',
         Organism != 'Offshore')

pdf('plots/fcm_sumDepletion.pdf', width = 15, height = 10)
fcm_23%>%
  filter(DayNight == 'Day',
         Organism != 'Influent',
         Organism != 'Water control',
         Organism != 'Offshore')%>%
  ggplot(aes(log10, cells_ul)) +
  geom_point(aes(color = Organism), size = 7) +
  geom_smooth(data = fcm23DayAlgae, method = 'lm', se = FALSE, linetype = 3, color = org_colors_no_water[[5]]) +
  geom_smooth(data = fcm23DayCorals, method = 'lm', color = org_colors_no_water[[3]]) +
  ylim(0.028,0.079) +
  xlim(7,10.5) +
  scale_color_manual(values = org_colors_no_water) +
  gen_theme() +
  labs(y = 'Specific Growth Rate', x = bquote(Labile ~Total ~Depletion ~(log[10] ~XIC)))

#Night

fcm_23%>%
  filter(DayNight == 'Night',
         Organism != 'Influent',
         Organism != 'Water control',
         Organism != 'Offshore')%>%
  ggplot(aes(log10, cells_ul)) +
  geom_point(aes(color = Organism), size = 7) +
  geom_smooth(method = 'lm', color = 'grey') +
  geom_smooth(data = fcm23NightAlgae, method = 'lm', se = FALSE, linetype = 3, color = org_colors_no_water[[5]]) +
  geom_smooth(data = fcm23NightCorals, method = 'lm', color = org_colors_no_water[[3]]) +
  ylim(0.028,0.079) +
  xlim(7,10.5) +
  scale_color_manual(values = org_colors_no_water) +
  gen_theme() +
  labs(y = 'Specific Growth Rate', x = bquote(Labile ~Total ~Depletion ~(log[10] ~XIC)))

dev.off()


# STATS -- Specific Growth Rate reduction ---------------------------------
fcmAverages <- fcm_23%>%
  group_by(Organism, DayNight)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()

sgrAverages <- fcmAverages%>%
  select(Organism, DayNight, cells_ul)%>%
  pivot_wider(names_from = `DayNight`, values_from = `cells_ul`)%>%
  mutate(percentDiff = (1-Day/Night)*100,
         Day24 = Day*24,
         Night24 = Night*24)

loadAverages <- fcmAverages%>%
  select(Organism, DayNight, final_cells)%>%
  pivot_wider(names_from = `DayNight`, values_from = `final_cells`)%>%
  mutate(percentDiff = (1-Day/Night)*100)


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


sgrAnova <- dorc_fcm_fdom%>%
  select(-c(1:4, 9:36, 38, 39))%>%
  spread(Timepoint, `Cells µL-1`)%>%
  mutate(hour = 24,
         final_cells = T4)%>%
  group_by(Organism, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         cells_ul = (log(final_cells - T0))/(hour))%>%
  select(-c(T0:TF, hour))%>%
  left_join(metabolitePool_labileChange,  by = c('Organism', 'DayNight', 'Replicate'))%>%
  mutate(log10 = log10(abs(difference)))%>%
  filter(Organism != 'Influent',
         Organism != 'Offshore')%>%
  # filter(!Organism %like any% c('Influent', 'Offshore'))%>%
  mutate(Organism = as.factor(Organism),
         Organism = fct_relevel(Organism, 'Water control'))%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(anova = map(data, ~aov(cells_ul ~ Organism, data = .x)%>%
                      tidy()),
         dunnet = map(data, ~DunnettTest(cells_ul ~ Organism, data = .x)))
  # mutate(test = 'sgr')

microbialGrowth_fdr <- microbialLoad_stats%>%
  bind_rows(sgr_stats)%>%
  mutate(fdr = p.adjust(data, method = "BH"))


# Figure S2 -- Labile heatmap -----------------------------------
depletion_hc <- feature_stats_wdf%>%
  select(-c(log10:asin))%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  filter(activity == 'labile')%>%
  spread(Timepoint, xic)%>%
  group_by(activity, DayNight, Organism, net_act)%>%
  nest()%>%
  mutate(data = map(data, ~filter(.x, !is.na(T0))%>%
                      group_by(ecoNetConsensus, Replicate)%>%
                      summarize_if(is.numeric, sum, na.rm = TRUE)%>%
                      ungroup()
  ))%>%
  unnest(data)%>%
  ungroup()%>%
  unite(sample, c('DayNight', 'Organism', 'Replicate'), sep = '_')%>%
  select(-c(network:numberOfNodes, activity, TF))%>%
  mutate(ecoNetConsensus = case_when(is.na(ecoNetConsensus) ~ '1 No Consensus',
                                     TRUE ~ ecoNetConsensus))%>%
  unite(netEcoNet, c('ecoNetConsensus', 'net_act'), sep = '_')%>%
  filter(!netEcoNet %like% '%Lignan%',
         !netEcoNet %like% '%No Con%')

depletionDfCraig <- depletion_hc%>%
  separate(sample, c('DayNight', 'Organism', 'Replicate'), sep = '_')%>%
  separate(netEcoNet, c('ecoNetConsensus', 'net_act'), sep = '_')%>%
  separate(ecoNetConsensus, c('superclass', 'class', 'subclass'), sep =';')

numberOfUnique <- depletionDfCraig%>%
  select(DayNight, Organism, subclass, T0)%>%
  group_by(Organism, DayNight, subclass)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  mutate(count = 1)%>%
  group_by(subclass, DayNight)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  filter(count < 2)
  

depletion_levels <- depletion_hc$netEcoNet%>% 
  as.factor()

organism_levels <- depletion_hc$sample%>%
  as.factor()
# fct_relevel(c('CCA', 'Turf', 'Dictyota'))

# pdf('plots/labileHeatmap.pdf', width = 20, height = 15)
depletion_hc%>%
  group_by(netEcoNet)%>%
  mutate(T0 = zscore(T0))%>%
  ggplot(aes(sample, netEcoNet, fill = T0)) +
  geom_tile() +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  scale_y_discrete(limits = rev(levels(depletion_levels))) +
  scale_x_discrete(limits = levels(organism_levels)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15))
# dev.off()

pdf('plots/labileViolinNets.pdf', width = 20, height = 30)
depletion_hc%>%
  filter(!netEcoNet %like% '%lignan%')%>%
  separate(sample, c('DayNight', 'Organism', 'Replicate'), sep = '_')%>%
  mutate(T0 = log10(T0))%>%
  ggplot(aes(netEcoNet, T0)) +
  geom_violin() +
  geom_point(aes(color = Organism), size = 4) + 
  scale_color_manual(values = org_colors_no_water) +
  # scale_color_manual(values = c('#F09837', 'grey')) +
  scale_x_discrete(limits = rev(levels(depletion_levels))) +
  # facet_wrap(~Organism, nrow = 1) +
  facet_wrap(~DayNight) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 15),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
        panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray")) +
  labs(y = bquote(Log[10] ~(XIC)))
dev.off()



# pdf('plots/recalcitrantViolinNets.pdf', width = 20, height = 50)
# depletion_hc%>%
#   filter(!netEcoNet %like% '%lignan%')%>%
#   separate(sample, c('DayNight', 'Organism', 'Replicate'), sep = '_')%>%
#   mutate(T0 = log10(T0))%>%
#   ggplot(aes(netEcoNet, T0)) +
#   geom_violin() +
#   geom_point(aes(color = Organism), size = 4) + 
#   scale_color_manual(values = org_colors_no_water) +
#   # scale_color_manual(values = c('#F09837', 'grey')) +
#   scale_x_discrete(limits = rev(levels(depletion_levels))) +
#   # facet_wrap(~Organism, nrow = 1) +
#   facet_wrap(~DayNight) +
#   coord_flip() +
#   theme(axis.text.x = element_text(size = 10),
#         legend.text = element_text(size = 25),
#         legend.title = element_text(size = 25),
#         panel.background = element_rect(fill = "transparent"), # bg of the panel
#         plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#         panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
#         panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray")) +
#   labs(y = bquote(Log[10] ~(XIC)))
# dev.off()

write_csv(depletion_hc%>% select(netEcoNet)%>% separate(netEcoNet, c('ecoNetConsensus', 'net_act'), sep = '_', remove = FALSE)%>% separate(ecoNetConsensus, c('Superclass', 'Class', 'Subclass'), sep = ';')%>% unique(), 'analysis/labileViolinColumn.csv')

# Figure 2b -- superclass percent pies------------------------------------
percentSuperclass <- org_lability%>%
  filter(!is.na(superclass_consensus))%>%
  group_by(Organism, DayNight)%>%
  mutate(percent = nxic/sum(nxic, na.rm = TRUE)*100,
         percentSd = std/sum(nxic, na.rm = TRUE)*100)
  # filter(DayNight == 'Day')%>%
  # ungroup()
  # group_by(Organism)%>%
  # nest()
# superclassColors <- c('#332288', '#88CCEE', '#DDCC77', '#CC6677', '#117733', '#44AA99', '#CC0033', 'black', 'grey')

pieChartSetup <- percentSuperclass%>%
  # filter(!superclass_consensus %like any% c('Lignans%', 'Alkaloids%', 'Phenylpropanoids%'))%>%
  mutate(colors = case_when(superclass_consensus == 'Benzenoids' ~ '#332288',
                            superclass_consensus %like% 'Lipids%' ~ '#88CCEE',
                            superclass_consensus %like% '%acids%' ~ '#DDCC77',
                            superclass_consensus %like% '%acids%' ~ '#CC6677',
                            superclass_consensus %like% '%nitrogen%' ~ '#117733',
                            superclass_consensus %like% '%oxygen%' ~ '#44AA99',
                            superclass_consensus %like% '%heterocyc%' ~ '#CC0033',
                            superclass_consensus %like% '%lkaloids%' ~ '#FFCCFF',
                            superclass_consensus %like% '%Lign%' ~ '#CCFFCC',
                            superclass_consensus %like% '%Phenyl%' ~ '#CCFFFF'))

colorsSuperclass <- pieChartSetup$colors
names(colorsSuperclass) <- pieChartSetup$superclass_consensus

pieCharts <- pieChartSetup%>%
  unite(title, c('Organism', 'DayNight'), sep= '_')%>%
  group_by(title)%>%
  mutate(total = sum(nxic, na.rm = TRUE),
         title = as.factor(title),
         title = fct_relevel(title, c('CCA_Day', 'Turf_Day', 'Dictyota_Day', 'Pocillopora verrucosa_Day', 'Porites lobata_Day',
                                      'CCA_Night', 'Turf_Night', 'Dictyota_Night', 'Pocillopora verrucosa_Night', 'Porites lobata_Night')))%>%
  # nest()%>%
  # mutate(data = map2(data, title, ~ 
  filter(!is.na(percent))
  



scales <- pieCharts$xic
# pos <- 0.5 * (cumsum(scales) + cumsum(c(0, scales[-length(scales)])))

pdf('plots/OrgCoxbomb.pdf', width = 12, height = 10)
pieCharts%>%
  group_by(title)%>%
  mutate(sizes = 0.5 * (cumsum(xic) + cumsum(c(0, xic[-length(xic)]))))%>%
  ungroup()%>%
  ggplot(aes(x = sizes)) +
  geom_bar(aes(y = xic, fill = superclass_consensus), color = "white", stat = "identity", width = scales) +
  facet_wrap(~title, nrow = 2) +
  labs(x= NULL, y= NULL) +
  scale_fill_manual(values = colorsSuperclass) +
  theme_classic() +
  coord_polar(theta = 'x') +
  scale_x_continuous(aes(labels = pieCharts$superclass_consensus),breaks = scales ) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 8, face = "bold"))
dev.off()




pdf('plots/OrgPieCharts.pdf', width = 12, height = 10)
pieCharts%>%
  ungroup()%>%
  group_by(title)%>%
  arrange(percent)%>%
  mutate(width = cumsum(percent))%>%
  ungroup()%>%
  mutate(size = log10(total/sum(total)*1000))%>%
  ggplot(aes(x = size/2, width = 0.8,  y = percent, fill = superclass_consensus)) +
  geom_bar(color = "white", stat = "identity") +
  facet_wrap(~title, nrow = 2) +
  labs(x= NULL, y= NULL) +
  scale_fill_manual(values = colorsSuperclass) +
  coord_polar(theta = "y", start = 0) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
dev.off()


# FIGURE xxxx -- superclass Cluster dendogram --------------------------------------------
superclassDendo <- feature_stats_wdf%>%
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
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  ungroup()%>%
  filter(activity == 'labile')%>%
  select(-c(network, numberOfNodes, ecoNetConsensusScore, ra, asin, net_act, log10))%>%
  group_by(superclass_consensus, DayNight, Replicate, Organism)%>%
  summarize_if(is.numeric, sum, na.rm  = TRUE)%>%
  ungroup()%>%
  mutate(superclass_consensus = ifelse(is.na(superclass_consensus), 'unclassified', superclass_consensus))%>%
  filter(superclass_consensus != 'Alkaloids and derivatives',
         superclass_consensus != 'Lignans, neolignans and related compounds',
         superclass_consensus != 'Phenylpropanoids and polyketides',
         superclass_consensus != 'Organic oxygen compounds')%>%
  group_by(superclass_consensus)%>%
  mutate(xic = log10(xic),
         xic = zscore(xic))%>%
  ungroup()

minXICSuperclass <- (superclassDendo%>%
  filter(xic > 0,
         !is.na(xic)))$xic%>%
  min()

pdf('~/Documents/GitHub/DORCIERR/data/plots/superclassCluster.pdf', width = 15, height = 10)
superclassDendo%>%
  # group_by(superclass_consensus, Organism, DayNight)%>%
  # summarise_if(is.numeric, mean)%>%
  # ungroup()%>%
  # unite(sample, c(Organism, DayNight), sep = '  ')%>%
  unite(sample, c(Organism, DayNight, Replicate), sep = '  ')%>%
  spread(sample, xic)%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  gather(sample, xic, 2:ncol(.))%>%
  mutate(sample = as.factor(sample),
         # sample = fct_relevel(sample, c('CCA  Day', 'CCA  Night', 'Turf  Day', 'Turf  Night','Dictyota  Day',  'Dictyota  Night')))%>%
         sample = fct_relevel(sample, c('CCA  Day  1', 
                                        'CCA  Day  2', 
                                        'CCA  Night  1', 
                                        'CCA  Night  2', 
                                        'Turf  Day  1', 
                                        'Turf  Day  2', 
                                        'Dictyota  Day  1', 
                                        'Dictyota  Day  2', 
                                        'Turf  Night  1', 
                                        'Turf  Night  2', 
                                        'Dictyota  Night  1',
                                        'Dictyota  Night  2',
                                        'Pocillopora verrucosa  Day  1',
                                        'Pocillopora verrucosa  Day  2',
                                        'Porites lobata  Day  1',
                                        'Porites lobata  Day  2'
                                        )))%>%
  spread(superclass_consensus, xic)%>%
  column_to_rownames(var = 'sample')%>%
  pheatmap::pheatmap(color = brewer.pal(n = 9, name = "Greys"), cluster_rows = FALSE)
dev.off()


# VIZUALIZATIONS -- PCOA of Consensus annotations, superclass permanova -------------------------
# Have to load in EmojiFont here otherwise it messes with unicode points
minValue <- (feature_stats_wdf%>% 
  filter(log10 > 0))$log10%>%
  min()

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
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  filter(superclass_consensus != 'Alkaloids and derivatives',
         superclass_consensus != 'Lignans, neolignans and related compounds',
         superclass_consensus != 'Phenylpropanoids and polyketides')%>%
  ungroup()%>%
  filter(activity == 'labile')%>%
  # left_join(log2_change_vals%>%
  #             ungroup(), by = c('feature_number', 'Replicate', 'Organism', 'DayNight'))%>%
  select(-c(network, numberOfNodes, ecoNetConsensusScore, ra, asin, net_act))%>%
  group_by(subclass_consensus, DayNight, Replicate, Organism)%>%
  summarize_if(is.numeric, sum, na.rm  = TRUE)%>%
  ungroup()%>%
  mutate(log10 = log10(xic +1))%>%
  select(-xic)%>%
  # filter(!is.na(subclass_consensus))%>%
  unite(sample, c('DayNight', 'Organism', 'Replicate'), sep = '  ')%>%
  spread(sample, log10)%>%
  # mutate_all(~replace(., is.na(.), 0))%>%
  gather(sample, log10, 2:ncol(.))%>%
  spread(subclass_consensus, log10)%>%
  column_to_rownames(var = "sample")%>%
  vegdist(distance = 'bray', na.rm = TRUE)%>%
  pcoa()




## Plot Eigenvalues
pcoa$values[1:10,]%>%
  as.data.frame()%>%
  rownames_to_column("Axis")%>%
  mutate(axis = as.numeric(Axis))%>%
  ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, color = "red", vjust = -0.5)

pdf("plots/annotationPcoaBray.pdf", width = 15, height = 13)
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
    axis.title = element_text(size = 25),
    legend.position = 'None') +
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
  filter(activity == 'labile')%>%
  select(-c(network, numberOfNodes, ecoNetConsensusScore, log10, ra, asin, net_act))%>%
  separate(ecoNetConsensus, c('superclass', 'class', 'subclass'), sep = ';')%>%
  group_by(subclass, DayNight, Replicate, Organism)%>%
  summarize_if(is.numeric, sum, na.rm  = TRUE)%>%
  ungroup()%>%
  mutate(xic = log10(xic + 1))%>%
  spread(subclass, xic)%>%
  mutate_all(~replace(., is.na(.), 0))

# This can check subclass or superclass, just ahve to change the grouping in the above section
permanova_organism <- permanova_df%>%
  mutate(orgGroup = case_when(Organism %like% 'Pocill%' | Organism %like% 'Porites%' ~ 'Coral',
                              TRUE ~ 'Algae'))%>%
  select(DayNight, Replicate, Organism, orgGroup, everything())%>%
  group_by(DayNight)%>%
  nest()%>%
  # mutate(data = map(data, ~ adonis2(.x[4:ncol(.x)] ~ orgGroup, .x, perm = 100000, method = "bray", na.rm = TRUE)))
  mutate(data = map(data, ~adonis2(.x[4:ncol(.x)] ~ Organism, .x, perm = 100000, method = "bray", na.rm = TRUE)))

permanovaOverall <- permanova_df%>%
  adonis2(.[4:ncol(.)] ~ Organism, ., perm = 100000, method = "bray", na.rm = TRUE)

# pca ---------------------------------------------------------------------

pc <- feature_stats_wdf%>%
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
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  ungroup()%>%
  filter(activity == 'labile')%>%
  # left_join(log2_change_vals%>%
  #             ungroup(), by = c('feature_number', 'Replicate', 'Organism', 'DayNight'))%>%
  select(-c(network, numberOfNodes, ecoNetConsensusScore, ra, asin, net_act))%>%
  group_by(superclass_consensus, DayNight, Replicate, Organism)%>%
  summarize_if(is.numeric, sum, na.rm  = TRUE)%>%
  ungroup()%>%
  mutate(xic = log10(xic))%>%
  select(-xic)%>%
  filter(!is.na(superclass_consensus))%>%
  unite(sample, c('DayNight', 'Organism', 'Replicate'), sep = '  ')%>%
  spread(sample, log10)%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  gather(sample, log10, 2:ncol(.))%>%
  spread(superclass_consensus, log10)%>%
  separate(sample, c('DayNight', 'Organism', 'Replicate'), sep = '  ')%>%
  select(-Replicate)
  # mutate(DayNight = ifelse(DayNight == 'Day', '\u2600', emoji('crescent_moon')))
  
# column_to_rownames(var = "sample")%>%
pcView <- prcomp(pc[,-c(1:2)], center = TRUE,
             scale. = TRUE)

pdf('plots/pcaSubclasses.pdf', width = 15, height = 10)
ggbiplot::ggbiplot(pcView,
              obs.scale = 1,
              var.scale = 1,
              var.axes = FALSE,
              point.size = 5,
              groups = pc$Organism,
              ellipse = TRUE,
              # circle = TRUE,
              ellipse.prob = 0.68,
              labels = pc$DayNight,
              labels.size = 5,
              ) +
  scale_color_manual(values = org_colors_no_water)+ 
  scale_fill_manual(values = org_colors_no_water) +
  theme(legend.direction = 'horizontal',
               legend.position = 'top') +
        # text = element_text(family = 'OpenSansEmoji')) +
  gen_theme() 
dev.off()


# VIZUALIZATIONS -- pcoa class groupings ----------------------------------
pcoaClass <- feature_stats_wdf%>%
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
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  ungroup()%>%
  filter(activity == 'labile')%>%
  select(-c(network, numberOfNodes, ecoNetConsensusScore, ra, asin, net_act))%>%
  group_by(superclass_consensus, subclass_consensus, DayNight, Replicate, Organism)%>%
  summarize_if(is.numeric, sum, na.rm  = TRUE)%>%
  ungroup()%>%
  unite(sample, c('DayNight', 'Organism', 'Replicate'), sep = '  ')%>%
  mutate(log10 = log10(xic + 1))%>%
  select(-xic)%>%
  spread(sample, log10)%>%
  # mutate_all(~replace(., is.na(.), minValue))%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  gather(sample, log10, 3:ncol(.))%>%
  # group_by(superclass_consensus, subclass_consensus)%>%
  # filter(sum(log10) > 0)%>%
  # ungroup()%>%
  group_by(superclass_consensus)%>%
  filter(superclass_consensus != 'Alkaloids and derivatives',
         # superclass_consensus != 'Benzenoids',
         superclass_consensus != 'Lignans, neolignans and related compounds',
         # superclass_consensus != 'Organoheterocyclic compounds',
         superclass_consensus != 'Phenylpropanoids and polyketides',
  # filter(superclass_consensus == 'Organic oxygen compounds',
  #        !sample %like% 'Day  Turf  2')%>%
         superclass_consensus != 0)%>%
  # mutate(superclass_consensus = ifelse(superclass_consensus == 0, 'Unclassified', superclass_consensus))%>%
  nest()%>%
  mutate(data = map(data, ~  spread(.x, subclass_consensus, log10)%>%
                      column_to_rownames(var = "sample")%>%
                      # metaMDS(distance = 'manhattan')),
                      vegdist(method = 'euclidian', na.rm = TRUE)%>%
                      pcoa()),
         plots = map2(data, superclass_consensus, ~ .x$vectors%>%
                       as.data.frame()%>%
                       rownames_to_column(var = "sample")%>%
                       separate(sample, c('DayNight', 'Organism', 'Replicate'), sep = '  ')%>%
                       mutate(shape = case_when(DayNight == 'Day' ~ '\u2600',
                                                TRUE ~ emoji('crescent_moon')))%>%
                       ggplot(., aes(x = Axis.1, y = Axis.2, color = Organism)) +
                       geom_text(aes(label = shape), cex = 13, family= 'OpenSansEmoji') +
                       scale_color_manual(values = org_colors_no_water) +
                       labs(title = .y) +
                       theme(
                         panel.background = element_rect(fill = "transparent"), # bg of the panel
                         plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                         panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
                         panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
                         legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                         legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
                         legend.text = element_text(face = "italic"),
                         axis.text = element_text(size = 25),
                         axis.title = element_text(size = 25),
                         legend.position = 'None') +
                       xlab(str_c("Axis 1", " (", round(.x$values$Relative_eig[1], digits = 4)*100, "%)", sep = "")) +
                       ylab(str_c("Axis 2", " (", round(.x$values$Relative_eig[2], digits = 4)*100, "%)", sep = ""))))

## Plot Eigenvalues
# pcoaClass$values[1:10,]%>%
#   as.data.frame()%>%
#   rownames_to_column("Axis")%>%
#   mutate(axis = as.numeric(Axis))%>%
#   ggplot(aes(reorder(Axis, axis), Relative_eig, label = round(Relative_eig, digits = 3))) +
#   geom_bar(stat = "identity") +
#   geom_text(size = 3, color = "red", vjust = -0.5)

pdf("plots/OrganicOxygenPcoaBrayCurtis.pdf", width = 15, height = 13)
pcoaClass$plots
dev.off()


# Fig 5c -- class/subclass View --------------------------------------
netViews <- lability_classes%>%
  group_by(feature_number, Organism, DayNight)%>%
  mutate(meanProd = mean(T0, na.rm = TRUE),
         val = TF - meanProd)%>%
  left_join(affinity%>%
              select(feature_number:DayNight, betterAffinity), by = c('feature_number', 'Organism', 'Replicate', 'DayNight'))%>%
  group_by(feature_number, Replicate)%>%
  mutate(count = 1,
         count = sum(count))


eclipseDfDeplete <- netViews%>%
  mutate(meanProd = log10(meanProd + 1))%>%
  unite(ellipse, c('Organism', 'DayNight'), sep = '_', remove = FALSE) %>%
  unite(facet, c('network', 'class_consensus', 'subclass_consensus'), sep = '_', remove = FALSE)

depleteLevels <- c('CCA_T0', "CCA_TF", "Turf_T0", "Turf_TF", "Dictyota_T0", "Dictyota_TF", 
                   "Pocillopora verrucosa_T0", "Pocillopora verrucosa_TF", "Porites lobata_T0", "Porites lobata_TF")

depleteBarPlots <- eclipseDfDeplete%>%
  filter(activity == 'labile')%>%
  # unite(facet, c(facet, DayNight), sep = '_')%>%
  group_by(Organism, Replicate, DayNight, superclass_consensus, class_consensus, subclass_consensus,network, net_act)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  select(-c(val:count, ecoNetConsensusScore, numberOfNodes))%>%
  gather(Timepoint, xic, T0:TF)%>%
  filter(Timepoint == 'T0' & Replicate == 1 |
         Timepoint == 'T0' & Replicate == 2 |
           Timepoint == 'TF')%>%
  spread(Organism, xic)%>%
  gather(Organism, xic, CCA:Turf)%>%
  spread(DayNight, xic)%>%
  gather(DayNight, xic, Day:Night)%>%
  group_by(Organism, Timepoint, DayNight, superclass_consensus, class_consensus, subclass_consensus,network, net_act)%>%
  # mutate(countReps = 1)%>%
  mutate(n = 1,
         n = sum(n),
         xic = log10(xic + 1),
         std = sd(xic, na.rm = TRUE)/sqrt(n))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  unite(xAxis, c(Organism, Timepoint), sep = '_', remove = FALSE)%>%
  unite(lineDeplete, c(Organism, DayNight), sep = '_', remove = FALSE)%>%
  mutate(xAxis = as.factor(xAxis),
         xAxis = fct_relevel(xAxis, depleteLevels))


plotDepletes <- function(x) {
  ggplot(x, aes(xAxis, xic, fill = Organism)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    geom_errorbar(aes(ymin = xic - std, ymax = xic + std), position = 'dodge', color = 'black', width = 0.5) +
    geom_line(aes(group = lineDeplete), size = 2, color = 'black') +
    facet_wrap(~DayNight, ncol = 2) +
    # scale_x_discrete(limits = levels(depleteLevels)) +
    scale_fill_manual(values = org_colors_no_water) + 
    scale_color_manual(values = c('white', 'white')) +
    gen_theme() +
    labs(x = '', y = '') +
    theme(legend.position = 'None',
          axis.text.x = element_blank(),
          strip.text = element_blank()) }

betterOptions <- c('1841', '21', '198', '3009', '343', '141', '2248', '3765', '3762', 
                   '1045', '17', '659', '199', '3352', '280', '1969')

finalOptions <- c('343', '141', '3765', '3762',  '199',  '280', '1969', '3585')

recalcitrantOptions <- c(
  '5297'
  # '1699', 
  # '2766', 
  # '3032', 
  # '-2510', 
  # '5673'
  )

depleteFacettedBarPlots <- depleteBarPlots%>%
  # filter(net_act %in% finalOptions)%>%
  # filter(net_act == 3585,
  #        DayNight == 'Night')%>%
  filter(net_act %in% recalcitrantOptions,
         DayNight == 'Day')%>%
  group_by(facet)%>%
  nest()%>%
  mutate(data = map(data, ~ plotDepletes(.x) +
                      ggtitle(facet) +
                      scale_y_continuous(limits = c(6,7), oob = rescale_none)))


# recalcitrantFacettedBarPlots <- lmerBarPlots%>%
#   filter(net_act %in% recalcitrantOptions,
#          DayNight == 'Night')%>%
#   group_by(facet)%>%
#   nest()%>%
#   mutate(data = map(data, ~ plotLMER(.x) +
#                       ggtitle(facet)))

pdf('plots/recalcitrantFinals.pdf', widt = 7.5, height = 10)
depleteFacettedBarPlots[[2]]
dev.off()

pdf('~/Documents/GitHub/DORCIERR/data/plots/Fig5FinalOptions.pdf', width = 7.5, height = 4)
depleteBarPlots%>%
  filter(network == '21')%>%
  plotDepletes() +
  ggtitle('Fatty acid esters') +
  scale_y_continuous(limits = c(1, 8.5), oob = rescale_none)

depleteBarPlots%>%
  filter(network == '343')%>%
  plotDepletes() +
  ggtitle('Glycerophosphoethanolamines') +
  scale_y_continuous(limits = c(3.6, 8.2), oob = rescale_none)
  

depleteBarPlots%>%
  filter(network == '4032')%>%
  plotDepletes() +
  ggtitle('Carboxylic acids') +
  scale_y_continuous(limits = c(6.3, 9.3), oob = rescale_none)

depleteBarPlots%>%
  filter(network == '3765')%>%
  plotDepletes() +
  ggtitle('Tricarboxylic acids and derivatives') +
  scale_y_continuous(limits = c(7.7, 8.55), oob = rescale_none)

depleteBarPlots%>%
  filter(network == '198')%>%
  plotDepletes() +
  ggtitle('Lineolic acids') +
  scale_y_continuous(limits = c(6.3, 8.1), oob = rescale_none)


#cluster 2
depleteBarPlots%>%
  filter(network == '2263')%>%
  plotDepletes() +
  ggtitle('Fatty amides') +
  scale_y_continuous(limits = c(7.1, 7.7), oob = rescale_none)

depleteBarPlots%>%
  filter(network == '103')%>%
  plotDepletes() +
  ggtitle('Fatty acyls') +
  scale_y_continuous(limits = c(5.5, 6.4), oob = rescale_none)

depleteBarPlots%>%
  filter(network == '5425')%>%
  plotDepletes() +
  ggtitle('Amino acids and peptides') +
  scale_y_continuous(limits = c(6.2, 8.4), oob = rescale_none)

depleteBarPlots%>%
  filter(network == '91')%>%
  plotDepletes() +
  ggtitle('Diterpenoids') +
  scale_y_continuous(limits = c(6.3, 9.7), oob = rescale_none)

#Cluster 3
depleteBarPlots%>%
  filter(network == '2447')%>%
  plotDepletes() +
  ggtitle('Prenol lipids') +
  scale_y_continuous(limits = c(6.1, 6.4), oob = rescale_none)

depleteBarPlots%>%
  filter(network == '280')%>%
  plotDepletes() +
  ggtitle('N-alkylindole') +
  scale_y_continuous(limits = c(6.9, 8.9), oob = rescale_none)

depleteBarPlots%>%
  filter(network == '225')%>%
  plotDepletes() +
  ggtitle('Amines') +
  scale_y_continuous(limits = c(5, 7.4), oob = rescale_none)

depleteBarPlots%>%
  filter(net_act == '-6536')%>%
  plotDepletes() +
  ggtitle('1,3,5 - triazines') +
  scale_y_continuous(limits = c(5.5, 7.6), oob = rescale_none)
dev.off()


# Figure 1c -- lmer example -----------------------------------------------
lmerExamples <- lability_classes%>%
  select(-c(net_act, ecoNetConsensus:matchSource, val))%>%
  unique()%>%
  gather(Timepoint, xic, T0:TF)%>%
  filter(Timepoint == 'T0' & Replicate == 1 |
           Timepoint == 'T0' & Replicate == 2 |
           Timepoint == 'TF')%>%
  pivot_wider(names_from = 'Organism', values_from = 'xic')%>%
  pivot_longer(7:ncol(.), names_to = 'Organism', values_to = 'xic')%>%
  group_by(Organism, Timepoint, DayNight, network, feature_number)%>%
  mutate(n = 1,
         n = sum(n),
         xic = log10(xic + 1),
         std = sd(xic, na.rm = TRUE)/sqrt(n))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  # unite(xAxis, c(Organism, Timepoint), sep = '_', remove = FALSE)%>%
  unite(lineDeplete, c(Organism, DayNight, feature_number), sep = '_', remove = FALSE)%>%
  mutate(Organism = as.factor(Organism),
         Organism = fct_relevel(Organism, depleteLevels))


plotF1c <- function(x) {
  ggplot(x, aes(xAxis, xic, color = Organism)) +
    geom_point(stat = 'identity') +
    # geom_bar(stat = 'identity', position = 'dodge') +
    geom_errorbar(aes(ymin = xic - std, ymax = xic + std), position = 'dodge', color = 'black') +
    geom_line(aes(group = lineDeplete), size = 2, color = 'black') +
    facet_wrap(~DayNight, ncol = 2) +
    # scale_x_discrete(limits = levels(depleteLevels)) +
    scale_color_manual(values = org_colors_no_water) +
    # scale_color_manual(values = c('white', 'white')) +
    gen_theme() +
    labs(x = '', y = '') +
    theme(legend.position = 'None',
          axis.text.x = element_blank(),
          strip.text = element_blank()) }


lmerExamplePlots <- lmerExamples%>%
  filter(network == 3585 & DayNight == 'Night' |
           network == 141 & DayNight == 'Night'|
           network == 5673 & DayNight == 'Day'|
           network == 5297 & DayNight == 'Day')%>%
  mutate(facet = network)%>%
  group_by(network)%>%
  nest()%>%
  mutate(data = map(data, ~  ggplot(.x, aes(Timepoint, xic, color = Organism)) +
                      # geom_bar(stat = 'identity', position = 'dodge') +
                      geom_errorbar(aes(ymin = xic - std, ymax = xic + std), position = 'dodge', color = 'black', width = 0.2) +
                      geom_point(stat = 'identity', size = 5) +
                      geom_line(aes(group = lineDeplete, color = Organism), size = 2) +
                      # ggtitle(.x$facet) +
                      # facet_wrap(~DayNight, ncol = 2) +
                      # scale_x_discrete(limits = levels(depleteLevels)) +
                      scale_color_manual(values = org_colors_no_water) +
                      gen_theme() +
                      labs(x = '', y = '') +
                      theme(legend.position = 'None',
                            axis.text.x = element_blank(),
                            strip.text = element_blank())))
  

pdf('~/Documents/GitHub/DORCIERR/data/plots/LmerExamples.pdf', width = 15, height = 10)
lmerExamplePlots[[2]]
dev.off()


# Figure 6 -- Superclass lability % change--------------------------------
labileChange <- feature_stats_wdf%>%
  ungroup()%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  filter(activity == 'labile')%>%
  # Timepoint == 'T0')%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  select(Organism, Replicate, Timepoint, DayNight, superclass_consensus, xic)%>%
  filter(Organism != 'CCA')%>%
  unique()%>%
  group_by(Organism, Replicate, Timepoint, DayNight)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'))%>%
  group_by(group, DayNight, Timepoint)%>%
  mutate(n = 1,
         n = sum(n),
         xic = log10(xic + 1),
         sterr = sd(xic)/sqrt(n))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  unite(colorFill, c('group', 'DayNight'), sep = '_', remove = FALSE)


#Depletion
labileDepletion <- metabolitePool_xic_change%>%
  select(-c(err:ra))%>%
  group_by(org_activity, Organism, activity, DayNight)%>%
  mutate(x_val = max(xic, na.rm = TRUE))%>%
  ungroup()%>%
  group_by(org_activity, Organism, Replicate, activity, Timepoint, DayNight)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  spread(Timepoint, xic)%>%
  group_by(activity, Organism, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         difference = T0-TF)%>%
  ungroup()%>%
  filter(activity == 'labile')%>%
  unite(group, c('Organism', 'DayNight'), sep = '_', remove = FALSE)%>%
  filter(Organism != 'CCA')%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'))%>%
  group_by(group, DayNight)%>%
  mutate(T0 = mean(T0, na.rm = TRUE),
         # depletion = T0-TF,
         n = 1,
         n = sum(n),
         depletion = ifelse(difference < 0, 0, difference),
         propDepletion = depletion/T0,
         depletion = log10(depletion + 1),
         propSterr = sd(propDepletion)/sqrt(n),
         sterr = sd(depletion)/sqrt(n))%>%
  # filter(depletion > 0)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  unite(colorFill, c('group', 'DayNight'), sep = '_', remove = FALSE)


# labileDepletion <- feature_stats_wdf%>%
#   ungroup()%>%
#   left_join(networking%>%
#               select(feature_number, network), by = 'feature_number')%>%
#   mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
#                              TRUE ~ network))%>%
#   left_join(all_activity, by = c('net_act', 'DayNight'))%>%
#   filter(activity == 'labile')%>%
#   # Timepoint == 'T0')%>%
#   select(feature_number, Organism, Replicate, Timepoint, DayNight, xic)%>%
#   # unique()%>%
#   filter(Organism != 'CCA')%>%
#   mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'))%>%
#   unique()%>%
#   group_by(Organism, group, Replicate, Timepoint, DayNight)%>%
#   summarize_if(is.numeric, sum)%>%
#   ungroup()%>%
#   spread(Timepoint, xic)%>%
#   group_by(group, DayNight)%>%
#   mutate(T0 = mean(T0, na.rm = TRUE),
#          depletion = T0-TF,
#          n = 1,
#          n = sum(n),
#          depletion = ifelse(depletion < 0, 0, depletion),
#          propDepletion = depletion/T0,
#          depletion = log10(depletion + 1),
#          propSterr = sd(propDepletion)/sqrt(n),
#          sterr = sd(depletion)/sqrt(n))%>%
#   filter(depletion > 0)%>%
#   summarize_if(is.numeric, mean, na.rm = TRUE)%>%
#   ungroup()%>%
#   unite(colorFill, c('group', 'DayNight'), sep = '_', remove = FALSE)

#Peak Load and SGR data frame
fcmGroup <- fcm_23%>%
  filter(Organism != 'CCA')%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'))%>%
  group_by(group, DayNight)%>%
  mutate(n = 1,
         n = sum(n),
         loadSterr = sd(final_cells)/sqrt(n),
         sgrSterr = sd(cells_ul)/sqrt(n))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  unite(colorFill, c('group', 'DayNight'), sep = '_', remove = FALSE)

#Superclasses
superclassChange <- feature_stats_wdf%>%
  ungroup()%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  filter(activity == 'labile')%>%
  # Timepoint == 'T0')%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  select(Organism, Replicate, Timepoint, DayNight, superclass_consensus, xic)%>%
  filter(Organism != 'CCA')%>%
  unique()%>%
  group_by(Organism, Replicate, Timepoint, DayNight, superclass_consensus)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  mutate(group = ifelse(Organism %like% c('Poc%', 'Por%'), 'Coral', 'Fleshy algae'))%>%
  group_by(group, DayNight, Timepoint, superclass_consensus)%>%
  mutate(n = 1,
         n = sum(n),
         xic = log10(xic + 1),
         sterr = sd(xic)/sqrt(n))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  unite(colorFill, c('group', 'DayNight'), sep = '_', remove = FALSE)%>%
  filter(superclass_consensus %like any% c('%nitrogen%' , 'Lipids%', '%heterocy%'))%>%
  group_by(superclass_consensus)%>%
  filter(Timepoint == 'T0')

figSixFills <- c('#EEB481', '#EF8C36', '#8DCA7D', '#67C94D')
figSixColors <- c('#EF8C36', '#C86C1D', '#67C94D', '#398A23')

# Making all of the plots
pdf('~/Documents/GitHub/DORCIERR/data/plots/figSixplots.pdf', width = 15, height = 10)
labileChange%>%
  filter(Timepoint == 'T0')%>%
  ggplot(aes(DayNight, xic, color = colorFill, fill = colorFill)) +
  geom_bar(stat = 'identity', size = 4) +
  geom_errorbar(aes(ymin = xic -sterr, ymax = xic + sterr), color = 'black', width = 0, size = 2) +
  geom_point(aes(y = xic +sterr), color = 'black', size = 4) +
  geom_point(aes(y = xic - sterr), color = 'black', size = 4) +
  geom_line(aes(group = group), color = 'black', size = 2, linetype = 'dashed') +
  facet_wrap(~group) +
  scale_color_manual(values = figSixColors) +
  scale_fill_manual(values = figSixFills) +
  gen_theme() +
  scale_y_continuous(limits = c(9.7, 10.4), oob = rescale_none) +
  labs(x = 'Diel cycle', y = 'Intensity (log10 XIC)') +
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank())

labileDepletion%>%
  ggplot(aes(DayNight, depletion, color = colorFill, fill = colorFill)) +
  geom_bar(stat = 'identity', size = 4) +
  geom_errorbar(aes(ymin = depletion -sterr, ymax = depletion + sterr), color = 'black', width = 0, size = 2) +
  geom_point(aes(y = depletion +sterr), color = 'black', size = 4) +
  geom_point(aes(y = depletion - sterr), color = 'black', size = 4) +
  geom_line(aes(group = group), color = 'black', size = 2, linetype = 'dashed') +
  facet_wrap(~group) +
  scale_color_manual(values = figSixColors) +
  scale_fill_manual(values = figSixFills) +
  gen_theme() +
  scale_y_continuous(limits = c(6.8, 10), oob = rescale_none) +
  labs(x = 'Diel cycle', y = 'Intensity (log10 XIC)')+
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank())

labileDepletion%>%
  ggplot(aes(DayNight, propDepletion, color = colorFill, fill = colorFill)) +
  geom_bar(stat = 'identity', size = 4) +
  geom_errorbar(aes(ymin = propDepletion -propSterr, ymax = propDepletion +propSterr), color = 'black', width = 0, size = 2) +
  geom_point(aes(y = propDepletion +propSterr), color = 'black', size = 4) +
  geom_point(aes(y = propDepletion -propSterr), color = 'black', size = 4) +
  geom_line(aes(group = group), color = 'black', size = 2, linetype = 'dashed') +
  facet_wrap(~group) +
  scale_color_manual(values = figSixColors) +
  scale_fill_manual(values = figSixFills) +
  gen_theme() +
  # scale_y_continuous(limits = c(6.9, 10), oob = rescale_none) +
  labs(x = 'Diel cycle', y = 'Intensity (log10 XIC)')+
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank())

fcmGroup%>%
  ggplot(aes(DayNight, final_cells, color = colorFill, fill = colorFill)) +
  geom_bar(stat = 'identity', size = 4) +
  geom_errorbar(aes(ymin = final_cells - loadSterr, ymax = final_cells + loadSterr), color = 'black', width = 0, size = 2) +
  geom_point(aes(y = final_cells +loadSterr), color = 'black', size = 4) +
  geom_point(aes(y = final_cells - loadSterr), color = 'black', size = 4) +
  geom_line(aes(group = group), color = 'black', size = 2, linetype = 'dashed') +
  facet_wrap(~group) +
  scale_color_manual(values = figSixColors) +
  scale_fill_manual(values = figSixFills) +
  gen_theme() +
  scale_y_continuous(limits = c(550, 820), oob = rescale_none) +
  labs(x = 'Diel cycle', y = 'Microbial Load (Cells ml-1)')+
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank())

fcmGroup%>%
  ggplot(aes(DayNight, cells_ul, color = colorFill, fill = colorFill)) +
  geom_bar(stat = 'identity', size = 4) +
  geom_errorbar(aes(ymin = cells_ul - sgrSterr, ymax = cells_ul + sgrSterr), color = 'black', width = 0, size = 2) +
  geom_point(aes(y = cells_ul +sgrSterr), color = 'black', size = 4) +
  geom_point(aes(y = cells_ul - sgrSterr), color = 'black', size = 4) +
  geom_line(aes(group = group), color = 'black', size = 2, linetype = 'dashed') +
  facet_wrap(~group) +
  scale_color_manual(values = figSixColors) +
  scale_fill_manual(values = figSixFills) +
  gen_theme() +
  scale_y_continuous(limits = c(0.05, 0.072), oob = rescale_none) +
  labs(x = 'Diel cycle', y = 'Specific Growth Rate (day-1)')+
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank())

superclassChange%>%
  filter(superclass_consensus %like% 'Lipids%')%>%
  ggplot(aes(DayNight, xic, color = colorFill, fill = colorFill)) +
  geom_bar(stat = 'identity', size = 4) +
  geom_errorbar(aes(ymin = xic -sterr, ymax = xic + sterr), color = 'black', width = 0, size = 2) +
  geom_point(aes(y = xic +sterr), color = 'black', size = 4) +
  geom_point(aes(y = xic - sterr), color = 'black', size = 4) +
  geom_line(aes(group = group), color = 'black', size = 2, linetype = 'dashed') +
  facet_wrap(~group) +
  scale_color_manual(values = figSixColors) +
  scale_fill_manual(values = figSixFills) +
  gen_theme() +
  scale_y_continuous(limits = c(8.5, 9.6), oob = rescale_none) +
  labs(x = 'Diel cycle', y = 'Lipids Intensity (log10 XIC)')+
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank())

superclassChange%>%
  filter(superclass_consensus %like% '%heter%')%>%
  ggplot(aes(DayNight, xic, color = colorFill, fill = colorFill)) +
  geom_bar(stat = 'identity', size = 4) +
  geom_errorbar(aes(ymin = xic -sterr, ymax = xic + sterr), color = 'black', width = 0, size = 2) +
  geom_point(aes(y = xic +sterr), color = 'black', size = 4) +
  geom_point(aes(y = xic - sterr), color = 'black', size = 4) +
  geom_line(aes(group = group), color = 'black', size = 2, linetype = 'dashed') +
  facet_wrap(~group) +
  scale_color_manual(values = figSixColors) +
  scale_fill_manual(values = figSixFills) +
  gen_theme() +
  scale_y_continuous(limits = c(6.8, 9), oob = rescale_none) +
  labs(x = 'Diel cycle', y = 'Organoheterocyclic Intensity (log10 XIC)')+
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank())

superclassChange%>%
  filter(superclass_consensus %like% '%nitrogen%')%>%
  ggplot(aes(DayNight, xic, color = colorFill, fill = colorFill)) +
  geom_bar(stat = 'identity', size = 4) +
  geom_errorbar(aes(ymin = xic -sterr, ymax = xic + sterr), color = 'black', width = 0, size = 2) +
  geom_point(aes(y = xic +sterr), color = 'black', size = 4) +
  geom_point(aes(y = xic - sterr), color = 'black', size = 4) +
  geom_line(aes(group = group), color = 'black', size = 2, linetype = 'dashed') +
  facet_wrap(~group) +
  scale_color_manual(values = figSixColors) +
  scale_fill_manual(values = figSixFills) +
  gen_theme() +
  scale_y_continuous(limits = c(6.8, 8.7), oob = rescale_none) +
  labs(x = 'Diel cycle', y = 'Nitrogen Intensity (log10 XIC)')+
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank())
dev.off()

# VIZUALIZATIONS -- Stoichiometry and labile networks ---------------------
stoich <- netViews%>%
  left_join(networking%>%
              select(feature_number, network, C:S),
            by = c('feature_number', 'network'))%>%
  # filter(activity == 'labile')%>%
  # filter(network %in% finalOptions)%>%
  group_by(Organism, Replicate, DayNight, activity)%>%
  mutate(nc = N/C,
         oc = O/C,
         pc = P/C,
         log10 = log10(T0 + 1),
         wnc = nc*T0/sum(T0),
         wpc = pc*T0/sum(T0))%>%
  filter(nc < 0.8)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  unite(ellipse, c('Organism', 'DayNight'), remove = FALSE, sep = '_')

stoich%>%
  filter(activity == 'labile')%>%
  ggplot(aes(wnc, wpc, color = Organism)) +
  geom_text(aes(label = DayNight), stat = 'identity') +
  ggalt::geom_encircle(aes(fill = Organism, color = Organism, group = Organism), alpha = 0.3, s_shape = 0.4) +
  # geom_text(aes(label = network)) +
  scale_color_manual(values = org_colors_no_water) +
  scale_fill_manual(values = org_colors_no_water) +
  gen_theme()
  # facet_wrap(~activity)


# pdf('plots/exampleLabileBarNetworksLipids.pdf', width = 30, height = 120)
# depleteBarPlots%>%
#   filter(superclass_consensus %like% c('%lipid%'))%>%
#   plotDepletes()
# dev.off()
# 
# pdf('plots/exampleLabileBarNetworksAcids.pdf', width = 30, height = 70)
# depleteBarPlots%>%
#   filter(superclass_consensus %like% c('%acids%'))%>%
#   plotDepletes()
# dev.off()
# 
# pdf('plots/exampleLabileBarNetworksAllElse.pdf', width = 30, height = 20)
# depleteBarPlots%>%
#   filter(superclass_consensus %like% c('%nitrogen%'))%>%
#   plotDepletes()
# 
# depleteBarPlots%>%
#   filter(superclass_consensus %like% c('%eterocyclic%'))%>%
#   plotDepletes()
# 
# depleteBarPlots%>%
#   filter(superclass_consensus %like% c('%Benzenoids%'))%>%
#   plotDepletes()
# 
# depleteBarPlots%>%
#   filter(superclass_consensus %like% c('%oxygen%'))%>%
#   plotDepletes()
# dev.off()

# pdf('plots/exampleLabileNetworks.pdf', width = 20, height = 30)
# eclipseDfDeplete%>%
#   ggplot(aes(meanProd, val)) +
#   ggalt::geom_encircle(aes(fill = Organism, color = Organism, group = ellipse), alpha = 0.3, s_shape = 0.2, expand = 0.1) +
#   geom_text(data = eclipseMeanDeplete, mapping = aes(label = shape, color = Organism), stat = 'identity', cex = 13, family= 'OpenSansEmoji') +
#   scale_fill_manual(values = org_colors_no_water) +
#   scale_color_manual(values = org_colors_no_water) +
#   # scale_shape_manual(values = c('\u2600', emoji('crescent_moon'))) +
#   facet_wrap(~facet, ncol = 2, scales = 'free_y') +
#   labs(x = 'T0 average production (XIC)', y = 'Utilization') +
#   gen_theme() +
#   theme(legend.position = 'None')
# dev.off()
# 
# pdf('plots/exampleRecalcitrantNetworks.pdf', width = 20, height = 20)
# eclipseDfRecalcitrant%>%
#   ggplot(aes(meanProd, val)) +
#   ggalt::geom_encircle(aes(fill = Organism, color = Organism, group = ellipse), alpha = 0.3, s_shape = 0.5) +
#   geom_text(data = eclipseMeanDeplete, mapping = aes(label = shape, color = Organism), stat = 'identity', cex = 13, family= 'OpenSansEmoji') +
#   scale_fill_manual(values = org_colors_no_water) +
#   scale_color_manual(values = org_colors_no_water) +
#   scale_shape_manual(values = c(16, 1)) +
#   facet_wrap(~facet, ncol = 2, scales = 'free_y') +
#   labs(x = 'T0 average production (XIC)', y = 'Utilization') +
#   gen_theme() +
#   theme(legend.position = 'None')
# dev.off()
# 
# names <- netViews%>%
#   ungroup()%>%
#   filter(network %in% c('199', '437', '1969', '1072', '21',
#                         '324', '2261',  '1907', '1758', '2023',
#                         '165', '814', '247', '147', '1633', '673',
#                         '1337', '3032', '2023', '2064'))%>%
#   select(network, ecoNetConsensus:matchSource)%>% 
#   unique()
# 
# libNets <- true_hits%>%
#   mutate(feature_number = as.character(feature_number))%>%
#   left_join(networking%>% 
#               select(feature_number, network), by = 'feature_number')%>%
#   filter(network %in% c('199', '437', '1969', '1072', '21',
#                         '324', '2261',  '1907', '1758', '2023',
#                         '165', '814', '247', '147', '1633', '673',
#                         '1337', '3032', '2023', '2064'))%>%
#   filter(!is.na(Compound_Name))%>%
#   select(network, Compound_Name)%>%
#   unique()

# VIZUALIZATIONS -- lability super, sub and class dataframes --------------------
classesIncluded <- c('Quinolines%', 'Lactones','Carboxylic%', 'Prenol%', 'Fatty%', 'Carbonyl%')
subclassesIncluded <- c('Carbohydrates%', 'Amines', 'Carbonyl%')


# Class lability subplots
class_lability <- lability_classes%>%
  filter(class_consensus %like any% classesIncluded,
         !is.na(superclass_consensus))%>%
  spread(activity, T0)%>%
  gather(activity, T0, labile:recalcitrant)%>%
  group_by(activity, DayNight, superclass_consensus, class_consensus, Replicate)%>%
  filter(!Replicate %in% c('3', '4'))%>%
  # filter(!is.na(val))%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, class_consensus, Replicate)%>%
  mutate(relativeT0Prod = T0/sum(T0, na.rm = TRUE),
         relativeT0Prod = case_when(is.na(relativeT0Prod) ~ 0,
                                    TRUE ~ relativeT0Prod),
         relativeT0Prod = relativeT0Prod*100)%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, class_consensus, activity)%>%
  mutate(std = sd(relativeT0Prod))%>%
  ungroup()%>%
  mutate(shape = case_when(DayNight == 'Day' ~ '\u2600',
                           TRUE ~ emoji('crescent_moon')))%>%
  mutate(class_consensus = case_when(is.na(class_consensus) ~ 'None',
                                     TRUE ~ class_consensus))%>%
  unite(facet, c('superclass_consensus', 'DayNight'), sep = '_', remove = FALSE)%>%
  unite(yAxis, c('superclass_consensus', 'class_consensus'), sep = '_', remove = FALSE)

#Subclass subplots
subclass_lability <- lability_classes%>%
  filter(subclass_consensus %like any% subclassesIncluded,
         !is.na(superclass_consensus))%>%
  spread(activity, T0)%>%
  gather(activity, T0, labile:recalcitrant)%>%
  group_by(activity, DayNight, superclass_consensus, class_consensus, subclass_consensus, Replicate)%>%
  # filter(!is.na(val))%>%
  filter(!Replicate %in% c('3', '4'))%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, class_consensus,subclass_consensus, Replicate)%>%
  mutate(relativeT0Prod = T0/sum(T0, na.rm = TRUE),
         relativeT0Prod = case_when(is.na(relativeT0Prod) ~ 0,
                                    TRUE ~ relativeT0Prod),
         relativeT0Prod = relativeT0Prod*100)%>%
  ungroup()%>%
  filter(!Replicate %in% c('3', '4'))%>%
  group_by(DayNight, superclass_consensus, class_consensus, subclass_consensus, activity)%>%
  mutate(std = sd(relativeT0Prod))%>%
  # filter(sum(T0) != 0)%>%
  ungroup()%>%
  mutate(shape = case_when(DayNight == 'Day' ~ '\u2600',
                           TRUE ~ emoji('crescent_moon')))%>%
  mutate(subclass_consensus = case_when(is.na(subclass_consensus) ~ 'None',
                                        TRUE ~ subclass_consensus))%>%
  unite(facet, c('superclass_consensus', 'DayNight'), sep = '_', remove = FALSE)%>%
  unite(yAxis, c('superclass_consensus', 'subclass_consensus'), sep = '_', remove = FALSE)

subclassClass_lability <- class_lability%>%
  bind_rows(subclass_lability)

# VIZUALIZATIONS -- plotting superclass lability --------------------------
pdf('plots/emojiSuperclassPercentProductionDayNight.pdf', width = 15, height = 10)
lability_classes%>%
  group_by(activity, DayNight, superclass_consensus, Replicate)%>%
  filter(!is.na(val),
         !superclass_consensus %like% 'Lignans%')%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, Replicate)%>%
  mutate(relativeT0Prod = T0/sum(T0, na.rm = TRUE),
         relativeT0Prod = case_when(is.na(relativeT0Prod) ~ 0,
                                    TRUE ~ relativeT0Prod))%>%
  ungroup()%>%
  filter(!Replicate %in% c('3', '4'))%>%
  group_by(DayNight, superclass_consensus, activity)%>%
  filter(sum(T0) != 0)%>%
  ungroup()%>%
  mutate(shape = case_when(DayNight == 'Day' ~ '\u2600',
                           TRUE ~ emoji('crescent_moon')))%>%
  ggplot(aes(superclass_consensus, relativeT0Prod*100, color = activity, shape = shape)) +
  geom_text(aes(label = shape), cex = 13, family= 'OpenSansEmoji') +
  coord_flip() +
  scale_color_manual(labels = c('Labile', 'Recalcitrant'), values = c('#EBCC2A', "#006658")) +
  # geom_bar(stat = 'identity') +
  # geom_errorbar(aes(ymax = total_production + std, ymin = total_production - std)) +
  # facet_wrap(~activity, labeller = labeller(activity = label_wrap_gen()), scales = 'free_y', nrow = 3) +
  labs(y = 'Percent production (%)', x = 'Putative Superclass Annotation', color = 'Network Reactivity') +
  facet_wrap(~DayNight) +
  # scale_y_log10() +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        axis.text.y = element_text(size = 14))
dev.off()


pdf('plots/SuperclassPercentBars.pdf', width = 15, height = 10)
lability_classes%>%
  group_by(activity, DayNight, superclass_consensus, Replicate)%>%
  filter(!is.na(val),
         !superclass_consensus %like any% c('Lignans%', 'Phenylprop%', 'Alkaloids%', 'Nucleosides%'))%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, Replicate)%>%
  mutate(relativeT0Prod = T0/sum(T0, na.rm = TRUE),
         relativeT0Prod = case_when(is.na(relativeT0Prod) ~ 0,
                                    TRUE ~ relativeT0Prod),
         relativeT0Prod = relativeT0Prod*100)%>%
  ungroup()%>%
  filter(!Replicate %in% c('3', '4'))%>%
  group_by(DayNight, superclass_consensus, activity)%>%
  # filter(sum(T0) != 0)%>%
  mutate(std = sd(relativeT0Prod))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  ggplot(aes(superclass_consensus, relativeT0Prod, fill = activity)) +
  geom_errorbar(aes(ymax = relativeT0Prod + std, ymin = relativeT0Prod - std), position = 'dodge') +
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  scale_fill_manual(labels = c('Labile', 'Recalcitrant'), values = c('#EBCC2A', "#006658")) +
  # geom_bar(stat = 'identity') +
  # facet_wrap(~activity, labeller = labeller(activity = label_wrap_gen()), scales = 'free_y', nrow = 3) +
  labs(y = 'Proportion of XIC produced', x = 'Putative Superclass Annotation', color = 'Network Reactivity') +
  facet_wrap(~DayNight) +
  # scale_y_log10() +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        axis.text.y = element_text(size = 14))

lability_classes%>%
  group_by(activity, DayNight, superclass_consensus, Replicate)%>%
  filter(!is.na(val),
         !superclass_consensus %like any% c('Lignans%', 'Phenylprop%', 'Alkaloids%', 'Nucleosides%'))%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, Replicate)%>%
  mutate(relativeT0Prod = T0/sum(T0, na.rm = TRUE),
         relativeT0Prod = case_when(is.na(relativeT0Prod) ~ 0,
                                    TRUE ~ relativeT0Prod),
         relativeT0Prod = relativeT0Prod*100)%>%
  ungroup()%>%
  filter(!Replicate %in% c('3', '4'))%>%
  group_by(DayNight, superclass_consensus, activity)%>%
  # filter(sum(T0) != 0)%>%
  mutate(std = sd(relativeT0Prod),
         std = case_when(activity == 'labile' ~ NA_real_,
                         TRUE ~ std))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  ggplot(aes(DayNight, relativeT0Prod, fill = activity))+
  geom_bar(stat = 'identity', position = 'stack') +
  # geom_errorbar(aes(ymax = relativeT0Prod + std, ymin = relativeT0Prod - std)) +
  # coord_flip() +
  scale_fill_manual(labels = c('Labile', 'Recalcitrant'), values = c('#EBCC2A', "#006658")) +
  # scale_color_manual(values = c('#F09837', 'grey')) +
  labs(y = 'Percent production (%)', x = 'Putative Superclass Annotation', color = 'Network Reactivity') +
  facet_wrap(~superclass_consensus, nrow = 1) +
  # scale_y_log10() +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        axis.text.y = element_text(size = 14))

lability_classes%>%
  group_by(activity, DayNight, superclass_consensus, Replicate)%>%
  filter(!is.na(val),
         !superclass_consensus %like any% c('Lignans%', 'Phenylprop%', 'Alkaloids%', 'Nucleosides%'))%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, Replicate)%>%
  mutate(relativeT0Prod = T0/sum(T0, na.rm = TRUE),
         relativeT0Prod = case_when(is.na(relativeT0Prod) ~ 0,
                                    TRUE ~ relativeT0Prod),
         relativeT0Prod = relativeT0Prod*100)%>%
  ungroup()%>%
  filter(!Replicate %in% c('3', '4'))%>%
  group_by(DayNight, superclass_consensus, activity)%>%
  # filter(sum(T0) != 0)%>%
  mutate(std = sd(relativeT0Prod),
         std = case_when(activity == 'labile' ~ NA_real_,
                         TRUE ~ std))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  ggplot(aes(superclass_consensus, relativeT0Prod, fill = activity))+
  geom_bar(stat = 'identity', position = 'stack') +
  # geom_errorbar(aes(ymax = relativeT0Prod + std, ymin = relativeT0Prod - std)) +
  # coord_flip() +
  scale_fill_manual(labels = c('Labile', 'Recalcitrant'), values = c('#EBCC2A', "#006658")) +
  # scale_color_manual(values = c('#F09837', 'grey')) +
  labs(y = 'Percent production (%)', x = 'Putative Superclass Annotation', color = 'Network Reactivity') +
  facet_wrap(~DayNight, nrow = 2) +
  # scale_y_log10() +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        axis.text.y = element_text(size = 14))

lability_classes%>%
  group_by(activity, DayNight, superclass_consensus, Replicate)%>%
  filter(!is.na(val),
         !superclass_consensus %like any% c('Lignans%', 'Phenylprop%', 'Alkaloids%', 'Nucleosides%'))%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, Replicate)%>%
  mutate(relativeT0Prod = T0/sum(T0, na.rm = TRUE),
         relativeT0Prod = case_when(is.na(relativeT0Prod) ~ 0,
                                    TRUE ~ relativeT0Prod),
         relativeT0Prod = relativeT0Prod*100)%>%
  ungroup()%>%
  filter(!Replicate %in% c('3', '4'))%>%
  group_by(DayNight, superclass_consensus, activity)%>%
  # filter(sum(T0) != 0)%>%
  mutate(std = sd(relativeT0Prod),
         std = case_when(activity == 'labile' ~ NA_real_,
                         TRUE ~ std))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  ggplot(aes(superclass_consensus, relativeT0Prod, fill = activity))+
  geom_bar(stat = 'identity', position = 'stack') +
  # geom_errorbar(aes(ymax = relativeT0Prod + std, ymin = relativeT0Prod - std)) +
  coord_flip() +
  scale_fill_manual(labels = c('Labile', 'Recalcitrant'), values = c('#EBCC2A', "#006658")) +
  # scale_color_manual(values = c('#F09837', 'grey')) +
  labs(y = 'Percent production (%)', x = 'Putative Superclass Annotation', color = 'Network Reactivity') +
  facet_wrap(~DayNight, nrow = 1) +
  # scale_y_log10() +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        axis.text.y = element_text(size = 14))
dev.off()


# VIZUALIZATIONS -- Plotting subclassClass lability -----------------------
pdf('plots/class_consensus.pdf', width = 15, height = 10)
subclassClass_lability%>%
  # filter(superclass_consensus %like any% c('Lipid%', '%hetero%', '%acids%'))%>%
  ggplot(aes(yAxis, relativeT0Prod, color = activity, shape = shape)) +
  geom_text(aes(label = shape), cex = 13, family= 'OpenSansEmoji') +
  facet_wrap(~DayNight) +
  coord_flip() +
  scale_color_manual(labels = c('Labile', 'Recalcitrant'), values = c('#EBCC2A', "#006658")) +
  # geom_bar(stat = 'identity') +
  # geom_errorbar(aes(ymax = total_production + std, ymin = total_production - std)) +
  # facet_wrap(~activity, labeller = labeller(activity = label_wrap_gen()), scales = 'free_y', nrow = 3) +
  labs(y = 'Percent production (%)', x = 'Putative Annotation', color = 'Network Reactivity') +
  # scale_y_log10() +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 5),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 5),
        legend.position = "None")

subclassClass_lability%>%
  ggplot(aes(yAxis, relativeT0Prod, fill = activity, shape = shape)) +
  geom_errorbar(aes(ymax = relativeT0Prod + std, ymin = relativeT0Prod - std), position = 'dodge') +
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_flip() +
  scale_fill_manual(labels = c('Labile', 'Recalcitrant'), values = c('#EBCC2A', "#006658")) +
  # geom_bar(stat = 'identity') +
  # facet_wrap(~activity, labeller = labeller(activity = label_wrap_gen()), scales = 'free_y', nrow = 3) +
  labs(y = 'Proportion of XIC produced', x = 'Putative Annotation', color = 'Network Reactivity') +
  facet_wrap(~DayNight) +
  scale_y_continuous(limits = c(0,100)) +
  # scale_y_log10() +
  gen_theme() +
  theme(strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13),
        axis.text.y = element_text(size = 14)) 
dev.off()

# VIZUALIZATIONS -- SPG depletion, together AND LINEAR MODEL -------------------------------
fcm_labile_lm_all <- fcm_23%>%
  lm(cells_ul ~ log10, data = .)

fca_p <- (fcm_labile_lm_all%>% 
            tidy()%>% 
            filter(term == 'log10'))$p.value

fca_f <- (fcm_labile_lm_all%>% 
            tidy()%>% 
            filter(term == 'log10'))$statistic

fca_slope <- fcm_labile_lm_all$coefficients["log10"]
fca_intercept <- fcm_labile_lm_all$coefficients["(Intercept)"]


fca_r2 <- summary(fcm_labile_lm_all)$adj.r.squared

pdf('plots/spgLabileChange.pdf', width = 15, height = 10)
fcm_23%>%
  group_by(Organism, DayNight)%>%
  mutate(xerr = plotrix::std.error(log10),
         yerr = plotrix::std.error(cells_ul))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  filter(Organism != 'Influent',
         Organism != 'Water control',
         Organism != 'Offshore')%>%
  mutate(shape = case_when(DayNight == 'Day' ~ '\u2600',
                           TRUE ~ emoji('crescent_moon')))%>%
  ggplot(aes(log10, cells_ul, color = Organism)) +
  geom_line(aes(group = Organism, color = Organism), size = 2) +
  geom_errorbarh(aes(xmin = log10-xerr, xmax = log10+xerr), size = 1, color = 'Grey') +
  geom_errorbar(aes(ymin = cells_ul - yerr, ymax = cells_ul + yerr), size = 1, color = 'Grey') +
  geom_text(aes(label = shape, color = Organism), cex = 13, family= 'OpenSansEmoji') +
  ylim(0.036,0.083) +
  scale_color_manual(values = org_colors_no_water) +
  gen_theme() +
  labs(y = 'Specific Growth Rate', x = bquote(Labile ~Total ~Depletion ~(log[10] ~XIC)))

fcm_23%>%
  group_by(Organism, DayNight)%>%
  mutate(xerr = sd(log10),
         yerr = sd(cells_ul))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  filter(Organism != 'Influent',
         Organism != 'Water control',
         Organism != 'Offshore')%>%
  mutate(shape = case_when(DayNight == 'Day' ~ '\u2600',
                           TRUE ~ emoji('crescent_moon')))%>%
  ggplot(aes(Organism, cells_ul, color = Organism)) +
  geom_line(aes(group = Organism, color = Organism), size = 2) +
  geom_errorbar(aes(ymin = cells_ul - yerr, ymax = cells_ul + yerr), size = 1, color = 'Grey', width = 0.5) +
  geom_text(aes(label = shape, color = Organism), cex = 13, family= 'OpenSansEmoji') +
  ylim(0.036,0.083) +
  scale_color_manual(values = org_colors_no_water) +
  gen_theme() +
  labs(y = 'Specific Growth Rate', x = 'Organism')

fcm_23%>%
  group_by(Organism, DayNight)%>%
  mutate(xerr = sd(log10))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  filter(Organism != 'Influent',
         Organism != 'Offshore')%>%
  mutate(shape = case_when(DayNight == 'Day' ~ '\u2600',
                           TRUE ~ emoji('crescent_moon')))%>%
  ggplot(aes(Organism, log10, color = Organism)) +
  geom_line(aes(group = Organism, color = Organism), size = 2) +
  geom_errorbar(aes(ymin = log10 - xerr, ymax = log10 + xerr), size = 1, color = 'Grey', width = 0.5) +
  geom_text(aes(label = shape, color = Organism), cex = 12, family= 'OpenSansEmoji') +
  scale_color_manual(values = org_colors_no_water) +
  gen_theme() +
  labs(y = bquote(Labile ~Utilization ~(log[10] ~XIC)), x = 'Organism')
dev.off()

prcntIncreaseSPG <-fcm_23%>%
  filter(Organism != 'CCA')%>%
  select(Organism, DayNight, Replicate, cells_ul)%>%
  spread(DayNight, cells_ul)%>%
  mutate(change = (Night-Day)/Day,
         orgType = case_when(Organism == 'Pocillopora verrucosa' ~ 'Coral',
                             Organism == 'Porites lobata' ~ 'Coral',
                             TRUE ~ 'Fleshyalgae'))%>%
  group_by(orgType)%>%
  mutate(sd = sd(change))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()

prcntIncreaseSPGOverall <-fcm_23%>%
  filter(Organism != 'CCA')%>%
  select(Organism, DayNight, Replicate, cells_ul)%>%
  spread(DayNight, cells_ul)%>%
  mutate(change = (Night-Day)/Day,
         orgType = case_when(Organism == 'Pocillopora verrucosa' ~ 'Coral',
                             Organism == 'Porites lobata' ~ 'Coral',
                             TRUE ~ 'Fleshyalgae'))%>%
  ungroup()%>%
  mutate(sd = sd(change))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()

prcntIncreaseMicLoad <-fcm_23%>%
  filter(Organism != 'CCA')%>%
  select(Organism, DayNight, Replicate, final_cells)%>%
  spread(DayNight, final_cells)%>%
  mutate(change = (Night-Day)/Day,
         orgType = case_when(Organism == 'Pocillopora verrucosa' ~ 'Coral',
                             Organism == 'Porites lobata' ~ 'Coral',
                             TRUE ~ 'Fleshyalgae'))%>%
  group_by(orgType)%>%
  mutate(sd = sd(change))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()

prcntIncreaseMicLoadOverall <-fcm_23%>%
  filter(Organism != 'CCA')%>%
  select(Organism, DayNight, Replicate, final_cells)%>%
  spread(DayNight, final_cells)%>%
  mutate(change = (Night-Day)/Day,
         orgType = case_when(Organism == 'Pocillopora verrucosa' ~ 'Coral',
                             Organism == 'Porites lobata' ~ 'Coral',
                             TRUE ~ 'Fleshyalgae'))%>%
  # group_by(orgType)%>%
  ungroup()%>%
  mutate(sd = sd(change))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()


# # VIZUALIZATIONS -- Cytoscape file ----------------------------------------
night_activity <- all_activity%>%
  filter(DayNight == 'Night')

day_activity <- all_activity%>%
  filter(DayNight == 'Day')

cytoDf<- log2_change_vals%>%
  inner_join(dom_stats_wdf%>%
               ungroup()%>%
               select(feature_number, Organism, DayNight)%>%
               unique(), 
             by = c('feature_number', 'Organism', 'DayNight'))%>%
  left_join(networking%>%
              select(feature_number, network),
            by = 'feature_number')%>%
  left_join(metadata%>%
              select(feature_number, `row m/z`), by = "feature_number")%>%
  filter(T0 != 0)%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  left_join(ecoNet%>% 
              rename(feature_number = 1)%>% 
              mutate(feature_number = as.character(feature_number)), 
            by = c('feature_number', 'network'))%>%
  select(feature_number, Organism, DayNight, net_act, network, ecoNetConsensus, log2_change)

cytoNets <- (cytoDf%>%
  select(network, feature_number)%>%
  unique()%>%
  filter(network != '-1')%>%
  group_by(network)%>%
  mutate(numberOfNodes = 1,
         numberOfNodes = sum(numberOfNodes))%>%
  ungroup()%>%
  select(network, numberOfNodes)%>%
  filter(numberOfNodes > 5))$network%>%
  unique()%>%
  as.vector()

# cytoExport <- cytoDf%>%
#   filter(network %in% cytoNets)%>%
#   group_by(feature_number, network, DayNight, ecoNetConsensus)%>%
#   summarize_if(is.numeric, mean, na.rm = TRUE)%>%
#   ungroup()%>%
#   group_by(feature_number, network, ecoNetConsensus)%>%
#   mutate(mostNeg = min(log2_change, na.rm = TRUE))%>%
#   spread(DayNight, log2_change)

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


Key_cyto_export <- feature_stats_wdf%>%
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
  select(feature_number, Organism, DayNight, net_act, ecoNetConsensus, xic)%>%
  group_by(feature_number, Organism, DayNight, ecoNetConsensus, net_act)%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  filter(net_act %in% finalOptions)%>%
  spread(DayNight, xic)%>%
  mutate(production = case_when(net_act %in% c('1969', '280', '199') ~ Day,
                                 TRUE ~ Night))%>%
  select(-c(Day, Night))%>%
  spread(Organism, production)%>%
  mutate(sum = rowSums(across(CCA:Turf), na.rm = TRUE))%>%
  left_join(night_activity, by = 'net_act')%>%
  left_join(day_activity, by = 'net_act', suffix = c('_night', '_day'))

write_csv(Key_cyto_export, 'analysis/KeycytoDiel.csv')

# 
# VIZUALIZATIONS -- SUPPLEMENT lability classes variables -----------------
lability_classes_supplements <- feature_stats_wdf%>%
  select(-c(log10, asin))%>%
  gather(responseVariable, val, xic:ra)%>%
  left_join(networking%>%
              select(feature_number, network), by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity, by = c('net_act', 'DayNight'))%>%
  filter(activity != 'accumolite')%>%
  # filter(DayNight == 'Day')%>%
  left_join(ecoNet%>%
              rename(feature_number = scan)%>%
              mutate(feature_number = as.character(feature_number)), by = c('feature_number', 'network'))%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  group_by(activity, DayNight, responseVariable)%>%
  nest()%>%
  mutate(data = map(data, ~ spread(.x, Timepoint, val)%>%
                      group_by(superclass_consensus, Replicate)%>%
                      summarize_if(is.numeric, sum, na.rm = TRUE)%>%
                      ungroup()))%>%
  unnest(data)%>%
  ungroup()%>%
  group_by(DayNight, superclass_consensus, Replicate, responseVariable)%>%
  mutate(relativeT0Prod = T0/sum(T0, na.rm = TRUE),
         relativeT0Prod = case_when(is.na(relativeT0Prod) ~ 0,
                                    TRUE ~ relativeT0Prod))%>%
  ungroup()%>%
  group_by(responseVariable)%>%
  nest()%>%
  mutate(plot = map(data, ~ 
                      ggplot(.x, aes(superclass_consensus, relativeT0Prod*100, color = activity)) +
                      geom_point(stat = 'identity', size = 4) +
                      coord_flip() +
                      ggtitle(responseVariable) +
                      scale_color_manual(labels = c('Labile', 'Recalcitrant'), values = c('#EBCC2A', "#006658")) +
                      # geom_bar(stat = 'identity') +
                      # geom_errorbar(aes(ymax = total_production + std, ymin = total_production - std)) +
                      # facet_wrap(~activity, labeller = labeller(activity = label_wrap_gen()), scales = 'free_y', nrow = 3) +
                      labs(y = 'Percent production (%)', x = 'Putative Superclass Annotation', color = 'Network Reactivity') +
                      facet_wrap(~DayNight) +
                      # scale_y_log10() +
                      gen_theme() +
                      theme(strip.background = element_rect(fill = "transparent"),
                            strip.text = element_text(size = 13),
                            axis.text.y = element_text(size = 14))),
         plotRaw = map(data, ~
                         ggplot(.x, aes(superclass_consensus, T0, color = activity)) +
                         geom_point(stat = 'identity', size = 4) +
                         coord_flip() +
                         ggtitle(responseVariable) +
                         scale_color_manual(labels = c('Labile', 'Recalcitrant'), values = c('#EBCC2A', "#006658")) +
                         # geom_bar(stat = 'identity') +
                         # geom_errorbar(aes(ymax = total_production + std, ymin = total_production - std)) +
                         # facet_wrap(~activity, labeller = labeller(activity = label_wrap_gen()), scales = 'free_y', nrow = 3) +
                         labs(y = 'raw value', x = 'Putative Superclass Annotation', color = 'Network Reactivity') +
                         facet_wrap(~DayNight) +
                         # scale_y_log10() +
                         gen_theme() +
                         theme(strip.background = element_rect(fill = "transparent"),
                               strip.text = element_text(size = 13),
                               axis.text.y = element_text(size = 14))))
# unnest(data)


pdf('plots/supplementalSuperclassPercentProductionDayNight.pdf', width = 15, height = 10)
lability_classes_supplements$plot
lability_classes_supplements$plotRaw
dev.off()




DictyotaAcids <- labileAcids%>%
  filter(Organism == 'Dictyota')%>%
  left_join(metadata%>%
              select(feature_number, C:NOSC), by = 'feature_number')

write_csv(DictyotaAcids, '~/Downloads/DictyotaAcids.csv')

pdf('~/Downloads/DictyotaAcids.pdf', width = 15, height = 10)
DictyotaAcids%>%
  filter(!Replicate %in% c(3,4))%>%
  mutate(network = as.character(network))%>%
  group_by(network, class_consensus, subclass_consensus, Replicate)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  group_by(network, class_consensus, subclass_consensus)%>%
  mutate(std = sd(T0, na.rm = TRUE))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  unite(consensusDisplay, c('class_consensus', 'subclass_consensus'), sep = '   ')%>%
  ggplot(aes(reorder(network, -T0), T0, fill = consensusDisplay)) +
  geom_errorbar(aes(ymin = T0 -std, ymax = T0 + std)) +
  geom_bar(stat = 'identity') +
  labs(y = 'Network Production (XIC)', x = 'Network') +
  gen_theme() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.text.x = element_text(size = 15))
dev.off()

eclipse199 <- netViews%>%
  mutate(meanProd = log10(meanProd + 1))%>%
  filter(network %in% c('199', '437', '627'),
         Organism == 'Dictyota'
         # filter(network %in% c('199', '437', '1969', '1072', '21',
         #                       '324', '2261',  '1907', '1758', '2023',
         #                       '280', '966')
         )%>%
  unite(ellipse, c('Organism', 'DayNight'), sep = '_', remove = FALSE) %>%
  unite(facet, c('network', 'activity', 'DayNight'), sep = '_', remove = FALSE)%>%
  group_by(network, activity, DayNight)%>%
  mutate(shape = case_when(DayNight == 'Day' ~ '\u2600',
                           TRUE ~ emoji('crescent_moon')),
         lowProp = case_when(meanProd < 1e6 ~ 1,
                             TRUE ~ 0),
         lowProp = lowProp/sum(lowProp),
         std = sd(meanProd))


eclipseMean199 <- eclipse199%>%
  group_by(ellipse, Organism, DayNight, facet)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  mutate(shape = case_when(DayNight == 'Day' ~ '\u2600',
                           TRUE ~ emoji('crescent_moon')))

pdf('~/Downloads/exampleLabileNetworks.pdf', width = 15, height = 10)
eclipse199%>%
  ggplot(aes(meanProd, val)) +
  geom_text(aes(label = shape, color = Organism), stat = 'identity', cex = 8, family = 'OpenSansEmoji') +
  # ggalt::geom_encircle(aes(fill = Organism, color = Organism, group = ellipse), alpha = 0.3, s_shape = 0.2, expand = 0.1) +
  # geom_text(data = eclipseMean199, mapping = aes(label = shape, color = Organism), stat = 'identity', cex = 13, family= 'OpenSansEmoji') +
  scale_fill_manual(values = org_colors_no_water) +
  scale_color_manual(values = org_colors_no_water) +
  # scale_shape_manual(values = c('\u2600', emoji('crescent_moon'))) +
  facet_wrap(~facet, ncol = 2, scales = 'free_y') +
  labs(x = 'T0 average production (XIC)', y = 'Utilization') +
  gen_theme() +
  theme(legend.position = 'None')
dev.off()





# VISUALIZATIONS -- conCISE plots -----------------------------------------
plotConcise <- function(x) {
  ggplot(x, aes(Organism, XIC, fill = Organism)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    geom_errorbar(aes(ymin = XIC - std, ymax = XIC + std), position = 'dodge', color = 'black') +
    geom_line(aes(group = lineDeplete), size = 2, color = 'black') +
    scale_fill_manual(values = org_colors_no_water) + 
    scale_color_manual(values = c('white', 'white')) +
    gen_theme() +
    labs(x = 'Organism', y = 'Intensity (XIC)') +
    theme(legend.position = 'None') }

conCISEFacettedBarPlots <- netViews%>%
  mutate(meanProd = log10(meanProd + 1))%>%
  unite(ellipse, c('Organism', 'DayNight'), sep = '_', remove = FALSE) %>%
  unite(facet, c('network', 'class_consensus', 'subclass_consensus'), sep = '_', remove = FALSE)%>%
  group_by(Organism, Replicate, DayNight, net_act, facet)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  select(-c(val:count, network, ecoNetConsensusScore, numberOfNodes))%>%
  gather(Timepoint, XIC, T0:TF)%>%
  filter(Timepoint == 'T0' & Replicate == 1 |
           Timepoint == 'T0' & Replicate == 2 |
           Timepoint == 'TF')%>%
  spread(Organism, XIC)%>%
  gather(Organism, XIC, CCA:Turf)%>%
  spread(DayNight, XIC)%>%
  gather(DayNight, XIC, Day:Night)%>%
  group_by(Organism, Timepoint, DayNight, net_act, facet)%>%
  mutate(std = sd(XIC, na.rm = TRUE))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  unite(xAxis, c(Organism, Timepoint), sep = '_', remove = FALSE)%>%
  unite(lineDeplete, c(Organism, DayNight), sep = '_', remove = FALSE)%>%
  filter(net_act %in% c('247', '2263', '131'),
         !xAxis %like% '%TF',
         DayNight == 'Day' | DayNight == 'Night' & net_act == '2263')%>%
  group_by(facet)%>%
  nest()%>%
  mutate(data = map(data, ~ plotConcise(.x) +
                      ggtitle(facet)))

pdf('~/Documents/SDSU_Scripps/ConCISE/investigations/Subnetworks/Dorc/production.pdf', width = 15, height = 10)
conCISEFacettedBarPlots[[2]]
dev.off()

conciseStackedBar <- netViews%>%
  mutate(meanProd = log10(meanProd + 1))%>%
  unite(ellipse, c('Organism', 'DayNight'), sep = '_', remove = FALSE) %>%
  unite(facet, c('network', 'class_consensus', 'subclass_consensus'), sep = '_', remove = FALSE)%>%
  group_by(Organism, Replicate, DayNight, net_act, facet, superclass_consensus)%>%
  summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  ungroup()%>%
  select(-c(val:count, network, ecoNetConsensusScore, numberOfNodes))%>%
  gather(Timepoint, XIC, T0:TF)%>%
  filter(Timepoint == 'T0' & Replicate == 1 |
           Timepoint == 'T0' & Replicate == 2 |
           Timepoint == 'TF')%>%
  spread(Organism, XIC)%>%
  gather(Organism, XIC, CCA:Turf)%>%
  spread(DayNight, XIC)%>%
  gather(DayNight, XIC, Day:Night)%>%
  group_by(Organism, Timepoint, DayNight, net_act, facet, superclass_consensus)%>%
  mutate(std = sd(XIC, na.rm = TRUE))%>%
  summarize_if(is.numeric, mean, na.rm = TRUE)%>%
  ungroup()%>%
  unite(xAxis, c(Organism, Timepoint), sep = '_', remove = FALSE)%>%
  unite(lineDeplete, c(Organism, DayNight), sep = '_', remove = FALSE)%>%
  filter(!xAxis %like% '%TF',
         DayNight == 'Day',
         # net_act != '199',
         !is.na(superclass_consensus),
         !is.na(XIC))%>%
  group_by(Organism)%>%
  mutate(percent = XIC/sum(XIC, na.rm = TRUE),
         superclassFig = case_when(superclass_consensus %like% 'Alka%' ~ 'Other',
                                   superclass_consensus %like% 'Phenyl%' ~ 'Other',
                                   superclass_consensus %like% 'Lignans%' ~ 'Other',
                                   TRUE ~ superclass_consensus))

superclassColors <- c('#332288', '#88CCEE', '#DDCC77', '#CC6677', '#117733', '#44AA99', '#CC0033', 'black', 'grey')



pdf('~/Documents/SDSU_Scripps/ConCISE/investigations/Subnetworks/Dorc/productionSTacked.pdf', width = 15, height = 10)
conciseStackedBar%>%
  ggplot(aes(Organism, percent, fill = superclassFig)) +
  geom_bar(stat = 'summary', fun.y = 'sum', position = 'stack') +
  gen_theme() +
  scale_fill_manual(values = superclassColors, na.translate = TRUE, na.value = 'grey') +
  labs(y = 'Percent of classified xic produced', x = 'Organism')
dev.off()



checkNet <- function(x) {
  net <- netViews%>% 
    filter(network == x)
  
  source <- net$matchSource%>% 
    unique()
  
  level <- net$ecoNetConsensusLevel%>%
    unique()
  
  num <- net$numberOfNodes%>%
    unique()
  
  class <- net$ecoNetConsensus%>%
    unique()
  
  score <- net$ecoNetConsensusScore%>%
    unique()
  
  paste0(x, ' is ', class, ' from ', source, '  ', level, '  ', num, ' nodes  ', score)
  
}


