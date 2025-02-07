# Microbial community metabolism of coral reef exometabolomes broadens the chemodiversity of labile dissolved organic matter
[![DOI](https://zenodo.org/badge/226174249.svg)](https://zenodo.org/badge/latestdoi/226174249) <br />
Zachary A. Quinlan, Craig E. Nelson, Irina Koester, Daniel Petras, Louis-Felix Nothias, Jacqueline Comstock, Brandie M. White, Lihini I. Aluwihare, Barbara A. Bailey, Craig A. Carlson, Pieter C. Dorrestein, Andreas F. Haas and *Linda Wegley Kelly
*Corresponding author

### Abstract:
Dissolved organic matter (DOM) comprises diverse compounds with variable bioavailability across aquatic ecosystems. The sources and quantities of DOM can influence microbial growth and community structure with effects on biogeochemical processes. To investigate the chemodiversity of labile DOM in tropical reef waters we tracked microbial utilization of over 3000 untargeted mass spectrometry ion features exuded from two coral and three algal species. Roughly half of these features clustered into over 500 biologically labile spectral subnetworks annotated to diverse structural superclasses, including benzenoids, lipids, organic acids, heterocyclics and phenylpropanoids, comprising on average one third of the ion richness and abundance within each chemical class. Distinct subsets of these labile compounds were exuded by algae and corals during the day and night, driving differential microbial growth and substrate utilization. This study expands the chemical diversity of labile marine DOM with implications for carbon cycling in coastal environments.

### Funding sources:
This research was sponsored by the US National Science Foundation Chemical Oceanography program, in collaboration with the Moorea Coral Reef Long-Term Ecological Research project (MCR-LTER), and the US National Science Foundation Graduate Research Fellowship Program.



## Repository Overview:
In this repository you will find the code used to analyze this data and all raw data <br /> <br />
We have tried to clearly index each script and file below but if you have any questions please email me: zquinlan@gmail.com <br /> <br />
If you're interested in using this data for otehr analyses, please reach out! There are so many questions that could be answered with the data and not enough time to answer them all...<br /><br />
All raw data is in the subdirectory data/raw/

*****
# Contents:
## code - all the code that we used to analyzed the data in the above manuscript
#### NewDOM.R
  - This code includes all the analyses we used to process the untargeted metabolomics data and flow cytometry data
#### MicrobialPipeline2021.R
  - This code includes all the analyses we used to process the 16S data

<br /><br />
## data - raw data, plots, and analyzed spreadsheets we exported during analysis
* not all of this was used in the publication <br />
### analysis:
- In here there are many exports from statistical analyses and some curated datasets although the majority that are used, remain in the scripts themselves.<br />

### plots:
- Plots and figures generated in R<br />

### raw -- Some additional dataframes are not identified in the list below. They were not helpful to us but we are leaving them in here in case others could utilize them
### Metabolomics 
#### sirius4_06012021
- All exports from SIRIUS 4. This includes canopus summaries and molecular formula <br />
#### Analog-Hits.csv
- analog matches from GNPS
#### Canopus_classes.csv
- CANOPUS classifications using SIRIUS 2 -- This is out of date as we have the SIRIUS 4 exports
#### Extracton_efficiency.csv
- This is the raw data of extraction efficiency of the PPL extracts measured directly on the LC MSMS
#### Library-Hits.csv
- library matches from GNPS
#### Mo'orea 2017 Mass spec sample codes - Sheet1.csv
- Metadata for each mass spec sample name
#### Morrea_feayures0Table_all_gap-filled5.csv
- Raw quantification table of feature ion intensity within each sample
#### Node_info.tsv
- Subnetwork numbers and other ancillary data such as m/z for each feature from GNPS
#### SIRIUS_Zodiac_converted.csv
- All sirius/zodiac molecular formula -- This is out of data as we have the sirius 4 exports
#### ecoNetConsensus.csv
- The export from [ConCISE](github.com/zquinlan/concise) using classyfire annotations from CANOPUS and GNPS library matches <br /><br /><br />
### Microbes
#### 052022_abundance_table_100.tsv
- 100% identity (ASV) table for the 16S data
#### 052022_abundance_table_97.tsv
- 97% identity (OTU) table for the 16S data
#### otu_repr_100.tre1.tsv
- distance matricies for ASV data
#### otu_repr_97.tre1.tsv
- distance matricies for OTU data

