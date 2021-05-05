var store = [{
        "title": "Contents",
        "excerpt":"Contents of Github Repository: DORCIERR README.md code canopus_model.R bicluster.py Microbial_pipeline.R metabolite_modeled_lability.R 12182019DORCIERR_working_pipeline.R metabolite_trends.R current_dorcierr_working_pipeline.R data test.csv mlr_df.csv analysis hc_depletolites_df.csv canopus_model_rough.txt multiple_reg_coefficients_whole_metabolome_aic_models.txt identify_nets_libid.csv lindaDepletolites.csv canopus_filtered_tidy.csv identify_nets_libid.numbers n_mulreg.csv Moorea2017_MolNetEnhancer.csv cyto_depletes.csv multiple_reg_coefficients_whole_metabolome_threemodels.txt lindaPieChartMetadata.csv model_rmse_r2.txt deplete_networks.csv org_random_residsmodels.txt multiple_reg_coefficients_whole_metabolome.txt canopus_df.csv aic_model_selection.txt dorcierr_day_benthic_exometabolite_features.csv tukey_log2.txt linear_model_dredged_variabled.txt dredge_model_selection.txt raw metabolomics SIRIUS_Zodiac_converted.csv Dorcierr_DIC_TopDepletolitesClass.csv Dorcierr_CCA_TopDepletolitesClassified.csv Extraction_efficiency.csv Morrea_Feayures-Table_all_Gap-Filled5.csv Analog-Hits.tsv Library-Hits.tsv Node_info.tsv...","categories": [],
        "tags": [],
        "url": "/%7B%7B%20site.baseurl%20%7D%7D/contents/",
        "teaser": null
      },{
        "title": "12182019dorcierr_working_pipeline.r",
        "excerpt":"## Script Written by Zach Quinlan 06/19/19 # Re-organization of DORCIERR_FCM_fDOM.R because it needs to be cleaner 07/15/2019 # Only working on daytime exudation and remineralization # Rewritten with changes to RR3 starting pipeline 10/11/2019 # 16s rRNA amplicon seque3nce data added in a upaded 10/18/2019 # Added in ClassyFire...","categories": [],
        "tags": [],
        "url": "/%7B%7B%20site.baseurl%20%7D%7D/contents/12182019DORCIERR_working_pipeline.R/",
        "teaser": null
      },{
        "title": "Microbial_pipeline.r",
        "excerpt":"## Script Written by Zach Quinlan 06/19/19 # Re-organization of DORCIERR_FCM_fDOM.R because it needs to be cleaner 07/15/2019 # Only working on daytime exudation and remineralization # Rewritten with changes to RR3 starting pipeline 10/11/2019 # 16s rRNA amplicon seque3nce data added in a upaded 10/18/2019 # Added in ClassyFire...","categories": [],
        "tags": [],
        "url": "/%7B%7B%20site.baseurl%20%7D%7D/contents/Microbial_pipeline.R/",
        "teaser": null
      },{
        "title": "Bicluster.py",
        "excerpt":"``` import seaborn as sns import matplotlib.pyplot as plt import pandas as pd import argparse #Arguments parser = argparse.ArgumentParser() parser.add_argument(‘infile’, help = ‘csv_dataframe’) parser.add_argument(‘sample_column’, help = ‘what is the name of the sample_column’) args = parser.parse_args() #Reading in file df = pd.read_csv(args.infile) #Setting font and color pallette #current_palette = sns.color_palette(“Greys”,...","categories": [],
        "tags": [],
        "url": "/%7B%7B%20site.baseurl%20%7D%7D/contents/bicluster.py/",
        "teaser": null
      },{
        "title": "Canopus_model.r",
        "excerpt":"angular_transform &lt;- function(x) { asin(sqrt(x)) } # cutoff_canl2 &lt;- canopus_chemonnt_tidy%&gt;% # group_by(canopus_annotation)%&gt;% # filter(max(canopus_probability) &gt;= 0.3)%&gt;% # ungroup()%&gt;% # separate(CLASS_STRING, paste('L', 1:5, sep = \"\"), sep = \";\")%&gt;% # group_by(L3)%&gt;% # filter(mean(canopus_probability) &gt;= 0.001)%&gt;% # group_by(L3)%&gt;% # mutate(max_prob = max(canopus_probability))%&gt;% # select(name, L3, max_prob)%&gt;% # unique()%&gt;% # spread(L3, max_prob) #...","categories": [],
        "tags": [],
        "url": "/%7B%7B%20site.baseurl%20%7D%7D/contents/canopus_model.R/",
        "teaser": null
      },{
        "title": "Current_dorcierr_working_pipeline.r",
        "excerpt":"## Script Written by Zach Quinlan 06/19/19 # Re-organization of DORCIERR_FCM_fDOM.R because it needs to be cleaner 07/15/2019 # Only working on daytime exudation and remineralization # Rewritten with changes to RR3 starting pipeline 10/11/2019 # 16s rRNA amplicon seque3nce data added in a upaded 10/18/2019 # Added in ClassyFire...","categories": [],
        "tags": [],
        "url": "/%7B%7B%20site.baseurl%20%7D%7D/contents/current_dorcierr_working_pipeline.R/",
        "teaser": null
      },{
        "title": "Metabolite_modeled_lability.r",
        "excerpt":"## Script Written by Zach Quinlan 06/19/19 # Re-organization of DORCIERR_FCM_fDOM.R because it needs to be cleaner 07/15/2019 # Only working on daytime exudation and remineralization # Rewritten with changes to RR3 starting pipeline 10/11/2019 # 16s rRNA amplicon seque3nce data added in a upaded 10/18/2019 # Added in ClassyFire...","categories": [],
        "tags": [],
        "url": "/%7B%7B%20site.baseurl%20%7D%7D/contents/metabolite_modeled_lability.R/",
        "teaser": null
      },{
        "title": "Metabolite_trends.r",
        "excerpt":"## Script Written by Zach Quinlan 06/19/19 # Re-organization of DORCIERR_FCM_fDOM.R because it needs to be cleaner 07/15/2019 # Only working on daytime exudation and remineralization # Rewritten with changes to RR3 starting pipeline 10/11/2019 # 16s rRNA amplicon seque3nce data added in a upaded 10/18/2019 # Added in ClassyFire...","categories": [],
        "tags": [],
        "url": "/%7B%7B%20site.baseurl%20%7D%7D/contents/metabolite_trends.R/",
        "teaser": null
      },]
