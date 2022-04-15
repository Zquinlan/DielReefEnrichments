craig <- dom_stats_wdf%>%
  ungroup()%>%
  left_join(networking%>%
              select(network, feature_number),
            by = 'feature_number')%>%
    mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                               TRUE ~ network))%>%
  left_join(all_activity,
            by = c('net_act', 'DayNight'))%>%
  left_join(ecoNet%>%
              rename(feature_number = 1)%>%
              mutate(feature_number = as.character(feature_number))%>%
              select(1:3), 
            by = c('feature_number', 'network'))%>%
  filter(Timepoint == 'T0')%>%
  select(net_act, network, DayNight, activity, ecoNetConsensus)%>%
  unique()
  
write_csv(craig, '~/Downloads/craig.csv')

zachLog2 <- log2_change_vals%>%
  select(feature_number, DayNight, log2_change)%>%
  left_join(networking%>%
             select(network, feature_number),
           by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  group_by(net_act, DayNight, log2_change)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()
  

zach <- dom_stats_wdf%>%
  ungroup()%>%
  left_join(networking%>%
              select(network, feature_number),
            by = 'feature_number')%>%
  mutate(net_act = case_when(network == -1 ~ -as.numeric(feature_number),
                             TRUE ~ network))%>%
  left_join(all_activity,
            by = c('net_act', 'DayNight'))%>%
  left_join(ecoNet%>%
              rename(feature_number = 1)%>%
              mutate(feature_number = as.character(feature_number))%>%
              select(1:3), 
            by = c('feature_number', 'network'))%>%
  inner_join(feature_stats_wdf%>%
               ungroup()%>%
               select(feature_number)%>%
               unique(),
             by = 'feature_number')%>%
  filter(Timepoint == 'T0')%>%
  select(net_act, network, DayNight, activity, ecoNetConsensus, Organism, Replicate, xic)%>%
  # group_by(net_act, DayNight, activity, Organism, ecoNetConsensus, Replicate)%>%
  # summarize_if(is.numeric, sum, na.rm = TRUE)%>%
  # ungroup()%>%
  group_by(net_act, DayNight, activity, ecoNetConsensus)%>%
  mutate(std = sd(xic))%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  separate(ecoNetConsensus, c('superclass_consensus', 'class_consensus', 'subclass_consensus'), remove = FALSE, sep = ';')%>%
  left_join(zachLog2, by = c('net_act', 'DayNight'))
  
pdf('~/Downloads/NetworkAveragedAverageProduction.pdf', width = 15, height = 10)
zach%>%
  ggplot(aes(DayNight, xic, color = activity)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_color_manual(labels = c('Depletolite', 'Recalcitrant'), values = c('#EBCC2A', "#006658")) +
  # scale_color_manual(values = c('#F09837', 'grey')) +
  facet_wrap(~superclass_consensus, scales = 'free_y') +
  gen_theme()

zach%>%
  ggplot(aes(xic, log2_change)) +
  geom_point(stat = 'identity', alpha = 0.1, color = 'grey') +
  geom_smooth(aes(color = activity), method = 'lm', size = 2) +
  scale_x_log10() +
  facet_wrap(~superclass_consensus) +
  scale_color_manual(labels = c('Depletolite', 'Recalcitrant'), values = c('#EBCC2A', "#006658")) +
  gen_theme()


dev.off()  

