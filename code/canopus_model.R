angular_transform <- function(x) {
  asin(sqrt(x))
}


# cutoff_canl2 <- canopus_chemonnt_tidy%>%
#   group_by(canopus_annotation)%>%
#   filter(max(canopus_probability) >= 0.3)%>%
#   ungroup()%>%
#   separate(CLASS_STRING, paste('L', 1:5, sep = ""), sep = ";")%>%
#   group_by(L3)%>%
#   filter(mean(canopus_probability) >= 0.001)%>%
#   group_by(L3)%>%
#   mutate(max_prob = max(canopus_probability))%>%
#   select(name, L3, max_prob)%>%
#   unique()%>%
#   spread(L3, max_prob)
#   
cutoff_canopus <- canopus_anotations%>%
  gather(classification, value, 2:ncol(.))%>%
  group_by(classification)%>%
  filter(max(value) >= 0.4)%>%
  spread(classification, value)

canopus_mulreg <- mul_reg%>%
  filter(NOSC < 0)%>%
  group_by(feature_number, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  select(feature_number, Organism, log2_change, `row m/z`, NOSC, logxic)%>%
  inner_join(cutoff_canopus%>%
               rename(feature_number = name)%>%
               mutate(feature_number = as.character(feature_number))%>%
               mutate_if(is.numeric, angular_transform), by = 'feature_number')%>%
  left_join(sirius_zodiac_anotations%>%
              select(feature_number, ZodiacScore)%>%
              mutate(feature_number = as.character(feature_number)), by = 'feature_number')%>%
  filter(ZodiacScore >= .98)%>%
  select(-ZodiacScore)%>%
  unite(row, c('feature_number', 'Organism'), sep = "_")%>%
  column_to_rownames('row')

recipe <- canopus_mulreg%>%
  recipe(log2_change ~ .)%>%
  step_center(all_predictors())%>%
  step_scale(all_predictors())

#checking gaussian
pdf("./plots/canopus_gaussian_check.pdf")
gaussian_check <- prep(recipe)%>%
  juice()%>%
  gather(variable, value)%>%
  inner_join(canopus_coeff%>% select(variable),
             by = 'variable')%>%
  group_by(variable)%>%
  nest()%>%
  mutate(non_transformed = map(data, ~ .x$value%>%
                                 na.omit()%>%
                                 car::qqPlot(ylab = variable, xlab = "Normal quantiles",
                                             main = paste("Non-transformed", variable))))
dev.off()

##canopus training set testing

canopus_split <- rsample::initial_split(canopus_mulreg, 9/10)

canopus_training <- rsample::training(canopus_split)
canopus_test <- rsample::testing(canopus_split)

canopus_x <- canopus_training[2:ncol(canopus_training)]%>%
  as.matrix()

canopus_y <- canopus_training['log2_change']%>%
  as.matrix()

canopus_test_x <-  canopus_test[2:ncol(canopus_test)]%>%
  as.matrix()

canopus_test_y <- canopus_test['log2_change']%>%
  as.matrix()

#running 
glmnet_canopus <- cv.glmnet(canopus_x, canopus_y, parallel = TRUE)

lambda_se <- glmnet_canopus$lambda.1se

#running model with parameters after lasso
canopus_model <- glmnet(canopus_x, canopus_y, intercept = TRUE, lambda = lambda_se, standarize = TRUE)

canopus_coeff <- coef(canopus_model)%>%
  as.matrix()%>%
  as.data.frame()%>%
  rename(coefficients = 1)%>%
  rownames_to_column("variable")%>%
  filter(coefficients != 0)

intercept <- canopus_coeff[[1,2]]

predicted_training <- predict(canopus_model, newx = canopus_x)
predicted <- predict(canopus_model, newx = canopus_test_x)

rootmse <- rmse_vec(as.vector(canopus_test_y), as.vector(predicted))
rootmse_training <- rmse_vec(as.vector(canopus_y), as.vector(predicted_training))

canopus_model$dev.ratio

predicted%>%
  as.data.frame()%>%
  rename(predicted = 1)%>%
  add_column(actuals = canopus_test_y)%>%
  rownames_to_column('sample')%>%
  ggplot() +
  geom_point(aes(sample, actuals), color = 'blue') +
  geom_point(aes(sample, predicted), color = 'red')



canopus_weighted_lability <- feature_stats_wdf%>%
  filter(Timepoint == 'T0',
         DayNight == 'Day')%>%
  inner_join(canopus_mulreg%>%
              rownames_to_column("row")%>%
              separate(row, c('feature_number', 'Organism'), sep = "_")%>%
              select(-log2_change)%>%
              gather(variable, value, 3:ncol(.))%>%
              inner_join(canopus_coeff, by = 'variable')%>%
              mutate(modeled_lab = value*coefficients)%>%
              select(feature_number, Organism, modeled_lab)%>%
              group_by(feature_number, Organism)%>%
              summarize_if(is.numeric, sum)%>%
              ungroup()%>%
              mutate(modeled_lab = modeled_lab + intercept),
            by = c('feature_number', 'Organism'))%>%
  group_by(Organism, Replicate)%>%
  mutate(weighted_lability = spatstat::weighted.median(modeled_lab, log10))%>%
  # mutate(weighted_lability = mean(modeled_lab))%>%
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
# 
# # dredge_n_select <- MuMIn::dredge(lm_test_canopus)
# # test <- AIC(lm_test_canopus)
# 

#Linear model
weight_lability_lm <- canopus_weighted_lability%>%
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

pdf('~/Documents/GitHub/DORCIERR/data/plots/canopus_weighted_lability.pdf', width = 15, height = 10)
canopus_weighted_lability%>%
  ggplot(aes(weighted_lability, cells_ul)) +
  geom_point(aes(color = Organism), stat = 'identity', size = 5) +
  geom_errorbar(aes(ymin = cells_ul - y_err, ymax = cells_ul + y_err)) +
  geom_errorbarh(aes(xmin = weighted_lability - x_err, xmax = weighted_lability + x_err))+
  scale_color_manual(values = org_colors_no_water) +
  geom_smooth(method = 'lm') +
  labs(y = bquote(Specific ~growth ~rate ~('Cells'~µL^-1 ~hr^-1)), x = "Metabolite pool modeled lability") +
  geom_text(aes(x = -1.4, y = 0.09,
                label = paste("p-value: ", wlab_p%>%
                                formatC(format = "e", digits = 2), sep = "")), size = 9) +
  geom_text(aes(x = -1.4, y = 0.086,
                label = paste("F statistic: ", wlab_f%>%
                                round(digits = 4), sep = "")), size = 9) +
  geom_text(aes(x = -1.4, y = 0.082,
                label = paste("r²: ", wlab_r2%>%
                                round(digits = 4), sep = "")), size = 9) +
  geom_text(aes(x = -1.4, y = 0.078,
                label = paste("Cells µL^-1", " = ", wlab_slope%>%
                                round(digits = 2), "*lability + ", wlab_intercept%>%
                                round(digits = 2), sep = "")), size = 9) +
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
