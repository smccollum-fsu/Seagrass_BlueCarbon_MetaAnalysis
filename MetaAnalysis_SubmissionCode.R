

library(readxl)
library(dplyr)
library(tidyr)
library(brms)
library(emmeans)
library(writexl)
library(bayestestR)
library(posterior)



df <- read_excel("./Desktop/MetaNewFolder/SeagrassBlueCarbon_MetaAnalysisDataset_Submission.xlsx") 

## FILTER THE DATASET

## Add observation numbers
df$Observation <- 1:nrow(df)

df_global <- df

## Double check that the NAs in the error is 0
df_global$SedCorgStockError[is.na(df_global$SedCorgStockError)] <- "0"

## Add a "Publication" column. Originally there was a "Paper" column with the author/date but we
    ## replaced it with this for the published dataset
df_global <- df_global %>% 
  group_by(Citation) %>% 
  mutate(PAPER = cur_group_id()) %>% 
  ungroup()


## Make some of the variables numeric
df_global <- df_global %>% 
  transform(SedCorgStock = as.numeric(SedCorgStock)) %>% 
  transform(SedCorgStockError = as.numeric(SedCorgStockError)) %>% 
  transform(SedCorgStockN = as.numeric(SedCorgStockN))

## Change some column names to those orignially used in the code.
    ## For the published dataset, we changed column names to be more up-to-date and descriptive.
    ## Here, we return them to the originally names to make it compatible with the code.

colnames(df_global)[4] <- "SpeciesSize"
colnames(df_global)[12] <- "Geomorphology"
colnames(df_global)[13] <- "Subtidal" ## Tidal_Class was originally a binary "Subtidal" with 1 (subtidal) and 0 (intertidal)


### MAKING DF FOR EACH OF THE MODELS

df_bioregion <- df_global %>% ## Only measurements where Bioregion  is reported
  filter(is.na(Bioregion) == FALSE)

df_sp <- df_bioregion %>%  ## Only measurements where species is reported, filter out poorly-reported species
  filter(
    is.na(Bioregion) == FALSE,
    Species != "Halodule emarginata",
    Species != "Halophila beccarii", 
    Species != "Halodule pinifolia", 
    Species != "Cymodocea rotundata", 
    Observation != 131)  ## A single observation of H. decipiencs in the TA # core = 1

## Remove observations with multiple species
df_sp <- df_sp[!grepl(",", df_sp$Species), ]
df_sp <- df_sp[!grepl("spp", df_sp$Species), ]

df_size_int <- df_sp %>% filter(Bioregion != "Temperate North Pacific" & Bioregion != "Temperate North Atlantic")
## Multiple sizes were not observed in these bioregions, thus their omission

df_bioregion_ss <- df_sp %>% ## Species size
  filter(is.na(Bioregion) == FALSE)


## NOTE: "Geomorph" or "Coastal" refer to hydrologic setting. These were the previous phrasings used in the code and manuscript.

df_geomorph <- df_sp %>%  ## Only measurements where Geomorphology AND species are present
  filter(is.na(Geomorphology) == FALSE)

df_geomorph_int <- df_geomorph %>% filter(Bioregion != "Mediterranean" & Bioregion != "Temperate Southern Ocean" & Bioregion != "Temperate North Atlantic")

df_td <- df_sp %>%   ## Only measurements where tidal depth AND species is reported
  filter(is.na(Subtidal) == FALSE)

df_td_int <- df_td %>% filter(Bioregion != "Mediterranean" & Bioregion != "Tropical Atlantic" & Bioregion != "Temperate Southern Ocean")

#### Establish priors ####

priors <- c(prior(student_t(3, 30, 200), lb = 0, class = b), ## prior assumption of the "slope"
            prior(student_t(3, 0, 200), class = sd)) ## prior assumption of the error around the mean


nullpriors<- c(prior(student_t(3, 30, 200), lb = 0, class = Intercept), ## prior assumption of the position of the intercept
               prior(student_t(3, 0, 200), class = sd)) ## prior assumption of the error around the mean





#### **** ####     MODELS     ##### **** ####


#### BIOREGION DATASET -- MODELLING THE INFLUENCE OF BIOREGION ####

mod_bioregion <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ Bioregion-1 + (1|PAPER),
                     df_bioregion,
                     prior = priors,
                     iter = 100000, cores = 4,
                     family = gaussian(link = "identity"),
                     save_pars = save_pars(all=TRUE)) ## Diagnostics clear

mod_bioregion_null <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ 1 + (1|PAPER),
                          df_bioregion,
                          family = gaussian(link = "identity"),
                          prior = nullpriors,
                          iter = 100000, cores = 4,
                          save_pars = save_pars(all=TRUE)) ## Diagnostics clear


#### Posterior model weights relative to null model

mcse(mod_bioregion)

post_prob(mod_bioregion, mod_bioregion_null) # 1.0

#### ####

mod_bioregion_em <- emmeans(mod_bioregion, "Bioregion", point.est = mean, mode = "hdi")

pairs(mod_bioregion_em)

bio_em_df <- as.data.frame(pairs(mod_bioregion_em))

write_xlsx(bio_em_df, "./Desktop/EMMEANSDataframes/BioregionPairs.xlsx")


###################################################

#### SPECIES DATASET -- USED FOR MODELLING SPECIES, SPECIES SIZE ####

## Null

sp.mod_null_sub <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ 1 + (1|PAPER),
                       df_sp,
                       prior = nullpriors,
                       iter = 100000, cores = 4,
                       family = gaussian(link = "identity"),
                       control = list(adapt_delta = .99),
                       save_pars = save_pars(all=TRUE)) ## Diagnostics clear

## Bioregion

sp.mod_bioregion <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~  Bioregion-1 + (1|PAPER),
                        df_sp,
                        prior = priors,
                        iter = 100000, cores = 4,
                        family = gaussian(link = "identity"),
                        control = list(adapt_delta = .99),
                        save_pars = save_pars(all=TRUE)) ## Diagnostics clear

#### Posterior model weights relative to null model

post_prob(sp.mod_bioregion, sp.mod_null) ## 1.0

#### ####


## Species Identity

sp.mod_sp <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE)  ~  Species-1 + (1|PAPER),
                 df_sp,
                 prior = priors,
                 iter = 100000, cores = 4,
                 family = gaussian(link = "identity"),
                 control = list(adapt_delta = .99),
                 save_pars = save_pars(all=TRUE)) ## Diagnostics clear

sp.mod_bioregion_sp_int <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~  (Bioregion:Species)-1 + (1|PAPER),
                               df_sp,
                               prior = priors,
                               iter = 100000, cores = 4,
                               family = gaussian(link = "identity"),
                               control = list(adapt_delta = .99),
                               save_pars = save_pars(all=TRUE)) ## Diagnostics clear


### MCSE
mcse(sp.mod_sp)
mcse(sp.mod_bioregion_sp_int)
###

#### Posterior model weights relative to null model

post_prob(sp.mod_sp, sp.mod_null) ## 1.0
post_prob(sp.mod_bioregion_sp_int, sp.mod_null) ## 1.0

####

#### EMMEANS !!!

mod_species_em <- emmeans(sp.mod_sp, "Species", point.est = mean, level = .95, mode = "hdi")

pairs(mod_species_em)

sp_em_df <- as.data.frame(pairs(mod_species_em))

write_xlsx(sp_em_df, "./Desktop/EMMEANSDataframes/SpeciesPairs.xlsx")


mod_species_int_em <- emmeans(sp.mod_bioregion_sp_int, "Bioregion", "Species", point.est = mean, level = .95, mode = "hdi")

pairs(mod_species_int_em)

sp_int_em_df <- as.data.frame(pairs(mod_species_int_em))

write_xlsx(sp_int_em_df, "./Desktop/EMMEANSDataframes/SpeciesInteractivePairs.xlsx")


mod_species_int_em_2 <- emmeans(sp.mod_bioregion_sp_int, "Species", "Bioregion", point.est = mean, level = .95, mode = "hdi")

pairs(mod_species_int_em_2)

sp_int_em_df_2 <- as.data.frame(pairs(mod_species_int_em_2))

write_xlsx(sp_int_em_df_2, "./Desktop/EMMEANSDataframes/SpeciesInteractivePairsInverted.xlsx")




## Size

sp.mod_size <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ SpeciesSize-1 + (1|PAPER),
                   df_sp,
                   prior = priors,
                   iter = 100000, cores = 4,
                   family = gaussian(link = "identity"),
                   control = list(adapt_delta = .99),
                   save_pars = save_pars(all=TRUE)) ## Diagnostics clear

sp.mod_size_int_null <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ 1 + (1|PAPER),
                            df_size_int,
                            prior = nullpriors,
                            iter = 100000, cores = 4,
                            family = gaussian(link = "identity"),
                            control = list(adapt_delta = .99),
                            save_pars = save_pars(all=TRUE)) ## Diagnostics clear

sp.mod_bioregion_size_int <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~  (Bioregion:SpeciesSize)-1 + (1|PAPER),
                                 df_size_int,
                                 prior = priors,
                                 iter = 100000, cores = 4,
                                 family = gaussian(link = "identity"),
                                 control = list(adapt_delta = .99),
                                 save_pars = save_pars(all=TRUE)) ## Diagnostics clear


##MCSE 

mcse(sp.mod_size)
mcse(sp.mod_bioregion_size_int)


#### Posterior model weights relative to null model

post_prob(sp.mod_size, sp.mod_null) ## .9999
post_prob(sp.mod_bioregion_size_int, sp.mod_size_int_null) ## 1.0


####################################################

#### EMMEANS !!!!

mod_size_em <- emmeans(sp.mod_size, "SpeciesSize", point.est = mean, level = .95, mode = "hdi")

pairs(mod_size_em)

size_em_df <- as.data.frame(pairs(mod_size_em))

write_xlsx(size_em_df, "./Desktop/EMMEANSDataframes/SizePairs.xlsx")



mod_size_int_em <- emmeans(sp.mod_bioregion_size_int, "SpeciesSize", "Bioregion", point.est = mean, level = .95, mode = "hdi")

pairs(mod_size_int_em)

size_em_df <- as.data.frame(pairs(mod_size_int_em))

write_xlsx(size_em_df, "./Desktop/EMMEANSDataframes/SizeInteractivePairs.xlsx")





#### TIDAL CLASS -- MODELLING TIDAL CLASS, SPECIES ID, AND SPECIES SIZE

## NULL

td.mod_null <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ 1 + (1|PAPER),
                   df_td,
                   prior = nullpriors,
                   iter = 100000, cores = 4,
                   family = gaussian(link = "identity"),
                   control = list(adapt_delta = .99),
                   save_pars = save_pars(all=TRUE)) ## Diagnostics clear

#### Bioregion

td.mod_bio <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ Bioregion-1 + (1| PAPER),
                  df_td,
                  prior = priors,
                  iter = 100000, cores = 4,
                  family = gaussian(link = "identity"),
                  control = list(adapt_delta = .99),
                  save_pars = save_pars(all=TRUE)) ## Diagnostics clear

#### Posterior model weights relative to null model

post_prob(td.mod_bio, td.mod_null) # 0.05542683

####

## Tidal class

td.mod_td <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ (Subtidal-1) + (1|PAPER),
                 df_td,
                 prior = priors,
                 iter = 100000, cores = 4,
                 family = gaussian(link = "identity"),
                 control = list(adapt_delta = .99),
                 save_pars = save_pars(all=TRUE)) ## Diagnostics clear

td.mod_td_int_null <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ 1 + (1|PAPER),
                          df_td_int,
                          prior = nullpriors,
                          iter = 100000, cores = 4,
                          family = gaussian(link = "identity"),
                          control = list(adapt_delta = .99),
                          save_pars = save_pars(all=TRUE)) ## Diagnostics clear

td.mod_td_bio_int <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ (Bioregion:Subtidal)-1 + (1|PAPER),
                         df_td_int,
                         prior = priors,
                         iter = 100000, cores = 4,
                         family = gaussian(link = "identity"),
                         control = list(adapt_delta = .99),
                         save_pars = save_pars(all=TRUE)) ## Diagnostics clear


## mcse
mcse(td.mod_td)
mcse(td.mod_td_bio_int)


#### Posterior model weights relative to null model

post_prob(td.mod_td, td.mod_null) ## 0.1421144
post_prob(td.mod_td_bio_int, td.mod_td_int_null) # 0.000037 


##### #####

#### EMMEANS !!!!!

mod_td_em <- emmeans(td.mod_td, "Subtidal", point.est = mean, level = .95, mode = "hdi")

pairs(mod_td_em)

td_em_df <- as.data.frame(pairs(mod_td_em))

write_xlsx(td_em_df, "./Desktop/EMMEANSDataframes/TidalClassPairs.xlsx")


mod_td_int_em <- emmeans(td.mod_td_bio_int, "Subtidal", "Bioregion", point.est = mean, level = .95, mode = "hdi")

pairs(mod_td_int_em)

td_int_em_df <- as.data.frame(pairs(mod_td_int_em))

write_xlsx(td_int_em_df, "./Desktop/EMMEANSDataframes/TidalClassInteractivePairs.xlsx")




####



## Species ID

td.mod_sp <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ (Species)-1 + (1|PAPER),
                 df_td,
                 prior = priors,
                 iter = 100000, cores = 4,
                 family = gaussian(link = "identity"),
                 control = list(adapt_delta = .99),
                 save_pars = save_pars(all=TRUE)) ## Diagnostics clear

td.mod_sp_bio <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ (Bioregion:Species)-1 + (1|PAPER),
                     df_td,
                     prior = priors,
                     iter = 100000, cores = 4,
                     family = gaussian(link = "identity"),
                     control = list(adapt_delta = .99),
                     save_pars = save_pars(all=TRUE)) ## Diagnostics clear


#### Posterior model weights relative to null model

post_prob(td.mod_sp, td.mod_null) ## 0.2588
post_prob(td.mod_sp_bio, td.mod_null) # 0.001279452

#### ####

## Species size

df_td_size_int <- df_td %>% filter(Bioregion != "Temperate North Atlantic" & Bioregion != "Temperate North Pacific")

td.mod_size <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ SpeciesSize-1 + (1|PAPER),
                   df_td,
                   prior = priors,
                   iter = 100000, cores = 4,
                   family = gaussian(link = "identity"),
                   control = list(adapt_delta = .99),
                   save_pars = save_pars(all=TRUE)) ## Diagnostics clear

td.mod_size_int_null <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ 1 + (1|PAPER),
                            df_td_size_int,
                            prior = nullpriors,
                            iter = 100000, cores = 4,
                            family = gaussian(link = "identity"),
                            control = list(adapt_delta = .99),
                            save_pars = save_pars(all=TRUE)) ## Diagnostics clear

td.mod_size_bio <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ (Bioregion:SpeciesSize)-1 + (1|PAPER),
                       df_td_size_int,
                       prior = priors,
                       iter = 100000, cores = 4,
                       family = gaussian(link = "identity"),
                       control = list(adapt_delta = .99),
                       save_pars = save_pars(all=TRUE)) ## Diagnostics clear

#### Posterior model weights relative to null model

post_prob(td.mod_size, td.mod_null) ## 0.0872
post_prob(td.mod_size_bio, td.mod_size_int_null) # 0.997 

#### #####

## Geomorphic setting

df_td_geo <- df_td %>% filter(is.na(Geomorphology) == FALSE)

td.mod_geomorph_null <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ 1 + (1|PAPER),
                            df_td_geo,
                            prior = nullpriors,
                            iter = 100000, cores = 4,
                            family = gaussian(link = "identity"),
                            control = list(adapt_delta = .99),
                            save_pars = save_pars(all=TRUE)) ## Diagnostics clear

td.mod_geomorph <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ Geomorphology-1 + (1|PAPER),
                       df_td_geo,
                       prior = priors,
                       iter = 100000, cores = 4,
                       family = gaussian(link = "identity"),
                       control = list(adapt_delta = .99),
                       save_pars = save_pars(all=TRUE)) ## Diagnostics clear

df_td_geo_int <- df_td_geo %>% filter(Bioregion != "Mediterranean" & Bioregion != "Temperate Southern Ocean" & Bioregion != "Temperate North Atlantic")

td.mod_geomorph_bio_null <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ 1 + (1|PAPER),
                                df_td_geo_int,
                                prior = nullpriors,
                                iter = 100000, cores = 4,
                                family = gaussian(link = "identity"),
                                control = list(adapt_delta = .99),
                                save_pars = save_pars(all=TRUE)) ## Diagnostics clear

td.mod_geomorph_bio <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ (Bioregion:Geomorphology)-1 + (1|PAPER),
                           df_td_geo_int,
                           prior = priors,
                           iter = 100000, cores = 4,
                           family = gaussian(link = "identity"),
                           control = list(adapt_delta = .99),
                           save_pars = save_pars(all=TRUE)) ## Diagnostics clear


#### Posterior model weights relative to null model

post_prob(td.mod_geomorph, td.mod_geomorph_null) #0.14178
post_prob(td.mod_geomorph_bio, td.mod_geomorph_bio_null) # 0.000019

###################################################


#### GEOMORPHIC SETTING -- MODELLING GEOMORPHIC SETTING, SPECIES, SPECIES SIZE

## Null
geo.mod_null <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ 1 + (1|PAPER),
                    df_geomorph,
                    prior = nullpriors,
                    iter = 100000, cores = 4,
                    family = gaussian(link = "identity"),
                    control = list(adapt_delta = .99),
                    save_pars = save_pars(all=TRUE)) ## Diagnostics clear

## Bioregion

geo.mod_bio <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ (Bioregion-1) + (1|PAPER),
                   df_geomorph,
                   prior = priors,
                   iter = 100000, cores = 4,
                   family =  gaussian(link = "identity"),
                   control = list(adapt_delta = .99),
                   save_pars = save_pars(all=TRUE)) ## Diagnostics clear

#### Posterior model weights relative to null model

post_prob(geo.mod_bio, geo.mod_null) # > 0.9999

#### #####



## Geomorphic setting

geo.mod_geo <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ (Geomorphology-1) + (1|PAPER),
                   df_geomorph,
                   prior = priors,
                   iter = 100000, cores = 4,
                   family =  gaussian(link = "identity"),
                   control = list(adapt_delta = .99),
                   save_pars = save_pars(all=TRUE)) ## Diagnostics clear

mcse(geo.mod_geo)

geo.mod_int_null <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ 1 + (1|PAPER),
                        df_geomorph_int,
                        prior = nullpriors,
                        iter = 100000, cores = 4,
                        family = gaussian(link = "identity"),
                        control = list(adapt_delta = .99),
                        save_pars = save_pars(all=TRUE)) ## Diagnostics clear

geo.mod_geo_bio <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ ((Bioregion:Geomorphology)-1) + (1|PAPER),
                       df_geomorph_int,
                       prior = priors,
                       iter = 100000, cores = 4,
                       family =  gaussian(link = "identity"),
                       control = list(adapt_delta = .99),
                       save_pars = save_pars(all=TRUE)) ## Diagnostics clear

#### Posterior model weights relative to null model

post_prob(geo.mod_geo, geo.mod_null) # 0.12596
post_prob(geo.mod_geo_bio, geo.mod_int_null) # > 0.00000081 

post_prob(geo.mod_geo_bio_sub, geo.mod_int_null_sub)

#### ####

#### EMMEANS !!!!

mcse(geo.mod_geo)

mod_geo_em <- emmeans(geo.mod_geo, "Geomorphology", point.est = mean, level = .95, mode = "hdi")

pairs(mod_geo_em)

geo_em_df <- as.data.frame(pairs(mod_geo_em))

write_xlsx(geo_em_df, "./Desktop/EMMEANSDataframes/GeomorphicSettingPairs.xlsx")

mcse(geo.mod_geo_bio)

mod_geo_int_em <- emmeans(geo.mod_geo_bio, "Geomorphology", "Bioregion", point.est = mean, level = .95, mode = "hdi")

pairs(mod_geo_int_em)

geo_int_em_df <- as.data.frame(pairs(mod_geo_int_em))

write_xlsx(geo_int_em_df, "./Desktop/EMMEANSDataframes/GeomorphicSettingInteractivePairs.xlsx")

#### Species

geo.mod_sp <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ (Species-1) + (1|PAPER),
                  df_geomorph,
                  prior = priors,
                  iter = 100000, cores = 4,
                  family = gaussian(link = "identity"),
                  control = list(adapt_delta = .99),
                  save_pars = save_pars(all=TRUE)) ## Diagnostics clear

geo.mod_sp_bio <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ ((Bioregion:Species-1)) + (1|PAPER),
                      df_geomorph,
                      prior = priors,
                      iter = 100000, cores = 4,
                      family = gaussian(link = "identity"),
                      control = list(adapt_delta = .99),
                      save_pars = save_pars(all=TRUE)) ## Diagnostics clear

#### Posterior model weights relative to null model

post_prob(geo.mod_sp, geo.mod_null) # > 0.9999
post_prob(geo.mod_sp_bio, geo.mod_null) # > 0.99

#### ####

#### Species Size

df_geomorph_size_int <- df_geomorph %>% filter(Bioregion != "Temperate North Atlantic" & Bioregion != "Temperate North Pacific")

geo.mod_size <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ (SpeciesSize-1) + (1|PAPER),
                    df_geomorph,
                    prior = priors,
                    iter = 100000, cores = 4,
                    family = gaussian(link = "identity"),
                    control = list(adapt_delta = .99),
                    save_pars = save_pars(all=TRUE)) ## Diagnostics clear

geo.mod_size_int_null <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ 1 + (1|PAPER),
                             df_geomorph_size_int,
                             prior = nullpriors,
                             iter = 100000, cores = 4,
                             family = gaussian(link = "identity"),
                             control = list(adapt_delta = .99),
                             save_pars = save_pars(all=TRUE)) ## Diagnostics clear

geo.mod_size_bio <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ ((Bioregion:SpeciesSize)-1) + (1|PAPER),
                        df_geomorph_size_int,
                        prior = priors,
                        iter = 100000, cores = 4,
                        family = gaussian(link = "identity"),
                        control = list(adapt_delta = .99),
                        save_pars = save_pars(all=TRUE)) ## Diagnostics clear

#### Posterior model weights relative to null model

post_prob(geo.mod_size, geo.mod_null) # 0.8316405
post_prob(geo.mod_size_bio, geo.mod_size_int_null) # .99999 

#### ####

df_td_geo <- df_geomorph %>% filter(Subtidal != "NA")

geo.mod_td <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ (Subtidal-1) + (1|PAPER),
                  df_td_geo,
                  prior = priors,
                  iter = 100000, cores = 4,
                  family = gaussian(link = "identity"),
                  control = list(adapt_delta = .99),
                  save_pars = save_pars(all=TRUE)) ## Diagnostics clear

df_td_geo_int <- df_td_geo %>% filter(Bioregion != "Tropical Atlantic", Bioregion != "Temperate Southern Ocean" & Bioregion != "Mediterranean")

geo.mod_td_int_null <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ 1 + (1|PAPER),
                           df_td_geo_int,
                           prior = nullpriors,
                           iter = 100000, cores = 4,
                           family = gaussian(link = "identity"),
                           control = list(adapt_delta = .99),
                           save_pars = save_pars(all=TRUE)) ## Diagnostics clear

geo.mod_td_bio <- brm(SedCorgStock | se(SedCorgStockError, sigma = TRUE) ~ ((Bioregion:Subtidal)-1) + (1|PAPER),
                      df_td_geo_int,
                      prior = priors,
                      iter = 100000, cores = 4,
                      family = gaussian(link = "identity"),
                      control = list(adapt_delta = .99),
                      save_pars = save_pars(all=TRUE)) ## Diagnostics clear

post_prob(geo.mod_td, td.mod_geomorph_null) # 0.0557
post_prob(geo.mod_td_bio, geo.mod_td_int_null) #  0.00000012 



### END OF CODE! You have produced models for each dataframe and variable. 



### Printing out the errors.

GlobalMCSE <- as.data.frame(mcse(mod_bioregion_null))

write_xlsx(as.data.frame(mcse(mod_bioregion_null)), "./Desktop/MCSE/GlobalMCSE.xlsx")

BioMCSE <- as.data.frame(mcse(mod_bioregion))

write_xlsx(as.data.frame(mcse(mod_bioregion)), "./Desktop/MCSE/BioMCSE.xlsx")

mcse(sp.mod_sp)

write_xlsx(as.data.frame(mcse(sp.mod_sp)), "./Desktop/MCSE/SpeciesMCSE.xlsx")

mcse(sp.mod_bioregion_sp_int)

write_xlsx(as.data.frame(mcse(sp.mod_bioregion_sp_int)), "./Desktop/MCSE/SpeciesInteractiveMCSE.xlsx")

mcse(sp.mod_size)

write_xlsx(as.data.frame(mcse(sp.mod_size)), "./Desktop/MCSE/SizeMCSE.xlsx")

mcse(sp.mod_bioregion_size_int)

write_xlsx(as.data.frame(mcse(sp.mod_bioregion_size_int)), "./Desktop/MCSE/SizeInteractiveMCSE.xlsx")

mcse(td.mod_td)

write_xlsx(as.data.frame(mcse(td.mod_td)), "./Desktop/MCSE/TidalClassMCSE.xlsx")

mcse(td.mod_td_bio_int)

write_xlsx(as.data.frame(mcse(td.mod_td_bio_int)), "./Desktop/MCSE/TidalClassInteractiveMCSE.xlsx")

mcse(geo.mod_geo)

write_xlsx(as.data.frame(mcse(geo.mod_geo)), "./Desktop/MCSE/GeomorphMCSE.xlsx")

mcse(geo.mod_geo_bio)

write_xlsx(as.data.frame(mcse(geo.mod_geo_bio)), "./Desktop/MCSE/GeomorphInteractiveMCSE.xlsx")







#### Calculating the parameter estimates

#### BIOREGION CI DATASET

post.bioregion <- as.data.frame(as_draws_matrix(mod_bioregion)) ## posterior_samples creates dataset of draws
post.bioregion <- post.bioregion[,c(1:6)] # Cut out everything but fixed effects columns (optional)

names(post.bioregion) <- c("Mediterranean", "Temperate North Atlantic", "Temperate North Pacific", "Temperate Southern Ocean", "Tropical Atlantic", "Tropical Indo-Pacific")
## Gather organizes the data in such a way that you can use ggplot intuitively
post.bioregion <- gather(post.bioregion, key = "Bioregion", value = "Est", 1:6 )

hdi_95 <- post.bioregion %>% 
  group_by(Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.95)) %>% 
  select(-c(CI)) %>% 
  rename(hdi5 = CI_low, hdi95 = CI_high)

hdi_90 <- post.bioregion %>% 
  group_by(Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.90)) %>% 
  select(-c(CI)) %>% 
  rename(hdi10 = CI_low, hdi90 = CI_high)

hdi_80 <- post.bioregion %>% 
  group_by(Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.80)) %>% 
  select(-c(CI)) %>% 
  rename(hdi20 = CI_low, hdi80 = CI_high)

mean <- post.bioregion %>% group_by(Bioregion) %>% summarize(mean = mean(Est))

bioregion_summary <- merge(merge(merge(hdi_80, hdi_90), hdi_95), mean)

##

post.global <- as.data.frame(as_draws_matrix(mod_bioregion_null))

post.global <- as.data.frame(post.global[,1])

names(post.global) <- "Global"
post.global <- gather(post.global, key = "Bioregion", value = "Est", 1)

hdi_95 <- post.global %>% 
  bayestestR::hdi(Bioregion, ci = 0.95, verbose = TRUE) %>% 
  select(-c(CI)) %>%
  rename(hdi5 = CI_low, hdi95 = CI_high)

hdi_90 <- post.global %>% 
  bayestestR::hdi(Bioregion, ci = 0.90, verbose = TRUE) %>% 
  select(-c(CI)) %>%
  rename(hdi10 = CI_low, hdi90 = CI_high)

hdi_80 <- post.global %>% 
  bayestestR::hdi(Bioregion, ci = 0.80, verbose = TRUE) %>% 
  select(-c(CI)) %>%
  rename(hdi20 = CI_low, hdi80 = CI_high)

mean <- post.global %>% summarize(mean = mean(Est))

global_summary <- merge(merge(merge(hdi_80, hdi_90), hdi_95), mean) %>% rename('Bioregion' = 'Parameter')

global_summary[,1] <- 'Global'

bioregion_global_summary <- rbind(global_summary, bioregion_summary)

bioregion_global_summary$N <- c(1230, 57, 233, 82, 488, 48, 322)

bioregion_global_summary$cores <- c(2174, 69, 378, 173, 628, 130, 796)

bioregion_global_summary$sig <- c(NA, 'a', 'bc', 'bd', 'be', 'ce', 'de')


write_xlsx(bioregion_global_summary, './Desktop/HDIntervalDataframes/BioregionCI.xlsx')







#### SPECIES X BIOREGION CI DATAFRAME

post.bioregion_sp <- as_draws_matrix(sp.mod_bioregion_sp_int)
post.bioregion_sp <- as.data.frame(post.bioregion_sp[,1:140])

post.bioregion_sp <- gather(post.bioregion_sp, key = "var", value = "Est", 1:ncol(post.bioregion_sp))
post.bioregion_sp <- post.bioregion_sp %>% 
  separate(col = "var", into = c("Bioregion", "Species"), sep = ":")

post.bioregion_sp_sub <- subset(
  x = post.bioregion_sp,
  Bioregion == "b_BioregionMediterranean" & Species == "SpeciesPosidoniaoceanica" | 
    Bioregion == "b_BioregionMediterranean" & Species == "SpeciesZosteramarina" | 
    Bioregion == "b_BioregionMediterranean" & Species == "SpeciesCymodoceanodosa" | 
    Bioregion == "b_BioregionMediterranean" & Species == "SpeciesZosteranoltii" |
    Bioregion == "b_BioregionTemperateNorthAtlantic" & Species == "SpeciesZosteranoltii" |
    Bioregion == "b_BioregionTemperateSouthernOcean" & Species == "SpeciesAmphibolisantarctica" | 
    Bioregion == "b_BioregionTemperateSouthernOcean" & Species == "SpeciesHalophilaovalis" | 
    Bioregion == "b_BioregionTemperateSouthernOcean" & Species == "SpeciesPosidoniaaustralis" | 
    Bioregion == "b_BioregionTemperateSouthernOcean" & Species == "SpeciesPosidoniasinuosa" | 
    Bioregion == "b_BioregionTropicalAtlantic" & Species == "SpeciesThalassiatestudinum" | 
    Bioregion == "b_BioregionTropicalAtlantic" & Species == "SpeciesHalodulewrightii" | 
    Bioregion == "b_BioregionTropicalIndoMPacific" & Species == "SpeciesCymodoceaserrulata" |
    Bioregion == "b_BioregionTropicalIndoMPacific" & Species == "SpeciesEnhalusacoroides" |
    Bioregion == "b_BioregionTropicalIndoMPacific" & Species == "SpeciesHaloduleuninervis" |
    Bioregion == "b_BioregionTropicalIndoMPacific" & Species == "SpeciesHalophiladecipiens" |
    Bioregion == "b_BioregionTropicalIndoMPacific" & Species == "SpeciesHalophilaovalis" |
    Bioregion == "b_BioregionTropicalIndoMPacific" & Species == "SpeciesHalophilastipulacea" |
    Bioregion == "b_BioregionTropicalIndoMPacific" & Species == "SpeciesSyringodiumisoetifolium" |
    Bioregion == "b_BioregionTropicalIndoMPacific" & Species == "SpeciesThalassiahemprichii" |
    Bioregion == "b_BioregionTropicalIndoMPacific" & Species == "SpeciesThalassodendronciliatum" |
    Bioregion == "b_BioregionTropicalIndoMPacific" & Species == "SpeciesZosteramuelleri" | 
    Bioregion == "b_BioregionTemperateNorthPacific" & Species == "SpeciesZosteramarina" |
    Bioregion == "b_BioregionTemperateNorthAtlantic" & Species == "SpeciesZosteramarina",
  select = c(Bioregion, Species, Est)
)

hdi_95 <- post.bioregion_sp_sub %>% 
  group_by(Species, Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.95)) %>% 
  select(-c(CI)) %>% 
  rename(hdi5 = CI_low, hdi95 = CI_high)

hdi_90 <- post.bioregion_sp_sub %>% 
  group_by(Species, Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.90)) %>% 
  select(-c(CI)) %>% 
  rename(hdi10 = CI_low, hdi90 = CI_high)

hdi_80 <- post.bioregion_sp_sub %>% 
  group_by(Species, Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.80)) %>% 
  select(-c(CI)) %>% 
  rename(hdi20 = CI_low, hdi80 = CI_high)

mean <- post.bioregion_sp_sub %>% group_by(Species, Bioregion) %>% summarize(mean = mean(Est))

bioxsp_summary <- merge(merge(merge(hdi_80, hdi_90), hdi_95), mean)

bioxsp_summary$sig <- c('cde', 
                        'b', 
                        'cde',
                        'ce', 
                        'cde',
                        'de',
                        'cde', 'cde', 'cde', 'cde', 'cde',
                        'a',
                        'cde',
                        'd',
                        'cde', 
                        'bcd', 
                        'cde',
                        'cde', 'cde', 'cde',
                        'bc', 
                        'cde', 'cde'
)
bioxsp_summary$cores <- c(59, 18, 20, 77, 57, 57,
                          4,
                          6, 24,
                          10, 52, 30, 47, 15, 36, 22, 23,
                          12, 372, 173,
                          33,
                          9, 6)

bioxsp_summary$N <- c(4, 6, 6, 33, 20, 17, 
                      4, 6, 9,
                      10, 8, 30, 6, 15, 35, 16, 23,
                      12, 227, 82,
                      4,
                      9, 6)


write_xlsx(bioxsp_summary, "./Desktop/HDIntervalDataframes/SpeciesXBioregionIntervals.xlsx")





## Tidal Depth

post.tidalglobal <- as.data.frame(as_draws_matrix(td.mod_td))
post.tidalglobal <- as.data.frame(post.tidalglobal[,1:2])

post.tidalglobal <- gather(post.tidalglobal, key = "TidalClass", value = "Est", 1:2 )

hdi_95 <- post.tidalglobal %>% 
  group_by(TidalClass) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.95)) %>% 
  select(-c(CI)) %>% 
  rename(hdi5 = CI_low, hdi95 = CI_high)

hdi_90 <- post.tidalglobal %>% 
  group_by(TidalClass) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.90)) %>% 
  select(-c(CI)) %>% 
  rename(hdi10 = CI_low, hdi90 = CI_high)

hdi_80 <- post.tidalglobal %>% 
  group_by(TidalClass) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.80)) %>% 
  select(-c(CI)) %>% 
  rename(hdi20 = CI_low, hdi80 = CI_high)

mean <- post.tidalglobal %>% group_by(TidalClass) %>% summarize(mean = mean(Est))

tidalglobal_summary <- merge(merge(merge(hdi_80, hdi_90), hdi_95), mean)

tidalglobal_summary$N <- c(99, 302)

tidalglobal_summary$cores <- c(245, 730)

write_xlsx(tidalglobal_summary, './Desktop/HDIntervalDataframes/GlobalTidalCI_Summary.xlsx')




## Tidal, Bioregional

post.bioregion_td <- as.data.frame(as_draws_matrix(td.mod_td_bio_int))
post.bioregion_td <- post.bioregion_td[,1:6]

post.bioregion_td <- gather(post.bioregion_td, key = "var", value = "Est", 1:ncol(post.bioregion_td))

post.bioregion_td <- post.bioregion_td %>% 
  separate(col = "var", into = c("Bioregion", "TidalClass"), sep = ":")

hdi_95 <- post.bioregion_td %>% 
  group_by(TidalClass, Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.95)) %>% 
  select(-c(CI)) %>% 
  rename(hdi5 = CI_low, hdi95 = CI_high)

hdi_90 <- post.bioregion_td %>% 
  group_by(TidalClass, Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.90)) %>% 
  select(-c(CI)) %>% 
  rename(hdi10 = CI_low, hdi90 = CI_high)

hdi_80 <- post.bioregion_td %>% 
  group_by(TidalClass, Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.80)) %>% 
  select(-c(CI)) %>% 
  rename(hdi20 = CI_low, hdi80 = CI_high)

mean <- post.bioregion_td %>% group_by(TidalClass, Bioregion) %>% summarize(mean = mean(Est))

tidaldepth_bioregion_summary <- merge(merge(merge(hdi_80, hdi_90), hdi_95), mean)

tidaldepth_bioregion_summary$N <- c(9, 54, 31, 209, 25, 5)

tidaldepth_bioregion_summary$cores <- c(9, 89, 115, 354, 81, 61)

write_xlsx(tidaldepth_bioregion_summary, './Desktop/HDIntervalDataframes/TidalXBioregionCI_Summary.xlsx')








### Hydrologic setting, global


post.geoglobal <- as.data.frame(as_draws_matrix(geo.mod_geo))
post.geoglobal <- as.data.frame(post.geoglobal[,1:2])

post.geoglobal <- gather(post.geoglobal, key = "Geomorphology", value = "Est", 1:2 )

hdi_95 <- post.geoglobal %>% 
  group_by(Geomorphology) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.95)) %>% 
  select(-c(CI)) %>% 
  rename(hdi5 = CI_low, hdi95 = CI_high)

hdi_90 <- post.geoglobal %>% 
  group_by(Geomorphology) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.90)) %>% 
  select(-c(CI)) %>% 
  rename(hdi10 = CI_low, hdi90 = CI_high)

hdi_80 <- post.geoglobal %>% 
  group_by(Geomorphology) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.80)) %>% 
  select(-c(CI)) %>% 
  rename(hdi20 = CI_low, hdi80 = CI_high)

mean <- post.geoglobal %>% group_by(Geomorphology) %>% summarize(mean = mean(Est))

geoglobal_summary <- merge(merge(merge(hdi_80, hdi_90), hdi_95), mean)

geoglobal_summary$N <- c(249, 63)

geoglobal_summary$cores <- c(605, 237)

write_xlsx(geoglobal_summary, './Desktop/HDIntervalDataframes/GlobalGeomorphCI_Summary.xlsx')


## Hydrologic setting, Bioregional

post.bioregion_geo <- as.data.frame(as_draws_matrix(geo.mod_geo_bio))
post.bioregion_geo <- post.bioregion_geo[,1:6]

post.bioregion_geo <- gather(post.bioregion_geo, key = "var", value = "Est", 1:ncol(post.bioregion_geo))

post.bioregion_geo <- post.bioregion_geo %>% 
  separate(col = "var", into = c("Bioregion", "Geomorphology"), sep = ":")

hdi_95 <- post.bioregion_geo %>% 
  group_by(Geomorphology, Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.95)) %>% 
  select(-c(CI)) %>% 
  rename(hdi5 = CI_low, hdi95 = CI_high)

hdi_90 <- post.bioregion_geo %>% 
  group_by(Geomorphology, Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.90)) %>% 
  select(-c(CI)) %>% 
  rename(hdi10 = CI_low, hdi90 = CI_high)

hdi_80 <- post.bioregion_geo %>% 
  group_by(Geomorphology, Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.80)) %>% 
  select(-c(CI)) %>% 
  rename(hdi20 = CI_low, hdi80 = CI_high)

mean <- post.bioregion_geo %>% group_by(Geomorphology, Bioregion) %>% summarize(mean = mean(Est))

geomorph_bioregion_summary <- merge(merge(merge(hdi_80, hdi_90), hdi_95), mean)

geomorph_bioregion_summary$N <- c(26, 14, 126, 21, 9, 28)

geomorph_bioregion_summary$cores <- c(76, 40, 176, 62, 29, 118)

write_xlsx(geomorph_bioregion_summary, './Desktop/HDIntervalDataframes/GeomorphXBioregionCI_Summary.xlsx')





### Species Size, Bioregional

post.bioregion_size <- as.data.frame(as_draws_matrix(sp.mod_bioregion_size_int))
post.bioregion_size <- post.bioregion_size[,1:8]

post.bioregion_size <- gather(post.bioregion_size, key = "var", value = "Est", 1:ncol(post.bioregion_size))

post.bioregion_size <- post.bioregion_size %>% 
  separate(col = "var", into = c("Bioregion", "SpeciesSize"), sep = ":")

hdi_95 <- post.bioregion_size %>% 
  group_by(SpeciesSize, Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.95)) %>% 
  select(-c(CI)) %>% 
  rename(hdi5 = CI_low, hdi95 = CI_high)

hdi_90 <- post.bioregion_size %>% 
  group_by(SpeciesSize, Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.90)) %>% 
  select(-c(CI)) %>% 
  rename(hdi10 = CI_low, hdi90 = CI_high)

hdi_80 <- post.bioregion_size %>% 
  group_by(SpeciesSize, Bioregion) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.80)) %>% 
  select(-c(CI)) %>% 
  rename(hdi20 = CI_low, hdi80 = CI_high)

mean <- post.bioregion_size %>% group_by(SpeciesSize, Bioregion) %>% summarize(mean = mean(Est))

size_bioregion_summary <- merge(merge(merge(hdi_80, hdi_90), hdi_95), mean)

size_bioregion_summary$N <- c(30, 18, 16, 91, 27, 6, 17, 68)

size_bioregion_summary$cores <- c(30, 158, 22, 136, 39, 6, 57, 163)

write_xlsx(size_bioregion_summary, './Desktop/HDIntervalDataframes/SizeXBioregionCI_Summary.xlsx')



## Species Size, global


post.sizeglobal <- as.data.frame(as_draws_matrix(sp.mod_size))
post.sizeglobal <- as.data.frame(post.sizeglobal[,1:2])

post.sizeglobal <- gather(post.sizeglobal, key = "SpeciesSize", value = "Est", 1:2 )

hdi_95 <- post.sizeglobal %>% 
  group_by(SpeciesSize) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.95)) %>% 
  select(-c(CI)) %>% 
  rename(hdi5 = CI_low, hdi95 = CI_high)

hdi_90 <- post.sizeglobal %>% 
  group_by(SpeciesSize) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.90)) %>% 
  select(-c(CI)) %>% 
  rename(hdi10 = CI_low, hdi90 = CI_high)

hdi_80 <- post.sizeglobal %>% 
  group_by(SpeciesSize) %>% 
  summarize(bayestestR::hdi(Est, ci = 0.80)) %>% 
  select(-c(CI)) %>% 
  rename(hdi20 = CI_low, hdi80 = CI_high)

mean <- post.sizeglobal %>% group_by(SpeciesSize) %>% summarize(mean = mean(Est))

sizeglobal_summary <- merge(merge(merge(hdi_80, hdi_90), hdi_95), mean)

sizeglobal_summary$N <- c(155, 433)

sizeglobal_summary$cores <- c(346, 816)

write_xlsx(sizeglobal_summary, './Desktop/HDIntervalDataframes/GlobalSizeCI_Summary.xlsx')





















