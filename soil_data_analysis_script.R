# soil data analysis script
# last updated 20180229
# contact: sheila saia (sms403@cornell.edu)

# ---- 1 load libraries ----

# clear ws
rm(list = ls())

# load packages
library(tidyverse)


# ---- 2 load data ----

# load icp data
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/icp_data/")
nahco3_icp_data = read_csv("icp_nahco3_calibrated.csv")
naoh_icp_data = read_csv("icp_naoh_calibrated.csv")
hcl_icp_data = read_csv("icp_hcl_calibrated.csv")
ox_icp_data = read_csv("icp_ox_calibrated.csv")

# load malachite data
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/malachite_data/")
nahco3_malachite_data = read_csv("malachite_nahco3_calibrated.csv")
naoh_malachite_data = read_csv("malachite_naoh_calibrated.csv")
hcl_malachite_data = read_csv("malachite_hcl_calibrated.csv")
ox_malachite_data = read_csv("malachite_ox_calibrated.csv")

# load isotope data
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/isotope_data/")
isotope_data = read_csv("isotope_data_reformatted_20171201.csv")

# load bonn sample, site info, and sample metadata
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/site_and_exp_design_data/")
bonn_sample_info = read_table2("bonn_sample_info.txt")
# site_info = read_table2("siteInfo.txt") # has lat/long
metadata = read_csv("all_joined_data.csv") # no lat/long


# ---- 3 reformat data ----

# merge all p pool data (nahco3, naoh, hcl, and ox)
p_pool_data = rbind(nahco3_icp_data, nahco3_malachite_data,
                    naoh_icp_data, naoh_malachite_data,
                    hcl_icp_data, hcl_malachite_data, 
                    ox_icp_data, ox_malachite_data) %>%
  select(-row_num) %>%
  mutate(row_num = seq(1:470))

# merge all hcl and ox p pool data into separate dataframes to go with isotope data
hcl_p_pool_data = rbind(hcl_icp_data, hcl_malachite_data)
ox_p_pool_data = rbind(ox_icp_data, ox_malachite_data)

# summarize select environmental metadata
metadata_sel = metadata %>% 
  select(SampleID, State, STI, texture, VMCCalibperc, pH, avgTempC, OMperc) %>%
  group_by(SampleID) %>% 
  summarize(State = unique(State), STI = unique(STI), texture = unique(texture), 
            avg_vmc_perc = mean(VMCCalibperc),
            avg_pH = mean(pH),
            avg_temp_c = mean(avgTempC),
            avg_om_perc = mean(OMperc))

# reformat isotope data and add in metadata
isotope_data_sel = isotope_data %>% 
  select(run_number, Bonn_SampleID:d18O) %>%
  filter(Extraction != "NA") %>%
  left_join(bonn_sample_info, by = "Bonn_SampleID") %>%
  select(-Rep, -Notes) %>% # pooled so Rep is not needed
  # left_join(site_info_sel, by = "SiteID") %>%
  left_join(metadata_sel, by = "SampleID")

# summarize p_pool_data and isotope data and combine
# hcl_isotope_data = isotope_data_sel %>%
#   filter(Extraction == "hcl") %>%
#   select(SampleID, d18O, SiteID:avg_om_perc, Extraction
#   
# 
# ox_isotope_data = isotope_data_sel %>%
#   filter(Extraction == "ox")  


# ---- 4 plot isotope data by extraction ----

# plot d18O value vs extraction (both states combined by points overlayed and colored by state)
ggplot(isotope_data_sel,aes(x=Extraction,y=d18O)) +
  geom_boxplot() +
  geom_point(size = 3, alpha = 0.5, aes(color = State)) +
  xlab("Extraction") +
  ylab("d18O") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# plot d18O value vs extraction by state
ggplot(isotope_data_sel,aes(x=Extraction,y=d18O,fill=State)) +
  geom_boxplot() +
  xlab("Extraction") +
  ylab("d18O") +
  scale_fill_manual(values = c("#f1a340", "#998ec3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
