# soil data analysis script
# last updated 20180229
# contact: sheila saia (sms403@cornell.edu)

# ---- 1 load libraries ----

# clear ws
rm(list = ls())

# load packages
library(tidyverse)
#library(car)


# ---- 2 load data ----

# load icp data
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/icp_data/")
nahco3_tp_data = read_csv("icp_nahco3_calibrated.csv")
naoh_tp_data = read_csv("icp_naoh_calibrated.csv")
hcl_tp_data = read_csv("icp_hcl_calibrated.csv")
ox_tp_data = read_csv("icp_ox_calibrated.csv")
totalp_tp_data = read_csv("icp_totalp_calibrated.csv")

# load malachite data
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/malachite_data/")
nahco3_dp_data = read_csv("malachite_nahco3_calibrated.csv")
naoh_dp_data = read_csv("malachite_naoh_calibrated.csv")
hcl_dp_data = read_csv("malachite_hcl_calibrated.csv")
ox_dp_data = read_csv("malachite_ox_calibrated.csv")
cacl2_dp_data = read_csv("molybdenum_cacl2_calibrated.csv")

# load isotope data
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/isotope_data/")
isotope_data = read_csv("isotope_data_reformatted_20171201.csv")

# load bonn sample, site info, and sample metadata
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/site_and_exp_design_data/")
bonn_sample_info = read_csv("bonn_sample_info.csv")
# site_info = read_table2("siteInfo.txt") # has lat/long
metadata = read_csv("all_joined_data.csv") # no lat/long


# ---- 3 reformat data ----

# merge all p pool data (nahco3, naoh, hcl, and ox)
p_pool_data = rbind(nahco3_tp_data, nahco3_dp_data,
                    naoh_tp_data, naoh_dp_data,
                    hcl_tp_data, hcl_dp_data, 
                    ox_tp_data, ox_dp_data,
                    totalp_tp_data, cacl2_dp_data) %>%
  select(-row_num) %>%
  mutate(row_num = seq(1:590))

# summarize select environmental metadata
metadata_sel = metadata %>% 
  select(SampleID, State, STI, texture, VMCCalibperc, pH, avgTempC, OMperc) %>%
  group_by(SampleID) %>% 
  summarize(State = unique(State), STI = unique(STI), texture = unique(texture), 
            avg_vmc_perc = mean(VMCCalibperc),
            avg_pH = mean(pH),
            avg_temp_c = mean(avgTempC),
            avg_om_perc = mean(OMperc))

# bonn sample info select
bonn_sample_info_sel = bonn_sample_info %>%
  filter(Notes == "pooled ox sample A through C") %>%
  select(-Rep, -Notes)
  
# reformat isotope data and add in metadata
isotope_data_sel = isotope_data %>% 
  select(run_number, Bonn_SampleID:d18O) %>%
  filter(Bonn_SampleID != "standard") %>%
  left_join(bonn_sample_info_sel, by = "Bonn_SampleID") %>%
  unique() %>% # for some reason joining doubles the length so use this to remove doubles
  left_join(metadata_sel, by = "SampleID")

length(isotope_data_sel$State[isotope_data_sel$State == "PA"]) # 21 samples
length(isotope_data_sel$State[isotope_data_sel$State == "NY"]) # 32 samples b/c some got re-analyzed

# summarize p_pool_data and isotope data and combine
# hcl_isotope_data = isotope_data_sel %>%
#   filter(Extraction == "hcl") %>%
#   select(SampleID, d18O, SiteID:avg_om_perc, Extraction
#   
# 
# ox_isotope_data = isotope_data_sel %>%
#   filter(Extraction == "ox")  


# ---- 4.1 stats for isotope vs extraction ----

# check for data normality
hist(isotope_data_sel$d18O)
qqPlot(isotope_data_sel$d18O)
# data is within bounds therefore is normally distributed

# anlysis of variance for differneces between extractions
iso_vs_ext_anova = aov(d18O ~ Extraction, data = isotope_data_sel)
summary(iso_vs_ext_anova)
# p-value > 0.5 so cannot reject null-hypothesis (= extractions are the same)

# analysis of variance for differences between extractions from different states$
iso_vs_ext_by_state_anova = aov(d18O ~ Extraction * State, data = isotope_data_sel)
summary(iso_vs_ext_by_state_anova)
# p-value for state is sign. so d18O are different between states

# ---- 4.2 plot isotope data by extraction ----

# plot d18O value vs extraction (both states combined by points overlayed and colored by state)
ggplot(isotope_data_sel,aes(x=Extraction,y=d18O)) +
  geom_boxplot(width = 0.5) +
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



# ---- 4.3 stats for isotope vs texture ----

unique(isotope_data_sel$texture[isotope_data_sel$State == "PA"])
unique(isotope_data_sel$texture[isotope_data_sel$State == "NY"])
# both have silt loam in common

# select only silt loam
isotope_data_sel_siltloam = isotope_data_sel %>%
  filter(texture == "silt_loam")

# check for data normality
hist(isotope_data_sel_siltloam$d18O)
qqPlot(isotope_data_sel_siltloam$d18O)
# data is within bounds therefore is normally distributed

# anlysis of variance for differneces between extractions
iso_vs_text_anova = aov(d18O ~ Extraction, data = isotope_data_sel_siltloam)
summary(iso_vs_text_anova)
# p-value > 0.5 so cannot reject null-hypothesis (= extractions are the same)

# analysis of variance for differences between extractions from different states$
iso_vs_text_by_state_anova = aov(d18O ~ Extraction * State, data = isotope_data_sel_siltloam)
summary(iso_vs_text_by_state_anova)
# p-value for state and extraction are sign. so extractions are different between states


# ---- 4.3 plot isotope data by texture ----

# plot d18O value vs extraction (both states combined by points overlayed and colored by state)
ggplot(isotope_data_sel,aes(x=Extraction,y=d18O)) +
  geom_point(size = 3, alpha = 0.5, aes(color = State, shape = texture)) +
  xlab("Extraction") +
  ylab("d18O") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# plot d18O value vs extraction (both states combined by points overlayed and colored by state)
ggplot(isotope_data_sel_siltloam,aes(x=Extraction,y=d18O)) +
  geom_boxplot(width = 0.5) +
  geom_point(size = 5, alpha = 0.5, aes(color = State, shape = texture)) +
  xlab("Extraction") +
  ylab("d18O") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  scale_shape_manual(values = c(18)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# ---- 4.4 stats for isotope data vs vmc ----

# check for data normality
hist(isotope_data_sel$d18O)
qqPlot(isotope_data_sel$d18O)
# data is within bounds therefore is normally distributed

# anlysis of variance for differneces between extractions
iso_vs_vmc_lm = lm(d18O ~ avg_vmc_perc, data = isotope_data_sel)
summary(iso_vs_vmc_lm)
# p-value < 0.5 for slope so there is a significant negative trend

# analysis of variance for differences between extractions from different states$
iso_vs_vmc_by_state_lm = lm(d18O ~ avg_vmc_perc * State, data = isotope_data_sel)
summary(iso_vs_vmc_by_state_lm)
# p-value for state and extraction are sign. so extractions are different between states


# ---- 4.4 plot isotope data vs vmc ----

# all data
ggplot(isotope_data_sel,aes(x=avg_vmc_perc,y=d18O, color = State)) +
  geom_point(size = 3, alpha = 0.5) +
  xlab("VMC (%)") +
  ylab("d18O") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  scale_shape_manual(values = c(18)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# just silt loams
ggplot(isotope_data_sel_siltloam,aes(x=avg_vmc_perc,y=d18O, color = State, shape = texture)) +
  geom_point(size = 3, alpha = 0.5) +
  xlab("VMC (%)") +
  ylab("d18O") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  scale_shape_manual(values = c(18)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# ---- 5.1 average calcs ----

# isotope summary data
isotope_data_sel_summary = isotope_data_sel %>%
  group_by(SampleID, Extraction) %>%
  summarise(SiteID = unique(SiteID),
            avg_d18O = mean(d18O, na.rm = TRUE),
            sd_d18O = sd(d18O, na.rm = TRUE),
            State = unique(State),
            STI = unique(STI),
            texture = unique(texture),
            avg_vmc_perc = unique(avg_vmc_perc),
            avg_pH = unique(avg_pH),
            avg_temp_c = unique(avg_temp_c),
            avg_om_perc = unique(avg_om_perc)) %>%
  na.omit() %>% # remove PA26 because only have ox results
  mutate(SampleIDext = paste0(SampleID,"_",Extraction))

# p pool summary data
p_pool_data_summary = p_pool_data %>%
  group_by(SampleID, Extraction, P_pool) %>%
  summarize(avg_P_mgperkg = mean(P_mgperkg, na.rm = TRUE),
            sd_P_mgperkg = sd(P_mgperkg, na.rm = TRUE)) %>%#,
            #avg_P_umol = mean(P_umol, na.rm = TRUE),
            #sd_P_umol = sd(P_umol, na.rm = TRUE)) %>%
  filter(SampleID != "samp_PA26_20150715") # remove because have no d18O data

# join isotope data and p pool data for hcl and ox samples
p_pool_isotope_data_summary = p_pool_data_summary %>%
  filter(Extraction == "hcl" | Extraction == "ox") %>%
  mutate(SampleIDext = paste0(SampleID,"_",Extraction)) %>%
  left_join(isotope_data_sel_summary, by = "SampleIDext") %>%
  select(SampleID = SampleID.x, Extraction = Extraction.x,
         P_pool:SampleIDext,SiteID:avg_om_perc) %>% # delete duplicate columns
  select(SampleIDext, SampleID, SiteID, State:avg_om_perc, 
         Extraction, P_pool:sd_P_mgperkg, avg_d18O, sd_d18O) # reorganized

# hedley data
hedley_p_pool_data = p_pool_data_summary %>%
  filter(Extraction == "hcl" | Extraction == "naoh" | Extraction == "nahco3")



# merge all hcl and ox p pool data into separate dataframes to go with isotope data
hcl_p_pool_data = rbind(hcl_tp_data, hcl_dp_data)
ox_p_pool_data = rbind(ox_tp_data, ox_dp_data)


# separate out DP and TP data
p_pool_isotope_dp_data_summary = p_pool_isotope_data_summary %>%
  filter(P_pool == "DP")
p_pool_isotope_tp_data_summary = p_pool_isotope_data_summary %>%
  filter(P_pool == "TP")

# plot dp data
ggplot(p_pool_isotope_dp_data_summary, aes(x = avg_P_mgperkg, y = avg_d18O, color = State, shape = Extraction)) +
  geom_point(size = 3, alpha = 0.75) +
  xlab("DP (mg/kg)") +
  ylab("d18O") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  #scale_shape_manual(values = c(18)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# plot tp data
ggplot(p_pool_isotope_tp_data_summary, aes(x = avg_P_mgperkg, y = avg_d18O, color = State, shape = Extraction)) +
  geom_point(size = 3, alpha = 0.75) +
  xlab("TP (mg/kg)") +
  ylab("d18O") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  #scale_shape_manual(values = c(18)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
