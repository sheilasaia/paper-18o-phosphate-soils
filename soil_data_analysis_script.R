# soil data analysis script
# last updated 20180229
# contact: sheila saia (sms403@cornell.edu)

# ---- 1 load libraries and funcitons ----

# clear ws
rm(list = ls())

# load packages
library(tidyverse)
#library(car)

# load home-made functions 
functions_path="/Users/ssaia/Documents/GitHub/yadkin-analysis/functions/"
source(paste0(functions_path,"multiplot.R")) # for creating plots with multiple figures


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

# ---- 5.1 summarize isotope and p pool data ----

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
  group_by(SampleID, Extraction, Method, P_pool) %>%
  summarize(avg_P_mgperkg = mean(P_mgperkg, na.rm = TRUE),
            sd_P_mgperkg = sd(P_mgperkg, na.rm = TRUE)) %>%#,
            #avg_P_umol = mean(P_umol, na.rm = TRUE),
            #sd_P_umol = sd(P_umol, na.rm = TRUE)) %>%
  filter(SampleID != "samp_PA26_20150715") # remove because have no d18O data


# ---- 5.2 join summary calcs ----

# join isotope and p pool data for hcl and ox extractions
hcl_ox_p_pool_isotope_data_summary = p_pool_data_summary %>%
  filter(Extraction == "hcl" | Extraction == "ox") %>%
  mutate(SampleIDext = paste0(SampleID,"_",Extraction)) %>%
  left_join(isotope_data_sel_summary, by = "SampleIDext") %>%
  select(SampleID = SampleID.x, Extraction = Extraction.x, Method,
         P_pool:SampleIDext,SiteID:avg_om_perc) %>% # delete duplicate columns
  select(SampleIDext, SampleID, SiteID, State:avg_om_perc,
         Extraction, Method, P_pool:sd_P_mgperkg, avg_d18O, sd_d18O) # reorganized

# join isotope and p pool data for cacl2 extractions
cacl2_p_pool_isotope_data_summary = p_pool_data_summary %>%
  filter(Extraction == "cacl2") %>%
  left_join(isotope_data_sel_summary, by = "SampleID") %>%
  select(SampleID = SampleID.x, Extraction = Extraction.x, Method,
         P_pool:SampleIDext,SiteID:avg_om_perc) %>% # delete duplicate columns
  select(SampleIDext, SampleID, SiteID, State:avg_om_perc,
         Extraction, Method, P_pool:sd_P_mgperkg, avg_d18O, sd_d18O) # reorganized


# ---- 5.3 plot d18O vs hcl and ox p pools by state ----

# select only dp data
hcl_ox_p_pool_isotope_data_summary_sel = p_pool_isotope_data_summary %>%
  filter(Extraction == "hcl" | Extraction == "ox") %>%
  filter(Method == "malachite") # until i figure out what's going on with bonn vs cornell icp ox data

# plot
ggplot(hcl_ox_p_pool_isotope_data_summary_sel, aes(x = avg_P_mgperkg, y = avg_d18O, color = State, shape = Extraction)) +
  geom_point(size = 3, alpha = 0.75) +
  xlab("Mean DP (mg/kg)") +
  ylab("Mean d18O") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  #scale_shape_manual(values = c(18)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# select only silt loams
hcl_ox_p_pool_isotope_data_summary_sel2 = hcl_ox_p_pool_isotope_data_summary_sel %>%
  filter(texture == "silt_loam")

# plot
ggplot(hcl_ox_p_pool_isotope_data_summary_sel2, aes(x = avg_P_mgperkg, y = avg_d18O, color = State, shape = Extraction)) +
  geom_point(size = 3, alpha = 0.75) +
  xlab("Mean DP Silt Loams (mg/kg)") +
  ylab("Mean d18O") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  #scale_shape_manual(values = c(18)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


# ---- 5.4 plot cacl2 pools vs isotope data ----

ggplot(cacl2_p_pool_isotope_data_summary, aes(x = avg_P_mgperkg, y = avg_d18O, color = State, shape = Extraction.y)) +
  geom_point(size = 3, alpha = 0.75) +
  xlab("Mean CaCl2-DP (mg/kg)") +
  ylab("Mean d18O") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  #scale_shape_manual(values = c(18)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


# ---- 6.1 join p pool and metadata select ----

p_pool_metadata_join = p_pool_data %>%
  left_join(metadata_sel, by = "SampleID") %>%
  select(-Bonn_SampleID, -row_num)

# hedley dp
hedley_dp_p_pool_metadata_join = p_pool_metadata_join %>%
  filter(Extraction == "nahco3" | Extraction == "naoh" | Extraction == "hcl") %>%
  filter(P_pool == "DP")
hedley_dp_p_pool_metadata_join$Extraction = factor(hedley_dp_p_pool_metadata_join$Extraction, levels = c("nahco3", "naoh", "hcl"))

# hedley tp
hedley_tp_p_pool_metadata_join = p_pool_metadata_join %>%
  filter(Extraction == "nahco3" | Extraction == "naoh" | Extraction == "hcl") %>%
  filter(P_pool == "TP")
hedley_tp_p_pool_metadata_join$Extraction = factor(hedley_tp_p_pool_metadata_join$Extraction, levels = c("nahco3", "naoh", "hcl"))

# cornell (all)
cornell_p_pool_metadata_join = p_pool_metadata_join %>%
  filter(Extraction == "cacl2" | Extraction == "ox" | Extraction == "totalp") %>%
  filter(Method != "malachite") # DP was from Bonn analysis
cornell_p_pool_metadata_join$Extraction = factor(cornell_p_pool_metadata_join$Extraction, levels = c("cacl2", "ox", "totalp"))

# bonn dp ox and cornell tp ox
ox_pool_metadata_join = p_pool_metadata_join %>%
  filter(Extraction == "ox")
ox_pool_metadata_join$Extraction = factor(ox_pool_metadata_join$Extraction, levels = c("cacl2", "ox", "totalp"))

# ---- 6.2 plot boxplots p pool vs extraction by state ----

# make plot object
my_p_vs_p_pool_by_state_plots = list()

# hedley dp
my_p_vs_p_pool_by_state_plots[[1]] = ggplot(hedley_dp_p_pool_metadata_join, aes(x = Extraction, y = P_mgperkg, color = State)) +
  geom_boxplot() +
  xlab("Extraction") +
  ylab("DP (mg/kg)") + 
  ylim(0, 950) +
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# hedley tp
my_p_vs_p_pool_by_state_plots[[2]] = ggplot(hedley_tp_p_pool_metadata_join, aes(x = Extraction, y = P_mgperkg, color = State)) +
  geom_boxplot() +
  xlab("Extraction") +
  ylab("TP (mg/kg)") + 
  ylim(0, 950) +
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# cornell (all)
my_p_vs_p_pool_by_state_plots[[3]] = ggplot(cornell_p_pool_metadata_join, aes(x = Extraction, y = P_mgperkg, color = State)) +
  geom_boxplot() +
  xlab("Extraction") +
  ylab("P (mg/kg)") + 
  ylim(0, 950) +
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

setwd("/Users/ssaia/Desktop")
cairo_pdf("p_vs_p_pool_by_state.pdf", width = 18, height = 6, pointsize = 18)
multiplot(plotlist = my_p_vs_p_pool_by_state_plots, cols = 3)
dev.off()

# ox comparision
ggplot(ox_pool_metadata_join, aes(x = P_pool, y = P_mgperkg, color = State)) +
  geom_boxplot() +
  xlab("P pool") +
  ylab("P (mg/kg)") + 
  ylim(0, 950) +
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


# ---- 6.3 plot dp and tp p pools (barchart) ----

hedley_p_pool_metadata_summary = p_pool_metadata_join %>% 
  filter(Extraction == "nahco3" | Extraction == "naoh" | Extraction == "hcl") %>%
  group_by(SampleID, Extraction, P_pool) %>%
  summarize(State = unique(State), STI = unique(STI), texture = unique(texture),
            avg_vmc_perc = unique(avg_vmc_perc), avg_pH = unique(avg_pH), avg_temp_c = unique(avg_temp_c), avg_om_perc = unique(avg_om_perc),
            avg_P_mgperkg = mean(P_mgperkg, na.rm = TRUE),
            sd_P_mgperkg = sd(P_mgperkg, na.rm = TRUE)) %>% #,
  #avg_P_umol = mean(P_umol, na.rm = TRUE),
  #sd_P_umol = sd(P_umol, na.rm = TRUE)) %>%
  filter(SampleID != "samp_PA26_20150715")  %>% # remove because have no d18O data
  mutate(Extraction_num = ifelse(Extraction == "hcl", 3,
                                 ifelse(Extraction == "naoh", 2, 1))) %>%
  mutate(SampleIDext = paste0(Extraction_num,"_",str_split(SampleID, "_")[[1]][2], "_", Extraction))

# set order of extractions
hedley_p_pool_metadata_summary$P_pool = factor(hedley_p_pool_metadata_summary$P_pool, levels = c("DP", "TP"))

# plot
ggplot(hedley_p_pool_metadata_summary, aes(x = reorder(SampleIDext, avg_vmc_perc), y = avg_P_mgperkg, fill = P_pool)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Sample_Extraction (ordered by increasing VMC)") +
  ylab("Mean P (mg/kg)") + 
  scale_fill_manual(values = c("grey75", "black")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))




# ---- 6.4 plots hedley pools vs environmental data ----

# make plot object
my_envir_plots = list()

# pH
my_envir_plots[[1]] = ggplot(hedley_p_pool_metadata_summary, aes(x = State, y = avg_pH, color = State)) +
  geom_boxplot(width = 0.5) +
  xlab("State") +
  ylab("Average pH") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# om
my_envir_plots[[2]] = ggplot(hedley_p_pool_metadata_summary, aes(x = State, y = avg_om_perc, color = State)) +
  geom_boxplot(width = 0.5) +
  xlab("State") +
  ylab("Average OM (%)") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# temp
my_envir_plots[[3]] = ggplot(hedley_p_pool_metadata_summary, aes(x = State, y = avg_temp_c, color = State)) +
  geom_boxplot(width = 0.5) +
  xlab("State") +
  ylab("Average Temperature (dec C)") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# vmc
my_envir_plots[[4]] = ggplot(hedley_p_pool_metadata_summary, aes(x = State, y = avg_vmc_perc, color = State)) +
  geom_boxplot(width = 0.5) +
  xlab("State") +
  ylab("Average VMC (%)") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

setwd("/Users/ssaia/Desktop")
cairo_pdf("environ_vars_by_state.pdf", width = 10, height = 10, pointsize = 18)
multiplot(plotlist = my_envir_plots, cols = 2)
dev.off()

# ---- extra ----

# merge all hcl and ox p pool data into separate dataframes to go with isotope data
# hcl_p_pool_data = rbind(hcl_tp_data, hcl_dp_data)
# ox_p_pool_data = rbind(ox_tp_data, ox_dp_data)


# separate out DP and TP data
p_pool_isotope_dp_data_summary = p_pool_isotope_data_summary %>%
  filter(P_pool == "DP")
p_pool_isotope_tp_data_summary = p_pool_isotope_data_summary %>%
  filter(P_pool == "TP")

# something is wierd weird with dp and tp pools...check that tp is more than dp

# dp data
my_iso_vs_p_plots = list()
my_iso_vs_p_plots[[1]] = ggplot(p_pool_isotope_dp_data_summary, aes(x = avg_P_mgperkg, y = avg_d18O, color = State, shape = Extraction)) +
  geom_point(size = 3, alpha = 0.75) +
  xlab("DP (mg/kg)") +
  ylab("d18O") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  #scale_shape_manual(values = c(18)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

# tp data
my_iso_vs_p_plots[[2]] = ggplot(p_pool_isotope_tp_data_summary, aes(x = avg_P_mgperkg, y = avg_d18O, color = State, shape = Extraction)) +
  geom_point(size = 3, alpha = 0.75) +
  xlab("TP (mg/kg)") +
  ylab("d18O") + 
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  #scale_shape_manual(values = c(18)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

setwd("/Users/ssaia/Desktop")
cairo_pdf("iso_vs_hedley_p.pdf", width = 15, height = 8.5, pointsize = 18)
multiplot(plotlist = my_iso_vs_p_plots, cols = 2)
dev.off()

