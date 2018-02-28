# icp reformatting and processing script
# last updated 20180229
# contact: sheila saia (sms403@cornell.edu)

# ---- 1 load libraries ----

library(tidyverse)

# ---- 2 load data ----

# set working directory for soil extraction data
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/icp_data/raw_icp_data/")

# load soil extraction data
nahco3_data_raw = read_csv("raw_reformatted_icp_nahco3.csv")
naoh_data_raw = read_csv("raw_reformatted_icp_naoh.csv")
hcl_data_raw = read_csv("raw_reformatted_icp_hcl.csv")

# set working directory for soil extraction weights
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/site_and_exp_design_data/")

# load in soil weight data
soil_wts = read_table2("bonn_hedley1_soil_weights.txt")

# select only what we need
soil_wts_sel = soil_wts %>%
  select(Bonn_SampleID,SampleID:Month,dry_soil_added_g)


# ---- 3 nahco3 data reformatting and processing ----

# select only what we need
nahco3_data_raw_blk = nahco3_data_raw %>% 
  filter(Bonn_SampleID == "blank") %>% 
  select(RowNum, Bonn_SampleID:P_RSDperc)
nahco3_data_raw_sel = nahco3_data_raw %>%
  filter(Extraction == "nahco3" & Bonn_SampleID != "blank") %>%
  select(RowNum, Bonn_SampleID:P_RSDperc)

# join soil_wts_sel and nahco3_data_sel
nahco3_raw_data_join = left_join(nahco3_data_raw_sel, soil_wts_sel, by = "Bonn_SampleID")

# calculate the average amount of P in the method blanks
avg_blk_nahco3 = mean(nahco3_data_raw_blk$P_mgperL)

# subtract average amount of P in method blank from samples and calculate TP in mg/kg and umol
vol_soln_added = 30 # volume of solution added for extraction (ml)
phosphate_molec_wt=95 # 95 g of phosphate = 1 mol
nahco3_data = nahco3_raw_data_join %>%
  mutate(Pfix_mgperL = P_mgperL - avg_blk_nahco3) %>%
  mutate(TP_mgperkg = Pfix_mgperL * (vol_soln_added / dry_soil_added_g) * Dilution) %>%
  mutate(TP_umol = (TP_mgperkg * dry_soil_added_g) / phosphate_molec_wt) %>%
  select(-P_mgperL, -P_RSDperc, -Pfix_mgperL)

# for 18O isolation, need at least 30 umol of TP

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/icp_data/")

# export
write_csv(nahco3_data, "icp_nahco3_calibrated.csv")


# ---- 4 naoh data reformatting and processing ----

# select only what we need
naoh_data_raw_blk = naoh_data_raw %>% 
  filter(Bonn_SampleID == "blank") %>% 
  select(RowNum, Bonn_SampleID:P_RSDperc)
naoh_data_raw_sel = naoh_data_raw %>%
  filter(Extraction == "naoh" & Bonn_SampleID != "blank") %>%
  select(RowNum, Bonn_SampleID:P_RSDperc)

# join soil_wts_sel and naoh_data_sel
naoh_raw_data_join = left_join(naoh_data_raw_sel, soil_wts_sel, by = "Bonn_SampleID")

# calculate the average amount of P in the method blanks
avg_blk_naoh = mean(naoh_data_raw_blk$P_mgperL)

# subtract average amount of P in method blank from samples and calculate TP in mg/kg and umol
vol_soln_added = 30 # volume of solution added for extraction (ml)
phosphate_molec_wt=95 # 95 g of phosphate = 1 mol
naoh_data = naoh_raw_data_join %>%
  mutate(Pfix_mgperL = P_mgperL - avg_blk_naoh) %>%
  mutate(TP_mgperkg = Pfix_mgperL * (vol_soln_added / dry_soil_added_g) * Dilution) %>%
  mutate(TP_umol = (TP_mgperkg * dry_soil_added_g) / phosphate_molec_wt) %>%
  select(-P_mgperL, -P_RSDperc, -Pfix_mgperL)

# for 18O isolation, need at least 30 umol of TP

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/icp_data/")

# export
write_csv(naoh_data, "icp_naoh_calibrated.csv")


# ---- 5 hcl data reformatting and processing ----

# select only what we need
hcl_data_raw_blk = hcl_data_raw %>% 
  filter(Bonn_SampleID == "blank") %>% 
  select(RowNum, Bonn_SampleID:P_RSDperc)
hcl_data_raw_sel = hcl_data_raw %>%
  filter(Extraction == "hcl" & Bonn_SampleID != "blank") %>%
  select(RowNum, Bonn_SampleID:P_RSDperc)

# join soil_wts_sel and hcl_data_sel
hcl_raw_data_join = left_join(hcl_data_raw_sel, soil_wts_sel, by = "Bonn_SampleID")

# calculate the average amount of P in the method blanks
avg_blk_hcl = mean(hcl_data_raw_blk$P_mgperL)

# subtract average amount of P in method blank from samples and calculate TP in mg/kg and umol
vol_soln_added = 30 # volume of solution added for extraction (ml)
phosphate_molec_wt=95 # 95 g of phosphate = 1 mol
hcl_data = hcl_raw_data_join %>%
  mutate(Pfix_mgperL = P_mgperL - avg_blk_hcl) %>%
  mutate(TP_mgperkg = Pfix_mgperL * (vol_soln_added / dry_soil_added_g) * Dilution) %>%
  mutate(TP_umol = (TP_mgperkg * dry_soil_added_g) / phosphate_molec_wt) %>%
  select(-P_mgperL, -P_RSDperc, -Pfix_mgperL)

# for 18O isolation, need at least 30 umol of TP

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/icp_data/")

# export
write_csv(hcl_data, "icp_hcl_calibrated.csv")

