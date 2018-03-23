# icp calibration script
# last updated 20180229
# contact: sheila saia (sms403@cornell.edu)

# ---- 1 load libraries ----

# clear ws
rm(list = ls())

# load packages
library(tidyverse)

# ---- 2 load data ----

# load soil extraction data
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/icp_data/raw_icp_data/")
nahco3_data_raw = read_csv("raw_reformatted_icp_nahco3.csv")
naoh_data_raw = read_csv("raw_reformatted_icp_naoh.csv")
hcl_data_raw = read_csv("raw_reformatted_icp_hcl.csv")

# load ox and total p data (did analysis at cornell)
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/site_and_exp_design_data")
cornell_data = read_csv("all_joined_data.csv")

# load in soil weight data
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/site_and_exp_design_data/")
soil_wts = read_table2("bonn_hedley1_soil_weights.txt")

# select only what we need
soil_wts_sel = soil_wts %>%
  select(Bonn_SampleID,SampleID:Month,dry_soil_added_g)


# ---- 3 nahco3 data reformatting and processing ----

# select only what we need
nahco3_data_raw_blk = nahco3_data_raw %>% 
  filter(Bonn_SampleID == "blank") %>% 
  select(Bonn_SampleID:P_RSDperc)
nahco3_data_raw_sel = nahco3_data_raw %>%
  filter(Extraction == "nahco3" & Bonn_SampleID != "blank") %>%
  select(Bonn_SampleID:P_RSDperc)

# join soil_wts_sel and nahco3_data_sel
nahco3_raw_data_join = left_join(nahco3_data_raw_sel, soil_wts_sel, by = "Bonn_SampleID")

# calculate the average amount of P in the method blanks
avg_blk_nahco3 = mean(nahco3_data_raw_blk$P_mgperL)

# subtract average amount of P in method blank from samples and calculate TP in mg/kg and umol
nahco3_vol_soln_added = 30 # volume of solution added for extraction (ml)
phosphate_molec_wt = 95 # 95 g of phosphate = 1 mol
nahco3_data = nahco3_raw_data_join %>%
  mutate(Pfix_mgperL = P_mgperL - avg_blk_nahco3) %>%
  mutate(P_mgperkg = Pfix_mgperL * (nahco3_vol_soln_added / dry_soil_added_g) * Dilution) %>%
  mutate(P_umol = (P_mgperkg * dry_soil_added_g) / phosphate_molec_wt) %>%
  mutate(row_num = seq(1:dim(nahco3_raw_data_join)[1])) %>%
  mutate(P_pool = "TP") %>%
  select(row_num, Bonn_SampleID, SampleID, Rep, Month, Extraction, Method, P_pool, P_mgperkg, P_umol)

# for 18O isolation, need at least 30 umol of TP

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/icp_data/")

# export
write_csv(nahco3_data, "icp_nahco3_calibrated.csv")


# ---- 4 naoh data reformatting and processing ----

# select only what we need
naoh_data_raw_blk = naoh_data_raw %>% 
  filter(Bonn_SampleID == "blank") %>% 
  select(Bonn_SampleID:P_RSDperc)
naoh_data_raw_sel = naoh_data_raw %>%
  filter(Extraction == "naoh" & Bonn_SampleID != "blank") %>%
  select(Bonn_SampleID:P_RSDperc)

# join soil_wts_sel and naoh_data_sel
naoh_raw_data_join = left_join(naoh_data_raw_sel, soil_wts_sel, by = "Bonn_SampleID")

# calculate the average amount of P in the method blanks
avg_blk_naoh = mean(naoh_data_raw_blk$P_mgperL)

# subtract average amount of P in method blank from samples and calculate TP in mg/kg and umol
naoh_vol_soln_added = 30 # volume of solution added for extraction (ml)
phosphate_molec_wt=95 # 95 g of phosphate = 1 mol
naoh_data = naoh_raw_data_join %>%
  mutate(Pfix_mgperL = P_mgperL - avg_blk_naoh) %>%
  mutate(P_mgperkg = Pfix_mgperL * (naoh_vol_soln_added / dry_soil_added_g) * Dilution) %>%
  mutate(P_umol = (P_mgperkg * dry_soil_added_g) / phosphate_molec_wt) %>%
  mutate(row_num = seq(1:dim(naoh_raw_data_join)[1])) %>%
  mutate(P_pool = "TP") %>%
  select(row_num, Bonn_SampleID, SampleID, Rep, Month, Extraction, Method, P_pool, P_mgperkg, P_umol)


# for 18O isolation, need at least 30 umol of TP

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/icp_data/")

# export
write_csv(naoh_data, "icp_naoh_calibrated.csv")


# ---- 5 hcl data reformatting and processing ----

# select only what we need
hcl_data_raw_blk = hcl_data_raw %>% 
  filter(Bonn_SampleID == "blank") %>% 
  select(Bonn_SampleID:P_RSDperc)
hcl_data_raw_sel = hcl_data_raw %>%
  filter(Extraction == "hcl" & Bonn_SampleID != "blank") %>%
  select(Bonn_SampleID:P_RSDperc)

# join soil_wts_sel and hcl_data_sel
hcl_raw_data_join = left_join(hcl_data_raw_sel, soil_wts_sel, by = "Bonn_SampleID")

# calculate the average amount of P in the method blanks
avg_blk_hcl = mean(hcl_data_raw_blk$P_mgperL)

# subtract average amount of P in method blank from samples and calculate TP in mg/kg and umol
hcl_vol_soln_added = 30 # volume of solution added for extraction (ml)
phosphate_molec_wt=95 # 95 g of phosphate = 1 mol
hcl_data = hcl_raw_data_join %>%
  mutate(Pfix_mgperL = P_mgperL - avg_blk_hcl) %>%
  mutate(P_mgperkg = Pfix_mgperL * (hcl_vol_soln_added / dry_soil_added_g) * Dilution) %>%
  mutate(P_umol = (P_mgperkg * dry_soil_added_g) / phosphate_molec_wt) %>%
  mutate(row_num = seq(1:dim(hcl_raw_data_join)[1])) %>%
  mutate(P_pool = "TP") %>%
  select(row_num, Bonn_SampleID, SampleID, Rep, Month, Extraction, Method, P_pool, P_mgperkg, P_umol)

# for 18O isolation, need at least 30 umol of TP

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/icp_data/")

# export
write_csv(hcl_data, "icp_hcl_calibrated.csv")


# ---- 6 ox data reformatting ----

# select only what you need from cornell_data
cornell_ox_data = cornell_data %>%
  mutate(Extraction = "ox", Method = "icp", P_pool = "TP") %>%
  select(SampleID, Rep = RepOx, Month, Extraction:P_pool, P_mgperkg = OxPmgkg) #, P_umol = OxPumol)

# define bonn sample ids
bonn_sample_ids_list = soil_wts_sel %>%
  select(Bonn_SampleID, SampleID) %>%
  na.omit() # drops blanks

# join bonn sample ids with cornell_ox_data
ox_data = left_join(bonn_sample_ids_list, cornell_ox_data, by = "SampleID") %>%
  mutate(row_num = seq(1:60)) %>%
  select(row_num, Bonn_SampleID:P_mgperkg) %>%
  mutate(P_umol = NA) # placeholder for now

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/icp_data/")

# export
write_csv(ox_data, "icp_ox_calibrated.csv")


# ---- 7 totalp data reformatting ----

# select only what you need from cornell_data
cornell_totalp_data = cornell_data %>%
  mutate(Extraction = "totalp", Method = "icp", P_pool = "TP") %>%
  select(SampleID, Rep = RepTot, Month, Extraction:P_pool, P_mgperkg = TotPmgkg) #, P_umol = TotPumol)

# define bonn sample ids
bonn_sample_ids_list = soil_wts_sel %>%
  select(Bonn_SampleID, SampleID) %>%
  na.omit() # drops blanks

# join bonn sample ids with cornell_cacl2_data
totalp_data = left_join(bonn_sample_ids_list, cornell_totalp_data, by = "SampleID") %>%
  mutate(row_num = seq(1:60)) %>%
  select(row_num, Bonn_SampleID:P_mgperkg) %>%
  mutate(P_umol = NA) # placeholder for now

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/icp_data/")

# export
write_csv(totalp_data, "icp_totalp_calibrated.csv")
