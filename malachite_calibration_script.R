# malachite calibration script
# last updated 20180229
# contact: sheila saia (sms403@cornell.edu)

# ---- 1 load libraries ----

# clear ws
rm(list = ls())

# load packages
library(tidyverse)

# ---- 2 load data ----

# load soil extraction data
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/malachite_data/raw_malachite_data/")
nahco3_data_raw = read_csv("raw_malachite_nahco3.csv", col_names = FALSE)
naoh_data_raw = read_csv("raw_malachite_naoh.csv", col_names = FALSE)
hcl_data_raw = read_csv("raw_malachite_hcl.csv", col_names = FALSE)
ox_data_raw = read_csv("raw_malachite_ox.csv", col_names = FALSE)

# load in soil weight and cacl2 p data
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/site_and_exp_design_data/")
soil_wts = read_table2("bonn_hedley1_soil_weights.txt")
cornell_data = read_csv("all_joined_data.csv")

# select only what we need
soil_wts_sel = soil_wts %>%
  select(-Notes)


# ---- 3 nahco3 data reformatting and processing ----

# standards
nahco3_standards = data.frame(std_abs = t(nahco3_data_raw[1,1:8]), # standard absorbances from plate reader data
                       vol_added_ul = c(0,20,40,60,80,100,120,200),
                       P_ppm=c(0,0.1,0.2,0.3,0.4,0.5,0.6,1))

# fit standards
nahco3_lm = lm(P_ppm ~ std_abs, data = nahco3_standards)

# print summary to check
summary(nahco3_lm)

# save variables for plotting
nahco3_intercept = summary(nahco3_lm)$coefficients[1,1]
nahco3_slope = summary(nahco3_lm)$coefficients[2,1]
nahco3_rsquared = summary(nahco3_lm)$r.squared

# plot line and observations for standard curve
plot(P_ppm ~ std_abs, data = nahco3_standards, pch = 16, xlab = "Absorbance", ylab = "Known PO4 as P Concentration (ppm)",
     main = paste("y =", round(nahco3_slope, 4), "x + ", round(nahco3_intercept, 4), ", r-squared =", round(nahco3_rsquared, 4), sep=" "))
abline(a = nahco3_intercept, b = nahco3_slope, col = "red", lwd = 3)

# make new dataframe to add sample abs readings
samp_slots = 96-12 #96-12 because first row are standards
nahco3_samples = data.frame(AssayID = rep(0, samp_slots), abs = rep(0, samp_slots))

# for loop to save abs data to sample table
nahco3_data_transposed = as.data.frame(t(nahco3_data_raw)) #transpose dataframe if samples are across rows
strt = 1 #start iterator outside loop to keep track of samples
for (i in 2:8) {
  end = strt+11
  nahco3_samples$abs[strt:end] = nahco3_data_transposed[,i]
  strt = strt + 12 #add 12 to iterator
}

# define AssayID for sample data
num_samples = 30
assay_id_list = seq(1:num_samples)
num_reps = 3
num_reps_per_sample = c(rep(c(3, 3, 2), num_samples/num_reps)) #number of reDPtitions you analyze each sample (A and B have 3 reps, C has 2 reps)
strt = 1
for (i in 1:num_samples) {
  end = strt + num_reps_per_sample[i] - 1
  nahco3_samples$AssayID[strt:end] = rep(assay_id_list[i], num_reps_per_sample[i])
  strt = strt + num_reps_per_sample[i]
}

# calculate raw and calibrated concentration
vol_soln_added = 30 #volume of solution added for extraction (ml)
nahco3_dilution = 200/20 #dilution multiple i.e. 10x dilution is represented as 10
phosphate_molec_wt=95 #95 g of phosphate = 1 mol
nahco3_data = nahco3_samples %>%
  mutate(rawP_ppm = abs * nahco3_slope + nahco3_intercept) %>%
  left_join(soil_wts_sel, by = "AssayID") %>% # join soil_wts_data
  mutate(row_num = seq(1:dim(nahco3_samples)[1])) %>%
  mutate(Extraction = "nahco3", Method = "malachite", P_pool = "DP") %>%
  mutate(P_mgperkg = rawP_ppm * (vol_soln_added / dry_soil_added_g) * nahco3_dilution) %>% # inorganic P pool = DP for dissolved P
  mutate(P_umol = (P_mgperkg * dry_soil_added_g) / phosphate_molec_wt) %>% # DP in micromoles
  select(row_num, Bonn_SampleID, SampleID, Month, Rep, Extraction:P_umol) %>%
  na.omit()

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/malachite_data/")

# export
write_csv(nahco3_data, "malachite_nahco3_calibrated.csv")


# ---- 4 naoh data reformatting and processing ----

# standards
naoh_standards = data.frame(std_abs = t(naoh_data_raw[1,1:8]), # standard absorbances from plate reader data
                              vol_added_ul = c(0,20,40,60,80,100,120,200),
                              P_ppm=c(0,0.1,0.2,0.3,0.4,0.5,0.6,1))

# fit standards
naoh_lm = lm(P_ppm ~ std_abs, data = naoh_standards)

# print summary to check
summary(naoh_lm)

# save variables for plotting
naoh_intercept = summary(naoh_lm)$coefficients[1,1]
naoh_slope = summary(naoh_lm)$coefficients[2,1]
naoh_rsquared = summary(naoh_lm)$r.squared

# plot line and observations for standard curve
plot(P_ppm ~ std_abs, data = naoh_standards, pch = 16, xlab = "Absorbance", ylab = "Known PO4 as P Concentration (ppm)",
     main = paste("y =", round(naoh_slope, 4), "x + ", round(naoh_intercept, 4), ", r-squared =", round(naoh_rsquared, 4), sep=" "))
abline(a = naoh_intercept, b = naoh_slope, col = "red", lwd = 3)

# make new dataframe to add sample abs readings
samp_slots = 96-12 #96-12 because first row are standards
naoh_samples = data.frame(AssayID = rep(0, samp_slots), abs = rep(0, samp_slots))

# for loop to save abs data to sample table
naoh_data_transposed = as.data.frame(t(naoh_data_raw)) #transpose dataframe if samples are across rows
strt = 1 #start iterator outside loop to keep track of samples
for (i in 2:8) {
  end = strt+11
  naoh_samples$abs[strt:end] = naoh_data_transposed[,i]
  strt = strt + 12 #add 12 to iterator
}

# define AssayID for sample data
num_samples = 30
assay_id_list = seq(1:num_samples)
num_reps = 3
num_reps_per_sample = c(rep(c(3, 3, 2), num_samples/num_reps)) #number of reDPtitions you analyze each sample (A and B have 3 reps, C has 2 reps)
strt = 1
for (i in 1:num_samples) {
  end = strt + num_reps_per_sample[i] - 1
  naoh_samples$AssayID[strt:end] = rep(assay_id_list[i], num_reps_per_sample[i])
  strt = strt + num_reps_per_sample[i]
}

# calculate raw and calibrated concentration
vol_soln_added = 30 #volume of solution added for extraction (ml)
naoh_dilution = 200/5 #dilution multiple i.e. 40x dilution is represented as 40
phosphate_molec_wt=95 #95 g of phosphate = 1 mol
naoh_data = naoh_samples %>%
  mutate(rawP_ppm = abs * naoh_slope + naoh_intercept) %>%
  left_join(soil_wts_sel, by = "AssayID") %>% # join soil_wts_data
  mutate(row_num = seq(1:dim(naoh_samples)[1])) %>%
  mutate(Extraction = "naoh", Method = "malachite", P_pool = "DP") %>%
  mutate(P_mgperkg = rawP_ppm * (vol_soln_added / dry_soil_added_g) * naoh_dilution) %>% # inorganic P pool = DP for dissolved P
  mutate(P_umol = (P_mgperkg * dry_soil_added_g) / phosphate_molec_wt) %>% # DP in micromoles
  select(row_num, Bonn_SampleID, SampleID, Month, Rep, Extraction:P_umol) %>%
  na.omit()

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/malachite_data/")

# export
write_csv(naoh_data, "malachite_naoh_calibrated.csv")


# ---- 5 hcl data reformatting and processing ----

# standards
hcl_standards = data.frame(std_abs = t(hcl_data_raw[1,1:8]), # standard absorbances from plate reader data
                            vol_added_ul = c(0,20,40,60,80,100,120,200),
                            P_ppm=c(0,0.1,0.2,0.3,0.4,0.5,0.6,1))

# fit standards
hcl_lm = lm(P_ppm ~ std_abs, data = hcl_standards)

# print summary to check
summary(hcl_lm)

# save variables for plotting
hcl_intercept = summary(hcl_lm)$coefficients[1,1]
hcl_slope = summary(hcl_lm)$coefficients[2,1]
hcl_rsquared = summary(hcl_lm)$r.squared

# plot line and observations for standard curve
plot(P_ppm ~ std_abs, data = hcl_standards, pch = 16, xlab = "Absorbance", ylab = "Known PO4 as P Concentration (ppm)",
     main = paste("y =", round(hcl_slope, 4), "x + ", round(hcl_intercept, 4), ", r-squared =", round(hcl_rsquared, 4), sep=" "))
abline(a = hcl_intercept, b = hcl_slope, col = "red", lwd = 3)

# make new dataframe to add sample abs readings
samp_slots = 96-12 #96-12 because first row are standards
hcl_samples = data.frame(AssayID = rep(0, samp_slots), abs = rep(0, samp_slots))

# for loop to save abs data to sample table
hcl_data_transposed = as.data.frame(t(hcl_data_raw)) #transpose dataframe if samples are across rows
strt = 1 #start iterator outside loop to keep track of samples
for (i in 2:8) {
  end = strt+11
  hcl_samples$abs[strt:end] = hcl_data_transposed[,i]
  strt = strt + 12 #add 12 to iterator
}

# define AssayID for sample data
num_samples = 30
assay_id_list = seq(1:num_samples)
num_reps = 3
num_reps_per_sample = c(rep(c(3, 3, 2), num_samples/num_reps)) #number of reDPtitions you analyze each sample (A and B have 3 reps, C has 2 reps)
strt = 1
for (i in 1:num_samples) {
  end = strt + num_reps_per_sample[i] - 1
  hcl_samples$AssayID[strt:end] = rep(assay_id_list[i], num_reps_per_sample[i])
  strt = strt + num_reps_per_sample[i]
}

# calculate raw and calibrated concentration
vol_soln_added = 30 #volume of solution added for extraction (ml)
hcl_dilution = 200/5 #dilution multiple i.e. 40x dilution is represented as 40
phosphate_molec_wt=95 #95 g of phosphate = 1 mol
hcl_data = hcl_samples %>%
  mutate(rawP_ppm = abs * hcl_slope + hcl_intercept) %>%
  left_join(soil_wts_sel, by = "AssayID") %>% # join soil_wts_data
  mutate(row_num = seq(1:dim(hcl_samples)[1])) %>%
  mutate(Extraction = "hcl", Method = "malachite", P_pool = "DP") %>%
  mutate(P_mgperkg = rawP_ppm * (vol_soln_added / dry_soil_added_g) * hcl_dilution) %>% # inorganic P pool = DP for dissolved P
  mutate(P_umol = (P_mgperkg * dry_soil_added_g) / phosphate_molec_wt) %>% # DP in micromoles
  select(row_num, Bonn_SampleID, SampleID, Month, Rep, Extraction:P_umol) %>%
  na.omit()

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/malachite_data/")

# export
write_csv(hcl_data, "malachite_hcl_calibrated.csv")


# ---- 6 ox data reformatting and processing ----

# standards
ox_standards = data.frame(std_abs = t(ox_data_raw[1,1:8]), # standard absorbances from plate reader data
                           vol_added_ul = c(0,20,40,60,80,100,120,200),
                           P_ppm=c(0,0.1,0.2,0.3,0.4,0.5,0.6,1))

# fit standards
ox_lm = lm(P_ppm ~ std_abs, data = ox_standards)

# print summary to check
summary(ox_lm)

# save variables for plotting
ox_intercept = summary(ox_lm)$coefficients[1,1]
ox_slope = summary(ox_lm)$coefficients[2,1]
ox_rsquared = summary(ox_lm)$r.squared

# plot line and observations for standard curve
plot(P_ppm ~ std_abs, data = ox_standards, pch = 16, xlab = "Absorbance", ylab = "Known PO4 as P Concentration (ppm)",
     main = paste("y =", round(ox_slope, 4), "x + ", round(ox_intercept, 4), ", r-squared =", round(ox_rsquared, 4), sep=" "))
abline(a = ox_intercept, b = ox_slope, col = "red", lwd = 3)

# make new dataframe to add sample abs readings
samp_slots = 96-12 #96-12 because first row are standards
ox_samples = data.frame(AssayID = rep(0, samp_slots), abs = rep(0, samp_slots))

# for loop to save abs data to sample table
ox_data_transposed = as.data.frame(t(ox_data_raw)) #transpose dataframe if samples are across rows
strt = 1 #start iterator outside loop to keep track of samples
for (i in 2:8) {
  end = strt+11
  ox_samples$abs[strt:end] = ox_data_transposed[,i]
  strt = strt + 12 #add 12 to iterator
}

# define AssayID for sample data
num_samples = 30
assay_id_list = seq(1:num_samples)
num_reps = 3
num_reps_per_sample = c(rep(c(3, 3, 2), num_samples/num_reps)) #number of reDPtitions you analyze each sample (A and B have 3 reps, C has 2 reps)
strt = 1
for (i in 1:num_samples) {
  end = strt + num_reps_per_sample[i] - 1
  ox_samples$AssayID[strt:end] = rep(assay_id_list[i], num_reps_per_sample[i])
  strt = strt + num_reps_per_sample[i]
}

# calculate raw and calibrated concentration
vol_soln_added = 30 #volume of solution added for extraction (ml)
ox_dilution = 200/1 #dilution multiple i.e. 200x dilution is represented as 200
phosphate_molec_wt=95 #95 g of phosphate = 1 mol
ox_data = ox_samples %>%
  mutate(rawP_ppm = abs * ox_slope + ox_intercept) %>%
  left_join(soil_wts_sel, by = "AssayID") %>% # join soil_wts_data
  mutate(row_num = seq(1:dim(ox_samples)[1])) %>%
  mutate(Extraction = "ox", Method = "malachite", P_pool = "DP") %>%
  mutate(P_mgperkg = rawP_ppm * (vol_soln_added / dry_soil_added_g) * ox_dilution) %>% # inorganic P pool = DP for dissolved P
  mutate(P_umol = (P_mgperkg * dry_soil_added_g) / phosphate_molec_wt) %>% # DP in micromoles
  select(row_num, Bonn_SampleID, SampleID, Month, Rep, Extraction:P_umol) %>%
  na.omit()

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/malachite_data/")

# export
write_csv(ox_data, "malachite_ox_calibrated.csv")

# ---- 7 cacl2 data reformatting ----

# select only what you need from cornell_data
cornell_cacl2_data = cornell_data %>%
  mutate(Extraction = "cacl2", Method = "molybdenum", P_pool = "DP") %>%
  select(SampleID, Rep = RepCacl2, Month, Extraction:P_pool, P_mgperkg = SRPmgkg) #, P_umol = SRPumol)

# define bonn sample ids
bonn_sample_ids_list = soil_wts_sel %>%
  select(Bonn_SampleID, SampleID) %>%
  na.omit() # drops blanks

# join bonn sample ids with cornell_cacl2_data
cacl2_data = left_join(bonn_sample_ids_list, cornell_cacl2_data, by = "SampleID") %>%
  mutate(row_num = seq(1:60)) %>%
  select(row_num, Bonn_SampleID:P_mgperkg) %>%
  mutate(P_umol = NA) # placeholder for now

# set working directory
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/malachite_data/")

# export
write_csv(cacl2_data, "molybdenum_cacl2_calibrated.csv")

