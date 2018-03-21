# malachite calibration script
# last updated 20180229
# contact: sheila saia (sms403@cornell.edu)

# ---- 1 load libraries ----

library(tidyverse)

# ---- 2 load data ----

# set working directory for soil extraction data
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/malachite_data/raw_malachite_data/")

# load soil extraction data
nahco3_data_raw = read_csv("raw_malachite_nahco3.csv", col_names = FALSE)
naoh_data_raw = read_csv("raw_malachite_naoh.csv", col_names = FALSE)
hcl_data_raw = read_csv("raw_malachite_hcl.csv", col_names = FALSE)

# set working directory for soil extraction weights
setwd("/Users/ssaia/Documents/phd/oxygen_isotopes/bonn_soil_analysis/data/site_and_exp_design_data/")

# load in soil weight data
soil_wts = read_table2("bonn_hedley1_soil_weights.txt")

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
#nahco3 dulted 200/20 = 10x, naoh diluted 200/5 = 40x, hcl diluted 200/5 = 40x, ox diluted 200/1 = 200x
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



