#============================================================================.
# Project: RBODI for assessing organ damage in SLE based on national register data
# Author: Álvaro Gómez and Annica Dominicus 
# Date: 2025-02-24, last update 2025-06-10
# Purpose: Derive Register-based Organ Damage Index (RBODI) for SLE
# Note: Code based on original script from 2023-12-17
#============================================================================.

# Program settings ---
# =============================================================================.

# Load packages
library(haven)
library(lubridate)
library(tidyverse)

# Definitions used
def_dpyr <- 365.25 # number of days per year
def_lookback <- 150*365 # time period before indexdate to identify prevalent cases for each item included in RBODI - here including diagnoses over the life span

#============================================================================.
# Function filter_by_code ----
#============================================================================.
# Select codes for medications (ATC) or diagnoses (ICD-10) including all lower level codes
#============================================================================.
# Arguments:
# df = dataset
# code_vector = vector of codes to extract
# code_column_name = name of variable containing the codes to be extracted
#============================================================================.
filter_by_code <- function(data, code_vector, code_column_name = "icd") {
  # Create pattern for string matching
  code_patterns <- paste0("^(", paste(code_vector, collapse = "|"), ")(\\d*|\\w*)$")
  
  # Filter the dataset where ICD codes match the pattern
  filtered_data <- data[grepl(code_patterns, data[[code_column_name]]), ]
  
  return(filtered_data)
}

#============================================================================.
# Function select_1st_diag_date ----
#============================================================================.
# Select first diagnosis date allowing for empty (null) dataset
#============================================================================.
# Arguments:
# df = dataset
#============================================================================.
select_1st_diag_date <- function(df){
  if(!is.null(df)){
    df %>%
      arrange(lopnr, diag_date) %>%
      group_by(lopnr) %>%
      filter(row_number() == 1)
  }
}

#============================================================================.
# Function add_score_post_pre ----
#============================================================================.
# Derive score post and pre indexdate
# Arguments:
# df = dataset
#============================================================================.
add_score_post_pre <- function(df){
  if(!is.null(df)){
    df %>%
      mutate(tmp_score=ifelse(incident==1 & bl_prevalent!=1, 1,0),
             tmp_score_bl=ifelse(bl_prevalent==1, 1,0))
  }
}

#============================================================================.
# Function post_pre ----
#============================================================================.
# Calculate first diagnosis date after indexdate (diag_date)
# Calculate last diagnosis date before indexdate (diag_date_pre)
# Indicators for whether incident (variable incident) or prevalent at indexdate (variable bl_prevalent)
#============================================================================.
# Arguments:
# df = dataset
# lookback = number of days to use as lookback period prior to indexdate to 
#            look for diagnoses to identify prevalent cases
#============================================================================.
post_pre <- function(df, lookback = NULL){
  
  df_all <- NULL
  
  if(!is.null(df) & "diag_date" %in% names(df)){
    df_post <- df %>% 
      filter(diag_date >= indexdate, diag_date <= exitdate) %>%
      arrange(lopnr, diag_date) %>%
      group_by(lopnr) %>%
      filter(row_number() == 1) 
    
    if(is.null(lookback)){
      df_pre <- df %>%
        filter(diag_date < indexdate) 
    } else {
      df_pre <- df %>%
        filter(diag_date >= indexdate - lookback, diag_date < indexdate)
    }
    
    df_pre <- df_pre %>% 
      arrange(lopnr, desc(diag_date)) %>%
      group_by(lopnr) %>%
      filter(row_number() == 1) %>%
      rename(diag_date_pre = diag_date)
    
    df_all <- safe_join(df_post, df_pre, by=c("lopnr", "indexdate", "exitdate"), join_type = "full") %>%
      mutate(incident = ifelse(!is.na(diag_date) & is.na(diag_date_pre),1,0),
             bl_prevalent = ifelse(!is.na(diag_date_pre),1,0))
    
  }
           
  return(df_all) 
}

#============================================================================.
# Function select_post_pre ----
#============================================================================.
# Calculate first diagnosis date after indexdate (diag_date)
# Calculate last diagnosis date before indexdate (diag_date_pre)
# Indicators for whether incident (incident) or prevalent at indexdate (bl_prevalent)
#============================================================================.
# Arguments:
# df = dataset
# lookback = number of days to use as lookback period prior to indexdate to 
#            look for diagnoses to identify prevalent cases
# set_score = score to be given for diagnosis to be extracted for RBODI calculation
# diagnosis_code = ICD-10 code
# diagnosis = name of diagnosis
#============================================================================.
select_post_pre <- function(df, lookback = NULL, set_score = NULL, diagnosis_code, diagnosis){
  
  df_all <- NULL
  
  if(!is.null(df) & "diag_date" %in% names(df)){
    df_post <- df %>% 
      filter(diag_date >= indexdate, diag_date <= exitdate) %>%
      arrange(lopnr, diag_date) %>%
      group_by(lopnr) %>%
      filter(row_number() == 1) 
    
    if(is.null(lookback)){
      df_pre <- df %>%
        filter(diag_date < indexdate) 
    } else {
      df_pre <- df %>%
        filter(diag_date >= indexdate - lookback, diag_date < indexdate)
    }
    
    df_pre <- df_pre %>% 
      arrange(lopnr, desc(diag_date)) %>%
      group_by(lopnr) %>%
      filter(row_number() == 1) %>%
      rename(diag_date_pre = diag_date)
    
    df_all <- safe_join(df_post, df_pre, by=c("lopnr", "indexdate", "exitdate"), join_type = "full") %>%
      mutate(incident = ifelse(!is.na(diag_date) & is.na(diag_date_pre),1,0),
             bl_prevalent = ifelse(!is.na(diag_date_pre),1,0),
             diagnosis_code = diagnosis_code,
             diagnosis = diagnosis)
  }
  
  if(!is.null(set_score)){
    df_all <- df_all %>%
      mutate(score=ifelse(incident==1 & bl_prevalent!=1, set_score,0),
             score_bl=ifelse(bl_prevalent==1, set_score,0))
  }
           
  return(df_all) 
}

#============================================================================.
# Function diag_date ----
#============================================================================.
# Extract all selected diagnosis dates from outpatient and inpatient care register
#============================================================================.
# Arguments:
# cohort = cohort dataset
# oppen = dataset with outpatient care visits
# sluten = dataset with inpatient care visits
# codes = ICD codes to be extracted
#============================================================================.
diag_date <- function(cohort, oppen=NULL, sluten=NULL, codes = NULL){
  
  # Outpatient
  if(!is.null(oppen)){
    names(oppen) <- tolower(names(oppen))
    oppen <- rename(oppen, diag_date = indatum) 
    
    if(!is.null(codes)){
      # Selecting from all diagnoses
      all_vars <- names(oppen)
      vars <- grep("^dia", all_vars, value = TRUE)
      vars <- vars[vars != "diag_date"]
      
      df1 <- select(cohort, lopnr) %>%
        inner_join(oppen) %>%
        group_by(lopnr, diag_date) %>%
        pivot_longer(cols = c("hdia", vars),
                     names_to = "diagnosis",
                     values_to = "icd") %>%
        filter_by_code(codes, "icd") %>%
        select(lopnr, diag_date) %>%
        distinct()
    } else{
      df1 <- select(cohort, lopnr) %>%
        inner_join(oppen) %>%
        select(lopnr, diag_date) %>%
        distinct()
    }
  }
  
  # Inpatient
  if(!is.null(sluten)){
    names(sluten) <- tolower(names(sluten))
    sluten <- rename(sluten, diag_date = utdatum)
    
    if(!is.null(codes)){
      # Selecting from all diagnoses
      all_vars <- names(sluten)
      vars <- grep("^dia", all_vars, value = TRUE)
      vars <- vars[vars != "diag_date"]
      
      df2 <- select(cohort, lopnr) %>%
        inner_join(sluten) %>%
        group_by(lopnr, diag_date) %>%
        pivot_longer(cols = c("hdia", vars),
                     names_to = "diagnosis",
                     values_to = "icd") %>%
        filter_by_code(codes, "icd") %>%
        select(lopnr, diag_date) %>%
        distinct()
    } else {
      df2 <- select(cohort, lopnr) %>%
        inner_join(sluten) %>%
        select(lopnr, diag_date) %>%
        distinct()
    }
    
    
  }
  
  # Check if first dataset is NULL, empty or missing
  df1_empty <- is.null(oppen) || is.null(df1) || length(df1) == 0 || (is.data.frame(df1) && nrow(df1) == 0)
  
  # Check if second dataset is NULL, empty or missing
  df2_empty <- is.null(sluten) || is.null(df2) || length(df2) == 0 || (is.data.frame(df2) && nrow(df2) == 0)
  
  # Handle all possible cases
  if (df1_empty && df2_empty) {
    return(NULL)  # Return empty data frame if both are empty
  } else if (df1_empty) {
    return(df2)
  } else if (df2_empty) {
    return(df1)
  } else {
    # If both datasets exist, bind them together
    return(distinct(rbind(df1, df2)))
  }
}

#============================================================================.
# Function proc_date ----
#============================================================================.
# Extract all selected procedure dates from outpatient and inpatient care register
#============================================================================.
# Arguments:
# cohort = cohort dataset
# oppen = dataset with outpatient care visits
# sluten = dataset with inpatient care visits
# codes = KVÅ codes to be extracted
#============================================================================.
proc_date <- function(cohort, oppen=NULL, sluten=NULL, codes = NULL){
  
  # Outpatient
  if(!is.null(oppen)){
    names(oppen) <- tolower(names(oppen))
    oppen <- rename(oppen, diag_date = indatum) 
    
    if(!is.null(codes)){
      df1 <- select(cohort, lopnr) %>%
        inner_join(oppen) %>%
        group_by(lopnr, diag_date) %>%
        pivot_longer(cols = c(starts_with("op")),
                     names_to = "procedure",
                     values_to = "kva") %>%
        filter_by_code(codes, "kva") %>%
        select(lopnr, diag_date) %>%
        distinct()
    } else {
      df1 <- select(cohort, lopnr) %>%
      inner_join(oppen) %>%
      select(lopnr, diag_date) %>%
      distinct()
    }
  }
  
  # Inpatient
  if(!is.null(sluten)){
    names(sluten) <- tolower(names(sluten))
    sluten <- rename(sluten, diag_date = utdatum) 
    
    if(!is.null(codes)){
      df2 <- select(cohort, lopnr) %>%
        inner_join(sluten) %>%
        group_by(lopnr, diag_date) %>%
        pivot_longer(cols = c(starts_with("op")),
                     names_to = "procedure",
                     values_to = "kva") %>%
        filter_by_code(codes, "kva") %>%
        select(lopnr, diag_date) %>%
        distinct()
    } else {
      df2 <- select(cohort, lopnr) %>%
        inner_join(sluten) %>%
        select(lopnr, diag_date) %>%
        distinct()
    }
  }
  
  # Check if first dataset is NULL, empty or missing
  df1_empty <- is.null(oppen) || is.null(df1) || length(df1) == 0 || (is.data.frame(df1) && nrow(df1) == 0)
  
  # Check if second dataset is NULL, empty or missing
  df2_empty <- is.null(sluten) || is.null(df2) || length(df2) == 0 || (is.data.frame(df2) && nrow(df2) == 0)
  
  # Handle all possible cases
  if (df1_empty && df2_empty) {
    return(NULL)  # Return empty data frame if both are empty
  } else if (df1_empty) {
    return(df2)
  } else if (df2_empty) {
    return(df1)
  } else {
    # If both datasets exist, bind them together
    return(distinct(rbind(df1, df2)))
  }
}

#============================================================================.
# Function drg_date ----
#============================================================================.
# Extract all selected procedure dates from outpatient and inpatient care register
#============================================================================.
# Arguments:
# cohort = cohort dataset
# oppen = dataset with outpatient care visits
# sluten = dataset with inpatient care visits
# codes = DRG codes to be extracted
#============================================================================.
drg_date <- function(cohort, oppen=NULL, sluten=NULL, codes = NULL){
  
  # Outpatient
  if(!is.null(oppen)){
    names(oppen) <- tolower(names(oppen))
    oppen <- rename(oppen, diag_date = indatum)
    
    if(!is.null(codes)){
      df1 <- select(cohort, lopnr) %>%
        inner_join(oppen) %>%
        filter_by_code(codes, "drg") %>%
        select(lopnr, diag_date) %>%
        distinct()
    } else {
      df1 <- select(cohort, lopnr) %>%
        inner_join(oppen) %>%
        select(lopnr, diag_date) %>%
        distinct()
    }
  }
  
  # Inpatient
  if(!is.null(sluten)){
    names(sluten) <- tolower(names(sluten))
    sluten <- rename(sluten, diag_date = utdatum) %>%
      select(lopnr, diag_date, drg)
    
    if(!is.null(codes)){
     df2 <- select(cohort, lopnr) %>%
        inner_join(sluten) %>%
        filter_by_code(codes, "drg") %>%
        select(lopnr, diag_date) %>%
        distinct()
    } else {
      df1 <- select(cohort, lopnr) %>%
        inner_join(sluten) %>%
        select(lopnr, diag_date) %>%
        distinct()
    }
  }
  
  # Check if first dataset is NULL, empty or missing
  df1_empty <- is.null(oppen) || is.null(df1) || length(df1) == 0 || (is.data.frame(df1) && nrow(df1) == 0)
  
  # Check if second dataset is NULL, empty or missing
  df2_empty <- is.null(sluten) || is.null(df2) || length(df2) == 0 || (is.data.frame(df2) && nrow(df2) == 0)
  
  # Handle all possible cases
  if (df1_empty && df2_empty) {
    return(NULL)  # Return empty data frame if both are empty
  } else if (df1_empty) {
    return(df2)
  } else if (df2_empty) {
    return(df1)
  } else {
    # If both datasets exist, bind them together
    return(rbind(df1, df2))
  }
}

#============================================================================.
# Function med_date ----
#============================================================================.
# Extract all dates of dispensations for selected medications from the Prescribed Drug Register
#============================================================================.
# Arguments:
# cohort = cohort dataset
# lm = dataset with dispensations
# codes = ATC codes to be extracted
# gap = number of days between two dispensations to be interpreted as a new period of treatment
# span = number of days required to be classified as repeated dispensations (within same period of treatment)
#============================================================================.
med_date <- function(cohort, lm, codes = NULL, gap, span){
  
    names(lm) <- tolower(names(lm))
    
    if(!is.null(codes)){
      lm <- lm %>%
        filter_by_code(codes, "atc")
    }
        
    df <- select(cohort, lopnr) %>%
      left_join(select(lm, lopnr, edatum, atc)) %>%
      rename(med_date = edatum) %>%  
      arrange(lopnr, med_date) %>%
      group_by(lopnr) %>%
      mutate(last2 = difftime(med_date, lag(med_date), units = "days"),
             med_gap = ifelse(last2 > gap, 1, 0),
             newseq = ifelse(med_gap ==0 | is.na(med_gap),0,1),
             seq = cumsum(newseq)) %>%
      group_by(lopnr, seq) %>%
      mutate(med_date_1st = as.Date(ifelse(is.infinite(suppressWarnings(min(med_date, na.rm = TRUE))), 
                                           NA, suppressWarnings(min(med_date, na.rm = TRUE))), origin = "1970-01-01"),
             med_span = difftime(med_date, med_date_1st, units = "days"),
             med_span_seq = ifelse(is.infinite(suppressWarnings(max(med_span, na.rm=TRUE))), 
                                   NA, suppressWarnings(max(med_span, na.rm=TRUE))), 
             med_seq_sel = ifelse(med_span_seq >= span, 1, 0)
      ) %>%
      filter(med_seq_sel == 1) %>%
      ungroup() %>%
      select(lopnr, med_date_1st) %>%
      rename(med_date = med_date_1st) %>%
      distinct()
  
    if(is.null(df) | nrow(df)==0) df <- NULL
    
  return(df)
}

#============================================================================.
# Function safe_join ----
#============================================================================.
# Combine datasets allowing for some datasets may be missing 
#============================================================================.
# Arguments:
# by = patient identifier
# join_type = type of join of datasets
#============================================================================.
safe_join <- function(..., by = "lopnr", join_type = "full") {
  datasets <- list(...)
  non_null_datasets <- Filter(Negate(is.null), datasets)
  
  if (length(non_null_datasets) == 0) {
    return(NULL)
  }
  
  result <- non_null_datasets[[1]]
  
  if (length(non_null_datasets) >1) {
    for (i in 2:length(non_null_datasets)) {
      if (join_type == "left") {
        result <- left_join(result, non_null_datasets[[i]], by = by)
      } else if (join_type == "inner") {
        result <- inner_join(result, non_null_datasets[[i]], by = by)
      } else if (join_type == "full") {
        result <- full_join(result, non_null_datasets[[i]], by = by)
      } else if (join_type == "right") {
        result <- right_join(result, non_null_datasets[[i]], by = by)
      }
    }
  }

  return(result)
}

#============================================================================.
# Function combine_diag_med ----
#============================================================================.
# Combine information on selected diagnoses and medications
#============================================================================.
# Arguments:
# diag = dataset with selected diagnoses
# med = dataset with selected medications
# diffdays = number of days within which diagnoses and dispensations should occur to be classified as an event
# method = "or" if one of diagnosis/medication required, "and" if both required
#============================================================================.
combine_diag_med <- function(diag, med, diffdays=NULL, method="or"){
  
  if(method=="or"){
    med <- rename(med, diag_date=med_date)
    
    data <- rbind(diag, med) %>%
      arrange(lopnr, diag_date) %>%
      group_by(lopnr) %>%
      filter(row_number() == 1)
  } else if(method=="and" & !is.null(diag) & !is.null(med)){
    data <- full_join(diag, med, relationship = "many-to-many") %>%
      mutate(diag_date_comb = case_when(abs(as.numeric(difftime(med_date, diag_date, units = "days"))) <= diffdays 
                                    ~ pmin(diag_date, med_date))) %>%
      filter(!is.na(diag_date_comb)) %>%
      ungroup() %>%
      select(lopnr, diag_date_comb) %>%
      rename(diag_date=diag_date_comb) %>% 
      arrange(lopnr, diag_date) %>%
      distinct()
  } else {
    data <- NULL
  }
  return(data)
}

#============================================================================.
# Read data ----
#============================================================================.
cohort <- readRDS(file.path("data", "derived", "010_incident_cases.rds"))
cohort$exitdate <- cohort$fu_end

# Read in all diagnoses once and save data for study cohort
oppen_all <- read_sas("K:/EA/SLE Linkage 2022/1. Raw data/Socialstyrelsen update VT2024/ut_r_par_ov_fk_12340_2024.sas7bdat")
names(oppen_all) <- tolower(names(oppen_all))
oppen <- inner_join(select(cohort, lopnr), oppen_all)
saveRDS(oppen, file.path(datapath_derived, paste0(prefix, "_oppen_all_cohort.rds")))
rm(oppen_all)

sluten_all <- read_sas("K:/EA/SLE Linkage 2022/1. Raw data/Socialstyrelsen update VT2024/ut_r_par_sv_fk_12340_2024.sas7bdat")
names(sluten_all) <- tolower(names(sluten_all))
sluten <- inner_join(select(cohort, lopnr), sluten_all)
saveRDS(sluten, file.path(datapath_derived, paste0(prefix, "_sluten_all_cohort.rds")))
rm(sluten_all)

# Read in all medications once and save data for study cohort and selected medications
sel_atc <- c("N05AH02", "N03", "B01AC23", "C04AD03", "A09AA02", "G03XC01",
             "H05AA02", "M05BA", "M05BB", "M05BX03", "M05BX04",
             "G03CA03", "G03CA57", "G03CX01", "G03FA01", "G03FA12", "G03FA15",
             "G03FA17", "G03FB05", "G03FB06", "G03FB09", "A10")

lm_all <- read_sas("K:/EA/SLE Linkage 2022/1. Raw data/Socialstyrelsen update VT2024/ut_r_lmed_fk_12340_2024.sas7bdat")
names(lm_all) <- tolower(names(lm_all))
lm_co <- inner_join(select(cohort, lopnr), lm_all)
lm <- filter_by_code(lm_co, sel_atc, code_column_name = "atc")

saveRDS(lm_co, file.path(datapath_derived, paste0(prefix, "_lm_cohort.rds")))
saveRDS(lm, file.path(datapath_derived, paste0(prefix, "_lm_selected_cohort.rds")))
rm(lm_all, lm_co, lm)

# Read data for study cohort
oppen <- readRDS(file.path(datapath_derived, "012_oppen_all_cohort.rds"))
sluten <- readRDS(file.path(datapath_derived, "012_sluten_all_cohort.rds"))
lm <- readRDS(file.path(datapath_derived, "012_lm_selected_cohort.rds"))

#### 1. Ocular domain ####
#============================================================================.
#### 1.1 Cataract ####
# Diagnoses
cataract_icd <- c("H25", "H26", "H28")
cataract_diag <- diag_date(cohort, oppen, sluten, cataract_icd) 

# Procedures - Currently don't have these codes  
cataract_kva <- c("CJC", "CJD", "CJE")
cataract_proc <- proc_date(cohort, oppen, sluten, cataract_kva)
  
# Combined
cataract <- safe_join(cataract_diag, cataract_proc) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=1.1, diagnosis="Cataract")
  
#### 1.2 Retinal change or atrophy ####
retina_icd <- c("H32", "H34", "H350", "H352", "H353", "H354", "H356", "H357", "H358", "H359", "H36", "H472", "H480", "H481")
retina <- diag_date(cohort, oppen, sluten, retina_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=1.2, diagnosis="Retina")

# Ocular domain ----
ocular <- rbind(cataract, retina) %>%
  arrange(lopnr, desc(incident), diag_date) %>%
  group_by(lopnr) %>%
  mutate(domain_code = 1,
         domain = "Ocular")

#### 2. Neuropsychiatric domain ####
#============================================================================.
#### 2.1 Cognitive impairment or psychosis ####

# Intellectual disability (exclusion)
intellectual_dis_excl_icd <- paste0("F", seq(70,79))
intellectual_dis_excl <- diag_date(cohort, oppen, sluten, intellectual_dis_excl_icd) %>%
  rename(diag_date_excl = diag_date)

# 2.1.1 Cognitive impairment ----
cognitive_imp_icd <- c("F00", "F01", "F02", "F03", "F04", "F051", "G30", "G311", "G311", "G318")
cognitive_imp_diag <- diag_date(cohort, oppen, sluten, cognitive_imp_icd)

# Exclude if any intellectual disability at any time
cognitive_imp <- filter(cognitive_imp_diag, !(lopnr %in% intellectual_dis_excl$lopnr))

# 2.1.2 Psychosis ---- 
# Definition requires clinical diagnosis or prescription
psychosis_icd <- c("F20", "F22", "F25")
psychosis_diag <- diag_date(cohort, oppen, sluten, psychosis_icd) 

# Exclude if any intellectual disability at any time
psychosis_use <- filter(psychosis_diag, !(lopnr %in% intellectual_dis_excl$lopnr))

# Medication
psychosis_atc <- c("N05AH02")
psychosis_med <- med_date(cohort, lm, psychosis_atc, span = 180, gap = 2*365) %>%
  rename(diag_date=med_date) %>%
  select(lopnr, diag_date)

# Combine: Diagnosis or medication
psychosis <- full_join(psychosis_use, psychosis_med)

# 2.1 Cognitive impairment OR psychosis ----
cogn_imp_psychosis <- full_join(cognitive_imp, psychosis) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=2.1, diagnosis="Cognitive impairment or psychosis")

#### 2.2 Epilepsy/Seizures ####
# Definition requires clinical diagnosis and prescription
epilepsy_icd <- c("G40", "G41")
epilepsy_diag <- diag_date(cohort, oppen, sluten, epilepsy_icd)
epilepsy_atc <- c("N03")
epilepsy_med <- med_date(cohort, lm, epilepsy_atc, span = 180, gap = 2*365)
epilepsy <- combine_diag_med(epilepsy_diag, epilepsy_med, diffdays=2*365, method="and") %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=2.2, diagnosis="Epilepsy")

#### 2.3 CVA ####
stroke_icd <- c("I60", "I61", "I63", "I64", "G45")
stroke <- diag_date(cohort, oppen, sluten, stroke_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=2.3, diagnosis="Stroke")

#### 2.4 Cranial or peripheral neuropathy ####
neuropathy_icd <- c("E114", "G50", "G51", "G52", "G53", "G572", "G573", "G574", "G587", "G588", "G589", "G59", "G61", "G62", "G63", "G64")
neuropathy <- diag_date(cohort, oppen, sluten, neuropathy_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=2.4, diagnosis="Neuropathy")

#### 2.5 Transverse myelitis ####
transverse_myelitis_icd <- c("G373")
transverse_myelitis <- diag_date(cohort, oppen, sluten, transverse_myelitis_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=2.5, diagnosis="Transverse myelitis")

# Neuropsychiatric domain ----
neuropsychiatric <- rbind(cogn_imp_psychosis, epilepsy, stroke, neuropathy, transverse_myelitis) %>%
  arrange(lopnr, desc(incident), diag_date) %>%
  group_by(lopnr) %>%
  mutate(domain_code = 2,
         domain = "Neuropsychiatric")

#### 3. Renal domain ####
#============================================================================.
#### 3.1 GFR <50% ####
gfr_lt50_icd <- c("N184")
gfr_lt50 <- diag_date(cohort, oppen, sluten, gfr_lt50_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = 1, diagnosis_code=3.1, diagnosis="gfr lt50")

#### 3.2 Proteinuria 24 h >3.5 g OR End-stage renal disease ####
# 3.2.1 ESRD 
esrd_icd <- c("N185")
esrd <- diag_date(cohort, oppen, sluten, esrd_icd)

# 3.2.2 Renal dialysis
# Diagnosis
renal_dialysis_icd <- c("Z992")
renal_dialysis_diag <- diag_date(cohort, oppen, sluten, renal_dialysis_icd)

# Procedure
renal_dialysis_kva <- c("DR016", "DR024", "DR061")
renal_dialysis_proc <- proc_date(cohort, oppen, sluten, renal_dialysis_kva)

# Combine
renal_dialysis <- safe_join(renal_dialysis_diag, renal_dialysis_proc, by="lopnr", join_type="full") %>%
  arrange(lopnr, diag_date) %>%
  group_by(lopnr) %>%
  filter(row_number() == 1)

# Renal tx (KVA)
renal_tx_icd <- c("Z940")
renal_tx_diag <- diag_date(cohort, oppen, sluten, renal_tx_icd)
renal_tx_kva <- c("KAS00", "KAS10", "KAS20")
renal_tx_proc <- proc_date(cohort, oppen, sluten, renal_tx_kva)
renal_tx <- safe_join(renal_tx_diag, renal_tx_proc) %>%
  select_1st_diag_date()

# Score set to 3 as in SDI
prot_renal <- rbind(esrd, renal_dialysis, renal_tx) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = 3, diagnosis_code=3.2, diagnosis="Proteinuria or end-stage renal disease")

# Renal domain ----
renal <- rbind(gfr_lt50, prot_renal) %>%
  arrange(lopnr, desc(incident), diag_date) %>%
  group_by(lopnr) %>%
  mutate(domain_code = 3,
         domain = "Renal")

# 4. Pulmonary domain ----
#============================================================================.
#### 4.1 Pulmonary hypertension ####
pah_icd <- c("I270", "I272", "I279")
pah <- diag_date(cohort, oppen, sluten, pah_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=4.1, diagnosis="Pulmonary hypertension")

#### 4.2 Pulmonary fibrosis ####
lung_fibrosis_icd <- c("J701", "J703", "J841")
lung_fibrosis <- diag_date(cohort, oppen, sluten, lung_fibrosis_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=4.2, diagnosis="Pulmonary fibrosis")

#### 4.3 Shrinking lung ####
# Requires X-ray 
shrink_lung <- NULL

#### 4.4 Pleural fibrosis ####
pleural_fib_icd <- c("J92", "J941")
pleural_fib <- diag_date(cohort, oppen, sluten, pleural_fib_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=4.4, diagnosis="Pleural fibrosis")

#### 4.5 Pulmonary infarction ####
pinf_icd <- c("I26")
pinf_diag <- diag_date(cohort, oppen, sluten, pinf_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=4.5, diagnosis="Pulmonary infarction")

# Missing pinf DRG codes 
pinf_drg <- c("GDB")
pinf_d <- drg_date(cohort, oppen, sluten, pinf_drg)
pinf <- rbind(pinf_diag, pinf_d)

# Pulmonary domain ----
pulmonary <- rbind(pah, lung_fibrosis, shrink_lung, pleural_fib, pinf) %>%
  arrange(lopnr, desc(incident), diag_date) %>%
  group_by(lopnr) %>%
  mutate(domain_code = 4,
         domain = "Pulmonary")

# 5. Cardiovascular domain ----
#============================================================================.
#### 5.1 Angina or CABG ####

# 5.1.1 Angina
angina_icd1 <- c("I25") # Chronic
angina_icd2 <- c("I20", "I24") # [Two or more codes, at least six months apart]

angina1 <- diag_date(cohort, oppen, sluten, angina_icd1) 
angina2 <- diag_date(cohort, oppen, sluten, angina_icd2) %>%
  group_by(lopnr) %>% 
  arrange(lopnr, diag_date) %>% 
  mutate(last2 = difftime(diag_date, lag(diag_date), units = "days"),
         diag_gap = ifelse(last2 >= 2*365, 1, 0),
         newseq = ifelse(diag_gap ==0 | is.na(diag_gap),0,1),
         seq = cumsum(newseq)) %>%
  group_by(lopnr, seq) %>%
  mutate(diag_date_1st = as.Date(ifelse(is.infinite(suppressWarnings(min(diag_date, na.rm = TRUE))), 
                                        NA, suppressWarnings(min(diag_date, na.rm = TRUE))), origin = "1970-01-01"),
         diag_span = difftime(diag_date, diag_date_1st, units = "days"),
         diag_span_seq = ifelse(is.infinite(suppressWarnings(max(diag_span, na.rm=TRUE))), 
                                NA, suppressWarnings(max(diag_span, na.rm=TRUE))), 
         seq_sel = ifelse(diag_span_seq >= 180, 1, 0)
  ) %>%
  filter(seq_sel == 1) %>%
  ungroup() %>%
  select(lopnr, diag_date_1st) %>%
  rename(diag_date = diag_date_1st) %>%
  distinct()

angina <- rbind(angina1, angina2)

# 5.1.2 CABG
cabg_icd <- c("Z951", "Z955")
cabg_diag <- diag_date(cohort, oppen, sluten, cabg_icd) 

cabg_drg <- c("FNA", "FNB", "FNC", "FND", "FNE", "FNEG")
cabg_d <- drg_date(cohort, oppen, sluten, cabg_drg) 

# CABG not detected
cabg <- rbind(cabg_diag, cabg_d)

# Angina or CABG
angina_cabg <- rbind(angina, cabg) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=5.1, diagnosis="Angina or CABG")

#### 5.2 MI ####
# Myocardial infarction ever (score 2 if > 1) 

# 5.2.1 1st MI episode
mi1_icd <- c("I21")
mi1 <- diag_date(cohort, oppen, sluten, mi1_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate)) %>%
  post_pre(lookback = def_lookback) %>%
  add_score_post_pre()

# 5.2.2 MI 2nd episode ####
mi2_icd <- c("I22")
mi2 <- diag_date(cohort, oppen, sluten, mi2_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate)) %>%
  post_pre(lookback = def_lookback) %>%
  add_score_post_pre()

# Combine and set score 2 if > 1, keeping all dates when increased score post baseline, and single summary score at bl
mi <- rbind(mi1, mi2) %>%
  inner_join(select(cohort, lopnr, indexdate, exitdate)) %>%
  arrange(lopnr, diag_date) %>%
  group_by(lopnr) %>%
  mutate(tmp1_score=cumsum(tmp_score),
         score = ifelse(tmp1_score >1,2,tmp1_score),
         tmp1_score_bl=cumsum(tmp_score_bl),
         score_bl = ifelse(tmp1_score_bl >1,2,tmp1_score_bl)) %>%
  select(-starts_with("tmp")) %>%
  filter(row_number() <= 2) %>%
  mutate(diagnosis_code = 5.2,
         diagnosis = "Myocardial infarction")

#### 5.3 Cardiomyopathy ####
cardiomyopathy_icd <- c("I42") # Exclude I255 from original definition since part of angina definition
cardiomyopathy <- diag_date(cohort, oppen, sluten, cardiomyopathy_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=5.3, diagnosis="Cardiomyopathy")

#### 5.4 Valvular disease ####
valvulopathy_icd <- c("I05", "I06", "I07", "I08", "I34", "I35", "I36", "I37", "Z952", "Z953", "Z954")
valvulopathy_drg <- c("FG", "FK", "FM", "FJE", "FJF")
valvulopathy <- diag_date(cohort, oppen, sluten, valvulopathy_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=5.4, diagnosis="Valvulopathy")

#### 5.5 Pericarditis ####
pericarditis_icd <- c("I092", "I310", "I311")
pericarditis_drg <- c("FEB", "FEF")
pericarditis <- diag_date(cohort, oppen, sluten, pericarditis_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=5.5, diagnosis="Pericarditis")

# Cardiovascular domain ----
cardiovascular <- rbind(angina_cabg, mi, cardiomyopathy, valvulopathy, pericarditis) %>%
  arrange(lopnr, desc(incident), diag_date) %>%
  group_by(lopnr) %>%
  mutate(domain_code = 5,
         domain = "Cardiovascular")

# 6. Peripheral vascular ----
#============================================================================.
#### 6.1 Claudication ####
claudication_icd <- c("I739B")
claudication_diag <- diag_date(cohort, oppen, sluten, claudication_icd) 
claudication_atc <- c("B01AC23", "C04AD03")
claudication_med <- med_date(cohort, lm, claudication_atc, span = 180, gap = 2*365) %>%
  rename(diag_date=med_date) %>%
  select(lopnr, diag_date)

# Combine: Diagnosis or medication
claudication <- full_join(claudication_diag, claudication_med) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=6.1, diagnosis="Claudication")

#### 6.2 Minor tissue loss ####
tisloss_minor_icd <- c("I739C")
tisloss_minor <- diag_date(cohort, oppen, sluten, tisloss_minor_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=6.2, diagnosis="Minor tissue loss")

#### 6.3 Significant tissue loss ####
# Significant tissue loss ever (e.g., loss of digit or limb, resection) (score 2 if > 1)

tisloss_major_kva <- c("NHQ")
tisloss_major_proc <- proc_date(cohort, oppen, sluten, tisloss_major_kva) 
tisloss_major_drg <- c("E20", "113", "114", "E21", "H09", "L01")
tisloss_major_d <- drg_date(cohort, oppen, sluten, tisloss_major_drg)

# Combine
tisloss_major_all <- safe_join(tisloss_major_proc, tisloss_major_d) %>%
  select(lopnr, diag_date) %>%
  distinct() %>%
  inner_join(select(cohort, lopnr, indexdate, exitdate))

# To set score 2 if > 1 - select first two diagnoses after indexdate
df_post <- tisloss_major_all %>% 
    filter(diag_date >= indexdate, diag_date <= exitdate) %>%
    group_by(lopnr) %>%
  filter(row_number() <= 2)

df_pre <- tisloss_major_all %>% 
  filter(diag_date >= indexdate - def_lookback, diag_date < indexdate) %>%
  arrange(lopnr, desc(diag_date)) %>%
  mutate(tmp_score=1) %>%
  group_by(lopnr) %>%
    mutate(tmp1_score=sum(tmp_score),
           score_bl = ifelse(tmp1_score >1,2,tmp1_score)) %>%
    select(-starts_with("tmp")) %>%
    filter(row_number() == 1) %>%
  rename(diag_date_pre = diag_date)


tisloss_major <- safe_join(df_post, df_pre, by=c("lopnr", "indexdate", "exitdate"), join_type = "full") %>%
  mutate(incident = ifelse(!is.na(diag_date) & is.na(diag_date_pre),1,0),
         bl_prevalent = ifelse(!is.na(diag_date_pre),1,0),
         diagnosis_code = 6.3 + 0.01*row_number(),
         diagnosis = "Significant tissue loss",
         tmp_score = ifelse(incident==1,1,0)) %>%
  group_by(lopnr) %>%
  mutate(tmp1_score=cumsum(tmp_score),
         score = ifelse(tmp1_score >1,2,tmp1_score)) %>%
  select(-starts_with("tmp"))


#### 6.4 DVT ####
vt_icd <- c("I80", "I81", "I82")
vt <- diag_date(cohort, oppen, sluten, vt_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=6.4, diagnosis="Venous thrombosis")

# Peripheral vascular domain ----
peripheral_vascular <- rbind(claudication, tisloss_minor, tisloss_major, vt) %>%
  arrange(lopnr, desc(incident), diag_date) %>%
  group_by(lopnr) %>%
  mutate(domain_code = 6,
         domain = "Peripheral vascular")

# 7. Gastrointestinal ----
#============================================================================.
#### 7.1 Infarction or resection bowel, spleen, lever or gallbladder ####
# 7.1.1 Diagnoses
gi_bowel_icd <- c("K550")
gi_bowel_diag <- diag_date(cohort, oppen, sluten, gi_bowel_icd)
gi_spleen_icd <- c("D735")
gi_spleen_diag <- diag_date(cohort, oppen, sluten, gi_spleen_icd)
gi_liver_icd <- c("K763")
gi_liver_diag <- diag_date(cohort, oppen, sluten, gi_liver_icd)

# 7.1.2 Procedures
gi_bowel_kva <- c("JFB", "JFH")
gi_bowel_proc <- proc_date(cohort, oppen, sluten, gi_bowel_kva)
gi_spleen_kva <- c("JMA")
gi_spleen_proc <- proc_date(cohort, oppen, sluten, gi_spleen_kva)
gi_liver_kva <- c("JJB")
gi_liver_proc <- proc_date(cohort, oppen, sluten, gi_liver_kva)
gi_gallbladder_kva <- c("JKA20", "JKA21")
gi_gallbladder_proc <- proc_date(cohort, oppen, sluten, gi_gallbladder_kva)

# 7.1.3 DRG codes
gi_bowel_drg <- c("F01", "146", "147", "F05", "148", "149", "F09", "152", "153", "F30", "166", "167")
gi_bowel_d <- drg_date(cohort, oppen, sluten, gi_bowel_drg)
gi_spleen_drg <- c("R01N", "392", "393O")
gi_spleen_d <- drg_date(cohort, oppen, sluten, gi_spleen_drg)
gi_liver_drg <- c("G01N", "480")
gi_liver_d <- drg_date(cohort, oppen, sluten, gi_liver_drg)
gi_gallbladder_drg <- c("G10", "196", "196", "G11", "197", "198", "G12", "493", "494", "G3O", "494O")
gi_gallbladder_d <- drg_date(cohort, oppen, sluten, gi_gallbladder_drg)

# Bowel
gi_bowel <- rbind(gi_bowel_diag, gi_bowel_proc, gi_bowel_d) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type="right") %>%
  post_pre(lookback = def_lookback) %>%
  add_score_post_pre()

# Spleen
gi_spleen <- rbind(gi_spleen_diag, gi_spleen_proc, gi_spleen_d) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate)) %>%
  post_pre(lookback = def_lookback) %>%
  add_score_post_pre()

# Liver
gi_liver <- rbind(gi_liver_diag, gi_liver_proc, gi_liver_d) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate)) %>%
  post_pre(lookback = def_lookback) %>%
  add_score_post_pre()

# Gallbladder
gi_gallbladder <- rbind(gi_gallbladder_proc, gi_gallbladder_d) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate)) %>%
  post_pre(lookback = def_lookback) %>%
  add_score_post_pre()

# Combine
gi <- rbind(gi_bowel, gi_spleen, gi_liver, gi_gallbladder) %>%
  inner_join(select(cohort, lopnr, indexdate, exitdate)) %>%
  arrange(lopnr, diag_date) %>%
  group_by(lopnr) %>%
  mutate(tmp1_score=cumsum(tmp_score),
         score = ifelse(tmp1_score >1,2,tmp1_score),
         tmp1_score_bl=cumsum(tmp_score_bl),
         score_bl = ifelse(tmp1_score_bl >1,2,tmp1_score_bl)) %>%
  select(-starts_with("tmp")) %>%
  filter(row_number() <= 2) %>%
  mutate(diagnosis_code = 7.1,
         diagnosis = "GI infarction or resection")

#### 7.2 Mesenteric insufficiency ####
mesen_icd  <- c("K551")
mesen <- diag_date(cohort, oppen, sluten, mesen_icd)

#### 7.3 Chronic peritonitis ####
peritonitis_icd <- c("K658", "N734")
peritonitis <- diag_date(cohort, oppen, sluten, peritonitis_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=7.3, diagnosis="Chronic peritonitis")

#### 7.4 Stricture or upper GIT surgery ####
uppergi_icd <- c("K222")
uppergi_diag <- diag_date(cohort, oppen, sluten, uppergi_icd) 
uppergi_kva <- c("JC", "JD")
uppergi_proc <- proc_date(cohort, oppen, sluten, uppergi_kva)
uppergi_drg <- c("L08", "288", "F11", "154", "155", "156", "F12", "F13O", "F35", "170", "171")
uppergi_d <- drg_date(cohort, oppen, sluten, uppergi_drg)
uppergi <- safe_join(uppergi_diag, uppergi_proc, uppergi_d) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=7.4, diagnosis="Stricture or upper GIT surgery")

#### 7.5 Pancreatic insufficiency ####
# Definition requires clinical diagnosis or prescription + K86.8, K86.9, K90.3, E16.9
# Diagnosis only
pancreatic_icd <- c("K860", "K861", "K863", "K868", "K869", "E169")
pancreatic_diag <- diag_date(cohort, oppen, sluten, pancreatic_icd)

# Prescription + any of selected diagnoses
pancreatic_atc <- c("A09AA02")
pancreatic_med <-  filter(lm, atc %in% pancreatic_atc) %>%
  select(lopnr, edatum) %>%
  rename(med_date = edatum) %>%
  distinct() %>%
  inner_join(select(cohort, lopnr))

pancreatic_sub_icd <- c("K868", "K869", "K903", "E169")
pancreatic_diag_sub <- diag_date(cohort, oppen, sluten, pancreatic_sub_icd)
pancreatic_comb <- combine_diag_med(pancreatic_diag_sub, pancreatic_med, diffdays=NULL, method="and") 

pancreas <- safe_join(pancreatic_diag, pancreatic_comb, by="lopnr", join_type="full") %>%
  select_1st_diag_date()

if(!is.null(pancreas)){
  pancreas <- pancreas %>%
    mutate(domain_code = 7.5, domain = "Pancreas")
}

# Gastrointestinal domain ----
gastrointestinal <- rbind(gi, mesen, peritonitis, uppergi, pancreas) %>%
  arrange(lopnr, desc(incident), diag_date) %>%
  group_by(lopnr) %>%
  mutate(domain_code = 7,
         domain = "Gastrointestinal")

# 8. Musculoskeletal ----
#============================================================================.
#### 8.1 Atrophy or weakness ####
msk_atrophy_icd <- c("M625")
msk_atrophy <- diag_date(cohort, oppen, sluten, msk_atrophy_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=8.1, diagnosis="Atrophy or weakness")

#### 8.2 Deforming or erosive arthritis ####
arthritis_icd <- c("M120", "M058L", "M058M", "M058N", "M068L", "M068M", "M068N")
arthritis <- diag_date(cohort, oppen, sluten, arthritis_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=8.2, diagnosis="Deforming or erosive arthritis")

#### 8.3 Osteoporosis with fracture or vertebral collapse ####
# Definition requires a code for osteoporosis and at least one code for
# non-vertebral fracture or vertebral collapse.
# ICD code M80-M82 or atc + fracture 

# Osteoporosis
osteoporosis_icd <- c("M80", "M81", "M82")
osteoporosis_diag <- diag_date(cohort, oppen, sluten, osteoporosis_icd) 
osteoporosis_atc <- c("G03XC01", "H05AA02", "M05BA", "M05BB", "M05BX03", "M05BX04")
osteoporosis_med <- med_date(cohort, lm, osteoporosis_atc, span = 180, gap = 2*365)
osteoporosis <- combine_diag_med(osteoporosis_diag, osteoporosis_med, diffdays=NULL, method="or") 

# Fracture
fx_icd <- c("S12", "S22", "S32", "S42", "S52", "S62", "S72", "S82", "S92", "T02", "T08", "T10", "T12", "T142")
fx <- diag_date(cohort, oppen, sluten, fx_icd)

# Vertebral collapse
vertebral_col_icd <- c("M485", "S120", "S121", "S122", "S220", "S221", "S320")
vertebral_col <- diag_date(cohort, oppen, sluten, vertebral_col_icd)

# Merge fracture with vertebral collapse
fxvercol <- rbind(fx, vertebral_col) %>%
  arrange(lopnr, diag_date) %>%
  filter(row_number() == 1) %>%
  rename(frac_date = diag_date)

# Merge osteoporosis with fracture or vertebral collapse
ost_frac <- safe_join(osteoporosis, fxvercol) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  filter(!is.na(diag_date) & !is.na(frac_date)) %>%
  rename(ost_date=diag_date) %>%
  group_by(lopnr) %>%
  mutate(diag_date = case_when(frac_date >= indexdate ~ min(frac_date),
                               frac_date < indexdate ~ max(frac_date))) %>% 
  select(-ost_date, -frac_date) %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=8.3, diagnosis="Osteoporosis with fracture or vertebral collapse")

#### 8.4 Avascular necrosis ####
# Not possible to score 2 based on codes we have
avascular_nec_icd <- c("M870", "M873", "M878", "M879")
avascular_nec <- diag_date(cohort, oppen, sluten, avascular_nec_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=8.4, diagnosis="Avascular necrosis")

#### 8.5 Osteomyelitis ####
osteomyelitis_icd <- c("M86")
osteomyelitis <- diag_date(cohort, oppen, sluten, osteomyelitis_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=8.5, diagnosis="Osteomyelitis")

#### 8.6 Ruptured tendons ####
ruptured_tendons_icd <- c("M662", "M663", "M664", "M665")
ruptured_tendons <- diag_date(cohort, oppen, sluten, ruptured_tendons_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=8.6, diagnosis="Ruptured tendons")

# Musculoskeletal domain ----
musculoskeletal <- rbind(msk_atrophy, arthritis, ost_frac, avascular_nec, osteomyelitis, ruptured_tendons) %>%
  arrange(lopnr, desc(incident), diag_date) %>%
  group_by(lopnr) %>%
  mutate(domain_code = 8,
         domain = "Musculoskeletal")

# 9. Skin ----
#============================================================================.
#### 9.1 Alopecia ####
alopecia_icd <- c("L66")
alopecia <- diag_date(cohort, oppen, sluten, alopecia_icd) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=9.1, diagnosis="Alopecia")

#### 9.2 Scarring or paniculum ####
scarring_icd <- c("L905")
scarring <- diag_date(cohort, oppen, sluten, scarring_icd)

#### 9.3 Skin ulceration ####
skin_ulcer_icd1 <- c("L984")
skin_ulcer_icd2 <- c("E106D", "E116D", "I702C", "I830", "I832", "L89", "L97") # [Two or more codes, at least six months apart]

skin_ulcer1 <- diag_date(cohort, oppen, sluten, skin_ulcer_icd1) 
skin_ulcer2 <- diag_date(cohort, oppen, sluten, skin_ulcer_icd2) %>%
  group_by(lopnr) %>% 
  arrange(lopnr, diag_date) %>% 
  mutate(last2 = difftime(diag_date, lag(diag_date), units = "days"),
       diag_gap = ifelse(last2 >= 2*365, 1, 0),
       newseq = ifelse(diag_gap ==0 | is.na(diag_gap),0,1),
       seq = cumsum(newseq)) %>%
  group_by(lopnr, seq) %>%
  mutate(diag_date_1st = as.Date(ifelse(is.infinite(suppressWarnings(min(diag_date, na.rm = TRUE))), 
                                       NA, suppressWarnings(min(diag_date, na.rm = TRUE))), origin = "1970-01-01"),
         diag_span = difftime(diag_date, diag_date_1st, units = "days"),
         diag_span_seq = ifelse(is.infinite(suppressWarnings(max(diag_span, na.rm=TRUE))), 
                               NA, suppressWarnings(max(diag_span, na.rm=TRUE))), 
         seq_sel = ifelse(diag_span_seq >= 180, 1, 0)
  ) %>%
  filter(seq_sel == 1) %>%
  ungroup() %>%
  select(lopnr, diag_date_1st) %>%
  rename(diag_date = diag_date_1st) %>%
  distinct()

skin_ulcer <- rbind(skin_ulcer1, skin_ulcer2) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=9.3, diagnosis="Skin ulceration")
  
# Skin domain ----
skin <- rbind(alopecia, scarring, skin_ulcer) %>%
  arrange(lopnr, desc(incident), diag_date) %>%
  group_by(lopnr) %>%
  mutate(domain_code = 9,
         domain = "Skin")

# 10. Premature gonadal failure ----
#============================================================================.
#### 10.1 Premature gonadal failure ####
# Maximum age 40 at first prescription of continuous prescriptions:
# [Two or more prescriptions, at least six months apart]

poi_icd <- c("E283", "E894")
poi_diag <- diag_date(cohort, oppen, sluten, poi_icd)
poi_atc <- c("G03CA03", "G03CA57", "G03CX01", "G03FA01", "G03FA12", "G03FA15", "G03FA17", "G03FB05", "G03FB06", "G03FB09")
poi_med <- med_date(cohort, lm, poi_atc, span = 180, gap = 2*365) %>%
  inner_join(select(cohort, lopnr, birthdate)) %>%
  filter(as.numeric(difftime(med_date, birthdate, units = "weeks"))/52 < 40) %>%
  select(-birthdate)

# Premature gonadal failure domain ----
poi <- combine_diag_med(poi_diag, poi_med, diffdays=365, method="or") %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=10.1, diagnosis="Premature gonadal failure") %>%
  arrange(lopnr, desc(incident), diag_date) %>%
  group_by(lopnr) %>%
  mutate(domain_code = 10,
         domain = "Premature gonadal failure")

# 11. Diabetes ----
#============================================================================.
#### 11.1 Diabetes ####
# [Two or more prescriptions, at least six months apart]
diabetes_icd <- c("E10", "E11", "E12", "E13", "E14")
diabetes_diag <- diag_date(cohort, oppen, sluten, diabetes_icd)%>%
  arrange(lopnr, diag_date)
diabetes_atc <- c("A10")
diabetes_med <- med_date(cohort, lm, diabetes_atc, span = 180, gap = 2*365) %>%
  arrange(lopnr, med_date)

# Diabetes domain ----
diabetes <- combine_diag_med(diabetes_diag, diabetes_med, diffdays=NULL, method="or") %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=11.1, diagnosis="Diabetes") %>%
  arrange(lopnr, desc(incident), diag_date) %>%
  group_by(lopnr) %>%
  mutate(domain_code = 11,
         domain = "Diabetes")

# 12. Malignancy ----
#============================================================================.
#### 12.1 Malignancy ####
## Malignant tumor
malign_icd <- c(sprintf("C%02d", 0:43), sprintf("C%02d", 45:97), sprintf("D%02d", 0:9))
malign <- diag_date(cohort, oppen, sluten, malign_icd) 

## Carcinoma in situ
cainsitu_icd <- sprintf("D%02d", 0:09)
cainsitu <- diag_date(cohort, oppen, sluten, cainsitu_icd)

## Combine
malignancy <- rbind(malign, cainsitu) %>%
  safe_join(select(cohort, lopnr, indexdate, exitdate), join_type = "inner") %>%
  select_post_pre(lookback = def_lookback, set_score = NULL, diagnosis_code=12.1, diagnosis="Malignancy") %>%
  arrange(lopnr, desc(incident), diag_date) %>%
  group_by(lopnr) %>%
  mutate(domain_code = 12,
         domain = "Malignancy")

# Combine all domains and set score if not previously set ----
#============================================================================.
df_all <- list(ocular, neuropsychiatric, renal, pulmonary, cardiovascular, 
                peripheral_vascular, gastrointestinal, musculoskeletal, 
                skin, poi, diabetes, malignancy) %>% 
  reduce(function(x, y) full_join(x, y, by = intersect(names(x), names(y)))) %>%
  ungroup() %>%
  arrange(lopnr, diagnosis_code)

if(!"score" %in% names(df_all)) df_all$score <- as.numeric(NA)

df_all2 <- df_all %>%
  ungroup() %>%
  right_join(select(cohort, lopnr)) %>%
  mutate(score2 = case_when(is.na(score) & incident==1 & bl_prevalent==0 ~ 1,
                            !is.na(score) & incident==1 & bl_prevalent==0 ~ score,
                            TRUE ~ 0),
         score2_bl = case_when(is.na(score_bl) & bl_prevalent==1 ~ 1,
                           !is.na(score_bl) ~ score_bl,
                           TRUE ~ 0)) %>%
  group_by(lopnr, domain_code) %>%
  arrange(lopnr, domain_code, diag_date) %>%
  mutate(domain_score = cumsum(score2),
         domain_score_bl = max(score2_bl, na.rm=TRUE),
         # Renal domain can be maximum score=3
         domain_score = ifelse(domain_score>3,3, domain_score), 
         domain_score_bl = ifelse(domain_score_bl>3,3, domain_score_bl)) 

saveRDS(df_all2, file.path(datapath_derived, paste0(prefix, "_rbodi_all_diagnoses_scores.rds")))

df_all3 <- df_all2[,!grepl("score", names(df_all2))]
saveRDS(df_all3, file.path(datapath_derived, paste0(prefix, "_rbodi_all_diagnoses.rds")))

# Overall RBODI score over time and baseline damage score ----
#============================================================================.

# Time points when increase in RBODI (after indexdate, RBODI at indexdate is 0)
df_rbodi <- df_all2 %>%
  filter(incident==1) %>%
  group_by(lopnr, domain_code, domain_score) %>%
  filter(row_number()==1) %>%
  group_by(lopnr) %>%
  arrange(lopnr, diag_date) %>%
  mutate(rbodi = cumsum(domain_score)) %>%
  filter(rbodi>0) %>%
  select(lopnr, diag_date, rbodi)

# RBODI=0 at indexdate
dfi <- select(cohort, lopnr, indexdate) %>%
  rename(diag_date = indexdate) %>%
  mutate(rbodi = 0)

# RBODI values from indexdate and onwards
df_rbodi2 <- rbind(df_rbodi, dfi)

# Comorbidity/organ damage score at baseline (based on RBODI scoring)
df_bl <- df_all2 %>%
  group_by(lopnr, domain_code) %>%
  arrange(lopnr, domain_code, diag_date) %>%
  filter(row_number()==1) %>%
  group_by(lopnr) %>%
  summarise(rbodi_bl= sum(domain_score_bl))

# Time points of RBODI changes and score at baseline
df_rbodi3 <- full_join(select(cohort, lopnr, indexdate, fu_end), df_bl) %>%
  full_join(df_rbodi2, by = "lopnr") %>%
  mutate(rbodi = ifelse(is.na(rbodi),0,rbodi),
         rbodi_bl = ifelse(is.na(rbodi_bl),0,rbodi_bl),
         time = as.numeric(difftime(diag_date, indexdate, units = "days"))/def_dpyr) %>%
  arrange(lopnr, diag_date)

saveRDS(df_rbodi3, file.path(datapath_derived, paste0(prefix, "_rbodi.rds")))
