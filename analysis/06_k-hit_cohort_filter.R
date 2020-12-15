# restrict case cohort and their matching controls if they have at least k schizo entries
# now that all potential case icds are saved on O2

start_time <- Sys.time()

library(dplyr)
library(tidyr)
library(stringr)
library(vroom)
library(lubridate)

k=3

## load into dataframes
# update read.csv functions to vroom versions to make it more efficient 
df_match <- vroom('XXX/schizo/data/schizo_match_BirthyearGenderZip.csv',col_select=c('case_id','schizo_onset','case_XXXXXFirstDate', 'case_XXXXXLastDate','ctrl_id','ctrl_XXXXXFirstDate', 'ctrl_XXXXXLastDate'))

# get all case icds
case_icds <- vroom('XXX/schizo/data/retrieved_data/case_ctrl_icds/case_icds_full.txt')

# obtained from phecode in 295 or 295.1, and exclude "V11.0 history of schizophrenia" 
schizo_icd_codes <- c('293.81', '293.82', '295', '295.0', '295.00', '295.01', '295.02', '295.03', '295.04', '295.05', '295.1', '295.10', '295.11', '295.12', '295.13', '295.14', '295.15', '295.2', '295.20', '295.21', '295.22', '295.23', '295.24', '295.25', '295.3', '295.30', '295.31', '295.32', '295.33', '295.34', '295.35', '295.4', '295.40', '295.41', '295.42', '295.43', '295.44', '295.45', '295.5', '295.50', '295.51', '295.52', '295.53', '295.54', '295.55', '295.6', '295.60', '295.61', '295.62', '295.63', '295.64', '295.65', '295.7', '295.70', '295.71', '295.72', '295.73', '295.74', '295.75', '295.8', '295.80', '295.81', '295.82', '295.83', '295.84', '295.85', '295.9', '295.90', '295.91', '295.92', '295.93', '295.94', '295.95', '298', '298.0', '298.1', '298.2', '298.8', '298.9', '298.90', 'F25', 'F25.9', 'F25.0', 'F25.1', 'F20', 'F20.9', 'F20.2', 'F20.5', 'F20.81', 'F20.89', 'F20.1', 'F20.0', 'F44.89', 'F25.8', 'F20.8', 'F20.3')

# filter out only the rows of schizo mentions 
case_schizo_icds <- filter(case_icds, Icd %in% schizo_icd_codes)
case_schizo_k <- case_schizo_icds %>%
    group_by(case_id) %>%
    filter(n() >= k) %>%
    distinct(case_id)

# only keep the pairs where case had at least three schizo icds
df_match_k <- merge(case_schizo_k, df_match, by='case_id')
write.csv(df_match_k, paste('XXX/schizo/data/schizo_match_BirthyearGenderZip_',k,'hit.csv',sep=''))

# get denominator - see who were enrolled at all the times (only need once for all phecodes)
case_all <- df_match_k[c('case_id','schizo_onset','case_XXXXXFirstDate')]
for(i in 0:8){
    print(i)
    case_all[ , paste0("time-", i)] <- (case_all$schizo_onset%m-%months(6*i)) >= case_all$case_XXXXXFirstDate
}

case_enrollment <- case_all[,4:12] %>%
    summarise_all(sum)

ctrl_all <- df_match_k[c('ctrl_id','schizo_onset','ctrl_XXXXXFirstDate')]
for(i in 0:8){
    print(i)
    ctrl_all[ , paste0("time-", i)] <- (ctrl_all$schizo_onset%m-%months(6*i)) >= ctrl_all$ctrl_XXXXXFirstDate
}

ctrl_enrollment <- ctrl_all[,4:12] %>%
    summarise_all(sum)
