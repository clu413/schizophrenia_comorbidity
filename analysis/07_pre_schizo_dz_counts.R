# for each phecode disease group, count the number of patients who had this disease in the schizo group vs their non-schizo matched controls before the schizo index dates
# only use 2-hit patients' icds

start_time <- Sys.time()

library(dplyr)
library(tidyr)
library(stringr)
library(vroom)
library(lubridate)

args <- commandArgs()
print(args)
phe=args[6]

phecode_df <- read.csv('XXX/schizo/data/phecode_update.txt',colClasses = c('PheWASCode'='character'))
df_icd <- phecode_df[(phecode_df$PheWASCode == phe & !str_detect(phecode_df$IcdString, 'istory')),]
icd_list <- unique(df_icd$IcdCode)
remove(phecode_df,df_icd)

## load into dataframes
case_icds_full <- vroom('XXX/schizo/data/retrieved_data/case_ctrl_icds/case_icds_3hit.txt')
print(paste('Finished loading all the cases!', Sys.time()))
# subset only the entries in this phecode 
case_icds_phe <- case_icds_full %>%
    filter(Icd %in% icd_list)
remove(case_icds_full)

ctrl_icds_full <- vroom('XXX/schizo/data/retrieved_data/case_ctrl_icds/ctrl_icds_3hit.txt')
print(paste('Finished loading all the ctrls!', Sys.time()))
ctrl_icds_phe <- ctrl_icds_full %>%
    filter(Icd %in% icd_list)
remove(ctrl_icds_full)

# read in patient info for the ones that meet the 3-hit schizo criteria
df_match <- vroom('XXX/schizo/data/schizo_match_BirthyearGenderZip_3hit.csv',col_select=c('case_id','schizo_onset','case_XXXXXFirstDate', 'case_XXXXXLastDate','ctrl_id','ctrl_XXXXXFirstDate', 'ctrl_XXXXXLastDate'))

# merge the dataframes 
case_df <- merge(df_match[c('case_id','schizo_onset','case_XXXXXFirstDate','case_XXXXXLastDate')], case_icds_phe, by='case_id')
ctrl_df <- merge(df_match[c('ctrl_id','schizo_onset','ctrl_XXXXXFirstDate','ctrl_XXXXXLastDate')], ctrl_icds_phe, by='ctrl_id')

# subset only the entries that were before or on the same day as schizo onset 
case_pre <- case_df %>%
    filter(DateServiceStarted <= schizo_onset)
ctrl_pre <- ctrl_df %>%
    filter(DateServiceStarted <= schizo_onset)
remove(df_match)

# filter out patients who only had one or two mention of dz before schizo onset
case_pre_3hit <- case_pre %>%
    group_by(case_id) %>%
    filter(n() >= 3) %>%
    ungroup()

# get patients who had at least 3 mentions of dz before schizo onset
case_pre_pt_counts_3hit <- length(unique(case_pre_3hit$case_id))

for(i in 0:7){
    case_pre_3hit[ , paste0("time-", i)] <-
        # icd entry occured after six months before timestamp index date
        ((case_pre_3hit$schizo_onset%m-%months(6*(i+1))) < case_pre_3hit$DateServiceStarted) &
        # icd entry occured before or on timestamp index date
        (case_pre_3hit$DateServiceStarted <= (case_pre_3hit$schizo_onset%m-%months(6*i))) #
}

case_dz_tally_3hit <- case_pre_3hit %>%
    group_by(case_id) %>%
    summarise(across('time-0':'time-7', any))

case_counts_3hit = colSums(case_dz_tally_3hit[,2:9])
print(paste('For phecode',phe,'that had three hits, at each six month period, we have'))
print(case_counts_3hit)

### repeat the steps for ctrls
ctrl_pre_3hit <- ctrl_pre %>%
    group_by(ctrl_id) %>%
    filter(n() >= 3) %>%
    ungroup()

# get patients who had at least 3 mentions of dz before schizo onset
ctrl_pre_pt_counts_3hit <- length(unique(ctrl_pre_3hit$ctrl_id))

for(i in 0:7){
    ctrl_pre_3hit[ , paste0("time-", i)] <-
        ((ctrl_pre_3hit$schizo_onset%m-%months(6*(i+1))) < ctrl_pre_3hit$DateServiceStarted) &
        (ctrl_pre_3hit$DateServiceStarted <= (ctrl_pre_3hit$schizo_onset%m-%months(6*i))) #
}

ctrl_dz_tally_3hit <- ctrl_pre_3hit %>%
    group_by(ctrl_id) %>%
    summarise(across('time-0':'time-7', any))

ctrl_counts_3hit = colSums(ctrl_dz_tally_3hit[,2:9])
print(paste('For phecode',phe,'that had three hits, at each six month period, we have'))
print(ctrl_counts_3hit)

## figure 1a: get number of patients who had at least three dz mentions before schizo in both case and control cohorts 
tab <- matrix(c(case_pre_pt_counts_3hit, ctrl_pre_pt_counts_3hit, 63133-case_pre_pt_counts_3hit, 63133-ctrl_pre_pt_counts_3hit), nrow=2)
fisher_result <- fisher.test(tab)
# columns: phecode, case_pre_counts, ctrl_pre_counts, OR, pval, ci_low, ci_high
line <- paste(phe, case_pre_pt_counts_3hit, ctrl_pre_pt_counts_3hit, fisher_result$estimate, fisher_result$p.value, fisher_result$conf.int[1], fisher_result$conf.int[2], sep = ',')
write(line,file="XXX/schizo/results/analysis_on_local_icds/pre_schizo_dz_counts_3hit.csv",append=TRUE)

## figure 4. columns: phecode, time-0, time-1, time-2, time-3, time-4, time-5, time-6, time-7
line_fg4_case <- paste(phe, paste(case_counts_3hit, collapse=','), sep=',')
line_fg4_ctrl <- paste(phe, paste(ctrl_counts_3hit, collapse=','), sep=',')
write(line_fg4_case, file="XXX/schizo/results/analysis_on_local_icds/pre_schizo_dz_counts_by_time_case_3hit.csv", append=TRUE)
write(line_fg4_ctrl, file="XXX/schizo/results/analysis_on_local_icds/pre_schizo_dz_counts_by_time_ctrl_3hit.csv", append=TRUE)

print(paste('phecode', phe, 'finished counting at', Sys.time(), 'Time used:', as.numeric(Sys.time()-start_time, units='mins'), 'mins'))
