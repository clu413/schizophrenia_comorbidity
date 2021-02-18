# time-to-dz analysis between schizoaffective cases and other schizophrenia patients
# use 3-hit criteria

start_time <- Sys.time()

library(vroom)
library(survival)
library(dplyr)

args <- commandArgs()
print(args)
phe=args[6]

# identify all the icd codes for this phecode
phecode_df <- read.csv('XXX/schizo/data/phecode_update.txt',colClasses = c('PheWASCode'='character'))
df <- phecode_df[phecode_df$PheWASCode == phe,]
#icd_list <- paste(shQuote(df$IcdCode, type="sh"), collapse=", ")
icd_list <- list(unique(df$IcdCode))[[1]]
remove(phecode_df, df)

# get all icds of schizo patients with 3-hit schizo
case_icds_3hit <- vroom('XXX/schizo/data/retrieved_data/case_ctrl_icds/case_icds_3hit.txt')
pairs <- vroom('XXXschizo/results/analysis_on_local_icds/schizoaffective_agegender_match_3hit_noHx_15yo_and_above.csv')
print(paste('Finished loading all the cases!', Sys.time()))

# subset only the entries in this phecode
icds_phe <- case_icds_3hit %>%
    filter(Icd %in% icd_list)
remove(case_icds_3hit)

### make sure neither case nor control in the pair had the condition 12 months before 
aff_existing <- pairs %>%
    merge(icds_phe, by.x='aff_id', by.y='case_id') %>%
    group_by(aff_id) %>%
    filter(!any(DateServiceStarted < aff_schizo_onset & DateServiceStarted >= (aff_schizo_onset - years(1)))) %>%
    ungroup() %>%
    distinct(aff_id) %>%
    unlist(use.names = FALSE)

# repeat in controls
other_existing <- pairs %>%
    merge(icds_phe, by.x='other_id', by.y='case_id') %>%
    group_by(other_id) %>%
    filter(!any(DateServiceStarted < other_schizo_onset & DateServiceStarted >= (other_schizo_onset - years(1)))) %>%
    ungroup() %>%
    distinct(other_id) %>%
    unlist(use.names = FALSE)

# remove any pair that had either case or contrl with existing condition
pairs <- pairs %>%
    filter(!aff_id %in% aff_existing) %>%
    filter(!other_id %in% other_existing)

# turn pairs table from wide to long, with column indicating aff(1) or other(0)
aff_dates <- pairs[2:4] %>%
    mutate(schizo_type =1)
colnames(aff_dates) <- c('case_id','schizo_onset','XXXXX_last_date','schizo_type')
other_dates <- pairs[5:7] %>%
    mutate(schizo_type =0)
colnames(other_dates) <- c('case_id','schizo_onset','XXXXX_last_date','schizo_type')
pairs_long <- rbind(aff_dates,other_dates)
# only get icd entries of all cohorts
icds_phe_cohort <- icds_phe %>%
    filter(case_id %in% pairs_long$case_id)
# filter the icd entries after schizo onset and make sure patient had at least three mentions
icds_phe_cohort_post <- merge(icds_phe_cohort, pairs_long, by='case_id', all.x=T) %>%
    filter( DateServiceStarted > schizo_onset) %>%
    group_by(case_id) %>%
    filter(n()>=3)

# for each patient, get their first date of dz, time to dz and time to last date
dz_df <- icds_phe_cohort_post %>%
    group_by(case_id) %>%
    arrange(DateServiceStarted) %>%
    slice(1L) %>%
    mutate(day_to_dz=DateServiceStarted-schizo_onset)

rest_days_to_last <- pairs_long %>%
    filter(!case_id %in% dz_df$case_id) %>%
    mutate(day_to_last=XXXXX_last_date-schizo_onset)

## ----------- time-to-dz -------------
df_life <- data.frame(matrix(nrow=0,ncol=3))
# status: 1 = has dz, 0 = censored
# schizo: 1 = aff, 0 = other 
colnames(df_life) <- c('day','status','schizo')

# go through the dz patients first, because they are fewer than rest
if (nrow(dz_df)!=0){
    for (row in 1:dim(dz_df)[1]){
        df_life[nrow(df_life)+1,] <- c(dz_df[row,]$day_to_dz, 1, dz_df[row,]$schizo_type)
    }
}

for (row in 1:dim(rest_days_to_last)[1]){
    df_life[nrow(df_life)+1,] <- c(rest_days_to_last[row,]$day_to_last, 0, rest_days_to_last[row,]$schizo_type)
}

sd <- survdiff(Surv(day, status) ~ schizo, data = df_life) # get output
cox_mod <- coxph(Surv(day, status) ~ schizo, data = df_life)
HR <- coef(summary(cox_mod))[2]
pval <- coef(summary(cox_mod))[5]

line <- paste(phe, sd$obs[1], sd$exp[1], sd$obs[2], sd$exp[2], HR, pval, sep = ',')
write(line,file="XXX/schizo/results/analysis_on_local_icds/time_to_dz_schizoaffective_3hit.csv",append=TRUE)
end_time <- Sys.time()
print(paste('Phewas code',phe,'started at',start_time,'and finished at', end_time,'. Time difference is', end_time-start_time))

