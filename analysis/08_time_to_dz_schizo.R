# time-to-dz analysis between schizophrenia cases and their non-schizophrenia controls
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

# get all icds of patients with 3-hit schizo and their matching controls
pairs <- vroom('XXX/schizo/data/schizo_match_BirthyearGenderZip_3hit_noHx_15yo_and_above.csv')
case_icds_3hit <- vroom('XXX/schizo/data/retrieved_data/case_ctrl_icds/case_icds_3hit.txt')
print(paste('Finished loading all the cases!', Sys.time()))

ctrl_icds_3hit <- vroom('XXX/schizo/data/retrieved_data/case_ctrl_icds/ctrl_icds_3hit.txt')
print(paste('Finished loading all the ctrls!', Sys.time()))

# subset only the entries in this phecode
case_icds_phe <- case_icds_3hit %>%
    filter(Icd %in% icd_list)
remove(case_icds_3hit)

ctrl_icds_phe <- ctrl_icds_3hit %>%
    filter(Icd %in% icd_list)
remove(ctrl_icds_3hit)

### make sure neither case nor control in the pair had the condition 12 months before 
case_existing <- pairs %>%
    merge(case_icds_phe) %>%
    group_by(case_id) %>%
    filter(!any(DateServiceStarted < schizo_onset & DateServiceStarted >= (schizo_onset - years(1)))) %>%
    ungroup() %>%
    distinct(case_id) %>%
    unlist(use.names = FALSE)

# repeat in controls
ctrl_existing <- pairs %>%
    merge(ctrl_icds_phe) %>%
    group_by(ctrl_id) %>%
    filter(!any(DateServiceStarted < schizo_onset & DateServiceStarted >= (schizo_onset - years(1)))) %>%
    ungroup() %>%
    distinct(ctrl_id) %>%
    unlist(use.names = FALSE)

# remove any pair that had either case or contrl with existing condition
pairs <- pairs %>%
    filter(!case_id %in% case_existing) %>%
    filter(!ctrl_id %in% ctrl_existing)

# turn pairs table from wide to long, with column indicating aff(schizo) or control(0)
case_dates <- pairs[,c("case_id","schizo_onset","case_XXXXXLastDate")] %>%
    mutate(schizo_status =1)
colnames(case_dates) <- c('pt_id','schizo_onset','XXXXX_last_date','schizo_status')
ctrl_dates <- pairs[,c("ctrl_id","schizo_onset","ctrl_XXXXXLastDate")] %>%
    mutate(schizo_status =0)
colnames(ctrl_dates) <- c('pt_id','schizo_onset','XXXXX_last_date','schizo_status')
pairs_long <- rbind(case_dates,ctrl_dates)

# only get icd entries of all 3-hit pts
# combine the icds tables for downstream filtering
colnames(case_icds_phe)[1] <- 'pt_id'
colnames(ctrl_icds_phe) <- colnames(case_icds_phe)
icds_phe <- rbind(case_icds_phe, ctrl_icds_phe)

# filter the icd entries after schizo onset and make sure patient had at least three mentions
icds_phe_post <- merge(icds_phe, pairs_long, by='pt_id', all.x=T) %>%
    filter(DateServiceStarted > schizo_onset) %>%
    group_by(pt_id) %>%
    filter(n()>=3)

# for each patient, get their first date of dz, time to dz and time to last date
dz_df <- icds_phe_post %>%
    group_by(pt_id) %>%
    arrange(DateServiceStarted) %>%
    slice(1L) %>%
    mutate(day_to_dz=DateServiceStarted-schizo_onset)

rest_days_to_last <- pairs_long %>%
    filter(!pt_id %in% dz_df$pt_id) %>%
    mutate(day_to_last=XXXXX_last_date-schizo_onset)

## ----------- time-to-dz -------------
# status: 1 = has dz, 0 = censored
# schizo: 1 = schizo, 0 = control

df_life <- rest_days_to_last[,c('day_to_last' ,'schizo_status')]
df_life$status <- 0 # does not have dz
df_life <- df_life[c(1,3,2)]
colnames(df_life) <- c('day','status','schizo')

# go through the dz patients
if (nrow(dz_df)!=0){
    dz_df <- dz_df[c('day_to_dz','schizo_status')]
    dz_df$status <- 1 # has dz
    dz_df <- dz_df[c(1,3,2)] # reorder the columns
    colnames(dz_df) <- c('day','status','schizo')
    df_life <- rbind(dz_df, df_life)
}

sd <- survdiff(Surv(day, status) ~ schizo, data = df_life) # get output
cox_mod <- coxph(Surv(day, status) ~ schizo, data = df_life)
HR <- coef(summary(cox_mod))[2]
pval <- coef(summary(cox_mod))[5]

line <- paste(phe, sd$obs[1], sd$exp[1], sd$obs[2], sd$exp[2], HR, pval, sep = ',')
write(line,file="XXX/schizo/results/analysis_on_local_icds/time_to_dz_schizo_3hit.csv",append=TRUE)
end_time <- Sys.time()

print(paste('Phewas code',phe,'started at',start_time,'and finished at', end_time,'. Time difference is', end_time-start_time))
