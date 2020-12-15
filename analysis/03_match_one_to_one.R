## matches case patients with a control patient in the same strata
## make sure the control stayed in the dataset after the case's first schizophrenia record

library(dplyr)
library(SqlServerJtds)
library(SqlTools)
library(FactToCube)
library(stringr)

args <- commandArgs()
strata=args[6]
print(paste('-------------------------processing',strata,'-----------------------'))
birthyear <- str_split(strata, "_" )[[1]][1]
gender <- str_split(strata, "_" )[[1]][2]
zipcode <- str_split(strata, "_" )[[1]][3]

print(paste('start time is',Sys.time()))

# XXXXX database
my.username= XXXXXX
my.password= "YYYYYYY"

cn = connect.sql.server(
  database="ZZZZZZZ",
  domain="MED",
  server.address="ZZZZZZZ",
  user=my.username,
  password=my.password)

DAT_DIR <- 'XXX/schizo/data/pairs/'
case_file <- paste(DAT_DIR, strata, '_case.txt', sep ='')
ctrl_file <- paste(DAT_DIR, strata, '_ctrl.txt', sep ='')

case_info <- read.table(case_file, colClasses = 'character', header = TRUE, sep=',')
ctrl_info <- read.table(ctrl_file, colClasses = 'character', header = TRUE, sep=',')

df_match <- data.frame(matrix(nrow=length(dim(case_info)[2]),ncol=7))
colnames(df_match) <- c('case_id','schizo_onset','case_XXXXXFirstDate','case_XXXXXLastDate','ctrl_id','ctrl_XXXXXFirstDate','ctrl_XXXXXLastDate')

j <- 0
for (case_index in 1:nrow(case_info)){
    print(j)
    case_row <- case_info[case_index,]
    case <- case_row['case_id']
    schizo_onset <- as.Date(case_row['Schizo_Onset'][1,1])
    case_XXXXXfirstdate <- as.Date(case_row['case_XXXXXFirstDate'][1,1])
    case_XXXXXlastdate <- as.Date(case_row['case_XXXXXLastDate'][1,1])
    year_before_schizo <- schizo_onset - 365

    ## reorder the control list so that some earlier members in the list do not get preference
    ctrl_info_reorder <- ctrl_info[sample(nrow(ctrl_info)),]
    match_ctrl_id <- NULL

    ## keep searching until we find a control that could match with the case
    for (ctrl_index in 1:dim(ctrl_info_reorder)[1]){
        ctrl_row <- ctrl_info_reorder[ctrl_index,]
        ctrl <- ctrl_row['ctrl_id'][1,1]
        ctrl_XXXXXfirstdate <- as.Date(ctrl_row['ctrl_XXXXXFirstDate'][1,1])
        ctrl_XXXXXlastdate <- as.Date(ctrl_row['ctrl_XXXXXLastDate'][1,1])
        if (ctrl_XXXXXlastdate > schizo_onset & ctrl_XXXXXirstdate < year_before_schizo){
            print('enrolled in XXXXX 12 months prior to schizo onset and stayed after schizo onset')
            #assign as matching control
            match_ctrl_id <- ctrl
            print(paste('found a match',match_ctrl_id))
            break
          }
        }
    j <- j + 1
    #remove matched control from potential controls for the next case patient
    if (!is.null(match_ctrl_id)){
        print(paste('Removing control',ctrl,'from all potential controls because they are already matched to a case. Current time is',Sys.time()))
        ctrl_info <- filter(ctrl_info, ctrl_id != ctrl)
        df_match[j,] <- c(case, as.character(schizo_onset), as.character(case_XXXXXfirstdate), as.character(case_XXXXXlastdate), match_ctrl_id, as.character(ctrl_aetnafirstdate), as.character(ctrl_aetnalastdate))
    } else{
        print(paste('Cannot find a matching control for case',case, '. Current time is', Sys.time()))
        df_match[j,] <- c(case, as.character(schizo_onset), as.character(case_XXXXXfirstdate), as.character(case_XXXXXlastdate), NA, NA, NA)
    }
}
print(paste('Finished matching strata',strata,'at', Sys.time()))
write.csv(df_match, paste(DAT_DIR,strata,'_match.csv',sep=''))
