library(vroom)
library(dplyr)

# load all cases's member Ids
aff <- read.csv('XXX/schizo/data/retrieved-data/schizo_affective_members_3hit.csv', header=T, colClasses = c('x'='character'))$x
other <- read.csv('XXX/schizo/data/retrieved-data/schizo_others_members_3hit.csv', header=T, colClasses = c('x'='character'))$x

members <- vroom('XXX/schizo/data/retrieved_data/schizoFinalCohortwIcd.txt')
df_aff <- members %>%
    filter(MemberId %in% aff) %>%
    select(MemberId:IcdCounts12mo)
df_other <- members %>%
    filter(MemberId %in% other) %>%
    select(MemberId:IcdCounts12mo)

df_match <- data.frame(aff_id=character(), aff_schizo_onset=as.Date(character()), aff_last_date=as.Date(character()), other_id=character(), other_schizo_onset=as.Date(character()), other_last_date=as.Date(character()), stringsAsFactors=FALSE)
for (row in 1:nrow(df_aff)){
    aff_id <- df_aff[row, 'MemberId'][[1]]
    birthyear <- df_aff[row, 'BirthYear'][[1]]
    gender <- df_aff[row,'Gender'][[1]]
    zip <- substr(df_aff[row, 'Zipcode'],0,3)
    aff_schizo_onset <- df_aff[row,'SchizoOnset'][[1]]
    aff_last_date <- df_aff[row,'XXXXXLastDate'][[1]]
    pot_ctrl <- df_other %>%
        filter(BirthYear == birthyear, Gender==gender) # only matching on age and sex
    if(dim(pot_ctrl)[1]!=0){ # at least one potential match
        other_sample <- pot_ctrl[sample(nrow(pot_ctrl),1),] # randomly pick one 
        other_id <- other_sample$MemberId
        other_schizo_onset <- other_sample$SchizoOnset
        other_last_date <- other_sample$XXXXXLastDate
        df_match[nrow(df_match)+1,] <- c(aff_id, as.character(aff_schizo_onset), as.character(aff_last_date), other_id, as.character(other_schizo_onset), as.character(other_last_date))
        df_other <- df_other[!(df_other$MemberId==other_id),] # remove this control so it will not be matched to other cases
    }
}


write.csv(df_match, 'XXX/schizo/results/analysis_on_local_icds/schizoaffective_agegender_match_3hit.csv')
