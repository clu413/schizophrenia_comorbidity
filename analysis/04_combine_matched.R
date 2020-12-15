## This script will combine all one-to-one matching files
## and filter out the cases that do not have a matching control 

#file_list <- list.files('XXX/schizo/data/pairs_new',pattern="*_match.csv$")
file_list <- list.files('XXX/schizo/data/pairs/',pattern="*_match.csv$")

#rm(dataset)
print(paste('There are',length(file_list),'files present'))

na_row_count <- 0
i <- 0
for (file in file_list){
  # create
  if (!exists("dataset")){
    temp_dataset <- read.table(paste('XXX/schizo/data/pairs/',file, sep=''), header=TRUE, sep=",",colClasses = c(rep('factor',6)))
    dataset <- temp_dataset[!is.na(temp_dataset$ctrl_id),]
    na_row_count <- na_row_count + sum(is.na(temp_dataset$ctrl_id))
    rm(temp_dataset)
  }
  # append 
  if (exists("dataset")){
    temp_dataset <- read.table(paste('XXX/schizo/data/pairs/',file,sep=''), header=TRUE, sep=",",colClasses = c(rep('factor',6)))
    subset <- temp_dataset[!is.na(temp_dataset$ctrl_id),]
    dataset <- rbind(dataset, subset)
    na_row_count <- na_row_count + sum(is.na(temp_dataset$ctrl_id))
    rm(temp_dataset,subset)
    i <- i + 1
    if(i %% 1000==0) {
        print(paste('processed',i,'files'))
    }
  }
}

write.csv(dataset, 'XXX/schizo/data/schizo_match_BirthyearGenderZip.csv')
print(paste(na_row_count,'cases had no controls, and',dim(dataset)[1],'cases did have controls. Unmatching rate is', round(na_row_count/(na_row_count+dim(dataset)[1]),2)))
