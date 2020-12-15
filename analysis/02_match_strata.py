## bin schizophrenia cases and controls based on birthyear, gender, and first three digits of zipcode (defined as strata)
## save the patients in stratas for potential future one-to-one matching

import sys
import csv
import re
import os
import io
import random
from datetime import datetime
import numpy as np
import pandas as pd
import bisect

DATA_DIR = 'YYYYY/schizo/data/retrieved_data/'
OUTPUT_DIR = 'YYYYY/schizo/data/pairs/'
caseInputFullInfo =  os.path.join(DATA_DIR,"schizoFinalCohortwIcd.txt")
ctrlInputFullInfo = os.path.join(DATA_DIR, "NonSchizoPatientsAgeGenderZipEnrollment.txt")

## map the case cohort into birth year, gender, zip groups
dzMemberMap = {}

with open(caseInputFullInfo,'r') as InputFile:
    next(InputFile)
    for line in InputFile:
        try:
            memberId, birthyear, gender, zipcode, schizomentions, schizoonset, XXXXXfirstdate, XXXXXlastdate, icdcounts12mo  = line.strip().split(',')
            memberId = memberId.strip('"')
            gender = gender.strip('"')
            birthyear = int(birthyear)
            zipcode = zipcode.strip('"')[:-2]
            schizoonset = schizoonset.strip('"')
            XXXXXfirstdate = XXXXXfirstdate.strip('"')
            XXXXXlastdate = XXXXXlastdate.strip('"')
        except:
            pass
        else:
            if birthyear >= 1900 and birthyear <=2017 and gender in ['M', 'F'] and len(zipcode) == 3:
                if not (birthyear, gender, zipcode) in dzMemberMap:
                    dzMemberMap[birthyear, gender, zipcode] = []
                dzMemberMap[birthyear, gender, zipcode].append([memberId, schizomentions, schizoonset, XXXXXfirstdate, XXXXXlastdate, icdcounts12mo])

print('finished mapping cases!')
# then map the control cohort into birth year, gender, zipcode groups
ctrlMemberMap={}
with open(ctrlInputFullInfo,'r') as InputFile:
    next(InputFile)
    for line in InputFile:
        try:
            memberId, birthyear, gender, zipcode, XXXXXfirstdate, XXXXXlastdate  = line.strip().split(',')
            memberId = memberId.strip('"')
            birthyear = int(birthyear[0:4])
            gender = gender.strip('"')
            zipcode = zipcode.strip('"')[:-2]
            XXXXXfirstdate = XXXXXfirstdate.strip('"')
            XXXXXlastdate = XXXXXlastdate.strip('"')
        except:
            pass
        else:
            if not (birthyear, gender, zipcode) in ctrlMemberMap:
                ctrlMemberMap[birthyear, gender, zipcode] = []
            ctrlMemberMap[birthyear, gender, zipcode].append([memberId,XXXXXfirstdate,XXXXXlastdate])

print('finished mapping ctrls!')

# some stats
tot_case = 0
tot_ctrl = 0
Unmatched_birthGenderZip = []

for birthyear, gender, zipcode in dzMemberMap:
    strata_name = str(birthyear)+'_'+gender+'_'+zipcode
    case_file = strata_name+'_case.txt'
    ctrl_file = strata_name+'_ctrl.txt'
    case_df = pd.DataFrame(data = dzMemberMap[birthyear, gender, zipcode])
    case_df.columns = ['case_id','SchizoMentions','Schizo_Onset','case_XXXXXFirstDate','case_XXXXXLastDate','IcdCounts12mo']
    case_df.to_csv(os.path.join(OUTPUT_DIR,case_file), sep=',', index=False)
    #if ctrlMemberMap.get((birthyear, gender, zipcode)) is not None:
    if (birthyear, gender, zipcode) in ctrlMemberMap:
        ctrl_df = pd.DataFrame(data = ctrlMemberMap[birthyear, gender, zipcode])
        ctrl_df.columns = ['ctrl_id','ctrl_XXXXXFirstDate','ctrl_XXXXXLastDate']
        ctrl_df.to_csv(os.path.join(OUTPUT_DIR,ctrl_file), sep=',',index=False)
        tot_ctrl += len(ctrlMemberMap[birthyear, gender, zipcode])
    else:
        Unmatched_birthGenderZip.append([birthyear, gender, zipcode])
    tot_case += len(dzMemberMap[birthyear, gender, zipcode])

print('tot_case: {}, tot_ctrl: {}, num_unmatched_strata: {}'.format(tot_case, tot_ctrl, len(Unmatched_birthGenderZip)))
print(Unmatched_birthGenderZip)
