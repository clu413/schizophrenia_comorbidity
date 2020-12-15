## Total number of unique members in database: 85972617 
SELECT COUNT(*)
FROM XXXXX.dbo.Members

## get unique memberIds of patients who have schizophrenia entries
SELECT 
  F.MemberId AS MemberID,
  F.DateServiceStarted, 
  F.Icd,
  F.IcdVersion
INTO 
  YYY.dbo.SchizoEntries
FROM
  XXXXX.dbo.FactIcd F 
  INNER JOIN PheWAS.dbo.PhewasCodeIcdCrosswalk I ON
  F.Icd=I.IcdCode
WHERE
  I.PheWASCode IN ('295', '295.1') AND Icd != 'V11.0' # V11.0 Icd code means personal history of schizophrenia - not reliable
GROUP BY
  F.MemberId, F.DateServiceStarted, F.Icd, F.IcdVersion
	
## Get all schizophrenia patients first recorded date of schizophrenia diagnosis: 406740
SELECT 
  MemberID, 
  COUNT(MemberID) as SchizoMentions, 
  MIN(DateServiceStarted) as SchizoOnset
INTO
  YYY.dbo.SchizoPatients
FROM 
  YYY.dbo.SchizoEntries
GROUP BY 
  MemberID
ORDER BY 
  SchizoMentions DESC

## Add age and gender information of schizophrenia patients: 398208 patients still left
SELECT 
  M.MemberId, 
  M.Gender, 
  M.BirthYear
INTO 
  YYY.dbo.SchizoPatientAgeGender
FROM 
  XXXXX.dbo.Members M
INNER JOIN 
  YYY.dbo.SchizoPatients P
ON 
  P.MemberId = M.MemberId
GROUP BY 
  M.MemberId, M.Gender, M.BirthYear

## Add zip codes of schizophrenia patients (Get the most recent zipcode before first report of schizo)
WITH cte as (
  SELECT 
    p.MemberId, 
    p.SchizoMentions, 
    p.SchizoOnset, 
    e.EffectiveDate, 
    e.Zipcode, 
    RN = row_number() over (partition by p.MemberId order by EffectiveDate DESC)
  FROM 
    XXXXX.dbo.Enrollment e
  INNER JOIN 
    YYY.dbo.SchizoPatients p
  ON 
    p.MemberId = e.MemberId
  AND 
    e.EffectiveDate <= p.SchizoOnset
  )
SELECT 
  MemberId, 
  Zipcode, 
  SchizoMentions, 
  SchizoOnset
INTO 
  YYY.dbo.SchizoPatientZip
FROM 
  cte 
WHERE 
  RN  = 1

## Combine information of schizophrenia patients so far
SELECT 
  z.MemberId, 
  a.BirthYear, 
  a.Gender, 
  z.Zipcode, 
  z.SchizoMentions, 
  z.SchizoOnset 
INTO 
  YYY.dbo.SchizoPatientAgeGenderZipMentionFirst
FROM 
  YYY.dbo.SchizoPatientZip z
INNER JOIN 
  YYY.dbo.SchizoPatientAgeGender a
ON 
  z.MemberId = a.MemberId

## Get all XXXXX Members enrollment periods 
SELECT 
  M.MemberId, 
  MIN(EffectiveDate) as XXXXXFirstDate, 
  MAX(EOMONTH(EffectiveDate)) as XXXXXLastDate 
INTO 
  YYY.dbo.AllMemberEnrollDates
FROM 
  XXXXX.dbo.Members M 
INNER JOIN 
  XXXXX.dbo.Enrollment E 
ON 
  M.MemberId=E.MemberId 
GROUP BY 
  M.MemberId

## Filter out schizophrenia patients who were not enrolled at least 12 months prior to schizo onset: 235979 
SELECT 
  m.*, 
  e.XXXXXFirstDate, 
  e.XXXXXLastDate
INTO 
  YYY.dbo.schizoFinalCohort
FROM 
  YYY.dbo.AllMemberEnrollDates e
INNER JOIN 
  YYY.dbo.SchizoPatientAgeGenderZipMentionFirst m
ON 
  m.MemberId = e.MemberId
WHERE 
  DATEADD(year, -1,  m.SchizoOnset) > e.XXXXXFirstDate 
  AND e.XXXXXLastDate > m.SchizoOnset

## Add icd counts between schizophrenia onset and 12 months prior: 235979 patients
SELECT 
  p.*,COUNT(Icd) as IcdCounts12mo
INTO 
  YYY.dbo.schizoFinalCohortwIcd
FROM 
  XXXXX.dbo.FactIcd f
INNER JOIN 
  YYY.dbo.schizoFinalCohort p
ON 
  f.MemberId = p.MemberId
WHERE 
  f.DateServiceStarted BETWEEN DATEADD(year, -1,  p.SchizoOnset) 
  AND p.SchizoOnset	
GROUP BY 
  p.MemberId, p.BirthYear, p.Gender, p.Zipcode, p.SchizoMentions, p.SchizoOnset, p.XXXXXFirstDate, p.XXXXXLastDate
