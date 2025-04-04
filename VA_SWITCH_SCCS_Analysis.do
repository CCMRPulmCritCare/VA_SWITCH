version 18
set more off
cap log close
clear all
set linesize 80

cd ""

local c_date = c(current_date)
local date = subinstr("`c_date'", " ", "", .)

log using "Logs\va_switch_sccs_`date'.log", replace

********************************************************************************
* Clinical Outcomes Following a National Inhaler Formulary Change
*	- SCCS Primary Analysis	
* Author: Sarah Seelye
*
* Date Created: 2024 Apr 19		
* Last updated: 2025 Apr 4
********************************************************************************

****************************************************************************************************
* Using Dataset with 0% Grace Period for Table 1 and Sensitivity Analyses in Supplemental Table S4 *
****************************************************************************************************

*------------------------------------------------------
* Merge in datasets with new variables for Table 1
*------------------------------------------------------

** Combat Veteran Status **
use "switch_military_service_20240423", clear

* create new combat veteran variable
tab combatvet
gen combatvet_numeric = combatvet=="Y"
tab combatvet_numeric combatvet, m
drop combatvet 
rename combatvet_numeric combatvet

* identify patients who ever served in combat
bysort patienticn (deployment_era): egen ever_combatvet = max(combatvet)

* drop variables no longer needed - only need ever_combatvet
drop deployment_era combatvet 

* keep one row per patient 
duplicates drop

* create temp file for subsequent merge
tempfile combat 
save `combat'

** Smoking Status ** 
use "switch_cohort_smoking_status_20240325", clear

* create temp file for subsequent merge 
tempfile smoking 
save `smoking'


** PFT **
use "switch_pft_data_with_qual_only_added_20240429", clear

* two rows of data contain two values for fev1s; recode using first value listed 
gen pft_date2 = date(pft_date, "YMD")	
format pft_date2 %td 
order pft_date2, after(pft_date)
rename pft_date pft_date_str
rename pft_date2 pft_date

order pft_date, after(patientsid)

* some pft_dates have missing data; correct these, as they're in clock format 
count if missing(pft_date) 
gen double pft_date3 = clock(pft_date_str, "YMDhms")
format pft_date3 %tc
gen pft_date4 = dofc(pft_date3)
format pft_date4 %td

replace pft_date = pft_date4 if missing(pft_date) 
drop pft_date_str pft_date3 pft_date4

* encode obstruction
tab obstruction, m
encode obstruction, gen(obstruction2)
tab obstruction obstruction2, m
drop obstruction
rename obstruction2 obstruction

tab obstruction, m //27 missing

* create new catetagorical variable for fev1_severity 
tab fev1_severity
gen fev1_severity2 = .
replace fev1_severity2 = 1 if fev1_severity=="Very Severe"
replace fev1_severity2 = 2 if fev1_severity=="Severe"
replace fev1_severity2 = 3 if fev1_severity=="Moderately Severe"
replace fev1_severity2 = 4 if fev1_severity=="Moderate"
replace fev1_severity2 = 5 if fev1_severity=="Mild"
replace fev1_severity2 = 6 if fev1_severity=="Normal"
lab def fev1_severity2 1 "Very Severe" 2 "Severe" 3 "Moderately Severe" ///
					   4 "Moderate" 5 "Mild" 6 "Normal"
lab val fev1_severity2 fev1_severity2

tab fev1_severity fev1_severity2
rename fev1_severity fev1_severity_str
rename fev1_severity2 fev1_severity
order fev1_severity, after(pft_date)

* confirm fev1 severity
	* Normal = >=80
	* Mild = 70-79
	* Moderate = 60-69
	* Moderately severe = 50-59
	* Severe = 35-49
	* Very severe = <35

bysort fev1_severity: sum fev1_perc_pred 

* create new fev1 severity cutpoints (<50, >=50)
gen fev1_lt50 = inrange(fev1_severity, 1, 2)
gen fev1_50plus = inrange(fev1_severity, 3, 6)

tab fev1_severity fev1_lt50
tab fev1_severity fev1_50plus

* obstruction 
tab obstruction, m

* indicator for whether patient has a pft 
gen pft_ind = 1

* keep only variables we need 
keep patienticn pft_date fev1_severity fev1_perc_pred fev1_fvc obstruction pft_ind

* keep the first pft for each patient between Jan 1 2018-Dec 31 2022
bysort patienticn (pft_date): gen pft_num = _n
keep if pft_num==1
drop pft_num

count 

* save for merge 
tempfile pft
save `pft'


** MERGE ** 

* open strict dataset 
use "switch_20241030", clear

* merge with combat file
merge m:1 patienticn using `combat' 
tab patrecnum if patrecnum==1 & _merge==1

* drop patients not in analytic dataset 
drop if _merge==2

drop _merge

* merge with smoking file 
merge m:1 patienticn using `smoking'

* drop patients not in analytic data set 
drop if _merge==2
drop _merge

* merge with pft 
merge m:1 patienticn using `pft'
drop if _merge==2 
drop _merge


*--------------------------------------
* Prep dataset for SCCS analysis
*--------------------------------------

* create exposure groups 
tab inhaler_type
tab inhaler_type, nol 

gen no_inhaler_exgrp = inhaler_type==0 
gen symbicort_exgrp = inhaler_type==1 
gen wixela_exgrp = inhaler_type==2 
gen other_exgrp = inhaler_type==3 

tab inhaler_type no_inhaler_exgrp, m
tab inhaler_type symbicort_exgrp, m
tab inhaler_type wixela_exgrp, m
tab inhaler_type other_exgrp, m 

* create loglength 
gen loglength = log(length)
tab length if loglength==.
recode loglength .=0 

count 
tab inhaler_type prednisone, ro
bysort inhaler_type: sum albuterol_equiv, de

sort patienticn start

* check number of prednisone scripts that patients receive 
bysort patienticn (patrecnum): egen tot_prednisone = total(prednisone)
tab tot_prednisone if patrecnum==1

* check number of intervals for each patient 
bysort patienticn (patrecnum): egen max_num_intervals = max(patrecnum)
tab max_num_intervals if patrecnum==1
drop if max_num_intervals==1
drop max_num_intervals

* recode pft to replace missing with 0 
tab pft_ind, m
recode pft_ind (.=0)

* creater never & former smoker category 
gen never_smoker = ever_smoker==0
gen former_smoker = ever_smoker==1 & current_smoker==0

* identify age at inhaler switch 
gsort patienticn -wixela_exgrp start
by patienticn: gen mkg_age_at_switch = age(dob, inhaler_releasedate_new)  if wixela_exgrp==1 & _n==1 & cohort5==1 //use age at switch for those who switch
by patienticn: replace mkg_age_at_switch = age(dob, td(01jul2021)) if cohort5==0 //use 7/01/2021 for those who don't switch
by patienticn: egen age_at_switch = max(mkg_age_at_switch)
drop mkg_age_at_switch
order age_at_switch, after(age)

* identify census region at inhaler switch
gsort patienticn -wixela_exgrp start

	* use region at switch date for those who switch 
	by patienticn: gen mkg_region_at_switch = address_census_region if wixela_exgrp==1 & _n==1 & cohort5==1 //use region at switch for those who switch

	* use region at 7/01/2021 for those who don't switch 
	gen preswitchdate_str = "01jul2021"
	gen preswitchdate = date(preswitchdate_str, "DMY")
	format preswitchdate %td
	gen preswitchdate_ind = (preswitchdate>=start & preswitchdate<=end) 
	
	by patienticn: replace mkg_region_at_switch = address_census_region if cohort5==0 & preswitchdate_ind==1
	
	* for patients who don't switch and enter the overall cohort after 7/01/2021, use the first region listed
	bysort patienticn (start): replace mkg_region_at_switch = address_census_region if preswitchdate<start & _n==1 & cohort5==0
	
	* create patient-level variable
	by patienticn: egen region_at_switch = max(mkg_region_at_switch)
	tab region_at_switch, m

	order region_at_switch, after(address_census_region)

	lab def region_at_switch 1 "Northeast" 2 "Midwest" 3 "South" 4 "West"
	lab val region_at_switch region_at_switch
	
	drop mkg_region_at_switch preswitchdate_str preswitchdate_ind

* identify date of switch to wixela. if patient doesn't have a switch date, 
* use 7/01/2021
gsort patienticn -wixela_exgrp start
by patienticn: gen patientswitchdate = inhaler_releasedate_new if wixela_exgrp==1 & _n==1 & cohort5==1
by patienticn: replace patientswitchdate = patientswitchdate[_n-1] if cohort5==1 & missing(patientswitchdate)
format patientswitchdate %td
replace patientswitchdate = preswitchdate if patientswitchdate==.

sort patienticn start

* identify 1-year pre-switch date
gen year_preswitch_date = patientswitchdate-365
format year_preswitch_date %td

* identify 1-year post-switch date 
gen year_postswitch_date = patientswitchdate+365
format year_postswitch_date %td

drop preswitchdate

* save tempfile  
tempfile intervals 
save `intervals'

*--------------------------------------------------------------
* count number of each outcome in the year prior switch
*--------------------------------------------------------------

* keep only dates of patient switch and year prior to patient switch to 
* merge with outcomes dataset 
keep if patrecnum==1
keep patienticn patientswitchdate year_preswitch_date

* merge patients' switch and preswitch dates with outcomes dataset 
merge 1:m patienticn using "switch_cohort_outcomes_20241010.dta"
drop if _merge==2
drop _merge 
sort patienticn datevalue

* create an indicator to identify records that fall within one year of the 
* preswitchdate; and keep only records that fall within the preswitch year  
gen year_preswitch_ind = inrange(datevalue, year_preswitch_date, patientswitchdate)
keep if year_preswitch_ind==1

* create new respiratory indicators for primary copd, asthma, pneumonia 
gen respiratory_primary_ed = inlist(1, asthma_primary_ed, copd_primary_ed, pneumonia_primary_ed)
gen respiratory_primary_hosp = inlist(1, asthma_primary_hosp, copd_primary_hosp, pneumonia_primary_hosp)

* keep only outcomes we need 
keep patienticn year_preswitch_date patientswitchdate datevalue 			///
	 year_preswitch_ind albuterol_inhaler_equiv_daily prednisone_daily 		///
	 any_cause_ed respiratory_primary_ed pneumonia_primary_ed				///
	 any_cause_hosp respiratory_primary_hosp pneumonia_primary_hosp

* count number of each outcome during patients' pre-switch year	 
bysort patienticn (datevalue): egen prednisone_yr_preswitch = sum(prednisone_daily) 
bysort patienticn (datevalue): egen albuterol_yr_preswitch = sum(albuterol_inhaler_equiv_daily) 
bysort patienticn (datevalue): egen edanycause_yr_preswitch = sum(any_cause_ed) 
bysort patienticn (datevalue): egen edresp_yr_preswitch = sum(respiratory_primary_ed) 
bysort patienticn (datevalue): egen edpneu_yr_preswitch = sum(pneumonia_primary_ed) 
bysort patienticn (datevalue): egen hospanycause_yr_preswitch = sum(any_cause_hosp) 
bysort patienticn (datevalue): egen hospresp_yr_preswitch = sum(respiratory_primary_hosp) 
bysort patienticn (datevalue): egen hosppneu_yr_preswitch = sum(pneumonia_primary_hosp) 

* keep only variables needed to merge back to intervals data set; 
* keep only one record per patient
keep patienticn prednisone_yr_preswitch albuterol_yr_preswitch 				///
	 edanycause_yr_preswitch edresp_yr_preswitch edpneu_yr_preswitch 		///
	 hospanycause_yr_preswitch hospresp_yr_preswitch hosppneu_yr_preswitch

duplicates drop 	 

* count number of patients who do not have any records during their 1-year 
* preswitch period
distinct patienticn
di 347486-341967
di 5519/347486 //1.6% of full cohort
	 
* save tempfile and merge back to intervals data set	 
tempfile yrpreswitch 
save `yrpreswitch'   

*--------------------------------------------------
* Merging data from year pre and post switch 
*--------------------------------------------------
use `intervals'
merge m:1 patienticn using `yrpreswitch' 
drop _merge 

*-------------------------------------------------------------------------------
* Create binary indicators for each of the outcomes in the year pre switch 
* and for the duration of the study
*-------------------------------------------------------------------------------

foreach var in prednisone edanycause edresp edpneu hospanycause   ///
			   hospresp hosppneu		{
	tab `var'_yr_preswitch cohort5 if patrecnum==1, m
} 

sum albuterol_yr_preswitch if patrecnum==1 & cohort5==0 
sum albuterol_yr_preswitch if patrecnum==1 & cohort5==1

* year pre-switch
foreach var in prednisone albuterol edanycause edresp edpneu hospanycause   ///
			   hospresp hosppneu		{
	gen `var'_yr_preswitch_ind = `var'_yr_preswitch>0 if !missing(`var'_yr_preswitch)
}

foreach var in prednisone albuterol edanycause edresp edpneu hospanycause   ///
			   hospresp hosppneu		{
	bysort `var'_yr_preswitch_ind: sum `var'_yr_preswitch
}			   

foreach var in prednisone albuterol edanycause edresp edpneu hospanycause   ///
			   hospresp hosppneu		{
	tab `var'_yr_preswitch_ind cohort5 if patrecnum==1
}			   


*----------------- 
* Table 1 
*-----------------

* Overall Cohort
tab male if patrecnum==1
sum age_at_switch if patrecnum==1, de
tab copd if patrecnum==1
tab asthma if patrecnum==1
tab current_smoker if patrecnum==1
tab former_smoker if patrecnum==1
tab never_smoker if patrecnum==1
tab ever_combatvet if patrecnum==1
tab region_at_switch if patrecnum==1

foreach var in albuterol prednisone edanycause edresp edpneu hospanycause   ///
			   hospresp hosppneu		{
	tab `var'_yr_preswitch_ind if patrecnum==1
}			   

* Excluded Patients from Cohort 5 (e.g. Overall Cohort = Excluded Patients + Cohort 5)
tab male cohort5 if patrecnum==1 , chi co
sum age_at_switch if patrecnum==1 & cohort5==0 , de
ttest age_at_switch if patrecnum==1, by(cohort5)
tab copd cohort5 if patrecnum==1 , chi co
tab asthma cohort5 if patrecnum==1 , chi co
tab current_smoker cohort5 if patrecnum==1 , chi co
tab former_smoker cohort5 if patrecnum==1 , chi co
tab never_smoker cohort5 if patrecnum==1 , chi co
tab ever_combatvet cohort5 if patrecnum==1 , chi co
tab region_at_switch cohort5 if patrecnum==1 , chi co

foreach var in albuterol prednisone edanycause edresp edpneu hospanycause   ///
			   hospresp hosppneu		{
	tab `var'_yr_preswitch_ind cohort5 if patrecnum==1 , chi co
}			   

sum albuterol_yr_preswitch if patrecnum==1 & cohort5==0, de
sum albuterol_yr_preswitch if patrecnum==1 & cohort5==1, de
ttest albuterol_yr_preswitch if patrecnum==1, by(cohort5)

* Keep cohort 5 for SCCS analysis 
tab cohort5
keep if cohort5==1
distinct patienticn

* Total in cohort 5
tab cohort5 if patrecnum==1

* Cohort 5 
tab male if patrecnum==1
sum age_at_switch if patrecnum==1, de
tab copd if patrecnum==1
tab asthma if patrecnum==1
tab current_smoker if patrecnum==1
tab former_smoker if patrecnum==1
tab never_smoker if patrecnum==1
tab ever_combatvet if patrecnum==1
tab region_at_switch if patrecnum==1

tab obstruction if patrecnum==1, m
tab fev1_severity if patrecnum==1, m


foreach var in albuterol prednisone edanycause edresp edpneu hospanycause   ///
			   hospresp hosppneu		{
	tab `var'_yr_preswitch_ind if patrecnum==1
}	


*-------------------------------------------------------------------------------
* 0% Grace Period Dataset - no exclusions; Cohort 5; Supplemental Table S4
*-------------------------------------------------------------------------------

* Albuterol inhaler equivalents	
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr

* Prednisone		
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr			

* Hospitalizations - All Cause
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* Hospitalizations - Respiratory 
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* Hospitalizations - Pneumonia 
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* ED - All Cause
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* ED - Respiratory 
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* ED - Pneumonia 
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr				
				
*-------------------------------------------------------------------------------
* Strict dataset - subset of patients with limited 'no controller periods' 
* (<3 months total)		
*-------------------------------------------------------------------------------

* create indicator for interval gap 
gen inhaler_gap = inhaler_type==0

* calculate number of days during all inhaler gaps 
gen gap_days = datediff(start, end, "day") if inhaler_gap==1
bysort patienticn (patrecnum): egen gap_days_tot = total(gap_days)

* indentify patients with more than 3 months of inhaler gaps 
gen inhaler_gap_lt3m = gap_days_tot<=91

tab inhaler_gap_lt3m 
tab inhaler_gap_lt3m if patrecnum==1 

* identify the last record for each patient 
gsort patienticn -patrecnum 
by patienticn: gen lastrec = _n 
replace lastrec = 0 if lastrec>1
sort patienticn patrecnum 

	* look at the number of gap days for the last record of a patient w/ gaps 
	* to confirm that everyone has 91 or fewer gap days
	sum gap_days if lastrec==1 & inhaler_gap==1
	
* inspect relationship between prednisone and inhaler type among people with
* <3m gap 
tab inhaler_type prednisone if inhaler_gap_lt3m==1, ro

bysort patienticn (patrecnum): egen ever_prednisone = max(prednisone)
tab inhaler_type prednisone if inhaler_gap_lt3m==1 & ever_prednisone>0, ro

tab ever_prednisone if patrecnum==1 & inhaler_gap_lt3m==1 
tab ever_prednisone if inhaler_gap_lt3m==1 

* SCCS analysis
preserve 

	keep if inhaler_gap_lt3m==1
	count 
	tab patrecnum if patrecnum==1 
	
		*---------------
		* Cohort 5
		*---------------
		
			* albuterol inhaler equivalents 
			xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp   ///
				c.age ib14.quarter##i.address_census_division , ///
				fe i(patienticn) offset(loglength) irr
			
			* Prednisone		
			xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp   ///
				c.age ib14.quarter##i.address_census_division , ///
				fe i(patienticn) offset(loglength) irr		

			* Hospitalizations - All Cause
			xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
					c.age ib14.quarter##i.address_census_division , ///
					fe i(patienticn) offset(loglength) irr		

			* Hospitalizations - Respiratory 
			xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
					c.age ib14.quarter##i.address_census_division , ///
					fe i(patienticn) offset(loglength) irr		

			* Hospitalizations - Pneumonia 
			xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
					c.age ib14.quarter##i.address_census_division , ///
					fe i(patienticn) offset(loglength) irr		

			* ED - All Cause
			xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
					c.age ib14.quarter##i.address_census_division , ///
					fe i(patienticn) offset(loglength) irr		

			* ED - Respiratory 
			xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
					c.age ib14.quarter##i.address_census_division , ///
					fe i(patienticn) offset(loglength) irr		

			* ED - Pneumonia 
			xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
					c.age ib14.quarter##i.address_census_division , ///
					fe i(patienticn) offset(loglength) irr				
						
restore	

						
********************************************************************************
* Using 33% Grace Period Dataset for Primary Analysis and Supplemental Table S4*
********************************************************************************

** Combat Veteran Status **
use "switch_military_service_20240423", clear

* create new combat vet variable
tab combatvet
gen combatvet_numeric = combatvet=="Y"
tab combatvet_numeric combatvet, m
drop combatvet 
rename combatvet_numeric combatvet

* identify patients who ever served in combat
bysort patienticn (deployment_era): egen ever_combatvet = max(combatvet)

* drop variables no longer needed - only need ever_combatvet
drop deployment_era combatvet 

* keep one row per patient 
duplicates drop

* create temp file to merge with strict dataset
tempfile combat 
save `combat'


** PFT **
use "switch_pft_data_with_qual_only_added_20240429", clear

* two rows of data contain two values for fev1s; recode using first value listed 
gen pft_date2 = date(pft_date, "YMD")	
format pft_date2 %td 
order pft_date2, after(pft_date)
rename pft_date pft_date_str
rename pft_date2 pft_date

order pft_date, after(patientsid)

* some pft_dates have missing data; correct these, as they're in clock format 
count if missing(pft_date) 
gen double pft_date3 = clock(pft_date_str, "YMDhms")
format pft_date3 %tc
gen pft_date4 = dofc(pft_date3)
format pft_date4 %td

replace pft_date = pft_date4 if missing(pft_date) 
drop pft_date_str pft_date3 pft_date4

* encode obstruction
tab obstruction, m
encode obstruction, gen(obstruction2)
tab obstruction obstruction2, m
drop obstruction
rename obstruction2 obstruction

tab obstruction, m 

* create new catetagorical variable for fev1_severity 
tab fev1_severity
gen fev1_severity2 = .
replace fev1_severity2 = 1 if fev1_severity=="Very Severe"
replace fev1_severity2 = 2 if fev1_severity=="Severe"
replace fev1_severity2 = 3 if fev1_severity=="Moderately Severe"
replace fev1_severity2 = 4 if fev1_severity=="Moderate"
replace fev1_severity2 = 5 if fev1_severity=="Mild"
replace fev1_severity2 = 6 if fev1_severity=="Normal"
lab def fev1_severity2 1 "Very Severe" 2 "Severe" 3 "Moderately Severe" ///
					   4 "Moderate" 5 "Mild" 6 "Normal"
lab val fev1_severity2 fev1_severity2

tab fev1_severity fev1_severity2
rename fev1_severity fev1_severity_str
rename fev1_severity2 fev1_severity
order fev1_severity, after(pft_date)

* confirm fev1 severity
	* Normal = >=80
	* Mild = 70-79
	* Moderate = 60-69
	* Moderately severe = 50-59
	* Severe = 35-49
	* Very severe = <35

bysort fev1_severity: sum fev1_perc_pred 

* obstruction 
tab obstruction, m

* indicator for whether patient has a pft 
gen pft_ind = 1

* keep only variables we need 
keep patienticn pft_date fev1_severity fev1_perc_pred fev1_fvc obstruction pft_ind

* keep the first pft for each patient between Jan 1 2018-Dec 31 2022
bysort patienticn (pft_date): gen pft_num = _n
keep if pft_num==1
drop pft_num

* save for merge 
tempfile pft
save `pft'

*** Next, open 33% leniency dataset ***
	
* open dataset 
use "switch_alt33_20241010", clear

* merge in smoking 
merge m:1 patienticn using "switch_cohort_smoking_status_20240325"
drop if _merge==2 
drop _merge

* merge in pft 
merge m:1 patienticn using `pft'
drop if _merge==2 
drop _merge

* merge with combat file
merge m:1 patienticn using `combat' 

* recode 4 patients with missing combat status as not a combat veteran
tab patrecnum if patrecnum==1 & _merge==1
replace ever_combatvet=0 if _merge==1
drop _merge

* create exposure groups 
tab inhaler_type
tab inhaler_type, nol 

gen no_inhaler_exgrp = inhaler_type==0 
gen symbicort_exgrp = inhaler_type==1 
gen wixela_exgrp = inhaler_type==2 
gen other_exgrp = inhaler_type==3 

tab inhaler_type no_inhaler_exgrp, m
tab inhaler_type symbicort_exgrp, m
tab inhaler_type wixela_exgrp, m
tab inhaler_type other_exgrp, m 

* create new fev1 severity cutpoints (<50, >=50)
gen fev1_lt50 = inrange(fev1_severity, 1, 2)
gen fev1_50plus = inrange(fev1_severity, 3, 6)

* create loglength 
gen loglength = log(length)
tab length if loglength==.
recode loglength .=0 

sort patienticn start

* check number of intervals for each patient 
bysort patienticn (patrecnum): egen max_num_intervals = max(patrecnum)
tab max_num_intervals if patrecnum==1
drop if max_num_intervals==1
drop max_num_intervals

* creater never & former smoker category 
gen never_smoker = ever_smoker==0
gen former_smoker = ever_smoker==1 & current_smoker==0

* identify preswitch date 
gen preswitchdate_str = "01jul2021"
gen preswitchdate = date(preswitchdate_str, "DMY")
format preswitchdate %td

* identify date of switch to wixela. if patient doesn't have a switch date, 
* use 7/01/2021
gsort patienticn -wixela_exgrp start
by patienticn: gen patientswitchdate = inhaler_releasedate_new if wixela_exgrp==1 & _n==1 & cohort5==1
by patienticn: replace patientswitchdate = patientswitchdate[_n-1] if cohort5==1 & missing(patientswitchdate)
format patientswitchdate %td
replace patientswitchdate = preswitchdate if patientswitchdate==.

sort patienticn start

* identify 1-year pre-switch date
gen year_preswitch_date = patientswitchdate-365
format year_preswitch_date %td

* save tempfile  
tempfile intervals 
save `intervals'

*--------------------------------------------------------------
* count number of each outcome in the year prior switch
*--------------------------------------------------------------

* keep only dates of patient switch and year prior to patient switch to 
* merge with outcomes dataset 
keep if patrecnum==1
keep patienticn patientswitchdate year_preswitch_date

* merge patients' switch and preswitch dates with outcomes dataset 
merge 1:m patienticn using "switch_cohort_outcomes_20241010.dta"
drop if _merge==2
drop _merge 
sort patienticn datevalue

* create an indicator to identify records that fall within one year of the 
* preswitchdate; and keep only records that fall within the preswitch year  
gen year_preswitch_ind = inrange(datevalue, year_preswitch_date, patientswitchdate)
keep if year_preswitch_ind==1

* create new respiratory indicators for primary copd, asthma, pneumonia 
gen respiratory_primary_ed = inlist(1, asthma_primary_ed, copd_primary_ed, pneumonia_primary_ed)
gen respiratory_primary_hosp = inlist(1, asthma_primary_hosp, copd_primary_hosp, pneumonia_primary_hosp)

* keep only outcomes we need 
keep patienticn year_preswitch_date patientswitchdate datevalue 			///
	 year_preswitch_ind albuterol_inhaler_equiv_daily prednisone_daily 		///
	 any_cause_ed respiratory_primary_ed pneumonia_primary_ed				///
	 any_cause_hosp respiratory_primary_hosp pneumonia_primary_hosp

* count number of each outcome during patients' pre-switch year	 
bysort patienticn (datevalue): egen prednisone_yr_preswitch = sum(prednisone_daily) 
bysort patienticn (datevalue): egen albuterol_yr_preswitch = sum(albuterol_inhaler_equiv_daily) 
bysort patienticn (datevalue): egen edanycause_yr_preswitch = sum(any_cause_ed) 
bysort patienticn (datevalue): egen edresp_yr_preswitch = sum(respiratory_primary_ed) 
bysort patienticn (datevalue): egen edpneu_yr_preswitch = sum(pneumonia_primary_ed) 
bysort patienticn (datevalue): egen hospanycause_yr_preswitch = sum(any_cause_hosp) 
bysort patienticn (datevalue): egen hospresp_yr_preswitch = sum(respiratory_primary_hosp) 
bysort patienticn (datevalue): egen hosppneu_yr_preswitch = sum(pneumonia_primary_hosp) 

* keep only variables needed to merge back to intervals data set; 
* keep only one record per patient
keep patienticn prednisone_yr_preswitch albuterol_yr_preswitch 				///
	 edanycause_yr_preswitch edresp_yr_preswitch edpneu_yr_preswitch 		///
	 hospanycause_yr_preswitch hospresp_yr_preswitch hosppneu_yr_preswitch

duplicates drop 	 

* count number of patients who do not have any records during their 1-year 
* preswitch period
distinct patienticn
di 347486-341991
di 5495/347486 //1.6% of full cohort
	 
* save tempfile and merge back to intervals data set	 
tempfile yrpreswitch 
save `yrpreswitch'

*--------------------------------------------------
* Merging data from year pre and post switch 
*--------------------------------------------------
use `intervals'
merge m:1 patienticn using `yrpreswitch' 
drop _merge 

*-------------------------------------------------------------------------------
* Create binary indicators for each of the outcomes in the year pre switch 
* and for the duration of the study
*-------------------------------------------------------------------------------

foreach var in prednisone edanycause edresp edpneu hospanycause   ///
			   hospresp hosppneu		{
	tab `var'_yr_preswitch cohort5 if patrecnum==1, m
} 

sum albuterol_yr_preswitch if patrecnum==1 & cohort5==0 
sum albuterol_yr_preswitch if patrecnum==1 & cohort5==1

* year pre-switch
foreach var in prednisone albuterol edanycause edresp edpneu hospanycause   ///
			   hospresp hosppneu		{
	gen `var'_yr_preswitch_ind = `var'_yr_preswitch>0 if !missing(`var'_yr_preswitch)
}

foreach var in prednisone albuterol edanycause edresp edpneu hospanycause   ///
			   hospresp hosppneu		{
	bysort `var'_yr_preswitch_ind: sum `var'_yr_preswitch
}			   

foreach var in prednisone albuterol edanycause edresp edpneu hospanycause   ///
			   hospresp hosppneu		{
	tab `var'_yr_preswitch_ind cohort5 if patrecnum==1, co
}	
		   
sum albuterol_yr_preswitch if patrecnum==1 & cohort5==0 , de
sum albuterol_yr_preswitch if patrecnum==1 & cohort5==1 , de

* create categories of albuterol use in the year pre-switch 
gen albuterol_yr_preswitch_rnd = round(albuterol_yr_preswitch)
gen albuterol_yr_preswitch_0or1 = albuterol_yr_preswitch_rnd<=1 if !missing(albuterol_yr_preswitch_rnd)
gen albuterol_yr_preswitch_2 = albuterol_yr_preswitch_rnd==2 if !missing(albuterol_yr_preswitch_rnd)
gen albuterol_yr_preswitch_3plus = albuterol_yr_preswitch_rnd>=3 if !missing(albuterol_yr_preswitch_rnd)

* Analysis will be on cohort 5 
tab cohort5
keep if cohort5==1
distinct patienticn

*--------------
* Cohort 5
*--------------

* Albuterol inhaler equivalents	
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr

* Prednisone		
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* Hospitalizations - All Cause
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* Hospitalizations - Respiratory 
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* Hospitalizations - Pneumonia 
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* ED - All Cause
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* ED - Respiratory 
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* ED - Pneumonia 
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr					
	
*----------------------------------------------------------		
* Subgroup analses for albuterol inhaler equivalents
*----------------------------------------------------------		

* male
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==1, ///
		fe i(patienticn) offset(loglength) irr		

* female		
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==0, ///
		fe i(patienticn) offset(loglength) irr			

* asthma 
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if asthma==1, ///
		fe i(patienticn) offset(loglength) irr	

* copd 
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if copd==1, ///
		fe i(patienticn) offset(loglength) irr	

* current smoker
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if current_smoker==1, ///
		fe i(patienticn) offset(loglength) irr	
		
* former smoker
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if former_smoker==1, ///
		fe i(patienticn) offset(loglength) irr			

* never smoker 
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if never_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* combat veteran
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if ever_combatvet==1, ///
		fe i(patienticn) offset(loglength) irr	

* albuterol fills in year pre switch 
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_0or1==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_2==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_3plus==1, ///
		fe i(patienticn) offset(loglength) irr			

* prednisone fills in year pre switch
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if prednisone_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr		

* any cause ED in year preswitch 
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if edanycause_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr	

* any cause hospitalizations in year preswitch
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if hospanycause_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr	

*----------------------------------------------------------		
* Subgroup analses for prednisone inhaler equivalents
*----------------------------------------------------------		

* male
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==1, ///
		fe i(patienticn) offset(loglength) irr		

* female		
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==0, ///
		fe i(patienticn) offset(loglength) irr			

* asthma 
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if asthma==1, ///
		fe i(patienticn) offset(loglength) irr	

* copd 
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if copd==1, ///
		fe i(patienticn) offset(loglength) irr	

* current smoker
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if current_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* never smoker 
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if never_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* former smoker
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if former_smoker==1, ///
		fe i(patienticn) offset(loglength) irr				

* combat veteran
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if ever_combatvet==1, ///
		fe i(patienticn) offset(loglength) irr		

* albuterol fills in year pre switch 
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_0or1==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_2==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_3plus==1, ///
		fe i(patienticn) offset(loglength) irr			

* prednisone fills in year pre switch
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if prednisone_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr		

* any cause ED in year preswitch 
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if edanycause_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr	

* any cause hospitalizations in year preswitch
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if hospanycause_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr	
		
*----------------------------------------------------------		
* Subgroup analyses for any cause ED visits
*----------------------------------------------------------		

* male
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==1, ///
		fe i(patienticn) offset(loglength) irr		

* female		
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==0, ///
		fe i(patienticn) offset(loglength) irr			

* asthma 
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if asthma==1, ///
		fe i(patienticn) offset(loglength) irr	

* copd 
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if copd==1, ///
		fe i(patienticn) offset(loglength) irr	

* current smoker
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if current_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* never smoker 
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if never_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* former smoker
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if former_smoker==1, ///
		fe i(patienticn) offset(loglength) irr				

* combat veteran
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if ever_combatvet==1, ///
		fe i(patienticn) offset(loglength) irr			

* albuterol fills in year pre switch 
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_0or1==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_2==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_3plus==1, ///
		fe i(patienticn) offset(loglength) irr			

* prednisone fills in year pre switch
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if prednisone_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr		

* any cause hospitalizations visits in year preswitch
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if hospanycause_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr		
	
	
*----------------------------------------------------------		
* Subgroup analyses for respiratory ED visits
*----------------------------------------------------------		

* male
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==1, ///
		fe i(patienticn) offset(loglength) irr		

* female		
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==0, ///
		fe i(patienticn) offset(loglength) irr			

* asthma 
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if asthma==1, ///
		fe i(patienticn) offset(loglength) irr	

* copd 
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if copd==1, ///
		fe i(patienticn) offset(loglength) irr	

* current smoker
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if current_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* never smoker 
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if never_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* former smoker
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if former_smoker==1, ///
		fe i(patienticn) offset(loglength) irr				

* combat veteran
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if ever_combatvet==1, ///
		fe i(patienticn) offset(loglength) irr			
	
* albuterol fills in year pre switch 
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_0or1==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_2==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_3plus==1, ///
		fe i(patienticn) offset(loglength) irr			

* prednisone fills in year pre switch
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if prednisone_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr		

* any cause hospitalizations in year preswitch
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if hospanycause_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr	
	
*----------------------------------------------------------		
* Subgroup analyses for pneumonia ED visits
*----------------------------------------------------------		

* male
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==1, ///
		fe i(patienticn) offset(loglength) irr		

* female		
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==0, ///
		fe i(patienticn) offset(loglength) irr			

* asthma 
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if asthma==1, ///
		fe i(patienticn) offset(loglength) irr	

* copd 
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if copd==1, ///
		fe i(patienticn) offset(loglength) irr	

* current smoker
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if current_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* never smoker 
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if never_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* former smoker
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if former_smoker==1, ///
		fe i(patienticn) offset(loglength) irr				

* combat veteran
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if ever_combatvet==1, ///
		fe i(patienticn) offset(loglength) irr			

* albuterol fills in year pre switch 
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_0or1==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_2==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_3plus==1, ///
		fe i(patienticn) offset(loglength) irr			

* prednisone fills in year pre switch
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if prednisone_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr		

* any cause hospitalizations in year preswitch
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if hospanycause_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr	
	
*----------------------------------------------------------		
* Subgroup analyses for any cause hospitalizations
*----------------------------------------------------------		

* male
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==1, ///
		fe i(patienticn) offset(loglength) irr		

* female		
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==0, ///
		fe i(patienticn) offset(loglength) irr			

* asthma 
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if asthma==1, ///
		fe i(patienticn) offset(loglength) irr	

* copd 
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if copd==1, ///
		fe i(patienticn) offset(loglength) irr	

* current smoker
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if current_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* never smoker 
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if never_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* former smoker
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if former_smoker==1, ///
		fe i(patienticn) offset(loglength) irr				

* combat veteran
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if ever_combatvet==1, ///
		fe i(patienticn) offset(loglength) irr		

* albuterol fills in year pre switch 
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_0or1==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_2==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_3plus==1, ///
		fe i(patienticn) offset(loglength) irr			

* prednisone fills in year pre switch
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if prednisone_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr		

* any cause ED in year preswitch 
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if edanycause_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr	
		
		
*----------------------------------------------------------		
* Subgroup analyses for respiratory hospitalizations
*----------------------------------------------------------		

* male
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==1, ///
		fe i(patienticn) offset(loglength) irr		

* female		
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==0, ///
		fe i(patienticn) offset(loglength) irr			

* asthma 
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if asthma==1, ///
		fe i(patienticn) offset(loglength) irr	

* copd 
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if copd==1, ///
		fe i(patienticn) offset(loglength) irr	

* current smoker
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if current_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* never smoker 
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if never_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* former smoker
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if former_smoker==1, ///
		fe i(patienticn) offset(loglength) irr				

* combat veteran
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if ever_combatvet==1, ///
		fe i(patienticn) offset(loglength) irr		

* albuterol fills in year pre switch 
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_0or1==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_2==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_3plus==1, ///
		fe i(patienticn) offset(loglength) irr			

* prednisone fills in year pre switch
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if prednisone_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr		

* any cause ED in year preswitch 
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if edanycause_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr	

*----------------------------------------------------------		
* Subgroup analyses for pneumonia hospitalizations
*----------------------------------------------------------		

* male
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==1, ///
		fe i(patienticn) offset(loglength) irr		

* female		
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if male==0, ///
		fe i(patienticn) offset(loglength) irr			

* asthma 
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if asthma==1, ///
		fe i(patienticn) offset(loglength) irr	

* copd 
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if copd==1, ///
		fe i(patienticn) offset(loglength) irr	

* current smoker
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if current_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* never smoker 
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if never_smoker==1, ///
		fe i(patienticn) offset(loglength) irr		

* former smoker
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if former_smoker==1, ///
		fe i(patienticn) offset(loglength) irr				
		
* combat veteran
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if ever_combatvet==1, ///
		fe i(patienticn) offset(loglength) irr			

* albuterol fills in year pre switch 
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_0or1==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_2==1, ///
		fe i(patienticn) offset(loglength) irr	
		
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if albuterol_yr_preswitch_3plus==1, ///
		fe i(patienticn) offset(loglength) irr			

* prednisone fills in year pre switch
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if prednisone_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr		

* any cause ED in year preswitch 
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if edanycause_yr_preswitch_ind==1, ///
		fe i(patienticn) offset(loglength) irr	 
				
		
*---------------------------------		
* FEV1 Subgroups by Outcome		
*---------------------------------

* fev1 >=50
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_50plus==1, ///
		fe i(patienticn) offset(loglength) irr 	
		
* fev1<50
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_lt50==1, ///
		fe i(patienticn) offset(loglength) irr		
		
* fev1 >=50
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_50plus==1, ///
		fe i(patienticn) offset(loglength) irr		

* fev1<50
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_lt50==1, ///
		fe i(patienticn) offset(loglength) irr		
		
* fev1 >=50
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_50plus==1, ///
		fe i(patienticn) offset(loglength) irr		

* fev1<50
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_lt50==1, ///
		fe i(patienticn) offset(loglength) irr		
				
* fev1 >=50
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_50plus==1, ///
		fe i(patienticn) offset(loglength) irr		

* fev1<50
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_lt50==1, ///
		fe i(patienticn) offset(loglength) irr		
				
* fev1 >=50
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_50plus==1, ///
		fe i(patienticn) offset(loglength) irr		

* fev1<50
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_lt50==1, ///
		fe i(patienticn) offset(loglength) irr		
				
* fev1 >=50
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_50plus==1, ///
		fe i(patienticn) offset(loglength) irr		

* fev1<50
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_lt50==1, ///
		fe i(patienticn) offset(loglength) irr		
				
* fev1 >=50
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_50plus==1, ///
		fe i(patienticn) offset(loglength) irr		

* fev1<50
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_lt50==1, ///
		fe i(patienticn) offset(loglength) irr		
				
* fev1 >=50
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_50plus==1, ///
		fe i(patienticn) offset(loglength) irr		

* fev1<50
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if fev1_lt50==1, ///
		fe i(patienticn) offset(loglength) irr		


*-------------------------------------------------------------------------------
* Create a subset of patients with limited 'no controller periods' 
* (<3 months total)	for Sensitivity Analyses in Supplemental Table S4	
*-------------------------------------------------------------------------------

* create indicator for interval gap 
gen inhaler_gap = inhaler_type==0

* calculate number of days during all inhaler gaps 
gen gap_days = datediff(start, end, "day") if inhaler_gap==1
bysort patienticn (patrecnum): egen gap_days_tot = total(gap_days)

* indentify patients with more than 3 months of inhaler gaps 
gen inhaler_gap_lt3m = gap_days_tot<=91

tab inhaler_gap_lt3m 
tab inhaler_gap_lt3m if patrecnum==1 

* identify the last record for each patient 
gsort patienticn -patrecnum 
by patienticn: gen lastrec = _n 
replace lastrec = 0 if lastrec>1
sort patienticn patrecnum 

	* look at the number of gap days for the last record of a patient w/ gaps 
	* to confirm that everyone has 91 or fewer gap days
	sum gap_days if lastrec==1 & inhaler_gap==1
	
* inspect relationship between prednisone and inhaler type among people with
* <3m gap 
tab inhaler_type prednisone if inhaler_gap_lt3m==1, ro

bysort patienticn (patrecnum): egen ever_prednisone = max(prednisone)
tab inhaler_type prednisone if inhaler_gap_lt3m==1 & ever_prednisone>0, ro

tab ever_prednisone if patrecnum==1 & inhaler_gap_lt3m==1 
tab ever_prednisone if inhaler_gap_lt3m==1 

* SCCS analysis
preserve 

	keep if inhaler_gap_lt3m==1
	count 
	tab patrecnum if patrecnum==1 
	
	
	*---------------
	* Cohort 5
	*---------------
	
	* albuterol inhaler equivalents 
	xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp   ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr
	
	* Prednisone		
	xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp   ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		
	
	* Hospitalizations - All Cause
	xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr		

	* Hospitalizations - Respiratory 
	xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr		

	* Hospitalizations - Pneumonia 
	xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr		

	* ED - All Cause
	xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr		

	* ED - Respiratory 
	xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr		

	* ED - Pneumonia 
	xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr				
			
	
restore	


*------------------------------------------------------------
* Total exposure and outcomes for each inhaler type
*------------------------------------------------------------

preserve 
	keep patrecnum patienticn start end inhaler_type albuterol_equiv  ///
		 prednisone any_cause_ed_total any_cause_hosp_total
	gen days_elapsed = end - start
	collapse (sum) days_elapsed (sum) albuterol_equiv (sum) prednisone ///
			 (sum) any_cause_ed_total (sum) any_cause_hosp_total	  ///
			 , by(inhaler_type)
	format days_elapsed %12.0g
	gen person_years = days_elapsed/365.25
	gen inhaler_type2 = inhaler_type 
	recode inhaler_type2 0=3
	label define inhaler_type2 1 "Symbicort" 2 "Wixela" 3 "None or Other"
	lab val inhaler_type2 inhaler_type2
	collapse (sum) person_years (sum) albuterol_equiv (sum) prednisone ///
			 (sum) any_cause_ed_total (sum) any_cause_hosp_total	  ///
			 , by(inhaler_type2)
	list inhaler_type2 person_years albuterol_equiv prednisone 	///
		 any_cause_ed_total any_cause_hosp_total
restore 


*--------------------------------------------------
* Pneumonia ED Visits and PNA Hospitalizations
*--------------------------------------------------

* Venn diagram 

preserve 

	* create a venn diagram for pneumonia ed visits and pneumonia hosps
	tab pneumonia_primary_ed_total pneumonia_primary_hosp_total

	gen pneumonia_primary_ed = pneumonia_primary_ed_total>0
	gen pneumonia_primary_hosp = pneumonia_primary_hosp_total>0

	tab pneumonia_primary_ed pneumonia_primary_hosp
	
	tab pneumonia_primary_ed pneumonia_primary_hosp if ///
		pneumonia_primary_ed==1 | pneumonia_primary_hosp==1, ro
	
	tab pneumonia_primary_ed pneumonia_primary_hosp if ///
		pneumonia_primary_ed==1 | pneumonia_primary_hosp==1, co
		
	tab pneumonia_primary_ed pneumonia_primary_hosp if ///
		pneumonia_primary_ed==1 | pneumonia_primary_hosp==1, cell	

	venndiag pneumonia_primary_ed pneumonia_primary_hosp

restore 


* ED visits occurring the day of or the day prior to 
	* (1) ANY hospitalization 
	* (2) PNA hospitalization


preserve 

	* use outcomes dataset
	use "switch_cohort_outcomes.dta", clear
	tab any_cause_ed
	tab any_cause_hosp
	tab any_cause_ed any_cause_hosp, m
	keep patienticn datevalue ed_arrival_date hosp_admitdate_from_ed admit_from_ed any_cause_ed pneumonia_primary_ed any_cause_hosp pneumonia_primary_hosp admitdate
	drop if (missing(ed_arrival_date) & missing(admitdate))
	count

	* (1) PNA ED visit on same day/day prior to any hospitalization 

	gen pnaed_sameday_anyhosp = 0
	replace pnaed_sameday_anyhosp = 1 if ((ed_arrival_date==admitdate) &  pneumonia_primary_ed==1)

	gen pnaed_dayprior_anyhosp = 0 
	gsort patienticn -datevalue
	by patienticn: replace pnaed_dayprior_anyhosp = 1 if ((ed_arrival_date-admitdate[_n-1]==-1) & pneumonia_primary_ed==1)
	sort patienticn datevalue

	* (2) PNA ED visit on same day/day prior to PNA hospitalization
	
	gen pnaed_sameday_pnahosp = 0
	replace pnaed_sameday_pnahosp = 1 if ((ed_arrival_date==admitdate) &  pneumonia_primary_ed==1 & pneumonia_primary_hosp==1)

	gen pnaed_dayprior_pnahosp = 0 
	gsort patienticn -datevalue
	by patienticn: replace pnaed_dayprior_pnahosp = 1 if ((ed_arrival_date-admitdate[_n-1]==-1) & pneumonia_primary_ed==1 & pneumonia_primary_hosp[_n-1]==1)
	sort patienticn datevalue

	* save dataset to merge with switch_alt33_analytic dataset
	keep patienticn datevalue pnaed_sameday_anyhosp pnaed_dayprior_anyhosp pnaed_sameday_pnahosp pnaed_dayprior_pnahosp
	save switch_alt33_pnaED, replace

	* because the current dataset includes discrete dates that must be merged 
	* between a range of dates in the switch_alt33_analytic dataset, 
	* given merging limitation in stata, export dataset to code in SAS using 
	* proc sql 
	use switch_alt33_analytic, clear
	keep patrecnum patienticn start end 
	save switch_alt33_cohort, replace
	
	* open the merged dataset 
	use pna_cohort_merge, clear
	
	collapse (sum) pnaed_dayprior_anyhosp (sum) pnaed_dayprior_pnahosp  ///
			 (sum) pnaed_sameday_anyhosp (sum) pnaed_sameday_pnahosp ///
			 , by (patienticn start end) 
	
	* create new variable that is same day or prior day 
		* any hosp 
		gen pnaed_sameprior_anyhosp = 0
		replace pnaed_sameprior_anyhosp = 1 if pnaed_dayprior_anyhosp==1 | pnaed_sameday_anyhosp==1
		tab pnaed_sameprior_anyhosp
	
		* pna hosp 
		gen pnaed_sameprior_pnahosp = 0
		replace pnaed_sameprior_pnahosp = 1 if pnaed_dayprior_pnahosp==1 | pnaed_sameday_pnahosp==1
		tab pnaed_sameprior_pnahosp
	
	drop pnaed_dayprior_anyhosp pnaed_dayprior_pnahosp pnaed_sameday_anyhosp pnaed_sameday_pnahosp
	
	* merge with 33% dataset 
	merge 1:1 patienticn start end using switch_alt33_analytic
	
	* number of pneumonia ed visits with any hospitalization on the same 
	* day or day prior 
	tab pneumonia_primary_ed_total
	gen pneumonia_primary_ed = pneumonia_primary_ed_total>0	
	tab pneumonia_primary_ed pnaed_sameprior_anyhosp, row
	
	* number of pneumonia ed visits with a pneumonia hospitalization on the same 
	* day or day prior 
	tab pneumonia_primary_ed pnaed_sameprior_pnahosp, row

	* check number of pneumonia hospitalizations 
	gen pneumonia_primary_hosp = pneumonia_primary_hosp_total>0
	tab pneumonia_primary_hosp
		
	* how many patients have a pneumonia ed visit and any cause/pneumonia
	* hospitalization?
	collapse (sum) pneumonia_primary_ed (sum) pneumonia_primary_hosp ///
			 (sum) pnaed_sameprior_anyhosp (sum) pnaed_sameprior_pnahosp ///
			 , by(patienticn)
	
	foreach var in 	pneumonia_primary_ed pneumonia_primary_hosp ///
					pnaed_sameprior_anyhosp pnaed_sameprior_pnahosp {
			gen `var'_ind = `var'>0
	}
	
	foreach var in 	pneumonia_primary_ed pneumonia_primary_hosp ///
					pnaed_sameprior_anyhosp pnaed_sameprior_pnahosp {
			tab `var'_ind
	}
	
	tab pnaed_sameprior_anyhosp_ind pneumonia_primary_ed_ind, co
	tab pnaed_sameprior_pnahosp_ind pneumonia_primary_ed_ind, co

	tab pneumonia_primary_ed_ind pneumonia_primary_hosp_ind
	
	* keep only variables we need to idenfity patients to exclude for the 
	* sensitivity analysis below and save dataset to merge with analytic 
	* dataset 
	keep patienticn pnaed_sameprior_anyhosp_ind pnaed_sameprior_pnahosp_ind
	save pna_sensitivity_analysis, replace
		
restore 


* SCCS Sensitivity Analysis: 
* For pneumonia ED outcome, remove those with PNA ED visits and PNA Hosps
use switch_alt33_analytic, clear 
merge m:1 patienticn using pna_sensitivity_analysis, nogen

* Pneumonia ED 
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division if pnaed_sameprior_pnahosp_ind==0, ///
		fe i(patienticn) offset(loglength) irr					

			
		
*************************************************************************************
* Using 100% Grace Period Dataset for Sensitivity Analyses in Supplemental Table S4 *
*************************************************************************************
	
* open dataset 
use "switch_alt100_20241010", clear

* create exposure groups 
tab inhaler_type
tab inhaler_type, nol 

gen no_inhaler_exgrp = inhaler_type==0 
gen symbicort_exgrp = inhaler_type==1 
gen wixela_exgrp = inhaler_type==2 
gen other_exgrp = inhaler_type==3 

tab inhaler_type no_inhaler_exgrp, m
tab inhaler_type symbicort_exgrp, m
tab inhaler_type wixela_exgrp, m
tab inhaler_type other_exgrp, m 

* create loglength 
gen loglength = log(length)
tab length if loglength==.
recode loglength .=0 

sort patienticn start

* check number of intervals for each patient 
bysort patienticn (patrecnum): egen max_num_intervals = max(patrecnum)
tab max_num_intervals if patrecnum==1
drop if max_num_intervals==1
drop max_num_intervals

* Analysis will be on cohort 5 
tab cohort5
keep if cohort5==1
distinct patienticn


*--------------
* Cohort 5
*--------------

* Albuterol inhaler equivalents	
xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr

* Prednisone		
xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* Hospitalizations - All Cause
xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* Hospitalizations - Respiratory 
xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* Hospitalizations - Pneumonia 
xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* ED - All Cause
xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* ED - Respiratory 
xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr		

* ED - Pneumonia 
xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
		c.age ib14.quarter##i.address_census_division , ///
		fe i(patienticn) offset(loglength) irr				
					
*-------------------------------------------------------------------------------
* Create a subset of patients with limited 'no controller periods' 
* (<3 months total)		
*-------------------------------------------------------------------------------

* create indicator for interval gap 
gen inhaler_gap = inhaler_type==0

* calculate number of days during all inhaler gaps 
gen gap_days = datediff(start, end, "day") if inhaler_gap==1
bysort patienticn (patrecnum): egen gap_days_tot = total(gap_days)

* indentify patients with more than 3 months of inhaler gaps 
gen inhaler_gap_lt3m = gap_days_tot<=91

tab inhaler_gap_lt3m 
tab inhaler_gap_lt3m if patrecnum==1 

* identify the last record for each patient 
gsort patienticn -patrecnum 
by patienticn: gen lastrec = _n 
replace lastrec = 0 if lastrec>1
sort patienticn patrecnum 

	* look at the number of gap days for the last record of a patient w/ gaps 
	* to confirm that everyone has 91 or fewer gap days
	sum gap_days if lastrec==1 & inhaler_gap==1
	
* inspect relationship between prednisone and inhaler type among people with
* <3m gap 
tab inhaler_type prednisone if inhaler_gap_lt3m==1, ro

bysort patienticn (patrecnum): egen ever_prednisone = max(prednisone)
tab inhaler_type prednisone if inhaler_gap_lt3m==1 & ever_prednisone==1, ro

tab ever_prednisone if patrecnum==1 & inhaler_gap_lt3m==1 
tab ever_prednisone if inhaler_gap_lt3m==1 

* SCCS analysis
preserve 

	keep if inhaler_gap_lt3m==1
	count 
	tab patrecnum if patrecnum==1 
	
	*--------------
	* Cohort 5
	*--------------

	* Albuterol inhaler equivalents	
	xi: xtpoisson albuterol_equiv i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr

	* Prednisone		
	xi: xtpoisson prednisone i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr		
	
	* Hospitalizations - All Cause
	xi: xtpoisson any_cause_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr		

	* Hospitalizations - Respiratory 
	xi: xtpoisson respiratory_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr		

	* Hospitalizations - Pneumonia 
	xi: xtpoisson pneumonia_primary_hosp_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr		

	* ED - All Cause
	xi: xtpoisson any_cause_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr		

	* ED - Respiratory 
	xi: xtpoisson respiratory_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr		

	* ED - Pneumonia 
	xi: xtpoisson pneumonia_primary_ed_total i.wixela_exgrp i.other_exgrp i.no_inhaler_exgrp i.symbicort_exgrp  ///
			c.age ib14.quarter##i.address_census_division , ///
			fe i(patienticn) offset(loglength) irr				
			
restore	
		
		
log close		
				