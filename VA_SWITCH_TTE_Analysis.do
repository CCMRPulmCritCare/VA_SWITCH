version 18
set more off
cap log close
clear all
set linesize 80

cd ""

local c_date = c(current_date)
local date = subinstr("`c_date'", " ", "", .)

log using "Logs\va_switch_tte_`date'.log", replace

********************************************************************************
* Clinical Outcomes Following a National Inhaler Formulary Change
*	- TTE Sensitivity Analysis
* Author: Sarah Seelye
*
* Date Created: 2024 Dec 10		
* Last updated: 2025 Apr 4
********************************************************************************


***************************************
* Construct Data Set for TTE Analysis *
***************************************

*------------------------------------------------------------
* Merge in datasets to include necessary variables
*------------------------------------------------------------

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

** Smoking Status ** 
use "switch_cohort_smoking_status_20240325", clear

* create temp file to merge with strict dataset 
tempfile smoking 
save `smoking'


** MERGE ** 

* open strict dataset 
use "switch_20241030", clear

* merge with combat file
merge m:1 patienticn using `combat' 
tab patrecnum if patrecnum==1 & _merge==1
list patienticn if patrecnum==1 & _merge==1 

* recode 4 patients with missing combat status as not a combat veteran
replace ever_combatvet=0 if _merge==1

* drop patients not in analytic dataset 
drop if _merge==2

drop _merge

* merge with smoking file 
merge m:1 patienticn using `smoking'

* drop patients not in analytic data set 
drop if _merge==2
drop _merge

*--------------------------------------
* Prep dataset for TTE analysis
*--------------------------------------

* only keep patients with more than 1 interval 
bysort patienticn (patrecnum): egen max_num_intervals = max(patrecnum)
tab max_num_intervals if patrecnum==1
drop if max_num_intervals==1
drop max_num_intervals

*--------------------------------------
* Identify eligibility criteria
*--------------------------------------
	* dispensed controller inhaler during 180 days before & including 6/30/2021
	* dispensed controller inhaler during 180 days after & including 7/01/2021 
gen preswitchdate_str = "01jul2021"
gen preswitchdate = date(preswitchdate_str, "DMY")
format preswitchdate %td

gen days_before_preswitch = preswitchdate-datevalue

gen preperiod = inrange(days_before_preswitch, 1, 181) 
gen postperiod = inrange(days_before_preswitch, -180, 0) 

bysort patienticn (datevalue): egen preperiod_pat = max(preperiod)
bysort patienticn (datevalue): egen postperiod_pat = max(postperiod)

	* only keep patients who are in the pre and the post periods 
	* (180 days before and after July 1, 2021)
	count if patrecnum==1 & (preperiod_pat==0 | postperiod_pat==0)
	count if patrecnum==1 & (preperiod_pat==1 | postperiod_pat==1)
	
	drop if (preperiod_pat==0 | postperiod_pat==0) //255243
	count if patrecnum==1 & (preperiod_pat==1 | postperiod_pat==1) //327580
	
	* identify patients who had a controller inhaler during each of the time periods 
	gen controller_inhaler = inhaler_type>0
	gen preperiod_inhaler = 0
	replace preperiod_inhaler = 1 if controller_inhaler==1 & preperiod==1
	
	gen postperiod_inhaler = 0
	replace postperiod_inhaler = 1 if controller_inhaler==1 & postperiod==1
	
	bysort patienticn (datevalue): egen preperiod_inhaler_pat = max(preperiod_inhaler)
	bysort patienticn (datevalue): egen postperiod_inhaler_pat = max(postperiod_inhaler)

	* only keep patients who had a controller inhaler during both pre and post 
	* periods 
	count if patrecnum==1 & (preperiod_inhaler_pat==0 | postperiod_inhaler_pat==0)
	count if patrecnum==1 & (preperiod_inhaler_pat==1 | postperiod_inhaler_pat==1)
	
	drop if (preperiod_inhaler_pat==0 | postperiod_inhaler_pat==0) //1,608,087
	count if patrecnum==1 & (preperiod_inhaler_pat==1 | postperiod_inhaler_pat==1) //276369

*-------------------------------------	
* identify treatment strategies 
*-------------------------------------
	* 1. Transition to Wixela: dispensed Wixela on or after July 1, 2021 and 
	*		for the next 180 days 
	* 2. Continuation of non-Wixela inhaler: dispensed a non-Wixela controller 
	*		during switch/post-period (7/1/2021 + 180) AND not dispensed Wixela
	
gen mkg_tx1 = 0
replace mkg_tx1 = 1 if inhaler_type==2 & postperiod==1
bysort patienticn (datevalue): egen switch_to_wixela = max(mkg_tx1)	

gen mkg_tx2 = 0
replace mkg_tx2 = 1 if inlist(inhaler_type, 1, 3) & postperiod==1
replace mkg_tx2 = 0 if switch_to_wixela==1 
bysort patienticn (datevalue): egen stay_on_nonwixela = max(mkg_tx2)	


*----------------------
* Enrollment date
*----------------------

* Identify the enrollment date for each treatment group (date of first inhaler 
* dispensed during the 180 day post-period for each group)

preserve 
	keep if patrecnum==1
	keep patienticn switch_to_wixela stay_on_nonwixela
	count //276369
	
	* merge with full inhaler prescription dataset 
	merge 1:m patienticn using "P:\ORD_Prescott_202306018D\Cohort Build\Sarah\Data\Stata\final_rx_cohort_20240105.dta"
	drop if _merge==2 
	drop _merge
	distinct patienticn //276369
	
	* identifying inhaler releasedatetimes
	sort patienticn releasedatetime
	gen inhaler_releasedatetime = releasedatetime 
	local period strpos(inhaler_releasedatetime, ".")
	gen releasedatetime2 = trim(cond(`period', substr(inhaler_releasedatetime, 1, `period' -1), inhaler_releasedatetime))

	gen double releasedatetime3 = clock(releasedatetime2, "YMDhms")
	format releasedatetime3 %tc

	gen double releasedate = dofc(releasedatetime3)
	format releasedate %td

	drop releasedatetime2 releasedatetime3

	rename releasedate inhaler_releasedate
	order inhaler_releasedate, after(inhaler_releasedatetime)
	sort patienticn inhaler_releasedate
	
	* identify switch post period 
	gen preswitchdate_str = "01jul2021"
	gen preswitchdate = date(preswitchdate_str, "DMY")
	format preswitchdate %td

	gen days_before_preswitch = preswitchdate-inhaler_releasedate

	gen postperiod = inrange(days_before_preswitch, -180, 0) //range is 0 to -180 becuase postswitch period includes july 1 and 180 days afterwards (neg number b/c of calculation of 'days_before_preswitch')

	* only keep records during the post period to identify enrollment date 
	keep if postperiod == 1
	distinct patienticn //263079
	
	* identify the earliest releasedate for wixela 
	gsort patienticn -switch_to_wixela inhaler_releasedate
	by patienticn: gen wixela_postperiod_num = _n if switch_to_wixela==1
	gen wixela_enrollment_date = inhaler_releasedate if wixela_postperiod_num==1
	format wixela_enrollment_date %td
	replace wixela_enrollment_date = wixela_enrollment_date[_n-1] if wixela_enrollment_date==. & switch_to_wixela==1
	format wixela_enrollment_date %td

	* identify the earliest releasedate for non-Wixela 
	gsort patienticn -stay_on_nonwixela inhaler_releasedate
	by patienticn: gen nonwixela_postperiod_num = _n if stay_on_nonwixela==1
	gen nonwixela_enrollment_date = inhaler_releasedate if nonwixela_postperiod_num==1
	format nonwixela_enrollment_date %td
	replace nonwixela_enrollment_date = nonwixela_enrollment_date[_n-1] if nonwixela_enrollment_date==. & stay_on_nonwixela==1
	format nonwixela_enrollment_date %td
	
	* identify the first inhaler post period 
	gen first_controller_postperiod = inlist(1, wixela_postperiod_num, nonwixela_postperiod_num)
	tab first_controller_postperiod //263079
	keep if first_controller_postperiod==1
	
	* identify enrollment date 
	gen enrollment_date = wixela_enrollment_date
	replace enrollment_date = nonwixela_enrollment_date if missing(enrollment_date)
	format enrollment_date %td
	
	* identify enrollment inhaler 
	gen enrollment_inhaler = .
	replace enrollment_inhaler = 1 if !missing(nonwixela_enrollment_date)
	replace enrollment_inhaler = 2 if !missing(wixela_enrollment_date)
	lab def enrollment_inhaler 1 "Non-Wixela" 2 "Wixela"
	lab val enrollment_inhaler enrollment_inhaler
	tab enrollment_inhaler
	tab enrollment_inhaler switch_to_wixela, m
	tab enrollment_inhaler stay_on_nonwixela, m

	* keep variables needed for merge 
	keep patienticn enrollment_inhaler enrollment_date
	order patienticn enrollment_inhaler enrollment_date
	
	* save tempfile to merge back with dataset 
	tempfile enrollment 
	save `enrollment'
	
restore 

merge m:1 patienticn using `enrollment'	
drop _merge
distinct patienticn //276369

* drop if enrollment inhaler date is missing - these patients neither switched 
* nor stayed on symbicort for the 180 day period; they should not be in either 
* treatment group
tab enrollment_inhaler switch_to_wixela
tab enrollment_inhaler stay_on_nonwixela

count if missing(enrollment_date)
count if missing(enrollment_inhaler)
	
drop if enrollment_inhaler == .

*--------------------------	
* Identify Outcomes 
*--------------------------
	* At 90- and 180-days post-enrollment (date of first controller inhaler dispensed, on or after 7/01/2021)
		*	1. Mortality
		*	2. Albuterol fills (count)
		*	3. Prednisone fills (count)
		* 	4. ED visits (all cause, resp, pna)
		* 	5. Hosps (all cause, resp, pna)


sort patienticn patrecnum

gen enrollment_plus90_date = enrollment_date + 90
format enrollment_plus90_date %td

gen enrollment_plus180_date = enrollment_date + 180
format enrollment_plus180_date %td

* mortality
gen mort90_enroll = inrange(dod, enrollment_date, enrollment_plus90_date)
gen mort180_enroll = inrange(dod, enrollment_date, enrollment_plus180_date)

tab mort90_enroll
tab mort180_enroll

* for the remaining outcomes, use the outcomes dataset.
* save a temporary file for merging back in later 
tempfile mkgtte 
save `mkgtte'

* keep patienticns for merging with outcomes 
keep if patrecnum==1
count 

* merge with outcomes
keep patienticn enrollment_date enrollment_plus90_date enrollment_plus180_date
merge 1:m patienticn using "switch_cohort_outcomes_20241010.dta"

* drop patienticns not in tte cohort 
drop if _merge==2
drop _merge

* create 90-day and 180-day enrollment indicators around the datevalue of each row of data
gen dateofservice_enroll90_ind = inrange(datevalue, enrollment_date, enrollment_plus90_date) 
gen dateofservice_enroll180_ind = inrange(datevalue, enrollment_date, enrollment_plus180_date) 

* create new respiratory indicators for primary copd, asthma, pneumonia 
gen respiratory_primary_ed = inlist(1, asthma_primary_ed, copd_primary_ed, pneumonia_primary_ed)
gen respiratory_primary_hosp = inlist(1, asthma_primary_hosp, copd_primary_hosp, pneumonia_primary_hosp)

* albuterol fills 
gsort patienticn -dateofservice_enroll90_ind datevalue
by patienticn: egen albuterol_enroll90_tot = total(albuterol_inhaler_equiv_daily) if dateofservice_enroll90_ind==1
by patienticn: egen albuterol_enroll180_tot = total(albuterol_inhaler_equiv_daily) if dateofservice_enroll180_ind==1

* prednisone fills 
by patienticn: egen prednisone_enroll90_tot = total(prednisone_daily) if dateofservice_enroll90_ind==1
by patienticn: egen prednisone_enroll180_tot = total(prednisone_daily) if dateofservice_enroll180_ind==1

* any cause ed
by patienticn: egen edanycause_enroll90_tot = total(any_cause_ed) if dateofservice_enroll90_ind==1
by patienticn: egen edanycause_enroll180_tot = total(any_cause_ed) if dateofservice_enroll180_ind==1

* respiratory ed 
by patienticn: egen edresp_enroll90_tot = total(respiratory_primary_ed) if dateofservice_enroll90_ind==1
by patienticn: egen edresp_enroll180_tot = total(respiratory_primary_ed) if dateofservice_enroll180_ind==1

* pneumonia ed 
by patienticn: egen edpneu_enroll90_tot = total(pneumonia_primary_ed) if dateofservice_enroll90_ind==1
by patienticn: egen edpneu_enroll180_tot = total(pneumonia_primary_ed) if dateofservice_enroll180_ind==1

* any cause hospitalization 
by patienticn: egen hospanycause_enroll90_tot = total(any_cause_hosp) if dateofservice_enroll90_ind==1
by patienticn: egen hospanycause_enroll180_tot = total(any_cause_hosp) if dateofservice_enroll180_ind==1

* respiratory hospitalization
by patienticn: egen hospresp_enroll90_tot = total(respiratory_primary_hosp) if dateofservice_enroll90_ind==1
by patienticn: egen hospresp_enroll180_tot = total(respiratory_primary_hosp) if dateofservice_enroll180_ind==1

* pneumonia hospitalization 
by patienticn: egen hosppneu_enroll90_tot = total(pneumonia_primary_hosp) if dateofservice_enroll90_ind==1
by patienticn: egen hosppneu_enroll180_tot = total(pneumonia_primary_hosp) if dateofservice_enroll180_ind==1

* create patient-level indicators and merge back to main dataset 
foreach var in 	albuterol_enroll90 albuterol_enroll180 ///
				prednisone_enroll90 prednisone_enroll180 ///
				edanycause_enroll90 edanycause_enroll180 ///
				edresp_enroll90 edresp_enroll180 ///
				edpneu_enroll90 edpneu_enroll180 ///
				hospanycause_enroll90 hospanycause_enroll180 /// 
				hospresp_enroll90 hospresp_enroll180 ///
				hosppneu_enroll90 hosppneu_enroll180	{
	bysort patienticn (datevalue): egen `var'_pat = max(`var'_tot)
}

keep patienticn albuterol_enroll90_pat albuterol_enroll180_pat 	///
				prednisone_enroll90_pat prednisone_enroll180_pat ///
				edanycause_enroll90_pat edanycause_enroll180_pat ///
				edresp_enroll90_pat edresp_enroll180_pat ///
				edpneu_enroll90_pat edpneu_enroll180_pat ///
				hospanycause_enroll90_pat hospanycause_enroll180_pat ///
				hospresp_enroll90_pat hospresp_enroll180_pat ///
				hosppneu_enroll90_pat hosppneu_enroll180_pat

foreach var in 	albuterol_enroll90_pat albuterol_enroll180_pat 	///
				prednisone_enroll90_pat prednisone_enroll180_pat ///
				edanycause_enroll90_pat edanycause_enroll180_pat ///
				edresp_enroll90_pat edresp_enroll180_pat ///
				edpneu_enroll90_pat edpneu_enroll180_pat ///
				hospanycause_enroll90_pat hospanycause_enroll180_pat ///
				hospresp_enroll90_pat hospresp_enroll180_pat ///
				hosppneu_enroll90_pat hosppneu_enroll180_pat  {
	count if missing(`var')
	recode `var' .=0			
}
				
duplicates drop				
count  

tempfile outcomes 
save `outcomes'

use `mkgtte', clear 
merge m:1 patienticn using `outcomes'
drop _merge
				
*------------------------------------
* Create covariates for models
*------------------------------------
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
* preswitchdate;  keep only records that fall within the preswitch year  
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

* creater never & former smoker category 
gen never_smoker = ever_smoker==0
gen former_smoker = ever_smoker==1 & current_smoker==0


*--------------------------------------------
* Set up data set for weighting with CEM 
*--------------------------------------------

* keep only variables needed 
keep patienticn male age_at_switch race_3cat ever_combatvet ///
	 current_smoker ever_smoker never_smoker former_smoker ///
	 region_at_switch asthma copd ///
	 switch_to_wixela stay_on_nonwixela /// 
	 mort90_enroll mort180_enroll ///
	 albuterol_enroll90_pat albuterol_enroll180_pat ///
	 prednisone_enroll90_pat prednisone_enroll180_pat ///
	 edanycause_enroll90_pat edanycause_enroll180_pat ///
	 edresp_enroll90_pat edresp_enroll180_pat ///
	 edpneu_enroll90_pat edpneu_enroll180_pat ///
	 hospanycause_enroll90_pat hospanycause_enroll180_pat ///
	 hospresp_enroll90_pat hospresp_enroll180_pat ///
	 hosppneu_enroll90_pat hosppneu_enroll180_pat ///
	 prednisone_yr_preswitch albuterol_yr_preswitch ///
	 edanycause_yr_preswitch edresp_yr_preswitch ///
	 edpneu_yr_preswitch hospanycause_yr_preswitch ///
	 hospresp_yr_preswitch hosppneu_yr_preswitch ///
	 prednisone_yr_preswitch_ind albuterol_yr_preswitch_ind ///
	 edanycause_yr_preswitch_ind edresp_yr_preswitch_ind ///
	 edpneu_yr_preswitch_ind hospanycause_yr_preswitch_ind ///
	 hospresp_yr_preswitch_ind hosppneu_yr_preswitch_ind 

duplicates drop
count  //276369

* create categories of albuterol use in the year pre-switch 
gen albuterol_yr_preswitch_rnd = round(albuterol_yr_preswitch)
gen albuterol_yr_preswitch_0or1 = albuterol_yr_preswitch_rnd<=1 if !missing(albuterol_yr_preswitch_rnd)
gen albuterol_yr_preswitch_2 = albuterol_yr_preswitch_rnd==2 if !missing(albuterol_yr_preswitch_rnd)
gen albuterol_yr_preswitch_3plus = albuterol_yr_preswitch_rnd>=3 if !missing(albuterol_yr_preswitch_rnd)

gen albuterol_yr_preswitch_cat = .
replace albuterol_yr_preswitch_cat = 1 if albuterol_yr_preswitch_0or1==1
replace albuterol_yr_preswitch_cat = 2 if albuterol_yr_preswitch_2==1
replace albuterol_yr_preswitch_cat = 3 if albuterol_yr_preswitch_3plus==1
label def albuterol_yr_preswitch_cat 1 "0-1 fills yr prior" ///
		2 "2 fills yr prior" 3 "3+ fills yr prior"
label val albuterol_yr_preswitch_cat albuterol_yr_preswitch_cat		
tab albuterol_yr_preswitch_cat

* create smoking status 
gen smoking_status = .
replace smoking_status = 1 if current_smoker==1
replace smoking_status = 2 if former_smoker==1
replace smoking_status = 3 if never_smoker==1
lab def smoking_status 1 "current" 2 "former" 3 "never"
lab val smoking_status smoking_status
tab smoking_status

* create a single variable for treatment 
tab switch_to_wixela stay_on_nonwixela, m
gen tx_switched = switch_to_wixela
tab tx_switched switch_to_wixela

* examine treatment variable 
tab tx_switched

* create binary indicators for outcome variables used in regressions 
foreach var in 	 edanycause_enroll90 edanycause_enroll180 ///
				 edresp_enroll90 edresp_enroll180 ///
				 edpneu_enroll90 edpneu_enroll180 ///
				 hospanycause_enroll90 hospanycause_enroll180 ///
				 hospresp_enroll90 hospresp_enroll180 ///
				 hosppneu_enroll90 hosppneu_enroll180 {
		tab `var'_pat, m
		gen `var'_ind = `var'_pat>0
		tab `var'_pat `var'_ind
	}
	
order edanycause_enroll90_ind- hosppneu_enroll180_ind, after(hosppneu_enroll180_pat)

* check global balance for the unmatched data 
imb age_at_switch male asthma copd smoking_status 	///
	region_at_switch hospanycause_yr_preswitch_ind  ///
	edanycause_yr_preswitch_ind albuterol_yr_preswitch_cat		  ///
	prednisone_yr_preswitch_ind, treatment(tx_switched)

* matching algorithm 
set seed 4893
cem age_at_switch male asthma copd smoking_status 	///
	region_at_switch hospanycause_yr_preswitch_ind  ///
	edanycause_yr_preswitch_ind albuterol_yr_preswitch_cat		  ///
	prednisone_yr_preswitch_ind, treatment(tx_switched)
		/* used cem's automated binning for continuous variables */

tab cem_matched tx_switched, co 
tab cem_matched if tx_switched==1 //95.9%
tab tx_switched if cem_weights>0

*-------------------------------------------
* eTable: Covariate balance and SMDs
*-------------------------------------------

* Note: Table completed in R in order to incorporate weights. stddiff doesn't 
* allow weights, and covbal won't produce accurate SMDs for categorical 
* variables with more than 2 categories 

* See R code saved here:
	* I:\VA-SWITCH\5. Identifiable Data\Sarah\Data\VINCI\Do\Identified\tte_covbal_smd_20250402

covbal tx_switched age_at_switch , wt(cem_weights)
covbal tx_switched male , wt(cem_weights)
covbal tx_switched asthma , wt(cem_weights)
covbal tx_switched copd , wt(cem_weights)
covbal tx_switched hospanycause_yr_preswitch_ind , wt(cem_weights)
covbal tx_switched edanycause_yr_preswitch_ind , wt(cem_weights)
*covbal tx_switched i.albuterol_yr_preswitch_cat , wt(cem_weights)
covbal tx_switched prednisone_yr_preswitch_ind , wt(cem_weights)
*covbal tx_switched i.smoking_status , wt(cem_weights)
*covbal tx_switched i.region_at_switch , wt(cem_weights)

* save dataset 
sort patienticn
save "switch_aim1_tte_20241218.dta", replace 


*-----------------------------------------		
* Unadjusted models using CEM weights		
*-----------------------------------------	

* Open data set for analysis
use Data\switch_aim1_tte_20241218, clear

* Numbers of those switched and not switched used in TTE analysis
tab tx_switched
tab tx_switched if cem_weights>0
			
* mortality 
logistic mort90_enroll i.tx_switched [pweight=cem_weights]
margins tx_switched
margins, dydx(tx_switched)

logistic mort180_enroll i.tx_switched [pweight=cem_weights]	
margins tx_switched 
margins, dydx(tx_switched)

* albuterol	
poisson albuterol_enroll90_pat i.tx_switched [pweight=cem_weights], irr
margins tx_switched
margins, dydx(tx_switched)

poisson albuterol_enroll180_pat i.tx_switched [pweight=cem_weights], irr	
margins tx_switched
margins, dydx(tx_switched)

* prednisone	
poisson prednisone_enroll90_pat i.tx_switched [pweight=cem_weights], irr
margins tx_switched

poisson prednisone_enroll180_pat i.tx_switched [pweight=cem_weights], irr
margins tx_switched
margins, dydx(tx_switched)

* ed any cause
logistic edanycause_enroll90_ind i.tx_switched [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

logistic edanycause_enroll180_ind i.tx_switched [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

* ed respiratory 
logistic edresp_enroll90_ind i.tx_switched [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

logistic edresp_enroll180_ind i.tx_switched [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

* ed pneumonia
logistic edpneu_enroll90_ind i.tx_switched [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

logistic edpneu_enroll180_ind i.tx_switched [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

* hospitalization any cause 
logistic hospanycause_enroll90_ind i.tx_switched [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)
 
logistic hospanycause_enroll180_ind i.tx_switched [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

* hosp respiratory
logistic hospresp_enroll90_ind i.tx_switched [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)
 
logistic hospresp_enroll180_ind i.tx_switched [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

* hosp pna 
logistic hosppneu_enroll90_ind i.tx_switched [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

logistic hosppneu_enroll180_ind i.tx_switched [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)


*-----------------------------------------		
* Adjusted models using CEM weights		
*-----------------------------------------	
local covar age_at_switch male asthma copd i.smoking_status ///
			i.region_at_switch hospanycause_yr_preswitch_ind  ///
			edanycause_yr_preswitch_ind albuterol_yr_preswitch_cat		  ///
			prednisone_yr_preswitch_ind


* mortality 
logistic mort90_enroll i.tx_switched `covar' [pweight=cem_weights]
margins tx_switched
margins, dydx(tx_switched)

logistic mort180_enroll i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched 
margins, dydx(tx_switched)
	
* albuterol	
poisson albuterol_enroll90_pat i.tx_switched `covar' [pweight=cem_weights], irr
margins tx_switched
margins, dydx(tx_switched)

poisson albuterol_enroll180_pat i.tx_switched `covar' [pweight=cem_weights], irr	
margins tx_switched
margins, dydx(tx_switched)

* prednisone	
poisson prednisone_enroll90_pat i.tx_switched `covar' [pweight=cem_weights], irr
margins tx_switched
margins, dydx(tx_switched)

poisson prednisone_enroll180_pat i.tx_switched `covar' [pweight=cem_weights], irr
margins tx_switched
margins, dydx(tx_switched)

* ed any cause
logistic edanycause_enroll90_ind i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)
 
logistic edanycause_enroll180_ind i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

* ed respiratory 
logistic edresp_enroll90_ind i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

logistic edresp_enroll180_ind i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

* ed pneumonia
logistic edpneu_enroll90_ind i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)
 
logistic edpneu_enroll180_ind i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

* hospitalization any cause 
logistic hospanycause_enroll90_ind i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)
 
logistic hospanycause_enroll180_ind i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

* hosp respiratory
logistic hospresp_enroll90_ind i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)
 
logistic hospresp_enroll180_ind i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

* hosp pna 
logistic hosppneu_enroll90_ind i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)

logistic hosppneu_enroll180_ind i.tx_switched `covar' [pweight=cem_weights]	
margins tx_switched
margins, dydx(tx_switched)


log close
