
			***** Wildfire smoke and Covid-19 case fatality ratios in the San Francisco Bay Area: a synthetic control analysis *****
			**** Version: June 10th, 2021 ****
			*** Code 1/3: Compliing data ***
		      


*===============================================================================
* Import the mortality data 
*===============================================================================

	* You can set the working directories in advance with global macros:
	
	//Anna's working directories
	global datain "C:\Users\annak\Dropbox\Projects\2020_San Diego\Wildfire smoke and COVID-19\Data"     
	global dataout "C:\Users\annak\Dropbox\Projects\2020_San Diego\Wildfire smoke and COVID-19\Data\Data_for_analysis"
	global results "C:\Users\annak\Dropbox\Projects\2020_San Diego\Wildfire smoke and COVID-19\Results"

	//Lara's working directories	
	global datain "F:\Lara\Projects\Wildfire smoke\Data"                                                 
	global dataout "F:\Lara\Projects\Wildfire smoke\Data_for_analysis"
	global results "F:\Lara\Projects\Wildfire smoke\Results"
				
	
	
	import delimited "$datain\covid\covid_deaths_usafacts.csv", clear encoding(UTF-8) 

		reshape long v, i(countyname state) j(day)  

		gen date= 21931 + day	
		gen date2= date
		format date2 %d
		 
		sort state countyname day
		 
		rename v total_d
		  
		by state countyname: gen death=total_d-total_d[_n-1]                    //daily covid deaths

		replace death=total_d if date==21936
		
		
		
		sum death if death < 0
		list state day countyfips date date2 total_d death if death < 0
		


	save "$datain\covid\Covid_deaths_by_county.dta", replace

*===============================================================================
* Import the case data 
*===============================================================================

	import delimited "$datain\covid\covid_confirmed_usafacts.csv", clear encoding(UTF-8) 

		reshape long v, i(countyname state) j(day)  

		gen date= 21931 + day
		gen date2= date
		format date2 %d
		 
		sort state countyname day
		 
		rename v total_cases
		 
		by state countyname: gen cases = total_cases-total_cases[_n-1]           //daily covid cases 

		replace cases=total_cases if date==21936
		
		
		sum cases if cases < 0
		list state day countyfips total_cases date date2 total_cases cases if cases < 0
		
		
		gen cases_neg=.
		replace cases_neg=1 if cases<0
		replace cases_neg=0 if cases>=0
		
		by state countyname: replace cases= cases[_n-1] if cases_neg==1
		by state countyname: replace cases= cases[_n-1] if cases_neg==1

		replace cases=. if cases_neg==1
			 
		by state countyname: gen cases_2wks = cases[_n-14]                       //daily covid cases lagged 2 weeks
		by state countyname: gen cases_avg7_14 = (cases[_n-14] +  ///              //rolling mean of daily covid cases (7 to 14 days earlier)
												  cases[_n-13] +  ///
												  cases[_n-12] +  ///
												  cases[_n-11] +  ///
												  cases[_n-10] +  ///
												  cases[_n-9 ] +  ///
												  cases[_n-8 ])/7 

	save "$datain\covid\Covid_cases_by_county.dta", replace


*===============================================================================
* Merge the covid deaths and hospitalization data 
*===============================================================================
	
	 use "$datain\covid\Covid_cases_by_county.dta", clear

		 replace countyname="Alexandria city" if countyname=="Alexandria City"
		 replace countyname="Broomfield County" if countyname=="Broomfield County and City"
		 replace countyname="Charlottesville city" if countyname=="Charlottesville City"
		 replace countyname="Chesapeake city" if countyname=="Chesapeake City"
		 replace countyname="Danville city" if countyname=="Danville City"
		 replace countyname="Fredericksburg city" if countyname=="Fredericksburg City"
		 replace countyname="Harrisonburg city" if countyname=="Harrisonburg City"
		 replace countyname="Lac qui Parle County" if countyname=="Lac Qui Parle County"
		 replace countyname="Manassas city" if countyname=="Manassas City"
		 replace countyname="Mathews County" if countyname=="Matthews County"
		 replace countyname="Norfolk city" if countyname=="Norfolk City"
		 replace countyname="Portsmouth city" if countyname=="Portsmouth City"
		 replace countyname="Richmond city" if countyname=="Richmond City"
		 replace countyname="Suffolk city" if countyname=="Suffolk City"
		 replace countyname="Virginia Beach city" if countyname=="Virginia Beach City"

		 merge 1:1 state countyname day using "$datain\covid\Covid_deaths_by_county.dta"
		 drop _merge
		 
	 save "$datain\covid\Covid_by_county.dta", replace

	 	 
*-------------------------Generate output measures------------------------------

		gen CFR_14dys=death/cases_2wks          		//daily deaths devided by diagnosed cases two weeks earlier
		
	
		sum CFR_14dys  if CFR_14dys<0
		list state day total_cases date2 cases total_cases cases death CFR_14dys if CFR_14dys < 0
		
		
		gen CFR_14dys_neg=.
		replace CFR_14dys_neg=1 if CFR_14dys<0
		replace CFR_14dys_neg=0 if CFR_14dys>=0
		replace CFR_14dys=. if CFR_14dys_neg==1      	
		
		gen CFR_same_day=death/cases                    //daily deaths devided by daily diagnosed cases
		gen CFR_same_day_neg=.
		replace CFR_same_day_neg=1 if CFR_same_day<0
		replace CFR_same_day_neg=0 if CFR_same_day>=0
		replace CFR_same_day=. if CFR_same_day_neg==1   //dropping negative ratios 

		gen CFR_7_14days=death/cases_avg7_14            //daily deaths devided by average daily diagnosed cases 7 to 14 days ago
		gen CFR_7_14days_neg=.
		replace CFR_7_14days_neg=1 if CFR_7_14days<0
		replace CFR_7_14days_neg=0 if CFR_7_14days>=0
		replace CFR_7_14days=. if CFR_7_14days_neg==1   //dropping negative ratios


	save "$datain\covid\Covid_by_county.dta", replace


	
*===============================================================================
* Prepare the smoke data (from R) 
*===============================================================================
	
	import delimited "$datain\Smoke_county_bin\smk_intsct_county_area_2020.csv", clear

		tostring dates, replace
		
		gen date2=date(dates, "YMD")
		gen date=date2
		format date %td
 
		rename geoid countyfips
		
		 save "$datain\Smoke\smoke_by_county", replace

		
		collapse date2, by( countyfips)
		replace date=22104
		
		save "$datain\Smoke\county_july8.dta", replace
		
		use "$datain\Smoke\smoke_by_county"
		
		append using "$datain\Smoke\county_july8.dta"

		sort countyfips date2
		
		by countyfips: replace medium_pct = (medium_pct[_n-1] + medium_pct[_n+1])/2 if medium_pct==.
		by countyfips: replace light_pct = (light_pct[_n-1] + light_pct[_n+1])/2 if light_pct==.
		by countyfips: replace heavy_pct = (heavy_pct[_n-1] + heavy_pct[_n+1])/2 if heavy_pct==.
		
		
		by countyfips: replace medium = (medium[_n-1] + medium[_n+1])/2 if medium==.
		by countyfips: replace light = (light[_n-1] + light[_n+1])/2 if light==.
		by countyfips: replace heavy = (heavy[_n-1] + heavy[_n+1])/2 if heavy==.


		replace date=date2 if date==.
  
   save "$datain\Smoke\smoke_by_county", replace

  *=============================================================================== 
* Merge the mobility data and save
*=============================================================================== 
	 
 import delimited "$datain\Mobility_data\Trips_by_Distance.csv" 
 keep if level=="County"
 
 gen date2=date( date, "YMD")
  drop if date2<21936  //restricting to same time period as other datasets Jan 22-Nov 10 2020
  drop if date2>22229
  
   rename date date3
   rename date2 date
   
   rename statepostalcode state

keep state countyname populationstayingathome populationnotstayingathome date
 
 
 save "$datain\Mobility_data\US_transportation_mobility.dta", replace
 
   

   
*=============================================================================== 
* Merge the health and smoke and mobility data and save
*=============================================================================== 
	
	use "$datain\covid\Covid_by_county.dta", clear

		drop if countyfips==0                        

		merge 1:1 countyfips date using "$datain\Smoke\smoke_by_county.dta"
		keep if _merge==3
		drop _merge
		
		merge 1:1 state countyname date using "F:\Lara\Projects\Wildfire smoke\Data\Mobility_data\US_transportation_mobility.dta"

		 drop if _merge==2
		 
		 drop _merge
		 
	save "$datain\Full_data\Covid_smoke_by_county.dta", replace
	
	
	
	*** Saved and used this data for further analysis (see R code "Diagnostics")
	saveold "$datain\Diagnostics\Covid_smoke_by_county_diagnostics.dta", replace version(12)


*=============================================================================== 
