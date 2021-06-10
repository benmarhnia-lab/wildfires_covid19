
 ***** Wildfire smoke and Covid-19 case fatality ratios in the San Francisco Bay Area: a synthetic control analysis *****
			**** Version: June 10th, 2021 ****
			*** Code 3/3: Data Analysis ***
		                         


ssc install synth, all
//cap ado uninstall synth_runner //in case already installed
net install synth_runner, from(https://raw.github.com/bquistorff/synth_runner/master/) replace
			
help synth
help synth_runner			

* Set working directories
	
//Anna's working directories
global datain "C:\Users\annak\Dropbox\Projects\2020_San Diego\Wildfire smoke and COVID-19\Data\Full_data"     
global dataout "C:\Users\annak\Dropbox\Projects\2020_San Diego\Wildfire smoke and COVID-19\Data\Data_for_analysis"
global results "C:\Users\annak\Dropbox\Projects\2020_San Diego\Wildfire smoke and COVID-19\Results"

//Lara's working directories
global datain "F:\Lara\Projects\Wildfire smoke\Data\Full_data"                                                 
global dataout "F:\Lara\Projects\Wildfire smoke\Data\Data_for_analysis"
global results "F:\Lara\Projects\Wildfire smoke\Results\Results_May"
				
			
*==============================================================================
* Prepare the data for analysis 
*============================================================================== 
	
	use "$datain\Covid_smoke_by_county_daily_all_revised.dta", clear
	
		rename heavy_smoke_p50 p_heavy50
		rename heavy_smoke_p70 p_heavy70
		rename heavy_smoke_p90	p_heavy90
		rename medium_smoke_p70 p_medium70
		
	save "$datain\Covid_smoke_by_county_daily_all_revised.dta", replace
		
		
		
		
	use "$datain\Covid_smoke_by_county_daily_all_revised.dta", clear
				
    	replace countyname = subinstr(countyname, " ", "", .)

		levelsof countyname if state=="CA" & eligible=="Yes" , local(Bay_area)  //create a list of eligible county names (to be used in the loop)        
		
		
		*** LOOPS START
		
		
				
		foreach exp in heavy70 heavy50 heavy90 medium70  {  
		
		foreach name of local Bay_area {              //loop through each eligible county                                                             
		
		foreach CFR in CFR14_lag7 CFR14_same_day {    //loop through each CFR measure                  
 
			use "$datain\Covid_smoke_by_county_daily_all_revised.dta", clear
			
 
				replace countyname = subinstr(countyname, " ", "", .)

				sort time
				gen day_smoke_treated = smoke_`exp'_day1 if countyname=="`name'" & state=="CA" & first_smoke_`exp'==1 //change based on percentile wanted
				replace day_smoke_treated = day_smoke_treated[_n-1] if missing(day_smoke_treated)
				gsort -time 
				replace day_smoke_treated = day_smoke_treated[_n-1] if missing(day_smoke_treated)
				sort id time

				gen X=.
				replace X=0 if time<day_smoke_treated
				replace X=1 if time>=day_smoke_treated
	 
				gen Y= `CFR'
	 			
				
				gen miss=.                                                      //remove donor counties with missing observations 
				replace miss=1 if missing(Y) 
				egen miss_by_county=mean(miss), by(countyname state)
				drop if miss_by_county==1 
				drop miss miss_by_county
				
				sort countyname state time
	 							
				gen miss=.                                                      //remove donor counties with missing mobility data
				rename populationstayingathome mobility
				replace miss=1 if missing(mobility) 
				egen miss_by_county=mean(miss), by(countyname state)
				drop if miss_by_county==1 
				drop miss miss_by_county
				
				gen outlier=.                                                   //remove donors with CFR>1
				replace outlier=1 if Y>1 & (countyname!="`name'" | state!="CA")
				egen outlier_by_county=mean(outlier), by(countyname state)
				drop if outlier_by_county==1 
				drop outlier outlier_by_county 
			
				gen t=day                                                         
				
				sort id time
	 
				gen treated=.
				replace treated=1 if countyname=="`name'" & state=="CA"
				replace treated=0 if (countyname!="`name'" | state!="CA")		

				egen avg_smoke=mean(p_`exp'), by(state countyname)              //remove donor counties which received more that 1% smoke (mostly CA counties)
				
				drop if avg_smoke>0.01 & treated!=1                             

				drop id
				egen id = group(countyname state)	
				
				egen day_one= min(t)                                            
																				
				gen countstate= countyname + state
				
				replace countstate = subinstr(countstate,"County", "", .)
				replace countstate = subinstr(countstate,".", "", .)
				replace countstate = subinstr(countstate,"'", "", .)																
																				
	 
				save "$dataout\Covid_smoke_by_county_daily_`exp'_`name'_`CFR'.dta", replace 
	
	
			}	
			
		    }
				
			}
			
			
  		*** LOOPS END

		
		
*==============================================================================
* Analysis 
*============================================================================== 

		*** LOOPS START    
		
		foreach exp in heavy90 medium70  heavy90 {
		
		foreach name in  AlamedaCounty   {   
		
		foreach CFR in CFR14_lag7 CFR14_same_day {   
 
 		macro drop Ypre week_first week_exp mobility* Y* day_first day_exp

			use "$dataout\Covid_smoke_by_county_daily_`exp'_`name'_`CFR'.dta", clear
			
							
				summarize day_one
				display r(min)
				global day_first = r(min)
				
				summarize day_smoke_treated
				display r(min)
				global day_exp = r(min)
				
				
								
				display $day_first
				display $day_exp
					
				forvalue i = $day_first / $day_exp {
						global Ypre $Ypre  Y(`i')  
						}
							
				forvalue i = $day_first / $day_exp {
						global mobility $mobility  mobility(`i')  
						}
			
				macro list		
				
				replace id=0 if treated==1
				
				replace Y= Y*1000
				

				tsset id t
						
				qui: synth_runner Y $Ypre $Xpre $mobility, trunit(0) trperiod($day_exp) pvals1s  gen_vars keep("$results\synth_runner\res_daily_`exp'_`name'_`CFR'") replace
				
				single_treatment_graphs
				graph save "$results\figures_gph\sing_treatm_daily_`exp'_`name'_`CFR'.gph", replace
				
				effect_graphs 
				graph save "$results\figures_gph\effect_graph_daily_`exp'_`name'_`CFR'.gph", replace
				
				
				outsheet using "$results\csv_files\results_daily_`exp'_`name'_`CFR'.csv" , comma replace

					
				
				*calculting confidence intervals
				use "$results\synth_runner\res_daily_`exp'_`name'_`CFR'.dta", clear

				mat pv = e(pvals_std_1s)
				mat pval = pv'
				mat list pval

				coefplot matrix(pv), vertical

				* Calculate 90% CIs / as in BMJ (2011) doi: https://doi.org/10.1136/bmj.d2090 

				mat effect = e(treat_control)
				svmat effect

				keep if t>=$day_exp & id==0

				svmat pval
				replace pval1=0.0001 if pval1==0
				gen lpval = ln(pval1)
				gen lpval1 =2.404*lpval
				gen lpval2 = 0.743 - lpval1
				gen lpval3 = sqrt(lpval2)
				gen zval = -0.862 + lpval3
				mkmat zval,nomissing

				gen tdiff = effect1-effect2
				gen tdiff2 = effect2-effect1
				svmat zval
				gen se = tdiff2/zval1
				gen upper_ci = tdiff - 1.645*se
				gen lower_ci = tdiff + 1.645*se

				gen ti = _n 
				tsset ti
				list upper_ci tdiff lower_ci ti

				

				save "$results\synth_runner\res_CI_daily_`exp'_`name'_`CFR'.dta", replace
				
				tsline upper_ci tdiff lower_ci, lcolor(grey black grey) lp(dash solid dash) title(`name') xtitle("Days after wildfire smoke event") ytitle("Difference in case fatality ratio")  legend(order (1 "90% CI" 2 "estimate")) 
				
				
				graph save  "$results\figures_gph\effect_CI_daily_`exp'_`name'_`CFR'.gph", replace
				
				
				
				clear
				
				cd "F:\Lara\Projects\Wildfire smoke\Results\Results_May\synth_control_units"


				
				
				use "$dataout\Covid_smoke_by_county_daily_`exp'_`name'_`CFR'.dta", clear
					
		
				
					
				display $day_first
				display $day_exp
					
				forvalue i = $day_first / $day_exp {
						global Ypre $Ypre  Y(`i')  
						}
							
				
				forvalue i = $day_first / $day_exp {
						global mobility $mobility  mobility(`i')  
							}
				
				
				macro list				
				
				replace id=0 if treated==1
				

				tsset id t
						
				
				synth Y $Ypre $Xpre $mobility , trunit(0) trperiod($day_exp) unitnames(countstate) keep("synth_daily_`exp'_`name'_`CFR'") replace

							
				
				
				** counterfactual fixed effect estimator
				
				use "$dataout\Covid_smoke_by_county_daily_`exp'_`name'_`CFR'.dta", clear
				
				display $day_first
				display $day_exp
					
				forvalue i = $day_first / $day_exp {
						global Ypre $Ypre  Y(`i')  
						}
								
				replace id=0 if treated==1
				
				replace Y= Y*1000
	
				tsset id t

				gen timing=t-$day_exp

				replace treated=0 if countyname=="`name'" & state=="CA" & timing<0


				fect Y, treat(treated) cov(mobility) unit(id) time(timing) method("fe") force("two-way")

				
				*effect_graphs
				graph save "$results\figures_gph\fece_daily_`exp'_`name'_`CFR'.gph", replace
				
				
			
		macro drop Ypre week_first week_exp
			
		
		}
		
		}
		
		}
		
  		*** LOOPS END
	
*==============================================================================
* Processing results
*==============================================================================				
foreach exp in  heavy70   {   
		
		foreach name in AlamedaCounty ContraCostaCounty MarinCounty SanFranciscoCounty SantaClaraCounty SonomaCounty  {   
		
		foreach CFR in CFR14_lag7 CFR14_same_day  {   

		use "$results\synth_control_units\synth_daily_`exp'_`name'_`CFR'"
		
		drop  _Y_treated _Y_synthetic _time
		
		keep if _W_Weight>0
		
		rename _Co_Number CountyState
		rename _W_Weight Weight
		
		save "$results\synth_control_units\weights_synth_daily_`exp'_`name'_`CFR'",  replace
		
		export excel using "$results\synth_weights\weights_daily_`exp'_`name'_`CFR'.xls" , replace
		}
		}
}				
				
*==============================================================================
* Combine all graphs
*==============================================================================				

* CFR_14days_new
				
		graph combine "$results\figures_gph\effect_graph_daily_heavy70_AlamedaCounty_CFR14_lag7"        ///
					  "$results\figures_gph\effect_graph_daily_heavy70_ContraCostaCounty_CFR14_lag7"    /// 	   
					  "$results\figures_gph\effect_graph_daily_heavy70_MarinCounty_CFR14_lag7"          /// 	  
					  "$results\figures_gph\effect_graph_daily_heavy70_SanFranciscoCounty_CFR14_lag7" /// 	  
					  "$results\figures_gph\effect_graph_daily_heavy70_SantaClaraCounty_CFR14_lag7"     /// 	   
					  "$results\figures_gph\effect_graph_daily_heavy70_SonomaCounty_CFR14_lag7",       ///
					  cols(6) title("") ycommon scheme(s1mono) commonscheme                   /// 	  
					  l1(CFR 7-14 days) b1(time) xsize(20) 
 
 graph combine    "$results\figures_gph\effect_graph_daily_heavy70_AlamedaCounty_CFR14_same_day"        ///
					  "$results\figures_gph\effect_graph_daily_heavy70_ContraCostaCounty_CFR14_same_day"    /// 	   
					  "$results\figures_gph\effect_graph_daily_heavy70_MarinCounty_CFR14_same_day"          /// 	  
					  "$results\figures_gph\effect_graph_daily_heavy70_SanFranciscoCounty_CFR14_same_day" /// 	  
					  "$results\figures_gph\effect_graph_daily_heavy70_SantaClaraCounty_CFR14_same_day"     /// 	   
					  "$results\figures_gph\effect_graph_daily_heavy70_SonomaCounty_CFR14_same_day",       ///
					  cols(6) title("") ycommon scheme(s1mono) commonscheme                   /// 	  
					  l1(CFR14 same day) b1(time) xsize(20) 
		
		
		
		graph export  "$results\figures_svg\CFR_14days_new_combined.svg", replace
		graph export  "$results\figures_png\CFR_14days_new_combined.png", replace			


graph combine "$results\figures_gph\effect_CI_daily_heavy70_AlamedaCounty_CFR14_lag7"        ///
					  "$results\figures_gph\effect_CI_daily_heavy70_ContraCostaCounty_CFR14_lag7"    /// 	   
					  "$results\figures_gph\effect_CI_daily_heavy70_SanFranciscoCounty_CFR14_lag7" /// 	  
					  "$results\figures_gph\effect_CI_daily_heavy70_SantaClaraCounty_CFR14_lag7"     /// 	   
					  "$results\figures_gph\effect_CI_daily_heavy70_SonomaCounty_CFR14_lag7",        ///
					  cols(6) title("COVID-19 Case Fatality Ratios (deaths/cases same day) after Wildfire Smoke in eligible Bay Area Counties") ycommon    scheme(s1color)                /// 	  
					  l1(CFR same day) b1(time) xsize(20) 
			

		graph combine "$results\figures_gph\fece_daily_heavy70_AlamedaCounty_CFR14_lag7"        ///
					  "$results\figures_gph\fece_daily_heavy70_ContraCostaCounty_CFR14_lag7"    /// 	   
					  "$results\figures_gph\fece_daily_heavy70_MarinCounty_CFR14_lag7"          /// 	  
					  "$results\figures_gph\fece_daily_heavy70_SanFranciscoCounty_CFR14_lag7" /// 	  
					  "$results\figures_gph\fece_daily_heavy70_SantaClaraCounty_CFR14_lag7"     /// 	   
					  "$results\figures_gph\fece_daily_heavy70_SonomaCounty_CFR14_lag7",       ///
					   title("") ycommon scheme(s1mono) commonscheme                   /// 	  
					  l1(CFR 7 day lag) b1(time) xsize(20) ysize(10)
					  
					  
			
*==============================================================================
*==============================================================================