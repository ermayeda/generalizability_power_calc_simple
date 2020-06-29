/*Create a local macro variable for "scalar_X = variable_X".
We will this local macro in the simulate command. This will pull in scalars  
fromthe data generation and analysis file and store each scalar as a variable
with the same name.*/
local simlist ""
foreach x in sample_n sample_pblack boot_mean_diff boot_SE_diff diff_observe  ///
	boot_diff_95CI_lb boot_diff_95CI_ub boot_diff_CI_inc_null { 
local simlist "`simlist' `x'=`x'"
}

simulate `simlist', ///
reps($B) seed(67208113): do "C:\Users\emayeda\Dropbox\ADRD_selection_grant\July_2020_resubmission\preliminary_analyses\aim3_multiple_do_files_exp\aim3powerv5_sample_exp.do"

foreach x in ///
sample_n sample_pblack boot_mean_diff boot_SE_diff diff_observe boot_diff_CI_inc_null {
summarize `x'
scalar mean_`x' = r(mean)
}


/*Save results. One row = one iteration of sample generation*/	   
*Excel file
export excel using sim_results_each_replication_p_black_$p_black, sheet("Results") sheetmodify firstrow(variables)

*Stata data file
save "sim_results_each_replication_p_black_$p_black.dta", replace


/**************************************************************************************************/
/***	summarize simulation results															***/
/***	results summarized across different Excel sheets										***/
/**************************************************************************************************/
tabstat sample_n sample_pblack boot_mean_diff boot_SE_diff diff_observe boot_diff_CI_inc_null, stat(n mean sd min max) save
matlist r(StatTotal)
matrix results = r(StatTotal)'
matlist results

putexcel set sim_results_summarized.xlsx, sheet(sheet1) replace
putexcel set sim_results_summarized.xlsx, sheet(p_black_$p_black) modify
putexcel A1 = matrix(results), names 



/*putexcel A1=("mean sample size") B1=("P Black)") C1=("Mean bootstrap diff") D1=("Bootstrap SE diff") E1=("Mean diff obs") F1=("P 95% CIs inc. null") ///
		 A2=(mean_sample_n) B2=(mean_sample_pblack) C2=(mean_boot_mean_diff) D2=(mean_boot_SE_diff) E2=(mean_diff_observe) F1=(mean_boot_diff_CI_inc_null) /// 
file sim_results_summarized.xlsx saved*/
