/*Program to sample from population and estimate population effect based on one sample
1. Draw sample of n people from population (copy 1 -- "real" observations). Note: we will vary the composition of the sample for power curves.
2. Estimate regression model in sample & store regression coefficients.
3. Draw bootstrap samples of size n from the sample and estimate regression model in each bootstrapped sample & store coefficients.
4. Load population data again. 
	a. Apply stored beta coefficients from the "observed" sample to the 
	population to estimate population effect of exp on outcome based on the 
	observed sample & store this pop estimate (this is the point estimate of
	the population effect in this sample). 
	b. Apply stored beta coefficients from the bootstrap samples to the 
	population to estimate population effect of exp on outcome based on each
	bootstrap sample & store each pop estimate. The SD of the pop estiamtes
	in the bootstrapped samples is the bootstrap Se of the pop estimaate based
	on the sample. Store the bootstrap SE. 

*Ref for writing bootstrap program: 
https://stats.idre.ucla.edu/stata/faq/how-do-i-write-my-own-bootstrap-program/

*Required package: ssc install randomselect
*/



*Create global variable for number of bootstrapped samples
global boot_reps = 250 //------------------------------------------------------------> Consider defining this global variable in a different step.


*Specify racial/ethnic composition of the sample ------------------------------------> We are going to vary this for power curve
*if n=1,000
* Sample 0.5:  97.5% white (n=975), 2.5% black (n=25)
* Sample 1:  95% white (n=950), 5%   black (n=50)
* Sample 2:  90% white (n=900), 10%  black (n=100)
* Sample 3:  85% white (n=850), 15%  black (n=150)
* Sample 4:  80% white (n=800), 20%  black (n=200)
* Sample 5:  75% white (n=750), 25%  black (n=250)
* Sample 6:  70% white (n=700), 30%  black (n=300)
* Sample 7:  65% white (n=650), 35%  black (n=350)
* Sample 8:  60% white (n=600), 40%  black (n=400)
* Sample 9:  55% white (n=550), 45%  black (n=550)
* Sample 10: 50% white (n=500), 50%  black (n=400)

*if n=600
* Sample 0.5:  97.5% white (n=585), 2.5% black (n=15)
* Sample 1:  95% white (n=570), 5%   black (n=30)
* Sample 2:  90% white (n=540), 10%  black (n=60)
* Sample 3:  85% white (n=510), 15%  black (n=90)
* Sample 4:  80% white (n=480), 20%  black (n=120)
* Sample 5:  75% white (n=450), 25%  black (n=150)
* Sample 6:  70% white (n=420), 30%  black (n=180)
* Sample 7:  65% white (n=390), 35%  black (n=210)
* Sample 8:  60% white (n=360), 40%  black (n=240)
* Sample 9:  55% white (n=330), 45%  black (n=270)
* Sample 10: 50% white (n=300), 50%  black (n=300)

local n_white = $p_white*600 
local n_black = $p_black*600 


*clear stored matrices
matrix drop _all


/*** Read in population data set ***/
use population_data_exp, replace


/*** Draw a sample of n=1,000 (based on "real" obs (copy=1) ***/
	
*local n_white = 950
*local n_black = 50
    
randomselect if copy==1 & black==0, gen(selectw) n(`n_white') 				
tab selectw copy
randomselect if copy==1 & black==1, gen(selectb) n(`n_black') 				
tab selectb copy
gen select = 0
replace select = 1 if copy==1 & selectw==1 | selectb==1
tab select copy
	
*store pblack
sum black if select==1
scalar sample_n=r(N)  		 //------------------------------------------------>store sample_n in sample in sim results	
scalar sample_pblack=r(mean) //------------------------------------------------>store sample_pblack in sample in sim results										


keep if select==1
/*Step 1: obtain estimated effect of exposure on amyloid in the sample and store 
results in matrix "observe"*/
	reg amyloid exp black blackexp /*if select==1 */
	predict pamyloid1 
		scalar b0 		= _b[_cons]
		scalar b_exp 	= _b[exp]
		scalar b_black 	= _b[black]
		scalar b_blackexp 	= _b[blackexp]
		matrix observe = (_b[_cons], _b[exp], _b[black], _b[blackexp])
		matrix colnames observe = b0_obs b_exp_obs b_black_obs b_blackexp_obs 		
	/*checking my work
		scalar list
		gen pamyloid1_hand_calc = b0 + exp*b_exp + black*b_black + blackexp*b_blackexp
		gen check = pamyloid1 - pamyloid1_hand_calc
		sum check*/

*To obtain SEs of population effect estimate, draw bootstrap samples from the sample of 1,000 people
	

/*Step 2: Write a program called myboot that samples the data with replacement 
and return the statistics of interest*/
	capture program drop myboot
	program define myboot, rclass 
	 preserve 

	 bsample 600 /*if copy==1*/ //bsample draws bootstrap samples (random samples with replacement) from the data in memory.
	  
		*estimate effect of exposure on amyloid in sample
		reg amyloid exp black blackexp 
		predict pamyloid_boot
			*store beta estimates as scalars
				return scalar b0_boot			= _b[_cons]
				return scalar b_exp_boot 		= _b[exp]
				return scalar b_black_boot 		= _b[black]
				return scalar b_blackexp_boot 	= _b[exp]
			*store beta estimates in a matrix
				matrix this_run = (_b[_cons], _b[exp], _b[black], _b[blackexp])
				matrix cumulative = nullmat(cumulative) \ this_run
				matrix colnames cumulative = b0_boot b_exp_boot b_black_boot b_blackexp_boot 
		/*store pblack in sample
		sum black 
		return scalar pblack=r(mean)*/
		
	 restore
	end


*Step 3
simulate b0_boot=r(b0_boot) b_exp_boot=r(b_exp_boot) ///
b_black_boot=r(b_black_boot) b_blackexp_boot=r(b_blackexp_boot), reps($boot_reps): myboot

*Step 4
bstat, stat(observe) n(1000)
/*bootstrap pblack=r(pblack) b0_boot=r(b0_boot) b_exp_boot=r(b_exp_boot) ///
b_black_boot=r(b_black_boot) b_blackexp_boot=r(b_blackexp_boot), ///
stat(observe) saving(myboot,replace) reps(100): myboot*/

*save bootstrap results in a dta file
save myboot, replace

matrix list observe
matrix list cumulative
matrix list this_run


/***Next, apply estimates from the sample to get the population effect estimates***/

*read in population dataset
use population_data_exp, clear
sort copy


*calculate observed estimated diff in sample using observed beta estimates
gen pamyloid_observe = b0 + exp*b_exp + black*b_black + blackexp*b_blackexp

*use effect estimate in sample to estimate expected value of amyloid in population if exposure set to 1 (copy 2)
sum pamyloid_observe if copy==2
scalar pamyloid_observe_exp_1=r(mean) 

*use effect estimate in sample to estimate expected value of amyloid in population if exposure set to 0 (copy 3)
sum pamyloid_observe if copy==3
scalar pamyloid_observe_exp_0=r(mean)

*estimated effect of exposure on amyloid in population -- estimated as difference in expected values if exp set to 1 vs. exp set to 0
scalar diff_observe = pamyloid_observe_exp_1 - pamyloid_observe_exp_0  //------------------------------------------------>store diff est observed in sample in sim results
	
scalar list diff_observe
	



*Apply estimates from each bootstrapped sample to obtain estimated SE of the population estimate
forvalues i= 1/$boot_reps {
	scalar b0_boot`i' = cumulative[`i',1]
	scalar b_exp_boot`i' = cumulative[`i',2]
	scalar b_black_boot`i' = cumulative[`i',3]
	scalar b_blackexp_boot`i' = cumulative[`i',4]

	*predicted amyloid for each observation in the dataset based on beta estimates
	gen pamyloid_`i' = b0_boot`i' + exp*b_exp_boot`i' + black*b_black_boot`i' + blackexp*b_blackexp_boot`i' 

	*use effect estimate in sample to estimate expected value of amyloid in population if exposure set to 1 (copy 2)
	quietly sum pamyloid_`i' if copy==2 //copy 2 = outcome if exposure set to 1
	scalar pamyloid_exp_1_`i'=r(mean) 

	*use effect estimate in sample to estimate expected value of amyloid in population if exposure set to 0 (copy 3)
	quietly sum pamyloid_`i' if copy==3 //copy 3 = outcome if exposure set to 0
	scalar pamyloid_exp_0_`i'=r(mean)

	*estimated effect of exposure on amyloid in population -- estimated as difference in expected values if exposure set to 1 vs. exposure set to 0
	scalar diff_`i' = pamyloid_exp_1_`i' - pamyloid_exp_0_`i'

	*store estimated effect of exposure on amyloid in population from each bootstrap sample in a matrix
	matrix diff_this_run = (diff_`i')
	matrix diff_cumulative = nullmat(diff_cumulative) \ diff_this_run
	matrix colnames diff_cumulative = diff 	
}

/*Look at results stored in matrices
matrix list diff_this_run
matrix list diff_cumulative
matrix list diff_this_run
*/


*create dataset with estimated effect of exposure on amyloid in population from each bootstrap sample in a matrix
clear 
set obs $boot_reps
gen boot_rep = _n
svmat diff_cumulative, name(diff)

sum diff
scalar boot_mean_diff = r(mean) //------------------------------------------------>store boot_mean_diff in sample in sim results	
scalar boot_SE_diff = r(sd)		//------------------------------------------------>store boot_SE_diff in sample in sim results	

*generate indicator variable for whether the 95% CI includes the null
scalar boot_diff_95CI_lb = diff_observe - 1.96*boot_SE_diff //------------------------------------------------>store boot_diff_95CI_lb in sample in sim results
scalar boot_diff_95CI_ub = diff_observe + 1.96*boot_SE_diff //------------------------------------------------>store boot_diff_95CI_ub in sample in sim results
scalar boot_diff_CI_inc_null = (boot_diff_95CI_lb < 0 & boot_diff_95CI_ub  > 0) //------------------------------------------------>store boot_diff_CI_inc_null in sample in sim results
