
set more off
clear all

di "$S_TIME"

timer clear 1

timer on 1

/*I'm not sure why, but when I try to specify 35% Black, 65% white, I'm only 
sampling 350 Blacks and 649 whites--so the bootstrap won't run since it's
999 people. But the code works fine for i = 8, 9, and 10 */
foreach i of numlist 0.5 1 3 5 6  {
	
	global p_white = 1-(`i'*0.05)
	global p_black = `i'*0.05
	
	global B = 200
	
	include aim3powerv5_run_sim_exp.do
	
	}
	
di "$S_TIME"

 timer off 1
 
 timer list 1

*Specify racial/ethnic composition of the sample 
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
