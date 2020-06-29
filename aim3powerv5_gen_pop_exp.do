/*Generate population data set --reflects US population of Black & white adults ages 60+
Includes 3 copies of each person: 
	1. copy 1 = "real" data
	2. copy 2 = exposure set to 1
	3. copy 3 = exposure set to 0 
*/


** set this up in a loop to rerun to get the variance estimates for the esimated effect of e4
** does race have an effect on amyloid?
clear all
set obs 500000 //2010 US Census pop Black & white adults ages 60+ n=50061267, but too computationally intensive to samle from 50M
gen id=_n

set seed 22346

*generate race as black as (0,1)
gen black=(id<=50000) // 10% of the pop. is black, 90% is white 	
*sum black				

*generate exposure as (0,1)
gen exp=runiform()<.403 if black==1  // prev htn (defined 140/90) 40.3% among Black adults
replace exp=runiform()<.278 if black==0  // prev htn = 27.8% among whites
/*Prevalnece of htn (140/20 mmHg) is from NHANES https://www.cdc.gov/nchs/data/databriefs/db289.pdf*/

*make three copies of the dataset 
expand 3
sort id
by id: gen copy=_n
* true population effect of exposure for blacks is 0.3								
* true population effect of exposure for whites is 0.2 (average of runiform())  	
* true population effect of exposure is 0.10*0.25 + 0.90*0.15 = 0.16

*generate amyloid as a function of exp and race (race has no main effect, effect is greater in Black vs. white people) 
gen amyloid=exp*0.15+black*exp*0.1+rnormal() if copy==1 //copy 1 = "real" data

*set exp to 1 in copy=2 and set exp to 0 in copy=3
replace exp=1 if copy==2 //copy 2 = outcome if exposure set to 1
replace exp=0 if copy==3 //copy 3 = outcome if exposure set to 0
sum amyloid
*histogram amyloid

*check that estimated effect of exposure on amyloid in "real" data aligns with data-generating rules (amyloid is only defined for copy 1)
gen blackexp=black*exp
reg amyloid exp black blackexp

*check work
reg amyloid i.black#i.exp if copy==1
margins exp


/*** Save population data set ***/
save population_data_exp, replace

/*
. *check work
. reg amyloid i.black#i.exp if copy==1

      Source |       SS           df       MS      Number of obs   =   500,000
-------------+----------------------------------   F(3, 499996)    =    924.80
       Model |  2773.65255         3  924.550851   Prob > F        =    0.0000
    Residual |  499860.697   499,996  .999729392   R-squared       =    0.0055
-------------+----------------------------------   Adj R-squared   =    0.0055
       Total |   502634.35   499,999  1.00527071   Root MSE        =    .99986

------------------------------------------------------------------------------
     amyloid |      Coef.   Std. Err.      t    P>|t|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
   black#exp |
        0 1  |   .1475538    .003331    44.30   0.000     .1410252    .1540825
        1 0  |   .0028461   .0060492     0.47   0.638    -.0090102    .0147024
        1 1  |   .2381153   .0072542    32.82   0.000     .2238974    .2523332
             |
       _cons |   .0038468   .0017528     2.19   0.028     .0004113    .0072822
------------------------------------------------------------------------------

. margins exp

Predictive margins                              Number of obs     =    500,000
Model VCE    : OLS

Expression   : Linear prediction, predict()

------------------------------------------------------------------------------
             |            Delta-method
             |     Margin   Std. Err.      t    P>|t|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
         exp |
          0  |   .0041314   .0016804     2.46   0.014     .0008378     .007425
          1  |   .1604568   .0026447    60.67   0.000     .1552733    .1656403
------------------------------------------------------------------------------
*/



