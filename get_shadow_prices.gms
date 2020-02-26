********************************************************************************************
*script to generate all FBA data for the 15 trials which will be used by the ED model which
*will be the input to makeKopticInputs.pl which will finally be the KOPTIC inputs
*By: Wheaton Schroeder
*Research Group: PI: Rajib Saha SSBIO University of Nebraska-Lincoln
*Latest Version: 11/21/2019
********************************************************************************************

$INLINECOM /*  */
$onlisting
$offdigit

options limrow = 1000
	optCR = 1E-9
	optCA = 1E-9
	iterlim = 1000000
	decimals = 8
	reslim = 1000000
	work = 50000000
    lp = cplex;
        

SETS
        jED   Reactions of k-ath model
$include "rxns.txt"

        iED   Metabolites of k-ath model
$include "mets.txt"

	reg_off(jED)	Reactions turned off (regulated)
$include "reg_rxns.txt"

	irrev_f(jED)	Reactions that are irreversible, forwards
$include "irrev_rxns_f.txt"

	irrev_b(jED)	reactions that are irreversible, backwards
$include "irrev_rxns_b.txt"
;

PARAMETERS
  UBED(jED) UB on rxns of ED 
  LBED(jED) LB on rxns of ED


  rxntype_ED(jED) Reaction type for rxns of ED 
$include "rxntype.txt"

	trial_count
	
	trial_count_temp
	
	conc(iED)			/*metabolite concentration*/
	
	density			/*density of E. dermatitidis*/
	
	epsilon		/*a very small number*/
	
;

epsilon = 1E-7;
density = 1;
trial_count = 1;
trial_count_temp = 1;

******************  Set the values of LB and UB ******************
SCALAR Vmax /10000/;
SCALAR BigM /9900/;

*** Set the experimental conditions

* irreversible reactions forward (rxntype = 1)
UBED(jED)$(rxntype_ED(jED) = 1) = Vmax;
LBED(jED)$(rxntype_ED(jED) = 1) = 0;

* For reversible reactions (rxntype = 0) 
UBED(jED)$(rxntype_ED(jED) = 0) = Vmax;
LBED(jED)$(rxntype_ED(jED) = 0) = -Vmax;

* irreversible reactions backwards (rxntype = -1)
UBED(jED)$(rxntype_ED(jED) = -1) = 0;
LBED(jED)$(rxntype_ED(jED) = -1) = -Vmax;


**********************************************trial 1************************************************
*               carbon limited, sucrose as only carbon source, LB(R350ex) = -1                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=-1; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
*uptake of a necessary resource
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

* ATPM
LBED('ATPM') = eps;
UBED('ATPM') = 30;

* Reactions turned off 
UBED(jED)$reg_off(jED) = 0;
LBED(jED)$reg_off(jED) = 0;

****************stoichiometry*****************
PARAMETERS
  SED(iED,jED) Stoichiometric matrix for ED 
$include "Sij.txt"
;

PARAMETERS
C(jED)
max(jED)
min(jED);
alias (j1ED, jED);

*set objective as the biomass function
c(jED) = 0;
*c('R00235[c]') = 1;
c('biomass_wt_37') = 1;

VARIABLES
	zprimal_ED     primal objective function (rxn fluxes in each loop) for ED
	zdual_ED     primal objective function (rxn fluxes in each loop) for ED
    vED(jED)       Flux of ED rxns
	lambda(iED)		dual variable associated with the mass balance
	
POSITIVE VARIABLES

	muLB(jED)		dual variable associated with the lower bound constraint
	muUB(jED)		dual variable associated with the upper bound constraint
	
;
	

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);


*********************EQUATIONS NAMES**************************
EQUATIONS
	primalobj_ED			primal objective function (rxn fluxes in each loop) for ED
	massbalance_ED(iED)  	Mass balance for ED
	dual_obj				dual objectie 
	dual_const_1(jED)			dual constraint 1
	dual_const_2(jED)			dual constraint 2
;

****************************** Primal Equations*******************************
primalobj_ED..		  zprimal_ED =e= sum(jED$(C(jED) eq 1),vED(jED));

massbalance_ED(iED)..     sum(jED,SED(iED,jED)*vED(jED)) =e= 0;
******************************************************************************

****************************** Dual Equations*********************************
dual_obj..		  	zdual_ED =e= 0 * sum(iED, lambda(iED)) + sum(jED, -LBED(jED) * muLB(jED)) + sum(jED, UBED(jED) * muUB(jED));
dual_const_1(jED)$(c(jED) = 0)..		sum(iED, SED(iED,jED) * lambda(iED)) - muLB(jED) + muUB(jED) =e= 0;
dual_const_2(jED)$(c(jED) = 1)..     	sum(iED, SED(iED,jED) * lambda(iED)) - muLB(jED) + muUB(jED) =e= 1;
******************************************************************************

model PRIMALED
/
	primalobj_ED
	massbalance_ED
/
;

model DUALED
/
	dual_obj
	dual_const_1
	dual_const_2
/
;

*no optimization file
PRIMALED.optfile=0;

*solve this scenario for maximum biomass
SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

*print output for the flux rates
FILE RXNOUT /rxn_rates_out.txt/;
PUT RXNOUT;
RXNOUT.lw = 25;
RXNOUT.pc = 4;

PUT "/"/;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*print output for the concentration
FILE METOUT /met_conc_out.txt/;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

PUT "/"/;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

*file for shadow prices
FILE SHADOW /shadow_price.csv/;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*write the header
PUT "Trial number","trial code";
PUT "growth rate";
PUT "X00001[c]";
PUT "C04033[c]";
PUT "C00083[c]";
PUT "C17937DHN[e]";
PUT "C05578[e]";
PUT "C05579[e]";
PUT "C01693[e]";
PUT "C05604[e]";
PUT "C00355[c]";
PUT "C00822[c]";
PUT "C00082[c]";
PUT "C17937Eu[e]";
PUT "C00544[c]";
PUT "X00002[c]";
PUT "C01179[c]";
PUT "C17937Pyo[e]";
PUT "C00353[c]";
PUT "C03427[c]";
PUT "C05421[c]";
PUT "C05414[c]";
PUT "C05430[c]";
PUT "C05431[c]";
PUT "C05432[c]";
PUT "C05435[c]";
PUT "C15867[c]";
PUT "C08613[c]";
PUT "C02094[c]";
PUT "C19892[c]";
PUT "C08607[c]"/;

*print the shadow price of just the pigments
PUT trial_count;
PUT "SCL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 2************************************************
*               carbon limited, sucrose as only carbon source, LB(R350ex) = -5                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=-5; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=-0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-0;
*glucose
LBED('R706ex')=-0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments and precursors
PUT trial_count,"SCM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 3************************************************
*               carbon limited, sucrose as only carbon source, LB(R350ex) = -10                     *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=-10; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-0;
*glucose
LBED('R706ex')=-0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"SCH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 4************************************************
*               carbon limited, ethanol as only carbon source, LB(R356ex) = -1                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=-1;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-0;
*glucose
LBED('R706ex')=-0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ECL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 5************************************************
*               carbon limited, ethanol as only carbon source, LB(R356ex) = -5                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=-5;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-0;
*glucose
LBED('R706ex')=-0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ECM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 6************************************************
*               carbon limited, ethanol as only carbon source, LB(R356ex) = -10                     *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=-10;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-0;
*glucose
LBED('R706ex')=-0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ECH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 7************************************************
*               carbon limited, acetate as only carbon source, LB(R705ex) = -1                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-1;
*glucose
LBED('R706ex')=-0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ACL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 8************************************************
*               carbon limited, acetate as only carbon source, LB(R705ex) = -5                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-5;
*glucose
LBED('R706ex')=-0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ACM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 9************************************************
*               carbon limited, acetate as only carbon source, LB(R705ex) = -10                     *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-10;
*glucose
LBED('R706ex')=-0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.lw = 25;
RXNOUT.pc = 4;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ACH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 10***********************************************
*               carbon limited, glucose as only carbon source, LB(R706ex) = -1                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-0;
*glucose
LBED('R706ex')=-1;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"GCL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 11***********************************************
*               carbon limited, glucose as only carbon source, LB(R706ex) = -5                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-0;
*glucose
LBED('R706ex')=-5;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"GCM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 12***********************************************
*               carbon limited, glucose as only carbon source, LB(R706ex) = -10                     *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-0;
*glucose
LBED('R706ex')=-10;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
MEtOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"GCH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 13***********************************************
*                  sucrose as carbon source, nitrogen limited, LB(R351ex) = -1                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=-1000; 
*ammonium
LBED('R351ex')=-1;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"SNL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 14***********************************************
*                  sucrose as carbon source, nitrogen limited, LB(R351ex) = -5                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=-1000; 
*ammonium
LBED('R351ex')=-5;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"SNM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 15***********************************************
*                 sucrose as carbon source, nitrogen limited, LB(R351ex) = -20                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=-1000; 
*ammonium
LBED('R351ex')=-20;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-0;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"SNH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 16***********************************************
*                  ethanol as carbon source, nitrogen limited, LB(R351ex) = -1                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=-1000;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-0;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ENL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 17***********************************************
*                  ethanol as carbon source, nitrogen limited, LB(R351ex) = -5                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-5;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=-1000;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ENM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 18***********************************************
*                 ethanol as carbon source, nitrogen limited, LB(R351ex) = -20                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-20;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=-1000;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ENH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 19***********************************************
*                  acetate as carbon source, nitrogen limited, LB(R351ex) = -1                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-1000;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ANL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 20***********************************************
*                  acetate as carbon source, nitrogen limited, LB(R351ex) = -5                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-5;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-1000;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ANM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 21***********************************************
*                 acetate as carbon source, nitrogen limited, LB(R351ex) = -20                      *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-20;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=-0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-1000;
*glucose
LBED('R706ex')=-0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ANH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 22***********************************************
*                 glucose as carbon source, nitrogen limited, LB(R351ex) = -1                       *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=-1000;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"GNL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 23***********************************************
*                 glucose as carbon source, nitrogen limited, LB(R351ex) = -5                       *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-5;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=-1000;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"GNM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 24***********************************************
*                glucose as carbon source, nitrogen limited, LB(R351ex) = -20                       *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-20;
*sulfate
LBED('R352ex')=-1000;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=-1000;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"GNH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 25***********************************************
*               sucrose as carbon source, sulfate limited, LB(R352ex) = -0.01                       *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=-1000; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-0.01;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"SSL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 26***********************************************
*               sucrose as carbon source, sulfate limited, LB(R352ex) = -0.1                        *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=-1000; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-0.1;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"SSM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 27***********************************************
*               sucrose as carbon source, sulfate limited, LB(R352ex) = -0.2                        *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=-1000; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-0.2;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"SSH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 28***********************************************
*              ethanol as carbon source, sulfate limited, LB(R352ex) = -0.01                        *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-0.01;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=-1000;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ESL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 29***********************************************
*              ethanol as carbon source, sulfate limited, LB(R352ex) = -0.1                         *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-0.1;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=-1000;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=0;
LBED('R707ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;
UBED('R707ex')=Vmax;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ESM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 30***********************************************
*              ethanol as carbon source, sulfate limited, LB(R352ex) = -0.2                         *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-0.2;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=-1000;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=0;
*glucose
LBED('R706ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ESH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 31***********************************************
*              acetate as carbon source, sulfate limited, LB(R352ex) = -0.01                        *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-0.01;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-1000;
*glucose
LBED('R706ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ASL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 32***********************************************
*              acetate as carbon source, sulfate limited, LB(R352ex) = -0.1                         *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-0.1;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-1000;
*glucose
LBED('R706ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ASM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 33***********************************************
*              acetate as carbon source, sulfate limited, LB(R352ex) = -0.2                         *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-0.2;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-1000;
*glucose
LBED('R706ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"ASH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 34***********************************************
*              acetate as carbon source, sulfate limited, LB(R352ex) = -0.01                        *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-0.01;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-1000;
*glucose
LBED('R706ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"GSL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 35***********************************************
*              acetate as carbon source, sulfate limited, LB(R352ex) = -0.1                         *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-0.1;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-1000;
*glucose
LBED('R706ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"GSM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 36***********************************************
*              acetate as carbon source, sulfate limited, LB(R352ex) = -0.2                         *
*****************************************************************************************************
*water
LBED('R348ex')=-1000; 
*phosphate 
LBED('R349ex')=-1000; 
*sucrose
LBED('R350ex')=0; 
*ammonium
LBED('R351ex')=-1000;
*sulfate
LBED('R352ex')=-0.2;
*protons
LBED('R353ex')=-1000;
*CO2
LBED('R354ex')=0;
*oxygen
LBED('R355ex')=-1000;
*ethanol
LBED('R356ex')=0;
LBED('R704ex')=-1000;
*acetate
LBED('R705ex')=-1000;
*glucose
LBED('R706ex')=0;

*set upper bounds for these exchanges
UBED('R348ex')=Vmax;
UBED('R349ex')=0;
UBED('R350ex')=0;
UBED('R351ex')=0;
UBED('R352ex')=0;
UBED('R353ex')=0;
UBED('R354ex')=Vmax;
UBED('R355ex')=0;
UBED('R356ex')=Vmax;
UBED('R704ex')=0;
UBED('R705ex')=0;
UBED('R706ex')=0;

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE PRIMALED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUT "/";

PUTCLOSE;

*add this trials results to the met output
METOUT.ap = 1;
PUT METOUT;
METOUT.pc = 4;
METOUT.lw = 25;

LOOP(iED,

	conc(iED) = 0.5 * sum(jED,abs(SED(iED,jED)*vED.l(jED))) * (1 / (vED.l('biomass_si') + epsilon)) * density * 1000;
	PUT trial_count,iED.tl,conc(iED)/;

);

PUT "/";

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*Solve the dual problem to calculate shadow price
SOLVE DUALED USING LP MINIMIZING zdual_ED;

*print the shadow price of just the pigments
PUT trial_count,"GSH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C05578[e]'):0:8;
PUT lambda.l('C05579[e]'):0:8;
PUT lambda.l('C01693[e]'):0:8;
PUT lambda.l('C05604[e]'):0:8;
PUT lambda.l('C00355[c]'):0:8;
PUT lambda.l('C00822[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C00544[c]'):0:8;
PUT lambda.l('X00002[c]'):0:8;
PUT lambda.l('C01179[c]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00353[c]'):0:8;
PUT lambda.l('C03427[c]'):0:8;
PUT lambda.l('C05421[c]'):0:8;
PUT lambda.l('C05414[c]'):0:8;
PUT lambda.l('C05430[c]'):0:8;
PUT lambda.l('C05431[c]'):0:8;
PUT lambda.l('C05432[c]'):0:8;
PUT lambda.l('C05435[c]'):0:8;
PUT lambda.l('C15867[c]'):0:8;
PUT lambda.l('C08613[c]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8/;
PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;
