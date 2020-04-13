********************************************************************************************
*script to generate all FBA data for the 15 trials which will be used by the ED model which
*will be the input to makeKopticInputs.pl which will finally be the KOPTIC inputs
*By: Wheaton Schroeder
*Research Group: PI: Rajib Saha SSBIO University of Nebraska-Lincoln
*Latest Version: 11/21/2019
********************************************************************************************

$INLINECOM /*  */
$EOLCOM !!
$onlisting
$offdigit
$onempty

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
	
	epsilon		/*a very small number*/
	
;

epsilon = 1E-7;
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
*               carbon limited, sucrose as only carbon source, LB(R003ex) = -0.25                   *
*****************************************************************************************************
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

LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=-0.25;	!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

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
c('biomass_si') = 1;

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
	primalobj_ED			primal problem objective function (rxn fluxes in each loop) for ED
	massbalance_ED(iED)  	Mass balance for ED
	dual_obj				dual problem objective 
	dual_const_1(jED)			dual constraint 1
	dual_const_2(jED)			dual constraint 2
	strong_duality				enforces strong duality
;

****************************** Primal Equations*******************************
primalobj_ED..		  zprimal_ED =e= sum(jED$(C(jED) eq 1),vED(jED));
massbalance_ED(iED)..     sum(jED,SED(iED,jED)*vED(jED)) =e= 0;
******************************************************************************

****************************** Dual Equations*********************************
dual_obj..		  	zdual_ED =e= sum(iED, 0 * lambda(iED)) + sum(jED, - LBED(jED) * muLB(jED)) + sum(jED, UBED(jED) * muUB(jED));
dual_const_1(jED)$(c(jED) = 0)..		sum(iED, SED(iED,jED) * lambda(iED)) - muLB(jED) + muUB(jED) =e= 0;
dual_const_2(jED)$(c(jED) = 1)..     	sum(iED, SED(iED,jED) * lambda(iED)) - muLB(jED) + muUB(jED) =e= 1;
******************************************************************************

strong_duality..		zdual_ED =e= zprimal_ED;

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

model STRONGED
/
	primalobj_ED
	massbalance_ED
	dual_obj
	dual_const_1
	dual_const_2
	strong_duality
/
;

*no optimization file
PRIMALED.optfile=1;
DUALED.optfile=1;
STRONGED.optfile=1;

*solve this scenario for maximum biomass
SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

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

*file for shadow prices
FILE SHADOW /shadow_price.csv/;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*write the header
PUT "Trial number","trial code";
PUT "growth rate";
PUT "X00001[c]";
PUT "C04033[c]";
PUT "C00779[c]";
PUT "C01173[c]";
PUT "C01624[c]";
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
PUT "C08607[c]";
PUT "C00097[c]";
PUT "Model Status"/;

*print the shadow price of just the pigments
PUT trial_count,"SCL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

*file for shadow prices related to malonyl-CoA
FILE SHADOWMCoA /shadow_price_MCoA.csv/;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;

*write the header
PUT "Trial number","trial code";
PUT "growth rate";
PUT "C00083[c]";
PUT "C00002[c]";
PUT "C00024[c]";
PUT "C00033[c]";
PUT "C00010[c]";
PUT "C00882[c]";
PUT "C01134[c]";
PUT "C04352[c]";
PUT "C00097[c]";
PUT "C03492[c]";
PUT "C00288[c]";
PUT "C00008[c]";
PUT "C00020[c]";
PUT "C00009[c]";
PUT "C00013[c]";
PUT "C00001[c]";
PUT "C00011[c]";
PUT "C00063[c]";
PUT "C00055[c]";
PUT "Model Status"/;

PUT trial_count,"SCL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

*file for shadow prices related to malonyl-CoA
FILE SHADOWBIO /shadow_price_biomass.csv/;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;

*make header
PUT "C00049[c]";
PUT "C00025[c]";
PUT "C00041[c]";
PUT "C00037[c]";
PUT "C00064[c]";
PUT "C00152[c]";
PUT "C00148[c]";
PUT "C00097[c]";
PUT "C00062[c]";
PUT "C00073[c]";
PUT "C00065[c]";
PUT "C00078[c]";
PUT "C00188[c]";
PUT "C00135[c]";
PUT "C00079[c]";
PUT "C00082[c]";
PUT "C00183[c]";
PUT "C00123[c]";
PUT "C00407[c]";
PUT "C00047[c]";
PUT "C00017[c]";
PUT "C00249[im]";
PUT "C01530[im]";
PUT "C06425[im]";
PUT "C00712[im]";
PUT "C01595[im]";
PUT "C01356[im]";
PUT "C00075[c]";
PUT "C00044[c]";
PUT "C00182[c]";
PUT "C00116[c]";
PUT "C00422[c]";
PUT "C00350[c]";
PUT "C02737[c]";
PUT "C00063[c]";
PUT "C01356[c]";
PUT "C00017[c";
PUT "C00286[c]";
PUT "C00131[c]";
PUT "C00458[c]";
PUT "C00459[c]";
PUT "cell[c]";
PUT "C00965[e]";
PUT "C01356[e]";
PUT "C17937DHN[e]";
PUT "C17937Eu[e]";
PUT "C17937Pyo[e]";
PUT "C00017[e]";
PUT "C00140[e]";
PUT "C00461[e]";
PUT "cell_wall[e]";
PUT "C02094[c]";
PUT "C19892[c]";
PUT "C08607[c]";
PUT "C00097[c]";
PUT "carotenoids[c]";
PUT "biomass"/;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 2************************************************
*               carbon limited, sucrose as only carbon source, LB(R003ex) = -0.5                    *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=-0.5;	!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;
SHADOW.pw = 30000;

*print the shadow price of just the pigments and precursors
PUT trial_count,"SCM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"SCM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 3************************************************
*               carbon limited, sucrose as only carbon source, LB(R003ex) = -1                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=-1;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"SCH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"SCH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 4************************************************
*               carbon limited, ethanol as only carbon source, LB(R009ex) = -1.5                    *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=-1.5;	!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!glucose

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ECL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ECL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 5************************************************
*               carbon limited, ethanol as only carbon source, LB(R009ex) = -3                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=-3;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ECM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ECM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 6************************************************
*               carbon limited, ethanol as only carbon source, LB(R009ex) = -6                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=-6;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ECH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ECH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 7************************************************
*               carbon limited, acetate as only carbon source, LB(R010ex) = -1.5                    *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=-1.5;	!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R011ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ACL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ACL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 8************************************************
*               carbon limited, acetate as only carbon source, LB(R010ex) = -3                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=-3;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ACM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ACM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 9************************************************
*               carbon limited, acetate as only carbon source, LB(R010ex) = -6                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=-6;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!glucose

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.lw = 25;
RXNOUT.pc = 4;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ACH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ACH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 10***********************************************
*               carbon limited, glucose as only carbon source, LB(R011ex) = -0.5                    *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=-0.5;	!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"GCL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"GCL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 11***********************************************
*               carbon limited, glucose as only carbon source, LB(R011ex) = -1                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=-1;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"GCM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"GCM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 12***********************************************
*               carbon limited, glucose as only carbon source, LB(R011ex) = -2                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=-2;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"GCH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"GCH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 13***********************************************
*                  sucrose as carbon source, nitrogen limited, LB(R004ex) = -0.05                   *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=-Vmax/12;	!!sucrose
LBED('R004ex')=-0.05;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;	!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"SNL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"SNL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 14***********************************************
*                  sucrose as carbon source, nitrogen limited, LB(R004ex) = -0.1                    *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=-Vmax/12;	!!sucrose
LBED('R004ex')=-0.1;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"SNM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"SNM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 15***********************************************
*                 sucrose as carbon source, nitrogen limited, LB(R004ex) = -0.2                     *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=-Vmax/12;	!!sucrose
LBED('R004ex')=-0.2;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"SNH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"SNH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 16***********************************************
*                  ethanol as carbon source, nitrogen limited, LB(R004ex) = -0.05                   *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-0.05;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=-Vmax/2;	!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ENL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ENL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 17***********************************************
*                  ethanol as carbon source, nitrogen limited, LB(R004ex) = -0.1                    *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-0.1;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=-Vmax/2;	!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;	!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ENM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ENM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 18***********************************************
*                 ethanol as carbon source, nitrogen limited, LB(R004ex) = -0.2                     *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-0.2;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=-Vmax/2;	!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate
vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ENH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"EHN";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 19***********************************************
*                  acetate as carbon source, nitrogen limited, LB(R004ex) = -0.05                   *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-0.05;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=-Vmax/2;	!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ANL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ANL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 20***********************************************
*                  acetate as carbon source, nitrogen limited, LB(R004ex) = -0.1                    *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-0.1;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=-Vmax/2;	!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ANM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ANM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 21***********************************************
*                 acetate as carbon source, nitrogen limited, LB(R004ex) = -0.2                     *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-0.2;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=-Vmax/2;	!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ANH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ANH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 22***********************************************
*                 glucose as carbon source, nitrogen limited, LB(R004ex) = -0.05                    *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-0.05;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=-Vmax/6;	!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"GNL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"GNL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 23***********************************************
*                 glucose as carbon source, nitrogen limited, LB(R004ex) = -0.1                     *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-0.1;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=-Vmax/6;	!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"GNM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"GNM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 24***********************************************
*                glucose as carbon source, nitrogen limited, LB(R004ex) = -0.2                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-0.2;	!!ammonia
LBED('R005ex')=-Vmax;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=-Vmax/6;	!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"GNH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"GNH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 25***********************************************
*               sucrose as carbon source, sulfate limited, LB(R005ex) = -0.0005                     *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=-Vmax/12;	!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-0.005;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"SSL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"SSL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 26***********************************************
*               sucrose as carbon source, sulfate limited, LB(R005ex) = -0.001                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=-Vmax/12;	!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-0.001;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate
vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"SSM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"SSM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 27***********************************************
*               sucrose as carbon source, sulfate limited, LB(R005ex) = -0.002                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=-Vmax/12;	!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-0.002;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"SSH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"SSH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 28***********************************************
*              ethanol as carbon source, sulfate limited, LB(R005ex) = -0.0005                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-0.0005;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=-Vmax/2;	!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;	!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ESL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ESL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 29***********************************************
*              ethanol as carbon source, sulfate limited, LB(R005ex) = -0.001                       *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-0.001;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=-Vmax/2;	!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ESM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ESM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 30***********************************************
*              ethanol as carbon source, sulfate limited, LB(R005ex) = -0.002                       *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-0.002;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=-Vmax/2;	!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ESH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ESH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 31***********************************************
*              acetate as carbon source, sulfate limited, LB(R005ex) = -0.0005                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-0.0005;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=-Vmax/2;	!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ASL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ASL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 32***********************************************
*              acetate as carbon source, sulfate limited, LB(R005ex) = -0.001                       *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-0.001;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=-Vmax/2;	!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ASM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ASM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 33***********************************************
*              acetate as carbon source, sulfate limited, LB(R005ex) = -0.002                       *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-0.002;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=-Vmax/2;	!!acetate
LBED('R011ex')=0;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"ASH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;

PUT trial_count,"ASH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 34***********************************************
*              glucose as carbon source, sulfate limited, LB(R005ex) = -0.0005                      *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-0.0005;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=-Vmax/6;	!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"GSL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;

PUT trial_count,"GSL";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 35***********************************************
*              glucose as carbon source, sulfate limited, LB(R005ex) = -0.001                       *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-0.001;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=-Vmax/6;	!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"GSM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;

PUT trial_count,"GSM";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;

**********************************************trial 36***********************************************
*              glucose as carbon source, sulfate limited, LB(R005ex) = -0.002                       *
*****************************************************************************************************
LBED('R001ex')=-Vmax;	!!water
LBED('R002ex')=-Vmax;	!!phosphate
LBED('R003ex')=0;		!!sucrose
LBED('R004ex')=-Vmax;	!!ammonia
LBED('R005ex')=-0.002;	!!sulfate
LBED('R006ex')=-Vmax;	!!H+
LBED('R007ex')=0;		!!Carbon dioxide
LBED('R008ex')=-Vmax;	!!oxygen
LBED('R009ex')=0;		!!ethanol
LBED('R010ex')=0;		!!acetate
LBED('R011ex')=-Vmax/6;	!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;	!!water
UBED('R002ex')=0;		!!phosphate
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;		!!ammonia
UBED('R005ex')=0;		!!sulfate
UBED('R006ex')=Vmax;	!!H+
UBED('R007ex')=Vmax;	!!Carbon dioxide
UBED('R008ex')=0;		!!oxygen
UBED('R009ex')=Vmax;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);

SOLVE STRONGED USING LP MAXIMIZING zprimal_ED;

RXNOUT.ap = 1;
PUT RXNOUT;
RXNOUT.pc = 4;
RXNOUT.lw = 25;

LOOP(jED,

	PUT trial_count,jED.tl,vED.l(jED):0:8/;

);

PUT "/";

PUTCLOSE;

SHADOW.ap = 1;
PUT SHADOW;
SHADOW.pc = 5;
SHADOW.lw = 25;
SHADOW.pw = 30000;

*print the shadow price of just the pigments
PUT trial_count,"GSH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('X00001[c]'):0:8;
PUT lambda.l('C04033[c]'):0:8;
PUT lambda.l('C00779[c]'):0:8;
PUT lambda.l('C01173[c]'):0:8;
PUT lambda.l('C01624[c]'):0:8;
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
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWMCoA.ap = 1;
PUT SHADOWMCoA;
SHADOWMCoA.pc = 5;
SHADOWMCoA.lw = 25;
SHADOWMCoA.pw = 30000;

PUT trial_count,"GSH";
PUT vED.l('biomass_wt_37'):0:8;
PUT lambda.l('C00083[c]'):0:8;
PUT lambda.l('C00002[c]'):0:8;
PUT lambda.l('C00024[c]'):0:8;
PUT lambda.l('C00033[c]'):0:8;
PUT lambda.l('C00010[c]'):0:8;
PUT lambda.l('C00882[c]'):0:8;
PUT lambda.l('C01134[c]'):0:8;
PUT lambda.l('C04352[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C03492[c]'):0:8;
PUT lambda.l('C00288[c]'):0:8;
PUT lambda.l('C00008[c]'):0:8;
PUT lambda.l('C00020[c]'):0:8;
PUT lambda.l('C00009[c]'):0:8;
PUT lambda.l('C00013[c]'):0:8;
PUT lambda.l('C00001[c]'):0:8;
PUT lambda.l('C00011[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C00055[c]'):0:8;
PUT STRONGED.ModelStat/;
PUTCLOSE;

SHADOWBIO.ap = 1;
PUT SHADOWBIO;
SHADOWBIO.pc = 5;
SHADOWBIO.lw = 25;
SHADOWBIO.pw = 30000;
SHADOWBIO.pw = 30000;

PUT lambda.l('C00049[c]'):0:8;
PUT lambda.l('C00025[c]'):0:8;
PUT lambda.l('C00041[c]'):0:8;
PUT lambda.l('C00037[c]'):0:8;
PUT lambda.l('C00064[c]'):0:8;
PUT lambda.l('C00152[c]'):0:8;
PUT lambda.l('C00148[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('C00062[c]'):0:8;
PUT lambda.l('C00073[c]'):0:8;
PUT lambda.l('C00065[c]'):0:8;
PUT lambda.l('C00078[c]'):0:8;
PUT lambda.l('C00188[c]'):0:8;
PUT lambda.l('C00135[c]'):0:8;
PUT lambda.l('C00079[c]'):0:8;
PUT lambda.l('C00082[c]'):0:8;
PUT lambda.l('C00183[c]'):0:8;
PUT lambda.l('C00123[c]'):0:8;
PUT lambda.l('C00407[c]'):0:8;
PUT lambda.l('C00047[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00249[im]'):0:8;
PUT lambda.l('C01530[im]'):0:8;
PUT lambda.l('C06425[im]'):0:8;
PUT lambda.l('C00712[im]'):0:8;
PUT lambda.l('C01595[im]'):0:8;
PUT lambda.l('C01356[im]'):0:8;
PUT lambda.l('C00075[c]'):0:8;
PUT lambda.l('C00044[c]'):0:8;
PUT lambda.l('C00182[c]'):0:8;
PUT lambda.l('C00116[c]'):0:8;
PUT lambda.l('C00422[c]'):0:8;
PUT lambda.l('C00350[c]'):0:8;
PUT lambda.l('C02737[c]'):0:8;
PUT lambda.l('C00063[c]'):0:8;
PUT lambda.l('C01356[c]'):0:8;
PUT lambda.l('C00017[c]'):0:8;
PUT lambda.l('C00286[c]'):0:8;
PUT lambda.l('C00131[c]'):0:8;
PUT lambda.l('C00458[c]'):0:8;
PUT lambda.l('C00459[c]'):0:8;
PUT lambda.l('cell[c]'):0:8;
PUT lambda.l('C00965[e]'):0:8;
PUT lambda.l('C01356[e]'):0:8;
PUT lambda.l('C17937DHN[e]'):0:8;
PUT lambda.l('C17937Eu[e]'):0:8;
PUT lambda.l('C17937Pyo[e]'):0:8;
PUT lambda.l('C00017[e]'):0:8;
PUT lambda.l('C00140[e]'):0:8;
PUT lambda.l('C00461[e]'):0:8;
PUT lambda.l('cell_wall[e]'):0:8;
PUT lambda.l('C02094[c]'):0:8;
PUT lambda.l('C19892[c]'):0:8;
PUT lambda.l('C08607[c]'):0:8;
PUT lambda.l('C00097[c]'):0:8;
PUT lambda.l('carotenoids[c]'):0:8;
PUT lambda.l('biomass'):0:8/;

PUTCLOSE;

*update the trial counter
trial_count_temp = trial_count + 1;
trial_count = trial_count_temp;
