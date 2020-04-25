*GAMS code for Flux Variability Analysis of Arabidopsis Thaliana seed model
*Original code by: Mohammad Mazharul Islam
*modified for k-ath by: Wheaton Schroeder

$INLINECOM /*  */
$EOLCOM !!
$ONLISTING
$ONEMPTY
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

LBED('R001ex')=-Vmax;
LBED('R002ex')=-Vmax;
LBED('R003ex')=-1;		!!sucrose
LBED('R004ex')=-Vmax;
LBED('R005ex')=-Vmax;
LBED('R006ex')=-Vmax;
LBED('R007ex')=0;
LBED('R008ex')=-Vmax;
LBED('R009ex')=-1;		!!ethanol
LBED('R010ex')=-1;		!!acetate
LBED('R011ex')=-1;		!!glucose
LBED('R012ex')=0;		!!urate

UBED('R001ex')=Vmax;
UBED('R002ex')=0;
UBED('R003ex')=0;		!!sucrose
UBED('R004ex')=0;
UBED('R005ex')=0;
UBED('R006ex')=Vmax;
UBED('R007ex')=Vmax;
UBED('R008ex')=0;
UBED('R009ex')=0;		!!ethanol
UBED('R010ex')=0;		!!acetate
UBED('R011ex')=0;		!!glucose
UBED('R012ex')=Vmax;	!!urate

*LBED('R708tex')=0;
*UBED('R708tex')=0;

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

VARIABLES
	zprimal_ED     primal objective function (rxn fluxes in each loop) for ED
    vED(jED)       Flux of ED rxns
;
	

vED.lo(jED)=LBED(jED);
vED.up(jED)=UBED(jED);


*********************EQUATIONS NAMES**************************
EQUATIONS
	primalobj_ED			primal objective function (rxn fluxes in each loop) for ED
	massbalance_ED(iED)  	Mass balance for ED
;

****************************** Equations*******************************
primalobj_ED..		  zprimal_ED =e= sum(jED$(C(jED) eq 1),vED(jED));

massbalance_ED(iED)..     sum(jED,SED(iED,jED)*vED(jED)) =e= 0;
*****************************************************************************

file result /FVA_result_ED.txt/;
put result;
result.lw = 20;
put 'FVA results for ED'/;
put 'Reaction','                 MIN','             MAX','             LB','          UB','       Model Status'/;
put '---------------------------------------------------------------------------------------------------'/;

model primalED
/
	primalobj_ED
	massbalance_ED
/
;

PUTCLOSE;

LOOP(j1ED,
	
	result.ap = 1;
	PUT result;
	
	C(jED) = 0;
	C(j1ED) = 1;
	
	primalED.optfile=1;
	Solve primalED using lp maximizing zprimal_ED;
	max(j1ED) = zprimal_ED.l;

	Solve primalED using lp minimizing zprimal_ED;
	min(j1ED) = zprimal_ED.l;
	put j1ED.tl,min(j1ED):0:8,'      ',max(j1ED):0:8,LBED(j1ED),UBED(j1ED),primalED.modelstat/;

	PUTCLOSE;

);

result.ap = 1;
PUT result;

put //;

put 'Fluxes hitting upper bounds'/;
put '------------------------------------------------------------------------'/;

LOOP(jED,
	if(max(jED)>= 9000, put jED.tl,
		'	Max flux =', max(jED)/;
	);
);

put 'Fluxes hitting lower bounds'/;
put '------------------------------------------------------------------------'/;
LOOP(jED,
	if(min(jED)<= -9000, put jED.tl,
		'	Min flux =', min(jED)/;
	);
);


put //;

put 'Blocked reactions'/;
put '------------------------------------------------------------------------'/;
LOOP(jED,
	if(bool_and(max(jED)= 0,min(jED)= 0), 
		put jED.tl/;
	);
);

PUT //;

put 'Possible problems |v|>=3000'/;
put '------------------------------------------------------------------------'/;
LOOP(jED,
	if(bool_or(max(jED)>=3000,min(jED)<=-3000), 
		put jED.tl,'	Max flux =', max(jED),'	Min flux =', min(jED)/;
	);
);

PUTCLOSE;