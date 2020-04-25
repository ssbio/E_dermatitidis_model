***************************************************************************************
*                  			FBA for the k-ath seed model				              *
*					Original Script for E-coli: Mohammad Islam						  *
*				Adaptation for Arabidopsis seed: Wheaton Schroeder					  *
***************************************************************************************

$INLINECOM /*  */
$onlisting

OPTIONS

        limrow = 1000
        optCR = 1E-9
        optCA = 1E-9
        iterlim = 100000
        decimals = 8
        reslim = 100000
        work = 50000000
		lp = cplex
;

SETS
        i *Set of metabolites
$include "mets.txt"

        j *Set of reactions
$include "rxns.txt"

        reversible(j) 		*Reversible rxns (other than exchange rxns)
$include "rev_rxns_no_ex.txt"

        irreversible_f(j)  	*Irreversible rxns proceeding forward
$include "irrev_rxns_f.txt"

		irreversible_b(j)	*Irreversible rxns proceding backwards
$include "irrev_rxns_b.txt"

        regulation(j) 		*Reactions that should off due to regulatory constraints
$include "reg_rxns.txt"

		medium(j)			*list of medium reactions
$include "medium.txt"
;

PARAMETERS
  UB(j) 	*Lowerbound on reaction fluxes

  LB(j) 	*Upperbound on reaction fluxes

  S(i,j) contains the Stoichiometric matrix of the metabolic model
$include "Sij.txt"

  pool		*pool variable for pool size iterations
  
    rxntype(j) Reaction type for rxns of ED 
$include "rxntype.txt"

	c(j)				*cost vector file lists objective reaction(s)

  ;

SCALAR Vmax /10000/;
c(j) = 0;
c('biomass_si') = 0;

***** Set the bounds *****
* irreversible reactions forward (rxntype = 1)
UB(j)$(rxntype(j) = 1) = Vmax;
LB(j)$(rxntype(j) = 1) = 0;

* For reversible reactions (rxntype = 0) 
UB(j)$(rxntype(j) = 0) = Vmax;
LB(j)$(rxntype(j) = 0) = -Vmax;

* irreversible reactions backwards (rxntype = -1)
UB(j)$(rxntype(j) = -1) = 0;
LB(j)$(rxntype(j) = -1) = -Vmax;

* Regulated Reactions
*LB(j)$regulation(j) = 0;
*UB(j)$regulation(j) = 0;

***** Set the conditions for the growth medium *****
LB(j)$medium(j) = -Vmax;
UB(j)$medium(j) = 0;

*set the limited substrate
LB('R350ex')=-10;

UB('R348ex')=Vmax;
UB('R349ex')=0;
UB('R350ex')=0;
UB('R351ex')=0;
UB('R352ex')=0;
UB('R353ex')=0;
UB('R354ex')=Vmax;
UB('R355ex')=0;
UB('R356ex')=Vmax;
UB('R704ex')=0;

VARIABLES

    v(j)      *Flux
    z         *primal objective variable
	lambda(i)	*dual variable associated with metabolites
	
NONNEGATIVE VARIABLES

	mu_LB(j)	*dual variable associated with LB constraints
	mu_UB(j)	*dual variable associated with UB constraints

;

*fix ATP maintenance
LB('ATPM') = 8.39;
UB('ATPM') = 8.39;

v.lo(j)=LB(j);
v.up(j)=UB(j);

EQUATIONS

	massbalance(i)   	/*Mass balance equations for each metabolite i*/
	obj		        	/*objective function garunteed optimality by formulation*/
	dual(j)				/*set of dual constrains associated with reactions*/	

;

obj..   			z =e= sum(j, c(j) * v(j) );
massbalance(i)..  	sum( j, S(i,j)*v(j) ) =e= 0;
dual(j)..			sum(i, lambda(i) * S(i,j)) + mu_UB(j) - mu_LB(j) =e= c(j);


************** Model definitions ********************
Model balance
/
  massbalance
  obj
  dual
/;

* First find out the maximum biomass

SOLVE balance USING LP MAXIMIZING z;

FILE FLUXRESULTS /ED_v_1.txt/;
PUT FLUXRESULTS;
FLUXRESULTS.pc=4;

PUT "/"/;

LOOP(j,
	PUT j.tl,v.l(j):0:8/;
);

PUT "/"/;


FILE CONCENTRATION /ED_pool_1.txt/;
PUT CONCENTRATION;
CONCENTRATION.pc=4;

PUT "/"/;


LOOP(i,
	pool = 0.5*sum(j,abs(S(i,j)*v.l(j)));
	IF(pool eq 0,
		PUT i.tl,1/;
	ELSE
		PUT i.tl,pool:0:8/;
	);
);


PUT "/"/;

FILE RESULTS /fba_results_1.txt/;
PUT RESULTS;
RESULTS.pc=6;

PUT 'Reaction','Rate','LB','UB','Model Status','mu_LB','mu_UB','.m value'/;
PUT '----------------------------------------------------------------------------------------'/;

LOOP(j,

	PUT j.tl,v.l(j),LB(j),UB(j),balance.modelstat,mu_LB.l(j):0:8,mu_UB.l(j):0:8,v.m(j):0:8/;

);

FILE SHADOW /shadow_price.txt/;
PUT SHADOW;
SHADOW.pc=6;

PUT 'Metabolite','lambda'/;
PUT '------------------------------------------------------------------------'/;
LOOP(i,

	PUT i.tl,lambda.l(i):0:8/;

);

