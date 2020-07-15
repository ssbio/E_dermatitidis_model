#!usr/bin/perl

#Written by: Wheaton Schroeder
#Latest Version: 07/14/2020

#Written to substitute KEGG IDs with BIGG IDs where applicable
#Note that BIGG has a way to convert to KEGG, KEGG does not have a way to convert to BIGG
#Therefore we have downloaded the BIGG metabolite and reaction databases and will transform these
#to hashes to convert to KEGG IDs to BIGG IDs to address a reviewer comment

#will not make these replacements in the original model, rather will create a copy with the 
#replacements made

use strict;

#make metabolite hash
my %kegg_to_bigg_mets = ( );

#make reaction hash
my %kegg_to_bigg_rxns = ( );

#read the BIGG metabolites file
open(BIGGMETS, "<bigg_mets.txt") or die "could not open BIGG mets file, reason: $!\n";
chomp(my @bigg_mets = <BIGGMETS>);

#for each line, determine if has a KEGG ID that it can map to, if so make hash entry
#otherwise do nothing
for(my $a = 0; $a <= $#bigg_mets; $a++) {
	
	#values are tab seperated
	my @line = split /\t/, $bigg_mets[$a];
	
	#check in the database_links column (4th column, so index of 3) for possible KEGG ID
	#also make sure the hash enry doesn't already exist, since each compartment gets its own BIGG
	#metabolite line
	if($line[3] =~ /KEGG Compound: http:\/\/identifiers.org\/kegg\.compound\/(C\d\d\d\d\d);/) {
		
		my $kegg_id = $1;
		
		if(not exists $kegg_to_bigg_mets{$kegg_id}) {
			
			#if KEGG ID is there, add hash entry
			#will reference the universal BIGG ID and will keep our own compartment info
			#just a bit different format
			$kegg_to_bigg_mets{$1} = $line[1];
			
		}
		
	}
	
}

#by this point metabolite hash should be done

#now make reaction hash

#read the BIGG metabolites file
open(BIGGRXNS, "<bigg_rxns.txt") or die "could not open BIGG rxns file, reason: $!\n";
chomp(my @bigg_rxns = <BIGGRXNS>);

#for each line, determine if has a KEGG ID that it can map to, if so make hash entry
#otherwise do nothing
for(my $b = 0; $b <= $#bigg_rxns; $b++) {
	
	#values are tab seperated
	my @line = split /\t/, $bigg_rxns[$b];
	
	#check in the database_links column (4th column, so index of 3) for possible KEGG ID
	#also make sure the hash enry doesn't already exist, since each compartment gets its own BIGG
	#metabolite line
	if($line[3] =~ /KEGG Reaction: http:\/\/identifiers.org\/kegg\.reaction\/(R\d\d\d\d\d);/) {
		
		my $kegg_id = $1;
		
		if(not exists $kegg_to_bigg_rxns{$kegg_id}) {
			
			#if KEGG ID is there, add hash entry
			#will reference the universal BIGG ID and will keep our own compartment info
			#just a bit different format
			$kegg_to_bigg_rxns{$1} = $line[0];;
			
			#note that in the BIGG reaction IDs there is no id without compartments
			#will remove compartment for id by removing last underscore and everything after it
			#using non-greeding regex
			$kegg_to_bigg_mets{$1} =~ s/_\w+?//;
			
		}
		
	}
	
}

#by this point reaction hash should be done

#now read the original model
open(ORIGMODEL, "<iEde2091.txt") or die "could not open the original iEde2091 model, reason: $!\n";
chomp(my @orig_model = <ORIGMODEL>);

#create a place to write the copy model with BIGG IDs
open(BIGGMODEL, ">iEde2091_BIGG.txt") or die "could not create/write to the BIGG version of iEde2091, reason: $!\n";

#for each line:
#get the reaction ID, convert to BIGG if in hash
#for each metabolite ID, convert to BIGG if in hash
#write new line to output file

#will keep count of fraction of reactions and metabolites substituted
my $rxns_subbed = 0;
my $total_rxns = 0;

#hash, key metabolite id, value is 1 if subbed, 0 otherwise
my %met_subbed = ( );

#array, stores all unique metabolite values
my @all_mets = ( );

for(my $c = 0; $c <= $#orig_model; $c++) {
	
	#skip formatting for comment lines, just parrot them
	if($orig_model[$c] =~ /^\#/) {
		
		printf BIGGMODEL "%s\n", $orig_model[$c];
		
	} else {
	
		#read the line and split into reaction id and stoichiometry
		(my $rxn_id, my $stoich) = split /\t/, $orig_model[$c];
		
		#for the reaction id, check if BIGG equivalent, replace if there is
		#first, get bareword id without compartments
		(my $bare_id, my $compartment) = split /\[/, $rxn_id;
		
		#my $bare_id =~ s/\[//;
		my $compartment = s/^/\[/;
		
		#if BIGG id exists
		if(exists $kegg_to_bigg_rxns{$bare_id}) {
			
			#replace with BIGG id
			$rxn_id = $kegg_to_bigg_rxns{$bare_id}.$compartment;
			$rxns_subbed++;
			$total_rxns++;
			
		} else {
			
			$total_rxns++;
			
		}
		
		#print the reaction id to the output
		printf BIGGMODEL "%s\t", $rxn_id;
		
		#next do the replacements for compounds
		#since have some propritary compounds, compounds will be identified by splitting the reaction into
		#products and reactants, then splitting by spaces, compounds are the odd numbered indeces
		(my $reactants, my $products) = split /\-/, $stoich;
		
		#also get the arrow
		my $arrow = "";
		
		#reversible arrow is the most exclusive case
		if($stoich =~ /\<\-\>/) {
			
			$arrow = "<->";
			
		} elsif ($stoich =~ /\<\-/) {
			
			$arrow = "<-";
			
		} else {
			
			#if not reversible or backwards must be forwards
			$arrow = "->";
			
		}
		
		#remove the arrowheads from both
		$reactants =~ s/\<//;
		$products =~ s/\>//;
		
		#split at spaces now
		my @react_arr = split /\s+/, $reactants;
		my @prod_arr = split /\s+/, $products;
		
		#for each compound, check if has BIGG id
		#start with 1 as 1 will be the first compound, increment by 3 (skips coefficients and "+" signs)
		for(my $d = 1; $d <= $#react_arr; $d+=3) {
			
			#get bareword compound
			(my $bare_comp, my $compart) = split /\[/, $react_arr[$d];
			
			#remove what got split and fix compartment
			$bare_comp =~ s/\[//;
			$compart =~ s/^/\[/;
			
			if(exists $kegg_to_bigg_mets{$bare_comp}) {
				
				$react_arr[$d] = $kegg_to_bigg_mets{$bare_comp}.$compart;
				
				#if entry in subbed does not exist, add it 
				if(not exists $met_subbed{$react_arr[$d]}) {
					
					$met_subbed{$react_arr[$d]} = 1;
					
					#add to array of all metabolites
					push @all_mets, $react_arr[$d];
					
				} #otherwise do nothing, should be recorded correctly first time
				
			} else {
				
				#if entry in subbed does not exist, add it 
				if(not exists $met_subbed{$react_arr[$d]}) {
					
					$met_subbed{$react_arr[$d]} = 0;
					
					#add to array of all metabolites
					push @all_mets, $react_arr[$d];
					
				} #otherwise do nothing, should be recorded correctly first time
				
			}
			
		}
		
		#need to shift the first emelement of the product array out, is nothing
		shift @prod_arr;
		
		#do the same for the products
		for(my $e = 1; $e <= $#prod_arr; $e+=3) {
			
			#get bareword compound
			(my $bare_comp, my $compart) = split /\[/, $prod_arr[$e];
			
			#remove what got split and fix compartment
			$bare_comp =~ s/\[//;
			$compart =~ s/^/\[/;
			
			if(exists $kegg_to_bigg_mets{$bare_comp}) {
				
				$prod_arr[$e] = $kegg_to_bigg_mets{$bare_comp}.$compart;
				
				#if entry in subbed does not exist, add it 
				if(not exists $met_subbed{$prod_arr[$e]}) {
					
					$met_subbed{$prod_arr[$e]} = 1;
					
					#add to array of all metabolites
					push @all_mets, $prod_arr[$e];
					
				} #otherwise do nothing, should be recorded correctly first time
				
			} else {
				
				#if entry in subbed does not exist, add it 
				if(not exists $met_subbed{$prod_arr[$e]}) {
					
					$met_subbed{$prod_arr[$e]} = 0;
					
					#add to array of all metabolites
					push @all_mets, $prod_arr[$e];
					
				} #otherwise do nothing, should be recorded correctly first time
				
			}
			
		}
		
		#make reactants array hold everthing
		push @react_arr, $arrow;
		push @react_arr, @prod_arr;
		
		#combine the reactants array into a single string
		my $subbed_stoich = join " ", @react_arr;
		
		#write to output file
		printf BIGGMODEL "%s\n", $subbed_stoich;
		
		#now repeat for each reaction
		
	}
	
}

my $frac_rxn_sub = ($rxns_subbed / $total_rxns) * 100;

printf "SUCCESS\nPercent reactions substituted: %s% (%s of %s)\n", $frac_rxn_sub,$rxns_subbed, $total_rxns;

my $total_mets;
my $subbed_mets;

for(my $f = 0; $f <= $#all_mets; $f++) {
	
	$subbed_mets+=$met_subbed{$all_mets[$f]};
	$total_mets++;
	
}

my $frac_met_sub = ($subbed_mets / $total_mets) * 100;

printf "Percent metabolites substituted: %s% (%s of %s)\n", $frac_met_sub, $subbed_mets, $total_mets;


