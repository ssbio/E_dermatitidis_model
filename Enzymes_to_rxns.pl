#usr/bin/perl -w

#Written by: Wheaton Schroeder
#Latest version: 04/24/2020

#Written to take a list of enzyme commission numbers with appended compartmentalization
#and turn it into a list of reaction stoichiometeries with included compartmentalization

use strict;
use LWP::UserAgent;

#read the input enzyme file
open(ECLIST, "<EC_list_Ede.txt") or die "could not open the list of enzyme commission numbers, reason: $!\n";
chomp(my @ec_list = <ECLIST>);

#prepare output file
open(RXNS, ">rxns_out.txt") or die "could not create/write the rxns_out.txt file, reason: $!\n";

#define user agent for doing internet things
my $UserAgent = LWP::UserAgent->new;

#array to store full list of reactions with compartment
my @all_rxns = ( );

#array to store all fully processed stoichiometery strings
my @all_stoich = ( );

#for each enzyme in the list
for(my $a = 0; $a <= $#ec_list; $a++) {

	#read the line
	my $temp = $ec_list[$a];
	
	#split the line into compartment and enzyme commission numbers
	(my $ec_num, my $comp) = split /\[/, $temp;
	
	#add the bracket back into the compartment
	$comp = "[".$comp;
	
	printf "EC number: %s\n", $ec_num;
	printf "compartment: %s\n", $comp;
	
	#look up which reactions are cross-linked with the EC number in KEGG using KEGG's restful API
	my $url = "http://rest.kegg.jp/link/reaction/".$ec_num;
	my $LinkPage = $UserAgent->get($url);
	my $LinkPage = $LinkPage->content;
	
	#should look like a list of enzyme classification numbers and reaction numbers, so need to 
	#get a list of reaction numbers
	my @rxns = split /\n/, $LinkPage;
	
	for(my $b = 0; $b <= $#rxns; $b++) {
	
		#remove the enzyme commission number and "rn:" string rom the entry so the raw reaction label is all that is left
		$rxns[$b] =~ s/^ec:$ec_num\s+rn://;
		
		printf "identified reaction: %s\n", $rxns[$b];
		
		#go to the web page to get stoichiometery only if isn't aleady in the array
		my $rxn_total = $rxns[$b].$comp;
		
		if(!($rxn_total ~~ @all_rxns)) {
		
			#add this reaction to the list of all reactions then
			push @all_rxns, $rxn_total;
			
			#now that we have the barework reaction, we can get the stoichiometery from KEGG
			my $url2 = "http://rest.kegg.jp/get/".$rxns[$b];
			my $RxnPage = $UserAgent->get($url2);
			my $RxnPage = $RxnPage->content;
			
			#holds the stoichiometery match for now
			my $stoich = "";
			
			#look for the stoichiometery
			while($RxnPage =~ m/EQUATION\s+(.+?)\n/g) {
			
				#if here then the stoichiometery was found
				$stoich = $1;
				
				#now we need to properly format the raw stoichiometery
				#replace the arrow format
				$stoich =~ s/=/-/;
				
				#next, add in explicit "1"s where there are implicit "1"s
				$stoich =~ s/\+\sC/+ 1 C/g;
				$stoich =~ s/\>\sC/> 1 C/g;
				$stoich =~ s/^C/1 C/g;
				
				#next, add compartmentalization to the stoichiometery
				$stoich =~ s/(C\d\d\d\d\d)/$1$comp/g;
				
				printf "formatted stoichiometery: %s\n\n", $stoich;
				
				#stoichiometery should be propery formatted by now, add to list of output lines
				push @all_stoich, $stoich;
				
			
			}
		
		}		
		
	}

}

#by this point, should have all reaction labels and stoichiometries
for(my $d = 0; $d <= $#all_rxns; $d++) {

	printf RXNS "%s\t%s\n", $all_rxns[$d], $all_stoich[$d];

}
