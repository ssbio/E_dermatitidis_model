#usr/bin/perl -w

#Created by Wheaton Schroeder to use protein names from table downloaded from NCBI database
#will use protein name as a search string in the Brenda Database. Takes an NCBI downloaded
#table of protein details (for E. dermatitidis the URL is https://www.ncbi.nlm.nih.gov/genome/proteins/2962?genome_assembly_id=34221&gi=-1)
#which has the following transformations:
#1) replace all commas with semicolons (necessary for csv formatting not to be thrown off)
#2)replace all tabs (\t) with commas (easily done in notepad ++)
#3) save the resultant file as a .csv
#this allows for easy of manual reading (as it may be opened in Excel) and ease of programmatic
#reading
#As this program cobers thousands of ORFs and all which are non-hypothetical are searched for in Brenda,
#the runtime of this program is by necessity long

#allows the use of searchBrendaProtein subroutine in the text common_function.pl
use LWP::UserAgent;
require 'common_functions.pl'; 

$headerLine = " ";
@linesToWrite = ( );
$ProtNameCol = 9; #column number where the protein name is. indexing begins at 0

#opens the file which contains the protein data from NCBI
if ( ! open PROTDATA, "<NCBIProteinDetails.csv") {
	die "Could not find/access NCBIProteinDetails.csv, reason: $!"; #fatal error
}

#opens the file to write new protein data
if ( ! open PROTNEW, ">NCBIProtECno.csv") {
	die "Could not create NCBIProtECno.csv, reason: $!"; #fatal error
}

#read the lines from the file
chomp(@lines = <PROTDATA>); #chomp removes white spaces at ends of each array element

#read in first line as the header
$headerLine = $lines[0];

#foreach line after the header
for (my $i = 0; $i <= $#lines - 1; $i++) {
	#search Brenda for the protein in each line
	$ECno = searchBrendaProtein($lines[$i+1], $ProtNameCol);
	$linesToWrite[$i] = $lines[$i+1].",".$ECno;
}

#begin writing protein information to NCBIProtECno.csv
select PROTNEW;
printf "%s, Brenda EC#s\n", $headerLine;

for ($j = 0; $j <= $#linesToWrite; $j++) { #for each protein line
	printf "%s\n", $linesToWrite[$j]; #write that line to the document
}

#after all lines are written, close the file	
close PROTNEW;