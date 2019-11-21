#usr/bin/perl -w

#Created by Wheaton Schroeder to use protein names from table downloaded from UniProt database
#will use protein name as a search string in the Brenda Database. Takes a UniProt downloaded
#table of protein details (for E. dermatitidis the URL is 
#http://www.uniprot.org/uniprot/?query=taxonomy%3A%22Exophiala+dermatitidis+%28strain+ATCC+34100+%2F+CBS+525.76+%2F+NIH%2FUT8656%29+%28Black+yeast%29+%28Wangiella+dermatitidis%29+%5B858893%5D%22&sort=score
#which has the following transformations:
#1) before downloading table, ensure that a column is present for EC number (can be added by selecting
#	the "columns" button above the table)
#2)Once downloaded, replace all commas with semicolons (necessary for csv formatting not to be thrown off)
#3)replace all tabs (\t) with commas (easily done in notepad ++)
#4) save the resultant file as a .csv
#this allows for easy of manual reading (as it may be opened in Excel) and ease of programmatic
#reading
#As this program cobers thousands of ORFs and all which are non-hypothetical are searched for in Brenda,
#the runtime of this program is by necessity long

#allows the use of searchBrendaProtein subroutine in the text common_function.pl
use LWP::UserAgent;
require 'common_functions.pl'; 

$headerLine = " ";
@linesToWrite = ( );
$ProtNameCol = 3; #column number where the protein name is. indexing begins at 0

#opens the file which contains the protein data from NCBI
if ( ! open PROTDATA, "<Uniprot_E_dermatitidis.csv") {
	die "Could not find/access Uniprot_E_dermatitidis.csv, reason: $!"; #fatal error
}

#opens the file to write new protein data
if ( ! open PROTNEW, ">UniProtECno.csv") {
	die "Could not create UniProtECno.csv, reason: $!"; #fatal error
}

#read the lines from the file
chomp(@lines = <PROTDATA>); #chomp removes white spaces at ends of each array element

#read in first line as the header
$headerLine = $lines[0];

#foreach line after the header
for ($i = 0; $i <= $#lines - 1; $i++) {
	#first check if the line already has an EC number given, it may then be safely ignored
	#this is considered because UniProt contains EC number information when available
	#if at least one EC number is there, we can assue that it can be ignored
	#need to count EC numbers that use hyphens instead of normal numbers
	if ($lines[$i+1] =~ /[0-9]+\.([0-9]+|-)\.([0-9]+|-)\.([0-9]+|-)/ig) {
		#if the line has an enzyme, then there is no need to search for it
		$linesToWrite[$i] = $lines[$i+1];
	} else {
		#The lines are formatted such that there are a lot of unneccesary commas at the 
		#end of each csv line, denoting blank cells. I note that in lines with content
		#therefore I will remove all commas from the ends of the csv lines, esssentially
		#deleting the blank cells

		#while there is a comma at the end of the line
		while ($lines[$i+1] =~ /,$/g) {
			#remove it
			$lines[$i+1] =~ s/,$//g;
		}
		#search Brenda for the protein in each line
		$ECno = searchBrendaProtein($lines[$i+1], $ProtNameCol);
		$linesToWrite[$i] = $lines[$i+1].",".$ECno;
	}
	
}

#begin writing protein information to NCBIProtECno.csv
select PROTNEW;
printf "%s\n", $headerLine;

for ($j = 0; $j <= $#linesToWrite; $j++) { #for each protein line
	printf "%s\n", $linesToWrite[$j]; #write that line to the document
}

#after all lines are written, close the file	
close PROTNEW;