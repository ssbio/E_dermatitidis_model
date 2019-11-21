#! usr/bin/perl -w
#
#Created to process a minimalist XML model and turn it into
#a csv document which will be easily converted to an excel
#file, the result should be useful for early metabolic models
#i.e. circa 2008
#All output files will be written to to same directory in which this
#script is contained.
#-w allows for warnings in first line, remove to remove warnings
#using file metModel.xml, runtime for the program was measured as
#9 minutes 55 seconds, this provides a runtime baseline
#most runtime comes from interacting with the KEGG webpages
#Written by: Wheaton Schroeder

#Input File: metModel.xml
#	should be formatted as a normal xml file output for a metabolic
#	model. Metabolites should be labeled as "species", reactions 
#	Should be labeled as "reactions", etc. See a member of the lab 
#	group of Rajib Saha to verify formatting. The input model must be
# 	saved with this file name for the code to function. Ideally this
#	will be the defaul format comming from COBRA

#Output File One: metabolites.csv
#	A comma-seperated values (CSV) document which can be easily opened
#	using microsoft excel. Contains each species, species ID (internal),
#	species ID consistent with KEGG database, species compartment, 
#	charge, and number and kinds of atoms for each metabolite. Once 
#	The file is created, open in excel, select "Save As", and select
#	current file extension for excel for full functionality

#Output File Two: reactions.csv
#	Contains reaction ID (internal), EC ID, stoichiometry, involved
#	metabolites, associated gene, and reaction bounds. As above can
#	"Save As" for full excel functionality after this program has 
#	been run. For ease, the workbooks could be combined.

#Output File Three: S.csv
#	Creates the "S" matrix for a flux balance analysis in a .csv fire
#	again can easily save as excel document and combine workbooks

#This allows for use of the package LWP, which is the library package
#for the worldwide web for perl. Specifically, this part of the package
#allows this program to get the content from the webpage
use LWP::UserAgent;

#open the xml file of the metabolic model or kill the program
if ( ! open MODEL, "<metModel.xml") {
	die "Could not find/access metModel.xml, reason: $!"; #fatal error
}

#create a new document metabolites.csv to write to. If this document
#exists aready, will replace the file! I do not expect this error unless
#the folder is designated "read only" or is open currently, error is fatal.
if ( ! open METABOLITES, ">metabolites.csv") {
	die "Could not create metabolites.csv, reason: $!"; #fatal error
}

#create a new document reactions.csv to write to. If this document
#exists aready, will replace the file! I do not expect this error unless
#the folder is designated "read only"
if ( ! open REACTIONS, ">reactions.csv") {
	die "Could not create reactions.csv, reason: $!"; #fatal error
}

#create a new document S.csv to write to. If this document
#exists aready, will replace the file! I do not expect this error unless
#the folder is designated "read only". error is fatal.
if ( ! open S, ">S.csv") {
	die "Could not create S.csv, reason: $!"; #fatal error
}

#stores model name
my $modelID = " ";
my $modelName = " ";
#stores metabolite information, each element a row on the metabolites.csv file
#this is a hash, which will allow descriptive strings to be accessed using internal
#metabolite ID. A hash is particularly useful for accessing a metabolites information
#from an ID given in reaction parameters
my %MetInfo = ( ); #empty list
#however, for S, a sense of order is important (hash is more like a barrel than a linear array)
#therefore a number to hash key is stored in this array
my @MetKeys = ( );
#stores reaction information, each element a row on the reactions.csv file
#this is a hash, which allows descriptive strings about the reaction to be accessed using
#internal Reaction ID. A hash is particularly useful for easy access
my %RxnInfo = ( ); #empty list
#however, for S, a sense of order is important (hash is more like a barrel than a linear array)
#therefore a number to hash key is stored in this array
my @RxnKeys = ( );
#stores unit info in an array
my $Units = "";
#stores top row for metabolies.csv file (column labels), also defines the
#order for data to be organized in @MetInfo
my @MetLabels = ( "ID (Internal)", "name", "Compartment", "Charge", "Has Formula?", "Chemical Formula", "KEGG ID (if program searched for it)");
#stores top row for reactions.csv (column labels), also defines order for
#data to be input
my @RxnLabels = ( "Reaction Name", "Reaction ID", "EC ID(s)", "Reaction", "Reversible?");
#stores atomic symbols order, the key is the atomic symbol, the value is the "place" of that symbol
#in the column header order, for instance if the colums were C,H,N,O,P,S, then giving "O" to this hash
#would result in "3" being returned (counting starts at 0)
my %AtomSym = ( );



#reads each line of input from metModel.xml, places each line in an array
chomp(@lines = <MODEL>); #chomp removes white spaces at ends of each array element
$MetCount = 0; #count to give order to hash keys when needed, e.g. in S
$RxnCount = 0; #count to give order to hash keys when needed, e.g. in S

for ($i = 0; $i <= $#lines; $i++) { #for each element in the @lines array, array indexing starts at 0
	#particularly, we are seeking info on compartments, species, genes,
	#and reactions, so search for those key words in the @Line array
	#This first block of code deals with the model name
	$_ = $lines[$i]; #set default variable to the line we are looking at now
	
	if(/id=\"/ && /model/) { #this block gets the model ID and Name
		#if the line has both strings, it is in the the model name and ID line
		$modelID = &getDelimitedText($_, "id=\"", "\""); #pulls substring from 4 to end of ID, not including quotation mark
		printf "\nfound model ID: %s \n", $modelID;
		$modelName = &getDelimitedText($_, "name=\"", "\""); #pulls model name
		printf "found model Name: %s \n", $modelName;
		
	} elsif(/listOfUnitDefinitions/) { #this block gets the units
		#in the list of unit definition blocks. Capture units until at end of listOfUnitDefinitions block
		#keeps a number for what line is being read
		my $j = $i; #starting with current index position
		#This is an infinite loop that will be broken once the line </listOfUnitDefinitions> is reached 
		LINE: while(1) {
			$_ = $lines[$j];
			last LINE if /\/listOfUnitDefinitions/;
			#we will pull the unit ID, more than that should not be required
			if (/id=\"/) {
				$Units = &getDelimitedText($_, "id=\"", "\""); #pulls substring from 4 to end of ID, not including quotation mark
			}
			$j = $j + 1;
		}
		$i = $j; #skip the lines we took care of or purposefully ignored above
	} elsif((/listofCompartments/) && !(/listofCompartments\//)) { #get compartment information
		#if there are compartments, find the names and information
		#the matching makes a point that the first string is matched but the second string is not
		#because the second string will be all that is present if there are no compartments
		#not sure if anything needs to be done with the compartments information
	} elsif((/listOfParameters/) && !(/listOfParameters\//)) {
		#again, may not be present, and not quite sure what to do with them at present. This
		#will probably be a future sophistication
	} elsif((/fbc:listOfObjectives/) && !(/fbc:listOfObjectives\//)) {
		#again, may not be present, and not quite sure what to do with them at present. This
		#will probably be a future sophistication
	} elsif (/listOfSpecies/) {
		#now we are at one of the big cunks of data, the species information
		#this will add strings into the @MetInfo array to, the strings will be formatted as would be
		#a line in metabolites.csv so that it can be written to one array element (line) at a time
		$TempMetInfo = "";
		$TempMetHandle = "";
		my $j = $i; #starting with current index position
		#This is an infinite loop that will be broken once the line </listOfSpecies> is reached 
		LINE: while(1) {
			$_ = $lines[$j];
			$name = " "; #initialize name for use in searching KEGG later
			last LINE if /\/listOfSpecies/;
			#we will pull the species information for all species
			$TempMetHandle = $MetCount.",";
			if (/id=\"/) { #get species ID
				$TempMetHandle = &getDelimitedText($_, "id=\"", "\""); #pulls substring from 4 to end of ID, not including quotation mark
				$TempMetInfo = $TempMetHandle.","; #add ID
			} 
			
			if (/name=\"/) { #add name information to met info string
				$name = &getDelimitedText($_, "name=\"", "\"");
				$name =~ s/,/;/g; #replace present commas in molecule names with semicolons to avoid them being taken apart into different columns
				$TempMetInfo = $TempMetInfo.$name.","; #pulls substring from 4 to end of ID, not including quotation mark
			} else {
				#no name given, leave name column blank
				$TempMetInfo = $TempMetInfo." ,";
			}
			if (/compartment=\"/) { #adds compartment info to met info string, if present
				$TempMetInfo = $TempMetInfo.&getDelimitedText($_, "compartment=\"", "\"").","; #pulls substring from 4 to end of ID, not including quotation mark
			} else {
				#no compartment given, leave compartment blank
				$TempMetInfo = $TempMetInfo." ,";
			}
			if (/fbc:charge=\"/) { #adds charge info to met info string, if present
				$TempMetInfo = $TempMetInfo.&getDelimitedText($_, "fbc:charge=\"", "\"").","; #pulls substring from beginning to end of ID, not including quotation mark
			} else {
				#no charge given, leave column blank
				$TempMetInfo = $TempMetInfo." ,";
			}
			if (/fbc:chemicalFormula/) { #adds chemical formula to met info string and quantity of each atom
				#get the chemical formula string
				$CF = &getDelimitedText($_, "fbc:chemicalFormula=\"", "\"");
				#adds TRUE string to Formula Given? colunmn, presence of "TRUE" or "FALSE" string will 
				#later be used to determine if screen scraping is necessary for chemical formula off KEGG database
				$TempMetInfo = $TempMetInfo."TRUE,".$CF.", NOT LOOKED UP,"; #pulls substring from beginning to end of ID, not including quotation mark
				#break chemical formula down to get types and counts of each element
				my ($AtomicSymsRef, $SubscriptsRef) = &breakdownCF($CF); #get the returned references
				@AtomicSyms = @$AtomicSymsRef; #now get the arrays that are referenced
				@Subscripts = @$SubscriptsRef; 
				#now we have atomic symbols and subscripts, and we need to add new symbols to the header @MetLabels, and make sure we add subscripts
				#in the correct order. We will create a new sub for that.
				my ($MetLabelsRef, $AtomSymRef, $formattedSub) = &updateSymsAndSubs(\@MetLabels, \%AtomSym, \@AtomicSyms, \@Subscripts);
				@MetLabels = @$MetLabelsRef; #update metabolic labels
				%AtomSym = %$AtomSymRef; #update atomic symbols hash
				$TempMetInfo = $TempMetInfo.$formattedSub." ,"; 
				
				#I AM WORKING HERE!
				
				
			} else {
				#no chemical formula is given, place false in Formula given? and append a blank column.
				#time to webcrawl to get the information
				#first, pull the model name, as that is what will be used the search the KEGG compound database
				#name is the second element in the CSV line $TempMetInfo, so pull the second line
				$_ = $name;
				$KEGGID = ""; #stores KEGG ID number for the compound
				$CF = ""; #stores chemical formula for chemical species
				s/,/;/g; #replace present commas in molecule names with semicolons to avoid them being taken apart into different columns
				#substitute underscores for spaces in the name, this way search may be successful
				s/_/ /g;
				chomp($_); #remove white spaces at beginning and end
				#if name is given, should contain a series of numbers, letters, parentheses, and hyphens
				#so, will try to capture a nonempty string of one or more of the those
				#ensure there is an alphanumeric character at the beginning of the string
				#then match anything else. If name not give, will not match an alphanumeric character
				#at the beginning of the string. KEGG does not provide charge. Manual curation required
				if (/^(|\()?\w+/) {
					#alright so a name is given. Let's find a chemical formula
					#to begin, remove qalifiers such as input/output, of, Artificial, Metabolite
					#also, replace semi-colons with commas for successful searching
					s/input\/output//ig; #case insensitive substitution
					s/ of //ig;
					s/artificial//ig;
					s/metabolite//ig;
					s/;/,/g;
					s/\.//g;
					#someone decided to use square brackets as well, so replacing square brackets as they mess up regular expressions
					s/\[//g;
					s/\]//g;
					#for some reason here, chomp didn't do the work of removing space characters from the end of the string, so do with substitution
					s/\s+$//g;
					s/^\s+//g;
					#replace all spaces with the string %20, since the name string is directly input to the URL
					$searchName = $_;
					$searchNameNoSpace = $searchName;
					$searchNameNoSpace = ($searchNameNoSpace =~ s/\s+//g);
					s/\s+/%20/g;
					#now the string is ready to search, store search string as $name
					$URLname = $searchName;
					#The below line is to present the user with a stream of output from running the script so they don'table
					#get impatient and terminate the program early. The bars on either side of the string ensure that
					#
					printf "Searching KEGG for chemical species |%s|\n", $searchName;
					#reset line
					$_ = $lines[$j];
					#name of search page dependent only on name of compound at one place, so can navigate directly
					#to the search results. Effectively the same as going to the compound home page, typing in 
					$searchResPage = "http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=compound&keywords=".$URLname."&page=1";
					$UserAgent = LWP::UserAgent->new; #creates new user agent to get webpage contents
					$searchResponse = $UserAgent->get($searchResPage); #get the webpage content
					$searchResponse = $searchResponse->content;
					$searchName =~ s/\(/\\\(/g;
					$searchName =~ s/\)/\\\)/g;
					$searchNameNoSpace =~ s/\(/\\\(/g;
					$searchNameNoSpace =~ s/\)/\\\)/g;
					if ($searchResponse =~ /(;|>)\s$searchName(;|<)/i || /(;|>)\s$searchNameNoSpace(;|<)/i) { #search for the name of the species bound between > <
						#between ; <, between > ; or between ; ; This is necessary because of common synonyms for the same 
						#chemical species may be listed. Keeping it to just the name bound up inside the symbols because
						#it is common that say "phosphate" will be after the name because of phosphorolation
						#now that the name has been found in the document, search for the instance of href"/dbget closest
						#to the name while still being in front of it. This contains a partial URL which when the string
						#http://www.genome.jp is put in front of it constitutes the URL for the compound page in KEGG
						$nameIndex = index($searchResponse,$&);#index of start of match
						#have to do three rear index searchers, to avoid hitting the Jmol or KegDraw URLs
						$partialURLindex = rindex($searchResponse, "href=\"", $nameIndex); #rear index search, gives last occurance of string before name
						$partialURLindex = rindex($searchResponse, "href=\"", $partialURLindex - 1); #rear index search, gives last occurance of string before name
						$partialURLindex = rindex($searchResponse, "href=\"", $partialURLindex - 1); #rear index search, gives last occurance of string before name
						$URLregionLength = $nameIndex - $partialURLindex;
						#substring containing from the beginning of the href match to beginning of name. This section will contain only one
						#correctly formatted partial URL, which will be captured shortly
						$containsURL = substr($searchResponse, $partialURLindex, $URLregionLength);
						$compoundURL = "";
						if ($containsURL =~ /href=\"\/dbget.+?\">/) {
							#captures compound URL, if has name will have URL
							$compoundURL = "http://www.genome.jp".&getDelimitedText($&, "href=\"", "\">");
							printf "Compound URL: %s\n", $compoundURL;
							if ($compoundURL =~ /C\d+$/) {
								#find KEGG ID in URL
								$KEGGID = $&;
							}
							$compoundResults = $UserAgent->get($compoundURL);
							$compoundResults = $compoundResults->content;
							if ($compoundResults =~ /Formula/i) {
							
								$formulaIndex = index($compoundResults, "Formula");
								$hasCF = substr($compoundResults, $formulaIndex, 10000); #capture full output from formula index to end
								if ($hasCF =~ />\s*([A-Z][a-zA-z0-9]+)\s*<br>/) {
									
									my ($AtomicSymsRef, $SubscriptsRef) = &breakdownCF($CF); #get the returned references
									@AtomicSyms = @$AtomicSymsRef; #now get the arrays that are referenced
									@Subscripts = @$SubscriptsRef; 
									#now we have atomic symbols and subscripts, and we need to add new symbols to the header @MetLabels, and make sure we add subscripts
									#in the correct order. We will create a new sub for that.
									
									my ($MetLabelsRef, $AtomSymRef, $formattedSub) = &updateSymsAndSubs(\@MetLabels, \%AtomSym, \@AtomicSyms, \@Subscripts);
									@MetLabels = @$MetLabelsRef; #update metabolic labels
									%AtomSym = %$AtomSymRef; #update atomic symbols hash
									printf "adding this to temp met info: %s\n", $formattedSub;
									$TempMetInfo = $TempMetInfo.$formattedSub." ,"; 
								} else {
									#no chemical formula given in KEGG
									$TempMetInfo = $TempMetInfo."FALSE, UNKOWN, KEGG PAGE DOES NOT LIST A FORMULA,";
								}
							}
						} else {
							#no link to follow for KEGG compound database page
							$TempMetInfo = $TempMetInfo."FALSE, NOT GIVEN, NO COMPOUND PAGE,";
						}
						
					} else { #search KEGG with no spaces, make sure that is not the problem
						$URLname = $searchNameNoSpace;
						$URLname =~ s/\\\(/\(/g;
						$URLname =~ s/\\\)/\)/g;
						$searchResPage = "http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=compound&keywords=".$URLname."&page=1";
						$UserAgent = LWP::UserAgent->new; #creates new user agent to get webpage contents
						$searchResponse = $UserAgent->get($searchResPage); #get the webpage content
						$searchResponse = $searchResponse->content;
						if ($searchResponse =~ /(;|>)\s$searchName(;|<)/i || /(;|>)\s$searchNameNoSpace(;|<)/i) { #search for the name of the species bound between > <
							#between ; <, between > ; or between ; ; This is necessary because of common synonyms for the same 
							#chemical species may be listed. Keeping it to just the name bound up inside the symbols because
							#it is common that say "phosphate" will be after the name because of phosphorolation
							#now that the name has been found in the document, search for the instance of href"/dbget closest
							#to the name while still being in front of it. This contains a partial URL which when the string
							#http://www.genome.jp is put in front of it constitutes the URL for the compound page in KEGG
							$nameIndex = index($searchResponse,$&);#index of start of match
							#have to do three rear index searchers, to avoid hitting the Jmol or KegDraw URLs
							$partialURLindex = rindex($searchResponse, "href=\"", $nameIndex); #rear index search, gives last occurance of string before name
							$partialURLindex = rindex($searchResponse, "href=\"", $partialURLindex - 1); #rear index search, gives last occurance of string before name
							$partialURLindex = rindex($searchResponse, "href=\"", $partialURLindex - 1); #rear index search, gives last occurance of string before name
							$URLregionLength = $nameIndex - $partialURLindex;
							#substring containing from the beginning of the href match to beginning of name. This section will contain only one
							#correctly formatted partial URL, which will be captured shortly
							$containsURL = substr($searchResponse, $partialURLindex, $URLregionLength);
							$compoundURL = "";
							if ($containsURL =~ /href=\"\/dbget.+?\">/) {
								#captures compound URL, if has name will have URL
								$compoundURL = "http://www.genome.jp".&getDelimitedText($&, "href=\"", "\">");
								if ($compoundURL =~ /C\d+$/) {
									#find KEGG ID in URL
									$KEGGID = $&;
								}
								$compoundResults = $UserAgent->get($compoundURL);
								$compoundResults = $compoundResults->content;
								if ($compoundResults =~ /Formula/i) {
									$formulaIndex = index($compoundResults, "Formula");
									$hasCF = substr($compoundResults, $formulaIndex, 10000); #capture full output from formula index to end
									if ($hasCF =~ />\s*([A-Z][a-zA-z0-9]+)\s*<br>/) {
										#find the chemical formula, save it as a match variable
										$CF = $1;
										#process chemical formula as normal
										$TempMetInfo = $TempMetInfo."FALSE,".$CF.",".$KEGGID.",";
										my ($AtomicSymsRef, $SubscriptsRef) = &breakdownCF($CF); #get the returned references
										@AtomicSyms = @$AtomicSymsRef; #now get the arrays that are referenced
										@Subscripts = @$SubscriptsRef; 
										#now we have atomic symbols and subscripts, and we need to add new symbols to the header @MetLabels, and make sure we add subscripts
										#in the correct order. We will create a new sub for that.
										my ($MetLabelsRef, $AtomSymRef, $formattedSub) = &updateSymsAndSubs(\@MetLabels, \%AtomSym, \@AtomicSyms, \@Subscripts);
										@MetLabels = @$MetLabelsRef; #update metabolic labels
										%AtomSym = %$AtomSymRef; #update atomic symbols hash
										printf "adding this to temp met info: %s\n", $formattedSub;
										$TempMetInfo = $TempMetInfo.$formattedSub." ,"; 
									
									} else {
										#no chemical formula given in KEGG
										$TempMetInfo = $TempMetInfo."FALSE, UNKOWN, KEGG PAGE DOES NOT LIST A FORMULA,";
									}
								}
							} else {
								#no link to follow for KEGG compound database page
								$TempMetInfo = $TempMetInfo."FALSE, NOT GIVEN, NO COMPOUND PAGE,";
							}
						
						} else {
							#compound has a name, but does not correspond to anything in the KEGG database
							$TempMetInfo = $TempMetInfo."FALSE, NOT GIVEN, NO SEARCH RESULTS,";
						}
					
					
					}
				
				} else {
					#no name is given, and since the ID take is arbitrary at best, can't look up compound
					$TempMetInfo = $TempMetInfo."FALSE, NOT GIVEN, NOT LOOKED UP,";
				}
				
			}
			$MetInfo{$TempMetHandle} = $TempMetInfo;
			$MetKeys[$MetCount] = $TempMetHandle;
			$MetCount = $MetCount + 1;
			$j = $j + 1;
		}
		$i = $j; #skip the lines we took care of or purposefully ignored above
	} elsif (/<reaction/) {
		#Temp variables declared out here since there is multiple reaction lines and multiple peices of 
		#information to grab
		$TempRxnInfo = " ";
		$TempRxnID = " ";
		$TempRxnName = " ";
		$TempRxnWritten = " ";
		$TempReversibility = "FALSE";
		$TempRxnUB = " ";
		$TempRxnLB = " ";
		@TempEC = ( );
		@TempReactants = ( );
		@TempReactantCoeff = ( );
		@TempProducts = ( );
		@TempProductCoeff = ( );
		#infinite loop to be terminated once </reaction> line is reached
		$j = $i; #create line counter, starting at current value of $i
		LINE: while (1) {
			$_ = $lines[$j];
			if (/<\/reaction>/) {
				#if have reached end of reaction block, pull everything together and update 
				#%ReactionInfo and @ReactionKeys
				$TempRxnInfo = $TempRxnName.","; #starts line with reaction name
				$TempRxnInfo = $TempRxnInfo.$TempRxnID.","; #adds reaction ID
				#list EC(s) if present
				if (exists $TempEC[0]) { #if an EC is listed
					for ($k = 0; $k <= $#TempEC; $k++) { #for each element in @TempEC
						$TempRxnInfo = $TempRxnInfo.$TempEC[$k]." "; #add the element seperated by a semicolon (so they stay in same cell)
					}
				} else { #no EC is given
					$TempRxnInfo = $TempRxnInfo."EC NOT GIVEN";
				}
				
				$TempRxnInfo = $TempRxnInfo.","; #end EC cell
				
				#Create a written reaction mechanism, start with reactants
				for ($k = 0; $k <= $#TempReactants; $k++) {
					$TempRxnWritten = $TempRxnWritten.$TempReactantCoeff[$k].$TempReactants[$k]." "; #adds coefficient, metabolite ID, and space
				}
				#now add the reaction arrow, if reversible double-headed, if not, forward
				if ($TempReversibility =~ /FALSE/) {
					$TempRxnWritten = $TempRxnWritten."=> ";
				} else {
					$TempRxnWritten = $TempRxnWritten."<=> ";
				}
				for ($k = 0; $k <= $#TempProducts; $k++) {
					$TempRxnWritten = $TempRxnWritten.$TempProductCoeff[$k].$TempProducts[$k]." "; #adds coefficent, metabolite ID, and space
				}
				#reaction is now written, add to temp info string, also add reversibility
				$TempRxnInfo = $TempRxnInfo.$TempRxnWritten.",".$TempReversibility.",";
				$RxnInfo{$TempRxnID} = $TempRxnInfo; #add info to the hash
				push @RxnKeys, $TempRxnID; #adds reaction ID to key order
				
				#reset temp variables
				$TempRxnInfo = " ";
				$TempRxnID = " ";
				$TempRxnName = " ";
				$TempRxnWritten = " ";
				$TempReversibility = "FALSE";
				$TempRxnUB = " ";
				$TempRxnLB = " ";
				$RxnCount = $RxnCount + 1; #update reaction count
				last LINE if /\/reaction/; #end while loop
			}
			
			if (/reaction.+id=\"/) { #if on a line that contains the ID
				$TempRxnID = &getDelimitedText($_, "id=\"", "\""); #store ID in temporary variable			
			} 
			
			if (/reaction.+name=\"/) {
				$TempRxnName = &getDelimitedText($_, "name=\"", "\""); #store ID in temporary variable
				#replace commas with semi-colons to keep name in same cell
				$TempRxnName =~ s/,/;/g;
			}
			
			if (/reversible=\"/) {
				$TempReversibility = &getDelimitedText($_, "reversible=\"", "\""); #store ID in temporary variable
			}
			
			if (/listOfReactants/) {
				$k = $j; #counter again, start at value of $j
				#infinite loop again
				LINE: while(1) {
					$_ = $lines[$k];
					if (/species=\"/) {
						#grab reactant species info
						$reactant = &getDelimitedText($_, "species=\"","\"");
						push @TempReactants, $reactant; #store ID in temporary variable
					}
					
					if (/stoichiometry=\"/) {
						#grab stoichiometric coefficient info
						$reacStoich = &getDelimitedText($_, "stoichiometry=\"", "\"");
						push @TempReactantCoeff, $reacStoich; #store ID in temporary variable
					}
					
					$j = $k;
					last LINE if /\/listOfReactants/; #break loop if reached end of reactant list
					$k = $k + 1;
				}
			}
			
			if (/listOfProducts/) {
				$k = $j; #counter again, start at value of $j
				#infinite loop again
				LINE: while(1) {
					$_ = $lines[$k];
					if (/species=\"/) {
						#grab product species info
						$product = &getDelimitedText($_, "species=\"", "\"");
						push @TempProducts, $product; #store ID in temporary variable
					}
					
					if (/stoichiometry=\"/) {
						#grab stoichiometric coefficient info
						$productStoich = &getDelimitedText($_, "stoichiometry=\"", "\"");
						push @TempProductCoeff, $productStoich; #store ID in temporary variable
					}
					
					$j = $k; #skip lines already processed/deliberately skipped
					last LINE if /\/listOfProducts/; #break loop if reached end of product list
					$k = $k + 1;
				}
			}
			
			if (/ENZYME_CLASSIFICATION:/) { #if enzyme classifications are listed in this line
				if (/[0-9]\.(-|[0-9])\.(-|[0-9])\.(-|[0-9]+)/g) { #captures correctly (or near correctly) formatted ec
					#allows for dashes instead of numbers for the last 3 digits as this appears to be common
					push @TempEC, $&;
					while(/[0-9]\.(-|[0-9])\.(-|[0-9])\.(-|[0-9]+)/g) { #searches for further ECs
						push @TempEC, $&;
					}
				} elsif (/Non enzymatic step/) { #if not an enzymatic reaction
					push @TempEC, "NON ENZYMATIC STEP";
				} elsif (/ASSOCIATION: <\/xhtmp:p>/) { #if nothing given for enzyme association
					push @TempEC, "INPUT STEP";
				} elsif (/Artificial/i) { #artificial reaction, ignore case because not consistent
					push @TempEC, "ARTIFICIAL REACTION";
				} elsif (/No EC/i) { #no EC ignore case because it is not consistent
					push @TempEC, "NO EC";
				} elsif (/ERG[0-9]+/i) { #for some reason some enzyme classifications are ERG then a string of digits
					push @TempEC, $&;
				} elsif (/SUR[0-9]+/i) { #for some reason some enzyme classifications are SUR then a string of digits
					push @TempEC, $&;
				} elsif (/:\s\w+</) { #one of the other random strings given for EC, should be the last statement
					if (/\w+/) { #now capture just the word
						push @TempEC, $&;
					}
				}
				#at this point should have all enzyme classification data
			}
			
			$j = $j + 1;
		}
	}
}

#before writing metabolites


#begin writing metabolite information to metabolites.csv
select METABOLITES;
my $header = join ",", @MetLabels;
printf "%s\n", $header;

for ($i = 0; $i <= $#MetKeys; $i++) { #for each metabolite key
	printf "%s\n", $MetInfo{$MetKeys[$i]}; #wriet that line to the document
}

#after all lines are written, close the file	
close METABOLITES;

select REACTIONS;
my $headerRxn = join ",",@RxnLabels;
printf "%s\n", $headerRxn;

for ($i = 0; $i <= $#RxnKeys; $i++) {
	printf "%s\n", $RxnInfo{$RxnKeys[$i]};
}

#afer all lines written, close the file
close REACTIONS;



#Below are the subroutines used in this code to save space and to allow other programs to use these subroutines
#Thes subroutines are used frequenty throughout the script above. Unfortunately in Perl there is no way to define the 
#number and type of inputs to a function, so it is not "fool-proof", so before using please read the comments and or 
#documentation.

#This subroutine pulls the text from between two delimiters, which is a very common operation in this code
#inputs:
#	$_[0] = string from which to pull the text
#	$_[1] = string which defines the beginning index (this string will not be pulled, will be skipped) beginning delimiter
#	$_[2] = string which defines the end index, as previous, will not be returned in pulled string
#outputs:
#	$soughtString = the string contained in $_[0] which lies between the delimiters $_[1] and $_[2]
sub getDelimitedText {
	$beginIndex = index($_[0], $_[1]) + length($_[1]); #gives index of beginning of string then adds # characters so not capturing delimiter
	$endIndex = index($_[0], $_[2], $beginIndex); #gives index of end delimiter
	$stringLength = $endIndex - $beginIndex;
	$soughtString = substr($_[0], $beginIndex, $stringLength); #pulls substring
	return $soughtString; #returns substring
}

#This subroutine takes in a string that is a chemical formula and breaks it down into individual atomic symbols
#and subscripts, where the subscripts and atomic symbols are in the same order. That is for the formula C6H12O6,
#$AtomicSym[1] = H, $Subscript[1]=12.
#inputs:
#	$_[0] = string that is the chemical formula, e.g. [C6H12O6]
#outputs:
#	\@AtomicSym = reference to store of atomic symbols found, e.g. [C, H, O]
#	\@Subscript = reference of store of substripts for found atomic symbols, if no subscript provided the default is 1, e.g. [6, 12, 6]
sub breakdownCF {
	#define arrays to work with
	my @AtomAndSub = ( );
	my @AtomicSym = ( );
	my @Subscript = ( );
	while($_[0] =~ /[A-z][a-z]?[0-9]*/g) { #match pattern to single atomic symbol, while optionally capturing subscript
		#for C6H12O6, the first interation should capture $& = "C6"
		push(@AtomAndSub, $&); #add the match to the atom and sub array
	}
	#now that the atomic symbol/subscript pairs have been split from the rest, split the symbol from the subscript
	for ($i = 0; $i <= $#AtomAndSub; $i++) { #for each item in @AtomAndSub
		while ($AtomAndSub[$i] =~ /[A-z][a-z]?/g) { #for each atomic symbol
			#add each symbol to the @AtomicSym array
			push @AtomicSym, $&;
		}
		$tempSub = 1; #set default subscript to 1
		#if there is a subscript other than one
		while ($AtomAndSub[$i] =~ /[0-9]+/g) {
			#then use that as the subscript value
			$tempSub = $&;
		}
		push @Subscript, $tempSub; #add substript to @Subscript
	}
	#now the atomic symbols and subscripts are broken up, return them as arrays
	return (\@AtomicSym, \@Subscript); #return the two arrays
	#this is actually a little complex in that it is not possible to return more than one array, so we are simply returning
	#an array of array references.
}

#This subroutine takes the atomic symbols and subscripts from breakdownCF and does two important things: 1) it adds new atomic symbols
#to the header for the metabolites document @MetLabels, and to the hash %AtomicSym and 2) formats the row of subscripts for the present species
#inputs:
#	\@MetLabels = reference to array of column labels for metabolite document. Will push new atomic symbols into it
#	\%AtomicSymbols = reference to hash of atomic symbols, key is symbol, value is number order. Will add new atomic symbols
#	\@AtomicSyms = array of atomic symbols present in the current species
# 	\@Subscripts = array of subscripts present in current species related to the atomic symbols
#outputs:
#	\@newMetLabels = reference to updated @MetLabels array
#	\%newAtomicSymbols = reference to updated %AtomicSymbols hash
#	$formattedSubs = string of subscripts formatted so that it is ready to append to the current metabolite row being generated
sub updateSymsAndSubs {
	$MetLabelsRef = $_[0]; #get the met labels reference
	$AtomicSymbolsRef = $_[1]; #get the atomic symbols reference
	$AtomicSymsRef = $_[2]; #get atomic syms for this species reference
	$SubscriptsRef = $_[3]; #get subscripts for this species reference
	my @newMetLabels = @$MetLabelsRef; #get array from reference
	my %newAtomicSymbols = %$AtomicSymbolsRef; #get symbols has from reference
	my @copyAtomicSyms = @$AtomicSymsRef;
	my @copySubscripts = @$SubscriptsRef;
	#creates an array in whic hto store what atoms are currently present
	my @PresentAtoms = ( );
	#output string defined here
	$formattedSubs = "";
	for ($i = 0; $i <= $#copyAtomicSyms; $i++) { #foreach atomic symbol in the species
		$TempAtomicSym = $copyAtomicSyms[$i];
		#1) check to see if it exists the the atomic symbol hash, if not add it to the hash, and the the labels
		if (! exists $newAtomicSymbols{$TempAtomicSym}) { #if the atomic symbol does not exist in the hash
			$numberofKeys = keys %newAtomicSymbols;
			$newAtomicSymbols{$TempAtomicSym} = $numberofKeys; #add new key to hash at the current number of keys
			#indexing begins at 0 as usual		
			push @newMetLabels, $copyAtomicSyms[$i]; #add new atomic symbol to end of metabolic labels
		}
		
		#add the subscript to the present symbols array at the index specified by the atomic symbol hash key
		#any species not present in the species will be undef, or in otherwords will not exist
		$PresentAtoms[$newAtomicSymbols{$copyAtomicSyms[$i]}] = $copySubscripts[$i]; 
	}
	
	
	#2) format the row for the metabolites document
	for ($i = 0; $i <= $#PresentAtoms; $i++) {
		#if the element in the present atoms array exists/is defined
		if (exists $PresentAtoms[$i]) {
			#then that atom is present and needs a subscript
			$formattedSubs = $formattedSubs.$PresentAtoms[$i].",";
		} else { #if it is not defined
			#it does not have any atoms of that symbol
			$formattedSubs = $formattedSubs." ,";
		}
	}

	#now the meat of the program is done, all that needs to be done is to return references to the updated hash and array
	#and to return the formatted string.
	return (\@newMetLabels, \%newAtomicSymbols, $formattedSubs);
}

#this subroutine is written to search the KEGG database for a specifica name. It also holds contingencies such as 
#a hash for amino acids (since they are commonly listed by their three letter designation), and a contingency for 
#protons (charge +1, formula hydrogen), proteins, polymers with a designated n-value, and finally a hash for common
#acronyms (as part of the amino acid hash)
#inputs:
#	$speciesName = chemical species name
#	$searchLine = line of input from which the species name is read
#outputs:
#	$KEGGID = KEGG ID+ for chemical compound, format C[0-9]+
#	$speciesCharge = formal charge of chemical species
#	$CF = chemical formula of the named compound
sub searchKEGGforCF {

	#first, create a hash for 3 letter protein codes that often appear in chemical species,
	#this will allow for easy lookup of amino acids. For most, a specific chirality is noted so that
	#search parameters may fit more exactly, since chirality does not affect the molecular formula
	#this hash gives KEGG ID, charge and molecular formula, allowing us to forgo a KEGG search
	#the return string may be split using the split function with delimeter :
	my %commonSpecies = (
		"Ala" => "L-Alanine:C00041:0:C2H7NO2",
		"Ile" => "L-Isoleucine:C00407:0:C6H13NO2",
		"Leu" => "L-Leucine:C00123:0:C6H12NO2",
		"Val" => "L-Valine:C00183:0:C5H11NO2",
		"Phe" => "L-Phenylalanine:C00079:0:C9H11NO2",
		"Trp" => "L-Tryptophan:C00078:0:C11H12N2O2",
		"Tyr" => "L-Tyrosine:C00082:0:C9H11NO3",
		"Asn" => "L-Asparagine:C00152:0:C4H8N203",
		"Cys" => "L-Cysteine:C00097:0:C3H7NO2S",
		"Gln" => "L-Glutamine:C00064:0:C5H10N2O3",
		"Met" => "L-Methionine:C00073:0:C5H11NO2S",
		"Ser" => "L-Serine:C00065:0:C3H7N03",
		"Thr" => "L-Threonine:c00188:0:C4H9NO3",
		"Asp" => "L-Aspartic Acid:C00049:0:C4H7NO4",
		"Glu" => "L-Glutamic Acid:C00025:0:C5H9NO4",
		"Arg" => "L-Arginine:C00062:0:C6H14N4O2",
		"His" => "L-Histidine:C00135:0:C6H9N3O2",
		"Lys" => "L-Lysine:C00047:0:C6H14N2O2",
		"Gly" => "Glycine:C00037:0:C2H5N02",
		"Pro" => "L-Proline:C00148:0:C5H9NO2",
		"ATP" => "Adenosine 5'-triphosphate:C00002:0:C10H16N5O13P3",
		"ADP" => "Adenosine 5'-diphosphate:C00008:0:C10H15N5O10P2",
		"AMP" => "Adenosine 5'-monophosphate:C00020:0:C10H14N5O7P1",
	);
	
	chomp($_[0]);
	$speciesName = $_[0];
	$searchLine = $_[1];
	$KEGGID = "";
	$speciesCharge = "";
	$CF = "";
	
	if (! ($speciesName =~ /^(|\()?\w+/)) { #if a word is not passed to the subroutine
		#allowing for chirality and alphanumeric characters given
		return (-1, -1); #return negative one's, no word is given so can't search
	}
	
	#otherwise, a word must have been passed, let's search for it
	#but first, we have to clear out "junk" words that are given in the name as descripttors
	#but will hinder the search in KEGG. Will also create two search strings, one with spaces and one
	#without spaces, since it does make a difference
	#to begin, remove qalifiers such as input/output, of, Artificial, Metabolite
	#also, replace semi-colons with commas for successful searching
	$speciesName =~ s/input//ig; #case insensitive substitution, remove input description
	$speciesName =~ s/output//ig; #case insensitive substitution, remove output description
	$speciesName =~ s/input//ig; #case insensitive substitution, remove input descriptions/ of //ig;
	$speciesName =~ s/artificial//ig;
	$speciesName =~ s/metabolite//ig;
	$speciesName =~ s/;/,/g; #replace semicolons because for CSV file writingsemicolons replace commas
	$speciesName =~ s/\.//g; #remove periods
	#someone decided to use square brackets as well, so replacing square brackets as they mess up regular expressions
	$speciesName =~ s/\[//g; #remove square brackets that interfere with KEGG searches
	$speciesName =~ s/\]//g;
	#for some reason here, chomp didn't do the work of removing space characters from the end of the string, so do with substitution
	$speciesName =~ s/\s+$//g;
	$speciesName =~ s/^\s+//g;
	#replace all spaces with the string %20, since the name string is directly input to the URL
	$searchName = $speciesName; #set a variable for searching the normal name
	$searchNameNoSpace = $searchName; #set a variable for searching name without spaces
	$searchNameNoSpace =~ s/\s+//g; #remove spaces to have a no spaces option for searching
	$searchName =~ s/\s+/%20/g; #replaces spaces with %20 which is the html character for a space, need for URL
	#give an output so the user does not become impatient and terminate early
	printf "Searching KEGG for chemical species |%s|\n", $speciesName;
	
	#now before we get too involved, let us see if we can forgo a KEGG search, by avoiding searches on
	#the provided chemical species by way of it being a proton, protein, or in common species hash
	if ($searchName =~ /proton/i) {
		#if the species string comtains the word proton, we can know the formula, and the charge
		#I looked up the below information manually
		$KEGGID = "C00080";
		$speciesCharge = "+1";
		$CF = "H";
		return ($KEGGID, $speciesCharge, $CF);
	} elsif ($searchName =~ /protein/i) { #if the species is named as a protien
		#then it will not be able to be found in KEGG compound, or even have a meaningful molecular
		#formula, as it may well vary on a species to species basis
		$KEGGID = "PROTEIN";
		$speciesCharge = "0";
		$CF = "AA POLYMER";
		return ($KEGGID, $speciesCharge, $CF);
	} else { #check if in common species hash
		#for each element in common species hash
		foreach $key (sort keys %commonSpecies) {
			#search to see if the name contains the hash key, but need to be careful so we just
			#pick up those instances where a protien or a form of ATP was intended
			#ATP/AMP/ADP will stand alone, however, it is possible that the two bound
			#amino acids may have the form "Ala-Leu" for example
			#although if stand alone will likely be full name such as "L-Leucine"
			#there are more sophisticated names with alternations of AA codes and molecule names,
			#but I have to be realistic and stop looking for specialized occurances at some point
			$hashOut = $commonSpecies{$key};
			#if species contains a three letter code anchored at either end
			if ($searchName =~ /-$key$/ || $searchName =~ /^$key-/ ) {
				#a hyphonated amino acid pair is present, call this sub recursively to find the other
				
				#get chemical formula for first
				
				#get chemical formula for second
				
				#add the chemical formulas together
				
				#now we have the information form both species, but they have a peptide bond, so remove
				#two hydrogens and one oxygen from combined chemical formula
				
				#finally, join the array into a cohesive string to return the chemical formula
				
				#return what we found, give KEGG ID for both amino acids
				
			} elsif ($hashOut =~ /$searchName/i) { #if the hash value for the key contains the search name
				#then the name fed must have been a full amino acid name, return the hash information
				($speciesName, $KEGGID, $speciesCharge, $CF) = split /:/, $hashOut;
				return ($KEGGID, $speciesCharge, $CF);
			}
			
			#otherwise, doesn't look like any common species are in here
		
		}
	
	}
	
	#alright, since we've gotten here, we haven't got it anywhere convenient, so we are looking it up 
	#in the KEGG database
	$searchResPage = "http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=compound&keywords=".$searchName."&page=1";
	#create a new user agent which pulls the html source code for the KEGG search page
	$UserAgent = LWP::UserAgent->new; #creates new user agent
	$searchResponse = $UserAgent->get($searchResPage);
	$searchResponse = $searchResponse->content; 
	#format search strings to avoid matching delimeters
	$searchName =~ s/\(/\\\(/g;
	$searchName =~ s/\)/\\\)/g;
	
	#if the html code for the webpage contains the search name surrounded by the proper symbols (to ensure not in search box)
	#between ; <, between > ; or between ; ; This is necessary because of common synonyms for the same 
	#chemical species may be listed. Keeping it to just the name bound up inside the symbols because
	#it is common that say "phosphate" will be after the name because of phosphorolation
	#now that the name has been found in the document, search for the instance of href"/dbget closest
	#to the name while still being in front of it. This contains a partial URL which when the string
	#http://www.genome.jp is put in front of it constitutes the URL for the compound page in KEGG
	if ($searchResponse =~ /(;|>)\s?$searchName\s?(;|<)/i || /(;|>)\s?$searchNameNoSpace\s?(;|<)/i) {
		$nameIndex = index($searchResponse,$&);#index of start of match
		#have to do three rear index searchers, to avoid hitting the Jmol or KegDraw URLs
		$partialURLindex = rindex($searchResponse, "href=\"", $nameIndex); #rear index search, gives last occurance of string before name, Jmol link index
		$partialURLindex = rindex($searchResponse, "href=\"", $partialURLindex - 1); #rear index search, gives last occurance of string before name, KegDraw link index
		$partialURLindex = rindex($searchResponse, "href=\"", $partialURLindex - 1); #rear index search, gives last occurance of string before name, compound page link index
		$URLregionLength = $nameIndex - $partialURLindex;
		#substring containing from the beginning of the href match to beginning of name. This section will contain only one
		#correctly formatted partial URL, which will be captured shortly
		$containsURL = substr($searchResponse, $partialURLindex, $URLregionLength);
		$compoundURL = "";
		if ($containsURL =~ /href=\"\/dbget.+?\">/) { #actually this isn't really an if, but helps for capturing
			#captures compound URL, if has name will have URL
			$compoundURL = "http://www.genome.jp".&getDelimitedText($&, "href=\"", "\">");
			printf "Compound URL: %s\n", $compoundURL;
			if ($compoundURL =~ /C\d+$/) { #captures compound ID at the end of the URL
				#find KEGG ID in URL
				$KEGGID = $&;
			}
			#create a new user agent to go to follow the link to the compound page
			$compoundResults = $UserAgent->get($compoundURL);
			$compoundResults = $compoundResults->content;
			#search 
			
			
			#I AM WORKING HERE!!!!!!!!!!!!!!!!!
			
			#check results for chemical formula
			if ($compoundResults =~ /Formula/i) {
							
								$formulaIndex = index($compoundResults, "Formula");
								$hasCF = substr($compoundResults, $formulaIndex, 10000); #capture full output from formula index to end
								if ($hasCF =~ />\s*([A-Z][a-zA-z0-9]+)\s*<br>/) {
									
									my ($AtomicSymsRef, $SubscriptsRef) = &breakdownCF($CF); #get the returned references
									@AtomicSyms = @$AtomicSymsRef; #now get the arrays that are referenced
									@Subscripts = @$SubscriptsRef; 
									#now we have atomic symbols and subscripts, and we need to add new symbols to the header @MetLabels, and make sure we add subscripts
									#in the correct order. We will create a new sub for that.
									
									my ($MetLabelsRef, $AtomSymRef, $formattedSub) = &updateSymsAndSubs(\@MetLabels, \%AtomSym, \@AtomicSyms, \@Subscripts);
									@MetLabels = @$MetLabelsRef; #update metabolic labels
									%AtomSym = %$AtomSymRef; #update atomic symbols hash
									printf "adding this to temp met info: %s\n", $formattedSub;
									$TempMetInfo = $TempMetInfo.$formattedSub." ,"; 
								} else {
									#no chemical formula given in KEGG
									$TempMetInfo = $TempMetInfo."FALSE, UNKOWN, KEGG PAGE DOES NOT LIST A FORMULA,";
								}
							}
						}
	
	}
	
	$searchResPageNS = "http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=compound&keywords=".$searchNameNoSpace."&page=1";
	$searchNameNoSpace =~ s/\(/\\\(/g;
	$searchNameNoSpace =~ s/\)/\\\)/g;
	

	
					 #search for the name of the species bound between > <
						#between ; <, between > ; or between ; ; This is necessary because of common synonyms for the same 
						#chemical species may be listed. Keeping it to just the name bound up inside the symbols because
						#it is common that say "phosphate" will be after the name because of phosphorolation
						#now that the name has been found in the document, search for the instance of href"/dbget closest
						#to the name while still being in front of it. This contains a partial URL which when the string
						#http://www.genome.jp is put in front of it constitutes the URL for the compound page in KEGG
						 else {
							#no link to follow for KEGG compound database page
							$TempMetInfo = $TempMetInfo."FALSE, NOT GIVEN, NO COMPOUND PAGE,";
						}
						
					}




}