############################
# Written by: Brett Pickett at JCVI (May 2016)
# Modified by: Jacob Porter at UVa (Summer 2021)
# Input 1: multiple sequence alignment in fasta format and
# Input 2 TSV file with ID (matching fasta header) in first column and each metadata category in additional columns
# Output 1: meta-CATS Goodness of Fit test results for each metadata category
# Output 2: meta-CATS Test of Independence results for each metadata category
# Output 3: ID-metadata-sequence input file required for meta-CATS (for each metadata category)
# Output 4: Summary table showing number of significant sequence positions for each metadata category and for the combination of all categories
############################
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use IPC::Run qw(run);
# use Statistics::R;

my ($seqFile, $metaDataFile, $seqType, $pvalue, $outputDir) = @ARGV;

# my $metaDataFile = 'metadata.tsv';
# my $seqFile = 'input.afa'; #aligned sequences


my %categories;
my $lineNum = 0;
my $numCategories = 0;
my @headers;
my @data;
my @fileNames;

#load metadata into hash
open (METADATA, "<$metaDataFile") || die "$metaDataFile: $!\n";  #check to see if file exists
print "Loading Metadata...\n";

while (my $line = <METADATA>) {     #read all lines of the file
	chop $line;
	if($lineNum ==1){ # read metadata lines and automatically populate 2D-hash`
		@data=split /\t/,$line;
		for (my $i = 1; $i < $numCategories; $i++){
			$categories{$data[0]}{$headers[$i]} = $data[$i]; #1D = seqID, 2D = metadata category, value = attribute
		}
	}
	elsif ($lineNum == 0){ #read header line and count number of metadata categories
		chomp $line;
		$lineNum = 1;
		@headers=split /\t/,$line;
		$numCategories = @headers;
	}
}
close (METADATA);

# read in sequence data and store in same hash as metadata
open (SEQS, "<$seqFile") || die "$seqFile: $!\n"; 	#check to see if file exists
my $tempName;
print "Loading Sequence Data...\n";

while (my $line1 = <SEQS>) {		#read all lines of the file
	chomp $line1;
	$line1
	if ($line1 =~m/^>(\S+)/){
		$tempName = $1;
		#print "$tempName\n";
	}
	else {
		$categories{$tempName}{Seq} .= $line1;
	}
}
close (SEQS);

#add 1 space between each sequence letter to format for downstream analysis
foreach my $outer (keys %categories) {
	if (exists $categories{$outer}{Seq}){
		$categories{$outer}{Seq} = join(" ",split(//,$categories{$outer}{Seq}));
		my $space = ' ';
		$categories{$outer}{Seq} = $space.$categories{$outer}{Seq};
	}
}

# calculate number of different fields for each category
my %categoryCounts;
for (my $k = 1; $k < $numCategories; $k++){
	foreach my $outer (keys %categories) {
		$categoryCounts{$headers[$k]}{$categories{$outer}{$headers[$k]}}++;
	}
}


#assign group numbers in new hash
my %assignments;
my @sortedOuter = sort keys %categoryCounts;
my $outerLength = @sortedOuter;
my @sortedInner;
my $innerLength;
my $number = 1;
foreach my $outer2 (@sortedOuter){
	@sortedInner = sort keys %{$categoryCounts{$outer2}};
	$innerLength = @sortedInner;
	$number = 1;
	foreach my $inner2 (@sortedInner) {
		$assignments{$outer2}{$inner2}= $number;
		$number++;
	}
}

# construct separate output files for each category and sequence
foreach my $category1 (@sortedOuter){#metadata category
	push(@fileNames, "$outputDir/$category1.txt");
	open (OUTFILE, ">$outputDir/$category1.txt") || die "$category1.txt: $!\n";
	my @seqLabels = sort keys %categories;
	foreach my $seq1 (@seqLabels){
		my $attribute = $categories{$seq1}{$category1}; #points to numerical assignment representing metadata attribute
		#Only output the label, category, and sequence for labels that have BOTH category AND sequence data
		if(exists $categories{$seq1}{Seq}){
		#	#print "present: $seq1\n";
			print OUTFILE "$seq1\t$assignments{$category1}{$attribute}\t$categories{$seq1}{Seq}\n";
		}
		else{
		#	#print "absent: $seq1\n";
			next;
		}
	}
	close (OUTFILE);
}

#Create a communication bridge with R and start R to run meta-CATS code (in R)
# my @results;
# my $R = Statistics::R->new();
print "Running Statistics...\n";
foreach my $file1 (@fileNames){#metadata category
	open (MSAFILE, "<$file1") || die "$file1: $!\n";
	my $boolean1 = 0;
	my $firstloop = 0;
	my $oldCategoryNum = 0;
	while (my $line1 = <MSAFILE>) {		#read all lines of the file
		chomp $line1;
		if ($line1 =~m/^.+?\t(.+)?\t.+/){
			my $newCategoryNum = $1;
			if ($firstloop == 0){
				$oldCategoryNum = $newCategoryNum;
				$firstloop = 1;
			}
			elsif ($oldCategoryNum != $newCategoryNum){
				$boolean1 = 1;
			}
			elsif ($oldCategoryNum == $newCategoryNum){
				next
			}
		}
	}
	close (MSAFILE);
	if($boolean1 == 1){
		print "\t$file1\n";
		my $base = basename($file1);
		$base = substr($base, 0, length($base) - 4);
		my $ok = run("Meta-CATS", $file1, $seqType, $pvalue, $outputDir, $base);
		if (!$ok)
		{
			die "Running Meta-CATS.R failed.";
		}
		copy_to_tsv($outputDir, "$base-chisqTable", $pvalue);
    	copy_to_tsv($outputDir, "$base-mcTable", $pvalue);
		unlink("$outputDir$base-chisqTable.txt") or warn "Unable to unlink $!";
		unlink("$outputDir$base-mcTable.txt") or warn "Unable to unlink $!";
		# $R->set('inFilename', $file1);
		# my $Rrun1 = $R->run_from_file("chisq_ViPR-perl_BEP.R");
		# my $returnedFilename = $R->get('outfilename1');
		# push(@results, "$returnedFilename");
	}
	elsif($boolean1 == 0){
		next;
	}
}

sub copy_to_tsv {
    my($work_dir, $basename, $pvalue) = @_;
    my $filename = "$work_dir$basename.txt";
    my $path = "$work_dir$basename.tsv";
    open my $fh, '<', $filename or die "Cannot open $filename: $!";
    open(IN, '>', $path) or die "Cannot open $path: $!";
    my $sel = 0; # Location of p-value.
    my $stuff = "\t"; # Delimiter for final columns.
    if (substr($basename, length($basename) - 7) eq "mcTable") {
        $sel = 1;
        $stuff = ",";
        print IN "Position\tMultiple_comparison_p-value\tSignificant\tGroups\n";
    } else {
        $sel = 2;
        print IN "Position\tChi-square_value\tP-value\tSignificant\tDegrees_of_freedom\tFewer_5\tResidue_Diversity\n";
    }
    my $count = 0;
    while ( my $line = <$fh> ) {
        next if(length($line) <= 1);
        if (substr($line, 0, 1) eq "\"") {
            next;
        }
        chomp $line;
        $count = $count + 1;
        my @columns = split(/\t/, substr($line, 5));
        # Remove whitespace and some formatting from columns.
        for ( my $i = 0; $i < scalar(@columns); $i++ ) {
            my $thing = $columns[$i];
            $thing=~ s/^\s+|\s+$//g;
            if ($i == scalar(@columns) - 1) {
                $thing =~ s/_/\x20/g;
            }
            @columns[$i] = $thing;
        }
        # Add a column if the position is statistically significant or not.
        my $obs_p_value = $columns[$sel];
        my $sig = "N";
        if ($obs_p_value ne "NA" and $obs_p_value < $pvalue) {
            $sig = "Y";
        }
        # Join the columns together and print the results in a TSV file.
        my $final = join($stuff, map { defined ? $_ : '' } @columns[$sel+1..scalar(@columns)]);
        $final =~ s/[,|\t]$//;
        print IN join("\t", @columns[0..$sel])."\t$sig\t".$final."\n";
    }
    close(IN);
    close($fh);
    return $count;
}

# #Parse all results files and save summaries of those that have highest number of significant positions
# print "Parsing Results...\n";
# open (SUMMARY, ">AA_Summary_Results-$metaDataFile.txt") || die "AA_Summary_Results-$metaDataFile.txt: $!\n";
# open (RESULTSALL, ">1AA_Summary_Results_ALL-$metaDataFile.txt") || die "AA_Summary_Results_ALL-$metaDataFile.txt: $!\n";
# print SUMMARY "Category Counts\nMetadata Category\t# Significant Positions\t# Total Positions\t% for each Category (# Sig / # Total )\n";
# print RESULTSALL "Position	Chi-Square Score	P-Value	Degrees of Freedom	Sparse Table (i.e. <5 of any residue)	Residue Diversity Between Groups\tResults_Filename\n";
# my %resultsStats;
# my @resultsArray;
# my $resultsLength = @results;
# my $numSigResults =0;
# for (my $p = 0; $p < $resultsLength; $p++){#metadata category
# 	open (RESULTSDATA, "<$results[$p]") || die "$results[$p]: $!\n";
# 	print "\t$results[$p]...\n";
# 	while (my $line4 = <RESULTSDATA>) {     #read all lines of the file
# 		chomp $line4;
# 		if ($line4=~m/^ .+/){
# 			@resultsArray=split /\t/,$line4;
# 			$resultsArray[0] =~ s/^\s+|\s+$//g;
# 			$resultsArray[2] =~ s/^\s+|\s+$//g;
# 			if($resultsArray[2] <= 0.05){
# 				$resultsStats{$results[$p]}{sig}++;
# 				print RESULTSALL "$line4\t$results[$p]\n";
# 			}
# 		}
# 		else{
# 			next;
# 		}
# 		$resultsStats{$results[$p]}{total}++;
# 	}
# 	my $quotient;
# 	if(defined $resultsStats{$results[$p]}{sig}){
# 		$quotient = 100*($resultsStats{$results[$p]}{sig}/$resultsStats{$results[$p]}{total});
# 		$numSigResults += $resultsStats{$results[$p]}{sig};
# 		print SUMMARY "$sortedOuter[$p]\t$resultsStats{$results[$p]}{sig}\t$resultsStats{$results[$p]}{total}\t";
# 		printf SUMMARY ("%.1f",$quotient);
# 		print SUMMARY "%\n";
# 	}
# 	else{
# 		print SUMMARY "$sortedOuter[$p]\tNone\tNA\n";
# 	}
# 	close (RESULTSDATA);
# }
# my @percValues;
# #add all significant sites for each category and calculate percentages for each
# print SUMMARY "\nComprehensive Counts:\nMetadata Category\t# Significant Positions in Single Category\t# Significant Positions from ALL Categories\tPercentage for ALL Categories\n";
# for (my $q = 0; $q < $resultsLength; $q++){#metadata category
# 	if(defined $resultsStats{$results[$q]}{sig}){
# 		print SUMMARY "$sortedOuter[$q]\t$resultsStats{$results[$q]}{sig}\t$numSigResults\t";
# 		$resultsStats{$results[$q]}{totalPerc} = 100*($resultsStats{$results[$q]}{sig}/$numSigResults);
# 		$resultsStats{$results[$q]}{totalPerc} = sprintf "%.1f", $resultsStats{$results[$q]}{totalPerc};
# 		print SUMMARY "$resultsStats{$results[$q]}{totalPerc} %\n";
# 		push(@percValues, "$resultsStats{$results[$q]}{totalPerc}");
# 	}
# }
# close (SUMMARY);
# $R->set('lbls', \@sortedOuter);
# $R->set('slices', \@percValues);
# $R->set('inFilename3', $metaDataFile);
# # my $Rrun2 = $R->run_from_file("pieChart.R");
# $R->stop();

# print "Complete\n";
# print "Number of significant results for each metadata category can be found in:\n\"AA_Summary_Results-$metaDataFile\"\n\n";
# # print "Pie chart of results can be found in:\n\"AA-$metaDataFile-PieChart.pdf\"\n";