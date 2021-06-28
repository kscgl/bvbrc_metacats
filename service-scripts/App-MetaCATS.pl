#The MSA application with variance analysis.

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig;

use strict;
use P3DataAPI;
use Data::Dumper;
use File::Basename;
use File::Slurp;
use LWP::UserAgent;
use JSON::XS;
use JSON;
use IPC::Run qw(run);
use Cwd;
use Clone;
use URI::Escape;

my $script = Bio::KBase::AppService::AppScript->new(\&process_metacats);
my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;

my $rc = $script->run(\@ARGV);

exit $rc;


sub process_metacats
{
    my($app, $app_def, $raw_params, $params) = @_;

    print "Proc Metacats ", Dumper($app_def, $raw_params, $params);
    my $token = $app->token();
    my $data_api_module = P3DataAPI->new($data_api, $token);
    my $output_folder = $app->result_folder();

    #
    # Create an output directory under the current dir. App service is meant to invoke
    # the app script in a working directory; we create a folder here to encapsulate
    # the job output.
    #
    # We also create a staging directory for the input files from the workspace.
    #

    my $cwd = getcwd();
    my $work_dir = "$cwd/work";
    my $stage_dir = "$cwd/stage";

    -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
    -d $stage_dir or mkdir $stage_dir or die "Cannot mkdir $stage_dir: $!";

    my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;
    my $dat = { data_api => $data_api };
    my $sstring = encode_json($dat);

    #
    # Read parameters and discover input files that need to be staged.
    #
    # Make a clone so we can maintain a list of refs to the paths to be
    # rewritten.
    #
    my %in_files;
    my $params_to_app = Clone::clone($params);
    #
    # Count the number of files.
    #
    my $prefix = $params_to_app->{output_file};
    my $alignment_type = "na";
    if ($params_to_app->{alignment_type} eq "aligned_protein_fasta") {
        $alignment_type = "aa";
    }
    my $p_value = $params_to_app->{p_value};
    my $seqFile = basename($params_to_app->{alignment_file});
    my $metaDataFile = basename($params_to_app->{group_file});
    #
    # Write files to the staging directory.
    #
    my @to_stage;
    push(@to_stage, $params_to_app->{alignment_file});
    push(@to_stage, $params_to_app->{group_file});
    my $staged = {};
    if (@to_stage)
    {
        warn Dumper(\@to_stage);
        $staged = $app->stage_in(\@to_stage, $stage_dir, 1);
        # while (my($orig, $staged_file) = each %$staged)
        # {
        #     my $path_ref = $in_files{$orig};
        #     $$path_ref = $staged_file;
        # }
    }
    # Replace leading and trailing gaps with the '#' symbol.
    open my $fh, '<', "$stage_dir/$seqFile" or die "Cannot open $stage_dir/$seqFile: $!";
    open(OUT, '>', "$work_dir/$seqFile") or die "Cannot open $work_dir/$seqFile: $!";
    my $seq_string = "";
    my $header = "";
    while ( my $line = <$fh> ) {
        chomp $line;
        if (substr($line, 0, 1) eq ">") {
            if ($header) {
                print OUT "$header\n";
            }
            $header = $line;
            if ($seq_string) {
                $seq_string = replace_gaps($seq_string);
                print OUT "$seq_string\n"

            }
            $seq_string = "";
        } else {
            $seq_string = $seq_string.$line;
        }
    }
    if ($header) {
        print OUT "$header\n";
    }
    if ($seq_string) {
        $seq_string = replace_gaps($seq_string);
        print OUT "$seq_string\n"

    }
    # Run the analysis.
    my @cmd = ("metadata_parser", "$work_dir/$seqFile", "$stage_dir/$metaDataFile", $alignment_type, "$p_value", "$work_dir/");
    run_cmd(\@cmd);

    copy_to_tsv($work_dir, "chisqTable", $p_value);
    copy_to_tsv($work_dir, "mcTable", $p_value);
    my @output_suffixes = (
        [qr/Table\.tsv$/, "tsv"],
        );
    opendir(D, $work_dir) or die "Cannot opendir $work_dir: $!";
    my @files = sort { $a cmp $b } grep { -f "$work_dir/$_" } readdir(D);
    my $output=1;
    for my $file (@files)
    {
	for my $suf (@output_suffixes)
	{
	    if ($file =~ $suf->[0])
	    {
 	    	$output=0;
		my $path = "$output_folder/$file";
		my $type = $suf->[1];
		$app->workspace->save_file_to_file("$work_dir/$file", {}, "$output_folder/$file", $type, 1,
					       (-s "$work_dir/$file" > 10_000 ? 1 : 0), # use shock for larger files
					       $token);
	    }
	}
    }
    #
    # Clean up staged input files.
    #
    while (my($orig, $staged_file) = each %$staged)
    {
	unlink($staged_file) or warn "Unable to unlink $staged_file: $!";
    }
    return $output;
}

sub replace_gaps {
    my($line) = @_;
    my $start_len = 0;
    my $end_len = 0;
    if ($line =~ /^(-+)/) {
        $start_len = length($1);
    }
    if ($line =~ /(-+)$/) {
       $end_len = length($1);
    }
    substr($line, 0, $start_len) = '#' x $start_len;
    substr($line, length($line) - $end_len, $end_len) = '#' x $end_len;
    return $line;
}

sub copy_to_tsv {
    my($work_dir, $basename, $p_value) = @_;
    my $work_dir = $_[0];
    my $basename = $_[1];
    my $filename = "$work_dir/$basename.txt";
    my $path = "$work_dir/$basename.tsv";
    open my $fh, '<', $filename or die "Cannot open $filename: $!";
    open(IN, '>', $path) or die "Cannot open $path: $!";
    my $sel = 0; # Location of p-value.
    my $stuff = ""; # Delimiter for final columns.
    if ($basename eq "chisqTable") {
        $sel = 2;
        $stuff = "\t";
        print IN "Position\tChi-square_value\tP-value\tSignificant\tDegrees_of_freedom\tFewer_5\tResidue_Diversity\n";
    } else {
        $sel = 1;
        $stuff = ",";
        print IN "Position\tMultiple_comparison_p-value\tSignificant\tGroups\n";
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
            my $thing = @columns[$i];
            $thing=~ s/^\s+|\s+$//g;
            if ($i == scalar(@columns) - 1) {
                $thing =~ s/,_/,/g;
            }
            @columns[$i] = $thing;
        }
        # Add a column if the position is statistically significant or not.
        my $obs_p_value = $columns[$sel];
        my $sig = "N";
        if ($obs_p_value < $p_value) {
            $sig = "Y";
        }
        # Join the columns together and print the results in a TSV file.
        my $final = join($stuff, @columns[$sel+1..scalar(@columns)]);
        $final =~ s/[,|\t]$//;
        print IN join("\t", @columns[0..$sel])."\t$sig\t".$final."\n";
    }
    close(IN);
    close($fh);
    return $count;
}

sub run_cmd() {
    my $cmd = $_[0];
    my $ok = run(@$cmd);
    if (!$ok)
    {
        die "Command failed: @$cmd\n";
    }
}