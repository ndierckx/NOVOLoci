#!/usr/bin/env perl
######################################################
#         SOFTWARE COPYRIGHT NOTICE AGREEMENT        #
#  Copyright (C) {2022-2024}  {Nicolas Dierckxsens}  #
#              All Rights Reserved                   #
#         See file LICENSE for details.              #
######################################################
#           NOVOLoci - Haplotype-aware assembly
#           nicolasdierckxsens@hotmail.com

use strict;

my @modules = ("Getopt::Long", "MCE::Child", "MCE::Channel", "Parallel::ForkManager");

foreach my $module (@modules) 
{
	eval "require $module";
    if ($@) 
	{
        print "Module $module is not installed: $@\n";
		system('perl', '-MCPAN', '-e', "install $module") == 0
        or die "Failed to install $module: $?\n";
    }
}

require Getopt::Long;
Getopt::Long->import();
require MCE::Child;
MCE::Child->import();
require MCE::Channel;
MCE::Channel->import();
require Parallel::ForkManager;
Parallel::ForkManager->import();
require MCE;
MCE->import();

use MCE::Channel;
use MCE::Child;
use Parallel::ForkManager;
use IO::Handle;

my $sequencing_depth_NP = '30';
my $sequencing_depth_PB = '30';
my $y = '1';
my $y0 = '1';
my %y;
undef %y;
my $iterations = "10000000000";
my $overlap = "12";
my $find_haps_in_seed = "";
my $keep_read_ids = "";
my $config = "";
my $project = "";
my $batch_file = "";
my $assembly_length_max = "";
my $type = "";
my $subsample = "";
my $save_reads = "";
my $seed_input0 = "";
my $ploidy = "2";
my $reference = "";
my $variance_detection = "";
my $NP_reads = "";
my $NP_reads_support = "";
my $NP_reads_support_SNR = "";
my $PB_reads = "";
my $output_path = "";
my $reverse_seed = "yes";
my $maxProcs = '5';
my $minimum_read_length_NP = '500';
my $minimum_read_length_PB = '500';

my $input_reads_DB_folder_NP = "";
my $input_reads_DB_folder_PB = "";

my $SNR_read_ahead = "";
my $SNR_read_back_ahead = "";
my $seed_batch = "";
my $seed_input = "";
my $use_quality_scores_NP = "";
my %quality_scores_NP;
undef %quality_scores_NP;

my $read = "";
my $read_back = "";
my $seed_id2 = "";
my $id = "";
my $position = '0';
my $position_back = '0';
my $compare_haps = "";
my $compare_haps_stop = "";
my %no_hap_track;
my $skip_hap = "";
my $nuc_unique_for_ONT = "";
my $reduce_last_600_PB = '1000';
my $reduce_last_600_NP = "";
my $allow_multi_match = "yes";
my %last_non_complex_region;
undef %last_non_complex_region;
my %original_seed_length;
undef %original_seed_length;
my $min_seed_length_PB = '100';

my %filehandle;
my %filehandle3;
my %filehandle4;
my %save_reads;
my %seed;
my %position;
my %position_back;

my %split_positions;
undef %split_positions;
my %split_positions_extra;
undef %split_positions_extra;
my %split_positions_DUP;
undef %split_positions_DUP;
my %split_positions_back;
undef %split_positions_back;
my %quality_scores;
undef %quality_scores;
my %quality_scores_gap;
undef %quality_scores_gap;
my %assembled_reads;
undef %assembled_reads;
my %printed_reads_NP;
undef %printed_reads_NP;
my %printed_reads_PB;
undef %printed_reads_PB;
my %cut_repeat_seq;
undef %cut_repeat_seq;
#my %position_correction;
my %hap_compare_pos;
my $hap_compare_mismatch_extend = '0';
my %PB_split_nucs;
my %PB_split_ids;
my $PB_extension = "";

my @seed_list_sorted;
undef @seed_list_sorted;
my $print_sep = '1';
my $repetitive_detect1 = "";
my $repetitive_detect2 = "";

my $full_reset_NP = "";
my $unresolvable_NP = "";
my $split_contigs_NP = "";
my $split_contigs_PB = "";
my $full_reset_PB = "";
my $unresolvable_PB = "";
my $TMP_directory = "";
my $max_file_count = '9990';
my $retry_NP = "";
my $retry_PB = "";

my %save_alignment_data_NP;
undef %save_alignment_data_NP;
my %rejected_alignment_data_NP;
undef %rejected_alignment_data_NP;
my %rejected_alignment_data_PB;
undef %rejected_alignment_data_PB;
my %save_alignment_data_NP_back;
undef %save_alignment_data_NP_back;
my %save_alignment_data_PB;
undef %save_alignment_data_PB;
my %save_alignment_data_PB_back;
undef %save_alignment_data_PB_back;
my %trace_back_split_NP;
undef %trace_back_split_NP;
my %track_split_NP;
undef %track_split_NP;
my %trace_back_split_PB;
undef %trace_back_split_PB;
my %track_split_PB;
undef %track_split_PB;
my %prev_position_hap_compare;
undef %prev_position_hap_compare;
my %first_back_assembly;
undef %first_back_assembly;

#Read the config file----------------------------------------------------------------------------------------------

GetOptions (
            "c=s" => \$config,
            ) or die "Incorrect usage!\n";

open(CONFIG, $config) or die "Error:Can't open the configuration file, please check the manual!\n\nUsage: perl NOVOLoci1.0.pl -c config.txt\n";

while (my $line = <CONFIG>)
{
    $line =~ tr/\r//d;
    $line =~ s/\R/\012/;
	$line =~ s/[ \t\xA0]+$//;                                                       #CHECK THISSSSSSSSSSSSSSSSSSSSSSSSSSs
    if ($line =~ m/.*Project name\s+\=\s+(.*?)(Assembly length.*)*$/)
    {
        $project = $1;
        chomp $project;
        my $project_tmp = $project;
        my $ggg;
        if ($project =~ m/batch\:(.*)/)
        {
            my $batch_file_tmp = $1;
            if ($batch_file eq "")
            {
                $batch_file = $batch_file_tmp;
                print "Batch file detected...\n\n";
                open(BATCH, $batch_file) or die "Error: $!\nCan't open the batch file, please check the manual!\n";
                $ggg = "yes";
            }
            while (my $line = <BATCH>)
            {
                $project = $line;
                chomp $project;
                last;
            }
            if ($project eq $project_tmp || $project eq "")
            {
                goto EXIT;
            }
            elsif ($ggg ne "yes")
            {
                print "\n\n------------------------------\n------------------------------\n";
                print "        NEXT SAMPLE\n";
                print "------------------------------\n------------------------------\n\n\n";
            }
        }
    }
    if ($line =~ m/.*Assembly length\s+\=\s+(.*?)(Subsample.*)*$/)
    {
        $assembly_length_max = $1;
        chomp $assembly_length_max;
        if ($assembly_length_max eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $assembly_length_max = $line;
                chomp $assembly_length_max;
                last;
            }
        }    
    }     
    if ($line =~ m/.*Subsample\s+\=\s+(.*?)(Save assembled reads.*)*$/)
    {
        $subsample = $1;
        chomp $subsample;
        if ($subsample eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $subsample = $line;
                chomp $subsample;
                last;
            }
        }
    }
    if ($line =~ m/.*Save assembled reads\s+\=\s+(.*?)(Seed Input.*)*$/)
    {
        $save_reads = $1;
        chomp $save_reads;
        if ($save_reads eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $save_reads = $line;
                chomp $save_reads;
                last;
            }
        }
    }
    if ($line =~ m/.*Seed Input\s+\=\s+(.*?)(Ploidy.*)*$/)
    {
        $seed_input0 = $1;
        chomp $seed_input0;
        if ($seed_input0 eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $seed_input0 = $line;
                chomp $seed_input0;
                last;
            }
        }
    }
    if ($line =~ m/.*Ploidy\s+\=\s+(.*?)(Reference sequence.*)*$/)
    {
        $ploidy = $1;
        chomp $ploidy;
        if ($ploidy eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $ploidy = $line;
                chomp $ploidy;
                last;
            }
        }
    }
    if ($line =~ m/.*Reference sequence\s+\=\s+(.*?)(Variance detection.*)*$/)
    {
        $reference = $1;
        chomp $reference;
        if ($reference eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $reference = $line;
                chomp $reference;
                last;
            }
        }
    }
    if ($line =~ m/.*Variance detection\s+\=\s+(.*?)(Cores.*)*$/)
    {
        $variance_detection = $1;
        chomp $variance_detection;
        if ($variance_detection eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $variance_detection = $line;
                chomp $variance_detection;
                last;
            }
        }
    }
    if ($line =~ m/.*Cores\s+\=\s+(.*?)(Output path.*)*$/)
    {
        $maxProcs = $1;
        chomp $maxProcs;
        if ($maxProcs eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $maxProcs = $line;
                chomp $maxProcs;
                last;
            }
        }
    }
    if ($line =~ m/.*Output path\s+\=\s+(.*?)(TMP path.*)*$/)
    {
        $output_path = $1;
        chomp $output_path;
        if ($output_path eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $output_path = $line;
                chomp $output_path;
                last;
            }
        }
    }
	if ($line =~ m/.*TMP path\s+\=\s+(.*?)(Nanopore reads:.*)*$/)
    {
        $TMP_directory = $1;
        chomp $TMP_directory;
        if ($TMP_directory eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $TMP_directory = $line;
                chomp $TMP_directory;
                last;
            }
        }
    }
	
    if ($line =~ m/.*Nanopore reads\s+\=\s+(.*?)(Local DB and NP reads.*)*$/)
    {
        $NP_reads = $1;
        chomp $NP_reads;
        if ($NP_reads eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $NP_reads = $line;
                chomp $NP_reads;
                last;
            }
        }
    }
    if ($line =~ m/.*Local DB and NP reads\s+\=\s+(.*?)(Sequencing depth NP.*)*$/)
    {
        $input_reads_DB_folder_NP = $1;
        chomp $input_reads_DB_folder_NP;
        if ($input_reads_DB_folder_NP eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $input_reads_DB_folder_NP = $line;
                chomp $input_reads_DB_folder_NP;
                last;
            }
        }
    }
    if ($line =~ m/.*Sequencing depth NP\s+\=\s+(.*?)(Min read length NP.*)*$/)
    {
        $sequencing_depth_NP = $1;
        chomp $sequencing_depth_NP;
        if ($sequencing_depth_NP eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $sequencing_depth_NP = $line;
                chomp $sequencing_depth_NP;
                last;
            }
        }    
    }
	if ($line =~ m/.*Min read length NP\s+\=\s+(.*?)(Use Quality scores.*)*$/)
    {
        $minimum_read_length_NP = $1;
        chomp $minimum_read_length_NP;
        if ($minimum_read_length_NP eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $minimum_read_length_NP = $line;
                chomp $minimum_read_length_NP;
                last;
            }
        }    
    }
	if ($line =~ m/.*Use Quality scores\s+\=\s+(.*?)(PacBio reads:.*)*$/)
    {
        $use_quality_scores_NP = $1;
        chomp $sequencing_depth_NP;
        if ($use_quality_scores_NP eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $use_quality_scores_NP = $line;
                chomp $use_quality_scores_NP;
                last;
            }
        }    
    }
    
    if ($line =~ m/.*PacBio reads\s+\=\s+(.*?)(Local DB and PB reads.*)*$/)
    {
        $PB_reads = $1;
        chomp $PB_reads;
        if ($PB_reads eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $PB_reads = $line;
                chomp $PB_reads;
                last;
            }
        }
    }
    if ($line =~ m/.*Local DB and PB reads\s+\=\s+(.*?)(Sequencing depth PB.*)*$/)
    {
        $input_reads_DB_folder_PB = $1;
        chomp $input_reads_DB_folder_PB;
        if ($input_reads_DB_folder_PB eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $input_reads_DB_folder_PB = $line;
                chomp $input_reads_DB_folder_PB;
                last;
            }
        }
    }
    if ($line =~ m/.*Sequencing depth PB\s+\=\s+(.*?)(Min read length PB.*)*$/)
    {
        $sequencing_depth_PB = $1;
        chomp $sequencing_depth_PB;
        if ($sequencing_depth_PB eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $sequencing_depth_PB = $line;
                chomp $sequencing_depth_PB;
                last;
            }
        }    
    }
	if ($line =~ m/.*Min read length PB\s+\=\s+(.*?)$/)
    {
        $minimum_read_length_PB = $1;
        chomp $minimum_read_length_PB;
        if ($minimum_read_length_PB eq "batch" && $batch_file ne "")
        {
            while (my $line = <BATCH>)
            {
                $minimum_read_length_PB = $line;
                chomp $minimum_read_length_PB;
                last;
            }
        }    
    }
}

close CONFIG;

if ($variance_detection eq "no")
{
    $variance_detection = "";
}
$project =~ tr/: /__/;

if ($output_path eq "")
{
    die "\nPlease give an output path in the config file: $output_path, $!\n";
}
my $output_path_test = substr $output_path, -1, 1;
my $output_path_test2 = substr $output_path, 0, 1;

if ($output_path_test ne "\\" && $output_path_test ne "/" && $output_path_test2 ne "\\" && $output_path_test2 ne "/" )
{
    die "\nCan not recognize the output path, it should end with a directory separator (/ or \): $output_path, $!\n";                         
}
elsif ($output_path ne "")
{
    if ($output_path_test ne "\\" && $output_path_test ne "/")
	{
		$output_path .= $output_path_test2;
		$output_path_test = $output_path_test2;
	}
	unless (-d $output_path)
    {
        mkdir $output_path;
    }
}
my $TMP_path_test = substr $TMP_directory, -1 , 1;
my $TMP_path_test2 = substr $TMP_directory, 0 , 1;

if ($TMP_directory eq "")
{
	$TMP_directory = $output_path;
	print "\nNo TMP directory was given, therefore the output path will be used as TMP directory, which could slow down the assembly: $TMP_directory, $!\n";
}
elsif ($TMP_path_test ne "\\" && $TMP_path_test ne "/" && $TMP_path_test2 ne "\\" && $TMP_path_test2 ne "/")
{
    die "\nCan not recognize the TMP path, it should end with a directory separator: $TMP_directory, $!\n";                         
}
elsif ($TMP_directory ne "")
{
    if ($TMP_path_test ne "\\" && $TMP_path_test ne "/")
	{
		$TMP_directory .= $TMP_path_test2;
	}
	unless (-d $TMP_directory)
    {
        mkdir $TMP_directory;
    }
}

my $output_file4  = $output_path."log_".$project.".txt";
open(OUTPUT4, ">" .$output_file4) or die "\nCan't open file $output_file4, $!\n";

print "\n\n-----------------------------------------------";
print "\nNOVOLoci\n";
print "Version 1.0\n";
print "Author: Nicolas Dierckxsens, (c) 2022-2024\n";
print "-----------------------------------------------\n\n";

print "\nInput parameters from the configuration file:   *** Verify if everything is correct ***\n\n";
print "Project:\n";
print "-----------------------\n";
print "Project name          = ".$project."\n";
print "Assembly length       = ".$assembly_length_max."\n";
print "Subsample             = ".$subsample."\n";
print "Save assembled reads  = ".$save_reads."\n";
print "Seed Input            = ".$seed_input0."\n";
print "Ploidy                = ".$ploidy."\n";
print "Reference sequence    = ".$reference."\n";
print "Variance detection    = ".$variance_detection."\n";
print "Cores                 = ".$maxProcs."\n";
print "Output path           = ".$output_path."\n";
print "TMP path              = ".$TMP_directory."\n\n";

print "Nanopore reads:\n";
print "-----------------------\n";
print "Nanopore reads        = ".$NP_reads."\n";
print "Local DB and NP reads = ".$input_reads_DB_folder_NP."\n";
print "Sequencing depth NP   = ".$sequencing_depth_NP."\n";
print "Min read length NP    = ".$minimum_read_length_NP."\n";
print "Use Quality scores    = ".$use_quality_scores_NP."\n\n";

print "PacBio reads:\n";
print "-----------------------\n";
print "PacBio reads          = ".$PB_reads."\n";
print "Local DB and PB reads = ".$input_reads_DB_folder_PB."\n";
print "Sequencing depth PB   = ".$sequencing_depth_PB."\n";
print "Min read length PB    = ".$minimum_read_length_PB."\n\n";



print OUTPUT4 "\n\n-----------------------------------------------";
print OUTPUT4 "\nNOVOLoci\n";
print OUTPUT4 "Version 1.0\n";
print OUTPUT4 "Author: Nicolas Dierckxsens, (c) 2022-2024\n";
print OUTPUT4 "-----------------------------------------------\n\n";

print OUTPUT4 "\nInput parameters from the configuration file:   *** Verify if everything is correct ***\n\n";
print OUTPUT4 "Project:\n";
print OUTPUT4 "-----------------------\n";
print OUTPUT4 "Project name          = ".$project."\n";
print OUTPUT4 "Assembly length       = ".$assembly_length_max."\n";
print OUTPUT4 "Subsample             = ".$subsample."\n";
print OUTPUT4 "Save assembled reads  = ".$save_reads."\n";
print OUTPUT4 "Seed Input            = ".$seed_input0."\n";
print OUTPUT4 "Ploidy                = ".$ploidy."\n";
print OUTPUT4 "Reference sequence    = ".$reference."\n";
print OUTPUT4 "Variance detection    = ".$variance_detection."\n";
print OUTPUT4 "Cores                 = ".$maxProcs."\n";
print OUTPUT4 "Output path           = ".$output_path."\n";
print OUTPUT4 "TMP path              = ".$TMP_directory."\n\n";

print OUTPUT4 "Nanopore reads:\n";
print OUTPUT4 "-----------------------\n";
print OUTPUT4 "Nanopore reads        = ".$NP_reads."\n";
print OUTPUT4 "Local DB and NP reads = ".$input_reads_DB_folder_NP."\n";
print OUTPUT4 "Sequencing depth NP   = ".$sequencing_depth_NP."\n";
print OUTPUT4 "Min read length NP    = ".$minimum_read_length_NP."\n";
print OUTPUT4 "Use Quality scores    = ".$use_quality_scores_NP."\n\n";

print OUTPUT4 "PacBio reads:\n";
print OUTPUT4 "-----------------------\n";
print OUTPUT4 "PacBio reads          = ".$PB_reads."\n";
print OUTPUT4 "Local DB and PB reads = ".$input_reads_DB_folder_PB."\n";
print OUTPUT4 "Sequencing depth PB   = ".$sequencing_depth_PB."\n";
print OUTPUT4 "Min read length PB    = ".$minimum_read_length_PB."\n\n";


#Warning messages-----------------------------------------------------------------------------------------------------------------------------

if (($sequencing_depth_NP < 1 || $sequencing_depth_NP eq "") && ($NP_reads ne "" || $input_reads_DB_folder_NP ne ""))
{
	die "\n'Please give an estimation of the Nanopore sequencing depth in the config file!\n";
}
if (($sequencing_depth_PB < 1 || $sequencing_depth_PB eq "") && ($PB_reads ne "" || $input_reads_DB_folder_PB ne ""))
{
	die "\n'Please give an estimation of the PacBio sequencing depth in the config file!\n";
}
if ($ploidy ne "1" && $ploidy ne "2")
{
    die "\n'The Ploidy option has to be '1' or '2'! Polyploid (>2) assemblies will be supported in the future\n";
}

$sequencing_depth_NP *= 0.7;
$sequencing_depth_PB *= 0.7;
$sequencing_depth_NP /= $ploidy;
$sequencing_depth_PB /= $ploidy;

if ($variance_detection eq "yes" && $reference eq "")
{
    #die "\nWhen variance detection is on, you must give a reference sequence, please check the configuration file!\n";
}
if ($save_reads ne "yes" && $save_reads ne "1" && $save_reads ne "2" && $save_reads ne ""  && $save_reads ne "no")
{
    die "\n'Save assembled reads' has to be '1', '2' or empty, please check the configuration file!\n";
}
if ($subsample =~ m/^\d*\.*\d*$/)
{
}
else
{
    die "\n'Subsample' can only be an integer!\n";
}
if ($save_reads eq "no")
{
    $save_reads = "";
}

my $USAGE = "\nUsage: perl NOVOLoci1.0.pl -c config_example.txt";

if (($NP_reads ne "" || $input_reads_DB_folder_NP ne "") && ($PB_reads ne "" || $input_reads_DB_folder_PB ne ""))
{
    $NP_reads_support = "yes";
}

chomp($maxProcs);
if ($maxProcs > 0)
{}
else
{
    die "\nMax Threads: '$maxProcs' has to be a value above 1!\n";
}

my $tmp_sequences_directory_NP = $output_path."tmp_sequences_NP";
my $tmp_sequences_directory_PB = $output_path."tmp_sequences_PB";

if ($input_reads_DB_folder_NP eq "" && $NP_reads ne "")
{
    mkdir $tmp_sequences_directory_NP;
}
if ($input_reads_DB_folder_PB eq "" && $PB_reads ne "")
{
    mkdir $tmp_sequences_directory_PB;
}

if ($input_reads_DB_folder_NP ne "")
{
    $tmp_sequences_directory_NP = $input_reads_DB_folder_NP."tmp_sequences_NP";
}
if ($input_reads_DB_folder_PB ne "")
{
    $tmp_sequences_directory_PB = $input_reads_DB_folder_PB."tmp_sequences_PB";
}

my $directory_DB_NP = $output_path."BLAST_DB_NP";
my $directory_DB_PB = $output_path."BLAST_DB_PB";

if ($input_reads_DB_folder_NP eq "" && $NP_reads ne "")
{
    mkdir $directory_DB_NP;
}
else
{
    $directory_DB_NP = $input_reads_DB_folder_NP."BLAST_DB_NP";
}

if ($input_reads_DB_folder_PB eq "" && $PB_reads ne "")
{
    mkdir $directory_DB_PB;
}
else
{
    $directory_DB_PB = $input_reads_DB_folder_PB."BLAST_DB_PB";
}

if ($use_quality_scores_NP ne "" && $use_quality_scores_NP ne "yes")
{
	open(INPUT7a, $use_quality_scores_NP) or die "\nCan't open quality score file, if you have a quality score file, you should give the path to that file,
	if it is your first run on this dataset, you should write \"yes\". If you don't want to use quality scores, leave it blank.\n";
	close INPUT7a;
}

if ($minimum_read_length_NP eq "")
{
	$minimum_read_length_NP = '1000';
}
if ($minimum_read_length_PB eq "")
{
	$minimum_read_length_PB = '500';
}

if ($PB_reads ne "" || $input_reads_DB_folder_PB ne "")
{
	die "\n'This is a beta version that only works for Nanopore, Pacbio and Hybrid assembly will be included with the next update!\n";
}

sub AT_rich_test
{
    my @str = @_;
    my $region_to_check = $str[0];
    $region_to_check =~ tr/N|K|R|Y|S|W|M|B|D|H|V/\./;
    my $extra = $str[1];
    my $AT_rich = "";

    my $A_rich_test = $region_to_check =~ tr/A/A/;
    my $T_rich_test = $region_to_check =~ tr/T/T/;
    my $G_rich_test = $region_to_check =~ tr/G/G/;
    my $C_rich_test = $region_to_check =~ tr/C/C/;
    my $dot_rich_test3 = $region_to_check =~ tr/\./\./;
    
    if ($A_rich_test+$dot_rich_test3 > length($region_to_check)-$extra || $T_rich_test+$dot_rich_test3 > length($region_to_check)-$extra
        || $G_rich_test+$dot_rich_test3 > length($region_to_check)-$extra || $C_rich_test+$dot_rich_test3 > length($region_to_check)-$extra
        || $A_rich_test+$T_rich_test > length($region_to_check)-$extra || $A_rich_test+$C_rich_test > length($region_to_check)-$extra 
        || $A_rich_test+$G_rich_test > length($region_to_check)-$extra || $C_rich_test+$T_rich_test > length($region_to_check)-$extra 
        || $C_rich_test+$G_rich_test > length($region_to_check)-$extra || $T_rich_test+$G_rich_test > length($region_to_check)-$extra)
    {
        $AT_rich = "yes";
    }
    return $AT_rich;
}
#Check in the extensions if an SNR is ahead---------------------------------------------------------------------------------------------
sub SNR_ahead
{
    my @str = @_;
    my %extensions_tmp = %{$str[0]};
    my $skip = $str[1];
    my $back = $str[2];
    my $SNR_ahead = '0';
    my %total;
    undef %total;
AHEAD:foreach my $extension_id_tmp (keys %extensions_tmp)
    {
        my $extension_tmp = $extensions_tmp{$extension_id_tmp};
        if (length($extension_tmp) > 10)
        {
            my $first11 = substr $extension_tmp, $skip, 12;
            #$first11 =~ tr/1234/ACTG/;
            my $A = $first11 =~ tr/1A/1A/;
            my $C = $first11 =~ tr/2C/2C/;
            my $T = $first11 =~ tr/3T/3T/;
            my $G = $first11 =~ tr/4G/4G/;
            my $N = $first11 =~ tr/N/N/;
            if ($A > 7 || $C > 7 || $T > 7 || $G > 7)
            {
                $SNR_ahead++;
                $total{'A'} = $A+$total{'A'};
                $total{'C'} = $C+$total{'C'};
                $total{'T'} = $T+$total{'T'};
                $total{'G'} = $A+$total{'G'};
            }
            elsif ($N < 2)
            {
                $SNR_ahead--;
            }
        }
        if ($SNR_ahead > 10 || $SNR_ahead < -2)
        {
            last AHEAD;
        }
    }
    if ($SNR_ahead > 7)
    {
        $SNR_ahead = "yes";
        my $highest = '0';
        foreach my $total (keys %total)
        {
            if ($total{$total} > $highest)
            {
                $highest = $total{$total};
                if ($back eq "yes")
                {
                    $SNR_read_back_ahead = $total;
                }
                else
                {
                    $SNR_read_ahead = $total;
                }        
            }  
        }
    }
    else
    {
        $SNR_ahead = "";
    }
    return $SNR_ahead;
}

my $output_file5  = $output_path."log_extended_".$project.".txt";
my $output_file10 = $output_path."Assembled_reads_PB_".$project.".fasta";
my $output_file11 = $output_path."Assembled_reads_NP_".$project.".fasta";
my $output_file12 = $output_path."structural_variation_".$project.".vcf";
my $output_file13 = $output_path."quality_scores_".$project.".txt";

if ($seed_input0 ne "")
{
    open(INPUT3, $seed_input0) or die "Can't open the seed file, $!\n";
}

if ($save_reads ne "" && ($PB_reads ne "" || $input_reads_DB_folder_PB ne ""))
{
    open(OUTPUT10, ">".$output_file10) or die "Can't open saved reads file $output_file10, $!\n";
}
if ($save_reads ne "" && ($NP_reads ne "" || $input_reads_DB_folder_NP ne ""))
{
    open(OUTPUT11, ">".$output_file11) or die "Can't open saved reads file $output_file11, $!\n";
}
if ($variance_detection eq "yes")
{
    open(OUTPUT12, ">" .$output_file12) or die "Can't open variation file $output_file12, $!\n";
}

#Open variance vcf file----------------------------------------------------------------------------------------------------------------------------------
if ($variance_detection eq "yes")
{
    my ($wday, $mon, $mday, $hour, $min, $sec, $year) = localtime;
    my @localtime = split / /, localtime;
    my %mon2num = qw(
    Jan 01  Feb 02  Mar 03  Apr 04  May 05  Jun 06
    Jul 07  Aug 08  Sep 09  Oct 10 Nov 11 Dec 12
    );
    my $month = $localtime[1];
    if (exists($mon2num{$localtime[1]}))
    {
       $month = $mon2num{$localtime[1]};
    }
    print OUTPUT12 "##fileformat=VCFv4.0\n";
    print OUTPUT12 "##fileDate=".$localtime[4].$month.$localtime[2]."\n";
    print OUTPUT12 "##reference=".$reference."\n";
    print OUTPUT12 "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n";
    print OUTPUT12 "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">\n";
    #print OUTPUT12 "##INFO=<ID=FR,Number=1,Type=Flag,Description=\"Detected on the forward(F) and/or reverse(R) strand\">\n";
    print OUTPUT12 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}

#Make hash of reference----------------------------------------------------------------------------------------------------------------------------------

my %hashref;
my %hashref2;
undef %hashref;
undef %hashref2;

if ($reference ne "")
{
    select(STDERR);
    $| = 1;
    select(STDOUT); # default
    $| = 1;
    print "\nScan reference sequence...";
    open(INPUT5, $reference) or die "\n\nCan't open reference file $reference, $!\n";
    my $ff2 = '0';
    my $value_ref2 = "";
    
    while (my $line = <INPUT5>)
    {
        if ($ff2 < 1)
        {
            $ff2++;
            next;
        }
        chomp $line;    
        $line =~ tr/actgn/ACTGN/;
        my $first = substr $line, 0, 1;
        my $line3;
        if ($first eq '>' || $first eq '@')
        {
            $line3 = $value_ref2."NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        }
        else
        {
            $line3 = $value_ref2.$line;
        }
        
        while (length($line3) > ((30*3)-1))
        {
            my $value_ref2b = substr $line3, 0, 30;
            my $line2 = $line3;
            $line3 = substr $line2, 1;
            
            if (exists($hashref{$value_ref2b}))
            {
                $hashref{$value_ref2b} .= ",".$ff2;
            }
            else
            {
                $hashref{$value_ref2b} = $ff2;
            }
            if (exists($hashref2{$ff2}))
            {
                $hashref2{$ff2} .= ",".$value_ref2b;
            }
            else
            {
                $hashref2{$ff2} = $value_ref2b;
            }
            $ff2++;
        }
        $value_ref2 = $line3;
    }
    while (length($value_ref2) > 1)
    {
        my $value_ref2b = substr $value_ref2, 0, 30;
        my $value_ref2bc = $value_ref2;
        $value_ref2 = substr $value_ref2bc, 1;
        
        if (exists($hashref{$value_ref2b}))
        {
            $hashref{$value_ref2b} .= ",".$ff2;
        }
        else
        {
            $hashref{$value_ref2b} = $ff2;
        }
        if (exists($hashref2{$ff2}))
        {
            $hashref2{$ff2} .= ",".$value_ref2b;
        }
        else
        {
            $hashref2{$ff2} = $value_ref2b;
        }
        $ff2++;
    }
    close INPUT5;
    print "...OK\n";
}

my %delete_input_files;
undef %delete_input_files;
my %original_input_files;
undef %original_input_files;
my @QUALITY_HASH;

my $maxProcs_tmp = $maxProcs-1;
if ($maxProcs_tmp < 1)
{
	$maxProcs_tmp = '1';
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------
#Make hash of ONT reads----------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------

my $read_count_NP = '0';
my $total_read_length_NP = '0';
my $max_read_length_NP = '0';
my $min_read_length_NP = '0';
my $removed_reads_min_NP = '0';
my $N50_NP = '0';
my @read_lengths_NP;

my $pm = new Parallel::ForkManager($maxProcs_tmp);
my %ret_data_NP;
undef %ret_data_NP;
my %ret_data2_NP;
undef %ret_data2_NP;

$pm->run_on_finish(
sub
{ 
    my @str = @_;
    my %hash_tmp = %{$str[5]};
    my $pid = $str[0];

    # retrieve data structure from child
    if (defined($pid))
    {  # children are not forced to send anything
        my $string = $hash_tmp{'1'};  # child passed a string reference
        my $string2 = $hash_tmp{'2'};
        my $string3 = $hash_tmp{'3'};
        my $string4 = $hash_tmp{'4'};
        my $string5 = $hash_tmp{'5'};
        my $string6 = $hash_tmp{'6'};
        my $string7 = $hash_tmp{'7'};
        my @string8 = @{$hash_tmp{'8'}};

        $ret_data_NP{$pid} = $string;
        $ret_data2_NP{$pid} = $string3;
        $read_count_NP += $string2;
        $total_read_length_NP += $string4;
        if ($string6 > $max_read_length_NP)
        {
            $max_read_length_NP = $string6;
        }
        if ($string5 < $min_read_length_NP || $min_read_length_NP eq '0')
        {
            $min_read_length_NP = $string5;
        }
        $removed_reads_min_NP += $string7;
        @read_lengths_NP = (@read_lengths_NP, @string8);
    }
    #else {  # problems occurring during storage or retrieval will throw a warning
  #print qq|No message received from child process $data_structure_reference!\n|;
#}
});

my $long_id_NP = '0';

if ($NP_reads ne "" || $input_reads_DB_folder_NP ne "")
{   
	my $NP_check_directory = substr $NP_reads, -1;
	if ($input_reads_DB_folder_NP eq "" && $NP_reads ne "")
	{
		if ($NP_check_directory eq $output_path_test)
		{
			$NP_check_directory = "yes";
		}
		elsif (-d $NP_reads)
		{
			$NP_reads .= $output_path_test;
			$NP_check_directory = "yes";
		}
		elsif (-e $NP_reads)
		{
		}
		else
		{
			die "\nThe Nanopore input doesn't not seems to be an existing file or folder!: ".$NP_reads."\n";
		}
	}
	
	select(STDERR);
    $| = 1;
    select(STDOUT); # default
    $| = 1;
    print "\nBuilding local databases for the Nanopore reads...";
    
    
    my %NP_reads;
    undef %NP_reads;
    my $count_files = '1';
    if ($NP_reads ne "")
    {
        my $disered_db_count = int($maxProcs_tmp/2);
        if ($NP_check_directory eq "yes")
        {
            opendir(DIR_NP, $NP_reads) or die "Could not open $NP_reads\n";
            
            for my $filename (sort readdir(DIR_NP))
            {
                my $last5 = substr $filename, -5;
                if ($last5 eq "fasta" || $last5 eq "fastq" || $last5 eq "tq.gz" || $last5 eq "ta.gz" || $last5 eq "q.bz2" || $last5 eq "a.bz2")
                {
                    my $NP_reads_tmp = $NP_reads.$filename;
                    $NP_reads{$NP_reads_tmp} = $count_files."a";
                    $original_input_files{$NP_reads_tmp} = undef;
                    $count_files++;
                }
                elsif ($filename ne "." && $filename ne "..")
                {
                    print "\n".$filename.": File extension not recognized!\n";
                    print OUTPUT4 "\n".$filename.": File extension not recognized!\n";
                }
            }
            closedir DIR_NP;
        }
        else
        {
            $NP_reads{$NP_reads} = "1a";
            $original_input_files{$NP_reads} = undef;
            $count_files++;
        }

        #my $total_filesize = '0';
        my $total_lines = '0';
        #my %file_sizes;
        #undef %file_sizes;
        my %file_lines;
        undef %file_lines;    
        
        foreach my $NP_reads_tmp (keys %NP_reads)
        {
            #my $filesize = "";
            my $check_zip_long = substr $NP_reads_tmp, -2;
            my $check_zip_long2 = substr $NP_reads_tmp, -3;           
            
            if ($check_zip_long eq "gz")
            {
                #my $FILE_LONG;
                #open $FILE_LONG, '<:raw', $NP_reads_tmp or die "Could not open file '$NP_reads_tmp': $!";
                #seek $FILE_LONG, -4, 2;  # Seek to four bytes from the end of the file
                
                #my $buf;
                #read $FILE_LONG, $buf, 4;  # Read the last four bytes
                #my $uncompressed_size = unpack 'V', $buf;  # Unpack the four bytes as a little-endian unsigned integer
                #$total_filesize += $uncompressed_size;    
                #close $FILE_LONG;
                
                my $new_filename = substr $NP_reads_tmp, 0, -3;
                my @new_filename = split /$output_path_test/, $new_filename;
                my $g = @new_filename;
                $g--;
                $new_filename = $output_path.$new_filename[$g];
                my $return_value = system("gzip -c -q -k -d ".$NP_reads_tmp." > ".$new_filename);
                
                $delete_input_files{$new_filename} = undef;
                
                #$file_sizes{$new_filename}{'2'} = $uncompressed_size;
                $NP_reads{$new_filename} = $NP_reads{$NP_reads_tmp};
                delete $NP_reads{$NP_reads_tmp};
                my $count_lines = qx(wc -l $new_filename);
                chomp($count_lines);
                my @count_lines = split /\s+/, $count_lines;
                $total_lines += $count_lines[0];
                $file_lines{$new_filename} = $count_lines[0];
            }
            elsif ($check_zip_long2 eq "bz2")
            {
                die "Can't read bz2 files, pleas decompress: $NP_reads_tmp, $!\n";
            }
            elsif ($check_zip_long2 eq "zip")
            {
                die "Can't read zip files, pleas decompress: $NP_reads_tmp, $!\n";
            }
            else
            {
                #$filesize = -s $NP_reads_tmp;
                #$total_filesize += $filesize;
                #$file_sizes{$NP_reads_tmp}{'2'} = $filesize;
                my $count_lines = qx(wc -l $NP_reads_tmp);
                chomp($count_lines);
                my @count_lines = split /\s+/, $count_lines;
                $total_lines += $count_lines[0];
                $file_lines{$NP_reads_tmp} = $count_lines[0];
            }        
        }
        #my $bytes_perDB = int($total_filesize/int(($maxProcs-1)/2))+100000;
        my $lines_perDB = int($total_lines/int($maxProcs_tmp/2))+8;
        my $adjusted_lines_perDB = (int($lines_perDB / 4) * 4)+4;
		                                          #print "\n".$lines_perDB.": ".$adjusted_lines_perDB." adjusted lines\n";
        
		my $NP_reads_DB = "";
NP_READS_DB: 
		
		my $count_total_files = keys %file_lines;
        my $count_files_tmp = '0';
        my $merge_files_check = "";
        my $merge_files_check_first = "";
        my %file_lines_merged;
        undef %file_lines_merged;
        my %NP_reads_tmp;
        undef %NP_reads_tmp;
     
        foreach my $NP_reads_tmp (sort keys %file_lines)
        {
            $count_files_tmp++;
#print $file_lines{$NP_reads_tmp}." FILE_LINES\n";
#print $count_files_tmp." COUNT_LINES\n";
#print $count_total_files." COUNT_LINES_TOTAL\n";
            if ($file_lines{$NP_reads_tmp} > $adjusted_lines_perDB*1.15)
            {
                my $command_split = "split --additional-suffix=".$NP_reads{$NP_reads_tmp}.".fasta -l ".$adjusted_lines_perDB." ".$NP_reads_tmp." ".$directory_DB_NP.$output_path_test;                   
                my $return_value = system($command_split);
                
                my $expected_files = $file_lines{$NP_reads_tmp}/$adjusted_lines_perDB;
                my $expected_files0 = $file_lines{$NP_reads_tmp}/$adjusted_lines_perDB;
                if ($expected_files0 > int($expected_files0))
                {
                    $expected_files = int($expected_files0)+1;
                }
                my $found_files = 0;

                while ($found_files < $expected_files)
                {
                    $found_files = 0;
                    
                    # Check for the existence of the expected files
                    for my $file (glob("$directory_DB_NP$output_path_test*"))
                    {
                        if (-s $file)
						{
							$found_files++;
						}
                    }           
                    # Sleep for a bit if not all files have been found
                    sleep(10) if $found_files < $expected_files;
                }
                delete $NP_reads{$NP_reads_tmp};
				delete $file_lines{$NP_reads_tmp};
            }
            elsif ($merge_files_check ne "")
            {
                my $combined_lines = $file_lines{$NP_reads_tmp} + $file_lines_merged{$merge_files_check};

                if ($combined_lines < $adjusted_lines_perDB*1.15)
                {  
                    my $output_tmp = $directory_DB_NP.$output_path_test.$count_files_tmp."m".$merge_files_check_first.".fasta";
                    my $command_merge = "cat ".$merge_files_check." ".$NP_reads_tmp." > ".$output_tmp;                   
                    my $return_value = system($command_merge);
                    delete $NP_reads{$NP_reads_tmp};
                    delete $NP_reads{$merge_files_check};
					delete $file_lines{$NP_reads_tmp};

                    my $time_before_check_tmp = time;
COUNT_LINES:                   
                    my $count_lines = qx(wc -l $output_tmp);
                    chomp($count_lines);
                    my @count_lines = split /\s+/, $count_lines;
                    my $count_lines2 = $count_lines[0];

                    if ($count_lines2 eq $combined_lines)
                    {      
                        unless (exists($original_input_files{$merge_files_check}))
                        {
                            unlink $merge_files_check;
                        }
                        unless (exists($original_input_files{$NP_reads_tmp}))
                        {
                            unlink $NP_reads_tmp;
                        } 
                    }
                    elsif (time > $time_before_check_tmp+600)
                    {
                        print "ERROR: Can't merge input files\n";
                        goto END1;
                    }
                    else
                    {
                        goto COUNT_LINES;
                    }
                    
                    $NP_reads_tmp{$output_tmp} = $count_files_tmp."m".$merge_files_check_first;
                    $merge_files_check = $output_tmp;
                    $file_lines_merged{$merge_files_check} = $combined_lines;
                    my $output_file_DB_tmp = $directory_DB_NP.$output_path_test.$NP_reads_tmp{$output_tmp}.".fasta"; 
                    $delete_input_files{$output_file_DB_tmp} = undef;
                    $delete_input_files{$output_tmp} = undef;               
                }
                elsif ($file_lines_merged{$merge_files_check} > $adjusted_lines_perDB/1.4)
                { 
                    $NP_reads{$merge_files_check} = $NP_reads_tmp{$merge_files_check}."a";                  
                    my $output_file_DB_tmp = $directory_DB_NP.$output_path_test.$NP_reads{$merge_files_check}.".fasta"; 
                    $delete_input_files{$output_file_DB_tmp} = undef;
                    $merge_files_check = "";
                    $merge_files_check_first = "";
                }
            } 
            if ($file_lines{$NP_reads_tmp} < $adjusted_lines_perDB/1.3 && $count_files_tmp < $count_total_files && $count_total_files > int($maxProcs_tmp/2) && $merge_files_check eq "")
            {
                $merge_files_check = $NP_reads_tmp;
                $merge_files_check_first = substr $NP_reads{$NP_reads_tmp}, 0, -1;
                $NP_reads_tmp{$NP_reads_tmp} = substr $NP_reads{$NP_reads_tmp}, 0, -1;
                $file_lines_merged{$NP_reads_tmp} = $file_lines{$NP_reads_tmp};
                next;
            }
            my $output_file_DB_tmp = $directory_DB_NP.$output_path_test.$NP_reads{$NP_reads_tmp}.".fasta"; 
            $delete_input_files{$output_file_DB_tmp} = undef;         
        }
        if ($merge_files_check ne "")
        {
            $NP_reads{$merge_files_check} = $NP_reads_tmp{$merge_files_check}."a";
            $delete_input_files{$merge_files_check} = undef;
            my $output_file_DB_tmp = $directory_DB_NP.$output_path_test.$NP_reads{$merge_files_check}.".fasta"; 
            $delete_input_files{$output_file_DB_tmp} = undef;
        }
		
		opendir(DIR_NP2, $directory_DB_NP.$output_path_test) or die "Could not open $directory_DB_NP.$output_path_test\n";
		my @files = sort readdir(DIR_NP2);
		foreach my $filename (@files)
		{
			my $last5 = substr $filename, -5;
			if ($last5 eq "fasta" || $last5 eq "fastq")
			{                       
				#my $FILE_tmp3;
				#open ($FILE_tmp3, '<:raw', $NP_reads_tmp) or die "\n\nCan't open long reads file $NP_reads_tmp, $!\n";
				#close $FILE_tmp3;
				my $NP_reads_tmp2 = $directory_DB_NP.$output_path_test.$filename;
				#substr $filename, -6, 6, "";
				$NP_reads{$NP_reads_tmp2} = $count_files."a";
				$count_files++;
				
				$delete_input_files{$NP_reads_tmp2} = undef;
				my $output_file_DB_tmp = $directory_DB_NP.$output_path_test.$NP_reads{$NP_reads_tmp2}.".fasta"; 
				$delete_input_files{$output_file_DB_tmp} = undef;
				
				my $count_lines = qx(wc -l $NP_reads_tmp2);
                chomp($count_lines);
                my @count_lines = split /\s+/, $count_lines;
                $total_lines += $count_lines[0];
                $file_lines{$NP_reads_tmp2} = $count_lines[0];
			}
			elsif ($filename ne "." && $filename ne "..")
			{
				print "\n".$filename.": File extension not recognized!\n";
				print OUTPUT4 "\n".$filename.": File extension not recognized!\n";
			}   
		}
		close DIR_NP2;
				
		if ($use_quality_scores_NP eq "yes")
        {
			foreach my $NP_reads_tmp (keys %NP_reads)
			{
				my $output_file16_tmp  = $output_path."QUALITY_SCORES_TMP_".$NP_reads{$NP_reads_tmp}.".txt";
				push @QUALITY_HASH, $output_file16_tmp;
			}
		}
		
		my $file_count = keys %NP_reads;
		
		if ($file_count > int($maxProcs_tmp/2) && $NP_reads_DB eq "")
		{
			$NP_reads_DB = "yes";	
			goto NP_READS_DB;
		}

#Build local BLAST databse and save reads to disk---------------------------------------------------------------------------
    
        foreach my $NP_reads_tmp (keys %NP_reads)
        {
            my $pid;
            if ($file_count > 1)
            {
                $pid = $pm->start and next;
            }
       srand();     
            my $check_zip_long = substr $NP_reads_tmp, -2;
            my $check_zip_long2 = substr $NP_reads_tmp, -3;
            my $firstLine_long;
            my $secondLine_long;
            my $thirdLine_long;        

            my $read_count_tmp = '0';
            my $total_read_length_tmp = '0';
            my $max_read_length_tmp = '0';
            my $min_read_length_tmp = '0';
            my $removed_reads_min_tmp = '0';
            my @read_lengths_tmp;
            
            my $FILE_tmp;
            if ($check_zip_long eq "gz")
            {
                open ($FILE_tmp, '-|', 'gzip', '-dc', $NP_reads_tmp) or die "Can't open file $NP_reads_tmp, $!\n";
            }
            elsif ($check_zip_long2 eq "bz2")
            {
                open ($FILE_tmp, '-|', 'bzip2', '-dc', $NP_reads_tmp) or die "Can't open file $NP_reads_tmp, $!\n";
            }
            else
            {
                open ($FILE_tmp, $NP_reads_tmp) or die "\n\nCan't open long reads file $NP_reads_tmp, $!\n";
            }
            $firstLine_long = <$FILE_tmp>;
            chomp $firstLine_long;
            $secondLine_long = <$FILE_tmp>;
            chomp $secondLine_long;
            $thirdLine_long = <$FILE_tmp>;
            chomp $thirdLine_long;
            close $FILE_tmp;
         
            my $no_quality_score_long_tmp = substr $thirdLine_long, 0, 1;
            my $quality_score_long = "";
            if ($thirdLine_long eq "+")
            {
                $quality_score_long = "yes";
            }
            
            my $read_limit = "10000000000000000000000000000000000000000000000000000";
            if ($subsample ne "")
            {
                my $command = "wc -l ".$NP_reads_tmp;
                my $lines = `$command`;
                chomp($lines);
                my @line_count = split / /, $lines;
                $read_limit = $subsample*($line_count[0]/4);
                if ($quality_score_long eq "")
                {
                    $read_limit = $subsample*($line_count[0]/2);
                }
            }
    
            my $FILE_LONG;
            if ($check_zip_long eq "gz")
            {
                open ($FILE_LONG, '-|', 'gzip', '-dc', $NP_reads_tmp) or die "Can't open file $NP_reads_tmp, $!\n";
            }
            elsif ($check_zip_long2 eq "bz2")
            {
                open ($FILE_LONG, '-|', 'bzip2', '-dc', $NP_reads_tmp) or die "Can't open file $NP_reads_tmp, $!\n";
            }
            else
            {
                open($FILE_LONG, $NP_reads_tmp) or die "\n\nCan't open long reads file $NP_reads_tmp, $!\n";
            }

            my $directory_tmp = $output_path."tmp_sequences_NP".$output_path_test.$NP_reads{$NP_reads_tmp};
            mkdir $directory_tmp;
            my $directory_DB_tmp = $directory_DB_NP.$output_path_test.$NP_reads{$NP_reads_tmp}.$output_path_test;
            mkdir $directory_DB_tmp;
            
            my $output_file_DB_tmp = $directory_DB_NP.$output_path_test.$NP_reads{$NP_reads_tmp}.".fasta";
            open(OUTPUT_DB_NP1, ">" .$output_file_DB_tmp) or die "\nCan't open file $output_file_DB_tmp, $!\n";
			
			my $FILE_HASH_TMP;
            
            if ($use_quality_scores_NP eq "yes")
            {
                my $output_file16_tmp  = $output_path."QUALITY_SCORES_TMP_".$NP_reads{$NP_reads_tmp}.".txt";
                open($FILE_HASH_TMP, ">".$output_file16_tmp) or die "Can't open file $output_file16_tmp, $!\n";
            }
            
            my $ww = '1';
            my $fail_or_pass = "";
            my $fail_or_pass2 = "";
            my $id_read = "";
			my $id_tmpi = "";
			my $removed_check = "";
            my %hash_NP_reads_kmer_tmp;
            undef %hash_NP_reads_kmer_tmp;
			my $count_files_in_folder = '0';
			my $subfolder_name = '1';
			my $directory_DB_sub_tmp = $output_path."tmp_sequences_NP".$output_path_test.$NP_reads{$NP_reads_tmp}.$output_path_test.$subfolder_name;
            mkdir $directory_DB_sub_tmp;
            
FILE_LONG:  while (my $line = <$FILE_LONG>)
            {
                if ($long_id_NP > $read_limit)
                {
                    last FILE_LONG;
                }
                chomp $line;    
                if ($fail_or_pass eq "fastq_fail" && $fail_or_pass2 > 0)
                {
                    $fail_or_pass2--;
                    next;
                }
                elsif ($ww eq '1')
                {
                    $fail_or_pass = substr $line, 0, 10;
                    if ($fail_or_pass eq "fastq_fail")
                    {
                        if ($quality_score_long eq "")
                        {
                            $fail_or_pass2 = '1';
                        }
                        else
                        {
                            $fail_or_pass2 = '3';
                        }
                        next;
                    }
                    if ($keep_read_ids eq "yes")
                    {
                        my @read_id = split / /, $line;
                        $id_read = $read_id[0];
                    }
                        
                    $ww++;
                    next;
                }
                elsif ($ww eq '2')
                {
                    #$line =~ tr/actgn/ACTGN/;
                    my $pos_tmp = '0';
                    $long_id_NP++;
                    
                    $id_tmpi = $NP_reads{$NP_reads_tmp}.$long_id_NP;
                    if ($keep_read_ids eq "yes")
                    {
                        $id_tmpi = $NP_reads{$NP_reads_tmp}.$id_read;
                    }
                    my $length_tmp = length($line);   
                    $total_read_length_tmp += $length_tmp;
                    
                    if ($length_tmp > $max_read_length_tmp)
                    {
                        $max_read_length_tmp = $length_tmp;
                    }
                    if ($length_tmp < $min_read_length_tmp || $min_read_length_tmp eq '0')
                    {
                        $min_read_length_tmp = $length_tmp;
                    }
                    $read_count_tmp++;
                    
                    push @read_lengths_tmp, $length_tmp;
                    
                    if ($length_tmp >= $minimum_read_length_NP)
                    {
                        #substr $line, 0, 40, "";
                        $count_files_in_folder++;
						if ($count_files_in_folder > $max_file_count)
						{
							$count_files_in_folder = '1';
							$subfolder_name++;
							$directory_DB_sub_tmp = $output_path."tmp_sequences_NP".$output_path_test.$NP_reads{$NP_reads_tmp}.$output_path_test.$subfolder_name;
							mkdir $directory_DB_sub_tmp;
						}
						$id_tmpi .= "a".$subfolder_name;
                        my $output_file_NP1  = $directory_tmp.$output_path_test.$subfolder_name.$output_path_test."sequence_tmp_NP_".$id_tmpi.".fasta";
                        unless (-e $output_file_NP1)
                        {
                            open(OUTPUT_LONG_NP1, ">" .$output_file_NP1) or die "\nCan't open file $output_file_NP1, $!\n";
                            print OUTPUT_LONG_NP1 $line;
                            close OUTPUT_LONG_NP1;
                        }
                        print OUTPUT_DB_NP1 ">".$id_tmpi."\n";
                        print OUTPUT_DB_NP1 $line."\n";
						$removed_check = "";
                    }
                    else
                    {
                        $removed_reads_min_tmp++;
						$removed_check = "yes";
                    }

                    if ($quality_score_long eq "")
                    {
                        $ww = '1';
                    }
                    else
                    {
                        $ww++
                    }
                }
                elsif ($ww eq '3')
                {
                    $ww++;
                    next;
                }
                elsif ($ww eq '4')
                {
                    $ww = '1';
					if ($use_quality_scores_NP eq "yes" && $removed_check eq "")
                    {
                        my $pos_tmp = '0';
						my $tmp_line = "";
						while ($pos_tmp < length($line)-1)
                        {
                            my $qscore = substr $line, $pos_tmp, 1;
                            $pos_tmp++;
							if ($qscore eq "!" || $qscore eq "\"" || $qscore eq "#" || $qscore eq "\$" || $qscore eq "%" || $qscore eq "&" || $qscore eq "'" || $qscore eq "(")
							{	
								if ($tmp_line eq "")
								{
									$tmp_line = $pos_tmp;
								}
								else
								{
									$tmp_line .= ",".$pos_tmp;
								}
							}
                        }
						print $FILE_HASH_TMP $id_tmpi."\n";
						print $FILE_HASH_TMP $tmp_line."\n";
                    }
                }
            }
            
            close $FILE_LONG;
			if ($use_quality_scores_NP eq "yes")
            {
				close $FILE_HASH_TMP;
			}
            close OUTPUT_DB_NP1;
            
            my $DB_direc_tmp = $directory_DB_tmp.$NP_reads{$NP_reads_tmp};
            my $DB_output_tmp = $output_path."DB_".$NP_reads{$NP_reads_tmp}."_tmp_file.txt";    
            my $command_make_DB = "makeblastdb -in ".$output_file_DB_tmp." -dbtype nucl -out ".$DB_direc_tmp." > ".$DB_output_tmp."";
            system($command_make_DB);
                      
            my %hash_tmp;
            if ($file_count > 1)
            {
                $hash_tmp{'1'} = $DB_direc_tmp;
                $hash_tmp{'2'} = $read_count_tmp;
                $hash_tmp{'3'} = $DB_output_tmp;
                $hash_tmp{'4'} = $total_read_length_tmp;
                $hash_tmp{'5'} = $min_read_length_tmp;
                $hash_tmp{'6'} = $max_read_length_tmp;
                $hash_tmp{'7'} = $removed_reads_min_tmp;
                $hash_tmp{'8'} = [@read_lengths_tmp];
                $pm->finish(0, \%hash_tmp);
            }
            else
            {
                $ret_data_NP{'1'} = $DB_direc_tmp;
                $ret_data2_NP{'1'} = $DB_output_tmp;
                $read_count_NP += $read_count_tmp;
                $total_read_length_NP += $total_read_length_tmp;
                $min_read_length_NP = $min_read_length_tmp;
                $max_read_length_NP = $max_read_length_tmp;
                $removed_reads_min_NP = $removed_reads_min_tmp;
                @read_lengths_NP = (@read_lengths_NP, @read_lengths_tmp);
            }
        }
        if ($file_count > 1)
        {
            $pm->wait_all_children;
        }
        #my $long1 = keys %hash_long_reads;
        #my $long2 = keys %hash_long_reads_kmer;
        #print $long1." LONG1\n";
        #print $long2." LONG2\n";
    }
    print "...OK\n";
}

if ($input_reads_DB_folder_NP ne "")
{
    my $directory_DB_tmp1 = $directory_DB_NP.$output_path_test;
    opendir(DIR2, $directory_DB_tmp1) or die "Could not open $directory_DB_tmp1\n";
    my $count_tmp = '1';
    while (my $filename = readdir(DIR2))
    {
        if (-d $directory_DB_tmp1.$filename && $filename ne "." && $filename ne "..")
        {
            my $directory_DB_tmp2 = $directory_DB_tmp1.$filename.$output_path_test.$filename;
            $ret_data_NP{$count_tmp} = $directory_DB_tmp2;
            $count_tmp++;
        }
    }
    closedir DIR2;
}

if ($read_count_NP > 0)
{
    my $total_read_length_NP_tmp = '0';
    
    foreach my $read_length_tmp (sort {$b <=> $a} @read_lengths_NP)
    {
        $total_read_length_NP_tmp += $read_length_tmp;
        if ($total_read_length_NP_tmp > $total_read_length_NP/2)
        {
            $N50_NP = $read_length_tmp;
            last;
        }
    }
    
    print "\n\nNanopore reads\n";
    print "--------------\n";
    print "Total count             : ".$read_count_NP."\n";
    print "N50                     : ".$N50_NP." bp\n";
    print "Average read length     : ".int($total_read_length_NP/$read_count_NP)." bp\n";
    print "Max read length         : ".$max_read_length_NP." bp\n";
    print "Min read length         : ".$min_read_length_NP." bp\n";
    print "Reads below min length  : ".$removed_reads_min_NP."\n\n";
    
    print OUTPUT4 "\n\nNanopore reads\n";
    print OUTPUT4 "--------------\n";
    print OUTPUT4 "Total count             : ".$read_count_NP."\n";
    print OUTPUT4 "N50                     : ".$N50_NP." bp\n";
    print OUTPUT4 "Average read length     : ".int($total_read_length_NP/$read_count_NP)." bp\n";
    print OUTPUT4 "Max read length         : ".$max_read_length_NP." bp\n";
    print OUTPUT4 "Min read length         : ".$min_read_length_NP." bp\n";
    print OUTPUT4 "Reads below min length  : ".$removed_reads_min_NP."\n\n";
}
#--------------------------------------------------------------------------------------------------------------------------------------------------------
#Make hash of PacBio reads----------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------

if ($PB_reads eq "")
{
	goto SKIP_PACBIO;
}
undef %original_input_files;
	
my $read_count_PB = '0';
my $total_read_length_PB = '0';
my $max_read_length_PB = '0';
my $min_read_length_PB = '0';
my $removed_reads_min_PB = '0';
my $N50_PB = '0';
my @read_lengths_PB;

my $pm_PB = new Parallel::ForkManager($maxProcs_tmp);
my %ret_data_PB;
undef %ret_data_PB;
my %ret_data2_PB;
undef %ret_data2_PB;

$pm_PB->run_on_finish(
sub
{ 
    my @str = @_;
    my %hash_tmp = %{$str[5]};
    my $pid_PB = $str[0];

    # retrieve data structure from child
    if (defined($pid_PB))
    {  # children are not forced to send anything
        my $string = $hash_tmp{'1'};  # child passed a string reference
        my $string2 = $hash_tmp{'2'};
        my $string3 = $hash_tmp{'3'};
        my $string4 = $hash_tmp{'4'};
        my $string5 = $hash_tmp{'5'};
        my $string6 = $hash_tmp{'6'};
        my $string7 = $hash_tmp{'7'};
        my @string8 = @{$hash_tmp{'8'}};

        $ret_data_PB{$pid_PB} = $string;
        $ret_data2_PB{$pid_PB} = $string3;
        $read_count_PB += $string2;
        $total_read_length_PB += $string4;
        if ($string6 > $max_read_length_PB)
        {
            $max_read_length_PB = $string6;
        }
        if ($string5 < $min_read_length_PB || $min_read_length_PB eq '0')
        {
            $min_read_length_PB = $string5;
        }
        $removed_reads_min_PB += $string7;
        @read_lengths_PB = (@read_lengths_PB, @string8);
    }
    #else {  # problems occurring during storage or retrieval will throw a warning
  #print qq|No message received from child process $data_structure_reference!\n|;
#}
});

my $long_id_PB = '0';

if ($PB_reads ne "" || $input_reads_DB_folder_PB ne "")
{
    my $PB_check_directory = substr $PB_reads, -1;

	if ($input_reads_DB_folder_PB eq "" && $PB_reads ne "")
	{
		if ($PB_check_directory eq $output_path_test)
		{
			$PB_check_directory = "yes";
		}
		elsif (-d $PB_reads)
		{
			$PB_reads .= $output_path_test;
			$PB_check_directory = "yes";
		}
		elsif (-e $PB_reads)
		{
		}
		else
		{
			die "\nThe PacBio input doesn't not seems to be an existing file or folder!: ".$PB_reads."\n";
		}
	}
	
	select(STDERR);
    $| = 1;
    select(STDOUT); # default
    $| = 1;
     print "\nBuilding local databases for the PacBio reads...";
    
    my %PB_reads;
    undef %PB_reads;
    my $count_files = '1';
   
    if ($PB_reads ne "")
    {
        my $disered_db_count = int(($maxProcs_tmp)/2);
        if ($PB_check_directory eq "yes")
        {
            opendir(DIR_PB, $PB_reads) or die "Could not open $PB_reads\n";
            
            for my $filename (sort readdir(DIR_PB))
            {
                my $last5 = substr $filename, -5;
                if ($last5 eq "fasta" || $last5 eq "fastq" || $last5 eq "tq.gz" || $last5 eq "ta.gz" || $last5 eq "q.bz2" || $last5 eq "a.bz2")
                {
                    my $PB_reads_tmp = $PB_reads.$filename;
                    $PB_reads{$PB_reads_tmp} = $count_files."a";
                    $original_input_files{$PB_reads_tmp} = undef;
                    $count_files++;
                }
                elsif ($filename ne "." && $filename ne "..")
                {
                    print "\n".$filename.": File extension not recognized!\n";
                    print OUTPUT4 "\n".$filename.": File extension not recognized!\n";
                }
            }
            closedir DIR_PB;
        }
        else
        {
            $PB_reads{$PB_reads} = "1a";
            $original_input_files{$PB_reads} = undef;
            $count_files++;
        }
        
        my $total_lines = '0';
        my %file_lines;
        undef %file_lines;    
        
        foreach my $PB_reads_tmp (keys %PB_reads)
        {
            my $check_zip_long = substr $PB_reads_tmp, -2;
            my $check_zip_long2 = substr $PB_reads_tmp, -3;           
            
            if ($check_zip_long eq "gz")
            {
                my $new_filename = substr $PB_reads_tmp, 0, -3;
                my @new_filename = split /$output_path_test/, $new_filename;
                my $g = @new_filename;
                $g--;
                $new_filename = $output_path.$new_filename[$g];
                my $return_value = system("gzip -c -q -k -d ".$PB_reads_tmp." > ".$new_filename);
                
                $delete_input_files{$new_filename} = undef;
                
                $PB_reads{$new_filename} = $PB_reads{$PB_reads_tmp};
                delete $PB_reads{$PB_reads_tmp};
                my $count_lines = qx(wc -l $new_filename);
                chomp($count_lines);
                my @count_lines = split /\s+/, $count_lines;
                $total_lines += $count_lines[0];
                $file_lines{$new_filename} = $count_lines[0];
            }
            elsif ($check_zip_long2 eq "bz2")
            {
                die "Can't read bz2 files, pleas decompress: $PB_reads_tmp, $!\n";
            }
            elsif ($check_zip_long2 eq "zip")
            {
                die "Can't read zip files, pleas decompress: $PB_reads_tmp, $!\n";
            }
            else
            {

                my $count_lines = qx(wc -l $PB_reads_tmp);
                chomp($count_lines);
                my @count_lines = split /\s+/, $count_lines;
                $total_lines += $count_lines[0];
                $file_lines{$PB_reads_tmp} = $count_lines[0];
            }        
        }
        my $lines_perDB = int($total_lines/int(($maxProcs_tmp)/2))+8;
        my $adjusted_lines_perDB = (int($lines_perDB / 4) * 4)+4;
		my $PB_reads_DB = "";
PB_READS_DB: 
		
		my $count_total_files = keys %file_lines;
        my $count_files_tmp = '0';
        my $merge_files_check = "";
        my $merge_files_check_first = "";
        my %file_lines_merged;
        undef %file_lines_merged;
        my %PB_reads_tmp;
        undef %PB_reads_tmp;
     
		foreach my $PB_reads_tmp (sort keys %file_lines)
        {
            $count_files_tmp++;

            if ($file_lines{$PB_reads_tmp} > $adjusted_lines_perDB*1.15)
            {
                my $command_split = "split --additional-suffix=".$PB_reads{$PB_reads_tmp}.".fasta -l ".$adjusted_lines_perDB." ".$PB_reads_tmp." ".$directory_DB_PB.$output_path_test;                   
                my $return_value = system($command_split);
                
                my $expected_files = $file_lines{$PB_reads_tmp}/$adjusted_lines_perDB;
                my $expected_files0 = $file_lines{$PB_reads_tmp}/$adjusted_lines_perDB;
				
                if ($expected_files0 > int($expected_files0))
                {
                    $expected_files = int($expected_files0)+1;
                }
                my $found_files = 0;

                while ($found_files < $expected_files)
                {
                    $found_files = 0;
                    
                    # Check for the existence of the expected files
                    for my $file (glob("$directory_DB_PB$output_path_test*"))
                    {
                        if (-s $file)
						{
							$found_files++;
						}
                    }           
                    # Sleep for a bit if not all files have been found
                    sleep(10) if $found_files < $expected_files;
                }
                delete $PB_reads{$PB_reads_tmp};
				delete $file_lines{$PB_reads_tmp};
            }
            elsif ($merge_files_check ne "")
            {
                my $combined_lines = $file_lines{$PB_reads_tmp} + $file_lines_merged{$merge_files_check};

                if ($combined_lines < $adjusted_lines_perDB*1.15)
                {  
                    my $output_tmp = $directory_DB_PB.$output_path_test.$count_files_tmp."m".$merge_files_check_first.".fasta";
                    my $command_merge = "cat ".$merge_files_check." ".$PB_reads_tmp." > ".$output_tmp;                   
                    my $return_value = system($command_merge);
                    delete $PB_reads{$PB_reads_tmp};
                    delete $PB_reads{$merge_files_check};
					delete $file_lines{$PB_reads_tmp};

                    my $time_before_check_tmp = time;
COUNT_LINES_PB:                   
                    my $count_lines = qx(wc -l $output_tmp);
                    chomp($count_lines);
                    my @count_lines = split /\s+/, $count_lines;
                    my $count_lines2 = $count_lines[0];

                    if ($count_lines2 eq $combined_lines)
                    {      
                        unless (exists($original_input_files{$merge_files_check}))
                        {
                            unlink $merge_files_check;
                        }
                        unless (exists($original_input_files{$PB_reads_tmp}))
                        {
                            unlink $PB_reads_tmp;
                        } 
                    }
                    elsif (time > $time_before_check_tmp+600)
                    {
                        print "ERROR: Can't merge input files\n";
                        goto END1;
                    }
                    else
                    {
                        goto COUNT_LINES_PB;
                    }
                    
                    $PB_reads_tmp{$output_tmp} = $count_files_tmp."m".$merge_files_check_first;
                    $merge_files_check = $output_tmp;
                    $file_lines_merged{$merge_files_check} = $combined_lines;
                    my $output_file_DB_tmp = $directory_DB_PB.$output_path_test.$PB_reads_tmp{$output_tmp}.".fasta"; 
                    $delete_input_files{$output_file_DB_tmp} = undef;
                    $delete_input_files{$output_tmp} = undef;               
                }
                elsif ($file_lines_merged{$merge_files_check} > $adjusted_lines_perDB/1.4)
                { 
                    $PB_reads{$merge_files_check} = $PB_reads_tmp{$merge_files_check}."a";                  
                    my $output_file_DB_tmp = $directory_DB_PB.$output_path_test.$PB_reads{$merge_files_check}.".fasta"; 
                    $delete_input_files{$output_file_DB_tmp} = undef;
                    $merge_files_check = "";
                    $merge_files_check_first = "";
                }
            } 
            if ($file_lines{$PB_reads_tmp} < $adjusted_lines_perDB/1.3 && $count_files_tmp < $count_total_files && $count_total_files > int($maxProcs_tmp/2) && $merge_files_check eq "")
            {
                $merge_files_check = $PB_reads_tmp;
                $merge_files_check_first = substr $PB_reads{$PB_reads_tmp}, 0, -1;
                $PB_reads_tmp{$PB_reads_tmp} = substr $PB_reads{$PB_reads_tmp}, 0, -1;
                $file_lines_merged{$PB_reads_tmp} = $file_lines{$PB_reads_tmp};
                next;
            }
            my $output_file_DB_tmp = $directory_DB_PB.$output_path_test.$PB_reads{$PB_reads_tmp}.".fasta"; 
            $delete_input_files{$output_file_DB_tmp} = undef;         
        }      
		if ($merge_files_check ne "")
        {
            $PB_reads{$merge_files_check} = $PB_reads_tmp{$merge_files_check}."a";
            $delete_input_files{$merge_files_check} = undef;
            my $output_file_DB_tmp = $directory_DB_PB.$output_path_test.$PB_reads{$merge_files_check}.".fasta"; 
            $delete_input_files{$output_file_DB_tmp} = undef;
        }
		
		opendir(DIR_PB2, $directory_DB_PB.$output_path_test) or die "Could not open $directory_DB_PB.$output_path_test\n";
		my @files = sort readdir(DIR_PB2);
		foreach my $filename (@files)
		{
			my $last5 = substr $filename, -5;
			if ($last5 eq "fasta" || $last5 eq "fastq")
			{                       
				my $PB_reads_tmp2 = $directory_DB_PB.$output_path_test.$filename;
				$PB_reads{$PB_reads_tmp2} = $count_files."a";
				$count_files++;
				
				$delete_input_files{$PB_reads_tmp2} = undef;
				my $output_file_DB_tmp = $directory_DB_PB.$output_path_test.$PB_reads{$PB_reads_tmp2}.".fasta"; 
				$delete_input_files{$output_file_DB_tmp} = undef;
				
				my $count_lines = qx(wc -l $PB_reads_tmp2);
                chomp($count_lines);
                my @count_lines = split /\s+/, $count_lines;
                $total_lines += $count_lines[0];
                $file_lines{$PB_reads_tmp2} = $count_lines[0];
			}
			elsif ($filename ne "." && $filename ne "..")
			{
				print "\n".$filename.": File extension not recognized!\n";
				print OUTPUT4 "\n".$filename.": File extension not recognized!\n";
			}   
		}
		close DIR_PB2;
				
		my $file_count = keys %PB_reads;
		
		if ($file_count > int($maxProcs_tmp/2) && $PB_reads_DB eq "")
		{
			$PB_reads_DB = "yes";	
			goto PB_READS_DB;
		}
      
#Build local BLAST databse and save reads to disk---------------------------------------------------------------------------
       
        foreach my $PB_reads_tmp (keys %PB_reads)
        {
            my $pid_PB;
            if ($file_count > 1)
            {
                $pid_PB = $pm_PB->start and next;
            }
       srand();     
            my $check_zip_long = substr $PB_reads_tmp, -2;
            my $check_zip_long2 = substr $PB_reads_tmp, -3;
            my $firstLine_long;
            my $secondLine_long;
            my $thirdLine_long;        

            my $read_count_tmp = '0';
            my $total_read_length_tmp = '0';
            my $max_read_length_tmp = '0';
            my $min_read_length_tmp = '0';
            my $removed_reads_min_tmp = '0';
            my @read_lengths_tmp;
            
            my $FILE_tmp;
            if ($check_zip_long eq "gz")
            {
                open ($FILE_tmp, '-|', 'gzip', '-dc', $PB_reads_tmp) or die "Can't open file $PB_reads_tmp, $!\n";
            }
            elsif ($check_zip_long2 eq "bz2")
            {
                open ($FILE_tmp, '-|', 'bzip2', '-dc', $PB_reads_tmp) or die "Can't open file $PB_reads_tmp, $!\n";
            }
            else
            {
                open ($FILE_tmp, $PB_reads_tmp) or die "\n\nCan't open long reads file $PB_reads_tmp, $!\n";
            }
            $firstLine_long = <$FILE_tmp>;
            chomp $firstLine_long;
            $secondLine_long = <$FILE_tmp>;
            chomp $secondLine_long;
            $thirdLine_long = <$FILE_tmp>;
            chomp $thirdLine_long;
            close $FILE_tmp;

            my $no_quality_score_long_tmp = substr $thirdLine_long, 0, 1;
            my $quality_score_long = "";
            if ($thirdLine_long eq "+")
            {
                $quality_score_long = "yes";
            }
			
			my $read_limit = "10000000000000000000000000000000000000000000000000000";
            if ($subsample ne "")
            {
                my $command = "wc -l ".$PB_reads_tmp;
                my $lines = `$command`;
                chomp($lines);
                my @line_count = split / /, $lines;
                $read_limit = $subsample*($line_count[0]/4);
                if ($quality_score_long eq "")
                {
                    $read_limit = $subsample*($line_count[0]/2);
                }
            }
            
            my $FILE_LONG;
            if ($check_zip_long eq "gz")
            {
                open ($FILE_LONG, '-|', 'gzip', '-dc', $PB_reads_tmp) or die "Can't open file $PB_reads_tmp, $!\n";
            }
            elsif ($check_zip_long2 eq "bz2")
            {
                open ($FILE_LONG, '-|', 'bzip2', '-dc', $PB_reads_tmp) or die "Can't open file $PB_reads_tmp, $!\n";
            }
            else
            {
                open($FILE_LONG, $PB_reads_tmp) or die "\n\nCan't open long reads file $PB_reads_tmp, $!\n";
            }

            my $directory_tmp = $output_path."tmp_sequences_PB".$output_path_test.$PB_reads{$PB_reads_tmp};
            mkdir $directory_tmp;
            my $directory_DB_tmp = $directory_DB_PB.$output_path_test.$PB_reads{$PB_reads_tmp}.$output_path_test;
            mkdir $directory_DB_tmp;
            
            my $output_file_DB_tmp = $directory_DB_PB.$output_path_test.$PB_reads{$PB_reads_tmp}.".fasta";
            open(OUTPUT_DB_PB1, ">" .$output_file_DB_tmp) or die "\nCan't open file $output_file_DB_tmp, $!\n";
            
            my $ww = '1';
            my $fail_or_pass = "";
            my $fail_or_pass2 = "";
            my $id_read = "";
			my $id_tmpi = "";
			my $removed_check = "";
            my %hash_PB_reads_kmer_tmp;
            undef %hash_PB_reads_kmer_tmp;
			my $count_files_in_folder = '0';
			my $subfolder_name = '1';
			my $directory_DB_sub_tmp = $output_path."tmp_sequences_PB".$output_path_test.$PB_reads{$PB_reads_tmp}.$output_path_test.$subfolder_name;
            mkdir $directory_DB_sub_tmp;
  
FILE_LONG_PB:while (my $line = <$FILE_LONG>)
            {
                if ($long_id_PB > $read_limit)
                {
                    last FILE_LONG_PB;
                }
                chomp $line;    
                if ($fail_or_pass eq "fastq_fail" && $fail_or_pass2 > 0)
                {
                    $fail_or_pass2--;
                    next;
                }
                elsif ($ww eq '1')
                {
                    $fail_or_pass = substr $line, 0, 10;
                    if ($fail_or_pass eq "fastq_fail")
                    {
                        if ($quality_score_long eq "")
                        {
                            $fail_or_pass2 = '1';
                        }
                        else
                        {
                            $fail_or_pass2 = '3';
                        }
                        next;
                    }
                    if ($keep_read_ids eq "yes")
                    {
                        my @read_id = split / /, $line;
                        $id_read = $read_id[0];
                    }
                        
                    $ww++;
                    next;
                }
                elsif ($ww eq '2')
                {
                    #$line =~ tr/actgn/ACTGN/;
                    my $pos_tmp = '0';
                    $long_id_PB++;
                    
                    my $id_tmpi = $PB_reads{$PB_reads_tmp}.$long_id_PB;
                    if ($keep_read_ids eq "yes")
                    {
                        $id_tmpi = $PB_reads{$PB_reads_tmp}.$id_read;
                    }
                    my $length_tmp = length($line);   
                    $total_read_length_tmp += $length_tmp;
                    
                    if ($length_tmp > $max_read_length_tmp)
                    {
                        $max_read_length_tmp = $length_tmp;
                    }
                    if ($length_tmp < $min_read_length_tmp || $min_read_length_tmp eq '0')
                    {
                        $min_read_length_tmp = $length_tmp;
                    }
                    $read_count_tmp++;
                    
                    push @read_lengths_tmp, $length_tmp;
					
					if ($length_tmp >= $minimum_read_length_PB)
                    {
                        #substr $line, 0, 40, "";
                        $count_files_in_folder++;
						if ($count_files_in_folder > $max_file_count)
						{
							$count_files_in_folder = '1';
							$subfolder_name++;
							$directory_DB_sub_tmp = $output_path."tmp_sequences_PB".$output_path_test.$PB_reads{$PB_reads_tmp}.$output_path_test.$subfolder_name;
							mkdir $directory_DB_sub_tmp;
						}
						$id_tmpi .= "a".$subfolder_name;
                        my $output_file_PB1  = $directory_tmp.$output_path_test.$subfolder_name.$output_path_test."sequence_tmp_PB_".$id_tmpi.".fasta";
                        unless (-e $output_file_PB1)
                        {
                            open(OUTPUT_LONG_PB1, ">" .$output_file_PB1) or die "\nCan't open file $output_file_PB1, $!\n";
                            print OUTPUT_LONG_PB1 $line;
                            close OUTPUT_LONG_PB1;
                        }
                        print OUTPUT_DB_PB1 ">".$id_tmpi."\n";
                        print OUTPUT_DB_PB1 $line."\n";
						$removed_check = "";
                    }
                    else
                    {
                        $removed_reads_min_tmp++;
						$removed_check = "yes";
                    }

                    if ($quality_score_long eq "")
                    {
                        $ww = '1';
                    }
                    else
                    {
                        $ww++
                    }
                }
                elsif ($ww eq '3')
                {
                    $ww++;
                    next;
                }
                elsif ($ww eq '4')
                {
                    $ww = '1';
                    next;
                }
            }
            
            close $FILE_LONG;
            close OUTPUT_DB_PB1;
            
            my $DB_direc_tmp = $directory_DB_tmp.$PB_reads{$PB_reads_tmp};
            my $DB_output_tmp = $output_path."DB_".$PB_reads{$PB_reads_tmp}."_tmp_file.txt";    
            my $command_make_DB = "makeblastdb -in ".$output_file_DB_tmp." -dbtype nucl -out ".$DB_direc_tmp." > ".$DB_output_tmp."";
            system($command_make_DB);
                      
            my %hash_tmp;
            if ($file_count > 1)
            {
                $hash_tmp{'1'} = $DB_direc_tmp;
                $hash_tmp{'2'} = $read_count_tmp;
                $hash_tmp{'3'} = $DB_output_tmp;
                $hash_tmp{'4'} = $total_read_length_tmp;
                $hash_tmp{'5'} = $min_read_length_tmp;
                $hash_tmp{'6'} = $max_read_length_tmp;
                $hash_tmp{'7'} = $removed_reads_min_tmp;
                $hash_tmp{'8'} = [@read_lengths_tmp];
                $pm_PB->finish(0, \%hash_tmp);
            }
            else
            {
                $ret_data_PB{$pid_PB} = $DB_direc_tmp;
                $ret_data2_PB{$pid_PB} = $DB_output_tmp;
                $read_count_PB += $read_count_tmp;
                $total_read_length_PB += $total_read_length_tmp;
                $min_read_length_PB = $min_read_length_tmp;
                $max_read_length_PB = $max_read_length_tmp;
                $removed_reads_min_PB = $removed_reads_min_tmp;
                @read_lengths_PB = (@read_lengths_PB, @read_lengths_tmp);
            }
        }
        if ($file_count > 1)
        {
            $pm_PB->wait_all_children;
        }
    }
    print "...OK\n";
}
SKIP_PACBIO:
if ($input_reads_DB_folder_PB ne "")
{
    my $directory_DB_tmp1 = $directory_DB_PB.$output_path_test;
    opendir(DIR_PB3, $directory_DB_tmp1) or die "Could not open $directory_DB_tmp1\n";
    my $count_tmp = '1';
    while (my $filename = readdir(DIR_PB3))
    {
        if (-d $directory_DB_tmp1.$filename && $filename ne "." && $filename ne "..")
        {
            my $directory_DB_tmp2 = $directory_DB_tmp1.$filename.$output_path_test.$filename;
            $ret_data_PB{$count_tmp} = $directory_DB_tmp2;
            $count_tmp++;
        }
    }
    closedir DIR_PB3;
}

if ($read_count_PB > 0)
{
    my $total_read_length_PB_tmp = '0';
    
    foreach my $read_length_tmp (sort {$b <=> $a} @read_lengths_PB)
    {
        $total_read_length_PB_tmp += $read_length_tmp;
        if ($total_read_length_PB_tmp > $total_read_length_PB/2)
        {
            $N50_PB = $read_length_tmp;
            last;
        }
    }
    
    print "\n\nPacBio reads\n";
    print "--------------\n";
    print "Total count             : ".$read_count_PB."\n";
    print "N50                     : ".$N50_PB." bp\n";
    print "Average read length     : ".int($total_read_length_PB/$read_count_PB)." bp\n";
    print "Max read length         : ".$max_read_length_PB." bp\n";
    print "Min read length         : ".$min_read_length_PB." bp\n";
    print "Reads below min length  : ".$removed_reads_min_PB."\n\n";
    
    print OUTPUT4 "\n\nPacBio reads\n";
    print OUTPUT4 "--------------\n";
    print OUTPUT4 "Total count             : ".$read_count_PB."\n";
    print OUTPUT4 "N50                     : ".$N50_PB." bp\n";
    print OUTPUT4 "Average read length     : ".int($total_read_length_PB/$read_count_PB)." bp\n";
    print OUTPUT4 "Max read length         : ".$max_read_length_PB." bp\n";
    print OUTPUT4 "Min read length         : ".$min_read_length_PB." bp\n";
    print OUTPUT4 "Reads below min length  : ".$removed_reads_min_PB."\n\n";
}  

#spin up worker early before creating big hash---------

my $chnl;
$chnl = MCE::Channel->new( impl => 'Simple' );

mce_child
{
    local $SIG{__WARN__} = sub {};
    while ( my ($cmd, @args) = $chnl->recv ) {
        local ($?, $!);
        system($cmd, @args);
        $chnl->send2($?, $!);
    }
};

sub syscmd {
    my $cmd = shift;
    return unless $cmd;

    $chnl->send($cmd, @_);
    my ($status, $errmsg) = $chnl->recv2;
    
    if ($status == -1) {
        print "SYSTEM: failed to execute ($cmd): $errmsg\n";
    }
    elsif ($status & 127) {
        printf "SYSTEM: $cmd died with signal %s, %s coredump\n",
            ($status & 127), ($status & 128) ? 'with' : 'without';
    }
    else {
        #printf "SYSTEM: $cmd exited with status %d\n", $status >> 8;
    }
}

#Load long reads quality hashes from file--------------------------------------------------
if ($use_quality_scores_NP eq "yes")
{
	foreach my $hash_tmp (@QUALITY_HASH)
	{
		my $merge_command2 = "cat ".$hash_tmp." >> ".$output_path."quality_scores_".$project.".txt";
		system($merge_command2);
		my $FILE_Q;
		open($FILE_Q, $hash_tmp) or die "Can't open variation file $hash_tmp, $!\n";
		my $id_tmp = "";
		while (my $line = <$FILE_Q>)
		{
			chomp($line);
			if ($id_tmp eq "")
			{
				$id_tmp = $line;
				next;
			}
			else
			{
				$quality_scores_NP{$id_tmp} = $line;
				$id_tmp = "";
			}   
		}
		close $FILE_Q;
	}
}
elsif ($use_quality_scores_NP ne "")
{
    my $count_kmer = '0';
    my $longest_kmer = '0';
    my $read_count = '0';

    open(INPUT7, $use_quality_scores_NP) or die "Can't open quality score file, $!\n";
    my $id_tmp = "";

    while (my $line = <INPUT7>)
    {
        chomp($line);
        if ($id_tmp eq "")
        {
            $id_tmp = $line;
            next;
        }
        else
        {
            $quality_scores_NP{$id_tmp} = $line;
			$id_tmp = "";
        }   
    }
    close INPUT7;
}

select(STDERR);
$| = 1;
select(STDOUT); # default
$| = 1;

print "\nPrepare seeds...";
print OUTPUT4 "\nPrepare seeds...\n";

#Retrieve first read from the given seed-----------------------------------------------------------------------------------------------------------------------------------

my $si = '0';
my $space_at_end2 = "";
my $id_line = "";
my %seeds_list;
my %seeds_list_sorted;
my $count_seeds0 = '0';

if ($assembly_length_max ne "WG")
{
	while (my $line = <INPUT3>)
	{
		chomp($line);
		$line =~ tr/\r//d;
		$line =~ s/\R/\012/;
		$line =~ s/[ \t\xA0]+$//; 
		if ($si eq 0)
		{
			my $last_character = substr $line, -1;       
			if ($last_character =~ m/\s|\t/g)
			{
				$space_at_end2 = "yes";
			}
			$id_line = $line;
		}
		if ($si > 0)
		{
			my $first_letter = substr $line, 0, 1;
			if ($first_letter eq ">")
			{
				$count_seeds0++;
				my $id_tmp = substr $id_line, 1;
				$seeds_list{$id_tmp} = $seed_input;
				$seeds_list_sorted{$count_seeds0} = $id_tmp;
				$y{$id_tmp} = '1';
				$seed_input = "";
				$seed_batch = "yes";
				$id_line = $line;
			}
			else
			{
				if ($space_at_end2 eq "yes")
				{
					chop($line);
				}
				my $seed_input_tmp = $seed_input;
				$seed_input = $seed_input_tmp.$line;
			}   
		}
		$si++;
	}
}
$count_seeds0++;
my $id_tmp = substr $id_line, 1;
$seeds_list{$id_tmp} = $seed_input;
$seeds_list_sorted{$count_seeds0} = $id_tmp;
$y{$id_tmp} = '1';

if ($ploidy > 1 && $count_seeds0 < 2)
{
    $find_haps_in_seed = "yes";
}

my %find_haps_in_seed;
my $found_haps_in_seed = "";
my $next_seed_print = "yes";
my $last_seed = "";
my $first_back_assembly = "";

#SEED select for WG mode----------------------------------------------------------------------------------------------------------------

my %split_contigs_reads;
undef %split_contigs_reads;
my %split_contigs_reads2;
undef %split_contigs_reads2;
my %split_contigs_ends;
undef %split_contigs_ends;
my %contig_connections;
undef %contig_connections;
my %reads_as_seeds;
undef %reads_as_seeds;
my $keep_track_of_reads_number = '0';
my $output_file20 = "";
if ($assembly_length_max eq "WG")
{
	my $output_file19  = $output_path."Assembled_read_ids_".$project.".txt";
	open(OUTPUT19, ">".$output_file19) or die "Can't open file $output_file19, $!\n";
	
	$output_file20  = $output_path."assemblies_all_".$project.".fasta";
	unless (-e $output_file20)
	{	
		open(OUTPUT20, ">".$output_file20) or die "Can't open file $output_file20, $!\n";
	}
	else
	{
		open(OUTPUT20, ">>".$output_file20) or die "Can't open file $output_file20, $!\n";
	}
	
	if ($seed_input0 ne "")
	{
		while (my $line = <INPUT3>)
		{
			chomp($line);
			$reads_as_seeds{$line} = undef;
			print OUTPUT19 $line."\n";
		}
	}
}

NEXT_SEED:

if ($assembly_length_max eq "WG")
{
	my $output_folder_tmp = $input_reads_DB_folder_NP;
	if ($input_reads_DB_folder_NP eq "")
	{
		$output_folder_tmp = $output_path;
	}
	my $DB_path = $output_folder_tmp."tmp_sequences_NP";
	
	if ($PB_reads ne "" || $input_reads_DB_folder_PB ne "")
	{
		$output_folder_tmp = $input_reads_DB_folder_PB;
		if ($input_reads_DB_folder_PB eq "")
		{
			$output_folder_tmp = $output_path;
		}
		$DB_path = $output_folder_tmp."tmp_sequences_PB";
	}
	
	opendir(DIR_DB_NP, $DB_path) or die "Could not open $DB_path\n";
	
	my $keep_track_of_reads_number_current = '0';
	
DIR_DB:	for my $filename (sort readdir(DIR_DB_NP))
	{
		if ($filename ne "." && $filename ne "..")
		{
			if (-d $DB_path.$output_path_test.$filename)
			{
				my $DB_path2 = $output_folder_tmp."tmp_sequences_NP".$output_path_test.$filename;
				
				if ($PB_reads ne "" || $input_reads_DB_folder_PB ne "")
				{
					$DB_path2 = $output_folder_tmp."tmp_sequences_PB".$output_path_test.$filename;
				}
	
				opendir(DIR_DB_NP2, $DB_path2) or die "Could not open $DB_path2\n";

				for my $filename2 (sort readdir(DIR_DB_NP2))
				{
					if ($filename2 ne "." && $filename2 ne "..")
					{
						if (-d $DB_path2.$output_path_test.$filename2)
						{
							my $DB_path3 = $output_folder_tmp."tmp_sequences_NP".$output_path_test.$filename.$output_path_test.$filename2;
				
							if ($PB_reads ne "" || $input_reads_DB_folder_PB ne "")
							{
								$DB_path3 = $output_folder_tmp."tmp_sequences_PB".$output_path_test.$filename.$output_path_test.$filename2;
							}
							
							opendir(DIR_DB_NP3, $DB_path3) or die "Could not open $DB_path3\n";

							for my $filename3 (sort readdir(DIR_DB_NP3))
							{
								if ($filename3 ne "." && $filename3 ne "..")
								{
									my $check_rev = substr $filename3, -9, 3;
									if ($check_rev ne "rev")
									{
										$keep_track_of_reads_number_current++;
										if ($keep_track_of_reads_number_current > $keep_track_of_reads_number)
										{
											substr $filename3, -6, 6, "";
											$keep_track_of_reads_number++;
											my @file_name = split /_/, $filename3;
											if (exists($reads_as_seeds{$file_name[3]}))
											{}
											else
											{
												my $sequence_tmp_file  = $DB_path3.$output_path_test.$filename3.".fasta";
												open(OUTPUT_SEQ_TMP, $sequence_tmp_file) or die "\nCan't open file $sequence_tmp_file, $!\n";
												$seed_input = <OUTPUT_SEQ_TMP>;
												$reads_as_seeds{$file_name[3]} = undef;
												print OUTPUT19 $file_name[3]."\n";

												if (length($seed_input) > 1500)
												{
													my $output_file25  = $TMP_directory."Sequence_blast_test_".$file_name[3].".fasta";
													open(OUTPUT25, ">".$output_file25) or die "Can't open file $output_file25, $!\n";
													OUTPUT25->autoflush(1);
													my $last1000 = substr $seed_input, -1000;
													print OUTPUT25 $last1000."\n";	
													close OUTPUT25;
													
													if (-e $output_file20 && -s $output_file25)
													{			
														my $file_tmp = $TMP_directory."blast_seed_test_".$file_name[3].".txt";
														my $command = "blastn -query ".$output_file25." -subject ".$output_file20." -out ".$file_tmp." -outfmt 7 -qcov_hsp_perc 40";
														syscmd($command);
														my $count_tmp = '0';
WG_SEED_ACCEPT:									
														if (-s $file_tmp)
														{
															open(SEED_TEST, $file_tmp) or die "\nCan't open file $file_tmp, $!\n";
															
															my $count_lines_tmp = '1';
															while (my $line_tmp = <SEED_TEST>)
															{
																chomp($line_tmp);
																if ($count_lines_tmp eq '4' && $line_tmp eq "# 0 hits found")
																{
																	close SEED_TEST;
																	goto WG_SEED_ACCEPT2;
																}
																elsif ($count_lines_tmp eq '5' && $line_tmp eq "# BLAST processed 1 queries")
																{
																	close SEED_TEST;
																	goto WG_SEED_ACCEPT;
																}
																elsif ($count_lines_tmp > 5)
																{
																	my @line_tmp = split /\t/, $line_tmp;
																	my $accuracy_tmp = $line_tmp[2];
							
																	if ($accuracy_tmp > 95 || ($PB_reads eq "" && $input_reads_DB_folder_PB eq "" && $accuracy_tmp > 70))
																	{
																		close SEED_TEST;
																		goto WG_SEED_ACCEPT3;
																	}  
																}
																$count_lines_tmp++;
															}
															close SEED_TEST;
														}
														elsif ($count_tmp < 100000)
														{
															$count_tmp++;
															goto WG_SEED_ACCEPT;
														}
													}
WG_SEED_ACCEPT2:									
													substr $seed_input, 0, 200, "";
													substr $seed_input, -200, 200, "";
				
													$seeds_list{$file_name[3]} = $seed_input;
													$seeds_list_sorted{$keep_track_of_reads_number} = $file_name[3];
													$y{$file_name[3]} = '1';
													
													close OUTPUT_SEQ_TMP;
													close DIR_DB_NP2;
													close DIR_DB_NP3;
													closedir DIR_DB_NP;
													print OUTPUT4 $file_name[3]." NEW_SEED\n";
													last DIR_DB;
												}
WG_SEED_ACCEPT3:								
												close OUTPUT_SEQ_TMP;
												close OUTPUT25;
											}
										}
									}
								}
							}
						}
					}
				}
				close DIR_DB_NP2;
			}
		}
	}
	closedir DIR_DB_NP;

	if ($seed_input eq "")
	{
		$last_seed = "yes";
		goto END1;
	}
}

print "...OK\n\n";
print OUTPUT4 "...OK\n\n";

FIRST_SEED:

my $low_coverage_check = '0';

foreach my $seed_input_tmp2 (sort {$a <=> $b} keys %seeds_list_sorted)
{
    my $seed_input_tmp = $seeds_list{$seeds_list_sorted{$seed_input_tmp2}};
    my $id_of_seed = $seeds_list_sorted{$seed_input_tmp2};
	if ($first_back_assembly{$seed_input_tmp2} ne "yes")
	{
		$seed_input_tmp =~ tr/actgn/ACTGN/;
		$original_seed_length{$id_of_seed} = length($seed_input_tmp);
	}
    if ($next_seed_print eq "yes")
    {
        push @seed_list_sorted, $seeds_list_sorted{$seed_input_tmp2};
    }
  
    if ($next_seed_print eq "yes" && $seed_batch eq "yes")
    {
        print "\nSeed ".$seed_input_tmp2.": ".$seeds_list_sorted{$seed_input_tmp2}."\n";
        print OUTPUT4 "\nSeed ".$seed_input_tmp2.": ".$seeds_list_sorted{$seed_input_tmp2}."\n";
        $next_seed_print = "";
        #push @seed_list_sorted, $seeds_list_sorted{$seed_input_tmp2};
    }
	my $fh;
	$filehandle{$id_of_seed} = $fh;
	$output_file5  = $output_path."log_extended_".$project."_".$id_of_seed.".txt";
	open($filehandle{$id_of_seed}, ">".$output_file5) or die "Can't open file $output_file5, $!\n";

	my $fh2;
	$filehandle4{$id_of_seed} = $fh2;
	$output_file13 = $output_path."quality_scores_".$project."_".$id_of_seed.".txt";
	open($filehandle4{$id_of_seed}, ">".$output_file13) or die "Can't open file $output_file13, $!\n";

	$seed{$id_of_seed} = $seed_input_tmp;
	$position{$id_of_seed} = length($seed{$id_of_seed});
	$position = length($seed{$id_of_seed});
	if ($find_haps_in_seed eq "yes")
	{
		$find_haps_in_seed{$id_of_seed} = "yes";
	}
    
    if ($find_haps_in_seed eq "")
    {        
        print "\nExtend the given seed directly: ".$id_of_seed."\n";
        print OUTPUT4 "\nExtend the given seed directly: ".$id_of_seed."\n";
        print {$filehandle{$id_of_seed}} "\nExtend the given seed directly: ".$id_of_seed."\n";

        #delete $seeds_list_sorted{$seed_input_tmp2};
        #delete $seeds_list{$seeds_list_sorted{$seed_input_tmp2}};
        $next_seed_print = "yes";
    }	
}

#my $output_file18  = $output_path."SNP_POS_".$project.".txt";
#open(OUTPUT18, ">".$output_file18) or die "Can't open file $output_file18, $!\n";

#----------------------------------------------------------------------------------------------------------------------------
#Prepare the haplotypes------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------

my %exclude_reads_hap1_NP;
my %exclude_reads_hap2_NP;
my %exclude_reads_hap1_PB;
my %exclude_reads_hap2_PB;
my %exclude_reads_hap1_NP_back;
my %exclude_reads_hap2_NP_back;
my %exclude_reads_hap1_PB_back;
my %exclude_reads_hap2_PB_back;
undef %exclude_reads_hap1_NP;
undef %exclude_reads_hap2_NP;
undef %exclude_reads_hap1_PB;
undef %exclude_reads_hap2_PB;
undef %exclude_reads_hap1_NP_back;
undef %exclude_reads_hap2_NP_back;
undef %exclude_reads_hap1_PB_back;
undef %exclude_reads_hap2_PB_back;

if ($find_haps_in_seed eq "" && $ploidy > 1)
{
    my $w = '0';
    my $input_file_seed_tmp = $TMP_directory."sequence_seed_tmp.fasta";
    my $output_file_seed_tmp = $TMP_directory."mafft_seed.txt";
    
    open(OUTPUT_LONG1, ">" .$input_file_seed_tmp) or die "\nCan't open file $input_file_seed_tmp, $!\n";
    
    foreach my $seed_input_tmp2 (sort {$a <=> $b} keys %seeds_list_sorted)
    {                                               
        print OUTPUT_LONG1 ">".$seed_input_tmp2."\n";
        print OUTPUT_LONG1 $seeds_list{$seeds_list_sorted{$seed_input_tmp2}}."\n";
        $hap_compare_pos{$seeds_list_sorted{$seed_input_tmp2}} = length($seeds_list{$seeds_list_sorted{$seed_input_tmp2}});
    }
    close OUTPUT_LONG1;
    chomp($input_file_seed_tmp);
        
    #my $command = "blastn -query ".$ref_seed_file." -subject ".$output_file_seed_tmp." -out blast_seed_tmp.txt -qcov_hsp_perc 60 -outfmt 4 -gapextend 1 -gapopen 2 &";
    my $cmd_seed = sprintf("mafft --thread 4 --op 0.4 --ep 1.06 --quiet --clustalout %s > ".$output_file_seed_tmp, $input_file_seed_tmp);
    syscmd($cmd_seed);
    
    
    open(INPUT_BLAST_VAR, $output_file_seed_tmp) or print "\n\nCan't open seed mafft file $output_file_seed_tmp, $!\n";                                 
    
#merge mafft lines-------------------------------------------------------------------------------             
    my $g = '0';
    my $query_line = "";
    my %subject_list_seed;
    undef %subject_list_seed;
    my $consensus_total_seed;
    
                    
INPUT_BLAST_VAR:while (my $line2 = <INPUT_BLAST_VAR>)
    {                                                     
        chomp($line2);
        if ($g > 2)
        {
            my @blast_result_tmp = split /\s+/, $line2;
            
            if (exists($seed{$seeds_list_sorted{$blast_result_tmp[0]}}))
            {
                my $subject_tmp = $subject_list_seed{$blast_result_tmp[0]};
                $subject_list_seed{$blast_result_tmp[0]} = $subject_tmp.$blast_result_tmp[1];
                $query_line = "yes";
            }  
            elsif ($query_line eq "yes")
            {
                my $consensus = substr $line2, 16, 60;
                $consensus_total_seed .= $consensus;
            }
        }
        $g++
    }
    close INPUT_BLAST_VAR;
#------------------------------------------------------------------------------------------------------
    
    chomp($consensus_total_seed);
    $consensus_total_seed =~ s/\s+$//;
    my @consensus_seed = split //, $consensus_total_seed;
    my $c = '0';                  
                    
PREPARE_HAP: foreach my $cons (@consensus_seed)
    {
        if ($cons eq "." || $cons eq " ")
        {
            foreach my $subject_tmp (keys %subject_list_seed)
            {
                my $seq_tmp = $subject_list_seed{$subject_tmp};
                my $nuc_match = substr $seq_tmp, $c, 1;
                $nuc_match =~ tr/actgn/ACTGN/;
                if ($nuc_match eq "N")
                {
                    $c++;
                    next PREPARE_HAP;
                }
            }
            foreach my $subject_tmp (keys %subject_list_seed)
            {
                my $seq_tmp = $subject_list_seed{$subject_tmp};
                my $nuc_match = substr $seq_tmp, $c, 1;
                $nuc_match =~ tr/actgn/ACTGN/;
                my $d = $c-$overlap;
                my $overlap_tmp = $overlap;
                if ($d < 0)
                {
                   $overlap_tmp += $d;
                   $d = '0';                         
                }
                my $read_end_tmpi = substr $seq_tmp, $d, $overlap_tmp;
                $read_end_tmpi =~ tr/-//d;
                $read_end_tmpi =~ tr/actgn/ACTGN/;
                my $gaps_correct_seq = substr $subject_list_seed{$subject_tmp}, 0, $c;
                my $gaps_correct = $gaps_correct_seq =~ tr/-/-/;
                my $pos = $c-$gaps_correct;
                if (exists($split_positions{$subject_tmp}{$pos}))
                {}
                else
                {
                    $split_positions{$seeds_list_sorted{$subject_tmp}}{$pos} = $read_end_tmpi.",".$nuc_match;
                    print OUTPUT18 $seeds_list_sorted{$subject_tmp}."\t".$pos."\t".$read_end_tmpi.",".$nuc_match."\n";
                }      
            }
        }
        $c++;
    }      
}

foreach my $seed_id_tmp (keys %seed)
{
    my $fh2;
    $filehandle3{$seed_id_tmp} = $fh2;
    my $output_file6  = $output_path."contigs_tmp_".$seed_id_tmp."_".$project.".fasta";
    open($filehandle3{$seed_id_tmp}, ">".$output_file6) or die "Can't open file $output_file6, $!\n";
}

#----------------------------------------------------------------------------------------------------------------------------------------------------
#Start assembly--------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------

close INPUT3;

my $sc = '0';
my $best_extension = "";
my $last_iteration_check = '0';

foreach (keys %seed)
{
    $sc++;    
}
if ($sc eq '1')
{
    print "\n\nStart Assembly...\n\n";
    print OUTPUT4 "\n\nStart Assembly...\n\n";
}
my $e = "0";
while ($y0 < $iterations)
{
    $e++;
SEED:
    my $shortest_seed = "no";
    my $seed_id = "";
    foreach my $seed_id_tmp (keys %seed)
    {    
	    if (((length($seed{$seed_id_tmp}) - $hap_compare_pos{$seed_id_tmp} < $shortest_seed) || $shortest_seed eq "no") && $skip_hap ne $seed_id_tmp)
        {
            $seed_id = $seed_id_tmp;
            $shortest_seed = length($seed{$seed_id_tmp}) - $hap_compare_pos{$seed_id_tmp};
        }
    }
    
    if ($e eq "20")
    {
        foreach my $seed_id (keys %seed)
        {  
            my $output_file6  = $output_path."contigs_tmp_".$seed_id."_".$project.".fasta";
            open($filehandle3{$seed_id}, ">".$output_file6) or die "Can't open file $output_file6, $!\n";
            print {$filehandle3{$seed_id}} ">".$seed_id."\n";
            my $m = '0';
            while (length($seed{$seed_id}) > $m)
            {
                my $value_ref2b = substr $seed{$seed_id}, $m, 150;
                $m += 150;
                print {$filehandle3{$seed_id}} $value_ref2b."\n";
            }
            
        }
        $e = '0';
    }
FULL_RESET: 
    if ($seed_id eq "")
    {
        goto END1;
    }
	if (exists($filehandle{$seed_id}))
    {
        $seed_id2 = $seed_id;
    }
    else
    {
        print $seed_id." SEED_NO_MATCH\n"; 
        $seed_id2 = $seed_id;      
        delete $seed{$seed_id};
        goto SEED;
    }
   
    if (exists($filehandle{$seed_id2}))
    {
    }
    else
    {
        print $seed_id2." SEED_ERROR2\n";   
        delete $seed{$seed_id};
        goto SEED;
    }
    
    $best_extension = "";
    if (($NP_reads ne "" || $input_reads_DB_folder_NP ne "") && ($PB_reads ne "" || $input_reads_DB_folder_PB ne "") && $NP_reads_support ne "yes2")
    {
        $NP_reads_support = "yes";
    }
	$y = $y{$seed_id};
    
    print {$filehandle{$seed_id2}} "\n".$y."\n\n";   
        
    if (exists($seed{$seed_id}))
    {     
        $read = $seed{$seed_id};
           
        print {$filehandle{$seed_id2}} "\n".$seed_id." SEED_exists\n\n";
        print {$filehandle{$seed_id2}} length($read)." READ_LENGTH\n";
        $id = $seed_id;
        
        if (exists($find_haps_in_seed{$id}))
        {
            $find_haps_in_seed = $find_haps_in_seed{$id};
        }
        else
        {
            $find_haps_in_seed = "";
        }
		$position = $position{$id};
		if ($position ne length($read))
		{
			print {$filehandle{$seed_id2}} $position." POS_ERROR\n";
		}
		
        if ($y eq '1' && $find_haps_in_seed eq "yes")
        {
            $position{$id} = '0';
            $position = '0';
            $position_back = '0';
            $position_back{$id} = '0';
        }
        
        my $hap_tag = substr $seed_id, -4;
		$first_back_assembly = $first_back_assembly{$id};
        
        my $time_START0 = time;
        
#Check positions of excluded and confirmed reads--------------------------------------------

        if (($NP_reads ne "" || $input_reads_DB_folder_NP ne "") && $find_haps_in_seed eq "")
        {
            if ($hap_tag eq "HAP1")
            {
                foreach my $confrim_id (keys %exclude_reads_hap1_NP)
                {
                    if ($position >= $exclude_reads_hap1_NP{$confrim_id})
                    {
                        delete $exclude_reads_hap1_NP{$confrim_id};
                    }
                }
                foreach my $confrim_id (keys %exclude_reads_hap1_NP_back)
                {
                    if ($position_back >= $exclude_reads_hap1_NP_back{$confrim_id})
                    {
                        delete $exclude_reads_hap1_NP_back{$confrim_id};
                    }
                }
            }
            elsif ($hap_tag eq "HAP2")
            {
                foreach my $confrim_id (keys %exclude_reads_hap2_NP)
                {
                    if ($position >= $exclude_reads_hap2_NP{$confrim_id})
                    {
                        delete $exclude_reads_hap2_NP{$confrim_id};
                    }
                }
                foreach my $confrim_id (keys %exclude_reads_hap2_NP_back)
                {
                    if ($position_back >= $exclude_reads_hap2_NP_back{$confrim_id})
                    {
                        delete $exclude_reads_hap2_NP_back{$confrim_id};
                    }
                }
			}
        }
		
        if (($PB_reads ne " "|| $input_reads_DB_folder_PB ne "") && $find_haps_in_seed eq "")
        {
            if ($hap_tag eq "HAP1")
            {
                foreach my $confrim_id (keys %exclude_reads_hap1_PB)
                {
                    if ($position >= $exclude_reads_hap1_PB{$confrim_id})
                    {
                        delete $exclude_reads_hap1_PB{$confrim_id};
                    }
                }
                foreach my $confrim_id (keys %exclude_reads_hap1_PB_back)
                {
                    if ($position_back >= $exclude_reads_hap1_PB_back{$confrim_id})
                    {
                        delete $exclude_reads_hap1_PB_back{$confrim_id};
                    }
                }
            }
            elsif ($hap_tag eq "HAP2")
            {
                foreach my $confrim_id (keys %exclude_reads_hap2_PB)
                {
                    if ($position >= $exclude_reads_hap2_PB{$confrim_id})
                    {
                        delete $exclude_reads_hap2_PB{$confrim_id};
                    }
                }
                foreach my $confrim_id (keys %exclude_reads_hap2_PB_back)
                {
                    if ($position_back >= $exclude_reads_hap2_PB_back{$confrim_id})
                    {
                        delete $exclude_reads_hap2_PB_back{$confrim_id};
                    }
                }
            }
        }
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------       
#Check var between haps----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------         
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------       
        my $count_seeds = keys %seed;
        my $check_var_between_hap = '0';
        my $extra_seq = "no";
        my $var_length = '3000';
        my $mafft_extra_quality_haps = "";
        
        if ($count_seeds > 1 && $ploidy > 1)
        {
            if ($PB_reads ne "")
            {
                $var_length = '1000';
            }
            if ($compare_haps ne "")
            {
                my %var_length;
                my $var_min = "";
                foreach my $id_tmp (keys %seed)
                {
                    my $var_length_tmp = length($seed{$id_tmp})-$hap_compare_pos{$id_tmp}-$hap_compare_mismatch_extend;
                    if ($var_length_tmp < $var_min || $var_min eq "")
                    {
                        $var_min = $var_length_tmp;
                    }
                }
                if ($var_min > 250)
                {
                    $var_length = $var_min;
                    print {$filehandle{$seed_id2}} $var_length." FORCED_HAP_COMPARE\n";
                }
                else
                {
                    $compare_haps = "";
                }
            }
            
            foreach my $id_tmp (keys %seed)
            {
                if (length($seed{$id_tmp}) > $hap_compare_pos{$id_tmp} + $var_length + $hap_compare_mismatch_extend)
                {
                    $check_var_between_hap++;
                    if (length($seed{$id_tmp})-($hap_compare_pos{$id_tmp} + $var_length + $hap_compare_mismatch_extend) < $extra_seq || $extra_seq eq "no")
                    {
                        $extra_seq = length($seed{$id_tmp})-($hap_compare_pos{$id_tmp} + $var_length + $hap_compare_mismatch_extend);
                    }
                }
            }
            if ($extra_seq < 0 || $extra_seq eq "no")
            {
                $extra_seq = 0;
            }
            
            if ($hap_compare_mismatch_extend > 5000)
            {
                $compare_haps = "";
                $check_var_between_hap = '0';
                $compare_haps_stop = "yes";
            }
#Find the lost haplotypes again-----------------------------------------------------------------------------------------------------------------

            if ($hap_compare_mismatch_extend > 5000 && (($position > $prev_position_hap_compare{$seed_id}+15000) || ($position > $prev_position_hap_compare{$seed_id}+5000 && ($PB_reads ne "" || $input_reads_DB_folder_PB ne ""))))
            {      
                print {$filehandle{$seed_id2}} $hap_compare_mismatch_extend." REALIGN HAPS\n";
				foreach my $id_tmp5 (keys %seed)
                {
                    my $last_1000 = substr $seed{$id_tmp5}, -1000;
					my $check_last_10000 = "";
					
                    my $ref_part = "";
                    foreach my $id_tmp2 (keys %seed)
                    {
                        if ($id_tmp2 ne $id_tmp5)
                        {
							$ref_part = substr $seed{$id_tmp2}, $hap_compare_pos{$id_tmp2};
                        }
                    }
CHECK_LAST_10000:                    
                    my $ref_file = $TMP_directory."hap_var_ref_tmp_".$id_tmp5.".fasta";
        
                    open(OUTPUT_HAP_VAR_REF, ">" .$ref_file) or die "\nCan't open file $ref_file, $!\n";
                    print OUTPUT_HAP_VAR_REF ">ref\n";
                    print OUTPUT_HAP_VAR_REF $ref_part;        
                    close OUTPUT_HAP_VAR_REF;
                    
                    my $last_1000_file = $TMP_directory."hap_var_last_1000_tmp_".$id_tmp5.".fasta";
        
                    open(OUTPUT_LAST1000, ">" .$last_1000_file) or die "\nCan't open file $last_1000_file, $!\n";
                    print OUTPUT_LAST1000 ">last1000\n";
                    print OUTPUT_LAST1000 $last_1000;        
                    close OUTPUT_LAST1000;
                
                    chomp($last_1000_file);
                    chomp($ref_file);
                    my $command = "blastn -query ".$last_1000_file." -subject ".$ref_file." -out ".$TMP_directory."blast_hap_var_tmp_".$id."_".$y."_".$id_tmp5.".txt -outfmt 7 -strand plus -qcov_hsp_perc 99 -word_size 60";
                    syscmd($command);
                    $compare_haps = "";
                    $check_var_between_hap = '0';
                    $compare_haps_stop = "yes";
					my $N_count = $last_1000 =~ tr/N/N/;
					my $extra_tmp = ($N_count/length($last_1000))*100;
					
					my $input_BLAST_tmp = $TMP_directory."blast_hap_var_tmp_".$id."_".$y."_".$id_tmp5.".txt";
					
					if (-s $input_BLAST_tmp)
					{
						open(BLAST_VAR_COMP, $input_BLAST_tmp) or die "\nCan't open file $input_BLAST_tmp, $!\n";
						my $count_lines_tmp = '1';
						
						while (my $line_tmp = <BLAST_VAR_COMP>)
                        {
                            chomp($line_tmp);
							if ($count_lines_tmp eq '4' && $line_tmp eq "# 0 hits found")
                            {
                                last;
                            }
							elsif ($count_lines_tmp eq '5' && $line_tmp ne "# 1 hits found" && length($ref_part) < 35000)
                            {
                                last;
                            }
                            elsif ($count_lines_tmp > 5)
                            {
                                my @line_tmp = split /\t/, $line_tmp;
                                my $accuracy_tmp = $line_tmp[2];
								#my $alignment_length = $line_tmp[3];
                                my $ref_pos_start_tmp = $line_tmp[8];
                                my $ref_pos_end_tmp = $line_tmp[9];

                                if ($accuracy_tmp >= 99.2-$extra_tmp)
                                {						
									if ($check_last_10000 eq "")
									{
										$check_last_10000 = "yes";
										if (length($ref_part) < 35000)
										{
											$last_1000 = substr $seed{$id_tmp5}, -length($ref_part);
										}
										else
										{
											$last_1000 = substr $seed{$id_tmp5}, -35000;
										}
										print {$filehandle{$seed_id2}} $ref_pos_end_tmp." HAPS_REFOUND\n";
										unlink $input_BLAST_tmp;
										
										goto CHECK_LAST_10000;
									}
									else
									{
										print {$filehandle{$seed_id2}} $ref_pos_end_tmp." HAPS_REFOUND2\n";
										undef %prev_position_hap_compare;
										$hap_compare_mismatch_extend = '0';
										$compare_haps_stop = "";
										foreach my $id_tmp2 (keys %seed)
										{
											if ($id_tmp2 eq $id_tmp5)
											{
												$hap_compare_pos{$id_tmp2} = length($seed{$id_tmp2}) - 35000;
											}
											else
											{
												my $current_pos_tmp = $hap_compare_pos{$id_tmp2};
												$hap_compare_pos{$id_tmp2} = $current_pos_tmp+$ref_pos_end_tmp - 35000;
											}
										}	
									}
                                }  
                            }
                            $count_lines_tmp++;
                        }
                        close BLAST_VAR_COMP;
					}
					unlink $input_BLAST_tmp;
                }
				$prev_position_hap_compare{$seed_id} = $position
            }
            
            if (($check_var_between_hap > 1 || $compare_haps ne "") && $find_haps_in_seed eq "")
            {  
                my %haps;
                undef %haps;
                my $haps = "";
                my $l = '200';
                foreach my $id_tmp (keys %seed)
                {
                    if ($hap_compare_pos{$id_tmp} < 200)
                    {
                        $l = 0;
                    }
                    my $seq = substr $seed{$id_tmp}, $hap_compare_pos{$id_tmp}-$l, $var_length+300 + $hap_compare_mismatch_extend + $extra_seq;
                    my $h = '80';
                    my $last_80 = substr $seq, -$h, 80;
                    my $AT_rich_tmp = AT_rich_test ($last_80, '15');
                    
                    while (length($seq) > $h+50 && $AT_rich_tmp eq "yes")
                    {
                        $h +=80;
                        $last_80 = substr $seq, -$h, 80;
                        $AT_rich_tmp = AT_rich_test ($last_80, '17');
                    }
                    
                    if ($h > 80)
                    {
                        substr $seq, -$h-100, $h+100, "";
                        print {$filehandle{$seed_id2}} $h." HAPS_AT_REMOVE\n"; 
                    }

                    if (length($seq) < 400)
                    {
                        goto SKIP_HAP_COMPARE;
                    }
                    
                    $haps{$id_tmp} = $seq;
                    $haps .= ">".$id_tmp."\n";
                    $haps .= $seq."\n";
                }
                
                my $output_file1  = $TMP_directory."haps.fasta";
                
                open(OUTPUT_LONG1, ">" .$output_file1) or die "\nCan't open file $output_file1, $!\n";
                print OUTPUT_LONG1 $haps;
                
                close OUTPUT_LONG1;
                
                chomp($output_file1);
                                
                #my $cmd = sprintf("blastn -query %s -subject %s -out blast_tmp3.txt -outfmt 4", $ref_file, $output_file1);
                #my $cmd = sprintf("blastn -query %s -subject %s -out blast_tmp3.txt -reward 1 -penalty -2 -gapopen 2 -gapextend 2 -outfmt 4", $ref_file, $output_file1);
                #my $cmd = sprintf("muscle -in %s -out blast_tmp3.txt -maxiters 1 -diags", $output_file1);
                my $cmd = "";
MAFFT_HAPS:                
                if ($mafft_extra_quality_haps eq "yes")
                {
                    $cmd = sprintf("mafft --thread 4 --op 0.1 --ep 0.1 --quiet --clustalout --maxiterate 100 --globalpair %s > ".$TMP_directory."hap_var_".$y.".txt ", $output_file1);
                }
                else
                {
                    $cmd = sprintf("mafft --thread 4 --op 0.1 --ep 0.1 --quiet --clustalout %s > ".$TMP_directory."hap_var_".$y.".txt ", $output_file1);
                }
                
                my $system_result = syscmd($cmd);                         
                
                my $input_BLAST_var  = $TMP_directory."hap_var_".$y.".txt";
                open(INPUT_BLAST_VAR, $input_BLAST_var) or print "\n\nCan't open blast haplotypes file $input_BLAST_var, $!\n";
                
#merge mafft lines-------------------------------------------------------------------------------             
                my $g = '0';
                my $query_line = "";
                my %subject_list;
                undef %subject_list;
                my $consensus_total = "";
                my %gaps_id;
                undef %gaps_id;
                my %length_id;
                undef %length_id;
                
INPUT_BLAST_VAR:while (my $line2 = <INPUT_BLAST_VAR>)
                {                                                     
                    chomp($line2);
                    if ($g > 2)
                    {
                        my @blast_result_tmp = split /\s+/, $line2;
                        
                        if (exists($seed{$blast_result_tmp[0]}))
                        {
                            my $subject_tmp = $subject_list{$blast_result_tmp[0]};
                            $subject_list{$blast_result_tmp[0]} = $subject_tmp.$blast_result_tmp[1];
                            $query_line = "yes";
                        }  
                        elsif ($query_line eq "yes")
                        {
                            my $consensus = substr $line2, 16, 60;
                            $consensus_total .= $consensus;
                        }
                    }
                    $g++
                }
                close INPUT_BLAST_VAR;
                foreach my $subject_id (keys %subject_list)
                {
                    my $gaps_tmp = $subject_list{$subject_id} =~ tr/-/-/;
                    $gaps_id{$subject_id} = $gaps_tmp;
                    $length_id{$subject_id} = length($subject_list{$subject_id});
                }
#------------------------------------------------------------------------------------------------------
                
                chomp($consensus_total);
                $consensus_total =~ s/\s+$//;
                my @consensus = split //, $consensus_total;
                
                my $count_stars_total = $consensus_total =~ tr/*/*/;
                my $non_stars_total = length($consensus_total)-$count_stars_total;
                print {$filehandle{$seed_id2}} $non_stars_total." NON_STARS\n";
                my $mm_count = '60';
                if ($PB_reads eq "")
                {
                    $mm_count += 0.05*length($consensus_total);
                }
                if ($non_stars_total > $mm_count && $mafft_extra_quality_haps eq "")
                {
                    $mafft_extra_quality_haps = "yes";
                    goto MAFFT_HAPS;
                }
                
                my $count_non_stars = '0';
                my $c = '0';
                my $skip_pos = '0';
                my $last_indel_pos = "";
                my $mismatch_detect = "";
                
                my $c2 = 100;
                my $last_100_cons = substr $consensus_total, -$c2, 100;
                my $count_stars_last100 = $last_100_cons =~ tr/*/*/;
                my $non_stars_last100 = length($last_100_cons)-$count_stars_last100;
                my $max_length_compare = 0;
                
                while ($non_stars_last100 > 15)
                { 
                    $c2 += 100;
                    $last_100_cons = substr $consensus_total, -$c2, 100;
                    $count_stars_last100 = $last_100_cons =~ tr/*/*/;
                    $non_stars_last100 = length($last_100_cons)-$count_stars_last100;
                    $max_length_compare += 100;     
                }
                foreach my $length_id_tmp (keys %length_id)
                {
                    if ($length_id{$length_id_tmp}-$max_length_compare < 500)
                    {
                        if ($PB_reads ne "")
                        {
                            $hap_compare_mismatch_extend += '2000';
                        }
                        else
                        {
                            $hap_compare_mismatch_extend += '3000';
                        }
                        goto SKIP_HAP_COMPARE;
                    }
                }
                if ($max_length_compare > 0)
                {
                    $max_length_compare += 200;
                }
                
CONSENSUS_HAP:  foreach my $cons (@consensus)
                {
					foreach my $gaps_id (keys %gaps_id)
                    {
                        if ($c > $length_id{$gaps_id}-$max_length_compare)
                        {       
							last CONSENSUS_HAP;
                        }
                    }
                    if ($skip_pos > 0)
                    {
                        $skip_pos--;
                        $c++;
                        next CONSENSUS_HAP;
                    }
                    
                    if ($cons ne "*")
                    {
                        $count_non_stars++;
                    }
                    
                    if ($count_non_stars > 10 && $cons ne "*" && $c > $last_indel_pos+50 && $compare_haps ne "yes2")
                    {
                        my $last_100 = substr $consensus_total, $c, 100;
                        my $count_star = $last_100 =~ tr/\*/\*/;
                        my $star_limit = '82';
                        if ($PB_reads ne "")
                        {
                            $star_limit = '90';
                        }
                        if ($count_star < $star_limit)
                        {
                            $mismatch_detect = "yes";
                            #goto END1;
                        }
                        else
                        {
                            $mismatch_detect = "";
                        }
                    }
#check SNR------------------------------                            
                    my $count_SNR = "";
                    foreach my $ext_line_tmp (sort {$a <=> $b} keys %subject_list)
                    {
                        my $next_10 = substr $subject_list{$ext_line_tmp}, $c+1, 10;
                        my $t_tmp = $next_10 =~ tr/t/t/;
                        my $a_tmp = $next_10 =~ tr/a/a/;
                        my $c_tmp = $next_10 =~ tr/c/c/;
                        my $g_tmp = $next_10 =~ tr/g/g/;
                        if ($t_tmp > 8 || $a_tmp > 8 || $c_tmp > 8 || $g_tmp > 8)
                        {
                            $count_SNR = "yes";
                            last;
                        }
                    }
#check SNR------------------------------  
                    my $indel_check0 = "";
                    foreach my $subject_tmp (keys %subject_list)
                    {
                        my $seq_tmp = $subject_list{$subject_tmp};
                        my $nuc_match = substr $seq_tmp, $c, 1;
                        if ($nuc_match eq "-")
                        {
                            $indel_check0 = "yes";
                        }
                    }
                    if ($cons eq " " && $indel_check0 eq "yes" && $PB_reads ne "" && $c > 4 && $count_SNR eq "")
                    {
                        my $c2 = $c;
                        while ($consensus[$c2] eq " ")
                        {
                            $c2++;
                        }
                        if ($c2-$c >= 30)
                        {
                            print {$filehandle{$seed_id2}} $c2-$c." C2-C\n";
                            foreach my $subject_tmp (sort keys %subject_list)
                            {
                                my $seq_tmp = $subject_list{$subject_tmp};
                                my $var_region = substr $seq_tmp, $c, $c2-$c;
                                my @var_region = split //, $var_region;
                                my $f = '0';
                                foreach my $var_nuc (@var_region)
                                {
                                    if ($var_nuc eq "-")
                                    {
                                        $f++;
                                    }
                                }
                                if ($f >= 30)
                                {     
                                    my $gaps_correct_seq = substr $subject_list{$subject_tmp}, 0, $c;
                                    my $gaps_correct = $gaps_correct_seq =~ tr/-/-/;
                                    my $pos_tmp = $hap_compare_pos{$subject_tmp}-$l-$gaps_correct+$c;
                                    print OUTPUT12 $subject_tmp."\t".$pos_tmp."\tDEL\t".$f."\n";
                                    $skip_pos = $f;
                                    $c++;
                                    $last_indel_pos = $c+$skip_pos;
                                    next CONSENSUS_HAP;
                                }
                                elsif ($f > 0)
                                {
                                    goto END1;
                                }
                            }
                        }
#detect and store indels----------------------------------------------------------------------------------------
                        elsif ($c > 4)
                        {
                            my $f_highest = '0';
                            my %read_ends;
                            my $shortest_read_end = "";
                            my $shortest_read_end_id = "";
                            foreach my $subject_tmp (sort keys %subject_list)
                            {
                                my $seq_tmp = $subject_list{$subject_tmp};
                                my $var_region = substr $seq_tmp, $c, $c2-$c;
                                my @var_region = split //, $var_region;
                                my $f = '0';
                                foreach my $var_nuc (@var_region)
                                {
                                    if ($var_nuc eq "-")
                                    {
                                        $f++;
                                    }
                                }
                                if ($f > $f_highest)
                                {
                                    $f_highest = $f;
                                }

                                my $read_end_tmpi = substr $seq_tmp, $c-3, $c2-$c+3+3;
                                $read_end_tmpi =~ tr/actgn/ACTGN/;
                                $read_end_tmpi =~ tr/-//d;
                                $read_ends{$subject_tmp} = $read_end_tmpi;
                                if (length($read_end_tmpi) < length($shortest_read_end) || $shortest_read_end eq "")
                                {
                                    $shortest_read_end = $read_end_tmpi;
                                    $shortest_read_end_id = $subject_tmp;
                                }
                            }
                            my $add_length_tmp = '0';
EXTEND_INDEL:
                            foreach my $subject_tmp (sort keys %read_ends)
                            {
                                if ($read_ends{$subject_tmp} ne $shortest_read_end)
                                {
                                    my $check_seq = substr $read_ends{$subject_tmp}, 0, length($shortest_read_end);
                                    if ($check_seq eq $shortest_read_end)
                                    {
                                        $add_length_tmp++;
                                        foreach my $subject_tmp (sort keys %subject_list)
                                        {
                                            my $read_end_tmpi = substr $subject_list{$subject_tmp}, $c-3, $c2-$c+3+3+$add_length_tmp;
                                            $read_end_tmpi =~ tr/actgn/ACTGN/;
                                            $read_end_tmpi =~ tr/-//d;
                                            $read_ends{$subject_tmp} = $read_end_tmpi;
                                            if ($subject_tmp eq $shortest_read_end_id)
                                            {
                                                $shortest_read_end = $read_end_tmpi;
                                            }
                                        }
                                        goto EXTEND_INDEL;
                                    }
                                    else
                                    {
                                        foreach my $subject_tmp (sort keys %subject_list)
                                        {
                                            my $read_end_tmpi = substr $subject_list{$subject_tmp}, $c-3, $c2-$c+3+3+$add_length_tmp+2;
                                            $read_end_tmpi =~ tr/actgn/ACTGN/;
                                            $read_end_tmpi =~ tr/-//d;
                                            $read_ends{$subject_tmp} = $read_end_tmpi;
                                            if ($subject_tmp eq $shortest_read_end_id)
                                            {
                                                $shortest_read_end = $read_end_tmpi;
                                            }
                                        }
                                    }
                                }
                            }
                            my $N_check = "";
                            foreach my $subject_tmp (sort keys %read_ends)
                            {
                                my $N_tmp = $read_ends{$subject_tmp} =~ tr/N/N/;
                                if ($N_tmp > 0)
                                {
                                    $N_check = "yes";
                                }
                            }
                            foreach my $subject_tmp (sort keys %read_ends)
                            {
                                my $gaps_correct_seq = substr $subject_list{$subject_tmp}, 0, $c;
                                my $gaps_correct = $gaps_correct_seq =~ tr/-/-/;
                                my $pos_tmp = $hap_compare_pos{$subject_tmp}-$l-$gaps_correct+$c;
                                my $read_end_tmpi = $read_ends{$subject_tmp};
                                
                                if (exists($split_positions{$subject_tmp}{$pos_tmp}))
                                {}
                                elsif ($N_check eq "")
                                {
                                    $split_positions{$subject_tmp}{$pos_tmp} = $read_end_tmpi.","."INDEL";
                                    print {$filehandle{$seed_id2}} $subject_tmp."\t".$pos_tmp."\t".$read_end_tmpi.",INDEL\n";
                                    my $seq_check = substr $subject_list{$subject_tmp}, 0, $c;
                                    print OUTPUT18 $seq_check." CHECK\n";
                                    print OUTPUT18 $subject_tmp."\t".$pos_tmp."\t".$read_end_tmpi.","."INDEL\n";
                                }
                                else
                                {
                                    print OUTPUT18 $subject_tmp."\t".$pos_tmp."\t".$read_end_tmpi.","."INDEL N!\n";
                                }
                            }
                            
                            $skip_pos = $f_highest;
                            $last_indel_pos = $c+$skip_pos;
                            $c++;
                            next CONSENSUS_HAP;
                        }
                    }
#Check SNPs-------------------------------------------------------------------------------------------------------------------------------------------
                    elsif (($cons eq "." || $cons eq " ") && $mismatch_detect eq "" && $c > 4 && (($consensus[$c-1] ne " " && $consensus[$c-2] ne " " && $consensus[$c-3] ne " " && $consensus[$c-4] ne " "
                        && $consensus[$c+1] ne " " && $consensus[$c+2] ne " ") || ($PB_reads ne "" && $consensus[$c-1] ne " ")))
                    {
                        my %pos_haps;
                        undef %pos_haps;

                        foreach my $subject_tmp (keys %subject_list)
                        {
                            my $seq_tmp = $subject_list{$subject_tmp};
                            my $nuc_match = substr $seq_tmp, $c, 1;
                            my $d = $c-$overlap;
                            my $overlap_tmp = $overlap;
                            if ($d < 0)
                            {
                               $overlap_tmp -= $d;
                               $d = '0';                         
                            }
                            my $read_end_tmpi = substr $seq_tmp, $d, $overlap_tmp;
                            my $check_del = $read_end_tmpi =~ tr/-/-/;

                            if ($nuc_match eq "n" || $nuc_match eq '-' || $check_del > '0')
                            {
                                #if ($check_del > '0' && $nuc_match ne "n")
                                #{
                                    #my $gaps_correct_seq = substr $subject_list{$subject_tmp}, 0, $c;
                                   # my $gaps_correct = $gaps_correct_seq =~ tr/-/-/;
                                   # my $pos = $hap_compare_pos{$subject_tmp}-200-$gaps_correct+$c;
                                    #print OUTPUT18 $subject_tmp."\t".$pos."\t".$read_end_tmpi.",".$nuc_match."\n";
                                #}
                                goto SKIP_CONSENSUS_HAP;
                            }
                            my $gaps_correct_seq = substr $subject_list{$subject_tmp}, 0, $c;
                            my $gaps_correct = $gaps_correct_seq =~ tr/-/-/;
                            my $pos = $hap_compare_pos{$subject_tmp}-$l-$gaps_correct+$c+1;
                            $pos_haps{$subject_tmp} = $pos;
                            my $one = $pos-1;
                            my $two = $pos+1;

                            if (exists($quality_scores{$subject_tmp}{$pos}) && $PB_reads eq "")
                            {                          
                                my @q_score_gap_tmp = split / /, $quality_scores_gap{$subject_tmp}{$pos}; 
								if (exists($quality_scores_gap{$subject_tmp}{$pos}) && $q_score_gap_tmp[0] < 0.7 && $compare_haps ne "yes2")
                                {
                                    print {$filehandle{$seed_id2}} $pos."\t".$quality_scores_gap{$subject_tmp}{$pos}." 1\n";
                                    goto SKIP_CONSENSUS_HAP;
                                }
								my @q_score_gap1_tmp = split / /, $quality_scores_gap{$subject_tmp}{$one}; 
                                if (exists($quality_scores_gap{$subject_tmp}{$one}) && $q_score_gap1_tmp[0] < 0.7)
                                {
                                    print {$filehandle{$seed_id2}} $pos."\t".$quality_scores_gap{$subject_tmp}{$one}." 2\n";
                                    goto SKIP_CONSENSUS_HAP;
                                }
								my @q_score_gap2_tmp = split / /, $quality_scores_gap{$subject_tmp}{$two}; 
                                if (exists($quality_scores_gap{$subject_tmp}{$two}) && $q_score_gap2_tmp[0] < 0.7)
                                {
                                    print {$filehandle{$seed_id2}} $pos."\t".$quality_scores_gap{$subject_tmp}{$two}." 3\n";
                                    goto SKIP_CONSENSUS_HAP;
                                }
								my @q_score_tmp = split / /, $quality_scores{$subject_tmp}{$pos}; 
                                if ($q_score_tmp[0] > 0.75 || $compare_haps eq "yes2")
                                {}
                                else
                                {                            
                                    goto SKIP_CONSENSUS_HAP;
                                }
                            }
                            elsif ($PB_reads eq "" && $compare_haps ne "yes2")
                            {
                                print {$filehandle{$seed_id2}} $subject_tmp." ".$pos." 4\n";
                                goto SKIP_CONSENSUS_HAP;
                            }
							if ($PB_reads eq "" && $compare_haps ne "yes2")
							{
								my $t_tmp = $read_end_tmpi =~ tr/TtnN/TtnN/;
                                my $a_tmp = $read_end_tmpi =~ tr/AanN/AanN/;
                                my $c_tmp = $read_end_tmpi =~ tr/CcnN/CcnN/;
                                my $g_tmp = $read_end_tmpi =~ tr/GgnN/GgnN/;
                                if ($t_tmp >= length($read_end_tmpi)-1 || $a_tmp >= length($read_end_tmpi)-1 || $c_tmp >= length($read_end_tmpi)-1 || $g_tmp >= length($read_end_tmpi)-1)
                                {
                                    goto SKIP_CONSENSUS_HAP;
                                }
							}
                        }
                        foreach my $subject_tmp (keys %subject_list)
                        {
                            my $seq_tmp = $subject_list{$subject_tmp};
                            my $nuc_match = substr $seq_tmp, $c, 1;
                            $nuc_match =~ tr/actgn/ACTGN/;
                            my $d = $c-$overlap;
                            my $overlap_tmp = $overlap;
                            if ($d < 0)
                            {
                               $overlap_tmp -= $d;
                               $d = '0';                         
                            }
                            my $read_end_tmpi = substr $seq_tmp, $d, $overlap_tmp;
                            $read_end_tmpi =~ tr/-//d;
                            $read_end_tmpi =~ tr/actgn/ACTGN/;
                            my $pos = $pos_haps{$subject_tmp};
                            
                            print {$filehandle{$seed_id2}} $pos." QUAL_POS\n";
                            my $one = $pos-1;
                            my $two = $pos+1;
                            my $three = $pos-2;
                            my $four = $pos+2;
                            print {$filehandle{$seed_id2}} $quality_scores{$subject_tmp}{$three}." QUAL-2\n";
                            print {$filehandle{$seed_id2}} $quality_scores{$subject_tmp}{$one}." QUAL-1\n";
                            print {$filehandle{$seed_id2}} $quality_scores{$subject_tmp}{$pos}." QUAL0\n";
                            print {$filehandle{$seed_id2}} $quality_scores{$subject_tmp}{$two}." QUAL+1\n";
                            print {$filehandle{$seed_id2}} $quality_scores{$subject_tmp}{$four}." QUAL+2\n";
     
                            if (exists($split_positions{$subject_tmp}{$pos}))
                            {}
                            else
                            {
                                $split_positions{$subject_tmp}{$pos} = $read_end_tmpi.",".$nuc_match;
                                print {$filehandle{$seed_id2}} $subject_tmp."\t".$pos."\t".$read_end_tmpi.",".$nuc_match."\n";
                                my $seq_check = substr $subject_list{$subject_tmp}, 0, $c;
                                print OUTPUT18 $seq_check." CHECK\n";
                                print OUTPUT18 $subject_tmp."\t".$pos."\t".$read_end_tmpi.",".$nuc_match."\n";
                            }           
                        }
                    }
SKIP_CONSENSUS_HAP:
                    if ($cons eq " " && $PB_reads eq "ddg")
                    {
                        my $n = "";
                        my $no_n = "";
                        my $new_nuc = "";
                        foreach my $subject_tmp (keys %subject_list)
                        {
                            my $seq_tmp = $subject_list{$subject_tmp};
                            my $nuc_match = substr $seq_tmp, $c, 1;
                            
                            if ($nuc_match eq "n")
                            {
                                $n = $subject_tmp;
                            }
                            if ($nuc_match eq "a" || $nuc_match eq "c" || $nuc_match eq "t" ||$nuc_match eq "g")
                            {
                                $no_n = $subject_tmp;
                                $new_nuc = $nuc_match;
                            }
                        }
                        if ($n ne "" && $no_n ne "")
                        {
                            my $gaps_correct_seq = substr $subject_list{$n}, 0, $c;
                            my $gaps_correct = $gaps_correct_seq =~ tr/-/-/;
                            my $pos = $hap_compare_pos{$n}-$l-$gaps_correct+$c;
                            my $seq_tmp = $seed{$n};
                            $new_nuc =~ tr/actg/ACTG/;         
                            substr $seq_tmp, $pos, 1, $new_nuc;
                            $seed{$n} = $seq_tmp;
                            if ($n eq $id)
                            {
                                $read = $seq_tmp;
                            }  
                        }
                    }
                    $c++;
                }

                $hap_compare_mismatch_extend = '0';
                foreach my $subject_list_tmp (keys %subject_list)
                {
                    my $pos_tmp = $hap_compare_pos{$subject_list_tmp};
                    my $gaps_correct_seq = substr $subject_list{$subject_list_tmp}, 0, $c;
                    my $gaps_correct = $gaps_correct_seq =~ tr/-/-/;
                    $hap_compare_pos{$subject_list_tmp} = $pos_tmp + $c - $gaps_correct;
                    print {$filehandle{$seed_id2}} $subject_list_tmp." POS_CORRECT_ID\n"; 
                    print {$filehandle{$seed_id2}} $hap_compare_pos{$subject_list_tmp}." POS_CORRECT\n"; 
                }  
            }
        }
        
#Check for repetitive region------------------------------------------------------  
			my $vl = '50000';
REP_TEST_NP0:		
		    my $end_repetitive = substr $read, -$vl;
            my $end_repetitiveb = substr $end_repetitive, -5000;
			
            $repetitive_detect2 = "";
            $repetitive_detect1 = "";
            my $u = '12';
REP_TEST_NP: while ($u < 250)
            { 
                my $repetitive_test2 = substr $end_repetitive, -$u, 12;
                my $repetitive_test2c = substr $end_repetitive, -$u-80, 92;
                my $SNR_checkA = $repetitive_test2 =~ tr/A/A/;
                my $SNR_checkC = $repetitive_test2 =~ tr/C/C/;
                my $SNR_checkT = $repetitive_test2 =~ tr/T/T/;
                my $SNR_checkG = $repetitive_test2 =~ tr/G/G/;
                my $SNR_checkN = $repetitive_test2 =~ tr/N/N/;
                if ($SNR_checkA > 8 || $SNR_checkC > 8 || $SNR_checkG > 8 || $SNR_checkT > 8 || $SNR_checkN > 0)
                {
                    $u += 50;
                    goto REP_TEST_NP;
                }
                else
                {
                    my $check_repetitive2 = $end_repetitiveb =~ s/$repetitive_test2/$repetitive_test2/g;
                    my $check_repetitive2c = $end_repetitive =~ s/$repetitive_test2c/$repetitive_test2c/g;
                    
                    if ($check_repetitive2 > 2)
                    {
                        $repetitive_detect1 = "yes";
                        print {$filehandle{$seed_id2}} $check_repetitive2." REP_DETECT1\n";
                    }

                    if ($check_repetitive2c > 1)
                    {
                        my @rep_split = split /$repetitive_test2c/, $end_repetitive;
                        if ($repetitive_detect2 eq "" || length($rep_split[0]) > $vl-$repetitive_detect2)
                        {
                            $repetitive_detect2 = $vl-length($rep_split[0]);
                        }
                        print {$filehandle{$seed_id2}} $repetitive_detect2." REP_DETECT3\n";
                    }

                    if ($repetitive_detect1 eq "" && $repetitive_detect2 eq "")
                    {
                        $u += 50;
                    }
                    elsif ($repetitive_detect2 ne "" && $u < 1000)
                    {
                        if ($check_repetitive2c > 1 && length($end_repetitive) < 200000 && length($read) > 200000)
						{
							$vl = '200000';
							goto REP_TEST_NP0;
						}
						elsif ($check_repetitive2c > 10 && length($end_repetitive) > 190000 && $repetitive_detect2 > 125000)
						{
							print {$filehandle{$seed_id2}} $repetitive_detect2." END_STUCK_IN_REP\n";
							$unresolvable_NP = "yes";
							$unresolvable_PB = "yes";
							goto AFTER_EXT;
						}
						last REP_TEST_NP;
                    }
                    else
                    {
                        last REP_TEST_NP;
                    }
                }
            }        
    
SKIP_HAP_COMPARE:                 
          
#Check if all databeses are created--------------------------------------------------------------------------------------------------------------------------------------------------------------      
                    
            if ($y0 eq '1')
            {
CHECK_DBs_NP:              
                my $time_before_check = time;
                if ($input_reads_DB_folder_NP eq "")
                {
                    my $count_DB_checked = '0';
                    foreach my $pid_tmp (keys %ret_data2_NP)
                    {
                        if (-e $ret_data2_NP{$pid_tmp})
                        {                 
                            open(DB_CHECK_NP, $ret_data2_NP{$pid_tmp}) or die "\nCan't open file $ret_data2_NP{$pid_tmp}, $!\n";

                            while (my $line = <DB_CHECK_NP>)
                            {
                                my $check_end = substr $line, 0, 34;
                                if ($check_end eq "Adding sequences from FASTA; added")
                                {
                                    $count_DB_checked++;
                                }                
                            }
                            close DB_CHECK_NP;
                        }      
                    }
                    if (time > $time_before_check+3600)
                    {
                        print "ERROR: Can't find all the databese output files\n";
                        goto END1;
                    }
                    elsif ($count_DB_checked < keys %ret_data2_NP)
                    {
                        goto CHECK_DBs_NP;
                    }
                }
                
CHECK_DBs_PB:              
                my $time_before_check2 = time;
                if ($input_reads_DB_folder_PB eq "")
                {
                    my $count_DB_checked = '0';
                    foreach my $pid_tmp (keys %ret_data2_PB)
                    {
                        if (-e $ret_data2_PB{$pid_tmp})
                        {                 
                            open(DB_CHECK_PB, $ret_data2_PB{$pid_tmp}) or die "\nCan't open file $ret_data2_PB{$pid_tmp}, $!\n";

                            while (my $line = <DB_CHECK_PB>)
                            {
                                my $check_end = substr $line, 0, 34;
                                if ($check_end eq "Adding sequences from FASTA; added")
                                {
                                    $count_DB_checked++;
                                }                
                            }
                            close DB_CHECK_PB;
                        }      
                    }
                    if (time > $time_before_check+3600)
                    {
                        print "ERROR: Can't find all the databese output files\n";
                        goto END1;
                    }
                    elsif ($count_DB_checked < keys %ret_data2_PB)
                    {
                        goto CHECK_DBs_PB;
                    }
                }
                
                foreach my $del_file_tmp (keys %delete_input_files)
                {
                    unlink  $del_file_tmp;
                }
				
				foreach my $hash_tmp (@QUALITY_HASH)
				{
					unlink $hash_tmp;
				}
            }
			
			my $time_START = time;
            my $time0b = $time_START - $time_START0;
            print {$filehandle{$seed_id2}} $time0b." TIME0\n";
PB_READS:
     
NP_READS:
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Check NP reads----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            if (($NP_reads ne "" || $input_reads_DB_folder_NP ne "") && length($read) > 300 && ($NP_reads_support eq "yes2" || ($PB_reads eq "" && $input_reads_DB_folder_PB eq "")))
            {     
                my $first_length_back = '900';
				my $flexible_seed = "";

                if ($first_length_back > length($read))
                {
                    $first_length_back = length($read)-50;
                }
                
                if ($find_haps_in_seed eq "yes")
                {
                    $first_length_back = length($read);
					if ($assembly_length_max eq "WG" && $first_length_back > 5100 && $first_back_assembly eq "")
					{
						$read = substr $read, 100, 5000;
						$seed{$id} = $read;
						$first_length_back = 5000;
					}
                }

                my $discontiguous_blast = "";
                my $last_hap_pos = "";
				my $last_hap_pos_DUP = "";
                my %hash_NP_reads_tmp;
                undef %hash_NP_reads_tmp;
				my %hash_NP_reads_tmp_back;
                undef %hash_NP_reads_tmp_back;
				my $query_coverage = '90';
				my $query_accuracy = '70';
				
SKIP_CONFIRMED_NP:
                               
                my %id_matches;
                undef %id_matches;
                my %double_matches;
                undef %double_matches;
                my %reverse_list;
                undef %reverse_list;
				my $only_saved_reads = "";
				my $new_contig_check = "";
				
				if ($y eq "1")
				{
					foreach my $id_read (keys %{$split_contigs_reads2{$id}})
					{
						$id_matches{$id_read} = undef;
						if ($split_contigs_reads2{$id}{$id_read} eq "yes")
						{
							$reverse_list{$id_read} = "yes";
						}
						$new_contig_check = "yes";
					}
					delete $split_contigs_reads2{$id};
				}

#Find matches in hash against last 600 bp--------------------------------------------------------------                           
	
LONGER_LAST_600_NP:

				my $last_600 = substr $read, -$first_length_back;
				my %first_pos_read;
                undef %first_pos_read;
                my $time_NP_find_matches = time;
				
				my %output_files_DB;
                undef %output_files_DB;
                my %output_files_DB2;
                undef %output_files_DB2;
				
				if ($new_contig_check eq "yes")
				{
					goto SKIP_TO_CONFIRMED_NP;
				}

				if (keys %{$trace_back_split_NP{$id}} > 0)
					{
							print {$filehandle{$seed_id2}} $last_600." LAST_600!!!!!!!!!!!!!!!!\n";         
					}
                
                my $query_file_DB = $TMP_directory."query_NP.fasta";
                open(INPUT_QUERY_NP, ">" .$query_file_DB) or die "\nCan't open file $query_file_DB, $!\n";
                #print INPUT_QUERY_NP ">ref\n";
                print INPUT_QUERY_NP $last_600;        
                close INPUT_QUERY_NP;       
                
                foreach my $pid_tmp (keys %ret_data_NP)
                {
                    my $file_tmp = $TMP_directory."blast_tmp_DB_NP_".$id."_".$y."_".$pid_tmp.".txt";
                    my $command_DB = "blastn -query ".$query_file_DB." -db ".$ret_data_NP{$pid_tmp}." -out ".$file_tmp." -outfmt 7 -qcov_hsp_perc 90 -num_threads 2 &";
                    syscmd($command_DB);             
                    $output_files_DB{$file_tmp} = undef;
                    $output_files_DB2{$file_tmp} = undef;
                }
				
#Add saved reads--------------------------------------------------------------                           						

SKIP_TO_CONFIRMED_NP:
				my %extensions;
                undef %extensions;
                my %extensions2;
                undef %extensions2;
                my %extensions2b;
                undef %extensions2b;
                my %extensions_nomatch;
                undef %extensions_nomatch;
                my %extensions_unknown;
                undef %extensions_unknown;
                my %extensions_nomatch2;
                undef %extensions_nomatch2;
                my %extensions_nomatch2b;
                undef %extensions_nomatch2b;                
                my %extensions_unknown2;
                undef %extensions_unknown2;
                my %store_mismatches_NP;
                undef %store_mismatches_NP;
                my %store_mismatches_all_NP;
                undef %store_mismatches_all_NP;
				my %store_mismatches_N_NP;
                undef %store_mismatches_N_NP;
				
				my $extensions_nomatch2b_count = '0';
                my $extensions_nomatch2b_count_saved = '0';
                my $total_matches_extra = '0';
                my $position_prev = "";
                my %alignment_length_saved;
                undef %alignment_length_saved;
                my %score_match_saved;
                undef %score_match_saved;
                my %score_no_match_saved;
                undef %score_no_match_saved;
                my %score_match_DUP_saved;
                undef %score_match_DUP_saved;
                my %score_no_match_DUP_saved;
                undef %score_no_match_DUP_saved;
                my %accuracy_saved;
                undef %accuracy_saved;               
                my %read_start_pos_rej;
                undef %read_start_pos_rej;
                my %read_start_pos_rej_saved;
                undef %read_start_pos_rej_saved;
                my %extensions_nomatch2b_saved;
                undef %extensions_nomatch2b_saved;
                
                foreach my $seed_id_tmp0 (keys %save_alignment_data_NP)
                {
                    if ($seed_id_tmp0 eq $seed_id)
                    {
                        foreach my $id_tmp7 (keys %{$save_alignment_data_NP{$seed_id_tmp0}})
                        {
                            my @alignment_data = split /_/, $save_alignment_data_NP{$seed_id_tmp0}{$id_tmp7};

                            if ($alignment_data[0] eq "yes" || $alignment_data[11] eq "")
                            { 
                                $alignment_length_saved{$id_tmp7} = $alignment_data[4];
                                $first_pos_read{$id_tmp7} = $alignment_data[2];
                                $position_prev = $alignment_data[3];
                                $score_match_saved{$id_tmp7} = $alignment_data[5];
                                $score_no_match_saved{$id_tmp7} = $alignment_data[6];
                                $score_match_DUP_saved{$id_tmp7} = $alignment_data[7];
                                $score_no_match_DUP_saved{$id_tmp7} = $alignment_data[8];
                                $accuracy_saved{$id_tmp7} = $alignment_data[9];
                            }
                           
                            if (($alignment_data[6] > $alignment_data[5] && $alignment_data[6] > 1) || ($alignment_data[8] > $alignment_data[7] && ($alignment_data[8] > 3
								|| $alignment_data[6]+$alignment_data[8] > $alignment_data[5]+$alignment_data[7]))
								|| ($hap_tag eq "HAP1" && exists($exclude_reads_hap1_NP{$id_tmp7})) || ($hap_tag eq "HAP2" && exists($exclude_reads_hap2_NP{$id_tmp7})))
                            {
								$extensions_nomatch2b_saved{$id_tmp7}{$alignment_data[5]}{$alignment_data[6]} = $alignment_data[3]+$alignment_data[10];
                                $extensions_nomatch2b_count_saved++;
								if ($alignment_data[1] eq "yes")
                                {
                                    $reverse_list{$id_tmp7} = undef;
                                }
                            }
							elsif ($alignment_data[0] ne "yes" && $alignment_data[11] ne "")
                            {
                                $read_start_pos_rej{$id_tmp7} = $alignment_data[11]+$position-$alignment_data[3];
                                $read_start_pos_rej_saved{$id_tmp7} = undef;
                                if ($alignment_data[1] eq "yes")
                                {
                                    $reverse_list{$id_tmp7} = undef;
                                }
                            }
                            else
                            {
								$id_matches{$id_tmp7} = undef;
                                $total_matches_extra++;
                                if ($alignment_data[1] eq "yes")
                                {
                                    $reverse_list{$id_tmp7} = undef;
                                }
                                elsif ($alignment_data[1] eq "yes2")
                                {
                                    $reverse_list{$id_tmp7} = undef;
                                    $double_matches{$id_tmp7} = undef;
                                }
								if ($alignment_data[12] ne "")
								{
									$hash_NP_reads_tmp{$id_tmp7} = $alignment_data[12];
								}
                            }				
                        }
                    }
                }                                          
                print {$filehandle{$seed_id2}} $total_matches_extra." LAST_1000_matches_EXTRA\n\n";
				
                
                my $database_count_tmp = '0';
DB_RESULTS_NP:  foreach my $blast_db_results_tmp (keys %output_files_DB)
                {
                    if (-s $blast_db_results_tmp)
                    {
						my $file_complete = "";
DB_RESULTS_NP1:                       
						open(BLAST_RESULTS_DB_NP, $blast_db_results_tmp) or die "\nCan't open file $blast_db_results_tmp, $!\n";
                        my $count_lines_tmp = '1';
						
                        while (my $line_tmp = <BLAST_RESULTS_DB_NP>)
                        {
                            chomp($line_tmp);
							if ($count_lines_tmp eq '4' && $line_tmp eq "# 0 hits found")
                            {
                                delete $output_files_DB{$blast_db_results_tmp};
                                next DB_RESULTS_NP;
                            }
							elsif ($count_lines_tmp eq '5' && $line_tmp eq "# BLAST processed 1 queries")
                            {
                                delete $output_files_DB{$blast_db_results_tmp};
                                next DB_RESULTS_NP;
                            }
                            elsif ($count_lines_tmp > 5 && $line_tmp eq "# BLAST processed 1 queries" && $file_complete eq "")
                            {
								$file_complete = "yes";
								delete $output_files_DB{$blast_db_results_tmp};
								goto DB_RESULTS_NP1;
                            }
                            elsif ($count_lines_tmp > 5 && $file_complete eq "yes")
                            {
                                my @line_tmp = split /\t/, $line_tmp;
                                my $id_tmp = $line_tmp[1];
                                my $accuracy_tmp = $line_tmp[2];
								my $alignment_length = $line_tmp[3];
                                my $read_pos_start_tmp = $line_tmp[8];
                                my $read_pos_end_tmp = $line_tmp[9];

                                if ($alignment_length < 0.8*$first_length_back && $first_length_back < 1100)
                                {
									next;
                                }
								if ($hap_tag eq "HAP1" && exists($exclude_reads_hap1_NP{$id_tmp}))
                                {
								    next;
                                }
                                if ($hap_tag eq "HAP2" && exists($exclude_reads_hap2_NP{$id_tmp}))
                                {
                                    next;
                                }
								if (exists($rejected_alignment_data_NP{$seed_id}{$id_tmp}))
                                {
									next;
                                }
                                if (exists($save_alignment_data_NP{$seed_id}{$id_tmp}))
                                {
									next;
                                }
                                if ($accuracy_tmp > $query_accuracy)
                                {
									if (exists($id_matches{$id_tmp}))
                                    {
                                        if ($read_pos_end_tmp < $read_pos_start_tmp) 
                                        {
                                            if (exists($reverse_list{$id_tmp}))
											{}
											else
											{
												$reverse_list{$id_tmp} = undef;		
											}	
                                        }
										$double_matches{$id_tmp} = undef;
                                        delete $first_pos_read{$id_tmp};	
                                    }
                                    else
                                    {
                                        $id_matches{$id_tmp} = undef;
										
										if ($assembly_length_max eq "WG")
										{
											$reads_as_seeds{$id_tmp} = undef;
											print OUTPUT19 $id_tmp."\n";
										}
						
                                        if ($read_pos_end_tmp < $read_pos_start_tmp) 
                                        {
                                            $reverse_list{$id_tmp} = undef;
                                            $first_pos_read{$id_tmp} = $read_pos_start_tmp;
                                        }
                                        else
                                        {
                                            $first_pos_read{$id_tmp} = $read_pos_end_tmp;
                                        }
                                    }        
                                }
                            }
                            $count_lines_tmp++;
                        }
                        close BLAST_RESULTS_DB_NP;
                    }
                }
                if (keys %output_files_DB > 0)
                {
                    goto DB_RESULTS_NP;
                }
                my $time_kmers = time;
                my $time1 = $time_kmers - $time_START;
                print {$filehandle{$seed_id2}} $time1." TIME1\n";
				
				my $total_matches = keys %id_matches;
				print {$filehandle{$seed_id2}} "\n".$total_matches." LAST_1000_matches\n";
				
				if ($assembly_length_max eq "WG" && $y eq "1" && ($total_matches > $sequencing_depth_NP*2 || ($total_matches < 10 && $total_matches < $sequencing_depth_NP*0.7)) && ($first_back_assembly eq "" || length($read) < 5000))
				{
					foreach my $blast_db_results_tmp (keys %output_files_DB2)
					{
						unlink $blast_db_results_tmp;
						delete $output_files_DB2{$blast_db_results_tmp};
					}
					$first_back_assembly = "yes";
					goto END1;
				}
				
				if ($total_matches < 5 && $y eq "1" && $flexible_seed eq "")
                {
					$query_coverage = '70';
					$query_accuracy = '60';
					undef %id_matches;
					undef %double_matches;
					undef %reverse_list;
					undef %hash_NP_reads_tmp;
					foreach my $blast_db_results_tmp (keys %output_files_DB2)
					{
						unlink $blast_db_results_tmp;
					}
					$flexible_seed = "yes";
                    goto SKIP_CONFIRMED_NP;
                }
				if ($total_matches < 10 && $total_matches < $sequencing_depth_NP && ($repetitive_detect1 ne "" || $repetitive_detect2 ne "") && $first_length_back < 6000)
				{
					$first_length_back += 1000;
					if ($total_matches > 1000)
					{
						$first_length_back += 1000;
					}
					$query_coverage = '60';
					undef %id_matches;
					undef %double_matches;
					undef %reverse_list;
					undef %hash_NP_reads_tmp;
					foreach my $blast_db_results_tmp (keys %output_files_DB2)
					{
						unlink $blast_db_results_tmp;
					}
					goto LONGER_LAST_600_NP;
				}
				
				if ($total_matches > 350 && $first_length_back < 2000)
				{
					$query_coverage = '90';
					$query_accuracy = '80';
					$first_length_back += 1500;
					if ($total_matches > 500)
					{
						$first_length_back += 1000;
					}
					if ($total_matches > 1000)
					{
						$first_length_back += 2000;
					}
					if ($total_matches > 2000)
					{
						$first_length_back += 3000;
					}
					undef %id_matches;
					undef %double_matches;
					undef %reverse_list;
					undef %hash_NP_reads_tmp;
					foreach my $blast_db_results_tmp (keys %output_files_DB2)
					{
						unlink $blast_db_results_tmp;
					}
					goto LONGER_LAST_600_NP;
				}
				
				if ($repetitive_detect2 ne "" && $repetitive_detect2 > $first_length_back && $total_matches > 1000 && $total_matches_extra > 9 && $only_saved_reads eq "")
				{
					$only_saved_reads = "yes";
					undef %id_matches;
					undef %double_matches;
					undef %reverse_list;
					undef %hash_NP_reads_tmp;
					undef %output_files_DB;
					undef %output_files_DB2;
					foreach my $blast_db_results_tmp (keys %output_files_DB2)
					{
						unlink $blast_db_results_tmp;
					}
					goto SKIP_TO_CONFIRMED_NP;
				}
                                                                                            
#--------------------------------------------------------------------------------------------------------------------------------------------------               
#Sort matches by alignment length------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------
		
                my $count_matches_with_high_scores = '0';
                my $longest_read = '0';              
                my %check_prev_splits;
                undef %check_prev_splits;
                my %accuracy;
                undef %accuracy;             
                my $count_limit1 = '0';         
               
                my %ref_part;
                undef %ref_part;              
                my %var_matches;
                undef %var_matches;
                my %var_matches_DUP;
                undef %var_matches_DUP;
                my %alignment_length_save;
                undef %alignment_length_save;
                my %long_read_end_pos_save;
                undef %long_read_end_pos_save;
				my %rejected_reads_save;
                undef %rejected_reads_save;
                my %multi_match;
                undef %multi_match;
                my %multi_match_pos;
                undef %multi_match_pos;
                my %span_complex_region;
                undef %span_complex_region;
                
                my %second_try;
                undef %second_try;
                my %length_ext_all;
                undef %length_ext_all;

                my $count_limit2 = '0';
                my $second_try = "";
                my $count_attempts = '0';
                
                my $ref_end_rejection = '0';
                my $read_start_pos_rejection = '0';
                my $no_hit_rejection = '0';
                my $almost_no_hit_rejection = '0';
                my $ext_too_short = '0';
                my $match_too_short = '0';
                my $overlap_too_short = '0';
                my $accuracy_rejection = '0';
                my $multi_match_rejection = '0';
                my $ref_end_rejection_threshold = '150';
                if ($discontiguous_blast ne "")
                {
                    $ref_end_rejection_threshold = '400';
                }
                my %score_matches_save;
                undef %score_matches_save;
                my $add_rejected_reads = "";
                my $add_no_match_reads = "";
                my %save_reads_for_next;
                undef %save_reads_for_next;          
                
ADD_REJ_POS_NP:
                my %id_by_length;
                undef %id_by_length;
                my %printed_refs;
                undef %printed_refs;
                my %printed_refs2;
                undef %printed_refs2;
                my %printed_blast;
                undef %printed_blast;
                my %add_rej_reads_extra;
                undef %add_rej_reads_extra;
				my %input_files_blast;
                undef %input_files_blast;
				my $count_total_matches = '0';
				my %input_length_done;
				undef %input_length_done;
				my %input_BLAST_done;
				undef %input_BLAST_done;
				my %double_match_check;
				undef %double_match_check;
				my %double_match_reject;
				undef %double_match_reject;
               
                foreach my $id_tmp (keys %id_matches)
                {        
					my $length_tmp = "";
					if (exists($first_pos_read{$id_tmp}))
					{
						$length_tmp = $first_pos_read{$id_tmp};
						if (exists($alignment_length_saved{$id_tmp}))
						{
							$length_tmp = $alignment_length_saved{$id_tmp}+$position-$position_prev;
						}
						elsif (exists($reverse_list{$id_tmp}) && exists($hash_NP_reads_tmp{$id_tmp}))
						{
							$length_tmp = length($hash_NP_reads_tmp{$id_tmp})-$first_pos_read{$id_tmp};
						}
					}
					elsif (exists($hash_NP_reads_tmp{$id_tmp}))
					{
						$length_tmp = length($hash_NP_reads_tmp{$id_tmp});
					}
					else
					{
						$length_tmp = 1;
					}
										
					if (exists($id_by_length{$length_tmp}))
					{
						while (exists($id_by_length{$length_tmp}))
						{
							$length_tmp += 1;
						}
						$id_by_length{$length_tmp} = $id_tmp;
					}
					else
					{
						$id_by_length{$length_tmp} = $id_tmp;
					}
					$count_total_matches++;
                }
#BLAST--------------------------------------------------------------------------------------------------------------------------
			
                foreach my $length_tmp (sort {$b <=> $a} keys %id_by_length)
                {                
                    if (exists($input_length_done{$length_tmp}))
					{
						next;
					}
					
					my $id_tmp2 = $id_by_length{$length_tmp};
					my $length_tmp2 = $length_tmp;
					
					if (exists($hash_NP_reads_tmp{$id_tmp2}))
					{}
					else
					{
						my @split_tmp = split /a/, $id_tmp2;
						my $output_file_NP1  = $tmp_sequences_directory_NP.$output_path_test.$split_tmp[0]."a".$output_path_test.$split_tmp[2].$output_path_test."sequence_tmp_NP_".$id_tmp2.".fasta";
						open(OUTPUT_LONG_NP4, $output_file_NP1) or die "\nCan't open file $output_file_NP1, $!\n";
						my $firstLine = <OUTPUT_LONG_NP4>;
						close OUTPUT_LONG_NP4;
						$hash_NP_reads_tmp{$id_tmp2} = $firstLine;
						$length_tmp2 = length($firstLine);
					}
					
					
                    my $double_match_tmp = "";                   
DOUBLE_MATCH1:                   
                    my $reverse_check = "";
                    my @split_tmp = split /a/, $id_tmp2;
                    my $output_file1  = $tmp_sequences_directory_NP.$output_path_test.$split_tmp[0]."a".$output_path_test.$split_tmp[2].$output_path_test."sequence_tmp_NP_".$id_tmp2.".fasta";
                    
                    if (exists($reverse_list{$id_tmp2}) && $double_match_tmp eq "")
                    {
                        $output_file1  = $tmp_sequences_directory_NP.$output_path_test.$split_tmp[0]."a".$output_path_test.$split_tmp[2].$output_path_test."sequence_tmp_NP_".$id_tmp2."_rev.fasta";
                        $reverse_check = "_rev";
                    }

                    if (exists($printed_reads_NP{$id_tmp2.$reverse_check}))
                    {
                    }
                    elsif (exists($reverse_list{$id_tmp2}) && $double_match_tmp eq "")
                    {                      
                        my $long_read_tmp = reverse($hash_NP_reads_tmp{$id_tmp2});
                        $long_read_tmp =~ tr/ACTG/TGAC/;
                        open(OUTPUT_LONG1, ">" .$output_file1) or die "\nCan't open file $output_file1, $!\n";
						OUTPUT_LONG1->autoflush(1);
                        #print OUTPUT_LONG1 ">".$id_tmp2."\n";
                        print OUTPUT_LONG1 $long_read_tmp;
                        close OUTPUT_LONG1;
                        $printed_reads_NP{$id_tmp2.$reverse_check} = undef;
                    }                               
                    
                    my $ref_part = "";
                    my $ref_file = "";                  
                    my $length_tmp_final = int(($length_tmp2*1.1)+50);
                    if ($length_tmp_final > length($read))
                    {
                        $length_tmp_final = length($read);
                    }
                    if (exists($printed_refs{$length_tmp_final}))
                    {
                        $ref_file = $printed_refs{$length_tmp_final};
                    }                
                    else
                    {
                        my $ref_part_found = "";
                        
                        if ($add_rejected_reads eq "")
                        {
                            foreach my $length_tmp3 (sort {$a <=> $b} keys %printed_refs)
                            {
                                if ($length_tmp3 >= $length_tmp_final-50 && $length_tmp3 < $length_tmp_final+1000)
                                {
                                    $ref_file = $printed_refs{$length_tmp3};
                                    $ref_part_found = "yes";
                                    $length_tmp_final = $length_tmp3;
                                    last;
                                }
                            }
                        }
                        
                        if ($ref_part_found eq "")
                        {
                            $ref_part = substr $read, -$length_tmp_final;
                            $ref_file = $TMP_directory."ref_tmp_".$id_tmp2."_".$y.".fasta";
                            $printed_refs{$length_tmp_final} = $ref_file;
                            $printed_refs2{$id_tmp2} = undef;
                            open(OUTPUT_LONG2, ">" .$ref_file) or die "\nCan't open file $ref_file, $!\n";
							OUTPUT_LONG2->autoflush(1);
                            #print OUTPUT_LONG2 ">ref\n";
                            print OUTPUT_LONG2 $ref_part;							
                            close OUTPUT_LONG2;
                        }
                    }

                    $ref_part{$id_tmp2} = $length_tmp_final;
                    $printed_blast{$id_tmp2} = $ref_file; 
                
                    chomp($output_file1);
                    chomp($ref_file);

                    $input_files_blast{$length_tmp}{$output_file1}{$ref_file}{$length_tmp_final} = $double_match_tmp;
					$input_length_done{$length_tmp} = undef;
					
					if (-s $output_file1)
					{
					}
					else
					{
						print "the file ".$output_file1." does not exist!\n";
					}
					if (-s $ref_file)
					{
					}
					else
					{
						print "the file ".$ref_file." does not exist!\n";
					}

					my $command = "blastn -query ".$ref_file." -subject ".$output_file1." -out ".$TMP_directory."blast_tmp_".$id."_".$y."_".$double_match_tmp."_".$id_tmp2.".txt -outfmt 5 -strand plus -culling_limit 1 &";
					#my $command = "blastn -task dc-megablast -query ".$ref_file." -subject ".$output_file1." -out ".$output_path."blast_tmp_".$id."_".$y."_".$double_match_tmp."_".$id_tmp2.".txt -outfmt 5 -strand plus -evalue 1e-50 -culling_limit 1 &";
					if ($discontiguous_blast eq "yes")
					{
						$command = "blastn -task dc-megablast -query ".$ref_file." -subject ".$output_file1." -out ".$TMP_directory."blast_tmp_".$id."_".$y."_".$double_match_tmp."_".$id_tmp2.".txt -outfmt 5 -strand plus -culling_limit 1 &";
					}
					elsif ($find_haps_in_seed ne "")
					{
						$command = "blastn -query ".$ref_file." -subject ".$output_file1." -out ".$TMP_directory."blast_tmp_".$id."_".$y."_".$double_match_tmp."_".$id_tmp2.".txt -outfmt 5 -strand plus -reward 2 -penalty -3 -evalue 1e-50 -culling_limit 1 &";
					}

					#my $command = sprintf("blastn -query %s -subject %s -out blast_tmp_".$id_tmp.".txt -outfmt 10 -gapextend 1 -gapopen 2", $output_file1, $ref_file);
					#my $command = "blastn -query ".$output_file1." -subject ".$ref_file." -out blast_tmp_".$id_tmp.$double_match_tmp.".txt -outfmt 10 -gapextend 1 -gapopen 2 -best_hit_score_edge 0.03 -best_hit_overhang 0.4 -evalue 1e-50";
					#my $command2 = 'cp sequence_tmp.fasta sequence_tmp_THIS_WORKS.fasta';
					#my $system_result = qx/$command2/;               
					syscmd($command);
					
					if ($save_reads eq "yes")
					{
						if (exists($assembled_reads{$id_tmp2}))
						{}
						else
						{
							$assembled_reads{$id_tmp2} = undef;
							if ($save_reads eq "yes")
							{
								print OUTPUT11 ">".$id_tmp2."_".$id."\n";
								print OUTPUT11 $hash_NP_reads_tmp{$id_tmp2}."\n";
							}  
						}
					}
					$input_BLAST_done{$length_tmp} = undef;
                    
                    if (exists($double_matches{$id_tmp2}) && $double_match_tmp eq "")
                    {
                        $double_match_tmp = "yes";
                        goto DOUBLE_MATCH1;
                    }
                }
				
				my $time_id_by_length = time;
                my $time_sort_matches = $time_id_by_length - $time_kmers;
                print {$filehandle{$seed_id2}} $time_sort_matches." TIME_SORT_MATCHES\n";
				
				my $double_matches_running = '0';

                my $count_BLAST_runnng = keys %input_BLAST_done;

#Read BLAST--------------------------------------------------------------------------------------------------------------------------

                foreach my $blast_db_results_tmp (keys %output_files_DB2)
                {
                    unlink $blast_db_results_tmp;
					delete $output_files_DB2{$blast_db_results_tmp};
                }
BLAST_RUN_NP:	
BLAST1_NP:      foreach my $length_tmp (sort {$a <=> $b} keys %input_files_blast)
                {                
                    my $id_tmp3 = $id_by_length{$length_tmp};
                    
                    foreach my $read_file_tmp (keys %{$input_files_blast{$length_tmp}})
                    {
						foreach my $ref_file_tmp (keys %{$input_files_blast{$length_tmp}{$read_file_tmp}})
						{
							foreach my $length_final_tmp (keys %{$input_files_blast{$length_tmp}{$read_file_tmp}{$ref_file_tmp}})
							{	
								my $double_match_tmp = $input_files_blast{$length_tmp}{$read_file_tmp}{$ref_file_tmp}{$length_final_tmp};              
								my $input_BLAST  =  $TMP_directory."blast_tmp_".$id."_".$y."_".$double_match_tmp."_".$id_tmp3.".txt";
								if ($second_try ne "")
								{
									$input_BLAST  = $TMP_directory."blast_tmp_".$id."_".$y."_".$double_match_tmp."_".$id_tmp3."_2.txt";
								}

								if (-s $input_BLAST)
								{
									my $prev_line = "";
									my $prev_line2 = "";
									open(INPUT_BLAST, $input_BLAST);
									while (my $line_tmp = <INPUT_BLAST>)
									{              
										$prev_line2 = $prev_line;
										$prev_line = $line_tmp;
									}
									close INPUT_BLAST;
									chomp($prev_line2);
									if ($prev_line2 eq "</BlastOutput>")
									{
										open(INPUT_BLAST, $input_BLAST);
									}
									else
									{
										next BLAST1_NP;
									}
								}
								else
								{
									next BLAST1_NP;
								}

								my $length_long_read_tmp = length($hash_NP_reads_tmp{$id_tmp3});
								#my $ref_part = substr $read, -$ref_part{$id_tmp3};
								#my $N_count_ref = $ref_part =~ tr/N/N/;
								if ($ref_part{$id_tmp3} eq '0' || $ref_part{$id_tmp3} eq "")
								{
									print {$filehandle{$seed_id2}} $length_tmp." 0ERROR\n";
									print {$filehandle{$seed_id2}} $id_tmp3." 0ERRORbbb\n";
									goto END1;
								}
								#my $N_correction = ($N_count_ref/$ref_part{$id_tmp3})*100;
								if (exists($double_matches{$id_tmp3}))
								{
									if (exists($double_match_check{$id_tmp3}))
									{
										$double_match_check{$id_tmp3} = "yes2";
									}
									else
									{
										$double_match_check{$id_tmp3} = "yes";
									}
								}							
								
#BLAST1&2-----------------------------------------------------------------------------------------------------
			
								my $j = '0';
								my $ref_pos = '0';
								my $score_match = '0';
								my $score_no_match = '0';
								my $score_match_DUP = '0';
								my $score_no_match_DUP = '0';
								my $bit_score_start = "";
								my $assembly_start_pos_start = "";
								my $assembly_end_pos_start = "";
								my $read_start_pos_start = "";
								my $read_end_pos_start = "";
								my $query_start = "";
								my $subject_start = "";
								my $mismatch_start = "";
								my $alignment_length2_start = "";
								my $pos_matches_start = "";
								my $query = "";
								my $subject = "";
								my $mismatch = "";
								my %query;
								undef %query;
								my %subject;
								undef %subject;
								my %mismatch;
								undef %mismatch;               
								my $bit_score = "";
								my $multiple_matches = "";
								my $multi_gap = '220';
			   
								my $alignment_length = "";
								my $alignment_length2 = '0';
								my $pos_matches = '0';
								my $long_read_start_pos = "";
								my $long_read_end_pos = "";
								my $assembly_end_pos = "";
								my $assembly_start_pos = "";
								my %long_read_start_pos;
								my %long_read_end_pos;
								my %assembly_end_pos;
								my %assembly_start_pos;
								undef %long_read_start_pos;
								undef %long_read_end_pos;
								undef %assembly_end_pos;
								undef %assembly_start_pos;
								my %multi_match_pos_tmp;
								undef %multi_match_pos_tmp;
								
								my $read_pass = "";
								my $accuracy = "";
				  
								my $largest_prev_split = "";
								my $largest_prev_split_ext = "";
								my $skip_id = "";
								my $count_multimatch = '1';
								
								my $assembly_start_pos_low = "1";
								my $long_read_start_pos_high = "1";
								my $assembly_end_pos_high = "1";
								my $long_read_end_pos_high = "1";
								my $match_confirmed = "";
			
INPUT_BLAST_NP:                 while (my $line2 = <INPUT_BLAST>)
								{                                                     
									chomp($line2);
	
									if ($j > 26)
									{
										my @blast_result_tmp;
										undef @blast_result_tmp;
										if ($j < $assembly_start_pos_start || $assembly_start_pos_start eq "")
										{
											@blast_result_tmp = split /\s+/, $line2;
										}
										
										if ($j eq '27' && $blast_result_tmp[0] ne "<Hit>")
										{
											#print {$filehandle{$seed_id2}} "\n".$id_tmp3." ID\n";
											#print {$filehandle{$seed_id2}} $blast_result_tmp[0] ." NO_HIT_REJECTION\n";
											$no_hit_rejection++;
											$skip_id = "no_hit_rejection";
											goto SKIP_BLAST1_NP;
										}
										elsif ($j eq '27')
										{
											#print {$filehandle{$seed_id2}} "\n".$id_tmp3." ID\n";
										}
										 
										if ($blast_result_tmp[0] eq "</BlastOutput>")
										{
											print {$filehandle{$seed_id2}} $assembly_start_pos." REF_START_POS  ".$long_read_start_pos." READ_START_POS  ".$long_read_end_pos." READ_END_POS\n";
											#if ($ref_pos < length($check_prev_splits2{$check_prev_splits_id}))
											#{
												#$read_pos += (length($check_prev_splits2{$check_prev_splits_id}) - $ref_pos);
											   # print {$filehandle{$seed_id2}} length($check_prev_splits2{$check_prev_splits_id}) - $ref_pos." READ_POS_CORRECTION\n";
										   # }
											last INPUT_BLAST_NP;
										}
										elsif ($bit_score_start ne "" && $bit_score_start eq $j)
										{
										}
										elsif ($assembly_start_pos_start ne "" && $assembly_start_pos_start eq $j)
										{
											my @blast_result_tmp2 = split /</, $line2;
											$assembly_start_pos = substr $blast_result_tmp2[1], 15;
											$assembly_start_pos{$count_multimatch} = $assembly_start_pos;
										}
										elsif ($assembly_end_pos_start ne "" && $assembly_end_pos_start eq $j)
										{
											my @blast_result_tmp2 = split /</, $line2;
											$assembly_end_pos = substr $blast_result_tmp2[1], 13;
											$assembly_end_pos{$count_multimatch} = $assembly_end_pos;
										}
										elsif ($read_start_pos_start ne "" && $read_start_pos_start eq $j)
										{
											my @blast_result_tmp2 = split /</, $line2;
											$long_read_start_pos = substr $blast_result_tmp2[1], 13;
											$long_read_start_pos{$count_multimatch} = $long_read_start_pos;     
										}
										elsif ($read_end_pos_start ne "" && $read_end_pos_start eq $j)
										{
											my @blast_result_tmp2 = split /</, $line2;
											$long_read_end_pos = substr $blast_result_tmp2[1], 11;
											$long_read_end_pos{$count_multimatch} = $long_read_end_pos;
										}
										elsif ($pos_matches_start ne "" && $pos_matches_start eq $j)
										{
											my @blast_result_tmp2 = split /</, $line2;
											my $pos_matches_tmp = substr $blast_result_tmp2[1], 13;
											$pos_matches += $pos_matches_tmp;
										}
										elsif ($alignment_length2_start ne "" && $alignment_length2_start eq $j)
										{   
											my @blast_result_tmp2 = split /</, $line2;
											my $alignment_length2_tmp = substr $blast_result_tmp2[1], 14;
											$alignment_length2 += $alignment_length2_tmp;
											
											#if ($multiple_matches eq "" && $alignment_length2_tmp < 0.3*$length_tmp && $alignment_length2_tmp < 0.3*length($ref_part))
											#{
												#$almost_no_hit_rejection++;
												#$skip_id = "yes";
												#close INPUT_BLAST;
												#goto SKIP_BLAST1_NP;
											#}
											
											if ($accuracy eq "")
											{
												$accuracy = $pos_matches/$alignment_length2;
											}
											else
											{
												my $accuracy_tmp = $pos_matches/$alignment_length2;
												my $accuracy_tmp2 = $accuracy;
												$accuracy = ($accuracy_tmp+$accuracy_tmp2)/2   
											}
											
											if ($accuracy < 0.7)
											{
												#print {$filehandle{$seed_id2}} $accuracy." ACCURACY_REJECTION\n";
												$accuracy_rejection++;
												$skip_id = "accuracy_rejection";
												goto SKIP_BLAST1_NP;
											}    
										}
										elsif ($query_start ne "" && $query_start eq $j)
										{
											my @blast_result_tmp2 = split /</, $line2;
											$query = substr $blast_result_tmp2[1], 9;
											$query{$count_multimatch} = $query;
											#print {$filehandle{$seed_id2}} $query." QUERY\n";    
										}
										elsif ($subject_start ne "" && $subject_start eq $j)
										{
											my @blast_result_tmp2 = split /</, $line2;
											$subject = substr $blast_result_tmp2[1], 9;
											$subject{$count_multimatch} = $subject;
											#print {$filehandle{$seed_id2}} $subject." SUBJECT\n";
										}
										elsif ($mismatch_start ne "" && $mismatch_start eq $j && $find_haps_in_seed eq "")
										{
											my @blast_result_tmp2 = split /</, $line2;
											$mismatch = substr $blast_result_tmp2[1], 12;
											$mismatch{$count_multimatch} = $mismatch;
											#print {$filehandle{$seed_id2}} $mismatch." MISMATCH\n";
										}
										elsif ($mismatch_start ne "" && $j eq $mismatch_start+2)
										{
											@blast_result_tmp = split /\s+/, $line2;
											my $multi_match_length_check = $assembly_end_pos-$assembly_start_pos;
											
											if (($blast_result_tmp[1] ne "<Hsp>" || (($long_read_start_pos <= 150 || $assembly_start_pos < 150) && $assembly_end_pos >= $ref_part{$id_tmp3}-200
												&& $long_read_end_pos <= $length_long_read_tmp-200)) && ($multiple_matches eq "" || ($multi_match_length_check > 1000 && $assembly_end_pos > $ref_part{$id_tmp3}-$ref_end_rejection_threshold && $long_read_end_pos > 400
												&& (($long_read_start_pos <= 150 || $assembly_start_pos < 150) && $assembly_end_pos >= $ref_part{$id_tmp3}-200 && $long_read_end_pos <= $length_long_read_tmp-200))))
											{                                                                            
												if ($assembly_end_pos < $ref_part{$id_tmp3}-$ref_end_rejection_threshold)
												{
													#print {$filehandle{$seed_id2}} $assembly_end_pos." REF_END_REJECTION\n";                                     
													$ref_end_rejection++;
													$skip_id = "ref_end_rejection";
													goto SKIP_BLAST1_NP;
												}                                
												elsif ($long_read_end_pos < 400)
												{
													$match_too_short++;
													$skip_id = "match_too_short";
													goto SKIP_BLAST1_NP;
												}
												elsif ($long_read_start_pos > 220 && $assembly_start_pos > 150)
												{
													#print {$filehandle{$seed_id2}} $assembly_start_pos." REF ".$long_read_start_pos." READ_START_POS_REJECTION\n";
													$read_start_pos_rejection++; 
													my $alignment_length_tmp = $assembly_end_pos - $assembly_start_pos + 1; 
													$read_start_pos_rej{$id_tmp3} = $alignment_length_tmp;
													$long_read_end_pos_save{$id_tmp3} = $long_read_end_pos;
													$skip_id = "read_start_pos_rejection";
													goto SKIP_BLAST1_NP;
												}
												elsif ($long_read_end_pos > $length_long_read_tmp-200)
												{
													#print {$filehandle{$seed_id2}} $long_read_end_pos." EXT_TOO_SHORT\n";
													$ext_too_short++;
													$skip_id = "ext_too_short";
													goto SKIP_BLAST1_NP;
												}
												
												$alignment_length = $assembly_end_pos - $assembly_start_pos + 1;                                  
												$alignment_length_save{$id_tmp3} = $alignment_length;
												$match_confirmed = "yes";
											}
#Verify multiple matches--------------------------------------------------------------------------                          
											elsif ($multiple_matches ne "")
											{                                                      
												my $last_multi_match = "";
												if ($multi_match_length_check < 400 && $multi_match_length_check > -400)
												{
													delete $assembly_start_pos{$count_multimatch};
													delete $assembly_end_pos{$count_multimatch};
													delete $long_read_start_pos{$count_multimatch};
													delete $long_read_end_pos{$count_multimatch};
													delete $mismatch{$count_multimatch};
													$last_multi_match = "yes";
												}
												
												my $end_pos_check = "";
												my $start_pos_check = "";
												my $rejection_tmp = "";
												my %end_pos_check;
												my %start_pos_check;
												undef %end_pos_check;
												undef %start_pos_check;
												
												foreach my $count_tmp (keys %assembly_end_pos)
												{
													if ($assembly_end_pos{$count_tmp} > $ref_part{$id_tmp3}-$ref_end_rejection_threshold && $long_read_end_pos{$count_tmp} > 400)
													{
														$end_pos_check = "yes";
														$end_pos_check{$count_tmp} = undef;
													}
												}
												if ($end_pos_check eq "yes")
												{
													foreach my $count_tmp (keys %assembly_start_pos)
													{
														if ($assembly_start_pos{$count_tmp} < 220 || $long_read_start_pos{$count_tmp} < 120)
														{
															$start_pos_check = "yes";
															$start_pos_check{$count_tmp} = undef;
														}
													}
													if ($start_pos_check ne "yes")
													{
														$rejection_tmp = "read_start";
													}
												}
												else
												{
													$rejection_tmp = "ref_end";
												}
												
												my $match_found_tmp = "";
												
												if ($rejection_tmp eq "")
												{
													my @alignment_order_assembly;                                    
													
ALIGNMENT_ORDER_NP0:                                foreach my $count_tmp (keys %start_pos_check)
													{
														undef @alignment_order_assembly;
														push @alignment_order_assembly, $count_tmp;
ALIGNMENT_ORDER_NP:                                           
														foreach my $count_tmp2 (keys %assembly_start_pos)
														{
															if ($assembly_start_pos{$count_tmp2} > $assembly_end_pos{$count_tmp}-$multi_gap && $assembly_start_pos{$count_tmp2} < $assembly_end_pos{$count_tmp}+$multi_gap
																&& $long_read_start_pos{$count_tmp2} > $long_read_end_pos{$count_tmp}-$multi_gap && $long_read_start_pos{$count_tmp2} < $long_read_end_pos{$count_tmp}+$multi_gap)
															{
																my $assembly_end_pos_tmp2 = $ref_part{$id_tmp3} - $assembly_end_pos{$count_tmp};
																while (exists($multi_match_pos{$assembly_end_pos_tmp2}))
																{
																	$assembly_end_pos_tmp2++;
																}
																$multi_match_pos_tmp{$assembly_end_pos_tmp2} = $id_tmp3;
																
																push @alignment_order_assembly, $count_tmp2;
																if (exists($end_pos_check{$count_tmp2}))
																{
																	$match_found_tmp = "yes";
																	$long_read_end_pos_high = $count_tmp2;
																	$assembly_end_pos_high = $count_tmp2;
																	$assembly_start_pos_low = $alignment_order_assembly[0];																
																	last ALIGNMENT_ORDER_NP0;
																}
																else
																{
																	$count_tmp = $count_tmp2;
																	goto ALIGNMENT_ORDER_NP;
																}
															}
														}
													}
												}
			  
												if (($blast_result_tmp[1] ne "<Hsp>" || $count_multimatch > 7 || $last_multi_match ne "") && $match_found_tmp eq "")
												{
													if ($end_pos_check eq "yes")
													{
														my $end_pos_rej_tmp = "";
														my $start_pos_rej_tmp = "";
														my $read_end_pos_rej_tmp = "";
														my $count_contigs = '0';
														foreach my $count_tmp (keys %end_pos_check)
														{
															my $lowest_count_tmp = $count_tmp;
															my $count_contigs_tmp = '1';
CHECK_MULTI_MATCH_REJ:																		
															foreach my $count_tmp2 (keys %assembly_start_pos)
															{
																if ($count_tmp2 ne $lowest_count_tmp && $assembly_end_pos{$count_tmp2} > $assembly_start_pos{$lowest_count_tmp}-$multi_gap && $assembly_end_pos{$count_tmp2} < $assembly_start_pos{$lowest_count_tmp}+$multi_gap
																	&& $long_read_end_pos{$count_tmp2} > $long_read_start_pos{$lowest_count_tmp}-$multi_gap && $long_read_end_pos{$count_tmp2} < $long_read_start_pos{$lowest_count_tmp}+$multi_gap
																	&& $assembly_start_pos{$count_tmp2} < $assembly_start_pos{$lowest_count_tmp})
																{
																	$lowest_count_tmp = $count_tmp2;
																	$count_contigs_tmp++;
																	goto CHECK_MULTI_MATCH_REJ;
																}
															}
															if ($end_pos_rej_tmp eq "" || (($assembly_end_pos{$count_tmp}-$assembly_start_pos{$lowest_count_tmp}) > ($end_pos_rej_tmp-$start_pos_rej_tmp)))
															{
																$end_pos_rej_tmp = $assembly_end_pos{$count_tmp};
																$start_pos_rej_tmp = $assembly_start_pos{$lowest_count_tmp};
																$read_end_pos_rej_tmp = $long_read_end_pos{$count_tmp};
																$count_contigs = $count_contigs_tmp;
															}
														}
														
														if ($end_pos_rej_tmp ne "" && $start_pos_rej_tmp ne "")
														{
															$assembly_end_pos = $end_pos_rej_tmp;
															$assembly_start_pos = $start_pos_rej_tmp;
															my $alignment_length_tmp = $assembly_end_pos - $assembly_start_pos + 1;
															$read_start_pos_rej{$id_tmp3} = $alignment_length_tmp;
															$long_read_end_pos_save{$id_tmp3} = $read_end_pos_rej_tmp;
	
															if ($count_contigs eq '1')
															{
																$read_start_pos_rejection++;
																$skip_id = "read_start_pos_rejection";
															}
															else
															{
																$skip_id = "multi_match_rejection";
																$multi_match_rejection++;
															}
															goto SKIP_BLAST1_NP;
														}
													}
													else
													{
														$skip_id = "ref_end_rejection";
														$ref_end_rejection++;
													}
													goto SKIP_BLAST1_NP;
												}
												elsif ($blast_result_tmp[1] eq "<Hsp>" && $count_multimatch < 8 && $match_found_tmp eq "")
												{
													$j++;
													$bit_score_start = $j+1;
													$assembly_start_pos_start = $j+4;
													$assembly_end_pos_start = $j+5;
													$read_start_pos_start = $j+6;
													$read_end_pos_start = $j+7;
													$pos_matches_start = $j+11;
													$alignment_length2_start = $j+13;
													$query_start = $j+14;
													$subject_start = $j+15;
													$mismatch_start = $j+16;
													$multiple_matches = "yes";
													$count_multimatch++;
													next INPUT_BLAST_NP;
												}
												elsif ($match_found_tmp ne "")
												{
													if ($long_read_end_pos{$assembly_end_pos_high} > $length_long_read_tmp-250)
													{
														$ext_too_short++;
														$skip_id = "ext_too_short";
														goto SKIP_BLAST1_NP;
													}
													elsif ($long_read_end_pos{$assembly_end_pos_high} < 500)
													{
														$match_too_short++;
														$skip_id = "match_too_short";
														goto SKIP_BLAST1_NP;
													}     
													
													foreach my $count_tmp3 (keys %mismatch)
													{
														if ($assembly_start_pos{$count_tmp3} < $assembly_start_pos{$assembly_start_pos_low}-100)
														{
															delete $mismatch{$count_tmp3};
														}
													}
																	
													$alignment_length = $assembly_end_pos{$assembly_end_pos_high} - $assembly_start_pos{$assembly_start_pos_low};
													$long_read_end_pos = $long_read_end_pos{$assembly_end_pos_high};
													print {$filehandle{$seed_id2}} "\n".$id_tmp3." MULTI_MATCH\n";
													$alignment_length_save{$id_tmp3} = $alignment_length;
													foreach my $multi_match_pos_tmp (keys %multi_match_pos_tmp)
													{
														$multi_match_pos{$multi_match_pos_tmp} = $multi_match_pos_tmp{$multi_match_pos_tmp};
														$multi_match{$id_tmp3}{$multi_match_pos_tmp} = undef;
													}
													$match_confirmed = "yes";
												}
												else
												{
													print {$filehandle{$seed_id2}} "\n".$id_tmp3." MULTI_MATCH_CHECK\n";
												}
											}
											elsif ($allow_multi_match eq "yes")
											{
												#print {$filehandle{$seed_id2}} "\n".$id_tmp3." MULTIPLE_MATCHES\n";
												$j++;
												$bit_score_start = $j+1;
												$assembly_start_pos_start = $j+4;
												$assembly_end_pos_start = $j+5;
												$read_start_pos_start = $j+6;
												$read_end_pos_start = $j+7;
												$pos_matches_start = $j+11;
												$alignment_length2_start = $j+13;
												$query_start = $j+14;
												$subject_start = $j+15;
												$mismatch_start = $j+16;
												$multiple_matches = "yes";
												$count_multimatch++;
												next INPUT_BLAST_NP;
											}
											else
											{
												$multi_match_rejection++;
												$skip_id = "multi_match_rejection";
												print {$filehandle{$seed_id2}} $assembly_end_pos." ".$ref_part{$id_tmp3}." ASS_END_POS_REJ\n";
												print {$filehandle{$seed_id2}} $assembly_end_pos." ".$assembly_start_pos." ASS_END_POS_REJ2\n";
												if ($assembly_end_pos > $ref_part{$id_tmp3}-$ref_end_rejection_threshold)
												{
													my $alignment_length_tmp = $assembly_end_pos - $assembly_start_pos + 1;
													$read_start_pos_rej{$id_tmp3} = $alignment_length_tmp;
													$long_read_end_pos_save{$id_tmp3} = $long_read_end_pos;
												}
												goto SKIP_BLAST1_NP;
											}
											
											if ($alignment_length < 550 && length($read) > 600)
											{
												#print {$filehandle{$seed_id2}} $alignment_length." OVERLAP_TOO_SHORT\n";
												$overlap_too_short++;
												$skip_id = "overlap_too_short";           
												goto SKIP_BLAST1_NP;
											}
											$count_matches_with_high_scores++;
											$long_read_end_pos_save{$id_tmp3} = $long_read_end_pos;
											
											my $pos_tmp = $assembly_start_pos-50;
											if ($pos_tmp < 0)
											{
												$pos_tmp = '0';
											}
											#my $ref_part_tmp = substr $ref_part, $pos_tmp;
											
											my $long_read_tmp = "";
			
											if (exists($reverse_list{$id_tmp3}) && $double_match_tmp eq "")
											{
												$long_read_tmp = reverse($hash_NP_reads_tmp{$id_tmp3});
												$long_read_tmp =~ tr/ACTG/TGAC/;
											}
											else
											{
												$long_read_tmp = $hash_NP_reads_tmp{$id_tmp3};
											}
											$check_prev_splits{$id_tmp3} = $long_read_tmp;
											
											$accuracy{$id_tmp3} = $accuracy;
								  
											my $track_length = '-1';
											
											if ($find_haps_in_seed eq "")
											{
												foreach my $id_split (keys %split_positions)
												{
													if ($id_split eq $id)
													{
SPLIT_NP:                                   		    foreach my $prev_splits (sort {$a <=> $b} keys %{$split_positions{$id_split}})
														{
															if ($last_hap_pos eq "")
															{
																$last_hap_pos = $position - $prev_splits;
															}
															$track_length = '-1';
			
															if ($alignment_length > $position - $prev_splits)
															{
																foreach my $count_mismatch_tmp (keys %mismatch)
																{ 
																	if ($assembly_end_pos{$assembly_end_pos_high}-$assembly_end_pos{$count_mismatch_tmp} < ($position - $prev_splits) &&
																		($assembly_end_pos{$assembly_end_pos_high}-$assembly_start_pos{$count_mismatch_tmp}) > ($position - $prev_splits))
																	{
																		$largest_prev_split = $ref_part{$id_tmp3}-$assembly_start_pos{$count_mismatch_tmp} - ($position - $prev_splits)+1; 
																		$largest_prev_split_ext = $split_positions{$id_split}{$prev_splits};
			
																		my $b = '0';
																		my $last12 = "";
																		my $last12_match = "";
																
																		while ($b < length($query{$count_mismatch_tmp}))
																		{
																			my $nuc_query = substr $query{$count_mismatch_tmp}, $b, 1;
																			if ($nuc_query ne "-")
																			{
																				$track_length += 1;
																				if ($track_length > $largest_prev_split-15 && $track_length < $largest_prev_split+5)
																				{
																					$last12 .= $nuc_query;
																					if (length($last12) > 12)
																					{
																						substr $last12, 0, 1, "";
																					}
																					if (length($last12) > 11)
																					{
																						my @largest_prev_split_ext = split /,/, $largest_prev_split_ext;
																						
																						if ($last12 eq $largest_prev_split_ext[0])
																						{
																							$last12_match = "yes";
																							$b++;
																							next;
																						}
																						if ($last12_match eq "yes")
																						{
																							$last12_match = "";
																							my $nuc_subject = substr $subject{$count_mismatch_tmp}, $b, 1;
																							my $nuc_mismatch_prev = substr $mismatch{$count_mismatch_tmp}, $b-1, 1;
																							my $nuc_mismatch_post = substr $mismatch{$count_mismatch_tmp}, $b+1, 1;
			
																							if ($nuc_query ne $largest_prev_split_ext[1])
																							{
																								print {$filehandle{$seed_id2}} $prev_splits." SPLIT_POS\n";
																				print {$filehandle{$seed_id2}} $largest_prev_split." SPLIT_POS2\n";
																				print {$filehandle{$seed_id2}} $largest_prev_split_ext." SPLIT_POS_ECT\n";
																				print {$filehandle{$seed_id2}} $assembly_start_pos{$count_mismatch_tmp}." START_POS\n";
																				print {$filehandle{$seed_id2}} $ref_part{$id_tmp3}." REF_PART\n";
																				print {$filehandle{$seed_id2}} $query{$count_mismatch_tmp}." QUERY\n";
																								print {$filehandle{$seed_id2}} $nuc_query." NUC_QUERY\n";
																								print {$filehandle{$seed_id2}} $nuc_subject." NUC_SUBJECT\n\n";
																								my $test_seq = substr $query{$count_mismatch_tmp}, $b-15, 30;
																								print {$filehandle{$seed_id2}} $test_seq." TEST_SEQ\n";
																							}
																							if ($nuc_mismatch_prev ne "|" && $nuc_mismatch_post ne "|")
																							{
																							}
																							elsif ($nuc_query eq $nuc_subject)
																							{
																								$score_match++;
																							}
																							elsif ($nuc_subject ne "-")
																							{
																								$score_no_match++;
																							}
																							next SPLIT_NP;
																						}
																					}
																				}
																			}
																			$b++;
																		}
																	}
																}
															}
														}
													}
												}
												
#Store mismatches of the aligned reads to split reads---------------------------------------------------------------------                              
												foreach my $count_mismatch_tmp (keys %mismatch)
												{
													if ($multiple_matches eq "yessss")
													{
									print {$filehandle{$seed_id2}} $id_tmp3." ".$ref_part{$id_tmp3}." ".$position." ".$assembly_start_pos{$count_mismatch_tmp}." ".$assembly_end_pos{$count_mismatch_tmp}." COUNT_MISMATCH\n";
									print {$filehandle{$seed_id2}} $query{$count_mismatch_tmp}." QUERY\n";
									print {$filehandle{$seed_id2}} $subject{$count_mismatch_tmp}." SUBJECT\n";
									print {$filehandle{$seed_id2}} $mismatch{$count_mismatch_tmp}." MISMATCH\n";
													}
													$track_length = '-1';
													my $track_query = $ref_part{$id_tmp3}-$assembly_start_pos{$count_mismatch_tmp}+1;
													my $track_query_full = $assembly_start_pos{$count_mismatch_tmp}+$position-$ref_part{$id_tmp3};
													my @blast_result_tmp3 = split /\s+/, $mismatch{$count_mismatch_tmp};
													
													foreach my $blast_tmp (@blast_result_tmp3)
													{
														$track_length += length($blast_tmp);
														$track_query -= length($blast_tmp);
														$track_query_full += length($blast_tmp);
						if ($multiple_matches eq "yesss")
													{
														print {$filehandle{$seed_id2}} $blast_tmp." BLAST_TMP ".$track_length." TRACK_LENGTH\n";
													}
														my $track_length_adjust = '0';
														my $nuc_query = substr $query{$count_mismatch_tmp}, $track_length+1, 1;
														my $nuc_query_prev = substr $query{$count_mismatch_tmp}, $track_length, 1;
														my $nuc_query_post = substr $query{$count_mismatch_tmp}, $track_length+2, 1;
														my $nuc_subject = substr $subject{$count_mismatch_tmp}, $track_length+1, 1;
														my $nuc_subject_prev = substr $subject{$count_mismatch_tmp}, $track_length, 1;
														my $nuc_subject_post = substr $subject{$count_mismatch_tmp}, $track_length+2, 1;
														while ($nuc_query ne $nuc_subject)
														{
															$track_length_adjust++;                                       
															if ($track_length < length($query{$count_mismatch_tmp})-1)
															{
																my $one = $track_query_full-1;
																my $two = $track_query_full+1;
																if (exists($quality_scores{$seed_id}{$track_query_full}))
																{
																	my @q_score_tmp = split / /, $quality_scores{$seed_id}{$track_query_full}; 
																	if ($q_score_tmp[0] < 0.8)
																	{
																		goto SKIP_STORE_MISMATCH_NP;
																	}												
																	my @q_score_gap_tmp = split / /, $quality_scores_gap{$seed_id}{$track_query_full}; 
																	if (exists($quality_scores_gap{$seed_id}{$track_query_full}) && $q_score_gap_tmp[0] < 0.7)
																	{
																		goto SKIP_STORE_MISMATCH_NP;
																	}
																	my @q_score_gap1_tmp = split / /, $quality_scores_gap{$seed_id}{$one}; 
																	if (exists($quality_scores_gap{$seed_id}{$one}) && $q_score_gap1_tmp[0] < 0.7)
																	{
																		goto SKIP_STORE_MISMATCH_NP;
																	}
																	my @q_score_gap2_tmp = split / /, $quality_scores_gap{$seed_id}{$two}; 
																	if (exists($quality_scores_gap{$seed_id}{$two}) && $q_score_gap2_tmp[0] < 0.7)
																	{
																		goto SKIP_STORE_MISMATCH_NP;
																	}
																}
																elsif ($PB_reads eq "" && $input_reads_DB_folder_PB eq "")
																{
																	if ($track_query_full > $original_seed_length{$id} && $assembly_length_max ne "WG")
																	{
																	print {$filehandle{$seed_id2}} $track_query_full." QS!!!\n";
																	}
																	goto SKIP_STORE_MISMATCH_NP;
																}
																
																if ($nuc_query ne "-" && $nuc_subject ne "-" && $nuc_query ne "N" && $nuc_query_prev ne "N" && $nuc_query_post ne "N"
																	&& $nuc_query_prev ne "-" && $nuc_query_post ne "-" && $nuc_subject_prev ne "-" && $nuc_subject_post ne "-")
																{
																	my $track_length_tmp = $track_length+$track_length_adjust;
																	my $seq_tmp = substr $read, -$track_query-30, 31;
																	
										if ($multiple_matches eq "yessssss")
													{
																	#print {$filehandle{$seed_id2}} $id_tmp3." ID\n";
																	print {$filehandle{$seed_id2}} $quality_scores{$seed_id}{$track_query_full}." QS\n";
																		print {$filehandle{$seed_id2}} $nuc_query." NQ\n";
														print {$filehandle{$seed_id2}} $nuc_subject." NS\n";
																		print {$filehandle{$seed_id2}} $track_query_full." TRACK_QUERY\n";
																		my $nuc_query_seq = substr $query{$count_mismatch_tmp}, $track_length_tmp-30, 31;
																		my $nuc_read_seq = substr $read, $track_query_full-31, 31;
																		print {$filehandle{$seed_id2}} $nuc_query_seq." SEQ1\n";
																	   print {$filehandle{$seed_id2}} $seq_tmp." SEQ2\n";
																	   print {$filehandle{$seed_id2}} $nuc_read_seq." READ_SEQ\n";
													}
													
																	$store_mismatches_NP{$id_tmp3}{$track_query_full} = $nuc_query.",".$nuc_subject.",".$seq_tmp;                 
																}																
SKIP_STORE_MISMATCH_NP:                                                    
																my $track_length_tmp = $track_length+$track_length_adjust;
																my $seq_tmp = substr $read, -$track_query-30, 31;
																if ($nuc_query ne "N")
																{
																	$store_mismatches_all_NP{$id_tmp3}{$track_query_full} = $nuc_query.",".$nuc_subject.",".$seq_tmp;
																}
																if ($nuc_query eq "N")
																{
																	$store_mismatches_N_NP{$id_tmp3}{$track_query_full} = $nuc_query.",".$nuc_subject.",".$seq_tmp;
																}
																if ($multiple_matches eq "yessssss")
													{
																print {$filehandle{$seed_id2}} $track_query_full." TRACK_QUERY_FULL\n";
													}
																#print {$filehandle{$seed_id2}} $track_length_adjust." TRACK_LENGTH_ADJUST\n";
																#print {$filehandle{$seed_id2}} $nuc_query." NUC_QUERY\n";
																#print {$filehandle{$seed_id2}} $nuc_subject." NUC_SUBJECT\n";
																
																if ($nuc_query ne "-")
																{
																	$track_query--;
																	$track_query_full++;
																}													
															}
															
															$nuc_query = substr $query{$count_mismatch_tmp}, $track_length+$track_length_adjust+1, 1;
															$nuc_subject = substr $subject{$count_mismatch_tmp}, $track_length+$track_length_adjust+1, 1;
														}                                                    
														$track_length += $track_length_adjust;  
													}
												}
												
#Check split positions derived from duplicated regions--------------------------------------------------------------------------------------------------------------------------
												foreach my $id_split (keys %split_positions_DUP)
												{
													if ($id_split eq $id)
													{
SPLIT_NP_DUP:                               	        foreach my $prev_splits (sort {$a <=> $b} keys %{$split_positions_DUP{$id_split}})
														{
															if ($last_hap_pos_DUP eq "")
															{
																$last_hap_pos_DUP = $position - $prev_splits;
															}
															$track_length = '-1';
			
															if ($alignment_length > $position - $prev_splits)
															{
																foreach my $count_mismatch_tmp (keys %mismatch)
																{ 
																	if ($assembly_end_pos{$assembly_end_pos_high}-$assembly_end_pos{$count_mismatch_tmp} < ($position - $prev_splits) &&
																		($assembly_end_pos{$assembly_end_pos_high}-$assembly_start_pos{$count_mismatch_tmp}) > ($position - $prev_splits))
																	{
																		$largest_prev_split = $ref_part{$id_tmp3}-$assembly_start_pos{$count_mismatch_tmp} - ($position - $prev_splits)+1; 
																		$largest_prev_split_ext = $split_positions_DUP{$id_split}{$prev_splits};
			
																		my $b = '0';
																		my $last12 = "";
																		my $last12_match = "";
																
																		while ($b < length($query{$count_mismatch_tmp}))
																		{
																			my $nuc_query = substr $query{$count_mismatch_tmp}, $b, 1;
																			if ($nuc_query ne "-")
																			{
																				$track_length += 1;
																				if ($track_length > $largest_prev_split-15 && $track_length < $largest_prev_split+5)
																				{
																					$last12 .= $nuc_query;
																					if (length($last12) > 12)
																					{
																						substr $last12, 0, 1, "";
																					}
																					if (length($last12) > 11)
																					{
																						my @largest_prev_split_ext = split /,/, $largest_prev_split_ext;
																						
																						if ($last12 eq $largest_prev_split_ext[0])
																						{
																							$last12_match = "yes";
																							$b++;
																							next;
																						}
																						if ($last12_match eq "yes")
																						{
																							$last12_match = "";
																							my $nuc_subject = substr $subject{$count_mismatch_tmp}, $b, 1;
																							my $nuc_mismatch_prev = substr $mismatch{$count_mismatch_tmp}, $b-1, 1;
																							my $nuc_mismatch_post = substr $mismatch{$count_mismatch_tmp}, $b+1, 1;
			
																							if ($nuc_subject ne $largest_prev_split_ext[1] && $nuc_query ne $nuc_subject && $nuc_subject ne "-")
																							{
																								if ($nuc_query eq "N")
																								{
																									delete $split_positions_DUP{$id_split}{$prev_splits};
																								}
																								print {$filehandle{$seed_id2}} $prev_splits." SPLIT_POS_DUP\n";
																				print {$filehandle{$seed_id2}} $largest_prev_split." SPLIT_POS2_DUP\n";
																				print {$filehandle{$seed_id2}} $largest_prev_split_ext." SPLIT_POS_ECT_DUP\n";
																				print {$filehandle{$seed_id2}} $assembly_start_pos{$count_mismatch_tmp}." START_POS_DUP\n";
																				print {$filehandle{$seed_id2}} $ref_part{$id_tmp3}." REF_PART_DUP\n";
																				#print {$filehandle{$seed_id2}} $query{$count_mismatch_tmp}." QUERY_DUP\n";
																								print {$filehandle{$seed_id2}} $nuc_query." NUC_QUERY_DUP\n";
																								print {$filehandle{$seed_id2}} $nuc_subject." NUC_SUBJECT_DUP\n\n";
																								my $test_seq = substr $query{$count_mismatch_tmp}, $b-15, 30;
																								print {$filehandle{$seed_id2}} $test_seq." TEST_SEQ_DUP\n";
																							}
																							if ($nuc_mismatch_prev ne "|" && $nuc_mismatch_post ne "|")
																							{
																							}
																							elsif ($nuc_query eq $nuc_subject)
																							{
																								$score_match_DUP++;
																							}
																							elsif ($nuc_subject ne "-")
																							{
																								$score_no_match_DUP++;
																							}
																							next SPLIT_NP_DUP;
																						}
																					}
																				}
																			}
																			$b++;
																		}
																	}
																}
															}
														}
													}
												}
											}
											if ($match_confirmed eq "yes")
											{
												goto SKIP_BLAST1_NP;
											}
										}
										elsif ($blast_result_tmp[1] eq "<Hsp>")
										{ 
											$j++;
											$bit_score_start = $j+1;
											$assembly_start_pos_start = $j+4;
											$assembly_end_pos_start = $j+5;
											$read_start_pos_start = $j+6;
											$read_end_pos_start = $j+7;
											$pos_matches_start = $j+11;
											$alignment_length2_start = $j+13;
											$query_start = $j+14;
											$subject_start = $j+15;
											$mismatch_start = $j+16;
											next INPUT_BLAST_NP;
										}
										else
										{
											$j++;
											next INPUT_BLAST_NP;
										}
									}
									$j++
								}
													
SKIP_BLAST1_NP:
								delete $input_files_blast{$length_tmp}{$read_file_tmp};
								if (keys %{$input_files_blast{$length_tmp}} < 1)
								{
									delete $input_files_blast{$length_tmp};
								}
								delete $input_BLAST_done{$length_tmp};
								if (exists($double_matches{$id_tmp3}) && $double_match_tmp eq "")
								{
									$double_matches_running--;
								}
								#print {$filehandle{$seed_id2}} $id_tmp3." FINISHED_IDs\n";
								#print {$filehandle{$seed_id2}} time." TIME\n";
								
								#delete $id_matches{$id_tmp3};
								#delete $id_by_length{$length_tmp};
								close INPUT_BLAST;
#ADD_REJ_READS-----------------------------------------------------------------------------------------------------------------               
								if (($add_rejected_reads ne "" && (exists($read_start_pos_rej_saved{$id_tmp3}) || exists($read_start_pos_rej_saved{$id_tmp3}))) || ($add_no_match_reads ne "" && exists($extensions_nomatch2b_saved{$id_tmp3})))
								{
									$add_rej_reads_extra{$id_tmp3} = undef;
									delete $save_reads_for_next{$id_tmp3};
								}
#-----------------------------------------------------------------------------------------------------------------
								
								if ($skip_id ne "")
								{
									$rejected_reads_save{$id_tmp3} = undef;
									if (exists($double_match_check{$id_tmp3}))
									{
										if ($double_match_check{$id_tmp3} eq "yes")
										{
											
										}
										elsif ($double_match_check{$id_tmp3} eq "yes2")
										{
											unless (exists($double_match_reject{$id_tmp3}))
											{
												delete $read_start_pos_rej{$id_tmp3};
												delete $rejected_reads_save{$id_tmp3};
												delete $add_rej_reads_extra{$id_tmp3};
												if ($double_match_tmp eq "")
												{
													delete $reverse_list{$id_tmp3};
												}
											}										
										}
										$double_match_reject{$id_tmp3} = undef;
									}
									
#Exclude reads------------------------
									if ($skip_id eq "multi_match_rejection" || $skip_id eq "ref_end_rejection")
									{
										my $exclude_pos = $long_read_end_pos;
										if (length($read) < $long_read_end_pos)
										{
											$exclude_pos = length($read);
										}
										if ($hap_tag eq "HAP1" && $repetitive_detect1 eq "" && $repetitive_detect2 eq "ff")
										{
											$exclude_reads_hap1_NP{$id_tmp3} = $position+$exclude_pos;
											#print {$filehandle{$seed_id2}} $id_tmp." EXCLUDE1\n";
										}
										elsif ($hap_tag eq "HAP2" && $repetitive_detect1 eq "" && $repetitive_detect2 eq "ff")
										{
											$exclude_reads_hap2_NP{$id_tmp3} = $position+$exclude_pos;
											#print {$filehandle{$seed_id2}} $id_tmp." EXCLUDE2\n";
										}
									}
#Exclude reads------------------------        
									next BLAST1_NP;  
								}
								else
								{
									if (exists($double_match_reject{$id_tmp3}))
									{
										delete $read_start_pos_rej{$id_tmp3};
										delete $add_rej_reads_extra{$id_tmp3};
										delete $rejected_reads_save{$id_tmp3};
										if ($double_match_tmp eq "yes")
										{
											delete $reverse_list{$id_tmp3};
										}
									}
									elsif ($double_match_check{$id_tmp3} eq "yes2")
									{
										print {$filehandle{$seed_id2}} $id_tmp3." DOUBLE_MATCH_ERROR\n";
									}
								}
							 
#Check for previous variations of the assembly in the long reads------------------------------------------------------------------------------------------                   
#---------------------------------------------------------------------------------------------------------------------------------------------------------
			
								print {$filehandle{$seed_id2}} "\n".$id_tmp3." ID\n";
								#$score_match += $score_match_saved{$id_tmp};
								#$score_no_match += $score_no_match_saved{$id_tmp};
								#$score_match_DUP += $score_match_DUP_saved{$id_tmp};
								#$score_no_match_DUP += $score_no_match_DUP_saved{$id_tmp};
								if ($score_match > 0)
								{
									print {$filehandle{$seed_id2}} $score_match." EXT_MATCH\n";
								}
								if ($score_no_match > 0)
								{
									print {$filehandle{$seed_id2}} $score_no_match." EXT_NO_MATCH\n";
								}
								if ($score_match_DUP ne 0)
								{
									print {$filehandle{$seed_id2}} $score_match_DUP." EXT_MATCH_DUP\n";
								}
								if ($score_no_match_DUP ne 0)
								{
									print {$filehandle{$seed_id2}} $score_no_match_DUP." EXT_NO_MATCH_DUP\n";
								}
								my $length_read = length($check_prev_splits{$id_tmp3});
								print {$filehandle{$seed_id2}} $alignment_length." LENGTH_MATCH\n";
								#print {$filehandle{$seed_id2}} $accuracy." ACCURACY\n";
								$score_matches_save{$id_tmp3} = $score_match."_".$score_no_match."_".$score_match_DUP."_".$score_no_match_DUP;   
			
								if ($score_no_match > $score_match || ($score_no_match_DUP > $score_match_DUP && ($score_no_match_DUP > 3 || $score_match eq '0' || $score_no_match+$score_no_match_DUP > $score_match+$score_match_DUP)))
								{
									my $ext = substr $check_prev_splits{$id_tmp3}, $long_read_end_pos-90;
									$length_ext_all{$id_tmp3} = length($ext);
									
									if (length($ext) > 200)
									{
										$extensions_nomatch{$ext} = $id_tmp3;
										$extensions_nomatch2{$id_tmp3} = $ext;
										if ($score_no_match > $score_match)
										{
											$extensions_nomatch2b{$id_tmp3}{$ext}{$score_no_match} = $score_match;
											$extensions_nomatch2b_count++;
										}
									}
									if ($score_no_match > $score_match && $score_no_match > 2 && $score_no_match > $score_match*2)
									{
#Exclude reads------------------------                        
										my $exclude_pos = $long_read_end_pos;
										if (length($read) < $long_read_end_pos)
										{
											$exclude_pos = length($read);
										}
										if ($hap_tag eq "HAP1" && $repetitive_detect1 eq "" && $repetitive_detect2 eq "dddddddd")
										{
											$exclude_reads_hap1_NP{$id_tmp3} = $position+$exclude_pos;
											#print {$filehandle{$seed_id2}} $id_tmp." EXCLUDE1\n";
										}
										elsif ($hap_tag eq "HAP2" && $repetitive_detect1 eq "" && $repetitive_detect2 eq "dddddddd")
										{
											$exclude_reads_hap2_NP{$id_tmp3} = $position+$exclude_pos;
											#print {$filehandle{$seed_id2}} $id_tmp." EXCLUDE2\n";
										}
#Exclude reads------------------------
									}
								}
								elsif ($score_match > $score_no_match*2 && $score_match > 1 && $score_match_DUP >= $score_no_match_DUP)
								{
									my $ext = substr $check_prev_splits{$id_tmp3}, $long_read_end_pos-90;
									$length_ext_all{$id_tmp3} = length($ext);
									
									if (length($ext) > 200)
									{
										$extensions{$ext} = $id_tmp3;
										$extensions2{$id_tmp3} = $ext;
										$extensions2b{$id_tmp3}{$ext}{$score_match} = $score_no_match;
			
										$var_matches{$id_tmp3} = ($score_match)-($score_no_match*1.7);
										$var_matches_DUP{$id_tmp3} = $score_match_DUP-($score_no_match_DUP*1.5);
									}
									else
									{
										print {$filehandle{$seed_id2}} length($ext)." LENGTH TOO SHORT\n";
									}
								}
								else
								{
									my $ext = substr $check_prev_splits{$id_tmp3}, $long_read_end_pos-90;
									$length_ext_all{$id_tmp3} = length($ext);
									
									if ($find_haps_in_seed eq "yes")
									{
										$ext = substr $check_prev_splits{$id_tmp3}, $long_read_start_pos, $long_read_end_pos-$long_read_start_pos;
									}
									
									if (length($ext) > 200)
									{                         
										$extensions_unknown{$ext} = $id_tmp3;
										$extensions_unknown2{$id_tmp3} = $ext;
										
										$var_matches{$id_tmp3} = ($score_match)-($score_no_match*1.7);
										$var_matches_DUP{$id_tmp3} = $score_match_DUP-($score_no_match_DUP*1.5);
									}
								}
								my $input_files_blast_tmp = keys %input_files_blast;
								if ($count_BLAST_runnng < $input_files_blast_tmp)
								{
									goto BLAST_RUN_NP;
								}
							}
						}
					}
                }
				
				
				my $count_remaining_ids = '0';
				foreach my $length_tmp (keys %input_files_blast)
				{
					my $count_remaining_ids2 = keys %{$input_files_blast{$length_tmp}};
					if ($count_remaining_ids2 > 0)
					{
						$count_remaining_ids++;
					}
					else
					{
						delete $input_files_blast{$length_tmp};
						delete $id_by_length{$length_tmp};
					}
				}

                my $limi = '1150';
                if ($count_remaining_ids > 0 && $count_limit1 < $limi)
                {
                    $count_limit1++;
					if ($count_limit1 > 1000)
					{
						sleep(1);
					}
                    goto BLAST_RUN_NP;
                }
                else
                {
                    foreach my $length_tmp (keys %input_files_blast)
					{
                        my $id_tmp3 = $id_by_length{$length_tmp};
						print {$filehandle{$seed_id2}} $id_tmp3." SKIPPED_IDs\n";
                    }
                }
                print {$filehandle{$seed_id2}} "\n".$ref_end_rejection." REF_END_REJECTION\n";
                print {$filehandle{$seed_id2}} $read_start_pos_rejection." READ_START_POS_REJECTION\n";
                print {$filehandle{$seed_id2}} $no_hit_rejection." NO_HIT_REJECTION\n";
                print {$filehandle{$seed_id2}} $almost_no_hit_rejection." ALMOST_NO_HIT_REJECTION\n";
                print {$filehandle{$seed_id2}} $accuracy_rejection." ACCURACY_REJECTION\n";
                print {$filehandle{$seed_id2}} $ext_too_short." EXT_TOO_SHORT\n";
                print {$filehandle{$seed_id2}} $overlap_too_short." OVERLAP_TOO_SHORT\n";
                print {$filehandle{$seed_id2}} $multi_match_rejection." MULTI_MATCH_REJECTION\n";
    
                print {$filehandle{$seed_id2}} "\n".$longest_read." LONGEST_READ\n";
                print {$filehandle{$seed_id2}} $count_matches_with_high_scores." HIGH_SCORE_matches\n\n";
		 
				if ($assembly_length_max eq "WG" && $y eq "1" && ($count_matches_with_high_scores > $sequencing_depth_NP*1.5 || ($total_matches < 10 && $total_matches < $sequencing_depth_NP*0.7)) && ($first_back_assembly eq "" || length($read) < 5000))
				{
					foreach my $id_tmp (keys %printed_refs2)
					{
						unlink $TMP_directory."ref_tmp_".$id_tmp."_".$y.".fasta";
					}
					foreach my $id_tmp (keys %printed_blast)
					{
						unlink $TMP_directory."blast_tmp_".$id."_".$y."__".$id_tmp.".txt";
						
						if (exists($double_matches{$id_tmp}))
						{
							unlink $TMP_directory."blast_tmp_".$id."_".$y."_yes_".$id_tmp.".txt";
						}
					}
					foreach my $id_tmp (keys %second_try)
					{
						unlink $TMP_directory."blast_tmp_".$id."_".$y."__".$id_tmp."_2.txt"
					}
					$first_back_assembly = "yes";
					goto END1;
				}  

                if ($second_try eq "" && $count_remaining_ids > 0)
                {
                    $second_try = "yes";
					foreach my $length_tmp (sort {$a <=> $b} keys %input_files_blast)
					{                
						my $id_tmp3 = $id_by_length{$length_tmp};

                        print {$filehandle{$seed_id2}} $id_tmp3." SECOND TRY\n";
                        my $ref_file = $printed_blast{$id_tmp3};
                        my @split_tmp = split /a/, $id_tmp3;
                        my $output_file1  = $tmp_sequences_directory_NP.$output_path_test.$split_tmp[0]."a".$output_path_test.$split_tmp[2].$output_path_test."sequence_tmp_NP_".$id_tmp3.".fasta";
                        if (exists($reverse_list{$id_tmp3}))
                        {
                            $output_file1  = $tmp_sequences_directory_NP.$output_path_test.$split_tmp[0]."a".$output_path_test.$split_tmp[2].$output_path_test."sequence_tmp_NP_".$id_tmp3."_rev.fasta";
                        }
                        $second_try{$id_tmp3} = undef;
                        my $command = "blastn -query ".$ref_file." -subject ".$output_file1." -out ".$TMP_directory."blast_tmp_".$id."_".$y."__".$id_tmp3."_2.txt -outfmt 5 -gapextend 1 -gapopen 2 -culling_limit 1 -strand plus -evalue 1e-50";
                        #my $command = "blastn -task dc-megablast -query ".$ref_file." -subject ".$output_file1." -out ".$output_path."blast_tmp_".$id."_".$y."__".$check_prev_splits_id."_2.txt -outfmt 5 -culling_limit 1 -strand plus -evalue 1e-50";
                    
                        syscmd($command);
                    }
                    goto BLAST_RUN_NP;
                }
                elsif ($count_remaining_ids > 0)
                {
                    foreach my $length_tmp (sort {$a <=> $b} keys %input_files_blast)
					{                
						my $id_tmp3 = $id_by_length{$length_tmp};
                        print {$filehandle{$seed_id2}} $id_tmp3." SKIPPED_IDs2\n";
                    }
                }

                foreach my $id_tmp (keys %printed_refs2)
                {
                    unlink $TMP_directory."ref_tmp_".$id_tmp."_".$y.".fasta";
                }
                foreach my $id_tmp (keys %printed_blast)
                {
                    unlink $TMP_directory."blast_tmp_".$id."_".$y."__".$id_tmp.".txt";
                    
                    if (exists($double_matches{$id_tmp}))
                    {
                        unlink $TMP_directory."blast_tmp_".$id."_".$y."_yes_".$id_tmp.".txt";
                    }
                }
                foreach my $id_tmp (keys %second_try)
                {
                    unlink $TMP_directory."blast_tmp_".$id."_".$y."__".$id_tmp."_2.txt"
                }
				
				if ($assembly_length_max eq "WG" && $y eq "1" && $read_start_pos_rejection > $count_matches_with_high_scores && ($first_back_assembly eq "" || length($read) < 5000))
				{
					$first_back_assembly = "yes";
					goto END1;
				}
#ADD_REJ_READS-----------------------------------------------------------------------------------------------------------------               
                if ($add_rejected_reads ne "" || $add_no_match_reads ne "")
                {
                    goto SELECT_LENGTH_NP2a;
                }        
#Remove multi match patterns-----------------------------------------------------------------------------------------------------------------                  
                my $gap_pos_prev = "";
                my $gap_id_prev = "";
                my $count_gap_pattern = '1';
                my $count_extra = '0';
                my %multi_match_remove;
                undef %multi_match_remove;
                foreach my $gap_pos_tmp (sort {$a <=> $b} keys %multi_match_pos)
                {             
                    foreach my $id_tmp (keys %alignment_length_save)
                    {
                        if (exists($multi_match{$id_tmp}))
						{}
						elsif ($alignment_length_save{$id_tmp} > $gap_pos_tmp+500)
                        {
                            $count_extra++;
                        }
                    }
                    if ($gap_pos_prev ne "" && $gap_pos_tmp < $gap_pos_prev+150 && $count_extra > 2)
                    {
                        $count_gap_pattern++;
                        $multi_match_remove{$gap_id_prev} = $gap_pos_tmp;
                        $multi_match_remove{$multi_match_pos{$gap_pos_tmp}} = $gap_pos_tmp;
                    }
                    elsif ($count_gap_pattern > 1 && $count_extra > $count_gap_pattern+1)
                    {
                        foreach my $id_tmp2 (keys %multi_match_remove)
                        {               
                            $extensions_nomatch{$extensions2{$id_tmp2}} = $id_tmp2;
							$extensions_nomatch2{$id_tmp2} = $extensions2{$id_tmp2};
							delete $extensions{$extensions2{$id_tmp2}};
                            delete $extensions2{$id_tmp2};
							delete $extensions2b{$id_tmp2};
							delete $extensions_unknown{$extensions2{$id_tmp2}};
							delete $extensions_unknown2{$id_tmp2};
                            delete $save_reads_for_next{$id_tmp2};
							$read_start_pos_rej{$id_tmp2} = $multi_match_remove{$id_tmp2};
							print {$filehandle{$seed_id2}} $id_tmp2." ".$multi_match_remove{$id_tmp2}." REMOVE_MULTI_MATCH\n";
							
                        }
                        print {$filehandle{$seed_id2}} $count_gap_pattern." REMOVE_MULTI_MATCH\n";
                        undef %multi_match_remove;
                        $count_gap_pattern = '1';
                        $count_extra = '0';
                    }
                    else
                    {
                        undef %multi_match_remove;
                        $count_gap_pattern = '1';
                        $count_extra = '0';
                    }

                    $gap_pos_prev = $gap_pos_tmp;
                    $gap_id_prev = $multi_match_pos{$gap_pos_tmp};
                }
                if ($count_gap_pattern > 1 && $count_extra > $count_gap_pattern+1)
                {
                    foreach my $id_tmp2 (keys %multi_match_remove)
                    {               
						$extensions_nomatch{$extensions2{$id_tmp2}} = $id_tmp2;
						$extensions_nomatch2{$id_tmp2} = $extensions2{$id_tmp2};
						delete $extensions{$extensions2{$id_tmp2}};
						delete $extensions2{$id_tmp2};
						delete $extensions2b{$id_tmp2};
						delete $extensions_unknown{$extensions2{$id_tmp2}};
						delete $extensions_unknown2{$id_tmp2};
						delete $save_reads_for_next{$id_tmp2};
                    }
                     print {$filehandle{$seed_id2}} $count_gap_pattern." REMOVE_MULTI_MATCH1\n";
                }

#-----------------------------------------------------------------------------------------------------------------  
                my $time_BLAST2 = time;
                my $time8 = $time_BLAST2 - $time_kmers;
                print {$filehandle{$seed_id2}} $time8." TIME_BLAST2\n\n";
                
                #if ($confirmed_reads_count_NP > 4 && $skip_confirmed eq "" && ($count_matches_with_high_scores < 4 || $count_matches_with_high_scores < $sequencing_depth_NP/4)
                    #&& $find_haps_in_seed eq "" && $only_confirmed eq "yes")
                #{
                    #$skip_confirmed = "yes";
                    #undef %id_matches;
                    #undef %reverse_list;
                    #undef %double_matches;
                   # goto SKIP_CONFIRMED_NP;
                #}
				if ($count_matches_with_high_scores < 5  && $full_reset_NP eq "")
				{
					if ($hap_tag eq "HAP1")
					{
						undef %exclude_reads_hap1_NP;
					}
					elsif ($hap_tag eq "HAP2")
					{
						undef %exclude_reads_hap2_NP;
					}
					undef %save_alignment_data_NP;
					undef %rejected_alignment_data_NP;
					$full_reset_NP = "A";
					print {$filehandle{$seed_id2}} $full_reset_NP." FULL_RESET0\n";
					$best_extension = "";
					$seed_id = $id;	
					$y++;
					$y{$id} = $y;
					goto FULL_RESET;
				}
                if ($count_matches_with_high_scores < 5 && $discontiguous_blast eq "" && $find_haps_in_seed eq "" && $sequencing_depth_NP < 80)
                {
                    print {$filehandle{$seed_id2}} "DISCONTIGUOUS_BLAST\n";
                    $discontiguous_blast = "yes";
                    undef %id_matches;
                    undef %reverse_list;
                    undef %double_matches;
                    goto SKIP_CONFIRMED_NP;
                }
				
				if ($count_matches_with_high_scores < $sequencing_depth_NP*1.5 && $position > 2000)
                {
                    $last_non_complex_region{$seed_id} = $position-2000;
                }
				
				if ($count_matches_with_high_scores < 5 && $count_matches_with_high_scores < $sequencing_depth_NP/2 && ($retry_NP eq "" || $retry_NP < $position-20000) && $position > 20000 && $PB_reads eq "" && $input_reads_DB_folder_PB eq "")
                {
                    $retry_NP = $position;
					if ($hap_tag eq "HAP1")
					{
						undef %exclude_reads_hap1_NP;
					}
					elsif ($hap_tag eq "HAP2")
					{
						undef %exclude_reads_hap2_NP;
					}
					undef %save_alignment_data_NP;
					undef %rejected_alignment_data_NP;
					$best_extension = "";
					$seed_id = $id;
					substr $read, -15000, 15000, "";
					$position -= 15000;
					$position{$id} = $position;
					$seed{$seed_id} = $read;
					$y++;
					$y{$id} = $y;
					goto FULL_RESET;
					print {$filehandle{$seed_id2}} "RETRY_NP\n";
                }				
                
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
#START reverse assembly----------------------------------------------------------------------------------------------------------------------------------------------------------                               
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
				
				if ($assembly_length_max eq "WdGsfs" && ($position_back > 0 || $position > $original_seed_length{$id}+10000))
				{
					if ($first_back_assembly eq "")
					{
						substr $read, 0, $original_seed_length{$id}, "";
						$first_back_assembly = "yes";
						$read_back = reverse($read);
						$read_back =~ tr/ACTG/TGAC/;
					}
		
					my $first_length_back = '900';
	
					if ($first_length_back > length($read_back))
					{
						$first_length_back = length($read_back)-50;
					}
	
					my $skip_confirmed_back = "";
					my $discontiguous_blast_back = "";
					my $last_hap_pos_back = "";
					my $last_hap_pos_DUP_back = "";
					my %hash_NP_reads_back_tmp;
					undef %hash_NP_reads_back_tmp;
					
SKIP_CONFIRMED_NP_BACK:                
					my %id_matches_back;
					undef %id_matches_back;
					my %double_matches_back;
					undef %double_matches_back;
					my %reverse_list_back;
					undef %reverse_list_back;
			  
	
#Find matches in hash against last 600 bp back--------------------------------------------------------------                           
						
					my $query_coverage_back = 90;
						
LONGER_LAST_600_NP_BACK:
	
					my $last_600 = substr $read_back, -$first_length_back;

					my %first_pos_read_back;
					undef %first_pos_read_back;
	
					if (keys %{$trace_back_split_NP{$id}} > 0)
						{
								#print {$filehandle{$seed_id2}} $last_600." LAST_600!!!!!!!!!!!!!!!!\n";         
						}
					
					my $query_file_DB = $TMP_directory."query_NP_back.fasta";
					open(INPUT_QUERY_NP, ">" .$query_file_DB) or die "\nCan't open file $query_file_DB, $!\n";
					#print INPUT_QUERY_NP ">ref\n";
					print INPUT_QUERY_NP $last_600;        
					close INPUT_QUERY_NP;
					
					my %output_files_DB_back;
					undef %output_files_DB_back;
					my %output_files_DB2_back;
					undef %output_files_DB2_back;
					
					foreach my $pid_tmp (keys %ret_data_NP)
					{
						my $file_tmp = $TMP_directory."blast_tmp_DB_NP_back_".$id."_".$y."_".$pid_tmp.".txt";
						my $command_DB = "blastn -query ".$query_file_DB." -db ".$ret_data_NP{$pid_tmp}." -out ".$file_tmp." -outfmt 7 -qcov_hsp_perc ".$query_coverage_back." -num_threads 2 &";
						syscmd($command_DB);             
						$output_files_DB_back{$file_tmp} = undef;
						$output_files_DB2_back{$file_tmp} = undef;
					}
		
#Add saved reads--------------------------------------------------------------                           						

SKIP_TO_CONFIRMED_NP_BACK:
					my %extensions_back;
					undef %extensions_back;
					my %extensions2_back;
					undef %extensions2_back;
					my %extensions2b_back;
					undef %extensions2b_back;
					my %extensions_nomatch_back;
					undef %extensions_nomatch_back;
					my %extensions_unknown_back;
					undef %extensions_unknown_back;
					my %extensions_nomatch2_back;
					undef %extensions_nomatch2_back;
					my %extensions_nomatch2b_back;
					undef %extensions_nomatch2b_back;                
					my %extensions_unknown2_back;
					undef %extensions_unknown2_back;
					my %store_mismatches_NP_back;
					undef %store_mismatches_NP_back;
					my %store_mismatches_all_NP_back;
					undef %store_mismatches_all_NP_back;
					my %store_mismatches_N_NP_back;
					undef %store_mismatches_N_NP_back;
					
					my $extensions_nomatch2b_count_back = '0';
					my $extensions_nomatch2b_count_saved_back = '0';
					my $total_matches_extra_back = '0';
					my $position_prev_back = "";
					my %alignment_length_saved_back;
					undef %alignment_length_saved_back;
					my %score_match_saved_back;
					undef %score_match_saved_back;
					my %score_no_match_saved_back;
					undef %score_no_match_saved_back;
					my %score_match_DUP_saved_back;
					undef %score_match_DUP_saved_back;
					my %score_no_match_DUP_saved_back;
					undef %score_no_match_DUP_saved_back;
					my %accuracy_saved_back;
					undef %accuracy_saved_back;               
					my %read_start_pos_rej_back;
					undef %read_start_pos_rej_back;
					my %read_start_pos_rej_saved_back;
					undef %read_start_pos_rej_saved_back;
					my %extensions_nomatch2b_saved_back;
					undef %extensions_nomatch2b_saved_back;
					
					foreach my $seed_id_tmp0 (keys %save_alignment_data_NP_back)
					{
						if ($seed_id_tmp0 eq $seed_id)
						{
							foreach my $id_tmp7 (keys %{$save_alignment_data_NP_back{$seed_id_tmp0}})
							{
								my @alignment_data = split /_/, $save_alignment_data_NP_back{$seed_id_tmp0}{$id_tmp7};
	
								if ($alignment_data[0] eq "yes" || $alignment_data[11] eq "")
								{ 
									$alignment_length_saved_back{$id_tmp7} = $alignment_data[4];
									$first_pos_read_back{$id_tmp7} = $alignment_data[2];
									$position_prev_back = $alignment_data[3];
									$score_match_saved_back{$id_tmp7} = $alignment_data[5];
									$score_no_match_saved_back{$id_tmp7} = $alignment_data[6];
									$score_match_DUP_saved_back{$id_tmp7} = $alignment_data[7];
									$score_no_match_DUP_saved_back{$id_tmp7} = $alignment_data[8];
									$accuracy_saved_back{$id_tmp7} = $alignment_data[9];
								}
							   
								if (($alignment_data[6] > $alignment_data[5] && $alignment_data[6] > 1) || ($alignment_data[8] > $alignment_data[7] && ($alignment_data[8] > 3
									|| $alignment_data[6]+$alignment_data[8] > $alignment_data[5]+$alignment_data[7])))
								{
									$extensions_nomatch2b_saved_back{$id_tmp7}{$alignment_data[5]}{$alignment_data[6]} = $alignment_data[3]+$alignment_data[10];
									$extensions_nomatch2b_count_saved_back++;
									if ($alignment_data[1] eq "yes")
									{
										$reverse_list_back{$id_tmp7} = undef;
									}
								}
								elsif ($alignment_data[0] ne "yes" && $alignment_data[11] ne "")
								{
									$read_start_pos_rej_back{$id_tmp7} = $alignment_data[11]+$position_back-$alignment_data[3];
									$read_start_pos_rej_saved_back{$id_tmp7} = undef;
									if ($alignment_data[1] eq "yes")
									{
										$reverse_list_back{$id_tmp7} = undef;
									}
								}
								else
								{
									$id_matches_back{$id_tmp7} = undef;
									$total_matches_extra_back++;
									if ($alignment_data[1] eq "yes")
									{
										$reverse_list_back{$id_tmp7} = undef;
									}
									elsif ($alignment_data[1] eq "yes2")
									{
										$reverse_list_back{$id_tmp7} = undef;
										$double_matches_back{$id_tmp7} = undef;
									}
									if ($alignment_data[12] ne "")
									{
										$hash_NP_reads_tmp_back{$id_tmp7} = $alignment_data[12];
									}
								}				
							}
						}
					}                                          
					print {$filehandle{$seed_id2}} $total_matches_extra_back." LAST_1000_matches_EXTRA_BACK\n\n";
					
					
					my $database_count_tmp = '0';
DB_RESULTS_NP_BACK: foreach my $blast_db_results_tmp (keys %output_files_DB_back)
					{
						if (-s $blast_db_results_tmp)
						{
							open(BLAST_RESULTS_DB_NP, $blast_db_results_tmp) or die "\nCan't open file $blast_db_results_tmp, $!\n";
							
							my $count_lines_tmp = '1';
							while (my $line_tmp = <BLAST_RESULTS_DB_NP>)
							{
								chomp($line_tmp);
								if ($count_lines_tmp eq '4' && $line_tmp eq "# 0 hits found")
								{
									delete $output_files_DB{$blast_db_results_tmp};
									next DB_RESULTS_NP;
								}
								elsif ($count_lines_tmp eq '5' && $line_tmp eq "# BLAST processed 1 queries")
								{
									delete $output_files_DB{$blast_db_results_tmp};
									next DB_RESULTS_NP;
								}
								elsif ($count_lines_tmp > 5 && $line_tmp eq "# BLAST processed 1 queries")
								{
									delete $output_files_DB{$blast_db_results_tmp};
									next DB_RESULTS_NP;
								}
								elsif ($count_lines_tmp > 5)
								{
									my @line_tmp = split /\t/, $line_tmp;
									my $id_tmp = $line_tmp[1];
									my $accuracy_tmp = $line_tmp[2];
									#my $alignment_length = $line_tmp[3];
									my $read_pos_start_tmp = $line_tmp[8];
									my $read_pos_end_tmp = $line_tmp[9];
									#if ($alignment_length < 0.9*$first_length_back)
									#{
										#next;
									#}
									if ($hap_tag eq "HAP1" && exists($exclude_reads_hap1_NP{$id_tmp}))
									{
										next;
									}
									if ($hap_tag eq "HAP2" && exists($exclude_reads_hap2_NP{$id_tmp}))
									{
										next;
									}
									if (exists($save_alignment_data_NP{$seed_id}{$id_tmp}))
									{
										next;
									}
									if ($accuracy_tmp > 70)
									{
										if (exists($id_matches{$id_tmp}))
										{
											if ($read_pos_end_tmp < $read_pos_start_tmp) 
											{
												$reverse_list{$id_tmp} = undef;
											}
											delete $first_pos_read{$id_tmp};
										}
										else
										{
											$id_matches{$id_tmp} = undef;
											if ($read_pos_end_tmp < $read_pos_start_tmp) 
											{
												$reverse_list{$id_tmp} = undef;
												$first_pos_read{$id_tmp} = $read_pos_start_tmp;
											}
											else
											{
												$first_pos_read{$id_tmp} = $read_pos_end_tmp;
											}
										}        
									}  
								}
								$count_lines_tmp++;
							}
							close BLAST_RESULTS_DB_NP;
						}
					}
					if (keys %output_files_DB > 0)
					{
						goto DB_RESULTS_NP;
					}
					my $time_kmers = time;
					my $time1 = $time_kmers - $time_START;
					print {$filehandle{$seed_id2}} $time1." TIME1\n";
					
					my $total_matches = keys %id_matches;
					print {$filehandle{$seed_id2}} "\n".$total_matches." LAST_1000_matches\n";
					
					if ($assembly_length_max eq "WG" && $y eq "1" && ($total_matches > $sequencing_depth_NP*5 || $total_matches < 5) && ($first_back_assembly eq "" || length($read) < 3000))
					{
						foreach my $blast_db_results_tmp (keys %output_files_DB2)
						{
							unlink $blast_db_results_tmp;
							delete $output_files_DB2{$blast_db_results_tmp};
						}
						$first_back_assembly = "yes";
						goto END1;
					}
					
					if ($total_matches < 10 && $total_matches < $sequencing_depth_NP && ($repetitive_detect1 ne "" || $repetitive_detect2 ne ""))
					{
						$first_length_back += 1000;
						if ($total_matches > 1000)
						{
							$first_length_back += 1000;
						}
						$query_coverage = '60';
						undef %id_matches;
						undef %double_matches;
						undef %reverse_list;
						undef %hash_NP_reads_tmp;
						foreach my $blast_db_results_tmp (keys %output_files_DB2)
						{
							unlink $blast_db_results_tmp;
						}
						goto LONGER_LAST_600_NP;
					}
					
					if ($total_matches > 350 && $first_length_back < 1000)
					{
						$first_length_back += 1000;
						if ($total_matches > 1000)
						{
							$first_length_back += 1000;
						}
						undef %id_matches;
						undef %double_matches;
						undef %reverse_list;
						undef %hash_NP_reads_tmp;
						foreach my $blast_db_results_tmp (keys %output_files_DB2)
						{
							unlink $blast_db_results_tmp;
						}
						goto LONGER_LAST_600_NP;
					}
																						   
					$total_matches += $total_matches_extra;
	
					if ($total_matches < $sequencing_depth_NP*2.3 && $position > 2000)
					{
						$last_non_complex_region{$seed_id} = $position-2000;
					}
				}
                
                                         
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
#Make consensus of extensions long reads NP-----------------------------------------------------------------------------------------------------------------------------------------                               
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
                my $mismatch_retry = '0';
                my $high_quality = "";
                my %split_positions_DUP_tmp;
                undef %split_positions_DUP_tmp;
                my $remove_reads_check = "";
                
MISMATCH_RETRY_NP:
                my $merged_ext = "";
                my %extensions2_tmp = %extensions2;
                my $ext2_count = keys %extensions2;
                print {$filehandle{$seed_id2}} $ext2_count." EXT_COUNT\n";
                
                #if ($skip_confirmed eq "" && ($ext2_count < 4 || $ext2_count < $sequencing_depth_NP/4)
                    #&& $find_haps_in_seed eq "" && $only_confirmed eq "yes")
                #{
                    #$skip_confirmed = "yes";
                    #undef %id_matches;
                    #undef %reverse_list;
                    #undef %double_matches;
                    #goto SKIP_CONFIRMED_NP;
                #}
                my $other_hap_read_count = keys %extensions_nomatch2;
                if ((($ext2_count < 3 || $ext2_count < $sequencing_depth_NP/2 || ($ext2_count < $sequencing_depth_NP && $count_matches_with_high_scores > $sequencing_depth_NP*3)
                      || ($ext2_count < $sequencing_depth_NP/1.5 && $ext2_count < 7) || ($ext2_count < ($count_matches_with_high_scores-$other_hap_read_count)/2))
                    && (keys %extensions_unknown > 0 || $other_hap_read_count > 0)) || ($last_hap_pos > 10000 && $ext2_count < $sequencing_depth_NP/2))
                {
                    %extensions = (%extensions, %extensions_unknown);
                    %extensions2_tmp = (%extensions2, %extensions_unknown2);
                    $merged_ext = "1";
                    print {$filehandle{$seed_id2}} " MERGE1\n";
                    if (keys %extensions < 3)
                    {
                        %extensions = (%extensions, %extensions_nomatch);
                        %extensions2_tmp = (%extensions2, %extensions_nomatch2);
                        $merged_ext = "2";
                        print {$filehandle{$seed_id2}} " MERGE2\n";
                    }
                }
                #if (keys %extensions < 3 && $confirmed_reads_count_NP > 4 && $skip_confirmed eq "" && $only_confirmed eq "yes")
                #{
                    #$skip_confirmed = "yes";
                    #$best_extension = "";
                    #goto SKIP_CONFIRMED_NP;
                #}
                
                my $longer_extension_for_repeat = "";
                my $length_extension = '1500';
				my $count_hap_matches = keys %extensions2_tmp;

#Reduce extensions-------------------------------------------------------------------------------------------------                            
                
REDUCE_EXTENSIONS_NP:

                my %scores;
                undef %scores;
                my %scores2;
                undef %scores2;
                
                my $acc_max = '0';
                my $acc_min = '100';
                foreach my $accuracy_id_tmp (sort {$a <=> $b} keys %accuracy)
                {
                    if (exists($extensions2_tmp{$accuracy_id_tmp}))
                    {
                        my $accuracy_tmp = $accuracy{$accuracy_id_tmp};
                        if ($accuracy_tmp > $acc_max)
                        {
                            $acc_max = $accuracy_tmp;
                        }
                        if ($accuracy_tmp < $acc_min)
                        {
                            $acc_min = $accuracy_tmp;
                        }
                    }
                }
                foreach my $accuracy_id_tmp (sort {$a <=> $b} keys %accuracy)
                {
                    if ((exists($extensions2_tmp{$accuracy_id_tmp})) && $acc_max-$acc_min > 0)
                    {
                        my $accuracy_tmp = $accuracy{$accuracy_id_tmp};
                        $scores{$accuracy_id_tmp} = (($accuracy_tmp-$acc_min)/($acc_max-$acc_min))/2;
                    }
                }
                
                my $overlap_max = '0';
                my $overlap_min = '0';
                foreach my $overlap_id_tmp (keys %alignment_length_save)
                {          
                    if (exists($multi_match{$overlap_id_tmp}))
                    {
                    }
                    elsif (exists($extensions2_tmp{$overlap_id_tmp}))
                    {
                        my $overlap_tmp = $alignment_length_save{$overlap_id_tmp};
                        if ($overlap_tmp > $overlap_max)
                        {
                            $overlap_max = $overlap_tmp;
                        }
                    }
                }
                foreach my $overlap_id_tmp (keys %alignment_length_save)
                {
                    if (exists($extensions2_tmp{$overlap_id_tmp}))
                    {
                        my $overlap_tmp = $alignment_length_save{$overlap_id_tmp};
						if (exists($multi_match{$overlap_id_tmp}))
                        {
                            $overlap_tmp *= 0.8;
                        }
                        my $overlap_score_tmp = '0';
                        if (($overlap_max-$overlap_min) ne '0')
                        {
                            $overlap_score_tmp = ($overlap_tmp-$overlap_min)/($overlap_max-$overlap_min);
                        }
                        my $score_tmp = $scores{$overlap_id_tmp};
                        
                        $scores{$overlap_id_tmp} = $overlap_score_tmp+$score_tmp;
                    }
                }
                
                my $var_max = '0';
                my $var_min = '100000000000000';
                
                foreach my $var_matches_id_tmp (sort {$a <=> $b} keys %var_matches)
                {     
                    my $var_matches_tmp = $var_matches{$var_matches_id_tmp};
                    if (exists($extensions2_tmp{$var_matches_id_tmp}))
                    {
                        if ($var_matches_tmp > $var_max)
                        {
                            $var_max = $var_matches_tmp;
                        }
                        if ($var_matches_tmp < $var_min)
                        {
                            $var_min = $var_matches_tmp;
                        }
                    }
                }
                foreach my $var_matches_id_tmp (sort {$a <=> $b} keys %var_matches)
                {     
                    my $var_matches_tmp = $var_matches{$var_matches_id_tmp};
                    if (exists($extensions2_tmp{$var_matches_id_tmp}))
                    {
                        my $var_score_tmp = '0';
                        if ($var_max-$var_min > 0)
                        {
                            $var_score_tmp = ($var_matches_tmp-$var_min)/($var_max-$var_min);
                        }
                        my $score_tmp = $scores{$var_matches_id_tmp};
                        $scores{$var_matches_id_tmp} = $score_tmp+$var_score_tmp;
                    }
                }
                
                my $var_max_DUP = '0';
                my $var_min_DUP = '100000000000000';
                
                foreach my $var_matches_id_tmp (sort {$a <=> $b} keys %var_matches_DUP)
                {     
                    my $var_matches_tmp = $var_matches{$var_matches_id_tmp};
                    if (exists($extensions2_tmp{$var_matches_id_tmp}))
                    {
                        if ($var_matches_tmp > $var_max_DUP)
                        {
                            $var_max_DUP = $var_matches_tmp;
                        }
                        if ($var_matches_tmp < $var_min_DUP)
                        {
                            $var_min_DUP = $var_matches_tmp;
                        }
                    }
                }
                foreach my $var_matches_id_tmp (sort {$a <=> $b} keys %var_matches_DUP)
                {     
                    my $var_matches_tmp = $var_matches{$var_matches_id_tmp};
                    if (exists($extensions2_tmp{$var_matches_id_tmp}))
                    {
                        my $var_score_tmp = '0';
                        if ($var_max_DUP-$var_min_DUP > 0)
                        {
                            $var_score_tmp = ($var_matches_tmp-$var_min_DUP)/($var_max_DUP-$var_min_DUP);
                        }
                        my $score_tmp = $scores{$var_matches_id_tmp};
                        $scores{$var_matches_id_tmp} = $score_tmp+$var_score_tmp;
                    }
                }
                      
                foreach my $score_id (keys %scores)
                {
                    if (exists($scores2{$scores{$score_id}}))
                    {
                        $scores2{$scores{$score_id}} .= ",$score_id"; 
                    }
                    else
                    {
                        $scores2{$scores{$score_id}} = $score_id;
                    }
                }

                my $c = '1';
                my $c2 = '0';
                my $limit2 = $sequencing_depth_NP;
				my $limit = 50;

                my $shortest_ext = '100000000000000';
                my $max_score = "";
				my $max_score3 = "";
                my $max_score_count = '0';
                my %length_ext0;
                undef %length_ext0;
				my %length_ext_10;
                undef %length_ext_10;
                
                foreach my $score_tmp (sort {$b <=> $a} keys %scores2)
                {
                    print {$filehandle{$seed_id2}} $scores2{$score_tmp}." ".$score_tmp." SCORE\n";             
                    my @ids_tmp = split /,/, $scores2{$score_tmp};
                    foreach my $ids_tmp (@ids_tmp)
                    {
						if (exists($extensions2_tmp{$ids_tmp}))
                        {
                            $max_score_count++;
                            if ($max_score eq "")
                            {
                                $max_score = $score_tmp;
                            }
							if ($max_score_count eq '3')
                            {
                                $max_score3 = $score_tmp;
                            }
                            
                            if ($max_score_count < 3)
                            {
                                $c2++;
                                $length_ext0{length($extensions2_tmp{$ids_tmp})} = $ids_tmp;
                            }
                            elsif (($c2 > $limit || $c2 > $limit2*1.8) && $c2 > 2 && $find_haps_in_seed eq "")
                            {
                                delete $extensions2_tmp{$ids_tmp};
                            }
							elsif ($c2 > $limit2*0.8 && $c2 > 10 && $find_haps_in_seed eq "" && $score_tmp < $max_score/3 && $score_tmp < 0.5 && $score_tmp < $max_score3*0.4)
                            {
                                delete $extensions2_tmp{$ids_tmp};
                            }
                            elsif ($c >= 5 && length($extensions2_tmp{$ids_tmp}) < ($shortest_ext*0.7) && length($extensions2_tmp{$ids_tmp}) < 1000 && $c2 > 2 && $find_haps_in_seed eq "")
                            {   
                                print {$filehandle{$seed_id2}} $ids_tmp." DEL\n";
                                print {$filehandle{$seed_id2}} length($extensions2_tmp{$ids_tmp})." LENGTH\n";
                                delete $extensions2_tmp{$ids_tmp};
                            }
                            elsif ($c < $limit2 || $score_tmp > 0.4*$max_score || $score_tmp > 0.6 || $c2 < 3 || ($count_matches_with_high_scores > $limit2 && $c2 < $sequencing_depth_NP*1.5))
                            {
                                $c2++;
                                $length_ext0{length($extensions2_tmp{$ids_tmp})} = $ids_tmp;
                            }
                            else
                            {
                                delete $extensions2_tmp{$ids_tmp};
                            }
                            
                            if ($c < 5 && length($extensions2_tmp{$ids_tmp}) < $shortest_ext && length($extensions2_tmp{$ids_tmp}) ne "")
                            {
                                $shortest_ext = length($extensions2_tmp{$ids_tmp});
                                print {$filehandle{$seed_id2}} $ids_tmp." ID\n";
                                print {$filehandle{$seed_id2}} $shortest_ext." SHORTEST_EXT\n";
                            }
							
							if ($c <= 10)
							{
								$length_ext_10{length($extensions2_tmp{$ids_tmp})} = $ids_tmp;
							}
                            $c++; 
                        }
                    } 
                }
#Increase consensus length-----------------------------------------------------------------------------------------------------            
                
SELECT_LENGTH_NP:               
                if ($PB_reads ne "" || $input_reads_DB_folder_PB ne "")
                {
                    $length_extension = '500';
                }

                if (exists($cut_repeat_seq{$seed_id}))
                {
                    $length_extension = '3900';
                }
                if ($longer_extension_for_repeat eq "" && $shortest_ext ne '100000000000000')
                {
                    $length_extension = $shortest_ext;
                }

                if ($length_extension < 800)
                {
                    $length_extension = '800';
                }
				if ($length_extension > 10000 && ($count_matches_with_high_scores < $sequencing_depth_NP*5 || $total_matches < $sequencing_depth_NP*10))
                {
                    $length_extension = '10000';
                }

                if ($find_haps_in_seed eq "yes")
                {
                    $length_extension = length($read)-25;
                }
                
                my %length_tmp;
                undef %length_tmp;
                my $total_tmp = keys %extensions2_tmp;
                foreach my $ids_tmp (keys %extensions2_tmp)
                {
                    $length_tmp{length($extensions2_tmp{$ids_tmp})} = undef;
                }
                my $count_tmp = '1';
                foreach my $length_ext_tmp (sort {$b <=> $a} keys %length_tmp)
                {
                    if (($count_tmp > 0.75*$total_tmp || $count_tmp > 25) && $count_tmp > 3)
                    {
                        if ($length_ext_tmp > $length_extension && $find_haps_in_seed eq "" && $length_ext_tmp < 10000)
                        {
                            $length_extension = $length_ext_tmp;
                        }
                        last;
                    }
                    $count_tmp++;
                }
				my $count_tmp2 = '1';
				foreach my $length_ext_tmp (sort {$b <=> $a} keys %length_ext_10)
                {
                    if ($count_tmp2 > 6)
                    {
                        if ($length_ext_tmp > $length_extension && $find_haps_in_seed eq "" && $length_ext_tmp < 10000)
                        {
                            $length_extension = $length_ext_tmp;
                        }
                        last;
                    }
                    $count_tmp2++;
                }
                $length_extension += 450;
                    
                if ($longer_extension_for_repeat ne "")
                {
                    my $count_tmp3 = '1';
                    foreach my $length_ext_tmp (sort {$b <=> $a} keys %length_tmp)
                    {
                        if ($count_tmp3 > 0.25*$total_tmp && $count_tmp3 > 3)
                        {
                            $length_extension = $length_ext_tmp;
                            last;
                        }
                        $count_tmp3++;
                    }
                    if ($length_extension < $longer_extension_for_repeat)
                    {                      
                    }
                    else
                    {
                        $length_extension = $longer_extension_for_repeat
                    }
                    
                }

                print {$filehandle{$seed_id2}} $length_extension." LENGTH_EXTENSION\n";
                
SELECT_LENGTH_NP2a:
#ADD_REJ_READS-----------------------------------------------------------------------------------------------------------------               
                if (($add_rejected_reads ne "" || $add_no_match_reads ne "") && $remove_reads_check eq "")
                {
                    foreach my $id_tmp1 (keys %add_rej_reads_extra)
                    {
                        my $long_read_tmp = "";
                        if (exists($reverse_list{$id_tmp1}))
                        {
                            $long_read_tmp = reverse($hash_NP_reads_tmp{$id_tmp1});
                            $long_read_tmp =~ tr/ACTG/TGAC/;
                        }
                        else
                        {
                            $long_read_tmp = $hash_NP_reads_tmp{$id_tmp1};
                        }
                        
                        my $long_read_end_pos_tmp = $long_read_end_pos_save{$id_tmp1};
                        print {$filehandle{$seed_id2}} $id_tmp1." IDD_ADD_REJ\n";
                        print {$filehandle{$seed_id2}} $long_read_end_pos_tmp." LONG_READ_END\n";
                        my $ext = substr $long_read_tmp, $long_read_end_pos_tmp-90, $length_extension;
                        if (exists($scores2{'0'}))
                        {
                            $scores2{'0'} .= ",$id_tmp1"; 
                        }
                        else
                        {
                            $scores2{'0'} = $id_tmp1;
                        }
                        $extensions2_tmp{$id_tmp1} = $ext;
                    }
                }        
#------------------------------------------------------------------------------------------------------------------------------
SELECT_LENGTH_NP2:
				my $extension_part_length = '700';
				my $devide_extension_count = int(($length_extension/$extension_part_length)+0.5);
				#if ($devide_extension_count > $maxProcs)
				#{
					#$devide_extension_count = $maxProcs;
				#}
				if ($devide_extension_count < 1)
				{
					$devide_extension_count = '1';
				}
				my $length_extension_part = int($length_extension/$devide_extension_count);
				my $length_extension_part_extra = '100';
				
				print {$filehandle{$seed_id2}} $length_extension_part." LENGTH_EXTENSION_PART\n";
				
                my $v = '1';
                my %length_ext;
                undef %length_ext;
                my %rank_to_id;
                undef %rank_to_id;
                my %id_to_rank;
                undef %id_to_rank;
                
                foreach my $score_tmp (sort {$b <=> $a} keys %scores2)
                {
                    my @ids_tmp = split /,/, $scores2{$score_tmp};
                    foreach my $ids_tmp (@ids_tmp)
                    {
                        if (exists($extensions2_tmp{$ids_tmp}))
                        {                     
                            my $ext_tmp2 = substr $extensions2_tmp{$ids_tmp}, 0, $length_extension;
							$rank_to_id{$v} = $ids_tmp;
                            $length_ext{$v} = length($ext_tmp2);                        
							$save_reads_for_next{$ids_tmp} = undef;   
                            $id_to_rank{$ids_tmp} = $v;
                            print {$filehandle{$seed_id2}} $v." V ".$ids_tmp." ID ".length($ext_tmp2)." LENGTH\n";
                            if (exists($alignment_length_save{$ids_tmp}))
                            {
                                if ($score_tmp eq '0')
                                {}
                                elsif ($alignment_length_save{$ids_tmp} > $position-$last_non_complex_region{$seed_id} && ($alignment_length_save{$ids_tmp} > $repetitive_detect2 || $repetitive_detect2 eq ""))
                                {
                                    $span_complex_region{$v} = $alignment_length_save{$ids_tmp}
                                }
                            }
                            $v++;
                        }
                    }
                }
				my $total_nuc_count_original = $v-1;
				
				my $clipped_ext = "";
				my $N = '0';
				my $N_resolved = '0';
				my $mafft_count = "1";
				my $cp_original = "";
				my %track_length_ext_total;
                undef %track_length_ext_total;
				my %SNP_patterns_prev;
                undef %SNP_patterns_prev;
				my %pos_pattern_list;
                undef %pos_pattern_list;
                my %read_patterns2;
                undef %read_patterns2;
                my %read_patterns_final;
                undef %read_patterns_final;
				my $first_split_pos = "";
                my %track_mismatch_ext;
                undef %track_mismatch_ext;
                my $track_mismatch_count = '0';
				my $post_pattern_match = "";
				my $post_pattern_match_count = '1';
				my $post_pattern_match_save = "";
				my $nuc_match = "";
                my $nuc_prev = "";
				my $total_count_prev_patterns = '0';
				my %total_count_prev_patterns;
				undef %total_count_prev_patterns;
				my %quality_scores_tmp;
                undef %quality_scores_tmp;
                my %quality_scores_gap_tmp;
                undef %quality_scores_gap_tmp;
				my %best_read_score;
                undef %best_read_score;
				
				my %mismatches_tmp;
                undef %mismatches_tmp;
                my %mismatches_tmp_all;
                undef %mismatches_tmp_all;
				my %reads_mismatch;
				undef %reads_mismatch;
				my %reads_mismatchb;
				undef %reads_mismatchb;
				my %reads_mismatch2;
				undef %reads_mismatch2;
				my %reads_mismatch2_tmp;
				undef %reads_mismatch2_tmp;
				my $mismatch_score = '0';                      
				my %reads_mismatch_all;
				undef %reads_mismatch_all;
				my %reads_mismatchb_all;
				undef %reads_mismatchb_all;
				my %reads_mismatch2_all;
				undef %reads_mismatch2_all;
				my %reads_mismatch2_tmp_all;
				undef %reads_mismatch2_tmp_all;
				my $mismatch_score_all = '0';
                my $mismatches_tmp_check = "";                     
                my $pos_pattern_list_check = "";
				my $time_BLAST3b = time;
				my %extensions_seed;
                undef %extensions_seed;
				my %haps_list;
                undef %haps_list;
                my $hap_position = "";
				my %clipped_ext_pos;
				undef %clipped_ext_pos;
MAFFT_NP:				
				my $other_seq = "";
				foreach my $rank_tmp7 (sort {$a <=> $b} keys %rank_to_id)
				{
					if (exists($extensions2_tmp{$rank_to_id{$rank_tmp7}}))
					{                     
						my $start_pos_tmp = 0;
						if ($mafft_count > 1)
						{
							$start_pos_tmp = $track_length_ext_total{$rank_tmp7};
						}
						if ($length_ext{$rank_tmp7} > $start_pos_tmp+500 || $mafft_count eq '1')
						{
							my $ext_tmp2 = substr $extensions2_tmp{$rank_to_id{$rank_tmp7}}, $start_pos_tmp, $length_extension_part+$length_extension_part_extra;
							$other_seq .= ">".$rank_tmp7."\n";
							$other_seq .= $ext_tmp2."\n";
						}
					}
				}
				my $time_mafft1 = time;
					
				my $output_file6  = $TMP_directory."sequence_tmp_".$mafft_count.".fasta";
			
				open(OUTPUT_LONG1, ">" .$output_file6) or die "\nCan't open file $output_file6, $!\n";
				print OUTPUT_LONG1 $other_seq;
				
				close OUTPUT_LONG1;
				
				chomp($output_file6);
				
				$print_sep = $y."_".$id;
			
				#my $cmd = sprintf("blastn -query %s -subject %s -out blast_tmp3.txt -outfmt 4", $ref_file, $output_file1);
				#my $cmd = sprintf("blastn -query %s -subject %s -out blast_tmp3.txt -reward 1 -penalty -2 -gapopen 2 -gapextend 2 -outfmt 4", $ref_file, $output_file1);
				#my $cmd = sprintf("muscle -in %s -out blast_tmp3.txt -maxiters 1 -diags", $output_file1);
				my $cmd = "";
				
				if ($high_quality eq "yes")
				{
					$high_quality = "yes2";
					print {$filehandle{$seed_id2}} "HIGH_QUALITY_MAFFT\n";
					$cmd = sprintf("mafft --op 0.2 --thread 4 --quiet --clustalout --maxiterate 100 --globalpair %s > ".$TMP_directory."mafft_tmp_NP_".$print_sep."_".$mismatch_retry."_".$high_quality."_".$mafft_count.".txt", $output_file6);
				}
				elsif (keys %extensions > 2)
				{
					$cmd = sprintf("mafft --op 0.5 --thread 4 --quiet --clustalout %s > ".$TMP_directory."mafft_tmp_NP_".$print_sep."_".$mismatch_retry."_".$high_quality."_".$mafft_count.".txt", $output_file6);
				}
				else
				{
					$cmd = sprintf("mafft --op 1.01 --ep 1.2 --thread 4 --quiet --clustalout %s > ".$TMP_directory."mafft_tmp_NP_".$print_sep."_".$mismatch_retry."_".$high_quality."_".$mafft_count.".txt", $output_file6);
				}
				
				my $system_result = syscmd($cmd);                         
				
				my $mafft_output_tmp = $TMP_directory."mafft_tmp_NP_".$print_sep."_".$mismatch_retry."_".$high_quality."_".$mafft_count.".txt";
				#$mafft_output{$mafft_count} = $output_tmp;	
				$mafft_count++;
				
				 my $time_maff3 = time;
                my $time_mafft = $time_maff3 - $time_mafft1;
                print {$filehandle{$seed_id2}} $time_mafft." TIME_MAFFT\n\n";
				unlink $output_file6;

#merge mafft lines-------------------------------------------------------------------------------             

				my $e = '0';
                my $query_line = "";
                my %subject_list;
                undef %subject_list;
                my %subject_list_original;
                undef %subject_list_original;
                my $consensus_total = "";
                my %gaps_id;
                undef %gaps_id;
                my %length_id;
                undef %length_id;
                my $cp = '20';
               
				open(INPUT_BLAST3, $mafft_output_tmp) or print "\n\nCan't open blast file $mafft_output_tmp, $!\n";

				while (my $line2 = <INPUT_BLAST3>)
				{                                                     
					chomp($line2);
					if ($e > 2)
					{
						my @blast_result_tmp = split /\s+/, $line2;
						
						my $read_id_tmp = $blast_result_tmp[0];
						
						if (exists($subject_list{$read_id_tmp}))
						{
							my $subject_tmp = $subject_list{$read_id_tmp};
							$subject_list{$read_id_tmp} = $subject_tmp.$blast_result_tmp[1];
						}
						elsif ($read_id_tmp =~ m/\d+/)
						{                                   
							$subject_list{$read_id_tmp} = $blast_result_tmp[1];
						}
						if ($query_line eq "yes")
						{
							my $consensus = substr $line2, 16, 60;
							$consensus_total .= $consensus;
						}
						elsif ($read_id_tmp eq $total_nuc_count_original)
						{
							$query_line = "yes";
						} 
					}
					$e++
				}
				close INPUT_BLAST3;
                               
                my %gaps_align;
                undef %gaps_align;
                foreach my $subject_id (keys %subject_list)
                {
                    my $read_tmp = $subject_list{$subject_id};
                    my $o = '0';
                    my $first_nuc = substr $read_tmp, $o, 1;
                    while ($first_nuc eq "-" && $o < 120+($total_nuc_count_original*2))
                    {
                        $o++;
                        $first_nuc = substr $read_tmp, $o, 1;
                    }
                    if (exists($gaps_align{$o}))
                    {
                        my $count_tmp = $gaps_align{$o};
                        $gaps_align{$o} = $count_tmp+1;
                    }
                    else
                    {
                        $gaps_align{$o} = '1';
                    } 
                }
                my $count_align = keys %subject_list;
                my $oo = '1';
                my $cut_gap = '1';
                foreach my $gap_length (sort {$a <=> $b} keys %gaps_align)
                {
                    $oo += $gaps_align{$gap_length};
                    if ($oo > 0.45*$count_align)
                    {
                        $cut_gap = $gap_length;
                        last;
                    }                 
                }
                if ($cut_gap > $cp)
                {
                    $cp = $cut_gap;
                }
				if ($mafft_count > 2)
				{
					$cp = '0';
				}
				else
				{
					$cp_original = $cp;
				}
                print {$filehandle{$seed_id2}} $cp." SKIP_LENGTH\n";
                
                my $too_much_gaps = '0';
                foreach my $subject_id (keys %subject_list)
                {
                    my $read_tmp = $subject_list{$subject_id};
                    my $last_20 = substr $read_tmp, -20;
                    my $gaps_20 = $last_20 =~ tr/-/-/;
                    while ($gaps_20 > 15)
                    {
                        substr $read_tmp, -20, 20, "";
                        $last_20 = substr $read_tmp, -20;
                        $gaps_20 = $last_20 =~ tr/-/-/;
                    }
                    
                    my $last_nuc = substr $read_tmp, -1;
                    while ($last_nuc eq "-")
                    {
                        chop($read_tmp);
                        $last_nuc = substr $read_tmp, -1;
                    }
                    substr $read_tmp, -80, 80, "";
                    $subject_list{$subject_id} = $read_tmp;
                    
                    my $gaps_tmp = $read_tmp =~ tr/-/-/;
                    $gaps_id{$subject_id} = $gaps_tmp;
                    $length_id{$subject_id} = length($read_tmp)-5;
                    
                    #if ($gaps_tmp > 0.10*length($read_tmp))
                    #{
                        #$too_much_gaps++;
                    #}
                    #if ($too_much_gaps > 0.2*$matches_count && $mafft_extra_quality eq "")
                    #{
                        #$mafft_extra_quality = "yes";
                        #print {$filehandle{$seed_id2}} "TOO_MANY_GAPS_MAFFT_EXTRA_QUALITY\n";
                        #goto MAFFT_PB;
                    #}
                }
				
				unlink $mafft_output_tmp;
				
#check quality scores--------------------------------------------------------------------------------------------------
                my %quality_scores_reads;
				undef %quality_scores_reads;
                if ($use_quality_scores_NP ne "")
                {
					foreach my $subject_rank (keys %subject_list)
					{
						if (exists($rank_to_id{$subject_rank}))
						{
							my $id_tmpie5 = $rank_to_id{$subject_rank};
							if (exists($quality_scores_NP{$id_tmpie5}))
							{
								my @list_tmp = split /,/, $quality_scores_NP{$id_tmpie5};
								foreach my $pos_tmp (@list_tmp)
								{
									$quality_scores_reads{$subject_rank}{$pos_tmp} = undef;
								}
							}
						}
					}
				}	
#check quality scores--------------------------------------------------------------------------------------------------
				
				foreach my $subject_rank (keys %subject_list)
                {
                    $subject_list_original{$subject_rank} = $subject_list{$subject_rank};
                }
                my %ignore_reads;
                undef %ignore_reads;
IGNORE_REMOVED_READS_NP:

                my $ignored_reads_count = keys %ignore_reads;

                if (keys %ignore_reads > 0)
                {
                    undef %subject_list;
                    print {$filehandle{$seed_id2}} $ignored_reads_count." IGNORE_COUNT\n";
                    foreach my $subject_rank (keys %subject_list_original)
                    {
                        if (exists($ignore_reads{$subject_rank}))
                        {
                            print {$filehandle{$seed_id2}} $subject_rank." REMOVE_READ\n";
                            delete $subject_list_original{$subject_rank};
                        }
                        else
                        {
                            $subject_list{$subject_rank} = $subject_list_original{$subject_rank};
                        }
                    }
                    $cp = $cp_original;
					if ($mafft_count > 2)
					{
						$cp = '0';
					}
                }
                
				if ($mafft_count < 3)
				{
					$clipped_ext = "";
					$N = '0';
				}

                
				my $best_extension_part = "";
                my %nucs;
                undef %nucs;
                my %nucs_rej;
                undef %nucs_rej;

                my $total_nuc_count = '0';
                my $total_nuc_count_rej = '0';
				my $local_pattern_matches2 = '0';
                my $nuc1 = "";
				$unresolvable_NP = "";

                my $only_2_reads = "";
                if (keys %extensions2_tmp eq '2' && keys %subject_list eq '2')
                {
                    $only_2_reads = "yes";
                    print {$filehandle{$seed_id2}} "ONLY_2_READS\n";
                }
				if (keys %extensions2_tmp < 4  && $full_reset_NP eq "")
				{
					if ($hap_tag eq "HAP1")
					{
						undef %exclude_reads_hap1_NP;
					}
					elsif ($hap_tag eq "HAP2")
					{
						undef %exclude_reads_hap2_NP;
					}
					undef %save_alignment_data_NP;
					undef %rejected_alignment_data_NP;
					$full_reset_NP = "A";
					print {$filehandle{$seed_id2}} $best_extension." FULL_RESET\n";
					$best_extension = "";
					$seed_id = $id;
					$y++;
					$y{$id} = $y;
					goto FULL_RESET;
				}
                my %track_length_ext;
                undef %track_length_ext;
				my $end_this_mafft_part = "";
					
				if ($mafft_count < 3)
				{
					foreach my $subject_rank (sort {$a <=> $b} keys %subject_list)
					{				
						my $seq_rank1c = substr $subject_list{$subject_rank}, 0, $cp;
						$seq_rank1c =~ tr/-//d;
						$track_length_ext{$subject_rank} = length($seq_rank1c);
					}
				}
                
                my $time_BLAST3 = time;
                my $time9 = $time_BLAST3 - $time_maff3;
                print {$filehandle{$seed_id2}} $time9." TIME_pre_CONS\n\n";
				my $loop_check = "";
				my $NP_reads_support_SNR2 = "";
                
INPUT_MAFFT3_NP: while ((keys %subject_list > 2 || $only_2_reads eq "yes") && ((keys %subject_list > 0.75*($total_nuc_count_original-$ignored_reads_count) && (keys %subject_list > 10
						|| keys %subject_list > $sequencing_depth_NP*1.4)) || keys %subject_list > ($total_nuc_count_original-1-$ignored_reads_count) || $clipped_ext ne "yes" 
                        || (keys %subject_list > 0.3*($total_nuc_count_original-$ignored_reads_count) && $longer_extension_for_repeat ne "" && keys %subject_list > 3))
                        && ($end_this_mafft_part eq "" || length($best_extension_part) < $length_extension_part-50))
                {
					$loop_check = "yes";
#clip_extension-------------------------------------------------------------------------------------------                                
                    if ($clipped_ext eq "")
					{
						$clipped_ext_pos{length($best_extension)} = $cp;
					}
					if (length($best_extension) > 14 && $clipped_ext eq "" && $find_haps_in_seed ne "yes" && $NP_reads_support_SNR2 eq "")
                    {
						my $check_start_assembly = substr $read, -length($best_extension), length($best_extension);
                        my $check_start_assembly2 = substr $read, -150;
						if ($NP_reads_support_SNR ne "")
						{
							$check_start_assembly2 = substr $read, -150, 146;
							$check_start_assembly = substr $read, -length($best_extension)-4, length($best_extension);
						}
                        my $check_start_assembly3 = $check_start_assembly2 =~ m/$check_start_assembly/;
                        if ($check_start_assembly3 > 1)
                        {
                            print {$filehandle{$seed_id2}} $check_start_assembly3." CLIP_ALERT\n";
                            goto SWITCH_CLIP_NP;
                        }
						if ($check_start_assembly eq $best_extension && $NP_reads_support_SNR ne "")
                        {
                            $NP_reads_support_SNR2 = $NP_reads_support_SNR;
							$NP_reads_support_SNR2 =~ tr/ACTGN/actgn/;
                        }
                        elsif ($check_start_assembly eq $best_extension)
                        {
                            $clipped_ext = "yes";
                        }
                        else
                        {
                            my $best_extension_tmp = $best_extension;
                            my $check_N = $best_extension_tmp =~ tr/N|\./\./;
                            my $check_N2 = $check_start_assembly =~ tr/N|\./\./;
                            if (($check_N > 0 || $check_N2 > 0) && $check_N < length($best_extension_tmp)*0.4 && $check_N2 < length($check_start_assembly)*0.4)
                            {
                                my $check_again = $check_start_assembly =~ s/$best_extension_tmp//;
                                my $check_again2 = $best_extension_tmp =~ s/$check_start_assembly//;
                                if ($check_again > 0 || $check_again2 > 0)
                                {
									if ($NP_reads_support_SNR ne "")
									{
										$NP_reads_support_SNR2 = $NP_reads_support_SNR;
										$NP_reads_support_SNR2 =~ tr/ACTGN/actgn/;
									}
									else
									{
										$clipped_ext = "yes";
									}
                                }
                            }
                        }
                        if ($clipped_ext eq "yes")
                        {
                            print {$filehandle{$seed_id2}} $best_extension." BEST_EXT_EXACT_MATCH\n";
							$best_extension = "";
                            $N = '0';
                            undef %SNP_patterns_prev;
                            undef %extensions_seed;
                            undef %quality_scores_tmp;
                            undef %quality_scores_gap_tmp;
                            undef %extensions_seed;
                            undef %best_read_score;
							undef %total_count_prev_patterns;
                            $total_count_prev_patterns = '0';
                        }
                    }
                    if (length($best_extension) > 120 && $clipped_ext eq "" && $find_haps_in_seed ne "yes" && $NP_reads_support_SNR2 eq "")
                    {
                        my $m = '15';
						my $best_extension_tmp_length = length($best_extension);
                        while ($m < length($best_extension)-10 && $clipped_ext eq "")
                        {
                            my $check_start_assembly = substr $read, -$m, 15;
                            $check_start_assembly =~ tr/N|\./\./;
                            my $best_extension_tmp = $best_extension;
							
                            my $check_15 = $best_extension_tmp =~ s/$check_start_assembly//;
                            if ($check_15 eq 1)
                            {
                                print {$filehandle{$seed_id2}} $best_extension." BEST_EXT_CLIP\n";
                                if ($best_extension =~ m/(.*$check_start_assembly)(.*)/)
                                {
                                    $best_extension = $2;
                                    my $best_extension_tmp = substr $best_extension, 0, $m-15+3;
                                    my $best_extension_tmp2 = substr $best_extension, $m-15+3;
                                    my $m2 = '14';
                                    my $matched = "";
                                    while ($m2 > 5 && $matched eq "" && $m > 15)
                                    {
                                        my $check_start_assembly2 = substr $read, -$m2, $m2;
										$check_start_assembly2 =~ tr/N|\./\./;
                                        if ($best_extension_tmp =~ m/(.*$check_start_assembly2)(.*)/)
                                        {
                                            $best_extension = $2;
                                            $best_extension .= $best_extension_tmp2;
                                            $matched = "yes";
                                            last
                                        }
										else
										{
											my $m3 = '0';
											while ($m3 < 8)
											{
												my $check_start_assembly3 = substr $best_extension_tmp, -$m2-$m3, $m2;
												$check_start_assembly3 =~ tr/N|\./\./;
												if ($check_start_assembly2 =~ m/$check_start_assembly3/)
												{
													$best_extension = "";
													if ($m3 > 0)
													{
														$best_extension = substr $best_extension_tmp, -$m3;
													}
													$best_extension .= $best_extension_tmp2;
													$matched = "yes";
													last
												}
												$m3++;
											}
										}
                                        $m2--;
                                    }
                                    if ($matched ne "yes" && $m > 15)
                                    {
                                        substr $best_extension, 0, $m-15, "";
										print {$filehandle{$seed_id2}} $m2." NO_EXACT_MATCH\n";
                                    }
                                    #my $M_seq = substr $best_extension, 0, $m-15, "";
                                    #$ext_remove = $1.$M_seq;
                                    #print {$filehandle{$seed_id2}} $ext_remove." EXT_REMOVE\n";   
                                }

                                print {$filehandle{$seed_id2}} $m." M\n";
                                print {$filehandle{$seed_id2}} $best_extension." BEST_EXT_CHOPPED\n";
								my $diff_tmp = $best_extension_tmp_length-length($best_extension);
                                $clipped_ext = "yes";
								$cp = $clipped_ext_pos{$diff_tmp};
								foreach my $subject_rank (sort {$a <=> $b} keys %subject_list)
								{				
									my $seq_rank1c = substr $subject_list{$subject_rank}, 0, $cp;
									$seq_rank1c =~ tr/-//d;
									$track_length_ext{$subject_rank} = length($seq_rank1c);
								}
                            }
                            $m++;
                        }
                        if ($clipped_ext eq "yes")
                        {
                            $best_extension = "";
                            $N = '0';
                            undef %SNP_patterns_prev;
                            undef %extensions_seed;
                            undef %quality_scores_tmp;
                            undef %quality_scores_gap_tmp;
                            undef %extensions_seed;
                            undef %best_read_score;
							undef %total_count_prev_patterns;
                            $total_count_prev_patterns = '0';
							if ($find_haps_in_seed ne "")
							{
								$extensions_seed{'HAP1'} = $best_extension;
								$extensions_seed{'HAP2'} = $best_extension;
							}
                        }
                    }				
SWITCH_CLIP_NP:                     
#-------------------------------------------------------------------------------------------------------
                    undef %nucs;
                    undef %nucs_rej;
                    $total_nuc_count = '0';
                    $total_nuc_count_rej = '0';
                    $nuc1 = "";
                    my %track_mismatch_ext0;
                    undef %track_mismatch_ext0;
                    $track_mismatch_count++;
					$nucs{'a'} = '0';
					$nucs{'c'} = '0';
					$nucs{'t'} = '0';
					$nucs{'g'} = '0';
					$nucs{"-"} = '0';
					my $nuc_top10 = "";
					my %nucs_by_rank;
					undef %nucs_by_rank;
					my %split_patterns_final_score;
                    undef %split_patterns_final_score;
					my $lowest_longest_match = "";

                    foreach my $subject_rank (sort {$a <=> $b} keys %subject_list)
                    {
                        my $nuc = substr $subject_list{$subject_rank}, $cp, 1;
						$nucs_by_rank{$subject_rank} = $nuc;
                        
                        if (($add_rejected_reads ne "" && $subject_rank > $add_rejected_reads) || ($add_no_match_reads ne "" && $subject_rank > $add_no_match_reads))
                        {
                            $nucs_rej{$nuc} += 1;
                            $total_nuc_count_rej++;
							
							if (exists($track_mismatch_ext0{$nuc}))
							{
								$track_mismatch_ext0{$nuc} .= ",".$subject_rank;
							}
							else
							{
								$track_mismatch_ext0{$nuc} = $subject_rank;
							}
                        }
                        else
                        {
                            $nucs{$nuc} += 1;
                            $total_nuc_count++;
							if ($nuc_top10 eq "")
							{
								$nuc_top10 = $nuc;
							}
							elsif ($subject_rank < 11 && $nuc_top10 ne $nuc)
							{
								$nuc_top10 = "no";
							}
                        }
                        
                        if ($nuc ne "-")
                        {
                            $track_length_ext{$subject_rank} += 1;
                        }
                        
                        if ($subject_rank eq "1")
                        {
							$nuc1 = $nuc;
                            my $length_rank1_tmp = $track_length_ext_total{$subject_rank}+$track_length_ext{$subject_rank};
                            my $posiie = length($best_extension)+$position;
                            if (exists($split_positions_DUP_tmp{$length_rank1_tmp}))
                            {
                                print {$filehandle{$seed_id2}} $seed_id."\t".$posiie."\t".$length_rank1_tmp."\t".$split_positions_DUP_tmp{$length_rank1_tmp}." DUP_POS_TMP\n";
                                $split_positions_DUP{$seed_id}{$posiie} = $split_positions_DUP_tmp{$length_rank1_tmp};
								my $read_end_tmpi = substr $best_extension, -$overlap;
								my @split_positions_DUP_tmp = split /,/, $split_positions_DUP_tmp{$length_rank1_tmp};
								my $selected_nuc_tmp = $split_positions_DUP_tmp[1];
								$split_positions_DUP{$seed_id}{$posiie} = $read_end_tmpi.",".$selected_nuc_tmp;
                            }
                        }
                        elsif ($nuc1 ne "" && (($add_rejected_reads eq "" && $add_no_match_reads eq "") || ($subject_rank < $add_rejected_reads && $subject_rank < $add_no_match_reads)))
                        {
                            if ($nuc eq $nuc1)
                            {
                                if (exists($best_read_score{$subject_rank}))
                                {
                                    my $score_tmp = $best_read_score{$subject_rank}+2;
                                    $best_read_score{$subject_rank} = $score_tmp;
                                }
                                else
                                {
                                    $best_read_score{$subject_rank} = '2';
                                }  
                            }
                            else
                            {
                                if (exists($best_read_score{$subject_rank}))
                                {
                                    my $score_tmp = $best_read_score{$subject_rank}-1;
                                    $best_read_score{$subject_rank} = $score_tmp;
                                }
                                else
                                {
                                    $best_read_score{$subject_rank} = '-1';
                                }
                            }
                        }                        
                    }
						
					if ($NP_reads_support_SNR2 ne "" && $nucs{$NP_reads_support_SNR2}+$nucs{"-"} < 0.8*$total_nuc_count)
					{
						$best_extension = "";
						$N = '0';
						undef %SNP_patterns_prev;
						undef %extensions_seed;
						undef %quality_scores_tmp;
						undef %quality_scores_gap_tmp;
						undef %extensions_seed;
						undef %best_read_score;
						undef %total_count_prev_patterns;
						$total_count_prev_patterns = '0';
						$NP_reads_support_SNR2 = "";
						$clipped_ext = "yes";
					}
					
#Remove rejected reads that don't align well----------------------------------------------------------------------------------------                    
                    if ($add_rejected_reads ne "" || $add_no_match_reads ne "")
                    #if ($add_rejected_reads eq "sghsh")
                    {
                        my $removed_reads = '0';
                        foreach my $nuc_5 (sort {$a <=> $b} keys %track_mismatch_ext0)
                        {
                            if ($nucs{$nuc_5} < $total_nuc_count*0.1)
                            {
                                my @track_mismatch_ext0 = split /,/, $track_mismatch_ext0{$nuc_5};
                                foreach my $rank_5 (@track_mismatch_ext0)
                                {
                                    #$track_mismatch_ext{$rank_5}{$track_mismatch_count} += 1;
                                    push @{ $track_mismatch_ext{$rank_5} }, $track_mismatch_count;
                                    my @rank_array = @{ $track_mismatch_ext{$rank_5} };
                                    my $array_count_tmp = @rank_array;
                                    my $threshold_tmp = 0.35*350;
                                    my $threshold_tmp_round = sprintf "%.0f", $threshold_tmp;
                                    
                                    if ($array_count_tmp > $threshold_tmp_round && $rank_array[$array_count_tmp-$threshold_tmp_round] > $track_mismatch_count-350)
                                    {
                                        foreach my $rank_6 (keys %track_mismatch_ext)
                                        {
                                            my @rank_array2 = @{ $track_mismatch_ext{$rank_6} };
                                            my $array_count_tmp2 = @rank_array2;
                                            my $threshold_tmp2 = 0.35*320;
                                            my $threshold_tmp_round2 = sprintf "%.0f", $threshold_tmp2;
                                            
                                            if ($array_count_tmp2 > $threshold_tmp_round2 && $rank_array2[$array_count_tmp2-$threshold_tmp_round2] > $track_mismatch_count-320)
                                            {
                                                $removed_reads++;
                                                $ignore_reads{$rank_6} = undef;
                                                my $id_tmp0 = $rank_to_id{$rank_6};
                                                delete $extensions2_tmp{$id_tmp0};
                                                delete $extensions{$extensions2{$id_tmp0}};
                                                delete $extensions2{$id_tmp0};
                                                delete $extensions_nomatch{$extensions2{$id_tmp0}};
                                                delete $extensions_nomatch2{$id_tmp0};
                                                delete $extensions_unknown{$extensions2{$id_tmp0}};
                                                delete $extensions_unknown2{$id_tmp0};
                                                delete $save_reads_for_next{$id_tmp0};
                                                delete $add_rej_reads_extra{$id_tmp0};
                                                delete $read_start_pos_rej{$id_tmp0};
												delete $rejected_reads_save{$id_tmp0};
                                                $ext2_count = keys %extensions2_tmp;              
                                            }
                                        }
                                    }
                                }
                            }                           
                        }
                        if ($removed_reads > 0)
                        {
                            $best_extension = "";
							$best_extension_part = "";
                            $mismatch_retry++;
                            undef %quality_scores_tmp;
                            my $ignore_count_tmp = keys %ignore_reads;
                            
                            foreach my $rank_7 (keys %ignore_reads)
                            {
                                my $id_tmp0 = $rank_to_id{$rank_7};
                                delete $extensions2_tmp{$id_tmp0};
                                delete $extensions{$extensions2{$id_tmp0}};
                                delete $extensions2{$id_tmp0};
                                delete $extensions_nomatch{$extensions2{$id_tmp0}};
                                delete $extensions_nomatch2{$id_tmp0};
                                delete $extensions_unknown{$extensions2{$id_tmp0}};
                                delete $extensions_unknown2{$id_tmp0};
                                delete $save_reads_for_next{$id_tmp0};
                                delete $add_rej_reads_extra{$id_tmp0};
                                delete $read_start_pos_rej{$id_tmp0};
								delete $rejected_reads_save{$id_tmp0};
                                $ext2_count = keys %extensions2_tmp;
                            }
                            if ($removed_reads eq '1000000000000000000000000000000' && $ignore_count_tmp eq '1')
                            {
                                print {$filehandle{$seed_id2}} "REMOVE_BAD_ALIGNMENTS_REJ\n";                     
                                goto IGNORE_REMOVED_READS_NP;
                            }
                            print {$filehandle{$seed_id2}} $removed_reads." REMOVE_BAD_ALIGNMENTS_REJ2\n";
                            goto MISMATCH_RETRY_NP;
                        }
                    }
#Check SNR ahead-------------------------------------------------------------------------------------------------------                   
                    $SNR_read_ahead = "";
					if ($nucs{"-"} > $total_nuc_count*0.15 && $nucs{"-"} < $total_nuc_count*0.9)
                    {
                        my $count_SNR = '0';
                        foreach my $subject_rank (sort {$a <=> $b} keys %subject_list)
                        {
                            my $p = '1';
                            my $p2 = '0';
							my %nucs_SNR;
							undef %nucs_SNR;
                            while ($p2 < 6 && $p < 30)
                            {
                                my $nuc2 = substr $subject_list{$subject_rank}, $cp+$p, 1;
                                if ($nuc2 ne "-")
                                {
                                    $nucs_SNR{$nuc2} += 1;
									$p2++;
                                }
                                $p++;
                            }
							foreach my $nuc_SNR_tmp (keys %nucs_SNR)
							{
								if ($nucs_SNR{$nuc_SNR_tmp} > 4)
								{
									$count_SNR++;
								}
							}
                        }
                        if ($count_SNR > $total_nuc_count*0.7)
                        {
                            $SNR_read_ahead = "yes";
                        }
                    }
                    if ($clipped_ext eq "yes")
                    {
                        #print {$filehandle{$seed_id2}} $nucs{"a"}." A ".$nucs{"c"}." C ".$nucs{"t"}." T ".$nucs{"g"}." G ".$nucs{"-"}." GAP\n";
                    }
#End NP assembly and return to PB---------------------------------------------------------------------------------------------------------------------------------------                                    
                    if ($NP_reads_support ne "" && $clipped_ext eq "yes" && length($best_extension) > $nuc_unique_for_ONT && (($nucs{"a"} < 0.8*$total_nuc_count &&
                      $nucs{"c"} < 0.8*$total_nuc_count && $nucs{"t"} < 0.8*$total_nuc_count && $nucs{"g"} < 0.8*$total_nuc_count && $nucs{"-"} < 0.8*$total_nuc_count) || $NP_reads_support eq "yes2"))
                    {
                        my $w = '1';
                        my $low_score = "";
                        if ($repetitive_detect1 eq "yes")
                        {
                            $w = '20';
                        }
                        foreach my $pos_tmp (sort {$a <=> $b} keys %quality_scores_tmp)
                        {
                            my @q_score_tmp = split / /, $quality_scores_tmp{$pos_tmp};
							if ($q_score_tmp[0] < 0.8 && $w > $nuc_unique_for_ONT)
                            {
                                $low_score = "yes";
                                last;
                            }
                            $w++;
                        }
                        if ($low_score eq "yes" || length($best_extension) <= 20)
                        {
                            my $best_extension_tmp = substr $best_extension, 0, $w;
                            $best_extension = $best_extension_tmp;
                            $NP_reads_support = "yes";
							$NP_reads_support_SNR = "";
                            print {$filehandle{$seed_id2}} $best_extension." LOW_SCORE_STOP\n";
                            my $one_match = "";
                            foreach my $nuc_tmp (keys %PB_split_nucs)
                            {
                                foreach my $seq_tmp (keys %{$PB_split_nucs{$nuc_tmp}})
								{
									print {$filehandle{$seed_id2}} $seq_tmp." SEQ_TMP\n";
									my $NP_seq = substr $best_extension, 0, length($seq_tmp);
									my $N_check = $NP_seq =~ tr/N/N/;
									if ($seq_tmp eq $NP_seq)
									{
										$one_match = "yes";
									}
									elsif ($NP_seq ne "N" && ($N_check < 1 || $N_check eq ""))
									{
										my @ids_PB = split /,/, $PB_split_nucs{$nuc_tmp}{$seq_tmp};
										print {$filehandle{$seed_id2}} $PB_split_nucs{$nuc_tmp}{$seq_tmp}." EXCLUDE_PB\n";
										foreach my $id_PB (@ids_PB)
										{
											if ($hap_tag eq "HAP1")
											{
												$exclude_reads_hap1_PB{$id_PB} = $position+10000;
											}
											elsif ($hap_tag eq "HAP2")
											{
												$exclude_reads_hap2_PB{$id_PB} = $position+10000;
											}
											delete $PB_split_ids{$id_PB};
										}
									}
                                }
                            }
                            if ($one_match ne "yes")
                            {
                                if ($hap_tag eq "HAP1")
                                {
                                    undef %exclude_reads_hap1_PB;
                                }
                                elsif ($hap_tag eq "HAP2")
                                {
                                    undef %exclude_reads_hap2_PB;
                                }
                            }
							else
							{
								$best_extension = "";
								substr $read, -length($PB_extension), length($PB_extension), "";                                                 
								$position -= length($PB_extension);
								$position{$id} = $position;
								$seed{$id} = $read;
								goto PB_READS;
							}
                            goto END_NP;
                        } 
                    }
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------                                    
         
                    if ($nucs{"-"} < $total_nuc_count*0.2)
                    {
                        #$total_nuc_count -= $nucs{"-"};
                    }
					my $q_check = "";
ADD_QUALITY_NP:                  
                    if (($nucs{"a"} > $total_nuc_count*0.76 || ($nucs{"a"} > $total_nuc_count*0.49 && $PB_reads eq "" && $input_reads_DB_folder_PB eq "" && $only_2_reads eq "yes" && ($nuc1 eq 'a' || $nuc1 eq "-")))
                       && ($nucs{"a"} > $total_nuc_count*0.83 || $total_nuc_count < 18 || $count_matches_with_high_scores < $sequencing_depth_NP*2.4 || $nuc_top10 ne "no"))
                    {
                        $best_extension .= "A";
						$best_extension_part .= "A";
                        if ($total_nuc_count > 2)
                        {
                            $quality_scores_tmp{length($best_extension)} = ($nucs{"a"}+($nucs{"-"}/2))/$total_nuc_count." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
                        }
                        #print {$filehandle{$seed_id2}} length($best_extension)." POS_QUAL ".$quality_scores_tmp{length($best_extension)}." QUAL ".$nucs{"a"}." A ".$nucs{"c"}." C ".$nucs{"t"}." T ".$nucs{"g"}." G ".$nucs{"-"}." GAP\n";
                        $nuc_match = "A";
                        goto SKIP_INPUT_BLAST3_NP;
                    }
                    elsif (($nucs{"c"} > $total_nuc_count*0.76 || ($nucs{"c"} > $total_nuc_count*0.49 && $PB_reads eq ""  && $input_reads_DB_folder_PB eq "" && $only_2_reads eq "yes" && ($nuc1 eq 'c' || $nuc1 eq "-")))
                        && ($nucs{"c"} > $total_nuc_count*0.83 || $total_nuc_count < 18 || $count_matches_with_high_scores < $sequencing_depth_NP*2.4 || $nuc_top10 ne "no"))
                    {
                        $best_extension .= "C";
						$best_extension_part .= "C";
                        if ($total_nuc_count > 2)
                        {
                            $quality_scores_tmp{length($best_extension)} = ($nucs{"c"}+($nucs{"-"}/2))/$total_nuc_count." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
                        }
                        #print {$filehandle{$seed_id2}} length($best_extension)." POS_QUAL ".$quality_scores_tmp{length($best_extension)}." QUAL ".$nucs{"a"}." A ".$nucs{"c"}." C ".$nucs{"t"}." T ".$nucs{"g"}." G ".$nucs{"-"}." GAP\n";
                        $nuc_match = "C";
                        goto SKIP_INPUT_BLAST3_NP;
                    }
                    elsif (($nucs{"t"} > $total_nuc_count*0.76 || ($nucs{"t"} > $total_nuc_count*0.49 && $PB_reads eq ""  && $input_reads_DB_folder_PB eq "" && $only_2_reads eq "yes" && ($nuc1 eq 't' || $nuc1 eq "-")))
                          && ($nucs{"t"} > $total_nuc_count*0.83 || $total_nuc_count < 18 || $count_matches_with_high_scores < $sequencing_depth_NP*2.4 || $nuc_top10 ne "no"))
                    {
                        $best_extension .= "T";
						$best_extension_part .= "T";
                        if ($total_nuc_count > 2)
                        {
                            $quality_scores_tmp{length($best_extension)} = ($nucs{"t"}+($nucs{"-"}/2))/$total_nuc_count." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
                        }
                        #print {$filehandle{$seed_id2}} length($best_extension)." POS_QUAL ".$quality_scores_tmp{length($best_extension)}." QUAL ".$nucs{"a"}." A ".$nucs{"c"}." C ".$nucs{"t"}." T ".$nucs{"g"}." G ".$nucs{"-"}." GAP\n";
                        $nuc_match = "T";
                        goto SKIP_INPUT_BLAST3_NP;
                    }
                    elsif (($nucs{"g"} > $total_nuc_count*0.76 || ($nucs{"g"} > $total_nuc_count*0.49 && $PB_reads eq "" && $input_reads_DB_folder_PB eq "" && $only_2_reads eq "yes" && ($nuc1 eq 'g' || $nuc1 eq "-")))
                        && ($nucs{"g"} > $total_nuc_count*0.83 || $total_nuc_count < 18 || $count_matches_with_high_scores < $sequencing_depth_NP*2.4 || $nuc_top10 ne "no"))
                    {
                        $best_extension .= "G";
						$best_extension_part .= "G";
                        if ($total_nuc_count > 2)
                        {
                            $quality_scores_tmp{length($best_extension)} = ($nucs{"g"}+($nucs{"-"}/2))/$total_nuc_count." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
                        }
                        #print {$filehandle{$seed_id2}} length($best_extension)." POS_QUAL ".$quality_scores_tmp{length($best_extension)}." QUAL ".$nucs{"a"}." A ".$nucs{"c"}." C ".$nucs{"t"}." T ".$nucs{"g"}." G ".$nucs{"-"}." GAP\n";
                        $nuc_match = "G";
                        goto SKIP_INPUT_BLAST3_NP;
                    }
                    elsif (($nucs{"-"} > $total_nuc_count*0.76 || ($nucs{"-"} > $total_nuc_count*0.49 && $nucs{"c"} < $total_nuc_count*0.25 && $nucs{"g"} < $total_nuc_count*0.25 &&
                         $nucs{"t"} < $total_nuc_count*0.25 && $nucs{"a"} < $total_nuc_count*0.25 && $PB_reads eq "" && $input_reads_DB_folder_PB eq "")) &&
						   ($nucs{"-"} > $total_nuc_count*0.83 || $total_nuc_count < 18 || $count_matches_with_high_scores < $sequencing_depth_NP*2.4 || $nuc_top10 ne "no"))
                    {
                        if ($nucs{"-"}/$total_nuc_count < 0.8)
						{
							$quality_scores_gap_tmp{length($best_extension)} = $nucs{"-"}/$total_nuc_count." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
						}
						
                        $nuc_match = "";
                        goto SKIP_INPUT_BLAST3_NP;
                    }
					elsif ($use_quality_scores_NP ne "" && $add_rejected_reads eq "" && $add_no_match_reads eq "" && $q_check eq "")
					{
						my %nucs_original;
						undef %nucs_original;
						foreach my $nucs_tmpi (keys %nucs)
						{
							$nucs_original{$nucs_tmpi} = $nucs{$nucs_tmpi};
						}
						undef %nucs;
						$total_nuc_count = '0';
						
						foreach my $rank_tmp (sort {$a <=> $b} keys %subject_list)
						{
							my $length_rank_tmp = $track_length_ext{$rank_tmp};
							if (exists($rank_to_id{$rank_tmp}))
							{
								my $id_tmp8 = $rank_to_id{$rank_tmp};
								if (exists($long_read_end_pos_save{$id_tmp8}))
								{
									my $total_length_rank_tmp = $track_length_ext{$rank_tmp}+$long_read_end_pos_save{$id_tmp8}-90;
									my $firstLine = $hash_NP_reads_tmp{$id_tmp8};
									if (exists($reverse_list{$id_tmp8}))
									{
										$total_length_rank_tmp = length($firstLine)-($track_length_ext{$rank_tmp}+$long_read_end_pos_save{$id_tmp8}-90);
									}
									my $nuc_tmp = substr $firstLine, $total_length_rank_tmp-10, 10;
									if (exists($reverse_list{$id_tmp8}))
									{
										$nuc_tmp = substr $firstLine, $total_length_rank_tmp, 10;
										my $nuc_tmp2 = reverse($nuc_tmp);
										$nuc_tmp2 =~ tr/ACTG/TGAC/;
										$nuc_tmp = $nuc_tmp2;
										#print {$filehandle{$seed_id2}} $rank_tmp." ID ".$nuc_tmp." NUC_REV\n";
									}
									else
									{
										#print {$filehandle{$seed_id2}} $rank_tmp." ID ".$nuc_tmp." NUC\n";
									}
									$total_length_rank_tmp++;
									
									if (exists($quality_scores_reads{$rank_tmp}{$total_length_rank_tmp}) || exists($quality_scores_reads{$rank_tmp}{$total_length_rank_tmp-1}) || exists($quality_scores_reads{$rank_tmp}{$total_length_rank_tmp+1}))
									{
										#print {$filehandle{$seed_id2}} $rank_tmp." ID ".$total_length_rank_tmp." POS\n";
										$q_check = "yes";
									}
									else
									{
										my $nuc_tmp2 = substr $subject_list{$rank_tmp}, $cp, 1;
										$nucs{$nuc_tmp2} += 1;
										$total_nuc_count++;
									}
								}
							}
						}
						my $check_error_pattern = "";
						foreach my $nucs_tmpi (keys %nucs)
						{
							if (exists($nucs_original{$nucs_tmpi}))
							{
								if ($nucs_original{$nucs_tmpi} > $total_nuc_count*0.3 && $nucs{$nucs_tmpi} > $nucs_original{$nucs_tmpi}*0.69)
								{
									$check_error_pattern = "yes";
								}
							}
						}
						if ($total_nuc_count < 5 || $check_error_pattern eq "")
						{
							undef %nucs;
							$total_nuc_count = '0';
							$q_check = "yes";
							
							foreach my $nucs_tmpi (keys %nucs_original)
							{
								$nucs{$nucs_tmpi} = $nucs_original{$nucs_tmpi};
								$total_nuc_count += $nucs_original{$nucs_tmpi};
							}
							goto ADD_QUALITY_NP;
						}
						if ($q_check eq "yes")
						{
							goto ADD_QUALITY_NP;
						}
					}
                    
#Resolve ambigious positions through other haplotype----------------------------------------------------------------------------------------------                                   
                    if ($find_haps_in_seed eq "" && $ext2_count > 0 && $merged_ext ne "")
                    {
                        my $nuc_highest_tmp = "";
                        my $count_highest_tmp = '0';
                        my $count_total_tmp = '0';
                        my %nuc_tmp2;
                        undef %nuc_tmp2;

                        foreach my $ranki_tmp (sort {$a <=> $b} keys %subject_list)
                        {
                            if (exists($rank_to_id{$ranki_tmp}))
                            {
                                if (exists($extensions2{$rank_to_id{$ranki_tmp}}))
                                {
                                    my $nuc_tmp = $nucs_by_rank{$ranki_tmp};
                                    $nuc_tmp2{$nuc_tmp} += 1;
                                    if ($nuc_tmp2{$nuc_tmp} > $count_highest_tmp)
                                    {
                                        $count_highest_tmp = $nuc_tmp2{$nuc_tmp};
                                        $nuc_highest_tmp = $nuc_tmp;
                                    }
                                    $count_total_tmp++;
                                }
                            }
                        }

                        if (($nuc_highest_tmp ne "" && $count_highest_tmp > $total_nuc_count*0.3) && $count_highest_tmp > 4 && $count_highest_tmp > $count_total_tmp-2 && $nuc_tmp2{$nuc_highest_tmp} > $count_total_tmp*0.8)
                        {
                            $nuc_match = $nuc_highest_tmp;
                            $nuc_match =~ tr/actgn/ACTGN/;
                            if ($nuc_highest_tmp eq "-")
                            {
                                $nuc_match = "";
                            }
                            $best_extension .= $nuc_match;
							$best_extension_part .= $nuc_match;
                            
                            if ($nuc_highest_tmp eq "-" && $nucs{"-"} > 0)
                            {
								if ($nucs{"-"}/$count_total_tmp < 0.8)
								{
									$quality_scores_gap_tmp{length($best_extension)} = $nucs{"-"}/$count_total_tmp." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
								}
                            }
                            else
                            {
                                $quality_scores_tmp{length($best_extension)} = $count_highest_tmp/$count_total_tmp." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
                            }           
                            
                            print {$filehandle{$seed_id2}} $nuc_match." N_CORRECTION\n";
							$total_count_prev_patterns{length($best_extension)} = undef;
                            goto SKIP_INPUT_BLAST3_NP
                        }
                    }     

                    
                    my $SNR_skip_check = "";
                    if ($nucs{"a"}+$nucs{"-"} < $total_nuc_count*0.9 && $nucs{"c"}+$nucs{"-"} < $total_nuc_count*0.9 && $nucs{"t"}+$nucs{"-"} < $total_nuc_count*0.9 && $nucs{"g"}+$nucs{"-"} < $total_nuc_count*0.9)
                    {
                        $SNR_skip_check = "yes";
                    }
					my $trace_back_check = "";
					my $find_haps_in_seed_check = "";
					if ($find_haps_in_seed ne "" && ($nucs{"a"}+$nucs{"-"} > $total_nuc_count*0.7 || $nucs{"c"}+$nucs{"-"} > $total_nuc_count*0.7 || $nucs{"t"}+$nucs{"-"} > $total_nuc_count*0.7
						|| $nucs{"g"}+$nucs{"-"} > $total_nuc_count*0.7) && ($nucs{"-"} < $total_nuc_count*0.4 || $nucs{"-"} > $total_nuc_count*0.75))
					{
						$find_haps_in_seed_check = "no";
					}
                                               
                    if (($clipped_ext ne "yes" || $best_extension eq "" || length($best_extension) < 1000 || $N < 5 || $N < length($best_extension)*0.08 ||
						 ($longer_extension_for_repeat ne "" && $N < length($best_extension)*0.1)) && $total_nuc_count > 5 && ($SNR_read_ahead eq "" || $SNR_skip_check eq "yes")
						 && $find_haps_in_seed_check eq "" && $found_haps_in_seed eq "")
                    {
#Split into haplotype groups----------------------------------------------------------------------                                                           
						my $time_split = time;		
                        my %separate_haps_NP;
                        undef %separate_haps_NP;

print {$filehandle{$seed_id2}} "\n".length($best_extension)." SEP_HAP\n";
print {$filehandle{$seed_id2}} $nucs{"a"}." A ".$nucs{"c"}." C ".$nucs{"t"}." T ".$nucs{"g"}." G ".$nucs{"-"}." GAP\n";                                    
                        my %SNP_patterns_now;
                        undef %SNP_patterns_now;
                        my %SNP_patterns_prev2;
                        undef %SNP_patterns_prev2;
						my %SNP_patterns_prev_match;
                        undef %SNP_patterns_prev_match;
                        my $current_pos_tmp = $position+length($best_extension);
                        my $selected_nuc = "";
						$total_count_prev_patterns{length($best_extension)} = undef;

                        my $first_rank = "";
                        foreach my $rank_tmp (sort {$a <=> $b} keys %subject_list)
                        {
                            my $nuc_tmp = $nucs_by_rank{$rank_tmp};
                            $SNP_patterns_prev{$current_pos_tmp}{$nuc_tmp}{$rank_tmp} = undef;
                            $SNP_patterns_now{$nuc_tmp}{$rank_tmp} = undef;
                            if ($first_rank eq "")
                            {
                                $first_rank = $nuc_tmp;
                            }
                            #print {$filehandle{$seed_id2}} $rank_tmp." ".$nuc_tmp."\n";
                            
                            if (exists($separate_haps_NP{$nuc_tmp}))
                            {                  
                                my $c_tmp = $separate_haps_NP{$nuc_tmp};
                                $separate_haps_NP{$nuc_tmp} = $c_tmp+1;      
                            }
                            else
                            {
                                $separate_haps_NP{$nuc_tmp} = 1;                          
                            }  
                        }
                        my $SNP_check = "";
                        if ($nucs{"-"} < $total_nuc_count*0.18 || ($NP_reads_support eq "yes2" && $nucs{"-"} < $total_nuc_count*0.35))
                        {
                            $SNP_check = "yes";
                            my $count_nucs = '0';
                            foreach my $nuc_tmp1 (keys %nucs)
                            {
                                if ($nuc_tmp1 ne "-" && $nucs{$nuc_tmp1} > 3 && $nucs{$nuc_tmp1} > $total_nuc_count*0.3)
                                {
                                    $count_nucs++;
                                }
                            }
                            if ($count_nucs > 1)
                            {
                                $SNP_check = "yes2";
                            }
                        }

#Compare previous split positions of this extension-----------------------------------------------------------------------------------------------------------------------                        
     
                        my %post_SNP_patterns;
                        undef %post_SNP_patterns;
						my %find_haps_SNPs;
                        undef %find_haps_SNPs;
                        my $prev_pos_count = '0';
						my $local_pattern_matches = '0';
						$local_pattern_matches2 = '0';
                        
POST_SNP_PATTERNS_TMP_NP: foreach my $pos_tmp (sort {$a <=> $b} keys %SNP_patterns_prev)
                        {
                            if ($pos_tmp ne $current_pos_tmp && $pos_tmp < $current_pos_tmp-1)
                            { 
                                my %post_SNP_patterns_tmp;
                                undef %post_SNP_patterns_tmp;
					#print {$filehandle{$seed_id2}} "\n".$pos_tmp." POS\n";   
                                foreach my $nuc_tmp (keys %{$SNP_patterns_prev{$pos_tmp}})
                                {
								    foreach my $rank_tmp (keys %{$SNP_patterns_prev{$pos_tmp}{$nuc_tmp}})
                                    {
										foreach my $nuc_now_tmp2 (keys %SNP_patterns_now)
                                        {              
                                            if (exists($SNP_patterns_now{$nuc_now_tmp2}{$rank_tmp}))
                                            {   
												if (exists($post_SNP_patterns_tmp{$nuc_now_tmp2}{$nuc_tmp}))
                                                {
                                                    my $count_tmp = $post_SNP_patterns_tmp{$nuc_now_tmp2}{$nuc_tmp}+1;
                                                    $post_SNP_patterns_tmp{$nuc_now_tmp2}{$nuc_tmp} = $count_tmp;
                                                }
                                                else
                                                {
                                                    $post_SNP_patterns_tmp{$nuc_now_tmp2}{$nuc_tmp} = '1';
                                                }
                                            }
                                        }
                                    }
                                }
                                
                                my $verified_tmp = "";
                                my %tmp_check_list;
                                undef %tmp_check_list;
								my $gap_count_tmp = keys %{$SNP_patterns_prev{$pos_tmp}{"-"}};
								
                                foreach my $nuc_now_tmp3 (keys %post_SNP_patterns_tmp)
                                {      
                                    my $count_nuc_now_tmp = keys %{$SNP_patterns_now{$nuc_now_tmp3}};
									
                #print {$filehandle{$seed_id2}} $nuc_now_tmp3." NUC_NOW\n";
				
                                    if (($count_nuc_now_tmp > 1 || $total_nuc_count < 11) && $count_nuc_now_tmp > $total_nuc_count*0.095)
                                    {
                                        my $pipi = 0.83;
										if ($count_nuc_now_tmp < 10 || $nuc_now_tmp3 eq "-")
										{
											$pipi = '0.85';
										}
										if ($count_nuc_now_tmp < 5)
										{
											$pipi = '0.95';
										}										
										
										foreach my $nuc_prev_tmp3 (keys %{$post_SNP_patterns_tmp{$nuc_now_tmp3}})
                                        {
                                            $count_nuc_now_tmp = keys %{$SNP_patterns_now{$nuc_now_tmp3}};
											my $count_nuc_prev_tmp = keys %{$SNP_patterns_prev{$pos_tmp}{$nuc_prev_tmp3}};
											if ($nuc_now_tmp3 ne "-" && $nuc_prev_tmp3 ne "-" && $find_haps_in_seed ne "" && $count_nuc_now_tmp > 9)
											{
												$pipi = '0.78';
											}
											
											if ($nuc_prev_tmp3 ne "-" && $gap_count_tmp < $total_nuc_count*0.15)
											{
												$count_nuc_now_tmp -= $post_SNP_patterns_tmp{$nuc_now_tmp3}{"-"};
												#print {$filehandle{$seed_id2}} $post_SNP_patterns_tmp{$nuc_now_tmp3}{"-"}." NUC_NOW_RANK_COUNT_ADJUST\n";
											}
											if ($nuc_now_tmp3 ne "-" && $nucs{"-"} < $total_nuc_count*0.15)
											{
												$count_nuc_prev_tmp -= $post_SNP_patterns_tmp{"-"}{$nuc_prev_tmp3};
												#print {$filehandle{$seed_id2}} $post_SNP_patterns_tmp{"-"}{$nuc_prev_tmp3}." NUC_PREV_RANK_COUNT_ADJUST\n";
											}
									#print {$filehandle{$seed_id2}} $nuc_prev_tmp3." NUC_PREV\n";
				#print {$filehandle{$seed_id2}} $count_nuc_prev_tmp." NUC_PREV_RANK_COUNT\n";
				#print {$filehandle{$seed_id2}} $count_nuc_now_tmp." NUC_NOW_RANK_COUNT\n";
				#print {$filehandle{$seed_id2}} $post_SNP_patterns_tmp{$nuc_now_tmp3}{$nuc_prev_tmp3}." MATCHES\n";  	
											if ($post_SNP_patterns_tmp{$nuc_now_tmp3}{$nuc_prev_tmp3} > $count_nuc_now_tmp*$pipi && $count_nuc_now_tmp > $post_SNP_patterns_tmp{$nuc_now_tmp3}{$nuc_prev_tmp3}*$pipi
												&& $count_nuc_now_tmp > $count_nuc_prev_tmp*$pipi && $count_nuc_prev_tmp > $count_nuc_now_tmp*$pipi)
                                            {
                                                if (exists($tmp_check_list{$nuc_prev_tmp3}))
                                                {
                                                    next POST_SNP_PATTERNS_TMP_NP;
                                                }
                                                $tmp_check_list{$nuc_prev_tmp3} = undef;
												if ($verified_tmp eq "yes")
												{
													$verified_tmp = "yes2";
												}
												else
												{
													$verified_tmp = "yes";
												}
                                            }
                                            elsif ($post_SNP_patterns_tmp{$nuc_now_tmp3}{$nuc_prev_tmp3} < $total_nuc_count*0.15)
                                            {              
                                            }
                                            else
                                            {
                                                next POST_SNP_PATTERNS_TMP_NP;
                                            }
                                        }
                                    }
                                }
                                if ($verified_tmp eq "yes2")
                                {
                                    print {$filehandle{$seed_id2}} $pos_tmp." POS_MATCH\n";
									if ($pos_tmp > $current_pos_tmp-150)
									{
										$local_pattern_matches++;
										if ($pos_tmp > $current_pos_tmp-25)
										{
											$local_pattern_matches2++;
										}	
									}		
									
									my $one_tmp = "";
									my $two_tmp = "";
									foreach my $nuc_tmp (keys %{$SNP_patterns_prev{$pos_tmp}})
									{
										foreach my $rank_tmp (keys %{$SNP_patterns_prev{$pos_tmp}{$nuc_tmp}})
										{
											if ($rank_tmp eq '1')
											{
												$one_tmp = $nuc_tmp;
											}
											if ($rank_tmp eq '2')
											{
												$two_tmp = $nuc_tmp;
											}
										}
									}
									if ($one_tmp eq $two_tmp)
									{
										$SNP_patterns_prev_match{$pos_tmp} = $one_tmp;
									}
									
                                    foreach my $nuc_now_tmp4 (keys %SNP_patterns_now)
                                    {
                                        my $highest_count_tmp2 = '0';
										my $highest_nuc_tmp2 = "";
										my $best_nuc_prev = "";
										foreach my $rank_now_tmp4 (keys %{$SNP_patterns_now{$nuc_now_tmp4}})
                                        {
                                            foreach my $nuc_prev_tmp4 (keys %{$SNP_patterns_prev{$pos_tmp}})
                                            {    
												foreach my $rank_prev_tmp4 (keys %{$SNP_patterns_prev{$pos_tmp}{$nuc_prev_tmp4}})
                                                {
                                                    if ($rank_prev_tmp4 eq $rank_now_tmp4)
                                                    {
                                                        my $highest_count_tmp = '0';
                                                        my $highest_nuc_tmp = "";
														
                                                        foreach my $nuc_now_tmp5 (keys %post_SNP_patterns_tmp)
                                                        {                                                          
                                                            foreach my $nuc_prev_tmp5 (keys %{$post_SNP_patterns_tmp{$nuc_now_tmp5}})
                                                            {
                                                                if ($nuc_prev_tmp5 eq $nuc_prev_tmp4)
                                                                {
                                                                    my $count_tmp = $post_SNP_patterns_tmp{$nuc_now_tmp5}{$nuc_prev_tmp5};
                                                                    if ($count_tmp > $highest_count_tmp)
                                                                    {
                                                                        $highest_count_tmp = $count_tmp;
                                                                        $highest_nuc_tmp = $nuc_now_tmp5;
                                                                    }
																	if ($count_tmp > $highest_count_tmp2)
                                                                    {
																		$highest_count_tmp2 = $count_tmp;
																		$highest_nuc_tmp2 = $nuc_now_tmp5;
																		$best_nuc_prev = $nuc_prev_tmp4;
                                                                    }   
                                                                }
                                                            }
                                                        }
                                                        if ($highest_nuc_tmp ne "")
                                                        {
                                                            $post_SNP_patterns{$rank_now_tmp4}{$highest_nuc_tmp}{$pos_tmp} = $nuc_prev_tmp4;
                                                        }
                                                    }
                                                }	
                                            }
                                        }
										if ($highest_nuc_tmp2 ne "")
										{
											$find_haps_SNPs{$highest_nuc_tmp2}{$pos_tmp} = $best_nuc_prev;
											#print {$filehandle{$seed_id2}} $highest_nuc_tmp2." HIGH_NUC\n";
											#print {$filehandle{$seed_id2}} $highest_count_tmp2." HIGH_NUC_COUNT\n";
											#print {$filehandle{$seed_id2}} $best_nuc_prev." NUC_PREV\n";
										}
                                    }
                                    $prev_pos_count++;
                                }
                            }
                        }
 #-----------------------------------------------------------------------------------------------------------------------                                    
                        my %split_patterns_final;
                        undef %split_patterns_final;
                        undef %split_patterns_final_score;

                        foreach my $rank_tmp6 (keys %post_SNP_patterns)
                        {                       
                            my $highest_count_tmp = '0';
                            my $highest_nuc_tmp = "";
                            foreach my $nuc_now_tmp6 (keys %{$post_SNP_patterns{$rank_tmp6}})
                            {
                                my $count_tmp = '0';
								foreach my $pos_tmp6 (keys %{$post_SNP_patterns{$rank_tmp6}{$nuc_now_tmp6}})
								{
									if ($post_SNP_patterns{$rank_tmp6}{$nuc_now_tmp6}{$pos_tmp6} eq "-" || $nuc_now_tmp6 eq "-")
									{									
										if ($find_haps_in_seed eq "")
										{
											$count_tmp += 0.5;
										}
										else
										{
											$count_tmp += 0.25;
										}
									}
									else
									{
										$count_tmp += 1.5;
									}
								}
								#my $count_tmp = keys %{$post_SNP_patterns{$rank_tmp6}{$nuc_now_tmp6}};
                                if ($count_tmp > $highest_count_tmp)
                                {
                                    $highest_count_tmp = $count_tmp;
                                    $highest_nuc_tmp = $nuc_now_tmp6;
                                }
                                elsif (exists($SNP_patterns_now{$nuc_now_tmp6}{$rank_tmp6}))
                                {
                                    if ($count_tmp eq $highest_count_tmp)
                                    {
                                        $highest_count_tmp = $count_tmp;
                                        $highest_nuc_tmp = $nuc_now_tmp6;
                                    }
                                }
                            }
                            $split_patterns_final{$rank_tmp6} = $highest_nuc_tmp;
                            $split_patterns_final_score{$highest_nuc_tmp}{$rank_tmp6} = $highest_count_tmp;
                        }
                        						
						if (keys %{$trace_back_split_NP{$id}} > 0 && $SNP_check ne "")
						{
							my $min_pos = $position+length($best_extension)-20;
							
							foreach my $min_pos_tmp (keys %{$trace_back_split_NP{$id}})
                            {
								if ($min_pos_tmp > $min_pos && $min_pos_tmp < $min_pos+40)
								{
									my $last_10 = substr $best_extension, -10, 10;
									my $last_10_prev = $trace_back_split_NP{$id}{$min_pos_tmp};
									my $N_check = $last_10_prev =~ tr/N/N/;
	
									print {$filehandle{$seed_id2}} $last_10_prev." ".$last_10." TRACE_BACK_POSITION!!!!!!!!!!!!!!!!\n";
									if ($last_10 eq $last_10_prev)
									{			
										$trace_back_check = "yes";
									}
									elsif ($N_check > 0)
									{
										$last_10_prev =~ tr/N/\./;
										if ($last_10 =~ m/$last_10_prev/)
										{
											$trace_back_check = "yes";
										}
									}
								}
							}
						}				

                        if (keys %split_patterns_final eq '0' && ($SNP_check ne "" || $trace_back_check eq "yes" || ($NP_reads_support eq "yes2" && $clipped_ext eq "yes" && $best_extension eq "")))
                        {
                            print {$filehandle{$seed_id2}} $SNP_check." SNP_CHECK\n";
                            foreach my $rank_tmp4 (sort {$a <=> $b} keys %subject_list)
                            {
								my $nuc_tmp = $nucs_by_rank{$rank_tmp4};
                                $split_patterns_final{$rank_tmp4} = $nuc_tmp;
                                $split_patterns_final_score{$nuc_tmp}{$rank_tmp4} = '1';
                            }
                        }
				
                  
                        my %reads_to_remove;
                        undef %reads_to_remove;
                        my $remove_reads = "";
                        my %average_rank_score;
                        undef %average_rank_score;
                                
                        $post_pattern_match = "";
                        my $post_pattern_match_extra = "";
                        my $post_pattern_match_average = "";
                        my %longest_match;
                        undef %longest_match;
						my %longest_match2;
                        undef %longest_match2;
                        my $longest_longest_match = "";
                        my $longest_longest_match_nuc = "";
						my $count_post_pattern_match3 = '0';
						my $find_haps_found = "";
                        
                        foreach my $rank_tmp13 (sort {$a <=> $b} keys %split_patterns_final)
                        {
                            print {$filehandle{$seed_id2}} $rank_tmp13." ".$nucs_by_rank{$rank_tmp13}." ".$split_patterns_final{$rank_tmp13}." NUC_TMP_PREV\n";
                        }
						
						if ($NP_reads_support eq "yes2" && $clipped_ext eq "yes" && length($best_extension) < 1 && $SNP_check ne "")
						{	
							foreach my $nuc_tmp13 (sort {$a <=> $b} keys %split_patterns_final_score)
							{
								if (exists($PB_split_nucs{$nuc_tmp13}))
								{}
								else
								{
									delete $split_patterns_final_score{$nuc_tmp13};
								}
							}
						}
						
                        foreach my $nuc_tmp13 (sort {$a <=> $b} keys %split_patterns_final_score)
                        {
                            my $total_count_tmp = '0';
                            my $total_score_count_tmp = '0';
                            my $count_ranks_tmp = keys %{$split_patterns_final_score{$nuc_tmp13}};
                            foreach my $rank_tmp13 (keys %{$split_patterns_final_score{$nuc_tmp13}})
                            {
                                $total_count_tmp++;
                                $total_score_count_tmp += $split_patterns_final_score{$nuc_tmp13}{$rank_tmp13};
                                if (exists($rank_to_id{$rank_tmp13}) && $count_ranks_tmp > 0.15*$total_nuc_count)
                                {
                                    my $id_tmp3 = $rank_to_id{$rank_tmp13};
                                    #print {$filehandle{$seed_id2}} $long_read_end_pos_save{$id_tmp3}." LONG_READ_END_TMP\n";
                                    if (exists($read_start_pos_rej{$id_tmp3}))
                                    {}
                                    elsif (exists($longest_match{$nuc_tmp13}))
                                    {
                                        if ($alignment_length_save{$id_tmp3} > $longest_match{$nuc_tmp13})
                                        {
                                            $longest_match{$nuc_tmp13} = $alignment_length_save{$id_tmp3};
                                        }
										$longest_match2{$nuc_tmp13}{$alignment_length_save{$id_tmp3}}{$id_tmp3} = undef;
                                    }
                                    else
                                    {
                                        $longest_match{$nuc_tmp13} = $alignment_length_save{$id_tmp3};
										$longest_match2{$nuc_tmp13}{$alignment_length_save{$id_tmp3}}{$id_tmp3} = undef;
                                    }
                                }
                            }
							$total_score_count_tmp *= ($count_ranks_tmp/$total_nuc_count)*2;
                            print {$filehandle{$seed_id2}} $nuc_tmp13." ".$total_score_count_tmp." SCORE ".$total_count_tmp." COUNT\n";
                            
                            #if ($confirmed_reads_count_NP > 4 && $skip_confirmed eq "" && $only_confirmed eq "yes" && $total_score_count_tmp/$total_count_tmp > 2 && $total_score_count_tmp/$total_count_tmp > $total_count_prev_patterns*0.2 && $total_count_tmp > 3)
                            #{
                                #$skip_confirmed = "yes";
                                #print {$filehandle{$seed_id2}} "SKIP_CONFIRMED2\n";
                                #$best_extension = "";
                                #goto SKIP_CONFIRMED_NP;
                            #}
#Split haplotypes in intial seed--------------------------------------------------------------------------------------------------						
							if ($find_haps_in_seed ne "" && $total_score_count_tmp/$total_count_tmp >= 1.5 && keys %find_haps_SNPs eq $ploidy)
							{
								if ($find_haps_found eq "")
								{
									$find_haps_found = "ye";
								}
								else
								{
									$found_haps_in_seed = "yes";
									my $second_hap = "";

									foreach my $nucie_now_tmp2 (sort {$a <=> $b} keys %find_haps_SNPs)
									{
										print {$filehandle{$seed_id2}} $nucie_now_tmp2." NUC_NOW\n";
										foreach my $pos_tmp8 (sort {$b <=> $a} keys %{$find_haps_SNPs{$nucie_now_tmp2}})
										{
											print {$filehandle{$seed_id2}} $pos_tmp8." CONNECT_POS\n";
											my $nucie_prev_tmp2 = $find_haps_SNPs{$nucie_now_tmp2}{$pos_tmp8};
											$nucie_prev_tmp2 =~ tr/actgn/ACTGN/;
											$nucie_prev_tmp2 =~ tr/-//d;
											
											if ($second_hap eq "")
											{
		print {$filehandle{$seed_id2}} $nucie_prev_tmp2." NUC_PREV\n";
												my $read_tmp = $extensions_seed{"HAP1"};
								print {$filehandle{$seed_id2}} $read_tmp." EXT\n";	
												substr $read_tmp, $pos_tmp8-$position, 1, $nucie_prev_tmp2;
									print {$filehandle{$seed_id2}} $read_tmp." EXT\n";		
												$extensions_seed{"HAP1"} = $read_tmp;
											}
											else
											{
												print {$filehandle{$seed_id2}} $nucie_prev_tmp2." NUC_PREV\n";
												my $read_tmp = $extensions_seed{"HAP2"};
								print {$filehandle{$seed_id2}} $read_tmp." EXT\n";	
												substr $read_tmp, $pos_tmp8-$position, 1, $nucie_prev_tmp2;
									print {$filehandle{$seed_id2}} $read_tmp." EXT\n";		
												$extensions_seed{"HAP2"} = $read_tmp;
											}
										}
										$nucie_now_tmp2 =~ tr/actgn/ACTGN/;
										$nucie_now_tmp2 =~ tr/-//d;
										if ($second_hap eq "")
										{
											$extensions_seed{"HAP1"} .= $nucie_now_tmp2;
										}
										else
										{
											$extensions_seed{"HAP2"} .= $nucie_now_tmp2;
										}
										
										$second_hap = "yes";
									}

                                    $hap_position = "yes";
									$best_extension .= "N";
									$best_extension_part .= "N";
									$quality_scores_tmp{length($best_extension)} = '0'." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
									$nuc_match = "N";
									goto SKIP_INPUT_BLAST3_NP;
								}
							}
#----------------------------------------------------------------------------------------------------------------------------------------------------			
							$total_count_prev_patterns = '0';
							foreach my $posie_tmpie (sort {$b <=> $a} keys %total_count_prev_patterns)
							{
								if ($posie_tmpie > length($best_extension)-150)
								{
									$total_count_prev_patterns++;
								}
							}
							my $total_count_prev_patterns_all = keys %total_count_prev_patterns;
             
                            if (((keys %split_patterns_final > 0.88*$total_nuc_count && ($total_score_count_tmp/$total_count_tmp > 1.5 || ($trace_back_check eq "yes" && $total_count_tmp > $total_nuc_count*0.2)) && $total_score_count_tmp > 0 && ($total_count_prev_patterns < 30
								|| ($total_count_prev_patterns < 20 && $total_nuc_count > 5) || ($total_score_count_tmp/$total_count_tmp > 10))) || $SNP_check eq "yes2")
								&& (length($best_extension) > 50 || $clipped_ext ne "") && $total_nuc_count > 5)
                            {            
                                print {$filehandle{$seed_id2}} $total_count_prev_patterns." TOTAL_COUNT_PREV_PATTERNS\n";
								print {$filehandle{$seed_id2}} $total_count_prev_patterns_all." TOTAL_COUNT_PREV_PATTERNS_ALL\n";
								
								if ($post_pattern_match eq "yes")
                                {
                                    $post_pattern_match = "yes2";
                                }
                                elsif ($post_pattern_match eq "")
                                {
                                    $post_pattern_match = "yes";
                                }
								my $vv = '0.3';
								my $ll = '0.15';
								if ($total_nuc_count > 8 && $count_ranks_tmp > 4)
								{
									$vv = '0.25';
									$ll = '0.1';
								}
								if ((($total_count_prev_patterns < 15 && $total_count_prev_patterns_all < 35) || ($total_score_count_tmp/$total_count_tmp > (($total_count_prev_patterns*$vv)+(10/$total_nuc_count)) &&
									 (($local_pattern_matches > 5 && ($count_ranks_tmp > $total_nuc_count*0.35 && $count_ranks_tmp < $total_nuc_count*0.65))
									  || $total_score_count_tmp/$total_count_tmp > (($total_count_prev_patterns_all*$ll)+(10/$total_nuc_count))))
									 || $total_score_count_tmp/$total_count_tmp > (33+((10/$total_nuc_count)*10)) || ($local_pattern_matches > 0.7*$total_count_prev_patterns && $local_pattern_matches > 10)
									 || $local_pattern_matches2 > 4) && ($total_count_tmp > 2 || $prev_pos_count > $total_count_prev_patterns_all*0.25))
                                {
									if ($post_pattern_match_average eq "yes")
									{
									}
									elsif ($post_pattern_match_average eq "ye")
									{
										$post_pattern_match_average = "yes";
									}
									else
									{
										$post_pattern_match_average = "ye";
										if ($total_score_count_tmp/$total_count_tmp > (40+((10/$total_nuc_count)*10)) && $count_ranks_tmp > 10 &&
										   (($local_pattern_matches > 0.68*$total_count_prev_patterns && $local_pattern_matches > 10) ||($count_ranks_tmp > $total_nuc_count*0.35 && $count_ranks_tmp < $total_nuc_count*0.65)))
										{
											$post_pattern_match_average = "yes";
										}
									}
                                }
                                if ((($total_score_count_tmp/$total_count_tmp > 1.6 && $SNP_check ne "") || $total_score_count_tmp/$total_count_tmp > 10)
									&& (($count_ranks_tmp > 3 && $count_ranks_tmp < $total_nuc_count-3) || $total_score_count_tmp/$total_count_tmp > 60) || $local_pattern_matches2 > 4)
                                {
                                    if ($post_pattern_match_extra eq "yes")
									{
									}
									elsif ($post_pattern_match_extra eq "ye")
									{
										$post_pattern_match_extra = "yes";
									}
									else
									{
										$post_pattern_match_extra = "ye";
										if ($total_score_count_tmp/$total_count_tmp > (40+((10/$total_nuc_count)*10)) && $count_ranks_tmp > 10 &&
										   (($local_pattern_matches > 0.68*$total_count_prev_patterns && $local_pattern_matches > 10) ||($count_ranks_tmp > $total_nuc_count*0.35 && $count_ranks_tmp < $total_nuc_count*0.65)))
										{
											$post_pattern_match_average = "yes";
										}
									}
									
									if ($first_split_pos eq "")
                                    {
                                        $first_split_pos = length($best_extension);
                                    }
                                }
                                print {$filehandle{$seed_id2}} $post_pattern_match_extra." POST_PATTERN_MATCH_EXTRA\n";
								print {$filehandle{$seed_id2}} $post_pattern_match_average." POST_PATTERN_MATCH_AVERAGE\n";
#Multi_match removal-------------------------------------------------------------------
								if ($post_pattern_match_extra eq "yes" && $post_pattern_match_average eq "yes" || ($NP_reads_support eq "yes2" && $clipped_ext eq "yes" && length($best_extension) < 2 && $SNP_check ne ""))
								{
									my $count_multi_match_tmp = '0';
									my $first_no_multi_match = "";
									my %gap_list_tmp;
									undef %gap_list_tmp;
									
									foreach my $rank_tmp5 (sort {$a <=> $b} keys %{$split_patterns_final_score{$nuc_tmp13}})
									{  
										if (exists($rank_to_id{$rank_tmp5}))
										{
											my $id_tmp = $rank_to_id{$rank_tmp5};
											if (exists($multi_match{$id_tmp}))
											{
												$count_multi_match_tmp++;
												foreach my $gap_pos_tmp (sort {$a <=> $b} keys %{$multi_match{$id_tmp}})
												{
													if (exists($gap_list_tmp{$gap_pos_tmp}))
													{
														$gap_list_tmp{$gap_pos_tmp} += 1;
													}
													else
													{
														$gap_list_tmp{$gap_pos_tmp} = '1';
													}  
												}
											}
											elsif ($first_no_multi_match eq "")
											{
												$first_no_multi_match = $rank_tmp5;
											}
											else
											{
												last;
											}
										}
									}
									my $gap_pos_prev = "";
									my $count_tmp = '0';
									my $max_count_tmp = '0';
									my $max_count_pos_tmp = '0';
									foreach my $gap_pos_tmp (sort {$a <=> $b} keys %gap_list_tmp)
									{
										if ($gap_pos_prev eq "")
										{
											$count_tmp = '1';
										}
										elsif ($gap_pos_tmp < $gap_pos_prev+350)
										{
											$count_tmp++;
											if ($count_tmp > $max_count_tmp)
											{
												$max_count_tmp = $count_tmp;
												$max_count_pos_tmp = $gap_pos_tmp;
											}
										}
										else
										{
											$count_tmp = '1';
										}
										$gap_pos_prev = $gap_pos_tmp;
									}
									
									if (($first_no_multi_match eq "" || $first_no_multi_match > 10 || ($first_no_multi_match > 5 && $count_multi_match_tmp > 2)) && $max_count_tmp > 1
										&& ($post_pattern_match_extra eq "yes" || ($max_count_tmp > 1 && $max_count_tmp > $count_multi_match_tmp*0.4)))
									{
										my $gap_region_tmp = substr $read, $max_count_pos_tmp-150, 900;
										my $N_count_tmp = $gap_region_tmp =~ tr/N/N/;
										
										my $first_rank_tmp = "";
										if ($N_count_tmp < 60)
										{
											foreach my $rank_tmp5 (sort {$a <=> $b} keys %{$split_patterns_final_score{$nuc_tmp13}})
											{
												if ($first_rank_tmp eq "")
												{
													$first_rank_tmp = $rank_tmp5;
												}
												if (exists($rank_to_id{$rank_tmp5}))
												{
													my $id_tmp = $rank_to_id{$rank_tmp5};
													my $exclude_pos = $alignment_length_save{$id_tmp};
													if (length($read) < $alignment_length_save{$id_tmp})
													{
														$exclude_pos = length($read);
													}
													if ($hap_tag eq "HAP1" && $repetitive_detect1 eq "" && $repetitive_detect2 eq "111")
													{
														$exclude_reads_hap1_NP{$id_tmp} = $position+$exclude_pos;
													}
													elsif ($hap_tag eq "HAP2" && $repetitive_detect1 eq "" && $repetitive_detect2 eq "111")
													{
														$exclude_reads_hap2_NP{$id_tmp} = $position+$exclude_pos;
													}
													$reads_to_remove{$first_rank_tmp}{$rank_tmp5} = undef;
													$remove_reads = "yes";
												}
											}
										}
										if ($remove_reads eq "yes")
										{
											print {$filehandle{$seed_id2}} $first_no_multi_match." MULTI_MATCH_PATTERN\n";
											goto REMOVE_READS_NP;
										}
									}
								}
                                print {$filehandle{$seed_id2}} $longest_match{$nuc_tmp13}." LONGEST_MATCH\n";
                                    
                                if (($longest_longest_match eq "" || ($longest_match{$nuc_tmp13} > $longest_longest_match && $total_count_tmp > $total_nuc_count*0.2)) && $nuc_tmp13 ne "" && $nucs{$nuc_tmp13} > $total_nuc_count*0.22)
                                {
                                    $longest_longest_match = $longest_match{$nuc_tmp13};
                                    $longest_longest_match_nuc = $nuc_tmp13;
                                }
#---------------------------------------------------------------------------------------------
								if ($post_pattern_match_extra ne "" && $post_pattern_match_save eq "yes" && $total_score_count_tmp/$total_count_tmp > 3 && $post_pattern_match_average ne "" &&
                                    ($total_score_count_tmp/$total_count_tmp > 10 || ($SNP_check eq "yes2" && $total_count_prev_patterns < 25)))
                                {
									$count_post_pattern_match3++;
								}
                                if ($post_pattern_match_extra eq "yes" && $post_pattern_match_save eq "yes" && $total_score_count_tmp/$total_count_tmp > 3 && $post_pattern_match_average eq "yes" &&
                                    ($total_score_count_tmp/$total_count_tmp > 10 || ($SNP_check eq "yes2" && $total_count_prev_patterns < 25)))
                                {                                    
                                    if ($post_pattern_match_save eq "yes" && $total_score_count_tmp/$total_count_tmp > 9 && $post_pattern_match_average eq "yes" && ($total_count_prev_patterns < 15 || $total_count_tmp > 6 || $total_score_count_tmp/$total_count_tmp > 25))
                                    {
                                        my $best_extension_part_tmp = substr $best_extension, -2000;
										my $CG = $best_extension_part_tmp =~ tr/CGN/CGN/;
										if ($CG > 0.52*(length($best_extension_part_tmp)) && $total_score_count_tmp/$total_count_tmp < (10+(10/$total_nuc_count)) && $total_count_prev_patterns > 15)
										{
											$post_pattern_match = "yes2";
										}
										else
										{
											$post_pattern_match = "yes3";
											$post_pattern_match_count++;
										}

										print {$filehandle{$seed_id2}} $post_pattern_match." POST_PATTERN_MATCH\n";
                                    }
                                }
                                $post_pattern_match_save = "yes";
                            }
                        }
						
						if ($count_post_pattern_match3 < 2 && $post_pattern_match eq "yes3")
						{
							$post_pattern_match = "yes2";
						}
						
						if ($trace_back_check eq "yes" || ($NP_reads_support eq "yes2" && $clipped_ext eq "yes" && length($best_extension) < 2 && $SNP_check ne ""))
						{	
							$post_pattern_match = "yes3";
							$post_pattern_match_average = "yes";
							$post_pattern_match_extra = "yes";
						}
						if ($local_pattern_matches2 > 8)
						{
							$post_pattern_match = "yes3";
						}

     my $time_split2 = time;
     my $testy_time = $time_split-$time_split2;
     if ($testy_time > 0.3)
     {
        print {$filehandle{$seed_id2}} $testy_time." TIME_CHECK_SPLIT\n";
     }        
                        
						my $min_pos = $position+length($best_extension)-20;
						my $track_check = "";
						if ($nucs{"-"} < 0.25*$total_nuc_count)
						{
							foreach my $min_pos_tmp (keys %{$track_split_NP{$id}})
							{
								if ($min_pos_tmp > $min_pos && $min_pos_tmp < $min_pos+40)
								{
									my $last_10 = substr $best_extension, -10, 10;
									my @last_10_prev = split /\+/, $track_split_NP{$id}{$min_pos_tmp};
									my $last_10_prev = $last_10_prev[0];
									my $N_check = $last_10_prev =~ tr/N/N/;
	
									print {$filehandle{$seed_id2}} $last_10_prev." ".$last_10." TRACK_POSITION!!!!!!!!!!!!!!!!\n";
									if ($last_10 eq $last_10_prev)
									{			
										$track_check = $last_10_prev[1];
									}
									elsif ($N_check > 0)
									{
										$last_10_prev =~ tr/N/\./;
										if ($last_10 =~ m/$last_10_prev/)
										{
											$track_check = $last_10_prev[1];
										}
									}
								}
							}
							if ($track_check ne "")
							{
								foreach my $nuc_tmp15 (sort {$a <=> $b} keys %split_patterns_final_score)
								{
									if ($nuc_tmp15 ne $track_check)
									{
										foreach my $rank_tmp15 (keys %{$split_patterns_final_score{$nuc_tmp15}})
										{
											if ($rank_tmp15 > 5)
											{
												$reads_to_remove{$rank_tmp15}{$rank_tmp15} = undef;
											}
										}
									}
								}
								$remove_reads = "yes";
								print {$filehandle{$seed_id2}} $track_check." REMOVE_BY_TRACK_BACK\n";
								goto REMOVE_READS_NP;
							}
						}
						
#--------------------------------------------------------------------------------------------------------------------------						
                        if ($mismatches_tmp_check eq "" && $find_haps_in_seed eq "")
                        {
                            $mismatches_tmp_check = "yes";
                        
							foreach my $rank_tmp (sort {$a <=> $b} keys %subject_list)
							{
								if (exists($rank_to_id{$rank_tmp}))
								{
									my $id_tmp = $rank_to_id{$rank_tmp};
									if (exists($store_mismatches_NP{$id_tmp}))
									{
										foreach my $posie (keys %{$store_mismatches_NP{$id_tmp}})
										{
											#if ($posie > 70)
											#{
												$mismatches_tmp{$posie}{$rank_tmp} = $store_mismatches_NP{$id_tmp}{$posie};
											#}
										}
									}
									if (exists($store_mismatches_all_NP{$id_tmp}))
									{
										foreach my $posie (keys %{$store_mismatches_all_NP{$id_tmp}})
										{
											#if ($posie > 70)
											#{
												$mismatches_tmp_all{$posie}{$rank_tmp} = $store_mismatches_all_NP{$id_tmp}{$posie};
											#}
										}
									}
								}
							}   

							foreach my $posie (sort {$a <=> $b} keys %mismatches_tmp)
							{          
								my $count_mm_tmp = keys %{$mismatches_tmp{$posie}};
	
								if (($count_mm_tmp > 2 || ($total_nuc_count < 10 && $count_mm_tmp > 1)) && $count_mm_tmp > $total_nuc_count*0.07)
								{                                    
									my $count_below_5 = '0';
									my $ii = '1';
									while ($ii < 5)
									{
										if (exists($mismatches_tmp{$posie}{$ii}))
										{
											$count_below_5++;
										}
										$ii++;
									}
									
									if ($count_below_5 < 5)
									{
										$mismatch_score++;
										foreach my $rank_tmp (keys %{$mismatches_tmp{$posie}})
										{
											$reads_mismatch{$rank_tmp}{$posie} = undef;
											$reads_mismatchb{$posie}{$rank_tmp} = undef;
										}
									}
								}
							}
                        
							foreach my $posie (sort {$a <=> $b} keys %mismatches_tmp_all)
							{          
								my $count_mm_tmp = keys %{$mismatches_tmp_all{$posie}};
	
								if (($count_mm_tmp > 2 || ($total_nuc_count < 10 && $count_mm_tmp > 1)) && $count_mm_tmp > $total_nuc_count*0.07)
								{                                    
									my $count_below_5 = '0';
									my $ii = '1';
									while ($ii < 5)
									{
										if (exists($mismatches_tmp_all{$posie}{$ii}))
										{
											$count_below_5++;
										}
										$ii++;
									}
									
									if ($count_below_5 < 5)
									{
										$mismatch_score_all++;
										foreach my $rank_tmp (keys %{$mismatches_tmp_all{$posie}})
										{
											$reads_mismatch_all{$rank_tmp}{$posie} = undef;
											$reads_mismatchb_all{$posie}{$rank_tmp} = undef;
										}
									}
								}
							}
                        
							foreach my $rank_tmp2 (sort {$a <=> $b} keys %reads_mismatch)
							{
								foreach my $posie_tmp2 (sort {$a <=> $b} keys %{$reads_mismatch{$rank_tmp2}})
								{
									foreach my $rank_tmp3 (sort {$a <=> $b} keys %{$reads_mismatchb{$posie_tmp2}})
									{
										if ($rank_tmp3 ne $rank_tmp2)
										{
											if (exists($reads_mismatch2_tmp{$rank_tmp2}{$rank_tmp3}))
											{
												$reads_mismatch2_tmp{$rank_tmp2}{$rank_tmp3}{$posie_tmp2} = undef;
											}
											elsif (exists($reads_mismatch2_tmp{$rank_tmp3}{$rank_tmp2}))
											{
												$reads_mismatch2_tmp{$rank_tmp3}{$rank_tmp2}{$posie_tmp2} =undef;
											}
											else
											{
												$reads_mismatch2_tmp{$rank_tmp2}{$rank_tmp3}{$posie_tmp2} = undef;
											}
										}
									}
								}
							}
							foreach my $rank_tmp2 (sort {$a <=> $b} keys %reads_mismatch2_tmp)
							{
								foreach my $rank_tmp3 (sort {$a <=> $b} keys %{$reads_mismatch2_tmp{$rank_tmp2}})
								{
									my $count_tmp = keys %{$reads_mismatch2_tmp{$rank_tmp2}{$rank_tmp3}};
									if ($count_tmp > 2 || ($total_nuc_count < 10 && $count_tmp > 1))
									{  
										foreach my $posie_tmp3 (sort {$a <=> $b} keys %{$reads_mismatch2_tmp{$rank_tmp2}{$rank_tmp3}})
										{
											$reads_mismatch2{$count_tmp}{$rank_tmp2}{$rank_tmp3}{$posie_tmp3} = undef;
										}
									}
									else
									{
										foreach my $posie_tmp3 (sort {$a <=> $b} keys %{$reads_mismatch2_tmp{$rank_tmp2}{$rank_tmp3}})
										{
											my $count_tmp2 = keys %{$reads_mismatchb{$posie_tmp3}};
											if ($count_tmp2 > 3 && $count_tmp2 > $total_nuc_count*0.15)
											{
												$reads_mismatch2{$count_tmp}{$rank_tmp2}{$rank_tmp3}{$posie_tmp3} = undef;
											} 
										}
									}
								}
							}
							
							foreach my $rank_tmp2 (sort {$a <=> $b} keys %reads_mismatch_all)
							{
								foreach my $posie_tmp2 (sort {$a <=> $b} keys %{$reads_mismatch_all{$rank_tmp2}})
								{
									foreach my $rank_tmp3 (sort {$a <=> $b} keys %{$reads_mismatchb_all{$posie_tmp2}})
									{
										if ($rank_tmp3 ne $rank_tmp2)
										{
											if (exists($reads_mismatch2_tmp_all{$rank_tmp2}{$rank_tmp3}))
											{
												$reads_mismatch2_tmp_all{$rank_tmp2}{$rank_tmp3}{$posie_tmp2} = undef;
											}
											elsif (exists($reads_mismatch2_tmp_all{$rank_tmp3}{$rank_tmp2}))
											{
												$reads_mismatch2_tmp_all{$rank_tmp3}{$rank_tmp2}{$posie_tmp2} =undef;
											}
											else
											{
												$reads_mismatch2_tmp_all{$rank_tmp2}{$rank_tmp3}{$posie_tmp2} = undef;
											}
										}
									}
								}
							}
							foreach my $rank_tmp2 (sort {$a <=> $b} keys %reads_mismatch2_tmp_all)
							{
								foreach my $rank_tmp3 (sort {$a <=> $b} keys %{$reads_mismatch2_tmp_all{$rank_tmp2}})
								{
									my $count_tmp = keys %{$reads_mismatch2_tmp_all{$rank_tmp2}{$rank_tmp3}};
									if ($count_tmp > 2 || ($total_nuc_count < 10 && $count_tmp > 1))
									{  
										foreach my $posie_tmp3 (sort {$a <=> $b} keys %{$reads_mismatch2_tmp_all{$rank_tmp2}{$rank_tmp3}})
										{
											$reads_mismatch2_all{$count_tmp}{$rank_tmp2}{$rank_tmp3}{$posie_tmp3} = undef;
										}
									}
								}
							}
							print {$filehandle{$seed_id2}} $mismatch_score_all." MISMATCH_SCORE_ALL\n";
						}
#--------------------------------------------------------------------------------------------------------------------------     
                        my %nucs_for_split_extra;
                        undef %nucs_for_split_extra;
                        my %pos_list_number;
                        undef %pos_list_number;
                        my $SNP_pattern_average = '0';
                        my $SNP_pattern_average_tmp = '0';
						my $high_score_save = "";
						my $highest_first_no_match = "";
                            
                        if ($mismatch_score > 0 && ($post_pattern_match_extra eq "yes" || $SNP_check eq "yes2" || $post_pattern_match eq "yes3") && $find_haps_in_seed eq "")
                        {                             
                            foreach my $rank_tmp (sort {$a <=> $b} keys %reads_mismatch)
                            {
                                my $pos_list = "";
                                my %post_SNP_matching_score;
                                undef %post_SNP_matching_score;
                                my $pos_list_number = '0';
                                
                                foreach my $posie_tmp (sort {$a <=> $b} keys %{$reads_mismatch{$rank_tmp}})
                                {                               
                                    if (exists($split_patterns_final{$rank_tmp}))
                                    {
                                        my $nuc_tmp18 = $split_patterns_final{$rank_tmp};
                                        if (exists($post_SNP_matching_score{$nuc_tmp18}))
                                        {
                                            my $count_tmp = $post_SNP_matching_score{$nuc_tmp18};
                                            $post_SNP_matching_score{$nuc_tmp18} = $count_tmp+1;
                                        }
                                        else
                                        {
                                            $post_SNP_matching_score{$nuc_tmp18} = '1';
                                        }                            
                                    }
                                    #print {$filehandle{$seed_id2}} $posie_tmp." POS ".$reads_mismatch{$rank_tmp}{$posie_tmp}."\n";
                                    if ($pos_list eq "")
                                    {
                                        $pos_list = $posie_tmp;
                                        $pos_list_number++;
                                    }
                                    else
                                    {
                                        $pos_list .= ";".$posie_tmp;
                                        $pos_list_number++;
                                    }
                                }
                                $pos_list_number{$rank_tmp} = $pos_list_number;
                                my $length_match_tmp = '1000';
                                if (exists($rank_to_id{$rank_tmp}))
                                {
                                    my $id_tmpi = $rank_to_id{$rank_tmp};
                                    $length_match_tmp = $alignment_length_save{$id_tmpi};
                                    if ($length_match_tmp > $position-$original_seed_length{$id})
                                    {
                                        $length_match_tmp -= $original_seed_length{$id} - ($position-$length_match_tmp);
                                    }
                                }
                                $SNP_pattern_average_tmp += $pos_list_number/$length_match_tmp;
                                
                                if ($pos_pattern_list_check eq "")
								{
									print {$filehandle{$seed_id2}} $rank_tmp." RANK\n";
									print {$filehandle{$seed_id2}} $pos_list."\n";
								}
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------                                                                
                                my $total_pattern_count = keys %{$reads_mismatch{$rank_tmp}};
                                foreach my $nuc_tmp19 (keys %post_SNP_matching_score)
                                {
                                    if ($post_SNP_matching_score{$nuc_tmp19} > 0.8*$total_pattern_count)
                                    {
                                        if (exists($SNP_patterns_prev2{$nuc_tmp19}))
                                        {
                                            my $count_tmp = $SNP_patterns_prev2{$nuc_tmp19};
                                            $SNP_patterns_prev2{$nuc_tmp19} = $count_tmp+1;
                                        }
                                        else
                                        {
                                            $SNP_patterns_prev2{$nuc_tmp19} = '1';
                                        }
										if ($pos_pattern_list_check eq "")
										{
											print {$filehandle{$seed_id2}} $nuc_tmp19." PATTERN_MATCH\n";
										}
                                    }
                                }               
                            }
                            
                            $SNP_pattern_average = $SNP_pattern_average_tmp/$total_nuc_count;                           
#------------------------------------------------------------------------------------------------------------------------
                            if ($pos_pattern_list_check eq "")
                            {
                                $pos_pattern_list_check = "yes";
								
								print {$filehandle{$seed_id2}} $SNP_pattern_average." SNP_PATTERN_AVERAGE\n";
								
                                my $count_patterns = '0';
								foreach my $count_tmp (sort {$b <=> $a} keys %reads_mismatch2)
								{
									my %pos_pattern_list_tmp;
									undef %pos_pattern_list_tmp;
									my %pos_pattern_list_out_of_range;
									undef %pos_pattern_list_out_of_range;
									
	 #print {$filehandle{$seed_id2}} "\n".$count_tmp." COUNT\n";
									foreach my $rank_tmp (sort {$a <=> $b} keys %{$reads_mismatch2{$count_tmp}})
									{
										my %read_patterns2_tmp;
										undef %read_patterns2_tmp;
										$read_patterns2_tmp{$rank_tmp} = undef;
						#print {$filehandle{$seed_id2}} $rank_tmp." RANK1\n";                
										foreach my $rank_tmp2 (sort {$a <=> $b} keys %{$reads_mismatch2{$count_tmp}{$rank_tmp}})
										{
											$read_patterns2_tmp{$rank_tmp2} = undef;
											
											foreach my $posie_tmp (sort {$a <=> $b} keys %{$reads_mismatch2{$count_tmp}{$rank_tmp}{$rank_tmp2}})
											{
												$pos_pattern_list_tmp{$posie_tmp} = '2';
												 #print {$filehandle{$seed_id2}} $posie_tmp." POSIE\n";
											}
									#print {$filehandle{$seed_id2}} $rank_tmp2." RANK2\n";
HIGHEST_RANK_NP:											
											my $highest_score_rank = '0';
											my $highest_rank_tmp = "";
											my %pos_pattern_list_tmp2;
											undef %pos_pattern_list_tmp2;
											my %pos_pattern_list_tmp3;
											undef %pos_pattern_list_tmp3;
											
											my $highest_pos_score_tmp = '0';
											
											foreach my $posie_tmp (sort {$a <=> $b} keys %pos_pattern_list_tmp)
											{
												if ($pos_pattern_list_tmp{$posie_tmp} > $highest_pos_score_tmp)
												{
													$highest_pos_score_tmp = $pos_pattern_list_tmp{$posie_tmp};
												}
											}
											
											foreach my $rank_tmp3 (sort {$b <=> $a} keys %reads_mismatch)
											{
												my $match_tmp = '0';
												my $match_tmp2 = '0';
												my $no_match_tmp = '0';
												my $count_tmp_new = '0';

												if (exists($read_patterns2_tmp{$rank_tmp3}))
												{
												}
												else
												{
													foreach my $posie_tmp (sort {$a <=> $b} keys %pos_pattern_list_tmp)
													{
														if ($pos_pattern_list_tmp{$posie_tmp} > $highest_pos_score_tmp*0.65)
														{
															if (exists($rank_to_id{$rank_tmp3}))
															{
																my $id_tmpi = $rank_to_id{$rank_tmp3};
																if (exists($alignment_length_save{$id_tmpi}))
																{
																	my $overlap_tmp = $alignment_length_save{$id_tmpi};
																	$count_tmp_new++;
																	if ($posie_tmp > $position-$overlap_tmp)
																	{                                
																		if (exists($store_mismatches_NP{$id_tmpi}{$posie_tmp}))
																		{
																			$match_tmp++;
																			$pos_pattern_list_tmp2{$rank_tmp3}{$posie_tmp} = '1';
																		}
																		elsif (exists($store_mismatches_all_NP{$id_tmpi}{$posie_tmp}))
																		{
																			$pos_pattern_list_tmp2{$rank_tmp3}{$posie_tmp} = '0.5';
																			$match_tmp2++;
																		}
																		else
																		{
																			$no_match_tmp++;
																			$pos_pattern_list_tmp2{$rank_tmp3}{$posie_tmp} = '-0.5';
																		}
																	}
																	else
																	{
																		$pos_pattern_list_tmp3{$rank_tmp3}{$posie_tmp} = '1';
																	}
																}
															}
														}
													}
												}
												
												if ($match_tmp > 0 && $no_match_tmp < ($match_tmp+$match_tmp2)*0.2 && $no_match_tmp < 0.2*$count_tmp_new &&
													(($match_tmp+$match_tmp2) > $count_tmp_new*0.5 || ($match_tmp > 3 && ($match_tmp+$match_tmp2) > $count_tmp_new*0.2) && $highest_pos_score_tmp < 3))
												{
													if ($match_tmp+$match_tmp2 > $highest_score_rank)
													{
														$highest_score_rank = $match_tmp+$match_tmp2;
														$highest_rank_tmp = $rank_tmp3;
													}
												}
											}
											if ($highest_rank_tmp ne "")
											{
												$read_patterns2_tmp{$highest_rank_tmp} = undef;
												foreach my $posie_tmp3 (keys %{$pos_pattern_list_tmp2{$highest_rank_tmp}})
												{
													$pos_pattern_list_tmp{$posie_tmp3} += $pos_pattern_list_tmp2{$highest_rank_tmp}{$posie_tmp3};
													
												}
												foreach my $posie_tmp3 (keys %{$pos_pattern_list_tmp3{$highest_rank_tmp}})
												{
													$pos_pattern_list_out_of_range{$posie_tmp3} += $pos_pattern_list_tmp3{$highest_rank_tmp}{$posie_tmp3};
												}
												
												foreach my $posie_tmp3 (keys %{$reads_mismatch{$highest_rank_tmp}})
												{	
													if (exists($pos_pattern_list_tmp2{$highest_rank_tmp}{$posie_tmp3}))
													{}
													else
													{
														$pos_pattern_list_tmp{$posie_tmp3} += 1;
													}
												}												
												goto HIGHEST_RANK_NP;
											}
											
											my $highest_pos_list = '0';
											foreach my $posie_tmp5 (keys %pos_pattern_list_tmp)
											{
												if ($pos_pattern_list_tmp{$posie_tmp5} > $highest_pos_list)
												{
													$highest_pos_list = $pos_pattern_list_tmp{$posie_tmp5};
												}
											}
											my $total_score_tmp = '0';
											foreach my $posie_tmp5 (keys %pos_pattern_list_tmp)
											{
												if ($pos_pattern_list_tmp{$posie_tmp5} > 1 && ($pos_pattern_list_tmp{$posie_tmp5} > $highest_pos_list*0.5 ||
												   ($pos_pattern_list_tmp{$posie_tmp5} > 3 && $pos_pattern_list_tmp{$posie_tmp5} > $total_nuc_count*0.1))
													&& $pos_pattern_list_tmp{$posie_tmp5}+$pos_pattern_list_out_of_range{$posie_tmp5} > $highest_pos_list*0.7)
												{
													$total_score_tmp += $pos_pattern_list_tmp{$posie_tmp5};
												}
												else
												{
													delete $pos_pattern_list_tmp{$posie_tmp5};
												}
											}
											
											my $final_score_tmp = $total_score_tmp/($total_nuc_count*0.25);
											
											if ($final_score_tmp > 1)
											{
												$count_patterns++;
												my $first_rank_tmp = "";
												foreach my $rank_tmp6 (sort {$a <=> $b} keys %read_patterns2_tmp)
												{
													$read_patterns2{$final_score_tmp}{$count_patterns}{$rank_tmp6} = $count_tmp;
													if ($first_rank_tmp eq "")
													{
														$first_rank_tmp = $rank_tmp6;
													}
												}
												foreach my $pos_tmp6 (sort {$a <=> $b} keys %pos_pattern_list_tmp)
												{
													$pos_pattern_list{$final_score_tmp}{$first_rank_tmp}{$pos_tmp6} = undef;
												}
											}
											
											undef %pos_pattern_list_tmp;
											undef %pos_pattern_list_out_of_range;
											undef %read_patterns2_tmp;
											$read_patterns2_tmp{$rank_tmp} = undef;
										}
									}                                     
								}
								foreach my $score_tmp5 (sort {$b <=> $a} keys %read_patterns2)
								{
PATTERN_NB_NP:                  	foreach my $pattern_number (sort {$a <=> $b} keys %{$read_patterns2{$score_tmp5}})
									{
										my $first = "";
										my $count_ranks = keys %{$read_patterns2{$score_tmp5}{$pattern_number}};
										if (($count_ranks < 3 && $total_nuc_count > 10 ) || $count_ranks < 2)
										{}
										else
										{
											foreach my $rank_tmp7 (sort {$a <=> $b} keys %{$read_patterns2{$score_tmp5}{$pattern_number}})
											{
												$first = $rank_tmp7;
												my $count_pos_tmp = keys %{$reads_mismatch{$rank_tmp7}};
												if ($count_ranks < 3 && $total_nuc_count > 10 && $read_patterns2{$score_tmp5}{$pattern_number}{$rank_tmp7} < 0.5*$count_pos_tmp)
												{
													next PATTERN_NB_NP;
												}
												last;
											}
											foreach my $rank_tmp8 (sort {$a <=> $b} keys %{$read_patterns2{$score_tmp5}{$pattern_number}})
											{
												if ($rank_tmp8 ne $first)
												{
													$read_patterns_final{$score_tmp5}{$first}{$rank_tmp8} = undef;
												}
											}
										}
									}
								}
							}
                            
#-----------------------------------------------------------------------------------------------------------------------------------------------
                           
                            my $highest_avg_score = "";
							my $highest_avg_score2 = "";
                            my $highest_rank_count = "";
                            my $highest_high_score_count = "";
							my $highest_count_matches = "";
                            my $current_score = "";
                            my $current_rank = "";
							my $selected_patterns_count = '0';
                            my %first_no_match;
                            undef %first_no_match;
                            my $score_diff = "";
                            my %highest_avg_score;
                            undef %highest_avg_score;
							my $score_gap = "";
							my %pattern_list;
							undef %pattern_list;
							my %rank_patterns;
							undef %rank_patterns;
							
							foreach my $score_tmp (sort {$b <=> $a} keys %read_patterns_final)
                            {
								foreach my $rank_tmp (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}})
                                {
									if (exists($subject_list{$rank_tmp}))
									{
									}
									else
									{
										my $first_rank_tmp = "";
										foreach my $rank_tmp2 (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}{$rank_tmp}})
										{
											$first_rank_tmp = $rank_tmp2;
										}
										foreach my $rank_tmp3 (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}{$rank_tmp}})
										{
											if ($rank_tmp3 ne $first_rank_tmp)
											{
												$read_patterns_final{$score_tmp}{$first_rank_tmp}{$rank_tmp3} = undef;
											}
										}
										delete $read_patterns_final{$score_tmp}{$rank_tmp};
									}
									foreach my $rank_tmp2 (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}{$rank_tmp}})
									{
										if (exists($subject_list{$rank_tmp2}))
										{
										}
										else
										{
											delete $read_patterns_final{$score_tmp}{$rank_tmp}{$rank_tmp2};
										}
									}  
								}
							}
                            
                            foreach my $score_tmp (sort {$b <=> $a} keys %read_patterns_final)
                            {       
								if ($highest_avg_score eq "")
                                {
                                    $highest_avg_score = $score_tmp;
                                    foreach my $rank_tmp (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}})
                                    {
                                        $highest_avg_score{$rank_tmp} = undef;
                                        foreach my $rank_tmp2 (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}{$rank_tmp}})
                                        {
                                            $highest_avg_score{$rank_tmp2} = undef;
                                        }   
                                    }
                                }
                                elsif ($score_tmp > 0.65*$highest_avg_score)
                                {
                                    foreach my $rank_tmp (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}})
                                    {
                                        if (exists($highest_avg_score{$rank_tmp}))
                                        {}
                                        else
                                        {
                                           $score_diff = "no"; 
                                        }
                                        foreach my $rank_tmp2 (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}{$rank_tmp}})
                                        {
                                            if (exists($highest_avg_score{$rank_tmp2}))
                                            {}
                                            else
                                            {
                                               $score_diff = "no"; 
                                            }
                                        }   
                                    }
                                }
                                
                                foreach my $rank_tmp (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}})
                                {
                                    my $pattern_print = $rank_tmp.";";
                                    my $rank_count_tmp = '0';
                                    my $high_score_count_tmp = '0';
                                    
                                    if ($rank_tmp <= 4)
                                    {
                                        $high_score_count_tmp++;
                                    }
                                    foreach my $rank_tmp2 (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}{$rank_tmp}})
                                    {
                                        $pattern_print .= $rank_tmp2.";";
                                        $rank_count_tmp++;
                                        if ($rank_tmp2 <= 4)
                                        {
                                            $high_score_count_tmp++;
                                        }
                                    }
                                    if ($highest_rank_count eq "")
                                    {
                                        $highest_rank_count = $rank_count_tmp;
                                    }
                                    if ($highest_high_score_count eq "")
                                    {
                                        $highest_high_score_count = $high_score_count_tmp;
                                    }
                                    
#match split pattern with SNP pattern------------------------------------------------------------------                                    
                                    my $count_matches_tmp = '0';
                                    my $count_total_tmp = '0';
                                    my $first_no_match = "";				
                                    my $nuc_match = $split_patterns_final{$rank_tmp};

                                    foreach my $rank_tmp13 (sort {$a <=> $b} keys %split_patterns_final)
                                    {
                                        if ($split_patterns_final{$rank_tmp13} eq $nuc_match)
                                        {
                                            if ($rank_tmp eq $rank_tmp13)
                                            {
                                                $count_matches_tmp++;
                                            }
                                            elsif (exists($read_patterns_final{$score_tmp}{$rank_tmp}{$rank_tmp13}))
                                            {
                                                $count_matches_tmp++;
                                            }
                                            else
                                            {
                                                $count_total_tmp++;
                                                if ($first_no_match eq "")
                                                {
                                                    $first_no_match = $rank_tmp13;
                                                }
                                            }
                                        }
                                        elsif (exists($read_patterns_final{$score_tmp}{$rank_tmp}{$rank_tmp13}))
                                        {
                                            if ($first_no_match eq "")
                                            {
                                                $first_no_match = $rank_tmp13;
                                            }
                                        }
                                    }
                                    if ($first_no_match eq "" && keys %split_patterns_final > 0.8*$total_nuc_count)
                                    {
                                        $first_no_match = $total_nuc_count;
                                    }
									if ($highest_count_matches eq "")
									{
										$highest_count_matches = $count_matches_tmp;
									}
									
									if (exists($pattern_list{$pattern_print}))
									{	
									}
									else
									{
										if ($highest_first_no_match eq "" || $score_tmp > $highest_avg_score*0.85)
										{
											$selected_patterns_count++;
										}
										
										print {$filehandle{$seed_id2}} $pattern_print." FINAL_READ_PATTERN_TEST ".$score_tmp." SCORE ".$first_no_match." FIRST_NO_MATCH\n";
										
										if (($count_matches_tmp > 2 || $score_tmp > $highest_avg_score*0.95) 
											&& (($score_tmp > $highest_avg_score*0.85)
											|| ($score_tmp > $highest_avg_score*0.5 && $score_tmp > 4 && $first_no_match > $highest_first_no_match && $rank_count_tmp > $highest_rank_count && $count_matches_tmp > $highest_count_matches)
											|| ($score_tmp > $highest_avg_score*0.7 && $post_pattern_match_extra eq "yes" && $rank_count_tmp > $highest_rank_count)
											|| ($post_pattern_match eq "yes3" && $score_tmp > 1 && $count_matches_tmp > $highest_count_matches && $first_no_match >= $highest_first_no_match && $score_tmp > $highest_avg_score*0.3))
											&& ($count_matches_tmp > 1 || $highest_first_no_match eq "" || ($score_tmp > 10 && $rank_count_tmp > 3))
											&& ($score_tmp > 3)
											&& ($first_no_match > $highest_first_no_match-3 || $highest_first_no_match eq "" || ($first_no_match >= $highest_first_no_match && $rank_count_tmp > $highest_rank_count)))
										{
											$highest_first_no_match = $first_no_match;
											$first_no_match{$score_tmp}{$rank_tmp}{$first_no_match} = undef;
											$rank_patterns{$score_tmp} = $rank_count_tmp." ".$first_no_match." ".$count_matches_tmp;
											print {$filehandle{$seed_id2}} "SELECT\n";
										}													

										if ($score_gap eq "" && $highest_avg_score2 ne "")
										{
											$score_gap = $highest_avg_score-$score_tmp;
										}
										if ($highest_avg_score2 eq "")
										{
											$highest_avg_score2 = $score_tmp;
										}
										$pattern_list{$pattern_print} = undef;
									}
                                }
                            }
							
							my %check_ranks;
                            undef %check_ranks;
                            foreach my $score_tmp2 (sort {$b <=> $a} keys %first_no_match)
                            {
                                if ($score_tmp2 > 5)
                                {
                                    foreach my $rank_tmp2 (keys %{$first_no_match{$score_tmp2}})
                                    {
                                        if ($rank_tmp2 < 11)
                                        {
                                            $check_ranks{$rank_tmp2} = undef;
                                        }
                                    }
                                }
                            }
                            if (keys %check_ranks > 7)
                            {
                                print {$filehandle{$seed_id2}} "\n\n".length($best_extension)." TERMINATE_CONFLICTING_MISMATCHES\n\n";
                                $best_extension = "";
                                goto AFTER_EXT;
                            }
							
							my $first_rank_tmp = "1";
							my %new_pattern_rank;
							undef %new_pattern_rank;
RANK_PATTERNS_NP:									
							my $high_score_tmp = "0";
							my $high_rank_count_tmp = "0";
							my $high_first_no_match_tmp = "0";
							my $high_count_matches_tmp = "0";
							
							foreach my $score_tmp2 (sort {$b <=> $a} keys %rank_patterns)
                            {
                                my @rank_patterns = split /\s/, $rank_patterns{$score_tmp2};
								if ($high_score_tmp eq "0"
								|| ($score_tmp2 > $high_score_tmp*0.85 && $rank_patterns[2] > $high_count_matches_tmp && $rank_patterns[1] >= $high_first_no_match_tmp)
								|| ($post_pattern_match_extra eq "yes" && $score_tmp2 > $high_score_tmp*0.7 && $rank_patterns[0] >= $high_rank_count_tmp && $rank_patterns[1] >= $high_first_no_match_tmp && $rank_patterns[2] >= $high_count_matches_tmp && ($rank_patterns[0] > $high_rank_count_tmp || $rank_patterns[1] > $high_first_no_match_tmp))
                                || ($post_pattern_match eq "yes3" && $score_tmp2 > $high_score_tmp*0.4 && $rank_patterns[1] > $high_first_no_match_tmp && $rank_patterns[2] > $high_count_matches_tmp)
                                || ($post_pattern_match ne "yes3" && $score_tmp2 > $high_score_tmp*0.8 && $rank_patterns[0] > $high_rank_count_tmp))
								{
									$high_score_tmp = $score_tmp2;
									$high_rank_count_tmp = $rank_patterns[0];
									$high_first_no_match_tmp = $rank_patterns[1];
									$high_count_matches_tmp = $rank_patterns[2];
                                }
                            }
							$new_pattern_rank{$first_rank_tmp}{$high_score_tmp} = undef;
							delete $rank_patterns{$high_score_tmp};
							
							if (keys %rank_patterns > 0)
							{
								$first_rank_tmp++;
								goto RANK_PATTERNS_NP;
							}
                            
                            #$split_patterns_final_score{$nuc_tmp12}{$rank_tmp9} = $score_tmp;
							
READ_PATTERN_FINAL_NP0:    	foreach my $rankie_tmp (sort {$a <=> $b} keys %new_pattern_rank)
                            {
								foreach my $score_tmp (sort {$b <=> $a} keys %{$new_pattern_rank{$rankie_tmp}})
								{
									if (exists($first_no_match{$score_tmp}))
									{
READ_PATTERN_FINAL_NP:              	foreach my $rank_tmp (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}})
										{                       
											if ($highest_first_no_match ne "")
											{
												if (exists($first_no_match{$score_tmp}{$rank_tmp}{$highest_first_no_match}))
												{   
												}
												else
												{
													next READ_PATTERN_FINAL_NP;
												}
											}
											my $pattern_print = $rank_tmp.";";
											my $count_ranks_tmp =  keys %{$read_patterns_final{$score_tmp}{$rank_tmp}};
																					
											foreach my $rank_tmp2 (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}{$rank_tmp}})
											{
												$pattern_print .= $rank_tmp2.";";
											}
		
											$reads_to_remove{$rank_tmp}{$rank_tmp} = $score_tmp;
											foreach my $rank_tmp2 (sort {$a <=> $b} keys %{$read_patterns_final{$score_tmp}{$rank_tmp}})
											{
												$reads_to_remove{$rank_tmp}{$rank_tmp2} = $score_tmp;
												$remove_reads = "yes";
												if (exists($nucs_for_split_extra{$split_patterns_final{$rank_tmp}}))
												{
													my $count_tmp = $nucs_for_split_extra{$split_patterns_final{$rank_tmp}};
													$nucs_for_split_extra{$split_patterns_final{$rank_tmp}} = $count_tmp+1;
												}
												else
												{
													$nucs_for_split_extra{$split_patterns_final{$rank_tmp}} = '1';
												}
											}
											print {$filehandle{$seed_id2}} $pattern_print." FINAL_READ_PATTERN ".$score_tmp." SCORE\n";
											
											if ($remove_reads eq "yes")
											{
												$current_score = $score_tmp;
												$current_rank = $rank_tmp;	
												$remove_reads = "";
												my $pos_count_tmp = keys %{$pos_pattern_list{$current_score}{$current_rank}};
												foreach my $rank_tmp0 (keys %reads_to_remove)
												{
													my $count_tmp = keys %{$reads_to_remove{$rank_tmp0}};
													my $SNP_pattern_average_remove = "";
													my $SNP_pattern_average_remove_tmp = '0';
													my $SNP_pattern_average_count = "";
													my $SNP_pattern_average_count_tmp = '0';
													my $cover_complete_assembly = "";
													
													if ($count_tmp > 0)
													{
														my $high_score_count = '0';
														my $high_score_count2 = '0';
														my $rank_one_check = "";
														my $score_diff_check = "";
														my $rank_tmp2 = "";
														
														foreach my $rank_tmp (sort {$a <=> $b} keys %{$reads_to_remove{$rank_tmp0}})
														{
															if ($rank_tmp <= 4)
															{
																$high_score_count++;
															}
															if ($rank_tmp <= 6)
															{
																$high_score_count2++;
															}
															if ($rank_tmp eq '1')
															{
																$rank_one_check = "yes";
																
																$rank_tmp2 = $rank_tmp;
																while (exists($reads_to_remove{$rank_tmp0}{$rank_tmp2}))
																{
																	$rank_tmp2++;
																}
																my $id_tmp01 = $rank_to_id{'1'};
																my $id_tmp02 = $rank_to_id{$rank_tmp2};
																my $score_tmp1 = "";
																my $score_tmp2 = "";
																foreach my $score_tmp (sort {$b <=> $a} keys %scores2)
																{
																	my @ids_tmp = split /,/, $scores2{$score_tmp};
																	foreach my $ids_tmp (@ids_tmp)
																	{
																		if ($ids_tmp eq $id_tmp01)
																		{
																			$score_tmp1 = $score_tmp;
																		}
																		if ($ids_tmp eq $id_tmp02)
																		{
																			$score_tmp2 = $score_tmp;
																		} 
																	}
																}
																if ($score_tmp2 > $score_tmp1*0.7)
																{
																	$score_diff_check = "yes";
																}
															}
	
															my $id_tmp = $rank_to_id{$rank_tmp};
															my $length_tmp = $alignment_length_save{$id_tmp};
															$SNP_pattern_average_remove_tmp += $pos_list_number{$rank_tmp}/$length_tmp;
															$SNP_pattern_average_count_tmp += $pos_list_number{$rank_tmp};										
				
															if (exists($rank_to_id{$rank_tmp}))
															{
																my $id_tmp = $rank_to_id{$rank_tmp};
																if (exists($extensions2{$id_tmp}))
																{
																	if ($alignment_length_save{$id_tmp} > length($read)-500)
																	{
																		$cover_complete_assembly = "yes";
																	}
																}
															}
														}
														
														$SNP_pattern_average_remove = $SNP_pattern_average_remove_tmp/$count_tmp;
														$SNP_pattern_average_count = $SNP_pattern_average_count_tmp/$count_tmp;
														
														print {$filehandle{$seed_id2}} $SNP_pattern_average_remove." SNP_PATTERN_AVERAGE_REMOVE\n";
														
														my $count_matches_tmp = '0';
														my $count_pos_no_matches_tmp = '0';
														my $count_no_matches_tmp = '0';
														my $first_no_match = "";
														my $first_no_match2 = '0';
														my $count_total_tmp = '0';
														my $nuc_match = $split_patterns_final{$rank_tmp0};
														
														foreach my $rank_tmp13 (sort {$a <=> $b} keys %split_patterns_final)
														{
															$count_total_tmp++;
															if ($split_patterns_final{$rank_tmp13} eq $nuc_match)
															{                                   
																if (exists($reads_to_remove{$rank_tmp0}{$rank_tmp13}))
																{
																	$count_matches_tmp++;
																	if ($first_no_match eq "")
																	{
																		$first_no_match2++;
																	}
																}
																else
																{
																	my $no_match = "";
																	my $no_match_tmp = '0';
																	my $match_tmp = '0';
																	if (exists($rank_to_id{$rank_tmp13}))
																	{
																		my $id_tmpi = $rank_to_id{$rank_tmp13};
																		if (exists($alignment_length_save{$id_tmpi}))
																		{
																			my $overlap_tmp = $alignment_length_save{$id_tmpi};
																			
																			foreach my $pos_snp_tmp (sort {$a <=> $b} keys %{$pos_pattern_list{$current_score}{$current_rank}})
																			{
																				if ($position-$pos_snp_tmp < $overlap_tmp)
																				{
																					if (exists($store_mismatches_all_NP{$id_tmpi}{$pos_snp_tmp}))
																					{
																						$match_tmp++;
																					}
																					elsif (exists($store_mismatches_NP{$id_tmpi}{$pos_snp_tmp}))
																					{
																						$match_tmp++;
																					}
																					else
																					{
																						$no_match_tmp++;
																					}
																				}
																			} 
																		}
																		elsif ($first_no_match eq "")
																		{
																			$first_no_match = $rank_tmp13;
																		}
																	}
																	if ($no_match_tmp > 0 && $no_match_tmp > ($match_tmp+$no_match_tmp)*0.3)
																	{
																		$no_match = "yes";
																	}
																	
																	if ($no_match eq "" && $post_pattern_match ne "" && $count_matches_tmp > 2 && $count_no_matches_tmp < 3)
																	{
																		$count_matches_tmp++;
																	}
																	elsif ($no_match ne "")
																	{
																		$count_no_matches_tmp++;
																	}
																	if ($first_no_match eq "" && ($no_match ne "" || $match_tmp < 3))
																	{
																		$first_no_match = $rank_tmp13;
																	}
																	elsif ($first_no_match eq "")
																	{
																		$first_no_match2++;
																	}
																}
															}
															elsif (exists($reads_to_remove{$rank_tmp0}{$rank_tmp13}))
															{
																if ($first_no_match eq "")
																{
																	$first_no_match = $rank_tmp13;
																}
																$count_no_matches_tmp++;
															}
															else
															{          
																if (exists($rank_to_id{$rank_tmp13}))
																{
																	my $id_tmpi = $rank_to_id{$rank_tmp13};
																	if (exists($alignment_length_save{$id_tmpi}))
																	{
																		my $overlap_tmp = $alignment_length_save{$id_tmpi};
																		my $no_pos_match_tmp = '0';
																		my $pos_match_tmp = '0';
																		
																		foreach my $pos_snp_tmp (sort {$a <=> $b} keys %{$pos_pattern_list{$current_score}{$current_rank}})
																		{
																			if ($position-$pos_snp_tmp < $overlap_tmp)
																			{
																				$pos_match_tmp++;
																			}
																			else
																			{
																				$no_pos_match_tmp++;
																			}
																		}
																		if ($pos_match_tmp > 0 && $pos_match_tmp > 0.5*($pos_match_tmp+$no_pos_match_tmp))
																		{
																			$count_pos_no_matches_tmp++;
																		}           
																	}
																}
															}
														}
														if ($first_no_match eq "" && keys %split_patterns_final > 0.8*$total_nuc_count)
														{
															$first_no_match = $total_nuc_count;
															if ($first_no_match < 10)
															{
																$first_no_match = '20';
															}
														}
														
														if ($count_pos_no_matches_tmp > 1 && $count_matches_tmp > 2 && $count_matches_tmp > 0.15*$count_total_tmp && $count_no_matches_tmp < 0.15*$count_matches_tmp
															&& ($first_no_match > 10 || ($first_no_match > 5 && $count_matches_tmp > 5)))
														{
															print {$filehandle{$seed_id2}} $count_matches_tmp." ADD_ALL_SNP_PATTERN\n";
															foreach my $rank_tmp13 (sort {$a <=> $b} keys %split_patterns_final)
															{
																if ($split_patterns_final{$rank_tmp13} eq $nuc_match)
																{
																	if (exists($reads_to_remove{$rank_tmp0}{$rank_tmp13}))
																	{                                                     
																	}
																	else
																	{
																		$reads_to_remove{$rank_tmp0}{$rank_tmp13} = undef;
																		if ($rank_tmp13 <= 4)
																		{
																			$high_score_count++;
																		}
																		if ($rank_tmp13 <= 6)
																		{
																			$high_score_count2++;
																		}
																	}
																}
															}
															if ($highest_first_no_match < 10)
															{
																$highest_first_no_match = '10';
															}
														}										
														
														if ($count_matches_tmp > 4 && $count_matches_tmp > 0.9*($count_no_matches_tmp+$count_matches_tmp))
														{
															if ($first_no_match < 10)
															{
																$first_no_match = '10';
																print {$filehandle{$seed_id2}} $first_no_match." CORRECT_FIRST_NO_MATCH\n";
															}
														}
														
														print {$filehandle{$seed_id2}} $count_no_matches_tmp." COUNT_NO_MATCHES\n";
														print {$filehandle{$seed_id2}} $count_matches_tmp." COUNT_MATCHES\n";
														if (($rank_one_check eq "yes" && $score_diff_check eq "" && $first_no_match < 10)
															|| ($high_score_count > 1 && $current_score < 7 && ($first_no_match < 10 || ($post_pattern_match ne "yes2" && $post_pattern_match ne "yes3") || $post_pattern_match_average eq ""))
															|| ($first_no_match < 10 && $total_nuc_count < 9)
															|| ($NP_reads_support eq "yes2" && ($current_score < 10 || $count_tmp < 4 || $first_no_match < 10 || $count_matches_tmp < 3))
															|| ($post_pattern_match eq "" && ($count_tmp < 4 || $current_score < 4) && ($current_score < 20 || $count_tmp < 3))
															|| ($post_pattern_match ne "yes3" && $count_tmp < 4 && ($current_score < 10 || $count_tmp < 3))
															|| ($first_no_match < 10 && $current_score < 3 && $count_tmp < 3)
															|| ($post_pattern_match_extra ne "yes" && ($count_tmp < 3 || $high_score_count > 1))
															|| ($first_no_match < 7 && $current_score < 4)
															|| ($first_no_match < 10 && $current_score < 40 && $post_pattern_match_extra ne "yes" && $total_nuc_count-$count_tmp < 11)
															|| ($post_pattern_match_extra ne "yes" && $current_score < 5)		
															|| ($count_matches_tmp/$count_tmp < 0.68 && $current_score < 15)
															|| (($first_no_match < 10 || $count_tmp < 3) && $score_gap < $current_score*0.5 && ($count_tmp < 5 || $current_score < 15) && $post_pattern_match_extra ne "yes")
															|| ($first_no_match < 5 && ($current_score < 10 || ($post_pattern_match_extra eq "" && $current_score < 30)))
															|| (($first_no_match < 10 || $count_matches_tmp < 3) && $current_score < 50 && $count_tmp < 5)
															|| ($first_no_match < 8 && $current_score < 20 && $count_tmp < 5 && $selected_patterns_count > 1)
															|| (($first_no_match < 7 || $first_no_match2 < 3) && $current_score < 5)
															|| ($post_pattern_match_extra ne "yes" && ($SNP_check eq "" || ($post_pattern_match_average ne "yes" && $count_tmp < 6)) && $current_score < 5)
															|| ($post_pattern_match_extra ne "yes" && ($SNP_check eq "" || ($count_tmp < 5 && $total_nuc_count > 15)) && $current_score < 10 && $first_no_match < 10)
															|| ($post_pattern_match_extra ne "yes" && $current_score < 20 && $count_tmp < 4 && ($cover_complete_assembly eq "yes" || $post_pattern_match_average ne "yes"))
															|| (($post_pattern_match_extra ne "yes" || $first_no_match < 10) && ($count_tmp < 3 || $rank_one_check eq "yes") && $current_score < 25)
															|| ($high_score_count2 > 4 && ($current_score < 50 || $total_nuc_count < 11 || $count_matches_with_high_scores < $sequencing_depth_NP))
															|| ($high_score_count > 1 && $count_tmp < 4 && $current_score < 11)
															|| ($SNP_pattern_average > 0.0005 && $SNP_pattern_average_remove > 0.0005 && $current_score < 3 && $count_tmp < 5)
															|| ($high_score_count2 > 0 && $count_tmp < 4 && $current_score < 15 && $first_no_match < 10)
															|| ($high_score_count > 1 && $SNP_check ne "yes2" && $current_score < 15 && $post_pattern_match ne "yes3")
															|| ($high_score_count > 0 && $current_score < 30 && $count_tmp < 4 && $total_nuc_count > 10 && $post_pattern_match_extra ne "yes")
															|| ($high_score_count > 0 && ($count_tmp < 3 || ($post_pattern_match_extra ne "yes" && ($first_no_match < 10 || $count_tmp < 5))) && $current_score < 10)
															|| ($high_score_count > 0 && ($first_no_match < 10 || $first_no_match2 < 3) && $current_score < 5)
															|| ($high_score_count > 2 && $first_no_match < 10 && $current_score < 15)
															|| ($high_score_count > 2 && $count_pos_no_matches_tmp < 2 && ($post_pattern_match_extra eq "yes" || $current_score < 10) && $post_pattern_match ne "yes3")
															|| ($first_no_match < 6 && $count_tmp < 4 && $post_pattern_match_extra ne "yes")
															|| ($high_score_count2 > 3 && $high_score_count > 2 && $current_score < 15 && $post_pattern_match ne "yes3")
															|| ($count_matches_tmp < 3 && $current_score < 5 && $count_tmp < 4 && ($post_pattern_match eq "yes" || $post_pattern_match eq "" || $post_pattern_match_extra ne "yes"))
															|| ($first_no_match < 10 && $post_pattern_match_extra ne "yes" && $SNP_check eq "" && (($high_score_count > 1 && $current_score < 10) || $rank_one_check eq "yes"))
															|| ($first_no_match < 10 && $count_matches_tmp < 3 && $current_score < 4)
															|| ($count_tmp < 4 && $total_nuc_count < 10 && ($first_no_match < 10 || $post_pattern_match eq "" || $post_pattern_match eq "yes" || $post_pattern_match_average eq ""))
															|| ($first_no_match < 9 && $count_matches_tmp < 3 && $current_score < 5 && $total_nuc_count < 13)
															|| ($SNP_pattern_average_remove < $SNP_pattern_average*1.6 && (($first_no_match < 10 && $current_score < 5) || ($count_tmp < 3 && ($post_pattern_match eq "" || $post_pattern_match eq "yes" || $first_no_match < 20))))
															|| ($first_no_match < 8 && $score_diff eq "no" && $current_score < 10 && ($count_tmp < 4 || $count_tmp < $total_nuc_count*0.15))
															|| ($SNP_pattern_average_remove < $SNP_pattern_average*1.6 && $high_score_count > 2 && $high_score_count2 > 3 && $pos_count_tmp < $SNP_pattern_average_count*0.3))
														{     
															undef %reads_to_remove;
															#delete $read_patterns_final{$current_score}{$current_rank};
															print {$filehandle{$seed_id2}} $high_score_count." SCORE ".$score_diff_check." CANCEL_REMOVE\n";
															$highest_first_no_match = "";
															next READ_PATTERN_FINAL_NP0;
														}
														else
														{
															$remove_reads = "yes";
															$high_score_save = $current_score;
															
#Exclude reads------------------------                                              
															if ($count_matches_with_high_scores > $sequencing_depth_NP && ($current_score > 20 || $first_no_match > 9))
															{
																my $id_tmpi = "";
																if (exists($rank_to_id{$rank_tmp0}))
																{
																	$id_tmpi = $rank_to_id{$rank_tmp0};
																}
																if ($hap_tag eq "HAP1" && $repetitive_detect1 eq "" && $repetitive_detect2 eq "")
																{
																	$exclude_reads_hap1_NP{$id_tmpi} = $position+$alignment_length_save{$id_tmpi};
																}
																elsif ($hap_tag eq "HAP2" && $repetitive_detect1 eq "" && $repetitive_detect2 eq "")
																{
																	$exclude_reads_hap2_NP{$id_tmpi} = $position+$alignment_length_save{$id_tmpi};
																}
																foreach my $rank_tmp5 (sort {$a <=> $b} keys %{$reads_to_remove{$rank_tmp0}})
																{
																	my $id_tmpi2 = "";
																	if (exists($rank_to_id{$rank_tmp5}))
																	{
																		$id_tmpi2 = $rank_to_id{$rank_tmp5};
																	}
																	if ($hap_tag eq "HAP1" && $repetitive_detect1 eq "" && $repetitive_detect2 eq "")
																	{
																		$exclude_reads_hap1_NP{$id_tmpi2} = $position+$alignment_length_save{$id_tmpi};
																	}
																	elsif ($hap_tag eq "HAP2" && $repetitive_detect1 eq "" && $repetitive_detect2 eq "")
																	{
																		$exclude_reads_hap2_NP{$id_tmpi2} = $position+$alignment_length_save{$id_tmpi};
																	}
																}
															}
															last READ_PATTERN_FINAL_NP0;
														}
													}
												}
											}
                                        }
                                    }                            
                                }  
                            }
#-------------------                            
                            my $count_pattern_lists = keys %reads_mismatch;
                            foreach my $nuc_tmp20 (keys %SNP_patterns_prev2)
                            {
                                if ($SNP_patterns_prev2{$nuc_tmp20} > $count_pattern_lists*0.25 && $SNP_patterns_prev2{$nuc_tmp20} > 1)
                                {
                                    print {$filehandle{$seed_id2}} $nuc_tmp20." ".$SNP_patterns_prev2{$nuc_tmp20}." PATTERN_MATCH_FINAL\n";
                                }
                            }
                        }
						
						if (($post_pattern_match eq "yes2" || $post_pattern_match eq "yes3" || $SNP_check eq "yes2") &&
                            $post_pattern_match_average eq "yes" && $post_pattern_match_extra eq "yes" && $remove_reads eq "" && $find_haps_in_seed eq "")
                        {
                            if ($assembly_length_max eq "WG" && $y eq "1" && ($first_back_assembly eq "" || length($read) < 5000))
							{
								$first_back_assembly = "yes";
								goto END1;
							}
						}
                        
#Select based on average score of each group---------------------------------------------------------------------------                                           
                        
                        if (((($post_pattern_match eq "yes2" || $post_pattern_match eq "yes3") && $post_pattern_match_average eq "yes" && $post_pattern_match_extra eq "yes")
							 || ($SNP_check eq "yes2" && $total_nuc_count > 13)) && $remove_reads eq "" && $find_haps_in_seed eq "" && $add_no_match_reads eq "" && $add_rejected_reads eq "")
                        {
							my %average_rank_score_first;
                            undef %average_rank_score_first;
                            my %average_rank_score_firstb;
                            undef %average_rank_score_firstb;
                            my $highest_nuc_tmp = "";
                            my $hihest_score_tmp = '0';
                            
                            foreach my $nuc_tmp (sort {$a <=> $b} keys %split_patterns_final_score)
                            {
                                my $rank_score = '0';
                                my $rank_count = '0';

                                foreach my $rank_tmp (sort {$a <=> $b} keys %{$split_patterns_final_score{$nuc_tmp}})
                                {
                                    my $score_tmp = '0';
                                    if (exists($rank_to_id{$rank_tmp}))
                                    {
                                        my $id_tmp = $rank_to_id{$rank_tmp};
                                        if (exists($scores{$id_tmp}))
                                        {
                                            $score_tmp = $scores{$id_tmp};
                                        }
                                    }
                                    if ($rank_score eq '0')
                                    {
                                        $average_rank_score_first{$nuc_tmp} = $score_tmp;
                                        $average_rank_score_firstb{$nuc_tmp} = $rank_tmp;
                                    }
                                    $rank_score += $score_tmp;
                                    $rank_count++;
                                }
                                my $average_rank = $rank_score/$rank_count;
                                if ($rank_count > 1 && $rank_count > $total_nuc_count*0.19)
                                {
                                    $average_rank_score{$nuc_tmp} = $average_rank;
                                    print {$filehandle{$seed_id2}} $nuc_tmp." ".$average_rank." AVERAGE_RANK\n";
                                    if ($average_rank > $hihest_score_tmp && $highest_nuc_tmp ne "-")
                                    {
                                        $hihest_score_tmp = $average_rank;
                                        $highest_nuc_tmp = $nuc_tmp;
                                    }
                                } 
                            }                          
                            if (exists($average_rank_score{$first_rank}))
                            {         
                                my $nuc_best = "";
                                foreach my $nuc_tmp (keys %average_rank_score)
                                {
                                    my $count_ranks = keys %{$SNP_patterns_now{$nuc_tmp}};
                                    if (($count_ranks > 2 && ($post_pattern_match eq "yes2" || $post_pattern_match eq "yes3" || $count_ranks > $total_nuc_count*0.3)) && $nuc_tmp ne $first_rank &&
                                        ($average_rank_score{$first_rank}*0.4 > $average_rank_score{$nuc_tmp} || ($average_rank_score_first{$first_rank} > 2*$average_rank_score_first{$nuc_tmp} && $average_rank_score_firstb{$nuc_tmp} > 4)
                                       || ($post_pattern_match eq "yes3" && $average_rank_score{$first_rank}*0.7 > $average_rank_score{$nuc_tmp}
                                        && $average_rank_score_firstb{$nuc_tmp}/$total_nuc_count > 0.3 && $average_rank_score_firstb{$nuc_tmp} > 3)
                                    || (($post_pattern_match_average ne "yes" || $post_pattern_match_extra ne "yes") && $average_rank_score{$first_rank}*0.65 > $average_rank_score{$nuc_tmp}
                                        && $average_rank_score_firstb{$nuc_tmp}/$total_nuc_count > 0.3)))
                                    {
                                        if (($post_pattern_match_extra eq "yes" && $post_pattern_match_average eq "yes") || ($post_pattern_match_extra eq "yes" && $SNP_check eq "yes2" && $total_nuc_count > 13))
                                        {
                                            my $count_tmpi = '0';
											print {$filehandle{$seed_id2}} $nuc_tmp." ".$average_rank_score{$nuc_tmp}." ".$average_rank_score{$first_rank}." AVERAGE_RANK_REMOVE\n";
                                            foreach my $rank_tmp (sort {$a <=> $b} keys %{$SNP_patterns_now{$nuc_tmp}})
                                            {
                                                if ($rank_tmp < 3)
                                                {
                                                    last;
                                                }
                                                if ($rank_tmp < 7)
                                                {
                                                    $count_tmpi++;
                                                    if ($count_tmpi > 1)
                                                    {
                                                        $remove_reads = "";
                                                        undef %reads_to_remove;
                                                        last;
                                                    }
                                                }
                                                $selected_nuc = $first_rank;
                                                $reads_to_remove{$rank_tmp}{$rank_tmp} = undef;
                                                $remove_reads = "yes";
                                            }
                                        }
										else
										{
											my $count_tmpi = '0';
                                            foreach my $rank_tmp (sort {$a <=> $b} keys %{$SNP_patterns_now{$nuc_tmp}})
                                            {
                                                if ($rank_tmp < 3)
                                                {
                                                    $nuc_best = "no";
                                                }
                                                if ($rank_tmp < 7)
                                                {
                                                    $count_tmpi++;
                                                    if ($count_tmpi > 1)
                                                    {
                                                        $nuc_best = "no";
                                                    }
                                                }
                                            }
										}
                                    }
                                    else
                                    {
                                        if ($nuc_best eq "")
                                        {
                                            $nuc_best = $nuc_tmp;
                                        }
                                        else
                                        {
                                            $nuc_best = "no";
                                        }
                                    }
                                }
                     
                                if ($remove_reads eq "yes")
                                {
                                    print {$filehandle{$seed_id2}} "REMOVE_BOTTOM_SCORES\n";
                                    goto REMOVE_READS_NP;
                                }
                            }

                            my $nuc_highest_tmp = "";
                            my $count_highest_tmp = '0';
                            my $count_total_tmp = '0';
                            my $count_total_ext2 = '0';
                            my %nuc_tmp2;
                            undef %nuc_tmp2;
    
                            foreach my $ranki_tmp (sort {$a <=> $b} keys %subject_list)
                            {
                                if (exists($rank_to_id{$ranki_tmp}))
                                {
                                    my $id_tmp = $rank_to_id{$ranki_tmp};
                                    if (exists($extensions2{$id_tmp}))
                                    {
                                        if ($alignment_length_save{$id_tmp} > length($read)-500)
                                        {
											my $nuc_tmp = $nucs_by_rank{$ranki_tmp};
                                            $nuc_tmp2{$nuc_tmp} += 1;
                                            if ($nuc_tmp2{$nuc_tmp} > $count_highest_tmp)
                                            {
                                                $count_highest_tmp = $nuc_tmp2{$nuc_tmp};
                                                $nuc_highest_tmp = $nuc_tmp;
                                            }
                                            $count_total_tmp++;
                                        }
                                    }
                                    $count_total_ext2++;
                                }
                            }
    
                            if ($nuc_highest_tmp ne "" && $count_highest_tmp > $total_nuc_count*0.1 && $count_highest_tmp > 1
                                && $count_highest_tmp eq $count_total_tmp && $nucs{$nuc_highest_tmp} > $total_nuc_count*0.3)
                            {
                                if ($post_pattern_match_average eq "yes" && $post_pattern_match_extra eq "yes")
                                {
                                    my $count_tmpi = '0';
                                    foreach my $nuc_tmp (keys %SNP_patterns_now)
                                    {
                                        if ($nuc_tmp ne $nuc_highest_tmp && $nuc_tmp ne "-" && $nucs{$nuc_tmp} > $total_nuc_count*0.2)
                                        {
                                            foreach my $rank_tmp (sort {$a <=> $b} keys %{$SNP_patterns_now{$nuc_tmp}})
                                            {
                                                if ($rank_tmp < 3)
                                                {
                                                    $remove_reads = "";
                                                    undef %reads_to_remove;
													last;
                                                }
                                                if ($rank_tmp < 6)
                                                {
                                                    $count_tmpi++;
                                                    if ($count_tmpi > 2)
                                                    {
                                                        $remove_reads = "";
                                                        undef %reads_to_remove;
                                                        last;
                                                    }
                                                }
                                                $selected_nuc = $first_rank;
                                                $reads_to_remove{$rank_tmp}{$rank_tmp} = undef;
                                                $remove_reads = "yes";
                                            }
                                        }
                                    }
                                    if ($remove_reads eq "yes")
                                    {
                                        print {$filehandle{$seed_id2}} "REMOVE_BOTTOM_SCORES2\n";
                                        goto REMOVE_READS_NP;
                                    }
                                }
                                else
                                {
                                    $nuc_match = $nuc_highest_tmp;
                                    $nuc_match =~ tr/actgn/ACTGN/;
                                    if ($nuc_highest_tmp eq "-")
                                    {
                                        $nuc_match = "";
                                    }
                                    $best_extension .= $nuc_match;
									$best_extension_part .= $nuc_match;
                                    
                                    if ($nuc_highest_tmp eq "-" && $nucs{"-"} > 0)
                                    {
										if ($nucs{"-"}/$count_total_ext2 < 0.8)
										{
											$quality_scores_gap_tmp{length($best_extension)} = $nucs{"-"}/$count_total_ext2." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
										}
                                    }
                                    else
                                    {
                                        $quality_scores_tmp{length($best_extension)} = $nucs{$nuc_highest_tmp}/$count_total_ext2." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
                                    }           
                                    
                                    print {$filehandle{$seed_id2}} $nuc_match." N_CORRECTION_LONG\n";
                                    goto SKIP_INPUT_BLAST3_NP
                                }
                            }
                        }
        
#Add no match reads to mafft-------------------------------------------------------------------------------------------------------------------------------                                        
                        
						my $count_ext2_counts = '0';
						if ($ext2_count > 0)
						{	
							foreach my $nuc_tmp (sort {$a <=> $b} keys %split_patterns_final_score)
                            {
								my $rank_count_tmp = keys %{$split_patterns_final_score{$nuc_tmp}};
								if ($rank_count_tmp > $total_nuc_count*0.2)
								{
									foreach my $rank_tmp (sort {$a <=> $b} keys %{$split_patterns_final_score{$nuc_tmp}})
									{
										if (exists($rank_to_id{$rank_tmp}))
										{
											my $id_tmp3 = $rank_to_id{$rank_tmp};
											if (exists($extensions2{$id_tmp3}))
											{
												$count_ext2_counts++;
												last;
											}
										}
									}
								}
							}
						}
						
						if ($add_no_match_reads eq "" && $add_rejected_reads eq "" && ($post_pattern_match eq "yes2" || $post_pattern_match eq "yes3" || $SNP_check eq "yes2") && $SNR_read_ahead eq "" && $post_pattern_match_extra eq "yes"
                            && $remove_reads eq "" && $extensions_nomatch2b_count+$extensions_nomatch2b_count_saved > 1 && $find_haps_in_seed eq "" && $count_ext2_counts < 2)
                        {          
							my %extensions_nomatch2b_tmp;
							undef %extensions_nomatch2b_tmp;
							my $score_score_no_match = '0';
							my $score_score_no_match1 = '0';
							
							foreach my $ext_id_tmp (keys %extensions_nomatch2b)
							{
								foreach my $ext_tmp (keys %{$extensions_nomatch2b{$ext_id_tmp}})
								{
									foreach my $score_no_match_tmp (keys %{$extensions_nomatch2b{$ext_id_tmp}{$ext_tmp}})
									{
										if ($score_no_match_tmp > 0 && $extensions_nomatch2b{$ext_id_tmp}{$ext_tmp}{$score_no_match_tmp} < $score_no_match_tmp*0.8 && length($ext_tmp) > length($best_extension))
										{
											if ($score_no_match_tmp > 1)
											{
												$score_score_no_match++;
											}
											elsif ($score_no_match_tmp eq 1)
											{
												$score_score_no_match1++;
											}
											$extensions_nomatch2b_tmp{$ext_id_tmp} = undef;
											print {$filehandle{$seed_id2}} $ext_id_tmp." ADD_NOMATCH\n";
										}
									}
								}
							}
							
							foreach my $ext_id_tmp (keys %extensions_nomatch2b_saved)
							{
								foreach my $score_no_match_tmp (keys %{$extensions_nomatch2b_saved{$ext_id_tmp}})
								{
									if ($score_no_match_tmp > 0)
									{
										foreach my $score_match_tmp (keys %{$extensions_nomatch2b_saved{$ext_id_tmp}{$score_no_match_tmp}})
										{
											if ($score_no_match_tmp*0.3 > $score_match_tmp && $extensions_nomatch2b_saved{$ext_id_tmp}{$score_no_match_tmp}{$score_match_tmp} > $position+length($best_extension))
											{
												if ($score_no_match_tmp > 1)
												{
													$score_score_no_match++;
												}
												elsif ($score_no_match_tmp eq 1)
												{
													$score_score_no_match1++;
												}   
											}
											$extensions_nomatch2b_tmp{$ext_id_tmp} = undef;
											print {$filehandle{$seed_id2}} $ext_id_tmp." ADD_NOMATCH_SAVED\n";
										}
									}
								}
							}
						
							print {$filehandle{$seed_id2}} $score_score_no_match." SCORE_NO_SCORE_MATCH\n";   
							if ((($score_score_no_match > 0 && $score_score_no_match1+$score_score_no_match > 0) || $score_score_no_match1 > 2) && $score_score_no_match < 20)
							{
								foreach my $id_tmp55 (keys %extensions_nomatch2b_tmp)
								{
									undef %id_matches;
									if (exists($extensions_nomatch2b_saved{$id_tmp55}))
									{
										$id_matches{$id_tmp55} = undef;
									}
									if (exists($extensions_nomatch2b{$id_tmp55}))
									{
										if (exists($scores2{'0'}))
										{
											$scores2{'0'} .= ",$id_tmp55"; 
										}
										else
										{
											$scores2{'0'} = $id_tmp55;
										}
										foreach my $ext_tmp (keys %{$extensions_nomatch2b{$id_tmp55}})
										{
											$extensions2_tmp{$id_tmp55} = $ext_tmp; 
										}
									}
								}
								foreach my $pos_tmp14 (keys %{$split_positions_DUP{$id}})
								{
									if ($pos_tmp14 > $position)
									{
										 delete $split_positions_DUP{$id}{$pos_tmp14};
									}
								}
								$add_no_match_reads = $total_nuc_count_original;
								print {$filehandle{$seed_id2}} $best_extension."\nHALLE0\n";
								$mismatch_retry++;
								$best_extension = "";
								$best_extension_part = "";
								$remove_reads_check = "";
								undef %quality_scores_tmp;
								undef %span_complex_region;
								
								if (keys %id_matches > 0)
								{
									goto ADD_REJ_POS_NP;
								}
								goto SELECT_LENGTH_NP2;
							}
                        }
                        elsif ($add_no_match_reads ne "" && ($post_pattern_match eq "yes3" || $SNP_check eq "yes2"))
                        {                            
                            my $score_score_match = '0';
                            my $score_score_match1 = '0';
                            my %score_matches_tmp;
                            undef %score_matches_tmp;
                            foreach my $ext_id_tmp (keys %extensions2b)
                            {
                                foreach my $ext_tmp (keys %{$extensions2b{$ext_id_tmp}})
                                {
                                    foreach my $score_match_tmp (keys %{$extensions2b{$ext_id_tmp}{$ext_tmp}})
                                    {
                                        my $rank_tmp0 = $id_to_rank{$ext_id_tmp};
                                        if (exists($split_patterns_final{$rank_tmp0}))
                                        {
                                            if ($score_match_tmp > 1)
                                            {
                                                $score_score_match++;
                                                $score_matches_tmp{$split_patterns_final{$rank_tmp0}} += 2;
                                            }
                                            elsif ($score_match_tmp eq 1)
                                            {
                                                $score_score_match1++;
                                                $score_matches_tmp{$split_patterns_final{$rank_tmp0}} += 1;
                                            }   
                                        }
                                    }
                                }
                            }
                            print {$filehandle{$seed_id2}} $score_score_match." SCORE_SCORE_MATCH\n";   
                            if (($score_score_match > 0 && $score_score_match1+$score_score_match > 0) || $score_score_match1 > 2)
                            {
                                foreach my $nuc_rej_tmp (keys %nucs_rej)
                                {
                                    if (exists($score_matches_tmp{$nuc_rej_tmp}))
                                    {
                                        if ($score_matches_tmp{$nuc_rej_tmp} > 1 && $nucs_rej{$nuc_rej_tmp} >= 0.1*$total_nuc_count_rej)
                                        {
                                            goto REMOVE_READS_NP;
                                        }
                                    }
                                }
                                foreach my $nuc_rej_tmp (keys %nucs_rej)
                                {
                                    if ($nucs_rej{$nuc_rej_tmp} > 1)
                                    {
                                        my $count_removed = '0';
                                        foreach my $rank_tmp8 (sort {$a <=> $b} keys %split_patterns_final)
                                        {
                                            if ($split_patterns_final{$rank_tmp8} eq $nuc_rej_tmp)
                                            {
                                                $count_removed++;
                                            }
                                        }
                                        if ($total_nuc_count-$count_removed < 6)
                                        {
                                            goto REMOVE_READS_NP;
                                        }
                                        
                                        $remove_reads = "yes";
                                        print {$filehandle{$seed_id2}} $add_no_match_reads." REMOVE_BY_NO_MATCH_READS0\n";
                                        print {$filehandle{$seed_id2}} $nuc_rej_tmp." REMOVE_BY_NO_MATCH_READS1\n";
                                        print {$filehandle{$seed_id2}} $total_nuc_count_rej." TOTAL_NUC_REJ\n";
                                        $add_no_match_reads = "";
                                        
                                        foreach my $rank_tmp8 (sort {$a <=> $b} keys %split_patterns_final)
                                        {
                                            if ($split_patterns_final{$rank_tmp8} eq $nuc_rej_tmp)
                                            {
                                                $reads_to_remove{$rank_tmp8}{$rank_tmp8} = undef;
                                                my $id_tmp0 = $rank_to_id{$rank_tmp8};
                                            }
                                        }
                                        foreach my $id_tmp5 (keys %extensions_nomatch2b)
                                        {
                                            delete $extensions2_tmp{$id_tmp5};
                                        }
                                        delete $scores2{'0'};
                                    }
                                }
                            }
                        }
#Add rejected reads to mafft-------------------------------------------------------------------------------------------------------------------------------                                        

                        $lowest_longest_match = "";
						my $lowest_longest_match_nuc = "";
						foreach my $nuc_tmp14 (keys %longest_match)
						{
							my $longest_match_tmp = $longest_match{$nuc_tmp14};
							if ($lowest_longest_match eq "" || $longest_match_tmp < $lowest_longest_match)
							{
								$lowest_longest_match = $longest_match_tmp;
								$lowest_longest_match_nuc = $nuc_tmp14;
							}
							if ($longest_longest_match eq "" || $longest_match_tmp > $longest_longest_match)
							{
								$longest_longest_match = $longest_match_tmp;
							}
						}
						
						
						if ($add_rejected_reads eq "" && $add_no_match_reads eq "" && $count_matches_with_high_scores > $sequencing_depth_NP*1.4 && $remove_reads eq "" &&
                            ($post_pattern_match eq "yes2" || $post_pattern_match eq "yes3" || $SNP_check eq "yes2") && $post_pattern_match_average eq "yes" && $post_pattern_match_extra eq "yes" && $find_haps_in_seed eq "")
                        {        
NEW_LONGEST_MATCH_NP:							
							print {$filehandle{$seed_id2}} $lowest_longest_match." ADD_REJECTED_READS\n";   
							foreach my $nuc_tmp14 (keys %longest_match2)
                            {
                                if ($nuc_tmp14 eq $lowest_longest_match_nuc)
								{
									my $check_tmp = "";
									my $next_lowest_longest_match = "";
									foreach my $length_tmp14 (sort {$b <=> $a} keys %{$longest_match2{$nuc_tmp14}})
									{
										if ($length_tmp14 eq $lowest_longest_match)
										{
											$check_tmp = "yes";
										}
										elsif ($check_tmp eq "yes")
										{
											$next_lowest_longest_match = $length_tmp14;
											$check_tmp = "";
										}
									}
									foreach my $length_tmp14 (sort {$b <=> $a} keys %{$longest_match2{$nuc_tmp14}})
									{
										if ($length_tmp14 eq $lowest_longest_match)
										{
											foreach my $id_tmp14 (keys %{$longest_match2{$nuc_tmp14}{$length_tmp14}})
											{
												if (exists($multi_match{$id_tmp14}))
												{
													foreach my $mm_pos_tmp14 (sort {$b <=> $a} keys %{$longest_match2{$id_tmp14}})
													{
														if ($mm_pos_tmp14 > $next_lowest_longest_match-5000)
														{
															$lowest_longest_match = $next_lowest_longest_match;
															goto NEW_LONGEST_MATCH_NP;
														}
													}
												}
											}
										}
									}
								}
                            }

                            my $count_tmp = '0';
                            my $count_tmp3 = '0';
                            my %read_start_pos_rej_tmp;
                            undef %read_start_pos_rej_tmp;
                            
                            foreach my $id_tmp5 (keys %read_start_pos_rej)
                            {                                   
								if ($read_start_pos_rej{$id_tmp5} > $lowest_longest_match && $read_start_pos_rej{$id_tmp5} < $longest_longest_match)
                                {
                                    $count_tmp++;
									print {$filehandle{$seed_id2}} $read_start_pos_rej{$id_tmp5}." REJ\n";  
                                    if (exists($read_start_pos_rej_tmp{$read_start_pos_rej{$id_tmp5}}))
                                    {
                                        $read_start_pos_rej_tmp{$read_start_pos_rej{$id_tmp5}}{$id_tmp5} = undef;
                                    }
                                    else
                                    {
                                        my $check_tmp = "";
                                        foreach my $pos_tmp1 (keys %read_start_pos_rej_tmp)
                                        {
                                            if ($pos_tmp1 > $read_start_pos_rej{$id_tmp5}-350 && $pos_tmp1 < $read_start_pos_rej{$id_tmp5}+350)
                                            {
                                                $read_start_pos_rej_tmp{$pos_tmp1}{$id_tmp5} = undef;
                                                $check_tmp = "yes";
                                                last;
                                            }
                                        }
                                        if ($check_tmp eq "")
                                        {
                                            $read_start_pos_rej_tmp{$read_start_pos_rej{$id_tmp5}}{$id_tmp5} = undef;
                                        }
                                    }
                                }
                            }
                            foreach my $pos_tmp5 (keys %read_start_pos_rej_tmp)
                            {
                                if (keys %{$read_start_pos_rej_tmp{$pos_tmp5}} > 1)
                                {
                                    $count_tmp3++;
                                }
                            }
							print {$filehandle{$seed_id2}} $count_tmp." COUNT_TMP\n";  
                            print {$filehandle{$seed_id2}} $count_tmp3." COUNT_TMP3\n";   
                            if ((($count_tmp > 1 && $count_tmp3 > 0 && $count_tmp3 < 3) || ($count_tmp > '0' && $post_pattern_match_extra eq "yes" && $post_pattern_match eq "yes3" && $post_pattern_match_average eq "yes"))
								&& $count_tmp < 20)
                            {
                                my $ff = keys %read_start_pos_rej_tmp;
								print {$filehandle{$seed_id2}} $ff." KEYS_REJ\n";  
								undef %id_matches;
                                foreach my $pos_tmp5 (keys %read_start_pos_rej_tmp)
                                {                                   
									my $rej_reads_count = keys %{$read_start_pos_rej_tmp{$pos_tmp5}};
									if ($count_tmp3 eq '0' || $rej_reads_count > 1)
									{
										foreach my $id_tmp5 (sort {$a <=> $b} keys %{$read_start_pos_rej_tmp{$pos_tmp5}})
										{ 
											print {$filehandle{$seed_id2}} $id_tmp5." ID ".$pos_tmp5." POS_REJ\n";
											
											if (exists($read_start_pos_rej_saved{$id_tmp5}))
											{
												$id_matches{$id_tmp5} = undef;
												print {$filehandle{$seed_id2}} $id_tmp5." ID ".$pos_tmp5." POS_REJ000\n";
												delete $save_alignment_data_NP{$seed_id}{$id_tmp5};
												delete $read_start_pos_rej{$id_tmp5};
											}
											else
											{
												my $long_read_tmp = "";
												if (exists($reverse_list{$id_tmp5}))
												{
													$long_read_tmp = reverse($hash_NP_reads_tmp{$id_tmp5});
													$long_read_tmp =~ tr/ACTG/TGAC/;
												}
												else
												{
													$long_read_tmp = $hash_NP_reads_tmp{$id_tmp5};
												}
												
												my $long_read_end_pos_tmp = $long_read_end_pos_save{$id_tmp5};
												my $ext = substr $long_read_tmp, $long_read_end_pos_tmp-90, $length_extension;
												if (exists($scores2{'0'}))
												{
													$scores2{'0'} .= ",$id_tmp5"; 
												}
												else
												{
													$scores2{'0'} = $id_tmp5;
												}
												$extensions2_tmp{$id_tmp5} = $ext;
											}
										}
									}
                                }
								foreach my $pos_tmp14 (keys %{$split_positions_DUP{$id}})
								{
									if ($pos_tmp14 > $position)
									{
										 delete $split_positions_DUP{$id}{$pos_tmp14};
									}
								}
								
								my $overlap_tmp = $overlap;
                                my $d = length($best_extension) - $overlap;
                                if ($d < 0)
                                {
                                   $overlap_tmp -= $d;
                                   $d = '0';                         
                                }
                                my $read_end_tmpi = substr $best_extension, $d, $overlap_tmp;
                                my $count_del = $read_end_tmpi =~ tr/-/-/;
                                my $h = '0';
                                while ($count_del > $h && $d > 0)
                                {
                                    $h = $count_del;
                                    $read_end_tmpi = substr $best_extension, $d-$count_del, $overlap_tmp+$count_del;
                                    $count_del = $read_end_tmpi =~ tr/-/-/;
                                }
                                $read_end_tmpi =~ tr/-//d;
                                $read_end_tmpi =~ tr/actgn/ACTGN/;
                                my $pos_tmp = $position+length($best_extension);
                                
                                if ($selected_nuc eq "")
                                {
                                    my $f = '0';
                                    foreach my $nuc_split_tmp (keys %nucs_for_split_extra)
                                    {
                                        if ($nucs_for_split_extra{$nuc_split_tmp} > $f)
                                        {
                                            $f = $nucs_for_split_extra{$nuc_split_tmp};
                                            $selected_nuc = $nuc_split_tmp;
                                        }
                                    }
                                }
								
								if ($SNP_check eq "yes2" && $trace_back_check eq "yes")
                                {
                                    my $pos_tmp2 = $track_length_ext_total{'1'};
                                    $split_positions_DUP_tmp{$pos_tmp2} = $read_end_tmpi.",".$selected_nuc;
									print {$filehandle{$seed_id2}} $seed_id."\t".$pos_tmp."\t".$read_end_tmpi.",".$selected_nuc." DUP_POS\n";
                                }
								
                                $add_rejected_reads = $total_nuc_count_original;
                                print {$filehandle{$seed_id2}} $best_extension."\nHALLE1\n";
                                $mismatch_retry++;
                                $best_extension = "";
								$best_extension_part = "";
                                undef %quality_scores_tmp;
                                undef %span_complex_region;
                                $remove_reads_check = "";
                                
                                if (keys %id_matches > 0)
                                {
                                    goto ADD_REJ_POS_NP;
                                }
                                goto SELECT_LENGTH_NP2;   
                            }
                        }
                        elsif ($add_rejected_reads ne "" && $longest_longest_match_nuc ne "" && ($post_pattern_match eq "yes3" || ($SNP_check eq "yes2" && $post_pattern_match_extra eq "yes" && $post_pattern_match_average eq "yes")))
                        {          
                            print {$filehandle{$seed_id2}} $longest_longest_match_nuc." LONGEST_NUC\n";
							foreach my $nuc_rej_tmp (keys %nucs_rej)
                            {
                                print {$filehandle{$seed_id2}} $nuc_rej_tmp." NUC_REJ\n";
								if ($nuc_rej_tmp eq $longest_longest_match_nuc && $nucs_rej{$nuc_rej_tmp} > $total_nuc_count_rej*0.04)
                                {
                                    goto REMOVE_READS_NP;
                                }
                            }
                            foreach my $nuc_rej_tmp (keys %nucs_rej)
                            {
                                my $count_tmp = keys %nucs_rej;
								my $count_tmp2 = "0";
								foreach my $rank_tmp8 (sort {$a <=> $b} keys %split_patterns_final)
								{
									if ($split_patterns_final{$rank_tmp8} eq $nuc_rej_tmp)
									{
										$count_tmp2++;
									}
								}
								if ($nuc_rej_tmp ne $longest_longest_match_nuc && ($nucs_rej{$nuc_rej_tmp} > 2 || ($nucs_rej{$nuc_rej_tmp} > 1 && $count_tmp eq "1")
									|| ($nucs_rej{$nuc_rej_tmp} > 0 && $count_tmp eq "1" && $post_pattern_match eq "yes3")) && $count_tmp2 > 0.15*$total_nuc_count)
                                {
                                    $remove_reads = "yes";
                                    print {$filehandle{$seed_id2}} $add_rejected_reads." REMOVE_BY_REJECTED_READS0\n";
                                    print {$filehandle{$seed_id2}} $nuc_rej_tmp." REMOVE_BY_REJECTED_READS1\n";
                                    print {$filehandle{$seed_id2}} $total_nuc_count_rej." TOTAL_NUC_REJ\n";
                                    $add_rejected_reads = "";
                                    
                                    foreach my $rank_tmp8 (sort {$a <=> $b} keys %split_patterns_final)
                                    {
                                        if ($split_patterns_final{$rank_tmp8} eq $nuc_rej_tmp)
                                        {
                                            $reads_to_remove{$rank_tmp8}{$rank_tmp8} = undef;
                                            my $id_tmp0 = $rank_to_id{$rank_tmp8};
                                            delete $read_start_pos_rej{$id_tmp0};
                                        }
                                    }
                                    foreach my $id_tmp5 (keys %read_start_pos_rej)
                                    {
                                        delete $extensions2_tmp{$id_tmp5};
                                    }
                                    delete $scores2{'0'};
									
									my $last_10 = substr $best_extension, -10, 10;
									#$track_split_NP{$id}{$position+length($best_extension)} = $last_10."+".$longest_longest_match_nuc;
									foreach my $posie_tmp (keys %SNP_patterns_prev_match)
									{
										my $last_10_tmp = substr $best_extension, $posie_tmp-$position-10, 10;
										#$track_split_NP{$id}{$posie_tmp} = $last_10_tmp."+".$SNP_patterns_prev_match{$posie_tmp};
									}
                                }
                            }
                        }					
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------                                          

REMOVE_READS_NP:                            
                        if (keys %reads_to_remove > 0)
						{}
						elsif ($remove_reads ne "")
						{
							$remove_reads = "";
							print {$filehandle{$seed_id2}} "REMOVE_ERROR\n";
							goto BASECALL2_NP;
						}
						if ($remove_reads eq "yes")
                        {
                            my $removed_reads_tmp = '0';
                            my $highest_rank = "";
                            my $add_rejected_reads_new = '0';
							
							my $reads_to_remove = '0';
							foreach my $rank_tmp0 (keys %reads_to_remove)
                            {
								foreach my $rank_tmp2 (keys %{$reads_to_remove{$rank_tmp0}})
								{
									$reads_to_remove++;
								}
							}
							if ($total_nuc_count - $reads_to_remove < 4 && $full_reset_NP ne "" && $post_pattern_match ne "yes3")
							{
								$remove_reads = "";
								undef %reads_to_remove;
								goto BASECALL2_NP;
							}
                            
                            foreach my $rank_tmp (sort {$a <=> $b} keys %subject_list)
                            {
                                my $removed = "";
                                foreach my $rank_tmp0 (keys %reads_to_remove)
                                {
                                    if (exists($reads_to_remove{$rank_tmp0}{$rank_tmp}))
                                    {
                                        $removed_reads_tmp++;
                                        $removed = "yes";
                                        my $id_tmp0 = $rank_to_id{$rank_tmp};
                                        delete $extensions{$extensions2{$id_tmp0}};
                                        delete $extensions2{$id_tmp0};
                                        delete $extensions2b{$id_tmp0};
                                        delete $extensions2_tmp{$id_tmp0};
                                        delete $extensions_nomatch{$extensions2{$id_tmp0}};
                                        delete $extensions_nomatch2{$id_tmp0};
                                        delete $extensions_nomatch2b{$id_tmp0};
                                        delete $extensions_unknown{$extensions2{$id_tmp0}};
                                        delete $extensions_unknown2{$id_tmp0};
                                        delete $save_reads_for_next{$id_tmp0};
                                        delete $add_rej_reads_extra{$id_tmp0};
                                        delete $read_start_pos_rej{$id_tmp0};
                                        if ($rank_tmp0 <= $add_rejected_reads && $add_rejected_reads ne "")
                                        {
                                            $add_rejected_reads_new++;
                                        }
                                    }
                                }
                                if ($highest_rank eq "" && $removed eq "")
                                {
                                    $highest_rank = $rank_tmp;
                                }
                            }
                            if ($removed_reads_tmp > 0 && keys %extensions2_tmp > 1)
                            {
                                my $overlap_tmp = $overlap;
                                my $d = length($best_extension) - $overlap;
                                if ($d < 0)
                                {
                                   $overlap_tmp -= $d;
                                   $d = '0';                         
                                }
                                my $read_end_tmpi = substr $best_extension, $d, $overlap_tmp;
                                my $count_del = $read_end_tmpi =~ tr/-/-/;
                                my $h = '0';
                                while ($count_del > $h && $d > 0)
                                {
                                    $h = $count_del;
                                    $read_end_tmpi = substr $best_extension, $d-$count_del, $overlap_tmp+$count_del;
                                    $count_del = $read_end_tmpi =~ tr/-/-/;
                                }
                                $read_end_tmpi =~ tr/-//d;
                                $read_end_tmpi =~ tr/actgn/ACTGN/;
                                my $pos_tmp = $position+length($best_extension);
                                
                                if ($selected_nuc eq "")
                                {
                                    my $f = '0';
                                    foreach my $nuc_split_tmp (keys %nucs_for_split_extra)
                                    {
                                        if ($nucs_for_split_extra{$nuc_split_tmp} > $f)
                                        {
                                            $f = $nucs_for_split_extra{$nuc_split_tmp};
                                            $selected_nuc = $nuc_split_tmp;
                                        }
                                    }
                                }
                                if ($SNP_check eq "yes2" && ((($post_pattern_match_extra eq "yes" || $high_score_save > 15) && $highest_first_no_match > 9) || $trace_back_check eq "yes"))
                                {
                                    my $pos_tmp2 = $track_length_ext_total{'1'};
                                    $split_positions_DUP_tmp{$pos_tmp2} = $read_end_tmpi.",".$selected_nuc;
									print {$filehandle{$seed_id2}} $seed_id."\t".$pos_tmp."\t".$read_end_tmpi.",".$selected_nuc." DUP_POS\n";
                                }

                                foreach my $pos_tmp14 (keys %{$split_positions_DUP{$id}})
								{
									if ($pos_tmp14 > $position)
									{
										 #delete $split_positions_DUP{$id}{$pos_tmp14};
									}
								}
                                print {$filehandle{$seed_id2}} $best_extension."\nHALLE2\n";
													
                                $mismatch_retry++;
                                $best_extension = "";
								$best_extension_part = "";
                                undef %quality_scores_tmp;
                                undef %span_complex_region;
								$remove_reads_check = "yes";
                                
                                if ($total_nuc_count_original > $add_rejected_reads && $add_rejected_reads ne "")
                                {
                                    $add_rejected_reads -= $add_rejected_reads_new;
                                    goto SELECT_LENGTH_NP;
                                }
                                if ($mismatch_retry > 500000000 && $total_nuc_count-$removed_reads_tmp > 7 && $high_quality eq "")
                                {
                                    $high_quality = "yes";
                                    goto SELECT_LENGTH_NP;
                                }
           
                                goto MISMATCH_RETRY_NP;
                            }
                        }
#Select reads based on length--------------------------------------------------------------------------------------  
                        elsif ($post_pattern_match eq "yes3" && $find_haps_in_seed eq "")
                        {
                            my $nuc_pattern = "";
                            my $extra_length = "";
                            my $extra_length2 = "";
                            my $count_tmp = '0';
                            foreach my $rank_tmp (sort {$a <=> $b} keys %span_complex_region)
                            {
                                print {$filehandle{$seed_id2}} $rank_tmp." SPAN_COMPLEX_REGION\n";
                                if (exists($split_patterns_final{$rank_tmp}))
                                {
                                    if ($nuc_pattern eq "")
                                    {
                                        $nuc_pattern = $split_patterns_final{$rank_tmp};
                                        $extra_length = $span_complex_region{$rank_tmp}-($position-$last_non_complex_region{$seed_id});
                                        $extra_length2 = $span_complex_region{$rank_tmp}-$repetitive_detect2;
                                        $count_tmp++;
                                    }
                                    elsif ($nuc_pattern ne $split_patterns_final{$rank_tmp})
                                    {
                                        goto ALL_MISMATCHES_NP;
                                    }
                                    else
                                    {
                                        $count_tmp++;
                                    }
                                }
                            }
                            if ($count_tmp > 1 && $extra_length > 5000 && $extra_length2 > 5000)
                            {
                                $remove_reads = "yes";
                                print {$filehandle{$seed_id2}} $count_tmp." REMOVE_BY_LENGTH\n";
                                
                                foreach my $rank_tmp (sort {$a <=> $b} keys %split_patterns_final)
                                {
                                    if ($split_patterns_final{$rank_tmp} ne $nuc_pattern)
                                    {
                                        $reads_to_remove{$rank_tmp}{$rank_tmp} = undef;
                                    }
                                }
                                goto REMOVE_READS_NP;
                            }
						}
ALL_MISMATCHES_NP:						
						if ($post_pattern_match_extra eq "yes" && (($post_pattern_match eq "yes2" && $count_matches_with_high_scores > $sequencing_depth_NP) || $post_pattern_match eq "yes3") && $post_pattern_match_average eq "yes" && $find_haps_in_seed eq "")
						{
#Select reads based on all mismatches--------------------------------------------------------------------------------------  
							my %mismatch_pattern_all;
							undef %mismatch_pattern_all;
							my %mismatch_pattern;
							undef %mismatch_pattern;
							my %mismatch_pattern_N;
							undef %mismatch_pattern_N;
							my $count_groups_tmp = '0';
							foreach my $nuc_tmp13 (sort {$a <=> $b} keys %split_patterns_final_score)
							{
								my $total_count_tmp = '0';
								my $total_score_count_tmp = '0';
								my $count_ranks_tmp = keys %{$split_patterns_final_score{$nuc_tmp13}};
								if ($count_ranks_tmp > 3)
								{
									$count_groups_tmp++;
									foreach my $rank_tmp13 (keys %{$split_patterns_final_score{$nuc_tmp13}})
									{
										if (exists($rank_to_id{$rank_tmp13}))
										{
											my $id_tmpi = $rank_to_id{$rank_tmp13};
											foreach my $pos_snp_tmp (sort {$a <=> $b} keys %{$store_mismatches_NP{$id_tmpi}})
											{
												$mismatch_pattern{$pos_snp_tmp}{$nuc_tmp13}{$id_tmpi} = undef;
											}
											foreach my $pos_snp_tmp (sort {$a <=> $b} keys %{$store_mismatches_all_NP{$id_tmpi}})
											{
												my @store_mismatches_all_NP = split /,/, $store_mismatches_all_NP{$id_tmpi}{$pos_snp_tmp};
												$mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp13}{$id_tmpi} = $store_mismatches_all_NP[0]."+".$store_mismatches_all_NP[1];
											}
											foreach my $pos_snp_tmp (sort {$a <=> $b} keys %{$store_mismatches_N_NP{$id_tmpi}})
											{
												my @store_mismatches_N_NP = split /,/, $store_mismatches_N_NP{$id_tmpi}{$pos_snp_tmp};
												$mismatch_pattern_N{$pos_snp_tmp}{$nuc_tmp13}{$store_mismatches_N_NP[1]}{$id_tmpi} = undef;
											}
										}
									}
								}	
							}
							if ($count_groups_tmp > 1)
							{
								my $count_tmp = keys %mismatch_pattern_all;
								print {$filehandle{$seed_id2}} $count_tmp." CHECK_ALL_MISMATCHES\n";
								my %pos_no_matches;
								undef %pos_no_matches;
								my %pos_no_matches1;
								undef %pos_no_matches1;
								my %rank_pattern;
								undef %rank_pattern;

								foreach my $pos_snp_tmp (sort {$b <=> $a} keys %mismatch_pattern_all)
								{
									my $one_nuc_match = "";
									my $one_nuc_no_match = "";
									
									if (exists($mismatch_pattern{$pos_snp_tmp}))
									{
										my $total_rank_count = '0';
										foreach my $nuc_tmp14 (keys %{$mismatch_pattern{$pos_snp_tmp}})
										{
											$total_rank_count += keys %{$mismatch_pattern{$pos_snp_tmp}{$nuc_tmp14}}
										}
										foreach my $nuc_tmp14 (keys %{$mismatch_pattern{$pos_snp_tmp}})
										{
											my $this_rank_count_tmp = keys %{$mismatch_pattern{$pos_snp_tmp}{$nuc_tmp14}};
											if ($this_rank_count_tmp < 2 && $this_rank_count_tmp < 0.15*$total_rank_count)
											{
												foreach my $id_tmpi14 (sort {$a <=> $b} keys %{$mismatch_pattern{$pos_snp_tmp}{$nuc_tmp14}})
												{
													delete $mismatch_pattern{$pos_snp_tmp}{$nuc_tmp14}{$id_tmpi14};
												}
												delete $mismatch_pattern{$pos_snp_tmp}{$nuc_tmp14};
											}
										}
										if (keys %{$mismatch_pattern{$pos_snp_tmp}} eq '1')
										{
											foreach my $nuc_tmp14 (keys %{$mismatch_pattern{$pos_snp_tmp}})
											{
												if (keys %{$mismatch_pattern{$pos_snp_tmp}{$nuc_tmp14}} > 3)
												{
													my $rank_pattern_tmp = "";
													foreach my $id_tmpi14 (sort {$a <=> $b} keys %{$mismatch_pattern{$pos_snp_tmp}{$nuc_tmp14}})
													{
														if ($rank_pattern_tmp eq "")
														{
															$rank_pattern_tmp = $id_tmpi14;
														}
														else
														{
															$rank_pattern_tmp .= "_".$id_tmpi14
														}
													}
													$rank_pattern{$rank_pattern_tmp}{$nuc_tmp14}{$pos_snp_tmp} = 3;
												}
											}
										}

										foreach my $nuc_tmp13 (sort {$a <=> $b} keys %split_patterns_final_score)
										{
											my $nuc_mismatches_count = '0';
											my $too_short_count = '0';
											my $match_count = '0';
											my $total_count_tmp = '0';
											
											foreach my $rank_tmp13 (keys %{$split_patterns_final_score{$nuc_tmp13}})
											{
												if (exists($rank_to_id{$rank_tmp13}))
												{
													my $id_tmpi = $rank_to_id{$rank_tmp13};
													if (exists($mismatch_pattern{$pos_snp_tmp}{$nuc_tmp13}{$id_tmpi}))
													{
														$nuc_mismatches_count++;
													}
													elsif (exists($alignment_length_save{$id_tmpi}))
													{
														my $overlap_tmp = $alignment_length_save{$id_tmpi};
														if ($position - $pos_snp_tmp > $overlap_tmp)
														{
															$too_short_count++;
														}
														else
														{
															$match_count++;
														}
													}
												}
												$total_count_tmp++;
											}
											if (($match_count < 2 && $nuc_mismatches_count > 3 && $nuc_mismatches_count > $total_count_tmp*0.3 && $nuc_mismatches_count > $match_count*4
												&& $match_count < $total_count_tmp*0.1) || ($match_count eq '0' && $nuc_mismatches_count > 2))
											{
												$one_nuc_no_match = $nuc_tmp13;
											}
											elsif (($nuc_mismatches_count < $total_count_tmp*0.1 && $match_count > 2 && $match_count > 3*$nuc_mismatches_count
												  ) || ($nuc_mismatches_count eq '0' && $match_count > 2))
											{
												$one_nuc_match = $nuc_tmp13;
											}
										}
									}
									
									if ($one_nuc_match ne "" && $one_nuc_no_match ne "")
									{
										$pos_no_matches{$one_nuc_no_match} += 2;
										$pos_no_matches1{$one_nuc_no_match} += 2;
										
										print {$filehandle{$seed_id2}} $pos_snp_tmp." FINAL_PATTERN_MATCH1\n";
										print {$filehandle{$seed_id2}} $one_nuc_no_match." FINAL_PATTERN_MATCH_NUC1\n";
									}
									else
									{
										$one_nuc_match = "";
										$one_nuc_no_match = "";
									
										my $total_rank_count = '0';
										foreach my $nuc_tmp14 (keys %{$mismatch_pattern_all{$pos_snp_tmp}})
										{
											$total_rank_count += keys %{$mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp14}}
										}
										foreach my $nuc_tmp14 (keys %{$mismatch_pattern_all{$pos_snp_tmp}})
										{
											my $this_rank_count_tmp = keys %{$mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp14}};
											if ($this_rank_count_tmp < 2 && $this_rank_count_tmp < 0.15*$total_rank_count)
											{
												foreach my $id_tmpi14 (sort {$a <=> $b} keys %{$mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp14}})
												{
													delete $mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp14}{$id_tmpi14};
												}
												delete $mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp14};
											}
										}
										
										if (keys %{$mismatch_pattern_all{$pos_snp_tmp}} eq '1')
										{
											foreach my $nuc_tmp14 (keys %{$mismatch_pattern_all{$pos_snp_tmp}})
											{
												if (keys %{$mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp14}} > 3)
												{
													my $rank_pattern_tmp = "";
													foreach my $id_tmpi14 (sort {$a <=> $b} keys %{$mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp14}})
													{
														if ($rank_pattern_tmp eq "")
														{
															$rank_pattern_tmp = $id_tmpi14;
														}
														else
														{
															$rank_pattern_tmp .= "_".$id_tmpi14
														}
													}
													$rank_pattern{$rank_pattern_tmp}{$nuc_tmp14}{$pos_snp_tmp} = 1;
												}
											}
										}
	
										foreach my $nuc_tmp13 (sort {$a <=> $b} keys %split_patterns_final_score)
										{
											my $nuc_mismatches_count = '0';
											my $too_short_count = '0';
											my $match_count = '0';
											my $total_count_tmp = '0';
											
											foreach my $rank_tmp13 (keys %{$split_patterns_final_score{$nuc_tmp13}})
											{
												if (exists($rank_to_id{$rank_tmp13}))
												{
													my $id_tmpi = $rank_to_id{$rank_tmp13};
													if (exists($mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp13}{$id_tmpi}))
													{
														$nuc_mismatches_count++;
													}
													elsif (exists($alignment_length_save{$id_tmpi}))
													{
														my $overlap_tmp = $alignment_length_save{$id_tmpi};
														if ($position - $pos_snp_tmp > $overlap_tmp)
														{
															$too_short_count++;
														}
														else
														{
															$match_count++;
														}
													}
												}
												$total_count_tmp++;
											}
											if (($match_count < 2 && $nuc_mismatches_count > 3 && $nuc_mismatches_count > $total_count_tmp*0.3 && $nuc_mismatches_count > $match_count*4
												&& $match_count < $total_count_tmp*0.1) || ($match_count eq '0' && $nuc_mismatches_count > 2))
											{
												$one_nuc_no_match = $nuc_tmp13;
											}
											elsif (($nuc_mismatches_count < $total_count_tmp*0.1 && $match_count > 2 && $match_count > 3*$nuc_mismatches_count
												  ) || ($nuc_mismatches_count eq '0' && $match_count > 2))
											{
												$one_nuc_match = $nuc_tmp13;
											}
										}
										if ($one_nuc_match ne "" && $one_nuc_no_match ne "")
										{
											$pos_no_matches{$one_nuc_no_match} += 1;
											
											print {$filehandle{$seed_id2}} $pos_snp_tmp." FINAL_PATTERN_MATCH2\n";
											print {$filehandle{$seed_id2}} $one_nuc_no_match." FINAL_PATTERN_MATCH_NUC2\n";
										}
									}
								}
	
								my $nuc_to_delete = '0';
								my $nuc_to_delete2 = '0';
								my $count_below_5 = '0';
								my $high_pos_count = "";
								my $high_pos_nuc = "";
								my $rank1_check = "";
								foreach my $nuc (keys %pos_no_matches)
								{
									my $pos_count = $pos_no_matches{$nuc};
									$nuc_to_delete++;
									print {$filehandle{$seed_id2}} $nuc." NUC\n";
									print {$filehandle{$seed_id2}} $pos_count." POS_COUNT\n";
									if ($pos_count > 5 || (keys %pos_no_matches < 2 && $pos_count > 4 && $post_pattern_match eq "yes3"))
									{
										$nuc_to_delete2++;
										foreach my $rank_tmp14 (keys %{$split_patterns_final_score{$nuc}})
										{
											$reads_to_remove{$rank_tmp14}{$rank_tmp14} = undef;
											if ($rank_tmp14 < 5)
											{
												$count_below_5++;
											}
											if ($rank_tmp14 eq "1")
											{
												$rank1_check = "yes";
												foreach my $length_tmp (keys %longest_match)
												{
													if ($length_tmp ne $longest_longest_match && $length_tmp > 0.85*$longest_longest_match)
													{
														$rank1_check = "";
													}
												}
											}
										}
										$high_pos_count = $pos_count;
										$high_pos_nuc = $nuc;
									}
									print {$filehandle{$seed_id2}} $count_below_5." COUNT_BELOW_5\n";
								}
								my $multi_match_check = "";
								foreach my $nuc (keys %pos_no_matches)
								{
									if ($high_pos_nuc ne $nuc)
									{
										my $pos_count = $pos_no_matches{$nuc};
										if ($pos_count > $high_pos_count*0.2)
										{
											$multi_match_check = "no";
										}
									}
								}
								
								my $pos_no_matches1_check = "";
								if (keys %pos_no_matches1 eq '1')
								{
									foreach my $nuci (keys %pos_no_matches1)
									{
										if ($pos_no_matches1{$nuci} > 4)
										{
											$pos_no_matches1_check = $nuci;
											$high_pos_nuc = $nuci;
										}
									}
								}
								if ($pos_no_matches1_check ne "" && $nuc_to_delete2 > 1)
								{
									undef %reads_to_remove;
									foreach my $rank_tmp14 (keys %{$split_patterns_final_score{$pos_no_matches1_check}})
									{
										$reads_to_remove{$rank_tmp14}{$rank_tmp14} = undef;
									}
								}
								if (((($nuc_to_delete eq '1' || $multi_match_check eq "") && $nuc_to_delete2 eq '1' && ($count_below_5 < 3 || ($high_pos_count > 20 && $count_below_5 < 4 && $total_nuc_count > 10)))
									 || ($pos_no_matches1_check ne "")) && $rank1_check eq "")
								{
									$remove_reads = "yes";
									print {$filehandle{$seed_id2}} $high_pos_nuc." REMOVE_BY_all_MISMATCH\n";
									
									my $one_nuc_match_tmp = "";
									foreach my $nuc_tmp13 (sort {$a <=> $b} keys %split_patterns_final_score)
									{
										my $count_ranks_tmp = keys %{$split_patterns_final_score{$nuc_tmp13}};
										if ($nuc_tmp13 ne $high_pos_nuc && $count_ranks_tmp > $total_nuc_count*0.25)
										{
											if ($one_nuc_match_tmp eq "")
											{
												$one_nuc_match_tmp = $nuc_tmp13
											}
											else
											{
												$one_nuc_match_tmp = "";
												last;
											}
										}
									}
									if ($one_nuc_match_tmp ne "")
									{
										my $last_10 = substr $best_extension, -10, 10;
										#$track_split_NP{$id}{$position+length($best_extension)} = $last_10."+".$one_nuc_match_tmp;
										foreach my $posie_tmp (keys %SNP_patterns_prev_match)
										{
											my $last_10_tmp = substr $best_extension, $posie_tmp-$position-10, 10;
											#$track_split_NP{$id}{$posie_tmp} = $last_10_tmp."+".$SNP_patterns_prev_match{$posie_tmp};
										}
									}
									
									goto REMOVE_READS_NP;
								}
#Go back to resolve previous SNP position-----------------------------------------------------------------------------------------------------------------
								elsif ($post_pattern_match eq "yes3")
								{
									my $found_pos = "";
									foreach my $pos_snp_tmp (sort {$a <=> $b} keys %mismatch_pattern_N)
									{
										my %N_matches;
										undef %N_matches;
											
										foreach my $nuc_split14 (keys %{$mismatch_pattern_N{$pos_snp_tmp}})
										{
											my $total_pattern_count_tmp = '0';
											my $split_patterns_count_tmp = keys %{$split_patterns_final_score{$nuc_split14}};
											foreach my $nuc_tmp14 (keys %{$mismatch_pattern_N{$pos_snp_tmp}{$nuc_split14}})
											{
												$total_pattern_count_tmp += keys %{$mismatch_pattern_N{$pos_snp_tmp}{$nuc_split14}{$nuc_tmp14}};
											}
											
											foreach my $nuc_tmp14 (keys %{$mismatch_pattern_N{$pos_snp_tmp}{$nuc_split14}})
											{
												my $nuc_count_tmp4 = keys %{$mismatch_pattern_N{$pos_snp_tmp}{$nuc_split14}{$nuc_tmp14}};
												if ($nuc_count_tmp4 > 0.8*$total_pattern_count_tmp && $nuc_count_tmp4 > 3 && ($nuc_count_tmp4 > $split_patterns_count_tmp*0.45 || $nuc_count_tmp4 > 5))
												{
													$N_matches{$nuc_tmp14} = $nuc_count_tmp4;
												}
											}
										}
										if (keys %N_matches > 1)
										{
											print {$filehandle{$seed_id2}} $pos_snp_tmp." FINAL_N_PATTERN_MATCH\n";
											my $indel_check_tmp = "";
											foreach my $nuc_tmp15 (keys %N_matches)
											{
												if ($nuc_tmp15 eq "-")
												{
													$indel_check_tmp = "yes";
												}
											}
N_MATCHES_NP:								foreach my $nuc_tmp15 (keys %N_matches)
											{
												print {$filehandle{$seed_id2}} $nuc_tmp15." FINAL_N_PATTERN_MATCH_NUC\n";
												print {$filehandle{$seed_id2}} $N_matches{$nuc_tmp15}." FINAL_N_PATTERN_MATCH_COUNT\n";
												
												foreach my $pos_tb_tmp (keys %{$trace_back_split_NP{$id}})
												{
													if ($pos_snp_tmp > $pos_tb_tmp-30 && $pos_snp_tmp < $pos_tb_tmp+30)
													{
														next N_MATCHES_NP;
													}
												}
												if (exists($trace_back_split_NP{$id}{$pos_snp_tmp}))
												{}
												elsif ($indel_check_tmp eq "")
												{
													my $last_10 = substr $read, $pos_snp_tmp-11, 10;
													$trace_back_split_NP{$id}{$pos_snp_tmp} = $last_10;
													print {$filehandle{$seed_id2}} $last_10." LAST10\n";

													if ($found_pos eq "")
													{
														$found_pos = $pos_snp_tmp;
													}
												}
											}		
										}
									}
									if ($found_pos ne "")
									{
										$read = substr $read, 0, $found_pos-800;
										$best_extension = "";
										if ($hap_tag eq "HAP1")
										{
											undef %exclude_reads_hap1_NP;
										}
										elsif ($hap_tag eq "HAP2")
										{
											undef %exclude_reads_hap2_NP;
										}
										undef %save_alignment_data_NP;
										$position = length($read);
										$position{$id} = $position;
										$seed{$id} = $read;
									}
									
									my $found_pos2 = "";
									my $yy = keys %mismatch_pattern_all;
									print {$filehandle{$seed_id2}} $yy." MISMATCH_PATTER_ALL\n";
MISMATCH_PATTER_ALL:				foreach my $pos_snp_tmp (sort {$a <=> $b} keys %mismatch_pattern_all)
									{
										my %all_matches;
										undef %all_matches;
										my %nuc_matches;
										undef %nuc_matches;
										my %nuc_matches2;
										undef %nuc_matches2;	
										my $total_pattern_count_tmp = '0';
										my $nucie_tmp = '0';
										my $nucie_tmp2 = '0';
					my $within_range2 = "";					
										foreach my $nuc_tmp14 (keys %{$mismatch_pattern_all{$pos_snp_tmp}})
										{
											my $tmp_count = keys %{$mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp14}};
											$total_pattern_count_tmp += $tmp_count;
											foreach my $id_tmp14 (keys %{$mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp14}})
											{
												my $nuc_tmp2 = $mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp14}{$id_tmp14};
												my @nuc_tmp2 = split /\+/, $nuc_tmp2;
												$nuc_matches{$nuc_tmp2[0]} += 1;
												$nuc_matches2{$nuc_tmp2[1]} += 1;
											}
											
											foreach my $nuc_tmp15 (keys %nuc_matches)
											{
												if ($nuc_matches{$nuc_tmp15} > 0.88*$tmp_count && $tmp_count > 3)
												{
													if ($nuc_matches{$nuc_tmp15} ne "-")
													{
														$nucie_tmp += 1;
													}
												}
											}
											foreach my $nuc_tmp15 (keys %nuc_matches2)
											{
												if ($nuc_matches2{$nuc_tmp15} > 0.88*$tmp_count && $tmp_count > 3)
												{
													if ($nuc_matches2{$nuc_tmp15} ne "-")
													{
														$nucie_tmp2 += 1;
													}
												}
											}
										}
										
										
										if ($nucie_tmp < 1 || $nucie_tmp2 < 1)
										{
											next MISMATCH_PATTER_ALL;
										}
										
										foreach my $nuc_tmp14 (keys %{$mismatch_pattern_all{$pos_snp_tmp}})
										{
											my $within_range = '0';
											my $total_count_tmp = '0';
										
											foreach my $rank_tmp13 (keys %{$split_patterns_final_score{$nuc_tmp14}})
											{
												if (exists($rank_to_id{$rank_tmp13}))
												{
													my $id_tmpi = $rank_to_id{$rank_tmp13};
													
													if (exists($alignment_length_save{$id_tmpi}))
													{
														my $overlap_tmp = $alignment_length_save{$id_tmpi};
														if ($position - $pos_snp_tmp > $overlap_tmp)
														{
														}
														else
														{
															$within_range++;
														}
													}
												}
												$total_count_tmp++;
											}
								if ($within_range > '0')
											{
												$within_range2 = "yes";

											}
											print {$filehandle{$seed_id2}} $nuc_tmp14." ".$within_range." NUCIE_TMPP\n";
											my $split_patterns_count_tmp = keys %{$split_patterns_final_score{$nuc_tmp14}};
											my $nuc_count_tmp4 = keys %{$mismatch_pattern_all{$pos_snp_tmp}{$nuc_tmp14}};
											print {$filehandle{$seed_id2}} $total_pattern_count_tmp." ".$split_patterns_count_tmp." ".$nuc_count_tmp4." NUCIE_TMPP2\n";
											if ($nuc_count_tmp4 > 0.89*$total_pattern_count_tmp && $nuc_count_tmp4 > 3 && $nuc_count_tmp4 > 0.89*$within_range && ($nuc_count_tmp4 > $split_patterns_count_tmp*0.49 || $nuc_count_tmp4 > 5))
											{
												$all_matches{$nuc_tmp14} = $nuc_count_tmp4;
												print {$filehandle{$seed_id2}} $nuc_tmp14." FINAL_ALL_PATTERN_MATCH_TEST\n";
											}
										}
										if ($within_range2 eq "")
											{
												print {$filehandle{$seed_id2}} $pos_snp_tmp." ".$total_pattern_count_tmp." ".$within_range2." FINAL_ALL_PATTERN_MATCH_TEST\n";
											}
										if (keys %all_matches eq 1)
										{
											print {$filehandle{$seed_id2}} $pos_snp_tmp." FINAL_ALL_PATTERN_MATCH\n";

ALL_MATCHES_NP:  							foreach my $nuc_tmp15 (keys %all_matches)
											{
												print {$filehandle{$seed_id2}} $nuc_tmp15." FINAL_ALL_PATTERN_MATCH_NUC\n";
												print {$filehandle{$seed_id2}} $all_matches{$nuc_tmp15}." FINAL_ALL_PATTERN_MATCH_COUNT\n";
												
												foreach my $pos_tb_tmp (keys %{$trace_back_split_NP{$id}})
												{
													if ($pos_snp_tmp > $pos_tb_tmp-50 && $pos_snp_tmp < $pos_tb_tmp+50)
													{
														next ALL_MATCHES_NP;
													}
												}
												if (exists($trace_back_split_NP{$id}{$pos_snp_tmp}))
												{}
												else
												{
													my $last_10 = substr $read, $pos_snp_tmp-11, 10;
													$trace_back_split_NP{$id}{$pos_snp_tmp} = $last_10;
													print {$filehandle{$seed_id2}} length($read)." READ\n";
													print {$filehandle{$seed_id2}} $pos_snp_tmp." POS\n";
													print {$filehandle{$seed_id2}} $last_10." LAST10\n";
													
													if ($found_pos2 eq "")
													{
														$found_pos2 = $pos_snp_tmp;
													}
												}
											}		
										}
									}
									if ($found_pos2 ne "")
									{
										$read = substr $read, 0, $found_pos2-800;
										$best_extension = "";
										if ($hap_tag eq "HAP1")
										{
											undef %exclude_reads_hap1_NP;
										}
										elsif ($hap_tag eq "HAP2")
										{
											undef %exclude_reads_hap2_NP;
										}
										undef %save_alignment_data_NP;
										$position = length($read);
										$position{$id} = $position;
										$seed{$id} = $read;
									}
									if ($best_extension eq "" && ($found_pos ne "" || $found_pos2 ne ""))
									{
										print {$filehandle{$seed_id2}} $position." HALLE5\n";
										goto NP_READS;
									}

									my %score_by_nuc;
									undef %score_by_nuc;
									foreach my $rank_pattern_tmp13 (sort {$a <=> $b} keys %rank_pattern)
									{
										if (keys %{$rank_pattern{$rank_pattern_tmp13}} eq '1')
										{
											foreach my $nuc_tmp14 (keys %{$rank_pattern{$rank_pattern_tmp13}})
											{
												if (keys %{$rank_pattern{$rank_pattern_tmp13}{$nuc_tmp14}} > 3)
												{
													my $count_tmp = '0';
													my $score_tmp = '0';
													foreach my $pos_snp_tmp (keys %{$rank_pattern{$rank_pattern_tmp13}{$nuc_tmp14}})
													{
														$count_tmp++;
														$score_tmp += $rank_pattern{$rank_pattern_tmp13}{$nuc_tmp14}{$pos_snp_tmp};
													}
													if ($count_tmp > 4 && $score_tmp > 5)
													{
														print {$filehandle{$seed_id2}} $rank_pattern_tmp13." RANK_PATTERN_TEST\n";
														print {$filehandle{$seed_id2}} $nuc_tmp14." NUC ".$count_tmp." COUNT ".$score_tmp." SCORE\n";
														$score_by_nuc{$nuc_tmp14} += $score_tmp;
													}
												}
											}
										}
									}
									my $nuc_to_remove = "";
									if (keys %score_by_nuc eq '1')
									{ 
										foreach my $nuc_tmp41 (keys %score_by_nuc)
										{
											if (($score_by_nuc{$nuc_tmp41} > 18 || ($N > 10 && $N_resolved > length($best_extension)*0.08 && $score_by_nuc{$nuc_tmp41} > 12)) && $split_patterns_final{'1'} ne $nuc_tmp41)
											{
												foreach my $rank_tmp8 (sort {$a <=> $b} keys %split_patterns_final)
												{
													if ($split_patterns_final{$rank_tmp8} eq $nuc_tmp41)
													{
														$reads_to_remove{$rank_tmp8}{$rank_tmp8} = undef;
														$remove_reads = "yes";
														$nuc_to_remove = $nuc_tmp41;
														if ($rank_tmp8 eq '1')
														{
															$remove_reads = "";
															undef %reads_to_remove;
															last;
														}
													}
												}
											}
										}
									}
									elsif ($post_pattern_match eq "yes3" && $N > 10 && $N_resolved > length($best_extension)*0.08)
									{
SCORE_BY_NUC_NP:     					foreach my $nuc_tmp41 (keys %score_by_nuc)
										{
											if ($score_by_nuc{$nuc_tmp41} > 15)
											{
												foreach my $nuc_tmp41b (keys %score_by_nuc)
												{
													if ($nuc_tmp41 ne $nuc_tmp41b && $score_by_nuc{$nuc_tmp41} > $score_by_nuc{$nuc_tmp41b}*5 && $nuc_to_remove eq "")
													{
														$nuc_to_remove = $nuc_tmp41;
													}
													elsif ($nuc_tmp41 ne $nuc_tmp41b && $score_by_nuc{$nuc_tmp41} > $score_by_nuc{$nuc_tmp41b}*5)
													{
														$nuc_to_remove = "";
														last SCORE_BY_NUC_NP;
													}
												}
											}
										}
										if ($split_patterns_final{'1'} ne $nuc_to_remove)
										{
											foreach my $rank_tmp8 (sort {$a <=> $b} keys %split_patterns_final)
											{
												if ($split_patterns_final{$rank_tmp8} eq $nuc_to_remove)
												{
													$reads_to_remove{$rank_tmp8}{$rank_tmp8} = undef;
													$remove_reads = "yes";
													if ($rank_tmp8 eq '1')
													{
														$remove_reads = "";
														undef %reads_to_remove;
														last;
													}
												}
											}
										}
									}
																			
									if ($remove_reads eq "yes")
									{
										my $reads_to_remove = '0';
										foreach my $rank_tmp0 (keys %reads_to_remove)
										{
											foreach my $rank_tmp2 (keys %{$reads_to_remove{$rank_tmp0}})
											{
												$reads_to_remove++;
											}
										}
										if ($reads_to_remove < $total_nuc_count*0.65 && ($total_nuc_count - $reads_to_remove >= 4 || $post_pattern_match ne "yes3"))
										{
											print {$filehandle{$seed_id2}} $nuc_to_remove." NUC ".$score_by_nuc{$nuc_to_remove}." REMOVE_BY_RANK_PATTERN\n";
											foreach my $rank_tmp0 (keys %reads_to_remove)
											{
												foreach my $rank_tmp2 (keys %{$reads_to_remove{$rank_tmp0}})
												{
													my $id_tmp0 = $rank_to_id{$rank_tmp2};
													delete $read_start_pos_rej{$id_tmp0};
												}
											}
											goto REMOVE_READS_NP;
										}
										else
										{
											$remove_reads = "";
											undef %reads_to_remove;
										}
									}
								}
							}
						}
                        
                        my $time_split3 = time;
     my $testy_time2 = $time_split2-$time_split3;
     if ($testy_time2 > 0)
    {
        print {$filehandle{$seed_id2}} $testy_time2." TIME_CHECK_SPLIT2\n";
    }
                    }
					if ($post_pattern_match eq "yes3" && $full_reset_NP eq "" && $find_haps_in_seed eq "")
                    {
                        if ($hap_tag eq "HAP1")
						{
							undef %exclude_reads_hap1_NP;
						}
						elsif ($hap_tag eq "HAP2")
						{
							undef %exclude_reads_hap2_NP;
						}
						undef %save_alignment_data_NP;
						undef %rejected_alignment_data_NP;
						$full_reset_NP = $best_extension;
						print {$filehandle{$seed_id2}} $best_extension." FULL_RESET\n";
						if ($NP_reads_support ne "")
						{
							$NP_reads_support = "yes";
							$NP_reads_support_SNR = "";
						}
						$best_extension = "";
						$seed_id = $id;
						$y++;
						$y{$id} = $y;
						goto FULL_RESET;
                    }
					elsif ($post_pattern_match eq "yes3" && ($N_resolved > length($best_extension)*0.08 || $local_pattern_matches2 > 8) && $trace_back_check eq "" && $post_pattern_match_count > 2 && $find_haps_in_seed eq "" && $PB_reads eq "" && $input_reads_DB_folder_PB)
                    {					
						print {$filehandle{$seed_id2}} $best_extension." LAST_EXTENSION\n";
						$unresolvable_NP = "yes";
					
						last INPUT_MAFFT3_NP;
                    }
					
					my $best_extension_part_tmp = substr $best_extension, -2000;
					my $CG = $best_extension_part_tmp =~ tr/CGN/CGN/;
					my $CG_rich = "";
					my $AF0 =0.55;
					my $AF2 = 0.7;
					if ($CG > 0.52*(length($best_extension_part_tmp)))
					{
						$CG_rich = "yes";
						$AF2 = 0.64;
					}
					if ($clipped_ext ne "yes")
					{
						$AF2 = 0.64;
					}
					if ($PB_reads ne "" || $input_reads_DB_folder_PB ne "")
					{
						$AF0 = 0.75;
					}
BASECALL2_NP:					
                    if ($found_haps_in_seed eq "" && ($clipped_ext ne "yes" || $best_extension eq "" || length($best_extension_part) < $length_extension_part*0.8 || ($SNR_read_ahead ne "" && $N < 4)
						|| $N < 5 || $N < length($best_extension)*0.1 || $NP_reads_support eq "yes2" || ($longer_extension_for_repeat ne "" && $N < length($best_extension)*0.15) || $total_nuc_count < 6 || ($CG_rich eq "yes" && $N < length($best_extension)*0.2)))
                    {
                        my $last3 = substr $best_extension, -3, 3;
						$N_resolved++;
						if (($nucs{"a"} > $total_nuc_count*$AF0 && ($nucs{"-"} > $total_nuc_count*0.35 || $nucs{"a"} > $total_nuc_count*$AF2
                         || ($nucs{"a"} > 2 && $nucs{"c"} < $total_nuc_count*0.25 && $nucs{"g"} < $total_nuc_count*0.25 && $nucs{"t"} < $total_nuc_count*0.25)))
                            || ($nucs{"a"} > $total_nuc_count*0.35 && $last3 eq "AAA" && $nucs{"-"} > $total_nuc_count*0.49)
							|| ($nucs{"a"} > $total_nuc_count*0.35 && (($last3 eq "CCC" && $nucs{"c"} > $total_nuc_count*0.35) || ($last3 eq "TTT" && $nucs{"t"} > $total_nuc_count*0.35) || ($last3 eq "GGG" && $nucs{"g"} > $total_nuc_count*0.35))))
                        {
                            $best_extension .= "A";
							$best_extension_part .= "A";
                            if ($total_nuc_count > 2)
                            {
                                $quality_scores_tmp{length($best_extension)} = ($nucs{"a"}+($nucs{"-"}/2))/$total_nuc_count." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
                            }
                            #print {$filehandle{$seed_id2}} length($best_extension)." POS_QUAL ".$quality_scores_tmp{length($best_extension)}." QUAL ".$nucs{"a"}." A ".$nucs{"c"}." C ".$nucs{"t"}." T ".$nucs{"g"}." G ".$nucs{"-"}." GAP\n";
                            $nuc_match = "A";
                        }
                        elsif (($nucs{"c"} > $total_nuc_count*$AF0 && ($nucs{"-"} > $total_nuc_count*0.35 || $nucs{"c"} > $total_nuc_count*$AF2 ||
                              ($nucs{"c"} > 2 && $nucs{"a"} < $total_nuc_count*0.25 && $nucs{"g"} < $total_nuc_count*0.25 && $nucs{"t"} < $total_nuc_count*0.25)))
                               || ($nucs{"c"} > $total_nuc_count*0.35 && $last3 eq "CCC" && $nucs{"-"} > $total_nuc_count*0.49)
								|| ($nucs{"c"} > $total_nuc_count*0.35 && (($last3 eq "AAA" && $nucs{"a"} > $total_nuc_count*0.35) || ($last3 eq "TTT" && $nucs{"t"} > $total_nuc_count*0.35) || ($last3 eq "GGG" && $nucs{"g"} > $total_nuc_count*0.35))))
                        {
                            $best_extension .= "C";
							$best_extension_part .= "C";
                            if ($total_nuc_count > 2)
                            {
                                $quality_scores_tmp{length($best_extension)} = ($nucs{"c"}+($nucs{"-"}/2))/$total_nuc_count." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
                            }
                            #print {$filehandle{$seed_id2}} length($best_extension)." POS_QUAL ".$quality_scores_tmp{length($best_extension)}." QUAL ".$nucs{"a"}." A ".$nucs{"c"}." C ".$nucs{"t"}." T ".$nucs{"g"}." G ".$nucs{"-"}." GAP\n";
                            $nuc_match = "C";
                        }
                        elsif (($nucs{"t"} > $total_nuc_count*$AF0 && ($nucs{"-"} > $total_nuc_count*0.35 || $nucs{"t"} > $total_nuc_count*$AF2  ||
                            ($nucs{"t"} > 2 && $nucs{"c"} < $total_nuc_count*0.25 && $nucs{"g"} < $total_nuc_count*0.25 && $nucs{"a"} < $total_nuc_count*0.25)))
                            || ($nucs{"t"} > $total_nuc_count*0.35 && $last3 eq "TTT" && $nucs{"-"} > $total_nuc_count*0.49)
							|| ($nucs{"t"} > $total_nuc_count*0.35 && (($last3 eq "CCC" && $nucs{"c"} > $total_nuc_count*0.35) || ($last3 eq "AAA" && $nucs{"a"} > $total_nuc_count*0.35) || ($last3 eq "GGG" && $nucs{"g"} > $total_nuc_count*0.35))))
                        {
                            $best_extension .= "T";
							$best_extension_part .= "T";
                            if ($total_nuc_count > 2)
                            {
                                $quality_scores_tmp{length($best_extension)} = ($nucs{"t"}+($nucs{"-"}/2))/$total_nuc_count." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
                            }
                            #print {$filehandle{$seed_id2}} length($best_extension)." POS_QUAL ".$quality_scores_tmp{length($best_extension)}." QUAL ".$nucs{"a"}." A ".$nucs{"c"}." C ".$nucs{"t"}." T ".$nucs{"g"}." G ".$nucs{"-"}." GAP\n";
                            $nuc_match = "T";
                        }
                        elsif (($nucs{"g"} > $total_nuc_count*$AF0 && ($nucs{"-"} > $total_nuc_count*0.35 || $nucs{"g"} > $total_nuc_count*$AF2 ||
                              ($nucs{"g"} > 2 && $nucs{"c"} < $total_nuc_count*0.25 && $nucs{"a"} < $total_nuc_count*0.25 && $nucs{"t"} < $total_nuc_count*0.25)))
                            || ($nucs{"g"} > $total_nuc_count*0.35 && $last3 eq "GGG" && $nucs{"-"} > $total_nuc_count*0.49)
							|| ($nucs{"g"} > $total_nuc_count*0.35 && (($last3 eq "CCC" && $nucs{"c"} > $total_nuc_count*0.35) || ($last3 eq "TTT" && $nucs{"t"} > $total_nuc_count*0.35) || ($last3 eq "AAA" && $nucs{"a"} > $total_nuc_count*0.35))))
                        {
                            $best_extension .= "G";
							$best_extension_part .= "G";
                            if ($total_nuc_count > 2)
                            {
                                $quality_scores_tmp{length($best_extension)} = ($nucs{"g"}+($nucs{"-"}/2))/$total_nuc_count." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
                            }
                            #print {$filehandle{$seed_id2}} length($best_extension)." POS_QUAL ".$quality_scores_tmp{length($best_extension)}." QUAL ".$nucs{"a"}." A ".$nucs{"c"}." C ".$nucs{"t"}." T ".$nucs{"g"}." G ".$nucs{"-"}." GAP\n";
                            $nuc_match = "G";
                        }
						elsif (($nucs{"-"} > $total_nuc_count*$AF0 || ($nucs{"-"} > $total_nuc_count_original*0.55 && $nucs{"c"} < $total_nuc_count*0.25 && $nucs{"g"} < $total_nuc_count*0.25 &&
                                $nucs{"t"} < $total_nuc_count*0.25 && $nucs{"a"} < $total_nuc_count*0.25)))
                        {
							if ($nucs{"-"}/$total_nuc_count < 0.8)
							{
								$quality_scores_gap_tmp{length($best_extension)} = $nucs{"-"}/$total_nuc_count." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
							}
                            $nuc_match = "";
                        }
                        else
                        {
                            $N++;
                            $best_extension .= "N";
							$best_extension_part .= "N";
                            $quality_scores_tmp{length($best_extension)} = '0'." ".$nucs{"a"}." ".$nucs{"c"}." ".$nucs{"t"}." ".$nucs{"g"}." ".$nucs{"-"};
                            #print {$filehandle{$seed_id2}} length($best_extension)." POS_QUAL ".$quality_scores_tmp{length($best_extension)}." QUAL ".$nucs{"a"}." A ".$nucs{"c"}." C ".$nucs{"t"}." T ".$nucs{"g"}." G ".$nucs{"-"}." GAP\n";
                            $nuc_match = "N";
                        }
					}
                    else
                    {
                        print {$filehandle{$seed_id2}} $nuc_match." TERMINATE EARLY2\n";
						goto AFTER_NEXT_MAFFT;
                    }
                                    
SKIP_INPUT_BLAST3_NP:                   
                    if ($hap_position eq "" && $find_haps_in_seed ne "")
                    {
                        $extensions_seed{"HAP1"} .= $nuc_match;
                        $extensions_seed{"HAP2"} .= $nuc_match;
                    }
                    $hap_position = "";

                    #if (($N > length($best_extension)*0.08 || ($N > length($best_extension)*0.045 && length($best_extension) > 1000)) && $N > 5 && $confirmed_reads_count_NP > 4 && $skip_confirmed eq "" && $only_confirmed eq "yes")
                   # {
                        #print {$filehandle{$seed_id2}} $best_extension." EXT_remove\n";
                        #$best_extension = "";
                        #$skip_confirmed = "yes";
                       # undef %quality_scores_tmp;
                        
                        #goto SKIP_CONFIRMED_NP;
                    #}
                    foreach my $subject_rank (keys %subject_list)
                    {
                        if (exists($length_ext{$subject_rank}))
                        {
                            if ($length_ext{$subject_rank} < $track_length_ext_total{$subject_rank}+$track_length_ext{$subject_rank}+100)
                            {
                                delete $subject_list{$subject_rank}
                            }
                            if ($length_ext{$subject_rank} < length($best_extension))
                            {
                                delete $subject_list{$subject_rank}
                            }
							
							if ($track_length_ext{$subject_rank} > $length_extension_part)
                            {
                                $end_this_mafft_part = "yes";
                            }
                        }
                    }
                    $nuc_prev = $nuc_match;
                    $cp++;
                }
#Go to the next mafft consensus------------------------------------------------------------------------------------------------------

				my $time_CONS3 = time;
                my $time10 = $time_CONS3 - $time_BLAST3;
                print {$filehandle{$seed_id2}} $time10." TIME_CONS\n\n";
				
				if (length($best_extension) < $length_extension-$extension_part_length && length($best_extension) > 0 && $unresolvable_NP eq "" && $loop_check eq "yes")
				{
					if ($mafft_count < 3)
					{
						foreach my $rank_tmp3 (keys %track_length_ext)
						{
							$track_length_ext_total{$rank_tmp3} = $track_length_ext{$rank_tmp3};
						}
					}
					else
					{
						foreach my $rank_tmp3 (keys %track_length_ext_total)
						{
							$track_length_ext_total{$rank_tmp3} += $track_length_ext{$rank_tmp3};
						}
					}
					print {$filehandle{$seed_id2}} length($best_extension)." NEXT_MAFFT\n";   
					goto MAFFT_NP;
				}
AFTER_NEXT_MAFFT:			
#remove bad alignments------------------------------------------------------------------------------------------               
                if ($find_haps_in_seed eq "" && keys %extensions > 5 && keys %extensions > $sequencing_depth_NP/2 && $unresolvable_NP eq "")
                {
                    print {$filehandle{$seed_id2}} length($best_extension)." EXT_LENGTH\n";                  
                    my $highest_score = '0';
                    my $count_remaining_reads = '1';
                    my %best_read_score2;
                    undef %best_read_score2;
                    
                    foreach my $read_numb_tmp (sort {$a <=> $b} keys %best_read_score)
                    {
                        my $length_align = $track_length_ext_total{$read_numb_tmp};

                        if ($length_align > 0)
                        {
                            my $score_by_length = $best_read_score{$read_numb_tmp}/$length_align;
                            print {$filehandle{$seed_id2}} $read_numb_tmp." = ".$score_by_length." ALIGNMENT_SCORES\n";
                            $best_read_score2{$score_by_length} = $read_numb_tmp;
                        }
                        $count_remaining_reads++;
                    }
                    my $gg = '0';
                    foreach my $score_tmp (sort {$b <=> $a} keys %best_read_score2)
                    {
                        $gg++;
                        if ($gg eq '2')
                        {
                            $highest_score = $score_tmp;
                        }
                    }
                    
                    my $removed_reads = '0';
                    foreach my $score_tmp (sort {$a <=> $b} keys %best_read_score2)
                    {
                        my $rank_tmp9 = $best_read_score2{$score_tmp};

                        if (($score_tmp < $highest_score*0.7 || $score_tmp eq "") && $count_remaining_reads > 2 && $count_remaining_reads > 0.5*($count_remaining_reads+$removed_reads)
							&& ($rank_tmp9 > 4 || ($score_tmp < $highest_score*0.7)))
                        {
                            $count_remaining_reads--;
                            $removed_reads++;
                            $ignore_reads{$rank_tmp9} = undef;
                            my $id_tmp0 = $rank_to_id{$rank_tmp9};
                            delete $extensions2_tmp{$id_tmp0};
                            delete $extensions{$extensions2{$id_tmp0}};
                            delete $extensions2{$id_tmp0};
                            delete $extensions_nomatch{$extensions2{$id_tmp0}};
                            delete $extensions_nomatch2{$id_tmp0};
                            delete $extensions_unknown{$extensions2{$id_tmp0}};
                            delete $extensions_unknown2{$id_tmp0};
                            delete $save_reads_for_next{$id_tmp0};
                            delete $add_rej_reads_extra{$id_tmp0};
                            $ext2_count = keys %extensions2_tmp;
                        }
                    }
                    if ($removed_reads > 0 && $removed_reads > 0.08*($count_remaining_reads+$removed_reads) && ($post_pattern_match ne "yes3" || $removed_reads < 3) && ($removed_reads < 0.3*($count_remaining_reads+$removed_reads) ||
                        ($ext2_count > 0 && $removed_reads < 0.5*($count_remaining_reads+$removed_reads)) || $removed_reads < 3 ||
                        ($count_remaining_reads > 4 && $count_matches_with_high_scores > $sequencing_depth_NP*5) || length($best_extension) < 500) && $count_remaining_reads > 4)
                    {		
                        #if ($confirmed_reads_count_NP > 4 && $skip_confirmed eq "" && $confirmed_reads_count_NP-$removed_reads < 4 && $sequencing_depth_NP > 12 && $only_confirmed eq "yes")
                        #{
                            #$skip_confirmed = "yes";
                            #print {$filehandle{$seed_id2}} "REMOVE_BAD_ALIGNMENTS SKIP_CONFIRMED\n";
                            #goto SKIP_CONFIRMED_NP;
                        #}
                        $mismatch_retry++;
                        undef %quality_scores_tmp; #CHECK THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        if ($removed_reads eq '1000000000')
                        {
                            print {$filehandle{$seed_id2}} "REMOVE_BAD_ALIGNMENTS0\n";
							if ($mafft_count < 3)
							{
								$best_extension = "";
							}
							else
							{
								substr $best_extension, -length($best_extension_part), length($best_extension_part), "";
							}
							$best_extension_part = "";
                            goto IGNORE_REMOVED_READS_NP;
                        }
						$best_extension = "";
                        print {$filehandle{$seed_id2}} $removed_reads." REMOVE_BAD_ALIGNMENTS\n";  
                        goto MISMATCH_RETRY_NP;
                    }
                    elsif ($first_split_pos > 400 && $post_pattern_match eq "yes3")
                    {
                    }
                    elsif (($post_pattern_match eq "yes3" || ($removed_reads > 0.3*($count_remaining_reads+$removed_reads) && $total_nuc_count > 5 && length($best_extension) < 350)) && $full_reset_NP eq "")
                    {
                        if ($hap_tag eq "HAP1")
						{
							undef %exclude_reads_hap1_NP;
						}
						elsif ($hap_tag eq "HAP2")
						{
							undef %exclude_reads_hap2_NP;
						}
						undef %save_alignment_data_NP;
						$full_reset_NP = $best_extension;
						$best_extension = "";
						print {$filehandle{$seed_id2}} $full_reset_NP." FULL_RESET2\n";  
						goto NP_READS;
                    }
                    elsif ($post_pattern_match eq "yes3" || ($removed_reads > 0.3*($count_remaining_reads+$removed_reads) && $total_nuc_count > 5 && length($best_extension) < 350))
                    {
                        $best_extension = "";                      
                        print {$filehandle{$seed_id2}} $last_non_complex_region{$seed_id}." LAST_NON_COMPLEX\n";  
                        goto AFTER_EXT;
                    }
                }
			
 #-------------------------------------------------------------------------------------------------------------------------

                if ($found_haps_in_seed ne "yes" && ($clipped_ext eq "yes" || ($find_haps_in_seed eq "yes" && $best_extension ne "")))
                {
					if ($unresolvable_NP ne "" && length($full_reset_NP) > length($best_extension)+1 && $split_contigs_NP eq "")
					{
						$best_extension = $full_reset_NP;
					}
					
					my $last_70b = substr $best_extension, -70, 70;
					my $N_test = $last_70b =~ tr/N/N/;
					while ($N_test > 5 && length($best_extension) > 500)
					{
						substr $best_extension, -100,100, "";  
						$last_70b = substr $best_extension, -70, 70;
						$N_test = $last_70b =~ tr/N/N/;
					}
                    print {$filehandle{$seed_id2}} $best_extension." EXT0\n";
                }
                elsif ($found_haps_in_seed eq "yes")
                {
                    print {$filehandle{$seed_id2}} $best_extension." BEST_EXT_FOUND_HAPS\n"; 
                    $best_extension = "";
                    delete $seed{$id};
                    $compare_haps = "yes2";
                    foreach my $id_tmp (keys %extensions_seed)
                    {
                        $seed{$id."_".$id_tmp} = $extensions_seed{$id_tmp};
						if ($find_haps_in_seed eq "yes2")
                        {
                            $seed{$id."_".$id_tmp} = $read.$extensions_seed{$id_tmp};
                        }
						$original_seed_length{$id."_".$id_tmp} = length($seed{$id."_".$id_tmp});
                        print {$filehandle{$seed_id2}} $extensions_seed{$id_tmp}." EXT_SEED\n"; 
                        
                        $hap_compare_pos{$id."_".$id_tmp} = 0;
                        delete $find_haps_in_seed{$id};
                        $position{$id."_".$id_tmp} = length($seed{$id."_".$id_tmp});

						my $fh;
						$filehandle{$id."_".$id_tmp} = $fh;
						$output_file5  = $output_path."log_extended_".$project."_".$id."_".$id_tmp.".txt";
						open($filehandle{$id."_".$id_tmp}, ">".$output_file5) or die "Can't open file $output_file5, $!\n";
						
						my $fh2;
						$filehandle4{$id."_".$id_tmp} = $fh2;
						$output_file13 = $output_path."quality_scores_".$project."_".$id."_".$id_tmp.".txt";
						open($filehandle4{$id."_".$id_tmp}, ">".$output_file13) or die "Can't open file $output_file13, $!\n";
                    }
                }
                elsif ($NP_reads_support eq "yes2" && $clipped_ext eq "")
                {        
                    $best_extension = "";
                    print {$filehandle{$seed_id2}} "BACK TO PB\n";
                    goto AFTER_EXT2;
                }
                else
                {
                    print {$filehandle{$seed_id2}} "NOCLIP\n";
                    print {$filehandle{$seed_id2}} $best_extension." EXT\n";
                    $best_extension = "";
                    delete $seed{$id};
                }
                
                my $cut_repeat_seq = "";
                
                if (length($best_extension) > 150 && ($clipped_ext eq "yes" || $found_haps_in_seed ne "yes") && $unresolvable_NP eq "")
                {
                    my $best_extension_original = length($best_extension);
CUT_AGAIN_NP:       my $h = '70';
                    my $last_70 = substr $best_extension, -$h, 70;
                    my $AT_rich_tmp = AT_rich_test ($last_70, '10');
                    
                    while (length($best_extension) > $h+90 && $AT_rich_tmp eq "yes")
                    {
                        $h += 70;
                        $last_70 = substr $best_extension, -$h, 70;
                        $AT_rich_tmp = AT_rich_test ($last_70, '10');
                    }
                    
                    if ($h > 70)
                    {
                        $h -= 60;
                        if (length($best_extension) > $h+150)
                        {
                            substr $best_extension, -$h, $h, "";
                            $cut_repeat_seq{$seed_id} = undef;
                            $cut_repeat_seq = "yes";
                            print {$filehandle{$seed_id2}} $best_extension." CUT_REPEAT_SEQ0\n";
                            if (length($best_extension) > 150)
                            {
                                goto CUT_AGAIN_NP;
                            }   
                        }
                        elsif ($longer_extension_for_repeat eq "" && $best_extension_original < 3900)
                        {
                            $longer_extension_for_repeat = $best_extension_original+3000;
                            $length_extension = 3900;
                            print {$filehandle{$seed_id2}} $longer_extension_for_repeat." LONGER_EXTENSIONS0\n";
                            $best_extension = "";
                            #$high_quality = "yes";
                            undef %quality_scores_tmp;
                            goto REDUCE_EXTENSIONS_NP;
                        }
                    }
                    my $h3 = '0';
                    my $h3_limit = '30';
                    my $check_rep15_full = "";
CUT_AGAIN_NP2:                   
                    my $h2 = '15';
                    my $last_15 = substr $best_extension, -$h2-$h3, 15;
                    my $last_250 = substr $best_extension, -250, 250;
                    my $check_rep15 = $last_250 =~ s/$last_15/$last_15/g;
                    my $N_check = $last_15 =~ tr/N/N/;
                    my $best_extension_tmp = $best_extension;
                    #my $check_rep15_full = $best_extension_tmp =~ s/$last_15/\+/g;
                 
                    my @check_rep15_full = split /$last_15/, $best_extension_tmp;
                    my %lengths_tmp;
                    undef %lengths_tmp;
                    
                    if (@check_rep15_full > 2)
                    {
                        foreach my $seq_tmp (@check_rep15_full)
                        {
                            if (exists($lengths_tmp{length($seq_tmp)}))
                            {
                                my $count_tmp = $lengths_tmp{length($seq_tmp)}+1;
                                $lengths_tmp{length($seq_tmp)} = $count_tmp;
                            }
                            else
                            {
                                 $lengths_tmp{length($seq_tmp)} = '1';
                            }                    
                        }
                        my $score_tmp = '0';
                        my $prev_length = "";
                        
                        foreach my $length_tmp (sort {$a <=> $a} keys %lengths_tmp)
                        {
                            if ($prev_length eq "" || $prev_length < $length_tmp-150 || $prev_length < $length_tmp*0.92)
                            {
                                $score_tmp = $lengths_tmp{$length_tmp};
                            }
                            else
                            {
                                $score_tmp += $lengths_tmp{$length_tmp};
                            }
                            
                            if ($score_tmp > 1)
                            {
                                $check_rep15_full = "yes";
                                if (length($check_rep15_full[0]) > $length_tmp)
                                {
                                    $check_rep15_full = "yes2";
                                }
                                print {$filehandle{$seed_id2}} $score_tmp." LONG_REP_TEST\n";
                            }
                        }
                    }
                    
                    while ($check_rep15 > 1 && length($best_extension) > $h2+$h3+50)
                    {
                        $h2 +=230;
                        $last_15 = substr $best_extension, -$h2-$h3, 15;
                        $last_250 = substr $best_extension, -250-$h2+15, 250;
                        $check_rep15 = $last_250 =~ s/$last_15/$last_15/g;
                    }
                    
                    if ($h2 > 15 || $check_rep15_full eq "yes")
                    {
                        $h2 -= 15;
                        
                        if ($longer_extension_for_repeat <= $length_extension)
                        {
                            $longer_extension_for_repeat += 3000;
                            print {$filehandle{$seed_id2}} $longer_extension_for_repeat." LONGER_EXTENSIONS1a\n";
                            $best_extension = "";
                            #$high_quality = "yes";
                            undef %quality_scores_tmp;
                            goto REDUCE_EXTENSIONS_NP;
                        }
                        elsif (length($best_extension) > $h2+$h3+250 && $h2 > 15)
                        {
                            substr $best_extension, -$h2-$h3-100, $h2+$h3+100, "";
                            $cut_repeat_seq{$seed_id} = undef;
                            $cut_repeat_seq = "yes";
                            print {$filehandle{$seed_id2}} $best_extension." CUT_REPEAT_SEQ1\n";
                            if (length($best_extension) > 150)
                            {
                                goto CUT_AGAIN_NP;
                            }                  
                        }
                        elsif ($check_rep15_full eq "yes2")
                        {
                            my $tmpie = length($check_rep15_full[0])-250;
                            if ($tmpie < 100)
                            {
                                $tmpie = '100';
                            }
                            substr $best_extension, $tmpie, length($best_extension)-$tmpie, "";
                            $cut_repeat_seq{$seed_id} = undef;
                            $cut_repeat_seq = "yes";
                            print {$filehandle{$seed_id2}} $best_extension." CUT_REPEAT_SEQ2\n";                 
                        }
                        elsif ($longer_extension_for_repeat eq "" && $best_extension_original < 3900)
                        {
                            $longer_extension_for_repeat = $best_extension_original+3000;
                            print {$filehandle{$seed_id2}} $length_extension." LONGER_EXTENSIONS1\n";
                            $best_extension = "";
                            #$high_quality = "yes";
                            undef %quality_scores_tmp;
                            goto REDUCE_EXTENSIONS_NP;
                        }
                    }
                    elsif ($N_check > 0)
                    {
                        $h3 += 5;
                        $h3_limit += 5;
                        goto CUT_AGAIN_NP2;
                    }
                    elsif ($h3 < $h3_limit)
                    {
                        $h3 += 5;
                        goto CUT_AGAIN_NP2;
                    }
                }
                if ($first_split_pos > 500 && $post_pattern_match eq "yes3" && $unresolvable_NP eq "")
                {            
                    my $ext_new_tmp = substr $best_extension, 0, $first_split_pos-100;
                    my $h = '70';
                    my $last_70 = substr $ext_new_tmp, -$h, 70;
                    my $AT_rich_tmp = AT_rich_test ($last_70, '10');
                    
                    if ($AT_rich_tmp eq "yes")
                    {
                        $best_extension = $ext_new_tmp;
                        print {$filehandle{$seed_id2}} $best_extension." EXT0a\n";
                    }
                }
END_NP:                
                
                my $time_CONSENSUS = time;
                my $time11 = $time_CONSENSUS - $time_BLAST3b;
                print {$filehandle{$seed_id2}} $time11." TIME_CONSENSUS\n";
                
                if ($cut_repeat_seq eq "")
                {
                    delete $cut_repeat_seq{$seed_id};
                }
                
                if ($clipped_ext eq "yes" && $best_extension ne "")
                {
                    foreach my $pos_tmp (sort {$a <=> $b} keys %quality_scores_tmp)
                    {
                        if ($pos_tmp <= length($best_extension))
                        {
                            #my $nuc_tmp = substr $best_extension, $pos_tmp-1, 1;
                            $quality_scores{$id}{$position+$pos_tmp} = $quality_scores_tmp{$pos_tmp};
                            #print {$filehandle{$seed_id2}} $position+$pos_tmp." POS_QUAL ".$quality_scores{$id}{$position+$pos_tmp}." QUAL\n";
                        }   
                    }
                    foreach my $pos_tmp (sort {$a <=> $b} keys %quality_scores_gap_tmp)
                    {
                        if ($pos_tmp <= length($best_extension))
                        {
                            $quality_scores_gap{$id}{$position+$pos_tmp} = $quality_scores_gap_tmp{$pos_tmp};
                            #print {$filehandle{$seed_id2}} $pos_tmp." POS_QUAL\n";
                            #print {$filehandle{$seed_id2}} $quality_scores{$id}{$position+$pos_tmp}." QUAL\n";
                        }   
                    }
                }
                close INPUT_BLAST3;
				
#Store alignment data for the next iteration---------------------------------------------------------------------------------------------------

                if ($find_haps_in_seed eq "" && $unresolvable_NP eq "")
                {
                    foreach my $seed_id_tmp5 (keys %save_alignment_data_NP)
                    {
                        if ($seed_id_tmp5 eq $seed_id)
                        {
                            foreach my $id_tmp4 (keys %{$save_alignment_data_NP{$seed_id_tmp5}})
                            {
                                my @alignment_data = split /_/, $save_alignment_data_NP{$seed_id_tmp5}{$id_tmp4};
                                if ($alignment_data[0] eq "yes" || (($alignment_data[3]+$alignment_data[10]) < $position))
                                {
                                    delete $save_alignment_data_NP{$seed_id_tmp5}{$id_tmp4};
                                }
                                elsif ($alignment_data[0] eq "no" && (($alignment_data[3]+9000 < $position) || (($alignment_data[4]+$alignment_data[3]+5000) < $position)))
                                {
                                    delete $save_alignment_data_NP{$seed_id_tmp5}{$id_tmp4};
                                }
                            }
                        }
                    }
					undef %rejected_alignment_data_NP;
                    foreach my $id_tmp5 (keys %long_read_end_pos_save)
                    {
                        my $assembled = "no";
                        my $reverse_tmp = "no";
                        my $alignment_length_tmp = '0';
                        my $score_matches_tmp = "___";
                        if (exists($save_reads_for_next{$id_tmp5}))
                        {
                            $assembled = "yes";
                        }
                        if (exists($reverse_list{$id_tmp5}))
                        {
                            $reverse_tmp = "yes";
                        }
                        if (exists($double_matches{$id_tmp5}))
                        {
                            $reverse_tmp = "yes2";
                        }
                        if (exists($alignment_length_save{$id_tmp5}))
                        {
                            $alignment_length_tmp += $alignment_length_save{$id_tmp5};
                        }
                        if (exists($score_matches_save{$id_tmp5}))
                        {
                            $score_matches_tmp = $score_matches_save{$id_tmp5};
                        }
                        $save_alignment_data_NP{$seed_id}{$id_tmp5} = $assembled."_".$reverse_tmp."_".$long_read_end_pos_save{$id_tmp5}."_".$position."_".$alignment_length_tmp
                        ."_".$score_matches_tmp."_".$accuracy{$id_tmp5}."_".$length_ext_all{$id_tmp5}."_".$read_start_pos_rej{$id_tmp5}."_".$hash_NP_reads_tmp{$id_tmp5};
    
                        #print {$filehandle{$seed_id2}} $id_tmp5." ".$save_alignment_data_NP{$seed_id}{$id_tmp5}." SAVVEEEE\n";
                    }
					foreach my $id_tmp5 (keys %rejected_reads_save)
					{
						$rejected_alignment_data_NP{$seed_id}{$id_tmp5} = undef;
					}
                }
#Post error correction -------------------------------------------------------------------------------------------------------------------------------------------------------------------------				
				my $read_count_N_check = keys %extensions2_tmp;
				if ($read_count_N_check > 5 && $NP_reads_support eq "")
				{
					my %N_mismatches_tmp;
					undef %N_mismatches_tmp;
					foreach my $id_tmp15 (keys %extensions2_tmp)
					{
						foreach my $pos_snp_tmp (sort {$a <=> $b} keys %{$store_mismatches_N_NP{$id_tmp15}})
						{
							my @store_mismatches_N_NP = split /,/, $store_mismatches_N_NP{$id_tmp15}{$pos_snp_tmp};
							$N_mismatches_tmp{$pos_snp_tmp}{$store_mismatches_N_NP[1]} += 1;
						} 
					}
					my $pos_assem_tmp = $position-1;
POST_ERROR_CORR:    while ($pos_assem_tmp > 0)
					{
						my $count_cov = '0';
						my $mismatch_count_tmp = '0';
						my $match_count_tmp = '0';
						my $N_count_tmp = '0';
						my %id_list_tmp;
						undef %id_list_tmp;
						my $corrected_check = "";
						
						foreach my $id_tmp15 (keys %extensions2_tmp)
						{
							my $overlap_pos_tmp = $position - $alignment_length_save{$id_tmp15};
							if ($pos_assem_tmp > $overlap_pos_tmp)
							{
								$count_cov++;
								if (exists($store_mismatches_all_NP{$id_tmp15}{$pos_assem_tmp}))
								{
									$mismatch_count_tmp++;
								}
								elsif (exists($store_mismatches_N_NP{$id_tmp15}{$pos_assem_tmp}))
								{
									$N_count_tmp++;
								}
								else
								{
									$match_count_tmp++;
								}
								$id_list_tmp{$id_tmp15} = undef;
							}
						}
						my $total_count_tmp = $mismatch_count_tmp+$N_count_tmp+$match_count_tmp;
						if ($total_count_tmp > 4 && $total_count_tmp > $read_count_N_check*0.3)
						{
							if (($mismatch_count_tmp > 4 && $mismatch_count_tmp > 0.9*$total_count_tmp) || ($mismatch_count_tmp > 5 && $mismatch_count_tmp > 0.75*$total_count_tmp))
							{
								my %all_mismatches_tmp;
								undef %all_mismatches_tmp;
								foreach my $id_tmp15 (keys %id_list_tmp)
								{
									if (exists($store_mismatches_all_NP{$id_tmp15}{$pos_assem_tmp}))
									{
										my @store_mismatches_all_NP = split /,/, $store_mismatches_all_NP{$id_tmp15}{$pos_assem_tmp};
										$all_mismatches_tmp{$store_mismatches_all_NP[1]} += 1;
										#print {$filehandle{$seed_id2}} $pos_assem_tmp." ".$store_mismatches_all_NP{$id_tmp15}{$pos_assem_tmp}." ALL_CORRECTED_TEST0\n";
									}
								}
								my $gap_check_count = "";
								foreach my $nuc_tmp15 (keys %all_mismatches_tmp)
								{
									if ($nuc_tmp15 eq "-" && $all_mismatches_tmp{$nuc_tmp15} > 0.25*$total_count_tmp)
									{
										$gap_check_count = "yes";
									}
								}
								foreach my $nuc_tmp15 (keys %all_mismatches_tmp)
								{
									print {$filehandle{$seed_id2}} $pos_assem_tmp." ".$nuc_tmp15." ".$all_mismatches_tmp{$nuc_tmp15}." ALL_CORRECTED_TEST0\n";
									if ($nuc_tmp15 ne "-" && (($all_mismatches_tmp{$nuc_tmp15} > 4 && $all_mismatches_tmp{$nuc_tmp15} > 0.85*$total_count_tmp) ||
										($all_mismatches_tmp{$nuc_tmp15} > 5 && $all_mismatches_tmp{$nuc_tmp15} > 0.74*$total_count_tmp && $all_mismatches_tmp{$nuc_tmp15} > $read_count_N_check*0.55)))
									{#ADDED $nuc_tmp15 ne "-" && 
										print {$filehandle{$seed_id2}} $pos_assem_tmp." ".$nuc_tmp15." ALL_CORRECTED\n";
										my $nuci_test = substr $read, $pos_assem_tmp-1, 1;
										$corrected_check = "yes";
										print {$filehandle{$seed_id2}} $all_mismatches_tmp{$nuc_tmp15}/$total_count_tmp." Nucie_Per\n";
										
										if ($nuci_test ne "-")
										{
											my $posie_tmp = $pos_assem_tmp-1;
											if ($nuc_tmp15 eq "-")
											{
												$nuc_tmp15 = "";
												my $posie_tmp2 = $posie_tmp;
												while ($posie_tmp2 < length($read))
												{
													if (exists($quality_scores{$seed_id}{$posie_tmp2+1}))
													{
														$quality_scores{$seed_id}{$posie_tmp2} = $quality_scores{$seed_id}{$posie_tmp2+1};
													}
													if (exists($quality_scores_gap{$seed_id}{$posie_tmp2+1}))
													{
														$quality_scores_gap{$seed_id}{$posie_tmp2} = $quality_scores_gap{$seed_id}{$posie_tmp2+1};
													}
													if (exists($split_positions{$seed_id}{$posie_tmp2}))
													{
														$split_positions{$seed_id}{$posie_tmp2-1} = $split_positions{$seed_id}{$posie_tmp2};
														delete $split_positions{$seed_id}{$posie_tmp2};
													}
													if (exists($split_positions_DUP{$seed_id}{$posie_tmp2}))
													{
														$split_positions_DUP{$seed_id}{$posie_tmp2-1} = $split_positions_DUP{$seed_id}{$posie_tmp2};
														delete $split_positions_DUP{$seed_id}{$posie_tmp2};
													}
													$posie_tmp2++;
												}
												$quality_scores_gap{$seed_id}{$posie_tmp} = $all_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15}/$total_count_tmp;
												$position--;
												$position{$id} = $position;
											}
											else
											{
												$quality_scores{$seed_id}{$posie_tmp} = ($all_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15}/$total_count_tmp)-0.05;
											}
											substr $read, $posie_tmp, 1, $nuc_tmp15;
										}
										else
										{
											print {$filehandle{$seed_id2}} $pos_assem_tmp." ".$nuc_tmp15." ALL_CORRECTED_ERROR\n";
											print {$filehandle{$seed_id2}} $nuci_test." Nucie_test\n";
										}
									}
								}
								print {$filehandle{$seed_id2}} $total_count_tmp." ".$pos_assem_tmp." MISMATCH_CHECK\n";
							}
							elsif ($N_count_tmp > 4)
							{
								if ($pos_assem_tmp > $position-11000)
								{
									my $q_score_line = "";
									my $line_part_tmp = "A";
									my $total_nuc_count_tmp = '0';
									my $highest_nuc_count_tmp = '0';
									my $highest_nuc_tmp = '0';
NEW_Q_LINE:											
									my $check_tmp = "";
									foreach my $nuc_tmp15 (keys %{$N_mismatches_tmp{$pos_assem_tmp}})
									{	
										if ($nuc_tmp15 eq $line_part_tmp)
										{
											$q_score_line .= " ".$N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15};
											$check_tmp = "yes";
											$total_nuc_count_tmp += $N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15};
											if ($N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15} > $highest_nuc_count_tmp)
											{
												$highest_nuc_count_tmp = $N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15};
												$highest_nuc_tmp = $nuc_tmp15;
											}
										}
									}
									
									if ($check_tmp eq "")
									{
										$q_score_line .= " 0"
									}
									if ($line_part_tmp eq "A")
									{
										$line_part_tmp = "C";
										goto NEW_Q_LINE;
									}
									if ($line_part_tmp eq "C")
									{
										$line_part_tmp = "T";
										goto NEW_Q_LINE;
									}
									if ($line_part_tmp eq "T")
									{
										$line_part_tmp = "G";
										goto NEW_Q_LINE;
									}
									if ($line_part_tmp eq "G")
									{
										$line_part_tmp = "-";
										goto NEW_Q_LINE;
									}
									if ($total_nuc_count_tmp > 0)
									{
										my $first_part = "";
										if ($highest_nuc_tmp eq "-")
										{
											$first_part = $highest_nuc_count_tmp/$total_nuc_count_tmp
										}
										else
										{
											$first_part = ($highest_nuc_count_tmp+(($N_mismatches_tmp{$pos_assem_tmp}{"-"})/2))/$total_nuc_count_tmp
										}
										$first_part .= $q_score_line;
										$quality_scores{$seed_id}{$pos_assem_tmp} = $first_part;
									}
									else
									{
										foreach my $nuc_tmp15 (keys %{$N_mismatches_tmp{$pos_assem_tmp}})
										{
											print {$filehandle{$seed_id2}} $pos_assem_tmp." ".$nuc_tmp15." ".$N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15}." N_CORRECTED0\n";
										}
									}
								}
							
								my $gap_check_count = "";
								foreach my $nuc_tmp15 (keys %{$N_mismatches_tmp{$pos_assem_tmp}})
								{
									if ($nuc_tmp15 eq "-" && $N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15} > 0.25*$total_count_tmp)
									{
										$gap_check_count = "yes";
									}
								}
								foreach my $nuc_tmp15 (keys %{$N_mismatches_tmp{$pos_assem_tmp}})
								{
									if ($nuc_tmp15 ne "-" && (($N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15} > 4 && $N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15} > 0.85*$total_count_tmp) ||
										($N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15} > 5 && $N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15} > 0.74*$total_count_tmp && $N_mismatches_tmp{$nuc_tmp15} > $read_count_N_check*0.55)))
									{#ADDED $nuc_tmp15 ne "-" && 
										print {$filehandle{$seed_id2}} $pos_assem_tmp." ".$nuc_tmp15." N_CORRECTED\n";
										my $nuci_test = substr $read, $pos_assem_tmp-1, 1;
										$corrected_check = "yes";
										print {$filehandle{$seed_id2}} $N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15}/$total_count_tmp." Nucie_Per\n";
										
										if ($nuci_test eq "N")
										{
											my $posie_tmp = $pos_assem_tmp-1;
											if ($nuc_tmp15 eq "-")
											{
												$nuc_tmp15 = "";
												my $posie_tmp2 = $posie_tmp;
												while ($posie_tmp2 < length($read))
												{
													if (exists($quality_scores{$seed_id}{$posie_tmp2+1}))
													{
														$quality_scores{$seed_id}{$posie_tmp2} = $quality_scores{$seed_id}{$posie_tmp2+1};
													}
													if (exists($quality_scores_gap{$seed_id}{$posie_tmp2+1}))
													{
														$quality_scores_gap{$seed_id}{$posie_tmp2} = $quality_scores_gap{$seed_id}{$posie_tmp2+1};
													}
													if (exists($split_positions{$seed_id}{$posie_tmp2}))
													{
														$split_positions{$seed_id}{$posie_tmp2-1} = $split_positions{$seed_id}{$posie_tmp2};
														delete $split_positions{$seed_id}{$posie_tmp2};
													}
													if (exists($split_positions_DUP{$seed_id}{$posie_tmp2}))
													{
														$split_positions_DUP{$seed_id}{$posie_tmp2-1} = $split_positions_DUP{$seed_id}{$posie_tmp2};
														delete $split_positions_DUP{$seed_id}{$posie_tmp2};
													}
													$posie_tmp2++;
												}
												$quality_scores_gap{$seed_id}{$posie_tmp} = $N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15}/$total_count_tmp;
												$position--;
												$position{$id} = $position;
											}
											else
											{
												$quality_scores{$seed_id}{$posie_tmp} = ($N_mismatches_tmp{$pos_assem_tmp}{$nuc_tmp15}/$total_count_tmp)-0.05;
											}
											substr $read, $posie_tmp, 1, $nuc_tmp15;		
										}
										else
										{
											print {$filehandle{$seed_id2}} $pos_assem_tmp." ".$nuc_tmp15." N_CORRECTED_ERROR\n";
											print {$filehandle{$seed_id2}} $nuci_test." Nucie_test\n";
										}
									}	
								}								
							}
							elsif (($match_count_tmp > 5 && $match_count_tmp > 0.9*$total_count_tmp) || ($match_count_tmp > 6 && $match_count_tmp > 0.75*$total_count_tmp))
							{
								if (exists($quality_scores{$seed_id}{$pos_assem_tmp}))
								{
									my $new_score = $match_count_tmp/$total_count_tmp;
									my @q_score_tmp = split / /, $quality_scores{$seed_id}{$pos_assem_tmp};  
									if ($new_score > $q_score_tmp[0] && $q_score_tmp[0] < 0.8)
									{
										#print {$filehandle{$seed_id2}} $quality_scores{$seed_id}{$pos_assem_tmp}." ".$new_score." NEW_SCORE\n";
										$quality_scores{$seed_id}{$pos_assem_tmp} = $new_score;
									}
								}
							}
						}
						else
						{
							last POST_ERROR_CORR;
						}
						$pos_assem_tmp--;
					}
				}	
            }
#---------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------

AFTER_EXT:
            my $time_END = time;
            my $time_ITERATION = $time_END-$time_START;
            print {$filehandle{$seed_id2}} $time_ITERATION." TIME_ITERATION\n";
            if ($time_ITERATION > 30)
            {
                print {$filehandle{$seed_id2}}"TIME_ALERT\n";
            }
            
            if (($best_extension eq "" || $unresolvable_NP eq "yes" || $unresolvable_PB eq "yes" || (length($read) > $assembly_length_max && $assembly_length_max ne "WG" && $assembly_length_max ne ""))
				&& $skip_hap ne $id && $found_haps_in_seed eq "" && $find_haps_in_seed ne "yes")
            { 
				my $output_file6  = $output_path."contigs_tmp_".$id."_".$project.".fasta";
				open($filehandle3{$id}, ">".$output_file6) or die "Can't open file $output_file6, $!\n";
				print {$filehandle3{$id}} ">".$id."\n";
				my $m = '0';
				while (length($read) > $m)
				{
					my $value_ref2b = substr $read, $m, 150;
					$m += 150;
					print {$filehandle3{$id}} $value_ref2b."\n";
				}
            	
				delete $seed{$id};
				foreach my $id_tmp (sort {$a <=> $b} keys %quality_scores)
				{
					if ($id eq $id_tmp)
					{
						foreach my $pos_qual (sort {$a <=> $b} keys %{$quality_scores{$id_tmp}})
						{
							print {$filehandle4{$seed_id2}} $pos_qual." ".$quality_scores{$id_tmp}{$pos_qual}."\n";
						}
					}
				}
				my $pos_in_ass = '1';
				my $read_polished = "";
				my $homopolymer_AF = 0.6;
				my $gap_AF = 1.3;
				if ($PB_reads ne "" || $input_reads_DB_folder_PB ne "")
				{
					$homopolymer_AF = 1;
					$gap_AF = 1.1;
				}
				while($pos_in_ass < length($read))
				{
					my $nuc_tmp = substr $read, $pos_in_ass-1, 1;
					
					if ($nuc_tmp eq "N" && exists($quality_scores{$id}{$pos_in_ass}))
					{
						my $before3 = substr $read, $pos_in_ass-4, 3;
						my $after3 = substr $read, $pos_in_ass, 3;
						my @q_score_tmp = split / /, $quality_scores{$id}{$pos_in_ass};
	print {$filehandle{$seed_id2}} $pos_in_ass." POS ".$quality_scores{$id}{$pos_in_ass}."\n";
						my $total_tmp = $q_score_tmp[1]+$q_score_tmp[2]+$q_score_tmp[3]+$q_score_tmp[4]+$q_score_tmp[5];
						
						if ($q_score_tmp[1] > $q_score_tmp[2]*1.2 && $q_score_tmp[1] > $q_score_tmp[3]*1.2 && $q_score_tmp[1] > $q_score_tmp[4]*1.2 &&
						   ($q_score_tmp[1] > $q_score_tmp[5]*0.8 || ($q_score_tmp[1] > $q_score_tmp[5]*$homopolymer_AF && ($before3 eq "AAA" || $after3 eq "AAA"))))
						{
							$read_polished .= "A";
							print {$filehandle{$seed_id2}} $pos_in_ass." A POLISHED\n";
						}
						elsif ($q_score_tmp[2] > $q_score_tmp[1]*1.2 && $q_score_tmp[2] > $q_score_tmp[3]*1.2 && $q_score_tmp[2] > $q_score_tmp[4]*1.2 && 
						      ($q_score_tmp[2] > $q_score_tmp[5]*0.8 || ($q_score_tmp[2] > $q_score_tmp[5]*$homopolymer_AF && ($before3 eq "CCC" || $after3 eq "CCC"))))
						{
							$read_polished .= "C";
							print {$filehandle{$seed_id2}} $pos_in_ass." C POLISHED\n";
						}
						elsif ($q_score_tmp[3] > $q_score_tmp[2]*1.2 && $q_score_tmp[3] > $q_score_tmp[1]*1.2 && $q_score_tmp[3] > $q_score_tmp[4]*1.2 &&
						      ($q_score_tmp[3] > $q_score_tmp[5]*0.8 || ($q_score_tmp[3] > $q_score_tmp[5]*$homopolymer_AF && ($before3 eq "TTT" || $after3 eq "TTT"))))
						{
							$read_polished .= "T";
							print {$filehandle{$seed_id2}} $pos_in_ass." T POLISHED\n";
						}
						elsif ($q_score_tmp[4] > $q_score_tmp[2]*1.2 && $q_score_tmp[4] > $q_score_tmp[3]*1.2 && $q_score_tmp[4] > $q_score_tmp[1]*1.2 && 
						      ($q_score_tmp[4] > $q_score_tmp[5]*0.8 || ($q_score_tmp[4] > $q_score_tmp[5]*$homopolymer_AF && ($before3 eq "GGG" || $after3 eq "GGG"))))
						{
							$read_polished .= "G";
							print {$filehandle{$seed_id2}} $pos_in_ass." G POLISHED\n";
						}
						elsif ($q_score_tmp[5] > $q_score_tmp[2]*$gap_AF && $q_score_tmp[5] > $q_score_tmp[3]*$gap_AF && $q_score_tmp[5] > $q_score_tmp[1]*$gap_AF && $q_score_tmp[5] > 0.55*$total_tmp)
						{
							print {$filehandle{$seed_id2}} $pos_in_ass." GAP POLISHED\n";
						}
						else
						{
							$read_polished .= "N";
						}
					}		
					else
					{
						$read_polished .= $nuc_tmp;
					}
					if (exists($quality_scores_gap{$id}{$pos_in_ass}))
					{
						my $before3 = substr $read, $pos_in_ass-4, 3;
						my $after3 = substr $read, $pos_in_ass, 3;
						my @q_score_tmp = split / /, $quality_scores_gap{$id}{$pos_in_ass};
						if ($q_score_tmp[0] < 0.75)
						{
							print {$filehandle{$seed_id2}} $pos_in_ass." POS ".$quality_scores_gap{$id}{$pos_in_ass}." GAP_SCORE\n";
							if 	($q_score_tmp[1] > $q_score_tmp[2]*1.2 && $q_score_tmp[1] > $q_score_tmp[3]*1.2 && $q_score_tmp[1] > $q_score_tmp[4]*1.2 &&
								($q_score_tmp[1] > $q_score_tmp[5]*0.8 || ($q_score_tmp[1] > $q_score_tmp[5]*0.6 && ($before3 eq "AAA" || $after3 eq "AAA"))))
							{
								$read_polished .= "A";
								print {$filehandle{$seed_id2}} $pos_in_ass." A POLISHED GAP\n";
							}
							elsif ($q_score_tmp[2] > $q_score_tmp[1]*1.2 && $q_score_tmp[2] > $q_score_tmp[3]*1.2 && $q_score_tmp[2] > $q_score_tmp[4]*1.2 && 
								  ($q_score_tmp[2] > $q_score_tmp[5]*0.8 || ($q_score_tmp[2] > $q_score_tmp[5]*0.6 && ($before3 eq "CCC" || $after3 eq "CCC"))))
							{
								$read_polished .= "C";
								print {$filehandle{$seed_id2}} $pos_in_ass." C POLISHED GAP\n";
							}
							elsif ($q_score_tmp[3] > $q_score_tmp[2]*1.2 && $q_score_tmp[3] > $q_score_tmp[1]*1.2 && $q_score_tmp[3] > $q_score_tmp[4]*1.2 &&
								  ($q_score_tmp[3] > $q_score_tmp[5]*0.8 || ($q_score_tmp[3] > $q_score_tmp[5]*0.6 && ($before3 eq "TTT" || $after3 eq "TTT"))))
							{
								$read_polished .= "T";
								print {$filehandle{$seed_id2}} $pos_in_ass." T POLISHED GAP\n";
							}
							elsif ($q_score_tmp[4] > $q_score_tmp[2]*1.2 && $q_score_tmp[4] > $q_score_tmp[3]*1.2 && $q_score_tmp[4] > $q_score_tmp[1]*1.2 && 
								  ($q_score_tmp[4] > $q_score_tmp[5]*0.8 || ($q_score_tmp[4] > $q_score_tmp[5]*0.6 && ($before3 eq "GGG" || $after3 eq "GGG"))))
							{
								$read_polished .= "G";
								print {$filehandle{$seed_id2}} $pos_in_ass." G POLISHED GAP\n";
							}
						}
					}
					$pos_in_ass++;
				}
				$read = $read_polished;
				
                if ($first_back_assembly eq "" && $y < 3 && $assembly_length_max eq "WG")
				{		
				}
				else
				{
					my $output_file7  = $output_path."assembly_".$project."_".$seed_id2.".fasta";               
					open(OUTPUT7, ">".$output_file7) or die "Can't open file $output_file7, $!\n";
					print OUTPUT7 ">".$id."\n";
					print OUTPUT20 ">".$id."\n";
					
					my $m = '0';
					while (length($read) > $m)
					{
						my $value_ref2b = substr $read, $m, 150;
						$m += 150;
						print OUTPUT7 $value_ref2b."\n";
						print OUTPUT20 $value_ref2b."\n";
					}
					close OUTPUT7;
				}
				
                $compare_haps_stop = "yes";
                $skip_hap = "";
                close $filehandle{$seed_id2};
            }
            elsif ($find_haps_in_seed eq "yes" && $found_haps_in_seed eq "")
            {
                $read = $best_extension;
                $position = length($best_extension);
                $position{$id} = $position;
                $best_extension = "";
                $seed{$id} = $read;
                $skip_hap = "";
                print {$filehandle{$seed_id2}} $position." POSSS\n";
            }
            elsif ($best_extension ne "" && $found_haps_in_seed eq "")
            {
                $read .= $best_extension;                                                 
                $position = length($read);
                $position{$id} = $position;
                $best_extension = "";
                $seed{$id} = $read;
                $skip_hap = "";
				$full_reset_NP = "";
				$full_reset_PB = "";
            }
            
            if ($find_haps_in_seed ne "")
            {
                if ($found_haps_in_seed eq "yes")
                {
                    $find_haps_in_seed = "";
                    $found_haps_in_seed = "";
                    $y = '1';
					$y{$id} = $y;
                    close $filehandle{$seed_id2};
                }
                else
                {
                    $find_haps_in_seed = "yes2";
                    $find_haps_in_seed{$id} = "yes2";
                }
            }
AFTER_EXT2:    
        }
    $y0++;
	$y++;
	$y{$id} = $y;
}

END1:
close OUTPUT18;
close $filehandle{$seed_id2};
close $filehandle3{$seed_id2};
close $filehandle4{$seed_id2};
if ($y < 4 && $first_back_assembly eq "")
{
	#unlink $output_file18;
	unlink $output_file5;
	unlink $output_file13;
	foreach my $seed_id_tmp (keys %seed)
	{
		my $contig_tmp_file = $output_path."contigs_tmp_".$seed_id_tmp."_".$project.".fasta";
		unlink $contig_tmp_file;
	}
	$first_back_assembly = "yes";
}

if (($assembly_length_max eq "WG" || $split_contigs_NP eq "yes" || $split_contigs_PB eq "yes") && $last_seed ne "yes")
{
	$y = '1';
	if ($ploidy > 1)
	{
		$find_haps_in_seed = "yes";
	}
	$seed_input = "";
	$seed_id2 = "";
	$position = '0';
	$position_back = '0';
	$compare_haps = "";
	$compare_haps_stop = "";
	undef %no_hap_track;
	$skip_hap = ""; #CHECKKKKKKKKKKKKKKKKKKK
	$nuc_unique_for_ONT = "";
	undef %last_non_complex_region;	
	undef %filehandle;
	undef %filehandle3;
	undef %filehandle4;
	undef %seed;
	undef %position;
	undef %position_back;	
	undef %split_positions;
	undef %split_positions_extra;
	undef %split_positions_DUP;
	undef %split_positions_back;
	undef %quality_scores;
	undef %quality_scores_gap;
	undef %cut_repeat_seq;
	undef %hap_compare_pos;
	$hap_compare_mismatch_extend = '0';
	undef %PB_split_nucs;
	undef %PB_split_ids;
	$PB_extension = "";
	undef @seed_list_sorted;
	$full_reset_NP = "";	
	undef %save_alignment_data_NP;
	undef %save_alignment_data_PB;
	undef %trace_back_split_NP;
	undef %seeds_list;
	undef %seeds_list_sorted;
	undef %prev_position_hap_compare;
	undef %exclude_reads_hap1_NP;
	undef %exclude_reads_hap2_NP;
	undef %exclude_reads_hap1_PB;
	undef %exclude_reads_hap2_PB;
	undef %exclude_reads_hap1_NP_back;
	undef %exclude_reads_hap2_NP_back;
	undef %exclude_reads_hap1_PB_back;
	undef %exclude_reads_hap2_PB_back;
	undef %trace_back_split_NP;

	if ($first_back_assembly eq "")
	{
		substr $read, 0, $original_seed_length{$id}, "";
		my $read_tmp = reverse($read);
		$read_tmp =~ tr/ACTG/TGAC/;
		$seed_input = $read_tmp;
		$seeds_list{$id."b"} = $seed_input;
		$first_back_assembly{$id."b"} = "yes";
		$seeds_list_sorted{$keep_track_of_reads_number} = $id."b";
		$y{$id."b"} = $y;
		$read = "";
		$id = "";
		goto FIRST_SEED;
	}
	elsif (exists($first_back_assembly{$id}))
	{
		if (exists($split_contigs_reads{$id}))
		{
			foreach my $length_match (keys %{$split_contigs_reads{$id}})
			{
				my $new_contig_seed = substr $read, -$length_match, $length_match;
				my $count_tmp = '1';
				$split_contigs_ends{$id} = $new_contig_seed;
				foreach my $nuc_tmp (keys %{$split_contigs_reads{$id}{$length_match}})
				{
					$seeds_list{$id."c".$count_tmp} = $new_contig_seed;
					$keep_track_of_reads_number++;
					$seeds_list_sorted{$keep_track_of_reads_number."c".$count_tmp} = $id."c".$count_tmp;
					$y{$id."c".$count_tmp} = $y;
					foreach my $id_read (keys %{$split_contigs_reads{$id}{$length_match}{$nuc_tmp}})
					{
						$split_contigs_reads2{$id."c".$count_tmp}{$id_read} = $split_contigs_reads{$id}{$length_match}{$nuc_tmp}{$id_read};
					}
					delete $split_contigs_reads{$id}{$length_match}{$nuc_tmp};
					$count_tmp++;
				}
			}		
			$read = "";
			$id = "";
			goto FIRST_SEED;
		}
	}
	elsif (exists($split_contigs_reads{$id}))
	{
		foreach my $length_match (keys %{$split_contigs_reads{$id}})
		{
			my $new_contig_seed = substr $read, -$length_match, $length_match;
			
			if (keys %split_contigs_ends > 0)
			{
				my $output_file29  = $TMP_directory."Sequence_new_contig_".$id.".fasta";
				open(OUTPUT29, ">".$output_file29) or die "Can't open file $output_file29, $!\n";
				OUTPUT29->autoflush(1);
				foreach my $id_contigs (keys %split_contigs_ends)
				{
					my $lenth_alignment = $length_match;
					if (length($split_contigs_ends{$id_contigs}) < $length_match)
					{
						$lenth_alignment = $split_contigs_ends{$id_contigs};
					}
					my $end_contig_seq = substr $split_contigs_ends{$id_contigs}, -$lenth_alignment;
					print OUTPUT29 ">".$id_contigs."\n";
					print OUTPUT29 $end_contig_seq."\n";	
					close OUTPUT29;
				}
	
				if (-s $output_file29)
				{			
					my $file_tmp = $TMP_directory."blast_seed_test_".$id.".txt";
					my $command = "blastn -query ".$output_file29." -subject ".$output_file20." -out ".$file_tmp." -outfmt 7 -qcov_hsp_perc 90";
					syscmd($command);
					my $count_tmp = '0';
WG_SEED_ACCEPT:									
					if (-s $file_tmp)
					{
						open(SEED_TEST, $file_tmp) or die "\nCan't open file $file_tmp, $!\n";
						
						my $count_lines_tmp = '1';
						while (my $line_tmp = <SEED_TEST>)
						{
							chomp($line_tmp);
							if ($count_lines_tmp eq '4' && $line_tmp eq "# 0 hits found")
							{
								close SEED_TEST;
								goto WG_SEED_ACCEPT2;
							}
							elsif ($count_lines_tmp eq '5' && $line_tmp eq "# BLAST processed 1 queries")
							{
								close SEED_TEST;
								goto WG_SEED_ACCEPT;
							}
							elsif ($count_lines_tmp > 5)
							{
								my @line_tmp = split /\t/, $line_tmp;
								my $accuracy_tmp = $line_tmp[2];
	
								if ($accuracy_tmp > 95 || ($PB_reads eq "" && $input_reads_DB_folder_PB eq "" && $accuracy_tmp > 70))
								{
									close SEED_TEST;
									print {$filehandle{$seed_id2}} $line_tmp." CONTIG_BREAK_MATCH\n";
									goto WG_SEED_ACCEPT3;
								}  
							}
							$count_lines_tmp++;
						}
						close SEED_TEST;
					}
					elsif ($count_tmp < 100000)
					{
						$count_tmp++;
						goto WG_SEED_ACCEPT;
					}
				}
WG_SEED_ACCEPT2:									
			}
			
			my $count_tmp = '1';
			$split_contigs_ends{$id} = $new_contig_seed;
			foreach my $nuc_tmp (keys %{$split_contigs_reads{$id}{$length_match}})
			{
				$seeds_list{$id."c".$count_tmp} = $new_contig_seed;
				$keep_track_of_reads_number++;
				$seeds_list_sorted{$keep_track_of_reads_number."c".$count_tmp} = $id."c".$count_tmp;
				$contig_connections{$id}{$id."c".$count_tmp} = undef;
				$y{$id."c".$count_tmp} = $y;
				foreach my $id_read (keys %{$split_contigs_reads{$id}{$length_match}{$nuc_tmp}})
				{
					$split_contigs_reads2{$id."c".$count_tmp}{$id_read} = $split_contigs_reads{$id}{$length_match}{$nuc_tmp}{$id_read};
				}
				delete $split_contigs_reads{$id}{$length_match}{$nuc_tmp};
				$count_tmp++;
			}
WG_SEED_ACCEPT3:
		}
		$read = "";
		$id = "";
		goto FIRST_SEED;
	}
	$read = "";
	$id = "";
	goto NEXT_SEED;
}
END_LAST:
close OUTPUT4;
close OUTPUT10;
close OUTPUT11;
close OUTPUT12;
close OUTPUT19;
close OUTPUT20;

sleep(10);
if ($NP_reads ne "" || $PB_reads ne "" || $input_reads_DB_folder_NP ne "" || $input_reads_DB_folder_PB)
{
    $chnl->end;
    MCE::Child->waitall;
}

print "\nThank you for using NOVOLoci!\n\n";
