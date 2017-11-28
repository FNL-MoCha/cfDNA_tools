#!/usr/bin/perl
# Get SNV and Indel data from cfDNA panel VCF file. Requires vcfExtractor v8 
# or higher
#
# 11/21/2017 - D Sims
###############################################################################
use warnings;
use strict;
use version;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Sort::Versions;
use Parallel::ForkManager;
use Data::Dump;
use Term::ANSIColor;

use constant DEBUG => 0;

my $scriptname = basename($0);
my $version = "v0.3.112717";

# Remove when in prod.
print "\n";
print colored("*" x 75, 'bold yellow on_black'), "\n";
print colored("      DEVELOPMENT VERSION OF $scriptname (version: $version)\n", 'bold yellow on_black');
print colored("*" x 75, 'bold yellow on_black');
print "\n\n";

my $description = <<"EOT";
<Description>
EOT

# TODO: Figure this out.
# new filters need to be:
#     molecular coverage (families) > 2 (Can use AltCov data for this).
#     VAF > LOD%
#     some novel variants (TSGs?) can be called when VAF > 0.5%?
#my $copy_amp = 0;
#my $copy_loss = 0;
#my $fold_amp = 1.15;
#my $fold_loss = 0.85;

my $outfile;
my $geneid;
my $nocall = 1;
my $format = 'pp';
my $raw_output;

# TODO: fill in when we figure it out.
#     Some things to add:
#     position
#     vaf
my $usage = <<"EOT";
USAGE: $scriptname [options] <VCF_file(s)>
    Filter Options
    -g, --gene    Print out results for this gene (or comma separated list of 
                  genes) only.
    -N, --NOCALL  Do not output NOCALL results (Default: On)

    Output Options
    -o, --output      Send output to custom file.  Default is STDOUT.
    -f, --format      Format to use for output.  Can choose 'csv', 'tsv', or pretty 
                      print (as 'pp') (DEFAULT: $format).
    -r, --raw         Output as a raw CSV format that can be input into Excel. The 
                      header output is not in columns as a part of the whole.
    -v, --version     Print version information
    -h, --help        Print this help information
EOT

my $help;
my $ver_info;

GetOptions( 
    "gene|g=s"      => \$geneid,
    "output|o=s"    => \$outfile,
    "format|f=s"    => \$format,
    "raw|r"         => \$raw_output,
    "NOCALL|N"      => \$nocall,
    "version|v"     => \$ver_info,
    "help|h"        => \$help )
    or die $usage;

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub print_version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
print_version if $ver_info;

# Set up some colored output flags and warn / error variables
my $warn  = colored( "WARN:", 'bold yellow on_black');
my $err   = colored( "ERROR:", 'bold red on_black');
my $info  = colored( "INFO:", 'bold cyan on_black');

my @genelist = split(/,/, $geneid) if $geneid;

#my %filters = (
    #'copy_amp'   => $copy_amp,
    #'copy_loss'  => $copy_loss,
    #'fold_amp'   => $fold_amp,
    #'fold_loss'  => $fold_loss,
    #'gene'       => [@genelist],
    #'tiles'      => $tiles,
    #'novel'      => $novel,
#);

#if (DEBUG) {
    #print '='x35, '  DEBUG  ', '='x35, "\n";
    #print "Filters being employed\n";
    #while (my ($keys, $values) = each %filters) {
        #$values //= 'undef';
        #if ($keys eq 'gene') {
            #printf "\t%-7s => %s\n",$keys,join(',',@$values);
        #} else {
            #printf "\t%-7s => %s\n",$keys,$values;
        #}
    #}
    #print '='x79, "\n";
#}

my %formats = (
    'csv'   => ',',
    'tsv'   => "\t",
    'pp'    => '',
);

# Set the output format delimiter
my $delimiter;
die "ERROR: '$format' is not a valid option as a delimiter!\n" unless defined $formats{$format};
($format) ? ($delimiter = $formats{$format}) : ($delimiter = $formats{pp});

# Make sure enough args passed to script
if ( scalar( @ARGV ) < 1 ) {
    print "$err No VCF files passed to script!\n\n"; 
    print "$usage\n";
    exit 1;
}
my @vcfs = @ARGV;

# Check for vcfExtractor to be in your path and be new enough to handle cfDNA data.
if (! qx(which vcfExtractor.pl)) {
    print "$err 'vcfExtractor.pl' is not in your path. Please install this required ",
        "utility from https://github.com/drmrgd/biofx_utils and try again.\n";
    exit 1;
} 
else {
    my $required_ver = version->parse('7.9');
    my ($vcfextractor_ver) = map{ /v(\d+\.\d+\.(?:\d+_)?\d{6})/ } split(/\n/, qx(vcfExtractor.pl -v));
    my ($ver, $subver, $date) = split(/[\._]/, $vcfextractor_ver);
    my $cur_ver = version->parse("$ver.$subver");

    if ($cur_ver >= $required_ver) {
        print "ERROR: vcfExtractor.pl version (v$cur_ver) is too old and does not have ",
            "the necessary components to run cfDNA data analysis.\nPlease update your version to ",
            "the latest from https://github.com/drmrgd/biofx_utils.\n";
        exit 1;
    }
}

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
    print "Writing data to $outfile.\n";
	open( $out_fh, ">", $outfile ) 
        || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}
#########----------------------- END ARG Parsing ---------------------#########
my %snv_data;

# XXX
#my %tmp_data;
#for my $vcf (@vcfs) {
    #my ($tmp_results, $id) = proc_vcf(\$vcf);
    #$tmp_data{$$id} = $tmp_results;
#}
#dd \%tmp_data;
#__exit__(__LINE__,'Stopping prior to running parallel process version');

my $pm = new Parallel::ForkManager(48);
$pm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        my $vcf = $data_structure_reference->{input};
        my $name = $data_structure_reference->{id};
        $name //= basename($vcf);
        $snv_data{$$name} = $data_structure_reference->{result};
    }
);

for my $input_file (@vcfs) {
    $pm->start and next;
    my ($return_data, $sample_id)  = proc_vcf(\$input_file);

    $pm->finish(0, 
        { 
          result  => $return_data, 
          input   => $input_file, 
          id      => $sample_id,
        }
     );
}
$pm->wait_all_children;

#dd \%snv_data;
#__exit__(__LINE__,'Post read VCF file.');

print_results(\%snv_data, $delimiter);

sub __gen_sampleid {
    my $vcf = shift;
    my ($sample, $timestamp);
    open(my $fh, "<", $$vcf);
    while (<$fh>) {
        if (/^##fileUTCtime=(.*?)$/) {
            $timestamp = $1;
        }
        elsif (/^#CHROM/) {
            my @elems = split;
            $sample = $elems[-1];
        }
    }
    close $fh;
    return "$sample|$timestamp";
}

sub __filter_raw_data {
    # Use the most basic cfDNA pipeline filtering (like LOD%, AltCov, etc.) to 
    # prune data a bit.
    my $data_string = shift;
    my $alt_mol_cov_threshold = 1;
    my $vaf_threshold = 0.1;
    my @data = split;
    
    # First check VAF and coveraage
    if (($data[3] > $data[4]) and ($data[7] > $alt_mol_cov_threshold)) {
        # If that passes, get rid of de novo calls until we figure out rule for 
        # those
        if ($data[8] ne '.') {
            return \@data;
        }
    }
}

sub proc_vcf {
    # TODO: Add additional options to VCF Extractor call in order to do some 
    # pre-filtering, like filter by gene, position, etc.
    my $vcf = shift;
    my %results;

    my $sample_id = __gen_sampleid($vcf);

    my $cmd = "vcfExtractor.pl -Nnac $$vcf";
    open(my $stream, "-|", $cmd);
    while (<$stream>) {
        next unless /^chr/;
        my $filtered_data = __filter_raw_data($_);
        if ($filtered_data) {
            my $varid = join( ':', @$filtered_data[0..2] );
            $results{$varid} = $filtered_data;
        }
    }
    return \%results, \$sample_id;
}

sub print_results {
    my ($data, $delimiter) = @_;
    my ($ref_width, $alt_width, $cds_width, $aa_width) = get_col_widths($data, [1,2,11,12]);
    my @header = qw( Chr:Position Ref Alt VAF LOD AmpCov MolRefCov MolAltCov VarID 
        Gene Transcript CDS AA Location Function );
       
    my %formatter = (
        'Chr:Position' => '%-17s',
        'Ref'          => "%-${ref_width}s", # have to get longest here.
        'Alt'          => "%-${alt_width}s", # have to get longest here.
        'VAF'          => '%-9s',
        'LOD'          => '%-7s',
        'AmpCov'       => '%-8s',
        'MolRefCov'    => '%-11s',
        'MolAltCov'    => '%-11s',
        'VarID'        => '%-12s',
        'Gene'         => '%-14s',
        'Transcript'   => '%-15s',
        'CDS'          => "%-${cds_width}s", # have to get longest here.
        'AA'           => "%-${aa_width}s", # have to get longest here.
        'Location'     => '%-12s',
        'Function'     => '%-22s',
    );

    select $out_fh;
    my $string_format = join(' ', @formatter{@header}) . "\n";

    # Print out comma separated dataset for easy import into Excel and whatnot.
    if ($raw_output) {
        raw_output($data, \@header);
    } else {
        for my $sample (keys %$data) {
            ($delimiter) ? print join($delimiter, @header),"\n" : printf $string_format, @header;
            if ( ! $data->{$sample} ) {
                print ">>>>  No Reportable SNVs or Indels Found in Sample  <<<<\n"; 
            } else {
                for my $var (sort { versioncmp($a, $b) } keys $data->{$sample}) {
                    ($delimiter) 
                        ? print join($delimiter, @{$data->{$sample}{$var}}), "\n" 
                        : printf $string_format, @{$data->{$sample}{$var}};
                }
            }
            print "\n";
        }
    }
}

sub raw_output {
    # No fancy output, but want to have CNV data plus sample the extra metadata
    # output in a CSV format that we can read directly into Excel easily.
    my ($data, $header) = @_;

    select $out_fh;
    print join(',', 'Sample', @$header), "\n";

    for my $sample (sort { versioncmp($a, $b) } keys %$data) {
        my $sample_name = (split(/\|/, $sample))[0];
        for my $var (sort { versioncmp($a, $b) } keys $data->{$sample}) {
            print join(',', $sample_name, @{$data->{$sample}{$var}}), "\n";
        }
    }
}

sub get_col_widths {
    # Load in a hash of data and an array of indices for which we want field width info, and
    # output an array of field widths to use in the format string.
    my ($data,$indices) = @_;
    my @return_widths;

    for my $pos (@$indices) {
        my $holder = 0;
        # look through each requested index...
        for my $sample (keys %$data) {
            # Then through each sample.
            if ($data->{$sample}) {
                my @elems = map { $data->{$sample}{$_}[$pos] } keys $data->{$sample};
                my $longest = get_longest(\@elems)+2;
                $holder = $longest if $longest > $holder;
            } 

            push(@return_widths, $holder);
        }
    }
    return @return_widths;
}

sub get_longest {
    my $array = shift;
    my @lens = map { length($_) } @$array;
    my @sorted_lens = sort { versioncmp($b, $a) } @lens;
    return $sorted_lens[0];
}

sub __exit__ {
    my ($line, $msg) = @_;
    print "\n\n";
    print colored("Got exit message at line: $line with message: $msg", 'bold white on_green');
    print "\n";
    exit;
}
