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
my $version = "v0.2.112117";

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
    print "ERROR: No VCF files passed to script!\n\n"; 
    print "$usage\n";
    exit 1;
}
my @vcfs = @ARGV;

# Check for required vcfExtractor to be in your path.
# TODO: Colorize error msg.
if (! qx(which vcfExtractor.pl)) {
    print "ERROR: vcfExtractor.pl is not in your path. Please install this required ",
        "utility from https://github.com/drmrgd/biofx_utils and try again.\n";
    exit 1;
} 
else {
    my $required_ver = version->parse('7.9');
    my ($vcfextractor_ver) = map{ /v(\d\.\d\.(?:\d_)?\d{6})/ } split(/\n/, qx(vcfExtractor.pl -v));
    my ($ver, $subver, $date) = split(/\./, $vcfextractor_ver);
    my $cur_ver = version->parse("$ver.$subver");

    if ($cur_ver <=> $required_ver) {
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

dd \%snv_data;
__exit__(__LINE__,'Post read VCF file.');

=cut
my $results = proc_cnv_data(\%cnv_data, \%filters);
#dd $results;
#__exit__(__LINE__,'Post result processing.');

print_results($results, $delimiter);

=cut

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

=cut
sub proc_cnv_data {
    my ($cnv_data, $filters) = @_;
    my %results;

    for my $sample ( keys %$cnv_data ) {
        $results{$sample} = [];
        my @outfields = qw( END LEN NUMTILES CN FD HS FUNC PVAL RMMDP MMDP );

        for my $cnv ( sort { versioncmp ( $a, $b ) } keys %{$$cnv_data{$sample}} ) {
            my %mapped_cnv_data;
            last if $cnv eq 'NONE';

            my ($chr, $start, $gene, undef) = split( /:/, $cnv );
            %mapped_cnv_data = map{ $_ => $cnv_data{$sample}->{$cnv}->{$_} } @outfields;
            @mapped_cnv_data{qw(chr start gene undef)} = split(/:/, $cnv);
            $mapped_cnv_data{HS} //= 'No';

            # Get OVAT Annot Data
            # XXX: Remove this for now since there is no OVAT annotation at the moment.  Will bring back 
            # if they do re-implement this.
            #my ($gene_class, $variant_class);
            #my $func = $mapped_cnv_data{FUNC};
            #if ( $func && $func =~ /oncomine/ ) {
                #my $json_annot = JSON->new->allow_singlequote->decode($func);
                #my $parsed_annot = $$json_annot[0];
                #$gene_class = $$parsed_annot{'oncomineGeneClass'};
                #$variant_class = $$parsed_annot{'oncomineVariantClass'};
            #} else {
                #$gene_class = $variant_class = '---';
            #}
            #$mapped_cnv_data{GC} = $gene_class;
            #$mapped_cnv_data{VC} = $variant_class;

            my @filtered_data = filter_results(\%mapped_cnv_data, $filters);
            push(@{$results{$sample}}, \@filtered_data) if @filtered_data;
        }
    }
    return \%results;
}

sub print_results {
    my ($data, $delimiter) = @_;
    my @header = qw( Chr Gene Start End Length Tiles CN FD p-val Med_Mol_Cov 
        Med_Read_Cov );
    #my $pp_format = "%-8s %-8s %-11s %-11s %-11s %-8s %-8s %-8s %-8s %-8s %-8s %-18s\n";

    my %formatter = (
        'Chr'          => '%-8s',
        'Gene'         => '%-8s',
        'Start'        => '%-11s', 
        'End'          => '%-11s',
        'Length'       => '%-11s',
        'Tiles'        => '%-8s',
        'CN'           => '%-8s',
        'FD'           => '%-8s',
        'p-val'        => '%-8s',
        'Med_Mol_Cov'  => '%-14s',
        'Med_Read_Cov' => '%-14s',
    );

    select $out_fh;
    my $string_format = join(' ', @formatter{@header}) . "\n";

    # Print out comma separated dataset for easy import into Excel and whatnot.
    if ($raw_output) {
        raw_output($data, \@header);
    } else {
        for my $sample (keys %$data) {
            my ($id, $gender, $mapd, $cellularity) = split( /:/, $sample );
            print "::: CNV Data For $id (Gender: $gender, Cellularity: $cellularity, MAPD: $mapd) :::\n";
            ($delimiter) ? print join($delimiter, @header),"\n" : printf $string_format, @header;
            if ( ! @{$$data{$sample}} ) {
                print ">>>>  No Reportable CNVs Found in Sample  <<<<\n"; 
            } else {
                for my $cnv (@{$$data{$sample}}) {
                    ($delimiter) ? print join($delimiter, @$cnv), "\n" : printf $string_format, @$cnv;
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
    print join(',', qw(Sample Gender MAPD Cellularity), @$header), "\n";

    for my $sample (sort { versioncmp($a, $b) } keys %$data) {
        my @elems = split(/:/, $sample);
        for my $cnv (@{$$data{$sample}}) {
            print join(',', @elems, @$cnv), "\n";
        }
    }
}

sub filter_results {
    # Filter out CNV data prior to printing it all out. Each call receives
    # a hash of data, and either None returned or a hash of specific fields
    # as defined in "return_data()"
    my ($data, $filters) = @_;
    my @cn_thresholds = @$filters{qw(copy_amp copy_loss fold_amp fold_loss)};

    # Gene level filter
    return if (@{$filters{gene}}) and ! grep {$$data{gene} eq $_} @{$filters{gene}};

    # We made it the whole way through; check for copy number thresholds
    (copy_number_filter($data, \@cn_thresholds)) ? return return_data($data) : return;
}

sub return_data {
    my $data = shift;
    my @fields = qw(chr gene start END LEN NUMTILES CN FD PVAL RMMDP MMDP);
    return @$data{@fields};
}

sub copy_number_filter {
    # Filter data based on either Fold Diff value or Copy Number value.
    my ($data, $threshold) = @_;
    my ($ca, $cl, $fa, $fl) = @$threshold;

    # XXX: Copy data is dependent on p-val. For now, let's just filter out anything
    # with a p-val > 10e-5 per TF QC doc.  Can tailor this a bit later.
    return 0 if $$data{'PVAL'} > 0.00005;

    # If we're using fold diff values, return 1 if FD > amp threshold or less 
    # than loss threshold.
    if ($fa and $fl) {
        return 1 if ($$data{'FD'} > $fa || $$data{'FD'} < $fl);
    }
    elsif ($ca and $cl) {
        return 1 if ($$data{'CN'} > $ca || $$data{'CN'} < $cl);
    }
    else {
        # Return everything if there are no filters.
        return 1;
    }

    # If we got here, the data does not meet threshold requirements and will be
    # filtered out.
    return 0;
}
=cut
sub __exit__ {
    my ($line, $msg) = @_;
    print "\n\n";
    print colored("Got exit message at line: $line with message: $msg", 'bold white on_green');
    print "\n";
    exit;
}
